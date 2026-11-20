use std::fs::{self, File};
use std::io;
use std::io::Write;
use std::io::{BufReader, BufWriter};

use anyhow::{Result, bail};
use codata as cc;
use fitsio::FitsFile;
use fitsio::hdu::FitsHdu;
use flate2::read::GzDecoder;
use ndarray::ArrayD;
use ndarray::Axis;
use ndarray::s;
use planetes_ext::types::{Matrix, RMatrix, RTensor, RTensor4, RVector, UVector};
use spectre::lamda::LAMDAData;

use crate::flits::config::Cell;

use super::config::{FLITSConfig, Grid, Molecule, RadTransition, SimState};

pub type CellMatrix = Matrix<Cell>;

#[derive(Debug, Clone, Copy)]
#[repr(i32)]
pub enum CellParameter {
    GasTemperature = 1,
    DustTemperature = 2,
    GasDensity = 3,
    Wavelength = 4,
    StellarSpectrum = 5,
    AbsOpacities = 7,
    ExtOpacities = 8,
    InternalField = 9,
    SourceFunction = 10,
    MolecularParticleDensities = 11,
    LineBroadening = 12,
    RelativeLevelPopulations = 13,
}

impl CellParameter {
    fn read_next(
        self,
        c: &mut CellMatrix,
        hdu: &FitsHdu,
        fptr: &mut FitsFile,
        grid: &mut Grid,
        flitsc: &mut FLITSConfig,
    ) -> Result<()> {
        match self {
            CellParameter::GasTemperature
            | CellParameter::DustTemperature
            | CellParameter::GasDensity
            | CellParameter::Wavelength
            | CellParameter::StellarSpectrum
            | CellParameter::AbsOpacities
            | CellParameter::ExtOpacities
            | CellParameter::InternalField
            | CellParameter::SourceFunction
            | CellParameter::MolecularParticleDensities
            | CellParameter::LineBroadening => {
                Self::read_hdu_generic(self, hdu, fptr, c, grid, flitsc)?;
            }

            CellParameter::RelativeLevelPopulations => {
                let nspec = grid.nspec;
                flitsc.npop0 = UVector::zeros(nspec);
                flitsc.mol_name0 = Vec::new();
                for imol in 0..nspec {
                    match fptr.iter().nth(13 + imol) {
                        Some(hdu) => {
                            Self::read_hdu_level_populations(imol, &hdu, fptr, c, grid, flitsc)?;
                        }
                        None => return Ok(()), // No more HDUs
                    }
                }
            }
        }

        Ok(())
    }

    #[allow(clippy::assigning_clones)]
    fn read_hdu_generic(
        kind: CellParameter,
        hdu: &FitsHdu,
        fptr: &mut FitsFile,
        c: &mut CellMatrix,
        grid: &mut Grid,
        flitsc: &mut FLITSConfig,
    ) -> Result<()> {
        match kind {
            CellParameter::GasTemperature => {
                let naxis = 2;
                let data = extract_data(fptr, hdu, naxis)?;
                for i in 1..grid.nr {
                    for j in 2..=grid.ntheta {
                        let val = data[[i - 1, grid.ntheta - j, 0, 0]];
                        let clamped_temp = if val < 1.0 { 1.0 } else { val };
                        c[[i, j]].tgas = clamped_temp;
                    }
                    c[[i, 0]].tgas = c[[i, 1]].tgas;
                }
            }
            CellParameter::DustTemperature => {
                let naxis = 2;
                let data = extract_data(fptr, hdu, naxis)?;
                for i in 1..grid.nr {
                    for j in 2..=grid.ntheta {
                        let val = data[[i - 1, grid.ntheta - j, 0, 0]];
                        let clamped_temp = if val < 1.0 { 1.0 } else { val };
                        c[[i, j]].tdust = clamped_temp;
                    }
                    c[[i, 0]].tdust = c[[i, 1]].tdust;
                }
            }
            CellParameter::GasDensity => {
                let naxis = 2;
                let data = extract_data(fptr, hdu, naxis)?;
                for i in 1..grid.nr {
                    for j in 2..=grid.ntheta {
                        let val = data[[i - 1, grid.ntheta - j, 0, 0]];
                        let clamped_density = if val < 1e-50 { 1e-60 } else { val };
                        c[[i, j]].dens = clamped_density;
                    }
                    c[[i, 0]].dens = c[[i, 1]].dens;
                }
            }
            CellParameter::Wavelength => {
                let naxis = 1;
                let data = extract_data(fptr, hdu, naxis)?;
                grid.nlam = data.shape()[0];
                for i in 0..=grid.nr {
                    for j in 0..=grid.ntheta {
                        c[[i, j]].kabs = RVector::zeros(grid.nlam);
                        c[[i, j]].albedo = RVector::zeros(grid.nlam);
                        c[[i, j]].kext = RVector::zeros(grid.nlam);
                        c[[i, j]].lrf = RVector::zeros(grid.nlam);
                        c[[i, j]].s = RVector::zeros(grid.nlam);
                    }
                }

                let lam_cont: RVector = &data.slice(s![..grid.nlam, 0, 0, 0]) * 1e4;
                if lam_cont[0] > flitsc.lmin {
                    if lam_cont[0] > flitsc.lmax {
                        bail!("Wavelength interval not in tracer file")
                    }
                    flitsc.lmin = lam_cont[0];
                    println!("Adjusting wavelength to {}", flitsc.lmin);
                }
                if lam_cont[grid.nlam - 1] < flitsc.lmax {
                    if lam_cont[grid.nlam - 1] < flitsc.lmin {
                        bail!("Wavelength interval not in tracer file")
                    }
                    flitsc.lmax = lam_cont[grid.nlam - 1];
                    println!("Adjusting wavelength to {}", flitsc.lmax);
                }
                grid.lam_cont = lam_cont;
            }
            CellParameter::StellarSpectrum => {
                let naxis = 1;
                let data = extract_data(fptr, hdu, naxis)?;
                let fstar = data.slice(s![..grid.nlam, 0, 0, 0]).to_owned();
                grid.fstar = fstar;
            }
            CellParameter::AbsOpacities => {
                let naxis = 4;
                let data = extract_data(fptr, hdu, naxis)?;
                for i in 1..grid.nr {
                    for j in 2..=grid.ntheta {
                        for k in 0..grid.nlam {
                            c[[i, j]].kabs[k] = data[[k, i - 1, grid.ntheta - j, 0]];
                        }
                    }
                    c[[i, 1]].kabs = &c[[i, 2]].kabs / 1e20;
                }
            }
            CellParameter::ExtOpacities => {
                let naxis = 4;
                let data = extract_data(fptr, hdu, naxis)?;
                for i in 1..grid.nr {
                    for j in 2..=grid.ntheta {
                        for k in 0..grid.nlam {
                            c[[i, j]].kext[k] = data[[k, i - 1, grid.ntheta - j, 0]];
                            c[[i, j]].albedo[k] = if c[[i, j]].kext[k] > 1e-150 {
                                (c[[i, j]].kext[k] - c[[i, j]].kabs[k]) / c[[i, j]].kext[k]
                            } else {
                                0.5
                            };
                        }
                    }
                    c[[i, 1]].kext = &c[[i, 2]].kext / 1e20;
                    c[[i, 1]].albedo = c[[i, 2]].albedo.clone();
                }
            }
            CellParameter::InternalField => {
                let naxis = 3;
                let data = extract_data(fptr, hdu, naxis)?;
                for i in 1..grid.nr {
                    for j in 2..=grid.ntheta {
                        for k in 0..grid.nlam {
                            c[[i, j]].lrf[k] = data[[k, i - 1, grid.ntheta - j, 0]];
                        }
                    }
                    c[[i, 1]].lrf = &c[[i, 2]].lrf / 1e20;
                }
            }
            CellParameter::SourceFunction => {
                let naxis = 3;
                let data = extract_data(fptr, hdu, naxis)?;
                for i in 1..grid.nr {
                    for j in 2..=grid.ntheta {
                        for k in 0..grid.nlam {
                            c[[i, j]].s[k] = data[[k, i - 1, grid.ntheta - j, 0]];
                        }
                    }
                    c[[i, 1]].s = &c[[i, 2]].s / 1e20;
                }
            }
            CellParameter::MolecularParticleDensities => {
                let naxis = 3;
                let data = extract_data(fptr, hdu, naxis)?;
                for i in 1..grid.nr {
                    for j in 2..=grid.ntheta {
                        for imol in 0..grid.nspec {
                            c[[i, j]].n0[imol] = data[[imol, i - 1, grid.ntheta - j, 0]];
                        }
                    }
                    c[[i, 1]].n0 = &c[[i, 2]].n0 / 1e20;
                }
            }
            CellParameter::LineBroadening => {
                let nmol = flitsc.line_files.len();
                let naxis = 3;
                let data = extract_data(fptr, hdu, naxis)?;
                for i in 1..grid.nr {
                    for j in 2..=grid.ntheta {
                        for imol in 0..nmol {
                            c[[i, j]].line_width0[imol] = data[[imol, i - 1, grid.ntheta - j, 0]];
                        }
                    }
                    c[[i, 1]].line_width0 = &c[[i, 2]].line_width0 / 1e20;
                }
            }
            CellParameter::RelativeLevelPopulations => {
                // This is processed separately
            }
        }
        Ok(())
    }

    #[allow(clippy::assigning_clones)]
    fn read_hdu_level_populations(
        imol: usize,
        hdu: &FitsHdu,
        fptr: &mut FitsFile,
        c: &mut CellMatrix,
        grid: &Grid,
        flitsc: &mut FLITSConfig,
    ) -> Result<()> {
        let nlevels = hdu.read_key::<f64>(fptr, "NLEV")? as usize;
        flitsc
            .mol_name0
            .push(hdu.read_key::<String>(fptr, "SPECIES")?);

        let data: ArrayD<f64> = hdu.read_image(fptr)?;
        let data: RTensor = data.into_dimensionality()?;
        let data = data.permuted_axes([2, 1, 0]).insert_axis(Axis(3));

        flitsc.npop0[imol] = nlevels;

        for i in 0..=grid.nr {
            for j in 0..=grid.ntheta {
                c[[i, j]].npop0[imol].n = RVector::zeros(nlevels);
            }
        }

        for i in 1..(grid.nr - 1) {
            for j in 2..grid.ntheta {
                for k in 0..nlevels {
                    c[[i, j]].npop0[imol].n[k] = data[[k, i - 1, grid.ntheta - j, 0]];
                }
            }
            c[[i, 1]].npop0[imol].n = c[[i, 2]].npop0[imol].n.clone();
        }

        Ok(())
    }
}

const CELL_PARAMS_ORDER: &[CellParameter] = &[
    CellParameter::GasTemperature,
    CellParameter::DustTemperature,
    CellParameter::GasDensity,
    CellParameter::Wavelength,
    CellParameter::StellarSpectrum,
    // Skipping the image HDU after StellarSpectrum
    CellParameter::AbsOpacities,
    CellParameter::ExtOpacities,
    CellParameter::InternalField,
    CellParameter::SourceFunction,
    CellParameter::MolecularParticleDensities,
    CellParameter::LineBroadening,
    CellParameter::RelativeLevelPopulations,
];

/// Reads the FLITS file, uncompressing if necessary, and extracts the grid and cell parameters into a [`SimState`] struct.
///
/// # Errors
/// If the FLITS file cannot be found, uncompressed, or read properly, an error is returned.
pub fn read_fits(flitsc: &mut FLITSConfig) -> Result<SimState> {
    flitsc.flits_file = flitsc
        .flits_file_compressed
        .strip_suffix(".gz")
        .map(ToString::to_string)
        .ok_or(anyhow::anyhow!("FLITS file name must end with .gz"))?;

    let fits_path = flitsc.config_dir.join(&flitsc.flits_file);

    if !fits_path.exists() {
        uncompress_fits_file(flitsc)?;
    }

    let mut fptr = FitsFile::open(&fits_path)?;

    let mut grid = read_grid(&mut fptr, flitsc)?;
    let c = read_cell_parameters(&mut fptr, flitsc, &mut grid)?;

    let sim_state = SimState {
        c,
        grid,
        ..Default::default()
    };

    Ok(sim_state)
}

fn extract_data(fptr: &mut FitsFile, hdu: &FitsHdu, naxis: usize) -> Result<RTensor4> {
    let data: ArrayD<f64> = hdu.read_image(fptr)?;
    match naxis {
        1 => {
            let data1: RVector = data.into_dimensionality()?;
            Ok(data1
                .insert_axis(Axis(1))
                .insert_axis(Axis(2))
                .insert_axis(Axis(3)))
        }
        2 => {
            let data2: RMatrix = data.into_dimensionality()?;
            let data_permuted = data2.permuted_axes([1, 0]);
            Ok(data_permuted.insert_axis(Axis(2)).insert_axis(Axis(3)))
        }
        3 | 4 => {
            let data3: RTensor = data.into_dimensionality()?;
            let data_permuted = data3.permuted_axes([2, 1, 0]);
            Ok(data_permuted.insert_axis(Axis(3)))
        }
        _ => Err(anyhow::anyhow!("Unsupported number of axes: {naxis}")),
    }
}

fn uncompress_fits_file(flitsc: &mut FLITSConfig) -> Result<()> {
    let gz_path = flitsc.config_dir.join(&flitsc.flits_file_compressed);
    let gz = File::open(&gz_path)?;
    let buffered = BufReader::new(gz);
    let mut decoder = GzDecoder::new(buffered);

    // Create output path in same directory as config_dir
    let output_path = flitsc.config_dir.join(&flitsc.flits_file);

    let output_file = File::create(&output_path)?;
    let mut writer = BufWriter::new(output_file);

    io::copy(&mut decoder, &mut writer)?;

    Ok(())
}

fn read_grid(fptr: &mut FitsFile, flitsc: &mut FLITSConfig) -> Result<Grid> {
    let hdu = fptr.primary_hdu()?;

    let mut grid = Grid::from_hdu(fptr, &hdu)?;

    flitsc.m_star = hdu.read_key::<f64>(fptr, "Mstar")?;
    flitsc.r_star = hdu.read_key::<f64>(fptr, "Rstar")?;
    flitsc.distance = hdu.read_key::<f64>(fptr, "distance")?;

    // Axis(0) --> r
    // Axis(1) --> theta
    // Axis(2) --> ???
    let data: ArrayD<f64> = hdu.read_image(fptr)?;
    let data3: RTensor = data.into_dimensionality()?;
    let array = data3.permuted_axes([2, 1, 0]).insert_axis(Axis(3));

    grid.initialize(&array, flitsc);

    Ok(grid)
}

fn read_cell_parameters(
    fptr: &mut FitsFile,
    flitsc: &mut FLITSConfig,
    grid: &mut Grid,
) -> Result<CellMatrix> {
    let mut c = CellMatrix::default(((grid.nr + 1), (grid.ntheta + 1)));

    for i in 0..=grid.nr {
        for j in 0..=grid.ntheta {
            c[[i, j]].allocate(grid.nspec);
        }
    }

    // Process rest of HDUs
    for param in CELL_PARAMS_ORDER {
        let hdu = fptr
            .iter()
            .nth(*param as usize)
            .ok_or_else(|| anyhow::anyhow!("No HDUs found in FITS file"))?;
        param.read_next(&mut c, &hdu, fptr, grid, flitsc)?;
    }

    set_boundary_conditions(&mut c, grid, flitsc);

    Ok(c)
}

fn set_boundary_conditions(c: &mut CellMatrix, grid: &Grid, flitsc: &FLITSConfig) {
    for i in 1..grid.nr {
        for j in 1..=grid.ntheta {
            let xx = grid.r_av[i];
            let zz = xx / grid.theta_av[[i, j]].tan();
            let rr = (xx * xx + zz * zz).sqrt();
            c[[i, j]].v = (cc::NEWTONIAN_CONSTANT_OF_GRAVITATION_CGS
                * flitsc.m_star
                * cc::astro::M_SUN_CGS
                * xx
                * xx
                / (rr * rr * rr))
                .sqrt();
        }
    }

    for j in 0..=grid.ntheta {
        c[[0, j]].kext.fill(1e-70);
        c[[grid.nr, j]].kext.fill(1e-70);
        c[[0, j]].kabs.fill(1e-70);
        c[[grid.nr, j]].kabs.fill(1e-70);
        let xx_0 = grid.r_av[0];
        let zz_0 = xx_0 / grid.theta_av[[0, j]].tan();
        let rr_0 = (xx_0 * xx_0 + zz_0 * zz_0).sqrt();
        c[[0, j]].v = (cc::NEWTONIAN_CONSTANT_OF_GRAVITATION_CGS
            * flitsc.m_star
            * cc::astro::M_SUN_CGS
            * xx_0
            * xx_0
            / (rr_0 * rr_0 * rr_0))
            .sqrt();
        c[[0, j]].dens = 1e-60;

        let xx_nr = grid.r_av[grid.nr];
        let zz_nr = xx_nr / grid.theta_av[[grid.nr, j]].tan();
        let rr_nr = (xx_nr * xx_nr + zz_nr * zz_nr).sqrt();
        c[[grid.nr, j]].v = (cc::NEWTONIAN_CONSTANT_OF_GRAVITATION_CGS
            * flitsc.m_star
            * cc::astro::M_SUN_CGS
            * xx_nr
            * xx_nr
            / (rr_nr * rr_nr * rr_nr))
            .sqrt();
        c[[grid.nr, j]].dens = 1e-60;
    }

    for i in 0..=grid.nr {
        c[[i, 0]].kext.fill(1e-70);
        c[[i, 0]].kabs.fill(1e-70);
        let xx = grid.r_av[i];
        let zz = xx / grid.theta_av[[i, 0]].tan();
        let rr = (xx * xx + zz * zz).sqrt();
        c[[i, 0]].v = (cc::NEWTONIAN_CONSTANT_OF_GRAVITATION_CGS
            * flitsc.m_star
            * cc::astro::M_SUN_CGS
            * xx
            * xx
            / (rr * rr * rr))
            .sqrt();
        c[[i, 0]].dens = 1e-60;
    }
}

/// Reads the LAMDA files specified in the configuration and populates the `Molecule` structs in the `FLITSConfig`.
///
/// # Errors
/// If any of the LAMDA files cannot be found or read properly, an error is returned.
pub fn read_lamda_files(flitsc: &FLITSConfig) -> Result<Vec<Molecule>> {
    let mut mol = Vec::with_capacity(flitsc.line_files.len());
    let two_h_over_c2 = 2.0 * cc::PLANCK_CONSTANT_CGS / cc::SPEED_OF_LIGHT_IN_VACUUM_CGS.powi(2);

    for line_file in &flitsc.line_files {
        let data = LAMDAData::from_path(flitsc.config_dir.join(line_file))?;
        println!("Reading line file for {}", data.name);
        let mut molecule = Molecule {
            name: data.name,
            weight: data.weight,
            LTE: flitsc.LTE,
            ..Default::default()
        };

        molecule.levels = data
            .levels
            .into_iter()
            .map(|mut level| {
                level.energy *= 11600.0 / 8065.54; // Convert energy from cm^-1 to K
                level
            })
            .collect();

        molecule.radset = data
            .radset
            .into_iter()
            .map(|rad| {
                let ilow = rad.low - 1; // zero-based index
                let iup = rad.up - 1;
                let freq = rad.freq * 1e9; // GHz to Hz
                let bul = rad.einst_a / (two_h_over_c2 * freq.powi(3));
                let blu = bul * molecule.levels[iup].weight / molecule.levels[ilow].weight;

                RadTransition {
                    id: rad.id,
                    up: rad.up,
                    low: rad.low,
                    einst_a: rad.einst_a,
                    einst_b_ul: bul,
                    einst_b_lu: blu,
                    freq,
                    lam: cc::SPEED_OF_LIGHT_IN_VACUUM_CGS * 1e4 / freq,
                    energy: rad.energy,
                }
            })
            .collect();

        mol.push(molecule);
    }
    Ok(mol)
}

/// Writes the image grid to a file in the output directory, scaling the values by the astronomical unit in CGS units.
///
/// # Errors
/// If the output directory cannot be created or the file cannot be written, an error is returned.
pub fn write_image_grid(imr: &[f64], r_sph_val: f64) -> Result<()> {
    fs::create_dir_all("output")?;

    let file = File::create("output/imagegrid.out")?;
    let mut writer = BufWriter::new(file);

    writeln!(writer, "{}", r_sph_val / cc::astro::AU_CGS)?;
    writeln!(writer)?; // empty line

    for val in imr {
        writeln!(writer, "{}", val / cc::astro::AU_CGS)?;
    }

    Ok(())
}
