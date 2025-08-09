use std::fs::File;
use std::io;

use std::io::{BufReader, BufWriter};

use anyhow::Result;
use fitsio::FitsFile;
use flate2::read::GzDecoder;
use ndarray::s;
use ndarray::{Array3, ArrayD, Axis};

use crate::constants as cc;
use crate::flits::config::{FLITSConfig, SimState};
use crate::types::{CellMatrix, RMatrix, RVector};

pub fn read_fits(flitsc: &mut FLITSConfig) -> Result<()> {
    flitsc.flits_file = flitsc
        .flits_file_compressed
        .strip_suffix(".gz")
        .map(|s| s.to_string())
        .ok_or(anyhow::anyhow!("FLITS file name must end with .gz"))?;

    let fits_path = flitsc.config_dir.join(&flitsc.flits_file);

    if !fits_path.exists() {
        uncompress_fits_file(flitsc)?;
    }

    let mut sim_state = SimState::default();

    let mut fptr = FitsFile::open(&fits_path)?;

    // Check dimensions
    let hdu = fptr.primary_hdu()?;
    let mut naxes = [0f64; 3];
    for (i, naxis) in naxes.iter_mut().enumerate() {
        let key_name = format!("NAXIS{}", i + 1);
        *naxis = hdu.read_key::<f64>(&mut fptr, &key_name)?;
    }

    flitsc.grid.rin = hdu.read_key::<f64>(&mut fptr, "Rin")?;
    flitsc.grid.rout = hdu.read_key::<f64>(&mut fptr, "Rout")?;
    flitsc.m_star = hdu.read_key::<f64>(&mut fptr, "Mstar")?;
    flitsc.r_star = hdu.read_key::<f64>(&mut fptr, "Rstar")?;
    flitsc.distance = hdu.read_key::<f64>(&mut fptr, "distance")?;
    flitsc.grid.nr = hdu.read_key::<f64>(&mut fptr, "NXX")? as usize;
    flitsc.grid.ntheta = hdu.read_key::<f64>(&mut fptr, "NZZ")? as usize;
    flitsc.grid.nspec = hdu.read_key::<f64>(&mut fptr, "NSPEC")? as usize;

    flitsc.grid.nr += 1;
    flitsc.grid.ntheta += 1;

    sim_state.c = CellMatrix::default(((flitsc.grid.nr + 1), (flitsc.grid.ntheta + 1)));
    flitsc.grid.r = RVector::zeros(flitsc.grid.nr + 2);
    flitsc.grid.theta = RMatrix::default(((flitsc.grid.nr + 1), (flitsc.grid.ntheta + 2)));

    flitsc.grid.r_av = RVector::zeros(flitsc.grid.nr + 1);
    flitsc.grid.theta_av = RMatrix::zeros(((flitsc.grid.nr + 1), (flitsc.grid.ntheta + 1)));

    // Axis(0) --> r
    // Axis(1) --> theta
    // Axis(2) --> ???
    let data: ArrayD<f64> = hdu.read_image(&mut fptr)?;
    let data3: Array3<f64> = data.into_dimensionality()?;
    let array = data3.permuted_axes([2, 1, 0]).insert_axis(Axis(3));

    flitsc.grid.rin = array[[0, 0, 0, 0]] / cc::AU_CGS;
    flitsc.grid.rout = array[[flitsc.grid.nr - 2, 0, 1, 0]] / cc::AU_CGS;

    flitsc.grid.r[0] = flitsc.r_star * cc::SOLAR_RADIUS_CGS;

    flitsc
        .grid
        .r
        .slice_mut(s![1..=flitsc.grid.nr - 1])
        .assign(&array.slice(s![.., 0, 0, 0]));

    flitsc.grid.r[flitsc.grid.nr] = array[[flitsc.grid.nr - 2, 0, 0, 0]];

    for i in 1..=(flitsc.grid.nr - 1) {
        flitsc.grid.r_av[i] = (flitsc.grid.r[i] + flitsc.grid.r[i + 1]).sqrt();
    }

    let mut theta_max = 0.0;

    for i in 1..=(flitsc.grid.nr - 1) {
        flitsc.grid.theta[[i, 0]] = 1.0;
        for j in 2..flitsc.grid.ntheta {
            flitsc.grid.theta[[i, j]] = array[[i - 1, flitsc.grid.ntheta - j, 3, 0]]
                / (flitsc.grid.r_av[i].powi(2)
                    + array[[i - 1, flitsc.grid.ntheta - j, 3, 0]].powi(2))
                .sqrt()
        }
        flitsc.grid.theta[[i, flitsc.grid.ntheta + 1]] = 0.0;
        flitsc.grid.theta_av[[i, 1]] =
            (flitsc.grid.theta[[i, 0]] + flitsc.grid.theta[[i, 2]]) / 2.0;

        if flitsc.cylindrical {
            flitsc.grid.theta_av[[i, 1]] = flitsc.grid.theta[[i, 1]].acos() / 2.0;
        } else {
            flitsc.grid.theta_av[[i, 1]] = 1.0;
            flitsc.grid.theta_av[[i, 0]] = 0.0;
        }

        for j in 1..=(flitsc.grid.ntheta) {
            flitsc.grid.theta_av[[i, j]] =
                ((flitsc.grid.theta[[i, j]] + flitsc.grid.theta[[i, j + 1]]) / 2.0).acos()
        }

        if flitsc.grid.theta[[i, 2]] > theta_max {
            theta_max = flitsc.grid.theta[[i, 2]];
        }
    }

    Ok(())
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
