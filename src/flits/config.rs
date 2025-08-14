//! Global configuration for the FLiTs program.
use std::io::BufRead;
use std::path::PathBuf;
use std::{fs::File, io::BufReader, path::Path as OSPath};

use anyhow::Result;
use ndarray::Array2;
use spectre::lamda::Level;

use crate::types::{BVector, PLVector, RMatrix, RTensor, RVector, UVector};

#[allow(non_snake_case)]
#[derive(Debug)]
pub struct FLITSConfig {
    // Stellar parameters
    /// Stellar mass (solar masses).
    pub m_star: f64,
    /// Stellar radius (solar radii).
    pub r_star: f64,
    /// Distance to the star (parsecs).
    pub distance: f64,

    // Wavelength grid and resolution
    /// Minimum wavelength (microns).
    pub lmin: f64,
    /// Maximum wavelength (microns).
    pub lmax: f64,
    pub rlines: f64,
    pub vresolution: f64,
    pub vres_mult: f64,
    pub tau_max: f64,
    pub vres_profile: f64,
    pub nv_profile: usize,
    pub accuracy: AccuracyLevel,
    pub do_blend: bool,
    pub cylindrical: bool,
    pub LTE: bool,

    // Files
    pub config_path: PathBuf,
    pub config_dir: PathBuf,
    pub line_files: Vec<String>,
    pub tracer_file: String,
    pub flits_file_compressed: String,
    pub flits_file: String,

    /// Model Structure
    pub model: Model,

    // Image grid
    pub image: Image,

    // Image cube
    pub imcube: ImageCube,

    /// Blackbodies
    pub maxt: usize,
    pub bb: RMatrix,

    /// Profiles
    pub profile: LineProfile,

    pub npop0: UVector,
    pub mol_name0: Vec<String>,

    pub inc: f64,
}

impl Default for FLITSConfig {
    fn default() -> Self {
        let nlam = 100;
        let maxt = 5000;
        FLITSConfig {
            m_star: 1.0,
            r_star: 1.0,
            lmin: 5.0,
            lmax: 50.0,
            distance: 1.0,
            config_path: PathBuf::new(),
            config_dir: PathBuf::new(),
            line_files: Vec::new(),
            flits_file_compressed: String::new(),
            flits_file: String::new(),
            maxt,
            model: Model::ProDiMo,
            npop0: UVector::default(maxt),
            mol_name0: Vec::new(),
            nv_profile: 10,
            profile: LineProfile::default(),
            rlines: 10.0,
            vres_mult: 1e5,
            tracer_file: String::new(),
            bb: RMatrix::default((nlam, maxt)),
            vresolution: 1e5,
            vres_profile: 1e4,
            tau_max: 15.0,
            accuracy: AccuracyLevel::DefaultSampling,
            LTE: false,
            cylindrical: true,
            do_blend: true,
            imcube: ImageCube::default(),
            image: Image::default(),
            inc: 30.0,
        }
    }
}

impl FLITSConfig {
    pub fn load_from_path<P: AsRef<OSPath>>(&mut self, path: P) -> Result<()> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);

        let lines: Vec<String> = reader.lines().map_while(Result::ok).collect();

        for line in lines {
            let (key, value) = line.split_once('=').ok_or_else(|| {
                anyhow::anyhow!("Invalid line format, expected 'key=value': {}", line)
            })?;
            // Remove quotes
            let value = value.trim_matches('\'');
            if key == "linefile" {
                self.line_files.push(value.to_string());
            }
            if key == "lmin" {
                let value = value.split_whitespace().next().ok_or_else(|| {
                    anyhow::anyhow!("Expected a value for lmin in line: {}", line)
                })?;
                self.lmin = value.parse::<f64>().expect("Failed to parse lmin");
            }
            if key == "lmax" {
                let value = value.split_whitespace().next().ok_or_else(|| {
                    anyhow::anyhow!("Expected a value for lmax in line: {}", line)
                })?;
                self.lmax = value.parse::<f64>().expect("Failed to parse lmax");
            }
            if key == "vres" {
                let value = value.split_whitespace().next().ok_or_else(|| {
                    anyhow::anyhow!("Expected a value for vres in line: {}", line)
                })?;
                self.vresolution = value.parse::<f64>().expect("Failed to parse vres");
            }
            if key == "vres_profile" {
                let value = value.split_whitespace().next().ok_or_else(|| {
                    anyhow::anyhow!("Expected a value for vres_profile in line: {}", line)
                })?;
                if value.contains('d') {
                    // replace with e
                    let new_value = value.replace('d', "e");
                    self.vres_profile = new_value
                        .parse::<f64>()
                        .expect("Failed to parse vres_profile");
                } else {
                    self.vres_profile = value.parse::<f64>().expect("Failed to parse vres_profile");
                }
            }
            if key == "LTE" {
                if value.contains("true") || value.contains("1") {
                    self.LTE = true;
                } else if value.contains("false") || value.contains("0") {
                    self.LTE = false;
                } else {
                    return Err(anyhow::anyhow!("Invalid value for LTE in line: {}", line));
                }
            }
            if key == "inc" {
                let value = value.trim();
                self.inc = value.parse::<f64>().expect("Failed to parse inc");
            }
            if key == "accuracy" {
                let value = value.trim();
                self.accuracy = match value {
                    "0" => AccuracyLevel::BasicSampling,
                    "1" => AccuracyLevel::DefaultSampling,
                    "2" => AccuracyLevel::RadialEnhanced,
                    "3" => AccuracyLevel::AzimuthalEnhanced,
                    "4" => AccuracyLevel::DualEnhanced,
                    "5" => AccuracyLevel::MaxSampling,
                    _ => {
                        return Err(anyhow::anyhow!(
                            "Invalid value for accuracy in line: {}",
                            line
                        ))
                    }
                };
            }
            if key == "flitsfile" {
                self.flits_file_compressed = value.to_string();
            }
        }

        if self.inc < 1.0 {
            self.inc = 1.0;
        }

        Ok(())
    }
}

#[derive(Default, Debug)]
pub struct Grid {
    pub r: RVector,
    pub theta: RMatrix,
    pub r_av: RVector,
    pub theta_av: RMatrix,
    pub r_sph: RVector,
    pub r_av_sph: RVector,
    pub nr: usize,
    pub ntheta: usize,
    pub nlam: usize,
    pub ilam1: usize,
    pub ilam2: usize,
    pub nlines: usize,
    pub nblends: usize,
    pub nspec: usize,
    pub rin: f64,
    pub rout: f64,
    pub lstar: f64,
    pub lam_cont: RVector,
    pub fstar: RVector,
}

#[derive(Default, Debug)]
pub struct Image {
    pub nimr: usize,
    pub ngrids: usize,
    pub npoints_tmp: usize,
    pub nimphi: UVector,
    pub npoints: UVector,
    pub vmax: f64,
}

#[derive(Default, Debug)]
pub struct ImageCube {
    pub xim: RTensor,
    pub yim: RTensor,
    pub data: RTensor,
    pub nint: usize,
    pub npix: usize,
    pub nvim: usize,
    pub imagecube: bool,
}

#[derive(Default, Debug)]
pub struct LineProfile {
    pub data: RVector,
    pub nz: BVector,
}

pub struct Line {
    pub jup: usize,
    pub jlow: usize,
    pub imol: usize,
    pub aul: f64,
    pub blu: f64,
    pub bul: f64,
    pub freq: f64,
    pub lam: f64,
    pub eup: f64,
}

pub struct Blend {
    pub lines: Vec<Line>,
    pub v: RVector,
    pub n: usize,
    pub nvmin: usize,
    pub nvmax: usize,
    pub next: Option<Box<Blend>>,
    pub ib0: UVector,
    pub nb0: UVector,
    pub lmin: f64,
    pub lmax: f64,
    pub lam: f64,
    pub f: f64,
    pub computing: bool,
    pub done: bool,
}

#[derive(Debug, Default)]
pub struct RadTransition {
    /// ID of the transition
    pub id: usize,
    /// Upper level ID
    pub up: usize,
    /// Lower level ID
    pub low: usize,
    /// Einstein A coefficient (s^-1)
    pub einst_a: f64,
    /// Einstein B coefficient (upper -> lower) (s^-1)
    pub einst_b_ul: f64,
    /// Einstein B coefficient (lower -> upper) (s^-1)
    pub einst_b_lu: f64,
    /// Frequency of the transition (GHz)
    pub freq: f64,
    /// Wavelength of the transition
    pub lam: f64,
    /// Energy of the transition in (K)
    pub energy: f64,
}

#[allow(non_snake_case)]
#[derive(Default, Debug)]
pub struct Molecule {
    pub name: String,
    pub weight: f64,
    pub levels: Vec<Level>,
    pub radset: Vec<RadTransition>,
    pub LTE: bool,
}

#[derive(Default)]
pub struct Cell {
    pub tdust: f64,
    pub tgas: f64,
    pub dens: f64,
    pub v: f64,
    pub line_width: RVector,
    pub n: RVector,
    pub n0: RVector,
    pub line_width0: RVector,
    pub npop0: PLVector,
    pub npop: PLVector,
    pub npopmax: RVector,
    pub line_emis: RVector,
    pub line_abs: RVector,
    pub kext_l: f64,
    pub therm_l: f64,
    pub scat_l: f64,
    pub it: usize,
    pub kabs: RVector,
    pub albedo: RVector,
    pub kext: RVector,
    pub lrf: RVector,
    pub s: RVector,
}

#[derive(Debug, Default)]
pub struct PopulationLevels {
    pub n: RVector,
}

#[derive(Default)]
pub struct Path {
    pub vmin: RVector,
    pub vmax: RVector,
    pub n: usize,
    pub npopmax: UVector,
    pub a: f64,
    pub r1: f64,
    pub r2: f64,
    pub phi1: f64,
    pub phi2: f64,
    pub r: f64,
    pub phi: f64,
    pub vx: f64,
    pub vy: f64,
    pub vz: f64,
    pub x: f64,
    pub y: f64,
    pub v: RVector,
    pub d: RVector,
    pub v1: RVector,
    pub v2: RVector,
    pub i: UVector,
    pub j: UVector,
    pub flux_cont: RVector,
    pub cont_contr: RVector,
    pub exptau_dust: RVector,
    pub s_dust: RVector,
}

pub struct Tracer {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub vx: f64,
    pub vy: f64,
    pub vz: f64,
    pub edge_nr: usize,
    pub i: usize,
    pub j: usize,
    pub on_edge: bool,
}

#[derive(Default)]
pub struct SimState {
    /// Matrix of cells (grid)
    pub c: Array2<Cell>,
    pub p: Array2<Path>,
    pub mol: Vec<Molecule>,
    pub line: Vec<Line>,
    pub blend: Vec<Blend>,
    pub grid: Grid,
}

#[derive(Debug)]
pub enum AccuracyLevel {
    BasicSampling,
    DefaultSampling,
    RadialEnhanced,
    AzimuthalEnhanced,
    DualEnhanced,
    MaxSampling,
}

#[derive(Debug)]
pub enum Model {
    ProDiMo,
    McMax,
}
