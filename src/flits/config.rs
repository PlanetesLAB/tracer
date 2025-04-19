//! Global configuration for the FLiTs program.

use ndarray::Array2;

use crate::types::{BVector, RMatrix, RTensor, RVector, UVector};

#[allow(non_snake_case)]
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
    pub line_file: String,
    pub tracer_file: String,

    /// Model Structure
    pub model: Model,
    // Grid setup
    pub grid: Grid,

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
            line_file: String::new(),
            maxt,
            model: Model::ProDiMo,
            npop0: UVector::default(maxt),
            nv_profile: 10,
            profile: LineProfile::default(),
            rlines: 10.0,
            vres_mult: 1e5,
            tracer_file: String::new(),
            bb: RMatrix::default((nlam, maxt)),
            grid: Grid::new(nlam),
            vresolution: 1e5,
            vres_profile: 1e4,
            tau_max: 15.0,
            accuracy: AccuracyLevel::Low,
            LTE: false,
            cylindrical: true,
            do_blend: true,
            imcube: ImageCube::default(),
            image: Image::default(),
        }
    }
}

impl FLITSConfig {
    fn init(&mut self) {
        todo!()
    }
    fn read(&mut self) {
        todo!()
    }
    fn prepare_structure(&mut self) {
        todo!()
    }
    fn setup_paths(&mut self) {
        todo!()
    }
    fn sort_lines(&mut self) {
        todo!()
    }
    fn write_line_list(&mut self) {
        todo!()
    }
    fn ray_trace_continuum(&mut self) {
        todo!()
    }
    fn ray_trace_lines(&mut self) {
        todo!()
    }
    pub fn run(&mut self) {
        self.init();
        self.read();
        self.prepare_structure();
        self.setup_paths();
        self.sort_lines();
        self.write_line_list();
        self.ray_trace_continuum();
        self.ray_trace_lines();
    }
}

#[derive(Default)]
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
    pub nmol: usize,
    pub nlines: usize,
    pub nblends: usize,
    pub nspec: usize,
    pub rin: f64,
    pub rout: f64,
    pub inc: f64,
    pub lstar: f64,
    pub lam_cont: RVector,
    pub fstar: RVector,
}

impl Grid {
    fn new(nlam: usize) -> Self {
        Self {
            nlam,
            inc: 30.0,
            ..Default::default()
        }
    }
}

#[derive(Default)]
pub struct Image {
    pub nimr: usize,
    pub ngrids: usize,
    pub npoints_tmp: usize,
    pub nimphi: UVector,
    pub npoints: UVector,
    pub vmax: f64,
}

#[derive(Default)]
pub struct ImageCube {
    pub xim: RTensor,
    pub yim: RTensor,
    pub data: RTensor,
    pub nint: usize,
    pub npix: usize,
    pub nvim: usize,
    pub imagecube: bool,
}

#[derive(Default)]
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

#[allow(non_snake_case)]
pub struct Molecule {
    pub energy: RVector,
    pub g: RVector,
    pub nlines: usize,
    pub nlevels: usize,
    pub mass: f64,
    pub lines: Vec<Line>,
    pub name: String,
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
    pub npop0: RMatrix,
    pub npop: RMatrix,
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
}

pub enum AccuracyLevel {
    Low,
    Medium,
    High,
}

pub enum Model {
    ProDiMo,
    McMax,
}
