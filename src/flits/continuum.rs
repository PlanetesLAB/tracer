//! Continuum ray-tracing

use extensions::types::{RMatrix, RVector};

use crate::utils::physics::planck;

use super::config::{FLITSConfig, Grid, Path};

pub fn ray_trace_continuum(flitsc: &mut FLITSConfig, grid: &Grid) {
    let mut bb = RMatrix::zeros((grid.nlam, flitsc.maxt));

    bb.indexed_iter_mut().for_each(|((ilam, i), val)| {
        let t = (i + 1) as f64;
        let lam = grid.lam_cont[ilam];
        *val = planck(t, lam)
    });

    let mut p: Vec<Vec<Path>> = Vec::with_capacity(flitsc.image.ngrids);
    for i in 0..flitsc.image.ngrids {
        for j in 0..flitsc.image.npoints[i] {
            let mut pth = Path::default();
            pth.flux_cont = RVector::zeros(grid.nlam);
            p.push(vec![pth]);
        }
    }
}
