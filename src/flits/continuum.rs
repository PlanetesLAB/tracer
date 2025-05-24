//! Continuum ray-tracing

use crate::types::RMatrix;

use super::{
    config::{FLITSConfig, Grid},
    utils::planck,
};

pub fn ray_trace_continuum(flitsc: &mut FLITSConfig, grid: &Grid) {
    flitsc.bb = RMatrix::zeros((grid.nlam, flitsc.maxt));

    flitsc.bb.indexed_iter_mut().for_each(|((ilam, i), val)| {
        let t = (i + 1) as f64;
        let lam = grid.lam_cont[ilam];
        *val = planck(t, lam)
    });
}
