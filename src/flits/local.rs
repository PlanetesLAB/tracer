//! Local Thermodynamic Equilibrium (LTE) calculations.

use super::config::{Grid, SimState};

#[allow(non_snake_case)]
pub fn compute_LTE(state: &mut SimState, grid: &Grid) {
    for i in 0..grid.nr {
        for j in 0..grid.ntheta {
            let t = state.c[[i, j]].tgas.max(3.0);
            for imol in 0..grid.nmol {
                if state.mol[imol].LTE {
                    let mut ut = 0.0;
                    for k in 0..state.mol[imol].nlevels {
                        ut += state.mol[imol].g[k] * (-state.mol[imol].energy[k] / t).exp();
                    }
                    for k in 0..state.mol[imol].nlevels {
                        state.c[[i, j]].npop[[imol, k]] =
                            state.mol[imol].g[k] * (-state.mol[imol].energy[k] / t).exp() / ut;
                    }
                }
            }
        }
    }
}
