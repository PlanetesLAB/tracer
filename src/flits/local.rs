//! Local Thermodynamic Equilibrium (LTE) calculations.

use super::config::SimState;

#[allow(non_snake_case)]
pub fn compute_LTE(state: &mut SimState) {
    let grid = &state.grid;
    let mol = &state.mol;

    let c = &mut state.c;

    for i in 0..=grid.nr {
        for j in 1..=grid.ntheta {
            let t = c[[i, j]].tgas.max(3.0);
            for (imol, molecule) in mol.iter().enumerate() {
                if molecule.LTE {
                    let mut ut = 0.0;
                    for k in 0..molecule.levels.len() {
                        ut += molecule.levels[k].weight * (-molecule.levels[k].energy / t).exp();
                    }
                    for k in 0..molecule.levels.len() {
                        c[[i, j]].npop[imol].n[k] =
                            molecule.levels[k].weight * (-molecule.levels[k].energy / t).exp() / ut;
                    }
                }
            }
        }
    }
}
