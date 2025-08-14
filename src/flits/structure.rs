use std::collections::HashMap;

use crate::flits::config::{FLITSConfig, SimState};

pub fn prepare_structure(sim_state: &mut SimState, flitsc: &FLITSConfig) {
    let c = &mut sim_state.c;
    let grid = &sim_state.grid;

    for i in 0..grid.nr {
        c[[i, 0]].dens = 1e-50;
    }

    for j in 0..grid.ntheta {
        c[[0, j]].dens = 1e-50;
        if flitsc.cylindrical {
            c[[grid.nr, j]].dens = 1e-50;
        }
    }

    let molecules = &mut sim_state.mol;
    let mol_names = &flitsc.mol_name0;

    let name_to_index: HashMap<&str, usize> = mol_names
        .iter()
        .enumerate()
        .map(|(i, name)| (name.trim(), i))
        .collect();

    let mut ispec = Vec::with_capacity(molecules.len());
    let mut maxlevels = 0;

    for mol in molecules.iter_mut() {
        let species_index = name_to_index.get(mol.name.trim()).copied();

        ispec.push(species_index);

        if let Some(idx) = species_index {
            if flitsc.npop0[idx] > maxlevels {
                maxlevels = flitsc.npop0[idx];
            }
        } else {
            println!("Species {} not found", mol.name.trim());
            println!("removing this species");
            mol.LTE = true;
        }

        if mol.levels.len() > maxlevels {
            maxlevels = mol.levels.len();
        }
    }
}
