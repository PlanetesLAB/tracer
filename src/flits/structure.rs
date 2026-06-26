use std::collections::HashMap;

use codata as cc;
use planetes_ext::types::{RVector, UVector};

use super::config::{FLITSConfig, PLVector, SimState};
use super::local::compute_LTE;

pub const MAXT: usize = 5000;

#[allow(clippy::assigning_clones)]
pub fn prepare_structure(state: &mut SimState, flitsc: &FLITSConfig) {
    let c = &mut state.c;
    let grid = &state.grid;

    for i in 0..grid.nr {
        c[[i, 0]].dens = 1e-50;
    }

    for j in 0..grid.ntheta {
        c[[0, j]].dens = 1e-50;
        if flitsc.cylindrical {
            c[[grid.nr, j]].dens = 1e-50;
        }
    }

    let molecules = &mut state.mol;
    let mol_names = &flitsc.mol_name0;

    let nmol = molecules.len();

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

    for i in 0..=grid.nr {
        for j in 0..=grid.ntheta {
            c[[i, j]].n = RVector::zeros(nmol);
            c[[i, j]].npop = PLVector::default(nmol);
            c[[i, j]].line_width = RVector::zeros(nmol);

            for imol in 0..nmol {
                let npop_n_size = match ispec[imol] {
                    Some(idx) => flitsc.npop0[idx].max(molecules[imol].levels.len()),
                    None => molecules[imol].levels.len(), // or handle differently if needed
                };
                c[[i, j]].npop[imol].n = RVector::zeros(npop_n_size);
                if let Some(idx) = ispec[imol] {
                    c[[i, j]].n[imol] = c[[i, j]].n0[idx];
                    c[[i, j]].line_width[imol] = c[[i, j]].line_width0[idx];
                    if c[[i, j]].line_width[imol] < flitsc.vres_profile * 3.0 {
                        c[[i, j]].line_width[imol] = flitsc.vres_profile * 3.0;
                    }
                    c[[i, j]].npop[imol].n = c[[i, j]].npop0[idx].n.clone();
                } else {
                    c[[i, j]].n[imol] = 1e-70;
                    c[[i, j]].line_width[imol] = 1e5;
                }
            }
        }
    }

    for i in 0..=grid.nr {
        for j in 0..=grid.ntheta {
            c[[i, j]].it = (c[[i, j]].tdust + 0.5) as usize;
            c[[i, j]].it = c[[i, j]].it.clamp(1, MAXT);
            c[[i, j]].npopmax = UVector::ones(nmol);
            for (imol, item) in ispec.iter().enumerate().take(nmol) {
                if let Some(idx) = *item {
                    // Find the first ipop in reverse where the condition is true
                    if let Some(ipop) = (0..flitsc.npop0[idx])
                        .rev()
                        .find(|&ipop| c[[i, j]].npop[imol].n[ipop] > 1e-150_f64)
                        && ipop > c[[i, j]].npopmax[imol]
                    {
                        c[[i, j]].npopmax[imol] = ipop;
                    }
                }
            }
        }
    }

    compute_LTE(state);

    let rlines = 1.0
        / (((1.0 + flitsc.vresolution / cc::SPEED_OF_LIGHT_IN_VACUUM_CGS)
            / (1.0 - flitsc.vresolution / cc::SPEED_OF_LIGHT_IN_VACUUM_CGS))
            .sqrt()
            - 1.0);

    println!("Mass of the star: {} solar masses", flitsc.m_star);
    println!("Minimum wavelength: {} microns", flitsc.lmin);
    println!("Maximum wavelength: {} microns", flitsc.lmax);
    println!("Resolution: {rlines} (dlam/lam)");
    println!("Velocity resolution: {} km/s", flitsc.vresolution / 1e5);
    println!("Inclination: {} degrees", flitsc.inc);
}
