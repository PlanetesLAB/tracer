use std::f64::consts::PI;

use anyhow::Result;

use crate::{
    flits::{
        config::{AccuracyLevel, FLITSConfig, Grid, SimState},
        io::write_image_grid,
    },
    types::UVector,
    utils::ArrayExtrema,
};

pub fn setup_paths(state: &SimState, flitsc: &FLITSConfig) -> Result<()> {
    let grid = &state.grid;
    let c = &state.c;

    let inc = flitsc.inc * PI / 180.0;
    let cos_inc = inc.cos();

    let mut ilam1 = 0;
    let mut ilam2 = grid.nlam;

    for ilam in 0..grid.nlam {
        if grid.lam_cont[ilam] < flitsc.lmin {
            ilam1 = ilam;
        }
        if grid.lam_cont[grid.nlam - 1 - ilam] > flitsc.lmax {
            ilam2 = grid.nlam - 1 - ilam;
        }
    }

    let (mut nr_reduce, res_inc, n_phi_min, n_phi_max) = match flitsc.accuracy {
        AccuracyLevel::MinSampling => (8, 1.0, 2, 10),
        AccuracyLevel::BasicSampling => (4, 1.0, 15, 60),
        AccuracyLevel::DefaultSampling => (2, 1.0, 30, 75),
        AccuracyLevel::RadialEnhanced => (1, 2.0, 30, 90),
        AccuracyLevel::AzimuthalEnhanced => (1, 4.0, 45, 90),
        AccuracyLevel::DualEnhanced => (1, 8.0, 45, 180),
        AccuracyLevel::MaxSampling => (1, 10.0, 120, 270),
    };

    while grid.nr / nr_reduce < 40 && nr_reduce > 1 {
        nr_reduce -= 1;
    }

    let mut imr = Vec::new();
    for i in (1..=grid.ntheta).rev() {
        let rr = trace_to_tau(i, state, flitsc);
        if rr > grid.r_sph[1] {
            imr.extend(vec![
                rr,
                rr * (inc - (PI / 2.0 - grid.theta_av[[grid.nr - 1, i]])).cos(),
                rr * (inc + (PI / 2.0 - grid.theta_av[[grid.nr - 1, i]])).cos(),
            ]);
        }
    }
    for i in (1..=grid.nr).step_by(nr_reduce) {
        imr.push(grid.r_av_sph[i]);
        if flitsc.cylindrical {
            imr.push(grid.r_av[i]);
        }
    }

    imr.sort_by(|a, b| a.total_cmp(b));

    let initial_len = imr.len();
    let buffer_len = (inc.sin().abs() * c[[1, grid.ntheta]].v / flitsc.vresolution * res_inc / 2.0)
        as usize
        + (grid.ntheta - 1) * 2
        + 100;
    let new_len = initial_len + buffer_len;

    imr.resize(new_len, 0.0);

    let mut cursor = initial_len - 1;
    while imr[cursor] > grid.r_sph[grid.nr + 1] {
        cursor -= 1;
    }

    for i in 2..=grid.ntheta {
        let rr = if flitsc.cylindrical {
            grid.r[1] / grid.theta_av[[grid.nr - 1, i]].sin()
        } else {
            grid.r_sph[1]
        };
        cursor += 1;

        let num_1 = rr * (inc - (PI / 2.0 - grid.theta_av[[grid.nr - 1, i]])).cos();
        let val_1 = num_1 / get_inclination_factor(rr, cos_inc, grid);
        let val_2 = num_1 / get_inclination_factor(val_1, cos_inc, grid);
        imr[cursor] = num_1 / get_inclination_factor(val_2, cos_inc, grid);
        cursor += 1;
        let num_2 = rr * (inc + (PI / 2.0 - grid.theta_av[[grid.nr - 1, i]])).cos();
        let val_1 = num_2 / get_inclination_factor(rr, cos_inc, grid);
        let val_2 = num_2 / get_inclination_factor(val_1, cos_inc, grid);
        imr[cursor] = num_2 / get_inclination_factor(val_2, cos_inc, grid);
    }

    let j = new_len - cursor - 1;

    for i in 0..j {
        cursor += 1;
        let r_log10 = grid.r[grid.nr].log10();
        let r_sph_log10 = (grid.r_sph[grid.nr + 1] / grid.r[grid.nr]).log10();
        imr[cursor] = 10.0f64.powf(r_log10 + r_sph_log10 * ((i + 1) as f64 - 0.1) / (j as f64));
    }

    imr.iter_mut().for_each(|val| {
        *val = val.abs();
    });

    imr.sort_by(|a, b| a.total_cmp(b));

    write_image_grid(imr, grid.r_sph[grid.nr + 1])?;

    let nim_phi = UVector::zeros(new_len);
    let start_phi = UVector::zeros(new_len);

    let mut n_imphi_max = 1;

    Ok(())
}

fn trace_to_tau(j: usize, state: &SimState, flitsc: &FLITSConfig) -> f64 {
    let grid = &state.grid;
    let c = &state.c;

    let mut tau = 0.0;
    for i in 1..=grid.nr {
        let tau_0 = if flitsc.cylindrical && j != 0 {
            (&c[[i, j]].kext * (grid.r[i + 1] - grid.r[i]) / grid.theta_av[[i, j]].sin())
                .maxval()
                .unwrap_or(0.0)
        } else {
            (&c[[i, j]].kext * (grid.r_sph[i + 1] - grid.r_sph[i]))
                .maxval()
                .unwrap_or(0.0)
        };
        if tau + tau_0 > 1.0 {
            if flitsc.cylindrical {
                let mut rr = grid.r[i] + (grid.r[i + 1] - grid.r[i]) * (1.0 - tau) / tau_0;
                rr /= grid.theta_av[[i, j]].sin();
                return rr;
            } else {
                let rr = grid.r_sph[i] + (grid.r_sph[i + 1] - grid.r_sph[i]) * (1.0 - tau) / tau_0;
                return rr;
            }
        }
        tau += tau_0;
    }
    0.0
}

fn get_inclination_factor(r: f64, cos_inc: f64, grid: &Grid) -> f64 {
    cos_inc
        + (1.0 - cos_inc)
            * ((r - grid.r_sph[1]) / grid.r_sph[grid.nr + 1] - grid.r_sph[1])
                .powi(2)
                .min(cos_inc)
}
