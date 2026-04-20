use std::f64::consts::PI;

use anyhow::Result;
use extensions::types::{BVector, Matrix, RMatrix, UMatrix, UVector};
use ndarray::array;
use rand::Rng;

use crate::{
    flits::{
        config::{AccuracyLevel, FLITSConfig, Grid, Path, Rotate, SimState},
        delaunay::Point,
        io::write_image_grid,
    },
    utils::ArrayExtrema,
};

pub type PathMatrix = Matrix<Path>;

/// Sets up the paths for ray tracing based on the provided simulation state and configuration.
///
/// # Errors
/// Returns an error if there is an issue writing the image grid to the output directory.
pub fn setup_paths(state: &SimState, flitsc: &FLITSConfig) -> Result<()> {
    const N_GRIDS: usize = 5;

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

    imr.sort_by(f64::total_cmp);

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

    for val in &mut imr {
        *val = val.abs();
    }

    imr.sort_by(f64::total_cmp);

    write_image_grid(&imr, grid.r_sph[grid.nr + 1])?;

    let mut nim_phi = UVector::zeros(new_len);
    let mut start_phi = BVector::default(new_len + 1);

    let mut n_imphi_max = 0;
    start_phi[0] = true;

    let mut n_pts_temp = 0;

    for i in 0..new_len {
        let j = (1..grid.nr)
            .find(|&j| imr[i] > grid.r[j] && imr[i] < grid.r[j + 1])
            .unwrap_or(0);
        nim_phi[i] =
            (inc.sin().abs() * c[[j, grid.ntheta]].v / flitsc.vresolution * res_inc) as usize;
        if nim_phi[i] < n_phi_min {
            nim_phi[i] = n_phi_min;
        }
        if nim_phi[i] > n_phi_max {
            nim_phi[i] = n_phi_max;
        }
        if start_phi[i] {
            nim_phi[i] += 1;
            start_phi[i + 1] = false;
        } else {
            start_phi[i + 1] = true;
        }
        if nim_phi[i] > n_imphi_max {
            n_imphi_max = nim_phi[i];
        }
        n_pts_temp += nim_phi[i];
    }

    let mut imphi = RMatrix::zeros((new_len, n_imphi_max));

    for i in 0..new_len {
        if start_phi[i] {
            for j in 0..nim_phi[i] {
                imphi[[i, j]] = PI * (j as f64 + 0.5) / (nim_phi[i] as f64);
            }
        } else {
            for j in 0..nim_phi[i] {
                imphi[[i, j]] = PI * (j as f64) / ((nim_phi[i] - 1) as f64);
            }
        }
    }

    let mut xy = Vec::new();
    let n_pts = UVector::zeros(N_GRIDS);

    let mat_len = 2 * n_pts_temp - 3;

    let mut p = PathMatrix::default((N_GRIDS, mat_len));
    let mut t_node = UMatrix::zeros((3, mat_len));
    let mut t_neigh = UMatrix::zeros((3, mat_len));
    let matrix = array![[1.0, 0.0], [0.0, 1.0]];

    let mut rng = rand::rng();

    let mut k = 0;
    for i in 0..N_GRIDS {
        for j in 0..n_pts_temp {
            let ir: usize = rng.random_range(0..=grid.nr);
            let it: usize = rng.random_range(0..=grid.ntheta);

            let rrr = rng.random::<f64>();
            let rr = (grid.r[ir].powi(2) * rrr + grid.r[ir + 1].powi(2) * (1.0 - rrr)).sqrt();

            let ttt = rng.random::<f64>();
            let tt = (grid.theta[[ir, it]] * ttt + grid.theta[[ir, it + 1]] * (1.0 - ttt)).acos();

            let ph = rng.random::<f64>() * 2.0 * PI;

            let mut im_vec = array![
                rr * ph.cos() * tt.sin(),
                rr * ph.sin() * tt.sin(),
                rr * tt.cos()
            ];
            let axis = array![0.0, 1.0, 0.0];
            let theta = -inc;
            im_vec.rotate(&axis, theta);

            xy.push(Point {
                x: im_vec[0],
                y: im_vec[1],
            });
        }
        // triangulate(&mut xy);
    }

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
            }
            let rr = grid.r_sph[i] + (grid.r_sph[i + 1] - grid.r_sph[i]) * (1.0 - tau) / tau_0;
            return rr;
        }
        tau += tau_0;
    }
    0.0
}

fn get_inclination_factor(r: f64, cos_inc: f64, grid: &Grid) -> f64 {
    cos_inc
        + (1.0 - cos_inc)
            * ((r - grid.r_sph[1]) / (grid.r_sph[grid.nr + 1] - grid.r_sph[1]))
                .powi(2)
                .min(cos_inc)
}
