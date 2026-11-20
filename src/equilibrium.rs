//! Optically thin statistical equilibrium in a prescribed radiation field.

use anyhow::Result;
use anyhow::{anyhow, ensure};
use nalgebra::{DMatrix, DVector};
use planetes_ext::types::{RMatrix, RVector};

use crate::collision::Level;
use crate::collision::check_levels;

/// Radiative transition with a prescribed background at its line frequency.
pub struct RadiativeTransition {
    pub upper: usize,
    pub lower: usize,
    /// Spontaneous emission rate in $\mathrm{s^{-1}}$.
    pub einstein_a: f64,
    /// Dimensionless photon occupation $c^2 J_\nu/(2h\nu^3)$.
    /// Set to zero for spontaneous emission into vacuum.
    pub photon_occupation: f64,
}

/// Solves for normalized level populations at fixed collision and radiation rates.
/// Supply the collision matrix in $\mathrm{s^{-1}}$, indexed by initial and
/// final level with zero diagonal. Returned fractions follow input level order.
///
/// # Errors
///
/// Requires valid levels, a matching finite nonnegative collision matrix,
/// and valid downward transitions with finite nonnegative rates and occupations.
/// Returns an error if normalization cannot uniquely fix the steady state or
/// the linear solve cannot represent its populations.
pub fn optically_thin_populations(
    levels: &[Level],
    collisions: &RMatrix,
    lines: &[RadiativeTransition],
) -> Result<RVector> {
    check_levels(levels)?;
    let n = levels.len();
    ensure!(
        collisions.dim() == (n, n),
        "collision matrix does not match level count"
    );
    ensure!(
        collisions.diag().iter().all(|&x| x == 0.0),
        "expected zero collision diagonal"
    );
    ensure!(
        collisions.iter().all(|x| x.is_finite() && *x >= 0.0),
        "invalid collision rate"
    );
    let mut rates = collisions.clone();
    for line in lines {
        let (u, l) = (line.upper, line.lower);
        ensure!(
            u < n && l < n && levels[u].wavenumber > levels[l].wavenumber,
            "invalid radiative transition {u} -> {l}"
        );
        let photons = line.photon_occupation;
        let a = line.einstein_a;
        ensure!(
            a.is_finite() && a >= 0.0 && photons.is_finite() && photons >= 0.0,
            "invalid radiation rate or photon occupation"
        );
        rates[[u, l]] += a * (1.0 + photons);
        rates[[l, u]] += a * photons * (levels[u].weight / levels[l].weight);
    }
    ensure!(
        rates.iter().all(|x| x.is_finite()),
        "transition rate overflow"
    );
    stationary_populations(&rates)
}

// Follow directed transitions. In reverse, this marks levels that can feed the seed.
fn mark_reachable(rates: &RMatrix, seed: usize, reverse: bool, visited: &mut [bool]) {
    let mut pending = vec![seed];
    visited[seed] = true;
    while let Some(i) = pending.pop() {
        let edges = if reverse {
            rates.column(i)
        } else {
            rates.row(i)
        };
        for (j, &rate) in edges.iter().enumerate() {
            if rate > 0.0 && !visited[j] {
                visited[j] = true;
                pending.push(j);
            }
        }
    }
}

fn stationary_populations(rates: &RMatrix) -> Result<RVector> {
    let n = rates.nrows();
    let mut visited = vec![false; n];
    let mut anchor = 0;
    // The last reverse-search root belongs to the unique closed class if one exists.
    for i in 0..n {
        if !visited[i] {
            anchor = i;
            mark_reachable(rates, i, true, &mut visited);
        }
    }
    visited.fill(false);
    mark_reachable(rates, anchor, true, &mut visited);
    ensure!(
        visited.iter().all(|&v| v),
        "multiple closed sets of levels; their relative populations are undetermined"
    );
    visited.fill(false);
    mark_reachable(rates, anchor, false, &mut visited);
    let active: Vec<_> = (0..n).filter(|&i| visited[i] && i != anchor).collect();
    let mut population = RVector::zeros(n);
    population[anchor] = 1.0;
    if active.is_empty() {
        return Ok(population);
    }

    // Solve inflow = outflow relative to the anchor population. Transient
    // levels have exactly zero population and never enter the factorization.
    let m = active.len();
    let mut balance = DMatrix::zeros(m, m);
    let mut inflow = DVector::zeros(m);
    for (row, &i) in active.iter().enumerate() {
        let loss = rates.row(i).sum();
        ensure!(loss.is_finite(), "total transition rate overflow");
        for (col, &j) in active.iter().enumerate() {
            balance[(row, col)] = if i == j { loss } else { -rates[[j, i]] };
        }
        inflow[row] = rates[[anchor, i]];
        let scale = balance.row(row).amax().max(inflow[row]);
        balance.row_mut(row).apply(|v| *v /= scale);
        inflow[row] /= scale;
    }
    let relative = balance
        .lu()
        .solve(&inflow)
        .ok_or_else(|| anyhow!("singular level balance matrix"))?;
    for (&i, &value) in active.iter().zip(relative.iter()) {
        ensure!(
            value.is_finite() && value >= 0.0,
            "unresolved equilibrium population at level {i}"
        );
        population[i] = value;
    }
    let peak = relative.amax().max(1.0);
    population /= peak;
    let total = population.sum();
    Ok(population / total)
}
