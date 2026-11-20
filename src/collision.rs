//! Inelastic collision rates for molecular excitation.

use anyhow::{Result, ensure};
use ndarray::Array2;

/// A molecular level, indexed by its position in the supplied level list.
pub struct Level {
    /// Term energy $E/(hc)$ in $\mathrm{cm^{-1}}$. An arbitrary common zero is allowed.
    pub wavenumber: f64,
    /// Statistical weight (degeneracy).
    pub weight: f64,
}

/// Downward rate coefficients ($\mathrm{cm^{3}\,s^{-1}}$) on the collider's temperature grid.
pub struct DownwardTransition<'a> {
    pub upper: usize,
    pub lower: usize,
    pub coefficients: &'a [f64],
}

/// One collision partner's number density ($\mathrm{cm^{-3}}$) and rate data.
/// Supply each downward transition once per partner.
pub struct Collider<'a> {
    pub density: f64,
    pub temperatures_k: &'a [f64],
    pub transitions: &'a [DownwardTransition<'a>],
}

/// Approximate thermal $\mathrm{H_2}$ densities (para, ortho), in the units of the total.
/// Use these for spin-resolved collision tables when only total $\mathrm{H_2}$ density is known.
/// The energy gap between $J=1$ and $J=0$ is $\Delta E/k_\mathrm{B}=170.6\,\mathrm{K}$.
/// Their statistical weights are $9$ and $1$,
/// including nuclear spin. The high-temperature ortho/para limit is $3$.
/// See <https://doi.org/10.1016/j.crhy.2005.03.002>.
///
/// # Errors
///
/// Returns an error for non-finite inputs, a negative density, or a
/// non-positive kinetic temperature.
pub fn thermal_h2_densities(total: f64, temperature_k: f64) -> Result<(f64, f64)> {
    ensure!(total.is_finite() && total >= 0.0, "invalid H2 density");
    check_temperature(temperature_k)?;
    let ratio = (9.0 * (-170.6 / temperature_k).exp()).min(3.0);
    let para = total / (1.0 + ratio);
    Ok((para, para * ratio))
}

fn check_temperature(temperature_k: f64) -> Result<()> {
    ensure!(
        temperature_k.is_finite() && temperature_k > 0.0,
        "invalid kinetic temperature"
    );
    Ok(())
}

/// Builds the collision matrix in $\mathrm{s^{-1}}$, indexed by `[initial_level, final_level]`.
/// Its diagonal is zero. A row sum gives the inverse collisional lifetime.
/// Downward coefficients interpolate linearly in temperature and hold their
/// endpoint values outside the grid. Upward rates follow detailed balance.
///
/// # Errors
///
/// Returns an error for non-finite inputs, a non-positive kinetic temperature,
/// an empty level list, or a non-positive statistical weight. Each collider
/// requires a non-negative density and a nonempty, strictly increasing grid
/// of positive temperatures. Each transition must connect distinct, valid
/// levels in descending or equal energy order and supply one non-negative
/// coefficient per grid point.
pub fn collision_rates(
    levels: &[Level],
    temperature_k: f64,
    colliders: &[Collider<'_>],
) -> Result<Array2<f64>> {
    check_temperature(temperature_k)?;
    ensure!(!levels.is_empty(), "empty level list");
    for level in levels {
        ensure!(
            level.wavenumber.is_finite() && level.weight.is_finite() && level.weight > 0.0,
            "invalid level energy or weight"
        );
    }

    let mut rates = Array2::<f64>::zeros((levels.len(), levels.len()));
    for collider in colliders {
        ensure!(
            collider.density.is_finite() && collider.density >= 0.0,
            "invalid collider density"
        );
        let grid = collider.temperatures_k;
        ensure!(
            !grid.is_empty()
                && grid.iter().all(|t| t.is_finite() && *t > 0.0)
                && grid.windows(2).all(|pair| pair[0] < pair[1]),
            "expected strictly increasing positive collision temperatures"
        );

        let hi = grid
            .partition_point(|&t| t < temperature_k)
            .min(grid.len() - 1);
        let lo = hi.saturating_sub(1);
        let fraction = if lo == hi {
            0.0
        } else {
            (temperature_k.min(grid[hi]) - grid[lo]) / (grid[hi] - grid[lo])
        };
        for transition in collider.transitions {
            let (u, l) = (transition.upper, transition.lower);
            ensure!(
                u < levels.len()
                    && l < levels.len()
                    && u != l
                    && levels[u].wavenumber >= levels[l].wavenumber,
                "invalid downward transition {u} -> {l}"
            );
            let k = transition.coefficients;
            ensure!(
                k.len() == grid.len() && k.iter().all(|x| x.is_finite() && *x >= 0.0),
                "invalid coefficients for transition {u} -> {l}"
            );
            let downward = collider.density * ((1.0 - fraction) * k[lo] + fraction * k[hi]);
            let gap_k = codata::SECOND_RADIATION_CONSTANT_CGS
                * (levels[u].wavenumber - levels[l].wavenumber);
            rates[[u, l]] += downward;
            rates[[l, u]] +=
                downward * (levels[u].weight / levels[l].weight) * (-gap_k / temperature_k).exp();
        }
    }
    Ok(rates)
}
