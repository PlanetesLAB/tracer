//! Local statistical equilibrium with line trapping in a homogeneous medium.

use anyhow::Result;
use anyhow::{bail, ensure};
use planetes_ext::types::{RMatrix, RVector};

use crate::collision::Level;
use crate::equilibrium::RadiativeTransition;
use crate::equilibrium::optically_thin_populations;
use crate::escape::EscapeGeometry;
use crate::line::DopplerWidth;
use crate::line::absorption_cross_section;

/// Converged fractional populations and line-centre depths, in input order.
pub struct TrappedEquilibrium {
    pub populations: RVector,
    pub optical_depths: RVector,
    /// Number of feedback solves after optically thin initialization.
    pub iterations: usize,
}

/// Solves local equilibrium with escape-modified spontaneous and stimulated rates.
/// Supply the molecular column in $\mathrm{cm^{-2}}$ and a common Gaussian width.
/// For a sphere the column spans its diameter. Each transition uses the frequency
/// implied by its level-energy difference.
///
/// The iteration retains 70% of the old population and averages successive
/// excitation temperatures. Convergence requires a population fixed point and a
/// mean relative excitation-temperature change below $10^{-6}$ for absorbing
/// lines with $\tau>0.01$. At a singular excitation temperature the
/// population test still applies.
///
/// # Errors
///
/// Returns an error for invalid equilibrium inputs, a negative or non-finite
/// column, an unsupported escape factor, or failure to converge in 9999 updates.
pub fn trapped_populations(
    levels: &[Level],
    collisions: &RMatrix,
    lines: &[RadiativeTransition],
    column: f64,
    width: DopplerWidth,
    geometry: EscapeGeometry,
) -> Result<TrappedEquilibrium> {
    ensure!(
        column.is_finite() && column >= 0.0,
        "invalid molecular column"
    );
    let mut populations = optically_thin_populations(levels, collisions, lines)?;
    let mut scales = Vec::with_capacity(lines.len());
    let mut temperatures = Vec::with_capacity(lines.len());
    for line in lines {
        let upper = &levels[line.upper];
        let lower = &levels[line.lower];
        let gap = upper.wavenumber - lower.wavenumber;
        // Opacity per unit population difference, integrated over the column.
        let opacity = column * absorption_cross_section(gap.recip(), line.einstein_a, width);
        let ratio = upper.weight / lower.weight;
        let energy_k = gap * codata::SECOND_RADIATION_CONSTANT_CGS;
        scales.push((opacity, ratio, energy_k));
        temperatures.push(excitation_temperature(line, &populations, ratio, energy_k));
    }
    let depths = |p: &RVector| {
        RVector::from_iter(
            lines
                .iter()
                .zip(&scales)
                .map(|(line, &(opacity, ratio, _))| {
                    opacity * (ratio * p[line.lower] - p[line.upper])
                }),
        )
    };
    let mut optical_depths = depths(&populations);
    let mut effective = lines.to_vec();
    for iteration in 1..=9999 {
        for ((line, adjusted), &tau) in lines.iter().zip(&mut effective).zip(&optical_depths) {
            let escape = geometry.probability(tau);
            ensure!(
                tau.is_finite() && escape.is_finite() && escape > 0.0,
                "invalid escape factor for transition {} -> {} at optical depth {tau}",
                line.upper,
                line.lower
            );
            adjusted.einstein_a = line.einstein_a * escape;
        }
        let candidate = optically_thin_populations(levels, collisions, &effective)?;
        let populations_stable = populations
            .iter()
            .zip(&candidate)
            .all(|(&old, &new)| (old - new).abs() <= 1e-6 * old.max(new).max(1e-20));
        let mut thick_count = 0_u32;
        let mut temperature_change = 0.0;
        for (i, line) in lines.iter().enumerate() {
            let (_, ratio, energy_k) = scales[i];
            let new = excitation_temperature(line, &candidate, ratio, energy_k);
            if let (Some(old), Some(current)) = (temperatures[i], new) {
                if optical_depths[i] > 0.01 {
                    thick_count += 1;
                    temperature_change += ((current - old) / current).abs();
                }
                temperatures[i] = Some(0.5 * old + 0.5 * current);
            } else {
                temperatures[i] = new;
            }
        }
        populations.zip_mut_with(&candidate, |old, &new| *old = 0.7 * *old + 0.3 * new);
        optical_depths = depths(&populations);
        if iteration >= 10
            && populations_stable
            && (thick_count == 0 || temperature_change < 1e-6 * f64::from(thick_count))
        {
            return Ok(TrappedEquilibrium {
                populations,
                optical_depths,
                iterations: iteration,
            });
        }
    }
    bail!("escape-probability iteration failed to converge in 9999 updates")
}

fn excitation_temperature(
    line: &RadiativeTransition,
    populations: &RVector,
    ratio: f64,
    energy_k: f64,
) -> Option<f64> {
    let upper = populations[line.upper];
    let lower = populations[line.lower];
    if upper == 0.0 || lower == 0.0 {
        return None;
    }
    let temperature = energy_k / (ratio.ln() + lower.ln() - upper.ln());
    temperature.is_finite().then_some(temperature)
}
