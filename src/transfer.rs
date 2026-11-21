//! Monochromatic radiative transfer along a ray.

use anyhow::Result;
use anyhow::ensure;

/// Mean transmission from uniformly distributed emission sites to the exit.
pub(crate) fn mean_transmission(tau: f64) -> f64 {
    if tau == 0.0 {
        1.0
    } else {
        -(-tau).exp_m1() / tau
    }
}

/// Advances specific intensity through one segment with constant coefficients.
/// Apply segments in propagation order, starting at the incident boundary.
/// Supply emissivity $j_\nu$ in $\mathrm{erg\,s^{-1}\,cm^{-3}\,Hz^{-1}\,sr^{-1}}$, signed net
/// absorption $\alpha_\nu$ in $\mathrm{cm^{-1}}$, and length in cm.
/// Incident and returned intensities are in $\mathrm{erg\,s^{-1}\,cm^{-2}\,Hz^{-1}\,sr^{-1}}$.
///
/// At zero opacity the intensity increases by $j_\nu L$. Negative opacity
/// gives unsaturated exponential amplification.
///
/// # Errors
///
/// Requires finite inputs, nonnegative incident intensity, emissivity and
/// length, and representable outgoing intensity.
pub fn constant_segment(
    incident: f64,
    emissivity: f64,
    absorption: f64,
    length: f64,
) -> Result<f64> {
    ensure!(
        absorption.is_finite()
            && [incident, emissivity, length]
                .iter()
                .all(|v| v.is_finite() && *v >= 0.0),
        "invalid constant-segment input"
    );
    if length == 0.0 || (incident == 0.0 && emissivity == 0.0) {
        return Ok(incident);
    }

    let tau = absorption * length;
    let transmitted = if incident == 0.0 {
        0.0
    } else if tau.abs() > 700.0 {
        // The transmission factor alone may overflow or underflow.
        (incident.ln() - tau).exp()
    } else {
        incident * (-tau).exp()
    };
    let emitted = if emissivity == 0.0 {
        0.0
    } else if tau < -700.0 {
        (emissivity.ln() - (-absorption).ln() - tau).exp()
    } else if tau > 1.0 {
        // In a thick cell j L may overflow even for a finite source function.
        let source = emissivity / absorption;
        let absorbed = -(-tau).exp_m1();
        if source.is_finite() {
            source * absorbed
        } else {
            (emissivity * absorbed) / absorption
        }
    } else {
        let unattenuated = emissivity * length;
        let mean = mean_transmission(tau);
        if unattenuated.is_normal() {
            unattenuated * mean
        } else {
            // Form the product in log space if j L alone lost its range.
            (emissivity.ln() + length.ln() + mean.ln()).exp()
        }
    };
    let outgoing = transmitted + emitted;
    ensure!(outgoing.is_finite(), "unrepresentable outgoing intensity");
    Ok(outgoing)
}
