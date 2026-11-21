//! Monochromatic radiative transfer along a ray.

use anyhow::Result;
use anyhow::ensure;

use crate::line::LineCoefficients;

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

/// Treatment of the interpolated emission in a two-endpoint segment.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EmissionLimit {
    /// Cap emission at the trapezoidal integral of unattenuated emissivity.
    Unattenuated,
    /// Use the linear source-function integral directly.
    Unbounded,
}

/// Advances intensity from the entrance to the exit of a segment using the
/// endpoint coefficients. Units follow [`constant_segment`]. The optical
/// depth is the trapezoidal integral of absorption. The source function
/// varies linearly with optical depth. When only one endpoint has positive
/// opacity, its source function supplies both endpoints.
///
/// For $\tau\le10^{-9}$ the emission is the trapezoidal emissivity integral.
/// For $\tau\le10^{-6}$ the incident transmission is $1-\tau$. This branch
/// approximates weak inversions and gives only first-order gain at large
/// negative depth. [`constant_segment`] gives exponential gain in a uniform cell.
///
/// # Errors
///
/// Requires finite inputs, nonnegative emissivities, incident intensity and
/// length, and representable optical depth and outgoing intensity.
pub fn two_endpoint_segment(
    incident: f64,
    entrance: &LineCoefficients,
    exit: &LineCoefficients,
    length: f64,
    limit: EmissionLimit,
) -> Result<f64> {
    ensure!(
        [incident, entrance.emissivity, exit.emissivity, length]
            .iter()
            .all(|v| v.is_finite() && *v >= 0.0)
            && entrance.absorption.is_finite()
            && exit.absorption.is_finite(),
        "invalid two-endpoint segment input"
    );
    if length == 0.0 {
        return Ok(incident);
    }
    let tau = entrance.absorption.midpoint(exit.absorption) * length;
    ensure!(tau.is_finite(), "unrepresentable segment optical depth");
    let unattenuated = entrance.emissivity.midpoint(exit.emissivity) * length;
    let transmission = if tau <= 1e-6 { 1.0 - tau } else { (-tau).exp() };
    let emitted = if tau <= 1e-9 {
        unattenuated
    } else {
        let (upstream, downstream) = if tau <= 1e-6 {
            (0.5 * tau, 0.5 * tau)
        } else {
            let mean = mean_transmission(tau);
            (mean - transmission, 1.0 - mean)
        };
        // Positive total depth guarantees at least one positive endpoint.
        let start = if entrance.absorption > 0.0 {
            entrance
        } else {
            exit
        };
        let end = if exit.absorption > 0.0 {
            exit
        } else {
            entrance
        };
        let emission = weighted_source(start, upstream) + weighted_source(end, downstream);
        match limit {
            EmissionLimit::Unattenuated => emission.min(unattenuated),
            EmissionLimit::Unbounded => emission,
        }
    };
    let outgoing = if tau > 700.0 && incident > 0.0 {
        (incident.ln() - tau).exp() + emitted
    } else {
        incident.mul_add(transmission, emitted)
    };
    ensure!(outgoing.is_finite(), "unrepresentable outgoing intensity");
    Ok(outgoing)
}

fn weighted_source(endpoint: &LineCoefficients, weight: f64) -> f64 {
    if weight == 0.0 || endpoint.emissivity == 0.0 {
        return 0.0;
    }
    let source = endpoint.emissivity / endpoint.absorption;
    if source.is_finite() {
        weight * source
    } else {
        (weight.ln() + endpoint.emissivity.ln() - endpoint.absorption.ln()).exp()
    }
}
