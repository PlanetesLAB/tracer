//! Doppler-shifted emission and absorption of an isolated Gaussian line.

use anyhow::Result;
use anyhow::ensure;

use crate::line::{DopplerWidth, LineCoefficients};
use crate::transfer::constant_segment;

/// Peak-normalized profile evaluated directly or sampled on a uniform grid.
#[derive(Default)]
pub enum GaussianProfile {
    #[default]
    Analytic,
    Tabulated(GaussianTable),
}

/// One symmetric half of a sampled Gaussian in units of the Doppler width.
pub struct GaussianTable {
    spacing: f64,
    samples: Vec<f64>,
}

impl GaussianProfile {
    /// Builds a lookup table with dimensionless spacing and an exclusive bin limit.
    /// For $x=(v_{\rm channel}-v_{\rm gas})/b$, use bin $\lfloor |x|/\Delta x\rfloor$.
    /// Bins beyond the limit and samples with $(k\Delta x)^2\ge10$ contribute zero.
    ///
    /// # Errors
    ///
    /// Requires positive finite spacing, at least one bin, and sufficient storage.
    pub fn tabulated(spacing: f64, bins: u32) -> Result<Self> {
        ensure!(
            spacing.is_finite() && spacing > 0.0 && bins > 0,
            "invalid Gaussian table grid"
        );
        let mut samples = Vec::new();
        for k in 0..bins {
            let squared = (f64::from(k) * spacing).powi(2);
            if squared >= 10.0 {
                break;
            }
            samples.try_reserve(1)?;
            samples.push((-squared).exp());
        }
        Ok(Self::Tabulated(GaussianTable { spacing, samples }))
    }

    /// Samples the line coefficients at one observer velocity. Both velocities
    /// are in $\mathrm{cm\,s^{-1}}$, positive toward the observer. For a
    /// frequency channel, $v_{\rm channel}=c(\nu-\nu_0)/\nu_0$ to first order.
    /// The input coefficients already include the line-centre normalization.
    ///
    /// # Errors
    ///
    /// Requires finite velocities and coefficients and nonnegative emissivity.
    pub fn coefficients(
        &self,
        centre: &LineCoefficients,
        width: DopplerWidth,
        channel_velocity: f64,
        projected_gas_velocity: f64,
    ) -> Result<LineCoefficients> {
        ensure!(
            channel_velocity.is_finite() && projected_gas_velocity.is_finite(),
            "invalid line velocity"
        );
        ensure!(
            centre.emissivity.is_finite()
                && centre.emissivity >= 0.0
                && centre.absorption.is_finite(),
            "invalid line-centre coefficients"
        );
        let offset = (channel_velocity - projected_gas_velocity) / width.b();
        let factor = match self {
            Self::Analytic => (-offset.powi(2)).exp(),
            Self::Tabulated(table) => table.sample(offset),
        };
        Ok(LineCoefficients {
            emissivity: factor * centre.emissivity,
            absorption: factor * centre.absorption,
        })
    }
}

impl GaussianTable {
    fn sample(&self, offset: f64) -> f64 {
        // A nonnegative float-to-index cast truncates toward line centre.
        // Out-of-range offsets saturate and the bounds check returns zero.
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let bin = (offset.abs() / self.spacing) as usize;
        self.samples.get(bin).copied().unwrap_or_default()
    }
}

/// Uniform segment of a ray through gas carrying a single line.
pub struct MovingLineSegment {
    pub centre: LineCoefficients,
    pub width: DopplerWidth,
    /// Projected gas velocity toward the observer in $\mathrm{cm\,s^{-1}}$.
    pub projected_velocity: f64,
    /// Path length through the segment in cm.
    pub length: f64,
}

/// Traces one isolated line at a single observer velocity, from the incident
/// boundary through the segments in propagation order. The returned specific
/// intensity has the same units as the incident intensity.
///
/// # Errors
///
/// Requires finite channel velocity and nonnegative finite incident intensity.
/// Propagates errors from profile sampling and segment transfer.
pub fn trace_isolated_line(
    incident: f64,
    channel_velocity: f64,
    segments: &[MovingLineSegment],
    profile: &GaussianProfile,
) -> Result<f64> {
    ensure!(
        channel_velocity.is_finite() && incident.is_finite() && incident >= 0.0,
        "invalid channel velocity or incident intensity"
    );
    segments.iter().try_fold(incident, |intensity, segment| {
        let line = profile.coefficients(
            &segment.centre,
            segment.width,
            channel_velocity,
            segment.projected_velocity,
        )?;
        constant_segment(intensity, line.emissivity, line.absorption, segment.length)
    })
}
