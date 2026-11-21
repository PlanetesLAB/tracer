//! Local molecular excitation and Gaussian line-centre coefficients.

use std::f64::consts::{LN_2, PI};

use anyhow::Result;
use anyhow::ensure;
use codata::aliases::{C_CGS, H_CGS};
use planetes_ext::physics::doppler_b;

use crate::collision::Level;
use crate::collision::{check_levels, check_temperature};

/// Fractional LTE populations, in input level order, at the supplied temperature.
/// The energy zero is arbitrary. The partition sum runs over the supplied levels.
///
/// # Errors
///
/// Requires a positive finite temperature, a nonempty level list, finite
/// wavenumbers, and positive finite statistical weights.
pub fn lte_populations(levels: &[Level], temperature_k: f64) -> Result<Vec<f64>> {
    check_temperature(temperature_k)?;
    check_levels(levels)?;
    let ground = levels
        .iter()
        .map(|level| level.wavenumber)
        .fold(f64::INFINITY, f64::min);
    let mut populations: Vec<_> = levels
        .iter()
        .map(|level| {
            let gap = level.wavenumber - ground;
            level.weight.ln() - gap * codata::SECOND_RADIATION_CONSTANT_CGS / temperature_k
        })
        .collect();
    let peak = populations
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    for log_weight in &mut populations {
        *log_weight = (*log_weight - peak).exp();
    }
    let sum: f64 = populations.iter().sum();
    for fraction in &mut populations {
        *fraction /= sum;
    }
    Ok(populations)
}

/// Gaussian velocity width $b$ for a profile proportional to $\exp[-(v/b)^2]$.
#[derive(Clone, Copy, Debug)]
pub struct DopplerWidth(f64);

impl DopplerWidth {
    /// Constructs a width from $b$ in $\mathrm{cm\,s^{-1}}$.
    ///
    /// # Errors
    ///
    /// Requires a positive finite width.
    pub fn from_b(b: f64) -> Result<Self> {
        ensure!(b.is_finite() && b > 0.0, "invalid Doppler width");
        Ok(Self(b))
    }

    /// Constructs a width from the velocity FWHM in $\mathrm{cm\,s^{-1}}$.
    ///
    /// # Errors
    ///
    /// Requires a positive finite Doppler width after conversion.
    pub fn from_fwhm(fwhm: f64) -> Result<Self> {
        Self::from_b(fwhm / (2.0 * LN_2.sqrt()))
    }

    /// Combines thermal and turbulent broadening, $b^2=2k_\mathrm{B}T/m+b_\mathrm{turb}^2$.
    /// Supply mass in grams and turbulent Doppler width in $\mathrm{cm\,s^{-1}}$.
    ///
    /// # Errors
    ///
    /// Requires positive finite temperature and mass, nonnegative finite turbulent
    /// width, and a representable positive result.
    pub fn from_temperature(temperature_k: f64, mass_g: f64, turbulent_b: f64) -> Result<Self> {
        check_temperature(temperature_k)?;
        ensure!(mass_g.is_finite() && mass_g > 0.0, "invalid particle mass");
        ensure!(
            turbulent_b.is_finite() && turbulent_b >= 0.0,
            "invalid turbulent width"
        );
        let thermal = doppler_b(temperature_k, mass_g);
        Self::from_b(thermal.hypot(turbulent_b))
    }

    /// Doppler width in $\mathrm{cm\,s^{-1}}$.
    #[must_use]
    pub fn b(self) -> f64 {
        self.0
    }

    /// Full width at half maximum in $\mathrm{cm\,s^{-1}}$.
    #[must_use]
    pub fn fwhm(self) -> f64 {
        self.0 * (2.0 * LN_2.sqrt())
    }
}

/// Frequency-specific emission and net absorption at line centre.
pub struct LineCoefficients {
    /// $j_{\nu_0}$ in $\mathrm{erg\,s^{-1}\,cm^{-3}\,Hz^{-1}\,sr^{-1}}$.
    pub emissivity: f64,
    /// $\alpha_{\nu_0}$ in $\mathrm{cm^{-1}}$, including stimulated emission.
    pub absorption: f64,
}

/// Line-centre spontaneous emissivity per upper-level number density.
pub(crate) fn emission_per_particle(einstein_a: f64, width: DopplerWidth) -> f64 {
    H_CGS * C_CGS * einstein_a / (4.0 * PI * PI.sqrt() * width.b())
}

/// Line-centre absorption cross section per unit population difference,
/// $A_{ul}\lambda^3/(8\pi^{3/2}b)$ in $\mathrm{cm^2}$.
pub(crate) fn absorption_cross_section(
    wavelength_cm: f64,
    einstein_a: f64,
    width: DopplerWidth,
) -> f64 {
    einstein_a * wavelength_cm.powi(3) / (8.0 * PI * PI.sqrt() * width.b())
}

/// Evaluates the Gaussian line centre for arbitrary level populations.
/// Supply the rest frequency in Hz, spontaneous rate in $\mathrm{s^{-1}}$,
/// upper-to-lower statistical weight ratio, and level number densities in
/// $\mathrm{cm^{-3}}$. An inverted transition has negative net absorption.
///
/// # Errors
///
/// Requires positive finite frequency and weight ratio, and nonnegative finite
/// spontaneous rate and number densities.
pub fn line_centre(
    frequency_hz: f64,
    einstein_a: f64,
    weight_ratio: f64,
    upper_density: f64,
    lower_density: f64,
    width: DopplerWidth,
) -> Result<LineCoefficients> {
    ensure!(
        frequency_hz.is_finite() && frequency_hz > 0.0,
        "invalid line frequency"
    );
    ensure!(
        weight_ratio.is_finite() && weight_ratio > 0.0,
        "invalid statistical weight ratio"
    );
    ensure!(
        [einstein_a, upper_density, lower_density]
            .iter()
            .all(|x| x.is_finite() && *x >= 0.0),
        "invalid spontaneous rate or level density"
    );
    let wavelength = C_CGS / frequency_hz;
    Ok(LineCoefficients {
        emissivity: emission_per_particle(einstein_a, width) * upper_density,
        absorption: absorption_cross_section(wavelength, einstein_a, width)
            * (weight_ratio * lower_density - upper_density),
    })
}
