//! Line observables along a ray through a uniform molecular column.

use std::f64::consts::PI;

use anyhow::Result;
use anyhow::ensure;
use codata::aliases::{C_CGS, H_CGS, K_CGS};

use crate::line::DopplerWidth;
use crate::line::{absorption_cross_section, emission_per_particle};
use crate::transfer::mean_transmission;

/// Homogeneous-source observables for one transition.
pub struct LineObservables {
    /// Excitation temperature in K. Equal populations per substate give infinity.
    /// Empty upper and lower levels have no defined excitation temperature.
    pub excitation_temperature: Option<f64>,
    /// Signed optical depth at line centre.
    pub optical_depth: f64,
    /// Background-subtracted line-centre specific intensity in
    /// $\mathrm{erg\,s^{-1}\,cm^{-2}\,Hz^{-1}\,sr^{-1}}$.
    pub excess_intensity: f64,
    /// Rayleigh-Jeans temperature of the excess intensity in K.
    pub radiation_temperature: f64,
    /// Gaussian velocity integral in $\mathrm{K\,km\,s^{-1}}$.
    pub integrated_temperature: f64,
    /// Isotropic-equivalent $4\pi\int\Delta I_\nu\,d\nu$, in $\mathrm{erg\,cm^{-2}\,s^{-1}}$.
    pub surface_flux: f64,
}

/// Evaluates one line using upper and lower level columns in $\mathrm{cm^{-2}}$.
/// Supply the rest frequency in Hz, spontaneous rate in $\mathrm{s^{-1}}$,
/// upper-to-lower statistical weight ratio, and incident specific intensity in
/// $\mathrm{erg\,s^{-1}\,cm^{-2}\,Hz^{-1}\,sr^{-1}}$.
/// Multiply the molecular column by each level fraction to obtain its column.
///
/// The integrals assume a Gaussian brightness profile with the supplied Doppler
/// width and the calculated peak. Opacity broadening can change the actual profile.
///
/// # Errors
///
/// Requires positive finite frequency and weight ratio, nonnegative finite
/// columns, spontaneous rate and background intensity, and representable
/// transfer results. Strong unsaturated amplification can overflow.
pub fn homogeneous_line(
    frequency_hz: f64,
    einstein_a: f64,
    weight_ratio: f64,
    upper_column: f64,
    lower_column: f64,
    width: DopplerWidth,
    background_intensity: f64,
) -> Result<LineObservables> {
    ensure!(
        frequency_hz.is_finite() && frequency_hz > 0.0,
        "invalid line frequency"
    );
    ensure!(
        weight_ratio.is_finite() && weight_ratio > 0.0,
        "invalid statistical weight ratio"
    );
    ensure!(
        [upper_column, lower_column, einstein_a, background_intensity]
            .iter()
            .all(|x| x.is_finite() && *x >= 0.0),
        "invalid column, spontaneous rate, or background intensity"
    );

    let wavelength = C_CGS / frequency_hz;
    let weighted_lower = weight_ratio * lower_column;
    let difference = weighted_lower - upper_column;
    let photon_temperature = H_CGS * frequency_hz / K_CGS;
    let excitation_temperature = if upper_column == 0.0 && lower_column == 0.0 {
        None
    } else {
        let log_ratio = if difference.abs() <= 0.5 * upper_column && upper_column > 0.0 {
            (difference / upper_column).ln_1p()
        } else {
            weight_ratio.ln() + lower_column.ln() - upper_column.ln()
        };
        Some(photon_temperature / log_ratio)
    };

    let cross_section = absorption_cross_section(wavelength, einstein_a, width);
    let optical_depth = cross_section * difference;
    let emission = emission_per_particle(einstein_a, width) * upper_column;
    // Integrate emissivity directly to keep the zero-opacity limit finite.
    let absorbed = -(-optical_depth).exp_m1();
    let excess_intensity = emission.mul_add(
        mean_transmission(optical_depth),
        -background_intensity * absorbed,
    );
    let radiation_temperature = excess_intensity * wavelength.powi(2) / (2.0 * K_CGS);
    let velocity_area = PI.sqrt() * width.b();
    let integrated_temperature = radiation_temperature * velocity_area / 1e5;
    let surface_flux = 4.0 * PI * excess_intensity * velocity_area / wavelength;
    ensure!(
        [
            optical_depth,
            excess_intensity,
            radiation_temperature,
            integrated_temperature,
            surface_flux
        ]
        .iter()
        .all(|x| x.is_finite()),
        "unrepresentable homogeneous line observable"
    );
    Ok(LineObservables {
        excitation_temperature,
        optical_depth,
        excess_intensity,
        radiation_temperature,
        integrated_temperature,
        surface_flux,
    })
}
