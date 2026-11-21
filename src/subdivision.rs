//! Uniform spatial sampling of a ray segment with a changing projected velocity.

use anyhow::Result;
use anyhow::ensure;

use crate::line::DopplerWidth;

/// Divides a ray segment into equal lengths using the projected velocity span,
/// the narrowest supplied Doppler width, and the velocity resolution of the
/// profile. All velocities are in $\mathrm{cm\,s^{-1}}$ and the length is in cm.
///
/// Returns pairs of (midpoint fraction, subsegment length) in propagation order.
/// Evaluate the local velocity field at each midpoint. The endpoint velocities
/// determine the count and need not describe a linear velocity field.
///
/// With $D=|v_{\rm exit}-v_{\rm entrance}|$, the width criterion gives
/// $\max(1,\lfloor 5D/b_{\min}\rfloor)$ samples, capped at
/// $1+\lfloor D/\Delta v_{\rm profile}\rfloor$. Equal endpoint velocities give
/// one sample, including for a field with an unresolved interior extremum.
///
/// # Errors
///
/// Requires finite endpoint velocities, nonnegative finite length, positive
/// finite profile resolution, at least one width, and a count that fits in a
/// 32-bit unsigned integer.
pub fn uniform_velocity_subdivision(
    length: f64,
    entrance_velocity: f64,
    exit_velocity: f64,
    widths: &[DopplerWidth],
    profile_resolution: f64,
) -> Result<impl Iterator<Item = (f64, f64)> + use<>> {
    ensure!(
        length.is_finite() && length >= 0.0,
        "invalid segment length"
    );
    ensure!(
        entrance_velocity.is_finite() && exit_velocity.is_finite(),
        "invalid endpoint velocity"
    );
    ensure!(
        profile_resolution.is_finite() && profile_resolution > 0.0,
        "invalid profile resolution"
    );
    ensure!(!widths.is_empty(), "missing Doppler width");
    let narrowest = widths.iter().map(|w| w.b()).fold(f64::INFINITY, f64::min);
    let span = (exit_velocity - entrance_velocity).abs();
    let scaled_span = span * 5.0;
    let width_steps = if scaled_span.is_finite() {
        scaled_span / narrowest
    } else {
        (span / narrowest) * 5.0
    };
    let requested = width_steps.floor().max(1.0);
    let cap = (span / profile_resolution).floor() + 1.0;
    let count = requested.min(cap);
    ensure!(
        count <= f64::from(u32::MAX),
        "velocity subdivision count exceeds supported range"
    );
    let step = length / count;
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let bins = count as u32;
    Ok((0..bins).map(move |i| ((f64::from(i) + 0.5) / count, step)))
}
