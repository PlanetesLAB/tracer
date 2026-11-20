//! Local escape-probability approximations for homogeneous media.

use std::f64::consts::PI;

/// Geometry assumed when estimating the fraction of line photons that escape.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EscapeGeometry {
    /// Homogeneous sphere (Osterbrock 1974, Appendix 2).
    UniformSphere,
    ExpandingSphere,
    /// Plane-parallel closure with effective depth 3τ (de Jong et al. 1975, ApJ 199, 69).
    PlaneParallelSlab,
}

impl EscapeGeometry {
    /// Evaluates the escape factor at dimensionless line-centre optical depth.
    /// For spheres, supply the depth across the diameter.
    ///
    /// The fits describe absorbing media. With signed depths they can return
    /// factors above unity, negative thick-limit values, or (for an expanding
    /// sphere at depth <= -14) NaN. No clamping is applied.
    #[must_use]
    pub fn probability(self, optical_depth: f64) -> f64 {
        let radius_depth = optical_depth / 2.0;
        match self {
            Self::UniformSphere if radius_depth.abs() < 0.1 => {
                1.0 - 0.75 * radius_depth + radius_depth.powi(2) / 2.5 - radius_depth.powi(3) / 6.0
                    + radius_depth.powi(4) / 17.5
            }
            Self::UniformSphere if radius_depth.abs() > 50.0 => 0.75 / radius_depth,
            Self::UniformSphere => {
                let correction = 0.5 / radius_depth.powi(2);
                0.75 / radius_depth
                    * (1.0 - correction
                        + (radius_depth.recip() + correction) * (-optical_depth).exp())
            }
            Self::ExpandingSphere if radius_depth.abs() < 0.01 => 1.0,
            Self::ExpandingSphere if radius_depth.abs() < 7.0 => {
                // De Jong, Boland & Dalgarno (1980) escape approximation.
                // Eq. A21 at https://doi.org/10.1088/0004-637X/693/2/1360 gives
                // the one-sided form. Doubling it gives a thin limit of unity.
                let scaled = 2.34 * radius_depth;
                -(-scaled).exp_m1() / scaled
            }
            Self::ExpandingSphere => 0.5 / (radius_depth * (radius_depth / PI.sqrt()).ln().sqrt()),
            Self::PlaneParallelSlab => {
                let scaled = 3.0 * optical_depth;
                if scaled.abs() < 0.1 {
                    1.0 - 1.5 * (optical_depth + optical_depth.powi(2))
                } else if scaled.abs() > 50.0 {
                    scaled.recip()
                } else {
                    -(-scaled).exp_m1() / scaled
                }
            }
        }
    }
}
