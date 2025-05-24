use crate::constants as cc;

pub fn planck(t: f64, lam: f64) -> f64 {
    let nu = cc::SPEED_OF_LIGHT_CGS / (lam * 1e-4);
    let x = cc::PLANCK_CGS * nu / (cc::BOLTZMANN_CGS * t);
    if x > 40.0 {
        2.0 * cc::PLANCK_CGS * nu.powi(3) / (cc::SPEED_OF_LIGHT_CGS.powi(2)) * (-x).exp()
    } else if x < 0.1 {
        2.0 * cc::PLANCK_CGS * nu.powi(3) / (cc::SPEED_OF_LIGHT_CGS.powi(2))
            * (-0.5 + 1.0 / x + x / 12.0 - x.powi(3) / 720.0)
    } else {
        2.0 * cc::PLANCK_CGS * nu.powi(3) / (cc::SPEED_OF_LIGHT_CGS.powi(2)) / (x.exp() - 1.0)
    }
}
