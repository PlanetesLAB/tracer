use std::cmp::PartialOrd;

use ndarray::{ArrayBase, Axis, OwnedRepr};
use ndarray::{Data, Dimension, RemoveAxis};

pub trait ArrayExtrema<T, D>
where
    D: Dimension,
{
    fn maxval(&self) -> Option<T>;
    fn minval(&self) -> Option<T>;

    fn maxval_along(&self, axis: Axis) -> Option<ArrayBase<OwnedRepr<T>, D::Smaller>>;
    fn minval_along(&self, axis: Axis) -> Option<ArrayBase<OwnedRepr<T>, D::Smaller>>;
}

impl<T, S, D> ArrayExtrema<T, D> for ArrayBase<S, D>
where
    T: PartialOrd + Copy,
    S: Data<Elem = T>,
    D: Dimension + RemoveAxis,
{
    fn maxval(&self) -> Option<T> {
        self.iter()
            .copied()
            .reduce(|a, b| if a > b { a } else { b })
    }

    fn minval(&self) -> Option<T> {
        self.iter()
            .copied()
            .reduce(|a, b| if a < b { a } else { b })
    }

    fn maxval_along(&self, axis: Axis) -> Option<ArrayBase<OwnedRepr<T>, D::Smaller>> {
        if self.is_empty() {
            None
        } else {
            Some(self.map_axis(axis, |subview| {
                subview
                    .iter()
                    .copied()
                    .reduce(|a, b| if a > b { a } else { b })
                    .unwrap()
            }))
        }
    }

    fn minval_along(&self, axis: Axis) -> Option<ArrayBase<OwnedRepr<T>, D::Smaller>> {
        if self.is_empty() {
            None
        } else {
            Some(self.map_axis(axis, |subview| {
                subview
                    .iter()
                    .copied()
                    .reduce(|a, b| if a < b { a } else { b })
                    .unwrap()
            }))
        }
    }
}

pub mod physics {
    use crate::constants as cc;

    #[must_use]
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
}
