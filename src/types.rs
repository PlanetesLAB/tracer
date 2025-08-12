use crate::flits::config::{Cell, PopulationLevels};

use ndarray::ArrayView1;
use ndarray::{Array1, Array2, Array3, Array4};

/// 1-dimensional real vector.
pub type RVector = Array1<f64>;

/// 1-dimensional real vector view.
pub type RVecView<'a> = ArrayView1<'a, f64>;

/// A real matrix (2-dimensional real array).
pub type RMatrix = Array2<f64>;

/// [`crate::flits::config::Cell`] matrix
pub type CellMatrix = Array2<Cell>;

/// [`crate::flits::config::PopulationLevels`] vector
pub type PLVector = Array1<PopulationLevels>;

/// A real tensor (3-dimensional real array).
pub type RTensor = Array3<f64>;

/// A real 4-tensor (4-dimensional real array).
pub type RTensor4 = Array4<f64>;

/// A 1-dimensional unsigned-integer vector.
pub type UVector = Array1<usize>;

/// A 1-dimensional boolean vector.
pub type BVector = Array1<bool>;
