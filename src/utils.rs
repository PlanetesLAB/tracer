use ndarray::{ArrayBase, Data, Dimension};
use std::cmp::PartialOrd;

pub trait ArrayExtrema<T> {
    fn maxval(&self) -> Option<T>;
    fn minval(&self) -> Option<T>;
}

impl<T, S, D> ArrayExtrema<T> for ArrayBase<S, D>
where
    T: PartialOrd + Copy,
    S: Data<Elem = T>,
    D: Dimension,
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
}
