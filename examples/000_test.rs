use ndarray::array;
use tracer::types::RMatrix;

fn main() {
    let points = array![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
    let p = points.as_slice();
    dbg!(p);
}
