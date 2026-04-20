use std::{cmp::Ordering, process::exit};

pub fn triangulate(mut points: Vec<Point>) {
    let tol = 100.0 * f64::EPSILON;
    let point_num = points.len();

    points.sort();

    for pair in points.windows(2) {
        let [p0, p1] = pair else { continue };

        assert_far_enough(p0.x, p1.x, tol);
        assert_far_enough(p0.y, p1.y, tol);
    }

    let mut m1 = 0;
    let mut m2 = 1;
    let mut j = 3;
    let mut m = j;

    let lr = loop {
        if point_num < j {
            exit(225);
        }

        m = j;
        let p = Point {
            x: points[m].x,
            y: points[m].y,
        };
        let line = Line {
            start: Point {
                x: points[m1].x,
                y: points[m1].y,
            },
            end: Point {
                x: points[m2].x,
                y: points[m2].y,
            },
        };
        let lr = line.relative_to(&p, 0.0);

        if lr != Position::On {
            break lr;
        }
        j += 1;
    };

    let tri_num = j - 2;

    let mut triangles = Vec::with_capacity(tri_num);
    // let mut triangle_neighbours = Vec::new();

    match lr {
        Position::Left => {
            triangles.push(Triangle {
                a: points[m1],
                b: points[m2],
                c: points[m],
            });
            // triangle_neighbours.push(todo!());
            // for i in 1..
        }
        _ => {}
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Point {
    pub x: f64,
    pub y: f64,
}

impl Eq for Point {}

impl PartialEq for Point {
    fn eq(&self, other: &Self) -> bool {
        self.x.total_cmp(&other.x) == Ordering::Equal
            && self.y.total_cmp(&other.y) == Ordering::Equal
    }
}

impl Ord for Point {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.x.total_cmp(&other.x) {
            Ordering::Equal => self.y.total_cmp(&other.y),
            ord => ord,
        }
    }
}

impl PartialOrd for Point {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other)) // <-- satisfies clippy
    }
}

pub struct Line {
    pub start: Point,
    pub end: Point,
}

pub struct Triangle {
    pub a: Point,
    pub b: Point,
    pub c: Point,
}

#[derive(Debug, PartialEq)]
pub enum Position {
    Left,
    Right,
    On,
}

impl Line {
    #[must_use]
    pub fn relative_to(&self, p: &Point, dv: f64) -> Position {
        let tol = 100.0 * f64::EPSILON;

        let dx = self.end.x - self.start.x;
        let dy = self.end.y - self.start.y;
        let dxu = p.x - self.start.x;
        let dyu = p.y - self.start.y;

        let tol_abs = tol
            * dx.abs()
                .max(dy.abs())
                .max(dxu.abs())
                .max(dyu.abs())
                .max(dv.abs());
        let t = dy * dxu - dx * dyu + dv * (dx * dx + dy * dy).sqrt();
        if tol_abs < t {
            Position::Right
        } else if -tol_abs <= t {
            Position::On
        } else {
            Position::Left
        }
    }
}

fn assert_far_enough(a: f64, b: f64, tol: f64) {
    let cmax = a.abs().max(b.abs());
    assert!(tol * (cmax + 1.0) < (a - b).abs());
}
