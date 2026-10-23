//! Points, directions and frames in space, and the knot insertion that
//! B-spline surfaces are clamped and cut into Bézier pieces with: the exact
//! geometry `surface3` is written in.

use std::ops::{Add, Mul, Neg, Sub};

#[derive(Clone, Copy, Debug, PartialEq, Default)]
pub struct V3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl V3 {
    pub const fn new(x: f64, y: f64, z: f64) -> V3 {
        V3 { x, y, z }
    }
    pub fn dot(self, o: V3) -> f64 {
        self.x * o.x + self.y * o.y + self.z * o.z
    }
    pub fn cross(self, o: V3) -> V3 {
        V3::new(
            self.y * o.z - self.z * o.y,
            self.z * o.x - self.x * o.z,
            self.x * o.y - self.y * o.x,
        )
    }
    pub fn len(self) -> f64 {
        self.dot(self).sqrt()
    }
    pub fn scale(self, s: f64) -> V3 {
        V3::new(self.x * s, self.y * s, self.z * s)
    }
    pub fn unit(self) -> V3 {
        let l = self.len();
        if l > 0.0 {
            self.scale(1.0 / l)
        } else {
            self
        }
    }
}

impl Add for V3 {
    type Output = V3;
    fn add(self, o: V3) -> V3 {
        V3::new(self.x + o.x, self.y + o.y, self.z + o.z)
    }
}
impl Sub for V3 {
    type Output = V3;
    fn sub(self, o: V3) -> V3 {
        V3::new(self.x - o.x, self.y - o.y, self.z - o.z)
    }
}
impl Mul<f64> for V3 {
    type Output = V3;
    fn mul(self, s: f64) -> V3 {
        V3::new(self.x * s, self.y * s, self.z * s)
    }
}
impl Neg for V3 {
    type Output = V3;
    fn neg(self) -> V3 {
        V3::new(-self.x, -self.y, -self.z)
    }
}

/// Right-handed frame: origin, z axis, x axis (y = z cross x).
#[derive(Clone, Copy, Debug)]
pub struct Frame {
    pub o: V3,
    pub z: V3,
    pub x: V3,
}

impl Frame {
    pub fn y(&self) -> V3 {
        self.z.cross(self.x)
    }
}

impl Frame {
    /// Coordinates of p in this frame.
    pub fn local(&self, p: V3) -> V3 {
        let d = p - self.o;
        V3::new(d.dot(self.x), d.dot(self.y()), d.dot(self.z))
    }
    /// World point from frame coordinates.
    pub fn world(&self, l: V3) -> V3 {
        self.o + self.x * l.x + self.y() * l.y + self.z * l.z
    }
    /// World direction from frame coordinates.
    pub fn dir(&self, l: V3) -> V3 {
        self.x * l.x + self.y() * l.y + self.z * l.z
    }
}

/// Inserts knot u once into a (local piece of a) B-spline: The NURBS Book
/// algorithm 5.1 for one insertion, in homogeneous coordinates. Point i has
/// support `knots[i] .. knots[i + p + 1]`.
pub fn insert_knot<const N: usize>(
    p: usize,
    knots: &mut Vec<f64>,
    pts: &mut Vec<[f64; N]>,
    u: f64,
) {
    let k = (0..knots.len() - 1)
        .rev()
        .find(|&i| knots[i] <= u && knots[i + 1] > u)
        .expect("knot inside the vector");
    let s = knots[..=k].iter().rev().take_while(|&&x| x == u).count();
    let mut q: Vec<[f64; N]> = Vec::with_capacity(pts.len() + 1);
    q.extend_from_slice(&pts[..=k - p]);
    for i in (k - p + 1)..=(k - s) {
        let a = (u - knots[i]) / (knots[i + p] - knots[i]);
        let mut r = [0.0; N];
        for c in 0..N {
            r[c] = a * pts[i][c] + (1.0 - a) * pts[i - 1][c];
        }
        q.push(r);
    }
    q.extend_from_slice(&pts[(k - s)..]);
    knots.insert(k + 1, u);
    *pts = q;
}

/// The same B-spline with clamped ends (first and last knots repeated degree + 1
/// times), exactly: insert each end of the domain `knots[p] .. knots[m - p]` until
/// it has multiplicity p, then drop the knots and points outside the domain.
/// A clamped input comes back unchanged.
pub fn clamp<const N: usize>(
    p: usize,
    knots: &[f64],
    pts: &[[f64; N]],
) -> (Vec<f64>, Vec<[f64; N]>) {
    let (mut kn, mut pt) = (knots.to_vec(), pts.to_vec());
    let clamped = |k: &[f64]| {
        k[..=p].iter().all(|&x| x == k[0])
            && k[k.len() - p - 1..].iter().all(|&x| x == k[k.len() - 1])
    };
    if clamped(&kn) {
        return (kn, pt);
    }
    // Start: raise knots[p] to multiplicity p, keep the points from last_a - p on.
    let a = kn[p];
    while kn.iter().filter(|&&x| x == a).count() < p {
        insert_knot(p, &mut kn, &mut pt, a);
    }
    let last_a = kn.iter().rposition(|&x| x == a).unwrap();
    let mut k2 = vec![a; p + 1];
    k2.extend_from_slice(&kn[last_a + 1..]);
    let (mut kn, mut pt) = (k2, pt[last_a - p..].to_vec());
    // End: the mirror image.
    let m = kn.len() - 1;
    let b = kn[m - p];
    while kn.iter().filter(|&&x| x == b).count() < p {
        insert_knot(p, &mut kn, &mut pt, b);
    }
    let first_b = kn.iter().position(|&x| x == b).unwrap();
    let mut k3 = kn[..first_b].to_vec();
    k3.extend(std::iter::repeat_n(b, p + 1));
    let n = k3.len() - p - 1;
    pt.truncate(n);
    (k3, pt)
}
