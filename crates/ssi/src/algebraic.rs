//! Implicit forms for B-spline patches: a polynomial F(x, y, z) whose zero set
//! holds the patch exactly, so that a NURBS cylinder, sphere or torus, or a
//! constructed free-form sheet of low implicit degree, gets a field like an
//! analytic surface's (`field::Field`).
//!
//! F is the null vector of the patch's composition matrix: column m holds the
//! Bernstein coefficients of monomial m composed with the patch's Bézier
//! pieces, built from their control points (no sampling), so F composed with
//! the patch is zero exactly when the matrix times F's coefficients is.
//! Degrees are tried lowest first, and one is taken only where the matrix's
//! least singular value is at the rounding (1e-13 of its largest): an exact
//! form, not a close one (a bicubic sheet's generic implicit degree is 18, but
//! the constructed sheets of the ssi-corpus have exact forms of degree 3 to
//! 10, their singular values falling 1e4-fold there and decaying smoothly
//! for a merely close fit). Then proved: composed with every piece, each
//! coefficient within 1e-11 of the terms' scale (`PROOF`). The fit proposes F; the
//! composition decides.

use crate::field::Bern;
use crate::patch::Patch;
use crate::piece::{norm, sub, Piece, V3};

/// An implicit polynomial of total degree `degree` in coordinates centered on
/// `center` and scaled by `length`: Σ coef[k] x^a y^b z^c over `monomials`.
#[derive(Clone, Debug)]
pub struct Algebraic {
    pub degree: usize,
    center: V3,
    length: f64,
    monomials: Vec<[usize; 3]>,
    coef: Vec<f64>,
}

/// The highest degree tried.
pub(crate) const MAX_DEGREE: usize = 10;

/// The least singular value of the composition matrix, over its largest, at
/// which a degree has an exact form: rounding, with the matrix's conditioning.
const EXACT: f64 = 1e-13;

/// The largest coefficient of F composed with the patch, over the terms'
/// scale, that proves the patch on F = 0. Exact forms from all the rows come
/// to 1e-15 or so, from a quarter of them to 1e-12 (1e-13 sent those back to
/// all the rows, half again a hard pair's time). The subset must keep four
/// rows a column: from two, a near form, not the exact one, passed this
/// bound, and missed a touch point told to the rounding (ssi-corpus
/// y23_fig1b_two_arcs_and_point).
const PROOF: f64 = 1e-11;

impl Algebraic {
    /// The patch's implicit form of lowest degree, when one up to `MAX_DEGREE` is proved.
    pub fn of(patch: &Patch, index: usize) -> Option<Algebraic> {
        let mut fitter = Fitter::new(patch, index)?;
        while fitter.open() {
            if let Some(form) = fitter.next() {
                return Some(form);
            }
        }
        None
    }

    /// g = W^k F(X / W) on the piece, a Bernstein polynomial of degree k times
    /// the piece's, and the scale of its terms.
    pub fn compose(&self, piece: &Piece) -> (Bern, f64) {
        let terms = terms(piece, self.center, self.length, self.degree);
        let weighted: Vec<(f64, &Bern)> = self.coef.iter().copied().zip(&terms).collect();
        let m = (0..4)
            .map(|d| coordinate(piece, self.center, self.length, d).max_abs())
            .fold(0.0_f64, f64::max);
        (Bern::sum(&weighted), m.powi(self.degree as i32))
    }

    /// F, its gradient (in space, not the scaled coordinates), and the size of
    /// its terms at p, Σ |c x^a y^b z^c|, for its rounding.
    fn value(&self, p: V3) -> (f64, V3, f64) {
        let x = [
            (p[0] - self.center[0]) / self.length,
            (p[1] - self.center[1]) / self.length,
            (p[2] - self.center[2]) / self.length,
        ];
        let power = |v: f64, e: usize| if e == 0 { 1.0 } else { v.powi(e as i32) };
        let mut f = 0.0;
        let mut size = 0.0;
        let mut g = [0.0; 3];
        for (m, c) in self.monomials.iter().zip(&self.coef) {
            let t = [power(x[0], m[0]), power(x[1], m[1]), power(x[2], m[2])];
            f += c * t[0] * t[1] * t[2];
            size += (c * t[0] * t[1] * t[2]).abs();
            for k in 0..3 {
                if m[k] > 0 {
                    let d = m[k] as f64 * power(x[k], m[k] - 1);
                    let others: f64 = (0..3).filter(|&j| j != k).map(|j| t[j]).product();
                    g[k] += c * d * others / self.length;
                }
            }
        }
        (f, g, size)
    }

    /// The point of F = 0 that Newton reaches from p along F's gradient: p's foot
    /// to first order, for the gap at a contact. A point where F is zero to its
    /// rounding is its own foot: stepping on from there, where the gradient is
    /// small (a shell's pole), walked 0.08 away (ssi-corpus y23_fig1g_three_points).
    pub fn foot(&self, p: V3) -> V3 {
        let mut q = p;
        for _ in 0..50 {
            let (f, g, size) = self.value(q);
            let g2 = g[0] * g[0] + g[1] * g[1] + g[2] * g[2];
            if !(g2 > 0.0) || f.abs() <= 64.0 * f64::EPSILON * size {
                break;
            }
            let step = f / g2;
            q = [q[0] - step * g[0], q[1] - step * g[1], q[2] - step * g[2]];
            if (step * g2.sqrt()).abs() <= 1e-16 * self.length {
                break;
            }
        }
        q
    }
}

/// The search for a patch's implicit form a degree at a time, lowest first
/// (`Algebraic::of`; several patches in step, `intersect_bodies`).
pub(crate) struct Fitter {
    pieces: Vec<Piece>,
    center: V3,
    length: f64,
    /// A tensor patch of degree (m, n) has an implicit form of degree 2 m n at most.
    cap: usize,
    powers: Vec<Powers>,
    degree: usize,
}

impl Fitter {
    pub fn new(patch: &Patch, index: usize) -> Option<Fitter> {
        let pieces = patch.bezier_pieces(index).ok()?;
        let points: Vec<V3> = pieces
            .iter()
            .flat_map(|piece| piece.points().to_vec())
            .collect();
        let mut center = [0.0; 3];
        for p in &points {
            for k in 0..3 {
                center[k] += p[k] / points.len() as f64;
            }
        }
        let length = points
            .iter()
            .fold(0.0_f64, |m, p| m.max(norm(sub(*p, center))));
        if !(length > 0.0) {
            return None;
        }
        let cap = MAX_DEGREE.min(2 * pieces[0].du * pieces[0].dv).max(1);
        let powers = pieces
            .iter()
            .map(|piece| Powers::new(piece, center, length))
            .collect();
        Some(Fitter {
            pieces,
            center,
            length,
            cap,
            powers,
            degree: 0,
        })
    }

    /// Whether degrees remain to try.
    pub fn open(&self) -> bool {
        self.degree < self.cap
    }

    /// The next degree's proved form, when it has one.
    pub fn next(&mut self) -> Option<Algebraic> {
        self.degree += 1;
        let degree = self.degree;
        let monomials = monomials(degree);
        // Column per monomial: its composition's coefficients over every piece.
        let mut columns = vec![Vec::new(); monomials.len()];
        for p in &mut self.powers {
            p.raise();
            for (column, term) in columns.iter_mut().zip(p.terms(&monomials)) {
                column.extend_from_slice(&term.c);
            }
        }
        let full = &columns;
        let proved = |columns: &[Vec<f64>]| -> Result<Algebraic, bool> {
            let Some((coef, ratio)) = null_vector(columns) else {
                return Err(false);
            };
            if ratio > EXACT {
                return Err(false);
            }
            // The proof: F composed with every piece, from the compositions of
            // the monomials already made (each piece's coefficients a run of
            // the columns), not composed again.
            let mut start = 0;
            for (piece, p) in self.pieces.iter().zip(&self.powers) {
                let n = (degree * piece.du + 1) * (degree * piece.dv + 1);
                let scale = p
                    .base
                    .iter()
                    .fold(0.0_f64, |m, b| m.max(b.max_abs()))
                    .powi(degree as i32);
                let largest = (start..start + n)
                    .map(|r| {
                        full.iter()
                            .zip(&coef)
                            .map(|(column, c)| c * column[r])
                            .sum::<f64>()
                            .abs()
                    })
                    .fold(0.0_f64, f64::max);
                if largest > PROOF * scale {
                    return Err(true);
                }
                start += n;
            }
            Ok(Algebraic {
                degree,
                center: self.center,
                length: self.length,
                monomials: monomials.clone(),
                coef,
            })
        };
        // The null vector from every stride-th row, about four a column; where
        // that one fails the proof, from them all (a sphere's pole pieces
        // repeat their coefficients, and a strided subset of them lost rank:
        // its null vector was not the sphere's, and the sphere times a plane
        // was proved at the next degree, ssi-corpus c4n_viviani_nurbs).
        let stride = (columns[0].len() / (4 * monomials.len())).max(1);
        if stride > 1 {
            let some: Vec<Vec<f64>> = columns
                .iter()
                .map(|c| c.iter().step_by(stride).copied().collect())
                .collect();
            match proved(&some) {
                Ok(form) => return Some(form),
                Err(false) => return None,
                Err(true) => {}
            }
        }
        proved(&columns).ok()
    }
}

/// The piece's homogeneous coordinate d (x, y, z scaled about the center, then W).
fn coordinate(piece: &Piece, center: V3, length: f64, d: usize) -> Bern {
    Bern {
        du: piece.du,
        dv: piece.dv,
        c: piece
            .c
            .iter()
            .map(|h| {
                if d == 3 {
                    h[3]
                } else {
                    (h[d] - center[d] * h[3]) / length
                }
            })
            .collect(),
    }
}

/// Each monomial of `monomials(degree)`, x^a y^b z^c, composed with the piece
/// and homogenized: X^a Y^b Z^c W^(degree - a - b - c). Built a coordinate at
/// a time, each monomial of one degree from one of the degree below (a
/// product per monomial: from powers, degree 10 took three large products
/// each, most of a second a piece).
fn terms(piece: &Piece, center: V3, length: f64, degree: usize) -> Vec<Bern> {
    let mut powers = Powers::new(piece, center, length);
    for _ in 0..degree {
        powers.raise();
    }
    powers
        .terms(&monomials(degree))
        .into_iter()
        .cloned()
        .collect()
}

/// The homogeneous monomials of one degree composed with a piece, raised a
/// degree at a time (fitting tries every degree up to the form's, and from
/// the start each time took five times as long).
struct Powers {
    base: Vec<Bern>,
    degree: usize,
    /// Exponents [a, b, c, w], each raised in coordinate order so each is built once.
    level: std::collections::HashMap<[usize; 4], Bern>,
}

impl Powers {
    fn new(piece: &Piece, center: V3, length: f64) -> Powers {
        let base = (0..4)
            .map(|d| coordinate(piece, center, length, d))
            .collect();
        Powers {
            base,
            degree: 0,
            level: std::collections::HashMap::from([(
                [0; 4],
                Bern {
                    du: 0,
                    dv: 0,
                    c: vec![1.0],
                },
            )]),
        }
    }

    fn raise(&mut self) {
        let mut next = std::collections::HashMap::with_capacity(self.level.len() * 2);
        for (e, b) in &self.level {
            let last = (0..4).rev().find(|&i| e[i] > 0).unwrap_or(0);
            for (i, coordinate) in self.base.iter().enumerate().skip(last) {
                let mut raised = *e;
                raised[i] += 1;
                next.insert(raised, b.mul(coordinate));
            }
        }
        self.level = next;
        self.degree += 1;
    }

    /// Each of `monomials` (of this degree), x^a y^b z^c, homogenized.
    fn terms(&self, monomials: &[[usize; 3]]) -> Vec<&Bern> {
        let k = self.degree;
        monomials
            .iter()
            .map(|&[a, b, c]| &self.level[&[a, b, c, k - a - b - c]])
            .collect()
    }
}

/// The exponents [a, b, c] with a + b + c at most `degree`.
fn monomials(degree: usize) -> Vec<[usize; 3]> {
    let mut out = Vec::new();
    for total in 0..=degree {
        for a in (0..=total).rev() {
            for b in (0..=total - a).rev() {
                out.push([a, b, total - a - b]);
            }
        }
    }
    out
}

/// The unit null vector of the matrix with these columns, and its least
/// singular value over its largest; None with fewer rows than columns (a null
/// vector there says nothing). Householder QR, whose R has the matrix's
/// singular values, then inverse iteration with R for the least (two
/// triangular solves a step, a zero pivot raised to the rounding) and power
/// iteration for the largest: a few times the QR's m n^2, where a Jacobi SVD
/// of a 961 x 286 matrix took seconds.
fn null_vector(columns: &[Vec<f64>]) -> Option<(Vec<f64>, f64)> {
    let n = columns.len();
    let m = columns.first()?.len();
    if m < n {
        return None;
    }
    let mut a = columns.to_vec();
    for j in 0..n {
        let norm = a[j][j..].iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm == 0.0 {
            continue;
        }
        let alpha = if a[j][j] > 0.0 { -norm } else { norm };
        let mut v = a[j][j..].to_vec();
        v[0] -= alpha;
        let v2: f64 = v.iter().map(|x| x * x).sum();
        if v2 == 0.0 {
            continue;
        }
        for column in a.iter_mut().skip(j) {
            let f = 2.0 * v.iter().zip(&column[j..]).map(|(x, y)| x * y).sum::<f64>() / v2;
            for (x, y) in column[j..].iter_mut().zip(&v) {
                *x -= f * y;
            }
        }
    }
    // R, upper triangular: r[i][j] = a[j][i] for i <= j.
    let r = |i: usize, j: usize| a[j][i];
    let largest_pivot = (0..n).fold(0.0_f64, |p, i| p.max(r(i, i).abs()));
    if !(largest_pivot > 0.0) {
        return None;
    }
    let floor = 1e-300_f64.max(f64::EPSILON * 1e-3 * largest_pivot);
    let pivot = |i: usize| {
        let d = r(i, i);
        if d.abs() < floor {
            floor.copysign(if d == 0.0 { 1.0 } else { d })
        } else {
            d
        }
    };
    let times = |x: &[f64]| -> Vec<f64> {
        (0..n)
            .map(|i| (i..n).map(|j| r(i, j) * x[j]).sum())
            .collect()
    };
    let unit = |mut x: Vec<f64>| {
        let l = x.iter().map(|v| v * v).sum::<f64>().sqrt();
        x.iter_mut().for_each(|v| *v /= l);
        x
    };
    let mut x = unit((0..n).map(|i| 1.0 + 0.1 * (i as f64).sin()).collect());
    for _ in 0..6 {
        // R^T y = x, then R z = y.
        let mut y = vec![0.0; n];
        for i in 0..n {
            let s: f64 = (0..i).map(|k| r(k, i) * y[k]).sum();
            y[i] = (x[i] - s) / pivot(i);
        }
        let mut z = vec![0.0; n];
        for i in (0..n).rev() {
            let s: f64 = (i + 1..n).map(|k| r(i, k) * z[k]).sum();
            z[i] = (y[i] - s) / pivot(i);
        }
        if !z.iter().all(|v| v.is_finite()) {
            return None;
        }
        x = unit(z);
    }
    let least = times(&x).iter().map(|v| v * v).sum::<f64>().sqrt();
    let mut big = unit(vec![1.0; n]);
    let mut most = 0.0;
    for _ in 0..30 {
        let y = times(&big);
        let back: Vec<f64> = (0..n)
            .map(|j| (0..=j).map(|i| r(i, j) * y[i]).sum())
            .collect();
        most = back.iter().map(|v| v * v).sum::<f64>().sqrt().sqrt();
        big = unit(back);
    }
    (most > 0.0).then(|| (x, least / most))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_null_vector_of_points_on_a_sphere_is_the_sphere() {
        let ms = monomials(2);
        let points: Vec<[f64; 3]> = (0..40)
            .map(|k| {
                let (a, b) = (k as f64 * 0.37, k as f64 * 0.91);
                [
                    0.8 * a.cos() * b.sin(),
                    0.8 * a.sin() * b.sin(),
                    0.8 * b.cos(),
                ]
            })
            .collect();
        let columns: Vec<Vec<f64>> = ms
            .iter()
            .map(|m| {
                points
                    .iter()
                    .map(|x| {
                        x[0].powi(m[0] as i32) * x[1].powi(m[1] as i32) * x[2].powi(m[2] as i32)
                    })
                    .collect()
            })
            .collect();
        let (v, ratio) = null_vector(&columns).expect("a sphere");
        assert!(ratio < 1e-14, "{ratio}");
        // x^2 + y^2 + z^2 - 0.64, up to scale.
        let at = |e: [usize; 3]| v[ms.iter().position(|m| *m == e).unwrap()];
        let unit = at([0, 0, 0]) / -0.64;
        for e in [[2, 0, 0], [0, 2, 0], [0, 0, 2]] {
            assert!((at(e) - unit).abs() < 1e-12, "{e:?}");
        }
        for e in [
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 1, 0],
            [1, 0, 1],
            [0, 1, 1],
        ] {
            assert!(at(e).abs() < 1e-12, "{e:?}");
        }
    }
}
