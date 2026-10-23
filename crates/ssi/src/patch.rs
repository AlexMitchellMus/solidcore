//! Rational B-spline patches, in the layout Form Studio's kernel emits, and
//! their rational Bézier pieces.
//!
//! A patch has `u_count × v_count` Cartesian poles, pole (i, j) at
//! `poles[i * v_count + j]` with i along u, optional weights (all 1 when
//! absent) and full knot vectors. Its domain is `[u_knots[p], u_knots[n + 1]]`
//! by the same in v, clamped or not.

use crate::surface3::{NurbsSurf, Surf};

use crate::piece::{Piece, H, MAX_DEGREE};
use crate::Error;

/// A tensor-product rational B-spline patch.
#[derive(Clone, Debug)]
pub struct Patch {
    pub u_degree: usize,
    pub v_degree: usize,
    pub u_count: usize,
    pub v_count: usize,
    pub poles: Vec<[f64; 3]>,
    pub weights: Option<Vec<f64>>,
    pub u_knots: Vec<f64>,
    pub v_knots: Vec<f64>,
    /// The analytic surface the patch lies on, when known (a plane, cylinder,
    /// cone, sphere or torus): its implicit equation then bounds, and proves,
    /// how other patches meet it (see `field`).
    pub analytic: Option<Surf>,
}

/// The span index of `t`: knots[k] <= t < knots[k + 1], within [p, n]; the domain's end goes to its last span.
fn span_of(knots: &[f64], p: usize, n: usize, t: f64) -> usize {
    if t >= knots[n + 1] {
        let mut k = n;
        while k > p && !(knots[k] < knots[k + 1]) {
            k -= 1;
        }
        return k;
    }
    let mut k = p;
    while k < n && knots[k + 1] <= t {
        k += 1;
    }
    k
}

/// The p + 1 nonzero basis functions at `t` in span `k` (The NURBS Book, A2.2).
fn basis_functions(knots: &[f64], p: usize, k: usize, t: f64) -> Vec<f64> {
    let mut n = vec![0.0; p + 1];
    let mut left = vec![0.0; p + 1];
    let mut right = vec![0.0; p + 1];
    n[0] = 1.0;
    for j in 1..=p {
        left[j] = t - knots[k + 1 - j];
        right[j] = knots[k + j] - t;
        let mut saved = 0.0;
        for r in 0..j {
            let temp = n[r] / (right[r + 1] + left[j - r]);
            n[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        n[j] = saved;
    }
    n
}

impl Patch {
    fn weight(&self, index: usize) -> f64 {
        self.weights.as_ref().map_or(1.0, |w| w[index])
    }

    /// The patch at (u, v) in its own knot domain.
    pub fn evaluate(&self, u: f64, v: f64) -> [f64; 3] {
        let (p, q) = (self.u_degree, self.v_degree);
        let ku = span_of(&self.u_knots, p, self.u_count - 1, u);
        let kv = span_of(&self.v_knots, q, self.v_count - 1, v);
        let nu = basis_functions(&self.u_knots, p, ku, u);
        let nv = basis_functions(&self.v_knots, q, kv, v);
        let mut h = [0.0; 4];
        for a in 0..=p {
            let i = ku - p + a;
            for b in 0..=q {
                let index = i * self.v_count + kv - q + b;
                let weight = nu[a] * nv[b] * self.weight(index);
                let pole = self.poles[index];
                h[0] += weight * pole[0];
                h[1] += weight * pole[1];
                h[2] += weight * pole[2];
                h[3] += weight;
            }
        }
        [h[0] / h[3], h[1] / h[3], h[2] / h[3]]
    }

    /// The patch's rational Bézier pieces, each tagged with `index` and its knot rectangle.
    pub fn bezier_pieces(&self, index: usize) -> Result<Vec<Piece>, Error> {
        let (p, q) = (self.u_degree, self.v_degree);
        if p == 0 || q == 0 || p > MAX_DEGREE || q > MAX_DEGREE {
            return Err(Error::BadPatch(format!(
                "patch {index} has degrees ({p}, {q}); 1 to {MAX_DEGREE} are taken"
            )));
        }
        if self.poles.len() != self.u_count * self.v_count
            || self.u_knots.len() != self.u_count + p + 1
            || self.v_knots.len() != self.v_count + q + 1
            || self
                .weights
                .as_ref()
                .is_some_and(|w| w.len() != self.poles.len())
        {
            return Err(Error::BadPatch(format!(
                "patch {index}: pole, weight and knot counts do not agree"
            )));
        }
        let mut homogeneous = Vec::with_capacity(self.poles.len());
        for (k, pole) in self.poles.iter().enumerate() {
            let w = self.weight(k);
            if !(w > 0.0) {
                return Err(Error::BadPatch(format!(
                    "patch {index} has a weight that is not positive: {w}"
                )));
            }
            homogeneous.push([w * pole[0], w * pole[1], w * pole[2], w]);
        }
        // Along v for every u row, then along u for every column of each v span.
        let mut rows = Vec::with_capacity(self.u_count);
        for i in 0..self.u_count {
            let row = homogeneous[i * self.v_count..(i + 1) * self.v_count].to_vec();
            rows.push(decompose_curve(&self.v_knots, q, row)?);
        }
        let closed = self.closed();
        let mut pieces = Vec::new();
        let v_spans = rows[0].len();
        for s in 0..v_spans {
            let mut columns = Vec::with_capacity(q + 1);
            for j in 0..=q {
                let column: Vec<H> = rows.iter().map(|row| row[s].points[j]).collect();
                columns.push(decompose_curve(&self.u_knots, p, column)?);
            }
            let u_spans = columns[0].len();
            for r in 0..u_spans {
                let mut c = Vec::with_capacity((p + 1) * (q + 1));
                for i in 0..=p {
                    for column in &columns {
                        c.push(column[r].points[i]);
                    }
                }
                let above = [r + 1 < u_spans || closed[0], s + 1 < v_spans || closed[1]];
                pieces.push(Piece::new(
                    index,
                    columns[0][r].t0,
                    columns[0][r].t1,
                    rows[0][s].t0,
                    rows[0][s].t1,
                    p,
                    q,
                    c,
                    above,
                ));
            }
        }
        Ok(pieces)
    }

    /// Whether one of the patch's four boundary edges is collapsed to a point
    /// (a pole, a cone's apex): all its poles one point.
    pub fn has_collapsed_edge(&self) -> bool {
        let (m, n) = (self.u_count, self.v_count);
        let same = |ids: &mut dyn Iterator<Item = usize>| {
            let pts: Vec<[f64; 3]> = ids.map(|k| self.poles[k]).collect();
            pts.windows(2)
                .all(|w| (0..3).all(|c| (w[0][c] - w[1][c]).abs() <= 1e-12 * (1.0 + w[0][c].abs())))
        };
        same(&mut (0..n))
            || same(&mut ((m - 1) * n..m * n))
            || same(&mut (0..m).map(|i| i * n))
            || same(&mut (0..m).map(|i| i * n + n - 1))
    }

    /// The patch's knot domain: [u_start, u_end, v_start, v_end].
    pub fn domain(&self) -> [f64; 4] {
        let (p, q) = (self.u_degree, self.v_degree);
        [
            self.u_knots[p],
            self.u_knots[self.u_count],
            self.v_knots[q],
            self.v_knots[self.v_count],
        ]
    }

    /// The patch as a `surface3` B-spline surface, for its certified closest point.
    pub fn nurbs(&self) -> Result<NurbsSurf, Error> {
        let n = self.v_count;
        let at = |i: usize, j: usize| {
            let p = self.poles[i * n + j];
            crate::geom3::V3::new(p[0], p[1], p[2])
        };
        let ctrl: Vec<Vec<_>> = (0..self.u_count)
            .map(|i| (0..n).map(|j| at(i, j)).collect())
            .collect();
        let weights: Option<Vec<Vec<f64>>> = self.weights.as_ref().map(|w| {
            (0..self.u_count)
                .map(|i| w[i * n..(i + 1) * n].to_vec())
                .collect()
        });
        NurbsSurf::new(
            self.u_degree,
            self.v_degree,
            &ctrl,
            weights.as_ref(),
            &self.u_knots,
            &self.v_knots,
        )
        .map_err(Error::BadPatch)
    }

    /// Whether the patch closes on itself in u and in v: its first and last
    /// rows of poles (and weights) coincide, so its two edges there are one seam.
    pub fn closed(&self) -> [bool; 2] {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for pole in &self.poles {
            for k in 0..3 {
                lo[k] = lo[k].min(pole[k]);
                hi[k] = hi[k].max(pole[k]);
            }
        }
        let tol = 1e-12 * (0..3).map(|k| hi[k] - lo[k]).fold(0.0, f64::max);
        let same = |i: usize, j: usize| {
            let (a, b) = (self.poles[i], self.poles[j]);
            (0..3).all(|k| (a[k] - b[k]).abs() <= tol)
                && (self.weight(i) - self.weight(j)).abs() <= 1e-12 * self.weight(i)
        };
        let (m, n) = (self.u_count, self.v_count);
        [
            (0..n).all(|j| same(j, (m - 1) * n + j)),
            (0..m).all(|i| same(i * n, i * n + n - 1)),
        ]
    }
}

/// One Boehm insertion of `t` into a curve of homogeneous points (The NURBS Book, A5.1).
fn insert_knot(knots: &[f64], p: usize, points: &[H], t: f64) -> (Vec<f64>, Vec<H>) {
    let n = points.len() - 1;
    let k = span_of(knots, p, n, t);
    let mut out = Vec::with_capacity(n + 2);
    out.extend_from_slice(&points[..=k - p]);
    for i in k - p + 1..=k {
        let alpha = (t - knots[i]) / (knots[i + p] - knots[i]);
        let (a, b) = (points[i - 1], points[i]);
        out.push([0, 1, 2, 3].map(|c| (1.0 - alpha) * a[c] + alpha * b[c]));
    }
    out.extend_from_slice(&points[k..]);
    let mut grown = knots[..=k].to_vec();
    grown.push(t);
    grown.extend_from_slice(&knots[k + 1..]);
    (grown, out)
}

struct Span {
    t0: f64,
    t1: f64,
    points: Vec<H>,
}

/// A curve's Bézier spans over its domain [knots[p], knots[n + 1]]: every knot
/// value there raised to multiplicity p, clamped or not, so span k's p + 1
/// control points are points[k - p ..= k].
fn decompose_curve(knots_in: &[f64], p: usize, points_in: Vec<H>) -> Result<Vec<Span>, Error> {
    let mut knots = knots_in.to_vec();
    let mut points = points_in;
    let n0 = points.len() - 1;
    let (start, end) = (knots[p], knots[n0 + 1]);
    if !(start < end) {
        return Err(Error::BadPatch(format!(
            "an empty knot domain [{start}, {end}]"
        )));
    }
    let mut values: Vec<f64> = knots
        .iter()
        .copied()
        .filter(|&t| t >= start && t <= end)
        .collect();
    values.sort_by(f64::total_cmp);
    values.dedup();
    for t in values {
        let mut multiplicity = knots.iter().filter(|&&knot| knot == t).count();
        while multiplicity < p {
            (knots, points) = insert_knot(&knots, p, &points, t);
            multiplicity += 1;
        }
    }
    let mut spans = Vec::new();
    for k in p..points.len() {
        if knots[k] < knots[k + 1] && knots[k] >= start && knots[k + 1] <= end {
            spans.push(Span {
                t0: knots[k],
                t1: knots[k + 1],
                points: points[k - p..=k].to_vec(),
            });
        }
    }
    Ok(spans)
}
