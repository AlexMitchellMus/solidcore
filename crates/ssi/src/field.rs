//! The field on a Bézier piece's parameter rectangle: an analytic surface's
//! implicit equation composed with the piece, as a tensor Bernstein polynomial.
//!
//! On a rational piece x(s, t) = X(s, t) / W(s, t), an implicit equation f of
//! degree k gives g = W^k f(X / W), a polynomial of degree (k du, k dv). W is
//! positive, so g has the sign of f at every point of the piece, and g's
//! Bernstein coefficients bound it there (Farouki 1986; Sederberg and Nishita
//! 1990). Every coefficient zero: the piece lies on the surface. Every
//! coefficient of one sign: the piece is apart from it. Both verdicts are read
//! off the coefficients, a proof up to rounding, and neither looks at a sampled
//! point's distance.

use crate::geom3::{Frame, V3};
use crate::surface3::Surf;

use crate::algebraic::Algebraic;
use crate::piece::Piece;

/// A tensor Bernstein polynomial of degree (du, dv), coefficient (i, j) at i * (dv + 1) + j.
#[derive(Clone, Debug)]
pub struct Bern {
    pub du: usize,
    pub dv: usize,
    pub c: Vec<f64>,
}

fn binomial(n: usize, k: usize) -> f64 {
    with_row(n, |row| row[k])
}

thread_local! {
    /// Pascal's triangle, grown a row at a time as far as asked: rebuilt for
    /// every product the binomial rows cost as much as the product itself, and
    /// built to row 255 up front they were most of a small case's time.
    static PASCAL: std::cell::RefCell<Vec<Vec<f64>>> = std::cell::RefCell::new(vec![vec![1.0]]);
}

/// `f` of row n of Pascal's triangle.
fn with_row<R>(n: usize, f: impl FnOnce(&[f64]) -> R) -> R {
    with_rows(n, |rows| f(&rows[n]))
}

/// `f` of Pascal's triangle, to row n at least.
fn with_rows<R>(n: usize, f: impl FnOnce(&[Vec<f64>]) -> R) -> R {
    PASCAL.with(|rows| {
        {
            let mut rows = rows.borrow_mut();
            while rows.len() <= n {
                let last = &rows[rows.len() - 1];
                let m = last.len();
                let row: Vec<f64> = (0..=m)
                    .map(|k| {
                        if k == 0 || k == m {
                            1.0
                        } else {
                            last[k - 1] + last[k]
                        }
                    })
                    .collect();
                rows.push(row);
            }
        }
        f(&rows.borrow())
    })
}

impl Bern {
    /// The product, by the Bernstein product rule: C(m, i) C(n, k) / C(m + n, i + k)
    /// in each direction, scaled in and out once instead of per term (binomial
    /// rows computed once per product; per term they were most of its cost).
    pub fn mul(&self, o: &Bern) -> Bern {
        let (du, dv) = (self.du + o.du, self.dv + o.dv);
        with_rows(du.max(dv), |rows| self.mul_with(o, rows))
    }

    /// `mul` with Pascal's triangle at hand (borrowed, not copied per product).
    fn mul_with(&self, o: &Bern, rows: &[Vec<f64>]) -> Bern {
        let (du, dv) = (self.du + o.du, self.dv + o.dv);
        let (bu_a, bv_a, bu_b, bv_b, bu, bv) = (
            &rows[self.du],
            &rows[self.dv],
            &rows[o.du],
            &rows[o.dv],
            &rows[du],
            &rows[dv],
        );
        // Scaled coefficients C(m, i) C(p, j) a_ij, multiplied as power-like terms, then divided back.
        let a: Vec<f64> = (0..=self.du)
            .flat_map(|i| (0..=self.dv).map(move |j| (i, j)))
            .map(|(i, j)| bu_a[i] * bv_a[j] * self.c[i * (self.dv + 1) + j])
            .collect();
        let b: Vec<f64> = (0..=o.du)
            .flat_map(|k| (0..=o.dv).map(move |l| (k, l)))
            .map(|(k, l)| bu_b[k] * bv_b[l] * o.c[k * (o.dv + 1) + l])
            .collect();
        let mut c = vec![0.0; (du + 1) * (dv + 1)];
        for i in 0..=self.du {
            for j in 0..=self.dv {
                let x = a[i * (self.dv + 1) + j];
                if x == 0.0 {
                    continue;
                }
                for k in 0..=o.du {
                    let base = (i + k) * (dv + 1) + j;
                    let row_b = &b[k * (o.dv + 1)..(k + 1) * (o.dv + 1)];
                    for (l, y) in row_b.iter().enumerate() {
                        c[base + l] += x * y;
                    }
                }
            }
        }
        for i in 0..=du {
            for j in 0..=dv {
                c[i * (dv + 1) + j] /= bu[i] * bv[j];
            }
        }
        Bern { du, dv, c }
    }

    /// Σ weight × term, over terms of one degree.
    pub(crate) fn sum(terms: &[(f64, &Bern)]) -> Bern {
        let first = terms[0].1;
        let mut c = vec![0.0; first.c.len()];
        for (w, t) in terms {
            for (k, x) in t.c.iter().enumerate() {
                c[k] += w * x;
            }
        }
        Bern {
            du: first.du,
            dv: first.dv,
            c,
        }
    }

    pub fn max_abs(&self) -> f64 {
        self.c.iter().fold(0.0, |m, x| m.max(x.abs()))
    }

    fn get(&self, i: usize, j: usize) -> f64 {
        self.c[i * (self.dv + 1) + j]
    }

    /// The value at (s, t) in [0, 1]^2, by de Casteljau along t then along s,
    /// on the stack up to degree 63 (a field point's root takes a dozen of
    /// these, and heap vectors for them were most of an emitted point's cost).
    pub fn at(&self, s: f64, t: f64) -> f64 {
        if self.du >= STACK || self.dv >= STACK {
            let column: Vec<f64> = (0..=self.du)
                .map(|i| casteljau(&self.c[i * (self.dv + 1)..(i + 1) * (self.dv + 1)], t))
                .collect();
            return casteljau(&column, s);
        }
        let mut column = [0.0; STACK];
        for (i, value) in column.iter_mut().enumerate().take(self.du + 1) {
            *value = casteljau(&self.c[i * (self.dv + 1)..(i + 1) * (self.dv + 1)], t);
        }
        casteljau(&column[..=self.du], s)
    }

    /// The two halves at s = x (dir 0) or t = x (dir 1), each over its own [0, 1].
    pub fn split(&self, dir: usize, x: f64) -> (Bern, Bern) {
        let (du, dv) = (self.du, self.dv);
        let (n, m) = if dir == 0 { (du, dv) } else { (dv, du) };
        let index = |i: usize, j: usize| {
            if dir == 0 {
                i * (dv + 1) + j
            } else {
                j * (dv + 1) + i
            }
        };
        let mut left = vec![0.0; self.c.len()];
        let mut right = vec![0.0; self.c.len()];
        let mut work = vec![0.0; n + 1];
        for j in 0..=m {
            for i in 0..=n {
                work[i] = self.c[index(i, j)];
            }
            for r in 0..=n {
                left[index(r, j)] = work[0];
                right[index(n - r, j)] = work[n - r];
                for i in 0..n - r {
                    work[i] = (1.0 - x) * work[i] + x * work[i + 1];
                }
            }
        }
        (Bern { du, dv, c: left }, Bern { du, dv, c: right })
    }

    /// The polynomial over [s0, s1] x [t0, t1], on its own [0, 1]^2.
    pub fn restrict(&self, s0: f64, s1: f64, t0: f64, t1: f64) -> Bern {
        let mut b = self.clone();
        for (dir, lo, hi) in [(0, s0, s1), (1, t0, t1)] {
            // Bounds past [0, 1] extrapolate: g's polynomial runs on past its square.
            if hi != 1.0 {
                b = b.split(dir, hi).0;
            }
            if lo != 0.0 {
                b = b.split(dir, lo / hi).1;
            }
        }
        b
    }

    /// The polynomial along the segment from (s0, t0) to (s1, t1) of its square,
    /// as a one-variable Bernstein polynomial of degree du + dv over the segment:
    /// the tensor's coefficients over the box with those corners (rows restricted
    /// along t, then columns along s, either way round), summed along the box's
    /// diagonal with weights C(du, i) C(dv, j) / C(du + dv, i + j).
    pub fn along(&self, (s0, t0): (f64, f64), (s1, t1): (f64, f64)) -> Vec<f64> {
        let (du, dv) = (self.du, self.dv);
        let mut q = self.c.clone();
        let mut work = Vec::with_capacity(du.max(dv) + 1);
        for row in q.chunks_mut(dv + 1) {
            restrict1_into(row, t0, t1, &mut work);
        }
        let mut column = vec![0.0; du + 1];
        for j in 0..=dv {
            for i in 0..=du {
                column[i] = q[i * (dv + 1) + j];
            }
            restrict1_into(&mut column, s0, s1, &mut work);
            for i in 0..=du {
                q[i * (dv + 1) + j] = column[i];
            }
        }
        let n = du + dv;
        let mut out = vec![0.0; n + 1];
        for i in 0..=du {
            for j in 0..=dv {
                out[i + j] += binomial(du, i) * binomial(dv, j) * q[i * (dv + 1) + j];
            }
        }
        for (k, x) in out.iter_mut().enumerate() {
            *x /= binomial(n, k);
        }
        out
    }

    /// h with g = x h, x the factor that vanishes on edge e (s, 1 - s, t or
    /// 1 - t for edges 0 to 3), for g zero on that edge: one degree lower in
    /// that direction, from x B_k^(n-1) = (k + 1) / n B_(k+1)^n and
    /// (1 - x) B_k^(n-1) = (n - k) / n B_k^n.
    pub fn divide_edge(&self, e: usize) -> Bern {
        let (du, dv) = (self.du, self.dv);
        let mut c = Vec::with_capacity(self.c.len());
        match e {
            0 | 1 => {
                for k in 0..du {
                    for j in 0..=dv {
                        c.push(if e == 0 {
                            self.get(k + 1, j) * du as f64 / (k + 1) as f64
                        } else {
                            self.get(k, j) * du as f64 / (du - k) as f64
                        });
                    }
                }
                Bern { du: du - 1, dv, c }
            }
            _ => {
                for i in 0..=du {
                    for k in 0..dv {
                        c.push(if e == 2 {
                            self.get(i, k + 1) * dv as f64 / (k + 1) as f64
                        } else {
                            self.get(i, k) * dv as f64 / (dv - k) as f64
                        });
                    }
                }
                Bern { du, dv: dv - 1, c }
            }
        }
    }

    /// The polynomial on the line s = x (dir 0, in t) or t = x (dir 1, in s).
    pub fn line(&self, dir: usize, x: f64) -> Vec<f64> {
        let (du, dv) = (self.du, self.dv);
        if dir == 0 {
            let mut column = vec![0.0; du + 1];
            (0..=dv)
                .map(|j| {
                    for (i, c) in column.iter_mut().enumerate() {
                        *c = self.get(i, j);
                    }
                    casteljau(&column, x)
                })
                .collect()
        } else {
            (0..=du)
                .map(|i| casteljau(&self.c[i * (dv + 1)..(i + 1) * (dv + 1)], x))
                .collect()
        }
    }

    /// The derivative along s (dir 0) or t (dir 1), a polynomial one degree lower there.
    pub fn derivative(&self, dir: usize) -> Bern {
        let (du, dv) = (self.du, self.dv);
        if (dir == 0 && du == 0) || (dir == 1 && dv == 0) {
            return Bern {
                du: 0,
                dv: 0,
                c: vec![0.0],
            };
        }
        let (nu, nv) = if dir == 0 { (du - 1, dv) } else { (du, dv - 1) };
        let mut c = Vec::with_capacity((nu + 1) * (nv + 1));
        for i in 0..=nu {
            for j in 0..=nv {
                c.push(if dir == 0 {
                    du as f64 * (self.get(i + 1, j) - self.get(i, j))
                } else {
                    dv as f64 * (self.get(i, j + 1) - self.get(i, j))
                });
            }
        }
        Bern { du: nu, dv: nv, c }
    }

    /// The edge s = 0, s = 1, t = 0 or t = 1 (edges 0 to 3) as a polynomial in the other parameter.
    pub fn edge(&self, e: usize) -> Vec<f64> {
        let (du, dv) = (self.du, self.dv);
        match e {
            0 => (0..=dv).map(|j| self.get(0, j)).collect(),
            1 => (0..=dv).map(|j| self.get(du, j)).collect(),
            2 => (0..=du).map(|i| self.get(i, 0)).collect(),
            _ => (0..=du).map(|i| self.get(i, dv)).collect(),
        }
    }

    /// Whether every coefficient is above `clear`, or every one below -`clear`.
    pub fn one_sign(&self, clear: f64) -> bool {
        self.c.iter().all(|&x| x > clear) || self.c.iter().all(|&x| x < -clear)
    }
}

/// Degrees below this are evaluated on the stack.
const STACK: usize = 64;

/// A one-variable Bernstein polynomial at t.
pub fn casteljau(c: &[f64], t: f64) -> f64 {
    let n = c.len();
    if n > STACK {
        let mut work = c.to_vec();
        for r in 1..n {
            for i in 0..n - r {
                work[i] = (1.0 - t) * work[i] + t * work[i + 1];
            }
        }
        return work[0];
    }
    let mut work = [0.0; STACK];
    work[..n].copy_from_slice(c);
    for r in 1..n {
        for i in 0..n - r {
            work[i] = (1.0 - t) * work[i] + t * work[i + 1];
        }
    }
    work[0]
}

/// The roots in [0, 1] of a one-variable Bernstein polynomial, each isolated by
/// subdivision until its coefficients change sign once (Descartes' rule for
/// Bernstein coefficients), then found by bisection held to that interval. A
/// polynomial that is zero to `clear` everywhere has no isolated roots: None.
pub fn roots(c: &[f64], clear: f64) -> Option<Vec<f64>> {
    if c.iter().all(|x| x.abs() <= clear) {
        return None;
    }
    let mut out = Vec::new();
    // Intervals on which the polynomial is zero to `clear` throughout: a double
    // root's rounding band, each run of them one root at its middle. Halving on
    // inside such a band keeps both halves of it alive at every level (2^52 of
    // them at the depth limit, a hang); it is stopped where the coefficients say so.
    let mut flat: Vec<(f64, f64)> = Vec::new();
    let mut stack = vec![(c.to_vec(), 0.0, 1.0, 0)];
    let mut budget = 20_000;
    while let Some((b, lo, hi, depth)) = stack.pop() {
        budget -= 1;
        if b.iter().all(|&x| x > clear) || b.iter().all(|&x| x < -clear) {
            continue;
        }
        if b.iter().all(|x| x.abs() <= clear) || depth > 52 || budget <= 0 {
            flat.push((lo, hi));
            continue;
        }
        let changes = b
            .windows(2)
            .filter(|w| (w[0] > 0.0) != (w[1] > 0.0))
            .count();
        let ends_differ = (b[0] > 0.0) != (b[b.len() - 1] > 0.0);
        if changes == 1 && ends_differ {
            // One simple root: bisection on the sign, to the last bit.
            let (mut a, mut z) = (0.0_f64, 1.0_f64);
            let fa = casteljau(&b, 0.0);
            for _ in 0..80 {
                let mid = 0.5 * (a + z);
                if (casteljau(&b, mid) > 0.0) == (fa > 0.0) {
                    a = mid;
                } else {
                    z = mid;
                }
            }
            out.push(lo + 0.5 * (a + z) * (hi - lo));
            continue;
        }
        let mid = 0.5 * (lo + hi);
        let (l, r) = split1(&b);
        stack.push((r, mid, hi, depth + 1));
        stack.push((l, lo, mid, depth + 1));
    }
    flat.sort_by(|x, y| x.0.total_cmp(&y.0));
    let mut run: Option<(f64, f64)> = None;
    for (lo, hi) in flat {
        run = match run {
            Some((a, z)) if lo <= z + 1e-12 => Some((a, z.max(hi))),
            Some((a, z)) => {
                out.push(0.5 * (a + z));
                Some((lo, hi))
            }
            None => Some((lo, hi)),
        };
    }
    if let Some((a, z)) = run {
        out.push(0.5 * (a + z));
    }
    out.sort_by(f64::total_cmp);
    out.dedup_by(|a, b| (*a - *b).abs() <= 1e-14);
    Some(out)
}

/// A one-variable Bernstein polynomial's coefficients over [a, b] of its [0, 1]
/// (b below a runs it backward): coefficient k is its blossom at a (n - k
/// times) and b (k times). In place, `work` a scratch buffer.
pub(crate) fn restrict1_into(c: &mut [f64], a: f64, b: f64, work: &mut Vec<f64>) {
    let n = c.len() - 1;
    let source = c.to_vec();
    for (k, out) in c.iter_mut().enumerate() {
        work.clear();
        work.extend_from_slice(&source);
        for r in 0..n {
            let x = if r < k { b } else { a };
            for i in 0..n - r {
                work[i] = (1.0 - x) * work[i] + x * work[i + 1];
            }
        }
        *out = work[0];
    }
}

/// A one-variable Bernstein polynomial raised to degree `to` (at least its own).
pub(crate) fn elevate1(c: &[f64], to: usize) -> Vec<f64> {
    let mut c = c.to_vec();
    while c.len() - 1 < to {
        let n = c.len() - 1;
        let mut up = vec![0.0; n + 2];
        up[0] = c[0];
        up[n + 1] = c[n];
        for i in 1..=n {
            let a = i as f64 / (n + 1) as f64;
            up[i] = a * c[i - 1] + (1.0 - a) * c[i];
        }
        c = up;
    }
    c
}

/// A one-variable Bernstein polynomial's halves at 1/2.
pub(crate) fn split1(c: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let n = c.len() - 1;
    let mut work = c.to_vec();
    let mut left = vec![0.0; n + 1];
    let mut right = vec![0.0; n + 1];
    for r in 0..=n {
        left[r] = work[0];
        right[n - r] = work[n - r];
        for i in 0..n - r {
            work[i] = 0.5 * (work[i] + work[i + 1]);
        }
    }
    (left, right)
}

/// How a piece stands to an analytic surface.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Relation {
    /// Every coefficient of g is zero (to 1e-9 of the terms' scale): the piece lies on the surface.
    On,
    /// Every coefficient has one sign, clear of rounding: the piece does not meet the surface.
    Apart,
    /// Neither: the piece may cross or touch the surface.
    Unknown,
}

/// The piece's homogeneous points in the frame: local x, y, z (each times W) and W, of the piece's degree.
fn local(piece: &Piece, f: &Frame) -> [Bern; 4] {
    let (y, z) = (f.y(), f.z);
    let axes = [f.x, y, z];
    let mut out = [(); 4].map(|_| Bern {
        du: piece.du,
        dv: piece.dv,
        c: Vec::with_capacity(piece.c.len()),
    });
    for h in &piece.c {
        let d = V3::new(
            h[0] - f.o.x * h[3],
            h[1] - f.o.y * h[3],
            h[2] - f.o.z * h[3],
        );
        for (k, axis) in axes.iter().enumerate() {
            out[k].c.push(d.dot(*axis));
        }
        out[3].c.push(h[3]);
    }
    out
}

/// The composed field g = W^k f(X / W) and the scale of its terms, or None for a B-spline surface.
pub fn compose(piece: &Piece, surf: &Surf) -> Option<(Bern, f64)> {
    let f = surf.frame()?;
    let [lx, ly, lz, w] = local(piece, f);
    // m bounds every term's size: the piece's coordinates and the surface's own lengths times W.
    let length = match surf {
        Surf::Plane { .. } | Surf::Nurbs(_) => 0.0,
        Surf::Cylinder { r, .. } | Surf::Sphere { r, .. } => *r,
        Surf::Cone { r, .. } => r.abs(),
        Surf::Torus { major, minor, .. } => major + minor,
    };
    let m = lx
        .max_abs()
        .max(ly.max_abs())
        .max(lz.max_abs())
        .max(length * w.max_abs());
    let square = |b: &Bern| b.mul(b);
    let radial = || Bern::sum(&[(1.0, &square(&lx)), (1.0, &square(&ly))]);
    let g = match surf {
        Surf::Plane { .. } => return Some((lz, m)),
        Surf::Sphere { r, .. } => {
            Bern::sum(&[(1.0, &radial()), (1.0, &square(&lz)), (-r * r, &square(&w))])
        }
        Surf::Cylinder { r, .. } => Bern::sum(&[(1.0, &radial()), (-r * r, &square(&w))]),
        Surf::Cone { r, a, .. } => {
            // STEP's cone: radius r + v tan a at height v.
            let rim = Bern::sum(&[(*r, &w), (a.tan(), &lz)]);
            Bern::sum(&[(1.0, &radial()), (-1.0, &square(&rim))])
        }
        Surf::Torus { major, minor, .. } => {
            let w2 = square(&w);
            let s = Bern::sum(&[
                (1.0, &radial()),
                (1.0, &square(&lz)),
                (major * major - minor * minor, &w2),
            ]);
            Bern::sum(&[
                (1.0, &square(&s)),
                (-4.0 * major * major, &w2.mul(&radial())),
            ])
        }
        Surf::Nurbs(_) => return None,
    };
    let k = if matches!(surf, Surf::Torus { .. }) {
        4
    } else {
        2
    };
    Some((g, m.powi(k)))
}

/// A surface with an implicit equation: analytic, or a B-spline patch's proved
/// algebraic form (`algebraic`).
#[derive(Clone, Copy, Debug)]
pub enum Field<'f> {
    Analytic(&'f Surf),
    Algebraic(&'f Algebraic),
}

impl Field<'_> {
    /// g = W^k f(X / W) on the piece and the scale of its terms, as `compose`.
    pub fn compose(&self, piece: &Piece) -> Option<(Bern, f64)> {
        match self {
            Field::Analytic(surf) => compose(piece, surf),
            Field::Algebraic(form) => Some(form.compose(piece)),
        }
    }

    /// The foot of p on the surface (to first order for an algebraic form).
    pub fn foot(&self, p: [f64; 3]) -> [f64; 3] {
        match self {
            Field::Analytic(surf) => {
                let f = surf.foot(V3::new(p[0], p[1], p[2]));
                [f.x, f.y, f.z]
            }
            Field::Algebraic(form) => form.foot(p),
        }
    }

    /// The implicit equation's degree: plane 1, quadrics 2, torus 4.
    pub fn degree(&self) -> usize {
        match self {
            Field::Analytic(Surf::Plane { .. }) => 1,
            Field::Analytic(Surf::Torus { .. }) => 4,
            Field::Analytic(_) => 2,
            Field::Algebraic(form) => form.degree,
        }
    }
}

/// How a piece stands to an analytic surface: on it, apart from it, or unknown.
pub fn relate(piece: &Piece, surf: &Surf) -> Relation {
    relate_field(piece, &Field::Analytic(surf))
}

/// `relate` for any field.
pub fn relate_field(piece: &Piece, field: &Field) -> Relation {
    let Some((g, scale)) = field.compose(piece) else {
        return Relation::Unknown;
    };
    if !(scale > 0.0) {
        return Relation::Unknown;
    }
    if g.max_abs() <= 1e-9 * scale {
        return Relation::On;
    }
    let clear = 1e-12 * scale;
    if g.c.iter().all(|&x| x > clear) || g.c.iter().all(|&x| x < -clear) {
        return Relation::Apart;
    }
    Relation::Unknown
}
