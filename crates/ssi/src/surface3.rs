//! Exact surfaces of B-rep faces: planes, circular cylinders and cones, spheres,
//! tori and rational B-spline surfaces, with the surface's own normal and the
//! closest point of an untrimmed surface.
//!
//! Analytic surfaces use the STEP (ISO 10303-42) parameterizations:
//!   cylinder  o + r e_r(u) + v z
//!   cone      o + (r + v tan a) e_r(u) + v z       (a = semi-angle)
//!   sphere    o + r cos v e_r(u) + r sin v z
//!   torus     o + (R + r cos v) e_r(u) + r sin v z
//! with e_r(u) = cos u x + sin u y in the surface frame. The surface normal is
//! the STEP one (outward from the axis, center or tube center; Su x Sv for
//! B-spline surfaces).
//!
//! Each analytic surface has exactly one local minimum of the distance to a
//! point (the foot), so the closest point of a trimmed face is its foot when
//! the foot lies in the face and on the face's boundary otherwise. B-spline
//! surfaces can have several; their closest point is a certified branch and
//! bound over Bezier patches.

use std::sync::Arc;

use crate::geom3::{clamp, insert_knot, Frame, V3};

#[derive(Clone, Debug)]
pub enum Surf {
    Plane { f: Frame },
    Cylinder { f: Frame, r: f64 },
    Cone { f: Frame, r: f64, a: f64 },
    Sphere { f: Frame, r: f64 },
    Torus { f: Frame, major: f64, minor: f64 },
    Nurbs(Arc<NurbsSurf>),
}

/// Axial coordinates of p in a frame: angle about z, distance from the axis,
/// height along z, and the radial and tangential unit vectors at that angle.
#[derive(Clone, Copy, Debug)]
pub struct Axial {
    pub theta: f64,
    pub rho: f64,
    pub h: f64,
    pub er: V3,
    pub et: V3,
}

pub fn axial(f: &Frame, p: V3) -> Axial {
    let l = f.local(p);
    let rho = (l.x * l.x + l.y * l.y).sqrt();
    let theta = if rho > 0.0 { l.y.atan2(l.x) } else { 0.0 };
    let (s, c) = theta.sin_cos();
    Axial {
        theta,
        rho,
        h: l.z,
        er: f.dir(V3::new(c, s, 0.0)),
        et: f.dir(V3::new(-s, c, 0.0)),
    }
}

impl Surf {
    pub fn name(&self) -> &'static str {
        match self {
            Surf::Plane { .. } => "plane",
            Surf::Cylinder { .. } => "cylinder",
            Surf::Cone { .. } => "cone",
            Surf::Sphere { .. } => "sphere",
            Surf::Torus { .. } => "torus",
            Surf::Nurbs(_) => "B-spline surface",
        }
    }

    pub fn frame(&self) -> Option<&Frame> {
        match self {
            Surf::Plane { f }
            | Surf::Cylinder { f, .. }
            | Surf::Cone { f, .. }
            | Surf::Sphere { f, .. }
            | Surf::Torus { f, .. } => Some(f),
            Surf::Nurbs(_) => None,
        }
    }

    /// Apex of a cone (height -r / tan a on its axis).
    pub fn apex(&self) -> Option<V3> {
        match self {
            Surf::Cone { f, r, a } => Some(f.world(V3::new(0.0, 0.0, -r / a.tan()))),
            _ => None,
        }
    }

    /// The foot of p on an analytic surface: its unique local minimum of distance.
    /// On a degenerate position (on the axis, at the center) any closest point.
    /// For a cone, the nappe r + v tan a >= 0 (the apex when p is behind it).
    pub fn foot(&self, p: V3) -> V3 {
        match self {
            Surf::Plane { f } => p - f.z * (p - f.o).dot(f.z),
            Surf::Cylinder { f, r } => {
                let ax = axial(f, p);
                f.o + ax.er * *r + f.z * ax.h
            }
            Surf::Cone { f, r, a } => {
                let ax = axial(f, p);
                let (sa, ca) = a.sin_cos();
                let ha = -r / a.tan();
                let s = ax.rho * sa + (ax.h - ha) * ca;
                if s <= 0.0 {
                    f.o + f.z * ha
                } else {
                    f.o + ax.er * (s * sa) + f.z * (ha + s * ca)
                }
            }
            Surf::Sphere { f, r } => {
                let d = p - f.o;
                let l = d.len();
                if l > 0.0 {
                    f.o + d * (r / l)
                } else {
                    f.o + f.x * *r
                }
            }
            Surf::Torus { f, major, minor } => {
                let ax = axial(f, p);
                let c = f.o + ax.er * *major;
                let d = p - c;
                let l = d.len();
                if l > 0.0 {
                    c + d * (minor / l)
                } else {
                    c + ax.er * *minor
                }
            }
            Surf::Nurbs(n) => n.closest(p, f64::INFINITY).map(|h| h.q).unwrap_or(p),
        }
    }

    /// The surface's own unit normal at the surface point nearest x.
    pub fn normal(&self, x: V3) -> V3 {
        match self {
            Surf::Plane { f } => f.z,
            Surf::Cylinder { f, .. } => axial(f, x).er,
            Surf::Cone { f, a, .. } => {
                let ax = axial(f, x);
                let (sa, ca) = a.sin_cos();
                ax.er * ca - f.z * sa
            }
            Surf::Sphere { f, .. } => (x - f.o).unit(),
            Surf::Torus { f, major, .. } => {
                let ax = axial(f, x);
                (x - (f.o + ax.er * *major)).unit()
            }
            Surf::Nurbs(n) => n
                .closest(x, f64::INFINITY)
                .map(|h| h.n)
                .unwrap_or(V3::new(0.0, 0.0, 1.0)),
        }
    }

    /// The foot of p and the surface's own unit normal there, from one search.
    pub fn foot_normal(&self, p: V3) -> (V3, V3) {
        match self {
            Surf::Nurbs(n) => n
                .closest(p, f64::INFINITY)
                .map(|h| (h.q, h.n))
                .unwrap_or((p, V3::new(0.0, 0.0, 1.0))),
            _ => {
                let q = self.foot(p);
                (q, self.normal(q))
            }
        }
    }

    /// Distance from x to the untrimmed surface (for tolerance checks).
    pub fn dist(&self, x: V3) -> f64 {
        (x - self.foot(x)).len()
    }
}

// ---------------------------------------------------------------- B-spline surfaces

/// Bernstein basis of degree n at t with first and second derivatives.
fn bernstein(n: usize, t: f64) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut levels: Vec<Vec<f64>> = vec![vec![1.0]];
    for k in 1..=n {
        let prev = &levels[k - 1];
        let mut b = vec![0.0; k + 1];
        for i in 0..=k {
            let left = if i > 0 { prev[i - 1] } else { 0.0 };
            let right = if i < k { prev[i] } else { 0.0 };
            b[i] = (1.0 - t) * right + t * left;
        }
        levels.push(b);
    }
    let get = |lv: &Vec<f64>, i: isize| {
        if i < 0 || i as usize >= lv.len() {
            0.0
        } else {
            lv[i as usize]
        }
    };
    let mut d1 = vec![0.0; n + 1];
    let mut d2 = vec![0.0; n + 1];
    for i in 0..=n {
        let ii = i as isize;
        if n >= 1 {
            let l = &levels[n - 1];
            d1[i] = n as f64 * (get(l, ii - 1) - get(l, ii));
        }
        if n >= 2 {
            let l = &levels[n - 2];
            d2[i] = (n * (n - 1)) as f64 * (get(l, ii - 2) - 2.0 * get(l, ii - 1) + get(l, ii));
        }
    }
    (levels.pop().unwrap(), d1, d2)
}

/// Newton iterations taken on patches (for profiling).
pub static NEWTON_ITERS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

/// Degrees below this evaluate on the stack.
const BMAX: usize = 16;

/// Bernstein basis of degree n < BMAX at t with first and second derivatives,
/// without allocating.
fn bernstein_small(n: usize, t: f64) -> ([f64; BMAX], [f64; BMAX], [f64; BMAX]) {
    // Levels n - 2 and n - 1 of the triangle are kept for the derivatives.
    let (mut cur, mut l1, mut l2) = ([0.0f64; BMAX], [0.0f64; BMAX], [0.0f64; BMAX]);
    cur[0] = 1.0;
    if n >= 1 && n - 1 == 0 {
        l1 = cur;
    }
    if n >= 2 && n - 2 == 0 {
        l2 = cur;
    }
    for k in 1..=n {
        for i in (0..=k).rev() {
            let left = if i > 0 { cur[i - 1] } else { 0.0 };
            let right = if i < k { cur[i] } else { 0.0 };
            cur[i] = (1.0 - t) * right + t * left;
        }
        if n >= 1 && k == n - 1 {
            l1 = cur;
        }
        if n >= 2 && k == n - 2 {
            l2 = cur;
        }
    }
    let get = |lv: &[f64; BMAX], len: usize, i: isize| {
        if i < 0 || i as usize >= len {
            0.0
        } else {
            lv[i as usize]
        }
    };
    let (mut d1, mut d2) = ([0.0f64; BMAX], [0.0f64; BMAX]);
    for i in 0..=n {
        let ii = i as isize;
        if n >= 1 {
            d1[i] = n as f64 * (get(&l1, n, ii - 1) - get(&l1, n, ii));
        }
        if n >= 2 {
            d2[i] = (n * (n - 1)) as f64
                * (get(&l2, n - 1, ii - 2) - 2.0 * get(&l2, n - 1, ii - 1) + get(&l2, n - 1, ii));
        }
    }
    (cur, d1, d2)
}

/// Point, first and second derivatives of a surface at local (s, t).
#[derive(Clone, Copy, Debug)]
pub struct SEval {
    pub p: V3,
    pub ps: V3,
    pub pt: V3,
    pub pss: V3,
    pub pst: V3,
    pub ptt: V3,
}

/// A rational Bezier patch in homogeneous form: `pw[i][j]`, i along u.
#[derive(Clone, Debug)]
pub struct Patch {
    pub pw: Vec<Vec<[f64; 4]>>,
    pub u0: f64,
    pub u1: f64,
    pub v0: f64,
    pub v1: f64,
}

fn v3h(h: [f64; 4]) -> V3 {
    V3::new(h[0], h[1], h[2])
}

fn casteljau_split(pts: &[[f64; 4]], s: f64) -> (Vec<[f64; 4]>, Vec<[f64; 4]>) {
    let n = pts.len() - 1;
    let mut b = pts.to_vec();
    let (mut left, mut right) = (vec![b[0]], vec![b[n]]);
    for level in 1..=n {
        for i in 0..=(n - level) {
            for c in 0..4 {
                b[i][c] = (1.0 - s) * b[i][c] + s * b[i + 1][c];
            }
        }
        left.push(b[0]);
        right.push(b[n - level]);
    }
    right.reverse();
    (left, right)
}

impl Patch {
    pub fn ctrl(&self) -> impl Iterator<Item = V3> + '_ {
        self.pw
            .iter()
            .flatten()
            .map(|h| V3::new(h[0] / h[3], h[1] / h[3], h[2] / h[3]))
    }
    pub fn bbox(&self) -> (V3, V3) {
        let mut lo = V3::new(f64::MAX, f64::MAX, f64::MAX);
        let mut hi = V3::new(f64::MIN, f64::MIN, f64::MIN);
        for q in self.ctrl() {
            lo = V3::new(lo.x.min(q.x), lo.y.min(q.y), lo.z.min(q.z));
            hi = V3::new(hi.x.max(q.x), hi.y.max(q.y), hi.z.max(q.z));
        }
        (lo, hi)
    }
    pub fn eval(&self, s: f64, t: f64) -> SEval {
        let (du, dv) = (self.pw.len() - 1, self.pw[0].len() - 1);
        if du < BMAX && dv < BMAX {
            return self.eval_small(s, t);
        }
        let (bu, bu1, bu2) = bernstein(du, s);
        let (bv, bv1, bv2) = bernstein(dv, t);
        let mut a = [[0.0f64; 4]; 6]; // A, As, At, Ass, Ast, Att
        for i in 0..=du {
            for j in 0..=dv {
                let h = self.pw[i][j];
                let w = [
                    bu[i] * bv[j],
                    bu1[i] * bv[j],
                    bu[i] * bv1[j],
                    bu2[i] * bv[j],
                    bu1[i] * bv1[j],
                    bu[i] * bv2[j],
                ];
                for k in 0..6 {
                    for c in 0..4 {
                        a[k][c] += w[k] * h[c];
                    }
                }
            }
        }
        let w = a[0][3];
        let p = v3h(a[0]) * (1.0 / w);
        let ps = (v3h(a[1]) - p * a[1][3]) * (1.0 / w);
        let pt = (v3h(a[2]) - p * a[2][3]) * (1.0 / w);
        let pss = (v3h(a[3]) - ps * (2.0 * a[1][3]) - p * a[3][3]) * (1.0 / w);
        let pst = (v3h(a[4]) - ps * a[2][3] - pt * a[1][3] - p * a[4][3]) * (1.0 / w);
        let ptt = (v3h(a[5]) - pt * (2.0 * a[2][3]) - p * a[5][3]) * (1.0 / w);
        SEval {
            p,
            ps,
            pt,
            pss,
            pst,
            ptt,
        }
    }
    /// eval for degrees below BMAX, with the basis on the stack.
    fn eval_small(&self, s: f64, t: f64) -> SEval {
        let (du, dv) = (self.pw.len() - 1, self.pw[0].len() - 1);
        let (bu, bu1, bu2) = bernstein_small(du, s);
        let (bv, bv1, bv2) = bernstein_small(dv, t);
        let mut a = [[0.0f64; 4]; 6]; // A, As, At, Ass, Ast, Att
        for i in 0..=du {
            for j in 0..=dv {
                let h = self.pw[i][j];
                let w = [
                    bu[i] * bv[j],
                    bu1[i] * bv[j],
                    bu[i] * bv1[j],
                    bu2[i] * bv[j],
                    bu1[i] * bv1[j],
                    bu[i] * bv2[j],
                ];
                for k in 0..6 {
                    for c in 0..4 {
                        a[k][c] += w[k] * h[c];
                    }
                }
            }
        }
        let w = a[0][3];
        let p = v3h(a[0]) * (1.0 / w);
        let ps = (v3h(a[1]) - p * a[1][3]) * (1.0 / w);
        let pt = (v3h(a[2]) - p * a[2][3]) * (1.0 / w);
        let pss = (v3h(a[3]) - ps * (2.0 * a[1][3]) - p * a[3][3]) * (1.0 / w);
        let pst = (v3h(a[4]) - ps * a[2][3] - pt * a[1][3] - p * a[4][3]) * (1.0 / w);
        let ptt = (v3h(a[5]) - pt * (2.0 * a[2][3]) - p * a[5][3]) * (1.0 / w);
        SEval {
            p,
            ps,
            pt,
            pss,
            pst,
            ptt,
        }
    }

    pub fn split_u(&self, s: f64) -> (Patch, Patch) {
        let (du, dv) = (self.pw.len() - 1, self.pw[0].len() - 1);
        let mut l = vec![vec![[0.0; 4]; dv + 1]; du + 1];
        let mut r = l.clone();
        for j in 0..=dv {
            let col: Vec<[f64; 4]> = (0..=du).map(|i| self.pw[i][j]).collect();
            let (a, b) = casteljau_split(&col, s);
            for i in 0..=du {
                l[i][j] = a[i];
                r[i][j] = b[i];
            }
        }
        let um = self.u0 + s * (self.u1 - self.u0);
        (
            Patch {
                pw: l,
                u1: um,
                ..self.clone()
            },
            Patch {
                pw: r,
                u0: um,
                ..self.clone()
            },
        )
    }
    pub fn split_v(&self, t: f64) -> (Patch, Patch) {
        let (mut l, mut r) = (Vec::new(), Vec::new());
        for row in &self.pw {
            let (a, b) = casteljau_split(row, t);
            l.push(a);
            r.push(b);
        }
        let vm = self.v0 + t * (self.v1 - self.v0);
        (
            Patch {
                pw: l,
                v1: vm,
                ..self.clone()
            },
            Patch {
                pw: r,
                v0: vm,
                ..self.clone()
            },
        )
    }
    /// Newton on |S - p|^2 over [0, 1]^2 from (s, t); Gauss-Newton where the
    /// Hessian is not positive definite.
    pub fn newton(&self, p: V3, s: f64, t: f64) -> (f64, f64) {
        self.newton_in(p, s, t, (0.0, 1.0), (0.0, 1.0))
    }

    /// Newton as above, with s and t kept in the given ranges; ranges past [0, 1]
    /// continue the patch's polynomial exactly.
    pub fn newton_in(&self, p: V3, s: f64, t: f64, sr: (f64, f64), tr: (f64, f64)) -> (f64, f64) {
        self.newton_iters(p, s, t, sr, tr, 50)
    }

    /// Newton over [0, 1]^2 from (s, t) with at most `iters` steps: enough for a
    /// branch and bound, which only needs a good point (its bounds certify).
    pub fn newton_capped(&self, p: V3, s: f64, t: f64, iters: usize) -> (f64, f64) {
        self.newton_iters(p, s, t, (0.0, 1.0), (0.0, 1.0), iters)
    }

    fn newton_iters(
        &self,
        p: V3,
        mut s: f64,
        mut t: f64,
        sr: (f64, f64),
        tr: (f64, f64),
        iters: usize,
    ) -> (f64, f64) {
        // A parameter sitting on its bound with its step pointing out is held
        // there, and Newton continues in the other one alone (the minimum along
        // that side of the patch); clamping both each step would bounce off the
        // bound until the iteration cap.
        let (mut hold_s, mut hold_t) = (false, false);
        for _ in 0..iters {
            NEWTON_ITERS.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            let e = self.eval(s, t);
            let r = e.p - p;
            let g = [e.ps.dot(r), e.pt.dot(r)];
            let (jss, jst, jtt) = (e.ps.dot(e.ps), e.ps.dot(e.pt), e.pt.dot(e.pt));
            let mut h = [jss + e.pss.dot(r), jst + e.pst.dot(r), jtt + e.ptt.dot(r)];
            let (ds, dt) = match (hold_s, hold_t) {
                (true, true) => break,
                (true, false) => (0.0, -g[1] / if h[2] > 0.0 { h[2] } else { jtt.max(1e-300) }),
                (false, true) => (-g[0] / if h[0] > 0.0 { h[0] } else { jss.max(1e-300) }, 0.0),
                (false, false) => {
                    if !(h[0] > 0.0 && h[0] * h[2] - h[1] * h[1] > 0.0) {
                        h = [jss, jst, jtt];
                    }
                    let det = h[0] * h[2] - h[1] * h[1];
                    if det.abs() < 1e-300 {
                        break;
                    }
                    (
                        -(h[2] * g[0] - h[1] * g[1]) / det,
                        -(h[0] * g[1] - h[1] * g[0]) / det,
                    )
                }
            };
            let at_bound =
                |x: f64, d: f64, r: (f64, f64)| (x <= r.0 && d < 0.0) || (x >= r.1 && d > 0.0);
            if !hold_s && at_bound(s, ds, sr) {
                hold_s = true;
                continue;
            }
            if !hold_t && at_bound(t, dt, tr) {
                hold_t = true;
                continue;
            }
            let (ns, nt) = ((s + ds).clamp(sr.0, sr.1), (t + dt).clamp(tr.0, tr.1));
            let moved = (ns - s).abs() + (nt - t).abs();
            s = ns;
            t = nt;
            if moved < 1e-13 {
                break;
            }
        }
        (s, t)
    }
}

/// Bezier points of the single knot span `knots[p] .. knots[p + 1]` of a local
/// B-spline piece (2p + 2 knots, p + 1 points): raise both span ends to
/// multiplicity p by knot insertion; the p + 1 points ending at the last copy of
/// the span's start are then its Bezier points.
pub fn span_bezier<const N: usize>(p: usize, knots: &[f64], pts: &[[f64; N]]) -> Vec<[f64; N]> {
    let (a, b) = (knots[p], knots[p + 1]);
    let (mut kn, mut pt) = (knots.to_vec(), pts.to_vec());
    while kn.iter().filter(|&&x| x == a).count() < p {
        insert_knot(p, &mut kn, &mut pt, a);
    }
    while kn.iter().filter(|&&x| x == b).count() < p {
        insert_knot(p, &mut kn, &mut pt, b);
    }
    let last_a = kn.iter().rposition(|&x| x == a).unwrap();
    pt[last_a - p..=last_a].to_vec()
}

/// Node of the quadtree over a B-spline surface's grid of knot spans.
#[derive(Clone, Debug)]
pub struct QNode {
    pub lo: V3,
    pub hi: V3,
    /// Span index ranges [i0, i1) along u and [j0, j1) along v.
    pub i0: u32,
    pub i1: u32,
    pub j0: u32,
    pub j1: u32,
    pub kids: Vec<u32>,
}

/// A clamped rational B-spline surface, kept as its control net. Bezier patches
/// are extracted one knot span at a time when a search reaches them (a surface
/// like the brain coral's, 363 thousand control points at degree 5, would need
/// hundreds of megabytes as patches). Each span is bounded by the box of the
/// control points that act on it (convex hull property, positive weights), and a
/// quadtree over the span grid carries the union boxes.
#[derive(Clone, Debug)]
pub struct NurbsSurf {
    pub du: usize,
    pub dv: usize,
    nv: usize,
    /// Homogeneous control net, point (i, j) at i * nv + j.
    pw: Vec<[f64; 4]>,
    ku: Vec<f64>,
    kv: Vec<f64>,
    /// Knot index s of each nonempty span ku[s] .. ku[s + 1] (and alike in v).
    su: Vec<usize>,
    sv: Vec<usize>,
    pub nodes: Vec<QNode>,
    pub u_range: (f64, f64),
    pub v_range: (f64, f64),
}

/// A point of a B-spline surface nearest a query.
#[derive(Clone, Copy, Debug)]
pub struct SurfHit {
    pub q: V3,
    pub dist: f64,
    pub u: f64,
    pub v: f64,
    /// Unit normal Su x Sv at the hit.
    pub n: V3,
}

/// What a branch and bound entry holds: a quadtree node, or an extracted patch.
pub enum Entry {
    Node(u32),
    Patch(Patch),
}

/// Min-heap entry for branch and bound over a B-spline surface.
pub struct PatchItem {
    pub lb: f64,
    pub entry: Entry,
    pub depth: u32,
    /// Known to lie wholly inside its trimmed face.
    pub inside: bool,
}
impl PartialEq for PatchItem {
    fn eq(&self, o: &Self) -> bool {
        self.lb == o.lb
    }
}
impl Eq for PatchItem {}
impl PartialOrd for PatchItem {
    fn partial_cmp(&self, o: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for PatchItem {
    fn cmp(&self, o: &Self) -> std::cmp::Ordering {
        o.lb.total_cmp(&self.lb)
    }
}

/// Newton steps per patch inside a branch and bound search.
pub const SEARCH_NEWTON: usize = 12;

/// Distances below this (mm) are treated as settled by the branch and bound searches.
pub const DIST_EPS: f64 = 1e-9;

pub fn box_dist(p: V3, b: &(V3, V3)) -> f64 {
    V3::new(
        (b.0.x - p.x).max(p.x - b.1.x).max(0.0),
        (b.0.y - p.y).max(p.y - b.1.y).max(0.0),
        (b.0.z - p.z).max(p.z - b.1.z).max(0.0),
    )
    .len()
}

/// Lower bound on the distance from p to a patch along unit direction n: the
/// patch lies in the hull of its control points, and every hull point y has
/// (p - y).n at least the least (p - c).n over the control points. Taken along
/// the direction to the patch's Newton point, it tightens quadratically as the
/// patch flattens, where a box bound only tightens linearly. It also holds for
/// every sub-patch (their control hulls lie in the parent's).
pub fn slab_bound(pa: &Patch, p: V3, n: V3) -> f64 {
    pa.ctrl()
        .map(|c| (p - c).dot(n))
        .fold(f64::INFINITY, f64::min)
}

impl Patch {
    /// The four quarters of the patch (split at the middle in u and in v).
    pub fn quarters(&self) -> [Patch; 4] {
        let (a, b) = self.split_u(0.5);
        let (a0, a1) = a.split_v(0.5);
        let (b0, b1) = b.split_v(0.5);
        [a0, a1, b0, b1]
    }
}

impl NurbsSurf {
    /// Rows `ctrl[i][j]` with i along u, weights alike, clamped knot vectors.
    pub fn new(
        du: usize,
        dv: usize,
        ctrl: &[Vec<V3>],
        weights: Option<&Vec<Vec<f64>>>,
        ku: &[f64],
        kv: &[f64],
    ) -> Result<NurbsSurf, String> {
        let nu = ctrl.len();
        let nv = ctrl.first().map_or(0, |r| r.len());
        if nu + du + 1 != ku.len() || nv + dv + 1 != kv.len() || ctrl.iter().any(|r| r.len() != nv)
        {
            return Err(format!("B-spline surface: {nu} x {nv} control points do not match degrees {du}, {dv} and {} + {} knots", ku.len(), kv.len()));
        }
        let w = |i: usize, j: usize| weights.map_or(1.0, |w| w[i][j]);
        if (0..nu).any(|i| (0..nv).any(|j| !(w(i, j) > 0.0))) {
            return Err("B-spline surface with a weight that is not positive".into());
        }
        // Homogeneous net, clamped exactly in each direction (knot insertion leaves
        // the surface unchanged; every row takes the same insertions).
        let mut net: Vec<Vec<[f64; 4]>> = (0..nu)
            .map(|i| {
                (0..nv)
                    .map(|j| {
                        let (c, wt) = (ctrl[i][j], w(i, j));
                        [c.x * wt, c.y * wt, c.z * wt, wt]
                    })
                    .collect()
            })
            .collect();
        let mut kv2 = kv.to_vec();
        for row in net.iter_mut() {
            let (k, r) = clamp(dv, kv, row);
            kv2 = k;
            *row = r;
        }
        let nv = kv2.len() - dv - 1;
        let mut ku2 = ku.to_vec();
        let mut cols: Vec<Vec<[f64; 4]>> = (0..nv)
            .map(|j| net.iter().map(|r| r[j]).collect())
            .collect();
        for col in cols.iter_mut() {
            let (k, c) = clamp(du, ku, col);
            ku2 = k;
            *col = c;
        }
        let nu = ku2.len() - du - 1;
        let (ku, kv) = (&ku2[..], &kv2[..]);
        let mut pw = Vec::with_capacity(nu * nv);
        for i in 0..nu {
            for j in 0..nv {
                pw.push(cols[j][i]);
            }
        }
        let spans = |k: &[f64], d: usize| -> Vec<usize> {
            (d..k.len() - d - 1).filter(|&s| k[s] < k[s + 1]).collect()
        };
        let (su, sv) = (spans(ku, du), spans(kv, dv));
        let mut ns = NurbsSurf {
            du,
            dv,
            nv,
            pw,
            ku: ku.to_vec(),
            kv: kv.to_vec(),
            su,
            sv,
            nodes: Vec::new(),
            u_range: (ku[du], ku[ku.len() - du - 1]),
            v_range: (kv[dv], kv[kv.len() - dv - 1]),
        };
        ns.build(0, ns.su.len() as u32, 0, ns.sv.len() as u32);
        Ok(ns)
    }

    pub fn span_count(&self) -> usize {
        self.su.len() * self.sv.len()
    }

    /// Spans along u and along v.
    pub fn spans_uv(&self) -> (usize, usize) {
        (self.su.len(), self.sv.len())
    }

    /// The span indices holding (u, v).
    pub fn span_at(&self, u: f64, v: f64) -> (usize, usize) {
        (
            Self::span_of(&self.su, &self.ku, u),
            Self::span_of(&self.sv, &self.kv, v),
        )
    }

    /// Box of the control points acting on span (a, b).
    fn span_box(&self, a: usize, b: usize) -> (V3, V3) {
        let (s, t) = (self.su[a], self.sv[b]);
        let mut lo = V3::new(f64::MAX, f64::MAX, f64::MAX);
        let mut hi = V3::new(f64::MIN, f64::MIN, f64::MIN);
        for i in s - self.du..=s {
            for j in t - self.dv..=t {
                let h = self.pw[i * self.nv + j];
                let q = V3::new(h[0] / h[3], h[1] / h[3], h[2] / h[3]);
                lo = V3::new(lo.x.min(q.x), lo.y.min(q.y), lo.z.min(q.z));
                hi = V3::new(hi.x.max(q.x), hi.y.max(q.y), hi.z.max(q.z));
            }
        }
        (lo, hi)
    }

    fn build(&mut self, i0: u32, i1: u32, j0: u32, j1: u32) -> u32 {
        let id = self.nodes.len() as u32;
        self.nodes.push(QNode {
            lo: V3::default(),
            hi: V3::default(),
            i0,
            i1,
            j0,
            j1,
            kids: Vec::new(),
        });
        let (lo, hi, kids) = if i1 - i0 == 1 && j1 - j0 == 1 {
            let (lo, hi) = self.span_box(i0 as usize, j0 as usize);
            (lo, hi, Vec::new())
        } else {
            let (im, jm) = ((i0 + i1) / 2, (j0 + j1) / 2);
            let mut kids = Vec::new();
            let iranges: Vec<(u32, u32)> = if i1 - i0 > 1 {
                vec![(i0, im), (im, i1)]
            } else {
                vec![(i0, i1)]
            };
            let jranges: Vec<(u32, u32)> = if j1 - j0 > 1 {
                vec![(j0, jm), (jm, j1)]
            } else {
                vec![(j0, j1)]
            };
            for &(a, b) in &iranges {
                for &(c, d) in &jranges {
                    kids.push(self.build(a, b, c, d));
                }
            }
            let mut lo = V3::new(f64::MAX, f64::MAX, f64::MAX);
            let mut hi = V3::new(f64::MIN, f64::MIN, f64::MIN);
            for &k in &kids {
                let n = &self.nodes[k as usize];
                lo = V3::new(lo.x.min(n.lo.x), lo.y.min(n.lo.y), lo.z.min(n.lo.z));
                hi = V3::new(hi.x.max(n.hi.x), hi.y.max(n.hi.y), hi.z.max(n.hi.z));
            }
            (lo, hi, kids)
        };
        let n = &mut self.nodes[id as usize];
        n.lo = lo;
        n.hi = hi;
        n.kids = kids;
        id
    }

    /// Parameter box of a node: (u0, u1, v0, v1).
    pub fn node_params(&self, n: &QNode) -> (f64, f64, f64, f64) {
        (
            self.ku[self.su[n.i0 as usize]],
            self.ku[self.su[n.i1 as usize - 1] + 1],
            self.kv[self.sv[n.j0 as usize]],
            self.kv[self.sv[n.j1 as usize - 1] + 1],
        )
    }

    /// The Bezier patch of span (a, b), by local knot insertion in u then v.
    pub fn patch(&self, a: usize, b: usize) -> Patch {
        let (s, t) = (self.su[a], self.sv[b]);
        let (du, dv) = (self.du, self.dv);
        let lu = &self.ku[s - du..=s + du + 1];
        let lv = &self.kv[t - dv..=t + dv + 1];
        // Along u for each local column j, then along v for each Bezier row i.
        let cols: Vec<Vec<[f64; 4]>> = (t - dv..=t)
            .map(|j| {
                span_bezier(
                    du,
                    lu,
                    &(s - du..=s)
                        .map(|i| self.pw[i * self.nv + j])
                        .collect::<Vec<_>>(),
                )
            })
            .collect();
        let pw: Vec<Vec<[f64; 4]>> = (0..=du)
            .map(|i| span_bezier(dv, lv, &cols.iter().map(|c| c[i]).collect::<Vec<_>>()))
            .collect();
        Patch {
            pw,
            u0: self.ku[s],
            u1: self.ku[s + 1],
            v0: self.kv[t],
            v1: self.kv[t + 1],
        }
    }

    /// The Bezier patch of the span holding (u, v).
    pub fn patch_at(&self, u: f64, v: f64) -> Patch {
        self.patch(
            Self::span_of(&self.su, &self.ku, u),
            Self::span_of(&self.sv, &self.kv, v),
        )
    }

    /// Index of the span holding parameter x (the last span for the far end).
    fn span_of(spans: &[usize], k: &[f64], x: f64) -> usize {
        spans.iter().rposition(|&s| k[s] <= x).unwrap_or_default()
    }

    /// Point and derivatives at global parameters (u, v); derivatives per unit u and v.
    pub fn eval(&self, u: f64, v: f64) -> SEval {
        let (a, b) = (
            Self::span_of(&self.su, &self.ku, u),
            Self::span_of(&self.sv, &self.kv, v),
        );
        let pa = self.patch(a, b);
        let (s, t) = ((u - pa.u0) / (pa.u1 - pa.u0), (v - pa.v0) / (pa.v1 - pa.v0));
        let e = pa.eval(s.clamp(0.0, 1.0), t.clamp(0.0, 1.0));
        let (x, y) = (1.0 / (pa.u1 - pa.u0), 1.0 / (pa.v1 - pa.v0));
        SEval {
            p: e.p,
            ps: e.ps * x,
            pt: e.pt * y,
            pss: e.pss * (x * x),
            pst: e.pst * (x * y),
            ptt: e.ptt * (y * y),
        }
    }

    /// Box holding the whole surface.
    pub fn bbox(&self) -> (V3, V3) {
        (self.nodes[0].lo, self.nodes[0].hi)
    }

    /// Certified closest point of the untrimmed surface: branch and bound down the
    /// span quadtree (box bounds) into extracted patches (box and support-plane
    /// bounds, Newton at each popped patch).
    pub fn closest(&self, p: V3, cutoff: f64) -> Option<SurfHit> {
        let mut best: Option<SurfHit> = None;
        let mut best_d = cutoff;
        let mut heap = std::collections::BinaryHeap::new();
        heap.push(PatchItem {
            lb: box_dist(p, &self.bbox()),
            entry: Entry::Node(0),
            depth: 0,
            inside: false,
        });
        let mut pops = 0;
        while let Some(PatchItem {
            lb, entry, depth, ..
        }) = heap.pop()
        {
            if lb >= best_d - DIST_EPS || best_d <= DIST_EPS {
                break;
            }
            pops += 1;
            if pops > 20000 {
                break;
            }
            let pa = match entry {
                Entry::Node(id) => {
                    let n = &self.nodes[id as usize];
                    if n.kids.is_empty() {
                        let pa = self.patch(n.i0 as usize, n.j0 as usize);
                        let lbp = box_dist(p, &pa.bbox()).max(lb);
                        if lbp < best_d - DIST_EPS {
                            heap.push(PatchItem {
                                lb: lbp,
                                entry: Entry::Patch(pa),
                                depth: 0,
                                inside: false,
                            });
                        }
                    } else {
                        for &k in &n.kids {
                            let c = &self.nodes[k as usize];
                            let lbk = box_dist(p, &(c.lo, c.hi));
                            if lbk < best_d - DIST_EPS {
                                heap.push(PatchItem {
                                    lb: lbk,
                                    entry: Entry::Node(k),
                                    depth: 0,
                                    inside: false,
                                });
                            }
                        }
                    }
                    continue;
                }
                Entry::Patch(pa) => pa,
            };
            let (s, t) = pa.newton_capped(p, 0.5, 0.5, SEARCH_NEWTON);
            let e = pa.eval(s, t);
            let d = (e.p - p).len();
            let n = e.ps.cross(e.pt).unit();
            if d < best_d {
                best_d = d;
                best = Some(SurfHit {
                    q: e.p,
                    dist: d,
                    u: pa.u0 + s * (pa.u1 - pa.u0),
                    v: pa.v0 + t * (pa.v1 - pa.v0),
                    n,
                });
            }
            let slab = slab_bound(&pa, p, if d > 1e-12 { (p - e.p) * (1.0 / d) } else { n });
            if slab >= best_d - DIST_EPS || depth >= 40 {
                continue;
            }
            for q in pa.quarters() {
                let lbq = box_dist(p, &q.bbox()).max(slab);
                if lbq < best_d - DIST_EPS {
                    heap.push(PatchItem {
                        lb: lbq,
                        entry: Entry::Patch(q),
                        depth: depth + 1,
                        inside: false,
                    });
                }
            }
        }
        best
    }

    /// Average span size (for sizing boundary boxes against patches).
    pub fn span_size(&self) -> f64 {
        let (lo, hi) = self.bbox();
        (hi - lo).len() / (self.span_count() as f64).sqrt().max(1.0)
    }

    /// True when the surface closes on itself in u (or v): its first and last
    /// iso-curves coincide.
    pub fn closed(&self, in_u: bool, tol: f64) -> bool {
        let (a, b) = if in_u { self.u_range } else { self.v_range };
        (0..=8).all(|k| {
            let s = k as f64 / 8.0;
            let (p, q) = if in_u {
                let v = self.v_range.0 + s * (self.v_range.1 - self.v_range.0);
                (self.eval(a, v).p, self.eval(b, v).p)
            } else {
                let u = self.u_range.0 + s * (self.u_range.1 - self.u_range.0);
                (self.eval(u, a).p, self.eval(u, b).p)
            };
            (p - q).len() < tol
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn frame() -> Frame {
        Frame {
            o: V3::new(1.0, 2.0, 3.0),
            z: V3::new(0.0, 0.6, 0.8),
            x: V3::new(1.0, 0.0, 0.0),
        }
    }

    /// Every analytic foot is a critical point of the distance (p - q parallel to the
    /// normal) and lies on its surface.
    #[test]
    fn analytic_feet_are_critical_points() {
        let f = frame();
        let surfs = [
            Surf::Plane { f },
            Surf::Cylinder { f, r: 2.0 },
            Surf::Cone { f, r: 2.0, a: 0.5 },
            Surf::Sphere { f, r: 2.0 },
            Surf::Torus {
                f,
                major: 5.0,
                minor: 1.5,
            },
        ];
        for s in &surfs {
            for k in 0..200 {
                let x = k as f64;
                let p = V3::new(
                    1.0 + 7.0 * (x * 0.37).sin(),
                    2.0 + 6.0 * (x * 0.73).cos(),
                    3.0 + 5.0 * (x * 1.31).sin(),
                );
                let q = s.foot(p);
                let n = s.normal(q);
                let d = p - q;
                if let Surf::Cone { .. } = s {
                    if (q - s.apex().unwrap()).len() < 1e-9 {
                        continue;
                    }
                }
                assert!(
                    d.cross(n).len() < 1e-9 * (1.0 + d.len()),
                    "{}: p - q not normal at {p:?}",
                    s.name()
                );
                assert!(s.dist(q) < 1e-12, "{}: foot off the surface", s.name());
            }
        }
    }

    /// A bicubic patch surface: the certified closest point is no farther than any
    /// point of a dense parameter sweep.
    #[test]
    fn nurbs_closest_beats_dense_sampling() {
        let mut ctrl = Vec::new();
        for i in 0..5 {
            let mut row = Vec::new();
            for j in 0..4 {
                let (x, y) = (i as f64 * 2.0, j as f64 * 2.0);
                row.push(V3::new(x, y, (x * 0.7).sin() * 2.0 + (y * 0.9).cos()));
            }
            ctrl.push(row);
        }
        let ku = [0.0, 0.0, 0.0, 0.0, 0.4, 1.0, 1.0, 1.0, 1.0];
        let kv = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
        let s = NurbsSurf::new(3, 3, &ctrl, None, &ku, &kv).unwrap();
        assert_eq!(s.span_count(), 2);
        for k in 0..40 {
            let x = k as f64;
            let p = V3::new(
                4.0 + 5.0 * (x * 0.37).sin(),
                3.0 + 4.0 * (x * 0.71).cos(),
                1.0 + 3.0 * (x * 1.3).sin(),
            );
            let h = s.closest(p, f64::INFINITY).unwrap();
            let mut dense = f64::INFINITY;
            for a in 0..=200 {
                for b in 0..=200 {
                    dense = dense.min((s.eval(a as f64 / 200.0, b as f64 / 200.0).p - p).len());
                }
            }
            assert!(
                h.dist <= dense + 1e-9,
                "closest {} vs dense {}",
                h.dist,
                dense
            );
            assert!(
                h.dist > dense - 0.05,
                "closest {} far below dense {}",
                h.dist,
                dense
            );
        }
    }
}
