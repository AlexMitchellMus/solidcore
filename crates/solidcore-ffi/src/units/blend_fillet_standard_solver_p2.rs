//! FFI wrappers for unit `blend.fillet_standard_solver_p2`.
//!
//! Shim: `harness/cpp/shim/blend_fillet_standard_solver_p2.cpp`
//! (`sc_blend_fillet_standard_solver_p2_eval`).
//! Upstream: `source/SMLib/src/SmFilletStandardSolver.cpp`.
//!
//! Encodings. A surface is `kind`: `0 n m p q r s`, `(n+1)(m+1)` control points `x y z w`
//! (v index fastest), `r+1` u knots, `s+1` v knots (`SmBSplineSurface(gw_SURFACE*, FALSE,
//! NULL)`), or `1 o[3] x[3] y[3] sx sy u0 v0 u1 v1` (`SmPlane(o, x, y, (sx, sy), domain,
//! ctx)`). A solver block is `cls same surfaceA [surfaceB] r0 r1 o0 o1 ext`: `cls` `0`
//! `SmConstantRadiusFS(ctx, 1e-6, 30deg, 2deg, r0, r1, A, B or A, o0, o1)`, `1`
//! `SmSurfaceSurfaceFS` (the same, UV guesses `(0.5, 0.5)`) or `2` an
//! `SmConstantRadiusAssistedFS` built on an `SmApiCreateBox` edge use whose rail `i`
//! offsets its base by `r_i` times `+-1` (`o_i`), with radii `r_i` and no edge uses; `ext`
//! is the extension factor. A surface selector is `1` `GetSurface(0)`, `2`
//! `GetSurface(1)`, `3` another offset surface of A, `4` an offset surface of a copy of A.
//! `kinds` is `nk` and `nk` type codes passed to `IsKindOf`.
//!
//! | [`Op`]   | upstream                                                  | parameters |
//! |----------|-----------------------------------------------------------|------------|
//! | `CrInit` | `SmConstantRadiusFS::LoadInitialValues`                   | `solver do2 rail sel1 s1_index dim` `dim` x `min max`, `np per[np] ng guess[ng]` |
//! | `Guess`  | `SmConstantRadiusAssistedFS::GuessPointOnPlaneSolve`      | `solver(cls 2) orig[3] normal[3] d1[4] d2[4] uv1[2] uv2[2]` |
//! | `Plane`  | `SmConstantRadiusAssistedFS::PointOnPlaneSolve`           | as `Guess` |
//! | `Ss`     | `SmSurfaceSurfaceFS` constructor (and `CalcFilletGeom`)   | `surfaceA surfaceB r0 r1 o0 o1 g0[2] g1[2] calc kinds` |
//! | `Types`  | `IsKindOf` of the assisted, surface/surface and variable-radius solvers | `kinds` |
//!
//! Output: the op status (`0` = `SM_SUCCESS`, else the `SmStatus`) and then: `CrInit` the
//! rail of `rpSrf2Offset` (`-1` `NULL`), `rlSurface2Index`, the intervals, `np' per[np']`,
//! `ng' guess[ng']` and both offset distances; `Guess` on success `found` and with
//! `found` both UV points, then both offset distances; `Plane` as `Guess`, after the UV points
//! per surface its point, derivatives and unit normal and then the curve point; `Ss` the
//! radii, orientations, offset distances, UV guesses, trim type, solver type, `GetType`
//! and the kinds (its status is that of `CalcFilletGeom` with `calc`, `-3` when the
//! shim did not run it, `0` without `calc`); `Types` per solver `GetType` and the kinds.
//!
//! [`eval`] accepts only what upstream handles without undefined behaviour: well-formed
//! clamped 3D B-spline surfaces (`1 <= p <= n`, at most [`MAX_INDEX`] control points per
//! direction, degrees at most [`MAX_DEGREE`], finite data, consistent rationality,
//! positive weights, non-decreasing knots with end multiplicity exactly `degree + 1` and
//! interior multiplicity at most the degree), planes with a non-reversed domain (an
//! offset surface's extension of a reversed one reads a half-built NLib surface, UB-187),
//! finite numbers, arrays of at most [`MAX_ARRAY`] entries, parameters within
//! [`MAX_PARAM`], a non-`NULL` first surface of `CrInit` and, with `do2`, the guess and
//! interval indices it reads and writes in range and rail 0 or 1; `Guess` and `Plane` only
//! on the assisted solver (`cls 2`). `CalcFilletGeom` past its parallel-normal exit would
//! step along the fillet and dereference the solver's `NULL` executive: the shim runs it
//! only when the offset normals at the guesses are parallel (its first exit).

use crate::{call_with_buffer, FfiError};

/// Largest control point index accepted by [`eval`].
pub const MAX_INDEX: usize = 64;
/// Largest NURB degree accepted by [`eval`].
pub const MAX_DEGREE: usize = 9;
/// Largest array length (guesses, intervals, periodicities, type codes) accepted by
/// [`eval`].
pub const MAX_ARRAY: usize = 32;
/// Largest `|parameter|` (guesses, UV points) accepted by [`eval`].
pub const MAX_PARAM: f64 = 1.0e6;
/// `NL_NOZ`.
const NOZ: f64 = f64::MAX;
/// `NL_NOW`.
const NOW: f64 = -f64::MAX;

/// The upstream routine run by the shim.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Op {
    /// `SmConstantRadiusFS::LoadInitialValues`.
    CrInit,
    /// `SmConstantRadiusAssistedFS::GuessPointOnPlaneSolve`.
    Guess,
    /// `SmConstantRadiusAssistedFS::PointOnPlaneSolve`.
    Plane,
    /// The `SmSurfaceSurfaceFS` constructor and `CalcFilletGeom`.
    Ss,
    /// `IsKindOf` of the assisted, surface/surface and variable-radius solvers.
    Types,
}

impl Op {
    fn code(self) -> i32 {
        match self {
            Op::CrInit => 0,
            Op::Guess => 1,
            Op::Plane => 2,
            Op::Ss => 3,
            Op::Types => 4,
        }
    }
}

mod sys {
    extern "C" {
        pub fn sc_blend_fillet_standard_solver_p2_eval(
            op: i32,
            input: *const f64,
            n_in: i64,
            out: *mut f64,
            cap: i64,
            out_len: *mut i64,
        ) -> i32;
    }
}

fn bad(msg: impl Into<String>) -> FfiError {
    FfiError::InvalidArgument(msg.into())
}

/// Sequential reader over the input.
struct Reader<'a> {
    input: &'a [f64],
    at: usize,
}

impl Reader<'_> {
    fn num(&mut self) -> Result<f64, FfiError> {
        let v = self
            .input
            .get(self.at)
            .copied()
            .ok_or_else(|| bad("input too short"))?;
        self.at += 1;
        Ok(v)
    }

    fn finite(&mut self) -> Result<f64, FfiError> {
        let v = self.num()?;
        if !v.is_finite() {
            return Err(bad("non-finite value"));
        }
        Ok(v)
    }

    fn param(&mut self) -> Result<f64, FfiError> {
        let v = self.finite()?;
        if v.abs() > MAX_PARAM {
            return Err(bad("parameter too large"));
        }
        Ok(v)
    }

    fn finites(&mut self, k: usize) -> Result<Vec<f64>, FfiError> {
        (0..k).map(|_| self.finite()).collect()
    }

    fn index(&mut self, max: usize) -> Result<usize, FfiError> {
        let v = self.num()?;
        if v.fract() != 0.0 || v < 0.0 || v > max as f64 {
            return Err(bad(format!("bad index/count {v}")));
        }
        Ok(v as usize)
    }

    fn flag(&mut self) -> Result<bool, FfiError> {
        Ok(self.index(1)? == 1)
    }

    fn done(&self) -> Result<(), FfiError> {
        if self.at != self.input.len() {
            return Err(bad("trailing input"));
        }
        Ok(())
    }
}

/// Checks a clamped knot vector of degree `p`.
fn check_knots(knots: &[f64], p: usize, what: &str) -> Result<(), FfiError> {
    let m = knots
        .len()
        .checked_sub(1)
        .ok_or_else(|| bad(format!("{what}: no knots")))?;
    if knots.windows(2).any(|w| w.first() > w.last()) {
        return Err(bad(format!("{what}: knots not non-decreasing")));
    }
    let (Some(&a), Some(&b)) = (knots.first(), knots.last()) else {
        return Err(bad(format!("{what}: empty knots")));
    };
    let short = || bad(format!("{what}: knots too short"));
    if m < 2 * p + 1 {
        return Err(short());
    }
    let head = knots.get(..=p).ok_or_else(short)?;
    let tail = knots.get(m - p..).ok_or_else(short)?;
    let inner = knots.get(p + 1..m - p).ok_or_else(short)?;
    if a >= b || head.iter().any(|&k| k != a) || tail.iter().any(|&k| k != b) {
        return Err(bad(format!("{what}: knots not clamped")));
    }
    if inner.iter().any(|&k| k == a || k == b) {
        return Err(bad(format!("{what}: end multiplicity above degree + 1")));
    }
    let mut run = 1;
    for w in inner.windows(2) {
        run = if w.first() == w.last() { run + 1 } else { 1 };
        if run > p {
            return Err(bad(format!(
                "{what}: interior multiplicity above the degree"
            )));
        }
    }
    Ok(())
}

/// Checks that all points are 3D and share their rational marker.
fn check_points(pts: &[[f64; 4]], what: &str) -> Result<(), FfiError> {
    let [_, _, _, w0] = *pts
        .first()
        .ok_or_else(|| bad(format!("{what}: no points")))?;
    for &[_, _, z, w] in pts {
        if z == NOZ {
            return Err(bad(format!("{what}: 2D data")));
        }
        if (w == NOW) != (w0 == NOW) {
            return Err(bad(format!("{what}: mixed rationality")));
        }
        if w != NOW && w <= 0.0 {
            return Err(bad(format!("{what}: non-positive weight")));
        }
    }
    Ok(())
}

fn points(rd: &mut Reader<'_>, k: usize) -> Result<Vec<[f64; 4]>, FfiError> {
    (0..k)
        .map(|_| Ok([rd.finite()?, rd.finite()?, rd.finite()?, rd.finite()?]))
        .collect()
}

/// A solver's base surface: a B-spline surface, or a plane whose domain is not reversed
/// (its NURB knots are the domain's corners, UB-187).
fn surface(rd: &mut Reader<'_>) -> Result<(), FfiError> {
    if rd.index(1)? == 1 {
        let d = rd.finites(15)?;
        if let [.., u0, v0, u1, v1] = d[..] {
            if u0 > u1 || v0 > v1 {
                return Err(bad("plane: reversed domain (decreasing NURB knots)"));
            }
        }
        return Ok(());
    }
    let n = rd.index(MAX_INDEX)?;
    let m = rd.index(MAX_INDEX)?;
    let p = rd.index(MAX_DEGREE)?;
    let q = rd.index(MAX_DEGREE)?;
    let r = rd.index(MAX_INDEX + MAX_DEGREE + 1)?;
    let s = rd.index(MAX_INDEX + MAX_DEGREE + 1)?;
    if p == 0 || p > n || q == 0 || q > m || r != n + p + 1 || s != m + q + 1 {
        return Err(bad(
            "surface: need 1 <= p <= n, 1 <= q <= m, r = n+p+1, s = m+q+1",
        ));
    }
    let pts = points(rd, (n + 1) * (m + 1))?;
    let u = rd.finites(r + 1)?;
    let v = rd.finites(s + 1)?;
    check_points(&pts, "surface")?;
    check_knots(&u, p, "surface u")?;
    check_knots(&v, q, "surface v")?;
    Ok(())
}

fn kinds(rd: &mut Reader<'_>) -> Result<(), FfiError> {
    let nk = rd.index(MAX_ARRAY)?;
    for _ in 0..nk {
        rd.index(1_000_000)?;
    }
    Ok(())
}

/// The solver block; returns its class.
fn solver_block(rd: &mut Reader<'_>) -> Result<usize, FfiError> {
    let cls = rd.index(2)?;
    let same = rd.flag()?;
    surface(rd)?;
    if !same {
        surface(rd)?;
    }
    rd.param()?;
    rd.param()?;
    rd.flag()?;
    rd.flag()?;
    rd.finite()?;
    Ok(cls)
}

fn validate_cr_init(rd: &mut Reader<'_>) -> Result<(), FfiError> {
    solver_block(rd)?;
    let do2 = rd.flag()?;
    let rail = rd.index(3)?;
    let sel1 = rd.index(4)?;
    if sel1 == 0 {
        return Err(bad("init: NULL first surface"));
    }
    let s1_index = rd.index(MAX_ARRAY)?;
    let dim = rd.index(MAX_ARRAY)?;
    rd.finites(2 * dim)?;
    let np = rd.index(MAX_ARRAY)?;
    for _ in 0..np {
        rd.flag()?;
    }
    let ng = rd.index(MAX_ARRAY)?;
    for _ in 0..ng {
        rd.param()?;
    }
    if do2 {
        if rail > 1 {
            return Err(bad("init: rail index above 1"));
        }
        // reads rGuessT[s1_index..+2], writes rIntervals[ng..+2]
        if s1_index + 1 >= ng || ng + 1 >= dim {
            return Err(bad("init: guess or interval index out of range"));
        }
    }
    Ok(())
}

fn validate_plane(rd: &mut Reader<'_>) -> Result<(), FfiError> {
    if solver_block(rd)? != 2 {
        return Err(bad("plane solve: needs the assisted solver (cls 2)"));
    }
    rd.finites(6)?;
    rd.finites(8)?;
    for _ in 0..4 {
        rd.param()?;
    }
    Ok(())
}

fn validate(op: Op, input: &[f64]) -> Result<(), FfiError> {
    let mut rd = Reader { input, at: 0 };
    let rd = &mut rd;
    match op {
        Op::CrInit => validate_cr_init(rd)?,
        Op::Guess | Op::Plane => validate_plane(rd)?,
        Op::Ss => {
            surface(rd)?;
            surface(rd)?;
            rd.param()?;
            rd.param()?;
            rd.flag()?;
            rd.flag()?;
            for _ in 0..4 {
                rd.param()?;
            }
            rd.flag()?;
            kinds(rd)?;
        }
        Op::Types => kinds(rd)?,
    }
    rd.done()
}

/// Runs `op` on `input` (layouts in the module docs); returns the shim's output record.
///
/// # Errors
///
/// [`FfiError::InvalidArgument`] for input the shim must not see; [`FfiError::Status`]
/// for a shim failure.
pub fn eval(op: Op, input: &[f64]) -> Result<Vec<f64>, FfiError> {
    validate(op, input)?;
    call_with_buffer(65536, |buf, len| {
        // SAFETY: input is valid for input.len() reads and holds exactly what the shim
        // reads for `op` (validated above: well-formed clamped 3D surfaces, non-reversed
        // plane domains, every index, count and flag the shim passes to upstream, the
        // array sizes upstream indexes); buf is valid for buf.len() writes and the shim
        // honours the capacity; len is valid.
        unsafe {
            sys::sc_blend_fillet_standard_solver_p2_eval(
                op.code(),
                input.as_ptr(),
                input.len() as i64,
                buf.as_mut_ptr(),
                buf.len() as i64,
                &raw mut *len,
            )
        }
    })
}
