//! FFI wrappers for unit `query.brep`.
//!
//! Shim: `harness/cpp/shim/query_brep.cpp` (`sc_query_brep_evaluate`,
//! `sc_query_brep_attribute`, `sc_query_brep_box`, `sc_query_brep_classify`).
//! Upstream: `source/SM_API/src/SmApiSurfaces.cpp`, `source/SMLib/inc/SmPolyDecimate.h`,
//! `source/SMLib/src/SmPolyDecimate.cpp`, `source/SMLib/src/SmBrep.cpp`.
//!
//! Notation as in [`super::surface_surface_cache_p3`]: `S` a surface `kind ..` (`0 B`:
//! B-spline; `1 o x y sx sy D`: bounded `SmPlane`; `2 B`: `SmCurveBoundedSurface` over the
//! B-spline), `V` a vector `x y z`.
//!
//! | [`EvalOp`]    | upstream                          |
//! |---------------|-----------------------------------|
//! | `Point`       | `SmApiEvaluateSurfacePoint`       |
//! | `Normal`      | `SmApiEvaluateSurfaceNormal`      |
//! | `Derivatives` | `SmApiEvaluateSurfaceDerivatives` |
//!
//! [`evaluate`] takes `has_surface [S] u v P Du Dv` (`has_surface = 0`: a `NULL` surface;
//! `P Du Dv` the values the outputs hold before the call) and returns the status (`0` =
//! `SM_SUCCESS`, else the `SmStatus`) followed by `P Du Dv` as the entry point leaves them.
//! The surface must be one [`super::surface_surface_cache_p3`] accepts; the parameter and
//! the output values may be anything (the entry points reject a parameter outside the
//! natural domain, NaN included, before evaluating).
//!
//! [`attribute`] takes `kind id has_behavior behavior set [data]`: `kind` an
//! [`AttrKind`], `behavior` an `SmAttributeBehaviorType` (`0..=4`, used when
//! `has_behavior`, else the constructor's default), and with `set` the data written
//! before the copy: `class distance error N sort_set sort` (`SmDecimateVertexAttr`),
//! `Q0 Q1 Q2 b c` (`SmQuadricVertexAttr`) or `P error sort_set sort`
//! (`SmQuadricEdgeAttr`); a sort index is written only when `sort_set`. An id or sort
//! index (`ULONG`) is its two 32-bit halves `hi lo`, so all 64 bits are exact. It
//! returns `id behavior type data` of the constructed attribute, then of its
//! `MakeCopy`, the data as `sort class distance error N`, `Q0 Q1 Q2 b c` or `P error
//! sort` (the quadric edge's sort index NaN NaN while unset). Classes are `0..=5`.
//!
//! [`box_steps`] builds a real `SmBrep` from unit `topology.edgeuse_p1`'s B-rep block
//! (`path` names the `.smb` file of a file block; validated and built once through that
//! unit's wrapper first) and runs the [`BoxStep`]s on it in order, each `kind a b c
//! defaults`: `Box` `CalculateBoundingBox(box, a, b, c)` (`bApprox bTight
//! bExpandByZoneTol3d`; with `defaults` the default arguments), `Tight`
//! `CalculateTightBoundingBox`, `Lock` `LockBoundingBox(a)` (with `defaults` the default
//! argument), `Unlock` `UnLockBoundingBox`. The output is, per step: for `Box` and
//! `Tight` the status and the box `min max` (seeded with `(-17 -19 -23) (29 31 37)`, so
//! what a failing call leaves is compared), for `Lock` the status, for `Unlock` `0`.
//! At most [`MAX_BOX_STEPS`] steps; a B-rep with a face without surface is rejected by
//! the shim (upstream dereferences it).
//!
//! [`classify`] builds the B-rep the same way and runs `Point3DClassify(P, tol, do_bnd,
//! pc)` with `pc` constructed with `tol`; it returns the status and `class kind index gap
//! src_tol presnap tparam u v param_set face_kind face_index` of `pc` (objects as `kind
//! index`: `1` face, `2` edge, `3` vertex, `4` region by traversal index, `0 -1` none).
//! The point and the tolerance must be finite (any sign or magnitude: upstream takes
//! them as given; a non-finite one could not end the tree descents).

use super::surface_surface_cache_p1::{bad, Reader};
use super::surface_surface_cache_p3::sized_surface;
use crate::{call_with_buffer, FfiError};

/// The entry point run by [`evaluate`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EvalOp {
    /// `SmApiEvaluateSurfacePoint`.
    Point,
    /// `SmApiEvaluateSurfaceNormal`.
    Normal,
    /// `SmApiEvaluateSurfaceDerivatives`.
    Derivatives,
}

impl EvalOp {
    fn code(self) -> i32 {
        match self {
            Self::Point => 0,
            Self::Normal => 1,
            Self::Derivatives => 2,
        }
    }
}

/// The decimation attribute class run by [`attribute`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AttrKind {
    /// `SmDecimateVertexAttr`.
    DecimateVertex,
    /// `SmQuadricVertexAttr`.
    QuadricVertex,
    /// `SmQuadricEdgeAttr`.
    QuadricEdge,
}

impl AttrKind {
    fn code(self) -> f64 {
        match self {
            Self::DecimateVertex => 0.0,
            Self::QuadricVertex => 1.0,
            Self::QuadricEdge => 2.0,
        }
    }

    /// The length of one `id behavior type data` output block.
    const fn block_len(self) -> usize {
        match self {
            Self::DecimateVertex => 12,
            Self::QuadricVertex => 17,
            Self::QuadricEdge => 10,
        }
    }
}

/// A step of [`box_steps`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoxStep {
    /// `CalculateBoundingBox(box, approx, tight, expand)`, or with the default arguments.
    Box {
        /// `bApprox`.
        approx: bool,
        /// `bTight`.
        tight: bool,
        /// `bExpandByZoneTol3d`.
        expand: bool,
        /// Call with the default arguments instead.
        defaults: bool,
    },
    /// `CalculateTightBoundingBox`.
    Tight,
    /// `LockBoundingBox(tight)`, or with the default argument.
    Lock {
        /// `bTight`.
        tight: bool,
        /// Call with the default argument instead.
        defaults: bool,
    },
    /// `UnLockBoundingBox`.
    Unlock,
}

impl BoxStep {
    fn record(self) -> [f64; 5] {
        let f = |b: bool| f64::from(u8::from(b));
        match self {
            Self::Box {
                approx,
                tight,
                expand,
                defaults,
            } => [0.0, f(approx), f(tight), f(expand), f(defaults)],
            Self::Tight => [1.0, 0.0, 0.0, 0.0, 0.0],
            Self::Lock { tight, defaults } => [2.0, f(tight), 0.0, 0.0, f(defaults)],
            Self::Unlock => [3.0, 0.0, 0.0, 0.0, 0.0],
        }
    }

    /// The length of the step's output.
    const fn output_len(self) -> usize {
        match self {
            Self::Box { .. } | Self::Tight => 7,
            Self::Lock { .. } | Self::Unlock => 1,
        }
    }
}

/// Largest number of [`box_steps`] steps.
pub const MAX_BOX_STEPS: usize = 8;


mod sys {
    extern "C" {
        pub fn sc_query_brep_classify(
            path: *const std::ffi::c_char,
            block: *const f64,
            n_block: i64,
            args: *const f64,
            n_args: i64,
            out: *mut f64,
            cap: i64,
            out_len: *mut i64,
        ) -> i32;
        pub fn sc_query_brep_box(
            path: *const std::ffi::c_char,
            block: *const f64,
            n_block: i64,
            args: *const f64,
            n_args: i64,
            out: *mut f64,
            cap: i64,
            out_len: *mut i64,
        ) -> i32;
        pub fn sc_query_brep_attribute(
            input: *const f64,
            n_in: i64,
            out: *mut f64,
            cap: i64,
            out_len: *mut i64,
        ) -> i32;
        pub fn sc_query_brep_evaluate(
            op: i32,
            input: *const f64,
            n_in: i64,
            out: *mut f64,
            cap: i64,
            out_len: *mut i64,
        ) -> i32;
    }
}

/// Checks an [`evaluate`] input record.
///
/// # Errors
///
/// [`FfiError::InvalidArgument`] for input the shim must not see.
pub fn validate_evaluate(input: &[f64]) -> Result<(), FfiError> {
    let mut rd = Reader { input, at: 0 };
    if rd.flag()? {
        sized_surface(&mut rd)?;
    }
    // u v P Du Dv: anything (see the module docs)
    for _ in 0..11 {
        rd.num()?;
    }
    rd.done()
}

/// Runs `op` on `input` (layout in the module docs); returns the status and `P Du Dv`.
///
/// # Errors
///
/// [`FfiError::InvalidArgument`] for input the shim must not see; [`FfiError::Status`]
/// for a shim failure.
pub fn evaluate(op: EvalOp, input: &[f64]) -> Result<Vec<f64>, FfiError> {
    validate_evaluate(input)?;
    let out = call_with_buffer(16, |buf, len| {
        // SAFETY: input is valid for input.len() reads and holds exactly what the shim
        // reads (validated above: a well-formed clamped surface or a plane with
        // independent axes, or none); buf is valid for buf.len() writes and the shim
        // honours the capacity; len is valid.
        unsafe {
            sys::sc_query_brep_evaluate(
                op.code(),
                input.as_ptr(),
                input.len() as i64,
                buf.as_mut_ptr(),
                buf.len() as i64,
                &raw mut *len,
            )
        }
    })?;
    if out.len() != 10 {
        return Err(bad("evaluate: unexpected output length"));
    }
    Ok(out)
}

/// An `ULONG` as its 32-bit halves `hi lo` (integers in `0..2^32`).
fn ulong(rd: &mut Reader<'_>) -> Result<(), FfiError> {
    for _ in 0..2 {
        let v = rd.num()?;
        if v.fract() != 0.0 || !(0.0..4_294_967_296.0).contains(&v) {
            return Err(bad(format!("bad ULONG half {v}")));
        }
    }
    Ok(())
}

/// Checks an [`attribute`] input record (`kind` replaces its first entry).
///
/// # Errors
///
/// [`FfiError::InvalidArgument`] for input the shim must not see.
pub fn validate_attribute(kind: AttrKind, data: &[f64]) -> Result<(), FfiError> {
    let mut rd = Reader { input: data, at: 0 };
    ulong(&mut rd)?; // id
    rd.flag()?; // has_behavior
    rd.index(4)?; // behavior
    if rd.flag()? {
        match kind {
            AttrKind::DecimateVertex => {
                rd.index(5)?;
                for _ in 0..5 {
                    rd.num()?;
                }
                rd.flag()?;
                ulong(&mut rd)?;
            }
            AttrKind::QuadricVertex => {
                for _ in 0..13 {
                    rd.num()?;
                }
            }
            AttrKind::QuadricEdge => {
                for _ in 0..4 {
                    rd.num()?;
                }
                rd.flag()?;
                ulong(&mut rd)?;
            }
        }
    }
    rd.done()
}

/// Constructs and copies a `kind` attribute from `data` (the record after `kind`, layout
/// in the module docs); returns both `id behavior type data` blocks.
///
/// # Errors
///
/// [`FfiError::InvalidArgument`] for input the shim must not see; [`FfiError::Status`]
/// for a shim failure.
pub fn attribute(kind: AttrKind, data: &[f64]) -> Result<Vec<f64>, FfiError> {
    validate_attribute(kind, data)?;
    let mut input = vec![kind.code()];
    input.extend_from_slice(data);
    let out = call_with_buffer(2 * kind.block_len(), |buf, len| {
        // SAFETY: input is valid for input.len() reads and holds exactly what the shim
        // reads (validated above: a known class and behaviour, a defined vertex class,
        // exact integers for the ULONG fields); buf is valid for buf.len() writes and the
        // shim honours the capacity; len is valid.
        unsafe {
            sys::sc_query_brep_attribute(
                input.as_ptr(),
                input.len() as i64,
                buf.as_mut_ptr(),
                buf.len() as i64,
                &raw mut *len,
            )
        }
    })?;
    if out.len() != 2 * kind.block_len() {
        return Err(bad("attribute: unexpected output length"));
    }
    Ok(out)
}

/// Runs `steps` on the B-rep of `block` (see the module docs); returns every step's
/// output in order.
///
/// # Errors
///
/// [`FfiError::InvalidArgument`] for a block the graph wrapper rejects, a B-rep that
/// cannot be built, no or too many steps; the shim's status (also for a face without
/// surface).
pub fn box_steps(path: &str, block: &[f64], steps: &[BoxStep]) -> Result<Vec<f64>, FfiError> {
    use super::topology_edgeuse_p1 as graph;
    if steps.is_empty() || steps.len() > MAX_BOX_STEPS {
        return Err(bad("box steps: count"));
    }
    let export = graph::eval(graph::Op::Export, path, block)?;
    if graph::counts(&export).is_none() {
        return Err(bad("native fixture construction failed"));
    }
    let mut args = vec![steps.len() as f64];
    for s in steps {
        args.extend(s.record());
    }
    let want: usize = steps.iter().map(|s| s.output_len()).sum();
    let path = std::ffi::CString::new(path).map_err(|_| bad("path contains NUL"))?;
    let n_block = i64::try_from(block.len()).map_err(|_| bad("block too large"))?;
    let out = call_with_buffer(want, |buf, len| {
        // SAFETY: the block was validated (and built) by the graph wrapper and the steps
        // are well-formed records; all buffers and the C string stay live for the call,
        // and the shim honours `cap`.
        unsafe {
            sys::sc_query_brep_box(
                path.as_ptr(),
                block.as_ptr(),
                n_block,
                args.as_ptr(),
                args.len() as i64,
                buf.as_mut_ptr(),
                buf.len() as i64,
                &raw mut *len,
            )
        }
    })?;
    if out.len() != want {
        return Err(bad("box steps: unexpected output length"));
    }
    Ok(out)
}

/// The length of a [`classify`] output.
pub const CLASSIFY_LEN: usize = 13;

/// Classifies `point` against the B-rep of `block` (see the module docs).
///
/// # Errors
///
/// [`FfiError::InvalidArgument`] for a block the graph wrapper rejects, a B-rep that
/// cannot be built, a non-finite point or tolerance; the shim's status (also for a face
/// without surface).
pub fn classify(
    path: &str,
    block: &[f64],
    point: [f64; 3],
    tol: f64,
    do_boundary: bool,
) -> Result<Vec<f64>, FfiError> {
    use super::topology_edgeuse_p1 as graph;
    if point.iter().any(|x| !x.is_finite()) {
        return Err(bad("classify: point"));
    }
    if !tol.is_finite() {
        return Err(bad("classify: tolerance"));
    }
    let export = graph::eval(graph::Op::Export, path, block)?;
    if graph::counts(&export).is_none() {
        return Err(bad("native fixture construction failed"));
    }
    let args = [
        point[0],
        point[1],
        point[2],
        tol,
        f64::from(u8::from(do_boundary)),
    ];
    let path = std::ffi::CString::new(path).map_err(|_| bad("path contains NUL"))?;
    let n_block = i64::try_from(block.len()).map_err(|_| bad("block too large"))?;
    let out = call_with_buffer(CLASSIFY_LEN, |buf, len| {
        // SAFETY: the block was validated (and built) by the graph wrapper and the
        // arguments are the five values the shim reads; all buffers and the C string stay
        // live for the call, and the shim honours `cap`.
        unsafe {
            sys::sc_query_brep_classify(
                path.as_ptr(),
                block.as_ptr(),
                n_block,
                args.as_ptr(),
                args.len() as i64,
                buf.as_mut_ptr(),
                buf.len() as i64,
                &raw mut *len,
            )
        }
    })?;
    if out.len() != CLASSIFY_LEN {
        return Err(bad("classify: unexpected output length"));
    }
    Ok(out)
}
