//! FFI wrappers for unit `query.brep`.
//!
//! Shim: `harness/cpp/shim/query_brep.cpp` (`sc_query_brep_evaluate`,
//! `sc_query_brep_attribute`).
//! Upstream: `source/SM_API/src/SmApiSurfaces.cpp`, `source/SMLib/inc/SmPolyDecimate.h`,
//! `source/SMLib/src/SmPolyDecimate.cpp`.
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
//! before the copy: `sort class distance error N` (`SmDecimateVertexAttr`), `Q0 Q1 Q2 b
//! c` (`SmQuadricVertexAttr`) or `P error sort_set sort` (`SmQuadricEdgeAttr`; the sort
//! index is written only when `sort_set`). It returns `id behavior type data` of the
//! constructed attribute, then of its `MakeCopy`; the data as written above (the
//! quadric edge's sort index NaN while unset). Ids and sort indices are integers in
//! `0..=2^53`, classes `0..=5`.

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
            Self::DecimateVertex => 10,
            Self::QuadricVertex => 16,
            Self::QuadricEdge => 8,
        }
    }
}

/// Largest id or sort index (every integer up to it is exact in an `f64`).
const MAX_ULONG: f64 = 9_007_199_254_740_992.0;

mod sys {
    extern "C" {
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

/// An integer in `0..=2^53` (an `ULONG` the shim converts exactly).
fn ulong(rd: &mut Reader<'_>) -> Result<(), FfiError> {
    let v = rd.num()?;
    if v.fract() != 0.0 || !(0.0..=MAX_ULONG).contains(&v) {
        return Err(bad(format!("bad ULONG {v}")));
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
                ulong(&mut rd)?;
                rd.index(5)?;
                for _ in 0..5 {
                    rd.num()?;
                }
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
