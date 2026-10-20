//! FFI wrappers for unit `query.brep`.
//!
//! Shim: `harness/cpp/shim/query_brep.cpp` (`sc_query_brep_evaluate`).
//! Upstream: `source/SM_API/src/SmApiSurfaces.cpp`.
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

mod sys {
    extern "C" {
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
