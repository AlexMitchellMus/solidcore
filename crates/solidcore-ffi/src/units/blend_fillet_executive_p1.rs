//! FFI wrappers for unit `blend.fillet_executive_p1` (the fillet laws). Shim:
//! `harness/cpp/shim/blend_fillet_executive_p1.cpp`. Upstream: `SmFilletExecutive.cpp`
//! (`SmFilletLaw`, `SmLinearFilletLaw`, `SmBSplineFilletLaw`: constructors, `Evaluate`,
//! `IsLinearLaw`, `GetType`, `IsKindOf`, `GetExtendedCurve`).
//!
//! Input layout (flat doubles, validated here):
//!
//! ```text
//! kind              0 linear `r0 r1`; 1 B-spline `n p m cpts(x y z w)... knots...`; 2 base
//! has_map e0 e1     edge map (e0 <= e1), used when has_map = 1 and kind != 2
//! ncalls            then per call `t i0 i1 reverse` (i0 <= i1)
//! nk kinds[nk]      type codes for IsKindOf
//! ```
use crate::units::blend_fillet_standard_solver_p1::{bad, bcurve, Reader, MAX_ARRAY, MAX_CALLS};
use crate::{call_with_buffer, FfiError};

mod sys {
    extern "C" {
        pub fn sc_blend_fillet_executive_p1_law_eval(
            input: *const f64,
            n_in: i64,
            out: *mut f64,
            cap: i64,
            out_len: *mut i64,
        ) -> i32;
    }
}

/// Reads an interval `a b` with `|a|, |b| <= MAX_PARAM` (`Reader::param`) and `a <= b`.
fn interval(rd: &mut Reader<'_>, what: &str) -> Result<(), FfiError> {
    let a = rd.param()?;
    let b = rd.param()?;
    if a > b {
        return Err(bad(format!("{what}: reversed interval")));
    }
    Ok(())
}

fn validate(input: &[f64]) -> Result<(), FfiError> {
    let mut rd = Reader { input, at: 0 };
    let rd = &mut rd;
    match rd.index(2)? {
        0 => {
            rd.finite()?;
            rd.finite()?;
        }
        1 => {
            bcurve(rd)?;
        }
        _ => {}
    }
    rd.flag()?;
    interval(rd, "edge map")?;
    let calls = rd.index(MAX_CALLS)?;
    for _ in 0..calls {
        rd.param()?;
        interval(rd, "edge interval")?;
        rd.flag()?;
    }
    let nk = rd.index(MAX_ARRAY)?;
    for _ in 0..nk {
        rd.index(i32::MAX as usize)?;
    }
    rd.done()
}

/// Builds the law described by `input` (layout in the module docs), evaluates it at every
/// call and returns the shim's output record.
///
/// # Errors
///
/// [`FfiError::InvalidArgument`] for input the shim must not see (a malformed curve,
/// non-finite or oversized values, a reversed interval, too many calls);
/// [`FfiError::Status`] for a shim failure.
pub fn law_eval(input: &[f64]) -> Result<Vec<f64>, FfiError> {
    validate(input)?;
    call_with_buffer(4096, |buf, len| {
        // SAFETY: input is valid for input.len() reads and holds exactly what the shim
        // reads (validated above: a well-formed clamped 3D curve, every count, flag and
        // index the shim passes to upstream); buf is valid for buf.len() writes and the
        // shim honours the capacity; len is valid.
        unsafe {
            sys::sc_blend_fillet_executive_p1_law_eval(
                input.as_ptr(),
                input.len() as i64,
                buf.as_mut_ptr(),
                buf.len() as i64,
                &raw mut *len,
            )
        }
    })
}
