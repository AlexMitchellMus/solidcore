//! FFI wrappers for unit `blend.fillet_corner_p1`.
//!
//! Shim: `harness/cpp/shim/blend_fillet_corner_p1.cpp` (`sc_blend_fillet_corner_p1_helper`),
//! which compiles `source/SMLib/src/SmFilletCorner.cpp` to reach its file-static helpers.
//!
//! | Op                     | Upstream                                  | Input after the op                 | Output |
//! |------------------------|-------------------------------------------|------------------------------------|--------|
//! | [`Op::Planar`]         | `sm_CreatePlanarSurface`                  | `n p(3).. tol` (`n <= 6`)          | `status` [surface block] |
//! | [`Op::Revolution`]     | `sm_CreateSrfOfRevolution`                | `n curve.. major tol` (`1..=8` 3D curves) | `status` [surface block] |
//! | [`Op::FilletRadius`]   | `sm_GetFilletRadius`                      | `p1 n1 p2 n2 orient`               | `radius span` |
//! | [`Op::SeamOut`]        | `sm_CheckSeamOut`                         | `surface u v tol closure`          | `ret u v` |
//! | [`Op::Closure`]        | `sm_SetSurfaceClosure`                    | `surface u0 v0 u1 v1`              | `closure` |
//! | [`Op::G1`]             | `sm_TestG1Continuity`                     | block, `v e1 e2 tol`               | `0 flag has_angle angle` |
//! | [`Op::FaceByTwo`]      | `sm_GetFaceByTwoEdgeuses`                 | block, `eu1 eu2 eu3`               | `0 face` |
//! | [`Op::Tangential`]     | `sm_CheckTangentialEdgeuse`               | block, `eu tol`                    | `0 flag` |
//! | [`Op::Concavity`]      | `sm_CheckFaceComplimentaryCornerConcavity`| block, `v eu1 eu2 ref`             | `0 flag` |
//!
//! Curves and surfaces are unit `topology.edgeuse_p1`'s curve and surface blocks
//! (`dim n p m cpts(x y z w) knots`; `1 O N` for an `SmPlane`, else the B-spline surface
//! block); a surface block in the output is the shim's `put_surface` (`1 O X Y scale
//! step_domain` for an `SmPlane`, `0 oob form n m p q r s cpts u v` for a B-spline
//! surface). `status` is `0` for `SM_SUCCESS`, else the `SmStatus`; `closure` is the
//! `SmSurfParamType` code.
//!
//! The topology ops ([`Op::is_brep`]) start with `topology.edgeuse_p1`'s B-rep block
//! (`path` is the `.smb` file of a file block) and address the built B-rep's vertices,
//! edges and edge uses by their index in that unit's export (`-1` = NULL where the
//! routine accepts it). Their output starts with `0`, or is `-1 status` when the B-rep
//! could not be built; a result of `-2` means "refused": the routine would dereference
//! null (an edge without a curve, an edge use without a radial, vertex use or edge).
//!
//! Inputs are validated before the call: the blocks, every index range, finite and
//! bounded coordinates, a non-empty curve array with the major arc inside it (upstream
//! indexes it unchecked) and closure codes `0..=4`.

use crate::{call_with_buffer, FfiError};

use super::topology_edgeuse_p1::{self as eu1, bad, brep_block, curve_up_to, surface, Reader};

/// The largest point count of [`Op::Planar`] (upstream handles 3 and 4; others fail).
pub const MAX_PLANAR_POINTS: usize = 6;
/// The largest curve count of [`Op::Revolution`].
pub const MAX_CURVES: usize = 8;
/// The largest control point index of an [`Op::Revolution`] curve.
pub const MAX_CURVE_INDEX: usize = 16;
/// The largest `SmSurfParamType` code (`SM_SP_NEITHER`).
pub const MAX_CLOSURE: usize = 4;

/// The helper run by [`eval`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Op {
    /// `sm_CreatePlanarSurface`.
    Planar = 0,
    /// `sm_CreateSrfOfRevolution`.
    Revolution = 1,
    /// `sm_GetFilletRadius`.
    FilletRadius = 2,
    /// `sm_CheckSeamOut`.
    SeamOut = 3,
    /// `sm_SetSurfaceClosure`.
    Closure = 4,
    /// `sm_TestG1Continuity`.
    G1 = 5,
    /// `sm_GetFaceByTwoEdgeuses`.
    FaceByTwo = 6,
    /// `sm_CheckTangentialEdgeuse`.
    Tangential = 7,
    /// `sm_CheckFaceComplimentaryCornerConcavity`.
    Concavity = 8,
}

impl Op {
    /// Every op, in code order.
    pub const ALL: [Self; 9] = [
        Self::Planar,
        Self::Revolution,
        Self::FilletRadius,
        Self::SeamOut,
        Self::Closure,
        Self::G1,
        Self::FaceByTwo,
        Self::Tangential,
        Self::Concavity,
    ];

    /// The shim's op code.
    #[must_use]
    pub fn code(self) -> i32 {
        self as i32
    }

    /// The op with shim code `code`, `None` outside `0..=8`.
    #[must_use]
    pub fn from_code(code: i64) -> Option<Self> {
        Self::ALL
            .into_iter()
            .find(|op| i64::from(op.code()) == code)
    }

    /// Whether the input starts with a B-rep block.
    #[must_use]
    pub fn is_brep(self) -> bool {
        matches!(
            self,
            Self::G1 | Self::FaceByTwo | Self::Tangential | Self::Concavity
        )
    }
}

mod sys {
    extern "C" {
        pub fn sc_blend_fillet_corner_p1_helper(
            op: i32,
            path: *const std::ffi::c_char,
            input: *const f64,
            n_in: i64,
            out: *mut f64,
            cap: i64,
            out_len: *mut i64,
        ) -> i32;
    }
}

fn raw(op: Op, path: &str, input: &[f64]) -> Result<Vec<f64>, FfiError> {
    let cpath = std::ffi::CString::new(path).map_err(|_| bad("path has a NUL byte"))?;
    call_with_buffer(1024, |buf, len| {
        // SAFETY: cpath is a valid NUL-terminated string; input is valid for
        // input.len() reads and holds what the shim reads for `op` (validated by
        // `eval`); buf is valid for buf.len() writes and the shim honours the capacity;
        // len is valid.
        unsafe {
            sys::sc_blend_fillet_corner_p1_helper(
                op.code(),
                cpath.as_ptr(),
                input.as_ptr(),
                input.len() as i64,
                buf.as_mut_ptr(),
                buf.len() as i64,
                &raw mut *len,
            )
        }
    })
}

/// The geometry ops' arguments.
fn geom_args(op: Op, rd: &mut Reader<'_>) -> Result<(), FfiError> {
    match op {
        Op::Planar => {
            let n = rd.index(MAX_PLANAR_POINTS)?;
            for _ in 0..n {
                rd.vec3()?;
            }
            rd.finite()?;
        }
        Op::Revolution => {
            let n = rd.index(MAX_CURVES)?;
            if n == 0 {
                return Err(bad(
                    "revolution: no curves (upstream indexes the empty array)",
                ));
            }
            for _ in 0..n {
                if rd.input.get(rd.at).copied() != Some(3.0) {
                    return Err(bad("revolution: curves must be 3D"));
                }
                curve_up_to(rd, MAX_CURVE_INDEX)?;
            }
            rd.below(n)?;
            rd.finite()?;
        }
        Op::FilletRadius => {
            for _ in 0..4 {
                rd.vec3()?;
            }
            rd.finite()?;
        }
        Op::SeamOut => {
            surface(rd)?;
            rd.finite()?;
            rd.finite()?;
            rd.finite()?;
            rd.index(MAX_CLOSURE)?;
        }
        Op::Closure => {
            surface(rd)?;
            rd.finites(4)?;
        }
        Op::G1 | Op::FaceByTwo | Op::Tangential | Op::Concavity => {
            return Err(bad("not a geometry op"));
        }
    }
    Ok(())
}

/// The topology ops' arguments; `c` holds the export's counts (`[.., edge uses (6),
/// edges (7), vertex uses (8), vertices (9)]`).
fn brep_args(op: Op, rd: &mut Reader<'_>, c: [usize; 10]) -> Result<(), FfiError> {
    let (n_eu, n_e, n_v) = (c[6], c[7], c[9]);
    match op {
        Op::G1 => {
            rd.below_or_none(n_v)?;
            rd.below_or_none(n_e)?;
            rd.below_or_none(n_e)?;
            rd.finite()?;
        }
        Op::FaceByTwo => {
            rd.below_or_none(n_eu)?;
            rd.below_or_none(n_eu)?;
            rd.below_or_none(n_eu)?;
        }
        Op::Tangential => {
            rd.below(n_eu)?;
            rd.finite()?;
        }
        Op::Concavity => {
            rd.below(n_v)?;
            rd.below(n_eu)?;
            rd.below(n_eu)?;
            rd.below_or_none(n_eu)?;
        }
        Op::Planar | Op::Revolution | Op::FilletRadius | Op::SeamOut | Op::Closure => {
            return Err(bad("not a topology op"));
        }
    }
    Ok(())
}

/// Validates `input` for `op` without running it.
///
/// # Errors
///
/// [`FfiError::InvalidArgument`] for a malformed input; for the topology ops, also the
/// errors of building and exporting the B-rep (which runs the shim of unit
/// `topology.edgeuse_p1`).
pub fn validate(op: Op, path: &str, input: &[f64]) -> Result<(), FfiError> {
    let mut rd = Reader { input, at: 0 };
    if !op.is_brep() {
        geom_args(op, &mut rd)?;
        return rd.done();
    }
    brep_block(&mut rd, path)?;
    let block = input.get(..rd.at).ok_or_else(|| bad("input too short"))?;
    let export = eu1::eval(eu1::Op::Export, path, block)?;
    let Some(c) = eu1::counts(&export) else {
        // the B-rep cannot be built: the shim reports it; the arguments are not checked
        return Ok(());
    };
    brep_args(op, &mut rd, c)?;
    rd.done()
}

/// Runs `op` on `input` (layouts in the module docs); `path` is the `.smb` file of a
/// file B-rep block (ignored otherwise).
///
/// # Errors
///
/// [`FfiError::InvalidArgument`] for a malformed input; [`FfiError::Status`] for a shim
/// failure.
pub fn eval(op: Op, path: &str, input: &[f64]) -> Result<Vec<f64>, FfiError> {
    validate(op, path, input)?;
    raw(op, path, input)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn op_codes_round_trip() {
        for op in Op::ALL {
            assert_eq!(Op::from_code(i64::from(op.code())), Some(op));
        }
        assert_eq!(Op::from_code(9), None);
        assert_eq!(Op::from_code(-1), None);
    }

    #[test]
    fn rejects_bad_geometry_inputs() {
        // too many planar points
        assert!(validate(Op::Planar, "", &[7.0]).is_err());
        // no curves
        assert!(validate(Op::Revolution, "", &[0.0, 0.0, 1e-6]).is_err());
        // trailing input
        assert!(validate(Op::FilletRadius, "", &[0.0; 14]).is_err());
        assert!(validate(Op::FilletRadius, "", &[0.0; 13]).is_ok());
        // closure code out of range
        let plane = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0];
        let mut v = plane.to_vec();
        v.extend([0.0, 0.0, 1e-6, 5.0]);
        assert!(validate(Op::SeamOut, "", &v).is_err());
    }
}
