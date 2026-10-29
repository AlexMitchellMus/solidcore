//! Calls from the `AssertValid` walk into checks other work packages port.
//!
//! Upstream: `source/SMLib/src/SmEdgeuse.cpp`, `SmBSplineCurve.cpp`, `SmLoopuse.cpp`,
//! `SmVertex.cpp`, `SmShell.cpp`, `SmFace.cpp`, `SmEdge.cpp`, `SmSurface.cpp` (and the
//! surface classes), `SmOffsetCurve.cpp`, `SmProjectedCurve.cpp`, `SmHermiteCurve.cpp`,
//! `SmCrvOnSurf.cpp`, `SmCrvInVolume.cpp`.
//!
//! | Rust (private seam)          | Upstream                                  |
//! |------------------------------|-------------------------------------------|
//! | `has_dog_leg_near_pole`      | `SmEdgeuse::HasDogLegNearPole`            |
//! | `assert_no_periodic_jump`    | `SmBSplineCurve::AssertNoPeriodicJump`    |
//! | `is_loopuse_side_in_surface` | `SmLoopuse::IsLoopuseSideInSurface`       |
//! | `vertex_check_pointers`      | `SmVertex::CheckPointers`                 |
//! | `is_inner_shell_of`          | `SmShell::IsInnerShellOf`                 |
//! | `assert_test_face37`         | `sm_AssertTestFace37`                     |
//! | `surface_assert_valid`       | `SmSurface::AssertValid` (virtual)        |
//! | `offset_curve_assert_valid`  | `SmOffsetCurve::AssertValid`              |
//! | `projected_curve_assert_valid` | `SmProjectedCurve::AssertValid`         |
//! | `hermite_curve_assert_valid` | `SmHermiteCurve::AssertValid`             |
//! | `crv_on_surf_assert_valid`   | `SmCrvOnSurf::AssertValid`                |
//! | `crv_in_volume_assert_valid` | `SmCrvInVolume::AssertValid`              |
//! | `assert_test_face32`         | `sm_AssertTestFace32`                     |
//! | `assert_test_face33`         | `sm_AssertTestFace33`                     |
//! | `assert_test_face38`         | `sm_AssertTestFace38`                     |
//! | `assert_test_edge13`         | `sm_AssertTestEdge13`                     |
//! | `assert_test_vertex1`        | `sm_AssertTestVertex1`                    |
//!
//! Deviations: one function per upstream callee, so each package replaces exactly one
//! body. Until a callee is ported its seam returns [`SolidError::Unported`]; a
//! difftest case reaching it is deferred, never passed. Upstream evaluates the level-2
//! seams (Face 32/33/38, Edge 13) only at `SM_LEVEL_2`; `sm_AssertTestVertex1` runs at
//! every level (only its report is level 2).
//!
//! Signatures are the contract with the packages porting the callees: keep them, and
//! replace a body with the call into the ported function.

use crate::curve::curve::Curve;
use crate::error::{Result, SolidError};
use crate::geom::assert_array::{AssertArray, ObjectKey};
use crate::geom::extent2d::Extent2d;
use crate::kernel::Kernel;
use crate::topology::{Brep, EdgeId, EdgeuseId, FaceId, LoopuseId, ShellId, VertexId};

fn unported(upstream: &'static str, unit: &'static str) -> SolidError {
    SolidError::Unported { upstream, unit }
}

// ---- level-0 helpers (WP2) ------------------------------------------------------------

/// `SmEdgeuse::HasDogLegNearPole()` of `eu`: `true` when the UV trim's control polygon
/// turns to follow a singular surface boundary (Edgeuse test 20 fails on `true`).
pub(super) fn has_dog_leg_near_pole(
    _k: &mut Kernel,
    _b: &mut Brep,
    _eu: EdgeuseId,
) -> Result<bool> {
    // UNPORTED(topology.edgeuse_p2): SmEdgeuse::HasDogLegNearPole
    Err(unported(
        "SmEdgeuse::HasDogLegNearPole",
        "topology.edgeuse_p2",
    ))
}

/// `m_pUVTrimCurve->AssertNoPeriodicJump(bClosedU, bClosedV, sUVDomain)`: `true` when
/// the UV trim `curve` of an edge use does not jump across a periodic boundary of the
/// surface domain `domain` (closed in U / V as given). Edgeuse test 7.
pub(super) fn assert_no_periodic_jump(
    _k: &Kernel,
    _curve: &Curve,
    _closed_u: bool,
    _closed_v: bool,
    _domain: &Extent2d,
) -> Result<bool> {
    // UNPORTED(nurbs.b_spline_curve_p3): SmBSplineCurve::AssertNoPeriodicJump
    Err(unported(
        "SmBSplineCurve::AssertNoPeriodicJump",
        "nurbs.b_spline_curve_p3",
    ))
}

/// `SmLoopuse::IsLoopuseSideInSurface()` of the loop's primary loop use `lu`. Loop
/// test 2.
pub(super) fn is_loopuse_side_in_surface(
    _k: &mut Kernel,
    _b: &mut Brep,
    _lu: LoopuseId,
) -> Result<bool> {
    // UNPORTED(topology.loopuse): SmLoopuse::IsLoopuseSideInSurface
    Err(unported(
        "SmLoopuse::IsLoopuseSideInSurface",
        "topology.loopuse",
    ))
}

/// `SmVertex::CheckPointers()` of `v`. Vertex test 2.
pub(super) fn vertex_check_pointers(_k: &mut Kernel, _b: &mut Brep, _v: VertexId) -> Result<bool> {
    // UNPORTED(topology.vertex): SmVertex::CheckPointers
    Err(unported("SmVertex::CheckPointers", "topology.vertex"))
}

/// `pInnerShell->IsInnerShellOf(pOuterShell)` (`== TRUE`). Brep test 12, evaluated only
/// for the inner shells of a region other than the infinite one.
pub(super) fn is_inner_shell_of(
    _k: &mut Kernel,
    _b: &mut Brep,
    _inner: ShellId,
    _outer: ShellId,
) -> Result<bool> {
    // UNPORTED(topology.shell): SmShell::IsInnerShellOf
    Err(unported("SmShell::IsInnerShellOf", "topology.shell"))
}

/// `sm_AssertTestFace37(this)` (`== TRUE`): the face's inner loops nest inside the
/// outer loop and outside each other (with `sm_GetLoopSampleUVPoint`). Face test 37.
pub(super) fn assert_test_face37(_k: &mut Kernel, _b: &mut Brep, _f: FaceId) -> Result<bool> {
    // UNPORTED(topology.face_p4): sm_AssertTestFace37
    Err(unported("sm_AssertTestFace37", "topology.face_p4"))
}

// ---- geometry checks (WP3a / WP3b) ---------------------------------------------------

/// `m_pSurface->AssertValid(pAList, eTestLevel, SM_NO_WALK)` for the surface of face
/// `f` (virtual: the surface's class check, with its base classes), reporting into
/// `list` about the surface object `key`. `level` is `SM_LEVEL_0`..`SM_LEVEL_2` (the
/// call is skipped at `SM_LEVEL_GIVEN`). Called only when the face has a surface.
pub(super) fn surface_assert_valid(
    _k: &mut Kernel,
    _b: &mut Brep,
    _f: FaceId,
    _list: Option<&mut AssertArray>,
    _level: u8,
    _key: ObjectKey,
) -> Result<bool> {
    // UNPORTED(surface.surface_assert_valid): SmSurface::AssertValid
    Err(unported(
        "SmSurface::AssertValid",
        "surface.surface_assert_valid",
    ))
}

/// `SmOffsetCurve::AssertValid` of `curve` (an offset curve), reporting about `key`.
pub(super) fn offset_curve_assert_valid(
    _k: &Kernel,
    _curve: &mut Curve,
    _list: Option<&mut AssertArray>,
    _level: u8,
    _key: ObjectKey,
) -> Result<bool> {
    // UNPORTED(curve.offset_curve): SmOffsetCurve::AssertValid
    Err(unported("SmOffsetCurve::AssertValid", "curve.offset_curve"))
}

/// `SmProjectedCurve::AssertValid` of `curve` (a projected curve), reporting about `key`.
pub(super) fn projected_curve_assert_valid(
    _k: &Kernel,
    _curve: &mut Curve,
    _list: Option<&mut AssertArray>,
    _level: u8,
    _key: ObjectKey,
) -> Result<bool> {
    // UNPORTED(curve.projected_curve): SmProjectedCurve::AssertValid
    Err(unported(
        "SmProjectedCurve::AssertValid",
        "curve.projected_curve",
    ))
}

/// `SmHermiteCurve::AssertValid` of `curve` (a Hermite curve), reporting about `key`.
pub(super) fn hermite_curve_assert_valid(
    _k: &Kernel,
    _curve: &mut Curve,
    _list: Option<&mut AssertArray>,
    _level: u8,
    _key: ObjectKey,
) -> Result<bool> {
    // UNPORTED(curve.hermite_curve): SmHermiteCurve::AssertValid
    Err(unported(
        "SmHermiteCurve::AssertValid",
        "curve.hermite_curve",
    ))
}

/// `SmCrvOnSurf::AssertValid` of `curve` (a curve on surface or a tangent field, which
/// inherits it), reporting about `key`.
pub(super) fn crv_on_surf_assert_valid(
    _k: &Kernel,
    _curve: &mut Curve,
    _list: Option<&mut AssertArray>,
    _level: u8,
    _key: ObjectKey,
) -> Result<bool> {
    // UNPORTED(curve.crv_on_surf): SmCrvOnSurf::AssertValid
    Err(unported("SmCrvOnSurf::AssertValid", "curve.crv_on_surf"))
}

/// `SmCrvInVolume::AssertValid` of `curve` (a curve in volume), reporting about `key`.
pub(super) fn crv_in_volume_assert_valid(
    _k: &Kernel,
    _curve: &mut Curve,
    _list: Option<&mut AssertArray>,
    _level: u8,
    _key: ObjectKey,
) -> Result<bool> {
    // UNPORTED(curve.crv_in_volume): SmCrvInVolume::AssertValid
    Err(unported(
        "SmCrvInVolume::AssertValid",
        "curve.crv_in_volume",
    ))
}

// ---- level-2 tests (WP4) --------------------------------------------------------------

/// `sm_AssertTestFace32(this)` (`== TRUE`): no coincident vertices
/// (`HasCoincidentVertices`). Face test 32, evaluated only at `SM_LEVEL_2` and above.
pub(super) fn assert_test_face32(_k: &mut Kernel, b: &mut Brep, f: FaceId) -> Result<bool> {
    crate::topology::face::assert_tests::assert_test_face32(b, f)
}

/// `sm_AssertTestFace33(this)` (`== TRUE`): no coincident edges (`HasCoincidentEdges`).
/// Face test 33, evaluated only at `SM_LEVEL_2` and above.
pub(super) fn assert_test_face33(k: &mut Kernel, b: &mut Brep, f: FaceId) -> Result<bool> {
    crate::topology::face::assert_tests::assert_test_face33(k, b, f)
}

/// `sm_AssertTestFace38(this)` (`== TRUE`): the face loops are disjoint in UV
/// (`sm_DoLoopsIntersectInUV`). Face test 38, evaluated only at `SM_LEVEL_2` and above
/// and only for a face with a surface.
pub(super) fn assert_test_face38(k: &mut Kernel, b: &mut Brep, f: FaceId) -> Result<bool> {
    crate::topology::face::assert_tests::assert_test_face38(k, b, f)
}

/// `sm_AssertTestEdge13(this, &dMaxGap)`: `(passed, max_gap)`, the edge curve within
/// `SmTol::GetXSectTol3d(this, this)` of every face surface (`IsOnFaceEdgeForm`).
/// Edge test 13, evaluated only at `SM_LEVEL_2` and above; `max_gap` is the report
/// value (upstream initialises it to `0.0`).
pub(super) fn assert_test_edge13(
    _k: &mut Kernel,
    _b: &mut Brep,
    _e: EdgeId,
) -> Result<(bool, f64)> {
    // UNPORTED(topology.edge): sm_AssertTestEdge13
    Err(unported("sm_AssertTestEdge13", "topology.edge"))
}

/// `sm_AssertTestVertex1(this, &dMaxGap3d, &pMaxGapFace)`: `(passed, max_gap_3d)`, the
/// vertex within the cross-section tolerance of its faces and on their corner
/// (`IsWithinXSectTol3dOfFaces`, `IsOnFaceCorner`). Called at every level; its report
/// (Vertex test 1) is level 2.
pub(super) fn assert_test_vertex1(
    _k: &mut Kernel,
    _b: &mut Brep,
    _v: VertexId,
) -> Result<(bool, f64)> {
    // UNPORTED(topology.vertex): sm_AssertTestVertex1
    Err(unported("sm_AssertTestVertex1", "topology.vertex"))
}
