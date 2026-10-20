//! Fillet corners, part 1: the file-static helpers of `SmFilletCorner.cpp` (edge
//! tangency at a vertex, the face shared by two edge uses, face tangency across an edge,
//! face-corner concavity, the planar and surface-of-revolution corner patches and the
//! variable-radius corner radius) and the `SmFilletCorner` methods this unit owns, as
//! methods of [`FilletCorner`] (the corner record of
//! [`crate::blend::fillet_corner::p3`]): data access, the predicates and setters, the
//! tangent-solver lookups, the extended-surface registry, the corner edge
//! re-initialisation and squeeze, the concave rail-rail test and the Nx1 degeneracy
//! check.
//!
//! Upstream: `source/SMLib/src/SmFilletCorner.cpp`, `source/SMLib/inc/SmFilletCorner.h`.
//!
//! | Rust                                         | Upstream                                    |
//! |----------------------------------------------|---------------------------------------------|
//! | [`test_g1_continuity`]                       | `sm_TestG1Continuity`                       |
//! | [`face_by_two_edgeuses`]                     | `sm_GetFaceByTwoEdgeuses`                   |
//! | [`check_tangential_edgeuse`]                 | `sm_CheckTangentialEdgeuse`                 |
//! | [`check_face_complimentary_corner_concavity`] | `sm_CheckFaceComplimentaryCornerConcavity` |
//! | [`create_srf_of_revolution`]                 | `sm_CreateSrfOfRevolution`                  |
//! | [`create_planar_surface`]                    | `sm_CreatePlanarSurface`                    |
//! | [`fillet_radius`]                            | `sm_GetFilletRadius`                        |
//! | [`check_seam_out`] (re-export)               | `sm_CheckSeamOut`                           |
//! | [`set_surface_closure`] (re-export)          | `sm_SetSurfaceClosure`                      |
//! | [`FilletCorner::this_approx_tol_3d`]         | `GetThisApproxTol3d` (and the `SmFilletCorner` branch of `SmTol::GetApproxTol3d(const SmObject*)`) |
//! | [`FilletCorner::fillet_edges`], [`FilletCorner::fillet_vertices`] | `GetFilletEdges`, `GetFilletVertices` |
//! | [`FilletCorner::filleted_vertex`], [`FilletCorner::set_back_dist`], [`FilletCorner::status`] | `GetFilletedVertex`, `GetSetBackDist`, `GetStatus` |
//! | [`FilletCorner::type_code`], [`FilletCorner::type_string`], [`FilletCorner::is_kind_of`] | `GetType`, `GetTypeString`, `IsKindOf` (`SM_COMMON`) |
//! | [`FilletCorner::add_fillet_vertex`]          | `AddFilletVertex`                           |
//! | [`FilletCorner::set_bevel`], [`FilletCorner::set_calc_as_const_rad`], [`FilletCorner::set_chamfer`] | `SetBevel`, `SetCalcAsConstRad`, `SmFilletConvexNxNCorner::SetChamfer` |
//! | [`FilletCorner::is_bevel_corner`], [`FilletCorner::is_blending_corner`], [`FilletCorner::is_degenerate_corner`], [`FilletCorner::is_extended_surface_patch`], [`FilletCorner::is_calc_as_const_rad`] | `IsBevelCorner`, `IsBlendingCorner`, `IsDegenerateCorner`, `IsExtendedSurfacePatch`, `IsCalcAsConstRad` |
//! | [`FilletCorner::is_tangent_corner`]          | `IsTangentCorner` (base, `SmFilletNx1ClosedCorner`, `SmFilletNx2Corner`) |
//! | [`FilletCorner::tangent_fillet_solver`]      | `GetTangentFilletSolver` (base, `SmFillet1x1Corner`, `SmFillet2x2Corner`) |
//! | [`FilletCorner::fillet_solvers`]             | `GetFilletSolvers`                          |
//! | [`FilletCorner::extended_surface`], [`FilletCorner::add_extended_original_surface`] | `GetExtendedSurface`, `AddExtendedOriginalSurface` |
//! | [`FilletCorner::re_initialize_corner_edges`] | `ReInitializeCornerEdges`                   |
//! | [`FilletCorner::squeeze_edge`]               | `SqueezeEdge`                               |
//! | [`FilletCorner::test_rail_rail_interpolation_curve`] | `SmFilletConcaveNxNCorner::TestRailRailInterpolationCurve` |
//! | [`check_nx1_degeneracy`]                     | `SmFilletNx1Corner::CheckNx1Degeneracy` (static) |
//! | `p3::FilletCorner::adjust_topology` / `adjust_tangent_corner` (`_` arms) | `SmFilletCorner::AdjustTopology`, `SmFillet4x3Corner::AdjustTopology`, `SmFilletNxNCorner::AdjustTopology`, `SmFilletCorner::AdjustTangentCorner` (all `return SM_SUCCESS`) |
//!
//! Deviations:
//! - **Shared helpers.** `sm_CheckSeamOut` and `sm_SetSurfaceClosure` are textually the
//!   same in `SmFilletCorner.cpp` and `SmFilletGeom.cpp` (each file has its own static
//!   copy): both are the ports in [`crate::blend::fillet_geom::p1`], re-exported here.
//! - **Out-params.** Results are return values: [`test_g1_continuity`] returns the
//!   flag and the optional angle (`pdAngDeg`, only set for non-degenerate tangents);
//!   [`fillet_radius`] returns `(radius, angle span)`. The `iDebugLevel` arguments
//!   (graphics only) are dropped.
//! - **Ignored statuses.** Where upstream ignores a callee's status (`SE`, or no check:
//!   `SmCurve::Evaluate`, `SmEdgeuse::EvaluateBinormal`, `SmVector3d::AngleBetween`,
//!   `Unitize`, `AdjustSTEPUVDomain`), the port continues with the values upstream's
//!   outputs keep: `SmVector3d` members start as `SM_UNDEF_DOUBLE`, `AngleBetween`
//!   reports `0`, `Unitize` stores `(0, 0, 1)`. Only [`SolidError::Unported`] propagates.
//!   `SmVector3d::CCWAngleBetween` (ignored in `sm_GetFilletRadius`) cannot fail there:
//!   its projected vectors pass a `100 * SM_EFF_ZERO_SQ` length test before
//!   `AngleBetween`, whose zero test is the smaller `SM_EFF_ZERO_SQ`.
//! - **Undefined inputs.** Where upstream dereferences null (an edge without a curve or
//!   edge uses, an edge use without a radial, vertex use or edge) or indexes past the
//!   curve array of `sm_CreateSrfOfRevolution`, the port returns
//!   [`SolidError::InvalidInput`]. A curve of `sm_CreateSrfOfRevolution` that is not an
//!   `SmBSplineCurve` (the array's element type) is [`SolidError::InvalidInput`].
//! - **Surfaces.** `SM_FILLETSURF_TYPE` is `SmSurface` (`SM_FILLET_SURFS` is defined):
//!   the new surfaces are [`Surface`] values (an `SmPlane` stays a plane).
//! - **Methods.** `GetContext` and `GetFilletExecutive` have no port: the [`Kernel`] and
//!   the executive's [`FilletExecutiveData`] are passed to the methods that need them
//!   (`m_pExecutive`'s `m_vExtendedSurfacesMap` is
//!   [`FilletExecutiveData::extended_surfaces_map`]). The virtual methods dispatch on
//!   [`FilletCorner::kind`]; the `SmFilletNx1ClosedCorner` / `SmFilletNx2Corner`
//!   overrides of `GetTangentFilletSolver` belong to unit `blend.fillet_corner_p2` and
//!   are [`SolidError::Unported`] there. `SqueezeEdge` removes the deleted edge and
//!   vertex from the corner's lists at once instead of leaving the `NULL` slots that
//!   `CalcCornerEdgeGeom` compresses. `ReInitializeCornerEdges` applies
//!   `SmEdge::SetCurve(NULL, TRUE)`'s notification effects to the pseudo B-rep's
//!   records directly (no spatial cache: the pseudo B-rep has none).
//!
//! No loops depend on input size beyond the given arrays; there is no recursion.

use crate::blend::fillet_corner::p3::{FilletCorner, FilletCornerKind, FilletExecutiveData};
use crate::blend::fillet_geom::p1::{FilStatus, FilletBrep, FilletEdgeType, FilletSolverId};
use crate::curve::curve::p1::GlobalPointSolveOptions;
use crate::curve::curve::{p2 as curve2, Curve, CurveKind};
use crate::error::{Result, SolidError};
use crate::geom::extent2d::Extent2d;
use crate::geom::global_solver::{SolutionRequested, SolverOperation};
use crate::geom::vector2d::Point2;
use crate::kernel::Kernel;
use crate::math::{rad_to_deg, Point3, Vec3, EFF_ZERO, EFF_ZERO_SQ, EFF_ZERO_SQRT, PI};
use crate::nurbs::b_spline_curve::BSplineCurve;
use crate::nurbs::b_spline_surface::BSplineSurface;
use crate::nurbs::b_spline_surface_n_lib::p2::OrientType;
use crate::surface::plane::Plane;
use crate::surface::surface::Surface;
use crate::topology::edgeuse::{p1 as eu1, p3 as eu3};
use crate::topology::{edge, vertexuse, Brep, EdgeId, EdgeuseId, FaceId, ObjectHandle, VertexId};

pub use crate::blend::fillet_geom::p1::{check_seam_out, set_surface_closure};

const OP_G1: &str = "blend::fillet_corner::p1::test_g1_continuity";
const OP_FACE: &str = "blend::fillet_corner::p1::face_by_two_edgeuses";
const OP_TANGENTIAL: &str = "blend::fillet_corner::p1::check_tangential_edgeuse";
const OP_CONCAVITY: &str = "blend::fillet_corner::p1::check_face_complimentary_corner_concavity";
const OP_REVOLUTION: &str = "blend::fillet_corner::p1::create_srf_of_revolution";
const OP_PLANAR: &str = "blend::fillet_corner::p1::create_planar_surface";

/// The reason of the inputs upstream dereferences null on.
const NULL_DEREF: &str = "upstream dereferences null";

/// The angle limit of `sm_CheckFaceComplimentaryCornerConcavity` ("Here we test against
/// 2 degrees").
const CONCAVITY_ANGLE_LIMIT_DEG: f64 = 2.0;

/// A status upstream ignores: `Unported` propagates, anything else gives `fallback`.
fn ignored<T>(r: Result<T>, fallback: T) -> Result<T> {
    match r {
        Ok(v) => Ok(v),
        Err(e @ SolidError::Unported { .. }) => Err(e),
        Err(_) => Ok(fallback),
    }
}

/// `pCrv->Evaluate(t, 1, TRUE, aPV)` with the status ignored: position and tangent,
/// `SM_UNDEF_DOUBLE` vectors when the evaluation fails.
fn evaluate_tangent(curve: &Curve, t: f64) -> Result<[Vec3; 2]> {
    let pv = ignored(curve.evaluate(t, 1, true, true).map(Some), None)?;
    Ok(match pv.as_deref() {
        Some([p, d, ..]) => [*p, *d],
        _ => [Vec3::UNDEF, Vec3::UNDEF],
    })
}

/// `pEU->EvaluateBinormal(t, FALSE, pnt, bin, pOptTangent, pOptNormal)` with the
/// status ignored (`SE`): binormal, tangent and face use normal, `SM_UNDEF_DOUBLE`
/// vectors when it fails.
fn evaluate_binormal(k: &Kernel, b: &mut Brep, eu: EdgeuseId, t: f64) -> Result<[Vec3; 3]> {
    let r = eu3::evaluate_binormal(k, b, eu, t, false, &eu3::BinormalOptions::default());
    let r = ignored(r.map(Some), None)?;
    Ok(match r {
        Some(e) => [e.binormal, e.tangent, e.face_normal],
        None => [Vec3::UNDEF, Vec3::UNDEF, Vec3::UNDEF],
    })
}

/// `pEU->GetRadial()`, dereferenced by every caller.
fn radial(b: &Brep, eu: EdgeuseId, op: &'static str) -> Result<EdgeuseId> {
    eu1::radial(b, eu)?.ok_or_else(|| SolidError::invalid_input(op, NULL_DEREF))
}

/// `pEU->GetEdge()`, dereferenced by every caller.
fn edge_of(b: &Brep, eu: EdgeuseId, op: &'static str) -> Result<EdgeId> {
    eu1::edge(b, eu)?.ok_or_else(|| SolidError::invalid_input(op, NULL_DEREF))
}

/// `sm_TestG1Continuity`: whether the end tangents of `edge1` and `edge2` at `vert`
/// are within `ang_tol_deg` degrees of parallel (either sense). Returns the flag and
/// the angle in degrees (`pdAngDeg`, set only when both tangents are longer than
/// `SM_EFF_ZERO`). `false` when any of the three is missing.
///
/// Edge 1 is evaluated at the end of its interval at `vert` (its start when `vert` is
/// its start vertex); edge 2 at the start of its interval when `vert` is its start
/// vertex and it is not edge 1, else at the end.
///
/// # Errors
///
/// [`SolidError::InvalidInput`] for an edge without a curve (upstream dereferences
/// null); [`SolidError::Unported`] from the curve evaluations.
pub fn test_g1_continuity(
    b: &Brep,
    vert: Option<VertexId>,
    edge1: Option<EdgeId>,
    edge2: Option<EdgeId>,
    ang_tol_deg: f64,
) -> Result<(bool, Option<f64>)> {
    // no work - not given 3 valid pointers
    let (Some(vert), Some(edge1), Some(edge2)) = (vert, edge1, edge2) else {
        return Ok((false, None));
    };

    let crv1 =
        edge::curve(b, edge1)?.ok_or_else(|| SolidError::invalid_input(OP_G1, NULL_DEREF))?;
    let ivl1 = edge::interval(b, edge1)?;
    let param1 = if edge::start_vertex(b, edge1)? == Some(vert) {
        ivl1.min()
    } else {
        ivl1.max()
    };

    // get Edge1 tangent at vertex
    let pv1 = evaluate_tangent(crv1, param1)?;

    let crv2 =
        edge::curve(b, edge2)?.ok_or_else(|| SolidError::invalid_input(OP_G1, NULL_DEREF))?;
    let ivl2 = edge::interval(b, edge2)?;
    let param2 = if edge2 == edge1 || edge::start_vertex(b, edge2)? != Some(vert) {
        ivl2.max()
    } else {
        ivl2.min()
    };

    // get Edge2 tangent at vertex
    let pv2 = evaluate_tangent(crv2, param2)?;

    if pv1[1].length_squared() > EFF_ZERO_SQ && pv2[1].length_squared() > EFF_ZERO_SQ {
        // get angle (degrees) between Edge1 and Edge2 endTangents (SE: 0 on failure)
        let mut ang = ignored(pv1[1].angle_between(pv2[1]), 0.0)?;
        if ang > PI / 2.0 {
            ang = PI - ang;
        }
        let ang_deg = rad_to_deg(ang);

        // check for in tolerance angle
        return Ok((ang_deg < ang_tol_deg, Some(ang_deg)));
    }

    // endTangents are not within tolerance
    Ok((false, None))
}

/// `sm_GetFaceByTwoEdgeuses`: the first face shared by the sectors of `eu1` and `eu2`
/// (each edge use's face and its radial's face, in that order), skipping the face of
/// `eu3` when given. `None` when either edge use is missing or no face is shared.
///
/// # Errors
///
/// [`SolidError::InvalidInput`] for an edge use without a radial (upstream
/// dereferences null).
pub fn face_by_two_edgeuses(
    b: &Brep,
    eu1: Option<EdgeuseId>,
    eu2: Option<EdgeuseId>,
    eu3: Option<EdgeuseId>,
) -> Result<Option<FaceId>> {
    let (Some(eu1), Some(eu2)) = (eu1, eu2) else {
        return Ok(None);
    };

    // let sFaces1 = two faces attached to pEU1->Sector
    let faces1 = [eu1::face(b, eu1)?, eu1::face(b, radial(b, eu1, OP_FACE)?)?];
    // let sFaces2 = two faces attached to pEU2->Sector
    let faces2 = [eu1::face(b, eu2)?, eu1::face(b, radial(b, eu2, OP_FACE)?)?];
    let skip = match eu3 {
        Some(eu3) => Some(eu1::face(b, eu3)?),
        None => None,
    };

    for face1 in faces1 {
        // when asked skip pEU3->Face
        if skip == Some(face1) {
            continue;
        }
        for face2 in faces2 {
            // return 1st found common face
            if face1 == face2 {
                return Ok(face2);
            }
        }
    }
    Ok(None)
}

/// `sm_CheckTangentialEdgeuse`: whether the face use normals of `eu` and its radial
/// at the middle of the edge's interval are parallel within `ang_tol_deg` degrees
/// (the edge separates two tangential faces).
///
/// # Errors
///
/// [`SolidError::InvalidInput`] for an edge use without an edge or radial (upstream
/// dereferences null); [`SolidError::Unported`] from `EvaluateBinormal`.
pub fn check_tangential_edgeuse(
    k: &Kernel,
    b: &mut Brep,
    eu: EdgeuseId,
    ang_tol_deg: f64,
) -> Result<bool> {
    let e = edge_of(b, eu, OP_TANGENTIAL)?;
    let t = edge::interval(b, e)?.evaluate(0.5);

    let [_, _, norm1] = evaluate_binormal(k, b, eu, t)?;
    let rad = radial(b, eu, OP_TANGENTIAL)?;
    let [_, _, norm2] = evaluate_binormal(k, b, rad, t)?;

    // Edge is adjacent to two tangential faces when the normals are parallel
    Ok(norm1.is_parallel_to(norm2, ang_tol_deg))
}

/// The binormal and the corner-outward tangent of `eu` at `corner` (the shared steps
/// of `sm_CheckFaceComplimentaryCornerConcavity` for each edge use).
fn corner_binormal(
    k: &Kernel,
    b: &mut Brep,
    corner: VertexId,
    eu: EdgeuseId,
) -> Result<(Vec3, Vec3)> {
    let vu = eu1::vertexuse(b, eu)?
        .ok_or_else(|| SolidError::invalid_input(OP_CONCAVITY, NULL_DEREF))?;
    let vert = vertexuse::vertex(b, vu)?;
    let orient = eu1::orientation(b, eu)?;

    // get parameter on pEU corresponding to CornerVert
    let at_corner = vert == Some(corner);
    let t = if (at_corner && orient == OrientType::Opposite)
        || (!at_corner && orient == OrientType::Same)
    {
        1.0
    } else {
        0.0
    };

    // get pEU tangent (pointing from corner into edge) and binormal at cornerVertex
    let e = edge_of(b, eu, OP_CONCAVITY)?;
    let ivl = edge::interval(b, e)?;
    let [bin, mut tangent, _] = evaluate_binormal(k, b, eu, ivl.evaluate(t))?;
    if !at_corner {
        // our corner vertex is the end of the edge, so flip the tangent to point away
        tangent = -tangent;
    }
    Ok((bin, tangent))
}

/// `sm_CheckFaceComplimentaryCornerConcavity`: whether the face corner at `corner`
/// bounded by `test_eu1` and `test_eu2` (both pointing out of the corner) is concave
/// (its inside angle runs from 180 degrees plus tolerance to 360). With `ref_eu` (an
/// edge use lying between them, which must separate tangential faces within 2 degrees,
/// else `false`) the face is `ref_eu`'s; without it the face the two edge uses share.
/// Each test edge use is replaced by its radial when the radial lies on that face; the
/// corner is concave when the dot product of each unit binormal with the other edge
/// use's unit tangent is below `-2 * PI / 180` (upstream compares it to the angle limit
/// in radians).
///
/// # Errors
///
/// [`SolidError::InvalidInput`] for an edge use without a radial, vertex use or edge
/// (upstream dereferences null); [`SolidError::Unported`] from `EvaluateBinormal`.
pub fn check_face_complimentary_corner_concavity(
    k: &Kernel,
    b: &mut Brep,
    corner: VertexId,
    test_eu1: EdgeuseId,
    test_eu2: EdgeuseId,
    ref_eu: Option<EdgeuseId>,
) -> Result<bool> {
    // If pOptRefEU is given, it needs to be a tangential edge
    if let Some(ref_eu) = ref_eu {
        if !check_tangential_edgeuse(k, b, ref_eu, CONCAVITY_ANGLE_LIMIT_DEG)? {
            return Ok(false);
        }
    }

    // get face spanning the edgeuse->Edges or any old face connected to the pOptRefEU->edge
    let common_face = match ref_eu {
        Some(ref_eu) => eu1::face(b, ref_eu)?,
        None => face_by_two_edgeuses(b, Some(test_eu1), Some(test_eu2), None)?,
    };

    // let pEU2 and pEU1 be Edgeuses to common face from their respective edges
    let mut eu_a = test_eu1;
    let mut eu_b = test_eu2;
    let radial2 = radial(b, test_eu2, OP_CONCAVITY)?;
    if common_face == eu1::face(b, radial2)? {
        eu_b = radial2;
    }
    let radial1 = radial(b, test_eu1, OP_CONCAVITY)?;
    if common_face == eu1::face(b, radial1)? {
        eu_a = radial1;
    }

    let (bin1, tangent1) = corner_binormal(k, b, corner, eu_a)?;
    let (bin2, tangent2) = corner_binormal(k, b, corner, eu_b)?;

    // Check if binormal of each edgeuse is pointing opposite way of the tangent of the
    // other edgeuse (Unitize statuses ignored)
    let bin1 = bin1.unitize_or_z();
    let bin2 = bin2.unitize_or_z();
    let tangent1 = tangent1.unitize_or_z();
    let tangent2 = tangent2.unitize_or_z();

    let val1 = bin1.dot(tangent2);
    let val2 = bin2.dot(tangent1);

    // SM_DEG2RAD(-dAngleLimitDegrees)
    let limit_rad = (-CONCAVITY_ANGLE_LIMIT_DEG) * (PI / 180.0);
    Ok(val1 < limit_rad && val2 < limit_rad)
}

/// The `SmBSplineCurve` part of `curve` (`SmTArray<SmBSplineCurve*>` elements: the
/// analytic curves derive from it upstream).
fn bspline_part(curve: &Curve) -> Result<&BSplineCurve> {
    match curve.kind() {
        CurveKind::BSpline(c) => Ok(c),
        CurveKind::Line(c) => Ok(&c.base),
        CurveKind::Ellipse(c) => Ok(c.base()),
        CurveKind::Hyperbola(c) => Ok(c.base()),
        CurveKind::Parabola(c) => Ok(c.base()),
        _ => Err(SolidError::invalid_input(
            OP_REVOLUTION,
            "boundary curve is not an SmBSplineCurve",
        )),
    }
}

/// `sm_CreateSrfOfRevolution`: the surface of revolution filling the counter-clockwise
/// loop `curves` (end to end), when `curves[major]` is an arc (`IsArc` with 5 samples
/// and tolerance `tol_3d`) whose axis and angle span define the revolution. With four
/// curves the opposite curve must also be an arc on a colinear axis with an angle span
/// within 2 degrees. The curve before the major arc is revolved by 1.1 times the arc's
/// span; the curve after it must lie on the result (`SimpleCoincidenceChecker` within
/// `tol_3d`).
///
/// # Errors
///
/// `SM_ERR` when a required curve is not an arc, the axes or spans disagree, or the
/// next curve does not lie on the surface; [`SolidError::InvalidInput`] for an empty
/// array, an index outside it (upstream reads past the array) or a curve that is not
/// an `SmBSplineCurve`; those of `CreateSurfOfRevolution` and the coincidence check.
pub fn create_srf_of_revolution(
    k: &Kernel,
    curves: &[Curve],
    major: usize,
    tol_3d: f64,
) -> Result<Surface> {
    let total = curves.len();
    let out_of_range = || SolidError::invalid_input(OP_REVOLUTION, "curve index outside the array");
    let major_arc = curves.get(major).ok_or_else(out_of_range)?;

    // no work - can't make SurfOfRevolution when MajorArc is not an Arc
    let Some((frame, _rad, start_ang, end_ang)) = major_arc.is_arc(5, tol_3d)? else {
        return Err(SolidError::sm_err(OP_REVOLUTION));
    };
    let origin = frame.origin();
    let ref_z_axis = frame.z_axis();
    let mut angle_span = end_ang - start_ang;

    // when given exactly 4 curves - make sure OppositeCurve->AxisOfRevolution ==
    // MajorCurve->AxisOfRevolution and have similar ArcAngle lengths
    if total == 4 {
        let other_arc = curves.get((major + 2) % 4).ok_or_else(out_of_range)?;
        let Some((frame1, _rad1, start_ang1, end_ang1)) = other_arc.is_arc(5, tol_3d)? else {
            return Err(SolidError::sm_err(OP_REVOLUTION));
        };
        let origin1 = frame1.origin();
        let ref_z_axis1 = frame1.z_axis();
        let angle_span1 = end_ang1 - start_ang1;

        // Quit: if two Arcs don't share Colinear AxesOfRotation and similar ArcAngle lengths
        if !ref_z_axis.is_colinear_with(ref_z_axis1, origin, origin1)
            || (angle_span - angle_span1).abs() > 2.0
        {
            return Err(SolidError::sm_err(OP_REVOLUTION));
        }
    }

    // Increase angle span by 10% to make sure it's big enough
    angle_span *= 1.1;

    // Let SurfOfRevolution GenCurve = Curve CW from MajorArcCurve
    let gen_index = (major + total - 1)
        .checked_rem(total)
        .ok_or_else(out_of_range)?;
    let gen_curve = bspline_part(curves.get(gen_index).ok_or_else(out_of_range)?)?;

    // Sweep GenCurve about MajorArcCurve AxisOfRevolution to make the fillet surface
    let surface = Surface::from(BSplineSurface::create_surf_of_revolution(
        k, gen_curve, origin, ref_z_axis, angle_span, false, None,
    )?);

    // let OtherCurve = Curve opposite to GenCurve
    let other_curve = curves.get((major + 1) % total).ok_or_else(out_of_range)?;

    // Check for OtherCurve/SurfOfRevolution coincidence
    let ivl = other_curve.natural_interval();
    let uv_domain = surface.natural_uv_domain();
    let found = surface.simple_coincidence_checker::<ObjectHandle>(
        k,
        &uv_domain,
        other_curve,
        &ivl,
        tol_3d,
        false,
        [None, None],
    )?;

    // when SurfOfRevolution is not coincident with OtherCurve - don't use it
    if found.is_none() {
        return Err(SolidError::sm_err(OP_REVOLUTION));
    }
    Ok(surface)
}

/// `sm_CreatePlanarSurface`: the planar corner patch through `plane_pnts`. Three points
/// give an `SmPlane` through their centroid with normal `(p1 - p0) x (p2 - p0)` and the
/// STEP domain `[-d, d]^2`, `d = |p1 - p0| + |p2 - p0|` (status ignored); four give the
/// bilinear B-spline surface with corners `p0, p1, p2, p3` in order around it.
/// `tol_3d` is unused upstream.
///
/// # Errors
///
/// `SM_ERR_INVALID_INPUT` (`Unitize`) for three collinear points; `SM_ERR` for any
/// other point count ("Not yet implemented"); those of `CreateBilinearSurface`.
pub fn create_planar_surface(k: &Kernel, plane_pnts: &[Point3], tol_3d: f64) -> Result<Surface> {
    let _ = tol_3d; // SM_REF1

    match plane_pnts {
        // when 3 points are given - we can always build an exactly interpolating plane
        [p0, p1, p2] => {
            // Make regular plane - centered on Point Centroid and oriented by given points
            let origin = (*p0 + *p1 + *p2) / 3.0;
            let x_axis = *p1 - *p0;
            let y_axis = *p2 - *p0;
            let z_axis = x_axis.cross(y_axis);

            // signal error for degenerate point set
            let z_axis = z_axis.unitize()?;

            // create infinite plane
            let mut plane = Surface::from(Plane::from_normal(origin, z_axis));

            // Adjust the newPlane->UVDomain to fit the given points (status ignored)
            let size = x_axis.length() + y_axis.length();
            let ext = Extent2d::from_points(Point2::new(-size, -size), Point2::new(size, size));
            ignored(plane.adjust_step_uv_domain(k, &ext), ())?;
            Ok(plane)
        }
        // Make NonRational Bilinear BSplineSurface interpolating given corners
        [p0, p1, p2, p3] => Ok(Surface::from(BSplineSurface::create_bilinear_surface(
            *p0, *p1, *p3, *p2,
        )?)),
        // other cases not yet implemented
        _ => Err(SolidError::sm_err(OP_PLANAR)),
    }
}

/// `sm_GetFilletRadius`: the radius and angle span (radians) of the circular arc
/// tangent to the two surfaces at `surf_pnt` / `surf_pnt2` with normals `surf_normal` /
/// `surf_normal2` (flipped when `offset_orientation > 0`): the counter-clockwise angle
/// between the normals about their cross product (folded to `PI - angle` above `PI`)
/// and `radius = |p2 - p1| / 2 / sin(angle / 2)`.
///
/// # Errors
///
/// None in practice (see the module docs on `CCWAngleBetween`).
pub fn fillet_radius(
    surf_pnt: Point3,
    surf_normal: Vec3,
    surf_pnt2: Point3,
    surf_normal2: Vec3,
    offset_orientation: f64,
) -> Result<(f64, f64)> {
    let (mut normal1, mut normal2) = (surf_normal, surf_normal2);
    if offset_orientation > 0.0 {
        normal1 = -normal1;
        normal2 = -normal2;
    }
    // Unitize status ignored
    let axis = normal1.cross(normal2).unitize_or_z();
    let mut ang_rad = axis.ccw_angle_between(normal1, normal2)?;
    if ang_rad > PI {
        ang_rad = PI - ang_rad;
    }
    let leng = surf_pnt.distance(surf_pnt2);
    let radius = leng / 2.0 / (ang_rad / 2.0).sin();
    Ok((radius, ang_rad))
}

// ---- SmFilletCorner methods: data access, predicates and the small helpers ---------------

const OP_SOLVERS: &str = "blend::fillet_corner::p1::fillet_solvers";
const OP_REINIT: &str = "blend::fillet_corner::p1::re_initialize_corner_edges";
const OP_SQUEEZE: &str = "blend::fillet_corner::p1::squeeze_edge";
const OP_RAIL_RAIL: &str = "blend::fillet_corner::p1::test_rail_rail_interpolation_curve";
const OP_NX1: &str = "blend::fillet_corner::p1::check_nx1_degeneracy";

/// `SmFilletCorner_TYPE` (`TOPO_BASE_TYPE + 36`).
pub const FILLET_CORNER_TYPE: i64 = 16036;
/// `SmFilletCorner::GetTypeString` (`SM_COMMON`).
pub const FILLET_CORNER_TYPE_STRING: &str = "SmFilletCorner_TYPE";

impl FilletCorner {
    /// `GetThisApproxTol3d` (its `SM_ASSERT_TOL` is debug-only). Also the
    /// `SmFilletCorner` branch of `SmTol::GetApproxTol3d(const SmObject*)`.
    #[must_use]
    pub const fn this_approx_tol_3d(&self) -> f64 {
        self.approx_tol_3d
    }

    /// `GetFilletEdges`: appends `m_vEdges` to `edges`.
    pub fn fillet_edges(&self, edges: &mut Vec<EdgeId>) {
        edges.extend_from_slice(&self.edges);
    }

    /// `GetFilletVertices`: appends `m_vVertices` to `verts`.
    pub fn fillet_vertices(&self, verts: &mut Vec<VertexId>) {
        verts.extend_from_slice(&self.vertices);
    }

    /// `GetFilletedVertex` (`m_cpVertex`).
    #[must_use]
    pub const fn filleted_vertex(&self) -> Option<VertexId> {
        self.vertex
    }

    /// `GetSetBackDist`.
    #[must_use]
    pub const fn set_back_dist(&self) -> f64 {
        self.set_back_dist
    }

    /// `GetStatus`.
    #[must_use]
    pub const fn status(&self) -> FilStatus {
        self.status
    }

    /// `GetType` (`SM_COMMON`): [`FILLET_CORNER_TYPE`] for every derived class.
    #[must_use]
    pub const fn type_code(&self) -> i64 {
        FILLET_CORNER_TYPE
    }

    /// `GetTypeString` (`SM_COMMON`).
    #[must_use]
    pub const fn type_string(&self) -> &'static str {
        FILLET_CORNER_TYPE_STRING
    }

    /// `IsKindOf`: an `SmFilletCorner` or what `SmObject::IsKindOf` accepts.
    #[must_use]
    pub const fn is_kind_of(&self, t: i64) -> bool {
        t == FILLET_CORNER_TYPE || crate::geom::object::is_kind_of(t)
    }

    /// `AddFilletVertex`.
    pub fn add_fillet_vertex(&mut self, v: VertexId) {
        self.vertices.push(v);
    }

    /// `SetBevel`.
    pub fn set_bevel(&mut self, is_bevel: bool) {
        self.bevel = is_bevel;
    }

    /// `SetCalcAsConstRad`.
    pub fn set_calc_as_const_rad(&mut self, is_const: bool) {
        self.calc_as_const_rad = is_const;
    }

    /// `SmFilletConvexNxNCorner::SetChamfer` (the member is `SmFilletCorner`'s; the
    /// setter is only declared on the convex `NxN` corner upstream).
    pub fn set_chamfer(&mut self, is_chamfer: bool) {
        self.chamfer = is_chamfer;
    }

    /// `IsBevelCorner`.
    #[must_use]
    pub const fn is_bevel_corner(&self) -> bool {
        self.bevel
    }

    /// `IsBlendingCorner`.
    #[must_use]
    pub const fn is_blending_corner(&self) -> bool {
        self.blending
    }

    /// `IsDegenerateCorner`.
    #[must_use]
    pub const fn is_degenerate_corner(&self) -> bool {
        self.degenerate
    }

    /// `IsExtendedSurfacePatch`.
    #[must_use]
    pub const fn is_extended_surface_patch(&self) -> bool {
        self.extended_surface_patch
    }

    /// `IsCalcAsConstRad`.
    #[must_use]
    pub const fn is_calc_as_const_rad(&self) -> bool {
        self.calc_as_const_rad
    }

    /// `IsTangentCorner` (virtual): `m_bG1Case` for the closed Nx1 and the Nx2 corner,
    /// `FALSE` for every other kind.
    #[must_use]
    pub const fn is_tangent_corner(&self) -> bool {
        match self.kind {
            FilletCornerKind::NByOneClosed | FilletCornerKind::NByTwo => self.g1_case,
            _ => false,
        }
    }

    /// `GetTangentFilletSolver(pFS)` (virtual): the solver of the other filleted edge
    /// use of a 2x2 corner, the solver itself for a 1x1 corner, `None` (`SM_ASSERT_ERR`)
    /// when the corner has another number of filleted edge uses or `fs` is not one of
    /// their solvers; `None` for the base class (`SM_REF1(pFS)`).
    ///
    /// # Errors
    ///
    /// [`SolidError::InvalidInput`] for a 1x1 or 2x2 corner without `m_pSolverEUs`
    /// (upstream dereferences it); [`SolidError::Unported`] for the closed Nx1 and Nx2
    /// corners (their overrides are unit `blend.fillet_corner_p2`).
    pub fn tangent_fillet_solver(
        &self,
        exec: &FilletExecutiveData,
        fs: FilletSolverId,
    ) -> Result<Option<FilletSolverId>> {
        const OP: &str = "blend::fillet_corner::p1::tangent_fillet_solver";
        match self.kind {
            FilletCornerKind::OneByOne => {
                let eus = self
                    .solver_eus
                    .as_deref()
                    .ok_or_else(|| SolidError::invalid_input(OP, NULL_DEREF))?;
                let [eu] = eus else {
                    return Ok(None);
                };
                let my_fs = exec.fillet_solver_of_edgeuse(*eu);
                if my_fs != Some(fs) {
                    return Ok(None);
                }
                Ok(my_fs)
            }
            FilletCornerKind::TwoByTwo => {
                let eus = self
                    .solver_eus
                    .as_deref()
                    .ok_or_else(|| SolidError::invalid_input(OP, NULL_DEREF))?;
                let [eu0, eu1_] = eus else {
                    return Ok(None);
                };
                let fs1 = exec.fillet_solver_of_edgeuse(*eu0);
                let fs2 = exec.fillet_solver_of_edgeuse(*eu1_);
                if fs1 == Some(fs) {
                    return Ok(fs2);
                }
                if fs2 == Some(fs) {
                    return Ok(fs1);
                }
                Ok(None)
            }
            FilletCornerKind::NByOneClosed => {
                // UNPORTED(blend.fillet_corner_p2): SmFilletNx1ClosedCorner::GetTangentFilletSolver
                Err(SolidError::unported(
                    "SmFilletNx1ClosedCorner::GetTangentFilletSolver",
                    "blend.fillet_corner_p2",
                ))
            }
            FilletCornerKind::NByTwo => {
                // UNPORTED(blend.fillet_corner_p2): SmFilletNx2Corner::GetTangentFilletSolver
                Err(SolidError::unported(
                    "SmFilletNx2Corner::GetTangentFilletSolver",
                    "blend.fillet_corner_p2",
                ))
            }
            _ => Ok(None),
        }
    }

    /// `GetFilletSolvers`: one entry per filleted edge use, each the solver of the
    /// FIRST filleted edge use (`(*m_pSolverEUs)[0]` inside the loop, as upstream).
    ///
    /// # Errors
    ///
    /// [`SolidError::InvalidInput`] without `m_pSolverEUs` (upstream dereferences it).
    pub fn fillet_solvers(
        &self,
        exec: &FilletExecutiveData,
    ) -> Result<Vec<Option<FilletSolverId>>> {
        let eus = self
            .solver_eus
            .as_deref()
            .ok_or_else(|| SolidError::invalid_input(OP_SOLVERS, NULL_DEREF))?;
        let mut solvers = Vec::new();
        for _ in eus {
            let first = eus
                .first()
                .ok_or_else(|| SolidError::invalid_input(OP_SOLVERS, "empty solver edge uses"))?;
            solvers.push(exec.fillet_solver_of_edgeuse(*first));
        }
        Ok(solvers)
    }

    /// `GetExtendedSurface(pOriginalFace)`: the extended surface registered in the
    /// executive for a face of the original B-rep, `None` when there is none.
    #[must_use]
    pub fn extended_surface<'a>(
        &self,
        exec: &'a FilletExecutiveData,
        original_face: FaceId,
    ) -> Option<&'a Surface> {
        exec.extended_surfaces_map.get(&original_face)
    }

    /// `AddExtendedOriginalSurface(pOriginalFace, pExtSurface)`: registers the pair in
    /// the executive's `m_vExtendedSurfacesMap`. When the face already has one, the
    /// surface with the larger natural domain is kept (the new one is dropped when its
    /// domain is contained in the old one's within `SM_EFF_ZERO`).
    pub fn add_extended_original_surface(
        &self,
        exec: &mut FilletExecutiveData,
        original_face: FaceId,
        ext_surface: Surface,
    ) {
        if let Some(old_surface) = exec.extended_surfaces_map.get(&original_face) {
            let new_domain = ext_surface.natural_uv_domain();
            let old_domain = old_surface.natural_uv_domain();
            if new_domain.is_contained_by(&old_domain, EFF_ZERO) {
                // Old domain was bigger, delete the new ext surface.
                return;
            }
            // else: the old surface is deleted (replaced below)
        }
        exec.extended_surfaces_map
            .insert(original_face, ext_surface);
    }

    /// `ReInitializeCornerEdges`: `SetCurve(NULL, TRUE)` on every corner fillet edge
    /// of the pseudo B-rep: the curve is dropped and, as the geometry-change
    /// notification does, the UV trim curves of its edge uses are deleted and the
    /// cached gaps of its edge uses, their mates and its vertex uses reset. An edge
    /// without a curve is left alone (`SmEdge::SetCurve` returns at once).
    ///
    /// # Errors
    ///
    /// [`SolidError::InvalidInput`] / `SM_ERR_NULL_POINTER` where the notification
    /// dereferences null (an edge use without mate or vertex use, a null vertex use);
    /// [`SolidError::NotFound`] for a dangling id.
    pub fn re_initialize_corner_edges(&self, pseudo: &mut FilletBrep) -> Result<()> {
        let b = pseudo.brep_mut();
        for &e in &self.edges {
            if b.edge(e)?.curve.is_none() {
                continue;
            }
            // side effect: delete all UVTrimCurves and the cached gaps (SmEdge::Notify,
            // SM_NO_CHANGE_GEOMETRY)
            for eu in edge::edgeuses(b, e, None, None)? {
                eu1::set_uv_trim_curve(b, eu, eu1::UvTrimCurveArg::Null, 0.0, true, None)?;
                let m = eu1::mate(b, eu)?.ok_or_else(|| {
                    SolidError::invalid_input(OP_REINIT, "Notify: edge use without mate")
                })?;
                b.edgeuse_mut(eu)?.max_edge_face_gap.reset();
                b.edgeuse_mut(m)?.max_edge_face_gap.reset();
            }
            for vu in edge::vertexuses(b, e, None)? {
                let vu = vu.ok_or_else(|| SolidError::null_pointer(OP_REINIT))?;
                let r = b.vertexuse_mut(vu)?;
                r.vertex_edge_gap.reset();
                r.vertex_face_gap.reset();
            }
            b.edge_mut(e)?.curve = None;
        }
        Ok(())
    }

    /// `SqueezeEdge(pEdgeToSqueeze)`: collapses a fillet edge of the pseudo B-rep onto
    /// its surviving vertex. The vertex marked `SM_FV_TO_BE_DELETED` is deleted: the
    /// vertex uses of its rail edges move to the survivor, its other edges (corner
    /// edges of this corner) are deleted, then the vertex itself.
    ///
    /// Deviation: upstream leaves `NULL` slots in `m_vEdges` / `m_vVertices` that
    /// `CalcCornerEdgeGeom` compresses afterwards; here the entries are removed at
    /// once (the relative order is the same).
    ///
    /// # Errors
    ///
    /// `SM_ERR` when neither vertex is to be deleted or an edge / the vertex is not
    /// this corner's; [`SolidError::InvalidInput`] where upstream dereferences null;
    /// those of the deletions.
    pub fn squeeze_edge(&mut self, pseudo: &mut FilletBrep, edge_to_squeeze: EdgeId) -> Result<()> {
        let null = || SolidError::invalid_input(OP_SQUEEZE, NULL_DEREF);
        // Get the vertex that will be deleted
        let mut surviving = edge::vertex(pseudo.brep(), edge_to_squeeze)?.ok_or_else(null)?;
        let mut delete = edge::other_vertex(pseudo.brep(), edge_to_squeeze, Some(surviving))?
            .ok_or_else(null)?;
        if pseudo.fillet_vertex(delete)?.status != FilStatus::ToBeDeleted {
            core::mem::swap(&mut surviving, &mut delete);
        }
        if pseudo.fillet_vertex(delete)?.status != FilStatus::ToBeDeleted {
            return Err(SolidError::sm_err(OP_SQUEEZE));
        }

        let vus = crate::topology::vertex::vertexuses(pseudo.brep(), delete)?;
        for vu in vus {
            let eu = eu1::vertexuse_edgeuse(pseudo.brep(), vu)?.ok_or_else(null)?;
            let fe = eu1::edge(pseudo.brep(), eu)?.ok_or_else(null)?;
            if pseudo.fillet_edge(fe)?.edge_type == FilletEdgeType::Rail {
                // pDeleteFV->Remove(pVU); pSurvivingFV->PostInsert(pVU) (the edge record's
                // vertex ends, derived from the vertex uses upstream, follow)
                let b = pseudo.brep_mut();
                b.vertex_mut(delete)?.vertexuses.retain(|&x| x != vu);
                b.vertexuse_mut(vu)?.vertex = Some(surviving);
                b.vertex_mut(surviving)?.vertexuses.push(vu);
                for end in &mut b.edge_mut(fe)?.vertices {
                    if *end == delete {
                        *end = surviving;
                    }
                }
            } else {
                let found = self
                    .edges
                    .iter()
                    .position(|&x| x == fe)
                    .ok_or_else(|| SolidError::sm_err(OP_SQUEEZE))?;
                self.edges.remove(found);
                pseudo.delete_fillet_edge(fe)?;
            }
        }
        let found = self
            .vertices
            .iter()
            .position(|&x| x == delete)
            .ok_or_else(|| SolidError::sm_err(OP_SQUEEZE))?;
        self.vertices.remove(found);
        pseudo.delete_fillet_vertex(delete)?;

        // all done
        Ok(())
    }

    /// `SmFilletConcaveNxNCorner::TestRailRailInterpolationCurve(pCurve)`: whether the
    /// corner vertex lies on the concave side of `curve` at its closest point (the
    /// vector from the closest point to the vertex opposes the curvature vector).
    ///
    /// # Errors
    ///
    /// [`SolidError::InvalidInput`] without `m_cpVertex` (upstream dereferences it);
    /// `SM_ERR` when the closest-point solve does not give exactly one solution; those
    /// of the solve, the evaluation and the unitisations (zero curvature or a vertex on
    /// the curve).
    pub fn test_rail_rail_interpolation_curve(
        &self,
        k: &Kernel,
        target: &Brep,
        curve: &Curve,
    ) -> Result<bool> {
        // Check to see if the original corner 'intersects' the curve
        // Drop corner vertex onto pCurve
        let vertex = self
            .vertex
            .ok_or_else(|| SolidError::invalid_input(OP_RAIL_RAIL, NULL_DEREF))?;
        let corner_pnt = target.vertex(vertex)?.point;
        let ivl = curve.natural_interval();
        let solutions = curve.global_point_solve::<ObjectHandle>(
            k,
            &ivl,
            SolverOperation::Minimize,
            corner_pnt,
            self.approx_tol_3d,
            SolutionRequested::Single,
            &GlobalPointSolveOptions::default(),
        )?;
        if solutions.len() != 1 {
            return Err(SolidError::sm_err(OP_RAIL_RAIL));
        }
        let t = solutions
            .get(0)
            .and_then(|s| s.start.parameters.first().copied())
            .ok_or_else(|| SolidError::sm_err(OP_RAIL_RAIL))?;
        // Find the 2nd deriv. of the curve at dT
        let geom = curve2::evaluate_geometric(curve, t, 2, true)?;
        let (Some(pos), Some(curvature)) = (geom.first(), geom.get(2)) else {
            return Err(SolidError::sm_err(OP_RAIL_RAIL));
        };
        let v = corner_pnt - *pos;
        let curvature = curvature.unitize()?;
        let v = v.unitize()?;
        Ok(v.dot(curvature) < 0.0)
    }
}

/// `SmFilletNx1Corner::CheckNx1Degeneracy(pCornerVert, pFilEdgeuse)` (static): whether
/// the two side edge uses next to the filleted edge use at the corner vertex (its
/// clockwise neighbour and its radial's counter-clockwise neighbour; the mate of the
/// first when they coincide, the 2x1 case) have parallel face use normals at the
/// corner end of their edges (within `SM_EFF_ZERO_SQRT`).
///
/// # Errors
///
/// [`SolidError::InvalidInput`] where upstream dereferences null (a missing neighbour,
/// radial, mate, vertex use or edge); those of the binormal evaluations
/// (`SmEdgeuse::EvaluateBinormal`, `SER`).
pub fn check_nx1_degeneracy(
    k: &Kernel,
    b: &mut Brep,
    corner_vert: VertexId,
    fil_edgeuse: EdgeuseId,
) -> Result<bool> {
    let null = || SolidError::invalid_input(OP_NX1, NULL_DEREF);
    // Get adjacent edgeuses of filleted edgeuse
    let side_eu1 = eu1::cw_edgeuse(b, fil_edgeuse)?.ok_or_else(null)?;
    let radial = eu1::radial(b, fil_edgeuse)?.ok_or_else(null)?;
    let mut side_eu2 = eu1::ccw_edgeuse(b, radial)?.ok_or_else(null)?;
    if side_eu1 == side_eu2 {
        // likely 2x1 case
        side_eu2 = eu1::mate(b, side_eu1)?.ok_or_else(null)?;
    }

    // the parameter of each side edge at the corner end
    let corner_param = |b: &Brep, eu: EdgeuseId| -> Result<f64> {
        let vu = eu1::vertexuse(b, eu)?.ok_or_else(null)?;
        let vert = b.vertexuse(vu)?.vertex;
        let e = eu1::edge(b, eu)?.ok_or_else(null)?;
        let ivl = edge::interval(b, e)?;
        let orientation = eu1::orientation(b, eu)?;
        let mut t = ivl.max();
        if (vert == Some(corner_vert) && orientation == OrientType::Same)
            || (vert != Some(corner_vert) && orientation == OrientType::Opposite)
        {
            t = ivl.min();
        }
        Ok(t)
    };
    let t1 = corner_param(b, side_eu1)?;
    let t2 = corner_param(b, side_eu2)?;
    let options = eu3::BinormalOptions::default();
    let faceuse_normal1 = eu3::evaluate_binormal(k, b, side_eu1, t1, false, &options)?.face_normal;
    let faceuse_normal2 = eu3::evaluate_binormal(k, b, side_eu2, t2, false, &options)?.face_normal;
    // Check if faceuse normals are parallel to each other
    Ok(faceuse_normal1.is_parallel_to(faceuse_normal2, EFF_ZERO_SQRT))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::SM_ERR;

    #[test]
    fn planar_surface_point_counts() {
        let k = Kernel::default();
        let p = [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ];
        assert!(create_planar_surface(&k, &p[..3], 0.0).is_ok());
        assert!(create_planar_surface(&k, &p, 0.0).is_ok());
        let e = create_planar_surface(&k, &p[..2], 0.0).unwrap_err();
        assert!(matches!(e, SolidError::Upstream { code: SM_ERR, .. }));
        let line = [p[0], p[1], Point3::new(2.0, 0.0, 0.0)];
        assert!(create_planar_surface(&k, &line, 0.0).is_err());
    }

    #[test]
    fn fillet_radius_quarter_turn() {
        // normals at right angles, points sqrt(2) apart: radius 1, span PI/2
        let (r, a) = fillet_radius(
            Point3::new(1.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            -1.0,
        )
        .unwrap();
        assert!((a - PI / 2.0).abs() < 1e-12);
        assert!((r - 1.0).abs() < 1e-12);
    }

    #[test]
    fn missing_inputs() {
        let b = Brep::default();
        assert_eq!(
            test_g1_continuity(&b, None, None, None, 1.0).unwrap(),
            (false, None)
        );
        assert_eq!(face_by_two_edgeuses(&b, None, None, None).unwrap(), None);
        let k = Kernel::default();
        assert!(create_srf_of_revolution(&k, &[], 0, 0.0).is_err());
    }
}
