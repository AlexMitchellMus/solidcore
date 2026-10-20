//! Fillet solver, part 2: analytic fillet geometry (cylinder, cone, torus), the
//! plane, rail/edge-use and rail/rail point solves, the cliff roll-over side
//! geometry, corner lookup and error recording.
//!
//! Upstream: `source/SMLib/src/SmFilletSolver.cpp`, `source/SMLib/inc/SmFilletSolver.h`.
//!
//! | Rust                                                     | Upstream                                     |
//! |----------------------------------------------------------|----------------------------------------------|
//! | [`FilletSolver::fillet_corner`]                          | `SmFilletSolver::GetFilletCorner`            |
//! | [`FilletSolver::calc_cylinder_fillet_geom`], [`FilletSolver::cylinder_fillet_geom_parts`] | `SmFilletSolver::CalcCylinderFilletGeom` |
//! | [`FilletSolver::calc_cone_fillet_geom`], [`FilletSolver::cone_fillet_geom_parts`]         | `SmFilletSolver::CalcConeFilletGeom`     |
//! | [`FilletSolver::calc_torus_fillet_geom`], [`FilletSolver::torus_fillet_geom_parts`]       | `SmFilletSolver::CalcTorusFilletGeom`    |
//! | [`FilletSolver::point_on_plane_solve`]                   | `SmFilletSolver::PointOnPlaneSolve`          |
//! | [`FilletSolver::rail_edgeuse_intersect`]                 | `SmFilletSolver::RailEdgeuseIntersect`       |
//! | [`FilletSolver::rail_rail_intersect`]                    | `SmFilletSolver::RailRailIntersect`          |
//! | [`find_rollover_side_edge`]                              | `sm_FindRolloverSideEdge` (static)           |
//! | [`FilletSolver::create_cliff_roll_over_side_fg`]         | `SmFilletSolver::CreateCliffRollOverSideFG`  |
//! | [`FilletSolver::record_fillet_error`]                    | `SmFilletSolver::RecordFilletError`          |
//!
//! Deviations:
//! - **Stacked on `blend.fillet_solver_p1`.** Everything here is a method of p1's
//!   [`FilletSolver`]; methods that call a virtual take the solver as
//!   `&mut dyn FilletSolverOps` (as p1's `ComputeSurfaceValues`). The p1 function
//!   objects are wrapped privately into [`EvalNFunctions`] for [`LocalSolveNd`].
//! - **Fillet geometries.** The ported `SmFilletGeom`, `SmFilletEdge` and
//!   `SmFilletVertex` (`blend.fillet_geom_p1`/`_p2`, partial) live in a `FilletBrep`,
//!   which this solver does not hold (`plan/findings.md`, "blend.fillet_solver_p2:
//!   fillet-geometry stubs need a FilletBrep"); those calls, the pseudo-B-rep and the
//!   executive return [`SolidError::Unported`] at the call site (tagged `UNPORTED`).
//!   `sm_SrfSrfIntersection` (a free function) is wired. The analytic
//!   `Calc*FilletGeom` methods build the fillet surface, rail curves and centre line
//!   first ([`FilletGeomParts`], `*_fillet_geom_parts`) and then store them into the
//!   fillet geometry (`SmFilletGeom::m_pFilletSurface`, `SmFilletEdge::SetCurve`, the
//!   unported step); upstream interleaves the two, which changes only the unported
//!   geometry's state on a failure.
//! - **Arguments.** Out-params become return values (`rbFoundSolution` /
//!   `rbFoundIntersection` the returned `bool` / `Option`); `SmContext` is the
//!   [`Kernel`]; `SmFilletErrorType` is its raw value; `SmFilletCorner*` is an opaque
//!   [`FilletCornerRef`].
//! - **Undefined inputs.** Where upstream dereferences `NULL` (missing edge use, edge,
//!   offset surface, surface generator, executive, UV trim curve), reads an
//!   uninitialised `lSurf2Index`, indexes past an array (`rTsectPnts` empty, a rail
//!   index above 1, `m_adUserDoubles` past `SM_MAX_USER_DOUBLES` in
//!   `PointOnPlaneSolve`, `pCC` past its end in `sm_FindRolloverSideEdge`) or wraps an
//!   unsigned size (`RailRailIntersect` with Jacobian sizes summing below 2), the port
//!   returns [`SolidError::InvalidInput`]. Writes of the solver intervals past the
//!   Jacobian size (upstream: `SM_ASSERT`, then into the fixed `SmExtentNd` storage)
//!   are dropped; the solve then rejects intervals not sized like the guess.
//! - **Iteration.** `sm_FindRolloverSideEdge`'s `while (TRUE)` advances its interval
//!   index on every pass that does not return, so it runs at most `pCC.len()` times;
//!   there is no recursion.

use crate::blend::fillet_custom_solver::OffsetSurface;
use crate::blend::fillet_geom::p2::{srf_srf_intersection, SrfSrfOptions};
use crate::blend::fillet_solver::p1::{
    check_edge_curve_domains, test_manifold_tangent_rollover, test_side_edge_convexity, FilStatus,
    FilletGeomRef, FilletSolver, FilletSolverOps, FilletSphereSolveEnfo,
    FilletSurfaceGeneratorType, RailRailIntersectEnfo, RailUvCurveIntersectEnfo,
};
use crate::curve::curve::p1::GlobalPointSolveOptions;
use crate::curve::curve::Curve;
use crate::curve::line::Line;
use crate::error::{ignore_status, Result, SolidError};
use crate::geom::axis2_placement::Axis2Placement;
use crate::geom::extent1d::{Extent1d, ExtentNd};
use crate::geom::extent2d::Extent2d;
use crate::geom::geom_utility::p2::ContinuityType;
use crate::geom::global_solver::{SolutionRequested, SolverOperation};
use crate::geom::local_solve_nd::{BoundaryHandler, EvalNFunctions, LocalSolveNd};
use crate::geom::matrix::Matrix;
use crate::geom::solution_array::{ObjectKind, SolutionObject};
use crate::geom::tol::ANG_TOL_RAD;
use crate::geom::vector2d::Point2;
use crate::intersect::curve_class::p1::CurveClassification;
use crate::intersect::point_class::ClassObject;
use crate::intersect::point_class::PointClassType;
use crate::intersect::surface_intersector::p1::{
    IntersectionPointType, SrfIndex, TsectPnt, MAX_USER_DOUBLES,
};
use crate::kernel::Kernel;
use crate::math::{rad_to_deg, Point3, Vec3, EFF_ZERO, EFF_ZERO_SQ, EFF_ZERO_SQRT, PI};
use crate::nurbs::b_spline_curve::p1::NurbCircleParam;
use crate::nurbs::b_spline_curve::BSplineCurve;
use crate::nurbs::b_spline_curve_n_lib::TrimMode;
use crate::nurbs::b_spline_surface::{BSplineSurface, SurfParamType};
use crate::nurbs::b_spline_surface_n_lib::p2::OrientType;
use crate::surface::plane::Plane;
use crate::surface::surf_of_revolution::SurfOfRevolution;
use crate::surface::surface::Surface;
use crate::topology::edgeuse::p1 as eu;
use crate::topology::edgeuse::p3::get_or_create_uv_trim_curve;
use crate::topology::face::p2 as face_p2;
use crate::topology::{edge, vertexuse, Brep, EdgeId, EdgeuseId, FaceId, VertexId};
use crate::topology::{vertex, VertexuseId};

const OP: &str = "blend::fillet_solver::p2";

/// `dAnalyticTol` of the `Calc*FilletGeom` methods.
const ANALYTIC_TOL: f64 = 1.0e-8;

/// `SM_FE_CROSS_SECTION` (`SmFilletEdgeType`, `SmFilletGeom.h`).
const FILLET_EDGE_CROSS_SECTION: i64 = 1;
/// `SM_FE_ON_EXTEND_EDGE`.
const FILLET_EDGE_ON_EXTEND_EDGE: i64 = 7;
/// `SM_FE_RAIL`.
const FILLET_EDGE_RAIL: i64 = 9;
/// `SM_FE_SPLIT_FACE`.
const FILLET_EDGE_SPLIT_FACE: i64 = 13;
/// `SM_FG_CLIFF_SIDE_PATCH` (`SmFilletGeomType`, `SmFilletGeom.h`).
const FILLET_GEOM_CLIFF_SIDE_PATCH: i64 = 4;

/// The curve and surface handles of `CreateCliffRollOverSideFG`'s
/// `GlobalCurveIntersect` (the cross-section curve and the side surface).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CliffHandle {
    /// `*pCurve`.
    Curve,
    /// `pSideSurface`.
    Surface,
}

impl SolutionObject for CliffHandle {
    fn kind(&self) -> ObjectKind {
        match self {
            Self::Curve => ObjectKind::Curve,
            Self::Surface => ObjectKind::Surface,
        }
    }
}

fn invalid(reason: &'static str) -> SolidError {
    SolidError::invalid_input(OP, reason)
}

fn sm_err() -> SolidError {
    SolidError::sm_err(OP)
}

/// An `SmFilletCorner*` (unit `blend.fillet_corner_p1`, not modelled).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FilletCornerRef;

/// The analytic fillet geometry a `Calc*FilletGeom` method builds before storing it into
/// the fillet geometry: the fillet surface (`m_pFilletSurface`), the rail curves
/// (`GetRail(i)->SetCurve`) with the interval each rail edge gets (`SetInterval`, `None`
/// where upstream sets none) and the centre line (`m_pCenterLineCurve`, `None` where
/// upstream sets none).
#[derive(Debug, Clone)]
pub struct FilletGeomParts {
    /// The fillet surface.
    pub surface: Surface,
    /// The rail curves.
    pub rails: [Curve; 2],
    /// The rail edge intervals.
    pub rail_intervals: [Option<Extent1d>; 2],
    /// The centre line.
    pub center_line: Option<Curve>,
}

/// The reference arguments of [`find_rollover_side_edge`] (`rlCurrIndx` and
/// `rpRolloverCornerEU` are read too; each path writes only some of them).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct RolloverSideEdge {
    /// `rlCurrIndx`.
    pub curr_index: usize,
    /// `rpSideEdge`.
    pub side_edge: Option<EdgeId>,
    /// `rbIsTangentialSideEdge`.
    pub is_tangential_side_edge: bool,
    /// `rbTrimSideFaces`.
    pub trim_side_faces: bool,
    /// `rbOppositeNormals`.
    pub opposite_normals: bool,
    /// `rpRolloverCornerEU`.
    pub rollover_corner_eu: Option<EdgeuseId>,
    /// `rpRolloverCornerV`.
    pub rollover_corner_v: Option<VertexId>,
}

// ---- function-object adapters ---------------------------------------------------------

struct SphereEval<'a>(FilletSphereSolveEnfo<'a>);

impl EvalNFunctions for SphereEval<'_> {
    fn evaluate(&mut self, x: &[f64], f: &mut [f64], j: Option<&mut Matrix>) -> Result<bool> {
        self.0.evaluate(x, f, j)
    }
}

struct RailUvEval<'a>(RailUvCurveIntersectEnfo<'a>);

impl EvalNFunctions for RailUvEval<'_> {
    fn evaluate(&mut self, x: &[f64], f: &mut [f64], j: Option<&mut Matrix>) -> Result<bool> {
        self.0.evaluate(x, f, j)
    }
}

struct RailRailEval<'a>(RailRailIntersectEnfo<'a>);

impl EvalNFunctions for RailRailEval<'_> {
    fn evaluate(&mut self, x: &[f64], f: &mut [f64], j: Option<&mut Matrix>) -> Result<bool> {
        self.0.evaluate(x, f, j)
    }
}

// ---- small helpers --------------------------------------------------------------------

/// `sIntervals[i] = e`; a write past the Jacobian size is dropped (module notes).
fn set_interval(intervals: &mut ExtentNd, i: usize, e: Extent1d) {
    if let Some(slot) = intervals.get_mut(i) {
        *slot = e;
    }
}

fn srf_of_rail(rail: usize) -> Result<SrfIndex> {
    match rail {
        0 => Ok(SrfIndex::First),
        1 => Ok(SrfIndex::Second),
        _ => Err(invalid("rail index out of range")),
    }
}

fn sol_at(x: &[f64], i: usize) -> Result<f64> {
    x.get(i)
        .copied()
        .ok_or_else(|| invalid("solution index out of range"))
}

/// A copy of the solver's offset surface of `rail` (`*GetSurface(rail)`), to hand to
/// `LoadInitialValues` alongside the solver; written back with [`put_surface`].
fn take_surface(s: &FilletSolver, rail: usize) -> Result<OffsetSurface> {
    s.surface(rail)?
        .map(|slot| slot.surface.clone())
        .ok_or_else(|| invalid("no offset surface"))
}

fn put_surface(s: &mut FilletSolver, rail: usize, surface: OffsetSurface) -> Result<()> {
    if let Some(slot) = s.surface_mut(rail)? {
        slot.surface = surface;
    }
    Ok(())
}

/// `GetSurface(rail)->EvaluatePoint(uv)` (keeping an extension it makes).
fn eval_offset(s: &mut FilletSolver, rail: usize, uv: Point2) -> Result<Point3> {
    s.surface_mut(rail)?
        .ok_or_else(|| invalid("no offset surface"))?
        .surface
        .surface
        .evaluate_point_extending(uv)
}

fn has_surface(s: &FilletSolver, rail: usize) -> Result<bool> {
    Ok(s.surface(rail)?.is_some())
}

/// `pCurve->Evaluate(t, 1, TRUE, sPV)`: point and first derivative.
fn curve_pv(c: &Curve, t: f64) -> Result<(Point3, Vec3)> {
    let pv = c.evaluate(t, 1, true, false)?;
    let p = pv
        .first()
        .copied()
        .ok_or_else(|| invalid("short curve evaluation"))?;
    let d = pv
        .get(1)
        .copied()
        .ok_or_else(|| invalid("short curve evaluation"))?;
    Ok((p, d))
}

/// The edge of the solver's first edge use (`GetEdgeuse(0)->GetEdge()`).
fn filleted_edge(s: &FilletSolver, b: &Brep) -> Result<EdgeId> {
    let eu0 = s.edgeuse(0)?.ok_or_else(|| invalid("no edge use"))?;
    eu::edge(b, eu0)?.ok_or_else(|| invalid("edge use without edge"))
}

/// The filleted edge's curve after `sm_CheckEdgeCurveDomains` (the trimmed copy, or the
/// edge's own curve).
fn checked_edge_curve(k: &Kernel, s: &FilletSolver, b: &Brep, e: EdgeId) -> Result<Curve> {
    if let Some(copy) = check_edge_curve_domains(k, b, e, s.extend_before, s.extend_after)? {
        return Ok(copy);
    }
    // NER(pCurve)
    edge::curve(b, e)?
        .cloned()
        .ok_or_else(|| SolidError::null_pointer(OP))
}

fn is_circular_fsg(s: &FilletSolver) -> Result<bool> {
    let fsg = s
        .fsg
        .as_ref()
        .ok_or_else(|| invalid("no fillet surface generator"))?;
    Ok(fsg.fillet_surface_generator_type() == FilletSurfaceGeneratorType::Circular)
}

/// `rTsectPnts[0]` and `rTsectPnts[GetSize() - 1]`.
fn first_last(tsects: &[TsectPnt]) -> Result<(&TsectPnt, &TsectPnt)> {
    let first = tsects.first().ok_or_else(|| invalid("no tsect points"))?;
    let last = tsects.last().ok_or_else(|| invalid("no tsect points"))?;
    Ok((first, last))
}

/// `SetupOffsetValues(rTsectPnt, &dRadius)`, returning the radius.
fn offset_radius(ops: &mut dyn FilletSolverOps, t: &TsectPnt) -> Result<f64> {
    let mut r = 0.0;
    ops.setup_offset_values(t, Some(&mut r))?;
    Ok(r)
}

/// Stores the parts into the fillet geometry (`pFilletGeom` or `m_vFilletGeoms[0]`).
fn store_parts(
    s: &FilletSolver,
    fillet_geom: Option<FilletGeomRef>,
    parts: FilletGeomParts,
) -> Result<()> {
    let fg = match fillet_geom {
        Some(fg) => fg,
        None => s.first_fillet_geom()?,
    };
    let _ = (fg, parts);
    // UNPORTED(blend.fillet_geom_p1): SmFilletGeom::m_pFilletSurface (member write)
    Err(SolidError::unported(
        "SmFilletGeom::m_pFilletSurface",
        "blend.fillet_geom_p1",
    ))
}

/// `GetOrCreateUVTrimCurve` hands back an `SmBSplineCurve*`: only the kinds deriving
/// from it (B-splines, lines, ellipses, hyperbolas, parabolas) can be a UV trim curve.
fn check_uv_bspline_kind(c: &Curve) -> Result<()> {
    use crate::curve::curve::CurveKind;
    match c.kind() {
        CurveKind::BSpline(_)
        | CurveKind::Line(_)
        | CurveKind::Ellipse(_)
        | CurveKind::Hyperbola(_)
        | CurveKind::Parabola(_) => Ok(()),
        _ => Err(invalid("UV trim curve is not an SmBSplineCurve")),
    }
}

/// `SmLine::CreateLineSegment(ctx, 3, a, b, pLine)` then `EditParameterization(ivl)`.
fn line_segment(a: Point3, b: Point3, ivl: &Extent1d) -> Result<Line> {
    let mut l = Line::create_line_segment(3, a, b, None)?;
    l.edit_parameterization(ivl)?;
    Ok(l)
}

impl FilletSolver {
    /// `GetFilletCorner(lIndex)`: the fillet corner at the start (0) or end (1) vertex
    /// of the filleted edge; `None` for another index.
    ///
    /// # Errors
    ///
    /// Those of [`FilletSolver::vertex`]; [`SolidError::InvalidInput`] without an
    /// executive (upstream dereferences `NULL`); [`SolidError::Unported`] for the corner
    /// lookup (`SmFilletExecutive::GetFilletCornerOfVertex`).
    pub fn fillet_corner(&self, b: &Brep, index: usize) -> Result<Option<FilletCornerRef>> {
        if index > 1 {
            return Ok(None);
        }
        let vtx = self.vertex(b, index)?;
        if self.executive.is_none() {
            return Err(invalid("no fillet executive"));
        }
        let _ = vtx;
        // UNPORTED(blend.fillet_executive_p1): SmFilletExecutive::GetFilletCornerOfVertex
        Err(SolidError::unported(
            "SmFilletExecutive::GetFilletCornerOfVertex",
            "blend.fillet_executive_p1",
        ))
    }

    /// `CalcCylinderFilletGeom(rTsectPnts, pFilletGeom)`: builds the cylinder of a
    /// constant-radius circular fillet along a straight edge
    /// ([`FilletSolver::cylinder_fillet_geom_parts`]) and stores it into `fillet_geom`
    /// (`m_vFilletGeoms[0]` for `None`).
    ///
    /// # Errors
    ///
    /// Those of [`FilletSolver::cylinder_fillet_geom_parts`]; then
    /// [`SolidError::Unported`] for the store (unit `blend.fillet_geom_p1`).
    pub fn calc_cylinder_fillet_geom(
        ops: &mut dyn FilletSolverOps,
        k: &Kernel,
        b: &Brep,
        tsects: &[TsectPnt],
        fillet_geom: Option<FilletGeomRef>,
    ) -> Result<()> {
        let parts = Self::cylinder_fillet_geom_parts(ops, k, b, tsects)?;
        store_parts(ops.solver(), fillet_geom, parts)
    }

    /// The geometry of [`FilletSolver::calc_cylinder_fillet_geom`]: the cylinder
    /// (`CreateConePatch` with equal radii, UV swapped, `CopyAndAddAnalytics`), the two
    /// rail lines and the centre line, each parameterised on the cylinder's u range.
    /// The ends are extended by the extension factor where the solver extends.
    ///
    /// # Errors
    ///
    /// `SM_ERR` when the generator is not circular, the edge curve is not a line, the
    /// radius changes, or the arc angle is degenerate or near `pi`; those of
    /// `SetupOffsetValues` (`SM_ERR` for the base solver), the evaluations and the
    /// constructors; [`SolidError::InvalidInput`] where upstream dereferences `NULL`
    /// or for no tsect points.
    pub fn cylinder_fillet_geom_parts(
        ops: &mut dyn FilletSolverOps,
        k: &Kernel,
        b: &Brep,
        tsects: &[TsectPnt],
    ) -> Result<FilletGeomParts> {
        let s = ops.solver();
        let e = filleted_edge(s, b)?;

        // Check for Edge's curve domain being bigger than the Edge domain.
        let curve = checked_edge_curve(k, s, b, e)?;

        // no work - not a circular cross-section on a straight line
        if !is_circular_fsg(s)? || curve.is_line(3, ANALYTIC_TOL)?.is_none() {
            return Err(sm_err());
        }

        // get start and end Offset Surface UV points
        let (t_first, t_last) = first_last(tsects)?;
        let start_uv0 = t_first.uv_pos(SrfIndex::First);
        let start_uv1 = t_first.uv_pos(SrfIndex::Second);
        let end_uv0 = t_last.uv_pos(SrfIndex::First);
        let end_uv1 = t_last.uv_pos(SrfIndex::Second);

        // set OffsetSurface offset = radius for begin point and get sCenter0
        let mut radius0 = offset_radius(ops, t_first)?;
        let mut center0 = eval_offset(ops.solver_mut(), 0, start_uv0)?;

        // set OffsetSurface offset = radius for end point and get sCenter1
        let mut radius1 = offset_radius(ops, t_last)?;
        let mut center1 = eval_offset(ops.solver_mut(), 0, end_uv0)?;

        // check state - cylinder case = constant radius fillets
        if (radius1 - radius0).abs() > EFF_ZERO {
            // Non-cylinder case
            return Err(sm_err());
        }

        // Evaluate points on base surfaces.  Use a zero-offset offset surface
        // instead of the base surface, in order to handle extended surfaces.
        let s = ops.solver_mut();
        s.set_surfaces_to_zero_offset()?;
        let mut start_pnt0 = eval_offset(s, 0, start_uv0)?;
        let mut start_pnt1 = eval_offset(s, 1, start_uv1)?;
        let mut end_pnt0 = eval_offset(s, 0, end_uv0)?;
        let mut end_pnt1 = eval_offset(s, 1, end_uv1)?;
        s.reload_surface_offsets()?;

        // get extension vectors along each rail curve and the center curve
        let f = s.surface_extension_factor;
        let ext0 = f * (end_pnt0 - start_pnt0);
        let ext1 = f * (end_pnt1 - start_pnt1);
        let center_ext = f * (center1 - center0);

        // when extending the start of the fillet
        if s.extend_before {
            start_pnt0 -= ext0;
            start_pnt1 -= ext1;
            center0 -= center_ext;
            radius0 = center0.distance(start_pnt0);
        }

        // when extending the end of the fillet
        if s.extend_after {
            end_pnt0 += ext0;
            end_pnt1 += ext1;
            center1 += center_ext;
            radius1 = center1.distance(end_pnt0);
        }
        let _ = radius1;

        // get the cylinder parameters
        let mut x_vec = start_pnt0 - center0;
        let mut y_vec = start_pnt1 - center0;
        let axis_of_revolution = center1 - center0;
        let height = axis_of_revolution.length();

        // Get angle of revolution
        let mut angle = axis_of_revolution.ccw_angle_between(x_vec, y_vec)?;

        // check state - avoid degenerate zero arc cylinders
        if angle.abs() < EFF_ZERO_SQRT {
            return Err(sm_err()); // Degenerate case here
        }

        // check state - avoid fillet arcs equal or greater than Pi
        if angle.abs() > PI - EFF_ZERO_SQRT {
            return Err(sm_err()); // Nearly 180 degree case here
        }

        // swap center to start vectors to force rotation angle to be positive
        if angle < 0.0 {
            std::mem::swap(&mut x_vec, &mut y_vec);
            angle = axis_of_revolution.ccw_angle_between(x_vec, y_vec)?;
        }

        // check state - avoid degenerate zero arc cylinders
        if angle < EFF_ZERO {
            return Err(sm_err());
        }

        // convert radians to degrees
        let angle = rad_to_deg(angle);

        // orthogonalize and normalize sXVec and sYVec
        let mut placement = Axis2Placement::default();
        y_vec = axis_of_revolution.cross(x_vec);
        let x_vec = x_vec.unitize()?;
        let y_vec = y_vec.unitize()?;

        // set up a coordinate system for cylinder center line
        placement.set_canonical(center0, x_vec, y_vec)?;

        // Create a temporary cylindrical cone patch
        let mut cone_patch = BSplineSurface::create_cone_patch(
            &placement,
            radius0,
            radius0,
            0.0,
            angle,
            height,
            NurbCircleParam::Quadratic,
        )?;

        // CreateConePatch() makes a STEP-parameterized cone,
        // we go the other way: u along the edge, v around the fillet surface.
        ignore_status(cone_patch.swap_uv())?;

        // Copy pConePatch, when possible as an analytic surface
        let cylinder = cone_patch.copy_and_add_analytics(k)?;

        // get cylinder u parameter range
        let uv_domain = cylinder.natural_uv_domain();
        let rail_crv_ivl = Extent1d::from_min_max(uv_domain.min.x, uv_domain.max.x);

        // Calc Rail0, Rail1 - set their extent = cylinder u parameter range
        let rail0 = line_segment(start_pnt0, end_pnt0, &rail_crv_ivl)?;
        let rail1 = line_segment(start_pnt1, end_pnt1, &rail_crv_ivl)?;

        // Calc and store centerline
        let center_line = line_segment(center0, center1, &rail_crv_ivl)?;

        Ok(FilletGeomParts {
            surface: cylinder,
            rails: [Curve::from(rail0), Curve::from(rail1)],
            rail_intervals: [None, None],
            center_line: Some(Curve::from(center_line)),
        })
    }

    /// `CalcConeFilletGeom(rTsectPnts, pFilletGeom)`: builds the cone of a circular
    /// fillet along a straight edge ([`FilletSolver::cone_fillet_geom_parts`]) and
    /// stores it into `fillet_geom` (`m_vFilletGeoms[0]` for `None`).
    ///
    /// # Errors
    ///
    /// Those of [`FilletSolver::cone_fillet_geom_parts`]; then
    /// [`SolidError::Unported`] for the store (unit `blend.fillet_geom_p1`).
    pub fn calc_cone_fillet_geom(
        ops: &mut dyn FilletSolverOps,
        k: &Kernel,
        b: &Brep,
        tsects: &[TsectPnt],
        fillet_geom: Option<FilletGeomRef>,
    ) -> Result<()> {
        let parts = Self::cone_fillet_geom_parts(ops, k, b, tsects)?;
        store_parts(ops.solver(), fillet_geom, parts)
    }

    /// The geometry of [`FilletSolver::calc_cone_fillet_geom`]: the first rail's line
    /// revolved about the centre line (`SmSurfOfRevolution::CreateCanonical`, its sweep
    /// limited to the fillet angle, UV swapped), its generator copy and its v-max iso
    /// curve as rails, and the centre line extended at both ends.
    ///
    /// # Errors
    ///
    /// `SM_ERR` when the generator is not circular, the edge curve is not a line, or the
    /// arc angle is degenerate or near `pi`; those of `SetupOffsetValues`, the
    /// evaluations and the constructors; `SM_ERR_NULL_POINTER` for a missing generator;
    /// [`SolidError::InvalidInput`] where upstream dereferences `NULL`.
    pub fn cone_fillet_geom_parts(
        ops: &mut dyn FilletSolverOps,
        k: &Kernel,
        b: &Brep,
        tsects: &[TsectPnt],
    ) -> Result<FilletGeomParts> {
        let s = ops.solver();
        let e = filleted_edge(s, b)?;

        // Check for Edge's curve domain being bigger than the Edge domain.
        let curve = checked_edge_curve(k, s, b, e)?;

        // no work - not a circular cross-section on a line
        if !is_circular_fsg(s)? || curve.is_line(3, ANALYTIC_TOL)?.is_none() {
            return Err(sm_err());
        }

        // Get some parameters for surf-of-revo creation
        let (t_first, t_last) = first_last(tsects)?;
        let start_uv0 = t_first.uv_pos(SrfIndex::First);
        let start_uv1 = t_first.uv_pos(SrfIndex::Second);
        let end_uv0 = t_last.uv_pos(SrfIndex::First);
        let end_uv1 = t_last.uv_pos(SrfIndex::Second);
        let _radius0 = offset_radius(ops, t_first)?;
        let mut center0 = eval_offset(ops.solver_mut(), 0, start_uv0)?;
        let _radius1 = offset_radius(ops, t_last)?;
        let mut center1 = eval_offset(ops.solver_mut(), 0, end_uv0)?;

        // Evaluate points on base surfaces.  Use a zero-offset offset surface
        // instead of the base surface, in order to handle extended surfaces.
        let s = ops.solver_mut();
        s.set_surfaces_to_zero_offset()?;
        let mut start_pnt0 = eval_offset(s, 0, start_uv0)?;
        let start_pnt1 = eval_offset(s, 1, start_uv1)?;
        let mut end_pnt0 = eval_offset(s, 0, end_uv0)?;
        let _end_pnt1 = eval_offset(s, 1, end_uv1)?;
        s.reload_surface_offsets()?;

        // Get origin & axis of revolution
        let v0 = start_pnt0 - center0;
        let v1 = start_pnt1 - center0;
        let mut axis_of_revolution = (center1 - center0).unitize()?;
        let angle = axis_of_revolution.ccw_angle_between(v0, v1)?;
        if angle.abs() < EFF_ZERO_SQRT {
            return Err(sm_err()); // Degenerate case here
        }
        if angle.abs() > PI - EFF_ZERO_SQRT {
            return Err(sm_err()); // Nearly 180 degree case here
        }
        if angle < 0.0 {
            axis_of_revolution = -axis_of_revolution;
        }
        let angle = axis_of_revolution.ccw_angle_between(v0, v1)?;
        if angle < EFF_ZERO {
            return Err(sm_err());
        }
        let angle = rad_to_deg(angle);

        let v0 = v0.unitize()?;
        let mut placement = Axis2Placement::default();
        let x_vec = v0;
        let y_vec = axis_of_revolution.cross(x_vec);
        let x_vec = x_vec.unitize()?;
        let y_vec = y_vec.unitize()?;
        placement.set_canonical(center0, x_vec, y_vec)?;
        let _ = placement;

        // Calc generator curve (Rail0)
        let f = s.surface_extension_factor;
        let dir = end_pnt0 - start_pnt0;
        let ext = dir * f;
        if s.extend_before {
            start_pnt0 -= ext;
        }
        if s.extend_after {
            end_pnt0 += ext;
        }
        let gen_curve = BSplineCurve::create_line_segment(3, start_pnt0, end_pnt0)?;

        // Now, create surface of revolution
        let mut cone = SurfOfRevolution::create_canonical(
            k,
            Curve::from(gen_curve),
            center0,
            axis_of_revolution,
            false,
        )?;
        let mut step_domain = cone.step_uv_domain();
        let mut uv_max = step_domain.max;
        uv_max.x = angle;
        step_domain = Extent2d::from_points(step_domain.min, uv_max);
        cone.adjust_step_uv_domain(k, &step_domain)?;
        let mut cone = Surface::from(cone);
        ignore_status(cone.swap_uv(k))?;

        let uv_domain = cone.natural_uv_domain();
        let ivl = Extent1d::from_min_max(uv_domain.min.x, uv_domain.max.x);

        // Calc rail0
        let gen = match cone.kind() {
            crate::surface::surface::SurfaceKind::SurfOfRevolution(r) => r.gen_curve(),
            _ => None,
        }
        .ok_or_else(|| invalid("surface of revolution without generator"))?;
        let mut rail_crv0 = gen.copy(k)?;
        rail_crv0.edit_parameterization(k, &ivl)?;
        let crv_ivl0 = rail_crv0.natural_interval();

        // Calc Rail1
        let rail_crv1 = cone
            .create_iso_parametric_curve(k, SurfParamType::V, uv_domain.max.y, 0.0, None, false)?
            .curve;
        let crv_ivl1 = rail_crv1.natural_interval();

        // Calc centerline
        let dir = center1 - center0;
        let ext = dir * f;
        center0 -= ext;
        center1 += ext;
        let mut center_line = BSplineCurve::create_line_segment(3, center0, center1)?;
        center_line.edit_parameterization(&ivl)?;

        Ok(FilletGeomParts {
            surface: cone,
            rails: [rail_crv0, rail_crv1],
            rail_intervals: [Some(crv_ivl0), Some(crv_ivl1)],
            center_line: Some(Curve::from(center_line)),
        })
    }

    /// `CalcTorusFilletGeom(rTsectPnts, pFilletGeom)`: builds the torus of a
    /// constant-radius circular fillet along a circular edge
    /// ([`FilletSolver::torus_fillet_geom_parts`]) and stores it into `fillet_geom`
    /// (`m_vFilletGeoms[0]` for `None`).
    ///
    /// # Errors
    ///
    /// Those of [`FilletSolver::torus_fillet_geom_parts`]; then
    /// [`SolidError::Unported`] for the store (unit `blend.fillet_geom_p1`).
    pub fn calc_torus_fillet_geom(
        ops: &mut dyn FilletSolverOps,
        k: &Kernel,
        b: &Brep,
        tsects: &[TsectPnt],
        fillet_geom: Option<FilletGeomRef>,
    ) -> Result<()> {
        let parts = Self::torus_fillet_geom_parts(ops, k, b, tsects, fillet_geom)?;
        let fg = match fillet_geom {
            Some(fg) => Some(fg),
            None => Some(ops.solver().first_fillet_geom()?),
        };
        store_parts(ops.solver(), fg, parts)
    }

    /// The geometry of [`FilletSolver::calc_torus_fillet_geom`]: the cross-section arc
    /// through the rail points revolved about the edge circle's axis
    /// (`SmSurfOfRevolution::CreateCanonical`, the sweep extended where the solver
    /// extends, at most 360 degrees, the seam moved by `pi` for a full-circle torus of an
    /// open edge, rotated back by the extension angle when extending before), and its
    /// v-min and v-max iso curves as rails (no centre line).
    ///
    /// # Errors
    ///
    /// `SM_ERR` when the generator is not circular, the edge curve is not an arc, the
    /// radius is tiny or changes, or the nearest point on the edge is not unique;
    /// [`SolidError::Unported`] for the rail end checks of an open edge that does not
    /// extend after (`SmFilletGeom::GetRail`, unit `blend.fillet_geom_p1`); those of
    /// `SetupOffsetValues`, the evaluations and the constructors;
    /// [`SolidError::InvalidInput`] where upstream dereferences `NULL`.
    pub fn torus_fillet_geom_parts(
        ops: &mut dyn FilletSolverOps,
        k: &Kernel,
        b: &Brep,
        tsects: &[TsectPnt],
        fillet_geom: Option<FilletGeomRef>,
    ) -> Result<FilletGeomParts> {
        let s = ops.solver();
        let e = filleted_edge(s, b)?;

        // Check for Edge's curve domain being bigger than the Edge domain.
        let curve = checked_edge_curve(k, s, b, e)?;

        // no work - not a circular cross-section on a circular arc
        // Side effect: get circular arc RefFrame, Radius, StartAng, and EndAng.
        if !is_circular_fsg(s)? {
            return Err(sm_err());
        }
        let Some((ref_frame, _arc_radius, start_ang, _end_ang)) = curve.is_arc(11, ANALYTIC_TOL)?
        else {
            return Err(sm_err());
        };

        // select new fillet surface storage spot
        let fillet_geom = match fillet_geom {
            Some(fg) => fg,
            None => s.first_fillet_geom()?,
        };

        // Get rotation point and axis for surf-of-revolution creation
        let mut axis_of_revolution = ref_frame.z_axis();
        let mut origin = ref_frame.origin();

        // Locate the generating curve(arc)
        let (t_first, t_last) = first_last(tsects)?;
        let start_uv0 = t_first.uv_pos(SrfIndex::First);
        let start_uv1 = t_first.uv_pos(SrfIndex::Second);

        // get centerPoint = offsetSurface0(uvPnt0)
        let arc_center = eval_offset(ops.solver_mut(), 0, start_uv0)?;

        // get cross-section radius size at begin and end fillet points
        let radius = offset_radius(ops, t_first)?;
        let radius1 = offset_radius(ops, t_last)?;
        if radius < EFF_ZERO || (radius1 - radius).abs() > EFF_ZERO {
            return Err(sm_err());
        }

        // Evaluate points on base surfaces.  Use a zero-offset offset surface
        // instead of the base surface, in order to handle extended surfaces.
        let s = ops.solver_mut();
        s.set_surfaces_to_zero_offset()?;
        let p0 = eval_offset(s, 0, start_uv0)?; // sP0 = rail1 point
        let p1 = eval_offset(s, 1, start_uv1)?; // sP1 = rail2 point
        s.reload_surface_offsets()?;

        // allocate the arc from the points [centerPoint, sP0, sP1]
        let mut arc =
            BSplineCurve::create_arc_from_points(arc_center, p0, p1, NurbCircleParam::Quadratic)?;

        // Determine if we need to reverse the axis of revolution
        let mut vec = p0 - origin;
        if vec.length() < ANALYTIC_TOL {
            vec = p1 - origin;
        }
        let dir = axis_of_revolution.cross(vec); // Start direction of rail curves

        // get circular to-be-filleted edge->curve point (pos,tang) nearest to arcCenter Point
        let ivl = edge::interval(b, e)?;
        let sols = curve.global_point_solve(
            k,
            &ivl,
            SolverOperation::Minimize,
            arc_center,
            ANALYTIC_TOL,
            SolutionRequested::Single,
            &GlobalPointSolveOptions::<crate::topology::ObjectHandle>::default(),
        )?;
        if sols.len() != 1 {
            return Err(sm_err());
        }
        let t_near = sols.get(0).ok_or_else(sm_err)?.start.param(0)?;
        let (_, near_tan) = curve_pv(&curve, t_near)?;

        // when to-be-filleted edge tangent nearest sweep curve opposes initial sweep direction
        if near_tan.dot(dir) < 0.0 {
            axis_of_revolution = -axis_of_revolution;
        }

        // Recalculate the origin of torus (see upstream notes).
        let tmp_vec = arc_center - origin;
        origin += axis_of_revolution * tmp_vec.dot(axis_of_revolution);

        // Evaluate points at the ends of the rail curves embedded within the base surfaces.
        let end_uv0 = t_last.uv_pos(SrfIndex::First);
        let end_uv1 = t_last.uv_pos(SrfIndex::Second);
        let s = ops.solver_mut();
        s.set_surfaces_to_zero_offset()?;
        let end_p0 = eval_offset(s, 0, end_uv0)?; // sEndP0 = rail1 point
        let end_p1 = eval_offset(s, 1, end_uv1)?; // sEndP1 = rail2 point
        s.reload_surface_offsets()?;

        let mut end_vec = end_p0 - origin;

        // If sEndVec is zero or parallel to sAxisOfRev, use sEndP1.
        let temp_cross = axis_of_revolution.cross(end_vec);
        if temp_cross.length() < ANALYTIC_TOL {
            end_vec = end_p1 - origin;
        }

        // CCWAngleBetween status ignored upstream: dSpanAng stays 0 on a failure
        let mut span_ang = match axis_of_revolution.ccw_angle_between(vec, end_vec) {
            Ok(a) => a,
            Err(e @ SolidError::Unported { .. }) => return Err(e),
            Err(_) => 0.0,
        };
        if span_ang < ANG_TOL_RAD {
            span_ang += 2.0 * PI;
        }

        let end_ang = rad_to_deg(span_ang);
        let mut angle_span = end_ang - start_ang;

        // Determine if we need to extend the span even when 'm_bExtendBefore' is FALSE
        if !curve.is_closed(&ivl, 0.0)? && !s.extend_after {
            // for both filletSector->edgeuses
            for i in 0..2 {
                let eu_i = s.edgeuse(i)?.ok_or_else(|| invalid("no edge use"))?;
                let vu = eu::vertexuse(b, eu_i)?
                    .ok_or_else(|| invalid("edge use without vertex use"))?;
                let _v = vertexuse::vertex(b, vu)?
                    .ok_or_else(|| invalid("vertex use without vertex"))?;
                let orient = eu::orientation(b, eu_i)?;
                if (orient == OrientType::Same && !s.extend_before)
                    || orient == OrientType::Opposite
                {
                    let _ = fillet_geom;
                    // UNPORTED(blend.fillet_geom_p1): SmFilletGeom::GetRail
                    return Err(SolidError::unported(
                        "SmFilletGeom::GetRail",
                        "blend.fillet_geom_p1",
                    ));
                }
            }
        }

        // extend angle span when requested - limit extent to 360 degrees
        let mut extend_flag = false;
        let extension_angle = (360.0 - angle_span) / 4.0;
        if angle_span < 350.0 {
            if s.extend_after {
                angle_span += extension_angle;
            }
            if s.extend_before {
                angle_span += extension_angle;
                extend_flag = true;
            }
        }
        if angle_span > 360.0 {
            angle_span = 360.0;
        }

        // Move the seam of a full-circle torus of an open edge out of the way [B466].
        let eu0 = s.edgeuse(0)?.ok_or_else(|| invalid("no edge use"))?;
        let e0 = eu::edge(b, eu0)?.ok_or_else(|| invalid("edge use without edge"))?;
        let e0_curve = edge::curve(b, e0)?.ok_or_else(|| invalid("edge without curve"))?;
        let move_seam = if e0_curve.is_closed(&edge::interval(b, e0)?, 0.0)? {
            false
        } else {
            angle_span > 360.0 - EFF_ZERO
        };

        if move_seam {
            // Revolve pArc around to the other side of the torus.  [B466]
            let mut cs = Axis2Placement::default();
            cs.init();
            cs.translate(-ref_frame.origin());
            cs.rotate_about_axis(PI, ref_frame.z_axis());
            cs.translate(ref_frame.origin());
            ignore_status(arc.transform(&cs, None))?;
        }

        // create surface of revolution
        let mut torus = SurfOfRevolution::create_canonical(
            k,
            Curve::from(arc),
            origin,
            axis_of_revolution,
            false,
        )?;

        // extend the torus uvDomain as needed
        let mut step_domain = torus.step_uv_domain();
        let mut uv_max = step_domain.max;
        uv_max.x = angle_span;
        step_domain = Extent2d::from_points(step_domain.min, uv_max);
        torus.adjust_step_uv_domain(k, &step_domain)?;

        // Z-Rotate the torus by dAngleExtendBefore
        if extend_flag {
            let mut trans = Axis2Placement::default();
            trans.translate(-origin);
            trans.rotate_about_axis((-extension_angle) * (PI / 180.0), axis_of_revolution);
            trans.translate(origin);
            torus.transform(k, &trans, None)?;
        }
        let torus = Surface::from(torus);

        // Calc rails: the min and max SM_SP_V isoParameter curves
        let uv_domain = torus.natural_uv_domain();
        let rail_crv0 = torus
            .create_iso_parametric_curve(k, SurfParamType::V, uv_domain.min.y, 0.0, None, false)?
            .curve;
        let crv_ivl0 = rail_crv0.natural_interval();
        let rail_crv1 = torus
            .create_iso_parametric_curve(k, SurfParamType::V, uv_domain.max.y, 0.0, None, false)?
            .curve;
        let crv_ivl1 = rail_crv1.natural_interval();

        Ok(FilletGeomParts {
            surface: torus,
            rails: [rail_crv0, rail_crv1],
            rail_intervals: [Some(crv_ivl0), Some(crv_ivl1)],
            center_line: None,
        })
    }

    /// `PointOnPlaneSolve(crPlaneOrig, crPlaneNormal, crUVDomain1, crUVDomain2, rUV1,
    /// rUV2, rbFoundSolution, rTsectPnt)`: the fillet point on the plane through
    /// `plane_orig` normal to `plane_normal`, solved with `ops`'s equations
    /// ([`FilletSphereSolveEnfo`]; a second try lets the plane adjust). Returns
    /// `rbFoundSolution`; when found, `tsect` gets both UV points and the further
    /// solution values as user doubles.
    ///
    /// # Errors
    ///
    /// Those of `GetJacobianSize`, `LoadInitialValues` (`SM_ERR` for the base solver)
    /// and the solves; [`SolidError::InvalidInput`] without a first offset surface or
    /// for more than [`MAX_USER_DOUBLES`] further values (upstream writes past
    /// `m_adUserDoubles`).
    pub fn point_on_plane_solve(
        ops: &mut dyn FilletSolverOps,
        k: &Kernel,
        brep: Option<&Brep>,
        plane_orig: Point3,
        plane_normal: Vec3,
        uv_domain1: &Extent2d,
        uv_domain2: &Extent2d,
        uv1: Point2,
        uv2: Point2,
        tsect: &mut TsectPnt,
    ) -> Result<bool> {
        // Now Set up Solver
        let jacobian_size = ops.jacobian_size()?;
        let mut intervals = ExtentNd::new(jacobian_size);
        set_interval(
            &mut intervals,
            0,
            Extent1d::from_min_max(uv_domain1.min.x, uv_domain1.max.x),
        );
        set_interval(
            &mut intervals,
            1,
            Extent1d::from_min_max(uv_domain1.min.y, uv_domain1.max.y),
        );
        set_interval(
            &mut intervals,
            2,
            Extent1d::from_min_max(uv_domain2.min.x, uv_domain2.max.x),
        );
        set_interval(
            &mut intervals,
            3,
            Extent1d::from_min_max(uv_domain2.min.y, uv_domain2.max.y),
        );

        let mut periodicities = vec![false, false, false, false];
        let mut guess_t = vec![uv1.x, uv1.y, uv2.x, uv2.y];

        let do_surf2_calcs = false; // Don't do surface2 calculations
        let surf1_index = 0; // UV values of Surface1
        let surf2_index = 2; // UV values of Surface2
        let rail_index = 0; // rail index for pSurface1

        // given current inputs - place a good guess for upcoming solve in sGuessT
        let mut surface1 = take_surface(ops.solver(), 0)?;
        let loaded = ops.load_initial_values(
            k,
            brep,
            do_surf2_calcs,
            rail_index,
            &mut surface1,
            surf1_index,
            &mut intervals,
            &mut periodicities,
            &mut guess_t,
        );
        put_surface(ops.solver_mut(), 0, surface1)?;
        // lSurf2Index is an in/out argument: a solver may move it
        let surf2_index = loaded?.map_or(surf2_index, |l| l.surface2_index);

        // get plane equation for given PlanePt and PlaneNormal
        let plane_d = -plane_orig.dot(plane_normal);
        let approx_tol = ops.solver().approx_tol_3d;

        // allocate NewtonRaphson evaluator - contains equation definitions being solved
        let mut eval = SphereEval(FilletSphereSolveEnfo::new(plane_normal, plane_d, ops));
        let mut ls = LocalSolveNd::new(&mut eval, Some(&intervals), Some(&periodicities));

        // set solve parameters
        ls.set_desired_accuracy(EFF_ZERO_SQ);
        ls.set_boundary_handler(BoundaryHandler::TotalBoundaryHits, 8);
        ls.set_maximum_iterations(30);

        // Solve equations with NewtonRaphson iteration
        let mut sol = ls.solve_it(&guess_t, approx_tol)?;

        // Try it again and allow the plane to adjust as convergance happens
        if !sol.found {
            ls.evaluator().0.set_plane_adjusts(true);
            sol = ls.solve_it(&guess_t, approx_tol)?;
        }

        // when a solution point was found
        if sol.found {
            let x = &sol.x;
            *tsect.uv_pos_mut(SrfIndex::First) =
                Point2::new(sol_at(x, surf1_index)?, sol_at(x, surf1_index + 1)?);
            *tsect.uv_pos_mut(SrfIndex::Second) =
                Point2::new(sol_at(x, surf2_index)?, sol_at(x, surf2_index + 1)?);
            for i in 4..x.len() {
                // upstream has no SM_MAX_USER_DOUBLES check here (UB-192)
                let slot = tsect
                    .user_doubles
                    .get_mut(i - 4)
                    .ok_or_else(|| invalid("more solution values than user doubles"))?;
                *slot = sol_at(x, i)?;
            }
        }
        Ok(sol.found)
    }

    /// `RailEdgeuseIntersect(lRailIndex, pSideEU, dGuessEdgeParameter,
    /// rbFoundIntersection, rTsectPnt, rdEdgeuseParameter)`: the fillet point where
    /// rail `rail_index` meets the side edge use's UV trim curve (built when missing),
    /// solved from the guess parameter on that curve with `ops`'s equations
    /// ([`RailUvCurveIntersectEnfo`]). Returns the edge use parameter
    /// (`rdEdgeuseParameter`) when found (`rbFoundIntersection`); `tsect` then gets both
    /// UV points, the parameter (`m_dCurveParameter`), a crossing type when the type was
    /// undefined or unknown, and the further solution values.
    ///
    /// # Errors
    ///
    /// `SM_ERR` without both offset surfaces or for more than [`MAX_USER_DOUBLES`]
    /// further values; those of `GetOrCreateUVTrimCurve`, the evaluations,
    /// `GetJacobianSize`, `LoadInitialValues` and the solve;
    /// [`SolidError::InvalidInput`] where upstream dereferences `NULL` (no UV trim
    /// curve, no side edge), for a UV trim curve that is not a B-spline, a rail index
    /// above 1, or no second surface index from `LoadInitialValues`.
    pub fn rail_edgeuse_intersect(
        ops: &mut dyn FilletSolverOps,
        k: &Kernel,
        b: &mut Brep,
        rail_index: usize,
        side_eu: EdgeuseId,
        guess_edge_parameter: f64,
        tsect: &mut TsectPnt,
    ) -> Result<Option<f64>> {
        // check state - need both offset surfaces
        if !has_surface(ops.solver(), 0)? || !has_surface(ops.solver(), 1)? {
            return Err(sm_err());
        }
        let rail_srf = srf_of_rail(rail_index)?;
        let other_srf = srf_of_rail(1 - rail_index)?;

        // get original Brep sideEdge->Curve to intersect
        let side_edge = eu::edge(b, side_eu)?;

        // let edgeuse->Surface UV guess = edgeuse->UVTrimCurve(GuessEdgeParameter)
        get_or_create_uv_trim_curve(k, b, side_eu, false, false, true)?.status?;
        let uv_trim = eu::uv_trim_curve(b, side_eu)?
            .0
            .cloned()
            .ok_or_else(|| invalid("no UV trim curve"))?;
        check_uv_bspline_kind(&uv_trim)?;
        let p2d = uv_trim.evaluate_point(guess_edge_parameter)?;
        let uv_guess = Point2::new(p2d.x, p2d.y);

        // Now Set up Solver
        let jacobian_size = ops.jacobian_size()? + 1;
        let mut intervals = ExtentNd::new(jacobian_size);
        set_interval(&mut intervals, 0, uv_trim.natural_interval());

        let base = ops
            .solver()
            .surface(rail_index)?
            .map(|slot| slot.surface.base_surface().clone())
            .ok_or_else(|| invalid("no offset surface"))?;
        let mut domain = base.natural_uv_domain();
        // virtual SetupOffsetExtension: dispatched on the solver's kind
        ops.setup_offset_extension(k, &base, &mut domain)?;
        set_interval(
            &mut intervals,
            1,
            Extent1d::from_min_max(domain.min.x, domain.max.x),
        );
        set_interval(
            &mut intervals,
            2,
            Extent1d::from_min_max(domain.min.y, domain.max.y),
        );

        let mut periodicities = vec![false, false, false];
        let mut guess_t = vec![guess_edge_parameter, uv_guess.x, uv_guess.y];

        let do_surf2_calcs = true; // Do surface2 calculations.
        let surf1_index = 1; // UV values of Surface1

        // given current inputs - place a good guess for upcoming solve in sGuessT
        let mut surface1 = take_surface(ops.solver(), rail_index)?;
        let loaded = ops.load_initial_values(
            k,
            Some(b),
            do_surf2_calcs,
            rail_index,
            &mut surface1,
            surf1_index,
            &mut intervals,
            &mut periodicities,
            &mut guess_t,
        );
        put_surface(ops.solver_mut(), rail_index, surface1)?;
        let loaded = loaded?;

        let side_tol = edge::tolerance(
            b,
            side_edge.ok_or_else(|| invalid("edge use without edge"))?,
        )?;
        let tol = ops.solver().approx_tol_3d + side_tol;

        // specify the rail/curve intersection evaluate function for the NewtonRaphson Solver
        let ivl = uv_trim.natural_interval();
        let mut eval = RailUvEval(RailUvCurveIntersectEnfo::new(
            &uv_trim, ivl, rail_index, ops,
        ));
        let mut ls = LocalSolveNd::new(&mut eval, Some(&intervals), Some(&periodicities));
        ls.set_desired_accuracy(EFF_ZERO_SQ);
        ls.set_boundary_handler(BoundaryHandler::TotalBoundaryHits, 8);
        ls.set_maximum_iterations(30);

        // Solve it
        let sol = ls.solve_it(&guess_t, tol)?;
        if !sol.found {
            return Ok(None);
        }

        // lSurf2Index is set by LoadInitialValues (uninitialised otherwise)
        let surf2_index = loaded
            .ok_or_else(|| invalid("LoadInitialValues set no second surface index"))?
            .surface2_index;
        let x = &sol.x;
        let edgeuse_parameter = sol_at(x, 0)?;

        // set intersection point
        *tsect.uv_pos_mut(rail_srf) =
            Point2::new(sol_at(x, surf1_index)?, sol_at(x, surf1_index + 1)?);
        *tsect.uv_pos_mut(other_srf) =
            Point2::new(sol_at(x, surf2_index)?, sol_at(x, surf2_index + 1)?);
        tsect.curve_parameter = edgeuse_parameter; // Just throw this into here for now.

        // Set this here, otherwise it will be ignored in LoadInitialValues().
        if matches!(
            tsect.point_type,
            IntersectionPointType::Undefined | IntersectionPointType::Unknown
        ) {
            tsect.point_type = IntersectionPointType::Crossing;
        }

        // Load other fillet surface variables being solved for into the fillet point.
        for (cnt, i) in ((surf2_index + 2)..x.len()).enumerate() {
            if cnt == MAX_USER_DOUBLES {
                return Err(sm_err());
            }
            if let Some(slot) = tsect.user_doubles.get_mut(cnt) {
                *slot = sol_at(x, i)?;
            }
        }
        Ok(Some(edgeuse_parameter))
    }

    /// `RailRailIntersect(crUVGuess, lRailIndex, pOtherFillet, lOtherRailIndex,
    /// rbFoundIntersection, rTsectPnt, rOtherTsectPnt)`: where rail `rail_index` of
    /// `ops` meets rail `other_rail_index` of `other` on their common surface, solved
    /// from `uv_guess` with both solvers' equations ([`RailRailIntersectEnfo`]).
    /// Returns `rbFoundIntersection`; when found, both tsect points get their UV points
    /// and further solution values.
    ///
    /// # Errors
    ///
    /// `SM_ERR` when a solver lacks an offset surface or for more than
    /// [`MAX_USER_DOUBLES`] further values; those of `GetJacobianSize`,
    /// `SetupOffsetExtension`, `LoadInitialValues` and the solve;
    /// [`SolidError::Unported`] for the temporary radius change
    /// (`SmTemporaryRadiusChange`, unit `blend.fillet_executive_p2`);
    /// [`SolidError::InvalidInput`] for Jacobian sizes summing below 2, a rail index
    /// above 1 or no second surface index from `LoadInitialValues`.
    pub fn rail_rail_intersect(
        ops: &mut dyn FilletSolverOps,
        other: &mut dyn FilletSolverOps,
        k: &Kernel,
        brep: Option<&Brep>,
        uv_guess: Point2,
        rail_index: usize,
        other_rail_index: usize,
        tsect: &mut TsectPnt,
        other_tsect: &mut TsectPnt,
    ) -> Result<bool> {
        // check inputs
        if !has_surface(ops.solver(), 0)?
            || !has_surface(ops.solver(), 1)?
            || !has_surface(other.solver(), 0)?
            || !has_surface(other.solver(), 1)?
        {
            return Err(sm_err());
        }

        // Set up Solver Size
        let jacobian_size = ops.jacobian_size()?;
        let other_jacobian_size = other.jacobian_size()?;
        let jacobian_size = (jacobian_size + other_jacobian_size)
            .checked_sub(2)
            .ok_or_else(|| invalid("Jacobian sizes sum below 2"))?;

        // compute and save extended domain size for TargetRail->Surface
        let base = ops
            .solver()
            .surface(rail_index)?
            .map(|slot| slot.surface.base_surface().clone())
            .ok_or_else(|| invalid("no offset surface"))?;
        let mut domain = base.natural_uv_domain();
        // virtual SetupOffsetExtension: dispatched on the solver's kind
        ops.setup_offset_extension(k, &base, &mut domain)?;
        let mut intervals = ExtentNd::new(jacobian_size);
        set_interval(
            &mut intervals,
            0,
            Extent1d::from_min_max(domain.min.x, domain.max.x),
        );
        set_interval(
            &mut intervals,
            1,
            Extent1d::from_min_max(domain.min.y, domain.max.y),
        );

        // mark periodicities for TargetRail->Surface
        let mut periodicities = vec![false, false];

        // Set Solver GuessPoint
        let mut guess_t = vec![uv_guess.x, uv_guess.y];

        // Temporarily replace the variable-radius functions with a constant. [Reg Iter:306]
        let this_eu = ops.solver().edgeuse(rail_index)?;
        let other_eu = other.solver().edgeuse(other_rail_index)?;
        let _ = (this_eu, other_eu);
        // UNPORTED(blend.fillet_executive_p2): SmTemporaryRadiusChange::SmTemporaryRadiusChange
        unported_temporary_radius_change()?;

        // Append a good lRailIndex->OtherSurface UVguessPoint into sGuessT
        let do_surf2_calcs = true;
        let surf1_index = 0;
        let mut surface1 = take_surface(ops.solver(), rail_index)?;
        let loaded = ops.load_initial_values(
            k,
            brep,
            do_surf2_calcs,
            rail_index,
            &mut surface1,
            surf1_index,
            &mut intervals,
            &mut periodicities,
            &mut guess_t,
        );
        put_surface(ops.solver_mut(), rail_index, surface1)?;
        let loaded = loaded?;

        // Append a good lOtherRailIndex->OtherSurface UVguessPoint into sGuessT
        let mut other_surface1 = take_surface(other.solver(), other_rail_index)?;
        let other_loaded = other.load_initial_values(
            k,
            brep,
            do_surf2_calcs,
            other_rail_index,
            &mut other_surface1,
            surf1_index,
            &mut intervals,
            &mut periodicities,
            &mut guess_t,
        );
        put_surface(other.solver_mut(), other_rail_index, other_surface1)?;
        let other_loaded = other_loaded?;

        let tol = ops.solver().approx_tol_3d + other.solver().approx_tol_3d;

        // Set Solver Manager
        let mut eval = RailRailEval(RailRailIntersectEnfo::new(
            rail_index,
            ops,
            other_rail_index,
            other,
        ));
        let mut ls = LocalSolveNd::new(&mut eval, Some(&intervals), Some(&periodicities));
        ls.set_desired_accuracy(EFF_ZERO_SQ);
        ls.set_boundary_handler(BoundaryHandler::TotalBoundaryHits, 8);
        ls.set_maximum_iterations(30);

        // Solve it
        let sol = ls.solve_it(&guess_t, tol)?;
        if !sol.found {
            return Ok(false);
        }

        let surf2_index = loaded
            .ok_or_else(|| invalid("LoadInitialValues set no second surface index"))?
            .surface2_index;
        let other_surf2_index = other_loaded
            .ok_or_else(|| invalid("LoadInitialValues set no second surface index"))?
            .surface2_index;
        let x = &sol.x;
        let rail_srf = srf_of_rail(rail_index)?;
        let rail_other_srf = srf_of_rail(1 - rail_index)?;
        *tsect.uv_pos_mut(rail_srf) =
            Point2::new(sol_at(x, surf1_index)?, sol_at(x, surf1_index + 1)?);
        *tsect.uv_pos_mut(rail_other_srf) =
            Point2::new(sol_at(x, surf2_index)?, sol_at(x, surf2_index + 1)?);

        // Load other fillet surface variables being solved for into the fillet point.
        for (cnt, i) in ((surf2_index + 2)..other_surf2_index).enumerate() {
            if cnt == MAX_USER_DOUBLES {
                return Err(sm_err());
            }
            if let Some(slot) = tsect.user_doubles.get_mut(cnt) {
                *slot = sol_at(x, i)?;
            }
        }

        let o_srf = srf_of_rail(other_rail_index)?;
        let o_other_srf = srf_of_rail(1 - other_rail_index)?;
        *other_tsect.uv_pos_mut(o_srf) =
            Point2::new(sol_at(x, surf1_index)?, sol_at(x, surf1_index + 1)?);
        *other_tsect.uv_pos_mut(o_other_srf) = Point2::new(
            sol_at(x, other_surf2_index)?,
            sol_at(x, other_surf2_index + 1)?,
        );

        // Load values into rOtherTsectPnt.
        for (cnt, i) in ((other_surf2_index + 2)..x.len()).enumerate() {
            if cnt == MAX_USER_DOUBLES {
                return Err(sm_err());
            }
            if let Some(slot) = other_tsect.user_doubles.get_mut(cnt) {
                *slot = sol_at(x, i)?;
            }
        }
        Ok(true)
    }

    /// `CreateCliffRollOverSideFG(pRailOnCliff, pSideEdge, pSideFace, bTrimSideFaces,
    /// dExtensionDist, rpSideFG)`: the side fillet geometry of a cliff roll-over. When
    /// trimming, the rail's original face becomes the side face; otherwise the side
    /// face's surface is extended by `extension_dist` (copied as an analytic surface
    /// where possible) and the rail loses its original face. At a fillet corner the
    /// corner vertex moves to where the cross-section edge meets that surface and the
    /// cross section is trimmed to it (and, without trimming and without a split-face
    /// edge, a split-face edge is made from the corner to the side edge). Without
    /// trimming a new cliff-side-patch geometry is made, with the rail and a new rail
    /// on the side edge; it is returned (`rpSideFG`), `None` when trimming.
    ///
    /// `b` is the B-rep being filleted (the side edge and face); `pseudo` is the
    /// executive's pseudo-B-rep holding `rail_on_cliff` and the fillet vertices.
    ///
    /// # Errors
    ///
    /// [`SolidError::InvalidInput`] without an executive or where upstream dereferences
    /// `NULL`; [`SolidError::Unported`] from the pseudo-B-rep lookup
    /// (`SmFilletExecutive::GetPseudoBrep`, unit `blend.fillet_executive_p1`) and the
    /// fillet vertex, edge and geometry calls after it; `SM_ERR_NULL_POINTER` (`NER`)
    /// for a side face without surface (when extending), no cross-section edge or
    /// curve, or a side edge without curve; `SM_ERR` when the cross section misses the
    /// side surface or the corner point does not drop onto the side edge; those of the
    /// extension, copy, intersections, evaluations and drops.
    pub fn create_cliff_roll_over_side_fg(
        &mut self,
        k: &Kernel,
        b: &Brep,
        pseudo: &mut Brep,
        rail_on_cliff: EdgeId,
        side_edge: EdgeId,
        side_face: FaceId,
        trim_side_faces: bool,
        extension_dist: f64,
    ) -> Result<Option<FilletGeomRef>> {
        const OP_C: &str = "blend::fillet_solver::p2::create_cliff_roll_over_side_fg";
        let nul = || SolidError::null_pointer(OP_C);

        // locals
        let rail_start_v = edge::vertex(pseudo, rail_on_cliff)?;
        let rail_end_v = edge::other_vertex(pseudo, rail_on_cliff, rail_start_v)?;
        if self.executive.is_none() {
            return Err(invalid("no fillet executive"));
        }
        // GetFilletExecutive()->GetPseudoBrep() is `pseudo`
        // UNPORTED(blend.fillet_executive_p1): SmFilletExecutive::GetPseudoBrep
        unported_on::<_, ()>(
            self.executive,
            "SmFilletExecutive::GetPseudoBrep",
            "blend.fillet_executive_p1",
        )?;
        let start_v = rail_start_v.ok_or_else(|| invalid("rail without start vertex"))?;
        // UNPORTED(blend.fillet_geom_p1): SmFilletVertex::GetFilletCorner
        let start_corner: Option<FilletCornerRef> = unported_on(
            start_v,
            "SmFilletVertex::GetFilletCorner",
            "blend.fillet_geom_p1",
        )?;
        let adj_corner_vert = if start_corner.is_some() {
            Some(start_v)
        } else {
            let end_v = rail_end_v.ok_or_else(|| invalid("rail without end vertex"))?;
            // UNPORTED(blend.fillet_geom_p1): SmFilletVertex::GetFilletCorner
            let end_corner: Option<FilletCornerRef> = unported_on(
                end_v,
                "SmFilletVertex::GetFilletCorner",
                "blend.fillet_geom_p1",
            )?;
            end_corner.map(|_| end_v)
        };

        // Pick the side surface to extend: the side face's surface or its copy
        // (SM_FILLETSURF_TYPE is SmSurface: the cast keeps any surface)
        let adj_surface = face_p2::surface(b, side_face)?.cloned();
        let side_surface = if trim_side_faces {
            // UNPORTED(blend.fillet_geom_p1): SmFilletEdge::SetOriginalFace
            unported_on::<_, ()>(
                (rail_on_cliff, Some(side_face)),
                "SmFilletEdge::SetOriginalFace",
                "blend.fillet_geom_p1",
            )?;
            adj_surface
        } else {
            let adj = adj_surface.ok_or_else(nul)?;
            let new_surface = adj.create_extended_surface(k, extension_dist, ContinuityType::G1)?;
            // Copy pNewSurface, when possible as an analytic surface
            let copy = new_surface.copy_and_add_analytics(k)?;
            // UNPORTED(blend.fillet_geom_p1): SmFilletEdge::SetOriginalFace
            unported_on::<_, ()>(
                (rail_on_cliff, None::<FaceId>),
                "SmFilletEdge::SetOriginalFace",
                "blend.fillet_geom_p1",
            )?;
            Some(copy)
        };

        // At a corner: move the corner vertex to the cross section / side surface
        // intersection and trim the cross section to run from there.
        if let Some(corner_v) = adj_corner_vert {
            let corner_edges = vertex::edges(pseudo, corner_v, None)?;

            // the last CROSS_SECTION and SPLIT_FACE fillet edges of the list
            let mut split_edge: Option<EdgeId> = None;
            let mut cross_section: Option<EdgeId> = None;
            for fe in &corner_edges {
                // UNPORTED(blend.fillet_geom_p1): SmFilletEdge::GetFilletEdgeType
                let fe_type: i64 = unported_on(
                    *fe,
                    "SmFilletEdge::GetFilletEdgeType",
                    "blend.fillet_geom_p1",
                )?;
                if fe_type == FILLET_EDGE_CROSS_SECTION {
                    cross_section = Some(*fe);
                } else if fe_type == FILLET_EDGE_SPLIT_FACE {
                    split_edge = Some(*fe);
                }
            }
            let cross_section = cross_section.ok_or_else(nul)?;

            // CrossSection curve locals
            let mut cs_curve = edge::curve(pseudo, cross_section)?
                .cloned()
                .ok_or_else(nul)?;
            let mut trim_ivl = cs_curve.natural_interval();

            // Intersect the cross section with the extended side surface
            let side = side_surface
                .as_ref()
                .ok_or_else(|| invalid("no side surface"))?;
            let solutions = side.global_curve_intersect(
                k,
                [CliffHandle::Curve, CliffHandle::Surface],
                &side.natural_uv_domain(),
                &cs_curve,
                &trim_ivl,
                self.approx_tol_3d,
            )?;
            let first = solutions.get(0).ok_or_else(sm_err)?;
            let param1 = first.start.param(0)?;
            let corner_pnt = cs_curve.evaluate_point(param1)?;

            // Update the corner vertex to the cross section / side face point
            // UNPORTED(topology.vertex): SmVertex::SetPoint (justified: findings.md "blend.fillet_solver_p1: pseudo-B-rep mutations")
            unported_on::<_, ()>(
                (corner_v, corner_pnt),
                "SmVertex::SetPoint",
                "topology.vertex",
            )?;

            // Drop the cross section's other vertex onto its curve
            let other_v = edge::other_vertex(pseudo, cross_section, Some(corner_v))?
                .ok_or_else(|| invalid("cross section without other vertex"))?;
            let dropped = cs_curve.drop_point(
                k,
                &trim_ivl,
                vertex::point(pseudo, other_v)?,
                None,
                self.approx_tol_3d,
                None,
                SolverOperation::Minimize,
            )?;
            if dropped.success {
                let param2 = dropped.parameter;
                // order the desired cross section interval
                trim_ivl = if param1 > param2 {
                    Extent1d::from_min_max(param2, param1)
                } else {
                    Extent1d::from_min_max(param1, param2)
                };

                // The cross section's curves belong to its edge: their trims read the
                // edge's tolerance as zone tolerance (as `sm_TrimRail`, p1 notes).
                let tk = Kernel::new(&crate::kernel::KernelOptions {
                    tolerances: crate::kernel::Tolerances {
                        zone_tol_3d: edge::tolerance(pseudo, cross_section)?,
                        ..k.tolerances
                    },
                    ..crate::kernel::KernelOptions::default()
                });
                // Set the cross section curve interval (status ignored upstream)
                ignore_status(cs_curve.trim(&tk, &mut trim_ivl, TrimMode::Notify))?;
                if let Some(c) = pseudo.edge_mut(cross_section)?.curve.as_mut() {
                    *c = cs_curve.clone();
                }

                // Set the cross section's UV trim curve intervals
                let prim = eu::primary_edgeuse(pseudo, cross_section)?
                    .ok_or_else(|| invalid("cross section without edge use"))?;
                let mate =
                    eu::mate(pseudo, prim)?.ok_or_else(|| invalid("edge use without mate"))?;
                for e_u in [prim, mate] {
                    if let Some(uv) = pseudo.edgeuse_mut(e_u)?.uv_trim_curve.as_mut() {
                        ignore_status(uv.trim(&tk, &mut trim_ivl, TrimMode::Notify))?;
                    }
                }
            }
            // else: SM_ASSERT_MSG only (the cross section interval stays)

            // Without a split-face edge and not trimming: make one
            if !trim_side_faces && split_edge.is_none() {
                // Create a split fillet edge between the corner point and its
                // projection onto the side curve
                // UNPORTED(blend.fillet_geom_p2): SmFilletVertex::SmFilletVertex (justified: findings.md "blend.fillet_solver_p2: fillet-geometry stubs need a FilletBrep")
                let end_vert: VertexId = unported_on(
                    &*pseudo,
                    "SmFilletVertex::SmFilletVertex",
                    "blend.fillet_geom_p2",
                )?;
                // UNPORTED(blend.fillet_geom_p5): SmFilletEdge::MakeFilletEdge
                let new_split: EdgeId = unported_on(
                    (end_vert, corner_v),
                    "SmFilletEdge::MakeFilletEdge",
                    "blend.fillet_geom_p5",
                )?;
                // UNPORTED(blend.fillet_geom_p1): SmFilletEdge::SetFilletEdgeType
                unported_on::<_, ()>(
                    (new_split, FILLET_EDGE_SPLIT_FACE),
                    "SmFilletEdge::SetFilletEdgeType",
                    "blend.fillet_geom_p1",
                )?;

                // drop the corner point to the side curve
                let side_curve = edge::curve(b, side_edge)?.ok_or_else(nul)?;
                let sols = side_curve.global_point_solve(
                    k,
                    &edge::interval(b, side_edge)?,
                    SolverOperation::Minimize,
                    corner_pnt,
                    EFF_ZERO_SQRT,
                    SolutionRequested::Single,
                    &GlobalPointSolveOptions::<crate::topology::ObjectHandle>::default(),
                )?;
                let sol = sols.get(0).ok_or_else(sm_err)?;
                let (end_pnt, end_tan) = curve_pv(side_curve, sol.start.param(0)?)?;
                // UNPORTED(topology.vertex): SmVertex::SetPoint (justified: findings.md "blend.fillet_solver_p1: pseudo-B-rep mutations")
                unported_on::<_, ()>((end_vert, end_pnt), "SmVertex::SetPoint", "topology.vertex")?;
                // UNPORTED(blend.fillet_geom_p1): SmFilletVertex::SetStatus
                unported_on::<_, ()>(
                    (end_vert, FilStatus::Processed),
                    "SmFilletVertex::SetStatus",
                    "blend.fillet_geom_p1",
                )?;

                // the plane through the end point perpendicular to the line
                // [corner, end] in the side curve's tangent plane
                let dir = corner_pnt - end_pnt;
                let plane_normal = end_tan.cross(dir).cross(dir).unitize()?;
                let mut plane = Plane::from_normal(end_pnt, plane_normal);
                plane.make_nurb()?;

                // the side surface / plane intersection curve (its UV curves are
                // deleted unused); the vertices' points are sEndPnt and the corner
                // point they were just set to
                let mut start_dir = dir;
                let plane_surface = Surface::from(plane);
                let curves = srf_srf_intersection(
                    k,
                    Some(end_pnt),
                    Some(corner_pnt),
                    Some(&mut start_dir),
                    None,
                    side,
                    &plane_surface,
                    self.approx_tol_3d,
                    30.0 * PI / 180.0,
                    SrfSrfOptions::default(),
                )?;
                // p3DCurve->SetOwner dereferences the curve
                let curve_3d = curves.curve_3d.ok_or_else(nul)?;
                // SetCurve(p3DCurve, TRUE) and SetOwner(pSplitEdge)
                // UNPORTED(topology.edge): SmEdge::SetCurve (justified: findings.md "blend.fillet_solver_p1: pseudo-B-rep mutations")
                unported_on::<_, ()>((new_split, curve_3d), "SmEdge::SetCurve", "topology.edge")?;
            }
        }

        // Trimming: done after setting the original face and the corner
        if trim_side_faces {
            return Ok(None);
        }

        // Create the new fillet geometry
        // UNPORTED(blend.fillet_geom_p5): SmFilletGeom::SmFilletGeom
        let side_fg: FilletGeomRef = unported_on(
            FILLET_GEOM_CLIFF_SIDE_PATCH,
            "SmFilletGeom::SmFilletGeom",
            "blend.fillet_geom_p5",
        )?;
        self.fillet_geoms.push(side_fg);

        let mut rail_eu = eu::primary_edgeuse(pseudo, rail_on_cliff)?
            .ok_or_else(|| invalid("rail without edge use"))?;
        if eu_vertex(pseudo, rail_eu)? != rail_start_v {
            rail_eu = eu::mate(pseudo, rail_eu)?.ok_or_else(|| invalid("edge use without mate"))?;
        }
        // UNPORTED(blend.fillet_geom_p1): SmFilletEdgeuse::SetFilletGeom
        unported_on::<_, ()>(
            (rail_eu, side_fg),
            "SmFilletEdgeuse::SetFilletGeom",
            "blend.fillet_geom_p1",
        )?;

        // Attach the surface
        // UNPORTED(blend.fillet_geom_p1): SmFilletGeom::SetFilletSurface
        unported_on::<_, ()>(
            (side_fg, side_surface),
            "SmFilletGeom::SetFilletSurface",
            "blend.fillet_geom_p1",
        )?;

        // Make some new topology objects (edge, edge uses, vertex uses)
        // UNPORTED(blend.fillet_geom_p1): SmFilletGeom::SetRail
        unported_on::<_, ()>(
            (side_fg, 0usize, rail_on_cliff),
            "SmFilletGeom::SetRail",
            "blend.fillet_geom_p1",
        )?;
        // UNPORTED(blend.fillet_geom_p2): SmFilletEdge::SmFilletEdge (justified: findings.md "blend.fillet_solver_p2: fillet-geometry stubs need a FilletBrep")
        let new_rail: EdgeId = unported_on(
            &*pseudo,
            "SmFilletEdge::SmFilletEdge",
            "blend.fillet_geom_p2",
        )?;
        // UNPORTED(blend.fillet_geom_p1): SmFilletGeom::SetRail
        unported_on::<_, ()>(
            (side_fg, 1usize, new_rail),
            "SmFilletGeom::SetRail",
            "blend.fillet_geom_p1",
        )?;
        // UNPORTED(blend.fillet_geom_p1): SmFilletEdge::SetFilletEdgeType
        unported_on::<_, ()>(
            (new_rail, FILLET_EDGE_RAIL),
            "SmFilletEdge::SetFilletEdgeType",
            "blend.fillet_geom_p1",
        )?;
        let mut new_eus: [Option<EdgeuseId>; 2] = [None, None];
        for slot in &mut new_eus {
            // UNPORTED(blend.fillet_geom_p1): SmFilletEdgeuse::SmFilletEdgeuse
            *slot = Some(unported_on(
                &*pseudo,
                "SmFilletEdgeuse::SmFilletEdgeuse",
                "blend.fillet_geom_p1",
            )?);
        }
        let [Some(new_eu1), Some(new_eu2)] = new_eus else {
            return Err(invalid("new edge uses missing"));
        };
        // UNPORTED(topology.proto_topology_p1): SmOwningTopology::PreInsert (justified: findings.md "blend.fillet_solver_p1: pseudo-B-rep mutations")
        unported_on::<_, ()>(
            (new_rail, new_eu1),
            "SmOwningTopology::PreInsert",
            "topology.proto_topology_p1",
        )?;
        // UNPORTED(topology.proto_topology_p1): SmOwningTopology::PostInsert (justified: findings.md "blend.fillet_solver_p1: pseudo-B-rep mutations")
        unported_on::<_, ()>(
            (new_rail, new_eu2),
            "SmOwningTopology::PostInsert",
            "topology.proto_topology_p1",
        )?;
        let mut new_vus: [Option<VertexuseId>; 2] = [None, None];
        for slot in &mut new_vus {
            // UNPORTED(blend.fillet_geom_p1): SmFilletVertexuse::SmFilletVertexuse
            *slot = Some(unported_on(
                &*pseudo,
                "SmFilletVertexuse::SmFilletVertexuse",
                "blend.fillet_geom_p1",
            )?);
        }
        let [Some(new_vu1), Some(new_vu2)] = new_vus else {
            return Err(invalid("new vertex uses missing"));
        };
        for (vu, e_u) in [(new_vu1, new_eu1), (new_vu2, new_eu2)] {
            // UNPORTED(blend.fillet_geom_p1): SmFilletVertexuse::SetProperty
            unported_on::<_, ()>(
                (vu, e_u),
                "SmFilletVertexuse::SetProperty",
                "blend.fillet_geom_p1",
            )?;
        }
        for (e_u, vu) in [(new_eu1, new_vu1), (new_eu2, new_vu2)] {
            // UNPORTED(blend.fillet_geom_p1): SmFilletEdgeuse::SetFilletVertexuse
            unported_on::<_, ()>(
                (e_u, vu),
                "SmFilletEdgeuse::SetFilletVertexuse",
                "blend.fillet_geom_p1",
            )?;
        }

        // Put each new vertex use into its rail end vertex, or into the other end of
        // a 'support' edge at that vertex (attaching the new geometry to it)
        for (rail_v, new_vu) in [(rail_start_v, new_vu1), (rail_end_v, new_vu2)] {
            let rail_v = rail_v.ok_or_else(|| invalid("rail without vertex"))?;
            let mut new_rail_v = rail_v;
            for fe in vertex::edges(pseudo, rail_v, None)? {
                // UNPORTED(blend.fillet_geom_p1): SmFilletEdge::GetFilletEdgeType
                let fe_type: i64 = unported_on(
                    fe,
                    "SmFilletEdge::GetFilletEdgeType",
                    "blend.fillet_geom_p1",
                )?;
                if fe_type == FILLET_EDGE_ON_EXTEND_EDGE || fe_type == FILLET_EDGE_SPLIT_FACE {
                    new_rail_v = edge::other_vertex(pseudo, fe, Some(rail_v))?
                        .ok_or_else(|| invalid("support edge without other vertex"))?;
                    // Attach the new fillet geometry to the support edge
                    let mut feu = eu::primary_edgeuse(pseudo, fe)?
                        .ok_or_else(|| invalid("support edge without edge use"))?;
                    // UNPORTED(blend.fillet_geom_p1): SmFilletEdgeuse::GetFilletGeom
                    let feu_fg: Option<FilletGeomRef> = unported_on(
                        feu,
                        "SmFilletEdgeuse::GetFilletGeom",
                        "blend.fillet_geom_p1",
                    )?;
                    if feu_fg.is_some() {
                        feu = eu::mate(pseudo, feu)?
                            .ok_or_else(|| invalid("edge use without mate"))?;
                    }
                    // UNPORTED(blend.fillet_geom_p1): SmFilletGeom::AddSideFilletEdgeuse
                    unported_on::<_, ()>(
                        (side_fg, feu),
                        "SmFilletGeom::AddSideFilletEdgeuse",
                        "blend.fillet_geom_p1",
                    )?;
                    // UNPORTED(blend.fillet_geom_p1): SmFilletEdgeuse::SetFilletGeom
                    unported_on::<_, ()>(
                        (feu, side_fg),
                        "SmFilletEdgeuse::SetFilletGeom",
                        "blend.fillet_geom_p1",
                    )?;
                    break;
                }
            }
            // UNPORTED(topology.proto_topology_p1): SmOwningTopology::PostInsert (justified: findings.md "blend.fillet_solver_p1: pseudo-B-rep mutations")
            unported_on::<_, ()>(
                (new_rail_v, new_vu),
                "SmOwningTopology::PostInsert",
                "topology.proto_topology_p1",
            )?;
        }

        // The new rail runs on the side edge
        // UNPORTED(blend.fillet_geom_p1): SmFilletEdge::SetOriginalEdge
        unported_on::<_, ()>(
            (new_rail, side_edge),
            "SmFilletEdge::SetOriginalEdge",
            "blend.fillet_geom_p1",
        )?;
        // UNPORTED(blend.fillet_geom_p1): SmFilletEdge::SetOriginalFace
        unported_on::<_, ()>(
            (new_rail, Some(side_face)),
            "SmFilletEdge::SetOriginalFace",
            "blend.fillet_geom_p1",
        )?;
        // UNPORTED(blend.fillet_geom_p6): SmFilletEdge::CalcEdgeOnEdge
        unported_on::<_, ()>(
            (new_rail, k),
            "SmFilletEdge::CalcEdgeOnEdge",
            "blend.fillet_geom_p6",
        )?;
        Ok(Some(side_fg))
    }

    /// `RecordFilletError(eType, rapTSPs, cComment)`: passes the error, with the
    /// filleted edge, the fillet surfaces, centre lines and rails built so far and the
    /// first tsect point's curve position, to the executive (`NoteFilletError`).
    ///
    /// # Errors
    ///
    /// [`SolidError::InvalidInput`] where upstream dereferences `NULL` (no edge use or
    /// executive); [`SolidError::Unported`] for the fillet geometry queries
    /// (`SmFilletGeom::GetFilletSurface`, unit `blend.fillet_geom_p1`) and the report
    /// (`SmFilletExecutive::NoteFilletError`, unit `blend.fillet_executive_p1`).
    pub fn record_fillet_error(
        &self,
        b: &Brep,
        error_type: i64,
        tsects: &[TsectPnt],
        comment: &str,
    ) -> Result<()> {
        let filleted = filleted_edge_opt(self, b)?;

        // Collect all existing fillet surfaces and rail curves (upstream adds rail 0
        // where it tests rail 1, WR-190)
        let geoms = self.fillet_geoms();
        if let Some(g) = geoms.first() {
            let _ = g;
            // UNPORTED(blend.fillet_geom_p1): SmFilletGeom::GetFilletSurface
            return Err(SolidError::unported(
                "SmFilletGeom::GetFilletSurface",
                "blend.fillet_geom_p1",
            ));
        }

        // Just look at one TsectPnt position for now
        let mut tsp_pos = Vec3::UNDEF;
        if let Some(t) = tsects.first() {
            if t.crv_pos().is_initialized() {
                tsp_pos = t.crv_pos();
            }
        }

        // pass the call along to FilletExecutive
        if self.executive.is_none() {
            return Err(invalid("no fillet executive"));
        }
        let _ = (error_type, self.status, filleted, tsp_pos, comment);
        // UNPORTED(blend.fillet_executive_p1): SmFilletExecutive::NoteFilletError
        Err(SolidError::unported(
            "SmFilletExecutive::NoteFilletError",
            "blend.fillet_executive_p1",
        ))
    }
}

/// `GetEdgeuse(0)->GetEdge()` (a missing edge use is a `NULL` dereference upstream).
fn filleted_edge_opt(s: &FilletSolver, b: &Brep) -> Result<Option<EdgeId>> {
    let eu0 = s.edgeuse(0)?.ok_or_else(|| invalid("no edge use"))?;
    eu::edge(b, eu0)
}

fn unported_temporary_radius_change() -> Result<()> {
    Err(SolidError::unported(
        "SmTemporaryRadiusChange::SmTemporaryRadiusChange",
        "blend.fillet_executive_p2",
    ))
}

/// A call into a unit that is not ported yet: `upstream` (unit `unit`) on `_args` (the
/// receiver and arguments of the C++ call, kept so each call site reads as upstream).
/// Every call site carries its `UNPORTED(<unit>)` tag.
fn unported_on<A, T>(_args: A, upstream: &'static str, unit: &'static str) -> Result<T> {
    Err(SolidError::unported(upstream, unit))
}

/// `(*pCC)[i]` (an index past the classification is out of bounds upstream).
fn interval_at<O: ClassObject>(
    cc: &CurveClassification<O>,
    i: usize,
) -> Result<&crate::intersect::curve_class::p1::CurveInterval<O>> {
    cc.intervals()
        .get(i)
        .ok_or_else(|| invalid("classification index out of range"))
}

/// `(SmEdge*)rPC.GetObject()` through the caller's handle map.
fn edge_object<O: ClassObject>(
    pc: &crate::intersect::point_class::PointClassification<O>,
    edge_of: &dyn Fn(O) -> Option<EdgeId>,
) -> Result<EdgeId> {
    pc.object()
        .and_then(edge_of)
        .ok_or_else(|| invalid("classification object is not an edge"))
}

/// `pEU->GetVertexuse()->GetVertex()`.
fn eu_vertex(b: &Brep, e: EdgeuseId) -> Result<Option<VertexId>> {
    let vu = eu::vertexuse(b, e)?.ok_or_else(|| invalid("edge use without vertex use"))?;
    vertexuse::vertex(b, vu)
}

/// `sm_FindRolloverSideEdge(pCC, rlCurrIndx, pRail, pOrigRailCurve, pOrigRailUVCurve,
/// pOrigFillet, bConvexFillet, rpSideEdge, rbIsTangentialSideEdge, rbTrimSideFaces,
/// rbOppositeNormals, rpRolloverCornerEU, rpRolloverCornerV, dTangencyTolRad,
/// iDebugLevel)`: walks the classified rail intervals from `io.curr_index` to the side
/// edge the rail runs into, telling a crash (trim the side faces) from a cliff (roll
/// over, possibly at a corner). `io` holds the reference arguments: upstream writes
/// only some of them on each path, the rest keep the caller's values. `pseudo` holds
/// the rail and its fillet vertices; `b` is the B-rep being filleted (the side edges,
/// their edge uses and faces; the tangency test may build UV trim curves on it);
/// `edge_of` maps a classification object to its edge (`(SmEdge*)GetObject()`).
/// `pOrigRailCurve` is unused upstream.
///
/// # Errors
///
/// [`SolidError::Unported`] from the rail's original face
/// (`SmFilletEdge::GetOriginalFace`, unit `blend.fillet_geom_p1`) on, and for
/// `SmFilletVertex::IsProcessed`; `SM_ERR_NULL_POINTER` (`NER`) without original face;
/// `SM_ERR` for an end class that is not an edge before the last interval, or a start
/// class that is not an edge on it; those of the evaluations, the angle and the
/// tangency and convexity tests; [`SolidError::InvalidInput`] where upstream
/// dereferences `NULL` or indexes past the classification.
pub fn find_rollover_side_edge<O: ClassObject>(
    k: &Kernel,
    pseudo: &Brep,
    b: &mut Brep,
    cc: &CurveClassification<O>,
    rail: EdgeId,
    orig_rail_uv_curve: &Curve,
    orig_fillet: &Surface,
    convex_fillet: bool,
    tangency_tol_rad: f64,
    edge_of: &dyn Fn(O) -> Option<EdgeId>,
    io: &mut RolloverSideEdge,
) -> Result<()> {
    const OP_R: &str = "blend::fillet_solver::p2::find_rollover_side_edge";
    let total_intervals = cc.intervals().len();
    let start_v = edge::vertex(pseudo, rail)?;
    let end_v = edge::other_vertex(pseudo, rail, start_v)?;
    // UNPORTED(blend.fillet_geom_p1): SmFilletEdge::GetOriginalFace
    let orig_face: Option<FaceId> = unported_on(
        rail,
        "SmFilletEdge::GetOriginalFace",
        "blend.fillet_geom_p1",
    )?;
    let orig_face = orig_face.ok_or_else(|| SolidError::null_pointer(OP_R))?;
    let is_processed = |v: Option<VertexId>| -> Result<bool> {
        // UNPORTED(blend.fillet_geom_p1): SmFilletVertex::IsProcessed
        unported_on(v, "SmFilletVertex::IsProcessed", "blend.fillet_geom_p1")
    };

    io.trim_side_faces = false;
    // Every pass that does not return advances `curr_index`, so the walk ends within
    // the classification's length (module notes).
    loop {
        let ivl = interval_at(cc, io.curr_index)?;

        if io.curr_index == 0
            && !is_processed(start_v)?
            && ivl.start.point_class() == PointClassType::Unknown
        {
            io.curr_index += 1; // Skip this interval, move on to next.
            continue;
        }

        let start_pc = &ivl.start;
        let end_pc = &ivl.end;
        if io.curr_index + 1 < total_intervals {
            if ivl.end.point_class() != PointClassType::Edge {
                return Err(sm_err()); // Unknown case
            }
            let side_edge = edge_object(&ivl.end, edge_of)?;

            // Check whether the side edge we hit is a tangent edge, and whether it's
            // convex.
            let (tangent, opposite) =
                test_manifold_tangent_rollover(k, b, side_edge, tangency_tol_rad)?;
            io.is_tangential_side_edge = tangent;
            io.opposite_normals = opposite;
            let convex_side_edge =
                test_side_edge_convexity(k, b, side_edge, io.is_tangential_side_edge)?;

            // Same convexity: the fillet crashes into the adjoining face (trim it);
            // otherwise it hangs over a cliff.
            if convex_fillet == convex_side_edge {
                // Trimming process
                io.trim_side_faces = true;
                io.curr_index += 1;
                continue;
            }

            if start_pc.point_class() == PointClassType::Edge
                && end_pc.point_class() == PointClassType::Edge
                && ivl.mid.point_class() == PointClassType::Unknown
            {
                let e1 = edge_object(start_pc, edge_of)?;
                let e2 = edge_object(end_pc, edge_of)?;
                let mut eu_on_orig_face =
                    eu::primary_edgeuse(b, e1)?.ok_or_else(|| invalid("edge without edge use"))?;
                if eu::face(b, eu_on_orig_face)? != Some(orig_face) {
                    eu_on_orig_face = eu::radial(b, eu_on_orig_face)?
                        .ok_or_else(|| invalid("edge use without radial"))?;
                }

                let start_param = start_pc.t_param();
                let end_param = end_pc.t_param();
                let e1_curve = edge::curve(b, e1)?.ok_or_else(|| invalid("edge without curve"))?;
                let (_, tangent_vec) = curve_pv(e1_curve, start_param)?;
                // Compute 'binormal' of rail curve
                let rail_param = ivl.interval.min();
                let pnt = orig_rail_uv_curve.evaluate_point(rail_param)?;
                let uv = Point2::new(pnt.x, pnt.y);
                let flags = crate::nurbs::b_spline_surface_n_lib::p2::EvaluateFlags {
                    only_upper_half: true,
                    ..Default::default()
                };
                let pvv = orig_fillet.evaluate(uv, 1, 1, flags)?;
                let mut bin_vec = pvv
                    .first()
                    .and_then(|row| row.get(1))
                    .copied()
                    .ok_or_else(|| invalid("short surface evaluation"))?;
                let center_uv = orig_fillet.natural_uv_domain().evaluate(0.5, 0.5);
                if uv.y > center_uv.y {
                    bin_vec = -bin_vec;
                }

                let angle = bin_vec.angle_between(tangent_vec)?;
                let eu_orient = if angle < PI / 2.0 {
                    OrientType::Same
                } else {
                    OrientType::Opposite
                };

                if let Some(corner_eu) = io.rollover_corner_eu {
                    let radial_eu = eu::radial(b, corner_eu)?
                        .ok_or_else(|| invalid("edge use without radial"))?;
                    let cw = eu::cw_edgeuse(b, radial_eu)?
                        .ok_or_else(|| invalid("edge use without CW edge use"))?;
                    let cw_edge = eu::edge(b, cw)?;
                    if io.is_tangential_side_edge {
                        let cw_edge = cw_edge.ok_or_else(|| invalid("edge use without edge"))?;
                        // rbOppositeNormals only; the tangency result is dropped
                        let (_, opposite) =
                            test_manifold_tangent_rollover(k, b, cw_edge, tangency_tol_rad)?;
                        io.opposite_normals = opposite;
                    }
                    if cw_edge != Some(e2) {
                        io.rollover_corner_v = eu_vertex(b, cw)?;
                        io.rollover_corner_eu = eu::cw_edgeuse(b, cw)?;
                        io.side_edge = cw_edge;
                        return Ok(());
                    }
                    io.rollover_corner_eu = None;
                    io.rollover_corner_v = None;
                } else {
                    if eu::orientation(b, eu_on_orig_face)? != eu_orient {
                        eu_on_orig_face = eu::mate(b, eu_on_orig_face)?
                            .ok_or_else(|| invalid("edge use without mate"))?;
                    }
                    let eu_r = eu::radial(b, eu_on_orig_face)?
                        .ok_or_else(|| invalid("edge use without radial"))?;
                    let mut rollover_corner = e1 != e2;
                    if e1 == e2
                        && edge::is_closed(b, e1)?
                        && ((start_param < end_param && eu_orient == OrientType::Opposite)
                            || (start_param > end_param && eu_orient == OrientType::Same))
                    {
                        rollover_corner = true;
                    }
                    if rollover_corner {
                        let eu_r_edge =
                            eu::edge(b, eu_r)?.ok_or_else(|| invalid("edge use without edge"))?;
                        let (tangent, opposite) =
                            test_manifold_tangent_rollover(k, b, eu_r_edge, tangency_tol_rad)?;
                        io.is_tangential_side_edge = tangent;
                        io.opposite_normals = opposite;
                        io.rollover_corner_v = eu_vertex(b, eu_r)?;
                        io.rollover_corner_eu = eu::cw_edgeuse(b, eu_r)?;
                        io.side_edge = if io.is_tangential_side_edge {
                            let corner_eu = io
                                .rollover_corner_eu
                                .ok_or_else(|| invalid("edge use without CW edge use"))?;
                            eu::edge(b, corner_eu)?
                        } else {
                            Some(eu_r_edge)
                        };
                        return Ok(());
                    }
                }
            }

            // If this is the one before the very last segment, do the following check
            if io.curr_index + 2 == total_intervals {
                let next_ivl = interval_at(cc, io.curr_index + 1)?;
                if next_ivl.end.point_class() == PointClassType::Unknown && !is_processed(end_v)? {
                    io.curr_index += 1;
                    return Ok(());
                }
            }
            io.rollover_corner_eu = None;
            io.side_edge = Some(side_edge);
            return Ok(());
        }

        // Process the very last segment
        if ivl.start.point_class() != PointClassType::Edge {
            return Err(sm_err()); // Unknown case
        }
        if ivl.end.point_class() == PointClassType::Unknown && !is_processed(end_v)? {
            return Ok(());
        }
        let side_edge = edge_object(&ivl.start, edge_of)?;

        // Check whether the side edge we hit is tangent, or concave/convex.
        let (tangent, opposite) =
            test_manifold_tangent_rollover(k, b, side_edge, tangency_tol_rad)?;
        io.is_tangential_side_edge = tangent;
        io.opposite_normals = opposite;
        let convex_side_edge =
            test_side_edge_convexity(k, b, side_edge, io.is_tangential_side_edge)?;

        // Cliff or crash (comment above).
        if convex_fillet == convex_side_edge {
            // Trimming process
            io.trim_side_faces = true;
            return Ok(());
        }
        io.side_edge = Some(side_edge);
        return Ok(());
    }
}
