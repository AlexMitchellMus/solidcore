//! The trimmed surface cache `SmTrimSrfCache`: a surface cache whose decomposition
//! tree also records the face's trim curves (edge uses), vertices and poly edges per
//! node, classifies its nodes in, on or out of the face, and tests and classifies
//! parameter points against the face boundary.
//!
//! Upstream: `source/SMLib/src/SmTrimSrfCache.cpp`, `source/SMLib/inc/SmTrimSrfCache.h`.
//!
//! | Rust                                      | Upstream                                   |
//! |-------------------------------------------|--------------------------------------------|
//! | [`TrimSrfCache`]                          | `SmTrimSrfCache`                           |
//! | [`TrimSrfCache::new`]                     | `SmTrimSrfCache::SmTrimSrfCache`           |
//! | `Drop` (implicit)                         | `SmTrimSrfCache::~SmTrimSrfCache`          |
//! | [`TYPE`], [`TYPE_STRING`], [`TrimSrfCache::type_code`], [`TrimSrfCache::type_string`] | `GetType`, `GetTypeString` (`SM_COMMON`) |
//! | [`is_kind_of`], [`TrimSrfCache::is_kind_of`] | `SmTrimSrfCache::IsKindOf`              |
//! | [`TrimSrfCache::point_test`]              | `SmTrimSrfCache::PointTest`                |
//! | [`TrimSrfCache::point_classify`]          | `SmTrimSrfCache::PointClassify`            |
//! | [`TrimSrfCache::build_tree`]              | `SmTrimSrfCache::BuildTree`                |
//! | [`TrimSrfCache::implant_edgeuse`]         | `SmTrimSrfCache::ImplantEdgeuse`           |
//! | [`TrimSrfCache::implant_poly_edge`], [`PolyEdgeInput`] | `SmTrimSrfCache::ImplantPolyEdge` |
//! | [`TrimSrfCache::implant_vertexuse`], [`VertexuseInput`] | `SmTrimSrfCache::ImplantVertexuse` |
//!
//! Deviations:
//! - **Inheritance.** `SmTrimSrfCache` derives from `SmSurfaceCache`; here it wraps a
//!   [`SurfaceCache`] ([`TrimSrfCache::base`], [`TrimSrfCache::base_mut`]) and adds
//!   the overriding methods. The constructor's `dAngleTolDeg` is handed to the base
//!   constructor's `dAngTolRad` unconverted, as upstream does: pass it as
//!   [`SurfaceCacheOptions::ang_tol_rad`].
//! - **Topology reads.** [`TrimSrfCache::build_tree_with_context`] walks the real
//!   B-rep through a borrowed [`FaceCacheContext`]. Caller handles are mapped from
//!   typed entities, never decoded from object keys. Exact UV B-spline copies are
//!   retained for subdivision to reimplant prior edge uses; context builds clear
//!   them on success and error. Public standalone implant calls retain their own
//!   supplied trim data without topology lookup. The final domain uses the owning face's
//!   surface, which need not be the generated surface being cached.
//! - **Node lists.** `SmTList::Prepend` onto a node's edge use, vertex and poly edge
//!   lists inserts at the front of the [`BezierAux2d`] `Vec`s (list order, head first).
//!   The list memory managers (`m_sVLMgr`, `m_sEdgeuseListMgr`) have no counterpart.
//! - **Temporary state.** `SmTemporaryChangeValue` (restored on every return) becomes
//!   explicit save / restore around the inner bodies; as upstream, the
//!   `m_bPointTestEnabled` set at the end of the face branch is overwritten by the
//!   restore when the branch ends.
//! - **Virtual dispatch.** `MarkFaceContainment` uses the existing base traversal
//!   with a prepared read-only face classifier and face-tolerance callback. Point
//!   tests on an unwrapped base retain trimmed dispatch and require an operation
//!   scoped [`FaceQueries`] service. A face-owned cache cannot be built without its
//!   context; opaque handles alone cannot reach topology. Poly-edge reimplantation
//!   still requires its separate poly-topology context.
//! - **Recursion.** Edge implantation uses an explicit subdivision stack. Its
//!   nested reimplantation always disables subdivision, so that call chain has
//!   at most one nested implantation. Geometric subdivision retains the base
//!   cache's documented bounds; construction uses its containment guard.
//! - **Undefined behaviour.** Upstream dereferences a null face in `PointTest` /
//!   `PointClassify` of a cache whose surface has no face; the port returns
//!   `SM_ERR_NULL_POINTER`. A node without surface data in this cache gives
//!   [`SolidError::InvalidInput`].
//! - **Not ported here.** `GetMemoryUsed` (the base class sum) and `Dump` (debug
//!   output).

use crate::curve::curve::Curve;
use crate::curve::curve_cache::{CurveCache, CurveCacheOptions};
use crate::error::{EntityRef, Result, SolidError};
use crate::geom::extent1d::Extent1d;
use crate::geom::extent2d::Extent2d;
use crate::geom::extent3d::Extent3d;
use crate::geom::tree::{NodeClass, NodeData, NodeId};
use crate::geom::vector2d::{Point2, Vec2};
use crate::geom::vector3d::{Point3, Vec3};
use crate::kernel::Kernel;
use crate::math::{max2, EFF_ZERO_SQRT};
use crate::nurbs::b_spline_curve::{BSplineCurve, BoundingBoxRequest};
use crate::nurbs::b_spline_surface::SurfParamType;
use crate::surface::surface::p2::PLANE_TYPE;
use crate::surface::surface::Surface;
use crate::surface::surface_cache::face_context::{
    FaceCacheContext, FaceQueries, PreparedFaceQueries,
};
use crate::surface::surface_cache::p1::{self as sc, FaceContainment, SurfaceCacheOptions};
use crate::surface::surface_cache::{BezierAux2d, SurfaceCache};

/// `SmTrimSrfCache_TYPE` (`SURF_BASE_TYPE + 120`).
pub const TYPE: i64 = 15120;

/// `SmTrimSrfCache::GetTypeString`.
pub const TYPE_STRING: &str = "SmTrimSrfCache_TYPE";

/// `SM_ZONE_TOL_3D`.
const ZONE_TOL_3D: f64 = crate::geom::tol::ZONE_TOL_3D;

/// `SM_SS_UMIN`.
const SS_UMIN: u64 = 1;
/// `SM_SS_VMIN`.
const SS_VMIN: u64 = 2;
/// `SM_SS_UMAX`.
const SS_UMAX: u64 = 4;
/// `SM_SS_VMAX`.
const SS_VMAX: u64 = 8;
/// `SM_SS_NONE`.
const SS_NONE: u64 = 0;

/// `SM_PC_UNKNOWN` (`SmPointClassificationType`).
pub const POINT_CLASS_UNKNOWN: i32 = 100;
/// `SM_PC_FACE`.
pub const POINT_CLASS_FACE: i32 = 102;
/// `SM_PC_EDGE`.
pub const POINT_CLASS_EDGE: i32 = 103;
/// `SM_PC_VERTEX`.
pub const POINT_CLASS_VERTEX: i32 = 104;

const OP_CLASSIFY: &str = "surface::trim_srf_cache::TrimSrfCache::point_classify";
const OP_EU: &str = "surface::trim_srf_cache::TrimSrfCache::implant_edgeuse";
const OP_VU: &str = "surface::trim_srf_cache::TrimSrfCache::implant_vertexuse";
const OP_BUILD: &str = "surface::trim_srf_cache::TrimSrfCache::build_tree_with_context";

/// `SmTrimSrfCache::IsKindOf`: the trimmed surface cache or a base type.
#[must_use]
pub const fn is_kind_of(t: i64) -> bool {
    t == TYPE || sc::is_kind_of(t)
}

fn invalid(operation: &'static str, reason: &'static str) -> SolidError {
    SolidError::InvalidInput { operation, reason }
}

/// What `ImplantVertexuse` reads from the vertex use (see the module deviations).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VertexuseInput<O> {
    /// The vertex use (stored in the node vertex lists).
    pub vertexuse: O,
    /// `SmVertexuse::ComputeUVPoint`.
    pub uv_point: Point2,
    /// `GetVertex()->GetTolerance()`.
    pub vertex_tolerance: f64,
}

/// What `ImplantPolyEdge` reads from the poly edge (see the module deviations).
#[derive(Debug, Clone, Copy)]
pub struct PolyEdgeInput<'a, O> {
    /// The poly edge (stored in the node poly edge lists).
    pub poly_edge: O,
    /// `GetStartPoint()` (its x and y are the parameter point).
    pub start: Point3,
    /// `GetEndPoint()`.
    pub end: Point3,
    /// `GetTolerance()`.
    pub tolerance: f64,
    /// `GetPolyFace()->GetOriginalFace()->GetSurface()`.
    pub surface: &'a Surface,
}

/// The trimmed surface cache (`SmTrimSrfCache`): a [`SurfaceCache`] whose tree
/// knows the face boundary.
#[derive(Debug, Clone)]
pub struct TrimSrfCache<O> {
    /// The `SmSurfaceCache` base.
    base: SurfaceCache<O>,
    /// Exact edge-use geometry for reimplantation, scoped to a context build.
    construction_edge_curves: Vec<(O, BSplineCurve)>,
}

impl<O: Copy + PartialEq> TrimSrfCache<O> {
    /// A trimmed cache of `surface` over its natural domain, not built
    /// (`SmTrimSrfCache::SmTrimSrfCache`): the state anticipates a face with trim
    /// curves (point testing on, boundary curves off) until [`TrimSrfCache::build_tree`]
    /// finds out.
    #[must_use]
    pub fn new(surface: Surface, opts: SurfaceCacheOptions) -> Self {
        let domain = surface.natural_uv_domain();
        let mut base = SurfaceCache::new(surface, domain, opts);

        // set state bits anticipating a Trim Surface cache with UVTrimCurves.
        // BuildTree() will change these values when there are no UVTrimCurves.
        base.face_was_modified = false;
        base.process_boundary_curves = false; // solvers skip naturalBoundaryCurve/target solutions
        base.point_test_enabled = true; // solvers only accept point solutions within trim boundaries
        base.have_t_surface_cache = false; // Set to true after adding UVTrimCurves to cache
        base.face_containment_done = FaceContainment::NotDone;
        base.trimmed_dispatch = true;
        Self {
            base,
            construction_edge_curves: Vec::new(),
        }
    }

    /// Wraps an existing cache as is (no state change).
    #[must_use]
    pub const fn from_base(base: SurfaceCache<O>) -> Self {
        let mut base = base;
        base.trimmed_dispatch = true;
        Self {
            base,
            construction_edge_curves: Vec::new(),
        }
    }

    /// The `SmSurfaceCache` base.
    #[must_use]
    pub const fn base(&self) -> &SurfaceCache<O> {
        &self.base
    }

    /// The `SmSurfaceCache` base, mutable.
    pub fn base_mut(&mut self) -> &mut SurfaceCache<O> {
        &mut self.base
    }

    /// Unwraps the `SmSurfaceCache` base.
    #[must_use]
    pub fn into_base(self) -> SurfaceCache<O> {
        self.base
    }

    /// `SmTrimSrfCache_TYPE` (`GetType`).
    #[must_use]
    pub const fn type_code(&self) -> i64 {
        TYPE
    }

    /// `"SmTrimSrfCache_TYPE"` (`GetTypeString`).
    #[must_use]
    pub const fn type_string(&self) -> &'static str {
        TYPE_STRING
    }

    /// [`is_kind_of`] (`IsKindOf`).
    #[must_use]
    pub const fn is_kind_of(&self, t: i64) -> bool {
        is_kind_of(t)
    }

    // ---- point tests ------------------------------------------------------------------

    /// Whether `uv_point` lies in or on the face boundary (`PointTest`): `(ok,
    /// object)`, `ok` true when point testing is off, the face was modified, or the
    /// point classifies to the face, one of its edges or vertices (`object`).
    /// `opt_zone_tol_3d` is the tolerance of the object the point comes from (else the
    /// face's).
    ///
    /// # Errors
    ///
    /// `SM_ERR_NULL_POINTER` when the surface has no face; those of the face
    /// classification; [`SolidError::Unported`] for its unported callees.
    pub fn point_test(
        &self,
        k: &Kernel,
        uv_point: Point2,
        opt_zone_tol_3d: Option<f64>,
    ) -> Result<(bool, Option<O>)> {
        self.base
            .point_test_with(k, uv_point, opt_zone_tol_3d, None)
    }

    /// Point test using an operation-scoped prepared face-query service.
    pub fn point_test_with(
        &self,
        k: &Kernel,
        uv_point: Point2,
        opt_zone_tol_3d: Option<f64>,
        queries: Option<&dyn FaceQueries<O>>,
    ) -> Result<(bool, Option<O>)> {
        self.base
            .point_test_with(k, uv_point, opt_zone_tol_3d, queries)
    }

    /// Classifies a point using the shared operation-scoped face-query service.
    pub fn point_classify_with(
        &self,
        k: &Kernel,
        uv_point: Point2,
        src_zone_tol_3d: f64,
        face: O,
        queries: &dyn FaceQueries<O>,
    ) -> Result<(i32, Option<O>)> {
        if self.base.face()? != Some(face) {
            return Err(invalid(
                OP_CLASSIFY,
                "classifier face does not match cache binding",
            ));
        }
        self.base
            .point_classify_with(k, uv_point, src_zone_tol_3d, queries)
    }

    /// Binds this cache to the caller's actual owning face without building its
    /// tree. Classification can then run while Boolean construction is delayed.
    ///
    /// # Errors
    /// A mismatched owner key, missing face surface, or failing handle mapping.
    pub fn bind_face_context(&mut self, context: &FaceCacheContext<'_, O>) -> Result<()> {
        let face_ref = EntityRef::Face(context.face);
        let expected = (context.key_of)(face_ref);
        if self.base.surface().face()? != Some(expected) {
            return Err(invalid(
                OP_BUILD,
                "surface owner does not match supplied face or B-rep",
            ));
        }
        crate::topology::face::p2::surface(context.brep, context.face)?
            .ok_or(SolidError::null_pointer(OP_BUILD))?;
        let handle = (context.handle_of)(face_ref)?;
        self.base.bind_trimmed_face(expected, handle)?;
        self.base
            .set_owner_zone_tol_3d(crate::topology::face::p2::tolerance(
                context.brep,
                context.face,
            )?);

        Ok(())
    }

    /// Builds this face's trimmed cache with a borrowed mutable B-rep context.
    ///
    /// Binding is checked against the surface's mapped owner key before construction.
    pub fn build_tree_with_context(
        &mut self,
        k: &Kernel,
        context: &mut FaceCacheContext<'_, O>,
    ) -> Result<()> {
        self.bind_face_context(context)?;
        let saved = self.base.face_containment_done;
        let fresh_tree = self.base.tree().is_none();
        self.base.face_containment_done = FaceContainment::UnderConstruction;
        let result = self.build_tree_context_body(k, context);
        self.finish_build(result, saved, fresh_tree)
    }

    /// Returns the decomposition tree, running the trimmed override with this
    /// owning-face context when construction was delayed or trims need loading.
    ///
    /// # Errors
    /// Those of [`Self::build_tree_with_context`].
    pub fn get_tree_with_context(
        &mut self,
        k: &Kernel,
        context: &mut FaceCacheContext<'_, O>,
    ) -> Result<&mut crate::geom::tree::Tree<O>> {
        if self.base.tree().is_none()
            || (!self.base.have_t_surface_cache
                && self.base.face_containment_done != FaceContainment::UnderConstruction)
        {
            self.build_tree_with_context(k, context)?;
        }
        self.base
            .tree_mut()
            .ok_or(SolidError::null_pointer(OP_BUILD))
    }

    /// The context build inside the containment guard, in upstream `BuildTree` order.
    fn build_tree_context_body(
        &mut self,
        k: &Kernel,
        context: &mut FaceCacheContext<'_, O>,
    ) -> Result<bool> {
        if self.base.tree().is_none() {
            surface_cache_build_tree(&mut self.base, k)?;
        }
        if context.brep.editing_enabled() {
            self.set_standalone_flags();
            return Ok(false);
        }

        let saved_point_test = self.base.point_test_enabled;
        let saved_editing = context.brep.editing_enabled();
        self.base.point_test_enabled = false;
        context.brep.set_editing_enabled(true);
        let result = self.implant_face_with_context(k, context);
        // No fallible operation is allowed between the body and these restores.
        context.brep.set_editing_enabled(saved_editing);
        self.base.point_test_enabled = saved_point_test;
        result
    }

    fn implant_face_with_context(
        &mut self,
        k: &Kernel,
        context: &mut FaceCacheContext<'_, O>,
    ) -> Result<bool> {
        let created = crate::topology::face::p3::create_uv_trim_curves(
            k,
            context.brep,
            context.face,
            true,
            None,
        );
        // Creation may raise tolerances even before reporting an error. Keep the
        // evaluation surface's owner tolerance in step with the actual face.
        self.base
            .set_owner_zone_tol_3d(crate::topology::face::p2::tolerance(
                context.brep,
                context.face,
            )?);
        created?;

        let loopuses =
            crate::topology::face::p4::upward_loopuses(context.brep, context.face, false, None)?;
        let mut face_domain = Extent2d::new();
        let do_subdivision = !self.base.surface().is_kind_of(PLANE_TYPE)?;

        // Preserve GetUpwardLoopuses / GetEdgeuses order and complete all edge
        // implants before computing or implanting any vertex use.
        for &lu in &loopuses {
            for eu in crate::topology::loopuse::edgeuses(context.brep, lu, None, None)? {
                let handle = (context.handle_of)(EntityRef::Edgeuse(eu))?;
                let curve = topology_uv_trim_curve(context.brep, eu)?;
                if let Some(curve) = &curve {
                    if self
                        .construction_edge_curves
                        .iter()
                        .any(|(h, _)| *h == handle)
                    {
                        return Err(invalid(OP_BUILD, "edge-use handles are not unique"));
                    }
                    self.construction_edge_curves
                        .push((handle, curve.clone_exact()));
                }
                self.implant_edgeuse(
                    k,
                    handle,
                    curve.as_ref(),
                    do_subdivision,
                    &mut face_domain,
                    None,
                )?;
            }
        }
        for &lu in &loopuses {
            for vu in crate::topology::loopuse::vertexuses(context.brep, lu, None)? {
                let vu = vu.ok_or(SolidError::null_pointer(OP_BUILD))?;
                let uv_point =
                    crate::topology::vertexuse::compute_uv_point(k, context.brep, vu, true)?;
                let vertex = crate::topology::vertexuse::vertex(context.brep, vu)?
                    .ok_or(SolidError::null_pointer(OP_BUILD))?;
                self.implant_vertexuse(
                    k,
                    VertexuseInput {
                        vertexuse: (context.handle_of)(EntityRef::Vertexuse(vu))?,
                        uv_point,
                        vertex_tolerance: crate::topology::vertex::tolerance(context.brep, vertex)?,
                    },
                    &mut face_domain,
                )?;
            }
        }

        let surface_domain = crate::topology::face::p2::surface(context.brep, context.face)?
            .ok_or(SolidError::null_pointer(OP_BUILD))?
            .natural_uv_domain();
        self.base
            .set_uv_domain(surface_domain.intersect(&face_domain)?);
        if loopuses.is_empty() {
            self.set_standalone_flags();
            return Ok(false);
        }

        self.base.process_boundary_curves = false;
        self.base.point_test_enabled = true;
        self.base.have_t_surface_cache = true;
        let queries = PreparedFaceQueries::new(context.brep, context.face, context.handle_of)?;
        let face = (context.handle_of)(EntityRef::Face(context.face))?;
        self.base.mark_face_containment_with(
            k,
            &|face| queries.face_zone_tol_3d(face),
            &|uv, src_tol| queries.face_point_classify(k, face, uv, src_tol),
        )?;
        Ok(true)
    }

    fn set_standalone_flags(&mut self) {
        self.base.process_boundary_curves = true;
        self.base.point_test_enabled = false;
        self.base.have_t_surface_cache = false;
    }

    /// Restores the containment guard and releases construction-only geometry,
    /// including on errors (no stale Bezier/subdivision references escape).
    fn finish_build(
        &mut self,
        result: Result<bool>,
        saved: FaceContainment,
        fresh_tree: bool,
    ) -> Result<()> {
        let cleanup = if fresh_tree {
            self.clean_up_build_data()
        } else {
            Ok(())
        };
        self.construction_edge_curves.clear();
        let result = result.and_then(|done| cleanup.map(|()| done));
        self.base.face_containment_done = if matches!(result, Ok(true)) {
            FaceContainment::Done
        } else {
            saved
        };
        if result.is_err() {
            self.base.have_t_surface_cache = false;
        }
        result.map(|_| ())
    }

    /// Classifies `uv_point` against the face (`PointClassify`): the
    /// `SmPointClassificationType` value and the face, edge or vertex the point is
    /// in or on (`None` for `SM_PC_UNKNOWN`). `src_zone_tol_3d` is the tolerance of
    /// the object the point comes from.
    ///
    /// # Errors
    ///
    /// `SM_ERR_NULL_POINTER` when the surface has no face; those of the face
    /// classification; [`SolidError::Unported`] for its unported callees.
    pub fn point_classify(
        &self,
        _k: &Kernel,
        _uv_point: Point2,
        _src_zone_tol_3d: f64,
    ) -> Result<(i32, Option<O>)> {
        self.base
            .face()?
            .ok_or(SolidError::null_pointer(OP_CLASSIFY))?;
        Err(SolidError::unported(
            "prepared face classifier for bound trim",
            "query.face_p2",
        ))
    }

    // ---- construction ------------------------------------------------------------------

    /// Builds a standalone surface's decomposition tree (`BuildTree`). A face-owned
    /// surface requires [`Self::build_tree_with_context`] to prepare and implant its
    /// actual trims; opaque handles alone cannot service that build.
    ///
    /// # Errors
    ///
    /// Those of tree construction; [`SolidError::Unported`] for a face-owned surface
    /// without a topology context.
    pub fn build_tree(&mut self, k: &Kernel) -> Result<()> {
        let saved = self.base.face_containment_done;
        let fresh_tree = self.base.tree().is_none();
        self.base.face_containment_done = FaceContainment::UnderConstruction;
        let result = (|| {
            if self.base.surface().face()?.is_some() {
                if self.base.tree().is_none() {
                    surface_cache_build_tree(&mut self.base, k)?;
                }
                // Preserve the original opaque-handle boundary after base-tree
                // construction; a mapped face still needs mutable build context.
                self.base.face()?;
                Err(SolidError::unported(
                    "TrimSrfCache::build_tree_with_context for a face-owned surface",
                    "surface.trim_srf_cache",
                ))
            } else {
                self.build_standalone_tree(k)
            }
        })();
        self.finish_build(result, saved, fresh_tree)
    }

    fn build_standalone_tree(&mut self, k: &Kernel) -> Result<bool> {
        if self.base.tree().is_none() {
            surface_cache_build_tree(&mut self.base, k)?;
        }
        self.set_standalone_flags();
        Ok(false)
    }

    /// Drops the construction data of a freshly built tree: the subdivision surfaces
    /// and the patches' Bezier surfaces; every node's surface becomes the cache
    /// surface.
    fn clean_up_build_data(&mut self) -> Result<()> {
        self.base.subdivision_surfaces_mut().clear();
        for p in self.base.patches_mut() {
            p.bezier = None;
        }
        let handle = self.base.surface_handle();
        if let Some(tree) = self.base.tree_mut() {
            let ids: Vec<NodeId> = (0..tree.nodes().len()).map(NodeId).collect();
            for id in ids {
                if let NodeData::Surface(a) = &mut tree.node_mut(id)?.data {
                    if a.slot.is_some() {
                        a.surface = handle;
                    }
                }
            }
        }
        Ok(())
    }

    // ---- implanting -----------------------------------------------------------------

    /// Puts `node` on the boundary (`SM_NC_ON_BOUNDARY`), marked when `mark`.
    fn set_on_boundary(&mut self, node: NodeId, mark: bool) -> Result<()> {
        let tree = self
            .base
            .tree_mut()
            .ok_or_else(|| invalid(OP_EU, "cache has no tree"))?;
        match &mut tree.node_mut(node)?.data {
            NodeData::Surface(a) => {
                a.node_class = NodeClass::OnBoundary;
                if mark {
                    a.is_marked = true;
                }
            }
            _ => return Err(invalid(OP_EU, "node has no surface data")),
        }
        Ok(())
    }

    /// `SmTreeNode::PropagateToParents` of `node`.
    fn propagate(&mut self, node: NodeId) -> Result<()> {
        self.base
            .tree_mut()
            .ok_or_else(|| invalid(OP_EU, "cache has no tree"))?
            .propagate_to_parents(node)
    }

    fn aux_mut(&mut self, node: NodeId) -> Result<&mut BezierAux2d<O>> {
        self.base.aux_mut(node)
    }

    /// Embeds the edge use `edgeuse`, whose UV trim curve is `uv_trim_curve`, into
    /// the decomposition (`ImplantEdgeuse`): every leaf near the curve (or near a pole
    /// side of the surface) is marked on the boundary and gets the edge use, the
    /// change is propagated to its ancestors; with `do_subdivision` large leaves are
    /// split first. Only `opt_tree_node` is visited when given, else every node whose
    /// domain meets the curve box. `face_extent` grows by the curve box.
    ///
    /// # Errors
    ///
    /// `SM_ERR_NULL_POINTER` without trim curve; those of the curve box, the curve
    /// cache and the tree; [`SolidError::Unported`] for `FindUVNodes` and
    /// `SubdivideNode`.
    #[allow(clippy::too_many_lines)]
    pub fn implant_edgeuse(
        &mut self,
        k: &Kernel,
        edgeuse: O,
        uv_trim_curve: Option<&BSplineCurve>,
        do_subdivision: bool,
        face_extent: &mut Extent2d,
        opt_tree_node: Option<NodeId>,
    ) -> Result<()> {
        // locals
        let eu_curve = uv_trim_curve.ok_or(SolidError::null_pointer(OP_EU))?;
        // Subdivision may need a curve implanted by an earlier call. Retain its
        // exact input data instead of trying to resolve an opaque edge-use handle.
        // Context builds discard these construction snapshots on every exit.
        if let Some((_, curve)) = self
            .construction_edge_curves
            .iter_mut()
            .find(|(handle, _)| *handle == edgeuse)
        {
            *curve = eu_curve.clone_exact();
        } else {
            self.construction_edge_curves
                .push((edgeuse, eu_curve.clone_exact()));
        }

        // set up for poles
        let (poles, u_ivl, v_ivl, u_tol, v_tol) = {
            let uv_domain = self.base.surface().natural_uv_domain();

            // save natural domain sizes
            let u_ivl = uv_domain.u_interval();
            let v_ivl = uv_domain.v_interval();
            let u_tol = u_ivl.length() / 100.0;
            let v_tol = v_ivl.length() / 100.0;

            // list m_sUVDomain boundaries near SrfNaturalDomain boundaries
            let d = self.base.uv_domain();
            let mut poles = u64::from((d.min().x - u_ivl.min()).abs() < u_tol)
                + (if (d.min().y - v_ivl.min()).abs() < v_tol {
                    2
                } else {
                    0
                })
                + (if (d.max().x - u_ivl.max()).abs() < u_tol {
                    4
                } else {
                    0
                })
                + (if (d.max().y - v_ivl.max()).abs() < v_tol {
                    8
                } else {
                    0
                });

            // eliminate boundaries which are not poles
            if poles != 0 {
                let srf_poles = surface_singularities(k, self.base.surface())?;
                poles &= srf_poles;
            }
            (poles, u_ivl, v_ivl, u_tol, v_tol)
        };

        // sCurveBBox = EUCurve 2D Bounding Box expanded by tolerance
        let bbox = eu_curve
            .calculate_bounding_box(
                &eu_curve.natural_interval(),
                BoundingBoxRequest {
                    normal: true,
                    ..BoundingBoxRequest::default()
                },
            )?
            .normal
            .ok_or_else(|| invalid(OP_EU, "no curve box"))?;
        let mut curve_bbox =
            Extent2d::from_bounds(bbox.min().x, bbox.min().y, bbox.max().x, bbox.max().y);
        let expand = (1.0 + curve_bbox.size().length()) * EFF_ZERO_SQRT;
        curve_bbox.expand_absolute(expand);

        // when EUCurve is a line - get its local values
        let line = eu_curve
            .is_line(4, ZONE_TOL_3D / 10.0)?
            .map(|(p, v)| (Point2::new(p.x, p.y), Vec2::new(v.x, v.y)));

        // Compute OffAxisTolerance - maximum size of non largest component
        //                            of the bounding box of the curve segment.
        let mut diff = self.base.uv_domain().size();
        let mut scale = 1.0;

        let mut clamp_node_growth = false;
        if diff.x > diff.y * 20.0 {
            clamp_node_growth = true;
            scale = diff.x / diff.y;
        }
        if diff.y > diff.x * 20.0 {
            clamp_node_growth = true;
            scale = diff.x / diff.y;
        }
        // Scale things to smaller size if one dimension of UV domain is
        // very large compared to the other.
        if diff.x > diff.y * 100.0 {
            diff.x /= 100.0;
        } else if diff.x > diff.y * 10.0 {
            diff.x = diff.y * 10.0;
        }

        if diff.y > diff.x * 100.0 {
            diff.y /= 100.0;
        } else if diff.y > diff.x * 10.0 {
            diff.y = diff.x * 10.0;
        }

        let uv_size = diff.length();
        let off_axis_tolerance = uv_size / self.base.trim_subdivision_factor;
        let mut min_node_size_squared = off_axis_tolerance * off_axis_tolerance;

        // increase dMinNodeSizeSquared when working on very long lines
        if line.is_some() {
            let line_size = curve_bbox.size().length() / 10.0;
            min_node_size_squared = max2(min_node_size_squared, line_size * line_size);
        }

        // sMinSize
        let mut min_size = Vec2::new(0.0, 0.0);
        if clamp_node_growth {
            let diff2 = self.base.uv_domain().size();
            min_size = diff2 / self.base.trim_subdivision_factor;
        }

        // find all tree nodes that intersect sCurveBBox - use pOptTreeNode if given
        let nodes = match opt_tree_node {
            None => find_uv_nodes(&mut self.base, k, &curve_bbox, true)?,
            Some(n) => vec![n],
        };

        // no work - no nodes
        if nodes.is_empty() {
            return Ok(());
        }

        // Create the curve cache for nonLinear EUCurves
        let ivl = eu_curve.natural_interval();
        let curve_cache = if line.is_none() {
            let mut cc: CurveCache<O> = CurveCache::new(
                Curve::from(eu_curve.clone_exact()),
                ivl,
                CurveCacheOptions {
                    make_full_tree: true,
                    chord_height_tolerance: 0.0,
                    angle_tol_deg: 20.0,
                    off_axis_tolerance,
                    min_segments: 0,
                    compute_aux_tree_data: true,
                    ..CurveCacheOptions::default()
                },
            );
            cc.tessellate(k)?;
            Some(cc)
        } else {
            None
        };

        // for every node that intersects the EUCurve bounding box
        let mut stack: Vec<NodeId> = Vec::new();
        for &start in &nodes {
            // init stack with current node
            stack.push(start);

            // while nodes exist in the stack
            while let Some(node) = stack.pop() {
                // skip non-leaf nodes
                let tree = self
                    .base
                    .tree()
                    .ok_or_else(|| invalid(OP_EU, "cache has no tree"))?;
                if tree.node(node)?.child1.is_some() {
                    continue;
                }

                // sNodeBBox = this node's 2D bounding box
                let node_dom = self.base.aux(node)?.uv_domain;
                let node_bbox =
                    Extent3d::from_points(Vec3::from(node_dom.min()), Vec3::from(node_dom.max()));

                // nodes by position
                let keep_for_position = match (&curve_cache, line) {
                    // when working with nonlinear curves - keep nodes that intersect the curve cache
                    (Some(cc), _) => cc
                        .tree()
                        .ok_or_else(|| invalid(OP_EU, "curve cache has no tree"))?
                        .intersects_box(&node_bbox)?,
                    // when working with lines - skip nodes too far from line
                    (None, Some((line_point, line_vec))) => {
                        // First do bounding box check (uv).
                        let mut keep =
                            !curve_bbox.are_disjoint(&Extent2d::from_extent3d(&node_bbox), 0.0);
                        // BBoxes intersect, check against line.
                        if keep {
                            keep = node_dom
                                .intersect_with_infinite_line(line_point, line_vec, None, None)
                                .is_some();
                        }
                        keep
                    }
                    (None, None) => false,
                };

                // also keep nodes near poles
                let mut keep_for_pole = false;
                if poles != 0 && !keep_for_position {
                    keep_for_pole = ((poles & SS_UMIN) != 0
                        && (node_dom.min().x - u_ivl.min()).abs() < u_tol)
                        || ((poles & SS_UMAX) != 0
                            && (node_dom.max().x - u_ivl.max()).abs() < u_tol)
                        || ((poles & SS_VMIN) != 0
                            && (node_dom.min().y - v_ivl.min()).abs() < v_tol)
                        || ((poles & SS_VMAX) != 0
                            && (node_dom.max().y - v_ivl.max()).abs() < v_tol);
                }

                // skip nodes too far away from this node's bounding box
                if !keep_for_position && !keep_for_pole {
                    continue;
                }

                // When not doing subdividing, when node is small enough
                let uv_diff = node_dom.size();
                if !do_subdivision
                    || uv_diff.length_squared() < min_node_size_squared
                    || uv_diff.x < min_size.x
                    || uv_diff.y < min_size.y
                {
                    // mark the AuxData as on the boundary, add the Edgeuse to the
                    // AuxData's curve list, propagate state change to ancestors
                    self.set_on_boundary(node, true)?;
                    self.aux_mut(node)?.edgeuse_list.insert(0, edgeuse);
                    self.propagate(node)?;
                } else {
                    // subdivide the node's largest UV direction
                    let dir = if uv_diff.x > uv_diff.y * scale {
                        SurfParamType::U
                    } else {
                        SurfParamType::V
                    };

                    // promote leaf node into a parent->Child(1,2) triplet
                    let (c1, c2) = self.subdivide_node(k, node, dir)?;

                    // put the children on the stack
                    stack.push(c1);
                    stack.push(c2);
                }
            }
        }

        // Update the domain of the face extents
        if face_extent.is_init() || face_extent.has_negative_area() {
            *face_extent = curve_bbox;
        } else {
            *face_extent = curve_bbox.union(face_extent);
        }

        // all done
        Ok(())
    }

    /// `SmSurfaceCache::SubdivideNode(node, dir, TRUE, NULL)` on this cache: the two new
    /// children. A parent on the boundary moves its edge uses into them through the
    /// trimmed surface cache's `ImplantEdgeuse(pEU, FALSE, sTempExt, pNode)` (virtual
    /// upstream), last first, and its poly edges through `ImplantPolyEdge`.
    ///
    /// # Errors
    ///
    /// Those of [`SurfaceCache::subdivide_node_split`]; [`SolidError::Unported`]:
    /// missing operation-scoped trim data for an edge use to move,
    /// `SmFace::GetSurface` (`topology.face_p2`) for a poly edge to move.
    fn subdivide_node(
        &mut self,
        k: &Kernel,
        node: NodeId,
        dir: SurfParamType,
    ) -> Result<(NodeId, NodeId)> {
        const OP: &str = "surface::trim_srf_cache::TrimSrfCache::subdivide_node";
        if let Some((node1, node2)) = self.base.subdivide_node_split(k, node, dir, true)? {
            if self.base.aux(node)?.poly_edge_list.is_empty() {
                // for every edge use in the parent node, last first
                while let Some(eu) = self.base.aux(node)?.edgeuse_list.last().copied() {
                    // implant the edge use in both children without subsequent
                    // subdivision (the status is ignored)
                    let curve = self
                        .construction_edge_curves
                        .iter()
                        .find(|(handle, _)| *handle == eu)
                        .map(|(_, curve)| curve.clone_exact())
                        .ok_or_else(|| {
                            SolidError::unported(
                                "operation-scoped edge-use trim for subdivision",
                                "surface.trim_srf_cache",
                            )
                        })?;
                    for child in [node1, node2] {
                        let mut temp_ext = Extent2d::new();
                        let r = self.implant_edgeuse(
                            k,
                            eu,
                            Some(&curve),
                            false,
                            &mut temp_ext,
                            Some(child),
                        );
                        if let Err(e) = r {
                            if e.is_unported() {
                                return Err(e);
                            }
                        }
                    }
                    self.base.aux_mut(node)?.edgeuse_list.pop();
                }
            } else {
                // UNPORTED(topology.face_p2): SmFace::GetSurface of the poly edge's
                // original face (the PolyEdgeInput of ImplantPolyEdge; the poly edge
                // points and tolerance are query::poly, but the cache keeps opaque
                // handles with no poly B-rep to read them from)
                return Err(SolidError::unported(
                    "SmFace::GetSurface",
                    "topology.face_p2",
                ));
            }
        }
        let n = self
            .base
            .tree()
            .ok_or(SolidError::InvalidInput {
                operation: OP,
                reason: "cache has no tree",
            })?
            .node(node)?;
        match (n.child1, n.child2) {
            (Some(c1), Some(c2)) => Ok((c1, c2)),
            _ => Err(SolidError::InvalidInput {
                operation: OP,
                reason: "split node has no children",
            }),
        }
    }

    /// Embeds a poly edge segment into `node` (`ImplantPolyEdge`): when the node's
    /// domain, enlarged by ten times the parametric tolerances, holds an end point,
    /// or the segment's line crosses the domain within the segment, the node is put
    /// on the boundary (neither marked nor propagated) and gets the poly edge.
    ///
    /// # Errors
    ///
    /// `InvalidInput` for a node without surface data; [`SolidError::Unported`]:
    /// `SmTol::MapTo2d` over a surface (`surface.surface_p1`).
    pub fn implant_poly_edge(
        &mut self,
        k: &Kernel,
        input: PolyEdgeInput<'_, O>,
        node: NodeId,
    ) -> Result<()> {
        // get segment endpoints
        let start = Point2::new(input.start.x, input.start.y);
        let end = Point2::new(input.end.x, input.end.y);

        // Tolerance: get uv space tols, expand the uv domain of the node by tols
        let tol_3d = input.tolerance;
        let mid_uv = (start + end) / 2.0;
        let tol_u = map_to_2d(k, tol_3d, mid_uv, Vec2::new(1.0, 0.0), input.surface)?;
        let tol_v = map_to_2d(k, tol_3d, mid_uv, Vec2::new(0.0, 1.0), input.surface)?;
        let tol_expansion = Vec2::new(tol_u, tol_v);

        let domain = self.base.aux(node)?.uv_domain;
        let mut domain_exp = domain;
        domain_exp.expand_absolute_vec(tol_expansion * 10.0);

        // locals for intersection test with node domain
        let line_vec = end - start;
        let is_degenerate = line_vec.x.abs() <= tol_u && line_vec.y.abs() <= tol_v;
        let unit = Extent1d::from_min_max(0.0, 1.0);

        // testing to see if node contains the start/end points of PolyEdge or
        // PolyEdge intersects the domain
        let mut embed = true;
        if !domain_exp.contains_point(start, 0.0) && !domain_exp.contains_point(end, 0.0) {
            if is_degenerate {
                embed = false;
            } else {
                match domain.intersect_with_infinite_line(start, line_vec, None, None) {
                    None => embed = false,
                    Some(ivl) => {
                        if ivl.intersect(&unit).is_err() {
                            embed = false;
                        }
                    }
                }
            }
        }

        if embed {
            // Add segment to this node's m_sPolyEdgeList
            self.set_on_boundary(node, false)?;
            self.aux_mut(node)?
                .poly_edge_list
                .insert(0, input.poly_edge);
        }
        Ok(())
    }

    /// Embeds a vertex use into the decomposition (`ImplantVertexuse`): every leaf
    /// meeting the point's parametric tolerance box (widened to a whole pole side) is
    /// marked on the boundary and gets the vertex use; `face_extent` grows by the box.
    ///
    /// # Errors
    ///
    /// Those of the tree, of `Evaluate1stDerivatives` and of
    /// [`SurfaceCache::find_uv_nodes`] and `GetSingularities`.
    pub fn implant_vertexuse(
        &mut self,
        k: &Kernel,
        input: VertexuseInput<O>,
        face_extent: &mut Extent2d,
    ) -> Result<()> {
        let _ = k;
        // locals: UVPoint, Vertex, 3D Surface Properties
        let uv = input.uv_point;
        let (_pnt, du, dv) = evaluate_1st_derivatives(k, self.base.surface(), uv)?;

        // get surface u and v speeds
        let mut du_len = du.length();
        let mut dv_len = dv.length();

        // keep the parameter space tolerance from getting too long and skinny
        if du_len < dv_len / 10.0 {
            du_len = dv_len;
        }
        if dv_len < du_len / 10.0 {
            dv_len = du_len;
        }

        // pick u and v tolerance
        let mut uv_tol = Vec2::new(
            input.vertex_tolerance / du_len,
            input.vertex_tolerance / dv_len,
        );

        // expand the UV tolerances when near surface poles
        let uv_domain = self.base.surface().natural_uv_domain();
        let u_ivl = uv_domain.u_interval();
        let v_ivl = uv_domain.v_interval();
        let u_tol = u_ivl.length() / 100.0;
        let v_tol = v_ivl.length() / 100.0;

        // Get boundaries close to sUV
        let mut poles = u64::from((uv.x - u_ivl.min()).abs() < u_tol)
            + (if (uv.y - v_ivl.min()).abs() < v_tol {
                2
            } else {
                0
            })
            + (if (uv.x - u_ivl.max()).abs() < u_tol {
                4
            } else {
                0
            })
            + (if (uv.y - v_ivl.max()).abs() < v_tol {
                8
            } else {
                0
            });

        // when sUV is close enough to boundaries for a pole check
        if poles != SS_NONE {
            // only check those boundaries which are poles
            poles &= surface_singularities(k, self.base.surface())?;

            // when UV point is in constant U singular side
            if (poles & SS_UMIN) != 0 || (poles & SS_UMAX) != 0 {
                uv_tol.y = v_ivl.length();
            }

            // when UV point is in constant V singular side
            if (poles & SS_VMIN) != 0 || (poles & SS_VMAX) != 0 {
                uv_tol.x = u_ivl.length();
            }
        }

        // set point extent to include the point and its tolerance
        let uv_ext = Extent2d::from_points(uv - uv_tol, uv + uv_tol);

        // get all nodes that intersect the point extent
        let nodes = find_uv_nodes(&mut self.base, k, &uv_ext, false)?;

        // for every intersecting node
        for node in nodes {
            // skip non-leaf nodes
            let tree = self
                .base
                .tree()
                .ok_or_else(|| invalid(OP_VU, "cache has no tree"))?;
            if tree.node(node)?.child1.is_some() {
                continue;
            }

            // mark node as on boundary, add vertexuse to node's vertex list,
            // propagate state change to ancestors
            self.set_on_boundary(node, true)?;
            self.aux_mut(node)?.vertex_list.insert(
                0,
                sc::VertexListItem {
                    vertexuse: Some(input.vertexuse),
                    uv_point: uv,
                    uv_tol,
                },
            );
            self.propagate(node)?;
        }

        // Update the domain of the face extents
        let size = face_extent.size();
        if size.x < 0.0 || size.y < 0.0 {
            *face_extent = uv_ext;
        } else {
            *face_extent = uv_ext.union(face_extent);
        }

        // all done
        Ok(())
    }
}

// ---- shared construction and geometry operations ------------------------------------

/// The actual `GetUVTrimCurve` B-spline part, without copy-constructor knot fix-up.
fn topology_uv_trim_curve(
    brep: &crate::topology::Brep,
    eu: crate::topology::EdgeuseId,
) -> Result<Option<BSplineCurve>> {
    use crate::curve::curve::CurveKind;

    let Some(curve) = crate::topology::edgeuse::p1::uv_trim_curve(brep, eu)?.0 else {
        return Ok(None);
    };
    let spline = match curve.kind() {
        CurveKind::BSpline(c) => c.clone_exact(),
        CurveKind::Line(c) => c.base_with_nurb()?.clone_exact(),
        CurveKind::Hyperbola(c) => c.base().clone_exact(),
        CurveKind::Ellipse(c) => c.base().clone_exact(),
        CurveKind::Parabola(c) => c.base().clone_exact(),
        _ => return Err(SolidError::null_pointer(OP_BUILD)),
    };
    Ok(Some(spline))
}

/// `SmSurfaceCache::BuildTree(&sBezierBlock)` (base class tree construction for a derived
/// cache, which finishes it).
fn surface_cache_build_tree<O: Copy + PartialEq>(
    cache: &mut SurfaceCache<O>,
    k: &Kernel,
) -> Result<()> {
    cache.build_tree(k, true)
}

/// `SmSurfaceCache::FindUVNodes`.
fn find_uv_nodes<O: Copy + PartialEq>(
    cache: &mut SurfaceCache<O>,
    k: &Kernel,
    uv_box: &Extent2d,
    near_poles: bool,
) -> Result<Vec<NodeId>> {
    cache.find_uv_nodes(k, uv_box, near_poles)
}

/// `SmSurface::GetSingularities`: the `SM_SS_*` sides that are poles.
fn surface_singularities(k: &Kernel, surface: &Surface) -> Result<u64> {
    Ok(u64::from(
        surface.get_singularities(k, false, false)?.singular,
    ))
}

/// `SmSurface::Evaluate1stDerivatives(uv, TRUE, TRUE, ..)`: point, `dS/du`, `dS/dv`.
///
/// When wiring, keep upstream's flag arguments exactly: `(uv, TRUE, TRUE, ..)`
/// (`SmTrimSrfCache.cpp:856`).
fn evaluate_1st_derivatives(
    k: &Kernel,
    surface: &Surface,
    uv: Point2,
) -> Result<(Point3, Vec3, Vec3)> {
    crate::surface::surface::p1::evaluate_1st_derivatives(
        k,
        surface,
        uv,
        crate::nurbs::knots::Side::Left,
        crate::nurbs::knots::Side::Left,
    )
}

/// `SmTol::MapTo2d(tol, uv, dir, surface)`.
fn map_to_2d(k: &Kernel, tol_3d: f64, uv: Point2, dir: Vec2, surface: &Surface) -> Result<f64> {
    crate::surface::surface::p1::map_to_2d_dir(k, tol_3d, uv, dir, surface)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::tree::{
        AuxDataType, NodeGeomType, SplitDir, SurfaceAux, Tree, TreeNode, TreeOwner,
    };
    use crate::nurbs::b_spline_surface::BSplineSurface;
    use crate::nurbs::crv_basic::Curve as NlCurve;
    use crate::nurbs::math::{CPoint, NOW, NOZ};
    use crate::nurbs::srf_basic::Surface as NlSurface;

    fn plane() -> Surface {
        let net = (0..2)
            .map(|i| {
                (0..2)
                    .map(|j| CPoint::new(f64::from(i), f64::from(j), 0.0, NOW))
                    .collect()
            })
            .collect();
        let s = NlSurface::new(
            net,
            1,
            1,
            vec![0.0, 0.0, 1.0, 1.0],
            vec![0.0, 0.0, 1.0, 1.0],
        )
        .expect("surface");
        Surface::from(BSplineSurface::from_nurb(s, false))
    }

    fn line_2d(a: [f64; 2], b: [f64; 2]) -> BSplineCurve {
        let pw = vec![
            CPoint::new(a[0], a[1], NOZ, NOW),
            CPoint::new(b[0], b[1], NOZ, NOW),
        ];
        BSplineCurve::from_nurb(
            NlCurve::new(pw, 1, vec![0.0, 0.0, 1.0, 1.0]).expect("curve"),
            2,
            false,
        )
    }

    /// A cache over [0.2, 0.8]^2 whose tree is the root and two halves split in u.
    fn cache() -> TrimSrfCache<i64> {
        let mut c = TrimSrfCache::new(plane(), SurfaceCacheOptions::default());
        c.base_mut()
            .set_uv_domain(Extent2d::from_bounds(0.2, 0.2, 0.8, 0.8));
        let doms = [
            Extent2d::from_bounds(0.2, 0.2, 0.8, 0.8),
            Extent2d::from_bounds(0.2, 0.2, 0.5, 0.8),
            Extent2d::from_bounds(0.5, 0.2, 0.8, 0.8),
        ];
        let links = [
            (None, Some(1), Some(2)),
            (Some(0), None, None),
            (Some(0), None, None),
        ];
        let mut tree: Tree<i64> = Tree::new(Some(TreeOwner::SurfaceCache));
        let mut auxs = Vec::new();
        for (k, (d, (p, c1, c2))) in doms.iter().zip(links).enumerate() {
            auxs.push(BezierAux2d {
                uv_domain: *d,
                owning_tree_node: Some(NodeId(k)),
                ..BezierAux2d::default()
            });
            tree.add_node(TreeNode {
                parent: p.map(NodeId),
                child1: c1.map(NodeId),
                child2: c2.map(NodeId),
                aux_type: AuxDataType::AuxData,
                data: NodeData::Surface(SurfaceAux {
                    split_dir: SplitDir::U,
                    is_marked: false,
                    node_class: NodeClass::Unknown,
                    face: None,
                    surface: None,
                    slot: Some(k),
                }),
                geom_type: NodeGeomType::Default,
                ..TreeNode::new()
            });
        }
        tree.set_top_node(Some(NodeId(0)));
        *c.base_mut().auxs_mut() = auxs;
        c.base_mut().set_tree(Some(tree));
        c
    }

    #[test]
    fn constructor_state() {
        let c: TrimSrfCache<i64> = TrimSrfCache::new(plane(), SurfaceCacheOptions::default());
        assert!(c.base().point_test_enabled);
        assert!(!c.base().process_boundary_curves);
        assert_eq!(c.type_code(), TYPE);
        assert!(c.is_kind_of(TYPE) && c.is_kind_of(sc::TYPE));
        assert!(!c.is_kind_of(16002));
    }

    #[test]
    fn point_test_off() {
        let k = Kernel::default();
        let mut c: TrimSrfCache<i64> = TrimSrfCache::new(plane(), SurfaceCacheOptions::default());
        c.base_mut().point_test_enabled = false;
        assert_eq!(
            c.point_test(&k, Point2::new(5.0, 5.0), None).unwrap(),
            (true, None)
        );
        c.base_mut().point_test_enabled = true;
        // no face to classify against (SM_ERR_NULL_POINTER)
        assert!(matches!(
            c.point_test(&k, Point2::new(5.0, 5.0), None),
            Err(SolidError::Upstream { code: 1006, .. })
        ));
    }

    #[test]
    fn implant_line_in_leaf() {
        let k = Kernel::default();
        let mut c = cache();
        let mut ext = Extent2d::new();
        let crv = line_2d([0.3, 0.3], [0.4, 0.7]);
        c.implant_edgeuse(&k, 7, Some(&crv), false, &mut ext, Some(NodeId(1)))
            .unwrap();
        assert_eq!(c.base().aux(NodeId(1)).unwrap().edgeuse_list, vec![7]);
        // the other leaf is not visited, the root is not a leaf
        c.implant_edgeuse(&k, 8, Some(&crv), false, &mut ext, Some(NodeId(2)))
            .unwrap();
        assert!(c.base().aux(NodeId(2)).unwrap().edgeuse_list.is_empty());
        assert!(ext.min().x < 0.3 && ext.max().y > 0.7);
        // without trim curve
        assert!(c
            .implant_edgeuse(&k, 9, None, false, &mut ext, Some(NodeId(1)))
            .is_err());
    }

    #[test]
    fn invalid_node_returns_error_without_changing_boundary() {
        let k = Kernel::default();
        let mut c = cache();
        let crv = line_2d([0.3, 0.3], [0.4, 0.7]);
        let mut ext = Extent2d::new();
        assert!(c
            .implant_edgeuse(&k, 7, Some(&crv), false, &mut ext, Some(NodeId(usize::MAX)))
            .is_err());
        assert!(c.base().aux(NodeId(1)).unwrap().edgeuse_list.is_empty());
        assert!(ext.is_init());
    }

    #[test]
    fn implant_vertexuse_away_from_poles_marks_the_leaf() {
        let k = Kernel::default();
        let mut c = cache();
        let mut ext = Extent2d::new();
        let vu = VertexuseInput {
            vertexuse: 7,
            uv_point: Point2::new(0.4, 0.4),
            vertex_tolerance: 0.001,
        };
        // away from the domain sides no pole check (GetSingularities) is made
        c.implant_vertexuse(&k, vu, &mut ext).unwrap();
        assert_eq!(ext, Extent2d::from_bounds(0.399, 0.399, 0.401, 0.401));
        let leaf = c.base().aux(NodeId(1)).unwrap();
        assert_eq!(leaf.vertex_list.len(), 1);
        assert_eq!(leaf.vertex_list[0].vertexuse, Some(7));
        assert!(c.base().aux(NodeId(2)).unwrap().vertex_list.is_empty());
    }

    #[test]
    fn implant_poly_edge_marks_the_leaf() {
        let k = &Kernel::default();
        let mut c = cache();
        let surface = plane();
        let edge = PolyEdgeInput {
            poly_edge: 8,
            start: Point3::new(0.3, 0.3, 0.0),
            end: Point3::new(0.4, 0.7, 0.0),
            tolerance: 0.001,
            surface: &surface,
        };
        // (SmTol::MapTo2d is ported with surface.surface_p1)
        c.implant_poly_edge(k, edge, NodeId(1)).unwrap();
        let aux = c.base().aux(NodeId(1)).unwrap();
        assert_eq!((aux.poly_edge_list.len(), aux.vertex_list.len()), (1, 0));
    }

    #[test]
    fn build_tree_without_face_restores_state() {
        let k = Kernel::default();
        let mut c: TrimSrfCache<i64> = TrimSrfCache::new(plane(), SurfaceCacheOptions::default());
        c.build_tree(&k).unwrap();
        assert_eq!(c.base().face_containment_done, FaceContainment::NotDone);
    }

    fn face_key(entity: EntityRef) -> crate::geom::object::ObjectRef {
        let EntityRef::Face(face) = entity else {
            panic!("test only maps faces");
        };
        crate::geom::object::ObjectRef {
            type_code: crate::topology::face::p2::TYPE,
            key: crate::geom::attribute::ObjectKey(10_000 + u64::from(face.index())),
        }
    }

    fn face_handle(entity: EntityRef) -> Result<i64> {
        let EntityRef::Face(face) = entity else {
            return Err(invalid(OP_BUILD, "test only maps faces"));
        };
        Ok(i64::from(face.index()) + 100)
    }

    fn bare_owned_face() -> (crate::topology::Brep, crate::topology::FaceId) {
        use crate::topology::{Brep, Face, ShellId};
        let mut brep = Brep::new();
        let surface = plane();
        let mut face = Face::new(ShellId::from_index(0));
        face.uv_domain = surface.natural_uv_domain();
        face.surface = Some(surface);
        face.tolerance = 2.5e-4;
        let face = brep.add_face_record(face).unwrap();
        brep.face_mut(face)
            .unwrap()
            .surface
            .as_mut()
            .unwrap()
            .set_owner(Some(face_key(EntityRef::Face(face))));
        brep.sync_face_surface_tolerance(face).unwrap();
        (brep, face)
    }

    #[derive(Debug)]
    struct Outside;

    impl FaceQueries<i64> for Outside {
        fn face_zone_tol_3d(&self, _face: i64) -> Result<f64> {
            Ok(2.5e-4)
        }

        fn face_point_classify(
            &self,
            _k: &Kernel,
            _face: i64,
            _uv: Point2,
            _tol: f64,
        ) -> Result<(i32, Option<i64>)> {
            Ok((POINT_CLASS_UNKNOWN, None))
        }
    }

    #[test]
    fn bound_trim_preserves_dispatch_mapping_and_face_tolerance() {
        let k = Kernel::default();
        let (mut brep, face) = bare_owned_face();
        brep.set_editing_enabled(true);
        let surface = brep
            .face(face)
            .unwrap()
            .surface
            .as_ref()
            .unwrap()
            .clone_for_evaluation()
            .unwrap();
        let mut cache = TrimSrfCache::new(surface, SurfaceCacheOptions::default());
        cache
            .build_tree_with_context(
                &k,
                &mut FaceCacheContext {
                    brep: &mut brep,
                    face,
                    key_of: &face_key,
                    handle_of: &face_handle,
                },
            )
            .unwrap();
        assert!(brep.editing_enabled());
        assert_eq!(cache.base().face().unwrap(), Some(100));
        assert_eq!(cache.base().surface().owner_zone_tol_3d(), Some(2.5e-4));

        let mut base = cache.into_base();
        base.point_test_enabled = true;
        let uv = Point2::new(0.5, 0.5);
        assert!(base.uv_domain().contains_point(uv, 0.0));
        assert!(base.point_test(uv).unwrap_err().is_unported());
        assert!(base
            .point_test_with(&k, uv, None, None)
            .unwrap_err()
            .is_unported());
        assert_eq!(
            base.point_test_with(&k, uv, None, Some(&Outside)).unwrap(),
            (false, None)
        );

        let mapped = base.map_handles(|handle| Ok(handle + 1_000)).unwrap();
        assert!(mapped.is_trimmed());
        assert_eq!(mapped.face().unwrap(), Some(1_100));
        assert_eq!(mapped.surface().owner(), base.surface().owner());
        assert_eq!(mapped.surface().owner_zone_tol_3d(), Some(2.5e-4));
        assert!(base
            .map_handles::<i64>(|_| Err(invalid(OP_BUILD, "unmappable")))
            .is_err());

        base.face_was_modified = true;
        assert_eq!(
            base.point_test_with(&k, uv, None, None).unwrap(),
            (true, None)
        );
        base.face_was_modified = false;
        base.point_test_enabled = false;
        assert_eq!(
            base.point_test_with(&k, uv, None, None).unwrap(),
            (true, None)
        );
    }

    #[test]
    fn context_build_restores_guards_when_real_topology_is_incomplete() {
        let k = Kernel::default();
        let (mut brep, face) = bare_owned_face();
        let surface = brep
            .face(face)
            .unwrap()
            .surface
            .as_ref()
            .unwrap()
            .clone_for_evaluation()
            .unwrap();
        let mut cache = TrimSrfCache::new(surface, SurfaceCacheOptions::default());
        assert!(cache
            .build_tree_with_context(
                &k,
                &mut FaceCacheContext {
                    brep: &mut brep,
                    face,
                    key_of: &face_key,
                    handle_of: &face_handle,
                }
            )
            .is_err());
        assert!(!brep.editing_enabled());
        assert!(cache.base().point_test_enabled);
        assert_eq!(cache.base().face_containment_done, FaceContainment::NotDone);
        assert!(!cache.base().have_t_surface_cache);
        assert!(cache.construction_edge_curves.is_empty());
    }

    #[test]
    fn context_rejects_a_different_owner_key_before_mutating_state() {
        let k = Kernel::default();
        let (mut brep, face) = bare_owned_face();
        let surface = brep
            .face(face)
            .unwrap()
            .surface
            .as_ref()
            .unwrap()
            .clone_for_evaluation()
            .unwrap();
        let mut cache = TrimSrfCache::new(surface, SurfaceCacheOptions::default());
        let wrong_key = |entity| {
            let mut key = face_key(entity);
            key.key.0 += 1;
            key
        };
        assert!(matches!(
            cache.build_tree_with_context(
                &k,
                &mut FaceCacheContext {
                    brep: &mut brep,
                    face,
                    key_of: &wrong_key,
                    handle_of: &face_handle,
                }
            ),
            Err(SolidError::InvalidInput { .. })
        ));
        assert!(!brep.editing_enabled());
        assert!(cache.base().tree().is_none());
        assert_eq!(cache.base().face_containment_done, FaceContainment::NotDone);
    }

    #[test]
    fn standalone_trimmed_cache_can_build_lazily_through_its_base() {
        let k = Kernel::default();
        let cache: TrimSrfCache<i64> = TrimSrfCache::new(plane(), SurfaceCacheOptions::default());
        let mut base = cache.into_base();
        assert!(base.tree().is_none());
        base.get_tree(&k).unwrap();
        assert!(base.is_trimmed());
        assert!(base.process_boundary_curves);
        assert!(!base.point_test_enabled);
        assert!(!base.have_t_surface_cache);
        assert_eq!(base.face_containment_done, FaceContainment::NotDone);
    }
}
