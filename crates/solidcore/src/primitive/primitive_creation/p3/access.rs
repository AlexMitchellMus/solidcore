//! External topology services of primitive construction.
//!
//! Upstream: `source/SMLib/src/SmPrimitiveCreation.cpp` (callee boundaries).
//!
//! | Rust | Upstream |
//! |------|----------|
//! | [`PrimitiveAccess`] | B-rep, face, classification and stitching services |
//!
//! Deviations: extends the existing statically dispatched sweep access interface.
//! Defaults defer only unavailable callees; partial lists remain caller-owned.

use crate::curve::composite_curve::CompositeCurve;
use crate::curve::curve::Curve;
use crate::error::{Result, SolidError};
use crate::geom::axis2_placement::Axis2Placement;
use crate::geom::extent2d::Extent2d;
use crate::geom::extent3d::Extent3d;
use crate::surface::surface::Surface;
use crate::sweep::topology_sweep::access::{FaceRef, RegionRef, SweepAccess};
use crate::sweep::topology_sweep::TopologySweep;
use crate::topology::{BrepId, EdgeId, FaceId};

/// Topology adapters used by primitive construction. Contexts, record ids and
/// resolver identities are shared with `SweepAccess`; construction tails never
/// replace upstream's topological edits with arena-only approximations.
pub trait PrimitiveAccess: SweepAccess {
    /// Resolve the caller's Kernel and mutable B-rep together, without copying
    /// either context or topology state.
    fn with_primitive_brep<T>(
        &mut self,
        brep: BrepId,
        body: impl FnOnce(&mut crate::kernel::Kernel, &mut crate::topology::Brep) -> Result<T>,
    ) -> Result<T>
    where
        Self: Sized;

    /// Create planar faces, retaining partial outputs and consumed curve slots.
    fn planar_faces(
        &mut self,
        region: RegionRef,
        curves: &mut Vec<Option<Curve>>,
        tolerance: f64,
        faces: &mut Vec<FaceId>,
    ) -> Result<()> {
        let _ = (region, curves, tolerance, faces);
        // UNPORTED(topology.brep_p8): SmBrep::CreatePlanarFacesWith3DCurves
        Err(SolidError::unported(
            "SmBrep::CreatePlanarFacesWith3DCurves",
            "topology.brep_p8",
        ))
    }
    /// Create one planar face from its boundary curves.
    fn planar_face(
        &mut self,
        region: RegionRef,
        curves: Vec<Curve>,
        tolerance: f64,
        face: &mut Option<FaceId>,
    ) -> Result<()> {
        let _ = (region, curves, tolerance, face);
        // UNPORTED(topology.brep_p8): SmBrep::CreatePlanarFaceWith3DCurves
        Err(SolidError::unported(
            "SmBrep::CreatePlanarFaceWith3DCurves",
            "topology.brep_p8",
        ))
    }
    /// Create a wire, consuming its curve on success.
    fn wire(&mut self, brep: BrepId, curve: &mut Option<Curve>) -> Result<EdgeId> {
        let _ = (brep, curve);
        // UNPORTED(topology.brep_p5): SmBrep::CreateWireEdgeFromCurve
        Err(SolidError::unported(
            "SmBrep::CreateWireEdgeFromCurve",
            "topology.brep_p5",
        ))
    }
    /// Split a surface at G1 discontinuities and append its faces.
    fn surface_faces(
        &mut self,
        brep: BrepId,
        surface: Surface,
        domain: Extent2d,
        faces: &mut Vec<FaceId>,
    ) -> Result<()> {
        let _ = (brep, surface, domain, faces);
        // UNPORTED(topology.brep_p8): SmBrep::CreateFacesFromSurface
        Err(SolidError::unported(
            "SmBrep::CreateFacesFromSurface",
            "topology.brep_p8",
        ))
    }
    /// Create a face, retaining the output slot on failure.
    fn surface_face(
        &mut self,
        brep: BrepId,
        surface: Surface,
        domain: Extent2d,
        face: &mut Option<FaceId>,
    ) -> Result<()> {
        let _ = (brep, surface, domain, face);
        // UNPORTED(topology.brep_p8): SmBrep::CreateFaceFromSurface
        Err(SolidError::unported(
            "SmBrep::CreateFaceFromSurface",
            "topology.brep_p8",
        ))
    }
    /// Delete a face without its connected edges/vertices; do region nesting.
    fn delete_primitive_face(&mut self, face: FaceRef) -> Result<()> {
        let _ = face;
        // UNPORTED(topology.brep_p3): SmBrep::DeleteFace
        Err(SolidError::unported(
            "SmBrep::DeleteFace",
            "topology.brep_p3",
        ))
    }
    /// Topology-intersector deletion, retaining the supplied face order.
    fn delete_primitive_faces(&mut self, brep: BrepId, faces: &[FaceId]) -> Result<()> {
        let _ = (brep, faces);
        // UNPORTED(boolean.topology_intersector_p2): SmTopologyIntersector::DeleteFaces
        Err(SolidError::unported(
            "SmTopologyIntersector::DeleteFaces",
            "boolean.topology_intersector_p2",
        ))
    }
    /// Stitch and orient the primitive's resulting B-rep.
    fn stitch_and_orient(&mut self, brep: BrepId) -> Result<()> {
        let _ = brep;
        // UNPORTED(topology.brep_p5): SmBrep::StitchAndOrient
        Err(SolidError::unported(
            "SmBrep::StitchAndOrient",
            "topology.brep_p5",
        ))
    }
    /// Transform the complete B-rep, including pre-existing topology.
    fn transform_primitive(&mut self, brep: BrepId, placement: &Axis2Placement) -> Result<()> {
        let _ = (brep, placement);
        // UNPORTED(topology.brep_p5): SmBrep::Transform
        Err(SolidError::unported(
            "SmBrep::Transform",
            "topology.brep_p5",
        ))
    }
    /// Create composites for all loops of an input face.
    fn face_composites(&mut self, face: FaceRef) -> Result<Vec<CompositeCurve>> {
        let _ = face;
        // UNPORTED(topology.face_p3): SmFace::CreateComposites
        Err(SolidError::unported(
            "SmFace::CreateComposites",
            "topology.face_p3",
        ))
    }
    /// Face's normal bounding box for the rotational axis check.
    fn primitive_face_box(&mut self, face: FaceRef) -> Result<Extent3d> {
        let _ = face;
        // UNPORTED(query.face_p1): SmFace::CalculateBoundingBox
        Err(SolidError::unported(
            "SmFace::CalculateBoundingBox",
            "query.face_p1",
        ))
    }
    /// Run `query::face::p2::curve_on_classify` with caller-resolved object handles;
    /// report whether any interval midpoint is `SM_PC_FACE`.
    fn axis_inside_face(&mut self, face: FaceRef, axis: &Curve) -> Result<bool>;
    /// Linear cap stitch: callback, fast edge compare, no vertex-edge splits,
    /// non-manifold mode, region nesting, default targets and tolerance.
    fn stitch_linear_caps(
        &mut self,
        brep: BrepId,
        sweep: &mut TopologySweep,
        tolerance: f64,
    ) -> Result<()> {
        let _ = (brep, sweep, tolerance);
        // UNPORTED(heal.stitch_p1): SmStitch::DoStitching
        Err(SolidError::unported(
            "SmStitch::DoStitching",
            "heal.stitch_p1",
        ))
    }
}
