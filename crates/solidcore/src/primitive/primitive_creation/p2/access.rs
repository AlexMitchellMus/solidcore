//! External B-rep services of the p2 primitive entry points.
//!
//! Upstream: `source/SMLib/src/SmPrimitiveCreation.cpp` (callee boundaries),
//! `source/SMLib/inc/SmTypes.h` (`SmTrimType`).
//!
//! | Rust | Upstream |
//! |------|----------|
//! | [`PrimitiveCreationAccess`] | B-rep, face, merge and stitching services |
//! | [`TrimType`] | `SmTrimType` (co-ported for `SmBrep::ProjectAndTrim`) |
//!
//! Deviations: extends p3's [`PrimitiveAccess`] as that trait extends
//! `SweepAccess`. Defaults defer only callees whose unit is not ported yet;
//! outputs no p2 caller reads (stitch counts and gaps) are not returned.

use crate::curve::composite_curve::OffsetMapAttribute;
use crate::error::{Result, SolidError};
use crate::math::{Point3, Vec3};
use crate::nurbs::b_spline_curve::p1::BSplineCurve;
use crate::primitive::primitive_creation::p3::PrimitiveAccess;
use crate::sweep::topology_sweep::access::FaceRef;
use crate::topology::{BrepId, FaceId};

/// Which side of a projected trimming curve `SmBrep::ProjectAndTrim` keeps
/// (`SmTrimType`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TrimType {
    /// `SM_TT_KEEP_POINT`: keep the side containing the reference point.
    KeepPoint,
    /// `SM_TT_DELETE_POINT`: remove the side containing the reference point.
    DeletePoint,
    /// `SM_TT_SPLIT`: keep both sides.
    Split,
}

/// Topology adapters used by the p2 primitives, on top of [`PrimitiveAccess`].
/// Construction tails never replace upstream's topological edits with arena-only
/// approximations.
pub trait PrimitiveCreationAccess: PrimitiveAccess {
    /// Copy `faces` of `source` into `destination` face by face, without the merge
    /// (`SmBrep::CopyFaces(faces, destination, NULL, FALSE)`; its result is unused
    /// here). Implementations resolve both B-reps for `topology::brep::p7::copy_faces`.
    fn copy_faces_into(
        &mut self,
        source: BrepId,
        faces: &[FaceId],
        destination: BrepId,
    ) -> Result<()>;
    /// Attach an offset map to the surface of `face` (`SmAObject::AddAttribute` of an
    /// `SM_AI_OFFSETMAP` attribute), in the caller's attribute store.
    fn attach_offset_map(&mut self, face: FaceRef, map: OffsetMapAttribute) -> Result<()>;
    /// The offset map of the surface of `face`, if any
    /// (`FindAttribute(SM_AI_OFFSETMAP)` on its surface).
    fn offset_map(&self, face: FaceRef) -> Result<Option<OffsetMapAttribute>>;

    /// Set the B-rep's zone tolerance (`SmBrep::SetTolerance(tol)`).
    fn set_brep_tolerance(&mut self, brep: BrepId, tolerance: f64) -> Result<()> {
        let _ = (brep, tolerance);
        // UNPORTED(topology.brep_p2): SmBrep::SetTolerance
        Err(SolidError::unported(
            "SmBrep::SetTolerance",
            "topology.brep_p2",
        ))
    }
    /// Shrink every face's surface to its UV trim curves and every edge's curve to
    /// the edge (`SmBrep::ShrinkGeometry`).
    fn shrink_geometry(&mut self, brep: BrepId) -> Result<()> {
        let _ = brep;
        // UNPORTED(topology.brep_p4): SmBrep::ShrinkGeometry
        Err(SolidError::unported(
            "SmBrep::ShrinkGeometry",
            "topology.brep_p4",
        ))
    }
    /// Whether two consecutive edges of the face's outer loop are lamina
    /// (`SmFace::HasConsecutiveLaminaEdges`).
    fn has_consecutive_lamina_edges(&mut self, face: FaceRef) -> Result<bool> {
        let _ = face;
        // UNPORTED(topology.face_p3): SmFace::HasConsecutiveLaminaEdges
        Err(SolidError::unported(
            "SmFace::HasConsecutiveLaminaEdges",
            "topology.face_p3",
        ))
    }
    /// Delete a face with the edges and vertices left unconnected, nesting regions
    /// (`SmBrep::DeleteFace(face, TRUE)`).
    fn delete_face_and_free_topology(&mut self, face: FaceRef) -> Result<()> {
        let _ = face;
        // UNPORTED(topology.brep_p3): SmBrep::DeleteFace
        Err(SolidError::unported(
            "SmBrep::DeleteFace",
            "topology.brep_p3",
        ))
    }
    /// Glue coincident boundary topology of the B-rep's faces within `tolerance`
    /// (`SmBrep::StitchFaces`; its counts and gaps are unused here).
    fn stitch_faces(&mut self, brep: BrepId, tolerance: f64) -> Result<()> {
        let _ = (brep, tolerance);
        // UNPORTED(topology.brep_p4): SmBrep::StitchFaces
        Err(SolidError::unported(
            "SmBrep::StitchFaces",
            "topology.brep_p4",
        ))
    }
    /// Project `curve` along `direction` onto the B-rep's faces and trim them
    /// (`SmBrep::ProjectAndTrim(curve, direction, reference, type, &tol, NULL, NULL)`).
    fn project_and_trim(
        &mut self,
        brep: BrepId,
        curve: &BSplineCurve,
        direction: Vec3,
        reference: Point3,
        trim_type: TrimType,
        approx_tol: f64,
    ) -> Result<()> {
        let _ = (brep, curve, direction, reference, trim_type, approx_tol);
        // UNPORTED(topology.brep_p10): SmBrep::ProjectAndTrim
        Err(SolidError::unported(
            "SmBrep::ProjectAndTrim",
            "topology.brep_p10",
        ))
    }
    /// Remove edges and vertices that only split faces or edges
    /// (`SmBrep::RemoveTopologicalEdgesAndVertices`).
    fn remove_topological_edges_and_vertices(&mut self, brep: BrepId) -> Result<()> {
        let _ = brep;
        // UNPORTED(topology.brep_p5): SmBrep::RemoveTopologicalEdgesAndVertices
        Err(SolidError::unported(
            "SmBrep::RemoveTopologicalEdgesAndVertices",
            "topology.brep_p5",
        ))
    }
    /// Whether the B-rep is a closed manifold solid (`SmBrep::IsManifoldSolid`).
    fn is_manifold_solid(&mut self, brep: BrepId) -> Result<bool> {
        let _ = brep;
        // UNPORTED(topology.brep_p2): SmBrep::IsManifoldSolid
        Err(SolidError::unported(
            "SmBrep::IsManifoldSolid",
            "topology.brep_p2",
        ))
    }
    /// Combine `breps` by the postfix tree `tree`, leaving the result first
    /// (`SmMerge::BooleanTrees`).
    fn boolean_trees(&mut self, breps: &mut Vec<BrepId>, tree: &[i64]) -> Result<()> {
        let _ = (breps, tree);
        // UNPORTED(boolean.merge_p3): SmMerge::BooleanTrees
        Err(SolidError::unported(
            "SmMerge::BooleanTrees",
            "boolean.merge_p3",
        ))
    }
}
