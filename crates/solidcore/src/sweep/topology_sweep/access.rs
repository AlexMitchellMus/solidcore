//! Typed access to the topology services called by the sweep algorithm.
//!
//! Upstream: `source/SMLib/src/SmTopologySweep.cpp` (calls to other units).
//!
//! | Rust | Upstream |
//! |------|----------|
//! | [`SweepAccess`] | the B-rep/context/topology services used by `SmTopologySweep` |
//!
//! Deviations: a statically dispatched access interface keeps B-rep and context
//! identities explicit. Defaults below defer only the named unavailable callee;
//! the sweep's remaining control flow is implemented by the caller of each service.

use super::{SweepEntity, TopologySweep};
use crate::curve::curve::Curve;
use crate::error::{Result, SolidError};
use crate::geom::attribute::TagData;
use crate::geom::axis2_placement::Axis2Placement;
use crate::geom::extent1d::Extent1d;
use crate::geom::extent2d::Extent2d;
use crate::geom::extent3d::Extent3d;
use crate::kernel::Kernel;
use crate::math::{Point3, Vec3};
use crate::nurbs::b_spline_surface_n_lib::p2::OrientType;
use crate::surface::surface::Surface;
use crate::topology::{Brep, BrepId, EdgeId, EdgeuseId, FaceId, RegionId, ShellId, VertexId};

/// A region and its owning B-rep.
pub type RegionRef = (BrepId, RegionId);
/// A shell and its owning B-rep.
pub type ShellRef = (BrepId, ShellId);
/// A face and its owning B-rep.
pub type FaceRef = (BrepId, FaceId);
/// An edge and its owning B-rep.
pub type EdgeRef = (BrepId, EdgeId);
/// A vertex and its owning B-rep.
pub type VertexRef = (BrepId, VertexId);
/// An edge use and its owning B-rep.
pub type EdgeuseRef = (BrepId, EdgeuseId);

/// Owned inputs to `SmBrep::MakeFaceWithCurves`.
#[derive(Debug)]
pub struct FaceInput {
    /// Target region.
    pub region: RegionRef,
    /// Edge count of each non-vertex loop.
    pub loop_counts: Vec<usize>,
    /// Trimming curves in loop order.
    pub curves: Vec<Curve>,
    /// Optional UV trimming curves.
    pub uv_curves: Option<Vec<Curve>>,
    /// Orientations in loop order.
    pub orientations: Vec<OrientType>,
    /// Locations of vertex loops.
    pub loop_points: Vec<Point3>,
    /// Surface being bounded.
    pub surface: Surface,
    /// Domain of that surface.
    pub domain: Extent2d,
    /// Face/surface orientation.
    pub orientation: OrientType,
}

/// Outputs of `MakeFaceWithCurves`, retaining each input curve's new owner.
#[derive(Debug)]
pub struct FaceBuild {
    /// Created face.
    pub face: FaceRef,
    /// Owner of each input curve, in input order; null curves have no owner.
    pub curve_edges: Vec<Option<EdgeRef>>,
}

/// Identifies a caller-owned `SmMerge` port instance for its scoped notifications.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SweepMergeId(pub u32);

/// Hidden-curve-removal scene. The access implementation retains the HCR session
/// until `delete_brep(brep)`, so silhouettes are queried after copying the scene.
#[derive(Debug)]
pub struct SweepScene {
    /// HCR scene B-rep, owned by this scoped session.
    pub brep: BrepId,
}

macro_rules! deferred {
    ($(#[$doc:meta])* fn $name:ident($($arg:ident: $ty:ty),* $(,)?) -> $out:ty; $(,)? $unit:literal, $symbol:literal $(,)?) => {
        $(#[$doc])*
        fn $name(&mut self, $($arg: $ty),*) -> Result<$out> {
            $(let _ = $arg;)*
            Err(SolidError::Unported { upstream: $symbol, unit: $unit })
        }
    };
}

/// B-rep resolution and external services for a topology sweep.
///
/// Implementations must resolve contexts and attribute keys using the caller's
/// object identities. Default dependency methods are explicit `Unported` calls;
/// they do not claim or reproduce code from open dependency PRs.
pub trait SweepAccess {
    /// Resolve a B-rep without changing its topology.
    fn brep(&self, id: BrepId) -> Result<&Brep>;
    /// Resolve a mutable B-rep.
    fn brep_mut(&mut self, id: BrepId) -> Result<&mut Brep>;
    /// Context of a B-rep, including its tolerances and caches.
    fn kernel(&self, id: BrepId) -> Result<&Kernel>;
    /// Mutable context; B-reps sharing a context must resolve to the same kernel.
    fn kernel_mut(&mut self, id: BrepId) -> Result<&mut Kernel>;
    /// Allocate an empty B-rep in `context_of`'s context.
    fn new_brep(&mut self, context_of: BrepId) -> Result<BrepId>;
    /// Release live temporary topology and geometry, including on error.
    /// A successful `union` already consumes its source, so it is never released twice.
    fn delete_brep(&mut self, id: BrepId);
    /// Find a tag by resolved object identity.
    fn tag(&self, entity: SweepEntity) -> Result<Option<TagData>>;
    /// Attach a tag using `SmAObject::AddAttribute` semantics.
    fn add_tag(&mut self, entity: SweepEntity, tag: TagData) -> Result<()>;
    /// Find the tag on an edge's original curve (null curve is an error).
    fn curve_tag(&self, edge: EdgeRef) -> Result<Option<TagData>>;
    /// Read/write the upstream editing flag. Return its previous value.
    fn set_editing(&mut self, brep: BrepId, editing: bool) -> Result<bool>;

    /// Release a scoped merge object and its notification registrations.
    fn delete_merge(&mut self, merge: SweepMergeId);

    // UNPORTED(boolean.merge_p2): SmMerge::SmMerge
    deferred!(
        /// Construct the merge notifier with both operands and caller context.
        fn new_merge(
            context_of: BrepId,
            destination: BrepId,
            source: BrepId,
            tolerance: f64,
            angle: f64,
        ) -> SweepMergeId;,
        "boolean.merge_p2",
        "SmMerge::SmMerge"
    );

    /// Keep a merge notifier alive for the complete upstream scope, releasing it
    /// after success or an early error.
    fn with_merge<T>(
        &mut self,
        context_of: BrepId,
        destination: BrepId,
        source: BrepId,
        tolerance: f64,
        angle: f64,
        body: impl FnOnce(&mut Self, SweepMergeId) -> Result<T>,
    ) -> Result<T>
    where
        Self: Sized,
    {
        let merge = self.new_merge(context_of, destination, source, tolerance, angle)?;
        let result = body(self, merge);
        self.delete_merge(merge);
        result
    }

    // UNPORTED(topology.brep_p1): SmBrep::GetInfiniteRegion
    deferred!(
        /// `SmBrep::GetInfiniteRegion`.
        fn infinite_region(brep: BrepId) -> RegionId;,
        "topology.brep_p1",
        "SmBrep::GetInfiniteRegion"
    );
    /// `SmFace::GetEdges`, using the existing faithful co-port.
    fn face_edges(&mut self, face: FaceRef) -> Result<Vec<EdgeId>> {
        crate::topology::topology_traverser::face_edges(self.brep(face.0)?, face.1)
    }
    /// `SmFace::GetVertices`, using the existing faithful co-port.
    fn face_vertices(&mut self, face: FaceRef) -> Result<Vec<VertexId>> {
        crate::query::face::p2::face_vertices(self.brep(face.0)?, face.1, false, None)
    }
    /// Acquire the caller context's unlocked mark and retain its lock until release.
    fn lock_mark(&mut self, brep: BrepId, mask: i32) -> Result<crate::geom::context::MarkType>;
    /// Release that mark lock (also on an early return).
    fn unlock_mark(&mut self, brep: BrepId, mark: crate::geom::context::MarkType) -> Result<()>;
    /// Query the persistent topology mark record using the owning context.
    fn is_marked(&self, entity: SweepEntity, mark: crate::geom::context::MarkType) -> Result<bool>;
    /// Set the persistent topology mark record using the owning context.
    fn mark(&mut self, entity: SweepEntity, mark: crate::geom::context::MarkType) -> Result<()>;
    /// `SmTopologyTraverser::CollectFaces` with the current lock and mark records.
    fn collect_faces(
        &mut self,
        face: FaceRef,
        mark: crate::geom::context::MarkType,
        faces: &mut Vec<FaceId>,
    ) -> Result<()>;

    // UNPORTED(topology.brep_p9): SmBrep::MakeShellVertex
    deferred!(
        /// `SmBrep::MakeShellVertex`.
        fn make_shell_vertex(region: RegionRef, point: Point3) -> (ShellRef, VertexRef);,
        "topology.brep_p9",
        "SmBrep::MakeShellVertex"
    );
    // UNPORTED(topology.brep_p3): SmBrep::MakeWireEdgeVertex
    deferred!(
        /// `SmBrep::MakeWireEdgeVertex`.
        fn make_wire_edge_vertex(
            region: RegionRef,
            start: VertexRef,
            curve: Curve,
            interval: Extent1d,
            orientation: OrientType,
            end: Point3,
        ) -> (EdgeRef, VertexRef);,
        "topology.brep_p3",
        "SmBrep::MakeWireEdgeVertex"
    );
    // UNPORTED(topology.brep_p3): SmBrep::MakeWireEdge
    deferred!(
        /// `SmBrep::MakeWireEdge`.
        fn make_wire_edge(
            region: RegionRef,
            start: VertexRef,
            end: VertexRef,
            curve: Curve,
            interval: Extent1d,
            orientation: OrientType,
        ) -> EdgeRef;,
        "topology.brep_p3",
        "SmBrep::MakeWireEdge"
    );
    // UNPORTED(topology.brep_p7): SmBrep::MakeFaceWithCurves
    deferred!(
        /// `SmBrep::MakeFaceWithCurves` and the resulting curve owners.
        fn make_face(input: FaceInput) -> FaceBuild;,
        "topology.brep_p7",
        "SmBrep::MakeFaceWithCurves"
    );
    // UNPORTED(topology.face_p3): SmFace::SetTolerance
    deferred!(
        /// `SmFace::SetTolerance` with upstream's default options.
        fn set_face_tolerance(face: FaceRef, tolerance: f64) -> ();,
        "topology.face_p3",
        "SmFace::SetTolerance"
    );
    // UNPORTED(topology.brep_p9): SmBrep::MakeVertexLoop
    deferred!(
        /// `SmBrep::MakeVertexLoop` (the new loop is not used by the sweep).
        fn make_vertex_loop(face: FaceRef, point: Point3) -> VertexRef;,
        "topology.brep_p9",
        "SmBrep::MakeVertexLoop"
    );
    // UNPORTED(topology.brep_p5): SmBrep::MakeEdgeInFace
    deferred!(
        /// `SmBrep::MakeEdgeInFace` (new loop/face outputs are unused).
        fn make_edge_in_face(
            face: FaceRef,
            start: VertexRef,
            end: VertexRef,
            curve: Curve,
            uv_curve: Option<Curve>,
            interval: Extent1d,
            orientation: OrientType,
            tolerance: f64,
        ) -> EdgeRef;,
        "topology.brep_p5",
        "SmBrep::MakeEdgeInFace"
    );
    // UNPORTED(topology.brep_p1): SmBrep::GetFaces
    deferred!(
        /// `SmBrep::GetFaces`.
        fn faces(brep: BrepId) -> Vec<FaceId>;,
        "topology.brep_p1",
        "SmBrep::GetFaces"
    );
    // UNPORTED(topology.shell): SmRegion::GetShells
    deferred!(
        /// `SmRegion::GetShells`.
        fn region_shells(region: RegionRef) -> Vec<ShellId>;,
        "topology.shell",
        "SmRegion::GetShells"
    );
    // UNPORTED(topology.shell): SmShell::GetFaceuses
    deferred!(
        /// Faces of the shell's face uses, retaining duplicates and order.
        fn shell_faceuse_faces(shell: ShellRef) -> Vec<FaceId>;,
        "topology.shell",
        "SmShell::GetFaceuses"
    );
    // UNPORTED(topology.brep_p5): SmBrep::RemoveFaces
    deferred!(
        /// `SmBrep::RemoveFaces`.
        fn remove_faces(brep: BrepId, faces: &[FaceId]) -> ();,
        "topology.brep_p5",
        "SmBrep::RemoveFaces"
    );
    // UNPORTED(blend.hcr): SmHCR::ComputeGlobalVisibility
    deferred!(
        /// HCR constructor, silhouette mode, visibility and scene ownership.
        fn hcr(
            brep: BrepId,
            inverse: Axis2Placement,
            tolerance: f64,
            angle: f64,
            seams: bool,
            smooth: bool,
            silhouettes: bool,
        ) -> SweepScene;,
        "blend.hcr",
        "SmHCR::ComputeGlobalVisibility"
    );
    // UNPORTED(blend.hcr): SmHCR::GetSilhouetteEdges
    deferred!(
        /// Query the retained HCR session at the original call site.
        fn silhouette_edges(scene: BrepId) -> Vec<EdgeId>;,
        "blend.hcr",
        "SmHCR::GetSilhouetteEdges"
    );
    // UNPORTED(topology.brep_p1): SmBrep::GetWireEdges
    deferred!(
        /// `SmBrep::GetWireEdges`.
        fn wire_edges(brep: BrepId) -> Vec<EdgeId>;,
        "topology.brep_p1",
        "SmBrep::GetWireEdges"
    );
    // UNPORTED(topology.brep_p2): SmBrep::GetShellVertices
    deferred!(
        /// `SmBrep::GetShellVertices`.
        fn shell_vertices(brep: BrepId) -> Vec<VertexId>;,
        "topology.brep_p2",
        "SmBrep::GetShellVertices"
    );
    // UNPORTED(topology.brep_p6): SmBrep::CopyEdges
    deferred!(
        /// `SmBrep::CopyEdges` with its output list supplied.
        fn copy_edges(source: BrepId, edges: &[EdgeId], destination: BrepId) -> Vec<EdgeId>;,
        "topology.brep_p6",
        "SmBrep::CopyEdges"
    );
    // UNPORTED(topology.brep_p3): SmBrep::CopyVertices
    deferred!(
        /// `SmBrep::CopyVertices` with its output list supplied.
        fn copy_vertices(
            source: BrepId,
            vertices: &[VertexId],
            destination: BrepId,
        ) -> Vec<VertexId>;,
        "topology.brep_p3",
        "SmBrep::CopyVertices"
    );
    // UNPORTED(topology.brep_p10): SmBrep::RemoveAllTopology
    deferred!(
        /// `SmBrep::RemoveAllTopology`.
        fn remove_all_topology(brep: BrepId) -> ();,
        "topology.brep_p10",
        "SmBrep::RemoveAllTopology"
    );
    // UNPORTED(topology.brep_p7): SmBrep::CopyFaces
    deferred!(
        /// `SmBrep::CopyFaces` with merge notification; returns scene edge mates.
        fn copy_faces(
            source: BrepId,
            faces: &[FaceId],
            destination: BrepId,
            use_merge: bool,
            merge: SweepMergeId,
            edge_mates: &mut BTreeEdgeMap,
        ) -> ();,
        "topology.brep_p7",
        "SmBrep::CopyFaces"
    );
    // UNPORTED(topology.brep_p5): SmBrep::Transform
    deferred!(
        /// `SmBrep::Transform`.
        fn transform(brep: BrepId, frame: Axis2Placement) -> ();,
        "topology.brep_p5",
        "SmBrep::Transform"
    );
    // UNPORTED(query.brep): SmBrep::CalculateBoundingBox
    deferred!(
        /// `SmBrep::CalculateBoundingBox`.
        fn bounding_box(brep: BrepId, output: &mut Extent3d) -> ();,
        "query.brep",
        "SmBrep::CalculateBoundingBox"
    );
    // UNPORTED(topology.edgeuse_p3): SmEdgeuse::EvaluateBinormal
    deferred!(
        /// `SmEdgeuse::EvaluateBinormal`: position and binormal.
        fn binormal(
            edgeuse: EdgeuseRef,
            parameter: f64,
            unitize: bool,
            point: &mut Point3,
            binormal: &mut Vec3,
        ) -> ();,
        "topology.edgeuse_p3",
        "SmEdgeuse::EvaluateBinormal"
    );
    // UNPORTED(topology.face_p3): SmFace::HasLaminaEdge
    deferred!(
        /// `SmFace::HasLaminaEdge`.
        fn has_lamina_edge(face: FaceRef) -> bool;,
        "topology.face_p3",
        "SmFace::HasLaminaEdge"
    );
    // UNPORTED(topology.brep_p3): SmBrep::MergeBrep
    deferred!(
        /// `SmBrep::MergeBrep`.
        fn merge_brep(destination: BrepId, source: BrepId) -> ();,
        "topology.brep_p3",
        "SmBrep::MergeBrep"
    );
    // UNPORTED(boolean.merge_p3): SmMerge::PiecewiseMerge
    deferred!(
        /// `SmMerge::PiecewiseMerge(NULL, TRUE, TRUE, ...)`.
        fn piecewise_merge(
            merge: SweepMergeId,
            use_existing_topology: bool,
            create_new_breps_for_faces: bool,
        ) -> BrepId;,
        "boolean.merge_p3",
        "SmMerge::PiecewiseMerge"
    );
    // UNPORTED(heal.stitch_p1): SmStitch::DoStitching
    deferred!(
        /// Sweep callback, default targets, non-manifold mode.
        fn stitch(
            brep: BrepId,
            sweep: &mut TopologySweep,
            tolerance: f64,
            making_manifold_solid: bool,
        ) -> ();,
        "heal.stitch_p1",
        "SmStitch::DoStitching"
    );
    // UNPORTED(topology.brep_p3): SmBrep::OrientTrimmedSurfaces
    deferred!(
        /// `SmBrep::OrientTrimmedSurfaces`.
        fn orient(brep: BrepId, outward: bool, maybe_open: bool, unify_only: bool) -> ();,
        "topology.brep_p3",
        "SmBrep::OrientTrimmedSurfaces"
    );
    // UNPORTED(boolean.merge_p4): SmMerge::ManifoldBoolean
    deferred!(
        /// `SmMerge::ManifoldBoolean(SM_BO_UNION, ...)`. Success consumes and
        /// unregisters the source B-rep; an error leaves it owned by the caller.
        fn union(merge: SweepMergeId) -> BrepId;,
        "boolean.merge_p4",
        "SmMerge::ManifoldBoolean"
    );
}

/// Correspondence from source to copied edges, filled by the merge notifier.
pub type BTreeEdgeMap = std::collections::BTreeMap<EdgeId, EdgeId>;
