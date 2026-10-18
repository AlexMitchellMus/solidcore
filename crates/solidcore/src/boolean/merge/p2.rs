//! The merge constructor, copying an other-operand edge or face use into the primary
//! B-rep, and merging freshly added topology into the rest of the B-rep.
//!
//! Upstream: `source/SMLib/src/SmMerge.cpp`, `source/SMLib/inc/SmMerge.h`.
//!
//! | Rust                               | Upstream                                       |
//! |------------------------------------|------------------------------------------------|
//! | [`MergeState::with_breps`]         | `SmMerge::SmMerge(crContext, pBrep, pOtherBrep, ...)` |
//! | [`merge_edge`]                     | `SmMerge::MergeEdge`                           |
//! | [`make_faceuse`], [`MadeFaceuse`]  | `SmMerge::MakeFaceuse`                         |
//! | [`TopologyCopy`]                   | the `Notify(SM_NO_COPY, ...)` identities       |
//! | `reverse_orientation` (private)    | `SM_REVERSE_ORIENTATION`                       |
//!
//! Deviations:
//! - Operands, mates and notification storage are borrowed explicitly, as in
//!   [`p1`](super::p1) and [`p5`](super::p5); arena ids replace pointers.
//!   `MergeEdge` and `MakeFaceuse` work on concrete [`Brep`] operands and the
//!   [`Relations`] mates, which gain face mates and edge/face `Relate` here.
//! - **Curve and surface mates.** An arena edge owns its curve and a face its
//!   surface, so neither can be shared by two edges or faces and neither has a mate
//!   of its own. `GetBrepMate(pCurve)` (`MergeEdge`) and `GetBrepMate(pOldSurface)`
//!   (`MakeFaceuse`) are reached only for an edge or face without a mate, so they
//!   are always null: the curve is always copied (`SmCurve::Copy`) and the surface
//!   always copied with `CopyAndAddAnalytics`. Their `Relate` calls record nothing
//!   observable and are not kept.
//! - Out-params become return values; on an error the outputs are not returned
//!   (upstream leaves `pOldRegion` null and `pNewRegion` as `MakeFace` left it).
//!   `MakeFaceuse`'s optional new-face list is cleared on entry and receives the new
//!   face on success, as upstream.
//! - Undefined behaviour upstream (a null mate faceuse, shell, edge, curve, vertex
//!   use or vertex dereferenced) is [`SolidError::InvalidInput`]; upstream's `NER`
//!   checks are `SM_ERR_NULL_POINTER`.
//! - `MakeFaceuse` treats a failed `FindRadialSector` as `SM_ERR`, as upstream,
//!   except [`SolidError::Unported`] and [`SolidError::InvalidInput`], which
//!   propagate.
//! - The infinite region `MakeFaceuse` clears around `MakeFace` is restored before
//!   `MakeFace`'s error propagates, as upstream.
//! - Callees not ported yet return [`SolidError::Unported`] at the call site,
//!   tagged `UNPORTED(<unit>)`: `SmBrep::GetMinEdgeBBoxSize` (`topology.brep_p2`),
//!   `SmBrep::GetTolerance` (`topology.brep_p1`), `SmBrep::MakeWireEdge`
//!   (`topology.brep_p3`) and `SmBrep::MakeFace` (`topology.brep_p7`). The code
//!   after each stub is ported in full.
//! - The `SM_DEBUG_CODE` / `SM_VALIDATE_TOPOLOGY` blocks are compiled out upstream
//!   and not ported. `SM_USE_OLDTOL` is the reference build's setting.
//! - **Recursion.** None.

use crate::boolean::merge::p1::{self, MergeOptions, MergeState, VertexCopy};
use crate::boolean::topology_intersector::p1::{Relations, TopologyIntersector};
use crate::curve::curve::p2::Curve;
use crate::error::{Result, SolidError};
use crate::geom::context::Context;
use crate::geom::extent1d::Extent1d;
use crate::geom::extent2d::Extent2d;
use crate::geom::object::{NotifyOperation, ObjectRef};
use crate::geom::solution_array::SolutionObject;
use crate::intersect::point_class::ClassObject;
use crate::kernel::Kernel;
use crate::math::{EFF_ZERO, PI};
use crate::nurbs::b_spline_surface_n_lib::p2::OrientType;
use crate::surface::surface::p2::Surface;
use crate::topology::edgeuse::{p1 as eu, p3 as eu3};
use crate::topology::topology as topo;
use crate::topology::{
    edge, loopuse, vertexuse, Brep, EdgeId, EdgeuseId, FaceId, FaceuseId, RegionId, ShellId,
    VertexId,
};

const OP: &str = "boolean::merge::p2";
const OP_EDGE: &str = "boolean::merge::p2::merge_edge";
const OP_FACEUSE: &str = "boolean::merge::p2::make_faceuse";

/// A dereference of a null pointer upstream (undefined behaviour).
const fn dangling(operation: &'static str, reason: &'static str) -> SolidError {
    SolidError::InvalidInput { operation, reason }
}

fn unavailable<T>(symbol: &'static str, unit: &'static str) -> Result<T> {
    Err(SolidError::unported(symbol, unit))
}

/// `SM_REVERSE_ORIENTATION(a)`: `SM_OT_SAME` becomes `SM_OT_OPPOSITE`, anything
/// else (including `SM_OT_UNKNOWN`) becomes `SM_OT_SAME`.
const fn reverse_orientation(o: OrientType) -> OrientType {
    match o {
        OrientType::Same => OrientType::Opposite,
        OrientType::Opposite | OrientType::Unknown => OrientType::Same,
    }
}

// ---- the constructor ---------------------------------------------------------------

fn min_edge_bbox_size(_brep: &Brep) -> Result<f64> {
    // UNPORTED(topology.brep_p2): SmBrep::GetMinEdgeBBoxSize
    unavailable("SmBrep::GetMinEdgeBBoxSize", "topology.brep_p2")
}

fn brep_tolerance(_brep: &Brep) -> Result<f64> {
    // UNPORTED(topology.brep_p1): SmBrep::GetTolerance
    unavailable("SmBrep::GetTolerance", "topology.brep_p1")
}

impl<O: ClassObject + SolutionObject + Ord> MergeState<TopologyIntersector<O>> {
    /// `SmMerge(crContext, pBrep, pOtherBrep, dApproxTol3d, dAngleTolRadians,
    /// pOptMergeOptions)`: the intersector over both B-reps with the given
    /// tolerances, the default controls, then the options (a null pointer keeps the
    /// defaults). A tolerance below `SM_EFF_ZERO` is replaced by
    /// `min(0.5 * min(brep tolerances), 0.33 * min edge bounding-box size)`; an angle
    /// tolerance of exactly `0.0` by 20 degrees. `brep` and `other` give each B-rep's
    /// arena and the handle the intersector refers to it by; `other` may be `None`
    /// (the stitch entry points) unless the tolerance must be computed.
    ///
    /// # Errors
    ///
    /// [`SolidError::InvalidInput`] for a computed tolerance without other B-rep
    /// (upstream dereferences null); [`SolidError::Unported`] from the B-rep queries.
    pub fn with_breps(
        brep: (&Brep, O),
        other: Option<(&Brep, O)>,
        approx_tol_3d: f64,
        ang_tol_rad: f64,
        options: Option<&MergeOptions>,
    ) -> Result<Self> {
        let intersector =
            TopologyIntersector::new(Some(brep.1), other.map(|o| o.1), approx_tol_3d, ang_tol_rad);
        // m_bManifoldBoolean(TRUE), m_bCookieCutter(FALSE), m_bImprinting(FALSE),
        // m_bBooleanPostProcess(TRUE), m_bIntersectionLoopClosed(TRUE),
        // m_bImprintAndClassifyFaces(FALSE), m_bKeepOtherBrep(FALSE), m_lTagID(0)
        let mut state = Self::new(intersector);
        state.set_merge_options(options);

        // The user can ask us to calculate a tolerance based on the geometry by
        // specifying a zero tolerance.
        if approx_tol_3d < EFF_ZERO {
            // Get the minimum Edge length.
            let edge1_size = min_edge_bbox_size(brep.0)?;
            let (other_brep, _) =
                other.ok_or(dangling(OP, "computed tolerance without other B-rep"))?;
            let edge2_size = min_edge_bbox_size(other_brep)?;
            let edge_size = edge1_size.min(edge2_size);

            // siWhichTol == 3: 1/2 of Min of two Brep's tols.
            let brep_tol = 0.5 * brep_tolerance(brep.0)?.min(brep_tolerance(other_brep)?);

            state
                .intersector
                .set_this_approx_tol_3d(brep_tol.min(edge_size * 0.33));
        }

        // Use 20 degrees as a default angle tolerance
        #[allow(clippy::float_cmp)] // upstream's exact `== 0.0`
        if ang_tol_rad == 0.0 {
            state.intersector.set_this_ang_tol_rad(20.0 * (PI / 180.0));
        }
        Ok(state)
    }
}

// ---- copying topology from the other operand -----------------------------------------

/// The notification storage and identities of the topology copied from the other
/// operand (`Notify(SM_NO_COPY, pCopy, ...)` of vertices, edges and faces). Keys must
/// be distinct across operands even when their arena ids coincide.
pub struct TopologyCopy<'a> {
    /// Records, attributes and vertex identities ([`p1::merge_vertex`]).
    pub vertices: VertexCopy<'a>,
    /// Notification identity of each other-operand edge.
    pub edge_source_key: &'a dyn Fn(EdgeId) -> ObjectRef,
    /// Notification identity of each primary edge, including freshly created ids.
    pub edge_destination_key: &'a dyn Fn(EdgeId) -> ObjectRef,
    /// Notification identity of each other-operand face.
    pub face_source_key: &'a dyn Fn(FaceId) -> ObjectRef,
    /// Notification identity of each primary face, including freshly created ids.
    pub face_destination_key: &'a dyn Fn(FaceId) -> ObjectRef,
}

impl TopologyCopy<'_> {
    /// `pSource->Notify(SM_NO_COPY, pCopy, SM_NO_GET_BREP(pCopy),
    /// SM_NO_GET_BREP(pSource))`: `SM_NO_GET_BREP` is null in the pinned source.
    fn notify(
        &mut self,
        context: Option<&mut Context>,
        source: ObjectRef,
        copy: ObjectRef,
    ) -> Result<()> {
        topo::notify(
            self.vertices.records,
            self.vertices.attributes,
            context,
            NotifyOperation::Copy,
            source,
            Some(copy),
            None,
            None,
        )
    }
}

#[allow(clippy::too_many_arguments)] // upstream's MakeWireEdge
fn make_wire_edge(
    brep: &mut Brep,
    region: Option<RegionId>,
    start: VertexId,
    end: VertexId,
    curve: Curve,
    interval: Extent1d,
    orientation: OrientType,
) -> Result<EdgeId> {
    let _ = (brep, region, start, end, curve, interval, orientation);
    // UNPORTED(topology.brep_p3): SmBrep::MakeWireEdge
    unavailable("SmBrep::MakeWireEdge", "topology.brep_p3")
}

/// Merge an edge of the other operand into the primary B-rep (`MergeEdge`). An edge
/// without a primary mate becomes a new primary wire edge in `region` (`None`: the
/// infinite region) over a copy of its curve, between the primary mates of its
/// vertices (merged in as shell vertices when they have none), and is related to it.
/// The primary edge then gets the larger tolerance and the source edge is notified
/// of the copy. A null edge use succeeds with no edge, as upstream.
///
/// # Errors
///
/// `SM_ERR_NULL_POINTER` when a vertex has no primary mate after merging;
/// [`SolidError::InvalidInput`] for a null edge, curve, vertex use or vertex of the
/// edge use (upstream dereferences them); those of the curve copy, `Relate` and
/// [`p1::merge_vertex`]; [`SolidError::Unported`] (module docs).
#[allow(clippy::too_many_arguments)]
pub fn merge_edge(
    k: &Kernel,
    brep: &mut Brep,
    other: &mut Brep,
    edgeuse: Option<EdgeuseId>,
    region: Option<RegionId>,
    relations: &mut Relations,
    copy: &mut TopologyCopy<'_>,
    mut context: Option<&mut Context>,
) -> Result<Option<EdgeId>> {
    // The caller knows better whether this is a problem.
    let Some(edgeuse) = edgeuse else {
        return Ok(None);
    };

    // get target edge
    let other_edge = eu::edge(other, edgeuse)?.ok_or(dangling(OP_EDGE, "edge use without edge"))?;

    // get Brep matching edge (probably NULL)
    let brep_edge = if let Some(brep_edge) = relations.edges.primary(other_edge) {
        brep_edge
    } else {
        // copy the curve (its Brep mate is always NULL; module docs)
        let source = other.edge(other_edge)?;
        let interval = source.interval;
        let new_curve = source
            .curve
            .as_ref()
            .ok_or(dangling(OP_EDGE, "edge without curve"))?
            .copy(k)?;

        // get OtherBrep edge->vertices
        let vu = eu::vertexuse(other, edgeuse)?
            .ok_or(dangling(OP_EDGE, "edge use without vertex use"))?;
        let mut v = vertexuse::vertex(other, vu)?;
        let mut v2 = edge::other_vertex(other, other_edge, v)?;

        // Swap vertices if the edgeuse is oriented opposite to the edge.
        if eu::orientation(other, edgeuse)? == OrientType::Opposite {
            core::mem::swap(&mut v, &mut v2);
        }

        // Get or make matching Brep vertices
        let mut mates = [None, None];
        for (slot, v) in mates.iter_mut().zip([v, v2]) {
            let mut mate = v.and_then(|v| relations.vertices.primary(v));
            if mate.is_none() {
                mate = p1::merge_vertex(
                    brep,
                    other,
                    v,
                    region,
                    relations,
                    &mut copy.vertices,
                    context.as_deref_mut(),
                )?;
            }
            *slot = Some(mate.ok_or(SolidError::null_pointer(OP_EDGE))?);
        }
        let [Some(bv), Some(bv2)] = mates else {
            return Err(SolidError::null_pointer(OP_EDGE));
        };

        // make the new Brep wire edge
        let brep_edge =
            make_wire_edge(brep, region, bv, bv2, new_curve, interval, OrientType::Same)?;

        // mate new edge with input otherBrep edge
        relations.relate_edges(brep, brep_edge, other, other_edge)?;
        brep_edge
    };

    // match tolerances for the mated edges
    let other_tolerance = other.edge(other_edge)?.tolerance;
    if other_tolerance > brep.edge(brep_edge)?.tolerance {
        edge::set_tolerance(brep, brep_edge, other_tolerance, true, true)?;
    }

    // notify OtherBrepEdge of the copy (so it can propagate attributes as needed)
    let (source, dest) = (
        (copy.edge_source_key)(other_edge),
        (copy.edge_destination_key)(brep_edge),
    );
    copy.notify(context, source, dest)?;

    Ok(Some(brep_edge))
}

/// What [`make_faceuse`] returns.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct MadeFaceuse {
    /// `pOldRegion`: the primary region receiving the new face (`None` when the face
    /// already had a mate).
    pub old_region: Option<RegionId>,
    /// `pNewRegion`: the region made when the new face splits `old_region`.
    pub new_region: Option<RegionId>,
}

fn shell_region_of(b: &Brep, fu: FaceuseId) -> Result<Option<RegionId>> {
    let shell: ShellId =
        eu::faceuse_shell(b, fu)?.ok_or(dangling(OP_FACEUSE, "face use without shell"))?;
    eu::shell_region(b, shell)
}

fn vertex_of_vertexuse(b: &Brep, vu: Option<crate::topology::VertexuseId>) -> Result<VertexId> {
    let vu = vu.ok_or(dangling(OP_FACEUSE, "null vertex use"))?;
    vertexuse::vertex(b, vu)?.ok_or(dangling(OP_FACEUSE, "vertex use without vertex"))
}

#[allow(clippy::too_many_arguments)] // upstream's MakeFace
fn make_face(
    brep: &mut Brep,
    region: RegionId,
    edge_loops: &[usize],
    edgeuses: &[EdgeuseId],
    loop_vertices: &[VertexId],
    surface: Surface,
    uv_domain: &Extent2d,
    orientation: OrientType,
) -> Result<(Option<RegionId>, Option<ShellId>, FaceId)> {
    // bVerifyRadialSectors = FALSE, bCalculateTrims = FALSE, cpOptUVCurves = NULL
    let _ = (
        brep,
        region,
        edge_loops,
        edgeuses,
        loop_vertices,
        surface,
        uv_domain,
        orientation,
    );
    // UNPORTED(topology.brep_p7): SmBrep::MakeFace
    unavailable("SmBrep::MakeFace", "topology.brep_p7")
}

/// The orientation of the primary edge use to look for in the radial sector of
/// `brep_edge` for the other-operand edge use `edgeuse` (opposite to it relative to
/// the vertices, since it will end up radial to it).
fn orientation_to_find(
    k: &Kernel,
    brep: &mut Brep,
    other: &mut Brep,
    edgeuse: EdgeuseId,
    other_edge: EdgeId,
    brep_edge: EdgeId,
    relations: &Relations,
) -> Result<OrientType> {
    let v = vertex_of_vertexuse(other, eu::vertexuse(other, edgeuse)?)?;
    let v2 = edge::other_vertex(other, other_edge, Some(v))?;
    let prim = eu::primary_edgeuse(brep, brep_edge)?
        .ok_or(dangling(OP_FACEUSE, "edge without primary edge use"))?;
    let vb = vertex_of_vertexuse(brep, eu::vertexuse(brep, prim)?)?;
    let vb2 = edge::other_vertex(brep, brep_edge, Some(vb))?;
    if Some(v) == v2 {
        // Here need to do additional work for closed loop stuff
        let tangent = |k: &Kernel, b: &mut Brep, use_: EdgeuseId, t: f64| -> Result<_> {
            eu3::normalized_evaluate(k, b, use_, t, false, true, true)?
                .first_deriv
                .ok_or(dangling(OP_FACEUSE, "no first derivative"))
        };
        let tan1 = tangent(k, other, edgeuse, 0.0)?;
        // Pointing from the vertex into the curve...
        let tan2 = -tangent(k, other, edgeuse, 1.0)?;
        let tan_ref = tangent(k, brep, prim, 0.0)?;
        let tan1 = tan1.unitize()?;
        let tan2 = tan2.unitize()?;
        let tan_ref = tan_ref.unitize()?;
        let orient = eu::orientation(brep, prim)?;
        if tan1.dot(tan_ref) > tan2.dot(tan_ref) {
            Ok(orient)
        } else {
            Ok(reverse_orientation(orient))
        }
    } else {
        let vb_rel = relations
            .vertices
            .primary(v)
            .ok_or(SolidError::null_pointer(OP_FACEUSE))?;
        if vb == vb_rel {
            eu::orientation(brep, prim)
        } else if vb2 == Some(vb_rel) {
            Ok(reverse_orientation(eu::orientation(brep, prim)?))
        } else {
            Err(SolidError::sm_err(OP_FACEUSE))
        }
    }
}

/// Make the face of an other-operand face use in the primary B-rep
/// (`MakeFaceuse`). The face use on the side of the other operand's infinite region
/// is used. A face that already has a primary mate is only notified of the copy.
/// Otherwise the containing primary region is `region` or found with
/// [`p1::find_brep_region_of_other_face`]; every vertex and edge of the face is
/// merged in ([`p1::merge_vertex`], [`merge_edge`]), the radial sector of each
/// primary edge chosen, the surface copied and the face made, related and notified.
/// `new_faces`, when given, is cleared and receives the new face.
///
/// # Errors
///
/// `SM_ERR` when no region contains the face, a vertex orientation cannot be matched
/// or no radial sector is found; `SM_ERR_NULL_POINTER` for a vertex without mate;
/// [`SolidError::InvalidInput`] for missing topology upstream dereferences or a zero
/// tangent (`Unitize`); those of the merges, evaluations, copies and `Relate`;
/// [`SolidError::Unported`] (module docs).
#[allow(clippy::too_many_arguments, clippy::too_many_lines)] // one upstream routine
pub fn make_faceuse(
    k: &Kernel,
    brep: &mut Brep,
    other: &mut Brep,
    faceuse: FaceuseId,
    region: Option<RegionId>,
    relations: &mut Relations,
    copy: &mut TopologyCopy<'_>,
    mut context: Option<&mut Context>,
    mut new_faces: Option<&mut Vec<FaceId>>,
) -> Result<MadeFaceuse> {
    // init output
    if let Some(list) = new_faces.as_deref_mut() {
        list.clear();
    }

    // If the Face of the input Faceuse bounds the other Brep's infinite region, then
    // set it up so that the face use is the one on the side of the infinite region.
    let mut faceuse = faceuse;
    let mate =
        eu::faceuse_mate(other, faceuse)?.ok_or(dangling(OP_FACEUSE, "face use without mate"))?;
    let other_inf = other.infinite_region();
    if shell_region_of(other, mate)? == other_inf && shell_region_of(other, faceuse)? != other_inf {
        faceuse = mate;
    }

    // Check to see if face is already mated to a Brep face
    let face_to_make =
        eu::faceuse_face(other, faceuse)?.ok_or(dangling(OP_FACEUSE, "face use without face"))?;
    if let Some(face) = relations.faces.primary(face_to_make) {
        // low work - notify of copy to propagate attributes and return
        let (source, dest) = (
            (copy.face_source_key)(face_to_make),
            (copy.face_destination_key)(face),
        );
        copy.notify(context, source, dest)?;
        return Ok(MadeFaceuse::default());
    }

    // when not given a Brep Region try to find the containing Brep Region
    let region = match region {
        Some(r) => Some(r),
        None => p1::find_brep_region_of_other_face(
            k,
            brep,
            other,
            faceuse,
            relations.edges.primary_map(),
        )?,
    };
    // error - unable to classify input OtherBrep faceuse into a Brep region
    let region = region.ok_or(SolidError::sm_err(OP_FACEUSE))?;

    // Now merge all unrelated input face->vertices into Brep's Region as shell
    // vertices, and pick up the radial sector edge uses.
    let mut face_edgeuses = Vec::new();
    let mut loop_vertices = Vec::new();
    let mut edge_loops = Vec::new();

    // for every Faceuse->Loopuse
    for lu in eu::faceuse_loopuses(other, faceuse, false, None)? {
        // get this loopuse's edgeuses
        let edgeuses = loopuse::edgeuses(other, lu, None, None)?;

        // load edge counts into sEdgeLoops array
        let vertex_loop = loopuse::is_vertex_loopuse(other, lu)?;
        if !vertex_loop {
            edge_loops.push(edgeuses.len());
        }

        // for every Loopuse->Vertexuse
        for (kk, vu) in loopuse::vertexuses(other, lu, None)?
            .into_iter()
            .enumerate()
        {
            let v = vertex_of_vertexuse(other, vu)?;

            // Merge the Vertex
            let brep_vertex = p1::merge_vertex(
                brep,
                other,
                Some(v),
                Some(region),
                relations,
                &mut copy.vertices,
                context.as_deref_mut(),
            )?;

            // load VertexLoopuses into sLoopVertices
            if kk == 0 && vertex_loop {
                loop_vertices.push(brep_vertex.ok_or(SolidError::null_pointer(OP_FACEUSE))?);
            }
        }

        // for every Loopuse->Edgeuse
        for eu_other in edgeuses {
            let other_edge =
                eu::edge(other, eu_other)?.ok_or(dangling(OP_FACEUSE, "edge use without edge"))?;

            // merge the edge
            let brep_edge = merge_edge(
                k,
                brep,
                other,
                Some(eu_other),
                Some(region),
                relations,
                copy,
                context.as_deref_mut(),
            )?
            .ok_or(dangling(OP_FACEUSE, "merged edge use without edge"))?;

            // Find the orientation of the Edgeuse to look for in the radial sector.
            let orient =
                orientation_to_find(k, brep, other, eu_other, other_edge, brep_edge, relations)?;

            // classify Edgeuse->Face against BrepEdge
            let found = match edge::find_radial_sector_between(
                k,
                brep,
                brep_edge,
                Some(&mut *other),
                eu_other,
                Some(region),
                orient,
                None,
            ) {
                Ok(sector) => sector.radial_edgeuse,
                Err(e @ (SolidError::Unported { .. } | SolidError::InvalidInput { .. })) => {
                    return Err(e)
                }
                Err(_) => None,
            };
            // add the classified radial sector's Edgeuse to the sFaceEdgeuses list
            face_edgeuses.push(found.ok_or(SolidError::sm_err(OP_FACEUSE))?);
        }
    }

    // Copy the surface, when possible as an analytic surface (its Brep mate is always
    // NULL; module docs)
    let face_record = other.face(face_to_make)?;
    let uv_domain = face_record.uv_domain;
    let new_surface = face_record
        .surface
        .as_ref()
        .ok_or(dangling(OP_FACEUSE, "face without surface"))?
        .copy_and_add_analytics(k)?;

    // If we have the case of a faceuse bordering an infinite region in the other
    // brep being merged into the infinite region in this brep, temporarily set Brep's
    // infinite region to NULL.
    let brep_inf = brep.infinite_region();
    let mate =
        eu::faceuse_mate(other, faceuse)?.ok_or(dangling(OP_FACEUSE, "face use without mate"))?;
    let saved_infinite_region = if Some(region) == brep_inf
        && shell_region_of(other, faceuse)? == other_inf
        && shell_region_of(other, mate)? != other_inf
    {
        brep.set_infinite_region_record(None);
        brep_inf
    } else {
        None
    };

    // make the new Brep Face within the classified Brep region
    let orientation = other.faceuse(faceuse)?.orientation;
    let made = make_face(
        brep,
        region,
        &edge_loops,
        &face_edgeuses,
        &loop_vertices,
        new_surface,
        &uv_domain,
        orientation,
    );
    if saved_infinite_region.is_some() {
        brep.set_infinite_region_record(saved_infinite_region);
    }
    let (new_region, _new_shell, new_face) = made?;

    relations.relate_faces(brep, new_face, other, face_to_make)?;
    let (source, dest) = (
        (copy.face_source_key)(face_to_make),
        (copy.face_destination_key)(new_face),
    );
    copy.notify(context, source, dest)?;

    // set output
    if let Some(list) = new_faces {
        list.push(new_face);
    }
    Ok(MadeFaceuse {
        old_region: Some(region),
        new_region,
    })
}
