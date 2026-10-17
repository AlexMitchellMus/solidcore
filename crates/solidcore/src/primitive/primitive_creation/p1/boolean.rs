//! Planar boolean classification, copying and cleanup after intersection insertion.
//!
//! Upstream: `source/SMLib/src/SmPrimitiveCreation.cpp` (`Boolean2D`).
//!
//! | Rust | Upstream |
//! |------|----------|
//! | `finish` (private) | `SmPrimitiveCreation::Boolean2D`, after `IntersectInsertRelate` |
//!
//! Deviations: marks are caller-owned records, separately for each arena. Both
//! operands use the passed Kernel. On success the result owns the primary arena;
//! the two consumed input slots become empty arenas. No recursion.

use super::primitives::{detailed_face, stitch_faces, Boolean2d};
use crate::error::{ignore_status, Result, SolidError};
use crate::geom::context::{MarkLock, MarkType};
use crate::kernel::Kernel;
use crate::nurbs::b_spline_surface_n_lib::p2::OrientType;
use crate::topology::edge;
use crate::topology::edgeuse::p1 as eu;
use crate::topology::topology_traverser::{Entity, Marks};
use crate::topology::{Brep, EdgeId, FaceId, VertexId};
use std::collections::BTreeMap;

#[derive(Default)]
struct PlanarIntersections {
    edges: Vec<(EdgeId, EdgeId)>,
    vertices: Vec<(VertexId, VertexId)>,
    to_other: BTreeMap<FaceId, FaceId>,
    to_primary: BTreeMap<FaceId, FaceId>,
}
impl PlanarIntersections {
    fn relate(&mut self, a: FaceId, b: FaceId) -> Result<()> {
        if self.to_other.get(&a).is_some_and(|v| *v != b)
            || self.to_primary.get(&b).is_some_and(|v| *v != a)
        {
            return Err(SolidError::sm_err(
                "primitive::primitive_creation::p1::boolean_2d::relate",
            ));
        }
        self.to_other.insert(a, b);
        self.to_primary.insert(b, a);
        Ok(())
    }
}
fn intersect_insert_relate(
    _k: &mut Kernel,
    _a: &mut Brep,
    _b: &mut Brep,
) -> Result<PlanarIntersections> {
    // UNPORTED(boolean.topology_intersector_p2): SmTopologyIntersector::IntersectInsertRelate
    Err(SolidError::unported(
        "SmTopologyIntersector::IntersectInsertRelate",
        "boolean.topology_intersector_p2",
    ))
}
fn delete_faces(
    _k: &Kernel,
    _a: &mut Brep,
    _intersections: &mut PlanarIntersections,
    _faces: &[FaceId],
) -> Result<()> {
    // UNPORTED(boolean.topology_intersector_p2): SmTopologyIntersector::DeleteFaces
    Err(SolidError::unported(
        "SmTopologyIntersector::DeleteFaces",
        "boolean.topology_intersector_p2",
    ))
}
fn remove_topology(_k: &Kernel, _a: &mut Brep, _marks: &mut Marks, _mark: MarkType) -> Result<()> {
    // UNPORTED(topology.brep_p5): SmBrep::RemoveTopologicalEdgesAndVertices
    Err(SolidError::unported(
        "SmBrep::RemoveTopologicalEdgesAndVertices",
        "topology.brep_p5",
    ))
}

pub(super) fn finish(
    k: &mut Kernel,
    a: &mut Brep,
    b: &mut Brep,
    first: FaceId,
    operation: Boolean2d,
    result: &mut Option<Brep>,
    marks: (&mut Marks, &mut Marks),
) -> Result<()> {
    let mut intersections = intersect_insert_relate(k, a, b)?;
    let mut lock = MarkLock::new_pair(core::slice::from_mut(&mut k.context), 0, 0)?;
    let mark_a = lock.mark_type1();
    let mark_b = lock.mark_type2();
    let (marks_a, marks_b) = marks;
    let outcome = (|| {
        for (primary, other) in &intersections.edges {
            marks_a.mark(k, Entity::Edge(*primary), mark_a);
            marks_b.mark(k, Entity::Edge(*other), mark_b);
        }
        for (primary, other) in &intersections.vertices {
            marks_a.mark(k, Entity::Vertex(*primary), mark_a);
            marks_b.mark(k, Entity::Vertex(*other), mark_b);
        }
        // The primary side propagates a failed radial-sector status; the other
        // side skips it. Mark each face BEFORE the radial query, as upstream.
        for (primary, other) in intersections.edges.clone() {
            for edgeuse in edge::edgeuses(a, primary, None, None)?
                .into_iter()
                .step_by(2)
            {
                let face = eu::face(a, edgeuse)?
                    .ok_or_else(|| SolidError::null_pointer("primitive::boolean_2d"))?;
                if marks_a.is_marked(k, Entity::Face(face), mark_a) {
                    continue;
                }
                marks_a.mark(k, Entity::Face(face), mark_a);
                let sector = edge::find_radial_sector_between(
                    k,
                    b,
                    other,
                    a,
                    edgeuse,
                    None,
                    OrientType::Same,
                    None,
                )?;
                if let Some(faceuse) = sector.coincident_faceuse {
                    let mate = eu::faceuse_face(b, faceuse)?
                        .ok_or_else(|| SolidError::null_pointer("primitive::boolean_2d"))?;
                    marks_b.mark(k, Entity::Face(mate), mark_b);
                    ignore_status(intersections.relate(face, mate))?;
                }
            }
            for edgeuse in edge::edgeuses(b, other, None, None)?.into_iter().step_by(2) {
                let face = eu::face(b, edgeuse)?
                    .ok_or_else(|| SolidError::null_pointer("primitive::boolean_2d"))?;
                if marks_b.is_marked(k, Entity::Face(face), mark_b) {
                    continue;
                }
                marks_b.mark(k, Entity::Face(face), mark_b);
                let sector = match edge::find_radial_sector_between(
                    k,
                    a,
                    primary,
                    b,
                    edgeuse,
                    None,
                    OrientType::Same,
                    None,
                ) {
                    Ok(value) => value,
                    Err(e) if e.is_unported() => return Err(e),
                    Err(_) => continue,
                };
                if let Some(faceuse) = sector.coincident_faceuse {
                    let mate = eu::faceuse_face(a, faceuse)?
                        .ok_or_else(|| SolidError::null_pointer("primitive::boolean_2d"))?;
                    marks_a.mark(k, Entity::Face(mate), mark_a);
                    intersections.relate(mate, face)?;
                }
            }
        }
        let mut delete = Vec::new();
        for face in crate::topology::brep::p1::faces(a, None)? {
            let common = intersections.to_other.contains_key(&face);
            if matches!(operation, Boolean2d::ExclusiveOr | Boolean2d::Difference) && common
                || operation == Boolean2d::Intersection && !common
            {
                delete.push(face);
            }
        }
        // Keep the original A faces until every B classification/copy is done.
        let mut covering = first;
        for face in crate::topology::brep::p1::faces(b, None)? {
            if matches!(
                operation,
                Boolean2d::ExclusiveOr | Boolean2d::Union | Boolean2d::Merge
            ) && !intersections.to_primary.contains_key(&face)
            {
                let boundary = crate::topology::face::p2::create_curves_from_face(
                    k,
                    b,
                    face,
                    OrientType::Same,
                    true,
                    false,
                )?;
                let original = a
                    .face(covering)?
                    .surface
                    .as_ref()
                    .ok_or_else(|| SolidError::null_pointer("primitive::boolean_2d"))?;
                let surface = if original.owner().is_some() {
                    original.copy()?
                } else {
                    original.clone_for_evaluation()?
                };
                let domain = surface.natural_uv_domain();
                covering = detailed_face(k, a, a.infinite_region(), boundary, surface, domain)?;
            }
        }
        ignore_status(delete_faces(k, a, &mut intersections, &delete))?;
        if operation != Boolean2d::Merge {
            stitch_faces(k, a, a.tolerance())?;
            remove_topology(k, a, marks_a, mark_a)?;
        }
        *result = Some(core::mem::take(a));
        *b = Brep::new();
        Ok(())
    })();
    lock.release(core::slice::from_mut(&mut k.context))?;
    outcome
}
