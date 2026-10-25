//! Original-face selection using the session's actual topology marks.

use super::{invalid, FilletBody, FilletObject, FilletSession, SessionObjects};
use crate::blend::fillet_executive::p3::{FilletExecutive, OriginalSurfaceTrimAccess};
use crate::blend::fillet_geom::p6::{
    TrimOriginalBody, TrimOriginalSurfacesEnv, TrimOriginalSurfacesMarks,
};
use crate::error::{EntityRef, Result};
use crate::geom::context::{MarkLock, MarkType, ALL_MARKS};
use crate::kernel::Kernel;
use crate::topology::{
    brep::{
        p1,
        p3::{seed_records, topology_ref},
    },
    topology_traverser::Entity,
    FaceId,
};

/// Face selections before the subsequent native deletion and merge operations.
#[derive(Debug, Default)]
pub struct OriginalSurfaceSelection {
    /// Intermediate faces selected for retention, in discovery order.
    pub kept: Vec<FaceId>,
    /// Target faces selected for deletion, in discovery order.
    pub deleted: Vec<FaceId>,
}

struct TrimMarks<'a> {
    records: &'a mut std::collections::BTreeMap<
        crate::geom::attribute::ObjectKey,
        crate::topology::topology::Topology,
    >,
    fillet_serial: u64,
    target_serial: u64,
}
impl TrimMarks<'_> {
    fn key(&self, body: TrimOriginalBody, entity: Entity) -> crate::geom::attribute::ObjectKey {
        let entity = match entity {
            Entity::Face(id) => EntityRef::Face(id),
            Entity::Edge(id) => EntityRef::Edge(id),
            Entity::Vertex(id) => EntityRef::Vertex(id),
            Entity::Faceuse(id) => EntityRef::Faceuse(id),
            Entity::Vertexuse(id) => EntityRef::Vertexuse(id),
            Entity::Shell(id) => EntityRef::Shell(id),
        };
        topology_ref(
            match body {
                TrimOriginalBody::Fillet => self.fillet_serial,
                TrimOriginalBody::Target => self.target_serial,
            },
            entity,
        )
        .key
    }
}
impl TrimOriginalSurfacesMarks for TrimMarks<'_> {
    fn is_marked(
        &self,
        k: &Kernel,
        body: TrimOriginalBody,
        entity: Entity,
        mark: MarkType,
    ) -> bool {
        self.records
            .get(&self.key(body, entity))
            .is_some_and(|topology| topology.is_marked(&k.context, mark))
    }
    fn mark(&mut self, k: &Kernel, body: TrimOriginalBody, entity: Entity, mark: MarkType) {
        let key = self.key(body, entity);
        self.records
            .entry(key)
            .or_default()
            .mark(Some(&k.context), mark);
    }
}

struct TrimAccess<'a, 'target> {
    objects: &'a mut SessionObjects<'target>,
    locks: &'a [MarkLock],
    fillet_serial: u64,
}
impl OriginalSurfaceTrimAccess<FilletObject, FaceId> for TrimAccess<'_, '_> {
    fn original_trim_geometries(&self, solver: FilletObject) -> Result<Vec<FilletObject>> {
        let FilletObject::Solver(id) = solver else {
            return Err(invalid("trim solver is not an actual solver"));
        };
        Ok(self
            .objects
            .solver(id)?
            .base()
            .fillet_geoms()
            .into_iter()
            .map(FilletObject::Geometry)
            .collect())
    }
    fn original_trim_geometry(
        &mut self,
        k: &mut Kernel,
        geom: FilletObject,
        kept: &mut Vec<FaceId>,
        deleted: &mut Vec<FaceId>,
        executive: &mut FilletExecutive<FilletObject>,
    ) -> Result<()> {
        let FilletObject::Geometry(geom) = geom else {
            return Err(invalid("trim geometry is not actual geometry"));
        };
        let id = self
            .objects
            .pseudo
            .geom(geom)?
            .fillet_solver()
            .ok_or_else(|| invalid("trim geometry has no actual solver"))?;
        let solver = self
            .objects
            .solvers
            .get(id.0 as usize)
            .and_then(Option::as_ref)
            .ok_or_else(|| invalid("trim solver is not live"))?
            .base();
        let fillet = self
            .objects
            .intermediate
            .as_mut()
            .ok_or_else(|| invalid("intermediate body is not allocated"))?;
        let mate =
            |edge| match executive
                .intersector
                .brep_mate(Some(FilletObject::IntermediateTopology(EntityRef::Edge(
                    edge,
                )))) {
                Some(FilletObject::Topology(FilletBody::Target, EntityRef::Edge(edge))) => {
                    Some(edge)
                }
                _ => None,
            };
        let key = |edge| topology_ref(self.fillet_serial, EntityRef::Edge(edge));
        let env = TrimOriginalSurfacesEnv {
            brep_mate: &mate,
            fillet_edge_key: &key,
            target_context: Some(0),
        };
        let mut marks = TrimMarks {
            records: &mut self.objects.notification_records,
            fillet_serial: self.fillet_serial,
            target_serial: self.objects.target_serial,
        };
        self.objects.pseudo.trim_original_surfaces_with_solver(
            k,
            geom,
            fillet,
            self.objects.target,
            &mut marks,
            kept,
            deleted,
            self.locks,
            &env,
            id,
            solver,
        )
    }
}

impl FilletSession<'_> {
    /// Select original faces using actual common-mate marks and retained solvers.
    /// This component ends its mark scope before returning; it does not delete faces
    /// or complete native trimming and merging.
    ///
    /// # Errors
    /// Missing actual body/context authority or the first selection failure.
    /// Both bodies must have identities issued by this same kernel context;
    /// a distinct target context needs the native second-lock capability.
    /// Previously written topology marks remain; locks and editing are restored.
    pub fn select_original_fillet_surfaces(
        &mut self,
        k: &mut Kernel,
        selection: &mut OriginalSurfaceSelection,
    ) -> Result<()> {
        let editing = self.objects.target.editing_enabled();
        self.objects.target.set_editing_enabled(true);
        let result = self.with_original_trim_marks(k, |session, k, locks, serial| {
            session.executive.select_original_surfaces(
                k,
                &mut TrimAccess {
                    objects: &mut session.objects,
                    locks,
                    fillet_serial: serial,
                },
                &mut selection.kept,
                &mut selection.deleted,
            )
        });
        self.objects.target.set_editing_enabled(editing);
        result
    }

    // The lock scope can also enclose later native trim continuation without
    // advancing marks or copying the authoritative record map between components.
    fn with_original_trim_marks<T>(
        &mut self,
        k: &mut Kernel,
        body: impl FnOnce(&mut Self, &mut Kernel, &[MarkLock], u64) -> Result<T>,
    ) -> Result<T> {
        let fillet = self
            .objects
            .intermediate
            .as_mut()
            .ok_or_else(|| invalid("intermediate body is not allocated"))?;
        let serial = k
            .context
            .issued_brep_serial(fillet.notify_identity_slot())
            .ok_or_else(|| invalid("intermediate body belongs to another context"))?;
        let target_serial = k
            .context
            .issued_brep_serial(self.objects.target.notify_identity_slot())
            .ok_or_else(|| invalid("target body belongs to another context"))?;
        if target_serial != self.objects.target_serial || serial == target_serial {
            return Err(invalid("trim body identity is not authoritative"));
        }
        seed_records(&mut self.objects.notification_records, fillet, &|entity| {
            topology_ref(serial, entity)
        });
        seed_records(
            &mut self.objects.notification_records,
            self.objects.target,
            &|entity| topology_ref(target_serial, entity),
        );
        let mut lock = MarkLock::new(std::slice::from_mut(&mut k.context), 0, ALL_MARKS)?;
        let result = (|| {
            self.mark_common_trim_topology(k, serial, lock.mark_type())?;
            body(self, k, std::slice::from_ref(&lock), serial)
        })();
        let released = lock.release(std::slice::from_mut(&mut k.context));
        match result {
            Err(error) => Err(error),
            Ok(value) => released.map(|()| value),
        }
    }

    fn mark_common_trim_topology(&mut self, k: &Kernel, serial: u64, mark: MarkType) -> Result<()> {
        let fillet = self
            .objects
            .intermediate
            .as_ref()
            .ok_or_else(|| invalid("intermediate body is not allocated"))?;
        let mut marks = TrimMarks {
            records: &mut self.objects.notification_records,
            fillet_serial: serial,
            target_serial: self.objects.target_serial,
        };
        for vertex in p1::vertices(self.objects.target, None)? {
            if let Some(mate) = self
                .executive
                .intersector
                .other_mate(Some(FilletObject::Topology(
                    FilletBody::Target,
                    EntityRef::Vertex(vertex),
                )))
            {
                let FilletObject::IntermediateTopology(EntityRef::Vertex(other)) = mate else {
                    return Err(invalid(
                        "common vertex mate is not actual intermediate topology",
                    ));
                };
                fillet.vertex(other)?;
                marks.mark(k, TrimOriginalBody::Fillet, Entity::Vertex(other), mark);
                marks.mark(k, TrimOriginalBody::Target, Entity::Vertex(vertex), mark);
            }
        }
        for edge in p1::edges(self.objects.target, None)? {
            if let Some(mate) = self
                .executive
                .intersector
                .other_mate(Some(FilletObject::Topology(
                    FilletBody::Target,
                    EntityRef::Edge(edge),
                )))
            {
                let FilletObject::IntermediateTopology(EntityRef::Edge(other)) = mate else {
                    return Err(invalid(
                        "common edge mate is not actual intermediate topology",
                    ));
                };
                fillet.edge(other)?;
                marks.mark(k, TrimOriginalBody::Target, Entity::Edge(edge), mark);
                marks.mark(k, TrimOriginalBody::Fillet, Entity::Edge(other), mark);
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::Point3;

    #[test]
    fn actual_box_selection_uses_retained_mates_and_releases_marks() {
        let mut k = Kernel::default();
        let mut target = crate::primitive::api_primitives::create_box(
            &mut k,
            Point3::new(0., 0., 0.),
            10.,
            20.,
            30.,
        )
        .unwrap();
        let edge = target.edges().ids().next().unwrap();
        let mut session = FilletSession::new(&mut k, &mut target).unwrap();
        session
            .prepare_constant_radius_edges(&mut k, &[edge], 0.5, 1, 1, 1.0)
            .unwrap();
        session.create_corners(&k, None).unwrap();
        session
            .with_filleting_context(&mut k, |session, k| {
                session.create_fillet_brep(k, None)?;
                session.shrink_fillet_faces(k)?;
                session.insert_fillet_intersection_topology(k)?;
                let editing = session.target().editing_enabled();
                let mut selected = OriginalSurfaceSelection::default();
                session.select_original_fillet_surfaces(k, &mut selected)?;
                assert_eq!(session.target().editing_enabled(), editing);
                assert!(!selected.kept.is_empty());
                for face in &selected.kept {
                    session.intermediate_brep().unwrap().face(*face)?;
                }
                for face in &selected.deleted {
                    session.target().face(*face)?;
                }
                for mark in [MarkType::Mark, MarkType::Mark2, MarkType::Mark3] {
                    assert!(!k.context.is_mark_locked(mark));
                }
                Ok(())
            })
            .unwrap();
    }

    #[test]
    fn scoped_error_retains_body_qualified_marks_and_unlocks() {
        let mut k = Kernel::default();
        let mut target = crate::primitive::api_primitives::create_box(
            &mut k,
            Point3::new(0., 0., 0.),
            1.,
            1.,
            1.,
        )
        .unwrap();
        let mut fillet = crate::primitive::api_primitives::create_box(
            &mut k,
            Point3::new(2., 0., 0.),
            1.,
            1.,
            1.,
        )
        .unwrap();
        crate::topology::brep::p3::notify_serial(&mut k.context, &mut fillet).unwrap();
        let mut session = FilletSession::new(&mut k, &mut target).unwrap();
        session.objects.intermediate = Some(fillet);
        let target_face = session.objects.target.faces().ids().next().unwrap();
        let existing_key =
            topology_ref(session.objects.target_serial, EntityRef::Face(target_face)).key;
        let existing = session
            .objects
            .notification_records
            .entry(existing_key)
            .or_default();
        existing.set_user_index1(73);
        existing.mark(Some(&k.context), MarkType::Mark2);
        let prior_mark = existing.get_mark(MarkType::Mark2);
        let mut captured = None;
        let result: Result<()> =
            session.with_original_trim_marks(&mut k, |session, k, locks, serial| {
                let mark = locks[0].mark_type();
                assert!(k.context.is_mark_locked(mark));
                let face = session.objects.target.faces().ids().next().unwrap();
                session.objects.intermediate.as_ref().unwrap().face(face)?;
                let mut marks = TrimMarks {
                    records: &mut session.objects.notification_records,
                    fillet_serial: serial,
                    target_serial: session.objects.target_serial,
                };
                marks.mark(k, TrimOriginalBody::Fillet, Entity::Face(face), mark);
                assert!(!marks.is_marked(k, TrimOriginalBody::Target, Entity::Face(face), mark));
                captured = Some((
                    marks.key(TrimOriginalBody::Fillet, Entity::Face(face)),
                    mark,
                ));
                Err(invalid("defined continuation failure"))
            });
        assert!(result.is_err());
        let (key, mark) = captured.unwrap();
        assert!(!k.context.is_mark_locked(mark));
        assert!(session.objects.notification_records[&key].is_marked(&k.context, mark));
        let retained = &session.objects.notification_records[&existing_key];
        assert_eq!(retained.user_index1(), 73);
        assert_eq!(retained.get_mark(MarkType::Mark2), prior_mark);
    }
    #[test]
    fn foreign_context_rejects_before_lock_or_mark_mutation() {
        let mut k = Kernel::default();
        let mut foreign = Kernel::default();
        let mut target = crate::primitive::api_primitives::create_box(
            &mut k,
            Point3::new(0., 0., 0.),
            1.,
            1.,
            1.,
        )
        .unwrap();
        let mut fillet = crate::primitive::api_primitives::create_box(
            &mut foreign,
            Point3::new(2., 0., 0.),
            1.,
            1.,
            1.,
        )
        .unwrap();
        crate::topology::brep::p3::notify_serial(&mut foreign.context, &mut fillet).unwrap();
        let mut session = FilletSession::new(&mut k, &mut target).unwrap();
        session.objects.intermediate = Some(fillet);
        let count = session.objects.notification_records.len();
        let mut selected = OriginalSurfaceSelection::default();
        assert!(session
            .select_original_fillet_surfaces(&mut k, &mut selected)
            .is_err());
        assert_eq!(session.objects.notification_records.len(), count);
        assert!(selected.kept.is_empty() && selected.deleted.is_empty());
        for mark in [MarkType::Mark, MarkType::Mark2, MarkType::Mark3] {
            assert!(!k.context.is_mark_locked(mark));
        }
    }
}
