//! Borrowed target/pseudo geometry for the executive's attached rail classifications.

use crate::blend::fillet_classification::{FilletBody, FilletClassObject};
use crate::curve::curve::Curve;
use crate::error::EntityRef;
use crate::geom::object::ObjectRef;
use crate::intersect::point_class::ClassGeometry;
use crate::math::Point3;
use crate::query::face::p2::CurveOnClassifyObjects;
use crate::surface::surface::Surface;
use crate::surface::volume::Volume;
use crate::topology::{Brep, BrepId, EdgeId, EdgeuseId, FaceId, LoopId, VertexId, VertexuseId};

/// Each body is borrowed from the same live session during the classification loan.
/// No indexed/legacy handles are accepted by this production provider.
pub(super) struct Objects<'a> {
    pub target: &'a Brep,
    pub pseudo: &'a Brep,
}

impl Objects<'_> {
    fn body(&self, body: FilletBody) -> &Brep {
        match body {
            FilletBody::Target => self.target,
            FilletBody::Pseudo => self.pseudo,
        }
    }

    fn body_of(&self, brep: &Brep) -> Option<FilletBody> {
        if std::ptr::eq(brep, self.target) {
            Some(FilletBody::Target)
        } else if std::ptr::eq(brep, self.pseudo) {
            Some(FilletBody::Pseudo)
        } else {
            None
        }
    }
}

impl ClassGeometry<FilletClassObject> for Objects<'_> {
    fn curve(&self, object: FilletClassObject) -> Option<&Curve> {
        match object {
            FilletClassObject::Curve(body, edge) => self.body(body).edge(edge).ok()?.curve.as_ref(),
            FilletClassObject::TrimCurve(body, eu) => {
                self.body(body).edgeuse(eu).ok()?.uv_trim_curve.as_ref()
            }
            _ => None,
        }
    }
    fn surface(&self, object: FilletClassObject) -> Option<&Surface> {
        let FilletClassObject::Surface(body, face) = object else {
            return None;
        };
        self.body(body).face(face).ok()?.surface.as_ref()
    }
    fn volume(&self, _: FilletClassObject) -> Option<&Volume> {
        None
    }
    fn point(&self, object: FilletClassObject) -> Option<Point3> {
        let FilletClassObject::Topology(body, EntityRef::Vertex(vertex)) = object else {
            return None;
        };
        Some(self.body(body).vertex(vertex).ok()?.point)
    }
    fn classification_vertex_point(&self, object: FilletClassObject) -> Option<Point3> {
        self.point(object)
    }
    fn topology_zone_tol_3d(&self, object: FilletClassObject) -> Option<f64> {
        let FilletClassObject::Topology(body, entity) = object else {
            return None;
        };
        let brep = self.body(body);
        match entity {
            EntityRef::Vertex(id) => Some(brep.vertex(id).ok()?.tolerance),
            EntityRef::Edge(id) => Some(brep.edge(id).ok()?.tolerance),
            EntityRef::Face(id) => Some(brep.face(id).ok()?.tolerance),
            EntityRef::Edgeuse(id) => Some(brep.edge(brep.edgeuse(id).ok()?.edge?).ok()?.tolerance),
            _ => None,
        }
    }
    fn topology_entity(&self, object: FilletClassObject) -> Option<(&Brep, EntityRef)> {
        let FilletClassObject::Topology(body, entity) = object else {
            return None;
        };
        Some((self.body(body), entity))
    }
    fn topology_handle(&self, entity: EntityRef) -> Option<FilletClassObject> {
        Some(FilletClassObject::Topology(FilletBody::Target, entity))
    }
    fn topology_handle_in(&self, brep: &Brep, entity: EntityRef) -> Option<FilletClassObject> {
        Some(FilletClassObject::Topology(self.body_of(brep)?, entity))
    }
    fn brep_id(&self, brep: &Brep) -> Option<BrepId> {
        self.body_of(brep)?;
        // These are the actual distinct serials issued by the kernel constructor,
        // not local topology indices or synthesized offsets between body arenas.
        Some(BrepId::from_index(
            u32::try_from(brep.notify_serial()?).ok()?,
        ))
    }
    fn owner_entity(&self, owner: ObjectRef) -> Option<(&Brep, EntityRef)> {
        for brep in [self.target, self.pseudo] {
            let Some(serial) = brep.notify_serial() else {
                continue;
            };
            let Some(entity) = crate::topology::brep::p3::topology_entity(serial, owner.key) else {
                continue;
            };
            if crate::topology::brep::p3::topology_ref(serial, entity) == owner {
                return Some((brep, entity));
            }
        }
        None
    }
}

impl CurveOnClassifyObjects<FilletClassObject> for Objects<'_> {
    fn face(&self, id: FaceId) -> FilletClassObject {
        FilletClassObject::Topology(FilletBody::Target, EntityRef::Face(id))
    }
    fn vertex(&self, id: VertexId) -> FilletClassObject {
        FilletClassObject::Topology(FilletBody::Target, EntityRef::Vertex(id))
    }
    fn edge(&self, id: EdgeId) -> FilletClassObject {
        FilletClassObject::Topology(FilletBody::Target, EntityRef::Edge(id))
    }
    fn vertexuse(&self, id: VertexuseId) -> FilletClassObject {
        FilletClassObject::Topology(FilletBody::Target, EntityRef::Vertexuse(id))
    }
    fn edgeuse(&self, id: EdgeuseId) -> FilletClassObject {
        FilletClassObject::Topology(FilletBody::Target, EntityRef::Edgeuse(id))
    }
    fn edge_curve(&self, id: EdgeId) -> FilletClassObject {
        FilletClassObject::Curve(FilletBody::Target, id)
    }
    fn trim_curve(&self, id: EdgeuseId) -> FilletClassObject {
        FilletClassObject::TrimCurve(FilletBody::Target, id)
    }
    fn face_surface(&self, id: FaceId) -> Option<FilletClassObject> {
        self.target.face(id).ok()?.surface.as_ref()?;
        Some(FilletClassObject::Surface(FilletBody::Target, id))
    }
    fn loop_id(&self, object: FilletClassObject) -> Option<LoopId> {
        match object {
            FilletClassObject::Topology(FilletBody::Target, EntityRef::Loop(id)) => Some(id),
            _ => None,
        }
    }
}

/// Live retained corner resolution; geometry borrows exist only inside `classify`.
pub(super) struct Context<'a> {
    pub corners: &'a [crate::blend::fillet_corner::p3::FilletCorner],
}

impl crate::blend::fillet_solver::p5::RailClassificationContext for Context<'_> {
    fn rail_curve(&self, edge: EdgeId) -> FilletClassObject {
        FilletClassObject::Curve(FilletBody::Pseudo, edge)
    }
    fn rail_mate_uv_curve(&self, pseudo: &Brep, edgeuse: EdgeuseId) -> Option<FilletClassObject> {
        pseudo.edgeuse(edgeuse).ok()?.uv_trim_curve.as_ref()?;
        Some(FilletClassObject::TrimCurve(FilletBody::Pseudo, edgeuse))
    }
    fn classify(
        &mut self,
        k: &crate::kernel::Kernel,
        target: &mut Brep,
        pseudo: &Brep,
        face: FaceId,
        cc: &mut crate::intersect::curve_class::p1::CurveClassification<FilletClassObject>,
    ) -> crate::error::Result<()> {
        let mut cache = crate::query::face::p2::CurveOnClassifyCache::default();
        let result = {
            let objects = Objects { target, pseudo };
            crate::query::face::p2::curve_on_classify_cached(
                k, target, face, &objects, true, cc, true, false, true, &mut cache,
            )
        };
        // Native const queries retain their lazy gap/UV edits, including mutations
        // made before a failed classification. Use the shared query synchronization.
        cache.apply_uv_trim_updates(target)?;
        result.map(|_| ())
    }
    fn corner_type(
        &self,
        id: crate::blend::fillet_geom::p1::FilletCornerId,
    ) -> crate::error::Result<crate::blend::fillet_geom::p1::FilletCornerType> {
        self.corners
            .get(id.0 as usize)
            .map(crate::blend::fillet_corner::p3::FilletCorner::corner_type)
            .ok_or_else(|| super::invalid("rollover corner is not live in this session"))
    }
}

/// Retained notification context for the actual installed recovery curve.
pub(super) struct Ownership<'a> {
    pub serial: u64,
    pub attributes: &'a mut crate::geom::attribute::AttributeStore,
    pub context: &'a mut crate::geom::context::Context,
}
impl<'a> Ownership<'a> {
    pub fn new(
        pseudo: &Brep,
        attributes: &'a mut crate::geom::attribute::AttributeStore,
        context: &'a mut crate::geom::context::Context,
    ) -> crate::error::Result<Self> {
        let serial = pseudo
            .notify_serial()
            .ok_or_else(|| super::invalid("pseudo B-rep has no kernel identity"))?;
        if context.issued_brep_serial(pseudo.notify_identity()) != Some(serial) {
            return Err(super::invalid(
                "rollover pseudo B-rep belongs to another context",
            ));
        }
        Ok(Self {
            serial,
            attributes,
            context,
        })
    }
}
impl crate::blend::fillet_solver::p1::RailOwnership for Ownership<'_> {
    fn edge_owner(&self, edge: EdgeId) -> crate::error::Result<(ObjectRef, Option<ObjectRef>)> {
        Ok((
            crate::topology::brep::p3::topology_ref(self.serial, EntityRef::Edge(edge)),
            // SM_NO_GET_BREP is NULL in the pinned build.
            None,
        ))
    }
    fn owner_changed(
        &mut self,
        _: crate::curve::curve::p2::OwnerChange,
    ) -> crate::error::Result<()> {
        // UNPORTED(blend.fillet_solver_p5): this legacy callback has no actual Curve identity.
        Err(crate::error::SolidError::unported(
            "SmFilletSolver::ProcessRollOver (actual curve notification context)",
            "blend.fillet_solver_p5",
        ))
    }
    fn set_curve_owner(
        &mut self,
        _edge: EdgeId,
        curve: &mut Curve,
        owner: ObjectRef,
        _: Option<ObjectRef>,
    ) -> crate::error::Result<()> {
        let old_owner = curve.owner();
        if old_owner != Some(owner) && curve.dim() == 3 {
            let type_code = curve.type_code();
            let key = self.context.geometry_key(curve.notify_identity_slot())?;
            let actual = ObjectRef { type_code, key };
            // SmCurve::SetOwner notifies while GetOwner still returns the old owner.
            crate::geom::object::report_owner_link(Some(self.context), actual.key, old_owner);
            crate::error::ignore_status(crate::geom::a_object::notify(
                self.attributes,
                Some(self.context),
                crate::geom::object::NotifyOperation::ChangeOwner,
                actual,
                Some(owner),
                None,
                old_owner,
            ))?;
        }
        curve.set_owner(Some(owner), None);
        if let Some(key) = curve.notify_key() {
            crate::geom::object::report_owner_link(Some(self.context), key, Some(owner));
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kernel::Kernel;
    use crate::primitive::api_primitives::create_box;

    #[test]
    fn equal_local_ids_keep_actual_body_and_geometry_identity() {
        let mut k = Kernel::default();
        let target = create_box(&mut k, Point3::new(0.0, 0.0, 0.0), 10.0, 20.0, 30.0).unwrap();
        let pseudo = create_box(&mut k, Point3::new(100.0, 200.0, 300.0), 4.0, 5.0, 6.0).unwrap();
        let target_edge = target.edges().ids().next().unwrap();
        let pseudo_edge = pseudo.edges().ids().next().unwrap();
        assert_eq!(target_edge, pseudo_edge);
        let objects = Objects {
            target: &target,
            pseudo: &pseudo,
        };
        let source = FilletClassObject::Curve(FilletBody::Pseudo, pseudo_edge);
        let boundary = objects.edge_curve(target_edge);
        assert_ne!(source, boundary);
        assert!(std::ptr::eq(
            objects.curve(source).unwrap(),
            pseudo.edge(pseudo_edge).unwrap().curve.as_ref().unwrap()
        ));
        assert!(std::ptr::eq(
            objects.curve(boundary).unwrap(),
            target.edge(target_edge).unwrap().curve.as_ref().unwrap()
        ));
        assert!(!std::ptr::eq(
            objects.curve(source).unwrap(),
            objects.curve(boundary).unwrap()
        ));
        assert!(objects
            .curve(FilletClassObject::Legacy(EntityRef::Edge(target_edge)))
            .is_none());
        let target_entity = EntityRef::Edge(target_edge);
        let scoped = objects.topology_handle_in(&pseudo, target_entity).unwrap();
        assert_eq!(
            scoped,
            FilletClassObject::Topology(FilletBody::Pseudo, target_entity)
        );
        assert!(std::ptr::eq(
            objects.topology_entity(scoped).unwrap().0,
            &raw const pseudo
        ));
        assert_ne!(objects.brep_id(&target), objects.brep_id(&pseudo));
        let actual_owner = target
            .edge(target_edge)
            .unwrap()
            .curve
            .as_ref()
            .unwrap()
            .owner()
            .unwrap();
        // Primitive construction retains its supplied producer key. This provider
        // does not infer a canonical body identity from that opaque key.
        assert!(objects.owner_entity(actual_owner).is_none());
        let canonical =
            crate::topology::brep::p3::topology_ref(target.notify_serial().unwrap(), target_entity);
        let (body, entity) = objects.owner_entity(canonical).unwrap();
        assert!(std::ptr::eq(body, &raw const target));
        assert_eq!(entity, target_entity);
        let foreign = create_box(&mut k, Point3::new(-10.0, 0.0, 0.0), 1.0, 1.0, 1.0).unwrap();
        let foreign_edge = foreign.edges().ids().next().unwrap();
        let foreign_owner = foreign
            .edge(foreign_edge)
            .unwrap()
            .curve
            .as_ref()
            .unwrap()
            .owner()
            .unwrap();
        assert!(objects.owner_entity(foreign_owner).is_none());
        let foreign_canonical = crate::topology::brep::p3::topology_ref(
            foreign.notify_serial().unwrap(),
            EntityRef::Edge(foreign_edge),
        );
        assert!(objects.owner_entity(foreign_canonical).is_none());
        assert!(objects
            .topology_handle_in(&foreign, EntityRef::Edge(foreign_edge))
            .is_none());
    }
}

#[cfg(test)]
mod rollover_notification_tests {
    use super::*;
    use crate::blend::fillet_geom::p1::{FilletBrep, FilletGeom, FilletGeomType, FilletVertex};
    use crate::blend::fillet_solver::p1::split_rail;
    use crate::geom::attribute::{AttributeStore, ObjectKey};
    use crate::geom::context::{Context as NotifyContext, NotifyCallback};
    use crate::geom::object::NotifyOperation;
    use crate::kernel::Kernel;
    use crate::math::Vec3;
    use std::collections::BTreeMap;
    use std::sync::{Arc, Mutex};

    #[derive(Debug, Clone, PartialEq)]
    struct Event {
        system: bool,
        caller: ObjectRef,
        data: [Option<ObjectRef>; 3],
        observed_owner: Option<ObjectRef>,
    }
    #[derive(Debug)]
    struct Callback {
        system: bool,
        events: Arc<Mutex<Vec<Event>>>,
        owners: BTreeMap<ObjectKey, Option<ObjectRef>>,
    }
    impl NotifyCallback for Callback {
        fn execute(
            &mut self,
            operation: NotifyOperation,
            caller: ObjectRef,
            data1: Option<ObjectRef>,
            data2: Option<ObjectRef>,
            data3: Option<ObjectRef>,
        ) -> crate::Result<()> {
            assert_eq!(operation, NotifyOperation::ChangeOwner);
            self.events.lock().unwrap().push(Event {
                system: self.system,
                caller,
                data: [data1, data2, data3],
                observed_owner: *self.owners.get(&caller.key).unwrap(),
            });
            // Native ignores callback failures, including dependency statuses.
            Err(crate::SolidError::unported("test callback", "test"))
        }
        fn owner_link(&mut self, object: ObjectKey, owner: Option<ObjectRef>) {
            self.owners.insert(object, owner);
        }
    }

    fn rail(
        k: &mut Kernel,
    ) -> (
        FilletBrep,
        crate::blend::fillet_geom::p1::FilletGeomId,
        EdgeId,
        ObjectRef,
    ) {
        let mut fb = FilletBrep::new();
        k.context
            .brep_serial(fb.brep_mut().notify_identity_slot())
            .unwrap();
        let vertices = [Vec3::ZERO, Vec3::X].map(|point| {
            let v = fb.brep_mut().add_vertex(point).unwrap();
            fb.insert_fillet_vertex(v, FilletVertex::new(k, None));
            v
        });
        let edge = fb
            .make_fillet_edge(Some(vertices[0]), Some(vertices[1]), None)
            .unwrap();
        let serial = fb.brep().notify_serial().unwrap();
        let owner = crate::topology::brep::p3::topology_ref(serial, EntityRef::Edge(edge));
        let mut curve = Curve::from(crate::curve::line::Line::from_points(
            Vec3::ZERO,
            Vec3::X,
            3,
        ));
        curve.set_owner(Some(owner), None);
        let original = ObjectRef {
            type_code: curve.type_code(),
            key: k
                .context
                .geometry_key(curve.notify_identity_slot())
                .unwrap(),
        };
        fb.brep_mut().edge_mut(edge).unwrap().curve = Some(curve);
        fb.brep_mut().edge_mut(edge).unwrap().tolerance = 0.125;
        let geom = fb.add_geom(FilletGeom::new(k, None, FilletGeomType::Default).unwrap());
        (fb, geom, edge, original)
    }

    #[test]
    fn rollover_split_notifies_actual_copy_before_owner_and_keeps_callback_failures_ignored() {
        let mut k = Kernel::default();
        let (mut fb, geom, old_edge, original) = rail(&mut k);
        let events = Arc::new(Mutex::new(Vec::new()));
        for system in [false, true] {
            let callback = Some(Box::new(Callback {
                system,
                events: events.clone(),
                owners: BTreeMap::new(),
            }) as Box<dyn NotifyCallback>);
            if system {
                k.context.set_sys_notify_callback(callback);
            } else {
                k.context.set_user_notify_callback(callback);
            }
        }
        let mut attributes = AttributeStore::default();
        let (_, new_edge) = NotifyContext::with_notification_context(&mut k, |k, context| {
            let mut ownership = Ownership::new(fb.brep(), &mut attributes, context).unwrap();
            split_rail(k, &mut fb, &mut ownership, geom, old_edge)
        })
        .unwrap()
        .unwrap();
        let installed = fb.brep().edge(new_edge).unwrap().curve.as_ref().unwrap();
        let actual = ObjectRef {
            type_code: installed.type_code(),
            key: installed.notify_key().unwrap(),
        };
        assert_ne!(actual, original);
        let owner = crate::topology::brep::p3::topology_ref(
            fb.brep().notify_serial().unwrap(),
            EntityRef::Edge(new_edge),
        );
        // Virtual Curve::copy constructs a fresh, unowned wrapper.
        let old_owner = None;
        assert_eq!(
            *events.lock().unwrap(),
            vec![
                Event {
                    system: false,
                    caller: actual,
                    data: [Some(owner), None, old_owner],
                    observed_owner: old_owner
                },
                Event {
                    system: true,
                    caller: actual,
                    data: [Some(owner), None, old_owner],
                    observed_owner: old_owner
                },
            ]
        );
        assert_eq!(installed.owner(), Some(owner));
        assert_eq!(fb.brep().edge(new_edge).unwrap().tolerance, 0.125);
        assert!(k.context.user_notify_callback().is_some());
        assert!(k.context.sys_notify_callback().is_some());
    }

    #[test]
    fn rollover_identity_refusal_retains_installed_copy_and_prior_topology_without_late_tolerance()
    {
        let mut k = Kernel::default();
        let (mut fb, geom, old_edge, _) = rail(&mut k);
        let old_vertices = fb.brep().vertices().len();
        let old_edges = fb.edge_list().len();
        let mut attributes = AttributeStore::default();
        let error = NotifyContext::with_notification_context(&mut k, |k, context| {
            let mut ownership = Ownership::new(fb.brep(), &mut attributes, context).unwrap();
            // A real active identity lease refuses issuance at the actual install.
            let lease = ownership.context.lend_notify().unwrap();
            let result = split_rail(k, &mut fb, &mut ownership, geom, old_edge);
            ownership.context.install_lease(lease);
            result
        })
        .unwrap()
        .unwrap_err();
        assert!(matches!(error, crate::SolidError::InvalidInput { .. }));
        assert_eq!(fb.brep().vertices().len(), old_vertices + 1);
        assert_eq!(fb.edge_list().len(), old_edges + 1);
        let new_edge = *fb.edge_list().last().unwrap();
        let installed = fb.brep().edge(new_edge).unwrap();
        assert!(installed.curve.is_some());
        assert_eq!(installed.tolerance, 0.0);
        assert_eq!(installed.curve.as_ref().unwrap().owner(), None);
    }

    #[test]
    fn rollover_refuses_foreign_domain_even_with_equal_numeric_serial() {
        let mut k = Kernel::default();
        let (fb, _, _, _) = rail(&mut k);
        let mut foreign = Kernel::default();
        let mut other = Brep::new();
        foreign
            .context
            .brep_serial(other.notify_identity_slot())
            .unwrap();
        assert_eq!(fb.brep().notify_serial(), other.notify_serial());
        let mut attributes = AttributeStore::default();
        assert!(Ownership::new(fb.brep(), &mut attributes, &mut foreign.context).is_err());
        assert!(attributes.attributes(ObjectKey(0)).is_empty());
    }
}
