//! Actual target/intermediate topology and pseudo-curve bindings for insertion.

use super::super::{FilletBody, FilletObject};
use super::curves::InsertionCurves;
use crate::blend::fillet_geom::p1::{FilletBrep, FilletGeomId, FilletSurfaceStorage};
use crate::curve::curve::Curve;
use crate::geom::attribute::AttributeStore;
use crate::geom::attribute::ObjectKey;
use crate::geom::context::Context;
use crate::geom::object::ObjectRef;
use crate::intersect::curve_class::p5::ClassGeometryMut;
use crate::intersect::point_class::ClassGeometry;
use crate::math::Point3;
use crate::surface::surface::Surface;
use crate::surface::volume::Volume;
use crate::topology::brep::p3::topology_ref;
use crate::topology::proto_topology::p1::{ProtoGeometry, ProtoTopology};
use crate::topology::topology::Topology;
use crate::topology::{Brep, BrepId, FaceId};
use crate::{EntityRef, Result, SolidError};
use std::collections::BTreeMap;

const OP: &str = "blend::fillet_executive::insertion_geometry";
fn invalid(reason: &'static str) -> SolidError {
    SolidError::invalid_input(OP, reason)
}

/// One actual body in the insertion operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum InsertionBody {
    /// Original target.
    Target,
    /// Retained pseudo topology and curves.
    Pseudo,
    /// Constructed intermediate fillet body.
    Intermediate,
}

/// An actual Surface instance, distinct from its owning Face slot.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct SurfaceBinding {
    role: InsertionBody,
    body: ObjectRef,
    face: FaceId,
    geometry: ObjectRef,
}
impl SurfaceBinding {
    /// Actual body role.
    #[must_use]
    pub const fn body_role(self) -> InsertionBody {
        self.role
    }
    /// Original actual Face slot.
    #[must_use]
    pub const fn face(self) -> FaceId {
        self.face
    }
    /// Issued identity of the actual geometry, not its Face.
    #[must_use]
    pub const fn identity(self) -> ObjectRef {
        self.geometry
    }
}

/// Mutable borrows of the three actual bodies; no topology or geometry is copied.
pub struct InsertionGeometry<'a> {
    target: &'a mut Brep,
    pseudo: &'a mut FilletBrep,
    intermediate: &'a mut Brep,
    curves: Option<&'a mut InsertionCurves>,
    notifications: Option<InsertionNotifications<'a>>,
}

/// Actual retained edit stores and the already lent notification context.
pub struct InsertionNotifications<'a> {
    /// Actual topology notification records.
    pub records: &'a mut BTreeMap<ObjectKey, Topology>,
    /// Actual attributes, including executive E and active face markers.
    pub attributes: &'a mut AttributeStore,
    /// Actual context lease, supplied by `with_lent_context`.
    pub context: &'a mut Context,
}
impl<'a> InsertionGeometry<'a> {
    /// Borrow the actual Session bodies in their established identity domains.
    /// # Errors
    /// Missing or colliding issued body identities.
    pub fn new(
        target: &'a mut Brep,
        pseudo: &'a mut FilletBrep,
        intermediate: &'a mut Brep,
    ) -> Result<Self> {
        let view = Self {
            target,
            pseudo,
            intermediate,
            curves: None,
            notifications: None,
        };
        let mut identities = Vec::new();
        for role in [
            InsertionBody::Target,
            InsertionBody::Pseudo,
            InsertionBody::Intermediate,
        ] {
            let identity = view.body_identity(role)?;
            if identities.contains(&identity) {
                return Err(invalid("insertion bodies have colliding identity domains"));
            }
            identities.push(identity);
        }
        Ok(view)
    }
    /// Supply actual operation-owned curves and retained notification stores.
    /// The context must be the real lease used by this Session's edits.
    #[must_use]
    pub fn with_notifications(
        mut self,
        curves: &'a mut InsertionCurves,
        notifications: InsertionNotifications<'a>,
    ) -> Self {
        self.curves = Some(curves);
        self.notifications = Some(notifications);
        self
    }
    /// Retain a real curve copy and return its actual issued operation handle.
    /// # Errors
    /// Missing actual edit stores or identity issuance failure.
    pub fn insert_curve_copy(&mut self, curve: Curve) -> Result<FilletObject> {
        let notifications = self
            .notifications
            .as_mut()
            .ok_or_else(|| invalid("missing insertion notification context"))?;
        let curves = self
            .curves
            .as_mut()
            .ok_or_else(|| invalid("missing insertion curve ownership"))?;
        Ok(FilletObject::InsertionCurve(
            curves.insert(notifications.context, curve)?,
        ))
    }
    /// Run one edit against the actual body and retained notification stores.
    /// # Errors
    /// Missing actual context or the original edit failure, with partial writes retained.
    pub fn with_edit<T>(
        &mut self,
        role: InsertionBody,
        edit: impl FnOnce(&mut Brep, &mut crate::topology::edge::NotifyEnv<'_>) -> Result<T>,
    ) -> Result<T> {
        let serial = self
            .body(role)
            .notify_serial()
            .ok_or_else(|| invalid("missing issued body identity"))?;
        let key_of = |entity| topology_ref(serial, entity);
        let body = match role {
            InsertionBody::Target => &mut *self.target,
            InsertionBody::Pseudo => self.pseudo.brep_mut(),
            InsertionBody::Intermediate => &mut *self.intermediate,
        };
        let notifications = self
            .notifications
            .as_mut()
            .ok_or_else(|| invalid("missing insertion notification context"))?;
        crate::topology::brep::p3::seed_records(notifications.records, body, &key_of);
        let mut env = crate::topology::edge::NotifyEnv {
            records: notifications.records,
            attributes: notifications.attributes,
            context: Some(notifications.context),
            key_of: &key_of,
        };
        edit(body, &mut env)
    }
    /// Native edge interval, including actual curve-domain fallback.
    #[must_use]
    pub fn edge_interval(&self, edge: FilletObject) -> Option<crate::geom::extent1d::Extent1d> {
        let (role, EntityRef::Edge(edge)) = Self::topology(edge)? else {
            return None;
        };
        crate::topology::edge::interval(self.body(role), edge).ok()
    }
    /// Current edge Curve getter; distinct from the Edge topology object.
    #[must_use]
    pub fn edge_curve(&self, edge: FilletObject) -> Option<Option<FilletObject>> {
        let (role, EntityRef::Edge(edge)) = Self::topology(edge)? else {
            return None;
        };
        Some(
            self.body(role)
                .edge(edge)
                .ok()?
                .curve
                .as_ref()
                .map(|_| match role {
                    InsertionBody::Target => FilletObject::Curve(FilletBody::Target, edge),
                    InsertionBody::Pseudo => FilletObject::Curve(FilletBody::Pseudo, edge),
                    InsertionBody::Intermediate => FilletObject::IntermediateCurve(edge),
                }),
        )
    }
    /// Current Face UV domain.
    #[must_use]
    pub fn face_uv_domain(&self, face: FilletObject) -> Option<crate::geom::extent2d::Extent2d> {
        let (role, EntityRef::Face(face)) = Self::topology(face)? else {
            return None;
        };
        crate::topology::face::p2::uv_domain(self.body(role), face).ok()
    }
    /// Write the actual Face UV domain through its existing setter.
    pub fn set_face_uv_domain(
        &mut self,
        face: FilletObject,
        domain: crate::geom::extent2d::Extent2d,
    ) -> Option<()> {
        let (role, EntityRef::Face(face)) = Self::topology(face)? else {
            return None;
        };
        crate::topology::face::p2::set_uv_domain(self.body_mut(role), face, &domain).ok()
    }
    /// Current Face Surface getter, distinct from a captured original binding.
    #[must_use]
    pub fn face_surface(&self, face: FilletObject) -> Option<Option<FilletObject>> {
        let (role, EntityRef::Face(face)) = Self::topology(face)? else {
            return None;
        };
        Some(
            self.body(role)
                .face(face)
                .ok()?
                .surface
                .as_ref()
                .map(|_| FilletObject::FaceSurface(role, face)),
        )
    }
    /// Find the actual edge use belonging to a Face in the same body.
    /// # Errors
    /// Invalid Edge/Face handles or topology walk failures.
    pub fn edge_edgeuse_of_face(
        &self,
        edge: FilletObject,
        face: FilletObject,
    ) -> Result<Option<FilletObject>> {
        let (role, EntityRef::Edge(edge)) =
            Self::topology(edge).ok_or_else(|| invalid("unknown Edge role"))?
        else {
            return Err(invalid("handle is not an Edge"));
        };
        let (face_role, EntityRef::Face(face)) =
            Self::topology(face).ok_or_else(|| invalid("unknown Face role"))?
        else {
            return Err(invalid("handle is not a Face"));
        };
        if role != face_role {
            return Ok(None);
        }
        Ok(
            crate::topology::edge::edgeuse_of_face(self.body(role), edge, Some(face))?
                .map(|eu| Self::handle(role, EntityRef::Edgeuse(eu))),
        )
    }
    /// Test connectivity through the actual Edge topology.
    /// # Errors
    /// Invalid handles or topology walk failures.
    pub fn edge_is_connected_to_face(
        &self,
        edge: FilletObject,
        face: FilletObject,
    ) -> Result<bool> {
        let (role, EntityRef::Edge(edge)) =
            Self::topology(edge).ok_or_else(|| invalid("unknown Edge role"))?
        else {
            return Err(invalid("handle is not an Edge"));
        };
        let (face_role, EntityRef::Face(face)) =
            Self::topology(face).ok_or_else(|| invalid("unknown Face role"))?
        else {
            return Err(invalid("handle is not a Face"));
        };
        if role != face_role {
            return Ok(false);
        }
        crate::topology::edge::is_connected_to_face(self.body(role), edge, Some(face))
    }
    /// Common edges of actual same-body vertices.
    /// # Errors
    /// Invalid/cross-body operands or topology walk failures.
    pub fn vertex_common_edges(
        &self,
        first: FilletObject,
        second: FilletObject,
    ) -> Result<Vec<FilletObject>> {
        let (role, EntityRef::Vertex(first)) =
            Self::topology(first).ok_or_else(|| invalid("unknown Vertex role"))?
        else {
            return Err(invalid("handle is not a Vertex"));
        };
        let (other, EntityRef::Vertex(second)) =
            Self::topology(second).ok_or_else(|| invalid("unknown Vertex role"))?
        else {
            return Err(invalid("handle is not a Vertex"));
        };
        if role != other {
            return Err(invalid("vertices belong to different bodies"));
        }
        Ok(
            crate::topology::vertex::common_edges(self.body(role), first, Some(second), None)?
                .into_iter()
                .map(|edge| Self::handle(role, EntityRef::Edge(edge)))
                .collect(),
        )
    }
    /// Calculate the actual (including partial) trim-derived UV domain.
    /// # Errors
    /// Invalid owning Face or a genuine missing lower capability.
    pub fn face_trim_uv_domain(
        &mut self,
        k: &crate::Kernel,
        face: FilletObject,
    ) -> Result<crate::geom::extent2d::Extent2d> {
        let (role, EntityRef::Face(face)) =
            Self::topology(face).ok_or_else(|| invalid("unknown Face role"))?
        else {
            return Err(invalid("handle is not a Face"));
        };
        let body = self.body_mut(role);
        body.face(face)?;
        let mut output = crate::topology::face::p1::UvDomainFromTrim {
            domain: crate::geom::extent2d::Extent2d::new(),
            max_uv_loop_gap: 0.0,
            second_max_uv_loop_gap: 0.0,
        };
        let result = crate::topology::face::p1::calculate_uv_domain_from_uv_trim_curves_into(
            k,
            body,
            face,
            true,
            false,
            false,
            &mut output,
        );
        if let Err(error) = result {
            if error.is_unported() {
                return Err(error);
            }
        }
        Ok(output.domain)
    }
    /// Classify UV on the actual Face using native ignored ordinary-status semantics.
    /// # Errors
    /// Invalid owning Face or a genuine missing lower capability.
    pub fn face_point_classify_uv(
        &mut self,
        k: &crate::Kernel,
        face: FilletObject,
        uv: crate::math::Point2,
        tolerance: f64,
    ) -> Result<crate::intersect::point_class::PointClassType> {
        let (role, EntityRef::Face(face)) =
            Self::topology(face).ok_or_else(|| invalid("unknown Face role"))?
        else {
            return Err(invalid("handle is not a Face"));
        };
        let body = self.body_mut(role);
        body.face(face)?;
        let mut output = crate::intersect::point_class::PointClassification::new(tolerance);
        let result = crate::query::face::p2::point_classify(
            k,
            body,
            face,
            uv,
            &mut output,
            crate::query::face::p2::PointClassifyOptions {
                src_zone_tol_3d: tolerance,
                do_boundary_intersections: true,
                use_uv_classification: true,
                ..crate::query::face::p2::PointClassifyOptions::default()
            },
        );
        if let Err(error) = result {
            if error.is_unported() {
                return Err(error);
            }
        }
        Ok(output.point_class())
    }
    fn body(&self, role: InsertionBody) -> &Brep {
        match role {
            InsertionBody::Target => self.target,
            InsertionBody::Pseudo => self.pseudo.brep(),
            InsertionBody::Intermediate => self.intermediate,
        }
    }
    fn body_mut(&mut self, role: InsertionBody) -> &mut Brep {
        match role {
            InsertionBody::Target => self.target,
            InsertionBody::Pseudo => self.pseudo.brep_mut(),
            InsertionBody::Intermediate => self.intermediate,
        }
    }
    fn body_identity(&self, role: InsertionBody) -> Result<ObjectRef> {
        let serial = self
            .body(role)
            .notify_serial()
            .ok_or_else(|| invalid("insertion body has no issued identity"))?;
        Ok(topology_ref(serial, EntityRef::Brep(BrepId::from_index(0))))
    }
    fn topology(object: FilletObject) -> Option<(InsertionBody, EntityRef)> {
        match object {
            FilletObject::Topology(FilletBody::Target, e) => Some((InsertionBody::Target, e)),
            FilletObject::Topology(FilletBody::Pseudo, e) => Some((InsertionBody::Pseudo, e)),
            FilletObject::IntermediateTopology(e) => Some((InsertionBody::Intermediate, e)),
            _ => None,
        }
    }
    fn handle(role: InsertionBody, entity: EntityRef) -> FilletObject {
        match role {
            InsertionBody::Target => FilletObject::Topology(FilletBody::Target, entity),
            InsertionBody::Pseudo => FilletObject::Topology(FilletBody::Pseudo, entity),
            InsertionBody::Intermediate => FilletObject::IntermediateTopology(entity),
        }
    }
    fn role_of(&self, brep: &Brep) -> Option<InsertionBody> {
        [
            InsertionBody::Target,
            InsertionBody::Pseudo,
            InsertionBody::Intermediate,
        ]
        .into_iter()
        .find(|&role| std::ptr::eq(self.body(role), brep))
    }
    /// Capture the original Surface through its real Context identity slot.
    /// # Errors
    /// Missing actual Surface or unavailable identity issuance.
    pub fn capture_surface(
        &mut self,
        context: &mut Context,
        role: InsertionBody,
        face: FaceId,
    ) -> Result<SurfaceBinding> {
        let body = self.body_identity(role)?;
        let surface = self
            .body_mut(role)
            .face_mut(face)?
            .surface
            .as_mut()
            .ok_or_else(|| SolidError::null_pointer(OP))?;
        let type_code = surface.type_code();
        let key = context.geometry_key(surface.notify_identity_slot())?;
        Ok(SurfaceBinding {
            role,
            body,
            face,
            geometry: ObjectRef { type_code, key },
        })
    }
    /// Capture the actual managed fillet Surface without a numerical copy.
    /// # Errors
    /// Foreign managed body or unavailable owned-Surface binding.
    pub fn capture_fillet_surface(
        &mut self,
        context: &mut Context,
        geometry: FilletGeomId,
    ) -> Result<Option<SurfaceBinding>> {
        let Some(storage) = self.pseudo.geom(geometry)?.fillet_surface.as_ref() else {
            return Ok(None);
        };
        let (body, face) = match storage {
            FilletSurfaceStorage::ManagedFace { body, face } => (*body, *face),
            FilletSurfaceStorage::Owned(_) => {
                // UNPORTED(blend.fillet_executive_p1): actual owned fillet Surface binding.
                return Err(SolidError::unported(
                    "SmFilletGeom::GetFilletSurface (owned insertion binding)",
                    "blend.fillet_executive_p1",
                ));
            }
        };
        if body != self.body_identity(InsertionBody::Intermediate)? {
            return Err(invalid("managed fillet Surface belongs to a foreign body"));
        }
        self.capture_surface(context, InsertionBody::Intermediate, face)
            .map(Some)
    }
    /// Resolve that same Surface instance; a replacement in the Face slot fails.
    /// # Errors
    /// Foreign body, deleted Face/Surface or replaced geometry identity.
    pub fn bound_surface(&self, binding: SurfaceBinding) -> Result<&Surface> {
        if binding.body != self.body_identity(binding.role)? {
            return Err(invalid("Surface binding belongs to a foreign body"));
        }
        let surface = self
            .body(binding.role)
            .face(binding.face)?
            .surface
            .as_ref()
            .ok_or_else(|| SolidError::null_pointer(OP))?;
        if surface.notify_key() != Some(binding.geometry.key)
            || surface.type_code() != binding.geometry.type_code
        {
            return Err(invalid("original insertion Surface was replaced"));
        }
        Ok(surface)
    }
    /// Actual owner Face, including the real Surface owner check.
    /// # Errors
    /// Missing original Surface or unresolved foreign owner.
    pub fn surface_face(&self, binding: SurfaceBinding) -> Result<Option<FilletObject>> {
        let Some(owner) = self.bound_surface(binding)?.owner() else {
            return Ok(None);
        };
        if let Some((role, entity @ EntityRef::Face(_))) = self.owner_binding(owner) {
            return Ok(Some(Self::handle(role, entity)));
        }
        // UNPORTED(blend.fillet_executive_p1): actual foreign Surface owner lookup.
        Err(SolidError::unported(
            "SmSurface::GetFace (foreign insertion owner)",
            "blend.fillet_executive_p1",
        ))
    }
    fn owner_binding(&self, owner: ObjectRef) -> Option<(InsertionBody, EntityRef)> {
        let mut found = None;
        for role in [
            InsertionBody::Target,
            InsertionBody::Pseudo,
            InsertionBody::Intermediate,
        ] {
            let body = self.body(role);
            let serial = body.notify_serial()?;
            let canonical = body
                .faces()
                .ids()
                .map(EntityRef::Face)
                .chain(body.edges().ids().map(EntityRef::Edge))
                .chain(body.vertices().ids().map(EntityRef::Vertex));
            let mut matches = Vec::new();
            for entity in canonical {
                if topology_ref(serial, entity) == owner {
                    matches.push(entity);
                }
            }
            for (&entity, &producer) in body.geometry_owner_keys() {
                let live = match entity {
                    EntityRef::Face(id) => body.face(id).is_ok(),
                    EntityRef::Edge(id) => body.edge(id).is_ok(),
                    _ => false,
                };
                if live && producer == owner && !matches.contains(&entity) {
                    matches.push(entity);
                }
            }
            for entity in matches {
                if found.is_some_and(|actual| actual != (role, entity)) {
                    return None;
                }
                found = Some((role, entity));
            }
        }
        found
    }
}

impl ClassGeometry<FilletObject> for InsertionGeometry<'_> {
    fn curve(&self, object: FilletObject) -> Option<&Curve> {
        match object {
            FilletObject::InsertionCurve(binding) => self.curves.as_ref()?.get(binding),
            FilletObject::Curve(body, edge) => self
                .body(match body {
                    FilletBody::Target => InsertionBody::Target,
                    FilletBody::Pseudo => InsertionBody::Pseudo,
                })
                .edge(edge)
                .ok()?
                .curve
                .as_ref(),
            FilletObject::IntermediateCurve(edge) => {
                self.intermediate.edge(edge).ok()?.curve.as_ref()
            }
            FilletObject::UvCurve(body, eu) => self
                .body(match body {
                    FilletBody::Target => InsertionBody::Target,
                    FilletBody::Pseudo => InsertionBody::Pseudo,
                })
                .edgeuse(eu)
                .ok()?
                .uv_trim_curve
                .as_ref(),
            FilletObject::IntermediateUvCurve(eu) => {
                self.intermediate.edgeuse(eu).ok()?.uv_trim_curve.as_ref()
            }
            _ => None,
        }
    }
    fn surface(&self, object: FilletObject) -> Option<&Surface> {
        match object {
            FilletObject::Surface(binding) => self.bound_surface(binding).ok(),
            FilletObject::FaceSurface(role, face) => {
                self.body(role).face(face).ok()?.surface.as_ref()
            }
            _ => None,
        }
    }
    fn volume(&self, _: FilletObject) -> Option<&Volume> {
        None
    }
    fn point(&self, object: FilletObject) -> Option<Point3> {
        let (role, EntityRef::Vertex(vertex)) = Self::topology(object)? else {
            return None;
        };
        Some(self.body(role).vertex(vertex).ok()?.point)
    }
    fn classification_vertex_point(&self, object: FilletObject) -> Option<Point3> {
        self.point(object)
    }
    fn topology_zone_tol_3d(&self, object: FilletObject) -> Option<f64> {
        let (role, entity) = Self::topology(object)?;
        let body = self.body(role);
        match entity {
            EntityRef::Vertex(id) => Some(body.vertex(id).ok()?.tolerance),
            EntityRef::Edge(id) => Some(body.edge(id).ok()?.tolerance),
            EntityRef::Face(id) => Some(body.face(id).ok()?.tolerance),
            EntityRef::Edgeuse(id) => Some(body.edge(body.edgeuse(id).ok()?.edge?).ok()?.tolerance),
            _ => None,
        }
    }
    fn topology_entity(&self, object: FilletObject) -> Option<(&Brep, EntityRef)> {
        let (role, entity) = Self::topology(object)?;
        Some((self.body(role), entity))
    }
    fn topology_handle_in(&self, brep: &Brep, entity: EntityRef) -> Option<FilletObject> {
        Some(Self::handle(self.role_of(brep)?, entity))
    }
    fn owner_entity(&self, owner: ObjectRef) -> Option<(&Brep, EntityRef)> {
        let (role, entity) = self.owner_binding(owner)?;
        Some((self.body(role), entity))
    }
    fn brep_id(&self, brep: &Brep) -> Option<BrepId> {
        self.role_of(brep)?;
        Some(BrepId::from_index(
            u32::try_from(brep.notify_serial()?).ok()?,
        ))
    }
}
impl ClassGeometryMut<FilletObject> for InsertionGeometry<'_> {
    fn curve_mut(&mut self, object: FilletObject) -> Option<&mut Curve> {
        match object {
            FilletObject::InsertionCurve(binding) => self.curves.as_mut()?.get_mut(binding),
            FilletObject::Curve(body, edge) => self
                .body_mut(match body {
                    FilletBody::Target => InsertionBody::Target,
                    FilletBody::Pseudo => InsertionBody::Pseudo,
                })
                .edge_mut(edge)
                .ok()?
                .curve
                .as_mut(),
            FilletObject::IntermediateCurve(edge) => {
                self.intermediate.edge_mut(edge).ok()?.curve.as_mut()
            }
            FilletObject::UvCurve(body, eu) => self
                .body_mut(match body {
                    FilletBody::Target => InsertionBody::Target,
                    FilletBody::Pseudo => InsertionBody::Pseudo,
                })
                .edgeuse_mut(eu)
                .ok()?
                .uv_trim_curve
                .as_mut(),
            FilletObject::IntermediateUvCurve(eu) => self
                .intermediate
                .edgeuse_mut(eu)
                .ok()?
                .uv_trim_curve
                .as_mut(),
            _ => None,
        }
    }
    fn topology_entity_mut(&mut self, object: FilletObject) -> Option<(&mut Brep, EntityRef)> {
        let (role, entity) = Self::topology(object)?;
        Some((self.body_mut(role), entity))
    }
}

impl ProtoGeometry<FilletObject> for InsertionGeometry<'_> {
    fn curve_mut(&mut self, object: FilletObject) -> Option<&mut Curve> {
        ClassGeometryMut::curve_mut(self, object)
    }
}

impl ProtoTopology<FilletObject> for InsertionGeometry<'_> {
    fn vertex_use_faces(&self, vertex: FilletObject) -> Result<Vec<Option<FilletObject>>> {
        let (role, entity) = Self::topology(vertex).ok_or(invalid("unknown topology role"))?;
        Ok(ProtoTopology::vertex_use_faces(self.body(role), entity)?
            .into_iter()
            .map(|value| value.map(|entity| Self::handle(role, entity)))
            .collect())
    }
    fn edge_use_faces(&self, edge: FilletObject) -> Result<Vec<Option<FilletObject>>> {
        let (role, entity) = Self::topology(edge).ok_or(invalid("unknown topology role"))?;
        Ok(ProtoTopology::edge_use_faces(self.body(role), entity)?
            .into_iter()
            .map(|value| value.map(|entity| Self::handle(role, entity)))
            .collect())
    }
    fn edge_start_vertex(&self, edge: FilletObject) -> Result<Option<FilletObject>> {
        let (role, entity) = Self::topology(edge).ok_or(invalid("unknown topology role"))?;
        Ok(ProtoTopology::edge_start_vertex(self.body(role), entity)?
            .map(|entity| Self::handle(role, entity)))
    }
    fn edge_other_vertex(
        &self,
        edge: FilletObject,
        vertex: Option<FilletObject>,
    ) -> Result<Option<FilletObject>> {
        let (role, entity) = Self::topology(edge).ok_or(invalid("unknown topology role"))?;
        let vertex = vertex
            .map(|vertex| {
                let (other_role, vertex) =
                    Self::topology(vertex).ok_or(invalid("unknown vertex role"))?;
                if role != other_role {
                    return Err(invalid("edge and vertex belong to different bodies"));
                }
                Ok(vertex)
            })
            .transpose()?;
        Ok(
            ProtoTopology::edge_other_vertex(self.body(role), entity, vertex)?
                .map(|entity| Self::handle(role, entity)),
        )
    }
    fn vertex_point(&self, vertex: FilletObject) -> Result<Point3> {
        let (role, entity) = Self::topology(vertex).ok_or(invalid("unknown topology role"))?;
        ProtoTopology::vertex_point(self.body(role), entity)
    }
    fn edge_brep(&self, edge: FilletObject) -> Result<Option<ObjectKey>> {
        let (role, entity) = Self::topology(edge).ok_or(invalid("unknown topology role"))?;
        let actual_body = ProtoTopology::edge_brep(self.body(role), entity)?;
        actual_body
            .map(|_| self.body_identity(role).map(|identity| identity.key))
            .transpose()
    }
    fn edge_faces(&self, edge: FilletObject) -> Result<Vec<FilletObject>> {
        let (role, entity) = Self::topology(edge).ok_or(invalid("unknown topology role"))?;
        Ok(ProtoTopology::edge_faces(self.body(role), entity)?
            .into_iter()
            .map(|entity| Self::handle(role, entity))
            .collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitive::api_primitives::create_box;
    use crate::Kernel;

    fn box_at(k: &mut Kernel, x: f64) -> Brep {
        let mut b = create_box(k, Point3::new(x, 0.0, 0.0), 10.0, 20.0, 30.0).unwrap();
        crate::topology::brep::p3::notify_serial(&mut k.context, &mut b).unwrap();
        b
    }

    #[test]
    fn actual_equal_local_ids_keep_body_and_geometry_identity() {
        let mut k = Kernel::default();
        let mut target = box_at(&mut k, 0.0);
        let mut intermediate = box_at(&mut k, 100.0);
        let mut pseudo = FilletBrep::new();
        *pseudo.brep_mut() = box_at(&mut k, 200.0);
        let mut view = InsertionGeometry::new(&mut target, &mut pseudo, &mut intermediate).unwrap();
        let face = FaceId::from_index(0);
        let a = view
            .capture_surface(&mut k.context, InsertionBody::Target, face)
            .unwrap();
        let b = view
            .capture_surface(&mut k.context, InsertionBody::Intermediate, face)
            .unwrap();
        assert_ne!(a, b);
        assert_ne!(a.identity(), b.identity());
        assert!(std::ptr::eq(
            view.bound_surface(a).unwrap(),
            view.target.face(face).unwrap().surface.as_ref().unwrap()
        ));
        let edge = view
            .target
            .edges()
            .ids()
            .find(|&id| view.target.edge(id).unwrap().curve.is_some())
            .unwrap();
        assert!(!std::ptr::eq(
            view.curve(FilletObject::Curve(FilletBody::Target, edge))
                .unwrap(),
            view.curve(FilletObject::IntermediateCurve(edge)).unwrap()
        ));
        let vertex = view.target.vertices().ids().next().unwrap();
        let target_edge = FilletObject::Topology(FilletBody::Target, EntityRef::Edge(edge));
        let start = view.edge_start_vertex(target_edge).unwrap().unwrap();
        assert!(matches!(
            start,
            FilletObject::Topology(FilletBody::Target, EntityRef::Vertex(_))
        ));
        assert!(view
            .edge_faces(target_edge)
            .unwrap()
            .iter()
            .all(|face| matches!(
                face,
                FilletObject::Topology(FilletBody::Target, EntityRef::Face(_))
            )));
        let faces = view.edge_faces(target_edge).unwrap();
        let actual_face = faces.first().copied().unwrap();
        assert!(view
            .edge_is_connected_to_face(target_edge, actual_face)
            .unwrap());
        assert!(view
            .edge_edgeuse_of_face(target_edge, actual_face)
            .unwrap()
            .is_some());
        let FilletObject::Topology(_, EntityRef::Face(face_id)) = actual_face else {
            panic!("actual Face expected")
        };
        let foreign_face = FilletObject::IntermediateTopology(EntityRef::Face(face_id));
        assert!(!view
            .edge_is_connected_to_face(target_edge, foreign_face)
            .unwrap());
        assert!(view
            .edge_edgeuse_of_face(target_edge, foreign_face)
            .unwrap()
            .is_none());
        let current = view.face_surface(actual_face).unwrap().unwrap();
        assert!(std::ptr::eq(
            view.surface(current).unwrap(),
            view.target.face(face_id).unwrap().surface.as_ref().unwrap()
        ));
        let foreign = FilletObject::Topology(FilletBody::Pseudo, EntityRef::Vertex(vertex));
        assert!(view.edge_other_vertex(target_edge, Some(foreign)).is_err());
        let before = view.target.vertex(vertex).unwrap().point;
        let object = FilletObject::Topology(FilletBody::Pseudo, EntityRef::Vertex(vertex));
        let (body, _) = view.topology_entity_mut(object).unwrap();
        body.vertex_mut(vertex).unwrap().point.x += 7.0;
        assert_eq!(view.target.vertex(vertex).unwrap().point, before);
        assert_ne!(view.pseudo.brep().vertex(vertex).unwrap().point, before);
    }

    #[test]
    fn actual_surface_replacement_cannot_reuse_captured_binding() {
        let mut k = Kernel::default();
        let mut target = box_at(&mut k, 0.0);
        let mut intermediate = box_at(&mut k, 100.0);
        let mut replacement = box_at(&mut k, 300.0);
        let mut pseudo = FilletBrep::new();
        *pseudo.brep_mut() = box_at(&mut k, 200.0);
        let mut view = InsertionGeometry::new(&mut target, &mut pseudo, &mut intermediate).unwrap();
        let face = FaceId::from_index(0);
        let binding = view
            .capture_surface(&mut k.context, InsertionBody::Intermediate, face)
            .unwrap();
        let actual = replacement.face_mut(face).unwrap().surface.take();
        view.intermediate.face_mut(face).unwrap().surface = actual;
        assert!(view.bound_surface(binding).is_err());
        let current = view
            .face_surface(FilletObject::IntermediateTopology(EntityRef::Face(face)))
            .unwrap()
            .unwrap();
        assert!(std::ptr::eq(
            view.surface(current).unwrap(),
            view.intermediate
                .face(face)
                .unwrap()
                .surface
                .as_ref()
                .unwrap()
        ));
    }
    #[test]
    fn retained_curves_and_lent_notifications_survive_partial_edit_failure() {
        use crate::geom::context::NotifyCallback;
        use crate::geom::object::NotifyOperation;
        use std::sync::{Arc, Mutex};
        #[derive(Debug)]
        struct Capture(Arc<Mutex<Vec<(NotifyOperation, ObjectRef)>>>);
        impl NotifyCallback for Capture {
            fn execute(
                &mut self,
                op: NotifyOperation,
                caller: ObjectRef,
                _one: Option<ObjectRef>,
                _two: Option<ObjectRef>,
                _three: Option<ObjectRef>,
            ) -> Result<()> {
                self.0.lock().unwrap().push((op, caller));
                Ok(())
            }
        }
        let mut k = Kernel::default();
        let mut target = box_at(&mut k, 0.0);
        let mut intermediate = box_at(&mut k, 100.0);
        let mut source = box_at(&mut k, 300.0);
        let edge = source
            .edges()
            .ids()
            .find(|&id| source.edge(id).unwrap().curve.is_some())
            .unwrap();
        let curve = source.edge_mut(edge).unwrap().curve.take().unwrap();
        let owner = curve.owner();
        let mut pseudo = FilletBrep::new();
        *pseudo.brep_mut() = box_at(&mut k, 200.0);
        let vertex = target.vertices().ids().next().unwrap();
        let serial = target.notify_serial().unwrap();
        let caller = topology_ref(serial, EntityRef::Vertex(vertex));
        let changed = Point3::new(17.0, 18.0, 19.0);
        let mut curves = InsertionCurves::default();
        let mut records = BTreeMap::new();
        let mut attributes = AttributeStore::default();
        let events = Arc::new(Mutex::new(Vec::new()));
        k.context
            .set_user_notify_callback(Some(Box::new(Capture(events.clone()))));
        let binding = crate::topology::brep::p3::with_lent_context(
            &mut k,
            |_, context| -> Result<FilletObject> {
                let mut view = InsertionGeometry::new(&mut target, &mut pseudo, &mut intermediate)?
                    .with_notifications(
                        &mut curves,
                        InsertionNotifications {
                            records: &mut records,
                            attributes: &mut attributes,
                            context,
                        },
                    );
                let binding = view.insert_curve_copy(curve)?;
                assert_eq!(view.curve(binding).unwrap().owner(), owner);
                assert!(ClassGeometryMut::curve_mut(&mut view, binding).is_some());
                let result: Result<()> = view.with_edit(InsertionBody::Target, |body, env| {
                    let vu_ref = |id| topology_ref(serial, EntityRef::Vertexuse(id));
                    crate::topology::vertex::set_point(
                        body,
                        vertex,
                        changed,
                        env.records,
                        env.attributes,
                        env.context.as_deref_mut(),
                        caller,
                        &vu_ref,
                    )?;
                    Err(invalid("deliberate failure after actual vertex edit"))
                });
                assert!(result.is_err());
                assert_eq!(view.target.vertex(vertex)?.point, changed);
                Ok(binding)
            },
        )
        .unwrap()
        .unwrap();
        let FilletObject::InsertionCurve(binding) = binding else {
            panic!("actual curve binding required")
        };
        assert_eq!(curves.checked(binding).unwrap().owner(), owner);
        assert!(InsertionCurves::default().get(binding).is_none());
        assert!(records.contains_key(&caller.key));
        assert!(events
            .lock()
            .unwrap()
            .iter()
            .any(|&(op, object)| op == NotifyOperation::ChangeGeometry && object == caller));
    }
}
