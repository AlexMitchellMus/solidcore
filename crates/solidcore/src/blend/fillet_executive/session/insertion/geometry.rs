//! Actual target/intermediate topology and pseudo-curve bindings for insertion.

use super::super::{FilletBody, FilletObject};
use crate::blend::fillet_geom::p1::{FilletBrep, FilletGeomId, FilletSurfaceStorage};
use crate::curve::curve::Curve;
use crate::geom::context::Context;
use crate::geom::object::ObjectRef;
use crate::intersect::curve_class::p5::ClassGeometryMut;
use crate::intersect::point_class::ClassGeometry;
use crate::math::Point3;
use crate::surface::surface::Surface;
use crate::surface::volume::Volume;
use crate::topology::brep::p3::topology_ref;
use crate::topology::{Brep, BrepId, FaceId};
use crate::{EntityRef, Result, SolidError};

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
        let FilletObject::Surface(binding) = object else {
            return None;
        };
        self.bound_surface(binding).ok()
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
    }
}
