//! Body-qualified identities retained by fillet rail classifications.

use crate::curve::curve::Curve;
use crate::error::{EntityRef, Result, SolidError};
use crate::geom::object::ObjectRef;
use crate::intersect::point_class::{ClassGeometry, ClassObject, ClassObjectKind};
use crate::math::Point3;
use crate::query::face::p2::CurveOnClassifyObjects;
use crate::surface::surface::Surface;
use crate::surface::volume::Volume;
use crate::topology::{Brep, BrepId, EdgeId, EdgeuseId, FaceId, LoopId, VertexId, VertexuseId};

/// Which actual body owns fillet geometry or topology.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum FilletBody {
    /// The caller's target B-rep.
    Target,
    /// The executive's pseudo-B-rep.
    Pseudo,
}

/// An actual classification object, or an explicitly indexed legacy handle.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum FilletClassObject {
    /// A handle supplied by an existing indexed classification caller.
    Legacy(EntityRef),
    /// Topology in the specified actual body.
    Topology(FilletBody, EntityRef),
    /// The retained 3D curve of this body's edge.
    Curve(FilletBody, EdgeId),
    /// The retained UV curve of this body's edgeuse.
    TrimCurve(FilletBody, EdgeuseId),
    /// The retained surface of this body's face.
    Surface(FilletBody, FaceId),
}

impl FilletClassObject {
    /// Resolve a target topology object without accepting pseudo or geometry handles.
    ///
    /// # Errors
    /// A handle outside the target topology domain.
    pub fn target_entity(self) -> Result<EntityRef> {
        match self {
            Self::Legacy(entity) | Self::Topology(FilletBody::Target, entity) => Ok(entity),
            _ => Err(SolidError::invalid_input(
                "blend::fillet_classification",
                "classification object is not target topology",
            )),
        }
    }

    /// Project only an explicitly legacy handle for an indexed provider.
    #[must_use]
    pub const fn legacy(self) -> Option<EntityRef> {
        match self {
            Self::Legacy(entity) => Some(entity),
            _ => None,
        }
    }
}

impl ClassObject for FilletClassObject {
    fn class_kind(&self) -> ClassObjectKind {
        match self {
            Self::Legacy(entity) | Self::Topology(_, entity) => entity.class_kind(),
            Self::Curve(..) | Self::TrimCurve(..) => ClassObjectKind::Curve,
            Self::Surface(..) => ClassObjectKind::Surface,
        }
    }
}

/// Checked namespace adapter for existing indexed geometry providers.
pub struct LegacyClassificationObjects<'a, T: ?Sized>(pub &'a T);

impl<T: ClassGeometry<EntityRef> + ?Sized> ClassGeometry<FilletClassObject>
    for LegacyClassificationObjects<'_, T>
{
    fn curve(&self, o: FilletClassObject) -> Option<&Curve> {
        self.0.curve(o.legacy()?)
    }
    fn surface(&self, o: FilletClassObject) -> Option<&Surface> {
        self.0.surface(o.legacy()?)
    }
    fn volume(&self, o: FilletClassObject) -> Option<&Volume> {
        self.0.volume(o.legacy()?)
    }
    fn point(&self, o: FilletClassObject) -> Option<Point3> {
        self.0.point(o.legacy()?)
    }
    fn classification_vertex_point(&self, o: FilletClassObject) -> Option<Point3> {
        self.0.classification_vertex_point(o.legacy()?)
    }
    fn topology_zone_tol_3d(&self, o: FilletClassObject) -> Option<f64> {
        self.0.topology_zone_tol_3d(o.legacy()?)
    }
    fn vertex_face_gap(
        &self,
        vertex: FilletClassObject,
        surface: &Surface,
        face: Option<FilletClassObject>,
    ) -> Option<f64> {
        let face = match face {
            Some(face) => Some(face.legacy()?),
            None => None,
        };
        self.0.vertex_face_gap(vertex.legacy()?, surface, face)
    }
    fn topology_entity(&self, o: FilletClassObject) -> Option<(&Brep, EntityRef)> {
        self.0.topology_entity(o.legacy()?)
    }
    fn topology_handle(&self, entity: EntityRef) -> Option<FilletClassObject> {
        self.0
            .topology_handle(entity)
            .map(FilletClassObject::Legacy)
    }
    fn topology_handle_in(&self, brep: &Brep, entity: EntityRef) -> Option<FilletClassObject> {
        self.0
            .topology_handle_in(brep, entity)
            .map(FilletClassObject::Legacy)
    }
    fn topology_intersector_mate(
        &self,
        ti: FilletClassObject,
        entity: FilletClassObject,
    ) -> Option<Option<FilletClassObject>> {
        self.0
            .topology_intersector_mate(ti.legacy()?, entity.legacy()?)
            .map(|v| v.map(FilletClassObject::Legacy))
    }
    fn brep_id(&self, brep: &Brep) -> Option<BrepId> {
        self.0.brep_id(brep)
    }
    fn owner_entity(&self, owner: ObjectRef) -> Option<(&Brep, EntityRef)> {
        self.0.owner_entity(owner)
    }
    fn owner_of(&self, owner: ObjectRef) -> Option<ObjectRef> {
        self.0.owner_of(owner)
    }
}

impl<T: CurveOnClassifyObjects<EntityRef> + ?Sized> CurveOnClassifyObjects<FilletClassObject>
    for LegacyClassificationObjects<'_, T>
{
    fn face(&self, f: FaceId) -> FilletClassObject {
        FilletClassObject::Legacy(self.0.face(f))
    }
    fn vertex(&self, v: VertexId) -> FilletClassObject {
        FilletClassObject::Legacy(self.0.vertex(v))
    }
    fn edge(&self, e: EdgeId) -> FilletClassObject {
        FilletClassObject::Legacy(self.0.edge(e))
    }
    fn vertexuse(&self, vu: VertexuseId) -> FilletClassObject {
        FilletClassObject::Legacy(self.0.vertexuse(vu))
    }
    fn edgeuse(&self, eu: EdgeuseId) -> FilletClassObject {
        FilletClassObject::Legacy(self.0.edgeuse(eu))
    }
    fn edge_curve(&self, e: EdgeId) -> FilletClassObject {
        FilletClassObject::Legacy(self.0.edge_curve(e))
    }
    fn trim_curve(&self, eu: EdgeuseId) -> FilletClassObject {
        FilletClassObject::Legacy(self.0.trim_curve(eu))
    }
    fn face_surface(&self, f: FaceId) -> Option<FilletClassObject> {
        self.0.face_surface(f).map(FilletClassObject::Legacy)
    }
    fn loop_id(&self, o: FilletClassObject) -> Option<LoopId> {
        self.0.loop_id(o.legacy()?)
    }
}
