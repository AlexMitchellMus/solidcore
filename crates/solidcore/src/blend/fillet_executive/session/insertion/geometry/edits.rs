//! Existing notified topology edits on the provider's actual bodies.

use super::{invalid, EntityRef, FilletObject, InsertionBody, InsertionGeometry, Result};
use crate::topology::proto_topology::p3::BrepGlue;
use crate::topology::{EdgeId, VertexId};
use crate::Kernel;

pub(super) fn body_role(object: FilletObject) -> Result<InsertionBody> {
    let (role, EntityRef::Brep(_)) =
        InsertionGeometry::topology(object).ok_or_else(|| invalid("unknown body role"))?
    else {
        return Err(invalid("handle is not a B-rep"));
    };
    Ok(role)
}
pub(super) fn vertex(role: InsertionBody, object: FilletObject) -> Result<VertexId> {
    let (actual, EntityRef::Vertex(id)) =
        InsertionGeometry::topology(object).ok_or_else(|| invalid("unknown Vertex role"))?
    else {
        return Err(invalid("handle is not a Vertex"));
    };
    if actual != role {
        return Err(invalid("Vertex is not in the actual B-rep"));
    }
    Ok(id)
}
fn edge(role: InsertionBody, object: FilletObject) -> Result<EdgeId> {
    let (actual, EntityRef::Edge(id)) =
        InsertionGeometry::topology(object).ok_or_else(|| invalid("unknown Edge role"))?
    else {
        return Err(invalid("handle is not an Edge"));
    };
    if actual != role {
        return Err(invalid("Edge is not in the actual B-rep"));
    }
    Ok(id)
}

impl InsertionGeometry<'_> {
    /// Squeeze an actual edge through the existing notified producer.
    /// # Errors
    /// Cross-body operands or the original squeeze failure.
    pub fn brep_squeeze_edge(
        &mut self,
        k: &mut Kernel,
        body: FilletObject,
        object: FilletObject,
        survivor: FilletObject,
    ) -> Result<()> {
        let role = body_role(body)?;
        let object = edge(role, object)?;
        let survivor = vertex(role, survivor)?;
        self.with_edit(role, |body, env| {
            crate::topology::brep::p3::squeeze_edge(k, body, env, object, survivor)
        })
    }
}

impl BrepGlue<FilletObject> for InsertionGeometry<'_> {
    fn glue_vertices(
        &mut self,
        k: &mut Kernel,
        body: FilletObject,
        first: FilletObject,
        second: FilletObject,
        wires: Option<&[FilletObject]>,
    ) -> Result<()> {
        let role = body_role(body)?;
        let first = vertex(role, first)?;
        let second = vertex(role, second)?;
        let wires = wires
            .map(|items| {
                items
                    .iter()
                    .map(|&object| edge(role, object))
                    .collect::<Result<Vec<_>>>()
            })
            .transpose()?;
        self.with_edit(role, |body, env| {
            let identity = (env.key_of)(EntityRef::Brep(crate::topology::BrepId::from_index(0)));
            crate::topology::brep::p5::glue_with_body_identity(
                k,
                body,
                Some(first),
                Some(second),
                wires.as_deref(),
                env,
                identity,
            )
        })
    }
    fn delete_edge(
        &mut self,
        k: &mut Kernel,
        body: FilletObject,
        object: FilletObject,
    ) -> Result<()> {
        let role = body_role(body)?;
        let id = edge(role, object)?;
        self.with_edit(role, |body, env| {
            let mut delete = crate::topology::brep::p4::DeleteEnv::from_notify(env);
            crate::topology::brep::p4::delete_edge(k, body, &mut delete, id, None)
        })
    }
}

impl crate::intersect::curve_class::p6::MergeIntervalBrep<FilletObject> for InsertionGeometry<'_> {
    fn merge_context(
        &mut self,
        _object: FilletObject,
    ) -> Option<crate::intersect::curve_class::p6::MergeContext<'_, FilletObject>> {
        // The legacy infallible callback cannot represent an actual new Curve identity.
        None
    }
    fn merge_context_with_curve_identity(
        &mut self,
        object: FilletObject,
    ) -> Option<crate::intersect::curve_class::p6::MergeContextWithIdentity<'_, FilletObject>> {
        let (role, target) = Self::topology(object)?;
        if !matches!(
            target,
            EntityRef::Vertex(_) | EntityRef::Edge(_) | EntityRef::Face(_) | EntityRef::Region(_)
        ) {
            return None;
        }
        let key_of = self.topology_keys.as_ref()?.key_of(role);
        let body = match role {
            InsertionBody::Target => &mut *self.target,
            InsertionBody::Pseudo => self.pseudo.brep_mut(),
            InsertionBody::Intermediate => &mut *self.intermediate,
        };
        let notifications = self.notifications.as_mut()?;
        crate::topology::brep::p3::seed_records(notifications.records, body, key_of);
        Some(
            crate::intersect::curve_class::p6::MergeContextWithIdentity {
                brep: body,
                target,
                handle: Box::new(move |entity| Self::handle(role, entity)),
                notifications: crate::topology::edge::NotifyEnv {
                    records: notifications.records,
                    attributes: notifications.attributes,
                    context: Some(notifications.context),
                    key_of,
                },
                curve_identity: crate::topology::brep::p9::CurveIdentity::Actual(
                    &crate::topology::brep::p9::issue_curve_identity,
                ),
            },
        )
    }
}
