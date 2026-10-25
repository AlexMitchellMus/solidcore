//! Existing notified topology edits on the provider's actual bodies.

use super::{invalid, EntityRef, FilletObject, InsertionBody, InsertionGeometry, Result};
use crate::topology::proto_topology::p3::BrepGlue;
use crate::topology::{EdgeId, VertexId};
use crate::Kernel;

fn body_role(object: FilletObject) -> Result<InsertionBody> {
    let (role, EntityRef::Brep(_)) =
        InsertionGeometry::topology(object).ok_or_else(|| invalid("unknown body role"))?
    else {
        return Err(invalid("handle is not a B-rep"));
    };
    Ok(role)
}
fn vertex(role: InsertionBody, object: FilletObject) -> Result<VertexId> {
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
            crate::topology::brep::p5::glue_with_notifications(
                k,
                body,
                Some(first),
                Some(second),
                wires.as_deref(),
                Some(env),
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
