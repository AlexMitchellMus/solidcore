//! Direct native helper proof; not a full ProcessRollOver fixture.
use super::{ffi, harness};
use crate::{HarnessError, Output};
use solidcore::blend::fillet_geom::p1::{FilletBrep, FilletGeom, FilletGeomType, FilletVertex};
use solidcore::blend::fillet_solver::p1::{split_rail, RailOwnership, RetainedRailOwnership};
use solidcore::curve::{curve::Curve, line::Line};
use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::geom::context::NotifyCallback;
use solidcore::geom::object::{NotifyOperation, ObjectRef};
use solidcore::math::Vec3;
use solidcore::topology::{brep::p3, edge};
use solidcore::EntityRef;
use solidcore::Kernel;
use std::collections::BTreeMap;
use std::sync::{Arc, Mutex};

#[derive(Debug)]
struct Event {
    system: bool,
    caller: ObjectRef,
    data: [Option<ObjectRef>; 3],
    owner: Option<ObjectRef>,
}
#[derive(Debug)]
struct Recorder {
    system: bool,
    events: Arc<Mutex<Vec<Event>>>,
    owners: BTreeMap<ObjectKey, Option<ObjectRef>>,
}
impl NotifyCallback for Recorder {
    fn execute(
        &mut self,
        operation: NotifyOperation,
        caller: ObjectRef,
        a: Option<ObjectRef>,
        b: Option<ObjectRef>,
        c: Option<ObjectRef>,
    ) -> solidcore::Result<()> {
        if operation == NotifyOperation::ChangeOwner {
            self.events.lock().unwrap().push(Event {
                system: self.system,
                caller,
                data: [a, b, c],
                owner: self.owners.get(&caller.key).copied().flatten(),
            });
        }
        Ok(())
    }
    fn owner_link(&mut self, object: ObjectKey, owner: Option<ObjectRef>) {
        self.owners.insert(object, owner);
    }
}
fn flag(value: bool) -> f64 {
    if value {
        1.
    } else {
        0.
    }
}
fn output(flat: &[f64]) -> Result<Output, HarnessError> {
    if flat.len() < 2 || flat[0] != 0.0 {
        return Err(harness("direct split did not complete"));
    }
    let count = flat[1] as usize;
    let boundary = 2usize
        .checked_add(
            count
                .checked_mul(7)
                .ok_or_else(|| harness("event overflow"))?,
        )
        .and_then(|n| n.checked_add(12))
        .ok_or_else(|| harness("record overflow"))?;
    if flat.len() != boundary + 15 {
        return Err(harness("invalid split record length"));
    }
    let mut exact = Vec::new();
    for &value in &flat[..boundary] {
        if !value.is_finite() || value.fract() != 0.0 {
            return Err(harness("noninteger split field"));
        }
        exact.push(value as i64);
    }
    let mut values = flat[boundary..boundary + 3].to_vec();
    for chunk in flat[boundary + 3..].chunks_exact(4) {
        if chunk[0] != 0.0 {
            return Err(harness("split curve sample failed"));
        }
        exact.push(chunk[0] as i64);
        values.extend_from_slice(&chunk[1..]);
    }
    Ok(Output::new()
        .ints("actual_identity_topology", exact)
        .values("actual_curve", values))
}
pub(super) fn cpp() -> Result<Output, HarnessError> {
    output(
        &ffi::eval_static(ffi::StaticOp::SplitRail, "", &[]).map_err(|e| harness(e.to_string()))?,
    )
}
pub(super) fn rust() -> Result<Output, HarnessError> {
    let run = || -> solidcore::Result<Vec<f64>> {
        let mut k = Kernel::default();
        let mut pseudo = FilletBrep::new();
        let mut vertices = Vec::new();
        for point in [Vec3::ZERO, Vec3::X] {
            let id = pseudo.brep_mut().add_vertex(point)?;
            pseudo.insert_fillet_vertex(id, FilletVertex::new(&k, None));
            vertices.push(id);
        }
        let rail = pseudo.make_fillet_edge(Some(vertices[0]), Some(vertices[1]), None)?;
        let line = Curve::from(Line::from_points(Vec3::ZERO, Vec3::X, 3));
        pseudo.brep_mut().edge_mut(rail)?.curve = Some(line);
        pseudo.brep_mut().edge_mut(rail)?.tolerance = 0.125;
        let mut attributes = AttributeStore::default();
        RetainedRailOwnership::with_context(
            &mut k,
            &mut pseudo,
            &mut attributes,
            |_, pseudo, ownership| {
                let (owner, brep) = ownership.edge_owner(rail)?;
                let curve = pseudo.brep_mut().edge_mut(rail)?.curve.as_mut().unwrap();
                ownership.set_curve_owner(rail, curve, owner, brep)
            },
        )?;
        let serial = pseudo.brep().notify_serial().unwrap();
        let old_key = pseudo
            .brep()
            .edge(rail)?
            .curve
            .as_ref()
            .unwrap()
            .notify_key()
            .unwrap();
        let geom = pseudo.add_geom(FilletGeom::new(&k, None, FilletGeomType::Default)?);
        let events = Arc::new(Mutex::new(Vec::new()));
        // This case owns a fresh Context, with no existing callback to displace.
        k.context.set_user_notify_callback(Some(Box::new(Recorder {
            system: false,
            events: events.clone(),
            owners: BTreeMap::new(),
        })));
        k.context.set_sys_notify_callback(Some(Box::new(Recorder {
            system: true,
            events: events.clone(),
            owners: BTreeMap::new(),
        })));
        let (middle, new_rail) = RetainedRailOwnership::with_context(
            &mut k,
            &mut pseudo,
            &mut attributes,
            |k, pseudo, ownership| split_rail(k, pseudo, ownership, geom, rail),
        )?;
        let copied = pseudo.brep().edge(new_rail)?.curve.as_ref().unwrap();
        let owner = p3::topology_ref(serial, EntityRef::Edge(new_rail));
        let actual = ObjectRef {
            key: copied.notify_key().unwrap(),
            type_code: copied.type_code(),
        };
        let events = events.lock().unwrap();
        let mut out = vec![0., events.len() as f64];
        for e in events.iter() {
            out.extend([
                flag(e.system),
                e.caller.type_code as f64,
                flag(e.caller == actual),
                flag(e.data[0] == Some(owner)),
                flag(e.data[1].is_none()),
                flag(e.data[2].is_none()),
                flag(e.owner.is_none()),
            ]);
        }
        out.extend([
            copied.type_code() as f64,
            flag(actual.key != old_key),
            flag(copied.owner() == Some(owner)),
            flag(edge::vertex(pseudo.brep(), rail)? == Some(vertices[0])),
            flag(edge::other_vertex(pseudo.brep(), rail, Some(vertices[0]))? == Some(middle)),
            flag(edge::vertex(pseudo.brep(), new_rail)? == Some(middle)),
            flag(edge::other_vertex(pseudo.brep(), new_rail, Some(middle))? == Some(vertices[1])),
            pseudo.brep().edge(rail)?.edgeuses.len() as f64,
            pseudo.brep().edge(new_rail)?.edgeuses.len() as f64,
            pseudo.brep().vertex(vertices[0])?.vertexuses.len() as f64,
            pseudo.brep().vertex(middle)?.vertexuses.len() as f64,
            pseudo.brep().vertex(vertices[1])?.vertexuses.len() as f64,
        ]);
        let domain = copied.natural_interval();
        out.extend([
            pseudo.brep().edge(new_rail)?.tolerance,
            domain.min,
            domain.max,
        ]);
        for f in [0., 0.5, 1.] {
            let p = copied.evaluate_point(&k, domain.evaluate(f))?;
            out.extend([0., p.x, p.y, p.z]);
        }
        Ok(out)
    };
    output(&run().map_err(|e| harness(e.to_string()))?)
}
