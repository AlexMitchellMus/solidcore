//! One real selection call, its ordered outputs, and passive retained mark records.
use super::{invalid, label, section, status, FilletSession, SolidResult};
use solidcore::blend::fillet_executive::session::{
    FilletBody, FilletObject, OriginalSurfaceSelection,
};
use solidcore::geom::context::MarkType;
use solidcore::kernel::Kernel;
use solidcore::topology::brep::{p1, p3::topology_ref};
use solidcore::EntityRef;

fn context(k: &Kernel, out: &mut Vec<f64>) {
    for mark in [
        MarkType::Mark,
        MarkType::Mark2,
        MarkType::Mark3,
        MarkType::MarkIo,
        MarkType::MarkAssert,
    ] {
        out.extend([
            f64::from(mark as i32),
            k.context.current_mark(mark) as f64,
            f64::from(k.context.is_mark_locked(mark)),
        ]);
    }
}

pub(super) fn run(
    k: &mut Kernel,
    session: &mut FilletSession<'_>,
) -> SolidResult<(Vec<f64>, bool)> {
    let mut out = vec![f64::from(session.executive().options.global_merge)];
    context(k, &mut out);
    let target_vertices = p1::vertices(session.target(), None)?;
    let target_edges = p1::edges(session.target(), None)?;
    let managed = session
        .intermediate_brep()
        .ok_or_else(|| invalid("selection requires actual managed body"))?;
    let managed_vertices = p1::vertices(managed, None)?;
    let managed_edges = p1::edges(managed, None)?;
    for entities in [
        target_vertices
            .iter()
            .copied()
            .map(EntityRef::Vertex)
            .collect::<Vec<_>>(),
        target_edges.iter().copied().map(EntityRef::Edge).collect(),
    ] {
        let mut mates = vec![entities.len() as f64];
        for (index, entity) in entities.into_iter().enumerate() {
            let mate = session
                .executive()
                .intersector
                .other_mate(Some(FilletObject::Topology(FilletBody::Target, entity)));
            let other = match (entity, mate) {
                (_, None) => -1.0,
                (
                    EntityRef::Vertex(_),
                    Some(FilletObject::IntermediateTopology(EntityRef::Vertex(id))),
                ) => label(&managed_vertices, Some(id)),
                (
                    EntityRef::Edge(_),
                    Some(FilletObject::IntermediateTopology(EntityRef::Edge(id))),
                ) => label(&managed_edges, Some(id)),
                _ => return Err(invalid("selection mate has unexpected actual body role")),
            };
            if mate.is_some() && other < 0.0 {
                return Err(invalid("selection mate is not in actual owning traversal"));
            }
            mates.extend([index as f64, f64::from(mate.is_some()), other]);
        }
        section(&mut out, &mates);
    }
    let mut selected = OriginalSurfaceSelection::default();
    let result = session.select_original_fillet_surfaces(k, &mut selected);
    let success = result.is_ok();
    out.push(match result {
        Ok(()) => 0.0,
        Err(error) => status(error)?,
    });
    context(k, &mut out);
    out.extend([
        f64::from(k.context.doing_boolean()),
        f64::from(session.target().editing_enabled()),
        f64::from(session.executive().intersector.doing_fillets),
    ]);
    let target_faces = p1::faces(session.target(), None)?;
    let managed_faces = p1::faces(
        session
            .intermediate_brep()
            .ok_or_else(|| invalid("selection lost actual managed body"))?,
        None,
    )?;
    for (selected, owning) in [
        (&selected.kept, &managed_faces),
        (&selected.deleted, &target_faces),
    ] {
        let mut labels = Vec::new();
        for id in selected {
            let id = label(owning, Some(*id));
            if id < 0.0 {
                return Err(invalid("selection output is outside actual owning body"));
            }
            labels.push(id);
        }
        section(&mut out, &labels);
    }
    let state = session.insertion_state()?;
    let managed_serial = state
        .intermediate
        .notify_serial()
        .ok_or_else(|| invalid("managed selection body lacks issued identity"))?;
    for (body, serial) in [
        (&*state.target, state.target_serial),
        (&*state.intermediate, managed_serial),
    ] {
        for entities in [
            p1::vertices(body, None)?
                .into_iter()
                .map(EntityRef::Vertex)
                .collect::<Vec<_>>(),
            p1::edges(body, None)?
                .into_iter()
                .map(EntityRef::Edge)
                .collect(),
            p1::faces(body, None)?
                .into_iter()
                .map(EntityRef::Face)
                .collect(),
        ] {
            let mut marks = vec![entities.len() as f64];
            for entity in entities {
                let topology = state
                    .notification_records
                    .get(&topology_ref(serial, entity).key)
                    .ok_or_else(|| invalid("selection actual topology record missing"))?;
                marks.extend([
                    topology.get_mark(MarkType::Mark) as f64,
                    f64::from(topology.is_marked(&k.context, MarkType::Mark)),
                ]);
            }
            section(&mut out, &marks);
        }
    }
    Ok((out, success))
}
