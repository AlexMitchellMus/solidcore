//! Actual intermediate body: geometry plus directed topology in owning-list order.
use super::super::{label, topology};
use super::{FilletSession, SolidResult};
use solidcore::kernel::Kernel;
use solidcore::topology::{brep::p1, edge, edgeuse::p1 as eu, face::p2 as face};

pub(super) fn snapshot(
    k: &mut Kernel,
    session: &FilletSession<'_>,
) -> SolidResult<(Vec<f64>, Vec<f64>)> {
    let Some(body) = session.intermediate_brep() else {
        return Ok((vec![0.0], Vec::new()));
    };
    let vertices = p1::vertices(body, None)?;
    let edges = p1::edges(body, None)?;
    let faces = p1::faces(body, None)?;
    let regions = p1::regions(body, None)?;
    let shells = p1::shells(k, body, None)?;
    let mut uses = Vec::new();
    for &id in &edges {
        uses.extend(edge::edgeuses(body, id, None, None)?);
    }
    let mut state = vec![
        1.0,
        vertices.len() as f64,
        edges.len() as f64,
        faces.len() as f64,
        regions.len() as f64,
        shells.len() as f64,
        f64::from(body.editing_enabled()),
    ];
    topology(body, &mut state)?;
    state.push(uses.len() as f64);
    for &id in &uses {
        let record = body.edgeuse(id)?;
        let vertex = record
            .vertexuse
            .map(|id| body.vertexuse(id).map(|vu| vu.vertex))
            .transpose()?
            .flatten();
        state.extend([
            label(&edges, record.edge),
            f64::from(record.orientation as i32),
            label(&vertices, vertex),
            label(&uses, eu::mate(body, id)?),
            label(&uses, eu::radial(body, id)?),
            label(&uses, eu::cw_edgeuse(body, id)?),
            label(&uses, eu::ccw_edgeuse(body, id)?),
            label(&faces, eu::face(body, id)?),
        ]);
    }
    for id in faces {
        let (a, b) = face::faceuses(body, id)?;
        for id in [a, b] {
            state.push(f64::from(id.is_some()));
            if let Some(id) = id {
                let use_ = body.faceuse(id)?;
                state.extend([
                    f64::from(use_.orientation as i32),
                    label(&shells, use_.shell),
                    label(&regions, eu::faceuse_region(body, id)?),
                ]);
            }
        }
    }
    for id in shells {
        state.push(label(&regions, body.shell(id)?.region));
    }
    for id in regions {
        state.push(f64::from(body.region(id)?.is_void));
    }
    let mut geometry = vec![body.tolerance()];
    crate::util::primitive::append(k, body, &mut geometry)?;
    Ok((state, geometry))
}
