//! Differential runner for unit `topology.brep_p4`: `SmBrep::DeleteEdge`,
//! `DeleteVertex`, `ShrinkGeometry`, `StitchFaces` and `~SmBrep` on real B-reps.
//!
//! The B-rep comes from the shared fixture (`sc_topology_fixture.h`): the C++ side
//! builds it, the Rust side imports the `topology.topology_traverser` export of the same
//! block record by record (with the vertex use lists), so the target indices name the
//! same objects on both sides. Compared: the status, and on success the ten record
//! counts of the B-rep afterwards (`StitchFaces` adds its outputs; the destructor its
//! context flag). Cases whose Rust path reaches an unported callee are deferred.
use super::topology_vertexuse::{self as shared, BrepIn};
use crate::{Case, HarnessError, Output, Rng, UnitDiff};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use solidcore::error::{EntityRef, SolidError};
use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::geom::object::ObjectRef;
use solidcore::kernel::Kernel;
use solidcore::topology::brep::p4::{self, DeleteEnv};
use solidcore::topology::topology::Topology;
use solidcore::topology::{Brep, EdgeId, FaceId, VertexId, VertexuseId};
use solidcore_ffi::units::topology_brep_p4::{self as ffi, Op, Request};
use solidcore_ffi::units::topology_vertexuse::counts;
use std::collections::BTreeMap;

/// Differential runner for `topology.brep_p4`.
pub struct Runner;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
enum OpIn {
    DeleteEdge,
    DeleteVertex,
    ShrinkGeometry,
    StitchFaces,
    Destroy,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct Inputs {
    brep: BrepIn,
    op: OpIn,
    #[serde(default = "yes")]
    editing: bool,
    /// The target edge or vertex index in walk order.
    #[serde(default)]
    target: i64,
    /// The face whose surface to keep (`DeleteEdge`) or to stitch, or `-1`.
    #[serde(default = "minus_one")]
    keep: i64,
    #[serde(default)]
    tol: f64,
    #[serde(default)]
    flag: bool,
}

fn yes() -> bool {
    true
}
fn minus_one() -> i64 {
    -1
}
fn err(e: &(impl ToString + ?Sized)) -> HarnessError {
    HarnessError::new(e.to_string())
}

/// The fixture block of `brep` (validate-only off, no cached gaps).
fn block(brep: &BrepIn) -> Result<(String, Vec<f64>), HarnessError> {
    let mut v = Vec::new();
    match brep {
        BrepIn::File { rebuild, .. } => v.extend([0.0, f64::from(u8::from(*rebuild))]),
        BrepIn::Wire { curve, t0, t1 } => {
            v.push(1.0);
            shared::push_curve(&mut v, curve)?;
            v.extend([*t0, *t1]);
        }
        BrepIn::Face { surface, dom } => {
            v.push(2.0);
            shared::push_surface(&mut v, surface)?;
            v.extend(dom);
        }
        BrepIn::Network => v.push(3.0),
    }
    v.extend([0.0, 0.0]);
    Ok((shared::path_of(brep)?, v))
}

fn request(i: &Inputs) -> Request {
    Request {
        op: match i.op {
            OpIn::DeleteEdge => Op::DeleteEdge,
            OpIn::DeleteVertex => Op::DeleteVertex,
            OpIn::ShrinkGeometry => Op::ShrinkGeometry,
            OpIn::StitchFaces => Op::StitchFaces,
            OpIn::Destroy => Op::Destroy,
        },
        editing: i.editing,
        target: i.target,
        keep: i.keep,
        tol: i.tol,
        flag: i.flag,
    }
}

/// The imported fixture with the vertex use lists in their raw order.
fn import(f: &solidcore_ffi::units::topology_topology_traverser::Fixture) -> Result<Brep, HarnessError> {
    let (mut b, _) = shared::import(&f.graph)?;
    for (i, list) in f.vertexuses.iter().enumerate() {
        b.vertex_mut(VertexId::from_index(i as u32))
            .map_err(|e| err(&e))?
            .vertexuses = list.iter().map(|&j| VertexuseId::from_index(j)).collect();
    }
    Ok(b)
}

/// The record counts of the shim's walk, in its kind order.
fn record_counts(b: &Brep) -> Vec<i64> {
    [
        b.regions().len(),
        b.shells().len(),
        b.faceuses().len(),
        b.faces().len(),
        b.loopuses().len(),
        b.loops().len(),
        b.edgeuses().len(),
        b.edges().len(),
        b.vertexuses().len(),
        b.vertices().len(),
    ]
    .iter()
    .map(|&n| n as i64)
    .collect()
}

fn status_output(status: f64) -> Output {
    Output::new().values("status", vec![status])
}

fn rust_status(e: &SolidError) -> Result<Output, HarnessError> {
    if let Some(o) = Output::unported(e) {
        return Ok(o);
    }
    match e {
        SolidError::Upstream { code, .. } => Ok(status_output(f64::from(*code))),
        other => Err(err(&format!("rust: {other}"))),
    }
}

fn key_of(e: EntityRef) -> ObjectRef {
    let (t, i) = match e {
        EntityRef::Edge(x) => (1u64, u64::from(x.index())),
        EntityRef::Vertex(x) => (2, u64::from(x.index())),
        EntityRef::Edgeuse(x) => (3, u64::from(x.index())),
        EntityRef::Vertexuse(x) => (4, u64::from(x.index())),
        EntityRef::Face(x) => (5, u64::from(x.index())),
        _ => (6, 0),
    };
    ObjectRef {
        type_code: 0,
        key: ObjectKey((t << 32) | i),
    }
}

impl UnitDiff for Runner {
    fn id(&self) -> &'static str {
        "topology.brep_p4"
    }

    fn run_cpp(&self, case: &Case) -> Result<Output, HarnessError> {
        let i: Inputs = case.parse()?;
        let (path, blk) = block(&i.brep)?;
        let out = ffi::eval(&path, &blk, request(&i)).map_err(|e| err(&e))?;
        let Some((&status, rest)) = out.split_first() else {
            return Err(err("empty shim output"));
        };
        if status == -1.0 {
            return Ok(Output::rejected("fixture B-rep did not build"));
        }
        if i.op == OpIn::Destroy {
            return Ok(status_output(status).values("doing_boolean", rest.to_vec()));
        }
        let mut o = status_output(status);
        if status == 0.0 {
            let (extra, cnt) = rest.split_at(rest.len().saturating_sub(10));
            if i.op == OpIn::StitchFaces {
                o = o.values("stitch", extra.to_vec());
            }
            o = o.ints("counts", cnt.iter().map(|&x| x as i64).collect::<Vec<_>>());
        }
        Ok(o)
    }

    fn run_rust(&self, case: &Case) -> Result<Output, HarnessError> {
        let i: Inputs = case.parse()?;
        let (path, blk) = block(&i.brep)?;
        let f = ffi::export(&path, &blk).map_err(|e| err(&e))?;
        if counts(&f.graph).is_none() {
            return Ok(Output::new());
        }
        let mut b = import(&f)?;
        b.set_editing_enabled(i.editing);
        let mut k = Kernel::default();
        let mut records: BTreeMap<ObjectKey, Topology> = BTreeMap::new();
        let mut attributes = AttributeStore::default();
        let mut env = DeleteEnv {
            records: &mut records,
            attributes: &mut attributes,
            key_of: &key_of,
        };
        let keep = u32::try_from(i.keep).ok().map(FaceId::from_index);
        let target = u32::try_from(i.target).unwrap_or(u32::MAX);
        let r = match i.op {
            OpIn::DeleteEdge => {
                p4::delete_edge(&mut k, &mut b, &mut env, EdgeId::from_index(target), keep).map(|()| None)
            }
            OpIn::DeleteVertex => {
                p4::delete_vertex(&mut k, &mut b, &mut env, VertexId::from_index(target)).map(|()| None)
            }
            OpIn::ShrinkGeometry => p4::shrink_geometry(&k, &mut b).map(|()| None),
            OpIn::StitchFaces => {
                let faces = [keep];
                p4::stitch_faces(&mut b, i.tol, keep.is_some().then_some(&faces[..]), i.flag).map(|s| {
                    Some(vec![
                        f64::from(s.stitched_edges),
                        f64::from(s.lamina_edges),
                        s.max_vertex_edge_gap,
                        s.max_edge_face_trim_curve_gap,
                    ])
                })
            }
            OpIn::Destroy => {
                return match p4::destroy(&mut k, &mut b) {
                    Ok(()) => Ok(status_output(0.0).values(
                        "doing_boolean",
                        vec![f64::from(u8::from(k.context.doing_boolean()))],
                    )),
                    Err(e) => rust_status(&e),
                };
            }
        };
        match r {
            Ok(extra) => {
                let mut o = status_output(0.0);
                if let Some(x) = extra {
                    o = o.values("stitch", x);
                }
                Ok(o.ints("counts", record_counts(&b)))
            }
            Err(e) => rust_status(&e),
        }
    }

    fn generate(&self, rng: &mut Rng) -> Option<Value> {
        let brep = shared::gen_brep(rng);
        let (path, blk) = block(&brep).ok()?;
        let f = ffi::export(&path, &blk).ok()?;
        let c = counts(&f.graph)?;
        let op = match rng.below(10) {
            0..=4 => OpIn::DeleteEdge,
            5..=6 => OpIn::DeleteVertex,
            7 => OpIn::ShrinkGeometry,
            8 => OpIn::StitchFaces,
            _ => OpIn::Destroy,
        };
        let n = match op {
            OpIn::DeleteEdge => c[7],
            OpIn::DeleteVertex => c[9],
            _ => 1,
        };
        if n == 0 {
            return None;
        }
        let target = rng.below(n as u64) as i64;
        let keep = if c[3] > 0 && rng.chance(0.4) {
            rng.below(c[3] as u64) as i64
        } else {
            -1
        };
        let tol = *rng.choose(&[1.0e-6, 1.0e-3, 0.1, 0.0])?;
        let i = Inputs {
            brep,
            op,
            editing: rng.below(8) != 0,
            target,
            keep,
            tol,
            flag: rng.chance(0.5),
        };
        serde_json::to_value(i).ok()
    }
}
