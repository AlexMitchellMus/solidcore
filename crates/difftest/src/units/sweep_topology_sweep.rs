//! Stateful Rust/C++ parity for topology sweep maps, the stitch callback, base
//! geometry dispatch and vertex helpers. Every script compares all four maps,
//! including null entries, reverse-map leftovers and both B-rep namespaces.
//! Vertex construction is deferred at the exact missing B-rep constructor;
//! whole-B-rep geometry construction and the tails behind those calls are not
//! claimed as compared coverage.
#![allow(clippy::indexing_slicing)]
use crate::{Case, CaseContext, HarnessError, Output, Rng, UnitDiff};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use solidcore::geom::{
    attribute::{ids, Attribute, AttributeStore, AttributeValue, Behavior, ObjectKey, TagData},
    context::MarkType,
};
use solidcore::kernel::Kernel;
use solidcore::math::Point3;
use solidcore::sweep::topology_sweep::access::{EdgeRef, FaceRef, SweepAccess};
use solidcore::sweep::topology_sweep::geometry::SweepGeometry;
use solidcore::sweep::topology_sweep::{SweepEntity as E, SweepGeometryCreation, TopologySweep};
use solidcore::sweep::{
    rotational_sweep_geometry::RotationalSweepGeometry,
    trans_sweep_geometry::TranslationalSweepGeometry,
};
use solidcore::topology::{Brep, BrepId, EdgeId, FaceId, Region, RegionId, Vertex, VertexId};
use solidcore::{Result, SolidError};
use std::collections::BTreeMap;

/// `sweep.topology_sweep` runner.
pub struct Runner;
#[derive(Debug, Clone, Serialize, Deserialize)]
struct Input {
    steps: Vec<[i64; 4]>,
    #[serde(default)]
    kind: u8,
    #[serde(default = "vector")]
    vector: [f64; 3],
}
fn vector() -> [f64; 3] {
    [1., 0., 0.]
}
fn entity(id: i64) -> Option<E> {
    if id == 64 {
        return None;
    }
    let b = BrepId::from_index(((id % 16) / 8) as u32);
    let i = (id % 8) as u32;
    Some(match id / 16 {
        0 => E::Vertex(b, VertexId::from_index(i)),
        1 => E::Edge(b, EdgeId::from_index(i)),
        2 => E::Face(b, FaceId::from_index(i)),
        _ => E::Region(b, RegionId::from_index(i)),
    })
}
fn entity_id(e: Option<E>) -> i64 {
    match e {
        None => 64,
        Some(E::Vertex(b, i)) => i64::from(b.index() * 8 + i.index()),
        Some(E::Edge(b, i)) => 16 + i64::from(b.index() * 8 + i.index()),
        Some(E::Face(b, i)) => 32 + i64::from(b.index() * 8 + i.index()),
        Some(E::Region(b, i)) => 48 + i64::from(b.index() * 8 + i.index()),
    }
}
fn status(e: &SolidError) -> i64 {
    match e {
        SolidError::Upstream { code, .. } => i64::from(*code),
        SolidError::InvalidInput { .. } => 1007,
        SolidError::Nlib { .. } => 1001,
        _ => -1,
    }
}
fn unsupported<T>() -> Result<T> {
    Err(SolidError::invalid_input(
        "runner",
        "fixture has no service for this path",
    ))
}
struct Model {
    k: Kernel,
    b: Vec<Brep>,
    attributes: AttributeStore,
    keys: BTreeMap<E, ObjectKey>,
    curve_keys: BTreeMap<EdgeRef, ObjectKey>,
}
impl Model {
    fn new(vector: [f64; 3]) -> Result<Self> {
        let mut b = vec![Brep::new(), Brep::new()];
        for brep in &mut b {
            brep.add_region_record(Region { in_brep: true })?;
            for i in 0..8 {
                brep.add_vertex_record(Vertex {
                    point: Point3::new(
                        vector[0] * f64::from(i),
                        vector[1] * f64::from(i),
                        vector[2] * f64::from(i),
                    ),
                    ..Vertex::default()
                })?;
            }
        }
        let mut attributes = AttributeStore::new();
        let mut keys = BTreeMap::new();
        let mut curve_keys = BTreeMap::new();
        for id in 0..64 {
            keys.insert(
                entity(id).ok_or_else(|| SolidError::null_pointer("runner"))?,
                ObjectKey(10_000 + id as u64 * 97),
            );
        }
        for (bid, brep) in b.iter_mut().enumerate() {
            for index in 0..8 {
                let bid = BrepId::from_index(bid as u32);
                let id = 16 + i64::from(bid.index()) * 8 + index;
                let edge = brep.add_edge_record(solidcore::topology::Edge {
                    vertices: [VertexId::from_index(0); 2],
                    edgeuses: Vec::new(),
                    interval: solidcore::geom::extent1d::Extent1d::from_min_max(0., 1.),
                    tolerance: 0.,
                    curve: Some(solidcore::curve::curve::Curve::from(
                        solidcore::curve::line::Line::create_line_segment(
                            3,
                            Point3::ZERO,
                            Point3::X,
                            None,
                        )?,
                    )),
                })?;
                let key = ObjectKey(40_000 + id as u64 * 59);
                curve_keys.insert((bid, edge), key);
                if id % 2 != 0 {
                    let h = attributes.insert(Attribute::tag(
                        ids::TAG,
                        TagData {
                            primary_id: 700 + id,
                            primary_ordinal: 12 + id,
                            ..TagData::default()
                        },
                        Behavior::Copy,
                    ))?;
                    attributes.add_attribute(key, h)?;
                }
            }
        }
        Ok(Self {
            k: Kernel::default(),
            b,
            attributes,
            keys,
            curve_keys,
        })
    }
}
impl SweepAccess for Model {
    fn brep(&self, id: BrepId) -> Result<&Brep> {
        self.b
            .get(id.index() as usize)
            .ok_or_else(|| SolidError::null_pointer("runner"))
    }
    fn brep_mut(&mut self, id: BrepId) -> Result<&mut Brep> {
        self.b
            .get_mut(id.index() as usize)
            .ok_or_else(|| SolidError::null_pointer("runner"))
    }
    fn kernel(&self, _id: BrepId) -> Result<&Kernel> {
        Ok(&self.k)
    }
    fn kernel_mut(&mut self, _id: BrepId) -> Result<&mut Kernel> {
        Ok(&mut self.k)
    }
    fn new_brep(&mut self, _id: BrepId) -> Result<BrepId> {
        unsupported()
    }
    fn delete_merge(&mut self, _id: solidcore::sweep::topology_sweep::access::SweepMergeId) {}
    fn delete_brep(&mut self, _id: BrepId) {}
    fn tag(&self, e: E) -> Result<Option<TagData>> {
        let key = self
            .keys
            .get(&e)
            .ok_or_else(|| SolidError::null_pointer("runner"))?;
        Ok(self
            .attributes
            .find_attribute(*key, ids::TAG)
            .and_then(|h| self.attributes.get(h))
            .and_then(|a| {
                if let AttributeValue::Tag(t) = a.value() {
                    Some(*t)
                } else {
                    None
                }
            }))
    }
    fn add_tag(&mut self, e: E, tag: TagData) -> Result<()> {
        let key = *self
            .keys
            .get(&e)
            .ok_or_else(|| SolidError::null_pointer("runner"))?;
        let h = self
            .attributes
            .insert(Attribute::tag(ids::TAG, tag, Behavior::Copy))?;
        self.attributes.add_attribute(key, h)
    }
    fn curve_tag(&self, edge: EdgeRef) -> Result<Option<TagData>> {
        let key = self
            .curve_keys
            .get(&edge)
            .ok_or_else(|| SolidError::null_pointer("runner"))?;
        Ok(self
            .attributes
            .find_attribute(*key, ids::TAG)
            .and_then(|h| self.attributes.get(h))
            .and_then(|a| {
                if let AttributeValue::Tag(t) = a.value() {
                    Some(*t)
                } else {
                    None
                }
            }))
    }
    fn set_editing(&mut self, _id: BrepId, _editing: bool) -> Result<bool> {
        unsupported()
    }
    fn lock_mark(&mut self, _id: BrepId, _mask: i32) -> Result<MarkType> {
        unsupported()
    }
    fn unlock_mark(&mut self, _id: BrepId, _mark: MarkType) -> Result<()> {
        unsupported()
    }
    fn mark(&mut self, _e: E, _mark: MarkType) -> Result<()> {
        unsupported()
    }
    fn is_marked(&self, _e: E, _mark: MarkType) -> Result<bool> {
        unsupported()
    }
    fn collect_faces(
        &mut self,
        _f: FaceRef,
        _mark: MarkType,
        _out: &mut Vec<FaceId>,
    ) -> Result<()> {
        unsupported()
    }
}
fn rust(i: &Input) -> Result<Vec<i64>> {
    let mut s = TopologySweep::default();
    let mut out = Vec::new();
    let mut model = Model::new(i.vector)?;
    let g = match i.kind {
        0 => SweepGeometry::Base(SweepGeometryCreation::default()),
        1 => SweepGeometry::Translation(TranslationalSweepGeometry::new(Point3::new(
            i.vector[0],
            i.vector[1],
            i.vector[2],
        ))),
        _ => SweepGeometry::Rotation(RotationalSweepGeometry::new(
            Point3::ZERO,
            Point3::Z,
            90.,
            1e-5,
        )),
    };
    for &[op, ai, bi, arg] in &i.steps {
        let a = entity(ai);
        let b = entity(bi);
        let mut keep = true;
        let mut distance = -123.;
        let result = match op {
            0 => {
                s.set_same(a, b);
                Ok(())
            }
            1 => {
                s.set_higher(a, b);
                Ok(())
            }
            2 => s.select_topology_to_keep(a, b, i.vector[0], &mut distance, &mut keep),
            3 => {
                s.remove_from_to_same(a);
                Ok(())
            }
            4 => {
                s.remove_from_from_same(a);
                Ok(())
            }
            5 => {
                let map = match arg {
                    0 => &mut s.to_same,
                    1 => &mut s.to_higher,
                    2 => &mut s.from_same,
                    _ => &mut s.from_lower,
                };
                map.insert(a, b);
                Ok(())
            }
            6 => {
                s.repetitions = arg as u64;
                Ok(())
            }
            7 => {
                // Both native defaults and the rotational override clear this
                // output before returning an error for too many turns.
                keep = false;
                g.is_sweep_closed(arg as u64).map(|v| keep = v)
            }
            8 => g.start_coord_system(Point3::ZERO).map(|_| ()),
            9 => g.sweep_transform().map(|_| ()),
            10 | 11 => {
                let Some(E::Vertex(b, v)) = a else {
                    return unsupported();
                };
                let region = (BrepId::from_index(0), RegionId::from_index(0));
                if op == 10 {
                    s.sweep_vertices_same(&mut model, &g, region, &[Some((b, v))])
                } else {
                    s.sweep_vertices_higher(&mut model, &g, region, &[Some((b, v))])
                }
            }
            15 => {
                s.do_stitching = false;
                s.do_merge = false;
                s.repetitions = 1;
                s.do_sweep(
                    &mut model,
                    &g,
                    Some(BrepId::from_index(0)),
                    Some((BrepId::from_index(0), RegionId::from_index(0))),
                    &mut solidcore::sweep::topology_sweep::engine::SweepTargets {
                        vertices_higher: Some(vec![Some((
                            BrepId::from_index(0),
                            VertexId::from_index(0),
                        ))]),
                        ..Default::default()
                    },
                )
            }
            14 => {
                let Some(E::Edge(ba, ia)) = a else {
                    return unsupported();
                };
                let Some(E::Edge(bb, ib)) = b else {
                    return unsupported();
                };
                s.tag_profile(
                    &mut model,
                    &solidcore::sweep::topology_sweep::engine::SweepTargets {
                        edges_higher: Some(vec![Some((ba, ia)), Some((bb, ib))]),
                        ..Default::default()
                    },
                )
            }
            13 => s.do_sweep(
                &mut model,
                &g,
                None,
                None,
                &mut solidcore::sweep::topology_sweep::engine::SweepTargets::default(),
            ),
            _ => {
                s.do_tagging = arg % 2 != 0;
                s.tag_id = arg;
                s.do_merge = arg % 2 != 0;
                s.do_stitching = (arg / 2) % 2 != 0;
                Ok(())
            }
        };
        let status = match result {
            Ok(()) => 0,
            Err(err @ SolidError::Unported { .. }) => return Err(err),
            Err(err) => status(&err),
        };
        out.extend([status, i64::from(keep), distance as i64]);
    }
    out.extend([
        i64::from(s.do_tagging),
        s.tag_id,
        s.repetitions as i64,
        i64::from(s.do_merge),
        i64::from(s.do_stitching),
    ]);
    for map in [&s.to_same, &s.to_higher, &s.from_same, &s.from_lower] {
        out.push(map.len() as i64);
        for id in 0..=64 {
            let e = entity(id);
            out.extend([
                i64::from(map.contains_key(&e)),
                entity_id(map.get(&e).copied().flatten()),
            ]);
        }
    }
    for id in 16..32 {
        let tag = model.tag(entity(id).ok_or_else(|| SolidError::null_pointer("runner"))?)?;
        out.push(i64::from(tag.is_some()));
        out.extend(tag.unwrap_or_default().to_array());
    }
    Ok(out)
}
impl UnitDiff for Runner {
    fn id(&self) -> &'static str {
        "sweep.topology_sweep"
    }
    fn run_cpp(&self, c: &Case) -> std::result::Result<Output, HarnessError> {
        let i: Input = c.parse()?;
        solidcore_ffi::units::sweep_topology_sweep::eval(&i.steps, i.vector, i.kind)
            .map(|v| Output::new().ints("state", v))
            .map_err(|e| HarnessError::new(e.to_string()))
    }
    fn run_rust(&self, c: &Case) -> std::result::Result<Output, HarnessError> {
        let i: Input = c.parse()?;
        Ok(match rust(&i) {
            Ok(v) => Output::new().ints("state", v),
            Err(e) => Output::unported(&e).unwrap_or_else(|| Output::err(e.to_string())),
        })
    }
    fn cases(&self, ctx: &CaseContext) -> std::result::Result<Vec<Case>, HarnessError> {
        let mut ctx = ctx.clone();
        ctx.seeds = ctx.seeds.max(512);
        crate::runner::default_cases(self, &ctx)
    }
    fn generate(&self, rng: &mut Rng) -> Option<Value> {
        let mut steps = Vec::new();
        for _ in 0..(4 + rng.below(24)) {
            let op = rng.below(10);
            let a = rng.below(65);
            let b = rng.below(65);
            let arg = if op == 5 {
                rng.below(4)
            } else {
                rng.below(1000)
            };
            steps.push([op as i64, a as i64, b as i64, arg as i64]);
        }
        if rng.below(3) == 0 {
            steps.push([14, 16 + rng.below(16) as i64, 16 + rng.below(16) as i64, 0]);
        }
        if rng.below(4) == 0 {
            steps.push([12, 0, 0, rng.below(1000) as i64]);
        }
        Some(
            json!({"steps":steps,"kind":rng.below(3),"vector":[rng.uniform(-10.,10.),rng.uniform(-10.,10.),rng.uniform(-10.,10.)]}),
        )
    }
}
