//! Runner for unit `blend.fillet_corner_p2`: `sm_FindConcaveEdgeuses` and the corner
//! methods ported so far (`AdjustConcaveCornerTopology`, `UpdateSplitFace`).
//!
//! C++: `solidcore_ffi::units::blend_fillet_corner_p2` (its shim is included in
//! `blend_fillet_corner_p1.cpp`, which compiles `source/SMLib/src/SmFilletCorner.cpp`).
//! Rust: `solidcore::blend::fillet_corner::p2`.
//!
//! **Corner cases** (`corner`) use unit `blend.fillet_corner_p3`'s corner harness (its
//! model, corner and shim, `sc_blend_fillet_corner_p3_corner`): an Nx2-shaped model
//! (`gen_nx2_model` there) with a corner whose method is
//! - `adjust_concave`: the virtual `AdjustTopology` of an Nx2 corner, which with three
//!   corner edges is `AdjustConcaveCornerTopology` (also other kinds and edge counts,
//!   which take the other branches). The shim's pre-check refuses where upstream
//!   dereferences null and does not run past the gap test (the first
//!   `SmFilletEdge::MakeFilletEdge`, `blend.fillet_geom_p5`: deferred);
//! - `update_split_face`: `UpdateSplitFace(orig, split)` on target faces, some fillet
//!   vertices classified on (or with their other face) a target face and some fillet
//!   edges with an original face. A fillet edge on `orig` reaches
//!   `SmFilletEdge::ShouldSwapFaces` (`blend.fillet_geom_p6`: deferred, not run in C++).
//!
//! Compared: the status, the corner's members and the model's whole export. 80 built-in
//! corner cases (a fixed seed) plus generated ones.
//!
//! **Cases.** `find_concave` on a real B-rep (unit `topology.edgeuse_p1`'s B-rep block:
//! the box-like fillet test files of `TestFiles/pt_TestFiles/Fillet` that unit
//! `blend.fillet_corner_p1` uses, with convex and concave corners, and that unit's
//! generated B-reps). The C++ side builds the B-rep, the Rust side imports the shim's
//! export. A vertex is selected; the edge uses are either the corner walk
//! `SmFilletCorner::Create` makes (from an edge use starting at the vertex: its clockwise
//! edge use's radial, skipping edges already seen, until the walk closes), or a selection
//! of the edge uses starting at the vertex (in export order, repeats allowed), sometimes
//! with edge uses of elsewhere. Built-in: every fillet file at its first 12 vertices, as
//! the walk and as the export ring, and the wire network of `prog_test` `my_test_topology`
//! at 4 vertices (its wire edge uses' radials). Compared: the concave pair (export
//! indices), or a count of 0 when no pair is concave.
//!
//! The shim refuses (`-2`) when a pair the loop reaches has an edge use upstream would
//! dereference null; the Rust side returns [`SolidError::InvalidInput`] there (whatever
//! its reason: a missing radial, or a missing edge as unit `topology.edgeuse_p1` reports
//! it) and both report `refused`. No B-rep the harness builds reaches it (every edge use
//! has a radial and an edge); the Rust half is covered by a unit test in
//! `fillet_corner::p2`. A case whose evaluated geometry the import holds only as its
//! B-spline part is deferred as unit `blend.fillet_corner_p1`'s runner defers it.

use serde::{Deserialize, Serialize};
use serde_json::Value;

use solidcore::blend::fillet_corner::p2 as corner;
use solidcore::error::SolidError;
use solidcore::kernel::Kernel;
use solidcore::topology::edgeuse::p1 as eu1;
use solidcore::topology::{Brep, EdgeuseId, VertexId};
use solidcore_ffi::units::blend_fillet_corner_p2 as ffi;
use solidcore_ffi::units::blend_fillet_corner_p3::Op as CornerOp;
use solidcore_ffi::units::topology_edgeuse_p1 as ffi1;

use super::blend_fillet_corner_p3::{
    corner_cpp_op, corner_rust_with, gen_nx2_model, CornerIn, EuSel,
};
use super::blend_fillet_geom_p1::Model;

use super::blend_fillet_corner_p1::{
    eu_vertex, gen_topo_brep, inexact, rust_err, Topo, CORNER_FILES,
};
use super::topology_edgeuse_p1::{brep_block_of, harness, import_with_types, path_of, sel, BrepIn};
use crate::case::CaseContext;
use crate::runner::default_cases;
use crate::{Case, HarnessError, Output, Rng, UnitDiff};

/// Registered runner.
pub struct Runner;

/// The upper bound of a corner walk (`SmFilletCorner::Create` stops when it closes).
const MAX_WALK: usize = 64;

/// Built-in corner cases (from [`CORNER_SEED`]).
const CORNER_CASES: usize = 80;

/// The seed of the built-in corner cases.
const CORNER_SEED: u64 = 0x00c0_47e2;

/// The corner method of a `corner` case.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
enum CornerMethod {
    /// The virtual `AdjustTopology` (`AdjustConcaveCornerTopology` for an Nx2 corner with
    /// three corner edges).
    AdjustConcave,
    /// `SmFilletCorner::UpdateSplitFace(orig, split)`.
    UpdateSplitFace,
}

/// One case.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "op", rename_all = "snake_case")]
enum Inputs {
    /// `sm_FindConcaveEdgeuses(v, all, ..)`. `walk`: `all` is the corner walk from the
    /// edge use `eus[0]` selects; else each of `eus` selects an edge use starting at the
    /// vertex (a negative one, `-1 - s`, any edge use of the B-rep), and no `eus` is all
    /// of them in export order.
    FindConcave {
        brep: BrepIn,
        v: i64,
        eus: Vec<i64>,
        walk: bool,
    },
    /// A corner method on unit `blend.fillet_corner_p3`'s corner harness; `split` are
    /// the target faces of `update_split_face`.
    Corner {
        method: CornerMethod,
        model: Model,
        corner: CornerIn,
        #[serde(default)]
        split: [usize; 2],
    },
}

/// The corner walk of `SmFilletCorner::Create` from `start`: the clockwise edge use's
/// radial, skipping an edge use whose edge the step started on or one already seen,
/// until it returns to `start`; `None` when a link is missing or it does not close.
fn corner_walk(b: &Brep, start: EdgeuseId) -> Option<Vec<EdgeuseId>> {
    let edge = |eu| eu1::edge(b, eu).ok().flatten();
    let mut all = vec![start];
    let mut eu = start;
    for _ in 0..MAX_WALK {
        let e0 = eu;
        let e1 = eu1::cw_edgeuse(b, e0).ok().flatten()?;
        eu = eu1::radial(b, e1).ok().flatten()?;
        if eu == start {
            return Some(all);
        }
        if edge(e0) == edge(e1) {
            continue;
        }
        if all.iter().any(|&p| edge(p) == edge(eu)) {
            continue;
        }
        all.push(eu);
    }
    None
}

/// Builds and imports the B-rep and resolves the vertex and edge uses: `args` is
/// `v eu..`. `Err(output)` when the case ends before the op, the same on both sides.
fn topo(i: &Inputs) -> Result<Result<Topo, Output>, HarnessError> {
    let Inputs::FindConcave { brep, v, eus, walk } = i else {
        return Err(harness("not a find_concave case"));
    };
    let path = path_of(Some(brep))?;
    let block = brep_block_of(Some(brep), false, false)?;
    let ex = ffi1::eval(ffi1::Op::Export, &path, &block)
        .map_err(|e| harness(format!("c++ export: {e}")))?;
    let Some(c) = ffi1::counts(&ex) else {
        let st = ex.get(1).copied().unwrap_or(0.0);
        return Ok(Err(Output::err(format!("brep_build_status_{st}"))));
    };
    let (n_eu, n_v) = (c[6], c[9]);
    let (b, types) = import_with_types(&ex)?;
    let vv = sel(*v, n_v);
    if vv < 0 {
        return Ok(Err(Output::err("no_vertex")));
    }
    let vid = VertexId::from_index(vv as u32);
    // the edge uses starting at the vertex (export order)
    let ring: Vec<EdgeuseId> = b
        .edgeuses()
        .iter()
        .map(|(id, _)| id)
        .filter(|&id| eu_vertex(&b, id) == Some(vid))
        .collect();
    if ring.is_empty() {
        return Ok(Err(Output::err("no_edge_use")));
    }
    let from_ring = |s: i64| ring.get(sel(s, ring.len()).max(0) as usize).copied();
    let chosen: Vec<EdgeuseId> = if *walk {
        let start = eus.first().copied().unwrap_or(0).rem_euclid(1 << 30);
        let Some(start) = from_ring(start) else {
            return Err(harness("ring selection"));
        };
        match corner_walk(&b, start) {
            Some(w) => w,
            None => return Ok(Err(Output::err("walk_open"))),
        }
    } else if eus.is_empty() {
        ring.clone()
    } else {
        eus.iter()
            .filter_map(|&s| {
                if s < 0 {
                    let r = sel(-1 - s, n_eu);
                    (r >= 0).then(|| EdgeuseId::from_index(r as u32))
                } else {
                    from_ring(s)
                }
            })
            .collect()
    };
    if chosen.is_empty() {
        return Ok(Err(Output::err("no_edge_use")));
    }
    let mut args = vec![vv];
    args.extend(chosen.iter().map(|e| i64::from(e.index())));
    Ok(Ok(Topo {
        path,
        block,
        brep: b,
        types,
        args,
    }))
}

fn corner_cpp(
    method: CornerMethod,
    model: &Model,
    c: &CornerIn,
    split: [usize; 2],
) -> Result<Output, HarnessError> {
    match method {
        CornerMethod::AdjustConcave => corner_cpp_op(CornerOp::AdjustTopology, &[], model, c),
        CornerMethod::UpdateSplitFace => corner_cpp_op(
            CornerOp::UpdateSplitFace,
            &[split[0] as f64, split[1] as f64],
            model,
            c,
        ),
    }
}

fn corner_rust(
    method: CornerMethod,
    model: &Model,
    c: &CornerIn,
    split: [usize; 2],
) -> Result<Output, HarnessError> {
    corner_rust_with(model, c, |c, env, faces| match method {
        CornerMethod::AdjustConcave => c.adjust_topology(env),
        CornerMethod::UpdateSplitFace => {
            let (Some(&orig), Some(&new)) = (faces.get(split[0]), faces.get(split[1])) else {
                return Err(SolidError::invalid_input("difftest", "split face index"));
            };
            c.update_split_face(env, orig, new)
        }
    })
}

fn run_cpp(i: &Inputs) -> Result<Output, HarnessError> {
    if let Inputs::Corner {
        method,
        model,
        corner,
        split,
    } = i
    {
        return corner_cpp(*method, model, corner, *split);
    }
    let t = match topo(i)? {
        Ok(t) => t,
        Err(o) => return Ok(o),
    };
    let mut v = t.block.clone();
    let (&vv, eus) = t
        .args
        .split_first()
        .ok_or_else(|| harness("missing vertex"))?;
    v.push(vv as f64);
    v.push(eus.len() as f64);
    v.extend(eus.iter().map(|&a| a as f64));
    let out = ffi::find_concave(&t.path, &v).map_err(|e| harness(format!("c++: {e}")))?;
    let Some((&head, rest)) = out.split_first() else {
        return Err(harness("empty output"));
    };
    if head != 0.0 {
        let st = rest.first().copied().unwrap_or(0.0);
        return Ok(Output::err(format!("brep_build_status_{st}")));
    }
    if rest.first().copied() == Some(-2.0) {
        return Ok(Output::err("refused"));
    }
    Ok(Output::new().values("concave", rest.to_vec()))
}

fn run_rust(i: &Inputs) -> Result<Output, HarnessError> {
    if let Inputs::Corner {
        method,
        model,
        corner,
        split,
    } = i
    {
        return corner_rust(*method, model, corner, *split);
    }
    let mut t = match topo(i)? {
        Ok(t) => t,
        Err(o) => return Ok(o),
    };
    let (&vv, eus) = t
        .args
        .split_first()
        .ok_or_else(|| harness("missing vertex"))?;
    let v = VertexId::from_index(u32::try_from(vv).map_err(|_| harness("index"))?);
    let eus = eus
        .iter()
        .map(|&a| u32::try_from(a).map(EdgeuseId::from_index))
        .collect::<Result<Vec<_>, _>>()
        .map_err(|_| harness("index"))?;
    if let Some(o) = inexact(&t, &[], &eus) {
        return Ok(o);
    }
    let k = Kernel::default();
    match corner::find_concave_edgeuses(&k, &mut t.brep, v, &eus) {
        Ok(c) => {
            let mut o = vec![c.len() as f64];
            o.extend(c.iter().map(|e| f64::from(e.index())));
            Ok(Output::new().values("concave", o))
        }
        // every InvalidInput of the op is a null dereference upstream (the empty list is
        // screened above): a missing radial (`NULL_DEREF`) or, from unit
        // `topology.edgeuse_p1`, a missing edge; the shim refuses both
        Err(SolidError::InvalidInput { .. }) => Ok(Output::err("refused")),
        Err(e) => Ok(rust_err(&e)),
    }
}

fn gen_find_concave(rng: &mut Rng) -> Inputs {
    let brep = gen_topo_brep(rng);
    let walk = rng.chance(0.5);
    let n = 1 + rng.below(6) as usize;
    let eus = (0..n)
        .map(|_| {
            if !walk && rng.chance(0.1) {
                -1 - rng.range_i64(0, 1 << 20)
            } else {
                rng.range_i64(0, 1 << 20)
            }
        })
        .collect();
    Inputs::FindConcave {
        brep,
        v: rng.range_i64(0, 1 << 20),
        eus,
        walk,
    }
}

/// A corner case on an Nx2-shaped model (`gen_nx2_model` of unit
/// `blend.fillet_corner_p3`): its target faces are `0` (plane x-z) and `1` (plane y-z),
/// its corner edges `4 5` (plus `2`, a rail, for the three-edge concave case), its
/// fillet vertices `0` (rails meet the side edge), `3`, `4`.
fn gen_corner(rng: &mut Rng) -> Inputs {
    let scale = *rng.choose(&[1.0, 1e-3, 1e3, 0.1, 1e5]).unwrap_or(&1.0);
    let tol = match rng.below(4) {
        0 => 1e-3 * scale,
        1 => 1e-6 * scale,
        2 => 0.0,
        _ => scale * 10f64.powf(rng.uniform(-9.0, -2.0)),
    };
    let (mut model, tsects) = gen_nx2_model(rng, scale);
    let eu = |e: usize, w: usize| -> EuSel { [e, w] };
    let n_all = if rng.chance(0.85) {
        4
    } else {
        rng.below(6) as usize
    };
    let all: Vec<EuSel> = (0..n_all).map(|i| eu(i % 4, usize::from(i >= 4))).collect();
    let solver = Some(vec![eu(1, 0), eu(2, rng.below(2) as usize)]);
    let mut verts: Vec<usize> = vec![0, 3, 4];
    if rng.chance(0.2) {
        verts.reverse();
    }
    if rng.chance(0.1) {
        verts = vec![3, 4];
    }
    let method = if rng.chance(0.5) {
        CornerMethod::AdjustConcave
    } else {
        CornerMethod::UpdateSplitFace
    };
    let (kind, edges, mut orig_faces, split) = match method {
        CornerMethod::AdjustConcave => {
            let kind = if rng.chance(0.9) {
                7
            } else {
                rng.below(11) as usize
            };
            let edges = if rng.chance(0.85) {
                vec![4, 5, 2]
            } else {
                vec![4, 5]
            };
            (kind, edges, Vec::new(), [0, 1])
        }
        CornerMethod::UpdateSplitFace => {
            // fillet vertices classified on a target face (or with it as their other face)
            for v in model.vertices.iter_mut().skip(1) {
                if rng.chance(0.6) {
                    v.pc_class = 102; // SM_PC_FACE
                    v.pc_kind = 3;
                    v.pc_index = rng.below(2) as i64;
                }
                if rng.chance(0.3) {
                    v.pc_face = rng.below(2) as i64;
                }
            }
            let orig = rng.below(2) as usize;
            let split = [orig, if rng.chance(0.9) { 1 - orig } else { orig }];
            let mut orig_faces = Vec::new();
            for e in [4usize, 5] {
                if rng.chance(0.3) {
                    orig_faces.push([e, rng.below(2) as usize]);
                }
            }
            (rng.below(11) as usize, vec![4, 5], orig_faces, split)
        }
    };
    orig_faces.dedup_by_key(|p| p[0]);
    let blending = rng.chance(0.25);
    Inputs::Corner {
        method,
        model,
        corner: CornerIn {
            kind,
            v_edge: if rng.chance(0.95) { 1 } else { -1 },
            v_sub: 0,
            solver,
            all: if rng.chance(0.95) { Some(all) } else { None },
            verts,
            edges,
            bevel: !blending,
            blending,
            chamfer: false,
            degenerate: false,
            ext: rng.chance(0.8),
            adjusted: false,
            setback: 0.0,
            tol,
            tan: rng.uniform(1e-4, 0.2),
            tsects,
            orig_faces,
            line_edges: Vec::new(),
            arc_edges: Vec::new(),
        },
        split,
    }
}

impl UnitDiff for Runner {
    fn id(&self) -> &'static str {
        "blend.fillet_corner_p2"
    }

    fn run_cpp(&self, case: &Case) -> Result<Output, HarnessError> {
        run_cpp(&case.parse::<Inputs>()?)
    }

    fn run_rust(&self, case: &Case) -> Result<Output, HarnessError> {
        run_rust(&case.parse::<Inputs>()?)
    }

    fn generate(&self, rng: &mut Rng) -> Option<Value> {
        let i = if rng.chance(0.5) {
            gen_find_concave(rng)
        } else {
            gen_corner(rng)
        };
        serde_json::to_value(i).ok()
    }

    fn cases(&self, ctx: &CaseContext) -> Result<Vec<Case>, HarnessError> {
        let mut v = default_cases(self, ctx)?;
        let push = |v: &mut Vec<Case>, name: String, upstream_ref: &str, i: Inputs| {
            let index = v.len();
            v.push(Case {
                unit: "blend.fillet_corner_p2".to_string(),
                name,
                source: crate::CaseSource::Fixture {
                    path: "built-in: blend.fillet_corner_p2 cases".to_string(),
                    index,
                    upstream_ref: Some(upstream_ref.to_string()),
                },
                inputs: serde_json::to_value(i).map_err(|e| harness(e.to_string()))?,
            });
            Ok::<(), HarnessError>(())
        };
        for file in CORNER_FILES {
            for vertex in 0..12 {
                for walk in [true, false] {
                    let i = Inputs::FindConcave {
                        brep: BrepIn::File {
                            file: file.to_string(),
                            rebuild: false,
                        },
                        v: vertex,
                        // the walk from the first edge use of the ring, or the whole ring
                        // in export order
                        eus: if walk { vec![0] } else { Vec::new() },
                        walk,
                    };
                    let name = format!(
                        "find_concave_{}_v{vertex}_{}",
                        file.trim_start_matches("Fillet/").trim_end_matches(".smb"),
                        if walk { "walk" } else { "ring" }
                    );
                    push(&mut v, name, "TestFiles/pt_TestFiles/Fillet", i)?;
                }
            }
        }
        // Wire edge uses (prog_test my_test_topology's network): their radial is the shell
        // edge use's next one, a branch the fillet files' loop edge uses never take.
        for vertex in 0..4 {
            let i = Inputs::FindConcave {
                brep: BrepIn::Network,
                v: vertex,
                eus: Vec::new(),
                walk: false,
            };
            let name = format!("find_concave_network_v{vertex}_ring");
            push(&mut v, name, "tests/prog_test (my_test_topology)", i)?;
        }
        let mut rng = Rng::new(CORNER_SEED);
        for n in 0..CORNER_CASES {
            let i = gen_corner(&mut rng);
            let Inputs::Corner { method, .. } = &i else {
                return Err(harness("corner case"));
            };
            let name = match method {
                CornerMethod::AdjustConcave => format!("corner_adjust_concave_{n:02}"),
                CornerMethod::UpdateSplitFace => format!("corner_update_split_face_{n:02}"),
            };
            push(
                &mut v,
                name,
                "SmFilletCorner.cpp (Nx2-shaped corner model)",
                i,
            )?;
        }
        Ok(v)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn inputs_round_trip() {
        let mut rng = Rng::new(7);
        for _ in 0..100 {
            let v = Runner.generate(&mut rng).unwrap();
            let back: Inputs =
                serde_json::from_value(v.clone()).unwrap_or_else(|e| panic!("{e}: {v}"));
            assert_eq!(serde_json::to_value(back).unwrap(), v);
        }
    }
}
