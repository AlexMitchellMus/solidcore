//! Runner for unit `blend.fillet_geom_p5`: the rail trims of a fillet geometry
//! (`TrimRailCurves`, `TrimCliffRails`), the cross-section and side-face curves of a
//! fillet edge, the fillet/fillet intersection vertex, the edge on an extended edge, the
//! split-face update and the pseudo-B-rep builders (`MakeFilletEdge`, `MakeRailEdge`,
//! `CreateEndCurve`).
//!
//! C++: `solidcore_ffi::units::blend_fillet_geom_p5` (upstream
//! `source/SMLib/src/SmFilletGeom.cpp`). Rust: `solidcore::blend::fillet_geom::p5`.
//!
//! A case is `blend.fillet_geom_p1`'s pseudo-B-rep model (built the same way on both
//! sides), `blend.fillet_geom_p4`'s base fillet solvers on its geometries (tolerances
//! only) and one query. Both sides produce the same flat record (layout in the FFI module
//! docs), split into keys: `status` (`0` = `SM_SUCCESS`, else the `SmStatus`; `-1` on the
//! Rust side for an error upstream cannot produce) and the query's values: for
//! `trim_rails` per rail its curve (`r{i}_c3d_*`), interval, primary and mate UV curves
//! and both vertex uses' UV points; for `cross_section` the edge curve and both UV curves
//! and every geometry's side-use count; for `side_face` the edge status and the side-use
//! counts; for `int_fillet` the vertex point and status; for the builders the edge count.
//!
//! Deferred on the Rust side (unported callees): every builder past its checks (the edge
//! use constructor needs `SmBrep()`'s infinite region, `topology.brep_p3`). The
//! split-face update and the surface/surface intersections (`blend.fillet_geom_p2`) are
//! wired; the split-face op is refused by the shim (`UpdateSplitFace(NULL, NULL)`
//! dereferences both faces in `ShouldSwapFaces` for every rail). The shim refuses
//! the inputs that would make upstream dereference what the harness cannot build (a real
//! fillet corner, a solver or rail it needs); the generator keeps only inputs both sides
//! accept. Every run has at least [`MIN_GENERATED`] generated cases.

use serde::{Deserialize, Serialize};
use serde_json::Value;

use solidcore::blend::fillet_geom::p1::{BoundaryTrimmingType, FilletSolverId};
use solidcore::blend::fillet_solver::p1::FilletSolver;
use solidcore::curve::curve::{Curve, CurveKind};
use solidcore::kernel::Kernel;
use solidcore::nurbs::b_spline_curve::BSplineCurve;
use solidcore::nurbs::math::NOW;
use solidcore::topology::{edge, edgeuse::p1 as eu1};
use solidcore::SolidError;
use solidcore_ffi::units::blend_fillet_geom_p5::{self as ffi, Op};

use super::blend_fillet_geom_p1::{
    curve_of, dump_bcurve, gen_model, gen_surface, push_model, scale_of, status, surface_of, Crv,
    EdgeIn, EuIn, GeomIn, Model, RustModel, Srf, VertexIn,
};
use crate::case::CaseContext;
use crate::runner::default_cases;
use crate::{Case, HarnessError, Output, Rng, UnitDiff};

/// Generated cases per run, whatever `--seeds` says (it can only raise this).
pub const MIN_GENERATED: usize = 300;

/// The runner.
pub struct Runner;

fn harness(msg: impl Into<String>) -> HarnessError {
    HarnessError::new(msg)
}

// ---- inputs ------------------------------------------------------------------------

/// A base fillet solver's tolerances.
#[derive(Debug, Clone, Serialize, Deserialize)]
struct SolverIn {
    approx: f64,
    ang: f64,
    tang: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "q", rename_all = "snake_case")]
enum Query {
    IntFillet { v: usize },
    CrossSection { e: usize },
    SideFace { e: usize },
    ExtEdge { e: usize },
    TrimRails { g: usize },
    SplitFace { g: usize },
    MakeEdge { vs: i64, ve: i64 },
    RailEdge { g: usize, rail: usize },
    EndCurve { g: usize, end: usize, trim: usize },
}

impl Query {
    const fn op(&self) -> Op {
        match self {
            Self::IntFillet { .. } => Op::IntFillet,
            Self::CrossSection { .. } => Op::CrossSection,
            Self::SideFace { .. } => Op::SideFace,
            Self::ExtEdge { .. } => Op::ExtEdge,
            Self::TrimRails { .. } => Op::TrimRails,
            Self::SplitFace { .. } => Op::SplitFace,
            Self::MakeEdge { .. } => Op::MakeEdge,
            Self::RailEdge { .. } => Op::RailEdge,
            Self::EndCurve { .. } => Op::EndCurve,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct Inputs {
    model: Model,
    #[serde(default)]
    solvers: Vec<SolverIn>,
    /// Per geometry: its solver index, `-1` (or missing) for none.
    #[serde(default)]
    geom_solvers: Vec<i64>,
    query: Query,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct CaseInputs {
    #[serde(flatten)]
    inputs: Inputs,
}

// ---- C++ side ------------------------------------------------------------------------

fn geom_solver(i: &Inputs, g: usize) -> i64 {
    i.geom_solvers.get(g).copied().unwrap_or(-1)
}

fn cpp_input(i: &Inputs) -> Result<Vec<f64>, HarnessError> {
    let mut o = Vec::new();
    push_model(&mut o, &i.model)?;
    o.push(i.solvers.len() as f64);
    for s in &i.solvers {
        o.extend([s.approx, s.ang, s.tang]);
    }
    for g in 0..i.model.geoms.len() {
        o.push(geom_solver(i, g) as f64);
    }
    match &i.query {
        Query::IntFillet { v } => o.push(*v as f64),
        Query::CrossSection { e } | Query::SideFace { e } | Query::ExtEdge { e } => {
            o.push(*e as f64);
        }
        Query::TrimRails { g } | Query::SplitFace { g } => o.push(*g as f64),
        Query::MakeEdge { vs, ve } => o.extend([*vs as f64, *ve as f64]),
        Query::RailEdge { g, rail } => o.extend([*g as f64, *rail as f64]),
        Query::EndCurve { g, end, trim } => o.extend([*g as f64, *end as f64, *trim as f64]),
    }
    Ok(o)
}

fn cpp(i: &Inputs) -> Result<Output, HarnessError> {
    let input = cpp_input(i)?;
    let flat = ffi::eval(i.query.op(), &input).map_err(|e| harness(e.to_string()))?;
    to_output(i, &flat)
}

// ---- output splitting ---------------------------------------------------------------

struct Split<'a> {
    flat: &'a [f64],
    at: usize,
    out: Output,
}

impl Split<'_> {
    fn take(&mut self, k: usize) -> Result<Vec<f64>, HarnessError> {
        let part = self
            .at
            .checked_add(k)
            .and_then(|end| self.flat.get(self.at..end))
            .ok_or_else(|| harness(format!("record too short ({} values)", self.flat.len())))?;
        self.at += k;
        Ok(part.to_vec())
    }

    fn put(&mut self, key: &str, v: Vec<f64>) {
        let out = std::mem::take(&mut self.out);
        self.out = out.values(key, v);
    }

    fn vals(&mut self, key: &str, k: usize) -> Result<(), HarnessError> {
        let v = self.take(k)?;
        self.put(key, v);
        Ok(())
    }

    /// A curve record `has [dim n p m cpts knots]` under `pfx`.
    fn curve(&mut self, pfx: &str) -> Result<(), HarnessError> {
        let has = self.take(1)?;
        let set = has.first().copied() == Some(1.0);
        self.put(&format!("{pfx}has"), has);
        if set {
            let dims = self.take(4)?;
            let cnt = |i: usize| dims.get(i).copied().unwrap_or(f64::NAN);
            let ok = |x: f64| x.fract() == 0.0 && (0.0..=1e6).contains(&x);
            let (n, m) = (cnt(1), cnt(3));
            if !(ok(n) && ok(m)) {
                return Err(harness(format!("bad curve dims {dims:?}")));
            }
            self.put(&format!("{pfx}dims"), dims);
            self.vals(&format!("{pfx}cpts"), 4 * (n as usize + 1))?;
            self.vals(&format!("{pfx}knots"), m as usize + 1)?;
        }
        Ok(())
    }
}

fn to_output(i: &Inputs, flat: &[f64]) -> Result<Output, HarnessError> {
    let mut s = Split {
        flat,
        at: 0,
        out: Output::new(),
    };
    let ng = i.model.geoms.len();
    s.vals("status", 1)?;
    match &i.query {
        Query::IntFillet { .. } => {
            s.vals("point", 3)?;
            s.vals("fil_status", 1)?;
        }
        Query::CrossSection { .. } => {
            s.curve("c3d_")?;
            s.curve("uv1_")?;
            s.curve("uv2_")?;
            s.vals("sides", ng)?;
        }
        Query::SideFace { .. } => {
            s.vals("edge_status", 1)?;
            s.vals("sides", ng)?;
        }
        Query::TrimRails { .. } => {
            for r in 0..2 {
                let has = s.take(1)?;
                let set = has.first().copied() == Some(1.0);
                s.put(&format!("r{r}_has"), has);
                if set {
                    s.curve(&format!("r{r}_c3d_"))?;
                    s.vals(&format!("r{r}_ivl"), 2)?;
                    s.curve(&format!("r{r}_uv1_"))?;
                    s.curve(&format!("r{r}_uv2_"))?;
                    s.vals(&format!("r{r}_vu_uv"), 4)?;
                }
            }
        }
        Query::MakeEdge { .. } | Query::RailEdge { .. } | Query::EndCurve { .. } => {
            s.vals("edges", 1)?;
        }
        Query::ExtEdge { .. } | Query::SplitFace { .. } => {}
    }
    if s.at != flat.len() {
        return Err(harness(format!(
            "record has {} trailing values",
            flat.len() - s.at
        )));
    }
    Ok(s.out)
}

// ---- Rust side -----------------------------------------------------------------------

/// A Rust result: the flat record, or the unported callee that stopped it.
enum Rec {
    Flat(Vec<f64>),
    Deferred(SolidError),
}

/// `0` for success, else the status of the error; an unported callee defers.
fn st<T>(r: &solidcore::Result<T>) -> Result<f64, SolidError> {
    match r {
        Ok(_) => Ok(0.0),
        Err(e) if e.is_unported() => Err(e.clone()),
        Err(e) => Ok(status(e)),
    }
}

/// The B-spline of a curve of a kind derived from `SmBSplineCurve` (what the shim's
/// `put_curve` exports).
fn b_spline_of(c: Option<&Curve>) -> Option<&BSplineCurve> {
    match c?.kind() {
        CurveKind::BSpline(b) => Some(b),
        CurveKind::Line(l) => Some(l.base()),
        CurveKind::Ellipse(e) => Some(e.base()),
        CurveKind::Hyperbola(h) => Some(h.base()),
        CurveKind::Parabola(p) => Some(p.base()),
        _ => None,
    }
}

const fn trim_type(code: usize) -> BoundaryTrimmingType {
    match code {
        0 => BoundaryTrimmingType::None,
        1 => BoundaryTrimmingType::Minimal,
        2 => BoundaryTrimmingType::Maximal,
        3 => BoundaryTrimmingType::Bevel,
        _ => BoundaryTrimmingType::Blend,
    }
}

#[allow(clippy::too_many_lines)]
fn rust(i: &Inputs) -> Result<Rec, HarnessError> {
    let h = |e: SolidError| harness(e.to_string());
    let mut rm = RustModel::build(&i.model)?;
    let solvers: Vec<FilletSolver> = i
        .solvers
        .iter()
        .map(|s| FilletSolver::new(s.approx, s.ang, s.tang))
        .collect();
    for g in 0..i.model.geoms.len() {
        if let Ok(si) = u32::try_from(geom_solver(i, g)) {
            let gid = rm.geom(g)?;
            rm.fb.geom_mut(gid).map_err(h)?.solver = Some(FilletSolverId(si));
        }
    }
    let k = Kernel::default();
    let target = rm.target.clone();
    let mut o = Vec::new();
    macro_rules! status_of {
        ($r:expr) => {
            match st(&$r) {
                Ok(x) => x,
                Err(e) => return Ok(Rec::Deferred(e)),
            }
        };
    }
    let sides = |rm: &RustModel, o: &mut Vec<f64>| -> Result<(), HarnessError> {
        for g in 0..i.model.geoms.len() {
            let gid = rm.geom(g)?;
            o.push(rm.fb.geom(gid).map_err(h)?.side_eus.len() as f64);
        }
        Ok(())
    };
    match &i.query {
        Query::IntFillet { v } => {
            let vid = rm.vertex(*v)?;
            let r = rm.fb.calc_fillet_int_fillet(&k, vid, &solvers, &target);
            o.push(status_of!(r));
            let p = rm.fb.vertex_point(vid).map_err(h)?;
            o.extend([p.x, p.y, p.z]);
            o.push(rm.fb.fillet_vertex(vid).map_err(h)?.status() as i64 as f64);
        }
        Query::CrossSection { e } => {
            let eid = rm.edge(*e)?;
            let r = rm.fb.calc_cross_section(&k, eid, None, &solvers, &target);
            o.push(status_of!(r));
            let b = rm.fb.brep();
            dump_bcurve(&mut o, b_spline_of(edge::curve(b, eid).map_err(h)?));
            let prim = eu1::primary_edgeuse(b, eid)
                .map_err(h)?
                .ok_or_else(|| harness("no primary edge use"))?;
            let mate = eu1::mate(b, prim)
                .map_err(h)?
                .ok_or_else(|| harness("no mate"))?;
            for eu in [prim, mate] {
                let c = b.edgeuse(eu).map_err(h)?.uv_trim_curve.as_ref();
                dump_bcurve(&mut o, b_spline_of(c));
            }
            sides(&rm, &mut o)?;
        }
        Query::SideFace { e } => {
            let eid = rm.edge(*e)?;
            let r = rm
                .fb
                .calc_fillet_int_side_face(&k, eid, None, &solvers, &target);
            o.push(status_of!(r));
            o.push(rm.fb.fillet_edge(eid).map_err(h)?.status() as i64 as f64);
            sides(&rm, &mut o)?;
        }
        Query::ExtEdge { e } => {
            let eid = rm.edge(*e)?;
            let r = rm.fb.calc_edge_on_extended_edge(&k, eid, None, &target);
            o.push(status_of!(r));
        }
        Query::TrimRails { g } => {
            let gid = rm.geom(*g)?;
            let r = rm.fb.trim_rail_curves(&k, gid, &solvers, &target);
            o.push(status_of!(r));
            let b = rm.fb.brep();
            for ri in 0..2 {
                let Some(rail) = rm.fb.geom(gid).map_err(h)?.rail(ri) else {
                    o.push(0.0);
                    continue;
                };
                o.push(1.0);
                dump_bcurve(&mut o, b_spline_of(edge::curve(b, rail).map_err(h)?));
                let ivl = edge::interval(b, rail).map_err(h)?;
                o.extend([ivl.min(), ivl.max()]);
                let prim = eu1::primary_edgeuse(b, rail)
                    .map_err(h)?
                    .ok_or_else(|| harness("no primary edge use"))?;
                let mate = eu1::mate(b, prim)
                    .map_err(h)?
                    .ok_or_else(|| harness("no mate"))?;
                for eu in [prim, mate] {
                    let c = b.edgeuse(eu).map_err(h)?.uv_trim_curve.as_ref();
                    dump_bcurve(&mut o, b_spline_of(c));
                }
                for eu in [prim, mate] {
                    let vu = b
                        .edgeuse(eu)
                        .map_err(h)?
                        .vertexuse
                        .ok_or_else(|| harness("no vertex use"))?;
                    let uv = rm.fb.fillet_vertexuse(vu).map_err(h)?.uv_point();
                    o.extend([uv.x, uv.y]);
                }
            }
        }
        Query::SplitFace { g } => {
            let gid = rm.geom(*g)?;
            let r = rm.fb.update_geom_split_face(gid, None, None);
            o.push(status_of!(r));
        }
        Query::MakeEdge { vs, ve } => {
            let v = |x: i64| -> Result<_, HarnessError> {
                match usize::try_from(x) {
                    Ok(i) => Ok(Some(rm.vertex(i)?)),
                    Err(_) => Ok(None),
                }
            };
            let (a, b) = (v(*vs)?, v(*ve)?);
            let r = rm.fb.make_fillet_edge(a, b, None);
            o.push(status_of!(r));
            o.push(rm.fb.edge_list().len() as f64);
        }
        Query::RailEdge { g, rail } => {
            let gid = rm.geom(*g)?;
            let r = rm.fb.make_rail_edge(&k, gid, *rail);
            o.push(status_of!(r));
            o.push(rm.fb.edge_list().len() as f64);
        }
        Query::EndCurve { g, end, trim } => {
            let gid = rm.geom(*g)?;
            let r = rm
                .fb
                .create_end_curve(&k, gid, *end, trim_type(*trim), &solvers);
            o.push(status_of!(r));
            o.push(rm.fb.edge_list().len() as f64);
        }
    }
    Ok(Rec::Flat(o))
}

// ---- generation ----------------------------------------------------------------------

/// A vertex at `p` with no classification (the model's defaults), `status` its fillet
/// status (`1` processed).
const fn vertex_at(p: [f64; 3], status: i64) -> VertexIn {
    VertexIn {
        point: p,
        vtype: 0,
        status,
        corner: -1,
        pc_class: 100,
        pc_kind: 0,
        pc_index: -1,
        pc_sub: 0,
        pc_face: -1,
        t_param: 0.0,
        uv: [0.0, 0.0],
        mates: [-1, -1],
    }
}

fn uv_line(a: [f64; 2], b: [f64; 2], ivl: [f64; 2]) -> Crv {
    Crv {
        dim: 2,
        cpts: vec![[a[0], a[1], 0.0, NOW], [b[0], b[1], 0.0, NOW]],
        p: 1,
        u: vec![ivl[0], ivl[0], ivl[1], ivl[1]],
    }
}

/// The iso-v boundary curve of `s` at its first (`last = false`) or last v row: the
/// curve of the row's control points over the u knots (the surface's own
/// parameterization).
fn row_curve(s: &Srf, last: bool) -> Crv {
    let cpts = s
        .net
        .iter()
        .filter_map(|row| if last { row.last() } else { row.first() })
        .copied()
        .collect();
    Crv {
        dim: 3,
        cpts,
        p: s.p,
        u: s.u.clone(),
    }
}

fn eval(c: &Crv, t: f64) -> Option<[f64; 3]> {
    let b = curve_of(c).ok()?;
    let p = Curve::from(b).evaluate_point(t).ok()?;
    Some([p.x, p.y, p.z])
}

fn surface_point(s: &Srf, u: f64, v: f64) -> Option<[f64; 3]> {
    let srf = surface_of(s).ok()?;
    let d = srf.natural_uv_domain();
    let uv = solidcore::math::Point2::new(
        d.min().x + u * (d.max().x - d.min().x),
        d.min().y + v * (d.max().y - d.min().y),
    );
    let p = srf.evaluate_point(uv).ok()?;
    Some([p.x, p.y, p.z])
}

fn solvers_for(rng: &mut Rng) -> Vec<SolverIn> {
    vec![SolverIn {
        approx: *rng.choose(&[1e-6, 1e-3, 1e-9, 1e-5]).unwrap_or(&1e-6),
        ang: rng.uniform(0.01, 0.5),
        tang: rng.uniform(0.001, 0.1),
    }]
}

/// Two rails on the v boundaries of a random fillet surface with their end vertices
/// (processed or not, on the rail or a little off it, sometimes coincident), the rails'
/// mate uses sometimes carrying a UV line; a default or cliff-rollover geometry.
#[allow(clippy::too_many_lines)]
fn gen_trim(rng: &mut Rng) -> Option<Inputs> {
    let scale = scale_of(rng);
    let srf = gen_surface(rng, scale);
    let mut rails = [row_curve(&srf, false), row_curve(&srf, true)];
    if rng.chance(0.15) {
        rails.swap(0, 1);
    }
    let mut vertices = Vec::new();
    let mut edges = Vec::new();
    for (r, crv) in rails.iter().enumerate() {
        let (a, b) = (*crv.u.first()?, *crv.u.last()?);
        let mut t0 = a + (b - a) * rng.uniform(0.0, 0.45);
        let mut t1 = a + (b - a) * rng.uniform(0.55, 1.0);
        if rng.chance(0.1) {
            std::mem::swap(&mut t0, &mut t1);
        }
        if rng.chance(0.1) {
            t1 = t0;
        }
        let noise = scale * *rng.choose(&[0.0, 0.0, 1e-7, 1e-5, 1e-3]).unwrap_or(&0.0);
        let mut p0 = eval(crv, t0)?;
        let mut p1 = eval(crv, t1)?;
        for c in p0.iter_mut().chain(p1.iter_mut()) {
            *c += rng.uniform(-noise, noise);
        }
        let s0 = i64::from(rng.chance(0.85));
        let s1 = i64::from(rng.chance(0.85));
        let base = vertices.len();
        vertices.push(vertex_at(p0, s0));
        vertices.push(vertex_at(p1, s1));
        let mate_uv = if rng.chance(0.4) {
            let y = rng.uniform(0.0, 1.0);
            Some(uv_line([a, y], [b, y + rng.uniform(-0.1, 0.1)], [a, b]))
        } else {
            None
        };
        edges.push(EdgeIn {
            vs: base,
            ve: base + 1,
            etype: if rng.chance(0.95) { 9 } else { 2 },
            status: 0,
            corner: -1,
            curve: Some(crv.clone()),
            interval: Some([a, b]),
            fbe: [-1, -1],
            eus: [
                EuIn {
                    uv_curve: None,
                    geom: 0,
                    vu_uv: [0.0, 0.0],
                },
                EuIn {
                    uv_curve: mate_uv,
                    geom: -1,
                    vu_uv: [0.0, 0.0],
                },
            ],
        });
        let _ = r;
    }
    let geom = GeomIn {
        gtype: if rng.chance(0.3) { 2 } else { 0 },
        status: 0,
        center: None,
        surface: Some(srf),
        rails: [0, 1],
        other1: Vec::new(),
        other2: Vec::new(),
        sides: Vec::new(),
    };
    Some(Inputs {
        model: Model {
            targets: Vec::new(),
            faces: Vec::new(),
            fillet_brep: Vec::new(),
            vertices,
            geoms: vec![geom],
            edges,
        },
        solvers: solvers_for(rng),
        geom_solvers: vec![if rng.chance(0.95) { 0 } else { -1 }],
        query: Query::TrimRails { g: 0 },
    })
}

/// A cross-section edge between two points of a random fillet surface (or of a second
/// geometry's surface on the mate use).
fn gen_cross(rng: &mut Rng) -> Option<Inputs> {
    let scale = scale_of(rng);
    let s1 = gen_surface(rng, scale);
    let u = rng.uniform(0.0, 1.0);
    let (v0, v1) = if rng.chance(0.5) {
        (0.0, 1.0)
    } else {
        (rng.uniform(0.0, 0.5), rng.uniform(0.5, 1.0))
    };
    let u1 = if rng.chance(0.6) {
        u
    } else {
        rng.uniform(0.0, 1.0)
    };
    let p0 = surface_point(&s1, u, v0)?;
    let p1 = surface_point(&s1, u1, v1)?;
    let two = rng.chance(0.3);
    let mut geoms = vec![GeomIn {
        gtype: 0,
        status: 0,
        center: None,
        surface: if rng.chance(0.95) {
            Some(s1.clone())
        } else {
            None
        },
        rails: [-1, -1],
        other1: Vec::new(),
        other2: Vec::new(),
        sides: Vec::new(),
    }];
    if two {
        geoms.push(GeomIn {
            surface: Some(s1),
            ..geoms.first()?.clone()
        });
    }
    let edges = vec![EdgeIn {
        vs: 0,
        ve: 1,
        etype: 1,
        status: 0,
        corner: -1,
        curve: None,
        interval: None,
        fbe: [-1, -1],
        eus: [
            EuIn {
                uv_curve: None,
                geom: 0,
                vu_uv: [0.0, 0.0],
            },
            EuIn {
                uv_curve: None,
                geom: if two { 1 } else { -1 },
                vu_uv: [0.0, 0.0],
            },
        ],
    }];
    let ng = geoms.len();
    Some(Inputs {
        model: Model {
            targets: Vec::new(),
            faces: Vec::new(),
            fillet_brep: Vec::new(),
            vertices: vec![vertex_at(p0, 1), vertex_at(p1, 1)],
            geoms,
            edges,
        },
        solvers: solvers_for(rng),
        geom_solvers: vec![0; ng],
        query: if rng.chance(0.8) {
            Query::CrossSection { e: 0 }
        } else {
            Query::SideFace { e: 0 }
        },
    })
}

/// A query on p1's random model, with base solvers on some geometries.
fn gen_topology(rng: &mut Rng) -> Inputs {
    let mut model = gen_model(rng);
    let nv = model.vertices.len().max(1);
    let ng = model.geoms.len().max(1);
    let ne = model.edges.len().max(1);
    let ns = rng.below(3) as usize;
    let solvers = (0..ns)
        .map(|_| solvers_for(rng).remove(0))
        .collect::<Vec<_>>();
    let geom_solvers = (0..model.geoms.len())
        .map(|_| {
            if ns > 0 && rng.chance(0.7) {
                rng.below(ns as u64) as i64
            } else {
                -1
            }
        })
        .collect::<Vec<_>>();
    for v in &mut model.vertices {
        if rng.chance(0.7) {
            v.corner = -1;
        }
    }
    let g = rng.below(ng as u64) as usize;
    let e = rng.below(ne as u64) as usize;
    let query = match rng.below(9) {
        0 => Query::IntFillet {
            v: rng.below(nv as u64) as usize,
        },
        1 => Query::CrossSection { e },
        2 => Query::SideFace { e },
        3 => Query::ExtEdge { e },
        4 => Query::TrimRails { g },
        5 => Query::SplitFace { g },
        6 => Query::MakeEdge {
            vs: rng.below(nv as u64) as i64 - i64::from(rng.chance(0.2)),
            ve: rng.below(nv as u64) as i64,
        },
        7 => Query::RailEdge {
            g,
            rail: rng.below(2) as usize,
        },
        _ => Query::EndCurve {
            g,
            end: rng.below(2) as usize,
            trim: rng.below(4) as usize,
        },
    };
    Inputs {
        model,
        solvers,
        geom_solvers,
        query,
    }
}

impl UnitDiff for Runner {
    fn id(&self) -> &'static str {
        "blend.fillet_geom_p5"
    }

    fn run_cpp(&self, case: &Case) -> Result<Output, HarnessError> {
        let c: CaseInputs = case.parse()?;
        cpp(&c.inputs)
    }

    fn run_rust(&self, case: &Case) -> Result<Output, HarnessError> {
        let i = case.parse::<CaseInputs>()?.inputs;
        // an input the FFI refuses is a harness error on both sides
        ffi::validate(i.query.op(), &cpp_input(&i)?).map_err(|e| harness(e.to_string()))?;
        match rust(&i)? {
            Rec::Flat(flat) => to_output(&i, &flat),
            Rec::Deferred(e) => {
                Ok(Output::unported(&e).unwrap_or_else(|| Output::err(e.to_string())))
            }
        }
    }

    fn generate(&self, rng: &mut Rng) -> Option<Value> {
        // keep only inputs the shim runs (it refuses those upstream would dereference a
        // corner, solver or rail the harness cannot provide for)
        for _ in 0..64 {
            let inputs = match rng.below(20) {
                0..=10 => match gen_trim(rng) {
                    Some(i) => i,
                    None => continue,
                },
                11..=15 => match gen_cross(rng) {
                    Some(i) => i,
                    None => continue,
                },
                _ => gen_topology(rng),
            };
            let Ok(flat) = cpp_input(&inputs) else {
                continue;
            };
            if ffi::eval(inputs.query.op(), &flat).is_ok() {
                return serde_json::to_value(inputs).ok();
            }
        }
        None
    }

    fn cases(&self, ctx: &CaseContext) -> Result<Vec<Case>, HarnessError> {
        let ctx = CaseContext {
            seeds: ctx.seeds.max(MIN_GENERATED),
            ..ctx.clone()
        };
        default_cases(self, &ctx)
    }
}
