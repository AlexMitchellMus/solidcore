//! Runner for unit `blend.fillet_corner_p1`: the file-static helpers of
//! `SmFilletCorner.cpp` (the corner methods wait for the `FilletCorner` record of #256).
//!
//! C++: `solidcore_ffi::units::blend_fillet_corner_p1` (the shim compiles
//! `source/SMLib/src/SmFilletCorner.cpp`). Rust: `solidcore::blend::fillet_corner::p1`.
//!
//! Two kinds of case:
//!
//! - **Geometry**: `planar` (`sm_CreatePlanarSurface`: 2 to 5 points, collinear,
//!   tiny and huge scales), `revolution` (`sm_CreateSrfOfRevolution`: boundary loops of a
//!   real surface of revolution built from exact rational arcs and a generating curve in a
//!   plane through the axis, with 3 or 4 curves, a pole, a perturbed opposite curve, an
//!   arc of a different span, a major curve that is no arc), `fillet_radius`
//!   (`sm_GetFilletRadius`), `seam_out` and `closure` (`sm_CheckSeamOut`,
//!   `sm_SetSurfaceClosure`, on planes, open and closed B-spline surfaces). Compared: the
//!   status and the whole surface block (the plane frame, scale and domain; the B-spline
//!   net, knots and form), or the values.
//! - **Topology**: `g1`, `face_by_two`, `tangential`, `concavity` on a real B-rep (unit
//!   `topology.edgeuse_p1`'s B-rep block: the box-like fillet test files of
//!   `TestFiles/pt_TestFiles/Fillet`, which have convex and concave corners, and that unit's
//!   generated B-reps). The C++ side builds the B-rep, the Rust side imports the shim's
//!   export. A vertex is selected; edges and edge uses are selected among those at it (the
//!   edges ending at it, the edge uses starting at it), so the inputs are real corners.
//!   Compared: the flags, the angle and the face (export index).
//!
//! The shim refuses (`-2`) the inputs where upstream dereferences null; the Rust side
//! returns [`SolidError::InvalidInput`] there and both report `refused` (a both-errored
//! case). A topology case whose evaluated geometry the import holds only as its B-spline
//! part (analytic surfaces and curves) is deferred as `topology.edgeuse_p2`'s runner
//! defers it.

use serde::{Deserialize, Serialize};
use serde_json::{json, Value};

use solidcore::blend::fillet_corner::p1 as corner;
use solidcore::curve::curve::Curve;
use solidcore::geom::extent2d::Extent2d;
use solidcore::geom::vector2d::Point2;
use solidcore::kernel::Kernel;
use solidcore::math::Vec3;
use solidcore::nurbs::b_spline_surface::p1::SurfParamType;
use solidcore::nurbs::math::NOW;
use solidcore::surface::surface::{Surface, SurfaceKind};
use solidcore::topology::edgeuse::p1 as eu1;
use solidcore::topology::{edge, Brep, EdgeId, EdgeuseId, VertexId};
use solidcore::SolidError;
use solidcore_ffi::units::blend_fillet_corner_p1::{self as ffi, Op};
use solidcore_ffi::units::topology_edgeuse_p1 as ffi1;

use super::topology_edgeuse_p1::{
    brep_block_of, bspline_of, domain, flag, gen_brep, gen_curve, gen_surface, harness,
    import_with_types, path_of, push_curve, push_surface, sel, surface_of, BrepIn, CrvIn,
    ImportTypes, SrfIn,
};
use crate::{Case, HarnessError, Output, Rng, UnitDiff};

/// Registered runner.
pub struct Runner;

/// The reason [`corner`] gives where upstream dereferences null.
pub(super) const NULL_DEREF: &str = "upstream dereferences null";

/// The box-like fillet test files (planar faces with convex and concave corners, one
/// all-NURBS box, one with a cylinder; filleted and B-spline parts whose tangential
/// edges reach the success paths of the tangency tests). Not `SmallBox.smb`: `SmBrep::ReadFromFile`
/// crashes on it (UB-1).
pub(super) const CORNER_FILES: [&str; 13] = [
    "Fillet/box1x1x2.smb",
    "Fillet/brick.smb",
    "Fillet/corner_box.smb",
    "Fillet/l_shaped_box.smb",
    "Fillet/6_sided_corner.smb",
    "Fillet/box_w_slot.smb",
    "Fillet/nurb_box.smb",
    "Fillet/box_on_box.smb",
    "Fillet/cyl_on_box.smb",
    "Fillet/filleted_box.smb",
    "Fillet/nurb_nurb_1.smb",
    "Fillet/nurb_plane_w_notch1.smb",
    "Fillet/brepForFillet1.smb",
];

/// One case.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "op", rename_all = "snake_case")]
enum Inputs {
    Planar {
        points: Vec<[f64; 3]>,
        tol: f64,
    },
    Revolution {
        curves: Vec<CrvIn>,
        major: usize,
        tol: f64,
    },
    FilletRadius {
        p1: [f64; 3],
        n1: [f64; 3],
        p2: [f64; 3],
        n2: [f64; 3],
        orient: f64,
    },
    SeamOut {
        surface: SrfIn,
        uv: [f64; 2],
        tol: f64,
        closure: i64,
    },
    Closure {
        surface: SrfIn,
        dom: [f64; 4],
    },
    G1 {
        brep: BrepIn,
        v: i64,
        e1: i64,
        e2: i64,
        tol: f64,
    },
    FaceByTwo {
        brep: BrepIn,
        v: i64,
        eu1: i64,
        eu2: i64,
        eu3: i64,
    },
    Tangential {
        brep: BrepIn,
        eu: i64,
        tol: f64,
    },
    Concavity {
        brep: BrepIn,
        v: i64,
        eu1: i64,
        eu2: i64,
        #[serde(rename = "ref")]
        ref_eu: i64,
        /// `eu2` is `eu1`'s corner mate at `v` (taken outward), as the edge uses around
        /// a corner are paired upstream; else `eu2` selects among those at `v`.
        #[serde(default)]
        adjacent: bool,
    },
}

// ---- shared -----------------------------------------------------------------------------

fn v3(a: [f64; 3]) -> Vec3 {
    Vec3::from_array(a)
}

/// The `SmStatus` of a Rust error: the upstream code, `SM_ERR_INVALID_INPUT` (1007) for
/// [`SolidError::InvalidInput`]; `None` for a null dereference (refused) or an error
/// with no upstream status.
fn status_of(e: &SolidError) -> Option<f64> {
    match e {
        SolidError::Upstream { code, .. } => Some(f64::from(*code)),
        SolidError::InvalidInput { reason, .. } if *reason != NULL_DEREF => Some(1007.0),
        _ => None,
    }
}

/// A Rust error as an output: deferred for `Unported`, `refused` for a null
/// dereference, else the status.
pub(super) fn rust_err(e: &SolidError) -> Output {
    if let Some(o) = Output::unported(e) {
        return o;
    }
    if matches!(e, SolidError::InvalidInput { reason, .. } if *reason == NULL_DEREF) {
        return Output::err("refused");
    }
    match status_of(e) {
        Some(st) => Output::new().values("status", vec![st]),
        None => Output::err(format!("rust: {e}")),
    }
}

/// A surface block, as the shim's `put_surface` (planes and B-spline surfaces).
fn put_surface(o: &mut Vec<f64>, s: &Surface) -> Result<(), HarnessError> {
    match s.kind() {
        SurfaceKind::Plane(p) => {
            o.push(1.0);
            let pos = p.position();
            for v in [pos.origin(), pos.x_axis(), pos.y_axis()] {
                o.extend([v.x, v.y, v.z]);
            }
            let sc = p.uv_scale();
            o.extend([sc.x, sc.y]);
            let d = p.step_uv_domain();
            o.extend([d.min().x, d.min().y, d.max().x, d.max().y]);
        }
        SurfaceKind::BSpline(bs) => {
            let nurb = bs
                .gw_nurb()
                .ok_or_else(|| harness("surface without NURBS data"))?;
            let (p, q) = nurb.degrees();
            o.extend([
                0.0,
                flag(bs.is_out_of_bounds_enabled()),
                bs.form() as u64 as f64,
                nurb.n() as f64,
                nurb.m() as f64,
                p as f64,
                q as f64,
                (nurb.u_knots().len() - 1) as f64,
                (nurb.v_knots().len() - 1) as f64,
            ]);
            for row in nurb.cpts() {
                for pt in row {
                    o.extend([pt.x, pt.y, pt.z, pt.w]);
                }
            }
            o.extend_from_slice(nurb.u_knots());
            o.extend_from_slice(nurb.v_knots());
        }
        _ => o.push(-2.0),
    }
    Ok(())
}

/// A C++ `status [surface block]` record as an output.
fn cpp_surface(out: &[f64]) -> Result<Output, HarnessError> {
    let st = *out.first().ok_or_else(|| harness("empty output"))?;
    let o = Output::new().values("status", vec![st]);
    Ok(if st == 0.0 {
        o.values("surface", out.get(1..).unwrap_or_default().to_vec())
    } else {
        o
    })
}

/// A Rust surface result as an output.
fn rust_surface(r: solidcore::Result<Surface>) -> Result<Output, HarnessError> {
    match r {
        Ok(s) => {
            let mut block = Vec::new();
            put_surface(&mut block, &s)?;
            Ok(Output::new()
                .values("status", vec![0.0])
                .values("surface", block))
        }
        Err(e) => Ok(rust_err(&e)),
    }
}

fn surface_args(s: &SrfIn) -> Result<Vec<f64>, HarnessError> {
    let mut v = Vec::new();
    push_surface(&mut v, s)?;
    Ok(v)
}

fn closure_of(code: i64) -> Result<SurfParamType, HarnessError> {
    Ok(match code {
        0 => SurfParamType::U,
        1 => SurfParamType::V,
        2 => SurfParamType::Unknown,
        3 => SurfParamType::Both,
        4 => SurfParamType::Neither,
        _ => return Err(harness(format!("closure code {code} outside 0..=4"))),
    })
}

// ---- topology setup ---------------------------------------------------------------------

/// The built B-rep of a topology case on both sides and the resolved graph indices of the
/// op's arguments (`-1` = NULL).
pub(super) struct Topo {
    pub(super) path: String,
    pub(super) block: Vec<f64>,
    pub(super) brep: Brep,
    pub(super) types: ImportTypes,
    pub(super) args: Vec<i64>,
}

pub(super) fn eu_vertex(b: &Brep, eu: EdgeuseId) -> Option<VertexId> {
    let vu = b.edgeuse(eu).ok()?.vertexuse?;
    b.vertexuse(vu).ok()?.vertex
}

fn idx_u32(i: i64) -> Result<u32, HarnessError> {
    u32::try_from(i).map_err(|_| harness("index"))
}

/// A selector into `ring` (export indices), falling back to all `n` when it is empty.
fn pick(s: i64, ring: &[i64], n: usize) -> i64 {
    if s < 0 {
        return -1;
    }
    let r = sel(s, ring.len());
    if r >= 0 {
        return ring.get(r as usize).copied().unwrap_or(-1);
    }
    sel(s, n)
}

/// Builds and imports the B-rep and resolves the selectors of `i`; `Err(output)` when the
/// case ends before the op (the B-rep cannot be built, a required entity is missing),
/// with the same output on both sides.
fn topo(i: &Inputs) -> Result<Result<Topo, Output>, HarnessError> {
    let (Inputs::G1 { brep, .. }
    | Inputs::FaceByTwo { brep, .. }
    | Inputs::Tangential { brep, .. }
    | Inputs::Concavity { brep, .. }) = i
    else {
        return Err(harness("not a topology case"));
    };
    let path = path_of(Some(brep))?;
    let block = brep_block_of(Some(brep), false, false)?;
    let ex = ffi1::eval(ffi1::Op::Export, &path, &block)
        .map_err(|e| harness(format!("c++ export: {e}")))?;
    let Some(c) = ffi1::counts(&ex) else {
        let st = ex.get(1).copied().unwrap_or(0.0);
        return Ok(Err(Output::err(format!("brep_build_status_{st}"))));
    };
    let (n_eu, n_e, n_v) = (c[6], c[7], c[9]);
    let (b, types) = import_with_types(&ex)?;

    // the edge uses starting at `v`, the edges ending at `v` (export order)
    let ring_eus = |v: i64| -> Vec<i64> {
        if v < 0 {
            return Vec::new();
        }
        let vid = VertexId::from_index(v as u32);
        b.edgeuses()
            .iter()
            .map(|(id, _)| id)
            .filter(|&id| eu_vertex(&b, id) == Some(vid))
            .map(|id| i64::from(id.index()))
            .collect()
    };
    let ring_edges = |v: i64| -> Vec<i64> {
        if v < 0 {
            return Vec::new();
        }
        let vid = Some(VertexId::from_index(v as u32));
        b.edges()
            .iter()
            .map(|(id, _)| id)
            .filter(|&id| {
                edge::start_vertex(&b, id).ok().flatten() == vid
                    || edge::end_vertex(&b, id).ok().flatten() == vid
            })
            .map(|id| i64::from(id.index()))
            .collect()
    };

    let args = match i {
        Inputs::G1 { v, e1, e2, .. } => {
            let vv = sel(*v, n_v);
            let ring = ring_edges(vv);
            vec![vv, pick(*e1, &ring, n_e), pick(*e2, &ring, n_e)]
        }
        Inputs::FaceByTwo {
            v, eu1, eu2, eu3, ..
        } => {
            let ring = ring_eus(sel(*v, n_v));
            vec![
                pick(*eu1, &ring, n_eu),
                pick(*eu2, &ring, n_eu),
                sel(*eu3, n_eu),
            ]
        }
        Inputs::Tangential { eu, .. } => {
            let e = sel(*eu, n_eu);
            if e < 0 {
                return Ok(Err(Output::err("no_edge_use")));
            }
            vec![e]
        }
        Inputs::Concavity {
            v,
            eu1,
            eu2,
            ref_eu,
            adjacent,
            ..
        } => {
            let vv = sel(*v, n_v);
            let ring = ring_eus(vv);
            if vv < 0 || ring.is_empty() {
                return Ok(Err(Output::err("no_edge_use")));
            }
            let a = pick(eu1.rem_euclid(1 << 30), &ring, n_eu);
            let mut bb = pick(eu2.rem_euclid(1 << 30), &ring, n_eu);
            if *adjacent {
                // the corner mate of `a` at the vertex, outward (itself or its mate)
                let vid = VertexId::from_index(idx_u32(vv)?);
                let eu_a = EdgeuseId::from_index(idx_u32(a)?);
                if let Ok(Some(cm)) = eu1::corner_mate_edgeuse(&b, eu_a, Some(vid)) {
                    let out = if eu_vertex(&b, cm) == Some(vid) {
                        Some(cm)
                    } else {
                        eu1::mate(&b, cm).ok().flatten()
                    };
                    if let Some(o) = out {
                        bb = i64::from(o.index());
                    }
                }
            }
            vec![vv, a, bb, sel(*ref_eu, n_eu)]
        }
        _ => return Err(harness("not a topology case")),
    };
    Ok(Ok(Topo {
        path,
        block,
        brep: b,
        types,
        args,
    }))
}

/// `SmBSplineSurface_TYPE`, `SmPlane_TYPE`; `SmBSplineCurve_TYPE`, `SmLine_TYPE`,
/// `SmCrvOnSurf_TYPE` (as `topology.edgeuse_p2`'s runner).
const EXACT_SURFACES: [i64; 2] = [15002, 15011];
const EXACT_CURVES: [i64; 3] = [14033, 14010, 14048];

/// The deferral label of cases whose geometry the shared B-rep export holds only as a
/// B-spline part (an analytic face surface or edge curve the import cannot rebuild
/// exactly): a harness limitation, not an unported unit (`topology.brep_p1`'s
/// convention).
const HARNESS_ANALYTIC_EXPORT: &str = "harness:analytic-geometry-export";

/// The deferral of a case that evaluates geometry the import holds only as its B-spline
/// part: the labels of `topology.edgeuse_p2`'s `inexact_geometry`.
pub(super) fn inexact(t: &Topo, edges: &[EdgeId], eus: &[EdgeuseId]) -> Option<Output> {
    let b = &t.brep;
    let mut all_edges = edges.to_vec();
    for &eu in eus {
        let mut both = vec![eu];
        if let Ok(Some(r)) = eu1::radial(b, eu) {
            both.push(r);
        }
        for e in both {
            let st = eu1::face(b, e)
                .ok()
                .flatten()
                .and_then(|f| t.types.face_surface_types.get(f.index() as usize).copied())
                .unwrap_or(-1);
            if st >= 0 && !EXACT_SURFACES.contains(&st) {
                let (unit, symbol) = match st {
                    15012 => ("surface.surface_p2", "SurfaceKind::Cone (SmCone dispatch)"),
                    15013 => (
                        "surface.surface_p2",
                        "SurfaceKind::Cone (SmCylinder dispatch)",
                    ),
                    15014 => (
                        "surface.surface_p2",
                        "SurfaceKind::Sphere (SmSphere dispatch)",
                    ),
                    15015 => (
                        "surface.surface_p2",
                        "SurfaceKind::Torus (SmTorus dispatch)",
                    ),
                    15030 => (
                        "surface.surf_of_extrusion_p2",
                        "SurfaceKind::SurfOfExtrusion",
                    ),
                    _ => (
                        HARNESS_ANALYTIC_EXPORT,
                        "face surface of an analytic kind exported as its B-spline part",
                    ),
                };
                return Some(Output::deferred(unit, symbol));
            }
            if let Ok(Some(ed)) = eu1::edge(b, e) {
                all_edges.push(ed);
            }
        }
    }
    for e in all_edges {
        let ct = t
            .types
            .edge_curve_types
            .get(e.index() as usize)
            .copied()
            .unwrap_or(-1);
        if ct >= 0 && !EXACT_CURVES.contains(&ct) {
            return Some(Output::deferred(
                HARNESS_ANALYTIC_EXPORT,
                "edge curve of an analytic kind exported as its B-spline part",
            ));
        }
    }
    None
}

fn opt<T>(i: i64, f: impl Fn(u32) -> T) -> Result<Option<T>, HarnessError> {
    Ok(if i < 0 { None } else { Some(f(idx_u32(i)?)) })
}

// ---- C++ side ---------------------------------------------------------------------------

fn run_cpp(i: &Inputs) -> Result<Output, HarnessError> {
    let call = |op: Op, path: &str, input: &[f64]| {
        ffi::eval(op, path, input).map_err(|e| harness(format!("c++: {e}")))
    };
    match i {
        Inputs::Planar { points, tol } => {
            let mut v = vec![points.len() as f64];
            for p in points {
                v.extend(p);
            }
            v.push(*tol);
            cpp_surface(&call(Op::Planar, "", &v)?)
        }
        Inputs::Revolution { curves, major, tol } => {
            let mut v = vec![curves.len() as f64];
            for c in curves {
                push_curve(&mut v, c)?;
            }
            v.extend([*major as f64, *tol]);
            cpp_surface(&call(Op::Revolution, "", &v)?)
        }
        Inputs::FilletRadius {
            p1,
            n1,
            p2,
            n2,
            orient,
        } => {
            let mut v = Vec::new();
            for a in [p1, n1, p2, n2] {
                v.extend(a);
            }
            v.push(*orient);
            Ok(Output::new().values("radius_span", call(Op::FilletRadius, "", &v)?))
        }
        Inputs::SeamOut {
            surface,
            uv,
            tol,
            closure,
        } => {
            let mut v = surface_args(surface)?;
            v.extend([uv[0], uv[1], *tol, *closure as f64]);
            Ok(Output::new().values("seam", call(Op::SeamOut, "", &v)?))
        }
        Inputs::Closure { surface, dom } => {
            let mut v = surface_args(surface)?;
            v.extend(dom);
            Ok(Output::new().values("closure", call(Op::Closure, "", &v)?))
        }
        _ => {
            let t = match topo(i)? {
                Ok(t) => t,
                Err(o) => return Ok(o),
            };
            let (op, tail): (Op, Vec<f64>) = match i {
                Inputs::G1 { tol, .. } => (Op::G1, vec![*tol]),
                Inputs::FaceByTwo { .. } => (Op::FaceByTwo, vec![]),
                Inputs::Tangential { tol, .. } => (Op::Tangential, vec![*tol]),
                _ => (Op::Concavity, vec![]),
            };
            let mut v = t.block.clone();
            v.extend(t.args.iter().map(|&a| a as f64));
            v.extend(tail);
            let out = call(op, &t.path, &v)?;
            let Some((&head, rest)) = out.split_first() else {
                return Err(harness("empty topology output"));
            };
            if head != 0.0 {
                let st = rest.first().copied().unwrap_or(0.0);
                return Ok(Output::err(format!("brep_build_status_{st}")));
            }
            if rest.first().copied() == Some(-2.0) {
                return Ok(Output::err("refused"));
            }
            let key = match op {
                Op::G1 => "g1",
                Op::FaceByTwo => "face",
                _ => "flag",
            };
            Ok(Output::new().values(key, rest.to_vec()))
        }
    }
}

// ---- Rust side --------------------------------------------------------------------------

fn run_rust(i: &Inputs) -> Result<Output, HarnessError> {
    let k = Kernel::default();
    match i {
        Inputs::Planar { points, tol } => {
            let pts: Vec<Vec3> = points.iter().map(|&p| v3(p)).collect();
            rust_surface(corner::create_planar_surface(&k, &pts, *tol))
        }
        Inputs::Revolution { curves, major, tol } => {
            let cs = curves
                .iter()
                .map(|c| bspline_of(c).map(Curve::from))
                .collect::<Result<Vec<_>, _>>()?;
            rust_surface(corner::create_srf_of_revolution(&k, &cs, *major, *tol))
        }
        Inputs::FilletRadius {
            p1,
            n1,
            p2,
            n2,
            orient,
        } => match corner::fillet_radius(v3(*p1), v3(*n1), v3(*p2), v3(*n2), *orient) {
            Ok((r, s)) => Ok(Output::new().values("radius_span", vec![r, s])),
            Err(e) => Ok(rust_err(&e)),
        },
        Inputs::SeamOut {
            surface,
            uv,
            tol,
            closure,
        } => {
            let s = surface_of(surface)?;
            let mut p = Point2::new(uv[0], uv[1]);
            let ret = corner::check_seam_out(&s, &mut p, *tol, closure_of(*closure)?);
            Ok(Output::new().values("seam", vec![flag(ret), p.x, p.y]))
        }
        Inputs::Closure { surface, dom } => {
            let s = surface_of(surface)?;
            let d = Extent2d::from_points(Point2::new(dom[0], dom[1]), Point2::new(dom[2], dom[3]));
            match corner::set_surface_closure(&s, &d) {
                Ok(c) => Ok(Output::new().values("closure", vec![c as u64 as f64])),
                Err(e) => Ok(rust_err(&e)),
            }
        }
        _ => {
            let mut t = match topo(i)? {
                Ok(t) => t,
                Err(o) => return Ok(o),
            };
            let args = t.args.clone();
            let a = |n: usize| -> Result<i64, HarnessError> {
                args.get(n)
                    .copied()
                    .ok_or_else(|| harness("missing argument"))
            };
            let eu = |x: i64| opt(x, EdgeuseId::from_index);
            let r: solidcore::Result<Vec<f64>> = match i {
                Inputs::G1 { tol, .. } => {
                    let (e1, e2) = (
                        opt(a(1)?, EdgeId::from_index)?,
                        opt(a(2)?, EdgeId::from_index)?,
                    );
                    let edges: Vec<EdgeId> = [e1, e2].into_iter().flatten().collect();
                    if let Some(o) = inexact(&t, &edges, &[]) {
                        return Ok(o);
                    }
                    let v = opt(a(0)?, VertexId::from_index)?;
                    corner::test_g1_continuity(&t.brep, v, e1, e2, *tol)
                        .map(|(f, ang)| vec![flag(f), flag(ang.is_some()), ang.unwrap_or(0.0)])
                }
                Inputs::FaceByTwo { .. } => {
                    corner::face_by_two_edgeuses(&t.brep, eu(a(0)?)?, eu(a(1)?)?, eu(a(2)?)?)
                        .map(|f| vec![f.map_or(-1.0, |f| f64::from(f.index()))])
                }
                Inputs::Tangential { tol, .. } => {
                    let e = EdgeuseId::from_index(idx_u32(a(0)?)?);
                    if let Some(o) = inexact(&t, &[], &[e]) {
                        return Ok(o);
                    }
                    corner::check_tangential_edgeuse(&k, &mut t.brep, e, *tol)
                        .map(|f| vec![flag(f)])
                }
                _ => {
                    let v = VertexId::from_index(idx_u32(a(0)?)?);
                    let (e1, e2) = (
                        EdgeuseId::from_index(idx_u32(a(1)?)?),
                        EdgeuseId::from_index(idx_u32(a(2)?)?),
                    );
                    let r = eu(a(3)?)?;
                    let eus: Vec<EdgeuseId> =
                        [Some(e1), Some(e2), r].into_iter().flatten().collect();
                    if let Some(o) = inexact(&t, &[], &eus) {
                        return Ok(o);
                    }
                    corner::check_face_complimentary_corner_concavity(&k, &mut t.brep, v, e1, e2, r)
                        .map(|f| vec![flag(f)])
                }
            };
            let key = match i {
                Inputs::G1 { .. } => "g1",
                Inputs::FaceByTwo { .. } => "face",
                _ => "flag",
            };
            match r {
                Ok(v) => Ok(Output::new().values(key, v)),
                Err(e) => Ok(rust_err(&e)),
            }
        }
    }
}

// ---- generation -------------------------------------------------------------------------

fn scale_of(rng: &mut Rng) -> f64 {
    *rng.choose(&[1e-4, 1.0, 1.0, 1.0, 1.0, 1e4]).unwrap_or(&1.0)
}

fn gen_planar(rng: &mut Rng) -> Inputs {
    let s = scale_of(rng);
    let n = *rng.choose(&[2usize, 3, 3, 3, 4, 4, 4, 5]).unwrap_or(&3);
    let mut points: Vec<[f64; 3]> = (0..n)
        .map(|_| rng.point3(-1.0, 1.0).map(|c| c * s))
        .collect();
    let collinear = rng.chance(0.15);
    let t = rng.uniform(-2.0, 2.0);
    if let [a, b, c] = points.as_mut_slice() {
        if collinear {
            // collinear (or coincident) points: Unitize fails
            *c = (v3(*a) + (v3(*b) - v3(*a)) * t).to_array();
        }
    }
    Inputs::Planar {
        points,
        tol: rng.uniform(0.0, 1e-3) * s,
    }
}

/// The exact rational quadratic arc of `angle` (radians, at most `PI / 2`) about the axis
/// through `c` along unit `z`, starting at `p` (homogeneous control points).
fn arc(c: Vec3, z: Vec3, p: Vec3, angle: f64) -> CrvIn {
    let center = c + z * (p - c).dot(z);
    let x = p - center;
    let y = z.cross(x);
    let at = |t: f64, scale: f64| center + (x * t.cos() + y * t.sin()) * scale;
    let half = angle / 2.0;
    let w = half.cos();
    let (p0, p1, p2) = (at(0.0, 1.0), at(half, 1.0 / w), at(angle, 1.0));
    CrvIn {
        dim: 3,
        cpts: vec![
            [p0.x, p0.y, p0.z, 1.0],
            [p1.x * w, p1.y * w, p1.z * w, w],
            [p2.x, p2.y, p2.z, 1.0],
        ],
        p: 2,
        u: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
    }
}

/// `p` rotated by `angle` about the z axis through the origin.
fn rot_z(p: [f64; 3], angle: f64) -> [f64; 3] {
    let [x, y, z] = p;
    let (s, c) = angle.sin_cos();
    [c * x - s * y, s * x + c * y, z]
}

/// A polynomial curve with control points `pts` (degree `min(3, n - 1)`, clamped
/// uniform knots).
fn poly(pts: &[[f64; 3]]) -> CrvIn {
    let n = pts.len().saturating_sub(1);
    let p = n.min(3);
    let mut u = vec![0.0; p + 1];
    for i in 1..(n + 1 - p) {
        u.push(i as f64 / (n + 1 - p) as f64);
    }
    u.extend(vec![1.0; p + 1]);
    CrvIn {
        dim: 3,
        cpts: pts.iter().map(|&[x, y, z]| [x, y, z, NOW]).collect(),
        p,
        u,
    }
}

/// The reversed curve (control points and knots).
fn reversed(c: &CrvIn) -> CrvIn {
    let (a, b) = (
        c.u.first().copied().unwrap_or(0.0),
        c.u.last().copied().unwrap_or(1.0),
    );
    CrvIn {
        dim: c.dim,
        cpts: c.cpts.iter().rev().copied().collect(),
        p: c.p,
        u: c.u.iter().rev().map(|t| a + b - t).collect(),
    }
}

/// A boundary loop of a surface of revolution about the z axis: the generating curve `G`
/// (in the xz plane), the arc from its end, `G` rotated (reversed), the arc back (absent
/// when `G` starts on the axis: a pole). The loop is `[G, arc_end, rot(G)^-1,
/// arc_start^-1]`, so the curve before each arc is a generating curve.
fn gen_revolution(rng: &mut Rng) -> Inputs {
    let s = scale_of(rng);
    let angle = rng.uniform(0.2, std::f64::consts::FRAC_PI_2);
    let pole = rng.chance(0.2);
    let n = 2 + rng.below(3) as usize;
    let gen: Vec<[f64; 3]> = (0..=n)
        .map(|i| {
            let t = i as f64 / n as f64;
            let r = if pole && i == 0 {
                0.0
            } else {
                rng.uniform(0.5, 2.0)
            };
            [r * s, 0.0, (t * 2.0 + rng.uniform(-0.1, 0.1)) * s]
        })
        .collect();
    let g = poly(&gen);
    let g_end = v3(gen.last().copied().unwrap_or_default());
    let g_start = v3(gen.first().copied().unwrap_or_default());
    let (o, z) = (Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0));
    let arc_end = arc(o, z, g_end, angle);
    let rot: Vec<[f64; 3]> = gen.iter().map(|&p| rot_z(p, angle)).collect();
    let g_rot = reversed(&poly(&rot));
    let mut curves = vec![g, arc_end, g_rot];
    if !pole {
        curves.push(reversed(&arc(o, z, g_start, angle)));
    }
    let total = curves.len();
    // the major arc: the arc at G's end (index 1) mostly
    let mut major = if rng.chance(0.8) {
        1
    } else {
        rng.below(total as u64) as usize
    };
    match rng.below(10) {
        // the opposite arc of a different span (4 curves)
        0 => {
            let shorter = reversed(&arc(o, z, g_start, angle * rng.uniform(0.5, 0.9)));
            if let Some(c) = curves.get_mut(3) {
                *c = shorter;
            }
        }
        // the curve after the major arc off the surface
        1 => {
            if let Some(c) = curves.get_mut((major + 1) % total) {
                for p in &mut c.cpts {
                    // homogeneous: shift x by 5% of the scale (weight 1 when non-rational)
                    let w = if p[3] == NOW { 1.0 } else { p[3] };
                    p[0] += 0.05 * s * w;
                }
            }
        }
        // a major curve that is no arc
        2 => {
            let other = gen_curve(rng, 3, s);
            if let Some(c) = curves.get_mut(major) {
                *c = other;
            }
        }
        // a random loop of curves
        3 => {
            let m = 1 + rng.below(5) as usize;
            curves = (0..m).map(|_| gen_curve(rng, 3, s)).collect();
            major = rng.below(m as u64) as usize;
        }
        _ => {}
    }
    Inputs::Revolution {
        curves,
        major,
        tol: *rng.choose(&[1e-9, 1e-6, 1e-4, 1e-2]).unwrap_or(&1e-6) * s,
    }
}

fn unit_vec(rng: &mut Rng) -> [f64; 3] {
    loop {
        let v = rng.point3(-1.0, 1.0);
        let l = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        if l > 0.1 {
            return v.map(|c| c / l);
        }
    }
}

fn gen_fillet_radius(rng: &mut Rng) -> Inputs {
    let s = scale_of(rng);
    let n1 = unit_vec(rng);
    let n2 = match rng.below(8) {
        0 => n1,
        1 => n1.map(|c| -c),
        _ => unit_vec(rng),
    };
    let p1 = rng.point3(-1.0, 1.0).map(|c| c * s);
    let p2 = if rng.chance(0.1) {
        p1
    } else {
        rng.point3(-1.0, 1.0).map(|c| c * s)
    };
    Inputs::FilletRadius {
        p1,
        n1: n1.map(|c| c * rng.uniform(0.5, 2.0)),
        p2,
        n2,
        orient: *rng.choose(&[-1.0, 1.0, 0.0]).unwrap_or(&1.0),
    }
}

/// A B-spline surface closed in u and/or v (the first and last control rows equal).
fn gen_closed_surface(rng: &mut Rng) -> SrfIn {
    let p = 1 + rng.below(2) as usize;
    let q = 1 + rng.below(2) as usize;
    let n = p + 1 + rng.below(3) as usize;
    let m = q + 1 + rng.below(3) as usize;
    let close_u = rng.chance(0.7);
    let close_v = !close_u || rng.chance(0.4);
    let mut net: Vec<Vec<[f64; 4]>> = (0..=n)
        .map(|i| {
            (0..=m)
                .map(|j| {
                    let a = std::f64::consts::TAU * i as f64 / n as f64;
                    let b = std::f64::consts::TAU * j as f64 / m as f64;
                    [
                        (2.0 + b.cos()) * a.cos(),
                        (2.0 + b.cos()) * a.sin(),
                        j as f64 + b.sin(),
                        NOW,
                    ]
                })
                .collect()
        })
        .collect();
    if close_u {
        let first = net.first().cloned();
        if let (Some(first), Some(last)) = (first, net.last_mut()) {
            *last = first;
        }
    }
    if close_v {
        for row in &mut net {
            let first = row.first().copied();
            if let (Some(first), Some(last)) = (first, row.last_mut()) {
                *last = first;
            }
        }
    }
    let knots = |deg: usize, cnt: usize| {
        let mut u = vec![0.0; deg + 1];
        for i in 1..(cnt + 1 - deg) {
            u.push(i as f64 / (cnt + 1 - deg) as f64);
        }
        u.extend(vec![1.0; deg + 1]);
        u
    };
    serde_json::from_value(json!({"bspline": {
        "net": net, "p": p, "q": q, "u": knots(p, n), "v": knots(q, m)
    }}))
    .unwrap_or(SrfIn::Plane {
        o: [0.0; 3],
        n: [0.0, 0.0, 1.0],
    })
}

fn gen_seam_surface(rng: &mut Rng) -> SrfIn {
    if rng.chance(0.5) {
        gen_closed_surface(rng)
    } else {
        gen_surface(rng)
    }
}

fn gen_seam_out(rng: &mut Rng) -> Inputs {
    let surface = gen_seam_surface(rng);
    let [u0, v0, u1, v1] = domain(&surface);
    let (du, dv) = (u1 - u0, v1 - v0);
    let uv = match rng.below(5) {
        0 => [
            u0 + rng.uniform(0.0, 1.0) * du,
            v0 + rng.uniform(0.0, 1.0) * dv,
        ],
        1 => [
            u0 - rng.uniform(0.0, 0.5) * du,
            v0 + rng.uniform(0.0, 1.0) * dv,
        ],
        2 => [
            u1 + rng.uniform(0.0, 0.5) * du,
            v1 + rng.uniform(0.0, 0.5) * dv,
        ],
        3 => [-1e300, v0],
        _ => [
            u0 + rng.uniform(-1.0, 2.0) * du,
            v0 + rng.uniform(-1.0, 2.0) * dv,
        ],
    };
    Inputs::SeamOut {
        surface,
        uv,
        tol: *rng.choose(&[0.0, 1e-9, 1e-6, 1e-2]).unwrap_or(&1e-6),
        closure: rng.below(5) as i64,
    }
}

fn gen_closure(rng: &mut Rng) -> Inputs {
    let surface = gen_seam_surface(rng);
    let [u0, v0, u1, v1] = domain(&surface);
    let dom = if rng.chance(0.7) {
        [u0, v0, u1, v1]
    } else {
        let a = rng.uniform(0.0, 0.4);
        let b = rng.uniform(0.6, 1.0);
        [
            u0 + a * (u1 - u0),
            v0 + a * (v1 - v0),
            u0 + b * (u1 - u0),
            v0 + b * (v1 - v0),
        ]
    };
    Inputs::Closure { surface, dom }
}

pub(super) fn gen_topo_brep(rng: &mut Rng) -> BrepIn {
    if rng.chance(0.8) {
        BrepIn::File {
            file: (*rng.choose(&CORNER_FILES).unwrap_or(&"Fillet/box1x1x2.smb")).to_string(),
            rebuild: false,
        }
    } else {
        gen_brep(rng)
    }
}

/// The fillet test files with concave corners (most vertices of the others are convex,
/// where the concavity test is `false`).
const CONCAVE_FILES: [&str; 5] = [
    "Fillet/l_shaped_box.smb",
    "Fillet/box_w_slot.smb",
    "Fillet/box_on_box.smb",
    "Fillet/6_sided_corner.smb",
    "Fillet/corner_box.smb",
];

fn concave_brep(rng: &mut Rng) -> BrepIn {
    BrepIn::File {
        file: (*rng
            .choose(&CONCAVE_FILES)
            .unwrap_or(&"Fillet/l_shaped_box.smb"))
        .to_string(),
        rebuild: false,
    }
}

fn gen_topology(rng: &mut Rng) -> Inputs {
    let brep = gen_topo_brep(rng);
    let s = |rng: &mut Rng| rng.range_i64(0, 1 << 20);
    let none_or = |rng: &mut Rng, p: f64| {
        if rng.chance(p) {
            -1
        } else {
            rng.range_i64(0, 1 << 20)
        }
    };
    match rng.below(4) {
        0 => Inputs::G1 {
            brep,
            v: none_or(rng, 0.05),
            e1: none_or(rng, 0.05),
            e2: none_or(rng, 0.05),
            tol: *rng.choose(&[0.5, 1.0, 5.0, 45.0, 90.0]).unwrap_or(&1.0),
        },
        1 => Inputs::FaceByTwo {
            brep,
            v: s(rng),
            eu1: none_or(rng, 0.05),
            eu2: none_or(rng, 0.05),
            eu3: none_or(rng, 0.6),
        },
        2 => Inputs::Tangential {
            brep,
            eu: s(rng),
            tol: *rng.choose(&[0.1, 2.0, 5.0, 45.0]).unwrap_or(&2.0),
        },
        _ => Inputs::Concavity {
            brep: if rng.chance(0.7) {
                concave_brep(rng)
            } else {
                brep
            },
            v: s(rng),
            eu1: s(rng),
            eu2: s(rng),
            ref_eu: none_or(rng, 0.8),
            adjacent: rng.chance(0.7),
        },
    }
}

impl UnitDiff for Runner {
    fn id(&self) -> &'static str {
        "blend.fillet_corner_p1"
    }

    fn run_cpp(&self, case: &Case) -> Result<Output, HarnessError> {
        run_cpp(&case.parse::<Inputs>()?)
    }

    fn run_rust(&self, case: &Case) -> Result<Output, HarnessError> {
        run_rust(&case.parse::<Inputs>()?)
    }

    fn generate(&self, rng: &mut Rng) -> Option<Value> {
        let i = match rng.below(12) {
            0 | 1 => gen_planar(rng),
            2..=4 => gen_revolution(rng),
            5 => gen_fillet_radius(rng),
            6 => gen_seam_out(rng),
            7 => gen_closure(rng),
            _ => gen_topology(rng),
        };
        serde_json::to_value(i).ok()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn inputs_round_trip() {
        let mut rng = Rng::new(11);
        for _ in 0..200 {
            let v = Runner.generate(&mut rng).unwrap();
            let back: Inputs =
                serde_json::from_value(v.clone()).unwrap_or_else(|e| panic!("{e}: {v}"));
            assert_eq!(serde_json::to_value(back).unwrap(), v);
        }
    }

    /// The op name of a case.
    fn op_name(case: &Case) -> String {
        case.inputs
            .get("op")
            .and_then(Value::as_str)
            .unwrap_or("?")
            .to_string()
    }

    /// With `CASE_CENSUS`: the outcome of every case of the unit run (fixtures plus
    /// `CASE_SEEDS` generated, default 16) by op: `deferred`, `refused` / other errors,
    /// `error_status` (a non-zero `SmStatus`) or the result (`success`, and for the
    /// flag ops `true` / `false`, for `face_by_two` `face` / `no_face`, for `seam_out`
    /// `moved` / `kept`); with `CASE_LIST=<kind>` as well, prints the name and inputs
    /// of every case of that kind (to pin generated cases as fixtures). Without it: runs
    /// the cases `CASE_INDICES` (or the first `CASE_PER_OP` of each op) on one side
    /// (`CASE_SIDE=rust`, else C++), for valgrind.
    #[test]
    #[ignore = "census and valgrind sample; run on demand"]
    fn case_sample() {
        let env = |name: &str| std::env::var(name).ok();
        let ctx = crate::case::CaseContext {
            root: crate::case::repo_root(),
            seeds: env("CASE_SEEDS").map_or(16, |s| s.parse().unwrap()),
        };
        let cases = Runner.cases(&ctx).unwrap();
        if env("CASE_CENSUS").is_some() {
            let mut tally = std::collections::BTreeMap::<(String, String), usize>::new();
            for case in &cases {
                let kind = match Runner.run_rust(case) {
                    Ok(o) if o.deferred.is_some() => "deferred".to_string(),
                    Ok(o) if o.error.is_some() => {
                        let e = o.error.unwrap_or_default();
                        e.split(|c: char| c == ':' || c.is_ascii_digit())
                            .next()
                            .unwrap_or("")
                            .to_string()
                    }
                    Ok(o) => {
                        let first = |k: &str| o.values.get(k).and_then(|v| v.first()).copied();
                        if first("status").is_some_and(|st| st != 0.0) {
                            "error_status".to_string()
                        } else if let Some(f) = first("flag") {
                            (if f == 0.0 { "false" } else { "true" }).to_string()
                        } else if let Some(f) = first("g1") {
                            (if f == 0.0 { "false" } else { "true" }).to_string()
                        } else if let Some(f) = first("face") {
                            (if f >= 0.0 { "face" } else { "no_face" }).to_string()
                        } else if let Some(f) = first("seam") {
                            (if f == 0.0 { "kept" } else { "moved" }).to_string()
                        } else {
                            "success".to_string()
                        }
                    }
                    Err(_) => "harness_error".to_string(),
                };
                if env("CASE_LIST").is_some_and(|k| k == kind) {
                    // the resolved arguments and, for a vertex op, the edge uses at the
                    // vertex (export order), so a seeded case can be pinned with small
                    // selectors
                    let mut resolved = String::new();
                    if let Ok(i) = serde_json::from_value::<Inputs>(case.inputs.clone()) {
                        if let Ok(Ok(t)) = topo(&i) {
                            let vid = t.args.first().and_then(|&v| u32::try_from(v).ok());
                            let ring: Vec<u32> = t
                                .brep
                                .edgeuses()
                                .iter()
                                .map(|(id, _)| id)
                                .filter(|&id| eu_vertex(&t.brep, id).map(VertexId::index) == vid)
                                .map(EdgeuseId::index)
                                .collect();
                            resolved = format!(
                                "args={:?} ring={ring:?} surfaces={:?} curves={:?}",
                                t.args, t.types.face_surface_types, t.types.edge_curve_types
                            );
                        }
                    }
                    eprintln!(
                        "case_item: {} {} {} {resolved}",
                        op_name(case),
                        case.name,
                        case.inputs
                    );
                }
                *tally.entry((op_name(case), kind)).or_default() += 1;
            }
            for ((op, kind), n) in tally {
                eprintln!("case_census: {op} {kind} {n}");
            }
            return;
        }
        let rust_side = env("CASE_SIDE").is_some_and(|s| s == "rust");
        let indices: Vec<usize> = if let Some(t) = env("CASE_INDICES") {
            t.split(',').map(|x| x.trim().parse().unwrap()).collect()
        } else {
            let per_op = env("CASE_PER_OP").map_or(10, |s| s.parse().unwrap());
            let mut taken = std::collections::BTreeMap::<String, usize>::new();
            (0..cases.len())
                .filter(|&i| {
                    let n = taken.entry(op_name(&cases[i])).or_default();
                    *n += 1;
                    *n <= per_op
                })
                .collect()
        };
        for i in indices {
            let case = &cases[i];
            let r = if rust_side {
                Runner.run_rust(case).map(|o| o.error)
            } else {
                Runner.run_cpp(case).map(|o| o.error)
            };
            eprintln!("case_sample: {i} {} {r:?}", case.name);
        }
    }

    #[test]
    fn arc_is_exact() {
        let c = arc(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(2.0, 0.0, 1.0),
            1.0,
        );
        let p = c.cpts[1];
        // the middle control point lies on the bisector at r / cos(a/2)
        let (x, y) = (p[0] / p[3], p[1] / p[3]);
        assert!(((x * x + y * y).sqrt() - 2.0 / 0.5f64.cos()).abs() < 1e-12);
    }
}
