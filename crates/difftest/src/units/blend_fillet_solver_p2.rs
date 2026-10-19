//! Runner for unit `blend.fillet_solver_p2`. C++: `solidcore_ffi::units::blend_fillet_solver_p2`.
//! Rust: `solidcore::blend::fillet_solver::p2`.
//!
//! The solves and the analytic fillet builders call the solver's virtuals
//! (`GetJacobianSize`, `LoadInitialValues`, `LoadJacobian`, `SetupOffsetValues`), which the
//! library's solvers inherit from units that are not ported. Both sides therefore run a
//! small rolling-ball solver: the shim's `P2FS` and its mirror [`TestSolver`] (offset
//! surface `rail` minus offset surface `1 - rail` as three equations in the two UV
//! points, plus `extra` linear unknowns), on the same offset surfaces of planes and
//! B-spline patches (blocks of unit `topology.edgeuse_p1`).
//!
//! Ops: `PointOnPlaneSolve` (`plane`: the rolling ball between a floor and a wall patch
//! cut by a plane across the edge, the plane adjusting on the second try; curved
//! patches; up to four user doubles), `RailEdgeuseIntersect` (`rail_eu`: the rail on a
//! face B-rep's surface meeting its side edges, either rail, missing offset surfaces
//! (`SM_ERR`), file B-reps; the side edge use's UV trim curve is built on demand),
//! `RailRailIntersect` (`rail_rail`: missing surfaces compare `SM_ERR`; with all four
//! surfaces the Rust side defers at `SmTemporaryRadiusChange`, unit
//! `blend.fillet_executive_p2`), the analytic builders `CalcCylinderFilletGeom`,
//! `CalcConeFilletGeom` and `CalcTorusFilletGeom` (`geom`: a circular fillet between a
//! floor and a wall patch along a straight wire edge, or between a floor and a tangent
//! wall along a rational circular arc edge (quarter to full circle); constant and
//! varying radius, extension before/after, swapped surfaces, a linear generator, the
//! wrong edge kind, jittered points; the C++ stores into a fresh `SmFilletGeom` with two
//! rails, the Rust side compares the built parts, [`FilletSolver::cylinder_fillet_geom_parts`]
//! etc.), `GetFilletCorner` above index 1 (`corner`) and `RecordFilletError`
//! (`record`: the shim skips it, upstream reaches the `NULL` executive; the Rust side
//! defers at `SmFilletGeom::GetFilletSurface`, unit `blend.fillet_geom_p1`).
//!
//! B-rep cases build the C++ B-rep, export it through unit `topology.edgeuse_p1` and
//! import the same records on the Rust side; entity selectors are resolved modulo the
//! export's counts.

use serde::{Deserialize, Serialize};
use serde_json::Value;

use solidcore::blend::fillet_custom_solver::{
    InitialValues, JacobianLoad, OffsetSurface, SurfaceId,
};
use solidcore::blend::fillet_solver::p1::{
    FilletSolver, FilletSolverOps, FilletSurfaceGenerator, OffsetSurfaceSlot,
};
use solidcore::blend::fillet_solver::p2::FilletGeomParts;
use solidcore::curve::curve::Curve;
use solidcore::geom::extent1d::{Extent1d, ExtentNd};
use solidcore::geom::extent2d::Extent2d;
use solidcore::geom::matrix::Matrix;
use solidcore::geom::vector2d::Point2;
use solidcore::geom::vector3d::Vec3;
use solidcore::intersect::surface_intersector::p1::{
    IntersectionPointType, SrfIndex, TsectPnt, MAX_USER_DOUBLES,
};
use solidcore::kernel::Kernel;
use solidcore::nurbs::knots::Side;
use solidcore::nurbs::math::NOW;
use solidcore::surface::offset_surface::SS_UNKNOWN;
use solidcore::surface::surface::p1::evaluate_1st_derivatives;
use solidcore::surface::surface::Surface;
use solidcore::topology::{Brep, EdgeuseId};
use solidcore::SolidError;
use solidcore_ffi::units::blend_fillet_solver_p2::{self as ffi, Op};
use solidcore_ffi::units::topology_edgeuse_p1 as ffi1;

use super::topology_edgeuse_p1::{
    brep_block_of, gen_surface, import, path_of, push_surface, sel, surface_of, BrepIn, BssIn,
    CrvIn, SrfIn, FILES,
};
use crate::{HarnessError, Output, Rng, UnitDiff};

/// The runner.
pub struct Runner;

fn harness(msg: impl Into<String>) -> HarnessError {
    HarnessError::new(msg)
}

/// The shim solver's parameters (`tol extra radius slope target gu gv`).
#[derive(Debug, Clone, Serialize, Deserialize)]
struct Cfg {
    tol: f64,
    extra: usize,
    radius: f64,
    slope: f64,
    target: f64,
    gu: f64,
    gv: f64,
}

/// Offset surface `i` of `base[i]` when `has[i]`, and the orientations.
#[derive(Debug, Clone, Serialize, Deserialize)]
struct Surfs {
    base: [SrfIn; 2],
    has: [bool; 2],
    off: [f64; 2],
    orient: [f64; 2],
}

/// A tsect point: both UV points and the curve parameter.
#[derive(Debug, Clone, Serialize, Deserialize)]
struct Tsect {
    uv0: [f64; 2],
    uv1: [f64; 2],
    cparam: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "op", rename_all = "snake_case")]
enum Inputs {
    Plane {
        cfg: Cfg,
        surfs: Surfs,
        orig: [f64; 3],
        normal: [f64; 3],
        dom1: [f64; 4],
        dom2: [f64; 4],
        uv1: [f64; 2],
        uv2: [f64; 2],
    },
    RailEu {
        cfg: Cfg,
        surfs: Surfs,
        brep: BrepIn,
        side: i64,
        frac: f64,
        rail: usize,
        ptype: u8,
    },
    RailRail {
        cfg: [Cfg; 2],
        surfs: [Surfs; 2],
        uv: [f64; 2],
        rail: usize,
        other: usize,
    },
    /// `which`: 0 cylinder, 1 cone, 2 torus; `fsg`: 0 linear, 1 circular.
    Geom {
        which: u8,
        cfg: Cfg,
        surfs: Surfs,
        fsg: u8,
        before: bool,
        after: bool,
        factor: f64,
        brep: BrepIn,
        eu: [i64; 2],
        tsects: Vec<Tsect>,
    },
    Corner {
        cfg: Cfg,
        index: u64,
    },
    Record {
        cfg: Cfg,
        brep: BrepIn,
        eu: i64,
    },
}

fn b(v: bool) -> f64 {
    if v {
        1.0
    } else {
        0.0
    }
}

fn push_cfg(o: &mut Vec<f64>, c: &Cfg) {
    o.extend([
        c.tol,
        c.extra as f64,
        c.radius,
        c.slope,
        c.target,
        c.gu,
        c.gv,
    ]);
}

fn push_surfs(o: &mut Vec<f64>, s: &Surfs) -> Result<(), HarnessError> {
    push_surface(o, &s.base[0])?;
    push_surface(o, &s.base[1])?;
    for (h, d) in s.has.iter().zip(&s.off) {
        o.push(b(*h));
        o.push(*d);
    }
    o.extend(s.orient);
    Ok(())
}

// ---- B-rep cases -------------------------------------------------------------------------

/// The case's B-rep: the C++ export's counts and the imported Rust B-rep.
struct BrepCase {
    path: String,
    block: Vec<f64>,
    counts: [usize; 10],
    brep: Brep,
}

/// Builds the case's B-rep on both sides; `Err(output)` when it cannot be built.
fn brep_case(brep: &BrepIn) -> Result<Result<BrepCase, Output>, HarnessError> {
    let path = path_of(Some(brep))?;
    let block = brep_block_of(Some(brep), false, false)?;
    let ex = ffi1::eval(ffi1::Op::Export, &path, &block)
        .map_err(|e| harness(format!("c++ export: {e}")))?;
    let Some(counts) = ffi1::counts(&ex) else {
        let st = ex.get(1).copied().unwrap_or(0.0);
        return Ok(Err(Output::err(format!("brep_build_status_{st}"))));
    };
    let brep = import(&ex)?;
    Ok(Ok(BrepCase {
        path,
        block,
        counts,
        brep,
    }))
}

impl Inputs {
    fn brep(&self) -> Option<&BrepIn> {
        match self {
            Inputs::RailEu { brep, .. }
            | Inputs::Geom { brep, .. }
            | Inputs::Record { brep, .. } => Some(brep),
            _ => None,
        }
    }

    fn op(&self) -> Op {
        match self {
            Inputs::Plane { .. } => Op::Plane,
            Inputs::RailEu { .. } => Op::RailEu,
            Inputs::RailRail { .. } => Op::RailRail,
            Inputs::Geom { .. } => Op::Geom,
            Inputs::Corner { .. } | Inputs::Record { .. } => Op::Misc,
        }
    }
}

/// The edge use counts of the case's export (`[n_eu]`).
fn n_eu(c: &BrepCase) -> usize {
    c.counts[6]
}

fn cpp_input(i: &Inputs, c: Option<&BrepCase>) -> Result<Vec<f64>, HarnessError> {
    let mut o = Vec::new();
    let block = || -> Result<&BrepCase, HarnessError> { c.ok_or_else(|| harness("no B-rep")) };
    match i {
        Inputs::Plane {
            cfg,
            surfs,
            orig,
            normal,
            dom1,
            dom2,
            uv1,
            uv2,
        } => {
            push_cfg(&mut o, cfg);
            push_surfs(&mut o, surfs)?;
            o.extend(orig);
            o.extend(normal);
            o.extend(dom1);
            o.extend(dom2);
            o.extend(uv1);
            o.extend(uv2);
        }
        Inputs::RailEu {
            cfg,
            surfs,
            side,
            frac,
            rail,
            ptype,
            ..
        } => {
            let c = block()?;
            push_cfg(&mut o, cfg);
            push_surfs(&mut o, surfs)?;
            o.extend_from_slice(&c.block);
            o.extend([
                sel(*side, n_eu(c)) as f64,
                *frac,
                *rail as f64,
                f64::from(*ptype),
            ]);
        }
        Inputs::RailRail {
            cfg,
            surfs,
            uv,
            rail,
            other,
        } => {
            for (c, s) in cfg.iter().zip(surfs) {
                push_cfg(&mut o, c);
                push_surfs(&mut o, s)?;
            }
            o.extend(uv);
            o.extend([*rail as f64, *other as f64]);
        }
        Inputs::Geom {
            which,
            cfg,
            surfs,
            fsg,
            before,
            after,
            factor,
            eu,
            tsects,
            ..
        } => {
            let c = block()?;
            o.push(f64::from(*which));
            push_cfg(&mut o, cfg);
            push_surfs(&mut o, surfs)?;
            o.extend([f64::from(*fsg), b(*before), b(*after), *factor]);
            o.extend_from_slice(&c.block);
            o.push(sel(eu[0], n_eu(c)) as f64);
            o.push(sel(eu[1], n_eu(c)) as f64);
            o.push(tsects.len() as f64);
            for t in tsects {
                o.extend(t.uv0);
                o.extend(t.uv1);
                o.push(t.cparam);
            }
        }
        Inputs::Corner { cfg, index } => {
            o.push(0.0);
            push_cfg(&mut o, cfg);
            o.push(*index as f64);
        }
        Inputs::Record { cfg, eu, .. } => {
            let c = block()?;
            o.push(1.0);
            push_cfg(&mut o, cfg);
            o.extend_from_slice(&c.block);
            o.push(sel(*eu, n_eu(c)) as f64);
        }
    }
    Ok(o)
}

fn cpp(i: &Inputs) -> Result<Output, HarnessError> {
    let case = match i.brep() {
        Some(brep) => match brep_case(brep)? {
            Ok(c) => Some(c),
            Err(out) => return Ok(out),
        },
        None => None,
    };
    let path = case.as_ref().map(|c| c.path.clone()).unwrap_or_default();
    let input = cpp_input(i, case.as_ref())?;
    let flat = ffi::eval(i.op(), &path, &input).map_err(|e| harness(e.to_string()))?;
    Ok(Output::new().values("flat", flat))
}

// ---- Rust side ------------------------------------------------------------------------

fn status(e: &SolidError) -> f64 {
    match e {
        SolidError::InvalidInput { .. } => 1007.0,
        SolidError::Upstream { code, .. } => f64::from(*code),
        _ => -1.0,
    }
}

/// The shim's `P2FS` (see its header): `GetJacobianSize` `4 + extra`; `SetupOffsetValues`
/// sets every present offset surface to `(radius + slope * curve_parameter) *
/// orientation`; `LoadInitialValues` appends the other surface's UV guess (the first
/// surface's plus `(gu, gv)`, extents its base domain) and the `extra` unknowns;
/// `LoadJacobian` loads `S_rail(uv1) - S_other(uv2)` and `x[k] - (target + i)`.
struct TestSolver {
    base: FilletSolver,
    extra: usize,
    radius: f64,
    slope: f64,
    target: f64,
    gu: f64,
    gv: f64,
}

fn slot_mut(s: &mut FilletSolver, rail: usize) -> Result<&mut OffsetSurfaceSlot, SolidError> {
    s.surfaces
        .get_mut(rail)
        .and_then(Option::as_mut)
        .ok_or_else(|| SolidError::invalid_input("difftest::TestSolver", "no offset surface"))
}

fn set_ext(e: &mut ExtentNd, i: usize, v: Extent1d) -> Result<(), SolidError> {
    let slot = e
        .get_mut(i)
        .ok_or_else(|| SolidError::invalid_input("difftest::TestSolver", "extent index"))?;
    *slot = v;
    Ok(())
}

fn at(x: &[f64], i: usize) -> Result<f64, SolidError> {
    x.get(i)
        .copied()
        .ok_or_else(|| SolidError::invalid_input("difftest::TestSolver", "index"))
}

impl FilletSolverOps for TestSolver {
    fn solver(&self) -> &FilletSolver {
        &self.base
    }

    fn solver_mut(&mut self) -> &mut FilletSolver {
        &mut self.base
    }

    fn jacobian_size(&self) -> Result<usize, SolidError> {
        Ok(4 + self.extra)
    }

    fn setup_offset_values(
        &mut self,
        tsect: &TsectPnt,
        offset_dist: Option<&mut f64>,
    ) -> Result<(), SolidError> {
        let r = self.radius + self.slope * tsect.curve_parameter;
        for i in 0..2 {
            let o = self.base.orientation(i)?;
            if let Some(slot) = self.base.surfaces.get_mut(i).and_then(Option::as_mut) {
                slot.surface.set_offset_distance(r * o);
            }
        }
        if let Some(d) = offset_dist {
            *d = r;
        }
        Ok(())
    }

    fn load_initial_values(
        &mut self,
        _k: &Kernel,
        _brep: Option<&Brep>,
        do_surf2_calcs: bool,
        rail_index: usize,
        _surface1: &mut OffsetSurface,
        surface1_index: usize,
        intervals: &mut ExtentNd,
        periodicities: &mut Vec<bool>,
        guess_t: &mut Vec<f64>,
    ) -> Result<Option<InitialValues>, SolidError> {
        let mut out = None;
        if do_surf2_calcs {
            let other = 1 - rail_index.min(1);
            let dom = slot_mut(&mut self.base, other)?
                .surface
                .base_surface()
                .natural_uv_domain();
            let idx2 = guess_t.len();
            set_ext(
                intervals,
                idx2,
                Extent1d::from_min_max(dom.min().x, dom.max().x),
            )?;
            set_ext(
                intervals,
                idx2 + 1,
                Extent1d::from_min_max(dom.min().y, dom.max().y),
            )?;
            periodicities.push(false);
            periodicities.push(false);
            let u = at(guess_t, surface1_index)? + self.gu;
            let v = at(guess_t, surface1_index + 1)? + self.gv;
            guess_t.push(u);
            guess_t.push(v);
            out = Some(InitialValues {
                surface2_rail: other,
                surface2_index: idx2,
            });
        }
        for i in 0..self.extra {
            set_ext(
                intervals,
                guess_t.len(),
                Extent1d::from_min_max(-100.0, 100.0),
            )?;
            periodicities.push(false);
            guess_t.push(0.5 * i as f64);
        }
        Ok(out)
    }

    fn load_jacobian(
        &mut self,
        num_equations: &mut usize,
        num_parameters: &mut usize,
        load: &mut JacobianLoad<'_>,
        f: &mut [f64],
        mut jacobian: Option<&mut Matrix>,
    ) -> Result<bool, SolidError> {
        let rail = load.rail_index.min(1);
        let other = 1 - rail;
        load.surface2 = Some(slot_mut(&mut self.base, other)?.surface.base_id);
        load.surf2_offset = *num_parameters;
        let (o1, o2) = (load.surf1_offset, load.surf2_offset);
        let x = load.x;
        let uv_a = Point2::new(at(x, o1)?, at(x, o1 + 1)?);
        let uv_b = Point2::new(at(x, o2)?, at(x, o2 + 1)?);
        let (pa, dua, dva) = slot_mut(&mut self.base, rail)?
            .surface
            .surface
            .with_evaluator(|e| evaluate_1st_derivatives(e, uv_a, Side::Left, Side::Left))?;
        let (pb, dub, dvb) = slot_mut(&mut self.base, other)?
            .surface
            .surface
            .with_evaluator(|e| evaluate_1st_derivatives(e, uv_b, Side::Left, Side::Left))?;
        let mut scale = pa.x.abs();
        if pa.y.abs() > scale {
            scale = pa.y.abs();
        }
        if pa.z.abs() > scale {
            scale = pa.z.abs();
        }
        let tol = self.base.conversion_tol * (1.0 + scale);
        let mut found = true;
        // per row: [S_a, dS_a/du, dS_a/dv, S_b, dS_b/du, dS_b/dv] components
        let rows = [
            [pa.x, dua.x, dva.x, pb.x, dub.x, dvb.x],
            [pa.y, dua.y, dva.y, pb.y, dub.y, dvb.y],
            [pa.z, dua.z, dva.z, pb.z, dub.z, dvb.z],
        ];
        let ne = *num_equations;
        for (c, [sa, ua, wa, sb, ub, wb]) in rows.into_iter().enumerate() {
            let v = sa - sb;
            *f.get_mut(ne + c)
                .ok_or_else(|| SolidError::invalid_input("difftest::TestSolver", "rF"))? = v;
            if v.abs() >= tol {
                found = false;
            }
            if let Some(j) = jacobian.as_deref_mut() {
                j.set_at(ne + c, o1, ua)?;
                j.set_at(ne + c, o1 + 1, wa)?;
                j.set_at(ne + c, o2, -ub)?;
                j.set_at(ne + c, o2 + 1, -wb)?;
            }
        }
        *num_equations += 3;
        *num_parameters += 2;
        for i in 0..self.extra {
            let (ne, np) = (*num_equations, *num_parameters);
            let v = at(x, np)? - (self.target + i as f64);
            *f.get_mut(ne)
                .ok_or_else(|| SolidError::invalid_input("difftest::TestSolver", "rF"))? = v;
            if v.abs() >= tol {
                found = false;
            }
            if let Some(j) = jacobian.as_deref_mut() {
                j.set_at(ne, np, 1.0)?;
            }
            *num_equations += 1;
            *num_parameters += 1;
        }
        Ok(found)
    }
}

/// The solver of `cfg` with the offset surfaces of `s` (distinct bases are distinct
/// objects).
fn solver_of(cfg: &Cfg, s: &Surfs) -> Result<TestSolver, HarnessError> {
    let k = Kernel::default();
    let err = |e: SolidError| harness(e.to_string());
    let mut fs = FilletSolver::new(cfg.tol, 0.1, 0.01);
    for (i, ((has, off), base)) in s.has.iter().zip(&s.off).zip(&s.base).enumerate() {
        if *has {
            let base = surface_of(base)?;
            let surface = OffsetSurface::new(&k, *off, base, SurfaceId::Free(i as u32), SS_UNKNOWN)
                .map_err(err)?;
            fs.set_surface(
                i,
                Some(OffsetSurfaceSlot {
                    key: i as u64,
                    surface,
                }),
            )
            .map_err(err)?;
        }
    }
    fs.orientations = s.orient;
    Ok(TestSolver {
        base: fs,
        extra: cfg.extra,
        radius: cfg.radius,
        slope: cfg.slope,
        target: cfg.target,
        gu: cfg.gu,
        gv: cfg.gv,
    })
}

fn put_offsets(o: &mut Vec<f64>, fs: &TestSolver) {
    for i in 0..2 {
        o.push(
            fs.base
                .surfaces
                .get(i)
                .and_then(Option::as_ref)
                .map_or(0.0, |s| s.surface.offset_distance()),
        );
    }
}

fn put_tsect(o: &mut Vec<f64>, t: &TsectPnt, n_user: usize) {
    for s in [SrfIndex::First, SrfIndex::Second] {
        let uv = t.uv_pos(s);
        o.extend([uv.x, uv.y]);
    }
    o.extend(t.user_doubles.iter().take(n_user.min(MAX_USER_DOUBLES)));
}

/// `put_crv` of the shim: `1 t0 t1 P(t0) P(mid) P(t1)`, `-1` for none.
fn put_crv(o: &mut Vec<f64>, c: Option<&Curve>) {
    let Some(c) = c else {
        o.push(-1.0);
        return;
    };
    o.push(1.0);
    let ivl = c.natural_interval();
    o.extend([ivl.min(), ivl.max()]);
    for t in [ivl.min(), ivl.mid(), ivl.max()] {
        let p = c.evaluate_point(t).unwrap_or(Vec3::UNDEF);
        o.extend([p.x, p.y, p.z]);
    }
}

/// `put_srf` of the shim: `1 u0 v0 u1 v1` and the 3 x 3 grid of points.
fn put_srf(o: &mut Vec<f64>, s: &Surface) {
    o.push(1.0);
    let d = s.natural_uv_domain();
    o.extend([d.min().x, d.min().y, d.max().x, d.max().y]);
    for a in [0.0, 0.5, 1.0] {
        for bb in [0.0, 0.5, 1.0] {
            let p = s.evaluate_point(d.evaluate(a, bb)).unwrap_or(Vec3::UNDEF);
            o.extend([p.x, p.y, p.z]);
        }
    }
}

fn put_parts(o: &mut Vec<f64>, which: u8, p: &FilletGeomParts) {
    put_srf(o, &p.surface);
    put_crv(o, Some(&p.rails[0]));
    put_crv(o, Some(&p.rails[1]));
    if which != 0 {
        for ivl in &p.rail_intervals {
            let ivl = ivl.unwrap_or_else(Extent1d::new);
            o.extend([ivl.min(), ivl.max()]);
        }
    }
    put_crv(o, p.center_line.as_ref());
}

fn eu_id(i: i64) -> Option<EdgeuseId> {
    u32::try_from(i).ok().map(EdgeuseId::from_index)
}

fn point_type(code: u8) -> IntersectionPointType {
    match code {
        0 => IntersectionPointType::Undefined,
        2 => IntersectionPointType::Crossing,
        3 => IntersectionPointType::TangentPoint,
        4 => IntersectionPointType::TangentCurve,
        5 => IntersectionPointType::Singularity,
        6 => IntersectionPointType::Coincidence,
        7 => IntersectionPointType::Touching,
        _ => IntersectionPointType::Unknown,
    }
}

fn point_type_code(t: IntersectionPointType) -> f64 {
    match t {
        IntersectionPointType::Undefined => 0.0,
        IntersectionPointType::Unknown => 1.0,
        IntersectionPointType::Crossing => 2.0,
        IntersectionPointType::TangentPoint => 3.0,
        IntersectionPointType::TangentCurve => 4.0,
        IntersectionPointType::Singularity => 5.0,
        IntersectionPointType::Coincidence => 6.0,
        IntersectionPointType::Touching => 7.0,
    }
}

fn dom2(d: &[f64; 4]) -> Extent2d {
    Extent2d::from_points(Point2::new(d[0], d[1]), Point2::new(d[2], d[3]))
}

/// Runs the case; `Err` is a status or a deferral (after the values pushed so far).
fn rust_run(i: &Inputs, o: &mut Vec<f64>) -> Result<Result<(), SolidError>, HarnessError> {
    let k = Kernel::default();
    Ok(match i {
        Inputs::Plane {
            cfg,
            surfs,
            orig,
            normal,
            dom1,
            dom2: d2,
            uv1,
            uv2,
        } => {
            let mut fs = solver_of(cfg, surfs)?;
            let mut t = TsectPnt::default();
            let r = FilletSolver::point_on_plane_solve(
                &mut fs,
                &k,
                None,
                Vec3::from_array(*orig),
                Vec3::from_array(*normal),
                &dom2(dom1),
                &dom2(d2),
                Point2::new(uv1[0], uv1[1]),
                Point2::new(uv2[0], uv2[1]),
                &mut t,
            );
            match r {
                Ok(found) => {
                    o.push(0.0);
                    put_offsets(o, &fs);
                    o.push(b(found));
                    if found {
                        put_tsect(o, &t, cfg.extra);
                    }
                    Ok(())
                }
                Err(e) => {
                    if let Some(_d) = Output::unported(&e) {
                        return Ok(Err(e));
                    }
                    o.push(status(&e));
                    put_offsets(o, &fs);
                    Ok(())
                }
            }
        }
        Inputs::RailEu {
            cfg,
            surfs,
            brep,
            side,
            frac,
            rail,
            ptype,
        } => {
            let Ok(mut c) = brep_case(brep)? else {
                return Err(harness("B-rep failed on the Rust side only"));
            };
            o.push(0.0);
            let side = eu_id(sel(*side, n_eu(&c))).ok_or_else(|| harness("side edge use"))?;
            let mut fs = solver_of(cfg, surfs)?;
            let mut t = TsectPnt {
                point_type: point_type(*ptype),
                ..TsectPnt::default()
            };
            let e = solidcore::topology::edgeuse::p1::edge(&c.brep, side)
                .map_err(|e| harness(e.to_string()))?
                .ok_or_else(|| harness("side edge use without edge"))?;
            let ivl = c.brep.edge(e).map_err(|e| harness(e.to_string()))?.interval;
            let guess = ivl.min() + (ivl.max() - ivl.min()) * *frac;
            let r = FilletSolver::rail_edgeuse_intersect(
                &mut fs,
                &k,
                &mut c.brep,
                *rail,
                side,
                guess,
                &mut t,
            );
            match r {
                Ok(param) => {
                    o.push(0.0);
                    put_offsets(o, &fs);
                    o.push(b(param.is_some()));
                    if let Some(p) = param {
                        o.push(p);
                        put_tsect(o, &t, cfg.extra);
                        o.push(t.curve_parameter);
                        o.push(point_type_code(t.point_type));
                    }
                    Ok(())
                }
                Err(e) => {
                    if Output::unported(&e).is_some() {
                        return Ok(Err(e));
                    }
                    o.push(status(&e));
                    put_offsets(o, &fs);
                    Ok(())
                }
            }
        }
        Inputs::RailRail {
            cfg,
            surfs,
            uv,
            rail,
            other,
        } => {
            let mut fa = solver_of(&cfg[0], &surfs[0])?;
            let mut fb = solver_of(&cfg[1], &surfs[1])?;
            let mut ta = TsectPnt::default();
            let mut tb = TsectPnt::default();
            let r = FilletSolver::rail_rail_intersect(
                &mut fa,
                &mut fb,
                &k,
                None,
                Point2::new(uv[0], uv[1]),
                *rail,
                *other,
                &mut ta,
                &mut tb,
            );
            match r {
                Ok(found) => {
                    o.push(0.0);
                    put_offsets(o, &fa);
                    put_offsets(o, &fb);
                    o.push(b(found));
                    if found {
                        put_tsect(o, &ta, cfg[0].extra);
                        put_tsect(o, &tb, cfg[1].extra);
                    }
                    Ok(())
                }
                Err(e) => {
                    if Output::unported(&e).is_some() {
                        return Ok(Err(e));
                    }
                    o.push(status(&e));
                    put_offsets(o, &fa);
                    put_offsets(o, &fb);
                    Ok(())
                }
            }
        }
        Inputs::Geom {
            which,
            cfg,
            surfs,
            fsg,
            before,
            after,
            factor,
            brep,
            eu,
            tsects,
        } => {
            let Ok(c) = brep_case(brep)? else {
                return Err(harness("B-rep failed on the Rust side only"));
            };
            o.push(0.0);
            let mut fs = solver_of(cfg, surfs)?;
            fs.base.set_extension_factor(*factor);
            fs.base.extend_before = *before;
            fs.base.extend_after = *after;
            fs.base.fsg = Some(if *fsg == 0 {
                FilletSurfaceGenerator::Linear
            } else {
                FilletSurfaceGenerator::circular_default()
            });
            fs.base.edgeuses = [eu_id(sel(eu[0], n_eu(&c))), eu_id(sel(eu[1], n_eu(&c)))];
            let pts: Vec<TsectPnt> = tsects
                .iter()
                .map(|t| {
                    let mut p = TsectPnt::default();
                    *p.uv_pos_mut(SrfIndex::First) = Point2::new(t.uv0[0], t.uv0[1]);
                    *p.uv_pos_mut(SrfIndex::Second) = Point2::new(t.uv1[0], t.uv1[1]);
                    p.curve_parameter = t.cparam;
                    p
                })
                .collect();
            let r = match which {
                0 => FilletSolver::cylinder_fillet_geom_parts(&mut fs, &k, &c.brep, &pts),
                1 => FilletSolver::cone_fillet_geom_parts(&mut fs, &k, &c.brep, &pts),
                _ => FilletSolver::torus_fillet_geom_parts(&mut fs, &k, &c.brep, &pts, None),
            };
            match r {
                Ok(parts) => {
                    o.push(0.0);
                    put_parts(o, *which, &parts);
                    Ok(())
                }
                Err(e) => Err(e),
            }
        }
        Inputs::Corner { cfg, index } => {
            let fs = solver_of(
                cfg,
                &Surfs {
                    base: [plane_z(), plane_z()],
                    has: [false, false],
                    off: [0.0, 0.0],
                    orient: [1.0, 1.0],
                },
            )?;
            let empty = Brep::new();
            let index = usize::try_from(*index).map_err(|_| harness("index"))?;
            match fs.base.fillet_corner(&empty, index) {
                Ok(c) => {
                    o.push(b(c.is_none()));
                    o.push(fs.jacobian_size().unwrap_or(0) as f64);
                    Ok(())
                }
                Err(e) => Err(e),
            }
        }
        Inputs::Record { cfg, brep, eu } => {
            let Ok(c) = brep_case(brep)? else {
                return Err(harness("B-rep failed on the Rust side only"));
            };
            o.push(0.0);
            let mut fs = solver_of(
                cfg,
                &Surfs {
                    base: [plane_z(), plane_z()],
                    has: [false, false],
                    off: [0.0, 0.0],
                    orient: [1.0, 1.0],
                },
            )?;
            fs.base.edgeuses = [eu_id(sel(*eu, n_eu(&c))), None];
            let pts = [TsectPnt::default()];
            match fs.base.record_fillet_error(&c.brep, 5, &pts, "difftest") {
                Ok(()) => {
                    o.push(-2.0);
                    Ok(())
                }
                Err(e) => Err(e),
            }
        }
    })
}

fn rust(i: &Inputs) -> Result<Output, HarnessError> {
    // the C++ side reports a B-rep that cannot be built
    if let Some(brep) = i.brep() {
        if let Err(out) = brep_case(brep)? {
            return Ok(out);
        }
    }
    let mut o = Vec::new();
    Ok(match rust_run(i, &mut o)? {
        Ok(()) => Output::new().values("flat", o),
        Err(e) => Output::unported(&e).unwrap_or_else(|| {
            o.push(status(&e));
            Output::new().values("flat", o)
        }),
    })
}

// ---- generators ------------------------------------------------------------------------

/// The plane `z = 0` (a filler base for solvers without surfaces).
fn plane_z() -> SrfIn {
    SrfIn::Plane {
        o: [0.0, 0.0, 0.0],
        n: [0.0, 0.0, 1.0],
    }
}

/// A bilinear patch `P(u, v) = f(u, v)` over `[-l, l]^2` (knots `-l -l l l`).
fn patch(l: f64, f: impl Fn(f64, f64) -> [f64; 3]) -> SrfIn {
    let net = [-l, l]
        .iter()
        .map(|&u| {
            [-l, l]
                .iter()
                .map(|&v| {
                    let [x, y, z] = f(u, v);
                    [x, y, z, NOW]
                })
                .collect()
        })
        .collect();
    SrfIn::Bspline(BssIn {
        oob: false,
        form: 0,
        net,
        p: 1,
        q: 1,
        u: vec![-l, -l, l, l],
        v: vec![-l, -l, l, l],
    })
}

const L: f64 = 10.0;

/// The floor `(u, v, 0)` (normal `+z`).
fn floor() -> SrfIn {
    patch(L, |u, v| [u, v, 0.0])
}

/// The wall `(v, 0, u)` (normal `+y`).
fn wall_y() -> SrfIn {
    patch(L, |u, v| [v, 0.0, u])
}

/// The wall `(x, u, v)` (normal `+x`).
fn wall_x(x: f64) -> SrfIn {
    patch(L, move |u, v| [x, u, v])
}

/// A straight degree-1 curve from `a` to `b`.
fn line(a: [f64; 3], bb: [f64; 3]) -> CrvIn {
    CrvIn {
        dim: 3,
        cpts: vec![[a[0], a[1], a[2], NOW], [bb[0], bb[1], bb[2], NOW]],
        p: 1,
        u: vec![0.0, 0.0, 1.0, 1.0],
    }
}

/// A rational quadratic arc of radius `r` about the origin in `z = 0`, from angle 0
/// through `quarters` quarter circles (1..=4).
fn arc(r: f64, quarters: usize) -> CrvIn {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let cs = [(1.0, 0.0), (0.0, 1.0), (-1.0, 0.0), (0.0, -1.0), (1.0, 0.0)];
    let mut cpts = Vec::new();
    let mut u = vec![0.0, 0.0, 0.0];
    for (q, pair) in cs.windows(2).take(quarters).enumerate() {
        let [(c0, s0), (c1, s1)] =
            [pair.first(), pair.last()].map(|p| p.copied().unwrap_or((1.0, 0.0)));
        if q == 0 {
            cpts.push([r * c0, r * s0, 0.0, 1.0]);
        }
        let (cx, cy) = (c0 + c1, s0 + s1);
        cpts.push([r * cx * w, r * cy * w, 0.0, w]);
        cpts.push([r * c1, r * s1, 0.0, 1.0]);
        let t = (q + 1) as f64 / quarters as f64;
        if q + 1 < quarters {
            u.extend([t, t]);
        }
    }
    u.extend([1.0, 1.0, 1.0]);
    CrvIn {
        dim: 3,
        cpts,
        p: 2,
        u,
    }
}

fn gen_cfg(rng: &mut Rng, extra_max: u64, radius: f64) -> Cfg {
    Cfg {
        tol: if rng.chance(0.8) {
            1.0e-6
        } else {
            rng.uniform(1e-8, 1e-3)
        },
        extra: rng.below(extra_max + 1) as usize,
        radius,
        slope: 0.0,
        target: rng.uniform(-3.0, 3.0),
        gu: 0.0,
        gv: 0.0,
    }
}

fn jitter(rng: &mut Rng, v: f64, s: f64) -> f64 {
    v + rng.uniform(-s, s)
}

/// The floor / `y = 0` wall surfaces with offsets `r`.
fn corner_surfs(r: f64) -> Surfs {
    Surfs {
        base: [floor(), wall_y()],
        has: [true, true],
        off: [r, r],
        orient: [1.0, 1.0],
    }
}

fn gen_plane(rng: &mut Rng) -> Inputs {
    let r = rng.uniform(0.2, 2.0);
    let x0 = rng.uniform(-5.0, 5.0);
    let mut surfs = corner_surfs(r);
    if rng.chance(0.2) {
        surfs.base = [gen_surface(rng), gen_surface(rng)];
    }
    let normal = if rng.chance(0.7) {
        [1.0, 0.0, 0.0]
    } else {
        [1.0, rng.uniform(-0.3, 0.3), rng.uniform(-0.3, 0.3)]
    };
    let s = if rng.chance(0.5) { 0.01 } else { 1.0 };
    let d = [-L, -L, L, L];
    Inputs::Plane {
        cfg: gen_cfg(rng, 4, r),
        surfs,
        orig: [x0, r, r],
        normal,
        dom1: d,
        dom2: d,
        uv1: [jitter(rng, x0, s), jitter(rng, r, s)],
        uv2: [jitter(rng, r, s), jitter(rng, x0, s)],
    }
}

fn gen_rail_eu(rng: &mut Rng) -> Inputs {
    let r = rng.uniform(0.2, 2.0);
    let mut surfs = corner_surfs(r);
    let brep = if rng.chance(0.15) {
        surfs.base = [gen_surface(rng), gen_surface(rng)];
        BrepIn::File {
            file: (*rng.choose(&FILES).unwrap_or(&"plane_square.smb")).to_string(),
            rebuild: false,
        }
    } else {
        let (a, bb) = (rng.uniform(-6.0, -1.0), rng.uniform(1.0, 6.0));
        let (c, d) = if rng.chance(0.8) {
            (rng.uniform(-3.0, r - 0.05), rng.uniform(r + 0.05, 5.0))
        } else {
            (rng.uniform(-3.0, 3.0), rng.uniform(3.5, 6.0))
        };
        // mostly a real SmPlane face (its UV trim curves come from the plane's own
        // projection), sometimes the floor patch
        let surface = if rng.chance(0.75) { plane_z() } else { floor() };
        surfs.base[0] = surface.clone();
        BrepIn::Face {
            surface,
            dom: [a, c, bb, d],
        }
    };
    if rng.chance(0.15) {
        if let Some(h) = surfs.has.get_mut(rng.below(2) as usize) {
            *h = false;
        }
    }
    let mut c = gen_cfg(rng, 5, r);
    if c.extra == 5 && rng.chance(0.7) {
        c.extra = rng.below(5) as usize;
    }
    Inputs::RailEu {
        cfg: c,
        surfs,
        brep,
        side: rng.below(16) as i64,
        frac: rng.uniform(0.0, 1.0),
        rail: usize::from(!rng.chance(0.8)),
        ptype: rng.below(8) as u8,
    }
}

fn gen_rail_rail(rng: &mut Rng) -> Inputs {
    let r = rng.uniform(0.2, 2.0);
    let mut surfs = [corner_surfs(r), corner_surfs(r)];
    surfs[1].base[1] = wall_x(rng.uniform(-3.0, 3.0));
    let mut any_missing = false;
    for s in &mut surfs {
        for h in &mut s.has {
            if rng.chance(0.3) {
                *h = false;
                any_missing = true;
            }
        }
    }
    if !any_missing && rng.chance(0.7) {
        let k = rng.below(2) as usize;
        let i = rng.below(2) as usize;
        if let Some(h) = surfs.get_mut(k).and_then(|s| s.has.get_mut(i)) {
            *h = false;
        }
    }
    Inputs::RailRail {
        cfg: [gen_cfg(rng, 2, r), gen_cfg(rng, 2, r)],
        surfs,
        uv: [rng.uniform(-3.0, 3.0), rng.uniform(-3.0, 3.0)],
        rail: rng.below(2) as usize,
        other: rng.below(2) as usize,
    }
}

fn tsect(uv0: [f64; 2], uv1: [f64; 2], cparam: f64) -> Tsect {
    Tsect { uv0, uv1, cparam }
}

/// A cylinder or cone case (`which` 0/1): the fillet between the floor and the `y = 0`
/// wall along the x axis.
fn gen_straight(rng: &mut Rng, which: u8) -> Inputs {
    let r = rng.uniform(0.2, 2.0);
    let (xs, xe) = (rng.uniform(-5.0, -1.0), rng.uniform(1.0, 5.0));
    let slope = if which == 1 || rng.chance(0.15) {
        if rng.chance(0.3) {
            0.0
        } else {
            rng.uniform(-0.5 * r, 0.5 * r)
        }
    } else {
        0.0
    };
    let r1 = r + slope;
    let mut surfs = corner_surfs(r);
    let mut ts = vec![tsect([xs, r], [r, xs], 0.0), tsect([xe, r1], [r1, xe], 1.0)];
    if rng.chance(0.2) {
        // the surfaces swapped: the arc runs the other way
        surfs.base = [wall_y(), floor()];
        for t in &mut ts {
            std::mem::swap(&mut t.uv0, &mut t.uv1);
        }
    }
    if rng.chance(0.2) {
        for t in &mut ts {
            for v in t.uv0.iter_mut().chain(t.uv1.iter_mut()) {
                *v = jitter(rng, *v, 0.05);
            }
        }
    }
    let curve = if rng.chance(0.1) {
        arc(2.0, 1)
    } else {
        line([xs - 0.5, 0.0, 0.0], [xe + 0.5, 0.0, 0.0])
    };
    let mut c = gen_cfg(rng, 0, r);
    c.slope = slope;
    Inputs::Geom {
        which,
        cfg: c,
        surfs,
        fsg: u8::from(!rng.chance(0.1)),
        before: rng.chance(0.5),
        after: rng.chance(0.5),
        factor: if rng.chance(0.7) {
            0.5
        } else {
            rng.uniform(0.0, 1.0)
        },
        brep: BrepIn::Wire {
            curve,
            t0: 0.0,
            t1: 1.0,
        },
        eu: [rng.below(2) as i64, rng.below(2) as i64],
        tsects: ts,
    }
}

/// A torus case: the inner fillet between the floor and the wall tangent to an arc edge
/// of radius `big` at angle 0.
fn gen_torus(rng: &mut Rng) -> Inputs {
    let big = rng.uniform(2.0, 5.0);
    let r = rng.uniform(0.2, 1.0);
    let rho = big - r;
    let quarters = 1 + rng.below(4) as usize;
    let span = quarters as f64 * std::f64::consts::FRAC_PI_2;
    let phi = span
        * if rng.chance(0.5) {
            1.0
        } else {
            rng.uniform(0.3, 0.95)
        };
    let slope = if rng.chance(0.1) { 0.2 * r } else { 0.0 };
    let surfs = Surfs {
        base: [floor(), wall_x(big)],
        has: [true, true],
        off: [r, -r],
        orient: [1.0, -1.0],
    };
    let ts = vec![
        tsect([rho, 0.0], [0.0, r], 0.0),
        tsect([rho * phi.cos(), rho * phi.sin()], [0.0, r], 1.0),
    ];
    let curve = if rng.chance(0.1) {
        line([0.0, 0.0, 0.0], [big, 0.0, 0.0])
    } else {
        arc(big, quarters)
    };
    let mut c = gen_cfg(rng, 0, r);
    c.slope = slope;
    Inputs::Geom {
        which: 2,
        cfg: c,
        surfs,
        fsg: u8::from(!rng.chance(0.1)),
        before: rng.chance(0.5),
        after: rng.chance(0.85),
        factor: 0.5,
        brep: BrepIn::Wire {
            curve,
            t0: 0.0,
            t1: 1.0,
        },
        eu: if rng.chance(0.5) { [0, 1] } else { [1, 0] },
        tsects: ts,
    }
}

fn gen_case(rng: &mut Rng) -> Inputs {
    match rng.below(20) {
        0..=3 => gen_plane(rng),
        4..=8 => gen_rail_eu(rng),
        9 | 10 => gen_rail_rail(rng),
        11..=13 => gen_straight(rng, 0),
        14 | 15 => gen_straight(rng, 1),
        16 | 17 => gen_torus(rng),
        18 => Inputs::Corner {
            cfg: gen_cfg(rng, 2, 1.0),
            index: 2 + rng.below(50),
        },
        _ => Inputs::Record {
            cfg: gen_cfg(rng, 0, 1.0),
            brep: if rng.chance(0.5) {
                BrepIn::Wire {
                    curve: line([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]),
                    t0: 0.0,
                    t1: 1.0,
                }
            } else {
                BrepIn::Face {
                    surface: floor(),
                    dom: [-1.0, -1.0, 1.0, 1.0],
                }
            },
            eu: rng.below(8) as i64,
        },
    }
}

impl UnitDiff for Runner {
    fn id(&self) -> &'static str {
        "blend.fillet_solver_p2"
    }

    fn generate(&self, rng: &mut Rng) -> Option<Value> {
        serde_json::to_value(gen_case(rng)).ok()
    }

    fn cases(&self, ctx: &crate::case::CaseContext) -> Result<Vec<crate::Case>, HarnessError> {
        let ctx = crate::case::CaseContext {
            seeds: ctx.seeds.max(300),
            ..ctx.clone()
        };
        crate::runner::default_cases(self, &ctx)
    }

    fn run_cpp(&self, case: &crate::Case) -> Result<Output, HarnessError> {
        cpp(&case.parse()?)
    }

    fn run_rust(&self, case: &crate::Case) -> Result<Output, HarnessError> {
        rust(&case.parse()?)
    }
}
