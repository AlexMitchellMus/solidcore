//! Diagnostic-only native test for review_cases.json#72, outside the corpus.
//!
//! Tags 1..=4 copy actual call arguments/results without extra evaluations.
//! Tags >=11 reconstruct the span AFTER the genuine operation has returned;
//! they must never be described as observations of its internal Hermite objects.
//! Per-side bit identity against the untouched runner is the observer safety gate.

use super::*;
use solidcore::blend::fillet_custom_solver::{InitialValues, JacobianLoad};
use solidcore::blend::fillet_intersector::{
    FilStatus, FilletSolverOps as IntersectorOps, HandlingAlgorithmType, StepSolve,
};
use solidcore::blend::fillet_solver::p1::{
    FilletSolver, FilletSolverOps as SolverOps, FilletSolverType,
};
use solidcore::blend::fillet_solver::p3::StepSolveOutcome;
use solidcore::curve::crv_on_surf::CrvOnSurf;
use solidcore::curve::hermite_curve::HermiteCurve;
use solidcore::geom::extent1d::ExtentNd;
use solidcore::geom::matrix::Matrix;
use solidcore::intersect::surface_intersector::p1::SrfIndex;
use solidcore::nurbs::b_spline_curve::BSplineCurve;
use solidcore::topology::Brep;
use solidcore::Result as KernelResult;
use std::io::Write;

#[derive(Debug)]
struct Event {
    tag: u32,
    values: Vec<f64>,
}

struct Observed {
    inner: ConstantRadiusAssistedFs,
    events: Vec<Event>,
    step: Option<f64>,
}

// All unobserved members delegate to the existing concrete solver. In particular,
// do NOT replace SolverOps::point_on_plane_solve with the similarly named
// IntersectorOps member: StepSolve currently dispatches through SolverOps.
macro_rules! delegate_mut {
    ($trait:path; $(fn $name:ident($($arg:ident: $ty:ty),*) -> $result:ty;)*) => {
        $(fn $name(&mut self, $($arg: $ty),*) -> $result {
            <ConstantRadiusAssistedFs as $trait>::$name(&mut self.inner, $($arg),*)
        })*
    };
}
macro_rules! delegate_ref {
    ($trait:path; $(fn $name:ident($($arg:ident: $ty:ty),*) -> $result:ty;)*) => {
        $(fn $name(&self, $($arg: $ty),*) -> $result {
            <ConstantRadiusAssistedFs as $trait>::$name(&self.inner, $($arg),*)
        })*
    };
}

fn result_status<T>(r: &KernelResult<T>) -> f64 {
    r.as_ref().err().map_or(0.0, status)
}

impl SolverOps for Observed {
    delegate_ref! { SolverOps;
        fn solver() -> &FilletSolver;
        fn jacobian_size() -> KernelResult<usize>;
        fn offset_radii_can_change() -> KernelResult<bool>;
        fn setup_offset_extension(k: &Kernel, surface: &Surface, domain: &mut Extent2d) -> KernelResult<()>;
        fn can_make_analytic() -> KernelResult<bool>;
        fn type_code() -> KernelResult<i64>;
        fn is_kind_of(t: i64) -> KernelResult<bool>;
        fn solver_type() -> KernelResult<FilletSolverType>;
        fn constant_distance() -> KernelResult<f64>;
    }
    delegate_mut! { SolverOps;
        fn solver_mut() -> &mut FilletSolver;
        fn setup_offset_values(point: &TsectPnt, radius: Option<&mut f64>) -> KernelResult<()>;
        fn set_extension_flags(brep: &Brep) -> KernelResult<()>;
        fn load_jacobian(k: &Kernel, equations: &mut usize, parameters: &mut usize, load: &mut JacobianLoad<'_>, f: &mut [f64], jacobian: Option<&mut Matrix>) -> KernelResult<bool>;
        fn load_initial_values(k: &Kernel, brep: Option<&Brep>, do_surf2_calcs: bool, rail_index: usize, surface1: &mut OffsetSurface, surface1_index: usize, intervals: &mut ExtentNd, periodicities: &mut Vec<bool>, guess_t: &mut Vec<f64>) -> KernelResult<Option<InitialValues>>;
    }

    fn point_on_plane_solve(
        &mut self,
        k: &Kernel,
        origin: Vec3,
        normal: Vec3,
        d0: &Extent2d,
        d1: &Extent2d,
        uv0: Vec2,
        uv1: Vec2,
        point: &mut TsectPnt,
    ) -> KernelResult<bool> {
        let mut before = Vec::new();
        put_vec(&mut before, origin);
        put_vec(&mut before, normal);
        put_uv(&mut before, uv0);
        put_uv(&mut before, uv1);
        for slot in &self.inner.solver.surfaces {
            before.push(slot.as_ref().unwrap().surface.offset_distance());
        }
        self.events.push(Event {
            tag: 1,
            values: before,
        });
        let result = SolverOps::point_on_plane_solve(
            &mut self.inner,
            k,
            origin,
            normal,
            d0,
            d1,
            uv0,
            uv1,
            point,
        );
        let mut after = vec![result_status(&result), flag(*result.as_ref().unwrap())];
        put_point(&mut after, point);
        self.events.push(Event {
            tag: 2,
            values: after,
        });
        result
    }
}

impl IntersectorOps for Observed {
    type FilletGeom = FilletGeomId;

    delegate_ref! { IntersectorOps;
        fn last_fillet_geom() -> KernelResult<FilletGeomId>;
        fn self_intersection_handler() -> KernelResult<Option<HandlingAlgorithmType>>;
        fn has_surface_generator() -> KernelResult<bool>;
    }
    delegate_mut! { IntersectorOps;
        fn compute_surface_values(k: &Kernel, uv: [Vec2; 2], point: &mut TsectPnt) -> KernelResult<()>;
        fn compute_uv_vectors(k: &Kernel, point: &mut TsectPnt) -> KernelResult<()>;
        fn point_on_plane_solve(k: &Kernel, origin: Vec3, normal: Vec3, d0: &Extent2d, d1: &Extent2d, uv0: Vec2, uv1: Vec2, point: &mut TsectPnt) -> KernelResult<bool>;
        fn set_surfaces_to_zero_offset() -> KernelResult<()>;
        fn reload_surface_offsets() -> KernelResult<()>;
        fn rail_surface_evaluate_1st_derivatives(k: &Kernel, rail: usize, uv: Vec2) -> KernelResult<(Vec3, Vec3, Vec3)>;
        fn set_status(value: FilStatus) -> KernelResult<()>;
        fn smooth_out_intersection(geom: &FilletGeomId, curves: &mut Vec<BSplineCurve>, uv0: &mut Vec<BSplineCurve>, uv1: &mut Vec<BSplineCurve>) -> KernelResult<()>;
        fn create_fillet_surface(k: &Kernel, center: &mut BSplineCurve, curve0: &mut BSplineCurve, curve1: &mut BSplineCurve, points: &[TsectPnt], uv0: &mut BSplineCurve, uv1: &mut BSplineCurve) -> KernelResult<Option<BSplineSurface>>;
        fn set_fillet_surface(pseudo: &mut FilletBrep, geom: &FilletGeomId, surface: Option<BSplineSurface>) -> KernelResult<()>;
        fn set_center_line_curve(pseudo: &mut FilletBrep, geom: &FilletGeomId, curve: BSplineCurve) -> KernelResult<()>;
        fn set_rail_curve(pseudo: &mut FilletBrep, geom: &FilletGeomId, rail: usize, curve: BSplineCurve) -> KernelResult<()>;
        fn set_rail_mate_uv_curve(pseudo: &mut FilletBrep, geom: &FilletGeomId, rail: usize, curve: BSplineCurve) -> KernelResult<()>;
    }

    fn step_solve(
        &mut self,
        k: &Kernel,
        current: &TsectPnt,
        step: &mut f64,
        d0: &Extent2d,
        d1: &Extent2d,
        flags: &mut StepSolve,
        next: &mut TsectPnt,
    ) -> KernelResult<()> {
        let mut out = StepSolveOutcome::default();
        let result = FilletSolver::step_solve(self, k, None, current, step, d0, d1, next, &mut out);
        flags.found_answer = out.found_good_point;
        flags.clipped = out.clipped;
        flags.boundary_hit = out.boundary_hit;
        self.step = Some(*step);
        let mut values = vec![*step];
        put_point(&mut values, next);
        self.events.push(Event { tag: 3, values });
        result
    }

    fn compute_point_values(
        &mut self,
        k: &Kernel,
        uv: [Vec2; 2],
        d0: &Extent2d,
        d1: &Extent2d,
        dir: f64,
        point: &mut TsectPnt,
        previous: Option<&TsectPnt>,
    ) -> KernelResult<()> {
        let result = IntersectorOps::compute_point_values(
            &mut self.inner,
            k,
            uv,
            d0,
            d1,
            dir,
            point,
            previous,
        );
        let mut values = vec![result_status(&result)];
        put_point(&mut values, point);
        self.events.push(Event { tag: 4, values });
        result
    }
}

fn reconstruct_span(
    k: &Kernel,
    surfaces: [&Surface; 2],
    domains: [Extent2d; 2],
    current: &TsectPnt,
    next: &TsectPnt,
    step: f64,
) -> Vec<Event> {
    let indices = [SrfIndex::First, SrfIndex::Second];
    let uv = indices.map(|s| {
        HermiteCurve::from_2d(
            current.uv_pos(s),
            current.uv_deriv(s) * step,
            next.uv_pos(s),
            next.uv_deriv(s) * step,
            2,
        )
    });
    let mut events = Vec::new();
    let rails = std::array::from_fn::<_, 2, _>(|i| {
        let map = CrvOnSurf::new(
            k,
            Curve::from(uv[i]),
            surfaces[i].clone_for_evaluation().unwrap(),
            Some(domains[i]),
            0,
        )
        .unwrap();
        let interval = map.natural_interval();
        let p0 = map.evaluate(k, interval.min(), 1, true, true).unwrap();
        let p1 = map.evaluate(k, interval.max(), 1, true, true).unwrap();
        let mut values = Vec::new();
        for p in [p0[0], p0[1], p1[0], p1[1]] {
            put_vec(&mut values, p);
        }
        events.push(Event {
            tag: 11 + i as u32,
            values,
        });
        HermiteCurve::new(p0[0], p0[1], p1[0], p1[1], 3)
    });
    let mut maximum = 0.0;
    for (i, t) in [0.5, 0.17, 0.83, 0.33, 0.67].into_iter().enumerate() {
        let raw = uv.map(|c| c.evaluate_point(t));
        let clamped = std::array::from_fn::<_, 2, _>(|j| {
            domains[j].clamp_point(Vec2::new(raw[j].x, raw[j].y))
        });
        let points =
            std::array::from_fn::<_, 2, _>(|j| surfaces[j].evaluate_point(k, clamped[j]).unwrap());
        let curve_points = rails.map(|c| c.evaluate_point(t));
        let ds0 = points[0].distance_squared(curve_points[0]);
        let ds1 = points[1].distance_squared(curve_points[1]);
        let distance = solidcore::math::max2(ds0, ds1);
        if distance > maximum {
            maximum = distance;
        }
        let mut values = vec![t];
        for p in raw {
            put_vec(&mut values, p);
        }
        for p in clamped {
            put_uv(&mut values, p);
        }
        for p in points {
            put_vec(&mut values, p);
        }
        for p in curve_points {
            put_vec(&mut values, p);
        }
        for (point, curve_point) in points.into_iter().zip(curve_points) {
            put_vec(&mut values, point - curve_point);
        }
        values.extend([ds0, ds1, distance, maximum]);
        events.push(Event {
            tag: 20 + i as u32,
            values,
        });
    }
    events.push(Event {
        tag: 30,
        values: vec![solidcore::geom::vector3d::sqrt(maximum)],
    });
    events
}

fn observed_rust(i: &Inputs) -> (Vec<f64>, Vec<Event>) {
    // Same construction as rust_record, kept here so the ordinary runner does
    // not acquire an instrumentation branch, callback, new record or state.
    let k = Kernel::default();
    let a = rust_surface(&i.a).unwrap();
    let b = i.b.as_ref().map(|s| rust_surface(s).unwrap());
    let bref = b.as_ref().unwrap_or(&a);
    let s = &i.solver;
    let orient = s.orient.map(|f| if f { 1.0 } else { -1.0 });
    let off0 = rail_offset(&k, s.radii[0] * orient[0], a.clone(), SurfaceId::Free(0))
        .unwrap_or_else(|_| panic!("first fixture offset construction failed"));
    let off1 = rail_offset(
        &k,
        s.radii[1] * orient[1],
        bref.clone(),
        SurfaceId::Free(u32::from(b.is_some())),
    )
    .unwrap_or_else(|_| panic!("second fixture offset construction failed"));
    let mut inner = ConstantRadiusAssistedFs::from_parts(
        s.tol,
        s.ang_tol,
        s.tangency_tol,
        s.radii,
        [Some(off0), Some(off1)],
        orient,
        [None, None],
    );
    for slot in inner.solver.surfaces.iter_mut().flatten() {
        slot.surface.surface.set_allow_extension(true);
    }
    let mut pseudo = FilletBrep::new();
    inner
        .solver
        .add_constructor_fillet_geom(&k, &mut pseudo, FilletSolverId(0))
        .unwrap();
    let mut solver = Observed {
        inner,
        events: Vec::new(),
        step: None,
    };
    let mut fi: Fi<'_> =
        FilletIntersector::new(&k, &a, &ext(i.dom1), bref, &ext(i.dom2), None).unwrap();
    fi.base.base.this_ang_tol_rad = i.ang_tol;
    fi.base.base.this_approx_tol_3d = i.approx_tol;
    fi.base.base.clipped_by_through_point = i.clipped;
    fi.base.base.curve_trace_direction = i.trace_dir;
    for p in &i.points {
        fi.base
            .base
            .alloc_point(rust_point(&point_values(p)).unwrap());
    }
    let OpIn::Next {
        a: from,
        step,
        b: to,
    } = i.op
    else {
        panic!("fixture must be Next")
    };
    assert!(
        !i.clipped,
        "reconstruction does not model finite-difference adjustment"
    );
    let out = fi
        .compute_next_point(&mut solver, &k, id(from), step, id(to))
        .unwrap();
    assert!(out.found_good_point && out.clipped && !out.boundary_hit);
    let mut record = vec![
        2.0,
        0.0,
        flag(out.found_good_point),
        out.deviation_found.unwrap(),
        out.angle_found_rad.unwrap(),
        flag(out.clipped),
        flag(out.boundary_hit),
    ];
    let si = &fi.base.base;
    put_dom(&mut record, &si.uv_domains[0]);
    put_dom(&mut record, &si.uv_domains[1]);
    record.extend([
        si.this_ang_tol_rad,
        si.this_approx_tol_3d,
        flag(si.do_boundary_point_on_curve_test),
        flag(fi.curr_fillet_touch_boundary()),
        flag(fi.self_intersect()),
        flag(si.clipped_by_through_point),
        si.curve_trace_direction,
        flag(fi.is_kind_of(TYPE)),
        flag(fi.is_kind_of(ADV_TYPE)),
        flag(fi.is_kind_of(7)),
    ]);
    for slot in &solver.inner.solver.surfaces {
        record.push(
            slot.as_ref()
                .map_or(ffi::UNSET, |s| s.surface.offset_distance()),
        );
    }
    record.extend([si.curve_points().len() as f64, si.points().len() as f64]);
    for p in si.points() {
        put_point(&mut record, p);
    }
    // Freeze the genuine operation's record before any extra evaluations.
    let reconstruction = reconstruct_span(
        &k,
        [&a, bref],
        si.uv_domains,
        si.point(id(from)).unwrap(),
        si.point(id(to)).unwrap(),
        solver.step.unwrap(),
    );
    solver.events.extend(reconstruction);
    (record, solver.events)
}

fn unpack_cpp(raw: &[f64]) -> (Vec<f64>, Vec<Event>) {
    let count = raw[0] as usize;
    assert_eq!(count, 189);
    let record = raw[1..=count].to_vec();
    let mut rest = &raw[count + 1..];
    let mut events = Vec::new();
    while !rest.is_empty() {
        assert!(rest.len() >= 2);
        let tag = rest[0] as u32;
        let len = rest[1] as usize;
        assert!(rest.len() >= len + 2);
        events.push(Event {
            tag,
            values: rest[2..len + 2].to_vec(),
        });
        rest = &rest[len + 2..];
    }
    (record, events)
}

fn bits(v: &[f64]) -> Vec<u64> {
    v.iter().map(|x| x.to_bits()).collect()
}

fn label(tag: u32) -> &'static str {
    match tag {
        1 => "observed.plane_input.origin_normal_uv_offsets",
        2 => "observed.plane_output.status_found_point80",
        3 => "observed.step_output.rescaled_step_point80",
        4 => "observed.point_values_output.status_point80",
        11 | 12 => "reconstruction.rail_endpoints.p0_d0_p1_d1",
        20..=24 => "reconstruction.sample.t_rawuv_clampeduv_surface_rail_delta_dist2_max",
        30 => "reconstruction.span_deviation",
        _ => panic!("unknown checkpoint tag {tag}"),
    }
}

fn difference(cpp: &Event, rust: &Event) -> Option<Value> {
    cpp.values
        .iter()
        .zip(&rust.values)
        .position(|(c, r)| c.to_bits() != r.to_bits())
        .map(|component| {
            json!({
                "tag": cpp.tag, "checkpoint": label(cpp.tag), "component": component,
                "cpp": cpp.values[component], "rust": rust.values[component],
                "cpp_bits": cpp.values[component].to_bits(),
                "rust_bits": rust.values[component].to_bits(),
            })
        })
}

#[test]
fn original_clipped_step_checkpoint_evidence() {
    let fixture: Value = serde_json::from_str(include_str!(
        "../../../../../fixtures/blend.fillet_intersector/review_cases.json"
    ))
    .unwrap();
    let case = &fixture["cases"][72];
    assert_eq!(case["name"], "rv_next_clipped_step");
    let input: Inputs = serde_json::from_value(case["inputs"].clone()).unwrap();
    let plain_cpp = ffi::run(&ffi_program(&input)).unwrap();
    let plain_rust =
        rust_record(&input).unwrap_or_else(|_| panic!("original Rust operation failed"));
    let (cpp_record, cpp_events) = unpack_cpp(&ffi::probe_next(&ffi_program(&input)).unwrap());
    let (rust_record, rust_events) = observed_rust(&input);
    assert_eq!(plain_cpp.len(), 189);
    assert_eq!(plain_rust.len(), 189);
    assert_eq!(
        bits(&plain_cpp),
        bits(&cpp_record),
        "C++ observation changed original record"
    );
    assert_eq!(
        bits(&plain_rust),
        bits(&rust_record),
        "Rust observation changed original record"
    );
    let expected = [
        (1, 12),
        (2, 82),
        (3, 81),
        (4, 81),
        (11, 12),
        (12, 12),
        (20, 33),
        (21, 33),
        (22, 33),
        (23, 33),
        (24, 33),
        (30, 1),
    ];
    for events in [&cpp_events, &rust_events] {
        assert_eq!(
            events
                .iter()
                .map(|e| (e.tag, e.values.len()))
                .collect::<Vec<_>>(),
            expected
        );
        assert!(events.iter().flat_map(|e| &e.values).all(|v| v.is_finite()));
    }
    let first_observed = cpp_events
        .iter()
        .zip(&rust_events)
        .filter(|(c, _)| c.tag <= 4)
        .find_map(|(c, r)| difference(c, r));
    let first_reconstructed = cpp_events
        .iter()
        .zip(&rust_events)
        .filter(|(c, _)| c.tag >= 11)
        .find_map(|(c, r)| difference(c, r));
    // Direct stdout intentionally survives libtest's passing-test capture, so
    // canonical CI's existing native-feature test job retains the evidence. No
    // new tracing option, environment knob, file or corpus comparison is added.
    let mut stdout = std::io::stdout().lock();
    writeln!(stdout, "FILLET_SPAN_PROBE {}", json!({
        "case": "review_cases.json#72 rv_next_clipped_step",
        "original_values": 189,
        "observer_records_bit_identical": true,
        "original_deviation": {"cpp": plain_cpp[3], "rust": plain_rust[3]},
        "first_observed_checkpoint_difference_not_internal_root_cause": first_observed,
        "first_post_operation_reconstruction_difference": first_reconstructed,
        "reconstruction_reproduces_original_deviation_bits": {
            "cpp": cpp_events.last().unwrap().values[0].to_bits() == plain_cpp[3].to_bits(),
            "rust": rust_events.last().unwrap().values[0].to_bits() == plain_rust[3].to_bits(),
        },
        "limitation": "No hook observes the internal LocalCurveIntersect or span locals; replay is explicitly reconstruction, never proof of their original cache state."
    })).unwrap();
    for (c, r) in cpp_events.iter().zip(&rust_events) {
        writeln!(
            stdout,
            "FILLET_SPAN_PROBE {}",
            json!({
                "tag": c.tag, "checkpoint": label(c.tag),
                "cpp": c.values, "rust": r.values,
                "cpp_bits": bits(&c.values), "rust_bits": bits(&r.values),
            })
        )
        .unwrap();
    }
}
