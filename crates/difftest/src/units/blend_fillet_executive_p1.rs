//! Runner for unit `blend.fillet_executive_p1` (the fillet laws). C++:
//! `solidcore_ffi::units::blend_fillet_executive_p1`. Rust:
//! `solidcore::blend::fillet_executive::p1`.
//!
//! One op, `law_eval`: build a law (base, linear or B-spline, with or without an edge
//! map), evaluate it at a sequence of parameters over edge intervals (forward or
//! reversed), then dump the extension the B-spline law cached and `IsKindOf` answers.
use serde::{Deserialize, Serialize};
use serde_json::Value;

use solidcore::blend::fillet_executive::p1::{
    FilletLaw, FilletLawKind, BSPLINE_FILLET_LAW_TYPE, FILLET_LAW_TYPE, LINEAR_FILLET_LAW_TYPE,
};
use solidcore::geom::extent1d::Extent1d;
use solidcore::geom::object::OBJECT_TYPE;
use solidcore::math::EFF_ZERO;
use solidcore::nurbs::math::NOW;
use solidcore_ffi::units::blend_fillet_executive_p1 as ffi;

use super::blend_fillet_standard_solver_p1::{
    bcurve_of, dump_curve, push_crv, status, uniform_knots, InlineCrv,
};
use crate::util::flag;
use crate::{Case, HarnessError, Output, Rng, UnitDiff};

/// `blend.fillet_executive_p1` runner.
pub struct Runner;

/// The law to build.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
enum LawIn {
    /// `SmFilletLaw()`.
    Base,
    /// `SmLinearFilletLaw(r0, r1, map)`.
    Linear { r0: f64, r1: f64 },
    /// `SmBSplineFilletLaw(curve, map)` (`[0, 0]`, i.e. none, without a map).
    Bspline { c: InlineCrv },
}

/// One `Evaluate(t, ivl, reverse)` call.
#[derive(Debug, Clone, Serialize, Deserialize)]
struct Call {
    t: f64,
    ivl: [f64; 2],
    reverse: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct Inputs {
    law: LawIn,
    /// The edge map (`pOptEdgeMap` / `crEdgeMap`); `None` = not given.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    map: Option<[f64; 2]>,
    calls: Vec<Call>,
    /// Type codes for `IsKindOf`.
    kinds: Vec<i64>,
}

fn cpp_input(i: &Inputs) -> Vec<f64> {
    let mut o = Vec::new();
    match &i.law {
        LawIn::Linear { r0, r1 } => o.extend([0.0, *r0, *r1]),
        LawIn::Bspline { c } => {
            o.push(1.0);
            push_crv(&mut o, c);
        }
        LawIn::Base => o.push(2.0),
    }
    match i.map {
        Some([a, b]) => o.extend([1.0, a, b]),
        None => o.extend([0.0, 0.0, 0.0]),
    }
    o.push(i.calls.len() as f64);
    for c in &i.calls {
        o.extend([c.t, c.ivl[0], c.ivl[1], flag(c.reverse)]);
    }
    o.push(i.kinds.len() as f64);
    o.extend(i.kinds.iter().map(|&k| k as f64));
    o
}

fn to_output(flat: &[f64]) -> Output {
    let (status, rest) = flat
        .split_first()
        .map_or((f64::NAN, &[][..]), |(s, r)| (*s, r));
    Output::new()
        .values("status", vec![status])
        .values("out", rest.to_vec())
}

/// The law the shim builds. The `SmLinearFilletLaw(r0, r1, pOptEdgeMap)` and
/// `SmBSplineFilletLaw(pCurve, crEdgeMap)` constructors belong to unit
/// `blend.fillet_executive_p2` (PR #267); until they land the harness mirrors their one
/// rule: an edge map shorter than `SM_EFF_ZERO` (or none) means `[0, 1]`.
fn law_of(i: &Inputs) -> Result<FilletLaw, HarnessError> {
    let kind = match &i.law {
        LawIn::Base => FilletLawKind::Base,
        LawIn::Linear { r0, r1 } => FilletLawKind::Linear {
            start_radius: *r0,
            end_radius: *r1,
        },
        LawIn::Bspline { c } => FilletLawKind::BSpline {
            law_curve: bcurve_of(c)?,
            extended_law: None,
        },
    };
    let mut law = FilletLaw::new(kind);
    if let (Some([a, b]), false) = (i.map, matches!(law.kind, FilletLawKind::Base)) {
        let map = Extent1d::from_min_max(a, b);
        if map.length() >= EFF_ZERO {
            law.edge_map = map;
        }
    }
    Ok(law)
}

fn rust(i: &Inputs) -> Result<Output, HarnessError> {
    let mut law = law_of(i)?;
    let mut o = vec![
        law.edge_map.min(),
        law.edge_map.max(),
        flag(law.is_linear_law().map_err(|e| HarnessError::new(e.to_string()))?),
        law.type_code() as f64,
    ];
    for c in &i.calls {
        let ivl = Extent1d::from_min_max(c.ivl[0], c.ivl[1]);
        match law.evaluate(c.t, &ivl, c.reverse) {
            Ok(v) => {
                o.push(0.0);
                o.extend(v);
            }
            Err(e) => {
                if let Some(out) = Output::unported(&e) {
                    return Ok(out);
                }
                o.push(status(&e));
            }
        }
    }
    dump_curve(&mut o, law.extended_curve());
    o.extend(i.kinds.iter().map(|&k| flag(law.is_kind_of(k))));
    Ok(Output::new().values("status", vec![0.0]).values("out", o))
}

// ---- generator ------------------------------------------------------------------------

/// A random clamped law curve: radii as x over `[a, b]`, degree 1..=5, sometimes rational
/// or with repeated interior knots.
fn gen_law_curve(rng: &mut Rng, scale: f64, a: f64, b: f64) -> InlineCrv {
    let p = 1 + rng.below(5) as usize;
    let n = p + rng.below(6) as usize;
    let rational = rng.chance(0.2);
    let cpts = (0..=n)
        .map(|_| {
            let w = if rational {
                rng.uniform(0.5, 2.0)
            } else {
                NOW
            };
            let x = scale * rng.uniform(0.2, 2.0);
            // the law reads x only; stray y/z must not matter
            let y = if rng.chance(0.3) { rng.uniform(-1.0, 1.0) } else { 0.0 };
            [x, y, 0.0, w]
        })
        .collect();
    let mut u = uniform_knots(n, p, a, b);
    if rng.chance(0.3) && n > p && p >= 2 {
        // repeat one interior knot (multiplicity 2 <= p)
        let i = p + 1 + rng.below((n - p) as u64) as usize;
        if let Some(&k) = u.get(i) {
            if let Some(next) = u.get_mut(i + 1) {
                if *next != b {
                    *next = k;
                }
            }
        }
    }
    InlineCrv { cpts, p, u }
}

fn gen_scale(rng: &mut Rng) -> f64 {
    match rng.below(5) {
        0 => 10f64.powf(rng.uniform(-6.0, -2.0)),
        1 => 10f64.powf(rng.uniform(2.0, 6.0)),
        _ => 10f64.powf(rng.uniform(-1.0, 1.0)),
    }
}

fn gen_calls(rng: &mut Rng, domain: [f64; 2]) -> Vec<Call> {
    let n = 1 + rng.below(5) as usize;
    let len = domain[1] - domain[0];
    (0..n)
        .map(|_| {
            let (i0, i1) = if rng.chance(0.08) {
                let v = rng.uniform(-3.0, 3.0) * len.max(1e-3);
                (v, v)
            } else {
                let i0 = rng.uniform(-2.0, 2.0) * len.max(1e-3);
                (i0, i0 + rng.uniform(0.1, 3.0) * len.max(1e-3))
            };
            let t = match rng.below(10) {
                0..=5 => i0 + rng.unit() * (i1 - i0),
                6 => i0,
                7 => i1,
                // beyond the interval: the law is extended (or snaps when within tolerance)
                8 => i1 + rng.uniform(0.0, 3.0) * (i1 - i0),
                _ => i0 - rng.uniform(0.0, 3.0) * (i1 - i0),
            };
            Call {
                t,
                ivl: [i0, i1],
                reverse: rng.chance(0.4),
            }
        })
        .collect()
}

impl UnitDiff for Runner {
    fn id(&self) -> &'static str {
        "blend.fillet_executive_p1"
    }

    fn run_cpp(&self, case: &Case) -> Result<Output, HarnessError> {
        let i: Inputs = case.parse()?;
        let flat = ffi::law_eval(&cpp_input(&i)).map_err(|e| HarnessError::new(e.to_string()))?;
        Ok(to_output(&flat))
    }

    fn run_rust(&self, case: &Case) -> Result<Output, HarnessError> {
        let i: Inputs = case.parse()?;
        rust(&i)
    }

    fn generate(&self, rng: &mut Rng) -> Option<Value> {
        let scale = gen_scale(rng);
        let a = rng.uniform(-5.0, 5.0) * 10f64.powf(rng.uniform(-2.0, 2.0));
        let b = a + rng.uniform(0.1, 10.0) * 10f64.powf(rng.uniform(-2.0, 2.0));
        let (law, domain) = match rng.below(10) {
            0 => (LawIn::Base, [0.0, 1.0]),
            1..=4 => (
                LawIn::Linear {
                    r0: scale * rng.uniform(0.0, 2.0),
                    r1: scale * rng.uniform(0.0, 2.0),
                },
                [0.0, 1.0],
            ),
            _ => (
                LawIn::Bspline {
                    c: gen_law_curve(rng, scale, a, b),
                },
                [a, b],
            ),
        };
        let map = match rng.below(10) {
            0..=2 => None,
            3..=4 => Some(domain),
            5 => {
                // zero length: the constructors fall back to [0, 1]
                let v = domain[0] + rng.unit() * (domain[1] - domain[0]);
                Some([v, v])
            }
            _ => {
                let len = domain[1] - domain[0];
                let m0 = domain[0] + rng.uniform(-0.5, 0.8) * len;
                Some([m0, m0 + rng.uniform(0.05, 1.5) * len])
            }
        };
        let edge0 = rng.uniform(-10.0, 10.0);
        let edge1 = edge0 + rng.uniform(0.1, 20.0);
        let calls = gen_calls(rng, [edge0, edge1]);
        let mut kinds = vec![
            FILLET_LAW_TYPE,
            LINEAR_FILLET_LAW_TYPE,
            BSPLINE_FILLET_LAW_TYPE,
            OBJECT_TYPE,
        ];
        kinds.push(rng.range_i64(0, 20000));
        serde_json::to_value(Inputs {
            law,
            map,
            calls,
            kinds,
        })
        .ok()
    }
}
