//! Runner for unit `query.brep`.
//!
//! C++: `solidcore_ffi::units::query_brep` (the public `SM_API` entry points). Rust:
//! `solidcore::query::brep`.
//!
//! An `evaluate` case is an entry point (`point`, `normal` or `derivatives`), an optional
//! surface (`null` for a `NULL` surface), a parameter `uv` (entries may be `"NaN"`, `"inf"`
//! or `"-inf"`) and the values the outputs hold before the call. Both sides produce the
//! status (`0` = `SM_SUCCESS`, else the `SmStatus`; `1007` for a Rust `InvalidInput`) and
//! `P Du Dv` as the entry point leaves them: the sentinels where it does not write them
//! (upstream evaluates into locals and writes the outputs only on success).

use serde::{Deserialize, Serialize};
use serde_json::Value;

use solidcore::math::{Point2, Vec3};
use solidcore::query::brep as qb;
use solidcore::surface::surface::Surface;
use solidcore_ffi::units::query_brep::{self as ffi, EvalOp};

use super::surface_surface_cache_p1::harness;
use super::surface_surface_cache_p3::{gen_surf, push_surf, surf_domain, surface_in, SurfIn};
use super::surface_surface_p7::{st, Fail};
use crate::case::CaseContext;
use crate::runner::default_cases;
use crate::{Case, HarnessError, Output, Rng, UnitDiff};

/// Generated cases per run, whatever `--seeds` says (it can only raise this).
pub const MIN_GENERATED: usize = 600;

/// `query.brep` runner.
pub struct Runner;

// ---- inputs -----------------------------------------------------------------------

/// The evaluation entry point of an `evaluate` case.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
enum What {
    Point,
    Normal,
    Derivatives,
}

impl What {
    fn ffi(self) -> EvalOp {
        match self {
            Self::Point => EvalOp::Point,
            Self::Normal => EvalOp::Normal,
            Self::Derivatives => EvalOp::Derivatives,
        }
    }
}

/// An `evaluate` case.
#[derive(Debug, Clone, Serialize, Deserialize)]
struct EvalIn {
    what: What,
    surf: Option<SurfIn>,
    #[serde(with = "non_finite2")]
    uv: [f64; 2],
    /// `P Du Dv` before the call.
    sentinel: [f64; 9],
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "op", rename_all = "snake_case")]
enum Inputs {
    Evaluate(EvalIn),
}

/// A `[f64; 2]` whose entries may be `"NaN"`, `"inf"` or `"-inf"` (JSON has no
/// non-finite numbers).
mod non_finite2 {
    use serde::de::Error;
    use serde::{Deserialize, Deserializer, Serialize, Serializer};
    use serde_json::Value;

    use crate::json::num;

    pub fn serialize<S: Serializer>(v: &[f64; 2], s: S) -> Result<S::Ok, S::Error> {
        v.map(num).serialize(s)
    }

    pub fn deserialize<'de, D: Deserializer<'de>>(d: D) -> Result<[f64; 2], D::Error> {
        let float = |v: &Value| match v {
            Value::Number(n) => n.as_f64().ok_or_else(|| D::Error::custom(format!("bad number {n}"))),
            Value::String(s) if s == "NaN" => Ok(f64::NAN),
            Value::String(s) if s == "inf" => Ok(f64::INFINITY),
            Value::String(s) if s == "-inf" => Ok(f64::NEG_INFINITY),
            other => Err(D::Error::custom(format!("bad float {other}"))),
        };
        let [a, b] = <[Value; 2]>::deserialize(d)?;
        Ok([float(&a)?, float(&b)?])
    }
}

// ---- C++ side -----------------------------------------------------------------------

fn eval_input(i: &EvalIn) -> Result<Vec<f64>, HarnessError> {
    let mut o = Vec::new();
    match &i.surf {
        Some(s) => {
            o.push(1.0);
            push_surf(&mut o, s)?;
        }
        None => o.push(0.0),
    }
    o.extend(i.uv);
    o.extend(i.sentinel);
    Ok(o)
}

fn eval_output(flat: &[f64]) -> Output {
    let (status, rest) = flat
        .split_first()
        .map_or((f64::NAN, &[][..]), |(s, r)| (*s, r));
    Output::new()
        .values("status", vec![status])
        .values("p_du_dv", rest.to_vec())
}

fn cpp_evaluate(i: &EvalIn) -> Result<Output, HarnessError> {
    let flat = ffi::evaluate(i.what.ffi(), &eval_input(i)?).map_err(|e| harness(e.to_string()))?;
    Ok(eval_output(&flat))
}

// ---- Rust side ----------------------------------------------------------------------

fn rust_evaluate(i: &EvalIn) -> Result<Output, Fail> {
    let surface: Option<Surface> = i.surf.as_ref().map(surface_in).transpose()?;
    let uv = Point2::new(i.uv[0], i.uv[1]);
    let mut out = i.sentinel.to_vec();
    let mut put = |at: usize, v: Vec3| {
        if let Some(slot) = out.get_mut(at..at + 3) {
            slot.copy_from_slice(&[v.x, v.y, v.z]);
        }
    };
    let code = match i.what {
        What::Point => st(qb::evaluate_surface_point(surface.as_ref(), uv))?.map(|p| put(0, p)),
        What::Normal => st(qb::evaluate_surface_normal(surface.as_ref(), uv))?.map(|n| put(0, n)),
        What::Derivatives => st(qb::evaluate_surface_derivatives(surface.as_ref(), uv))?.map(
            |(p, du, dv)| {
                put(0, p);
                put(3, du);
                put(6, dv);
            },
        ),
    }
    .map_or_else(|code| code, |()| 0.0);
    let mut flat = vec![code];
    flat.extend(out);
    Ok(eval_output(&flat))
}

// ---- generation ---------------------------------------------------------------------

fn gen_what(rng: &mut Rng) -> What {
    match rng.below(3) {
        0 => What::Point,
        1 => What::Normal,
        _ => What::Derivatives,
    }
}

/// A parameter for a surface over `d`: mostly inside, else on a corner or side, just
/// outside, or non-finite.
fn gen_uv(rng: &mut Rng, d: [f64; 4]) -> [f64; 2] {
    let [u0, v0, u1, v1] = d;
    let inside = |rng: &mut Rng| [rng.uniform(u0, u1), rng.uniform(v0, v1)];
    match rng.below(12) {
        0 => [
            if rng.chance(0.5) { u0 } else { u1 },
            if rng.chance(0.5) { v0 } else { v1 },
        ],
        1 => {
            let [u, v] = inside(rng);
            if rng.chance(0.5) {
                [if rng.chance(0.5) { u0 } else { u1 }, v]
            } else {
                [u, if rng.chance(0.5) { v0 } else { v1 }]
            }
        }
        2 => {
            // just outside one bound
            let [u, v] = inside(rng);
            let du = (u1 - u0).abs().max(1.0) * 1e-9;
            let dv = (v1 - v0).abs().max(1.0) * 1e-9;
            match rng.below(4) {
                0 => [u0 - du, v],
                1 => [u1 + du, v],
                2 => [u, v0 - dv],
                _ => [u, v1 + dv],
            }
        }
        3 => {
            let bad = match rng.below(3) {
                0 => f64::NAN,
                1 => f64::INFINITY,
                _ => f64::NEG_INFINITY,
            };
            let [u, v] = inside(rng);
            if rng.chance(0.5) {
                [bad, v]
            } else {
                [u, bad]
            }
        }
        _ => inside(rng),
    }
}

fn gen_eval(rng: &mut Rng) -> Result<EvalIn, HarnessError> {
    let what = gen_what(rng);
    let mut sentinel = [0.0; 9];
    for s in &mut sentinel {
        *s = rng.uniform(-1e3, 1e3);
    }
    if rng.chance(0.03) {
        return Ok(EvalIn {
            what,
            surf: None,
            uv: [rng.unit(), rng.unit()],
            sentinel,
        });
    }
    let surf = gen_surf(rng);
    let uv = gen_uv(rng, surf_domain(&surf)?);
    Ok(EvalIn {
        what,
        surf: Some(surf),
        uv,
        sentinel,
    })
}

// ---- runner -------------------------------------------------------------------------

impl UnitDiff for Runner {
    fn id(&self) -> &'static str {
        "query.brep"
    }

    fn generate(&self, rng: &mut Rng) -> Option<Value> {
        for _ in 0..64 {
            let Ok(i) = gen_eval(rng) else {
                continue;
            };
            let Ok(input) = eval_input(&i) else {
                continue;
            };
            if ffi::validate_evaluate(&input).is_ok() {
                return serde_json::to_value(Inputs::Evaluate(i)).ok();
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

    fn run_cpp(&self, case: &Case) -> Result<Output, HarnessError> {
        match case.parse::<Inputs>()? {
            Inputs::Evaluate(i) => cpp_evaluate(&i),
        }
    }

    fn run_rust(&self, case: &Case) -> Result<Output, HarnessError> {
        let r = match case.parse::<Inputs>()? {
            Inputs::Evaluate(i) => rust_evaluate(&i),
        };
        match r {
            Ok(o) => Ok(o),
            Err(Fail::Deferred(o)) => Ok(*o),
            Err(Fail::Harness(e)) => Err(e),
        }
    }
}

