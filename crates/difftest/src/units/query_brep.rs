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
//!
//! An `attribute` case is a decimation attribute class (`decimate_vertex`,
//! `quadric_vertex` or `quadric_edge`), its id, an optional behaviour (`null` for the
//! constructor's default) and optional data written before the copy, laid out as the
//! ffi module documents. Both sides produce `id behavior type data` of the constructed
//! attribute, then of its `MakeCopy` (the quadric edge's sort index NaN while unset).

use serde::{Deserialize, Serialize};
use serde_json::Value;

use solidcore::geom::attribute::{Attribute, AttributeValue, Behavior};
use solidcore::math::{Point2, Vec3};
use solidcore::query::brep::{self as qb, PolyVertexClass};
use solidcore::surface::surface::Surface;
use solidcore_ffi::units::query_brep::{self as ffi, AttrKind, EvalOp};

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

/// The class of an `attribute` case.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
enum Kind {
    DecimateVertex,
    QuadricVertex,
    QuadricEdge,
}

impl Kind {
    fn ffi(self) -> AttrKind {
        match self {
            Self::DecimateVertex => AttrKind::DecimateVertex,
            Self::QuadricVertex => AttrKind::QuadricVertex,
            Self::QuadricEdge => AttrKind::QuadricEdge,
        }
    }
}

/// An `attribute` case.
#[derive(Debug, Clone, Serialize, Deserialize)]
struct AttrIn {
    kind: Kind,
    id: u64,
    behavior: Option<u8>,
    /// The data record written before the copy (see the ffi module).
    data: Option<Vec<f64>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "op", rename_all = "snake_case")]
enum Inputs {
    Evaluate(EvalIn),
    Attribute(AttrIn),
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

fn attr_input(i: &AttrIn) -> Vec<f64> {
    let mut o = vec![
        i.id as f64,
        f64::from(u8::from(i.behavior.is_some())),
        f64::from(i.behavior.unwrap_or(0)),
        f64::from(u8::from(i.data.is_some())),
    ];
    o.extend(i.data.iter().flatten());
    o
}

fn cpp_attribute(i: &AttrIn) -> Result<Output, HarnessError> {
    let flat = ffi::attribute(i.kind.ffi(), &attr_input(i)).map_err(|e| harness(e.to_string()))?;
    Ok(Output::new().values("attributes", flat))
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

/// Reads the `at`th entry of a data record.
fn datum(d: &[f64], at: usize) -> Result<f64, HarnessError> {
    d.get(at)
        .copied()
        .ok_or_else(|| harness("attribute data too short"))
}

fn vec_at(d: &[f64], at: usize) -> Result<Vec3, HarnessError> {
    Ok(Vec3::new(datum(d, at)?, datum(d, at + 1)?, datum(d, at + 2)?))
}

/// Writes a data record into the attribute (as the shim's `set_data`).
fn set_data(a: &mut Attribute, d: &[f64]) -> Result<(), HarnessError> {
    match a.value_mut() {
        AttributeValue::DecimateVertex(v) => {
            v.sort_index = datum(d, 0)? as u64;
            let class = datum(d, 1)?;
            v.vertex_class = PolyVertexClass::ALL
                .get(class as usize)
                .copied()
                .ok_or_else(|| harness(format!("bad vertex class {class}")))?;
            v.vertex_distance = datum(d, 2)?;
            v.max_face_error = datum(d, 3)?;
            v.average_plane_normal = vec_at(d, 4)?;
        }
        AttributeValue::QuadricVertex(v) => {
            v.q = [vec_at(d, 0)?, vec_at(d, 3)?, vec_at(d, 6)?];
            v.b = vec_at(d, 9)?;
            v.c = datum(d, 12)?;
        }
        AttributeValue::QuadricEdge(v) => {
            v.new_point = vec_at(d, 0)?;
            v.error = datum(d, 3)?;
            if datum(d, 4)? != 0.0 {
                v.sort_index = Some(datum(d, 5)? as u64);
            }
        }
        _ => return Err(harness("not a decimation attribute")),
    }
    Ok(())
}

fn push_attribute(o: &mut Vec<f64>, a: &Attribute) -> Result<(), HarnessError> {
    let pv = |o: &mut Vec<f64>, v: Vec3| o.extend([v.x, v.y, v.z]);
    o.extend([
        a.attribute_id() as f64,
        f64::from(a.behavior().code()),
        a.type_code() as f64,
    ]);
    match a.value() {
        AttributeValue::DecimateVertex(v) => {
            o.extend([
                v.sort_index as f64,
                f64::from(v.vertex_class as u8),
                v.vertex_distance,
                v.max_face_error,
            ]);
            pv(o, v.average_plane_normal);
        }
        AttributeValue::QuadricVertex(v) => {
            for r in v.q {
                pv(o, r);
            }
            pv(o, v.b);
            o.push(v.c);
        }
        AttributeValue::QuadricEdge(v) => {
            pv(o, v.new_point);
            o.push(v.error);
            o.push(v.sort_index.map_or(f64::NAN, |s| s as f64));
        }
        _ => return Err(harness("not a decimation attribute")),
    }
    Ok(())
}

fn rust_attribute(i: &AttrIn) -> Result<Output, HarnessError> {
    let behavior = match i.behavior {
        Some(b) => Behavior::from_code(i32::from(b)).map_err(|e| harness(e.to_string()))?,
        // the constructors' default argument
        None => Behavior::Copy,
    };
    let mut a = match i.kind {
        Kind::DecimateVertex => qb::decimate_vertex_attribute(i.id, behavior),
        Kind::QuadricVertex => qb::quadric_vertex_attribute(i.id, behavior),
        Kind::QuadricEdge => qb::quadric_edge_attribute(i.id, behavior),
    };
    if let Some(d) = &i.data {
        set_data(&mut a, d)?;
    }
    let c = a.make_copy();
    let mut o = Vec::new();
    push_attribute(&mut o, &a)?;
    push_attribute(&mut o, &c)?;
    Ok(Output::new().values("attributes", o))
}

// ---- generation ---------------------------------------------------------------------

/// A value for a data field: mostly moderate, sometimes zero, huge, tiny or non-finite
/// (the attributes copy whatever they hold).
fn gen_value(rng: &mut Rng) -> f64 {
    match rng.below(16) {
        0 => 0.0,
        1 => -1e300,
        2 => 1e-300,
        3 => f64::MIN_POSITIVE / 4.0,
        _ => rng.uniform(-1e3, 1e3),
    }
}

/// An exact `ULONG` the ffi accepts.
fn gen_ulong(rng: &mut Rng) -> u64 {
    match rng.below(6) {
        0 => 0,
        1 => u64::from(u32::MAX),
        2 => 1 << 53,
        _ => rng.below(1 << 20),
    }
}

fn gen_attr(rng: &mut Rng) -> AttrIn {
    let kind = match rng.below(3) {
        0 => Kind::DecimateVertex,
        1 => Kind::QuadricVertex,
        _ => Kind::QuadricEdge,
    };
    let id = gen_ulong(rng);
    let behavior = rng.chance(0.7).then(|| rng.below(5) as u8);
    let data = rng.chance(0.75).then(|| match kind {
        Kind::DecimateVertex => {
            let mut d = vec![gen_ulong(rng) as f64, rng.below(6) as f64];
            d.extend((0..5).map(|_| gen_value(rng)));
            d
        }
        Kind::QuadricVertex => (0..13).map(|_| gen_value(rng)).collect(),
        Kind::QuadricEdge => {
            let mut d: Vec<f64> = (0..4).map(|_| gen_value(rng)).collect();
            d.push(f64::from(u8::from(rng.chance(0.5))));
            d.push(gen_ulong(rng) as f64);
            d
        }
    });
    AttrIn {
        kind,
        id,
        behavior,
        data,
    }
}


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
        if rng.chance(0.2) {
            let i = gen_attr(rng);
            return serde_json::to_value(Inputs::Attribute(i)).ok();
        }
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
            Inputs::Attribute(i) => cpp_attribute(&i),
        }
    }

    fn run_rust(&self, case: &Case) -> Result<Output, HarnessError> {
        let r = match case.parse::<Inputs>()? {
            Inputs::Evaluate(i) => rust_evaluate(&i),
            Inputs::Attribute(i) => return rust_attribute(&i),
        };
        match r {
            Ok(o) => Ok(o),
            Err(Fail::Deferred(o)) => Ok(*o),
            Err(Fail::Harness(e)) => Err(e),
        }
    }
}

