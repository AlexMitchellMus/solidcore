//! The full record of a UV trim curve, as the shims' `sc_trim_curve::put_trim_curve`
//! (`harness/cpp/shim/sc_trim_curve_record.h`, which documents the layout): the curve's
//! type, its NURB (degree, homogeneous control points, knots), its natural interval and
//! three samples whose evaluation failures are reported as their status.

use solidcore::curve::curve::p2::{Curve, CurveKind};
use solidcore::kernel::Kernel;
use solidcore::SolidError;

/// Appends the record of `c` (`-1` for none). A failed evaluation that is not an
/// upstream status (an unported callee) is returned.
pub(crate) fn put_trim_curve(
    k: &Kernel,
    o: &mut Vec<f64>,
    c: Option<&Curve>,
) -> Result<(), SolidError> {
    let Some(c) = c else {
        o.push(-1.0);
        return Ok(());
    };
    o.extend([1.0, c.type_code() as f64]);
    let base = match c.kind() {
        CurveKind::BSpline(x) => Some(x),
        CurveKind::Ellipse(x) => Some(x.base()),
        CurveKind::Hyperbola(x) => Some(x.base()),
        CurveKind::Parabola(x) => Some(x.base()),
        CurveKind::Line(x) => Some(x.base()),
        _ => None,
    };
    match base.and_then(|x| x.gw_nurb()) {
        Some(cur) => {
            o.extend([cur.n() as f64, cur.degree() as f64, cur.m() as f64]);
            for p in cur.cpts() {
                o.extend([p.x, p.y, p.z, p.w]);
            }
            o.extend_from_slice(cur.knots());
        }
        None => o.extend([-1.0; 3]),
    }
    let ivl = c.natural_interval();
    o.extend([ivl.min(), ivl.max()]);
    for i in 0..=2u8 {
        match c.evaluate_point(k, ivl.evaluate(f64::from(i) / 2.0)) {
            Ok(p) => o.extend([0.0, p.x, p.y, p.z]),
            Err(SolidError::Upstream { code, .. }) => o.push(f64::from(code)),
            Err(e) => return Err(e),
        }
    }
    Ok(())
}

/// The length of the record [`put_trim_curve`] wrote at the start of `r`, or `None` when
/// `r` is too short to hold one.
pub(crate) fn trim_curve_len(r: &[f64]) -> Option<usize> {
    if *r.first()? == -1.0 {
        return Some(1);
    }
    // 1 type n p m
    let n = *r.get(2)?;
    let m = *r.get(4)?;
    let mut at = 5;
    if n >= 0.0 {
        at += 4 * (n as usize + 1) + (m as usize + 1);
    }
    at += 2; // t0 t1
    for _ in 0..3 {
        // `0 x y z`, or a failed evaluation's status
        at += if *r.get(at)? == 0.0 { 4 } else { 1 };
    }
    (at <= r.len()).then_some(at)
}
