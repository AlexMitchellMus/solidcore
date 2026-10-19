//! Fillet laws: the radius (and its derivatives) along a filleted edge, constant,
//! linear or given by a B-spline (part 1 of `SmFilletExecutive`: the `SmFilletLaw`
//! family).
//!
//! Upstream: `source/SMLib/src/SmFilletExecutive.cpp` (lines 5126-6132),
//! `source/SMLib/inc/SmFilletExecutive.h` (lines 660-899).
//!
//! The law type itself ([`FilletLaw`], [`FilletLawKind`]) was CO-PORTED by unit
//! `blend.fillet_standard_solver_p1` into
//! [`crate::blend::fillet_standard_solver::p1`]; it is extended there (edge map, arc
//! length, cached extension) and re-exported here. This module holds its behaviour.
//!
//! | Rust                                   | Upstream                                              |
//! |----------------------------------------|-------------------------------------------------------|
//! | [`FilletLaw::new`] (co-port)           | `SmFilletLaw::SmFilletLaw`, `~SmFilletLaw`, `~SmBSplineFilletLaw` (`Drop`) |
//! | [`FilletLaw::evaluate`]                | `SmFilletLaw::Evaluate`, `SmLinearFilletLaw::Evaluate`, `SmBSplineFilletLaw::Evaluate` |
//! | [`FilletLaw::is_linear_law`] (co-port) | `SmFilletLaw::IsLinearLaw`, `SmLinearFilletLaw::IsLinearLaw`, `SmBSplineFilletLaw::IsLinearLaw` |
//! | [`FilletLaw::law_curve`]               | `SmBSplineFilletLaw::GetLawCurve`                     |
//! | [`FilletLaw::extended_curve`]          | `SmBSplineFilletLaw::GetExtendedCurve`                |
//! | field `FilletLawKind::BSpline::law_curve` / `extended_law` | `SetLawCurve`, `SetExtendedCurve`     |
//! | field `FilletLawKind::Linear::start_radius` / `end_radius` | `Get/SetStartRadius`, `Get/SetEndRadius` |
//! | field `FilletLaw::edge_arc_length`     | `SmFilletLaw::SetEdgeArcLength`                       |
//! | [`FilletLaw::extended_edge_curve`], [`FilletLaw::set_extended_edge_curve`] (co-port) | `GetExtendedEdgeCurve`, `SetExtendedEdgeCurve` |
//! | [`FilletLaw::type_code`], [`FilletLaw::type_string`] | `GetType`, `GetTypeString` of the three classes |
//! | [`FilletLaw::is_kind_of`]              | `IsKindOf` of the three classes                       |
//! | [`FILLET_LAW_TYPE`], [`LINEAR_FILLET_LAW_TYPE`], [`BSPLINE_FILLET_LAW_TYPE`] | `SmFilletLaw_TYPE`, `SmLinearFilletLaw_TYPE`, `SmBSplineFilletLaw_TYPE` |
//!
//! The `SmLinearFilletLaw` and `SmBSplineFilletLaw` constructors belong to unit
//! `blend.fillet_executive_p2`: a law is built with [`FilletLaw::new`] and its public
//! fields until they land.
//!
//! Deviations:
//! - The class hierarchy is the enum [`FilletLawKind`]; the trivial getters and setters
//!   are public fields.
//! - `Evaluate` is `const` upstream and caches the B-spline extension through a cast
//!   away from `const`; here [`FilletLaw::evaluate`] takes `&mut self` and stores it in
//!   `FilletLawKind::BSpline::extended_law`. The context passed to
//!   `CreateExtendedCurve` only allocates, so no `Kernel` is needed.
//! - The "unable to extend non-B-spline laws" branch is unreachable (both the law curve
//!   and its extension are B-splines) and not ported.
//! - A G1-closed law curve is not extended: `CreateExtendedCurve` hands back an
//!   unextended copy (NLib evaluates such a curve periodically) and the evaluation
//!   clamps, as upstream. `CreateExtendedCurve` returning no curve at all is impossible
//!   for a parameter outside the domain; that path (a `NULL` dereference upstream) is
//!   `SM_ERR_NULL_POINTER` ([`SolidError::null_pointer`]).
//! - The out-parameter `dValues[3]` is the returned `[f64; 3]`.
//! - No recursion.

pub use crate::blend::fillet_standard_solver::p1::{FilletLaw, FilletLawKind};
use crate::error::{Result, SolidError};
use crate::geom::extent1d::Extent1d;
use crate::geom::geom_utility::p2::ContinuityType;
use crate::geom::object;
use crate::math::EFF_ZERO;
use crate::nurbs::b_spline_curve::BSplineCurve;

const OP_EVAL: &str = "blend::fillet_executive::p1::FilletLaw::evaluate";

/// `TOPO_BASE_TYPE`.
const TOPO_BASE_TYPE: i64 = 16000;
/// `SmFilletLaw_TYPE` (`TOPO_BASE_TYPE + 39`).
pub const FILLET_LAW_TYPE: i64 = TOPO_BASE_TYPE + 39;
/// `SmLinearFilletLaw_TYPE` (`TOPO_BASE_TYPE + 40`).
pub const LINEAR_FILLET_LAW_TYPE: i64 = TOPO_BASE_TYPE + 40;
/// `SmBSplineFilletLaw_TYPE` (`TOPO_BASE_TYPE + 41`).
pub const BSPLINE_FILLET_LAW_TYPE: i64 = TOPO_BASE_TYPE + 41;

impl FilletLaw {
    /// `GetLawCurve`: the B-spline radius function (`None` for the other kinds).
    #[must_use]
    pub const fn law_curve(&self) -> Option<&BSplineCurve> {
        match &self.kind {
            FilletLawKind::BSpline { law_curve, .. } => Some(law_curve),
            _ => None,
        }
    }

    /// `GetExtendedCurve`: the extension of the law curve the last out-of-domain
    /// evaluation created (`None` before that, and for the other kinds).
    #[must_use]
    pub const fn extended_curve(&self) -> Option<&BSplineCurve> {
        match &self.kind {
            FilletLawKind::BSpline { extended_law, .. } => extended_law.as_ref(),
            _ => None,
        }
    }

    /// `GetType`: the run-time type code of the law's class.
    #[must_use]
    pub const fn type_code(&self) -> i64 {
        match &self.kind {
            FilletLawKind::Base => FILLET_LAW_TYPE,
            FilletLawKind::Linear { .. } => LINEAR_FILLET_LAW_TYPE,
            FilletLawKind::BSpline { .. } => BSPLINE_FILLET_LAW_TYPE,
        }
    }

    /// `GetTypeString`.
    #[must_use]
    pub const fn type_string(&self) -> &'static str {
        match &self.kind {
            FilletLawKind::Base => "SmFilletLaw_TYPE",
            FilletLawKind::Linear { .. } => "SmLinearFilletLaw_TYPE",
            FilletLawKind::BSpline { .. } => "SmBSplineFilletLaw_TYPE",
        }
    }

    /// `IsKindOf`: the law's own type, `SmFilletLaw_TYPE` or `SmObject_TYPE`.
    #[must_use]
    pub const fn is_kind_of(&self, t: i64) -> bool {
        t == self.type_code() || t == FILLET_LAW_TYPE || object::is_kind_of(t)
    }

    /// `Evaluate(dParameter, crCurveInterval, bReverseOrientation, dValues)`: the
    /// fillet radius at `parameter` of the edge interval `interval` and its first two
    /// derivatives with respect to that parameter (`reverse`: the law runs from the
    /// interval's end to its start).
    ///
    /// A linear law evaluates anywhere. A B-spline law maps the parameter through the
    /// edge map onto the law curve; a parameter within `SM_EFF_ZERO * (1 + |min| +
    /// |max|)` outside the curve's domain snaps to the end, a parameter further out
    /// extends the curve (`C-infinity`, to twice the overshoot and at least a tenth of
    /// the domain) and keeps the extension for later evaluations.
    ///
    /// # Errors
    ///
    /// `SM_ERR` for the base law and for a zero-length `interval` of a B-spline law;
    /// the errors of [`BSplineCurve::create_extended_curve`] and
    /// [`BSplineCurve::evaluate`].
    pub fn evaluate(
        &mut self,
        parameter: f64,
        interval: &Extent1d,
        reverse: bool,
    ) -> Result<[f64; 3]> {
        let edge_map = self.edge_map;
        match &mut self.kind {
            FilletLawKind::Base => Err(SolidError::sm_err(OP_EVAL)),
            FilletLawKind::Linear {
                start_radius,
                end_radius,
            } => Ok(linear_evaluate(
                *start_radius,
                *end_radius,
                &edge_map,
                parameter,
                interval,
                reverse,
            )),
            FilletLawKind::BSpline {
                law_curve,
                extended_law,
            } => bspline_evaluate(
                law_curve,
                extended_law,
                &edge_map,
                parameter,
                interval,
                reverse,
            ),
        }
    }
}

/// `SmLinearFilletLaw::Evaluate`.
fn linear_evaluate(
    start_value: f64,
    end_value: f64,
    edge_map: &Extent1d,
    parameter: f64,
    interval: &Extent1d,
    reverse: bool,
) -> [f64; 3] {
    // Map dParameter into domain of Law which is contained in the edge map.
    let d_min = interval.min();
    let d_max = interval.max();
    let t_norm = if reverse {
        (d_max - parameter) / (d_max - d_min)
    } else {
        (parameter - d_min) / (d_max - d_min)
    };
    let t = edge_map.min() + t_norm * (edge_map.max() - edge_map.min());

    // evaluate function, 1st, and 2nd derivatives at param = dParameter
    let value = start_value + t * (end_value - start_value);
    let end_map_value = start_value + edge_map.max() * (end_value - start_value);
    let start_map_value = start_value + edge_map.min() * (end_value - start_value);
    let mut deriv1 = (end_map_value - start_map_value) / interval.length();
    let deriv2 = 0.0;

    // negate 1st derivatives for reversed orientations
    if reverse {
        deriv1 = -deriv1;
    }

    [value, deriv1, deriv2]
}

/// `SmBSplineFilletLaw::Evaluate`.
fn bspline_evaluate(
    law_curve: &BSplineCurve,
    extended_law: &mut Option<BSplineCurve>,
    edge_map: &Extent1d,
    parameter: f64,
    interval: &Extent1d,
    reverse: bool,
) -> Result<[f64; 3]> {
    if interval.length() == 0.0 {
        return Err(SolidError::sm_err(OP_EVAL));
    }

    // select function to evaluate
    let law = extended_law.as_ref().unwrap_or(law_curve);

    // Compute the scale difference between the parameters of the map and the curve.
    let ivl = law.natural_interval();
    let crv_ivl_size = interval.length();
    let map_ivl_size = edge_map.length();
    let param_scale = crv_ivl_size / map_ivl_size;

    // Map parameter to unitized range accounting for bReverseOrientation
    let d_min = interval.min();
    let d_max = interval.max();

    let t = if reverse {
        (d_max - parameter) / (d_max - d_min)
    } else {
        (parameter - d_min) / (d_max - d_min)
    };

    // Map parameter back to interval of the curve.
    let mut law_t = edge_map.evaluate(t);

    // when parameter is beyond current range
    let mut new_extension = None;
    if !ivl.contains_value(law_t, 0.0) {
        // select param space tolerance
        let scaled_tol = EFF_ZERO * (1.0 + ivl.min().abs() + ivl.max().abs());

        // Try to snap to begin point.
        if law_t < ivl.min() && law_t + scaled_tol > ivl.min() {
            law_t = ivl.min();
        }

        // Try to snap to end point.
        if law_t > ivl.max() && law_t - scaled_tol < ivl.max() {
            law_t = ivl.max();
        }

        // when parameter is still beyond current range
        if !ivl.contains_value(law_t, 0.0) {
            // create an extension of the curve to the requested param point
            // But extend farther than requested, so that we won't have to do this too often.
            // This doubles the extension distance:
            let mut ext_param = if law_t < ivl.min() {
                law_t - (ivl.min() - law_t)
            } else {
                law_t + (law_t - ivl.max())
            };
            // But never do a really tiny extension: it might get snapped to the end,
            // and therefor not extended at all.  [Fillet 308]
            let len = ivl.length();
            if ext_param < ivl.min() {
                if ext_param > ivl.min() - len / 10.0 {
                    ext_param = ivl.min() - len / 10.0;
                }
            } else if ext_param < ivl.max() + len / 10.0 {
                ext_param = ivl.max() + len / 10.0;
            }

            // precise extension (otherwise it adds another length of sIvl);
            // C-infinity works better than G1_G2.  [B551]
            new_extension =
                Some(law.create_extended_curve(ext_param, ContinuityType::CInfinity, true)?);
        } // end need to extend fillet-radius function check after snapping to endPoints
    } // end need to extend fillet-radius function check

    if let Some(extension) = new_extension {
        // save and use the extended curve (the previous one is dropped)
        *extended_law = extension;
        if extended_law.is_none() {
            // unreachable for an out-of-domain parameter (upstream would dereference NULL)
            return Err(SolidError::null_pointer(OP_EVAL));
        }
    }
    let law = extended_law.as_ref().unwrap_or(law_curve);

    let pv = law.evaluate(law_t, 2, true, true)?;
    let x = |i: usize| {
        pv.get(i)
            .map(|v| v.x)
            .ok_or_else(|| SolidError::sm_err(OP_EVAL))
    };

    let value = x(0)?;
    let mut deriv1 = x(1)? / param_scale;
    if reverse {
        deriv1 = -deriv1;
    }
    let deriv2 = x(2)? / (param_scale * param_scale);

    Ok([value, deriv1, deriv2])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::nurbs::crv_basic::Curve as NurbCurve;
    use crate::nurbs::math::{CPoint, NOW};

    /// A clamped degree-`p` curve with the radii `xs` as x coordinates over `[a, b]`.
    fn radius_curve(xs: &[f64], p: usize, a: f64, b: f64) -> BSplineCurve {
        let pw = xs.iter().map(|&x| CPoint::new(x, 0.0, 0.0, NOW)).collect();
        let spans = xs.len() - p;
        let mut u = vec![a; p + 1];
        for i in 1..spans {
            u.push(a + (b - a) * i as f64 / spans as f64);
        }
        u.extend(std::iter::repeat_n(b, p + 1));
        BSplineCurve::from_nurb(NurbCurve::new(pw, p, u).expect("curve"), 3, false)
    }

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= 1e-12 * (1.0 + a.abs().max(b.abs()))
    }

    fn linear(r0: f64, r1: f64) -> FilletLaw {
        FilletLaw::new(FilletLawKind::Linear {
            start_radius: r0,
            end_radius: r1,
        })
    }

    fn bspline(law_curve: BSplineCurve, map: Extent1d) -> FilletLaw {
        let mut law = FilletLaw::new(FilletLawKind::BSpline {
            law_curve,
            extended_law: None,
        });
        law.edge_map = map;
        law
    }

    #[test]
    fn linear_law_values_and_reverse() {
        let mut law = linear(2.0, 0.0);
        assert_eq!(law.edge_map.min(), 0.0);
        assert_eq!(law.edge_map.max(), 1.0);
        assert_eq!(law.edge_arc_length, 0.0);
        let ivl = Extent1d::from_min_max(0.0, 10.0);
        assert_eq!(
            law.evaluate(2.5, &ivl, false).expect("eval"),
            [1.5, -0.2, 0.0]
        );
        assert_eq!(
            law.evaluate(2.5, &ivl, true).expect("eval"),
            [0.5, 0.2, 0.0]
        );
        // a sub-map: only [0.5, 1] of the function maps to the edge
        let mut law = linear(2.0, 0.0);
        law.edge_map = Extent1d::from_min_max(0.5, 1.0);
        assert_eq!(
            law.evaluate(0.0, &ivl, false).expect("eval"),
            [1.0, -0.1, 0.0]
        );
        // a zero-length map means none
        assert!(law.is_linear_law().expect("linear"));
        assert_eq!(law.type_code(), LINEAR_FILLET_LAW_TYPE);
        assert!(law.is_kind_of(FILLET_LAW_TYPE));
        assert!(law.is_kind_of(11001));
        assert!(!law.is_kind_of(BSPLINE_FILLET_LAW_TYPE));
    }

    #[test]
    fn base_law_errors() {
        let mut law = FilletLaw::new(FilletLawKind::Base);
        assert!(law
            .evaluate(0.5, &Extent1d::from_min_max(0.0, 1.0), false)
            .is_err());
        assert_eq!(law.type_string(), "SmFilletLaw_TYPE");
        assert!(law.is_kind_of(FILLET_LAW_TYPE));
        assert!(!law.is_linear_law().expect("linear"));
    }

    #[test]
    fn bspline_law_inside_snap_and_extension() {
        // a degree-1 law from 1 to 3 over [0, 2], mapped onto the edge interval [10, 14]
        let mut law = bspline(
            radius_curve(&[1.0, 3.0], 1, 0.0, 2.0),
            Extent1d::from_min_max(0.0, 2.0),
        );
        assert!(law.is_linear_law().expect("linear"));
        assert_eq!(law.type_code(), BSPLINE_FILLET_LAW_TYPE);
        let ivl = Extent1d::from_min_max(10.0, 14.0);
        let v = law.evaluate(11.0, &ivl, false).expect("eval");
        assert!(
            close(v[0], 1.5) && close(v[1], 0.5) && close(v[2], 0.0),
            "{v:?}"
        );
        let v = law.evaluate(11.0, &ivl, true).expect("eval");
        assert!(close(v[0], 2.5) && close(v[1], -0.5), "{v:?}");
        assert!(law.extended_curve().is_none());
        // just outside: snaps to the end
        let v = law.evaluate(14.0 + 1e-13, &ivl, false).expect("eval");
        assert!(close(v[0], 3.0), "{v:?}");
        assert!(law.extended_curve().is_none());
        // further outside: extends the line (to at least a tenth of the domain)
        let v = law.evaluate(15.0, &ivl, false).expect("eval");
        assert!(close(v[0], 3.5) && close(v[1], 0.5), "{v:?}");
        let ext = law.extended_curve().expect("extended").natural_interval();
        assert!(ext.max() >= 3.0 && ext.min() <= 0.0, "{ext:?}");
        // zero-length interval
        assert!(law
            .evaluate(1.0, &Extent1d::from_min_max(3.0, 3.0), false)
            .is_err());
    }

    #[test]
    fn closed_law_curve_extension_is_a_copy() {
        // smFilletTest.cpp:307: the cubic r, 1.2r, 0.8r, r is G1-closed, so NLib hands
        // back an unextended copy and the evaluation clamps to the end
        let mut law = bspline(
            radius_curve(&[1.0, 1.2, 0.8, 1.0], 3, 0.0, 10.0),
            Extent1d::from_min_max(0.0, 10.0),
        );
        let ivl = Extent1d::from_min_max(0.0, 10.0);
        let v = law.evaluate(5.0, &ivl, false).expect("eval");
        assert!(close(v[0], 1.0), "{v:?}");
        let v = law.evaluate(12.0, &ivl, false).expect("eval");
        assert!(close(v[0], 1.0), "{v:?}");
        let ext = law.extended_curve().expect("copy").natural_interval();
        assert_eq!((ext.min(), ext.max()), (0.0, 10.0));
    }
}
