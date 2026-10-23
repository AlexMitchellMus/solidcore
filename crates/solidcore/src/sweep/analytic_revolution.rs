//! The analytic surface a revolved profile curve sweeps: a cylinder, cone, sphere or
//! torus bounded to the curve and the turn (SolidCraft change).
//!
//! Upstream: none. SolidCraft fix tangency-labs/smlib 0b143f7 (branch
//! `port-revolve-analytic`, tangency-labs/solidcore#545). At `c9c5979`
//! `SmPrimitiveCreation::CreateRotationalSweep` makes analytic faces only for a
//! one-curve profile (`sm_CreateAnalyticRotationalSweep`). Each curve of a longer
//! profile gets an `SmSurfOfRevolution` (a whole turn) or the B-spline surface of
//! revolution of `SmRotationalSweepGeometry::EdgeSweepHigher` (a partial turn), whose
//! `SmBSplineSurface::CopyAndAddAnalytics` recognises no cylinder, cone, sphere or torus
//! the face's trims can stay on.
//!
//! | Rust                    | Upstream      |
//! |-------------------------|---------------|
//! | [`analytic_revolution`] | (SolidCraft)  |
//!
//! Deviations: not in upstream. The surface is laid out as `sm_CreateRotSrf`'s B-spline
//! surface of revolution: `u` the turn (`0` at the curve, rising with the rotation about
//! the axis), `v` the curve from its start to its end. It is only returned when its
//! corners, the curve's start tangent and the turn's sense agree with that layout, so a
//! caller can swap it in for the B-spline and keep the trims and loop the B-spline
//! layout gives. Any other curve, or a construction that fails, gives `None` and the
//! caller keeps the upstream surface. No recursion.

use crate::curve::curve::Curve;
use crate::error::{ignored, Result};
use crate::geom::axis2_placement::Axis2Placement;
use crate::geom::extent2d::Extent2d;
use crate::geom::geom_utility::p1::{line_closest_point, line_point_distance};
use crate::geom::vector3d::{Point3, Vec3};
use crate::kernel::Kernel;
use crate::math::{EFF_ZERO, PI};
use crate::nurbs::b_spline_surface_n_lib::p2::EvaluateFlags;
use crate::surface::cone::Cone;
use crate::surface::sphere::{Sphere, SphereOptions};
use crate::surface::surface::Surface;
use crate::surface::torus::{Torus, TorusOptions};

/// Curve samples tested for lying in the half-plane through the axis.
const PLANE_SAMPLES: u32 = 8;

/// The analytic surface swept by revolving `curve` (over its natural interval) by
/// `angle_deg` degrees about the axis through `base_pt` along `axis`: a cylinder for a
/// line parallel to the axis, a cone for another line that is not perpendicular to
/// it, a sphere for a circular arc centred on the axis and a torus for one off it, the
/// arc in a plane through the axis. Bounded to the curve and the turn, and laid out as
/// the B-spline surface of revolution (see the module header): when `reference` is
/// given (that B-spline), the result must reproduce its corners, its tangent along the
/// curve at the start and its sense of turn, else the layout the curve itself gives.
/// `tol` is the distance tolerance of those tests and of the shape recognition.
///
/// `None` for any other curve (a line perpendicular to the axis sweeps a plane), a
/// curve not in one half-plane through the axis, or a surface that does not pass.
///
/// # Errors
///
/// [`crate::SolidError::Unported`] from the curve queries or surface construction;
/// other failures give `None`.
pub fn analytic_revolution(
    k: &Kernel,
    curve: &Curve,
    base_pt: Point3,
    axis: Vec3,
    angle_deg: f64,
    tol: f64,
    reference: Option<&Surface>,
) -> Result<Option<Surface>> {
    let Some(axis) = ignored(axis.unitize())? else {
        return Ok(None);
    };
    // (written so that a NaN turn or tolerance is refused)
    let usable = angle_deg > 0.0 && angle_deg <= 360.0 + EFF_ZERO && tol > 0.0;
    if !usable {
        return Ok(None);
    }
    let Some(Some(candidate)) = ignored(candidate(k, curve, base_pt, axis, angle_deg, tol))? else {
        return Ok(None);
    };
    let expected = match reference {
        Some(surface) => ignored(Layout::of(k, surface, base_pt, axis))?,
        None => ignored(Layout::of_curve(k, curve, base_pt, axis, angle_deg))?,
    };
    let actual = ignored(Layout::of(k, &candidate, base_pt, axis))?;
    Ok(match (expected, actual) {
        (Some(expected), Some(actual)) if actual.matches(&expected, tol) => Some(candidate),
        _ => None,
    })
}

/// The unit direction from the axis to `point`, `None` within `tol` of the axis.
fn radial(base_pt: Point3, axis: Vec3, point: Point3, tol: f64) -> Option<Vec3> {
    let foot = base_pt + axis * line_closest_point(base_pt, axis, point);
    let offset = point - foot;
    let length = offset.length();
    (length > tol).then(|| offset / length)
}

/// An angle in degrees taken into `[0, 360)`.
fn turn_into_range(angle_deg: f64) -> f64 {
    let a = angle_deg % 360.0;
    if a < 0.0 {
        a + 360.0
    } else {
        a
    }
}

/// Whether the curve lies in the half-plane bounded by the axis and holding `x`
/// (`y = axis x x`), within `tol`, at its ends and `PLANE_SAMPLES - 1` inner points.
fn in_half_plane(
    k: &Kernel,
    curve: &Curve,
    base_pt: Point3,
    x: Vec3,
    y: Vec3,
    tol: f64,
) -> Result<bool> {
    let interval = curve.natural_interval();
    for i in 0..=PLANE_SAMPLES {
        let t = interval.evaluate(f64::from(i) / f64::from(PLANE_SAMPLES));
        let p = curve.evaluate_point(k, t)? - base_pt;
        if p.dot(y).abs() > tol || p.dot(x) < -tol {
            return Ok(false);
        }
    }
    Ok(true)
}

/// The analytic surface for the curve's shape, before the layout test.
fn candidate(
    k: &Kernel,
    curve: &Curve,
    base_pt: Point3,
    axis: Vec3,
    angle_deg: f64,
    tol: f64,
) -> Result<Option<Surface>> {
    let ends = curve.ends(k)?;
    let (start_dist, _) = line_point_distance(base_pt, axis, ends.start);
    let (end_dist, _) = line_point_distance(base_pt, axis, ends.end);

    if curve.is_linear(k, EFF_ZERO, false)?.linear {
        let start_h = line_closest_point(base_pt, axis, ends.start);
        let end_h = line_closest_point(base_pt, axis, ends.end);
        let height = end_h - start_h;
        // perpendicular to the axis: a plane, which the callers make
        if height.abs() <= tol {
            return Ok(None);
        }
        // from the end farther from the axis (the other may be the apex)
        let far = if start_dist >= end_dist {
            ends.start
        } else {
            ends.end
        };
        let Some(x) = radial(base_pt, axis, far, tol) else {
            return Ok(None);
        };
        let y = axis.cross(x);
        if !in_half_plane(k, curve, base_pt, x, y, tol)? {
            return Ok(None);
        }
        // the generator line rises along the axis from the origin; inside out runs
        // the NURB from the top, so that v always runs from the curve's start
        let (low_h, bot, top, inside_out) = if height > 0.0 {
            (start_h, start_dist, end_dist, false)
        } else {
            (end_h, end_dist, start_dist, true)
        };
        let origin = base_pt + axis * low_h;
        let mut cone = if (start_dist - end_dist).abs() <= tol {
            Cone::new_cylinder_bounded(
                origin,
                x,
                y,
                start_dist,
                0.0,
                angle_deg,
                height.abs(),
                false,
                inside_out,
                false,
            )
        } else {
            Cone::new_bounded(
                origin,
                x,
                y,
                bot,
                top,
                0.0,
                angle_deg,
                height.abs(),
                false,
                inside_out,
                false,
            )
        };
        cone.make_nurb(k)?;
        return Ok(Some(Surface::from(cone)));
    }

    let Some((frame, radius, arc_start, arc_end)) = curve.is_arc(k, 8, tol)? else {
        return Ok(None);
    };
    let sweep = arc_end - arc_start;
    let is_arc = sweep > 0.0 && sweep <= 360.0 * (1.0 + EFF_ZERO);
    if !is_arc {
        return Ok(None);
    }
    let centre = frame.origin();
    let (centre_dist, _) = line_point_distance(base_pt, axis, centre);
    let interval = curve.natural_interval();
    let mid = curve.evaluate_point(k, interval.evaluate(0.5))?;
    let to_deg = 180.0 / PI;

    if centre_dist <= tol {
        // a sphere: its meridian is the arc, in the plane of x and the axis
        let Some(x) = radial(base_pt, axis, mid, tol) else {
            return Ok(None);
        };
        let y = axis.cross(x);
        if !in_half_plane(k, curve, base_pt, x, y, tol)? {
            return Ok(None);
        }
        let latitude = |p: Point3| {
            let d = p - centre;
            d.dot(axis).atan2(d.dot(x)) * to_deg
        };
        let (start_lat, end_lat) = (latitude(ends.start), latitude(ends.end));
        // the meridian rises from x towards the axis (x cross axis is -y)
        let rising = frame.z_axis().dot(-y) > 0.0;
        let (low, high, inside_out) = if rising {
            (start_lat, end_lat, false)
        } else {
            (end_lat, start_lat, true)
        };
        if high - low <= EFF_ZERO {
            return Ok(None);
        }
        let mut sphere = Sphere::new(
            k,
            centre,
            x,
            y,
            Extent2d::from_bounds(0.0, low, angle_deg, high),
            radius,
            SphereOptions {
                swap_uv: false,
                inside_out,
                circle_nurb: None,
            },
        )?;
        sphere.make_nurb(k)?;
        return Ok(Some(Surface::from(sphere)));
    }

    // a torus: its minor circle is the arc's, centred off the axis
    let Some(x) = radial(base_pt, axis, centre, tol) else {
        return Ok(None);
    };
    let y = axis.cross(x);
    if !in_half_plane(k, curve, base_pt, x, y, tol)? {
        return Ok(None);
    }
    let origin = base_pt + axis * line_closest_point(base_pt, axis, centre);
    let minor_angle = |p: Point3| {
        let d = p - centre;
        d.dot(axis).atan2(d.dot(x)) * to_deg
    };
    let rising = frame.z_axis().dot(-y) > 0.0;
    // the tube angle range starts in [0, 360): one starting below 0 is refused or
    // built wrong upstream (tangency-labs/solidcore#386)
    let (low, inside_out) = if rising {
        (turn_into_range(minor_angle(ends.start)), false)
    } else {
        (turn_into_range(minor_angle(ends.end)), true)
    };
    let mut torus = Torus::new(
        k,
        origin,
        x,
        y,
        Extent2d::from_bounds(0.0, low, angle_deg, low + sweep),
        centre_dist,
        radius,
        TorusOptions {
            swap_uv: false,
            inside_out,
            circle_nurb: None,
        },
    )?;
    torus.make_nurb(k)?;
    Ok(Some(Surface::from(torus)))
}

/// What a surface of revolution's natural domain shows of its layout.
#[derive(Debug, Clone, Copy)]
struct Layout {
    /// The points at the corners `(umin, vmin)`, `(umin, vmax)`, `(umax, vmin)`,
    /// `(umax, vmax)`.
    corners: [Point3; 4],
    /// The tangent along `v` at `(umin, vmin)`.
    start_tangent: Vec3,
    /// The sign of the `u` tangent against the rotation about the axis, at the
    /// domain's centre.
    turn: f64,
}

impl Layout {
    /// The layout of `surface` about the axis through `base_pt` along unit `axis`.
    fn of(k: &Kernel, surface: &Surface, base_pt: Point3, axis: Vec3) -> Result<Self> {
        let d = surface.natural_uv_domain();
        let mut corners = [Point3::ZERO; 4];
        for (corner, (u, v)) in
            corners
                .iter_mut()
                .zip([(0.0, 0.0), (0.0, 1.0), (1.0, 0.0), (1.0, 1.0)])
        {
            *corner = surface.evaluate_point(k, d.evaluate(u, v))?;
        }
        let flags = EvaluateFlags::default();
        let start = surface.evaluate(k, d.evaluate(0.0, 0.0), 0, 1, flags)?;
        let start_tangent = start
            .first()
            .and_then(|row| row.get(1))
            .copied()
            .unwrap_or(Vec3::ZERO);
        let centre = surface.evaluate(k, d.evaluate(0.5, 0.5), 1, 0, flags)?;
        let point = centre.first().and_then(|row| row.first()).copied();
        let du = centre.get(1).and_then(|row| row.first()).copied();
        let turn = match (point, du) {
            (Some(p), Some(du)) => du.dot(axis.cross(p - base_pt)).signum(),
            _ => 0.0,
        };
        Ok(Self {
            corners,
            start_tangent,
            turn,
        })
    }

    /// The layout of the B-spline surface of revolution of `curve` by `angle_deg`
    /// degrees: `v` along the curve from its start, `u` turning positively.
    fn of_curve(
        k: &Kernel,
        curve: &Curve,
        base_pt: Point3,
        axis: Vec3,
        angle_deg: f64,
    ) -> Result<Self> {
        let interval = curve.natural_interval();
        let start = curve.evaluate(k, interval.min(), 1, false, true)?;
        let start_pt = start.first().copied().unwrap_or(Vec3::UNDEF);
        let start_tangent = start.get(1).copied().unwrap_or(Vec3::ZERO);
        let end_pt = curve.evaluate_point(k, interval.max())?;
        let mut turn = Axis2Placement::default();
        turn.init();
        if angle_deg < 360.0 - EFF_ZERO {
            turn.rotate_about_axis_at_point(angle_deg * (PI / 180.0), base_pt, axis);
        }
        Ok(Self {
            corners: [
                start_pt,
                end_pt,
                turn.transform_point(start_pt),
                turn.transform_point(end_pt),
            ],
            start_tangent,
            turn: 1.0,
        })
    }

    /// Whether this layout is `other`'s: the corners within `tol`, the start
    /// tangents the same way and the same, defined, sense of turn.
    fn matches(&self, other: &Self, tol: f64) -> bool {
        self.corners
            .iter()
            .zip(&other.corners)
            .all(|(a, b)| (*a - *b).length() <= tol)
            && self.start_tangent.dot(other.start_tangent) > 0.0
            && self.turn != 0.0
            && self.turn == other.turn
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::curve::curve::CurveKind;
    use crate::nurbs::b_spline_curve::{BSplineCurve, NurbCircleParam};
    use crate::nurbs::b_spline_surface::BSplineSurface;

    fn bspline_rev(k: &Kernel, curve: &Curve, angle: f64) -> Surface {
        let CurveKind::BSpline(b) = curve.kind() else {
            panic!("B-spline profile")
        };
        Surface::from(
            BSplineSurface::create_surf_of_revolution_with_line_vector(
                k,
                b,
                None,
                Vec3::ZERO,
                Vec3::Z,
                angle,
                false,
                None,
            )
            .unwrap(),
        )
    }

    fn line(a: Vec3, b: Vec3) -> Curve {
        Curve::from(BSplineCurve::create_line_segment(3, a, b).unwrap())
    }

    /// The arc of radius `r` about `centre` in the XZ plane from `start` to `end`
    /// degrees, measured from `x` towards `y`.
    fn arc(centre: Vec3, x: Vec3, y: Vec3, r: f64, start: f64, end: f64) -> Curve {
        Curve::from(
            BSplineCurve::create_circle_segment(
                3,
                &Axis2Placement::new(centre, x, y),
                r,
                start,
                end,
                NurbCircleParam::Quadratic,
            )
            .unwrap(),
        )
    }

    fn kind(s: &Surface) -> &'static str {
        use crate::surface::surface::SurfaceKind;
        match s.kind() {
            SurfaceKind::Cone(c) if c.is_cylinder_class() => "cylinder",
            SurfaceKind::Cone(_) => "cone",
            SurfaceKind::Sphere(_) => "sphere",
            SurfaceKind::Torus(_) => "torus",
            _ => "other",
        }
    }

    /// Every profile kind, both ways along the curve: (name, curve, expected kind).
    fn profiles() -> Vec<(&'static str, Curve, &'static str)> {
        let mut out = vec![
            (
                "cylinder up",
                line(Vec3::new(2.0, 0.0, 0.0), Vec3::new(2.0, 0.0, 3.0)),
                "cylinder",
            ),
            (
                "cylinder down",
                line(Vec3::new(2.0, 0.0, 3.0), Vec3::new(2.0, 0.0, 0.0)),
                "cylinder",
            ),
            (
                "cone up",
                line(Vec3::new(1.0, 0.0, 0.0), Vec3::new(2.0, 0.0, 1.0)),
                "cone",
            ),
            (
                "cone down",
                line(Vec3::new(2.0, 0.0, 1.0), Vec3::new(1.0, 0.0, 0.0)),
                "cone",
            ),
            (
                "drill point",
                line(Vec3::new(1.0, 0.0, 1.0), Vec3::new(0.0, 0.0, 0.0)),
                "cone",
            ),
            (
                "apex first",
                line(Vec3::new(0.0, 0.0, 2.0), Vec3::new(1.5, 0.0, 0.0)),
                "cone",
            ),
        ];
        // the partial sphere's half circle: x = -Z, z = Y, pole to pole
        out.push((
            "half circle",
            arc(Vec3::ZERO, -Vec3::Z, Vec3::X, 3.0, 0.0, 180.0),
            "sphere",
        ));
        out.push((
            "half circle down",
            arc(Vec3::ZERO, Vec3::Z, Vec3::X, 3.0, 0.0, 180.0),
            "sphere",
        ));
        out.push((
            "dome",
            arc(Vec3::ZERO, Vec3::X, Vec3::Z, 2.0, 0.0, 90.0),
            "sphere",
        ));
        out.push((
            "dome down",
            arc(Vec3::ZERO, Vec3::Z, Vec3::X, 2.0, 0.0, 90.0),
            "sphere",
        ));
        // the partial torus's minor circle: x = -Z, about Y, from the bottom
        out.push((
            "minor circle",
            arc(Vec3::new(5.0, 0.0, 0.0), -Vec3::Z, Vec3::X, 1.5, 0.0, 360.0),
            "torus",
        ));
        out.push((
            "minor circle reversed",
            arc(
                Vec3::new(5.0, 0.0, 0.0),
                -Vec3::Z,
                -Vec3::X,
                1.5,
                0.0,
                360.0,
            ),
            "torus",
        ));
        out.push((
            "horn minor circle",
            arc(
                Vec3::new(5.0, 0.0, 0.0),
                -Vec3::Z,
                Vec3::X,
                5.0,
                -180.0,
                180.0,
            ),
            "torus",
        ));
        out.push((
            "fillet arc",
            arc(
                Vec3::new(4.0, 0.0, 1.0),
                Vec3::X,
                Vec3::Z,
                1.0,
                270.0,
                360.0,
            ),
            "torus",
        ));
        out.push((
            "fillet arc reversed",
            arc(Vec3::new(4.0, 0.0, 1.0), Vec3::Z, Vec3::X, 1.0, 0.0, 90.0),
            "torus",
        ));
        out
    }

    #[test]
    fn each_profile_curve_sweeps_its_analytic_surface_in_the_b_spline_layout() {
        let k = Kernel::default();
        for angle in [90.0, 200.0, 360.0] {
            for (name, curve, expected) in profiles() {
                let reference = bspline_rev(&k, &curve, angle);
                let from_curve = Layout::of_curve(&k, &curve, Vec3::ZERO, Vec3::Z, angle).unwrap();
                let b_layout = Layout::of(&k, &reference, Vec3::ZERO, Vec3::Z).unwrap();
                assert!(
                    b_layout.matches(&from_curve, 1e-9),
                    "{name} {angle}: B-spline {b_layout:?} curve {from_curve:?}"
                );
                for r in [None, Some(&reference)] {
                    let s = analytic_revolution(&k, &curve, Vec3::ZERO, Vec3::Z, angle, 1e-6, r)
                        .unwrap()
                        .unwrap_or_else(|| panic!("{name} {angle}: no analytic surface"));
                    assert_eq!(kind(&s), expected, "{name} {angle}");
                }
            }
        }
    }

    /// The distance from `p`, turned back into the half-plane `y = 0, x >= 0`, to the
    /// line or circle the profile lies on.
    fn off_profile(curve: &Curve, p: Vec3) -> f64 {
        let k = Kernel::default();
        let q = Vec3::new((p.x * p.x + p.y * p.y).sqrt(), 0.0, p.z);
        let ends = curve.ends(&k).unwrap();
        if let Some((frame, r, _, _)) = curve.is_arc(&k, 8, 1e-9).unwrap() {
            return ((q - frame.origin()).length() - r).abs();
        }
        let d = (ends.end - ends.start).unitize().unwrap();
        (q - ends.start).cross(d).length()
    }

    #[test]
    fn the_analytic_surface_is_the_swept_one() {
        let k = Kernel::default();
        for (name, curve, _) in profiles() {
            let reference = bspline_rev(&k, &curve, 120.0);
            let s = analytic_revolution(
                &k,
                &curve,
                Vec3::ZERO,
                Vec3::Z,
                120.0,
                1e-6,
                Some(&reference),
            )
            .unwrap()
            .unwrap();
            let a = s.natural_uv_domain();
            for i in 0..=6 {
                for j in 0..=6 {
                    let (f, g) = (f64::from(i) / 6.0, f64::from(j) / 6.0);
                    let p = s.evaluate_point(&k, a.evaluate(f, g)).unwrap();
                    assert!(off_profile(&curve, p) < 1e-9, "{name} {f} {g}: {p:?}");
                    let azimuth = p.y.atan2(p.x).to_degrees();
                    if p.x.hypot(p.y) > 1e-3 {
                        assert!(
                            (-1e-9..=120.0 + 1e-9).contains(&azimuth),
                            "{name}: {azimuth}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn other_curves_keep_the_upstream_surface() {
        let k = Kernel::default();
        let none = |c: &Curve| {
            analytic_revolution(&k, c, Vec3::ZERO, Vec3::Z, 90.0, 1e-6, None)
                .unwrap()
                .is_none()
        };
        // perpendicular to the axis (a plane), on the axis, skew to the axis
        assert!(none(&line(
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(2.0, 0.0, 0.0)
        )));
        assert!(none(&line(Vec3::ZERO, Vec3::new(0.0, 0.0, 2.0))));
        assert!(none(&line(
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(1.0, 1.0, 1.0)
        )));
        // an arc not in a plane through the axis
        assert!(none(&arc(
            Vec3::new(4.0, 0.0, 0.0),
            Vec3::X,
            Vec3::Y,
            1.0,
            0.0,
            90.0
        )));
        // an arc in a plane through the axis, crossing to its other side
        assert!(none(&arc(Vec3::ZERO, Vec3::X, Vec3::Z, 2.0, 0.0, 270.0)));
        // no turn, no axis
        let c = line(Vec3::new(2.0, 0.0, 0.0), Vec3::new(2.0, 0.0, 3.0));
        assert!(
            analytic_revolution(&k, &c, Vec3::ZERO, Vec3::Z, 0.0, 1e-6, None)
                .unwrap()
                .is_none()
        );
        assert!(
            analytic_revolution(&k, &c, Vec3::ZERO, Vec3::ZERO, 90.0, 1e-6, None)
                .unwrap()
                .is_none()
        );
    }
}
