//! Regression for #1054 (WR-474): `SmCurve::CalculateTightBoundingBox` adds a B-spline's
//! tangent-discontinuous knots, and always its last knot, whatever the queried interval is,
//! so the tight box of an edge on part of a NURBS curve reaches points the edge never does.
//!
//! The fix is opt-in, `Kernel::fixes.curve_tight_box_interval_knots`: with it on, only the
//! knots strictly inside the interval are added; with it off (the default) the box keeps the
//! pinned upstream result. The first test is `PapaEcho13`'s from #1561, run both ways.

#![allow(clippy::unwrap_used)] // test helpers

use solidcore::curve::curve::Curve;
use solidcore::geom::axis2_placement::Axis2Placement;
use solidcore::geom::extent1d::Extent1d;
use solidcore::geom::extent3d::Extent3d;
use solidcore::nurbs::b_spline_curve::{BSplineCurve, NurbCircleParam};
use solidcore::{Kernel, Vec3};

fn arc(x_axis: [f64; 3], y_axis: [f64; 3], start_deg: f64, end_deg: f64) -> Curve {
    let frame = Axis2Placement::new(
        Vec3::from_array([0.0, 0.0, 0.0]),
        Vec3::from_array(x_axis),
        Vec3::from_array(y_axis),
    );
    Curve::from(
        BSplineCurve::create_circle_segment(
            3,
            &frame,
            50.0,
            start_deg,
            end_deg,
            NurbCircleParam::Quadratic,
        )
        .unwrap(),
    )
}

fn kernel(fixed: bool) -> Kernel {
    let mut k = Kernel::default();
    k.fixes.curve_tight_box_interval_knots = fixed;
    k
}

fn corners(bx: &Extent3d) -> [f64; 6] {
    [
        bx.u_min(),
        bx.v_min(),
        bx.w_min(),
        bx.u_max(),
        bx.v_max(),
        bx.w_max(),
    ]
}

fn assert_close(got: [f64; 6], want: [f64; 6]) {
    for (g, w) in got.iter().zip(want) {
        assert!((g - w).abs() < 1e-6, "box {got:?}, expected {want:?}");
    }
}

#[test]
fn half_quarter_circle_stops_at_its_own_45_degree_point() {
    // A quadratic NURBS quarter circle of radius 50 from (50, 0) to (0, 50), queried over
    // the first half of its interval (0 to 45 deg).
    let c = arc([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.0, 90.0);
    let ivl = c.natural_interval();
    let half = Extent1d::from_min_max(ivl.min(), 0.5 * (ivl.max() - ivl.min()) + ivl.min());
    let r = 50.0 / 2.0_f64.sqrt();

    let fixed = c
        .calculate_tight_bounding_box(&kernel(true), &half)
        .unwrap();
    assert_close(corners(&fixed), [r, 0.0, 0.0, 50.0, r, 0.0]);

    // Upstream (flag off) also adds the last knot, the 90 deg end (0, 50).
    let pinned = c
        .calculate_tight_bounding_box(&kernel(false), &half)
        .unwrap();
    assert_close(corners(&pinned), [0.0, 0.0, 0.0, 50.0, 50.0, 0.0]);
}

#[test]
fn sphere_cap_seam_box_stays_on_the_cap() {
    // The issue's seam: R 50, a cap from the pole down to polar angle 30 deg. The seam edge
    // is the first 30 deg of a quarter-circle NURBS in the xz-plane that runs from polar
    // angle -30 deg, (25, 0, 43.30), over the pole to polar angle 60 deg, (-43.30, 0, 25).
    let c = arc([0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], -30.0, 60.0);
    let k = kernel(false);
    let ivl = c.natural_interval();
    // The parameter of the pole (x = 0; x falls along the curve).
    let (mut lo, mut hi) = (ivl.min(), ivl.max());
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if c.evaluate_point(&k, mid).unwrap().x > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let pole = 0.5 * (lo + hi);
    assert!(
        (pole - ivl.min() - 0.341 * (ivl.max() - ivl.min())).abs() < 1e-3,
        "pole at {pole} of {ivl:?}"
    );
    let seam = Extent1d::from_min_max(ivl.min(), pole);
    let z30 = 50.0 * 30.0_f64.to_radians().cos();

    let fixed = c
        .calculate_tight_bounding_box(&kernel(true), &seam)
        .unwrap();
    assert_close(corners(&fixed), [0.0, 0.0, z30, 25.0, 0.0, 50.0]);

    // Upstream takes in the last knot, (-43.30, 0, 25), which the seam never reaches.
    let pinned = c.calculate_tight_bounding_box(&k, &seam).unwrap();
    assert_close(corners(&pinned), [-z30, 0.0, 25.0, 25.0, 0.0, 50.0]);
}
