//! #1105: whole-trim adaptive steps can miss features between UV spline knots.
//! The opt-in fix keeps the requested accuracy and integrates each knot span.
#![allow(clippy::unwrap_used, clippy::indexing_slicing)] // test fixtures

use solidcore::curve::{curve::Curve, iso_curve::IsoCurve, line::Line};
use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::geom::extent2d::Extent2d;
use solidcore::mass::integrator::properties::{
    boundary_properties, brep_properties, MassPropertyBehavior, PropertyAttributes,
};
use solidcore::math::{Point2, Vec2, Vec3};
use solidcore::nurbs::b_spline_curve::p1::BSplineCurveForm;
use solidcore::nurbs::b_spline_surface::p1::SurfParamType;
use solidcore::nurbs::b_spline_surface_n_lib::p2::OrientType;
use solidcore::primitive::api_curves::{
    create_canonical_curve, create_line_segment, CanonicalCurveParams,
};
use solidcore::primitive::api_primitives::create_sphere_no_pole_with_tolerance;
use solidcore::query::brep::{calculate_bounding_box, BoundingBoxOptions};
use solidcore::surface::{
    curve_bounded_surface::CurveBoundedSurface, plane::Plane, surface::Surface,
};
use solidcore::topology::face::p2::{PPF_ALL, PPF_AREA};
use solidcore::Kernel;

#[test]
fn narrow_trim_feature_is_missed_by_default_and_integrated_when_enabled() {
    let domain = Extent2d::from_points(Point2::new(0.0, 0.0), Point2::new(1.0, 1.1));
    let surface = Surface::from(Plane::new(
        Vec3::ZERO,
        Vec3::X,
        Vec3::Y,
        Vec2::new(1.0, 1.0),
        &domain,
    ));
    // A rectangle of area 0.1 with a narrow triangular feature of area 0.005.
    // The whole-curve Cash-Karp stages miss its support, 0.61 < t < 0.62.
    let points = [
        Vec3::new(1.0, 0.1, 0.0),
        Vec3::new(0.39, 0.1, 0.0),
        Vec3::new(0.385, 1.1, 0.0),
        Vec3::new(0.38, 0.1, 0.0),
        Vec3::new(0.0, 0.1, 0.0),
    ];
    let trims = [
        Curve::from(
            create_canonical_curve(&CanonicalCurveParams {
                points: &points,
                knots: &[0.0, 0.61, 0.615, 0.62, 1.0],
                knot_multiplicities: &[2, 1, 1, 1, 2],
                weights: None,
                degree: 1,
                form: BSplineCurveForm::Unspecified,
            })
            .unwrap(),
        ),
        Curve::from(create_line_segment(points[4], Vec3::ZERO).unwrap()),
        Curve::from(create_line_segment(Vec3::ZERO, Vec3::X).unwrap()),
        Curve::from(create_line_segment(Vec3::X, points[0]).unwrap()),
    ];
    let mut k = Kernel::default();
    assert!(!k.fixes.precise_properties_at_trim_knots);
    for accuracy in [1e-10, 1e-4] {
        for enabled in [false, true] {
            k.fixes.precise_properties_at_trim_knots = enabled;
            for flags in [PPF_AREA, PPF_ALL] {
                let p = boundary_properties(
                    &k,
                    &surface,
                    &trims,
                    &domain,
                    OrientType::Same,
                    accuracy,
                    Vec3::ZERO,
                    1.0,
                    0.0,
                    flags,
                )
                .unwrap();
                if enabled {
                    assert!((p.area - 0.105).abs() <= 1e-9 * 0.105, "{p:?}");
                } else {
                    assert!((p.area - 0.1).abs() <= 1e-9 * 0.1, "{p:?}");
                }
            }
        }
    }
}

#[test]
fn opt_in_integrates_whole_trim_when_knots_are_unavailable() {
    let domain = Extent2d::from_points(Point2::new(0.0, 0.0), Point2::new(1.0, 1.0));
    let surface = Surface::from(Plane::new(
        Vec3::ZERO,
        Vec3::X,
        Vec3::Y,
        Vec2::new(1.0, 1.0),
        &domain,
    ));
    let mut k = Kernel::default();
    // This surface evaluates through its plane, but its base knot query returns
    // SM_ERR. Its iso curve is the top edge of the square, from right to left.
    let uv_source = Surface::from(
        CurveBoundedSurface::new(
            Surface::from(Plane::new(
                Vec3::new(1.0, 1.0, 0.0),
                -Vec3::X,
                Vec3::Y,
                Vec2::new(1.0, 1.0),
                &domain,
            )),
            Vec::new(),
            Vec::new(),
            true,
        )
        .unwrap(),
    );
    let top = Curve::from(
        IsoCurve::from_surface(uv_source, SurfParamType::V, 0.0, false, false).unwrap(),
    );
    assert!(top.knots(&k, Some(&top.natural_interval())).is_err());
    let trims = [
        Curve::from(Line::from_points(Vec3::ZERO, Vec3::X, 3)),
        Curve::from(Line::from_points(Vec3::X, Vec3::new(1.0, 1.0, 0.0), 3)),
        top,
        Curve::from(Line::from_points(Vec3::Y, Vec3::ZERO, 3)),
    ];
    let base = boundary_properties(
        &k,
        &surface,
        &trims,
        &domain,
        OrientType::Same,
        1e-4,
        Vec3::ZERO,
        1.0,
        0.0,
        PPF_AREA,
    )
    .unwrap();
    k.fixes.precise_properties_at_trim_knots = true;
    let fixed = boundary_properties(
        &k,
        &surface,
        &trims,
        &domain,
        OrientType::Same,
        1e-4,
        Vec3::ZERO,
        1.0,
        0.0,
        PPF_AREA,
    )
    .unwrap();
    assert!((base.area - 1.0).abs() < 1e-9, "{base:?}");
    assert_eq!(fixed.area, base.area);
}

#[test]
fn pole_free_spheres_keep_analytic_area_and_volume_with_trim_span_integration() {
    for (center, radius) in [(Vec3::ZERO, 2.0), (Vec3::new(1.0, -2.0, 0.5), 5.0)] {
        let mut k = Kernel::default();
        k.fixes.precise_properties_at_trim_knots = true;
        let mut b = create_sphere_no_pole_with_tolerance(&mut k, center, radius, 1e-9).unwrap();
        let bbox = calculate_bounding_box(&k, &b, BoundingBoxOptions::default()).unwrap();
        let store = AttributeStore::default();
        let key_of = |_| ObjectKey(0);
        let p = brep_properties(
            &k,
            &mut b,
            1e-10,
            bbox.mid(),
            MassPropertyBehavior::default(),
            &PropertyAttributes {
                store: &store,
                key_of: &key_of,
            },
        )
        .unwrap();
        let area = 4.0 * std::f64::consts::PI * radius * radius;
        let volume = area * radius / 3.0;
        assert!(
            (p.area - area).abs() <= 1e-9 * area,
            "radius {radius}: {p:?}"
        );
        assert!(
            (p.volume - volume).abs() <= 1e-9 * volume,
            "radius {radius}: {p:?}"
        );
    }
}
