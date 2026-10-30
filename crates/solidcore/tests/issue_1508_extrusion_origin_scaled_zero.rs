//! Regression tests for tangency-labs/solidcore#1508: `SmSurfOfExtrusion::AssertValid`
//! test 8 ("origin is a point on the generating curve") point-solves the origin with an
//! unscaled `SM_EFF_ZERO`, so an origin that the solve resolves a few 1e-12 off the
//! curve is reported missing at coordinates of about 10.
//!
//! The pinned C++ keeps that unscaled solve (and so does the port by default). With
//! `Kernel::fixes.extrusion_origin_scaled_zero` the solve takes
//! `SM_EFF_ZERO * (1.0 + origin.max_dimension())`, as the check's other point tests
//! scale theirs (WR-511, #1508; tangency-labs/smlib#31).
//!
//! These tests go through the production `Surface::assert_valid`, not the arithmetic of
//! the tolerance: the default path must still report test 8, the opt-in path must not,
//! and an origin genuinely off the generator must fail either way. On main this override
//! returns `Unported`, so the off/on comparison is inside this newly completed port; it
//! is not a pre-existing main result.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers

use solidcore::curve::curve::Curve;
use solidcore::geom::assert_array::{AssertArray, ObjectKey};
use solidcore::geom::extent2d::Extent2d;
use solidcore::nurbs::b_spline_curve::BSplineCurve;
use solidcore::nurbs::crv_basic::Curve as NurbCurve;
use solidcore::nurbs::math::CPoint;
use solidcore::surface::surf_of_extrusion::SurfOfExtrusion;
use solidcore::surface::surface::assert_valid::LEVEL_GIVEN;
use solidcore::surface::surface::Surface;
use solidcore::{Kernel, Vec3};

/// A planar quadratic Bezier `(0, 0) - (1, 1) - (2, 0)` translated by `offset` in x and
/// y (the text-outline shape of the issue).
fn quadratic_bezier(offset: f64) -> Curve {
    let pw = vec![
        CPoint::from_point(Vec3::new(offset, offset, 0.0)),
        CPoint::from_point(Vec3::new(offset + 1.0, offset + 1.0, 0.0)),
        CPoint::from_point(Vec3::new(offset + 2.0, offset, 0.0)),
    ];
    let nurb = NurbCurve::new(pw, 2, vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]).unwrap();
    Curve::from(BSplineCurve::from_nurb(nurb, 3, false))
}

/// The extrusion of [`quadratic_bezier`] at `offset` along z, with its origin at the
/// curve point 45.678% along its interval plus `origin_delta` (the solve's residual).
fn extrusion(k: &Kernel, offset: f64, origin_delta: Vec3) -> Surface {
    let gen = quadratic_bezier(offset);
    let origin = gen.evaluate_point(k, 0.45678).unwrap() + origin_delta;
    Surface::from(
        SurfOfExtrusion::new(
            k,
            gen,
            origin,
            Vec3::new(0.0, 0.0, 10.0),
            Extent2d::from_bounds(0.0, 0.0, 1.0, 1.0),
            false,
        )
        .unwrap(),
    )
}

/// `Surface::assert_valid` at `level`: whether every test passed and the failed test
/// numbers.
fn check(k: &Kernel, s: &Surface, level: u8) -> (bool, Vec<u64>) {
    let mut list = AssertArray::new();
    let ok = s
        .assert_valid(k, Some(&mut list), level, ObjectKey(1))
        .unwrap();
    let failed = list
        .reports()
        .iter()
        .filter(|r| !r.ok)
        .map(|r| r.test_index)
        .collect();
    (ok, failed)
}

fn set_fix(k: &mut Kernel, on: bool) {
    k.fixes.extrusion_origin_scaled_zero = on;
}

#[test]
fn ordinary_origin_is_valid_with_and_without_the_fix() {
    let mut k = Kernel::default();
    let s = extrusion(&k, 0.0, Vec3::ZERO);
    for on in [false, true] {
        set_fix(&mut k, on);
        assert_eq!(check(&k, &s, 0), (true, vec![]), "fix {on}");
    }
}

#[test]
fn issue_scale_origin_at_the_solver_residual_needs_the_scaled_zero() {
    // The C++ solve resolves an origin evaluated on this quadratic to a few 1e-12 at
    // coordinates of about 10 (#1508); model that residual with a 5e-12 out-of-plane
    // offset, which the unscaled SM_EFF_ZERO rejects and the scaled zero accepts.
    let mut k = Kernel::default();
    let s = extrusion(&k, 10.0, Vec3::new(0.0, 0.0, 5.0e-12));

    set_fix(&mut k, false);
    let default = check(&k, &s, 0);
    assert_eq!(default, (false, vec![8]), "default");

    set_fix(&mut k, true);
    let fixed = check(&k, &s, 0);
    assert_eq!(fixed, (true, vec![]), "fix");
    eprintln!("issue_1508 origin residual: default {default:?}, fix {fixed:?}");
}

#[test]
fn translated_quadratic_origin_is_valid_with_the_fix() {
    // An origin evaluated exactly on the curve, far from the origin of coordinates:
    // the solve's absolute residual is above the unscaled zero here, so only the scaled
    // zero can accept it. (The default's residual depends on the platform's rounding;
    // the assertion below allows it to be empty, but it must not report any other test.)
    let mut k = Kernel::default();
    let s = extrusion(&k, 1.0e8, Vec3::ZERO);

    set_fix(&mut k, false);
    let (_, default_failed) = check(&k, &s, 0);
    assert!(
        default_failed.iter().all(|&test| test == 8),
        "default reports only test 8: {default_failed:?}"
    );

    set_fix(&mut k, true);
    let fixed = check(&k, &s, 0);
    assert_eq!(fixed, (true, vec![]), "fix");
    eprintln!("issue_1508 translated origin: default {default_failed:?}, fix {fixed:?}");
}

#[test]
fn an_origin_off_the_generator_still_fails_with_the_fix() {
    let mut k = Kernel::default();
    let s = extrusion(&k, 10.0, Vec3::new(0.0, 0.0, 1.0e-6));
    for on in [false, true] {
        set_fix(&mut k, on);
        assert_eq!(check(&k, &s, 0), (false, vec![8]), "fix {on}");
    }
}

#[test]
fn a_made_nurb_extrusion_passes_at_every_level() {
    // The complete override: with a NURB the base checks, tests 5/6 (STEP <-> NURB) and
    // test 25 (generator parameterization) run too.
    let mut k = Kernel::default();
    let ext = SurfOfExtrusion::create_canonical(&k, quadratic_bezier(0.0), Vec3::Z).unwrap();
    let s = Surface::from(ext);
    for on in [false, true] {
        set_fix(&mut k, on);
        for level in [0, 1, 2, LEVEL_GIVEN] {
            assert_eq!(
                check(&k, &s, level),
                (true, vec![]),
                "fix {on} level {level}"
            );
        }
    }
}
