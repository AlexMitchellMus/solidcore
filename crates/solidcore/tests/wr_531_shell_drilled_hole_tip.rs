//! Regression test for WR-531 (the fork's tangency-labs/smlib@d5dadaa, capi test
//! `shelled_box_with_a_drilled_hole`): a box with a blind drilled hole (r 3, a 118
//! degree tip) shelled 2 mm inward, open at a side. The hole's tip is a cone's apex
//! whose seam is its only edge: the extended offset builds no offset there, so the
//! cone's offset ends open on a circle and the shell fails. With
//! `offset_point_vertex_caps` the apex's offset is a sphere zone of radius 2 and the
//! shell is sound.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers

#[path = "support/offset_fixture.rs"]
mod offset_fixture;

use offset_fixture::{planar_faces_at, shell};
use solidcore::boolean::api_brep::boolean_difference;
use solidcore::curve::curve::Curve;
use solidcore::math::{Point3, Vec3};
use solidcore::primitive::api_curves::create_line_segment;
use solidcore::primitive::api_primitives::create_box;
use solidcore::surface::surface::SurfaceKind;
use solidcore::sweep::api_brep::create_rotational_sweep;
use solidcore::topology::brep::validation::{check_brep, CheckLevel};
use solidcore::{Brep, Kernel};

/// A 40 x 40 x 20 box with a hole r 3 drilled 10 deep from the top at (20, 20), its
/// tip a 118 degree cone.
fn drilled_box(k: &mut Kernel) -> Brep {
    let (r, top, depth) = (3.0, 20.0, 10.0);
    let tip = r / 59.0_f64.to_radians().tan();
    let bottom = top - depth;
    let at = |x: f64, z: f64| Point3::new(20.0 + x, 20.0, z);
    let profile = [
        (at(0.0, bottom - tip), at(r, bottom)),
        (at(r, bottom), at(r, top + 1.0)),
        (at(r, top + 1.0), at(0.0, top + 1.0)),
        (at(0.0, top + 1.0), at(0.0, bottom - tip)),
    ];
    let curves: Vec<Curve> = profile
        .iter()
        .map(|&(a, b)| Curve::from(create_line_segment(a, b).unwrap()))
        .collect();
    let curves: Vec<Option<&Curve>> = curves.iter().map(Some).collect();
    let drill = create_rotational_sweep(
        k,
        &curves,
        at(0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        360.0,
        true,
    )
    .unwrap();
    let block = create_box(k, Point3::new(0.0, 0.0, 0.0), 40.0, 40.0, 20.0).unwrap();
    boolean_difference(k, block, drill).unwrap()
}

fn shelled(fixed: bool) -> solidcore::Result<(Kernel, Brep)> {
    let mut k = Kernel::default();
    k.fixes.offset_point_vertex_caps = fixed;
    let part = drilled_box(&mut k);
    let side = planar_faces_at(&k, &part, 0, 0.0);
    assert_eq!(side.len(), 1, "the side x = 0");
    let result = shell(&mut k, part, -2.0, &side)?;
    Ok((k, result))
}

#[test]
fn wr_531_drilled_hole_shells_with_a_sphere_cap() {
    let (mut k, mut b) = shelled(true).unwrap();
    let c = check_brep(&mut k, &mut b, CheckLevel::Level2).unwrap();
    assert_eq!(c.failure_count(), 0, "{:?}", c.lines());
    let spheres = b
        .face_list()
        .iter()
        .filter(|&&f| {
            b.face(f)
                .unwrap()
                .surface
                .as_ref()
                .is_some_and(|s| matches!(s.kind(), SurfaceKind::Sphere(_)))
        })
        .count();
    assert_eq!(spheres, 1, "one cap on a sphere");
}

/// Without the opt-in the shell fails, as upstream does.
#[test]
fn wr_531_without_the_opt_in_the_shell_fails() {
    assert!(shelled(false).is_err());
}
