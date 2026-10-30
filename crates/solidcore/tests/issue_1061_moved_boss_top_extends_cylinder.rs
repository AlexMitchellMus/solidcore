//! tangency-labs/solidcore#1061: a local extended offset (`SmBrep::LocalOperation`, the
//! capi's `smlib_offset_faces`) that moves a boss's top face up extends the boss's wall
//! cylinder by a tangent (G1R) extension, which leaves the cylinder at once: the body
//! came out 2.6e-4 relative off its closed-form volume. With
//! `offset_extends_cones_as_themselves` (WR-532, the fork's tangency-labs/smlib@b283ea2)
//! the wall is extended as the cylinder it is.
//!
//! The part: a 60 x 40 x 10 plate with an r10 boss on it at (30, 20), its top at z = 20,
//! moved to z = 22: the plate plus a cylinder r10 from z = 10 to 22.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers

#[path = "support/offset_fixture.rs"]
mod offset_fixture;

use offset_fixture::{planar_faces_at, volume, Notifications};
use solidcore::boolean::api_brep::boolean_union;
use solidcore::math::{Point3, Vec3, PI};
use solidcore::primitive::api_primitives::{create_box, create_cylinder};
use solidcore::surface::plane::Plane;
use solidcore::surface::surface::Surface;
use solidcore::topology::brep::local_operation::local_operation;
use solidcore::Kernel;

/// The moved part's volume relative to its closed form, `(got - want) / want`.
fn moved_top(fixed: bool) -> f64 {
    let mut k = Kernel::default();
    k.fixes.offset_extends_cones_as_themselves = fixed;
    let plate = create_box(&mut k, Point3::new(0.0, 0.0, 0.0), 60.0, 40.0, 10.0).unwrap();
    let boss = create_cylinder(&mut k, Point3::new(30.0, 20.0, 5.0), 10.0, 15.0).unwrap();
    let mut part = boolean_union(&mut k, plate, boss).unwrap();
    let top = planar_faces_at(&k, &part, 2, 20.0);
    assert_eq!(top.len(), 1, "the boss's top face");
    let moved = Surface::from(Plane::from_normal(Vec3::new(30.0, 20.0, 22.0), Vec3::Z));
    let mut notifications = Notifications::new(&mut part);
    local_operation(
        &mut k,
        &mut part,
        &top,
        vec![moved],
        &mut notifications.context(),
    )
    .unwrap();
    solidcore::topology::integrity::check(&part).unwrap();
    let want = 60.0 * 40.0 * 10.0 + PI * 100.0 * 12.0;
    let got = volume(&k, &mut part).unwrap();
    eprintln!("fixed {fixed}: volume {got}, closed form {want}");
    (got - want) / want
}

#[test]
fn issue_1061_moved_boss_top_keeps_the_wall_a_cylinder() {
    let rel = moved_top(true);
    assert!(rel.abs() < 1e-9, "{rel}");
}

/// Without the opt-in the wall is extended approximately, as upstream does.
#[test]
fn issue_1061_without_the_opt_in_the_wall_is_extended_approximately() {
    let rel = moved_top(false);
    assert!(rel.abs() > 1e-7, "{rel}");
}
