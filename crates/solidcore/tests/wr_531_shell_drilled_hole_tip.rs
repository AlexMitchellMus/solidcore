//! Regression test for WR-531 (the fork's tangency-labs/smlib@d5dadaa, capi test
//! `shelled_box_with_a_drilled_hole`): a box with a blind drilled hole (r 3, a 118
//! degree tip) shelled 2 mm inward, open at a side. The hole's tip is a cone's apex
//! whose seam is its only edge: the extended offset builds no offset there, so the
//! cone's offset ends open on a circle and the offset is no closed solid. With
//! `offset_point_vertex_caps` the apex's offset is a sphere zone of radius 2 and the
//! shell is a sound closed solid.
//!
//! Without the opt-in the pinned `SmOffsetExecutive::ShellBrep` hands that open offset
//! back as a success (`SmOffsetExecutive.cpp:255-262`: not a manifold solid and no
//! offset solid asked for, `return SM_SUCCESS; // was SM_ERR: ??`). The fork's "the
//! shell fails" is that open offset refused by its WR-454 check
//! (`shell_empty_cavity_fails`, tangency-labs/smlib@fbe43bb), so that case runs with it.
//!
//! The result is checked with the geometry owners the body itself records: a shell
//! result names its edges' curves by the caller's merge and result identities while
//! being labelled `NamingScheme::Notify` (#1670, open), so `check_brep`'s owner
//! comparison in test 61 fails on every such edge of any merged shell. Every other
//! check, and the rest of test 61 (each edge's ends on its vertices), still runs.
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
use solidcore::geom::attribute::ObjectKey;
use solidcore::geom::object::ObjectRef;
use solidcore::math::{Point3, Vec3};
use solidcore::primitive::api_curves::create_line_segment;
use solidcore::primitive::api_primitives::create_box;
use solidcore::surface::surface::SurfaceKind;
use solidcore::sweep::api_brep::create_rotational_sweep;
use solidcore::topology::brep::p2::is_manifold_solid;
use solidcore::topology::brep::validation::{check_brep_with_keys, CheckLevel};
use solidcore::{Brep, EntityRef, Kernel};
use std::collections::BTreeMap;

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

/// The drilled box shelled 2 inward, open at x = 0, with WR-531's `caps` and WR-454's
/// `refuse_open`.
fn shelled(caps: bool, refuse_open: bool) -> solidcore::Result<(Kernel, Brep)> {
    let mut k = Kernel::default();
    k.fixes.offset_point_vertex_caps = caps;
    k.fixes.shell_empty_cavity_fails = refuse_open;
    let part = drilled_box(&mut k);
    let side = planar_faces_at(&k, &part, 0, 0.0);
    assert_eq!(side.len(), 1, "the side x = 0");
    let result = shell(&mut k, part, -2.0, &side)?;
    Ok((k, result))
}

/// The geometry owners `b` records: each edge's curve's, face's surface's and edge
/// use's UV trim curve's (the only entities the check names as owners).
fn recorded_owners(b: &Brep) -> BTreeMap<EntityRef, ObjectRef> {
    let mut owners = BTreeMap::new();
    for e in b.edges().ids() {
        if let Some(o) = b.edge(e).unwrap().curve.as_ref().and_then(Curve::owner) {
            owners.insert(EntityRef::Edge(e), o);
        }
    }
    for f in b.faces().ids() {
        if let Some(o) = b.face(f).unwrap().surface.as_ref().and_then(|s| s.owner()) {
            owners.insert(EntityRef::Face(f), o);
        }
    }
    for eu in b.edgeuses().ids() {
        if let Some(o) = b
            .edgeuse(eu)
            .unwrap()
            .uv_trim_curve
            .as_ref()
            .and_then(Curve::owner)
        {
            owners.insert(EntityRef::Edgeuse(eu), o);
        }
    }
    owners
}

/// The level-2 check's failure count and report lines under the owners the body
/// records (#1670).
fn check(k: &mut Kernel, b: &mut Brep) -> (usize, Vec<String>) {
    let owners = recorded_owners(b);
    let none = ObjectRef {
        type_code: 0,
        key: ObjectKey(0),
    };
    let keys = |e: EntityRef| owners.get(&e).copied().unwrap_or(none);
    let c = check_brep_with_keys(k, b, CheckLevel::Level2, &keys).unwrap();
    (c.failure_count(), c.lines())
}

fn spheres(b: &Brep) -> usize {
    b.face_list()
        .iter()
        .filter(|&&f| {
            b.face(f)
                .unwrap()
                .surface
                .as_ref()
                .is_some_and(|s| matches!(s.kind(), SurfaceKind::Sphere(_)))
        })
        .count()
}

#[test]
fn wr_531_drilled_hole_shells_with_a_sphere_cap() {
    for refuse_open in [false, true] {
        let (mut k, mut b) = shelled(true, refuse_open).unwrap();
        assert!(is_manifold_solid(&b).unwrap(), "refuse_open {refuse_open}");
        let (failed, lines) = check(&mut k, &mut b);
        assert_eq!(failed, 0, "refuse_open {refuse_open}: {lines:?}");
        assert_eq!(spheres(&b), 1, "one cap on a sphere");
    }
}

/// Without the opt-in the offset is left open at the apex: the pinned `ShellBrep` hands
/// it back as a success, and WR-454's check refuses it (the fork's failing shell).
#[test]
fn wr_531_without_the_opt_in_the_shell_is_left_open() {
    let (_, b) = shelled(false, false).unwrap();
    assert!(!is_manifold_solid(&b).unwrap(), "open at the apex");
    assert_eq!(spheres(&b), 0);
    assert!(shelled(false, true).is_err());
}
