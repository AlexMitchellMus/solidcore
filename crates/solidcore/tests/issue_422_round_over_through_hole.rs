//! tangency-labs/solidcore#422 (WR-212): an R4 round over a through hole returns success with
//! a body that is not a closed manifold solid (SolidCraft rig case fillets/324). The round
//! would have to be trimmed by the hole's wall, which SMLib does not do. The opt-in
//! `Kernel::fixes.fillet_refuse_non_manifold_result` refuses the result instead, as the
//! issue's "Expected" allows ("a refusal naming the edge that can't be blended"). The
//! fork's trim by the hole's wall (tangency-labs/smlib@3dbe54a) is not ported.
//!
//! The case is the rig's: a 40 x 20 x 10 block minus a through hole of radius 3 about
//! (20, 14.5), rounded R4 along its top long edge at y = 20. The round covers y = 16..20
//! and the hole reaches y = 17.5.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers

use solidcore::blend::fillet_executive::api::fillet_edges;
use solidcore::boolean::api_brep::boolean_difference;
use solidcore::math::Point3;
use solidcore::primitive::api_primitives::{create_box, create_cylinder};
use solidcore::topology::brep::p2::is_manifold_solid;
use solidcore::topology::{edge, vertex, EdgeId};
use solidcore::{Brep, EntityRef, Kernel, Result, SolidError};

/// The block with its hole, and its top long edge at y = 20.
fn block(k: &mut Kernel) -> (Brep, EdgeId) {
    let b = create_box(k, Point3::new(0.0, 0.0, 0.0), 40.0, 20.0, 10.0).unwrap();
    let hole = create_cylinder(k, Point3::new(20.0, 14.5, -1.0), 3.0, 12.0).unwrap();
    let b = boolean_difference(k, b, hole).unwrap();
    let top: Vec<EdgeId> = b
        .edges()
        .ids()
        .filter(|&e| {
            let (s, t) = edge::start_end_vertices(&b, e).unwrap();
            let (Some(s), Some(t)) = (s, t) else {
                return false;
            };
            let (p, q) = (vertex::point(&b, s).unwrap(), vertex::point(&b, t).unwrap());
            [p, q]
                .iter()
                .all(|v| (v.y - 20.0).abs() < 1e-9 && (v.z - 10.0).abs() < 1e-9)
                && (p.x - q.x).abs() > 39.0
        })
        .collect();
    assert_eq!(top.len(), 1);
    (b, top[0])
}

/// The round of radius `r` with the fix as `fix` says: its status and whether the block is
/// a closed manifold solid after it.
fn round(fix: bool, r: f64) -> (Result<()>, bool) {
    let mut k = Kernel::default();
    k.fixes.fillet_refuse_non_manifold_result = fix;
    let (mut b, top) = block(&mut k);
    let result = fillet_edges(&mut k, &mut b, &[top], r, 1, 1, 1.0);
    let solid = is_manifold_solid(&b).unwrap();
    (result, solid)
}

#[test]
fn issue_422_flag_is_off_by_default() {
    assert!(!Kernel::default().fixes.fillet_refuse_non_manifold_result);
}

/// Upstream: the R4 round reports success and leaves a body that is not a closed solid.
#[test]
fn issue_422_upstream_returns_an_open_body_as_success() {
    let (result, solid) = round(false, 4.0);
    result.unwrap();
    assert!(!solid);
}

/// With the fix the R4 round is refused with `SM_FILERR_BAD_MERGE` (1001), naming the
/// edge of the result that is not manifold. (The refusal comes after the merge, so the
/// target is left as the merge made it.)
#[test]
fn issue_422_round_into_the_hole_is_refused_with_the_fix() {
    let (result, _) = round(true, 4.0);
    match result {
        Err(SolidError::Upstream {
            code: 1001,
            entity: Some(EntityRef::Edge(_)),
            ..
        }) => {}
        other => panic!("not refused naming an edge: {other:?}"),
    }
}

/// A round clear of the hole (R2, the hole's wall at y = 17.5) still builds a solid.
#[test]
fn issue_422_round_clear_of_the_hole_is_kept() {
    for fix in [false, true] {
        let (result, solid) = round(fix, 2.0);
        result.unwrap();
        assert!(solid, "fix {fix}");
    }
}
