//! Regression test for #1473: `SmShell::IsInnerShellOf` on a box with a sealed internal
//! void (the issue's repro: a 20 mm cube minus a 10 mm cube at (5, 5, 5)).

use solidcore::boolean::api_brep as sm;
use solidcore::geom::a_object::BooleanOperation;
use solidcore::kernel::Kernel;
use solidcore::math::Vec3;
use solidcore::primitive::api_primitives::create_box;
use solidcore::topology::shell::{is_inner_shell_of, SM_FALSE, SM_TRUE};
use solidcore::topology::{brep::p1, Brep, ShellId};

/// The largest |coordinate| of the shell's face vertices: 20 for the outer wall, 15 for
/// the cavity's.
fn shell_extent(b: &Brep, s: ShellId) -> f64 {
    let mut m: f64 = 0.0;
    for &f in &b.shell(s).unwrap().faces {
        for &l in &b.face(f).unwrap().loops {
            for le in &b.loop_(l).unwrap().edges {
                for v in b.edge(le.edge).unwrap().vertices {
                    let p = b.vertex(v).unwrap().point;
                    m = m.max(p.x.abs()).max(p.y.abs()).max(p.z.abs());
                }
            }
        }
    }
    m
}

fn void_body(k: &mut Kernel) -> Brep {
    let outer = create_box(k, Vec3::new(0.0, 0.0, 0.0), 20.0, 20.0, 20.0).unwrap();
    let cavity = create_box(k, Vec3::new(5.0, 5.0, 5.0), 10.0, 10.0, 10.0).unwrap();
    sm::boolean(k, outer, cavity, BooleanOperation::Difference).unwrap().result
}

fn check(k: &mut Kernel) {
    let mut b = void_body(k);
    let regions = p1::regions(&b, None).unwrap();
    let infinite = b.infinite_region();
    let mut solid_checked = false;
    for r in regions {
        let shells = b.region(r).unwrap().shells.clone();
        if Some(r) == infinite {
            // the infinite region's one shell (the 20 mm cube's outside) is an inner
            // shell of the infinite region
            assert_eq!(shells.len(), 1);
            assert_eq!(is_inner_shell_of(k, &mut b, shells[0], None).unwrap(), SM_TRUE);
            continue;
        }
        if shells.len() == 1 {
            // the void: another shell contains it
            assert_eq!(shell_extent(&b, shells[0]), 15.0);
            assert_eq!(is_inner_shell_of(k, &mut b, shells[0], None).unwrap(), SM_FALSE);
            continue;
        }
        // the solid region: the outer wall encloses the cavity's wall
        assert_eq!(shells.len(), 2);
        let (outer, cavity) = if shell_extent(&b, shells[0]) > shell_extent(&b, shells[1]) {
            (shells[0], shells[1])
        } else {
            (shells[1], shells[0])
        };
        assert_eq!(shell_extent(&b, outer), 20.0);
        assert_eq!(is_inner_shell_of(k, &mut b, cavity, Some(outer)).unwrap(), SM_TRUE);
        assert_eq!(is_inner_shell_of(k, &mut b, outer, Some(cavity)).unwrap(), SM_FALSE);
        solid_checked = true;
    }
    assert!(solid_checked, "no region with two shells");
}

#[test]
fn issue_1473_box_with_sealed_void_upstream_statuses() {
    let mut k = Kernel::default();
    check(&mut k);
}

#[test]
fn issue_1473_box_with_sealed_void_fork_statuses() {
    let mut k = Kernel::default();
    k.fixes.shell_inner_shell_of_statuses = true;
    check(&mut k);
}
