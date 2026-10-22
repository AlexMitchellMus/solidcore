//! Builds primitive solids through the public API (pure Rust, no FFI) and checks that
//! each is a real, usable solid: a closed manifold, the Euler-Poincare characteristic,
//! volume and area against the analytic formulas, and its bounding box.
//!
//! Run with `cargo run -p solidcore --example primitives`. No tessellator is ported
//! yet, so no STL is written.
use std::f64::consts::PI;

use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::geom::vector3d::Point3;
use solidcore::kernel::Kernel;
use solidcore::mass::integrator::properties::{brep_area, brep_volume, PropertyAttributes};
use solidcore::primitive::api_primitives as api;
use solidcore::query::api_queries::{brep_bounding_box, brep_is_manifold_solid};
use solidcore::topology::Brep;
use solidcore::{Result, SolidError};

/// Integration accuracy handed to the mass-property integrator (its finest).
const ACCURACY: f64 = 1e-8;
/// Relative tolerance of the analytic volume and area comparisons.
const REL_TOL: f64 = 1e-9;
/// Relative tolerance (of the box size) of the bounding-box comparison.
const BOX_TOL: f64 = 1e-6;

/// What a primitive must measure.
#[derive(Debug)]
struct Expected {
    /// Display name.
    name: &'static str,
    /// `2 * (shells - genus)`: 2 for a ball, 0 for a solid torus.
    euler: i64,
    /// Analytic volume.
    volume: f64,
    /// Analytic surface area.
    area: f64,
    /// Analytic bounding box (min, max).
    bbox: ([f64; 3], [f64; 3]),
}

/// What a primitive measured.
#[derive(Debug)]
struct Report {
    manifold: bool,
    v: i64,
    e: i64,
    f: i64,
    l: i64,
    volume: f64,
    area: f64,
    bbox: ([f64; 3], [f64; 3]),
}

impl Report {
    /// `V - E + F`.
    fn vef(&self) -> i64 {
        self.v - self.e + self.f
    }

    /// Euler-Poincare `V - E + F - (L - F)` (equals `V - E + F` when every face has
    /// one loop).
    fn euler(&self) -> i64 {
        self.vef() - (self.l - self.f)
    }
}

fn rel(a: f64, b: f64) -> f64 {
    (a - b).abs() / a.abs().max(b.abs())
}

fn count(n: usize) -> i64 {
    i64::try_from(n).unwrap_or(i64::MAX)
}

fn measure(k: &Kernel, b: &Brep) -> Result<Report> {
    let store = AttributeStore::default();
    let attrs = PropertyAttributes {
        store: &store,
        key_of: &|_| ObjectKey(0),
    };
    let ext = brep_bounding_box(k, b, true)?;
    Ok(Report {
        manifold: brep_is_manifold_solid(b)?,
        v: count(b.vertices().len()),
        e: count(b.edges().len()),
        f: count(b.faces().len()),
        l: count(b.loops().len()),
        volume: brep_volume(k, b, ACCURACY, &attrs)?,
        area: brep_area(k, b, ACCURACY)?,
        bbox: (
            [ext.min_at(0), ext.min_at(1), ext.min_at(2)],
            [ext.max_at(0), ext.max_at(1), ext.max_at(2)],
        ),
    })
}

/// The failed checks of `r` against `x` (empty when all pass).
fn failures(r: &Report, x: &Expected) -> Vec<String> {
    let mut out = Vec::new();
    if !r.manifold {
        out.push("not a closed manifold solid".to_string());
    }
    if r.euler() != x.euler {
        out.push(format!("Euler-Poincare {} != {}", r.euler(), x.euler));
    }
    if rel(r.volume, x.volume) > REL_TOL {
        out.push(format!("volume rel error {:e}", rel(r.volume, x.volume)));
    }
    if rel(r.area, x.area) > REL_TOL {
        out.push(format!("area rel error {:e}", rel(r.area, x.area)));
    }
    let axes = r
        .bbox
        .0
        .iter()
        .zip(&r.bbox.1)
        .zip(x.bbox.0.iter().zip(&x.bbox.1));
    for (i, ((lo, hi), (xlo, xhi))) in axes.enumerate() {
        let scale = xhi - xlo;
        if (lo - xlo).abs() > BOX_TOL * scale || (hi - xhi).abs() > BOX_TOL * scale {
            out.push(format!("bbox axis {i}: [{lo}, {hi}] != [{xlo}, {xhi}]"));
        }
    }
    out
}

/// Sphere centre and radius used by the sphere cases.
const SPHERE: ([f64; 3], f64) = ([1.0, -2.0, 0.5], 3.0);
/// Torus centre, major and minor radius used by the torus cases.
const TORUS: ([f64; 3], f64, f64) = ([0.0, 1.0, -1.0], 5.0, 1.5);

fn pt(c: [f64; 3]) -> Point3 {
    Point3::new(c[0], c[1], c[2])
}

fn sphere_expected(name: &'static str) -> Expected {
    let (c, r) = SPHERE;
    Expected {
        name,
        euler: 2,
        volume: 4.0 / 3.0 * PI * r * r * r,
        area: 4.0 * PI * r * r,
        bbox: (
            [c[0] - r, c[1] - r, c[2] - r],
            [c[0] + r, c[1] + r, c[2] + r],
        ),
    }
}

fn sphere(k: &mut Kernel) -> Result<(Brep, Expected)> {
    let (c, r) = SPHERE;
    Ok((api::create_sphere(k, pt(c), r)?, sphere_expected("sphere")))
}

fn sphere_no_pole(k: &mut Kernel) -> Result<(Brep, Expected)> {
    let (c, r) = SPHERE;
    Ok((
        api::create_sphere_no_pole(k, pt(c), r)?,
        sphere_expected("sphere (no pole)"),
    ))
}

/// The 0..90 degree wedge of the sphere: a quarter of the ball, its outer area a
/// quarter of the sphere's plus two half discs.
fn partial_sphere(k: &mut Kernel) -> Result<(Brep, Expected)> {
    let (c, r) = SPHERE;
    let b = api::create_partial_sphere(k, pt(c), r, 0.0, 90.0)?;
    Ok((
        b,
        Expected {
            name: "partial sphere 0..90",
            euler: 2,
            volume: PI * r * r * r / 3.0,
            area: PI * r * r + PI * r * r,
            bbox: ([c[0], c[1], c[2] - r], [c[0] + r, c[1] + r, c[2] + r]),
        },
    ))
}

fn torus(k: &mut Kernel) -> Result<(Brep, Expected)> {
    let (c, big, small) = TORUS;
    let o = big + small;
    Ok((
        api::create_torus(k, pt(c), big, small)?,
        Expected {
            name: "torus",
            euler: 0,
            volume: 2.0 * PI * PI * big * small * small,
            area: 4.0 * PI * PI * big * small,
            bbox: (
                [c[0] - o, c[1] - o, c[2] - small],
                [c[0] + o, c[1] + o, c[2] + small],
            ),
        },
    ))
}

/// The 0..90 degree sector of the torus: a quarter of the ring plus two end discs.
fn partial_torus(k: &mut Kernel) -> Result<(Brep, Expected)> {
    let (c, big, small) = TORUS;
    let o = big + small;
    Ok((
        api::create_partial_torus(k, pt(c), big, small, 0.0, 90.0)?,
        Expected {
            name: "partial torus 0..90",
            euler: 2,
            volume: PI * PI * big * small * small / 2.0,
            area: PI * PI * big * small + 2.0 * PI * small * small,
            bbox: (
                [c[0], c[1], c[2] - small],
                [c[0] + o, c[1] + o, c[2] + small],
            ),
        },
    ))
}

type Builder = fn(&mut Kernel) -> Result<(Brep, Expected)>;

const BUILDERS: [(&str, Builder); 5] = [
    ("sphere", sphere),
    ("sphere_no_pole", sphere_no_pole),
    ("partial_sphere", partial_sphere),
    ("torus", torus),
    ("partial_torus", partial_torus),
];

/// Builds, measures and checks one primitive, printing what it found.
fn check(label: &str, build: Builder) -> std::result::Result<(), String> {
    let mut k = Kernel::default();
    let (b, x) = build(&mut k).map_err(|e| format!("{label}: build failed: {e}"))?;
    let r = measure(&k, &b).map_err(|e| format!("{label}: measure failed: {e}"))?;
    println!(
        "{}: closed manifold solid: {}\n  V={} E={} F={} L={}  V-E+F={}  V-E+F-(L-F)={} (expect {})\n  volume {:.15} analytic {:.15} rel {:.2e}\n  area   {:.15} analytic {:.15} rel {:.2e}\n  bbox min {:?} max {:?}",
        x.name, r.manifold, r.v, r.e, r.f, r.l, r.vef(), r.euler(), x.euler,
        r.volume, x.volume, rel(r.volume, x.volume),
        r.area, x.area, rel(r.area, x.area),
        r.bbox.0, r.bbox.1
    );
    let f = failures(&r, &x);
    if f.is_empty() {
        println!("  OK");
        Ok(())
    } else {
        Err(format!("{}: {}", x.name, f.join("; ")))
    }
}

fn main() -> Result<()> {
    let mut failed = 0;
    for (label, build) in BUILDERS {
        if let Err(e) = check(label, build) {
            println!("{e}");
            failed += 1;
        }
    }
    println!("STL: not written (no tessellator is ported in solidcore yet)");
    if failed == 0 {
        Ok(())
    } else {
        Err(SolidError::invalid_input(
            "primitives",
            "some primitive checks failed",
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_solid(build: Builder) {
        let r = check("test", build);
        assert!(r.is_ok(), "{r:?}");
    }

    #[test]
    #[ignore = "blocked: SmSurface::DropCurve (PR #463), then SmTol::GetZoneTol3d of a ray-fire solution object (geom.tol)"]
    fn sphere_is_a_closed_solid() {
        assert_solid(sphere);
    }

    #[test]
    #[ignore = "blocked: SmSurface::DropCurve (PR #463), then SmEdge::GetTolerance in point classification during stitching (topology.edge, intersect/point_class.rs topology_of)"]
    fn sphere_no_pole_is_a_closed_solid() {
        assert_solid(sphere_no_pole);
    }

    #[test]
    #[ignore = "blocked: SmFace::CurveOnClassify in PrimitiveModel::axis_inside_face (query.face_p2)"]
    fn partial_sphere_is_a_closed_solid() {
        assert_solid(partial_sphere);
    }

    #[test]
    #[ignore = "blocked: SmSurface::DropCurve (PR #463), then SmTol::GetZoneTol3d of a ray-fire solution object (geom.tol)"]
    fn torus_is_a_closed_solid() {
        assert_solid(torus);
    }

    #[test]
    #[ignore = "blocked: SmFace::CurveOnClassify in PrimitiveModel::axis_inside_face (query.face_p2)"]
    fn partial_torus_is_a_closed_solid() {
        assert_solid(partial_torus);
    }
}
