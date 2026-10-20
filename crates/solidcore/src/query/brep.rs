//! B-rep bounding boxes and point classification, the polygon decimation attributes,
//! and the checked surface evaluation entry points of `SM_API`.
//!
//! Upstream: `source/SM_API/src/SmApiSurfaces.cpp`.
//!
//! | Rust                               | Upstream                          |
//! |------------------------------------|-----------------------------------|
//! | [`evaluate_surface_point`]         | `SmApiEvaluateSurfacePoint`       |
//! | [`evaluate_surface_normal`]        | `SmApiEvaluateSurfaceNormal`      |
//! | [`evaluate_surface_derivatives`]   | `SmApiEvaluateSurfaceDerivatives` |
//!
//! Deviations:
//! - **Out-params.** The evaluated point, normal or point and derivatives are the
//!   return value; on failure there is none (upstream evaluates into locals and leaves
//!   the caller's outputs untouched).
//! - **Pointers.** The `NULL` surface the entry points test for is `None`.
//! - **Statuses.** `SM_ERR_INVALID_INPUT` is [`SolidError::InvalidInput`]; errors of
//!   the surface evaluations pass through unchanged.

use crate::error::{Result, SolidError};
use crate::math::{Point2, Point3, Vec3};
use crate::nurbs::knots::Side;
use crate::surface::surface::Surface;

const OP_EVAL_POINT: &str = "query::brep::evaluate_surface_point";
const OP_EVAL_NORMAL: &str = "query::brep::evaluate_surface_normal";
const OP_EVAL_DERIVS: &str = "query::brep::evaluate_surface_derivatives";

// ---- SM_API surface evaluation ------------------------------------------------------

/// The input checks the three evaluation entry points share: a surface, and a
/// parameter inside its natural domain (bounds included).
fn checked_surface<'s>(
    surface: Option<&'s Surface>,
    uv: Point2,
    op: &'static str,
) -> Result<&'s Surface> {
    let Some(surface) = surface else {
        return Err(SolidError::invalid_input(op, "no surface"));
    };

    // Negated in-range test so a NaN parameter is rejected too.
    let domain = surface.natural_uv_domain();
    if !(uv.x >= domain.min().x
        && uv.x <= domain.max().x
        && uv.y >= domain.min().y
        && uv.y <= domain.max().y)
    {
        return Err(SolidError::invalid_input(
            op,
            "parameter outside the natural domain",
        ));
    }
    Ok(surface)
}

/// `SmApiEvaluateSurfacePoint`: the point of `surface` at `uv`.
///
/// # Errors
///
/// [`SolidError::InvalidInput`] without a surface or for a parameter outside (or NaN
/// against) the surface's natural domain; those of [`Surface::evaluate_point`].
pub fn evaluate_surface_point(surface: Option<&Surface>, uv: Point2) -> Result<Point3> {
    checked_surface(surface, uv, OP_EVAL_POINT)?.evaluate_point(uv)
}

/// `SmApiEvaluateSurfaceNormal`: the unit normal of `surface` at `uv`, from the upper
/// interval in both directions on a knot (`EvaluateNormal(uv, TRUE, TRUE)`).
///
/// # Errors
///
/// Those of [`evaluate_surface_point`]'s input checks; those of
/// [`Surface::evaluate_normal`].
pub fn evaluate_surface_normal(surface: Option<&Surface>, uv: Point2) -> Result<Vec3> {
    checked_surface(surface, uv, OP_EVAL_NORMAL)?.evaluate_normal(uv, Side::Left, Side::Left)
}

/// `SmApiEvaluateSurfaceDerivatives`: the point and first derivatives `(P, Du, Dv)` of
/// `surface` at `uv`, from the upper interval in both directions on a knot
/// (`Evaluate1stDerivatives(uv, TRUE, TRUE)`).
///
/// # Errors
///
/// Those of [`evaluate_surface_point`]'s input checks; those of
/// [`Surface::evaluate_1st_derivatives`].
pub fn evaluate_surface_derivatives(
    surface: Option<&Surface>,
    uv: Point2,
) -> Result<(Point3, Vec3, Vec3)> {
    checked_surface(surface, uv, OP_EVAL_DERIVS)?.evaluate_1st_derivatives(uv, true, true)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::extent2d::Extent2d;
    use crate::geom::vector2d::Vec2;
    use crate::surface::plane::p1::Plane;

    fn unit_plane() -> Surface {
        Surface::from(Plane::new(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec2::new(1.0, 1.0),
            &Extent2d::from_bounds(0.0, 0.0, 1.0, 1.0),
        ))
    }

    #[test]
    fn rejects_missing_surface_and_out_of_domain_parameters() {
        let s = unit_plane();
        let bad = |r: Result<Point3>| matches!(r, Err(SolidError::InvalidInput { .. }));
        assert!(bad(evaluate_surface_point(None, Point2::new(0.5, 0.5))));
        assert!(bad(evaluate_surface_point(Some(&s), Point2::new(1.5, 0.5))));
        assert!(bad(evaluate_surface_point(
            Some(&s),
            Point2::new(f64::NAN, 0.5)
        )));
        // the domain bounds are inside
        assert!(evaluate_surface_point(Some(&s), Point2::new(1.0, 0.0)).is_ok());
    }

    #[test]
    fn evaluates_a_plane() {
        let s = unit_plane();
        let uv = Point2::new(0.25, 0.75);
        let p = evaluate_surface_point(Some(&s), uv).unwrap();
        assert_eq!((p.x, p.y, p.z), (0.25, 0.75, 0.0));
        let n = evaluate_surface_normal(Some(&s), uv).unwrap();
        assert_eq!((n.x, n.y, n.z), (0.0, 0.0, 1.0));
        let (q, du, dv) = evaluate_surface_derivatives(Some(&s), uv).unwrap();
        assert_eq!((q.x, q.y, q.z), (p.x, p.y, p.z));
        assert_eq!(
            (du.x, du.y, du.z, dv.x, dv.y, dv.z),
            (1.0, 0.0, 0.0, 0.0, 1.0, 0.0)
        );
    }
}
