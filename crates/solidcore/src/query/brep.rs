//! B-rep bounding boxes and point classification, the polygon decimation attributes,
//! and the checked surface evaluation entry points of `SM_API`.
//!
//! Upstream: `source/SM_API/src/SmApiSurfaces.cpp`, `source/SMLib/inc/SmPolyDecimate.h`,
//! `source/SMLib/src/SmPolyDecimate.cpp`.
//!
//! | Rust                                         | Upstream                                       |
//! |----------------------------------------------|------------------------------------------------|
//! | [`evaluate_surface_point`]                   | `SmApiEvaluateSurfacePoint`                    |
//! | [`evaluate_surface_normal`]                  | `SmApiEvaluateSurfaceNormal`                   |
//! | [`evaluate_surface_derivatives`]             | `SmApiEvaluateSurfaceDerivatives`              |
//! | [`PolyVertexClass`]                          | `SmPolyVertexClass` (`SmPoly.h`)               |
//! | [`decimate_vertex_attribute`], [`DecimateVertexData`] | `SmDecimateVertexAttr::SmDecimateVertexAttr` (`~` is `Drop`) |
//! | [`DecimateVertexData::make_copy`]            | `SmDecimateVertexAttr::MakeCopy`               |
//! | [`quadric_vertex_attribute`], [`QuadricVertexData`] | `SmQuadricVertexAttr::SmQuadricVertexAttr` (`~` is `Drop`) |
//! | [`Attribute::make_copy`] (data cloned)       | `SmQuadricVertexAttr::MakeCopy`                |
//! | [`quadric_edge_attribute`], [`QuadricEdgeData`] | `SmQuadricEdgeAttr::SmQuadricEdgeAttr` (`~` is `Drop`) |
//! | [`QuadricEdgeData::make_copy`]               | `SmQuadricEdgeAttr::MakeCopy`                  |
//!
//! Deviations:
//! - **Out-params.** The evaluated point, normal or point and derivatives are the
//!   return value; on failure there is none (upstream evaluates into locals and leaves
//!   the caller's outputs untouched).
//! - **Pointers.** The `NULL` surface the entry points test for is `None`.
//! - **Statuses.** `SM_ERR_INVALID_INPUT` is [`SolidError::InvalidInput`]; errors of
//!   the surface evaluations pass through unchanged.
//! - **Decimation attributes.** They are [`Attribute`]s whose [`AttributeValue`] is
//!   the class's data, as for the other attribute classes; `MakeCopy` is
//!   [`Attribute::make_copy`] (same id and behaviour, no users, upstream's `new
//!   (crContext)` has no counterpart). The destructors only overwrite fields of an
//!   object being freed (`~SmDecimateVertexAttr`); `Drop` frees it. `SmQuadricEdgeAttr`'s
//!   constructor leaves `m_lSortIndex` uninitialised: it is `None` until set
//!   ([`QuadricEdgeData::sort_index`]). `SmPolyVertexClass` (`SmPoly.h`) is defined
//!   here, the only place that uses it so far.

use crate::error::{Result, SolidError};
use crate::geom::attribute::{Attribute, AttributeValue, Behavior};
use crate::math::{Point2, Point3, Vec3, BIG_DOUBLE};
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

// ---- polygon decimation attributes ---------------------------------------------------

/// `SmPolyVertexClass`: the class of a polygon vertex during decimation; the
/// discriminant is the upstream value.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum PolyVertexClass {
    /// `SM_PVC_NOT_COMPUTED`: not yet computed.
    NotComputed = 0,
    /// `SM_PVC_UNKNOWN`: connected to a spine edge, or non-manifold.
    Unknown = 1,
    /// `SM_PVC_SIMPLE`: an interior vertex on no lamina, feature or silhouette edge.
    Simple = 2,
    /// `SM_PVC_BOUNDARY`: on the boundary of the mesh (has lamina edges).
    Boundary = 3,
    /// `SM_PVC_INTERIOR`: an interior point between two feature edges.
    Interior = 4,
    /// `SM_PVC_CORNER`: an interior point between three or more feature edges.
    Corner = 5,
}

impl PolyVertexClass {
    /// Every class, in discriminant order.
    pub const ALL: [Self; 6] = [
        Self::NotComputed,
        Self::Unknown,
        Self::Simple,
        Self::Boundary,
        Self::Interior,
        Self::Corner,
    ];
}

/// The data of an `SmDecimateVertexAttr` (attribute id `SM_AI_VERTEX_DECIMATION_VALUES`).
#[derive(Debug, Clone, PartialEq)]
pub struct DecimateVertexData {
    /// `m_lSortIndex`.
    pub sort_index: u64,
    /// `m_eVertClass`.
    pub vertex_class: PolyVertexClass,
    /// `m_dVertexDistance`: the weighted chord height from the vertex to the baseline
    /// between its feature or boundary edges.
    pub vertex_distance: f64,
    /// `m_dMaxFaceError`: the largest error of the faces around the vertex.
    pub max_face_error: f64,
    /// `m_vAveragePlaneNormal`: the area-weighted normal of the faces around the vertex.
    pub average_plane_normal: Vec3,
}

impl DecimateVertexData {
    /// `SmDecimateVertexAttr::MakeCopy`'s data: everything but the sort index, which
    /// keeps the constructor's `0`.
    #[must_use]
    pub fn make_copy(&self) -> Self {
        Self {
            sort_index: 0,
            ..self.clone()
        }
    }
}

/// The data of an `SmQuadricVertexAttr` (attribute id `SM_AI_QUADRIC_VERTEX_VALUES`): the
/// quadric error form `NewVert . Q . NewVert + 2 b . NewVert + c` of the planes of the
/// original polygons around a vertex.
#[derive(Debug, Clone, PartialEq)]
pub struct QuadricVertexData {
    /// `m_Q`: the rows of the quadratic error matrix.
    pub q: [Vec3; 3],
    /// `m_b`.
    pub b: Vec3,
    /// `m_c`.
    pub c: f64,
}

/// The data of an `SmQuadricEdgeAttr` (attribute id `SM_AI_QUADRIC_EDGE_VALUES`): the
/// vertex an edge contracts to and its error.
#[derive(Debug, Clone, PartialEq)]
pub struct QuadricEdgeData {
    /// `m_vNewPoint`.
    pub new_point: Point3,
    /// `m_dError`.
    pub error: f64,
    /// `m_lSortIndex`; `None` while uninitialised (as the constructor leaves it).
    pub sort_index: Option<u64>,
}

impl QuadricEdgeData {
    /// `SmQuadricEdgeAttr::MakeCopy`'s data: the point and error, sort index `0`.
    #[must_use]
    pub fn make_copy(&self) -> Self {
        Self {
            new_point: self.new_point,
            error: self.error,
            sort_index: Some(0),
        }
    }
}

/// `SmDecimateVertexAttr(lAttributeID, eBehavior)` (upstream default
/// [`Behavior::Copy`]): sort index `0`, class not computed, vertex distance
/// `-SM_BIG_DOUBLE`, no face error, zero average normal.
#[must_use]
pub fn decimate_vertex_attribute(attribute_id: u64, behavior: Behavior) -> Attribute {
    Attribute::with_value(
        attribute_id,
        behavior,
        AttributeValue::DecimateVertex(DecimateVertexData {
            sort_index: 0,
            vertex_class: PolyVertexClass::NotComputed,
            vertex_distance: -BIG_DOUBLE,
            max_face_error: 0.0,
            average_plane_normal: Vec3::new(0.0, 0.0, 0.0),
        }),
    )
}

/// `SmQuadricVertexAttr(lAttributeID, eBehavior)` (upstream default
/// [`Behavior::Copy`]): the zero error form.
#[must_use]
pub fn quadric_vertex_attribute(attribute_id: u64, behavior: Behavior) -> Attribute {
    let zero = Vec3::new(0.0, 0.0, 0.0);
    Attribute::with_value(
        attribute_id,
        behavior,
        AttributeValue::QuadricVertex(QuadricVertexData {
            q: [zero, zero, zero],
            b: zero,
            c: 0.0,
        }),
    )
}

/// `SmQuadricEdgeAttr(lAttributeID, eBehavior)` (upstream default [`Behavior::Copy`]):
/// the origin as new point, no error, no sort index.
#[must_use]
pub fn quadric_edge_attribute(attribute_id: u64, behavior: Behavior) -> Attribute {
    Attribute::with_value(
        attribute_id,
        behavior,
        AttributeValue::QuadricEdge(QuadricEdgeData {
            new_point: Vec3::new(0.0, 0.0, 0.0),
            error: 0.0,
            sort_index: None,
        }),
    )
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
    fn decimation_attribute_copies_reset_sort_indices() {
        let mut a = decimate_vertex_attribute(7, Behavior::Reference);
        if let AttributeValue::DecimateVertex(d) = a.value_mut() {
            d.sort_index = 9;
            d.vertex_class = PolyVertexClass::Corner;
            d.max_face_error = 0.5;
        }
        let c = a.make_copy();
        assert_eq!((c.attribute_id(), c.behavior()), (7, Behavior::Reference));
        let AttributeValue::DecimateVertex(d) = c.value() else {
            panic!("class lost");
        };
        assert_eq!(
            (d.sort_index, d.vertex_class, d.max_face_error),
            (0, PolyVertexClass::Corner, 0.5)
        );
        assert_eq!(d.vertex_distance, -BIG_DOUBLE);

        let e = quadric_edge_attribute(8, Behavior::Copy);
        let AttributeValue::QuadricEdge(d) = e.make_copy().value().clone() else {
            panic!("class lost");
        };
        assert_eq!(d.sort_index, Some(0));
        assert_eq!(e.type_code(), crate::geom::attribute::type_code::ATTRIBUTE);
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
