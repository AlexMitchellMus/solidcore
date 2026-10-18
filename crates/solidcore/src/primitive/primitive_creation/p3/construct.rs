//! Curve-sweep, swung, skinned and blended primitive orchestration.
//!
//! Upstream: `source/SMLib/src/SmPrimitiveCreation.cpp`.
//!
//! | Rust | Upstream |
//! |------|----------|
//! | [`PrimitiveCreation::create_curve_sweep`] | `CreateCurveSweep` |
//! | [`PrimitiveCreation::create_swung_primitive`] | `CreateSwungPrimitive` |
//! | [`PrimitiveCreation::create_skin_from_faces`] | `CreateSkinFromFaces` |
//! | [`PrimitiveCreation::create_blend_primitive`] | `CreateBlendPrimitive` |
//!
//! Deviations: typed service resolution replaces pointers. Output face lists
//! preserve partial results; owned geometry handles cleanup. No recursion.

use super::{
    geometry::{prepare_curve_sweep, GeometryHandle},
    PrimitiveAccess,
};
use crate::curve::composite_curve::CompositeCurve;
use crate::curve::curve::p2::CurvePropertyType;
use crate::curve::curve::p7::PropertyAnalysisOptions;
use crate::curve::curve::Curve;
use crate::error::{ignore_status, Result, SolidError};
use crate::math::Vec3;
use crate::nurbs::b_spline_curve::p1::BSplineCurve;
use crate::nurbs::b_spline_surface::p1::{BSplineSurface, SurfParamType};
use crate::nurbs::b_spline_surface_n_lib::p2::SkinnedSurfaceOptions;
use crate::nurbs::knots::Side;
use crate::primitive::primitive_creation::{p1::PrimitiveCreation, SweepOptions};
use crate::surface::surface::{
    p2::{BlendEndType, CurvatureType, DerivSurfDefinition},
    Surface,
};
use crate::sweep::topology_sweep::access::{FaceRef, RegionRef};
use crate::topology::{BrepId, EdgeuseId, FaceId};

const OP: &str = "primitive::primitive_creation::p3";

pub(super) fn target(primitive: &PrimitiveCreation, brep: BrepId) -> Result<RegionRef> {
    Ok((
        brep,
        primitive
            .region()
            .ok_or_else(|| SolidError::null_pointer(OP))?,
    ))
}
fn face_surface(a: &impl PrimitiveAccess, face: FaceRef) -> Result<&Surface> {
    crate::topology::face::p2::surface(a.brep(face.0)?, face.1)?
        .ok_or_else(|| SolidError::null_pointer(OP))
}

/// Parameters for `CreateBlendPrimitive`.
#[derive(Debug, Clone, Copy)]
pub struct BlendPrimitiveOptions {
    /// Requested blend continuity.
    pub curvature: CurvatureType,
    /// Smooth or cusped ends.
    pub ends: BlendEndType,
    /// 0 guesses directions; 1 matches starts; 2 matches start to end.
    pub curves_direction: u64,
    /// Starting takeoff length.
    pub start_tangent_length: f64,
    /// Ending takeoff length (divided by the start, including zero, upstream).
    pub end_tangent_length: f64,
}

impl PrimitiveCreation {
    /// Build swung surfaces for every B-spline-compatible path/profile pair.
    /// Only pair (0,0) initializes the scale, even if earlier slots were skipped.
    pub fn create_swung_primitive(
        &mut self,
        a: &mut impl PrimitiveAccess,
        brep: BrepId,
        paths: &[Option<Curve>],
        profiles: &[Option<Curve>],
        scale: f64,
        approx_tol: f64,
    ) -> Result<()> {
        self.init_entity_maps();
        let _ = target(self, brep)?;
        let mut curve_scale = 0.0;
        let mut faces = Vec::new();
        for (i, path) in paths.iter().enumerate() {
            let Some(path) = path.as_ref().and_then(super::geometry::bspline) else {
                continue;
            };
            for (j, profile) in profiles.iter().enumerate() {
                let Some(profile) = profile.as_ref().and_then(super::geometry::bspline) else {
                    continue;
                };
                if i == 0 && j == 0 {
                    let mut point = profile.evaluate_point(profile.natural_interval().min())?;
                    let normal = Vec3::Y;
                    let d = -point.dot(normal);
                    let curve = Curve::from(path.clone());
                    let property = curve.global_property_analysis(
                        a.kernel(brep)?,
                        &profile.natural_interval(),
                        CurvePropertyType::PlaneIntersection,
                        approx_tol,
                        &PropertyAnalysisOptions {
                            property_value: Some(d),
                            opt_vectors: &[normal],
                            handle: Some(GeometryHandle(curve.natural_interval())),
                            ..Default::default()
                        },
                    )?;
                    if let Some(solution) = property.solutions.get(0) {
                        point = path.evaluate_point(solution.start.param(0)?)?;
                    }
                    curve_scale = (1.0 / point.x).abs();
                }
                let surface =
                    BSplineSurface::create_swung_surface(&profile, &path, scale * curve_scale)?;
                let surface = Surface::from(surface);
                let domain = surface.natural_uv_domain();
                a.surface_faces(brep, surface, domain, &mut faces)?;
            }
        }
        Ok(())
    }

    /// Loft matching loops of faces, or the supplied sequence of profile curves.
    /// `new_faces` is reset before any failing call, as upstream.
    pub fn create_skin_from_faces(
        &mut self,
        a: &mut impl PrimitiveAccess,
        brep: BrepId,
        faces: &[FaceRef],
        profiles: &[BSplineCurve],
        approx_tol: f64,
        degree: usize,
        new_faces: &mut Vec<FaceId>,
    ) -> Result<()> {
        self.init_entity_maps();
        new_faces.clear();
        let _ = target(self, brep)?;
        let mut skin_curves = Vec::new();
        let mut layers = Vec::new();
        let mut layer_count = 0;
        for (i, &face) in faces.iter().enumerate() {
            let domain = crate::topology::face::p2::uv_domain(a.brep(face.0)?, face.1)?;
            face_surface(a, face)?.evaluate_normal(
                domain.evaluate(0.5, 0.5),
                Side::Left,
                Side::Left,
            )?;
            let composites = a.face_composites(face)?;
            if i == 0 {
                layer_count = composites.len();
            } else if layer_count != composites.len() {
                return Err(SolidError::sm_err(OP));
            }
            for (j, mut composite) in composites.into_iter().enumerate() {
                skin_curves.push(composite.make_composite_nurb()?);
                layers.push(j);
            }
        }
        if faces.is_empty() {
            skin_curves.extend_from_slice(profiles);
            layer_count = 1;
            layers.resize(skin_curves.len(), 0);
        }
        let mut made = Vec::new();
        for layer in 0..layer_count {
            let curves: Vec<_> = skin_curves
                .iter()
                .zip(&layers)
                .filter(|(_, l)| **l == layer)
                .map(|(c, _)| c.clone())
                .collect();
            let surface = BSplineSurface::create_skinned_surface(
                &curves,
                SkinnedSurfaceOptions {
                    approx_tol,
                    degree,
                    deriv_surfs: Some([None, None]),
                    ..Default::default()
                },
            )?;
            let surface = Surface::from(surface);
            let domain = surface.natural_uv_domain();
            ignore_status(a.surface_faces(brep, surface, domain, &mut made))?;
            new_faces.extend_from_slice(&made);
        }
        Ok(())
    }

    /// Construct a blend and add its face, preserving the ignored face status.
    pub fn create_blend_primitive(
        &mut self,
        a: &mut impl PrimitiveAccess,
        brep: BrepId,
        first: Option<EdgeuseId>,
        second: Option<EdgeuseId>,
        output: &mut Option<FaceId>,
        options: BlendPrimitiveOptions,
    ) -> Result<()> {
        self.init_entity_maps();
        let _ = target(self, brep)?;
        let flag = match options.curvature {
            CurvatureType::NoCurvature => 1,
            CurvatureType::G3FromSurface => 3,
            _ => 2,
        };
        let mut derivatives = DerivSurfDefinition {
            curvature_type: options.curvature,
            blend_end_type: options.ends,
            start_tan_length: options.start_tangent_length,
            end_tan_scale: options.end_tangent_length / options.start_tangent_length,
            ..Default::default()
        };
        let surface = crate::topology::brep::p1::create_edge_edge_blend(
            a.kernel(brep)?,
            a.brep(brep)?,
            flag,
            first,
            second,
            options.curves_direction,
            Some(&mut derivatives),
        )?;
        a.set_editing(brep, true)?;
        let domain = surface.natural_uv_domain();
        ignore_status(a.surface_face(brep, surface, domain, output))?;
        Ok(())
    }

    /// Sweep copies of the profiles along the path, optionally scaling and capping.
    /// Only `side_faces` is reset; cap output slots retain their upstream lifetime.
    pub fn create_curve_sweep(
        &mut self,
        a: &mut impl PrimitiveAccess,
        brep: BrepId,
        profiles: &[Option<Curve>],
        path: &BSplineCurve,
        scale_reference: Option<&BSplineCurve>,
        scale_curve: Option<&BSplineCurve>,
        approx_tol: f64,
        options: SweepOptions,
        test_continuity: bool,
        start_faces: &mut Vec<FaceId>,
        side_faces: &mut Vec<FaceId>,
        end_faces: &mut Vec<FaceId>,
    ) -> Result<()> {
        self.init_entity_maps();
        side_faces.clear();
        if profiles.is_empty() {
            return Err(SolidError::sm_err(OP));
        }
        let region = target(self, brep)?;
        let previous = a.kernel(brep)?.context.doing_boolean();
        a.kernel_mut(brep)?.context.set_doing_boolean(true);
        let result = (|| {
            let prepared = prepare_curve_sweep(
                a.kernel(brep)?,
                profiles,
                path,
                scale_reference,
                scale_curve,
                options,
                test_continuity,
            )?;
            let count = prepared.profiles.len();
            let curves = prepared.profiles.into_iter().map(Curve::from).collect();
            let composites = CompositeCurve::build_composites_from_curves(
                curves,
                a.brep(brep)?.tolerance(),
                0.0,
            )?;
            for mut composite in composites.composites {
                let profile = composite.make_composite_nurb()?;
                let surface = BSplineSurface::create_sweep_surface(
                    &profile,
                    &prepared.path,
                    prepared.scale.as_ref(),
                    options.translational_sweep,
                    approx_tol,
                )?;
                let surface = Surface::from(surface);
                let domain = surface.natural_uv_domain();
                let mut made = Vec::new();
                ignore_status(a.surface_faces(brep, surface, domain, &mut made))?;
                side_faces.extend(made);
            }
            ignore_status(a.transform_primitive(brep, &prepared.rotation_back))?;
            let placement = if options.move_to_profile {
                &prepared.profile_back
            } else {
                &prepared.path_back
            };
            ignore_status(a.transform_primitive(brep, placement))?;
            if options.cap_ends {
                let mut starts = Vec::new();
                let mut ends = Vec::new();
                let mut direction = SurfParamType::U;
                let total = side_faces.len();
                for (i, &face) in side_faces.iter().enumerate() {
                    let surface = face_surface(a, (brep, face))?;
                    let domain = surface.natural_uv_domain();
                    if i < count {
                        direction = if options.translational_sweep {
                            SurfParamType::V
                        } else {
                            SurfParamType::U
                        };
                        let param = if options.translational_sweep {
                            domain.min().y
                        } else {
                            domain.min().x
                        };
                        let curve = surface
                            .create_iso_parametric_curve(
                                a.kernel(brep)?,
                                direction,
                                param,
                                approx_tol / 10.0,
                                None,
                                false,
                            )?
                            .curve;
                        starts.push(Some(curve));
                    }
                    if count > total || i >= total - count {
                        let param = if options.translational_sweep {
                            direction = SurfParamType::V;
                            domain.max().y
                        } else {
                            domain.max().x
                        };
                        let curve = surface
                            .create_iso_parametric_curve(
                                a.kernel(brep)?,
                                direction,
                                param,
                                approx_tol / 10.0,
                                None,
                                false,
                            )?
                            .curve;
                        ends.push(Some(curve));
                    }
                }
                a.planar_faces(region, &mut starts, approx_tol, start_faces)?;
                a.planar_faces(region, &mut ends, approx_tol, end_faces)?;
            }
            ignore_status(a.stitch_and_orient(brep))?;
            Ok(())
        })();
        a.kernel_mut(brep)?.context.set_doing_boolean(previous);
        result
    }
}
