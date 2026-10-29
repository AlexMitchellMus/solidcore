//! Convert a SolidCore B-rep into the neutral [`step_import::Model`].
//!
//! This is the export counterpart of [`crate::build`]. It maps exact SolidCore
//! topology and analytic/B-spline geometry to the neutral model that
//! [`step_export`] writes, so an engine B-rep can be exported without being
//! tessellated.
//!
//! The initial bridge writes one `MANIFOLD_SOLID_BREP` per B-rep, with all live faces
//! in one closed shell. Curves and surfaces are written as their exact analytic or
//! NURBS entities. Pcurves are not yet emitted; the 3D edge curves carry the trim
//! geometry.

use std::collections::BTreeMap;

use solidcore::curve::curve::p2::CurveKind;
use solidcore::topology::brep::p2::is_sheet;
use solidcore::topology::edgeuse::p1 as eu;
use solidcore::topology::{Brep, VertexId};
use step_import::model::{
    Body, BodyKind, Bound, Loop, OrientedEdge, Shell, Surface as ModelSurface, Vertex,
    Curve as ModelCurve, Edge, Face, Model, Frame,
};

/// Why a B-rep could not be mapped to the neutral model.
#[derive(Debug, Clone)]
pub struct ExportModelError {
    /// Human-readable reason.
    pub message: String,
}

impl std::fmt::Display for ExportModelError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.message)
    }
}

impl std::error::Error for ExportModelError {}

fn err(message: impl Into<String>) -> ExportModelError {
    ExportModelError { message: message.into() }
}

fn point(p: solidcore::math::Point3) -> [f64; 3] {
    [p.x, p.y, p.z]
}

fn vec(v: solidcore::math::Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}

fn unit(v: solidcore::math::Vec3) -> Result<[f64; 3], ExportModelError> {
    let raw = vec(v);
    let length = (raw[0] * raw[0] + raw[1] * raw[1] + raw[2] * raw[2]).sqrt();
    if !(length.is_finite() && length > 0.0) {
        return Err(err("degenerate direction"));
    }
    Ok([raw[0] / length, raw[1] / length, raw[2] / length])
}

fn frame(position: &solidcore::geom::axis2_placement::Axis2Placement) -> Result<Frame, ExportModelError> {
    let z = unit(position.z_axis())?;
    let x = unit(position.x_axis())?;
    let y = [
        z[1] * x[2] - z[2] * x[1],
        z[2] * x[0] - z[0] * x[2],
        z[0] * x[1] - z[1] * x[0],
    ];
    Ok(Frame { origin: point(position.origin()), x, y, z })
}

fn curve(curve: &solidcore::curve::curve::Curve) -> Result<ModelCurve, ExportModelError> {
    match curve.kind() {
        CurveKind::Line(line) => Ok(ModelCurve::Line {
            origin: point(line.line_point()),
            dir: unit(line.line_vector())?,
            magnitude: line.line_scale().abs(),
        }),
        CurveKind::Ellipse(ellipse) => {
            let frame = frame(ellipse.position())?;
            let x = ellipse.x_radius();
            let y = ellipse.y_radius();
            if (x - y).abs() <= 1.0e-12 * x.abs().max(y.abs()).max(1.0) {
                Ok(ModelCurve::Circle { frame, radius: x })
            } else {
                Ok(ModelCurve::Ellipse { frame, semi1: x, semi2: y })
            }
        }
        CurveKind::BSpline(curve) => {
            let (points, weights) = curve.control_polygon().map_err(|e| err(e.to_string()))?;
            let (knots, mults) = curve.knots(None).map_err(|e| err(e.to_string()))?;
            Ok(ModelCurve::BSpline(step_import::model::BSplineCurve {
                degree: curve.degree().map_err(|e| err(e.to_string()))? as u32,
                points: points.into_iter().map(point).collect(),
                weights: (!weights.is_empty()).then_some(weights),
                knots,
                mults: mults.into_iter().map(|m| m as u32).collect(),
            }))
        }
        other => Err(err(format!("unsupported curve kind {other:?}"))),
    }
}

fn surface(surface: &solidcore::surface::surface::Surface) -> Result<ModelSurface, ExportModelError> {
    use solidcore::surface::surface::p2::SurfaceKind;
    match surface.kind() {
        SurfaceKind::Plane(plane) => Ok(ModelSurface::Plane { frame: frame(plane.position())? }),
        SurfaceKind::Sphere(sphere) => Ok(ModelSurface::Sphere {
            frame: frame(sphere.base().position())?,
            radius: sphere.radius(),
        }),
        SurfaceKind::Torus(torus) => Ok(ModelSurface::Torus {
            frame: frame(torus.base().position())?,
            major: torus.major_radius(),
            minor: torus.minor_radius(),
        }),
        SurfaceKind::Cone(cone) => Ok(ModelSurface::Cone {
            frame: frame(cone.rev().position())?,
            radius: cone.base_radius(),
            semi_angle_deg: cone.semi_angle_deg(),
        }),
        SurfaceKind::SurfOfRevolution(revolution) => Ok(ModelSurface::Revolution {
            curve: Box::new(curve(revolution.gen_curve().ok_or_else(|| err("revolution has no generator"))?)?),
            axis_origin: point(revolution.position().origin()),
            axis_dir: unit(revolution.position().z_axis())?,
        }),
        SurfaceKind::SurfOfExtrusion(extrusion) => {
            let vector = extrusion.extrusion_vector();
            let magnitude = (vector.x * vector.x + vector.y * vector.y + vector.z * vector.z).sqrt();
            Ok(ModelSurface::Extrusion {
                curve: Box::new(curve(extrusion.gen_curve().ok_or_else(|| err("extrusion has no generator"))?)?),
                dir: unit(vector)?,
                magnitude,
            })
        }
        SurfaceKind::BSpline(bspline) => {
            let net = bspline.control_point_net().map_err(|e| err(e.to_string()))?;
            let v_count = net.v_count.max(1);
            let points: Vec<Vec<[f64; 3]>> = net
                .points
                .chunks(v_count)
                .map(|row| row.iter().map(|p| point(*p)).collect())
                .collect();
            let weights = if net.weights.is_empty() {
                None
            } else {
                Some(net.weights.chunks(v_count).map(|row| row.to_vec()).collect())
            };
            let canonical = bspline.canonical().map_err(|e| err(e.to_string()))?;
            Ok(ModelSurface::BSpline(step_import::model::BSplineSurface {
                u_degree: canonical.u_degree as u32,
                v_degree: canonical.v_degree as u32,
                points,
                weights,
                u_knots: canonical.u_knots,
                u_mults: canonical.u_multiplicities.into_iter().map(|m| m as u32).collect(),
                v_knots: canonical.v_knots,
                v_mults: canonical.v_multiplicities.into_iter().map(|m| m as u32).collect(),
            }))
        }
        other => Err(err(format!("unsupported surface kind {other:?}"))),
    }
}

/// Maps one B-rep to the neutral STEP model.
pub fn brep_to_model(brep: &Brep) -> Result<Model, ExportModelError> {
    let vertex_ids = brep.vertex_list().to_vec();
    let vertex_index: BTreeMap<VertexId, usize> = vertex_ids.iter().enumerate().map(|(i, v)| (*v, i)).collect();
    let vertices = vertex_ids.iter().map(|id| {
        brep.vertex(*id).map(|v| Vertex { step_id: id.index() as u64 + 1, point: point(v.point) }).map_err(|e| err(e.to_string()))
    }).collect::<Result<Vec<_>, _>>()?;

    let mut edges = Vec::new();
    let mut edge_index = BTreeMap::new();
    for id in brep.edge_list() {
        let edge = brep.edge(*id).map_err(|e| err(e.to_string()))?;
        let (start, end) = solidcore::topology::edge::start_end_vertices(brep, *id).map_err(|e| err(e.to_string()))?;
        let start = start.and_then(|v| vertex_index.get(&v).copied()).ok_or_else(|| err("edge has no start vertex"))?;
        let end = end.and_then(|v| vertex_index.get(&v).copied()).ok_or_else(|| err("edge has no end vertex"))?;
        let edge_curve = edge.curve.as_ref().ok_or_else(|| err("edge has no curve"))?;
        let interval = edge.interval;
        edge_index.insert(*id, edges.len());
        edges.push(Edge { step_id: id.index() as u64 + 1, start, end, curve: curve(edge_curve)?, same_sense: true, t_start: interval.min(), t_end: interval.max() });
    }

    let mut faces = Vec::new();
    let mut face_indices = Vec::new();
    let mut pcurves = std::collections::BTreeMap::new();
    let mut pcurve_kernel = solidcore::kernel::Kernel::default();
    for (face_pos, id) in brep.face_list().iter().enumerate() {
        let face = brep.face(*id).map_err(|e| err(e.to_string()))?;
        let surf = face.surface.as_ref().ok_or_else(|| err("face has no surface"))?;
        let faceuse = face.faceuse.ok_or_else(|| err("face has no faceuse"))?;
        let loops = eu::faceuse_loopuses(brep, faceuse, false, None).map_err(|e| err(e.to_string()))?;
        let mut bounds = Vec::new();
        for (loop_index, loopuse) in loops.iter().enumerate() {
            let edgeuses = eu::loopuse_edgeuses(brep, *loopuse).map_err(|e| err(e.to_string()))?;
            let mut oriented = Vec::new();
            let planar_face = matches!(
                face.surface.as_ref().map(|surface| surface.kind()),
                Some(solidcore::surface::surface::p2::SurfaceKind::Plane(_))
            );
            for (use_index, edgeuse) in edgeuses.into_iter().enumerate() {
                if !planar_face && let Ok((Some(uv_curve), _)) = eu::uv_trim_curve(brep, edgeuse) {
                    let interval = uv_curve.natural_interval();
                    let (low, high) = (interval.min(), interval.max());
                    let mut uv = Vec::with_capacity(65);
                    for step in 0..=64 {
                        let t = low + (high - low) * step as f64 / 64.0;
                        if let Ok(p) = uv_curve.evaluate_point(&pcurve_kernel, t) {
                            uv.push([p.x, p.y]);
                        }
                    }
                    if uv.len() >= 2 {
                        let opposite = eu::orientation(brep, edgeuse)
                            .map(|orientation| orientation != solidcore::nurbs::b_spline_surface_n_lib::p2::OrientType::Same)
                            .unwrap_or(false);
                        if opposite {
                            uv.reverse();
                        }
                        pcurves.insert((id.index() as u64 + 1, loop_index, use_index), uv);
                    }
                }
                let edge = eu::edge(brep, edgeuse).map_err(|e| err(e.to_string()))?.ok_or_else(|| err("edgeuse has no edge"))?;
                let forward = eu::orientation(brep, edgeuse).map_err(|e| err(e.to_string()))? == solidcore::nurbs::b_spline_surface_n_lib::p2::OrientType::Same;
                oriented.push(OrientedEdge { edge: *edge_index.get(&edge).ok_or_else(|| err("edgeuse edge not indexed"))?, forward });
            }
            let orientation = solidcore::topology::loopuse::orientation(brep, *loopuse).map_err(|e| err(e.to_string()))? == solidcore::nurbs::b_spline_surface_n_lib::p2::OrientType::Same;
            bounds.push(Loop { outer: loop_index == 0, orientation, bound: Bound::Edges(oriented) });
        }
        let same_sense = brep.faceuse(faceuse).map_err(|e| err(e.to_string()))?.orientation == solidcore::nurbs::b_spline_surface_n_lib::p2::OrientType::Same;
        face_indices.push(faces.len());
        faces.push(Face { step_id: id.index() as u64 + 1, surface: surface(surf)?, same_sense, loops: bounds });
    }
    if faces.is_empty() {
        return Err(err("B-rep has no faces"));
    }
    let mut kernel = solidcore::kernel::Kernel::default();
    let sheet = is_sheet(&mut kernel, brep).map_err(|e| err(e.to_string()))?;
    let body = Body {
        step_id: 1,
        name: None,
        kind: if sheet { BodyKind::ShellModel } else { BodyKind::ManifoldSolid },
        unit_mm: 1.0,
        uncertainty_mm: None,
        shells: vec![Shell { step_id: 1, closed: !sheet, orientation: true, faces: face_indices }],
        faces,
        edges,
        vertices,
    };
    Ok(Model { bodies: vec![body], pcurves, ..Model::default() })
}

#[cfg(test)]
mod tests {
    use super::*;
    use solidcore::kernel::Kernel;
    use solidcore::math::Point3;
    use solidcore::primitive::api_primitives as api;

    #[test]
    fn maps_a_box_and_a_sphere_to_the_neutral_model() {
        let mut kernel = Kernel::default();
        for brep in [
            api::create_box(&mut kernel, Point3::new(-1.0, -1.0, -1.0), 2.0, 2.0, 2.0).unwrap(),
            api::create_sphere(&mut kernel, Point3::new(0.0, 0.0, 0.0), 2.0).unwrap(),
        ] {
            let model = brep_to_model(&brep).expect("brep to model");
            assert_eq!(model.bodies.len(), 1);
            let body = &model.bodies[0];
            assert!(!body.faces.is_empty());
            assert!(!body.edges.is_empty());
            assert!(!body.vertices.is_empty());
            assert!(body.faces.iter().any(|face| matches!(face.surface, ModelSurface::Plane { .. } | ModelSurface::Sphere { .. })));
        }
    }
}
