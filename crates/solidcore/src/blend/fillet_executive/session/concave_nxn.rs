//! Live concave `NxN` corner topology and mate vertex geometry over the session's
//! objects.
//!
//! Upstream: `source/SMLib/src/SmFilletCorner.cpp` (`SmFilletConcaveNxNCorner::MakeCornerTopology`,
//! 5196–5329), `source/SMLib/src/SmFilletGeom.cpp` (`SmFilletVertex::CalcMateGeom`).
//!
//! Deviations: each outward edge use's solver is resolved before the corner builds
//! anything (upstream looks them up inside its loop, without the mate fallback the
//! convex corner uses; the lookups have no effects). A setback corner goes to the
//! convex `NxN` session path, as upstream calls `SmFilletConvexNxNCorner::MakeCornerTopology`.
//! `CalcMateGeom` follows [`FilletBrep::calc_mate_geom`] with the main vertex
//! computed by the session's vertex dispatch; its recursion through nested mates is
//! bounded by [`MAX_MATE_DEPTH`].

use std::collections::BTreeMap;

use super::{
    invalid, live_corner_vert_geom, FilletBody, FilletObject, SessionObjects, SolverObject,
};
use crate::blend::fillet_corner::p3::FilletCorner;
use crate::blend::fillet_corner::topology::NxNTopologyEnv;
use crate::blend::fillet_executive::p3::FilletExecutive;
use crate::blend::fillet_geom::p1::{evaluate_mate, FilStatus, FilletBrep, FilletCornerId};
use crate::blend::fillet_geom::p3::MAX_MATE_DEPTH;
use crate::error::{EntityRef, Result, SolidError};
use crate::intersect::surface_intersector::p1::SrfIndex;
use crate::kernel::Kernel;
use crate::surface::surface::Surface;
use crate::topology::{Brep, FaceId, VertexId};

const OP_MATE: &str = "blend::fillet_executive::session::calc_mate_geom_live";

/// `SmFilletVertex::CalcMateGeom` over the live session: computes the main vertex
/// (`GetMate(0)`) first through the session's vertex dispatch (status ignored, as
/// upstream; [`SolidError::Unported`] and the nesting bound propagate), then places
/// the mate on the opposite rail of its first fillet geometry at the main vertex's
/// solver point on the mate's face.
///
/// # Errors
/// As [`FilletBrep::calc_mate_geom`].
#[allow(clippy::too_many_arguments)] // the session's borrowed object stores
pub(super) fn calc_mate_geom_live(
    k: &Kernel,
    target: &mut Brep,
    pseudo: &mut FilletBrep,
    solvers: &mut [Option<SolverObject>],
    corners: &[FilletCorner],
    extended_surfaces: &BTreeMap<FaceId, Surface>,
    v: VertexId,
    depth: usize,
) -> Result<()> {
    if depth >= MAX_MATE_DEPTH {
        return Err(SolidError::invalid_input(
            OP_MATE,
            "mate vertices nested deeper than MAX_MATE_DEPTH",
        ));
    }

    // our mate (which is the 'real' fillet vertex) must be done first.
    let main_fv = pseudo
        .fillet_vertex(v)?
        .mate(0)
        .ok_or_else(|| SolidError::null_pointer(OP_MATE))?;
    if !pseudo.fillet_vertex(main_fv)?.is_processed() {
        match live_corner_vert_geom(
            k,
            target,
            pseudo,
            solvers,
            corners,
            extended_surfaces,
            main_fv,
            depth + 1,
        ) {
            Err(e) if e.is_unported() || is_depth_error(&e) => return Err(e),
            _ => {} // status ignored
        }
    }
    if !pseudo.fillet_vertex(main_fv)?.is_processed() {
        return Ok(());
    }

    // Calculating Main's geom can cause us to be calculated also.
    if pseudo.fillet_vertex(v)?.is_processed() {
        return Ok(());
    }

    // Just work from the first FilletGeom.
    let geoms = pseudo.vertex_fillet_geoms(v)?;

    // check state - at least 1 railEdge must be connected to this vertex
    let fillet_geom = *geoms.first().ok_or_else(|| SolidError::sm_err(OP_MATE))?;

    // We have to figure out which rail of the FilletGeom we're on.
    let other_side_eu = match pseudo.fillet_vertex(main_fv)?.point_class_object() {
        None => return Ok(()),
        Some(EntityRef::Edgeuse(eu)) => eu,
        Some(_) => {
            return Err(SolidError::invalid_input(
                OP_MATE,
                "main vertex classification object is not an edge use",
            ))
        }
    };

    let solver_id = pseudo
        .geom(fillet_geom)?
        .fillet_solver()
        .ok_or_else(|| SolidError::null_pointer(OP_MATE))?;
    let fil_solver = solvers
        .get(solver_id.0 as usize)
        .and_then(Option::as_ref)
        .map(SolverObject::base)
        .ok_or_else(|| invalid("mate solver is not live"))?;
    let rail_index = fil_solver.find_index_of_rail_x_side_edgeuse(target, other_side_eu)?;
    let rail_index = 1 - rail_index; // We're on the opposite rail.

    let main_vtx_use = pseudo
        .vu_at_rail_end(main_fv, Some(fillet_geom), None, None)?
        .ok_or_else(|| {
            SolidError::invalid_input(OP_MATE, "main vertex without a rail-end vertex use")
        })?;
    let srf = if rail_index == 0 {
        SrfIndex::First
    } else {
        SrfIndex::Second
    };
    let uv = pseudo
        .fillet_vertexuse(main_vtx_use)?
        .tsect_pnt()
        .uv_pos(srf);

    if pseudo
        .vu_at_rail_end(v, Some(fillet_geom), None, None)?
        .is_none()
    {
        return Ok(());
    }

    let face = match pseudo.fillet_vertex(v)?.point_class_object() {
        Some(EntityRef::Face(f)) => f,
        None => return Err(SolidError::null_pointer(OP_MATE)),
        Some(_) => {
            return Err(SolidError::invalid_input(
                OP_MATE,
                "classification object is not a face",
            ))
        }
    };

    // compute sPnt = pFace->Surface(sUV), use ExtendedSurface as needed
    let pnt = evaluate_mate(k, target, face, uv)?;

    // set mate outputs
    pseudo.fillet_vertex_mut(v)?.set_original_uv(uv);
    pseudo.set_vertex_point(v, pnt)?;
    pseudo
        .fillet_vertex_mut(v)?
        .set_status(FilStatus::Processed);
    Ok(())
}

fn is_depth_error(e: &SolidError) -> bool {
    matches!(
        e,
        SolidError::InvalidInput { operation, reason }
            if *operation == OP_MATE && reason.contains("MAX_MATE_DEPTH")
    )
}

impl SessionObjects<'_> {
    pub(super) fn make_concave_nxn_corner_topology(
        &mut self,
        k: &Kernel,
        id: FilletCornerId,
        executive: &mut FilletExecutive<FilletObject>,
    ) -> Result<()> {
        const NATIVE: &str = "SmFilletConcaveNxNCorner::MakeCornerTopology";
        let corner = self
            .corners
            .get(id.0 as usize)
            .ok_or_else(|| invalid("corner is not live"))?;
        // Determine if this is a setback corner
        if corner.set_back_dist > corner.approx_tol_3d * 100.0 {
            // Let SmFilletConvexNxNCorner handle setbacks
            return self.make_convex_nxn_corner_topology(k, id, executive);
        }
        let edgeuses = corner
            .all_eus
            .as_deref()
            .ok_or_else(|| invalid("corner has no outward edge uses"))?;
        let mut solver_ids = Vec::new();
        for &eu in edgeuses {
            let object = FilletObject::Topology(FilletBody::Target, EntityRef::Edgeuse(eu));
            match executive.fillet_solver_of_edgeuse(object) {
                Some(FilletObject::Solver(solver)) => solver_ids.push(solver),
                None => return Err(SolidError::null_pointer(NATIVE)),
                Some(_) => return Err(invalid("edge use does not name an actual solver")),
            }
        }
        let SessionObjects {
            target,
            pseudo,
            solvers,
            corners,
            ..
        } = self;
        let solvers = solver_ids
            .into_iter()
            .map(|solver| {
                solvers
                    .get(solver.0 as usize)
                    .and_then(Option::as_ref)
                    .map(SolverObject::base)
                    .ok_or_else(|| invalid("corner solver is not live"))
            })
            .collect::<Result<Vec<_>>>()?;
        let corner = corners
            .get_mut(id.0 as usize)
            .ok_or_else(|| invalid("corner is not live"))?;
        corner.make_concave_nxn_corner_topology(NxNTopologyEnv {
            k,
            target,
            pseudo,
            corner: id,
            solvers,
        })
    }
}
