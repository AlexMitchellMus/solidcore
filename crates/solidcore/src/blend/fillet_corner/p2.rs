//! Fillet corners, part 2: the concave edge uses of a corner and the topology adjustments
//! of the `Nx2` corner (and, in later batches, the variable-radius and analytic corner
//! surfaces and the `Nx1` closed and `NxN` adjustments).
//!
//! Upstream: `source/SMLib/src/SmFilletCorner.cpp`, `source/SMLib/inc/SmFilletCorner.h`.
//!
//! | Rust                                              | Upstream                                         |
//! |---------------------------------------------------|--------------------------------------------------|
//! | [`find_concave_edgeuses`]                         | `sm_FindConcaveEdgeuses`                         |
//! | [`FilletCorner::adjust_concave_corner_topology`]  | `SmFilletNx2Corner::AdjustConcaveCornerTopology` |
//! | [`FilletCorner::update_split_face`]               | `SmFilletCorner::UpdateSplitFace`                |
//! | [`FilletCorner::tangent_fillet_solver_nx1_closed`] | `SmFilletNx1ClosedCorner::GetTangentFilletSolver` |
//! | [`FilletCorner::tangent_fillet_solver_nx2`]       | `SmFilletNx2Corner::GetTangentFilletSolver`      |
//!
//! Deviations:
//! - **Out-param.** The concave edge uses (`rConcaveEUs`, appended to) are returned; the
//!   `iDebugLevel` argument (graphics only) is dropped.
//! - **Undefined input.** Upstream reads `(*pAllEUs)[GetSize() - 1]` unchecked, past the
//!   array when it is empty; an empty `all_eus` is [`SolidError::InvalidInput`]. Its
//!   callers pass a corner's `m_pAllEUs`, which `SmFilletCorner::Create` never leaves
//!   empty.
//! - **Corner methods** are methods of `blend.fillet_corner_p3`'s [`FilletCorner`] taking
//!   its [`CornerEnv`] (that module's docs). Where upstream dereferences `NULL` or casts a
//!   classification object blindly, the port returns [`SolidError::InvalidInput`];
//!   `NER` is `SM_ERR_NULL_POINTER` and `SER` propagates. Evaluations whose status
//!   upstream ignores leave [`Point3::UNDEF`]. The new fillet vertices are created with
//!   their point (upstream `SetPoint`s them later; no vertex use is attached in between).
//! - **`UpdateSplitFace`** reads `m_vVertices` / `m_vEdges` directly (the inline
//!   `GetFilletVertices` / `GetFilletEdges` of unit `blend.fillet_corner_p1` append
//!   them). Its callees always return `SM_SUCCESS` upstream (ignored); their ports'
//!   errors (lookups, and the `NULL` dereferences of `SmFilletEdge::ShouldSwapFaces` and
//!   face classification) propagate. The faces are ids (upstream compares and
//!   dereferences the pointers); a fillet edge's corner is resolved in the caller's
//!   `corners` (upstream's `GetFilletCorner()` pointer).
//! - **`GetTangentFilletSolver`** overrides take the executive's edge-use-to-solver map
//!   ([`FilletExecutiveData`]) and the B-rep being filleted (for `sm_TestG1Continuity`)
//!   explicitly; `SM_ASSERT_ERR` (a debug assertion) is not reproduced, the `NULL` it
//!   returns is `None`.
//! - **Unported callees.** `SmFilletEdge::MakeFilletEdge` (`blend.fillet_geom_p5`)
//!   returns [`SolidError::Unported`]; the code after it is ported (including
//!   `AddExtendedOriginalSurface`, unit `blend.fillet_corner_p1`, on the executive's map
//!   in [`CornerEnv`]), so wiring it is the only change needed.
//!
//! Every loop runs over a given list once; there is no recursion.

use crate::blend::fillet_corner::p1::{
    check_face_complimentary_corner_concavity, test_g1_continuity,
};
use crate::blend::fillet_corner::p3::{
    move_vertexuse, primary_eu, set_vertex_point, unchecked_point, CornerEnv, FilletCorner,
    FilletExecutiveData,
};
use crate::blend::fillet_geom::p1::{
    FilStatus, FilletBrep, FilletEdgeType, FilletGeomId, FilletSolverId, FilletVertex,
    FilletVertexType,
};
use crate::curve::curve::Curve;
use crate::error::{EntityRef, Result, SolidError};
use crate::geom::geom_utility::p2::ContinuityType;
use crate::intersect::point_class::PointClassType;
use crate::intersect::surface_intersector::p1::{SrfIndex, TsectPnt};
use crate::kernel::Kernel;
use crate::math::{rad_to_deg, Point3, Vec3, BIG_DOUBLE, EFF_ZERO};
use crate::surface::surface::Surface;
use crate::topology::edgeuse::p1 as eu1;
use crate::topology::{edge, Brep, EdgeId, EdgeuseId, FaceId, VertexId, VertexuseId};

const OP_CONCAVE: &str = "blend::fillet_corner::p2::find_concave_edgeuses";
const OP_CONVEX_TOPO: &str = "blend::fillet_corner::p2::adjust_topology_convex_nxn";

/// `SM_BIG_DOUBLE`.
const SM_BIG_DOUBLE: f64 = BIG_DOUBLE;
const OP_CONCAVE_TOPO: &str = "blend::fillet_corner::p2::adjust_concave_corner_topology";
const OP_TANGENT_SOLVER: &str = "blend::fillet_corner::p2::tangent_fillet_solver";

/// `sm_FindConcaveEdgeuses`: the first pair of consecutive edge uses of `all_eus`
/// (cyclically, starting with the last and the first) that bound a concave face corner
/// at `vertex` (`sm_CheckFaceComplimentaryCornerConcavity` without a reference edge use),
/// as `[previous, current]`; empty when no pair does.
///
/// # Errors
///
/// [`SolidError::InvalidInput`] for an empty `all_eus` (upstream indexes past the
/// array); those of [`check_face_complimentary_corner_concavity`].
pub fn find_concave_edgeuses(
    k: &Kernel,
    b: &mut Brep,
    vertex: VertexId,
    all_eus: &[EdgeuseId],
) -> Result<Vec<EdgeuseId>> {
    let Some(&last) = all_eus.last() else {
        return Err(SolidError::invalid_input(
            OP_CONCAVE,
            "no edge uses (upstream indexes past the array)",
        ));
    };
    let mut concave = Vec::new();
    let mut prev = last;
    for &curr in all_eus {
        if check_face_complimentary_corner_concavity(k, b, vertex, prev, curr, None)? {
            concave.push(prev);
            concave.push(curr);
            return Ok(concave);
        }
        prev = curr;
    }
    Ok(concave)
}

impl FilletCorner {
    /// `SmFilletNx1ClosedCorner::GetTangentFilletSolver(pFS)`: the corner's one solver
    /// when it is `fs` and its edge meets itself with G1 continuity at the corner vertex
    /// (`sm_TestG1Continuity` within `m_dTangencyTolRadians`); `None` otherwise
    /// (`SM_ASSERT_ERR` when the corner does not have exactly one filleted edge use or
    /// `fs` is not its solver; not an error when it is not G1).
    ///
    /// # Errors
    ///
    /// [`SolidError::InvalidInput`] for an unset `m_pSolverEUs` (upstream dereferences
    /// it); those of [`test_g1_continuity`].
    pub fn tangent_fillet_solver_nx1_closed(
        &self,
        k: &Kernel,
        b: &Brep,
        exec: &FilletExecutiveData,
        fs: FilletSolverId,
    ) -> Result<Option<FilletSolverId>> {
        let solver_eus = self
            .solver_eus
            .as_deref()
            .ok_or_else(|| SolidError::invalid_input(OP_TANGENT_SOLVER, "m_pSolverEUs is NULL"))?;
        let [edge_use] = solver_eus else {
            return Ok(None); // SM_ASSERT_ERR
        };
        let fil_solver = exec.fillet_solver_of_edgeuse(*edge_use);

        if fil_solver != Some(fs) {
            return Ok(None); // SM_ASSERT_ERR
        }

        // We have to test its continuity.
        let edge = eu1::edge(b, *edge_use)?;
        let tangency_tol_deg = rad_to_deg(self.tangency_tol_radians);
        if test_g1_continuity(k, b, self.vertex, edge, edge, tangency_tol_deg)?.0 {
            return Ok(fil_solver);
        }

        Ok(None) // not an error.
    }

    /// `SmFilletNx2Corner::GetTangentFilletSolver(pFS)`: when the corner's two filleted
    /// edges meet with G1 continuity at the corner vertex, the solver of the other one
    /// than `fs`; `None` otherwise (`SM_ASSERT_ERR` when the corner does not have
    /// exactly two filleted edge uses or `fs` is neither solver).
    ///
    /// # Errors
    ///
    /// [`SolidError::InvalidInput`] for an unset `m_pSolverEUs` (upstream dereferences
    /// it); those of [`test_g1_continuity`].
    pub fn tangent_fillet_solver_nx2(
        &self,
        k: &Kernel,
        b: &Brep,
        exec: &FilletExecutiveData,
        fs: FilletSolverId,
    ) -> Result<Option<FilletSolverId>> {
        let solver_eus = self
            .solver_eus
            .as_deref()
            .ok_or_else(|| SolidError::invalid_input(OP_TANGENT_SOLVER, "m_pSolverEUs is NULL"))?;
        let [eu0, eu1_] = solver_eus else {
            return Ok(None); // SM_ASSERT_ERR
        };

        // We have to test its continuity.
        let edge0 = eu1::edge(b, *eu0)?;
        let edge1 = eu1::edge(b, *eu1_)?;

        let tangency_tol_deg = rad_to_deg(self.tangency_tol_radians);
        if !test_g1_continuity(k, b, self.vertex, edge0, edge1, tangency_tol_deg)?.0 {
            return Ok(None); // not an error.
        }

        let fs0 = exec.fillet_solver_of_edgeuse(*eu0);
        let fs1 = exec.fillet_solver_of_edgeuse(*eu1_);

        if fs0 == Some(fs) {
            return Ok(fs1);
        }
        if fs1 == Some(fs) {
            return Ok(fs0);
        }

        Ok(None) // SM_ASSERT_ERR
    }

    /// `SmFilletCorner::UpdateSplitFace`: the target face `orig_face` was split, making
    /// `split_face`; every fillet vertex of the corner (then every fillet edge) moves to
    /// the face it now lies in. `corners` resolves the fillet corner of a fillet edge
    /// (upstream's `SmFilletEdge::GetFilletCorner()` pointer), as
    /// [`FilletBrep::should_swap_faces`] does.
    ///
    /// # Errors
    ///
    /// Those of [`FilletBrep::update_split_face`] and
    /// [`FilletBrep::update_edge_split_face`].
    pub fn update_split_face(
        &self,
        env: &mut CornerEnv<'_>,
        orig_face: FaceId,
        split_face: FaceId,
        corners: &[FilletCorner],
    ) -> Result<()> {
        for &v in &self.vertices {
            env.pseudo
                .update_split_face(env.k, v, env.target, orig_face, split_face)?;
        }
        for &e in &self.edges {
            env.pseudo.update_edge_split_face(
                env.k,
                e,
                Some(orig_face),
                Some(split_face),
                env.target,
                corners,
            )?;
        }
        Ok(())
    }

    /// `SmFilletNx2Corner::AdjustConcaveCornerTopology`: for a mixed-convexity corner
    /// whose two fillets do not meet on the side edge (their rails hit it more than
    /// `m_dThisApproxTol3d` apart), splits the `SM_FV_RAIL_X_EDGEUSE` vertex: two new
    /// such vertices where the rails hit the side edge, the old one becomes
    /// `SM_FV_FILLET_X_FILLET`, and three new fillet edges fill the gap on the extended
    /// side face.
    ///
    /// # Errors
    ///
    /// `SM_ERR_NULL_POINTER` when no `SM_FV_RAIL_X_EDGEUSE` vertex has two fillet
    /// geometries (or a cross-section edge use and its mate have none);
    /// [`SolidError::InvalidInput`] where upstream dereferences `NULL` or casts a
    /// classification object that is not an edge use; [`SolidError::Unported`] (module
    /// docs).
    pub fn adjust_concave_corner_topology(&mut self, env: &mut CornerEnv<'_>) -> Result<()> {
        // Find the fillet vertex to be adjusted, pAdjustVert.
        // It's the only SM_FV_RAIL_X_EDGEUSE.
        let mut found = None;
        for &v in &self.vertices {
            if env.pseudo.fillet_vertex(v)?.fillet_vertex_type() != FilletVertexType::RailXEdgeuse {
                continue;
            }
            if let [g1, g2] = env.pseudo.vertex_fillet_geoms(v)?[..] {
                found = Some((v, g1, g2));
                break;
            }
        }
        let (adjust_vert, fillet_geom1, fillet_geom2) =
            found.ok_or_else(|| SolidError::null_pointer(OP_CONCAVE_TOPO))?;

        let Some(EntityRef::Edgeuse(side_eu1)) =
            env.pseudo.fillet_vertex(adjust_vert)?.point_class_object()
        else {
            return Err(null_deref("classification object is not an edge use"));
        };
        let vu1 = env
            .pseudo
            .vu_at_rail_end(adjust_vert, Some(fillet_geom1), None, None)?
            .ok_or_else(|| null_deref("no rail vertex use"))?;
        let vu2 = env
            .pseudo
            .vu_at_rail_end(adjust_vert, Some(fillet_geom2), Some(vu1), None)?
            .ok_or_else(|| null_deref("no rail vertex use"))?;
        let edge_param1 = env
            .pseudo
            .fillet_vertexuse(vu1)?
            .tsect_pnt()
            .curve_parameter;
        let edge_param2 = env
            .pseudo
            .fillet_vertexuse(vu2)?
            .tsect_pnt()
            .curve_parameter;
        let side_edge = eu1::edge(env.target, side_eu1)?
            .ok_or_else(|| null_deref("side edge use without edge"))?;
        let side_edge_curve = env
            .target
            .edge(side_edge)?
            .curve
            .as_ref()
            .ok_or_else(|| null_deref("side edge without curve"))?;
        let pnt1 = unchecked_point(side_edge_curve.evaluate_point(env.k, edge_param1))?;
        let pnt2 = unchecked_point(side_edge_curve.evaluate_point(env.k, edge_param2))?;

        // See if we have coincident geometry, i.e. two rails intersect at this vertex.
        let gap_dist = pnt1.distance(pnt2);
        if gap_dist < self.approx_tol_3d {
            // No adjustments are needed
            return Ok(());
        }

        // See which of the two rail-side edge intersections is closer to the filleted
        // vertex, to decide Large/Small (read before the new vertices are made: nothing
        // here depends on them).
        let vert_pnt = env
            .target
            .vertex(
                self.vertex
                    .ok_or_else(|| null_deref("m_cpVertex is NULL"))?,
            )?
            .point;
        let near1 = pnt1.distance(vert_pnt) < pnt2.distance(vert_pnt);
        let (pnt_l, t_l, pnt_s, t_s) = if near1 {
            (pnt2, edge_param2, pnt1, edge_param1)
        } else {
            (pnt1, edge_param1, pnt2, edge_param2)
        };
        let (large_fg, small_fg, vu_on_large_fg, vu_on_small_fg) = if near1 {
            (fillet_geom2, fillet_geom1, vu2, vu1)
        } else {
            (fillet_geom1, fillet_geom2, vu1, vu2)
        };
        let side_face = if near1 {
            eu1::face(env.target, side_eu1)?
        } else {
            let radial = eu1::radial(env.target, side_eu1)?
                .ok_or_else(|| null_deref("side edge use without radial partner"))?;
            eu1::face(env.target, radial)?
        };

        // Create two new FilletVertices where each rail hits the side edge.
        let new_vert_l = new_rail_x_edgeuse_vertex(env, pnt_l, t_l)?; // "Large"
        let new_vert_s = new_rail_x_edgeuse_vertex(env, pnt_s, t_s)?; // "Small"
        self.vertices.push(new_vert_l);
        self.vertices.push(new_vert_s);

        let fv = env.pseudo.fillet_vertex_mut(adjust_vert)?;
        fv.set_fillet_vertex_type(FilletVertexType::FilletXFillet);
        fv.set_status(FilStatus::Unprocessed);

        move_vertexuse(env.pseudo, vu_on_large_fg, adjust_vert, new_vert_l)?;
        move_vertexuse(env.pseudo, vu_on_small_fg, adjust_vert, new_vert_s)?;

        let new_edge = make_fillet_edge(env, adjust_vert, new_vert_l)?;
        let fe = env.pseudo.fillet_edge_mut(new_edge)?;
        fe.set_fillet_edge_type(FilletEdgeType::FilletXSideFace);
        fe.set_original_face(side_face);
        let prim_eu = primary_eu(env.pseudo, new_edge, OP_CONCAVE_TOPO)?;
        env.pseudo
            .set_edgeuse_fillet_geom(prim_eu, Some(large_fg))?;

        // Need to extend the side surface
        let side_face_id = side_face.ok_or_else(|| null_deref("no side face"))?;
        let side_surface = env
            .target
            .face(side_face_id)?
            .surface
            .as_ref()
            .ok_or_else(|| null_deref("side face without surface"))?;
        let new_surface =
            side_surface.create_extended_surface(env.k, 10.0 * gap_dist, ContinuityType::G1)?;

        // Copy pNewSurface, when possible as an analytic surface.
        let new_analytic_surface = new_surface.copy_and_add_analytics(env.k)?;
        self.add_extended_original_surface(env.exec, side_face_id, new_analytic_surface);

        let new_edge1 = make_fillet_edge(env, new_vert_l, new_vert_s)?;
        let fe1 = env.pseudo.fillet_edge_mut(new_edge1)?;
        fe1.set_fillet_edge_type(FilletEdgeType::OnEdge);
        fe1.set_original_edge(Some(side_edge));
        fe1.set_original_face(side_face);

        let new_edge2 = make_fillet_edge(env, new_vert_s, adjust_vert)?;
        let fe2 = env.pseudo.fillet_edge_mut(new_edge2)?;
        fe2.set_fillet_edge_type(FilletEdgeType::RailExtension);
        fe2.set_original_face(side_face);
        let prim_eu2 = primary_eu(env.pseudo, new_edge2, OP_CONCAVE_TOPO)?;
        env.pseudo
            .set_edgeuse_fillet_geom(prim_eu2, Some(small_fg))?;

        self.edges.push(new_edge);
        self.edges.push(new_edge1);
        self.edges.push(new_edge2);

        // pAdjustVert comes in with two mates, in the common face.
        // Find them -- they are at the other end of the two FilletEdges
        // of type SM_FE_CROSS_SECTION -- and reset them as mates of
        // the two FilletVertices we just created.
        let vus = env.pseudo.brep().vertex(adjust_vert)?.vertexuses.clone();
        for vu in vus {
            let eu = eu1::vertexuse_edgeuse(env.pseudo.brep(), vu)?
                .ok_or_else(|| null_deref("vertex use without edge use"))?;
            let fillet_e = env
                .pseudo
                .brep()
                .edgeuse(eu)?
                .edge
                .ok_or_else(|| null_deref("edge use without edge"))?;
            if env.pseudo.fillet_edge(fillet_e)?.fillet_edge_type() != FilletEdgeType::CrossSection
            {
                continue;
            }

            // FilletEdge of type SM_FE_CROSS_SECTION:
            let other_vert = edge::other_vertex(env.pseudo.brep(), fillet_e, Some(adjust_vert))?;
            let mut fg = env.pseudo.fillet_edgeuse(eu)?.fillet_geom();
            if fg.is_none() {
                let mate = eu1::mate(env.pseudo.brep(), eu)?
                    .ok_or_else(|| null_deref("edge use without mate"))?;
                fg = env.pseudo.fillet_edgeuse(mate)?.fillet_geom();
                if fg.is_none() {
                    return Err(SolidError::null_pointer(OP_CONCAVE_TOPO));
                }
            }
            // pOtherVert->GetMate(0) is read only when the geometry matches
            let has_mate = |pseudo: &FilletBrep| -> Result<bool> {
                let other =
                    other_vert.ok_or_else(|| null_deref("cross section without other vertex"))?;
                Ok(pseudo.fillet_vertex(other)?.mate(0).is_some())
            };
            let new_v = if fg == Some(small_fg) && has_mate(env.pseudo)? {
                Some(new_vert_s)
            } else if fg == Some(large_fg) && has_mate(env.pseudo)? {
                Some(new_vert_l)
            } else {
                None
            };
            if let (Some(new_v), Some(other)) = (new_v, other_vert) {
                env.pseudo
                    .fillet_vertex_mut(new_v)?
                    .set_mate(0, Some(other));
                env.pseudo
                    .fillet_vertex_mut(other)?
                    .set_mate(0, Some(new_v));
            }
        }
        let fv = env.pseudo.fillet_vertex_mut(adjust_vert)?;
        fv.set_mate(0, None);
        fv.set_mate(1, None);

        // Allow fillets to be extended at the corner
        self.blending = false;

        self.topology_adjusted = true;

        Ok(())
    }
}

impl FilletCorner {
    /// `SmFilletConvexNxNCorner::AdjustTopology` (also `SmFilletConcaveNxNCorner`'s): a
    /// three-sided corner patch that is too tight is stretched into a four-sided one.
    /// The `SM_FV_RAIL_X_RAIL` vertex whose two cross-section iso edges (estimated on
    /// its original face) are farthest apart, without crossing and of comparable length,
    /// is replaced by a rail-rail-interpolating edge between itself and a new
    /// `SM_FV_MATE` vertex, when that distance is at least a fifth of the shortest
    /// corner edge.
    ///
    /// # Errors
    ///
    /// `SM_ERR` for a corner vertex that is not `SM_FV_RAIL_X_RAIL`;
    /// `SM_ERR_NULL_POINTER` for a vertex without original face, a face without surface,
    /// a cross section or corner edge without fillet geometry, a filleted edge without
    /// curve or a new vertex without rail vertex use; those of the evaluations;
    /// [`SolidError::InvalidInput`] where upstream dereferences `NULL` or indexes past
    /// an array (module docs).
    #[allow(clippy::too_many_lines)]
    pub fn adjust_topology_convex_nxn(&mut self, env: &mut CornerEnv<'_>) -> Result<()> {
        const OP: &str = OP_CONVEX_TOPO;
        // Only non-chamfer & three-sided corners can be adjusted for now
        let total = self.vertices.len();
        if total > 3
        /* || m_bChamfer */
        {
            return Ok(());
        }

        let mut found_vert_index = None;
        let mut min_dist = SM_BIG_DOUBLE;
        let mut max_dist = 0.0;

        // Look for a vertex of the corner patch that can be 'stretched',
        // i.e. replacing the vertex by a rail-rail-interpolating edge
        for ii in 0..total {
            let curr_vert = vertex_at(&self.vertices, ii, OP)?;
            if env.pseudo.fillet_vertex(curr_vert)?.fillet_vertex_type()
                != FilletVertexType::RailXRail
            {
                return Err(SolidError::sm_err(OP)); // Unknown case
            }
            let next_vert = vertex_at(&self.vertices, (ii + 1) % total, OP)?;

            let vert_pnt = env.pseudo.vertex_point(curr_vert)?;
            let dist_to_next_vert = vert_pnt.distance(env.pseudo.vertex_point(next_vert)?);
            if min_dist > dist_to_next_vert {
                min_dist = dist_to_next_vert;
            }

            // Get the original face which pCurrVert is on
            let orig_face = point_class_face(env.pseudo, curr_vert)
                .ok_or_else(|| SolidError::null_pointer(OP))?;

            let mut is_candidate = true;
            let vertexuses = env.pseudo.brep().vertex(curr_vert)?.vertexuses.clone();
            let mut pnt = [Point3::UNDEF; 4];
            let mut vec = [Vec3::UNDEF; 2];
            let mut count = 0usize;
            // Try to determine the estimated distance between two adjacent fillets
            // if iso-edges are used in the stretched corner patch.
            for vu in vertexuses {
                let eu = eu1::vertexuse_edgeuse(env.pseudo.brep(), vu)?
                    .ok_or_else(|| invalid(OP, "vertex use without edge use"))?;
                let e = env
                    .pseudo
                    .brep()
                    .edgeuse(eu)?
                    .edge
                    .ok_or_else(|| invalid(OP, "edge use without edge"))?;
                if env.pseudo.fillet_edge(e)?.fillet_edge_type() != FilletEdgeType::CrossSection {
                    continue;
                }
                let fillet_geom = edgeuse_or_mate_geom(env.pseudo, eu, OP)?
                    .ok_or_else(|| SolidError::null_pointer(OP))?;
                let fil_solver = solver_of_geom(env, fillet_geom, OP)?;

                let other_v = edge::other_vertex(env.pseudo.brep(), e, Some(curr_vert))?
                    .ok_or_else(|| invalid(OP, "cross section without other vertex"))?;
                let rail_vu = env
                    .pseudo
                    .vu_at_rail_end(other_v, Some(fillet_geom), None, None)?
                    .ok_or_else(|| invalid(OP, "no rail vertex use"))?;
                let rail_index = env
                    .solvers
                    .get(fil_solver)
                    .ok_or_else(|| invalid(OP, "fillet solver id not in the solver table"))?
                    .find_index_of_rail_on_face(env.target, Some(orig_face))?;
                let uv = env
                    .pseudo
                    .fillet_vertexuse(rail_vu)?
                    .tsect_pnt()
                    .uv_pos(srf_index(rail_index));
                // (more than two cross sections: upstream writes past sVec)
                let (Some(p), true) = (pnt.get_mut(count), count < vec.len()) else {
                    return Err(invalid(
                        OP,
                        "more than two cross sections at a corner vertex",
                    ));
                };
                *p = face_surface(env.target, orig_face, OP)?.evaluate_point(env.k, uv)?;

                // The following is a simple way to determine if two iso-edges
                // are crossing each other.
                let v = *p - vert_pnt;
                if let Some(slot) = vec.get_mut(count) {
                    *slot = v;
                }
                let filleted_eu = env
                    .solvers
                    .get(fil_solver)
                    .ok_or_else(|| invalid(OP, "fillet solver id not in the solver table"))?
                    .edgeuse(0)?
                    .ok_or_else(|| invalid(OP, "fillet solver without edge use"))?;
                let filleted_edge = eu1::edge(env.target, filleted_eu)?
                    .ok_or_else(|| invalid(OP, "filleted edge use without edge"))?;
                let ivl = edge::interval(env.target, filleted_edge)?;
                let curve = env
                    .target
                    .edge(filleted_edge)?
                    .curve
                    .as_ref()
                    .ok_or_else(|| SolidError::null_pointer(OP))?;
                let tangent = if edge::start_vertex(env.target, filleted_edge)? == self.vertex {
                    tangent_at(env.k, curve, ivl.min(), OP)?
                } else {
                    -tangent_at(env.k, curve, ivl.max(), OP)?
                };
                if tangent.dot(v) < -EFF_ZERO {
                    // Reject this vertex if crossing of the iso-edges happened
                    is_candidate = false;
                    break;
                }
                count += 1;
            }

            if !is_candidate {
                continue;
            }

            // SM_ASSERT(lCount == 2); (fewer leave sVec at its SM_UNDEF_DOUBLE members)
            let [vec0, vec1] = vec;
            let vec_len1 = vec0.length();
            let vec_len2 = vec1.length();
            if vec_len1 < 0.2 * vec_len2 || vec_len2 < 0.2 * vec_len1 {
                continue;
            }

            let [pnt0, pnt1, ..] = pnt;
            let dist = pnt0.distance(pnt1);
            if dist > max_dist {
                max_dist = dist;
                found_vert_index = Some(ii);
            }
        }

        // Determine if we have found a good candidate. Let's use
        // 1/5th of the shortest corner edge as the threshold
        let Some(found_vert_index) = found_vert_index else {
            // Do not adjust this corner
            return Ok(());
        };
        if max_dist < min_dist / 5.0 {
            // Do not adjust this corner
            return Ok(());
        }

        // Insert an interpolating edge to substitute the found vertex
        let adjust_vert = vertex_at(&self.vertices, found_vert_index, OP)?;
        // (a vertex without face: upstream calls GetSurface on NULL)
        let orig_face = point_class_face(env.pseudo, adjust_vert)
            .ok_or_else(|| invalid(OP, "corner vertex without original face"))?;
        let orig_surf = env
            .target
            .face(orig_face)?
            .surface
            .clone()
            .ok_or_else(|| SolidError::null_pointer(OP))?;

        // Determine where to inset the new edge
        let mut found_edge_index = None;
        let mut next_e = None;
        let mut prev_e = *total
            .checked_sub(1)
            .and_then(|i| self.edges.get(i))
            .ok_or_else(|| invalid(OP, "fewer corner edges than vertices"))?;
        for ii in 0..total {
            let e = *self
                .edges
                .get(ii)
                .ok_or_else(|| invalid(OP, "fewer corner edges than vertices"))?;
            next_e = Some(e);
            if edge::vertex(env.pseudo.brep(), e)? == Some(adjust_vert) {
                found_edge_index = Some(ii);
                break;
            }
            prev_e = e;
        }
        // SM_ASSERT(lFoundEdgeIndex != 9999); (upstream inserts at 9999 otherwise)
        let found_edge_index =
            found_edge_index.ok_or_else(|| invalid(OP, "no corner edge starts at the vertex"))?;
        let next_e = next_e.ok_or_else(|| invalid(OP, "no corner edge"))?;

        // Set pAdjustVert as the SM_FV_MATE type and compute its geometry
        let prev_eu = primary_eu(env.pseudo, prev_e, OP)?;
        let prev_fg = env
            .pseudo
            .fillet_edgeuse(prev_eu)?
            .fillet_geom()
            .ok_or_else(|| SolidError::null_pointer(OP))?;
        let prev_v = vertex_at(&self.vertices, (found_vert_index + total - 1) % total, OP)?;
        let prev_vu = env
            .pseudo
            .vu_at_rail_end(prev_v, Some(prev_fg), None, None)?
            .ok_or_else(|| invalid(OP, "no rail vertex use"))?;
        let prev_fs = solver_of_geom(env, prev_fg, OP)?;
        let rail_index1 = env
            .solvers
            .get(prev_fs)
            .ok_or_else(|| invalid(OP, "fillet solver id not in the solver table"))?
            .find_index_of_rail_on_face(env.target, Some(orig_face))?;
        let tsect_pnt1 = env.pseudo.fillet_vertexuse(prev_vu)?.tsect_pnt().clone();
        let uv1 = tsect_pnt1.uv_pos(srf_index(rail_index1));
        let pnt1 = orig_surf.evaluate_point(env.k, uv1)?;
        env.pseudo
            .fillet_vertex_mut(adjust_vert)?
            .set_original_uv(uv1);
        set_vertex_point(env.pseudo, adjust_vert, pnt1)?;
        let fv = env.pseudo.fillet_vertex_mut(adjust_vert)?;
        fv.set_fillet_vertex_type(FilletVertexType::Mate);
        fv.set_mate(0, Some(prev_v));
        env.pseudo
            .fillet_vertex_mut(prev_v)?
            .set_mate(0, Some(adjust_vert));
        let vu1 = env
            .pseudo
            .vu_at_rail_end(adjust_vert, Some(prev_fg), None, None)?
            .ok_or_else(|| SolidError::null_pointer(OP))?;
        copy_uv_pos(env.pseudo, vu1, &tsect_pnt1)?;

        // Make a new vertex
        // Compute its geometry
        let next_eu = primary_eu(env.pseudo, next_e, OP)?;
        let next_fg = env
            .pseudo
            .fillet_edgeuse(next_eu)?
            .fillet_geom()
            .ok_or_else(|| SolidError::null_pointer(OP))?;
        let next_v = vertex_at(&self.vertices, (found_vert_index + 1) % total, OP)?;
        let next_vu = env
            .pseudo
            .vu_at_rail_end(next_v, Some(next_fg), None, None)?
            .ok_or_else(|| invalid(OP, "no rail vertex use"))?;
        let next_fs = solver_of_geom(env, next_fg, OP)?;
        let rail_index2 = env
            .solvers
            .get(next_fs)
            .ok_or_else(|| invalid(OP, "fillet solver id not in the solver table"))?
            .find_index_of_rail_on_face(env.target, Some(orig_face))?;
        let tsect_pnt2 = env.pseudo.fillet_vertexuse(next_vu)?.tsect_pnt().clone();
        let uv2 = tsect_pnt2.uv_pos(srf_index(rail_index2));
        let pnt2 = orig_surf.evaluate_point(env.k, uv2)?;
        let new_vert = env.pseudo.brep_mut().add_vertex(pnt2)?;
        let mut fv = FilletVertex::new(env.k, env.corner);
        fv.set_original_uv(uv2);
        fv.set_point_class(PointClassType::Face, Some(EntityRef::Face(orig_face)));
        fv.set_fillet_vertex_type(FilletVertexType::Mate);
        fv.set_mate(0, Some(next_v));
        fv.set_status(FilStatus::Processed);
        env.pseudo.insert_fillet_vertex(new_vert, fv);
        env.pseudo
            .fillet_vertex_mut(next_v)?
            .set_mate(0, Some(new_vert));
        let vu2 = env
            .pseudo
            .vu_at_rail_end(adjust_vert, Some(next_fg), None, None)?
            .ok_or_else(|| SolidError::null_pointer(OP))?;
        copy_uv_pos(env.pseudo, vu2, &tsect_pnt2)?;

        // Move pVU2 from pAdjustVert to pNewVert
        move_vertexuse(env.pseudo, vu2, adjust_vert, new_vert)?;

        let mut vu3 = env
            .pseudo
            .brep()
            .edgeuse(next_eu)?
            .vertexuse
            .ok_or_else(|| invalid(OP, "edge use without vertex use"))?;
        if env.pseudo.brep().vertexuse(vu3)?.vertex != Some(adjust_vert) {
            let mate = eu1::mate(env.pseudo.brep(), next_eu)?
                .ok_or_else(|| invalid(OP, "edge use without mate"))?;
            vu3 = env
                .pseudo
                .brep()
                .edgeuse(mate)?
                .vertexuse
                .ok_or_else(|| invalid(OP, "edge use without vertex use"))?;
        }
        // Move pVU3 from pAdjustVert to pNewVert
        move_vertexuse(env.pseudo, vu3, adjust_vert, new_vert)?;

        self.vertices.insert(found_vert_index, new_vert);

        let new_edge = make_fillet_edge(env, adjust_vert, new_vert)?;
        let fe = env.pseudo.fillet_edge_mut(new_edge)?;
        fe.set_fillet_edge_type(FilletEdgeType::RailRailInterpolation);
        fe.set_original_face(Some(orig_face));

        self.edges.insert(found_edge_index, new_edge);

        self.topology_adjusted = true;

        Ok(())
    }
}

fn invalid(op: &'static str, reason: &'static str) -> SolidError {
    SolidError::invalid_input(op, reason)
}

/// `m_vVertices[i]`.
fn vertex_at(vertices: &[VertexId], i: usize, op: &'static str) -> Result<VertexId> {
    vertices
        .get(i)
        .copied()
        .ok_or_else(|| SolidError::invalid_input(op, "vertex index"))
}

/// `SM_CAST_PTR(SmFace, pV->GetPointClassObject())`: the face a fillet vertex is
/// classified on, `None` for another kind or none.
fn point_class_face(pseudo: &FilletBrep, v: VertexId) -> Option<FaceId> {
    match pseudo.fillet_vertex(v).ok()?.point_class_object() {
        Some(EntityRef::Face(f)) => Some(f),
        _ => None,
    }
}

/// The surface of a face of the B-rep being filleted (`GetSurface()`, dereferenced).
fn face_surface<'b>(b: &'b Brep, f: FaceId, op: &'static str) -> Result<&'b Surface> {
    b.face(f)?
        .surface
        .as_ref()
        .ok_or_else(|| SolidError::invalid_input(op, "face without surface"))
}

/// `pEU->GetFilletGeom()`, or its mate's when it has none.
fn edgeuse_or_mate_geom(
    pseudo: &FilletBrep,
    eu: EdgeuseId,
    op: &'static str,
) -> Result<Option<FilletGeomId>> {
    let fg = pseudo.fillet_edgeuse(eu)?.fillet_geom();
    if fg.is_some() {
        return Ok(fg);
    }
    let mate = eu1::mate(pseudo.brep(), eu)?
        .ok_or_else(|| SolidError::invalid_input(op, "edge use without mate"))?;
    Ok(pseudo.fillet_edgeuse(mate)?.fillet_geom())
}

/// The index in [`CornerEnv::solvers`] of `pFilletGeom->GetFilletSolver()`.
fn solver_of_geom(env: &CornerEnv<'_>, g: FilletGeomId, op: &'static str) -> Result<usize> {
    let si = env
        .pseudo
        .geom(g)?
        .fillet_solver()
        .ok_or_else(|| SolidError::invalid_input(op, "fillet geometry without solver"))?;
    usize::try_from(si.0).map_err(|_| SolidError::invalid_input(op, "fillet solver id"))
}

/// The surface of a rail index (`UVPos(lRailIndex)`).
const fn srf_index(rail_index: usize) -> SrfIndex {
    if rail_index == 0 {
        SrfIndex::First
    } else {
        SrfIndex::Second
    }
}

/// `pCurve->Evaluate(t, 1, TRUE, sPV)` (`SER`): the first derivative.
fn tangent_at(k: &Kernel, curve: &Curve, t: f64, op: &'static str) -> Result<Vec3> {
    curve
        .evaluate(k, t, 1, true, true)?
        .get(1)
        .copied()
        .ok_or_else(|| SolidError::invalid_input(op, "evaluation without derivative"))
}

/// `rTsectPnt.UVPos(0) = rFrom.UVPos(0); rTsectPnt.UVPos(1) = rFrom.UVPos(1)` on the
/// intersection point of vertex use `vu`.
fn copy_uv_pos(pseudo: &mut FilletBrep, vu: VertexuseId, from: &TsectPnt) -> Result<()> {
    let t = pseudo.fillet_vertexuse_mut(vu)?.tsect_pnt_mut();
    *t.uv_pos_mut(SrfIndex::First) = from.uv_pos(SrfIndex::First);
    *t.uv_pos_mut(SrfIndex::Second) = from.uv_pos(SrfIndex::Second);
    Ok(())
}

fn null_deref(reason: &'static str) -> SolidError {
    SolidError::invalid_input(OP_CONCAVE_TOPO, reason)
}

/// `new (m_pPseudoBrep) SmFilletVertex(this)` of type `SM_FV_RAIL_X_EDGEUSE` with its
/// point, original curve parameter and `SM_FIL_PROCESSED` status.
fn new_rail_x_edgeuse_vertex(env: &mut CornerEnv<'_>, point: Point3, t: f64) -> Result<VertexId> {
    let v = env.pseudo.brep_mut().add_vertex(point)?;
    let mut fv = FilletVertex::new(env.k, env.corner);
    fv.set_fillet_vertex_type(FilletVertexType::RailXEdgeuse);
    fv.set_original_t_param(t);
    fv.set_status(FilStatus::Processed);
    env.pseudo.insert_fillet_vertex(v, fv);
    Ok(v)
}

/// `SmFilletEdge::MakeFilletEdge(m_pPseudoBrep, start, end, rpNewEdge, this)` (unit
/// `blend.fillet_geom_p5`).
pub(super) fn make_fillet_edge(
    env: &mut CornerEnv<'_>,
    start: VertexId,
    end: VertexId,
) -> Result<EdgeId> {
    env.pseudo
        .make_fillet_edge(Some(start), Some(end), env.corner)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::topology::Edgeuse;

    #[test]
    fn empty_edge_use_list_is_invalid_input() {
        let k = Kernel::default();
        let mut b = Brep::new();
        let v = VertexId::from_index(0);
        let r = find_concave_edgeuses(&k, &mut b, v, &[]);
        assert!(matches!(r, Err(SolidError::InvalidInput { .. })));
    }

    /// The difftest refusal arm: a pair edge use without edge (so without radial) is the
    /// `InvalidInput` of a null dereference upstream, which the runner reports as `refused`.
    #[test]
    fn edge_use_without_edge_is_invalid_input() {
        let k = Kernel::default();
        let mut b = Brep::new();
        let eu = b.add_edgeuse_record(Edgeuse::default()).unwrap();
        let v = VertexId::from_index(0);
        let r = find_concave_edgeuses(&k, &mut b, v, &[eu]);
        assert!(matches!(r, Err(SolidError::InvalidInput { .. })), "{r:?}");
    }
}
