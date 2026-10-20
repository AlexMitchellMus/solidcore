//! Fillet corners, part 2: the concave edge uses of a corner (and, once the corner record
//! of `blend.fillet_corner_p3` lands, the variable-radius and analytic corner surfaces and
//! the topology adjustments of the 1x1, `Nx1`, `Nx2` and `NxN` corners).
//!
//! Upstream: `source/SMLib/src/SmFilletCorner.cpp`, `source/SMLib/inc/SmFilletCorner.h`.
//!
//! | Rust                       | Upstream                 |
//! |----------------------------|--------------------------|
//! | [`find_concave_edgeuses`]  | `sm_FindConcaveEdgeuses` |
//!
//! Deviations:
//! - **Out-param.** The concave edge uses (`rConcaveEUs`, appended to) are returned; the
//!   `iDebugLevel` argument (graphics only) is dropped.
//! - **Undefined input.** Upstream reads `(*pAllEUs)[GetSize() - 1]` unchecked, past the
//!   array when it is empty; an empty `all_eus` is [`SolidError::InvalidInput`]. Its
//!   callers pass a corner's `m_pAllEUs`, which `SmFilletCorner::Create` never leaves
//!   empty.
//!
//! The loop runs over the given edge uses once; there is no recursion.

use crate::blend::fillet_corner::p1::check_face_complimentary_corner_concavity;
use crate::error::{Result, SolidError};
use crate::kernel::Kernel;
use crate::topology::{Brep, EdgeuseId, VertexId};

const OP_CONCAVE: &str = "blend::fillet_corner::p2::find_concave_edgeuses";

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_edge_use_list_is_invalid_input() {
        let k = Kernel::default();
        let mut b = Brep::new();
        let v = VertexId::from_index(0);
        let r = find_concave_edgeuses(&k, &mut b, v, &[]);
        assert!(matches!(r, Err(SolidError::InvalidInput { .. })));
    }
}
