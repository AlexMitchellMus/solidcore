Updates: plan/upstream-issues/UB-256--offset-piecewise-merge-deleted-face-uaf.md

#### UB-256 (`PiecewiseMerge` deleted-face UAF) (update, 2026-10-26). Fixed in v1.0.1

- **Fixed upstream in v1.0.1 (3fd5b95):** `SmOffsetExecutive::RemoveExcessFaces` now removes the deleted faces from `m_vFUToF`, `m_vEToF`, `m_vVToF` and `m_vBlendToFaces` (both directions for the blend relation) right after `RemoveFaces`. It also prunes `m_vEUToEU` of the deleted and dead edge uses. `PiecewiseMerge` therefore no longer iterates freed `SmFace` pointers (UB-258 addendum).
- **Consistent with our earlier C++ measurement:** shelling a 30x20x10 box with all twelve edges filleted at r3 crashed in 6 of 25 heap layouts at 1.0.0 and succeeds in all 25 at 1.0.1. That measurement did not identify which 1.0.1 change fixed it; this forgetting is the 1.0.1 change that removes the freed faces from the relations `PiecewiseMerge` reads.
- **Found by:** reading the c9c5979..3fd5b95 diff during the reference pin move.
- **Confirmed:** no. Read at 3fd5b95 `SmOffsetExecutive.cpp` (`RemoveExcessFaces`); case `tiny_tolerance_sheet_above_boundary` of `offset.offset_executive_p2` not yet rerun against the v1.0.1 reference.
- **Rust port:** follows v1.0.1 (`OffsetExecutive::remove_excess_faces`).
