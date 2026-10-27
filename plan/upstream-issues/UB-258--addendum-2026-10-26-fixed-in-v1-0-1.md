Updates: UB-258

#### UB-258 (update, 2026-10-26). Fixed in v1.0.1

- **Fixed upstream in v1.0.1 (3fd5b95), `SmOffsetExecutive::RemoveExcessFaces`.** It collects the faces to delete (`sGoneFaces`) and their own edge uses (`sGone`) before `pResult->RemoveFaces(sDeleteFaces)`. After the removal, when any face was deleted, it calls `sm_ForgetSeconds` on `m_vFUToF`, `m_vEToF`, `m_vVToF` and `m_vBlendToFaces`, and `sm_ForgetFirsts` on `m_vBlendToFaces`. `m_vEUToEU` also forgets every pair whose result edge use is in `sGone`, or is no edge use of a remaining result edge (`sLive`). These are the new static helpers `sm_ForgetSeconds` and `sm_ForgetFirsts`.
- **Differences from the fork fix the port carried (tangency-labs/smlib 5fe723c, #368):** upstream forgets after `RemoveFaces` rather than before, and it also prunes `m_vEUToEU`.
- **Found by:** reading the c9c5979..3fd5b95 diff during the reference pin move.
- **Confirmed:** no. Read at 3fd5b95 `source/SMLib/src/SmOffsetExecutive.cpp` (`RemoveExcessFaces`); no ASan rerun of the #368 trigger at v1.0.1.
- **Rust port:** `OffsetExecutive::remove_excess_faces` follows v1.0.1, unconditionally: `removed_faces_edgeuses`, then `remove_faces`, then `forget_seconds`/`forget_firsts`. See the #368 finding's addendum.
