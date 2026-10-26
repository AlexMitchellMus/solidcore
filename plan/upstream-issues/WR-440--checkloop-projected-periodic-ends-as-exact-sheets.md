#### WR-440. `CheckLoop` treats projected periodic ends as exact UV sheets (#530)
- **Severity:** Wrong result (a valid loop is rejected with `SM_ERR`, or a closed loop is classified as closed except at singularities)
- **Area:** Topology / trimming
- **Location:** `source/SMLib/src/SmTrimmingTools.cpp`, `SmTrimmingTools::CheckLoop` (`GetEdgeUVExtremes` around 3109; large-UV-gap predicates around 3220).
- **What happens:** the raw UV jump between consecutive edge-use ends is tested against a tenth of the natural domain even when one end comes from a vertex projection, which is known only modulo the period. The jump can then produce SM_ERR 1001 with a loop-list-3/test-1 report, or a seam-only closure-1 classification at a singularity.
- **Trigger:** an endpoint without an exact UV trim is projected to the opposite side of a closed natural domain from the adjoining endpoint. Its representative differs by one period although its surface point closes in 3D.
- **How found:** issue #530; adversarial/public-file fixtures in `fixtures/topology.trimming_tools/upstream.json`, `projected_seam_cone_eu0..3` (PR #709 by fkirill, filed there as WR-264; renumbered because #700 also uses WR-264).
- **Confirmed:** yes, unchanged pinned C++ reference: eu1/eu2 report U gap 15.30733729460359 against 1.5307337294603591 and return 1001; eu0/eu3 return success with closure 1. Exact UV trim jumps must still be rejected.
