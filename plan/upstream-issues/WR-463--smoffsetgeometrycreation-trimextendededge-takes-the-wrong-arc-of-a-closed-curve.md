#### WR-463. `SmOffsetGeometryCreation::TrimExtendedEdge` takes the wrong arc of a closed intersection curve, so a hemisphere's shell is no solid (#403)
- **Severity:** Wrong result
- **Area:** Offset and shell
- **Location:** `source/SMLib/src/SmOffsetGeometryCreation.cpp`, `SmOffsetGeometryCreation::TrimExtendedEdge`, called by `MateEdges`.
- **What happens:** the offset of an edge that lies on a closed curve (here the circle where the offset sphere meets the copied flat face) is trimmed between the two points where the base edge's ends map. On a closed curve two points bound two arcs, and the trim always takes the interval between them. The rim of a hemisphere is two arcs of one circle (270 and 90 degrees); both offset edges get the same 270 degree arc, so the flat face holds one arc twice and none for the other, and the cavity is no solid.
- **Trigger:** a sphere r10 at the origin intersected with a box (-20,-20,0)-(20,20,20), shelled 2 inward open at the flat face (z = 0).
- **How found:** SolidCraft's regression rig (tangency-labs/solidcore#403); the intermittency the rig saw is not reproduced (the other cause, pointer-hash order in `SmMapPtrToPtr`, is absent in the port).
- **Confirmed:** yes, by the fork's capi test in tangency-labs/smlib@fa7ff4a, and in the port by unit test `issue_403_hemisphere_rim_arcs_keep_their_own_trim` (no solid, volume 0, without the fix).
- **Fix:** tangency-labs/smlib@fa7ff4a takes, on a closed curve, the arc whose midpoint is nearer the base edge's midpoint where that arc is an interval of the curve (a trim point on the seam); where it would cross the seam the edge gets no trim curve.
- **Rust port:** `trim_extended_edge` behaves as upstream by default. With `Kernel::fixes.offset_closed_curve_trim_arc` it does the same; the hemisphere shells to a solid of volume 2/3 pi (10^3 - 8^3).
