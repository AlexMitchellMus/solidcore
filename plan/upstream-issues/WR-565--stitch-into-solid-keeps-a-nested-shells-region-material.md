#### WR-565. `SmStitch::StitchIntoSolid` keeps the region inside a nested closed shell as material
- **Severity:** Wrong result
- **Area:** Topology
- **Location:** `source/SMLib/src/SmStitch.cpp`, `SmStitch::StitchIntoSolid` (def. `:2659`), the `rbProducesASolid` block after the stitch loop (`:2713`).
- **What happens:** the stitch passes (`DoStitching` with region nesting) only find the infinite region and swap its void flag (`SM_VS_INFINITEREGION_SWAP`). When the faces close into two shells, one inside the other, the region inside the inner shell stays material and each shell stays in the region the gluing put it in. `StitchIntoSolid` still reports a solid. The volume counts the cavity as material, and SMLib's check can refuse the body ("Inner shell is not geometrically inside its containing outerShell").
- **Trigger:** twelve separate planar faces, the six of the cube -5..5 and the six of the cube -2..2, stitched with `StitchIntoSolid`: three regions with only the infinite region void, volume 1000 instead of 10^3 - 4^3 = 936. The SolidCraft fork's case is a profile with a hole swept round a closed path.
- **How found:** the SolidCraft SMLib fork (tangency-labs/smlib@351c57c), which after the stitch runs `FindAndSetInfiniteRegion(NULL, SM_VS_NESTEDSOLID)` with editing disabled.
- **Confirmed:** in the Rust port, which follows the pinned upstream by default (volume 1000, one void region).
- **Rust port:** `stitch_into_solid_with` (`crates/solidcore/src/heal/stitch/p2.rs`) applies the fork's step under the opt-in `Kernel::fixes.stitch_solid_is_nested_cavities_void` (off by default, so the default stays at C++ parity). Test: `crates/solidcore/tests/fork_351c57c_stitch_nested_solid.rs` (flag off: volume 1000; flag on: infinite region and cavity void, volume 936).
