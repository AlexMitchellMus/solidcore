#### WR-466. The fillet point solve's finite-difference step is a process-wide static
- **Severity:** Wrong result
- **Area:** Blends
- **Location:** `source/SMLib/src/SmFilletSolver.cpp`, `SmFilletSolver::ComputePointValues` (`:3878` in the fork, `:3604` at the pinned commit): `static double sdStaticStep = (TANGENT) ? 1.0e-3 : 1.0e-5;`.
- **What happens:** the static is initialised by the first point any blend of the process computes, and every later point takes that step. A tangent point steps 1e-5 after an earlier crossing point, and a blend's result depends on the blends before it. Same inputs gave a sound body in a fresh process and an unsound one (internal pole, repeated control points, cusp, missing seam) after earlier blends.
- **Trigger:** box 60 x 40 x 20, r5 on the top edge at y = 40, r3 on the top edge at x = 60, then r2.5 or r2.9 on the vertical edge at x = 60, y = 40, after earlier blends in the process (#1173). Same record family as ND-3.
- **How found:** the fork's investigation (gdb: `.data`/.bss diff at the r2.5 call; forcing the value flips the result in both orders).
- **Confirmed:** yes, by the fork (tangency-labs/smlib@0a05860).
- **Fix:** choose the step per point.
- **Rust port:** the `Kernel::fillet_point_step` latch reproduces upstream (per kernel). With `Kernel::fixes.fillet_point_step_per_point` the step is chosen per point and the latch is not used. Test `per_point_step_ignores_the_latch_with_the_fix`.
