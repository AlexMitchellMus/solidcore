# surface.cone_p2: recheck the residual cylinder fix against official 1.0.1 (2026-10-27)

Updates: plan/findings/surface.cone_p2--cylinder-nurb-rebuild-issue713--addendum-2026-10-26-upstreamfixes.md

PR: https://github.com/tangency-labs/solidcore/pull/724

The original issue #713 was closed by #1077, which added the opt-in mean-radius
correction. The later maintainer review on #724 confirms a remaining defect:
when the accepted radii have been equalized, their difference is zero and the
subsequent conditional angle update leaves an earlier nonzero semi-angle intact.
This PR completes that existing opt-in; it does not reimplement the closed issue.

The current project pin is official NVIDIA usd-brep 1.0.1,
`3fd5b9548d19aa8c12c29bc9637a96d42a06f8a5`, as declared in
`scripts/upstream-pin.env`. The entire
`SmCone::RebuildSTEPFromNURBParameters` function body in
`source/SMLib/src/SmCone.cpp` is identical to the previous official pin
`c9c5979f04cea8c38213004d4f04df0dcfb40d05` (compared directly in the two checkouts).
In particular, it still uses the strict `SM_EFF_ZERO * (1 + r)` radius equality
and changes the semi-angle only for a nonzero radius difference.

The patch quoted in #713 is from the tangency-labs/smlib fork, identified there
as `f06b48a` on `sync/1.0.1-capi`; it is distinct from the official 1.0.1 release.
The quoted patch adds cylinder radius reconciliation and taper rejection, not
an unconditional angle reset after accepting an edited cylinder. There is no
new official 1.0.1 rebuild implementation to port. Main's other 1.0.1 changes,
including its ellipse-trim code, were retained by merging main `8cb29a1fa`.

The stale-angle reset therefore remains inside the existing
`Kernel::fixes.cylinder_rebuild_mean_radius` accepted-cylinder branch. The default
continues to follow the project's official pin, now 1.0.1, and the taper error
still precedes all analytic-state changes. Four added library regressions plus
main's tests (14 cone-p2 tests in total) pass after the merge. The C++ libraries
and sparse upstream checkout used for revalidation both identify commit `3fd5b95`.
