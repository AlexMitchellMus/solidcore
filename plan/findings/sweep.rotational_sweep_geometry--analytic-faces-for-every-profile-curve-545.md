# sweep.rotational_sweep_geometry: analytic faces for every profile curve of a revolution, issue #545 (2026-10-24)

SolidCraft change (tangency-labs/smlib 0b143f7; that commit is not public, so the port
follows the issue's description). Upstream `WR-268` (`plan/upstream-issues/WR-268--createrotationalsweep-analytic-faces-only-for-one-curve.md`).

- `sweep::analytic_revolution` gives a line or circular arc in one half-plane through the
  axis its cylinder, cone, sphere or torus, bounded to the curve and the turn, in the
  layout of `sm_CreateRotSrf`'s B-spline (u the turn, v the curve). It is only used when
  its corners, start tangent and sense of turn match that layout (against the B-spline
  itself in `EdgeSweepHigher`), so the loop `EdgeSweepHigher` builds stays the same
  (test `analytic_revolutions_keep_the_b_spline_loop`: same roles, orientations,
  tolerances and normal; the UV trims lie on the analytic surface).
- **Switches.** `PrimitiveCreation::set_analytic_revolutions` and
  `RotationalSweepGeometry::set_analytic_revolutions` are off by default, so every
  unit-level difftest compares upstream's surfaces: `sweep.rotational_sweep_geometry`
  (surfaces compared), `primitive.primitive_creation_p3` (`partial_rotation_line`
  compares the face's surface at 3x3 points; a cylinder there would differ from the
  C++'s B-spline) and `p2`. The `SM_API` entry points (`PrimitiveModel::build_new`,
  `api_primitives::new_primitive`) turn it on.
- **API-level difftests with it on.** `sweep.api_surfaces` revolution cases revolve
  interpolated waves and 2D NURBs, which get no analytic surface: unchanged.
  `sweep.api_brep` compares status and topology counts only. A local Rust run of its
  generated profile kinds (circle, rectangle of lines, two half arcs about x; 45 to 360
  degrees; capped or not; five scales) gave the same status and counts with and without
  the switch in 130 of 150 runs. The other 20, **one full circle swept part of a turn
  without caps**, built a B-spline torus face (1 face, 3 edges, 2 vertices) and now
  defer at `surface.trim_srf_cache` (`TrimSrfCache::build_tree_with_face for a
  face-owned surface`): a later step queries the analytic torus face through its
  trimmed surface cache. Expect those `api_brep` cases to move from pass to deferred;
  they come back when the trimmed surface cache is ported.
- Whole turn: a drilled hole (wall and drill point) now gives a cylinder and a cone
  face through `SmApiCreateRotationalSweep` (was two `SurfOfRevolution`). A dome on its
  disc and a torus from two half circles stitch through the trimmed surface cache and
  defer there, as their `SurfOfRevolution` faces did. Partial turns of a closed
  profile defer as before (`topology.face_p4` or `surface.trim_srf_cache`), so no
  partial solid is certified end to end yet.
- A torus tube range always starts in [0, 360): issue #386 reports a range starting
  below 0 is refused or built wrong.
