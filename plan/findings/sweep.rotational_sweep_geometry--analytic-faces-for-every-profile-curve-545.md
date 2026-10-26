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
- **Recognition is tight.** The first CI run failed `sweep.api_surfaces` `gen0102`: a
  0.1 mm interpolated wave (scale 1e-5) revolved 180 degrees lay within the B-rep
  tolerance (1e-6) of an arc, so it became a torus 2.2e-7 off the C++'s B-spline
  (`revolution_points[0]`: cpp 0, rust 2.2e-7). The line/arc recognition and layout
  tests now use at most `1e-9 * (1 + size)` (`FIT_REL`), so only curves that are lines
  or arcs to rounding qualify (test `a_curve_only_near_an_arc_keeps_the_upstream_surface`).
  A real line or arc revolved through `SmApiCreateSurfaceRevolution` would still differ
  from the C++ inside the domain (the analytic surface is parameterised by angle, the
  B-spline rationally); no `api_surfaces` case revolves one.
- A torus tube range always starts in [0, 360): issue #386 reports a range starting
  below 0 is refused or built wrong.
- **Boundary check; full minor circles.** After merging main (`surface.trim_srf_cache`
  now builds face-owned caches), `create_partial_torus` failed with
  `SM_ERR_NULL_POINTER` in the topology sweep: its profile is a full minor circle from
  the bottom (tube angle -90), and the torus's iso circle on that side
  (`torus.rs` `iso_curve`, faithful) asks for `[270, 630]`, which
  `SmEllipse::AdjustSTEPInterval` collapses to one parameter (`WR-211`; the fork's
  74bed94 fixes it, the port does not). `MakeFaceWithCurves` then dropped those
  degenerate trims. `analytic_revolution` now also requires the analytic surface's
  natural boundary to have the B-spline's (the same sides degenerate, the others' lengths
  within 1e-3), so a full minor circle not starting at tube angle 0 keeps the B-spline
  (the partial torus, and `sweep.api_brep`'s full-circle profiles: those cases are back
  on main's path, so the deferral described above no longer happens). Tube ranges stay
  inside `[-360, 360]` (`tube_range`). With the merge, the whole-turn dome on its disc
  (sphere + plane) and the torus of two half circles (two tori) build closed manifold
  solids.
- `crates/solidcore/tests/primitives_solid.rs` `precise_sphere_properties_scale_and_translate`
  (from main) fails on Windows with or without the switches (`create_sphere` never
  reaches this code); main's Linux CI passes it.
- **Rebased on main after #899 (2026-10-26).** #899 (issue #826) measures the torus iso
  circle's end angles from the circle's centre, not the torus origin (`torus.rs`
  `iso_curve`, a deviation). The full minor circle from the bottom now has a whole iso
  circle, so its boundary matches the B-spline's and it becomes a torus as well: the
  partial torus (`create_partial_torus`, 90 degrees) now has an analytic torus face
  and is still a manifold solid (`api_primitives` test
  `complete_solids_construct_and_partial_solids_report_dependencies`). `sweep.api_brep`
  full-circle profiles that do not start at tube angle 0 may take the analytic path
  again; the PR's CI difftest shows their status and counts. The boundary check stays:
  it still keeps the B-spline wherever the analytic sides differ.
