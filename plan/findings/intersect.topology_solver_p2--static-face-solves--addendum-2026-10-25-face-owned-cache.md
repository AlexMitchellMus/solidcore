# intersect.topology_solver_p2: face-owned surfaces reach their face's cache in solver identical-solution tests, line intersections and drop traces (2026-10-25)

Updates: plan/findings.md, "## intersect.topology_solver_p2 static face solves: answering `GetNaturalUVDomain` exposes 9 mismatches (open, PR #599, 2026-10-23)"

PR #890 (`wire/solve-handle-domains`).

Closes the open entry "intersect.topology_solver_p2 static face solves: answering
`GetNaturalUVDomain` exposes 9 mismatches" (PR #599). Its suspected cause was right: the
9 cases disagreed because the Rust solves ran over a fresh plain cache of a standalone
plane, while C++ solves over the face's `SmTrimSrfCache` (`smsurf_GetSurfaceCache`).

- **topology_solver_p2 static adapter.** `StaticScene` now returns the faces' own
  (face-owned) planes and gives the solver a face service
  (`BrepFaceServices` over an identical second `CreateFaceFromSurface` build of the
  shim's faces; `Brep::clone` does not keep the surfaces' owners, so a clone cannot
  serve). `StaticHandle` answers `GetNaturalUVDomain` for the planes. With that, the 9
  former mismatches and the 2 `SmTLSVertexFace` maximize cases match C++, and the 8
  silhouette cases (operation 14) match once `sm_GetFaceSilhouettes` and
  `CreateSilhouetteCurves` share the face's cache (upstream reads the one cache of the
  surface in both).
- **`GlobalPointSolve_0`** (the maximize branch of `GlobalPointSolve`) takes the caller's
  face context for its cache and leaf trim tests, as the curve/surface solve already did.
- **`SmCSIGlobalSolver::IdenticalSolutions`** of a face-aware solve reads its operands'
  natural interval and UV domain (`SolveOperands`, as `SmCSGlobalSolver` did): the
  `GetNaturalUVDomain` deferrals of `boolean.topology_intersector_p1` (3) and
  `primitive.primitive_creation_p2` (2).
- **`sm_ClassifyCurveSurfaceXSectPoint`** reads `pCurveOwner->GetTolerance()` of an
  edge's curve through the new `FaceQueries::curve_edge_tolerance` (answered by the
  intersector's face-cache context), in the CSI leaf solve and in
  `sm_GlobalCurveIntersectCombineSolutions`' coincidence test.
- **`SmBrep::DeleteTopologicalVertex`**: the coincidence check's curve handles carry the
  curves' natural intervals (`SmCurve::GetNaturalInterval`, 3 cases).
- **Line intersections and drop traces.** `GlobalLineIntersect`'s coincident drop gets
  the face's cache from the context when the caller holds none; the `SmPlane` override's
  base solve runs on the owning surface (upstream calls `SmSurface::GlobalLineIntersect`
  on `this`; the port ran it on an ownerless copy, so a face's trimmed cache could not
  build its tree); `sm_GlobalCurveIntersectCombineSolutions`' `DropAndTrimCurve` traces on
  the face's cache; `SmSurfaceDropCurve`'s `DropCurve` runs with point testing off, as
  `DoTrace` leaves the cache (`SmSurfaceDropCurve.cpp:1629`), so its start drop keeps the
  fast-drop solutions unclassified.

Since merging main (#945 and the `faces_to_delete` LiveBrepObjects change): the CSI
operands come from #945 (the trees' owners on every solve, `csi_identical_solutions_in`),
and the 5 planar-path cases that deferred in `SmBrep::ProjectAndTrim`'s ray
classification (`planar_path_two_lines`, `gen0069`, `gen0071`, `gen0167`, `gen0184`) run:
`primitive.primitive_creation_p2` 300 pass, 0 fail, 0 deferred, 8 known mismatches.
