# TODO: restore skip_cone_seams for unbounded ST-ACIS cones

- **Area:** `step-brep` / `step-import`
- **Fixtures:** NIST/ST-ACIS `1797609in.stp`, `2827056.stp`
- **Symptom:** display mesh collapses into a folded “pretzel” because generated
  seams on unbounded non-cylindrical cones run out to SolidCore’s infinity
  marker (±1e6).
- **Previous solution:** `BuildOptions { skip_cone_seams: true }` and
  `build_without_cone_seams(body)`. The rewrite that removed the vendored STEP
  snapshots also removed this option.
- **Current SolidNode workaround:** detect infinity-marker coordinates in the
  display mesh, rebuild with the vendored `build_without_cone_seams`, validate,
  and remesh.
- **Required:** re-expose `skip_cone_seams` in upstream `BuildOptions` and
  `build_without_cone_seams`, or provide an equivalent fallback API. Then remove
  the vendored STEP snapshots from SolidNode.
- **Tests:** both fixtures must produce recognizable, bounded solids; assert a
  non-trivial triangle count (the current test only requires a non-empty mesh).
