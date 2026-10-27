# blend.fillet_executive_p3: sweep5 sphere_diff_3cylinders e29 vertex tolerance (2026-10-27)

Unit: `blend.fillet_executive_p3`. Fixture: `fixtures/blend.fillet_executive_p3/sweep_round5.json`.

Input: `sweep5-sphere_diff_3cylinders-e29-r1.0` (`Fillet/sphere_diff_3cylinders.smb`, edge 29,
circular G1, radius 1.0, thumbweight 1).

Values: everything agrees except `result_geometry[111]`, the tolerance of result vertex 27:
C++ 9.100516385135238e-5, Rust 9.304666970990486e-5 (rel 2.2e-2). The rail trims agree
(rail 0 (82.647291068524822, 102.85355699527092), rail 1 (82.647290569503795,
102.85355973819911), both sides). So do the per-stage UV-trim counts and lengths from
CreateFilletBrep to the end of merge.

Suspected cause: a tolerance raise from a UV trim drop (`SmEdgeuse::CreateUVTrimCurve` /
`SetUVTrimCurve`) that gets a different gap, as in
`plan/findings/blend.fillet_executive_p3--bottle-e29-r0.1-vertex-tolerance.md`. That one was
a `Surface::clone` standing in for upstream's shared surface object. The remaining
`clone_keeping_owner` callers that addendum lists are candidates. Not localised yet.
