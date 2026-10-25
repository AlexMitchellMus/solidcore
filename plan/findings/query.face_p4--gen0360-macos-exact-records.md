# query.face_p4: the importer's analytic approximation also flipped gen0360 on macOS; resolved by exact records (2026-10-25)

Follows "query.face_p4: two scale-1e-3 cylinder cases move with the importer's curve
approximation" above. On macOS (Apple Silicon, Apple clang 21) a third cylinder case
failed: `gen0360` (seed `4634847627957580463`, a seam-crossing band and a rectangle on an
`SmCylinder`, scale 1, swapped): `classify_loop_ios[8]` cpp 1 rust 2, `[10]` cpp 1 rust 0,
`[11]` cpp 1 rust 5, `looptree.len` cpp 17 rust 9. Linux CI passed it. The C++ side
evaluates the real `SmCylinder` and its circle edges (Apple libm `sin`/`cos`), the Rust
side imported their B-spline parts only, and on macOS the bisector landed on the other
side of a loop. It is not the sin/cos fusion of issue #357/#654: the case fails the same
with a C++ reference built with `-fno-builtin-sin -fno-builtin-cos`, and with a harness
`__sincos_stret` that calls `sin` and `cos` separately for both sides. (Apple's
`__sincos_stret` returns a different `sin` than `sin()` for about 0.5% of arguments in
`[-7, 7]`; `cos` agreed for all 10M probed arguments.)

Fix: the `query.face_p4` shim (and `query.face_p3`, which includes it) now puts the exact
database record before an analytic face surface (`sc_analytic_export.h`) and an ellipse
or circle edge curve (`sc_exact_curve_export.h`), and the shared `topology.face_p2`
importer rebuilds them with the ported readers. `gen0360`, `gen0157` and `gen0300` pass on
macOS (the latter two were known mismatches, entries removed) and are named regressions
`cyl_*` in `fixtures/query.face_p4/cases.json`. `query.face_p3` compares 197 more cases
with values (`no_values` 235 -> 38): faces and edges that were `inexact_geometry` before.
