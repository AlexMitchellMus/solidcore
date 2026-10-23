# curve.curve_p1: Range endpoint refinement (#527, 2026-10-23)

Unit: `curve.curve_p1`; production caller: `curve.curve_p8`; PR: #662.

The Rust port copied `sm_RefineCCSolutions`' endpoint acceptance defect from
pinned SMLib `c9c5979` (`SmCurve.cpp:193`). Both ends were compared to oriented
fraction zero. End refinement must instead compare to fraction one. The strict
`0.002` shift limit, orientation normalization, local solver, desired precision,
and deviation improvement test are unchanged.

Two library-direct regressions failed before the change and pass after it:

- A line `L(t)=(t,0,0)` and a quadratic bow `B(t)=(t,0.1*t*(1-t),0)`
  intersect at exactly 0 and 1. A candidate range `[0.0001,0.9999]` has endpoint
  deviations about `9.999e-6`. Previously only its start was polished; its end
  stayed at `0.9999`. Rust now polishes both ends below `1e-10` deviation, for
  all four independent forward/reversed curve parameterizations.
- A line `L(t)=(t,0,0)` and `C(t)=(t,0.1*t,0)` intersect only at zero.
  On a candidate range `[0.0001,0.9]`, both local solves converge to zero.
  Previously the end was accepted there too (about `5.55e-16`), collapsing the
  candidate range. Rust now accepts the small start correction and preserves
  the original end and its `0.09` deviation because the end shift is too large.

The five `refine_range_end_*` fixtures in the existing `edge_cases.json` exercise
these same inputs on both sides. Before the change they compare exactly because
both implementations share the defect. The C++ shim's file-static helper remains
an unchanged copy of the pinned upstream implementation; no reference formula or
tolerance is altered to hide the intentional Rust correction. The resulting
native disagreements are recorded explicitly in `crates/difftest/known_mismatches/curve.curve_p1.toml` and
`crates/difftest/known_mismatches/curve.curve_p8.toml`, and
`curve.curve_p1` is consequently partial rather than green until the reference
is corrected. This is an upstream defect, not a local-solver numerical mismatch.


The caller's generated `gen2160` (seed `7577089526937863746`) exposed the same
end-to-start collapse in a real `GlobalCurveIntersect`, now also saved as
`range_end_preserves_noncollapsed_tail` in the existing caller `edge.json`.
C++ returns two solutions, Rust three; their first 20 output values are identical.
The extra approximate range runs from `(6.11037890900769, 0.4336618890164721)`
to `(6.344696713196761, 0.5053108900962493)` in the two curve parameters.
An isolated native helper call on this tail confirms that C++ incorrectly moves
its end to `(6.110378909007681, 0.43366188901647057)`, at fraction approximately
zero rather than one. Its later classification discards the collapsed range.
Rust preserves the original end (deviation `0.3323583728011054`); the new
library-direct runner regression samples eleven corresponding points along the
tail and verifies they remain within the requested tolerance `0.5`.
Both the generated case and its stable fixture are explicit known mismatches;
`curve.curve_p8` is also partial. No comparison input is rejected or deferred,
and no tolerance changes. Other existing cases continue to compare exactly.


Validation on macOS arm64 / Rust 1.98: 2,072 kernel tests pass; the default
runner tests and the native-feature caller regression pass. Formatting,
workspace Clippy with and without `cpp` (`-D warnings`) and wasm32 check pass.
Native helper: 3,096 cases, 3,091 exact passes and five known mismatches;
caller: 3,069 cases, 3,067 exact passes and two known mismatches. Zero unexpected
failures or deferrals. The full native gate is not re-certified: the already
recorded macOS reference defects remain unrelated to this change.

## Review update (2026-10-26)

Main's opt-in `Kernel::fixes.refine_cc_range_ends` contains the entire correction:
there is no additional production guard. The production regression now asserts
two solutions with the default policy and three with the switch enabled, checking
the saved tail's residual at eleven samples. Five helper fixtures and the saved
production fixture remain in the differential suite with default C++ behavior.
The former known mismatches and status downgrades are no longer needed. The WR-260
record from this branch is replaced by the WR-383 addendum; earlier numerical
validation above belongs to the earlier implementation and bases.
