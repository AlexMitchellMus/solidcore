# Primitive creation p3 comparisons

`cases.json` calls the pinned native implementation, including the actual
file-local `sm_CreateRuledSurfaceFromOffsetCurve` (the shim includes the source
with renamed externally visible symbols). The generated cases use fixed harness
seeds and scales from 1e-6 through 1e6.

| Path | Evidence compared on both sides |
| --- | --- |
| Ordinary ruled surface, surface-intersection attribute transfer | 25 surface samples, three input curve states, status and flags |
| Corner to pivot, including rational weight matching | Resulting surface, offset degree/control count/rationality and samples; original generator states |
| Circular arc before/past the apex | Surface samples plus offset/generator mutations; past-apex offset is intentionally not translated |
| Missing map | Error status and all three curves after degree reduction |
| Already-used generator, pre-existing missing-face state | Error status, marriage flags, and transformed offset before the early return |
| Sliver rejection | Error and mutated inputs, with below/equal/above length-threshold fixtures |
| Invalid capped repetition counts, empty curve sweep, null blend edge uses | Entry-point status and edit/context/output-list state; these are early exits, not geometry coverage |
| Swung primitive with no profiles | Successful empty-loop return and unchanged state, not a constructed surface |

The eight nontrivial construction fixtures execute the real C++ public entry
points. Rust defers at `CreatePlanarFacesWith3DCurves`, `CreateFaceFromSurface`
or `CreateFacesFromSurface`; no topology or later transformations are compared
for those cases. All eight owned functions have full Rust bodies, including
control flow after dependency adapters. The production caller must supply a
`PrimitiveAccess` resolver extending the existing `SweepAccess` interface.
No concrete production construction adapter is supplied by this unit.

Curve-sweep scaling, profile orientation, face-based skinning, classification,
cap construction, stitching, map-copy tails and topology edits are implemented
but **not certified by these differential results**. A kernel regression checks
that preparation preserves the input curves and that the inverse placements
recover the original path. Another checks the transform/marriage mutations on
the missing-face early return. Neither replaces native comparison of the tails.

The empty skin fixture is rejected before C++ under existing **UB-78**, not
counted as a parity pass. The offset helper preflights the selected three-point
circle solve (**UB-57**) and rejects non-finite apex placement. Failed native
curve construction or non-finite fixture samples are construction rejections.
The Rust port still checks these inputs safely. No comparison tolerance or
absolute floor is added.

Legacy sweep APIs remain blocked at their earlier `SmBrep` constructor. Their
outdated p3 blocker labels now identify that actual missing dependency. Five
read-only sweep services resolve the integrated B-rep/region/shell accessors;
no arena-order replacement is used.


Review regressions (`review.json`) allow finite zero/negative rule tolerances and
zero/negative curve-sweep approximation tolerance through the checked FFI. Four
helper cases compare real values on both sides; two curve-sweep cases still defer
at topology construction. The updated gate has 293 cases: 284 pass, 8 deferred,
1 UB-78 rejection and 0 failures.

Swung scaling retains the original derived `Curve` for virtual point evaluation
and global property analysis. Two Rust regressions use imported analytic lines
with a cached NURB differing by two ULPs, checking profile evaluation without a
plane intersection and path evaluation after an intersection. Both distinguish
the correct analytic result from a sliced B-spline. They do not add a native
comparison of the deferred swung construction/topology tail.
