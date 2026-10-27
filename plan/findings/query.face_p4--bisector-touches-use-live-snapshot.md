# query.face_p4: bisector touches must use the live classification snapshot (2026-10-27)

Units: `query.face_p4`, caller `heal.face_props`. PR: #1334.

## Input and observed failure

Canonical main `fbd9b3f74` and PR #1334 head `eed0062ac` fail the same existing generated case, `heal.face_props/s3gen0112`, seed `9132544051722639607`. It constructs three semicircular loops on a swapped 180-degree surface of revolution and runs Stage 3 face properties followed by the loop snapshot. The input is retained unchanged in the existing generator.

The reports compare 166 values. Only `values.v[2]` differs: C++ `1` (`NoClassifyZeroCrosscnt`), Rust `2` (`NoClassifyOddCrosscnt`). The local matching upstream 1.0.1 reference reproduces this discrepancy. Canonical push CI `40565730025` fails its sweep on this case; child PR baseline checks consequently fail with `exact-base CI not successful`.

## First differing computation

Isolated diagnostic copies of `SmFace.cpp` and `SmCurveClass.cpp`, linked before the unchanged reference archives in a throwaway executable, and temporary Rust traces established that both sides choose identical bisector geometry and try the same six candidates. Neither has a UV guess initially. A stale UV-guess read is therefore not the cause of this input.

After the final candidate, both sides have four raw touches. Native classifies their directions as `Lamina / On / In / Lamina`; Rust used `Lamina / In / In / Lamina`. Native reduction removes the coincident run and leaves zero touches; Rust retains four and finds one target-loop crossing.

The native first vertex sector's end UV is `(0.99381936021794437, 1)`, with end tangent `(0, 1)`. Rust's stale topology view instead gives end UV `(0.7499298600583857, 0)` and tangent `(0, -1)`. This reverses the end arm's projected 3D tangent and loses the coincident departure.

`CurveOnClassify` leaves UV trim curves and tolerance changes in `CurveOnClassifyCache`'s shared B-rep snapshot. The subsequent bisector `GetLoopTouches` / `GetLoopCrossings` adapters previously resolved topology through `LoopObjectCurves` over the original B-rep. Native queries continue reading the same live topology they modified.

## Correction and limits

Wrap both bisector touch adapters with the existing `SnapshotObjects` resolver over the same classification snapshot. Vertex-sector topology and tolerances then observe the retained UV curves. Original classification-curve and surface identities stay unchanged. No additional B-rep clone, geometry approximation, status special case, fixture change, mismatch marker, skip, or comparison-tolerance change is introduced. Cached Stage 2 seam lists retain their existing handling.

The change does not alter the bisector's separate UV-guess or short-edge-culling reads; neither was demonstrated to cause this case. No broader classification cutover is claimed.

## Trace-free native verification

The final runner links the unchanged reference archives; diagnostic copies are not used for these gates. All temporary Rust tracing was removed.

```text
PASS heal.face_props cases=1873 passed=1873 failed=0 deferred=0 max_rel_err=0e0 with_values=1864 no_values=9 vacuous=0
PASS query.face_p4 cases=503 passed=503 failed=0 deferred=0 max_rel_err=0e0 with_values=476 no_values=27 vacuous=0
SENSITIVE heal.face_props perturbed=1 detected=1 unperturbed=0 deferred=0
```

The perturbation selection is the original failing `s3gen0112`. Existing non-value/status cases are not counted as geometry comparisons. These are focused local native gates, not a claim of canonical CI success or a repaired canonical-main baseline.
