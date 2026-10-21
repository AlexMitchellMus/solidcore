# Issue #294: private edge revolution transport, Phases 0–1

The isolated integration base is `c4241a48` (main `7b57d8de`, PR #356
`b338b918`, PR #304 `f7ba5d5e`). Both dependency histories and authorship are
preserved by merges. #356 is unchanged. Ownership was refreshed on 2026-10-21
at 07:51 UTC: all 86 open PR file lists, relevant PR discussions and the #76/#294
issue discussions were checked. No private edge transport implementation was
claimed in that inventory. #313's refreshed head was `6d59d6fa`; its shared
dispatch ownership, the analytic variant follow-ups and #136's DropCurve work
remain separate. Only #304 and #356 were needed for this integration.

The baseline report and both runner outputs were captured before the adapter
change. A diagnostic copy of the unchanged edge importer retained native type
codes before flattening, resolved the selected edge and enumerated its incident
faces. The 186-case import manifest, six separate DropCurve cases, full before
and after outputs, comparison script and logs are retained locally under
`bin/edge-coverage/`.

## Implementation

The private edge shim calls #304's `sc_analytic_export::put_exact_head` before
the regular B-spline block. The private importer uses `read_exact_surface` and
`check_exact_base`, without changing either shared helper. The native helper's
NURB-presence/non-mutation checks, native binary record round trip and polar
curve checks are retained. Analytic domains, generator, polar curve and the
surface/generator out-of-bounds flags travel through the existing format. The
reconstructed base must match the independent regular block bitwise. Existing
face ownership/tolerance metadata and fixture input layouts are preserved.

A second transport head where the regular block belongs is rejected without
recursive parsing. Native regressions cover the 72 revolution operations and
damaged bytes, truncated/trailing records, mismatched/missing bases and nested
heads, alongside the existing exact line/B-spline/plane/same-B-rep controls.

## Default differential results

| Measurement | Integrated baseline | Private adapter |
| --- | ---: | ---: |
| Cases | 1787 | 1787 |
| Passes with compared values | 1595 | 1667 |
| Empty passes / both errored | 0 | 0 |
| Import deferrals | 186 | 114 |
| DropCurve deferrals | 6 | 6 |
| Failures / rejections | 0 | 0 |
| Maximum relative error | 0 | 0 |

All 72 revolution cases now execute both operations and compare 120 output
values. Every one was formerly C++ `inexact_geometry` with empty outputs and a
Rust `topology.edge` import deferral; neither operation executed at baseline.
Every one now passes with nonempty outputs, no operation error or deferral, and
zero relative error. The other **1715 case records and both runner outputs are
unchanged**, including all 1595 former passes. Relative tolerance remains
`1e-9`; the absolute-floor list remains empty.

These operations return predicates, traversal indices and statuses. In particular,
the topological result is false and its distance remains the upstream `1e20`
sentinel; no newly computed nonzero geometric distance is claimed. Blend returns
the null edge-use index `-1`. The table below records the actual `state` arrays
from both engines, including these outcomes.

Existing controls still compare 4 values (`line_join_fw_tangent_start`), 18
(`line_join_fw_trimmed`), 18 (`net_trimmed_e0`), 1
(`plane_square_e0_tangent`) and 50 (`Mouse_e0_query`).

## Per-case before/after evidence

For **each** row, before = C++ refusal / Rust import deferral; after = genuine
pass with the identical C++ and Rust `state` arrays shown. `_seam` and `_seam_f`
are the two seam forms; `_rb_` uses rebuilt UV trim curves. The native selected
edge curves are B-splines (`14033`), and their incident face is an exact
`SmSurfOfRevolution` (`15020`), including the `torus.smb` fixtures.

| Case | After C++ state = Rust state | Values compared |
| --- | --- | ---: |
| `torus_e0_seam` | `[1.0]` | 1 |
| `torus_e0_seam_f` | `[1.0]` | 1 |
| `torus_e0_topological` | `[0.0, -1.0, 1e+20]` | 3 |
| `torus_e0_tangent` | `[1.0]` | 1 |
| `torus_e0_blend` | `[-1.0]` | 1 |
| `torus_e0_radial` | `[0.0, 2.0, -1.0]` | 3 |
| `torus_e2_seam` | `[1.0]` | 1 |
| `torus_e2_seam_f` | `[1.0]` | 1 |
| `torus_e2_topological` | `[0.0, -1.0, 1e+20]` | 3 |
| `torus_e2_tangent` | `[1.0]` | 1 |
| `torus_e2_blend` | `[-1.0]` | 1 |
| `torus_e2_radial` | `[0.0, 2.0, -1.0]` | 3 |
| `torus_e5_seam` | `[1.0]` | 1 |
| `torus_e5_seam_f` | `[1.0]` | 1 |
| `torus_e5_topological` | `[0.0, -1.0, 1e+20]` | 3 |
| `torus_e5_tangent` | `[1.0]` | 1 |
| `torus_e5_blend` | `[-1.0]` | 1 |
| `torus_e5_radial` | `[0.0, 1.0, -1.0]` | 3 |
| `torus_e9_seam` | `[1.0]` | 1 |
| `torus_e9_seam_f` | `[1.0]` | 1 |
| `torus_e9_topological` | `[0.0, -1.0, 1e+20]` | 3 |
| `torus_e9_tangent` | `[1.0]` | 1 |
| `torus_e9_blend` | `[-1.0]` | 1 |
| `torus_e9_radial` | `[0.0, 1.0, -1.0]` | 3 |
| `torus_e11_seam` | `[1.0]` | 1 |
| `torus_e11_seam_f` | `[1.0]` | 1 |
| `torus_e11_topological` | `[0.0, -1.0, 1e+20]` | 3 |
| `torus_e11_tangent` | `[1.0]` | 1 |
| `torus_e11_blend` | `[-1.0]` | 1 |
| `torus_e11_radial` | `[0.0, 1.0, -1.0]` | 3 |
| `torus_rb_e0_seam` | `[1.0]` | 1 |
| `torus_rb_e0_seam_f` | `[1.0]` | 1 |
| `torus_rb_e0_topological` | `[0.0, -1.0, 1e+20]` | 3 |
| `torus_rb_e0_tangent` | `[1.0]` | 1 |
| `torus_rb_e0_blend` | `[-1.0]` | 1 |
| `torus_rb_e0_radial` | `[0.0, 2.0, -1.0]` | 3 |
| `surf_rev_e0_seam` | `[0.0]` | 1 |
| `surf_rev_e0_seam_f` | `[0.0]` | 1 |
| `surf_rev_e0_topological` | `[0.0, -1.0, 1e+20]` | 3 |
| `surf_rev_e0_tangent` | `[0.0]` | 1 |
| `surf_rev_e0_blend` | `[-1.0]` | 1 |
| `surf_rev_e0_radial` | `[0.0, 0.0, -1.0]` | 3 |
| `surf_rev_e2_seam` | `[0.0]` | 1 |
| `surf_rev_e2_seam_f` | `[0.0]` | 1 |
| `surf_rev_e2_topological` | `[0.0, -1.0, 1e+20]` | 3 |
| `surf_rev_e2_tangent` | `[0.0]` | 1 |
| `surf_rev_e2_blend` | `[-1.0]` | 1 |
| `surf_rev_e2_radial` | `[0.0, 2.0, -1.0]` | 3 |
| `surf_rev_e5_seam` | `[0.0]` | 1 |
| `surf_rev_e5_seam_f` | `[0.0]` | 1 |
| `surf_rev_e5_topological` | `[0.0, -1.0, 1e+20]` | 3 |
| `surf_rev_e5_tangent` | `[0.0]` | 1 |
| `surf_rev_e5_blend` | `[-1.0]` | 1 |
| `surf_rev_e5_radial` | `[0.0, 1.0, -1.0]` | 3 |
| `surf_rev_e9_seam` | `[0.0]` | 1 |
| `surf_rev_e9_seam_f` | `[0.0]` | 1 |
| `surf_rev_e9_topological` | `[0.0, -1.0, 1e+20]` | 3 |
| `surf_rev_e9_tangent` | `[0.0]` | 1 |
| `surf_rev_e9_blend` | `[-1.0]` | 1 |
| `surf_rev_e9_radial` | `[0.0, 1.0, -1.0]` | 3 |
| `surf_rev_e11_seam` | `[0.0]` | 1 |
| `surf_rev_e11_seam_f` | `[0.0]` | 1 |
| `surf_rev_e11_topological` | `[0.0, -1.0, 1e+20]` | 3 |
| `surf_rev_e11_tangent` | `[0.0]` | 1 |
| `surf_rev_e11_blend` | `[-1.0]` | 1 |
| `surf_rev_e11_radial` | `[0.0, 3.0, -1.0]` | 3 |
| `surf_rev_rb_e0_seam` | `[0.0]` | 1 |
| `surf_rev_rb_e0_seam_f` | `[0.0]` | 1 |
| `surf_rev_rb_e0_topological` | `[0.0, -1.0, 1e+20]` | 3 |
| `surf_rev_rb_e0_tangent` | `[0.0]` | 1 |
| `surf_rev_rb_e0_blend` | `[-1.0]` | 1 |
| `surf_rev_rb_e0_radial` | `[0.0, 0.0, -1.0]` | 3 |

## Remaining deferrals and ownership

The 114 remaining import deferrals are unchanged: 42 cone cases (`cone.smb`
and six `box_with_large_fillet.smb` cases), 36 sphere cases (`sphere.smb`), and
36 extrusion cases (`Mouse.smb`). Their original unsupported incident-face
types are `15012`, `15014` and `15030` respectively. The existing native refusal
and Rust import deferral remain in effect; no B-spline stand-in is newly accepted.

The six unchanged `intersect.surface_drop_lib_p3` deferrals are
`Endcap_e0_uv_of_surface`, `Endcap_e2_uv_of_surface`,
`Endcap_e5_uv_of_surface`, `Endcap_e9_uv_of_surface`,
`Endcap_e11_uv_of_surface`, and `Endcap_rb_e0_uv_of_surface`.

No shared dispatch change or new numerical port was needed by these 72 cases.
The unit remains **partial** with 120 deferrals. Sphere/extrusion variants,
cone integration and DropCurve remain outside this bounded slice.

## Validation

The first C++ Clippy run exited 101 on `manual_let_else` in the new test helper.
That test-only destructuring was corrected; no lint was suppressed. The native
test build was deliberately interrupted (exit 143) to restart with that fix. The initial
lint failure and interrupted-build logs are retained separately from final checks.

All final checks exited **0**:

| Check | Final result |
| --- | --- |
| `cargo fmt --all --check` | Pass, including rerun after the test fix |
| `cargo clippy --workspace --all-targets -- -D warnings` | Pass |
| Same Clippy gate with `--features difftest/cpp` | Pass after the test-only fix |
| Native `topology_edge::tests`, difftest profile, C++ enabled, one test thread | 1 passed; covers all 72 cases, exact controls and damaged records |
| `cargo test -p solidcore` | 1551 passed; zero failed or ignored |
| `cargo check -p solidcore --target wasm32-unknown-unknown` | Pass |
| Default edge differential on the final source tree, 16 seeds | 1667 passes with values, 120 deferred, zero failures or empty passes |
| Edge `--perturb` | 1667 perturbed, 1667 detected, zero unperturbed, 120 deferred |
| Edge `--fuzz 2000` | 107 rejected, 235 malformed, 1658 ran, zero differences, crashes, panics or timeouts |
| Full compiled differential runner, `--green --jobs 4` | 236 units: 116 green, 120 partial, zero failed |

The final default report was regenerated after the test fix and agrees with the
measured after report. No numerical disagreement or operation failure was exposed
in the 72-case slice. Fuzz differences and all abnormal outcomes were counted
separately rather than treated as passes.

Build jobs were limited to two and differential workers to four; Cargo/C++
commands ran under `nice -n 10`. Existing native reference libraries were reused.
Only the edge report was regenerated in the versioned report directory; the full
gate summary and raw validation logs remain local under `bin/edge-coverage/`.
No tolerance or guard was relaxed, and nothing was published or pushed.
