# topology.brep_constructor: port, difftest and the STEP corpus built through SmBrepConstructor (2026-10-25)

PR #1000, which builds on #985 (`feat/step-import`).

## Unit

`SmBrepConstructor` (SmBrepData.cpp lines 78-1189) is in `topology::brep_constructor`. The
difftest cases are constructor scripts (`step_import::script`). Both kernels run them:
the C++ side in a forked child, so an upstream null dereference becomes a crash marker
that is compared with the port's constructor `InvalidInput`. Each case compares every
call's status and a full dump of the B-rep that was built.

Result: `PARTIAL topology.brep_constructor cases=158 passed=99 failed=0 deferred=59
max_rel_err=0e0 blocked_by=topology.brep_validate_pointers(59)`.
- With `--seeds 400` (522 cases): 410 passed, 0 failed, 112 deferred.
- `--perturb`: SENSITIVE 87/87.
- `--fuzz 600`: FUZZ-PASS.

All deferrals are `EndShell(TRUE)`. Upstream calls `SmBrep::ValidatePointers` first and
ignores its status, but that call creates the missing UV trim curves of every face that
already has one (`plan/upstream-issues/WR-268--...`). It is not ported; the new unit
`topology.brep_validate_pointers` tracks it.

## Corpus (`difftest --step-corpus [--constructor] --jobs 4`, 302 bodies)

| build | MATCH | HEAL-DIFF | DIFF | DEFERRED | CRASH | SKIPPED |
|---|---|---|---|---|---|---|
| faces + stitch (#985 recipe) | 177 | 101 | 8 | 13 | 3 | 4 |
| SmBrepConstructor | 126 | 141 | 30 | 5 | 0 | 4 |

Constructor path:
- Every construct, counts and manifold stage agrees between Rust and C++.
- All 30 DIFFs first differ at `area`; the stitch path has 3 area DIFFs.
- Deferrals: `SmEdge::GetTolerance` 4, `SmFace::GetTolerance` 1. They come from the
  existing curve/surface owner stubs, reached in the stages after construction.
- No C++ crash: the 3 stitch-path crashes are in `FindInfiniteRegion` during the
  stitch, and the constructor path does not stitch.

Both paths, C++:
- 262 bodies are manifold solids with the constructor, 244 with the stitch path.
- Of the 244 bodies that are solids on both paths, 129 have the same |volume| within
  1e-6, and 39 have a negative volume on one path (the stitch path leaves face
  orientation to the stitch).

The area DIFFs are not investigated yet. Constructor-built faces have no UV trim curves,
so the area integration has to create them, which the stitch path does not need.
