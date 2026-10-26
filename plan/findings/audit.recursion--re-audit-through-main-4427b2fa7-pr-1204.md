# audit.recursion: re-audit of the main advance `92454ca6..4427b2fa7` (PR #1204, 2026-10-26)

Re-audit of the whole crate, following the PR #273 audit and the #835 source review.
Range: base `92454ca6` (the main that #835's source review screened through) to
`4427b2fa7` (main tip at review time; CI run 40378725048 green there: 557,740 cases,
0 fail, 0 deferred). `git diff --name-only` lists 215 files under
`crates/solidcore/src`, +45729/−5424, merged through 273 PR merge commits
(`git log --merges --oneline 92454ca6..4427b2fa7`; highlights: #1196 the S1 B-rep
writer and part-file reader, #1159 the
`refine_cci_to_surfaces` port, #1195 sweeps batch, #1109 fillet Jacobian overflow
guards, #1167/#1155 primitive-p2 guard paths, and a large set of `fix`/`wire`/
difftest-corpus batches).

Question asked of the delta: do the merged PRs introduce an input-driven recursive or
mutually reentrant traversal, do they change the bounds of existing ones, and do the
existing explicit depth/work bounds still cover the new callers without changing the
visitation order?

Method: enumerate the merged PRs; screen every changed source file for self- and
mutual name matches, added `while` loops and reader/nesting lines; read and verify
every hit; verify the existing bounds at the range tip. 172 of the 215 files carried
any signal; 112 of those carried a recursion/loop/reader hit, and every such hit was
read and resolved (accounting below).

### Bounds history: no bound changed inside the range

All three nesting bounds were already 32 at the range start, and none moved inside it:

- `MAX_GEOMETRY_NESTING = 32` (`surface/surface/p2.rs:552` at both ends) — the 64→32
  reduction is `aa45b3aa0` (the PR #273 re-audit), an ancestor of `92454ca6`.
- `MAX_SURFACE_NESTING` (`surface/offset_surface.rs:1929` at the tip) is the alias
  `= MAX_GEOMETRY_NESTING` at both `92454ca6` (`:1836`) and `4427b2fa7` (`:1929`);
  its 64→32 reduction (`bfdd1bdad`) likewise predates the range.
- `MAX_READ_NESTING = 32` (`topology/database_io.rs:701`) predates the range (#308).

Earlier reviews must not re-report the 64→32 reductions as in-range changes; the
`git log -G` value-pattern pitfall (`git log -L` and patternless `-G` under-report
merge commits) made the `git show <commit>:<file>` check above the tie-breaker.

### New recursion in the range (three paths, all bounded)

**`refine_cci_to_surfaces`** (`intersect/curve_class/p2.rs`): the recursive call itself
is new in the range (`97a52695c`, #1159, the large-gap refinement). The bound and the
counter predate it: `MAX_REFINE_RECURSION = 5` and the `recursion_count` field were
already at `92454ca6` (`:102`), where the large-gap part was still
`UNPORTED(intersect.curve_class_p2)`.

- `MAX_REFINE_RECURSION = 5` (`:126`), an explicit `recursion_count` parameter bumped
  first (`:1053`) and checked at entry (`:1056`, returns the unchanged output past the
  bound) — mirroring upstream's `lRecursionCount`/`lMaxRecursion`, with the bound made
  explicit rather than left to the caller.
- The recursive call in the large-gap path (`:2063`, `refine_cci_to_surfaces_mut`) is a
  single tail call that shares the already-incremented counter; it is not inside a
  loop (it sits in the `if pt_was_refined1 || pt_was_refined2` block, `:2002`-`:2073`;
  the function's only `while` ends at `:1629`), so one input reaches at most 5 levels.

**`glue_with_notification_identity`** (`topology/brep/p5.rs:1136`, new in `312acd985`,
#947): without a notification environment it re-enters itself once through
`plain_edit` with `Some(env)` (`:1146`-`:1147`); that call takes the `Some` branch and
does not recurse again. Depth 1, not input-driven.

**`has_seam_problem_mut` ↔ `FaceProps::set_props_stage2_with`**
(`topology/face/p4.rs:1846`, `heal/face_props.rs:1886`, both new in `a24a2434b`, #958):
when the face props are below `Stage2A` (`p4.rs:1924`), `has_seam_problem_mut` runs
`set_props_stage2_with` (`:1932`), whose seam callback (`face_props.rs:2062`) calls
`has_seam_problem_mut` again (`p4.rs:1942`) with those props. `set_props_stage2_with`
sets `prop_stage = Stage2A` (`:1918`) before the callback runs, so the inner call's
`prop_stage < Stage2A` check (`p4.rs:1924`) skips the stage-2 run (and
`set_props_stage2_with`'s own guard, `face_props.rs:1896`, would return at once): depth
2, as upstream
("increment PropState here - prevent infinite loops"). The bound is that ordering, not
a counter (see Q4).

### The S1 writer/reader (#1196)

- **Writer** (`topology/brep_data_write.rs`, `topology/geometry_write.rs`):
  `BrepData::write_to_db` → `write_curve_to_db`/`write_surface_to_db` → the
  `geometry_write` entry points, which start at depth 0 and nest through
  `write_curve_at`/`write_surface_at`, each taking an explicit `depth`, calling
  `check_depth` at entry and returning `InvalidInput` past `MAX_GEOMETRY_NESTING`.
  The increment happens at the genuine nesting points (generator curves at `:256`
  and `:320`, the surface-of-revolution polar curve at `:263`, `depth + 1`);
  decompositions that are one level (cone → surface of revolution) pass the same
  depth through.
- **Part reader** (`topology/part_file.rs`, 928 lines): no recursion of its own.
  Nested curve reads go through `Curve::read_from_db`, which carries its own
  `enter_nested`/`leave_nested` (the composite branch at `brep_data.rs:1858` is
  covered by the read-nesting bound); every file count passes a `checked_count`
  (`MAX_PART_ITEMS = 1_000_000`, `MAX_PART_BREPS = 100_000`,
  `MAX_ATTRIBUTE_VALUES = 1_000_000`) and `ensure_available` before any allocation.

### The remaining new code, verified bounded

- `memo_pure` (`curve/p5.rs`): explicit stack that bails past
  `MAX_GEOMETRY_NESTING`.
- `geom/tree.rs` `side()`: explicit stack with a `depth > self.nodes.len()` cyclic
  guard.
- `read_from_db_with_curve_object` (`curve/ellipse/p1.rs`,
  `surface/surf_of_revolution/p3.rs`): the nested reads go through
  `Curve::read_from_db`'s own nesting; the stand-alone
  `PolarConversion::read_from_db` (`:489`) rejects the curve-owning case instead of
  recursing.

### New non-recursion (input/size) bounds added in the range, verified

- `MAX_DERIVS_ARRAY_BYTES` = 64 MiB (`nurbs/basis_adv.rs`)
- `BASIS_DERIVS_MAX_BYTES` = 8 MiB (`nurbs/basis_basic.rs`)
- `MAX_BOUNDARY_SAMPLES` = 65 536 (`surface/surface/p6.rs:106`, enforced at `:1212`
  before the `while ii <= num_samples` loop)
- `MAX_SHELL_SPLITS` = 1024 (`topology/brep/p3.rs:2581`)
- `MAX_WALK` = 100 000 (`topology/brep/validation/pointers.rs:37`, checked at
  `:101`/`:199`)
- the part-file count bounds above

### Bounds table at `4427b2fa7`

| Bound | Value | Location |
|-------|-------|----------|
| `MAX_GEOMETRY_NESTING` | 32 | `surface/surface/p2.rs:552` |
| `MAX_SURFACE_NESTING` | = `MAX_GEOMETRY_NESTING` | `surface/offset_surface.rs:1929` |
| `MAX_READ_NESTING` | 32 | `topology/database_io.rs:701` |
| `MAX_SOLVE_DEPTH` | 1 000 | `geom/global_solver.rs:117` |
| `MAX_SOLVE_BRANCHES` | 1<<20 | `geom/global_solver.rs:122` |
| `MAX_MATE_DEPTH` | 64 | `blend/fillet_geom/p3.rs:118` |
| `MAX_SPLIT_AT_SEAM_DEPTH` | 8 | `topology/face/p1.rs:178` (check at `:2113`) |
| `MAX_LOOP_NESTING` | 256 | `topology/brep/p1.rs:177` |
| `MAX_LOOP_MARK_CALLS` | 1<<20 | `topology/brep/p1.rs:1238` |
| `MAX_REFINE_RECURSION` | 5 | `intersect/curve_class/p2.rs:126` (new, #1159) |

The `MAX_SOLVE_DEPTH` per-call cap of #835's table was re-verified: `solve_branch`
(`global_solver.rs:1665`) bails with `InvalidInput` past `MAX_SOLVE_DEPTH + 1` nested
calls (depth increments at `:1671`, restores at `:1673`), and the unit test at
`:2424` builds a `MAX_SOLVE_DEPTH + 1`-level chain and asserts the error. The static
`check_solve_depth` (`:679`) rejects cyclic child links and over-deep trees before
any recursion, and `MAX_SOLVE_BRANCHES` bounds total work, not just depth.

### Scanner-hit triage (the 112 recursion/loop/reader hits, all resolved)

Every self/mutual name match, mutual pair and added `while` loop in the 112
signal-carrying files was read and resolved as either a same-name artifact (trait
declaration vs implementation,
inherent vs trait method, free function vs method, or two unrelated types sharing a
name), an existing bounded traversal, or the new code above. Notable clusters:

- `blend/fillet_executive/*` (`p1..p3`, `session*`, `topology_changes`): the
  `mapping_*`/`creation_smooth`/`smooth_target_brep`/`finish_*`/`tag_merge_*` matches
  are the session's trait-vs-inherent pairs; the stage `while let Some(x) = vec.get(index)`
  loops are the same index-increasing executive-owned-vector loops #835 already
  bounded, and `topology_changes.rs` loops over bounded record lists.
- `boolean/merge/p4/arena_pair.rs`, `boolean/topology_intersector/p1.rs`,
  `query/face/p2..p3.rs`, `topology/face/p2..p4.rs`, `topology/brep/p10.rs`,
  `offset/offset_executive/edgeuse.rs`, `heal/heal_data/p3.rs`, the
  `fillet_*` accessor clusters and `query/brep.rs`: same-name accessor /
  trait-declaration-vs-implementation matches across the entity-handle traits; no
  recursion.
- `topology/topology.rs` and `topology/topology_traverser.rs` MUTUAL pairs:
  trait declaration vs implementation.
- The added `while` loops: `intersect/curve_class/p2.rs` (the bounded refine loop
  above), `offset/offset_executive/crossing_loops.rs` (`while !*at(&done, tt)` marks
  `done` each iteration), `surface/surface/p2.rs` (`while found`, each join removes a
  loop), `surface/surface/p7.rs` (`while step < 0.1`, at most ×10 from 1e-8),
  `topology/brep/p3.rs` (`MAX_SHELL_SPLITS` inner, monotonic-index outer),
  `topology/brep/validation/spatial.rs` (explicit `pending` stacks),
  `nurbs/srf_adv/p3.rs` (test-only grid loop, `dim < 10_000` assert),
  `curve/curve/p5.rs` (the bounded explicit stack), `geom/tree.rs` (the bounded
  explicit stack) — all bounded or finite.

### Open questions (not certifications)

- **Q1 — the refine bound is a port decision.** Upstream's `lMaxRecursion` is the
  same 5; the port makes the counter explicit and returns the unchanged output past
  the bound instead of continuing. When `intersect.curve_class_p2` is next re-run,
  a fixture that drives `refine_cci_to_surfaces` to the bound (five levels of
  large-gap solutions) would pin the bail path.
- **Q2 — writer depth is writer-managed.** The S1 writer's bound is enforced at
  `write_curve_at`/`write_surface_at` entry from the depth they are given; the
  increment is at the nesting sites above. A future write path that composes a
  nested geometry without incrementing would write at the parent's depth rather
  than fail. The read side does not have this shape (its nesting is per-call in
  `read_from_db`).
- **Q3 — work, not stack.** `MAX_SOLVE_DEPTH`/`MAX_MATE_DEPTH` bound call depth;
  the total number of solver/mate evaluations per stage is bounded by
  `MAX_SOLVE_BRANCHES`/the stage's own vector sizes, not by the depth counters.
  Carried over from #835, unchanged by this delta.
- **Q4 — the seam/stage-2 recursion is bounded by ordering.** `has_seam_problem_mut` ↔
  `set_props_stage2_with` stops at depth 2 only because `prop_stage = Stage2A` is set
  before the seam callback runs. A change that moves that assignment after the callback,
  or a caller that passes props reset below `Stage2A` into the callback, would make the
  pair recurse without bound. No test pins the depth; a unit test asserting a single
  inner `has_seam_problem_mut` per stage-2 run would.

### Verification

Checks run on `4427b2fa7` (this PR is notes-only; no kernel change):

- `cargo fmt --all --check`: clean.
- `cargo clippy --workspace --all-targets -- -D warnings`: clean.
- `cargo test -p solidcore`: 2404 passed, 0 failed, 3 ignored.
- `cargo check -p solidcore --target wasm32-unknown-unknown`: clean.
- `python3 harness/plan/units.py check`: 262 units OK.
- Targeted difftest (local; the C++ reference built from source with Apple-Clang —
  the FFI shims needed `-Wno-c++11-narrowing`, a toolchain difference: Apple-Clang
  hard-errors the C++11 narrowing that the CI's GNU toolchain only warns — plus
  `CARGO_BUILD_JOBS=2`, `RUSTC_WRAPPER=sccache`, `nice -n 10`, load 1.8):
  - `surface.volume_p2` (the 64→32-rescaled fixtures): PASS, 11220/11220, 0 deferred.
  - `topology.brep_data` (S1 reader): PASS, 2379 cases, 1006 with values,
    1373 `skipped_platform`, 0 deferred.
  - `topology.part_file` (S1 part reader, out-of-plan corpus): PASS, 76/76,
    0 deferred.
  - `intersect.curve_class_p2` (the new `refine_cci_to_surfaces`): PASS,
    1671/1671, 0 deferred.
  - `intersect.curve_class_p5`: PASS, 1411/1411, 0 deferred.

  `difftest: 5 unit(s), 5 passed, 0 partial, 0 failed`; all five with
  `max_rel_err=0e0`. The CI's own run at `4427b2fa7` (40378725048) is the
  authoritative all-green baseline (557,740 cases, 0 fail, 0 deferred).

### Scope limits

- The 215 changed files were screens-first; the 112 recursion/loop/reader-signal
  files were read in full, the remaining 60 (function-list signal only, no recursion
  hit) where the changed call chains enter them. This is a source re-audit of
  the delta on top of the #273 whole-crate audit, not a fresh whole-crate audit.
- No worst-case wasm stack use is measured by this review; the wasm check (below)
  only establishes that the crate still compiles for `wasm32-unknown-unknown`.
- Files owned by other lanes were not modified by this PR; it changes no kernel
  code, no tolerance and no unit status (`audit.recursion` stays green, notes-only).
- Unmerged contributor branches were not treated as landed.

### Changed source files with signal (172; the 112 recursion/loop/reader-hit files all triaged)

`blend/fillet_corner/{p1,p2,p3,topology}.rs`, `blend/fillet_executive/{api,p1,p2,p3,
session,session/finish,session/insertion,session/insertion/corner,session/insertion/
curves,session/insertion/geometry,session/insertion/geometry/edits,session/insertion/
geometry/intersector,session/insertion/geometry/proto,session/insertion/geometry/
query,session/insertion/geometry/topology_objects,session/insertion/live,session/
original_merge,session/original_trim,topology_changes}.rs`,
`blend/fillet_geom/{p1,p2,p4,p6}.rs`, `blend/fillet_intersector.rs`,
`blend/fillet_solver/{p1,p1/rails,p2,p3,p4,p5/tests,surfaces}.rs`,
`blend/fillet_standard_solver/p2.rs` (fillet cluster),
`boolean/merge/{p1,p2,p3,p4,p4/arena_pair,p4/kernel_adapter,p5}.rs`,
`boolean/topology_intersector/{p1,p2,p3}.rs`, `curve/crv_on_surf.rs`,
`curve/curve/{p1,p2,p5,p6}.rs`, `curve/curve_cache.rs`, `curve/ellipse/{p1,p2}.rs`,
`curve/line.rs`, `geom/{assert_array,attribute,context,geom_utility/p1,hash_table,
local_solve1d,local_solve_nd,matrix,object,tol,tree}.rs`,
`heal/face_props.rs`, `heal/heal_data/p3.rs`, `heal/stitch/p1.rs`,
`intersect/adv_surface_intersector.rs`, `intersect/curve_class/{p1,p2,p3,p5,p6,p7,p9}.rs`,
`intersect/point_class.rs`, `intersect/ray_tracer.rs`, `intersect/surface_drop_curve.rs`,
`intersect/surface_intersector/p1.rs`, `intersect/topology_solver/p1.rs`, `kernel.rs`,
`nurbs/{b_spline_curve/{p1,p2,p3},b_spline_curve_n_lib,b_spline_surface/p1,
b_spline_surface_n_lib/p2,basis_adv,basis_basic,bezier,crv_basic,knots,nurbs_crv,
nurbs_srf,spiral,srf_adv/{p3,p4},srf_basic}.rs`,
`offset/offset_executive/{crossing_loops,edgeuse,merge,p1,p2}.rs`,
`offset/offset_geometry_creation.rs`, `primitive/primitive_creation/{p1/primitives,
p2/planar_path,p3/geometry,p3/sweeps}.rs`, `query/api_queries.rs`, `query/brep.rs`,
`query/brep/ray.rs`, `query/brep/ray/live.rs`, `query/face/{p1,p2,p3}.rs`,
`query/poly.rs`, `surface/{offset_surface,plane/p1,sphere,surf_of_revolution/{p1,p2,p3},
surface/{p1,p2,p3,p5,p6,p7,p8,p9},surface_cache/{face_context,p1,p3},
surface_point_solve/p2,torus,trim_srf_cache}.rs`, `sweep/topology_sweep/{store,
store/classification,store/merge}.rs`, `topology.rs`,
`topology/brep/{p1,p2,p3,p4,p5,p7,p9,p10,validation,validation/coincidence,
validation/pointers,validation/spatial}.rs`, `topology/brep_data.rs`,
`topology/brep_data_write.rs`, `topology/database_io.rs`, `topology/edge.rs`,
`topology/edgeuse/p1.rs`, `topology/face/{p1,p2,p3,p4}.rs`,
`topology/geometry_write.rs`, `topology/loopuse.rs`, `topology/part_file.rs`,
`topology/proto_topology/p1/proto.rs`, `topology/proto_topology/{p2,p3}.rs`,
`topology/topology.rs`, `topology/topology_traverser.rs`,
`topology/trimming_tools.rs`.
