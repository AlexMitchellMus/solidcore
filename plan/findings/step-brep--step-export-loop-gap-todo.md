# TODO: exact STEP export produces a B-rep its reader rejects

- **Area:** `step-brep` / `step-export` / `step-import` round-trip
- **Reported from:** SolidNode gear graph (sketch → loft → Boolean cut → exporter)
- **Observed:** `write_brep_exact` succeeds, but importing the written STEP fails:

```text
solid #74546: topology::brep::p7: invalid input: loop gap larger than tolerance
```

## Why this is not closed

Healing the evaluated B-rep in `brep_to_model` makes this particular file
importable, but it does not prove the writer preserved the input. It can change
the exported topology/geometry and is not an acceptable final answer.

## Investigation required

1. Reproduce with a minimized exact B-rep (two-section loft, then a Boolean cut).
2. Before writing, compare for every loop:
   - live edge start/end vertices against curve endpoints at the edge interval
   - the reported loop-gap distance against `zone_tol`
3. Determine whether the gap already exists in the evaluated B-rep or is
   introduced by `brep_to_model`/`step-export`.
4. Compare face/edge/vertex counts, bounds and volume before and after any repair.

## Possible fixes

- If the kernel result is open: heal/close where the Boolean or skin result is
  produced, not only in export.
- If the writer loses closure: preserve/emit the exact loop endpoints and
  tolerances in `brep_to_model`/`step-export`, or reject the export rather than
  write a file the reader rejects.
- Reader-side tolerance changes are only acceptable with an upstream-issue
  record and a minimized reproducer.

## Acceptance criteria

- Export → import returns the same topology and geometry (within documented
  tolerances) for the gear reproducer.
- No export-side healing is required for this case.
- Regression test covers the minimal B-rep and the full-file round-trip.
