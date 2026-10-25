# Findings

One file per finding: Rust-vs-C++ disagreements, upstream bugs, nondeterminism and other
facts later agents need. Every PR adds its own new file, so parallel PRs never conflict.
Never append to a shared file. Findings from before 2026-10-25 are in the frozen archive
[`plan/findings.md`](../findings.md).

- **Name:** `<unit id or general>--<short-slug>.md`, lowercase, e.g.
  `topology.brep_p3--turn-to-nurbs-uv-trim-notify.md`. Pick a slug nobody else will.
- **Content:** start with `# <unit>: <title> (<YYYY-MM-DD>)`, then the unit, inputs, both
  values, suspected cause and the PR. Later updates go in the same file, or a new file
  that links to it.
- **References** from code, `known_mismatches.toml` or unit notes give the file path,
  e.g. `plan/findings/topology.brep_p3--turn-to-nurbs-uv-trim-notify.md`.
