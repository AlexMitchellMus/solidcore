# TODO: missing import/export formats

- **Area:** format support / import-export
- **Current:** `.solidnode` load/save; STEP `.step`/`.stp` load and export
  (exact B-rep, faceted, supported assemblies).

## STEP representation gaps (do first)

These are already STEP files; they need reader/repair work rather than a new
container parser.

1. **Geometric-set / surface / sheet / wireframe models** — display-only import
   with an explicit non-exact warning. No exact/B-rep claim.
2. **Assemblies with open, no-loop, or loop-gap bodies** — skip or repair bad
   bodies per body instead of rejecting the whole assembly; expose the measured
   gap in diagnostics.
3. **AP242 tessellated / PMI / product-manufacturing semantics** — at least
   classify and report them precisely; full semantics later.
4. **AP214 / AP242 export** — the writer currently emits an AP203-style schema.

## File formats

| Format | Load | Save | Priority | Notes |
| --- | --- | --- | --- | --- |
| STL `.stl` | Missing | Missing | P1 | Mesh handoff / 3D printing |
| OBJ `.obj` | Missing | Missing | P1 | Mesh interchange and debugging |
| 3MF `.3mf` | Missing | Missing | P2 | Slicer / manufacturing |
| glTF / GLB | Missing | Missing | P2 | Web / AR viewer export |
| PLY `.ply` | Missing | Missing | P3 | Scan / mesh data |
| IGES `.igs` / `.iges` | Missing | Low priority | P3 | Legacy CAD import |
| USD / USDZ | Missing | Maybe | P3 | Visualization; not a CAD authority |
| Parasolid `.x_t` / `.x_b` | Missing | Missing | P4 | Read-only `parasolid-core` covers limited profiles; no open-source writer |
| ACIS `.sat` / `.sab` | Missing | Missing | P4 | Proprietary kernel |

## Parasolid note

`parasolid-core` (MIT/Apache-2.0) parses `X_T`/`X_B` for five verified schema
profiles into a Parasolid-native B-rep. It is read-only: `write_xb` returns the
retained original bytes, not an encoder. A SolidCore adapter would be required,
and supported coverage is limited. Saving Parasolid needs the licensed Siemens
SDK or a commercial translator.

## Acceptance criteria

- Every format has a manifest entry with a declared outcome: exact, display-only,
  assembly, or unsupported.
- Display-only imports are visibly distinct from exact solids and never feed
  exact Boolean/mass operations.
- No format claims success merely because a mesh is non-empty; bounds/volume
  checks apply where an independent reference exists.
- P0 STEP gaps are covered by regression fixtures before new container parsers
  are added.
