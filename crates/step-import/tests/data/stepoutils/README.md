# StepOutils corpus

Small, self-authored STEP cases for the same parser and builder branches the
STEP Tools corpus exercises. Unlike `corpus/steptools/`, these files are ours
and may be committed.

| File | Branch it exercises |
|---|---|
| `valid_sphere.stp` | baseline closed analytic solid |
| `dangling_reference.stp` | a shape representation naming a missing entity |
| `loopless_face.stp` | a face with no bounds |
| `multi_winding_loop.stp` | an edge loop that winds the same edges twice |
| `loose_surface_face.stp` | an open shell / sheet rather than a solid |

These are intentionally tiny. Add a case here when a corpus failure reduces to
one of these structural shapes; keep vendor-specific geometry out of the
committed tree.
