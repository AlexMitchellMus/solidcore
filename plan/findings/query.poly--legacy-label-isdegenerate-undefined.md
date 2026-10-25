# query.poly: IsDegenerate on a labelled fake SmFace is undefined, rejected by the shim (2026-10-25)

Unit: query.poly (also notes on intersect.curve_class_p1).

Inputs: `gen0378`, `gen0496` (22 more at `--seeds 6000`): `degenerate` op on an edge whose
poly face has an original-face label (not a real `SmFace`) while `OKBackPtrs` holds.
Upstream `SmPoly.cpp:21281-21283` then calls `pFace->GetSurface()` on the label: undefined
in C++. They were deferred (`topology.face_p2`, `SmFace::GetSurface`), but no port of
`GetSurface` can make a fake face comparable.

Handling: the shim decides from the built model (face label + back pointers), before any
call, and returns `STATUS_UNDEFINED_INPUT` (-1101); `run_cpp` maps that status to
`Output::rejected`. No C++ result or sentinel is inspected. The same applies to
`owner_degenerate` edges whose label no explicit owner resolves. Granularity is the whole
case (the Rust side deferred the whole case too). Real original faces stay compared by the
`owner_degenerate` cases.

intersect.curve_class_p1 (face point class guard): with the shim's zero-size-UV-domain guard
removed, a face whose UV domain size is (0,1), (1,0) or (0,0) segfaults in
`sm_ComputeFacePointClass`; so rejecting either-axis zero is exact (measured), not
over-broad. Exact-zero input domains are normalised by face construction; sizes down to 1e-100
on the other axis run. Cause not located.
