# blend.fillet_executive_p3: overhanging box union, all edges at r=0.5, fails in Rust only (2026-10-26)

Unit: `blend.fillet_executive_p3` (end-to-end `SmApiFilletEdges`, real-solid runner).

Input (not committed as a fixture, since it would fail the gate):

```json
{"real_solid": {"source": {"kind": "box_union",
  "a": {"origin": [0, 0, 0], "size": [20, 10, 10]},
  "b": {"origin": [5, -2, 10], "size": [10, 6, 10]}},
  "radius": 0.5, "xsect": 1, "continuity": 1, "thumbweight": 1, "all_edges": true}}
```

Values: C++ returns `SM_SUCCESS` with 58 vertices, 114 edges and 58 faces, volume
1303.4715432511623. Rust returns 1001 (`FILEXEC_ERR`) and leaves the 18/28/12 body,
volume 1360. The same source with `all_edges` at r=1.0 agrees with C++, and so do
its 21 single-edge fillets at r=0.5, 1 and 2 (`stubs3-overhang-union-*` in
`fixtures/blend.fillet_executive_p3/stubs_round3.json`). Those single-edge cases
reach the Nx1 multi-neighbour end-face gap checks; they fail with 1001 in both
kernels.

The same thing happens, Rust 1001 against C++ success, with `all_edges` at r=0.5
when `b` is `{"origin": [5, -2, 10], "size": [10, 14, 10]}` (C++ volume
1496.645519676933, Rust 1560) and when it is `{"origin": [-2, -2, 10], "size": [10, 6, 10]}`
(C++ 1319.4710975188434, Rust 1376). All three variants agree at r=1.0.

Suspected cause: not the Nx1 gap checks. A debug trace showed this all-edges r=0.5
run never reaches them. The divergence is somewhere else in the all-edges corner
set (several corners of the overhang meet at r=0.5) and has not been narrowed down.

Found while probing for fixtures that reach the gap checks (slice
`slice/p1-stubs-round3`).
