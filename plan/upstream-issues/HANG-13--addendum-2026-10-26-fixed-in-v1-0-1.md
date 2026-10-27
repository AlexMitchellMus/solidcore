Updates: HANG-13

#### HANG-13 (update, 2026-10-26). Fixed in v1.0.1

Fixed in v1.0.1 (3fd5b95): SmPrimitiveCreation::CreateLinearSweep (curve overload), the sweep vector is unitized and a negative distance becomes a positive one along the reversed direction before the analytic domain SmExtent2d(0, 0, 1, dSweepDist) is built, so the full-circle cylinder branch no longer gets a reversed V domain. Verification: fixtures/sweep.api_brep/hang13-negative-full-circle.json and the four primitive.primitive_creation_p1 linear_circle_cap*_distance-1 cases, now compared with the FFI and shim preflights removed, pending the pin101 difftest run.
