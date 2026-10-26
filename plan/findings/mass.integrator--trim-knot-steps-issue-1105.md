# mass.integrator: trim knot spans and Windows sphere precision (2026-10-25)

Units: `mass.integrator`, `topology.face_p2`. Issue #1105; PR #1106.

At base `caf670b7f214c7a51d0c77d451c170c7f9b11080`, the existing
`precise_sphere_properties_scale_and_translate` test fails on Windows GNU
Rust 1.98.0 and passes on Linux x86-64. The relevant kernel is unchanged from
the independently reproduced `27b234a99d646d010e50e0dd4ad9d6e1a88bb87d` and
`d674288e9f24758bd8a6fd45c98b64f1d5094ae1` baselines.

Input: pole-free sphere, radius 5, center `(1,-2,0.5)`, construction tolerance
`1e-9`, precise-property accuracy `1e-10`, computation origin `bbox.mid()`,
default mass behavior. Analytic area is `314.159265358979326` and volume is
`523.598775598298857`; the existing relative limit is `1e-9`.

| Engine / mode | Area | Volume |
| --- | --- | --- |
| Windows GNU, upstream default | 314.159263637494462 | 523.598774064344411 |
| Linux Rust, upstream default | 314.159265344620053 | 523.59877561869075 |
| Linux pinned native C++, upstream default | 314.159265344620053 | 523.59877561869075 |
| Windows GNU, trim-knot opt-in | 314.159265365340616 | 523.598775571151577 |

Almost all of the Windows area deficit comes from FaceId(3), whose default
area is `31.415924823658806`, versus approximately `31.4159265359`. The
underlying surface's sampled radial residual is below `1.78e-15`; sampled
stored trim/edge gaps are below `9.89e-10`, within construction tolerance.
The affected face's control data agrees bit for bit between hosts. These
measurements exclude a large surface or trim displacement as the explanation.

The cap's oriented cubic UV trims each have 297 unique knots. Upstream
`SmFace::ComputePreciseProperties` integrates each entire trim with adaptive
Cash-Karp steps; it retries a smaller initial step only on an integration
error. A successful local error estimate need not detect a feature between
sample points. Integrating individual trim knot spans at the same requested
accuracy removes this deficit on Windows. A focused planar regression also
demonstrates a narrow feature missed by whole-curve integration. Platform
rounding changes the adaptive schedule; the exact first divergent step has
not been isolated. Native Windows C++ confirmation remains unavailable.

`Kernel::fixes.precise_properties_at_trim_knots` enables the remedy. The
outer absolute error-scale floor is shared among spans; inner iso integration
and the caller's accuracy remain unchanged. Default mode preserves the
pinned upstream step and accumulation order. No expected value, tolerance,
skip, mismatch entry, or unit status changes.

The original default-mode Windows precision test still fails with this opt-in
field off. A passing opt-in regression or Linux run does not resolve that
default-mode gate. Keep #1105 open and PR #1106 draft until required validation
and the treatment of upstream defaults are settled. The independent partial
primitive coverage change belongs to #1070.
