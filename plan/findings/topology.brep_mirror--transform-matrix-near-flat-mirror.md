# topology.brep_mirror: matrix Transform of a near-flat revolution still differs after the ray-fire fix (2026-10-31)

Unit: `topology.brep_mirror` (`SmBrep::Transform(const double[4][4], SmBoolean)`, PR #1595).

Input: generated case `mat0076` (seed 17267955321482986155), `Solids/revolution.smb`,
matrix rows `(0.8058680190465212, -0.29026436029085756, 0, 0)`,
`(-0.10430573503780678, 0.141964763216248, 0.9843608686062126, 0)`,
`(5.828353536945563e-11, -7.932649432106016e-11, 1.7616378843798325e-11, 0)`,
`(43.871017221482376, 4.321056959395637, 5.973499327348293, 1)`. The transform takes
the `bMirrorV` and `bMirrorUt` branches and flattens the body to about 1e-10 along one
axis.

Fixed by #1600: the second `Mirror`'s `FindAndSetInfiniteRegion` returned `SM_ERR`
(C++ `SM_SUCCESS`). Its ray fire rebuilt the plain `SmSurfaceCache` fallback of the
degenerate face as a trimmed cache. A temporary C++ probe confirmed that, with editing
disabled as in `Mirror`, the C++ cache of that face is a plain `SmSurfaceCache`
(type 15100) and the others are `SmTrimSrfCache` (15120). Rust now matches that, and the
statuses agree.

Still open: 310 values differ. `geometry[120..135]` is one curve's five samples in
reverse order: C++ `(7.641457250313985, -2.3347845500435853, 3.9374434744248505)`,
Rust `(0, 0, 0)` at the first sample, and so on. `record[342..415]` are graph values
off by about 1e-6 relative (e.g. C++ `-1.0000014427543823e1`, Rust
`-1.0000016944616139e1`). In Rust the first ray of that search still fails
(`BSplineSurface::evaluate` SM_ERR, a zero-tangent repair inside the LSI leaf
`DropPointsToTangentPlane` on the degenerate face's plain cache). The second ray finds
the region. Whether the C++ first ray fails the same way is not known yet. The ray
fires create UV trim curves, so a different ray history is the suspected cause of the
reversed curve and the small graph differences.

Narrowed (2026-10-30): the sequence was replayed step by step, `Transform(V)`, `Mirror(V)`,
`Transform(Uᵀ)`, `Mirror(Uᵀ)`. Every curve, surface, UV trim and edge/face/vertex
tolerance is bit-identical through `Transform(Uᵀ)`. After `Mirror(Uᵀ)` only the UV trims
of edge uses on the flattened face differ (the `v` coordinate, about 1e-6 relative).
Skipping `FindAndSetInfiniteRegion` does not change them, so they come from `Mirror`'s
own `CreateUVTrimCurves`. For those edge uses all four `DropCurve`/`DropAndTrimCurve`
attempts fail on both sides alike (the standalone drops agree bit for bit in
`intersect.surface_drop_lib_p3`). The trims come from the 101-point brute-force drop in
`SmEdgeuse::CreateUVTrimCurve`.

Root: `SmSurface::GlobalPointSolve(domain, SM_SO_MINIMIZE, point, 1e-5, NULL,
SM_SR_ALL)` on that face's B-spline surface. Standalone `surface.surface_p4` `gps`
cases built from the brute-force inputs fail 88 of 202. For example, the first sample
point: C++ uv `(-6.217248937900877e-15, -2.6622142923926617e-05)`, Rust
`(-7.719268526242007e-15, -4.4668287143539916e-05)`. Both distances are about 1e-15:
the surface is flat to rounding along `v` there, so both answers have converged, and the
choice depends on the solver's iteration path. Normal evaluation on the same surface
agrees bit for bit. A difference in operation order inside the shared global point
solve (tree seeding, Newton or the stopping test) is therefore suspected. The generated
`gps` cases in `surface.surface_p4` only cover unit-scale surfaces (`gps_ok`), so this
regime is not exercised there. The reversed-curve values above were not re-checked
after the narrowing.
