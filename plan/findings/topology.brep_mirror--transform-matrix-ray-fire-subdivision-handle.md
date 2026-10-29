# topology.brep_mirror: matrix Transform's Mirror fails its infinite-region ray fire on a near-flat body (2026-10-30)

Unit: `topology.brep_mirror` (`SmBrep::Transform(const double[4][4], SmBoolean)`, PR #1595).

Input: generated case `mat0076` (seed 17267955321482986155), `Solids/revolution.smb`,
matrix rows `(0.8058680190465212, -0.29026436029085756, 0, 0)`,
`(-0.10430573503780678, 0.141964763216248, 0.9843608686062126, 0)`,
`(5.828353536945563e-11, -7.932649432106016e-11, 1.7616378843798325e-11, 0)`,
`(43.871017221482376, 4.321056959395637, 5.973499327348293, 1)`, `bScaleTol = TRUE`
(a rotation, a shear and a near-zero third singular value: the body is flattened to
about 1e-10 thickness). The same matrix without the translation fails the same way.

Values: C++ status `0` (SM_SUCCESS); Rust status `1001` (SM_ERR), and the edge-curve
samples from `geometry[120]` on differ (e.g. C++ `51.51247447179636`, Rust
`43.871017221482376`), as do a few graph values.

Both sides decompose the matrix the same way (`SVDecompose` is bit-identical in
`geom.matrix`), and both take the `bMirrorV` and `bMirrorUt` branches. The divergence is
in the second `Mirror` (about `Uᵀ`, after the `Uᵀ` transform scaled by
`(5.6e-11, 0.827, 1.025)`). There `FindAndSetInfiniteRegion` (multi-region B-rep) fires
its rays. The first ray, after `Mirror(V)`, succeeds. All three rays of the second search
fail in Rust with `InvalidInput { surface::surface_cache::p3::SurfaceCache::subdivide_node,
"cache cannot name new subdivision surfaces" }`. The face line solve's cache
(`topology::face::p4::FaceLineCache::build`) maps every handle to the one solve surface,
and it sets no subdivision-handle factory (`SurfaceCache::set_subdivision_handles`). So
when the line solve on the flattened face has to subdivide the tree further, it cannot
name the halves, where upstream allocates new `SmSurface` objects. The search then ends
with `SM_ERR` (`topology::brep::p9::find_infinite_region`).

Suspected cause: a port-side gap in the B-rep ray-fire path (face line cache handle
model; related to #1032), not in the matrix transform. Not fixed in #1595.
