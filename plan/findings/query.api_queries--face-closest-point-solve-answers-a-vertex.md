# query.api_queries: the face closest-point solve over `LiveBrepObjects` answers a vertex, not the face (2026-10-30)

- **Unit:** `query.api_queries` (`face_closest_point_in_shape`, the port of `SmApiFaceClosestPoint`), reached through the topology solver (`intersect.topology_solver`, `TopologySolver::shape_point_solve`) over a `query::brep::LiveBrepObjects` view.
- **What was run:** `solidcore-smlib-api`'s pick helper (`crates/solidcore-smlib-api/src/pick.rs`, the port of the SolidCraft fork's `smlib_pick::face_distance`). It was first written as:
  1. `Shape::new(k, &(BrepId 0, &b), &[], &[(BrepId 0, f)], &[], &[])`;
  2. `LiveBrepObjects::new(&mut b, &BrepIdentity { brep_id: 0, keys, owner_of: |_| None }, None)` (no B-rep cache);
  3. `face_closest_point_in_shape(k, &shape, &objects, (BrepId 0, f), p)`.

  The keys were first the arena identity. They were then each face's and edge's own geometry owner, which is what `Brep::owner_identity` gives for them. The owners resolved (surface owners type 16002, keys 17179869185 to 17179869190).
- **Inputs:** `create_box` at the origin, 10 × 20 × 30, run under `UpstreamFixes::solidcraft()`, at the centre of three of its faces.
- **Values (Rust, from the leader's diagnostic run of `tests/diag.rs`, 2026-10-30):**
  - At p = (5, 10, 30), the centre of the top face, every face answers a corner. Face 0 gives 32.0156 at (10, 0, 0). Faces 1, 3 and 5 give 11.1803 at (10, 20, 30), face 2 gives 11.1803 at (0, 20, 30), and face 4 gives 11.1803 at (10, 0, 30). The expected values are 0 for the top face and 30 for the bottom.
  - At p = (0, 10, 15), the centre of the x = 0 face, the answers are 18.0278 at box corners (and 20.6155 for face 3). The expected value is 0 for the x = 0 face.
  - At p = (5, 10, 0), the centre of the bottom face, the answers are 11.1803 or 32.0156 at corners. The expected value is 0.
  - The minimum over the face is never found. The solve answers one of the face's vertices.
- **C++ value:** not run on this input. `SmApiFaceClosestPoint` on a box face at its own centre is 0, which is the SolidCraft fork's pick rule: the C API's `capi/tests` pick faces by their centres.
- **Coverage:** `crates/difftest/src/units/query_api_queries.rs` case 19 compares `face_closest_point_in_shape` with C++. It goes through difftest's own adapter (`Objects`), not a bare `LiveBrepObjects`. That adapter builds a `BrepCache` with `build_trees` (`brep_cache`), keeps a snapshot of each B-rep for `face_entity`, and forwards the rest to `LiveBrepObjects`. The path a library caller takes, a bare `LiveBrepObjects` with no cache and no `face_entity`, is not covered by any difftest.
- **Suspected cause:** the bare `LiveBrepObjects` view differs from difftest's adapter. It has no `brep_cache` (`TopologyObjects::brep_cache` errors), and `face_entity` is the trait's default `None`. So the face's interior candidate is never accepted (a classification or local face solve that errors, its error dropped). Only the shape's vertex candidates remain, and the first is returned (`SolutionRequested::Single`). Not confirmed: needs a difftest case of this exact call on a box face centre against pinned C++, then a trace of the face's local solve.
- **Workaround in the crate:** `smlib-api`'s `face_distance` is the surface's closest point (`surface_closest_point`), counted when `face_classify_uv` puts it on the face, together with the face's edges' closest points (`edge_closest_point`). The nearest wins. That is the same minimum for a trimmed face. The pick helper uses it, and its unit test asserts distance 0 at the box's face centres.
- **Follow-up (queued, not taken):** a difftest case of `face_closest_point_in_shape` over a bare `LiveBrepObjects` on `create_box` face centres against pinned C++, then a fix of the adapter or the solve. Once that matches, `pick::face_distance` can go back to the solve.
