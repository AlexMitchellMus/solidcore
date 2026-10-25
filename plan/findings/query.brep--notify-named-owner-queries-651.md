# query.brep: remaining default-identity query and owner-cleanup sites resolved (#651) (2026-10-25)

- **PR:** #962 (closes #651)
- **Context:** PRs #840 and #905 added `NamingScheme`, `Brep::geometry_owner_keys`, and `BrepQuery::of_brep(_with_id)`, updating `boolean::merge::{p1,p4,p5}` and `topology::brep::p1::{surfaces,curves}`. Nine remaining library sites still hardcoded `BrepIdentity::default()`, `BrepQuery::default()`, `arena_identity(...)`, or `BrepQuery::new(identity)` (dropping `geometry_owner_keys` overrides on mixed-scheme stitches/splits).
- **Fix:**
  - Generalized `BrepQuery::with_proven_owners(identity, b)` for any lifetime `'a`, exposed `BrepQuery::effective_identity`, and added `BrepClassifier::of_brep` / `of_brep_with_id`.
  - Updated `query::api_queries` (`Target::new`), `blend::hcr` (`find_and_set_infinite_region` and `brep_surfaces`), `heal::heal_data::p1` (`brep_find_and_set_infinite_region`), `heal::stitch::p2` (`find_and_set_infinite_region_with_identity`), `topology::brep::p1` (`create_parallel_projection_curves`), `topology::brep::p3` (`destroy_face`, `squeeze_edge`, `turn_to_nurbs_in`, and `brep_curve_solve`), and `topology::brep::p6` (`split_closed`) to resolve geometry owners via the B-rep's naming scheme and proven owner keys.
  - Added regression tests verifying that `create_box` (`NamingScheme::Notify`) builds full spatial trees and classifies interior points into its finite solid region.
