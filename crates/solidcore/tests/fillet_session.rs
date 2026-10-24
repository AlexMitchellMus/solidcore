//! Actual-object lifetime controls for the production fillet session.

use solidcore::blend::fillet_corner::p3::FilletCornerKind;
use solidcore::blend::fillet_executive::p3::VertexEntry;
use solidcore::blend::fillet_executive::session::{FilletBody, FilletObject, FilletSession};
use solidcore::blend::fillet_solver::p1::FilletSurfaceGenerator;
use solidcore::geom::attribute::ids;
use solidcore::math::{Point3, PI};
use solidcore::primitive::api_primitives::create_box;
use solidcore::topology::edge;
use solidcore::{EntityRef, Kernel};

#[test]
fn actual_box_solver_stays_authoritative_through_corner_creation() {
    let mut k = Kernel::default();
    let mut target = create_box(&mut k, Point3::new(0.0, 0.0, 0.0), 10.0, 20.0, 30.0).unwrap();
    let edge = target.edges().ids().next().unwrap();
    let edgeuse = edge::blend_edgeuse(&k, &mut target, edge).unwrap().unwrap();
    let mut session = FilletSession::new(&mut k, &mut target).unwrap();
    let generator = session
        .add_surface_generator(FilletSurfaceGenerator::circular(false, 0.05, false, false))
        .unwrap();
    let solver = session
        .new_constant_radius_solver(&k, edgeuse, 0.5, 1e-4, 20.0 * PI / 180.0, 10.0 * PI / 180.0)
        .unwrap();
    session.set_solver_generator(solver, generator).unwrap();
    let actual_base = std::ptr::from_ref(session.solver(solver).unwrap().base());
    assert!(session
        .solver(solver)
        .unwrap()
        .base()
        .surface(0)
        .unwrap()
        .is_some());
    assert!(session
        .solver(solver)
        .unwrap()
        .base()
        .surface(1)
        .unwrap()
        .is_some());
    session.load_solver(&k, solver).unwrap();
    assert_eq!(
        session.executive().solvers,
        vec![FilletObject::Solver(solver)]
    );
    assert!(session
        .solver(solver)
        .unwrap()
        .base()
        .fillet_executive()
        .is_some());
    assert_eq!(session.executive().filleted_vertices.len(), 2);
    session.create_corners(&k, None).unwrap();
    assert_eq!(session.executive().corners.len(), 2);
    assert!(std::ptr::eq(
        actual_base,
        session.solver(solver).unwrap().base()
    ));
    for vertex in &session.executive().filleted_vertices {
        let Some(FilletObject::Topology(FilletBody::Target, EntityRef::Vertex(expected))) = vertex
        else {
            panic!("actual target vertex required");
        };
        let Some(VertexEntry::Corner(FilletObject::Corner(id))) =
            session.executive().vertex_corner_map.get(vertex)
        else {
            panic!("actual retained corner required");
        };
        assert_eq!(session.corner(*id).unwrap().vertex, Some(*expected));
        assert_eq!(session.corner(*id).unwrap().kind, FilletCornerKind::NByOne);
        session.target().vertex(*expected).unwrap();
    }
    let key = session.target_identity(EntityRef::Edge(edge));
    assert!(session
        .attributes()
        .find_attribute(key.key, ids::MASS_PROPERTIES)
        .is_none());
    assert_ne!(
        FilletObject::Topology(FilletBody::Target, EntityRef::Edge(edge)),
        FilletObject::Topology(FilletBody::Pseudo, EntityRef::Edge(edge))
    );
}
