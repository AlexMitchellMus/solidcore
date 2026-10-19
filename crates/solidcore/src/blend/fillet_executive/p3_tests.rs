//! Control-flow tests for executive p3. Dependency responses are deliberately
//! scripted here; these tests do not establish C++ geometric parity.
#![allow(clippy::unwrap_used, clippy::indexing_slicing)]

use super::*;
use crate::geom::solution_array::ObjectKind;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct Handle(u8);
impl ClassObject for Handle {
    fn class_kind(&self) -> ClassObjectKind {
        ClassObjectKind::Other
    }
}
impl SolutionObject for Handle {
    fn kind(&self) -> ObjectKind {
        ObjectKind::Other
    }
}

#[derive(Default)]
struct Access {
    process_error: Option<SolidError>,
    piecewise_error: Option<SolidError>,
    stitch_error: bool,
    editing: bool,
    fillet_editing: bool,
    events: Vec<&'static str>,
    stitch_tolerances: Vec<f64>,
    stitch_lamina: u64,
    corner_metadata: Option<CornerView<Handle>>,
    corner_point: Point3,
    constant_results: Vec<bool>,
    corner_fail: Option<&'static str>,
}
fn absent<T>() -> Result<T> {
    Err(invalid("test handle has no geometry"))
}
impl ExecutiveAccess<Handle> for Access {
    fn make_corner_topology(&mut self, _: &mut Kernel, _: Handle) -> Result<()> {
        self.events.push("corner_topology");
        Ok(())
    }
    fn corner(&self, _: Handle) -> Result<CornerView<Handle>> {
        self.corner_metadata
            .clone()
            .ok_or_else(|| invalid("no corner metadata"))
    }
    fn solver(&self, _: Handle) -> Result<SolverView<Handle>> {
        Ok(SolverView {
            geoms: vec![],
            edge: Handle(5),
            approx_tolerance: 0.125,
        })
    }
    fn geom_face(&self, _: Handle) -> Result<Option<Handle>> {
        absent()
    }
    fn topo_edge(&self, _: Handle) -> Result<TopoEdgeView<Handle>> {
        absent()
    }
    fn vertex_point(&self, _: Handle) -> Result<Point3> {
        Ok(self.corner_point)
    }
    fn edge_curve_ends(&self, _: &Kernel, _: Handle) -> Result<(Point3, Point3)> {
        Ok((Point3::new(0.0, 0.0, 0.0), Point3::new(10.0, 0.0, 0.0)))
    }
    fn set_calc_as_const_rad(&mut self, _: Handle, value: bool) -> Result<()> {
        self.constant_results.push(value);
        Ok(())
    }
    fn new_brep(&mut self, _: &Kernel) -> Result<Handle> {
        self.events.push("new");
        Ok(Handle(1))
    }
    fn delete_brep(&mut self, _: Handle) -> Result<()> {
        self.events.push("delete");
        Ok(())
    }
    fn tolerance(&self, _: &Kernel, _: Handle, _: bool) -> Result<f64> {
        Ok(0.125)
    }
    fn editing_enabled(&self, _: Handle) -> Result<bool> {
        Ok(self.editing)
    }
    fn set_editing_enabled(&mut self, brep: Handle, value: bool) -> Result<()> {
        if brep == Handle(0) {
            self.editing = value;
        } else {
            self.fillet_editing = value;
        }
        Ok(())
    }
    fn edges(&self, _: Handle) -> Result<Vec<Handle>> {
        Ok(vec![])
    }
    fn vertices(&self, _: Handle) -> Result<Vec<Handle>> {
        Ok(vec![])
    }
    fn faces(&self, _: Handle) -> Result<Vec<Handle>> {
        Ok(vec![Handle(2)])
    }
    fn add_attribute(&mut self, _: Handle, _: Option<Handle>) -> Result<()> {
        absent()
    }
    fn new_keep_attribute(&mut self, _: &Kernel) -> Result<Handle> {
        absent()
    }
    fn attribute_users(&self, _: Handle) -> Result<Vec<Handle>> {
        absent()
    }
    fn delete_attribute(&mut self, _: Handle) -> Result<()> {
        absent()
    }
    fn with_marks<T>(
        &mut self,
        k: &mut Kernel,
        _: Handle,
        _: Handle,
        _: &[(Handle, Handle)],
        body: impl FnOnce(&mut Self, &mut Kernel, &[ContextMark]) -> Result<T>,
    ) -> Result<T> {
        self.events.push("marks");
        body(self, k, &[])
    }
    fn shrink_face(&mut self, _: &mut Kernel, _: Handle) -> Result<()> {
        absent()
    }
    fn notify_face(&mut self, _: &mut Kernel, _: Handle, _: Handle) -> Result<()> {
        self.events.push("notify");
        Ok(())
    }
    fn boolean_post_process(
        &mut self,
        _: &mut Kernel,
        _: &mut FilletExecutive<Handle>,
        all: bool,
    ) -> Result<()> {
        assert!(all);
        self.events.push("post");
        Ok(())
    }
    fn process_smooth_edge_sequences(
        &mut self,
        k: &mut Kernel,
        _: &mut FilletExecutive<Handle>,
    ) -> Result<()> {
        self.events.push("smooth");
        // DoFilleting sets this flag before invoking CreateFilletBrep.
        let _ = k.context.doing_boolean();
        self.process_error.clone().map_or(Ok(()), Err)
    }
    fn set_brep_tolerance(
        &mut self,
        _: &mut Kernel,
        _: Handle,
        tol: f64,
        larger: bool,
        cascade: bool,
    ) -> Result<()> {
        assert_eq!(tol, 0.125);
        assert!(larger && cascade);
        self.events.push("tolerance");
        Ok(())
    }
    fn stitch_faces(
        &mut self,
        _: &mut Kernel,
        _: Handle,
        _: f64,
        faces: Option<&[Handle]>,
        spine: bool,
    ) -> Result<StitchAttempt> {
        assert!(faces.is_none());
        assert!(!spine);
        self.events.push("stitch_faces");
        Ok(StitchAttempt {
            status: None,
            stitched: 0,
            lamina: 0,
            max_vertex_gap: 0.0,
            max_edge_gap: 0.0,
        })
    }
    fn piecewise_merge(&mut self, k: &mut Kernel, _: &mut FilletExecutive<Handle>) -> Result<()> {
        assert!(k.context.doing_boolean());
        assert!(self.editing);
        self.events.push("piecewise");
        self.piecewise_error.clone().map_or(Ok(()), Err)
    }
    fn set_fillet_brep_edge_map(
        &mut self,
        _: &mut Kernel,
        _: &mut FilletExecutive<Handle>,
    ) -> Result<()> {
        self.events.push("edge_map");
        Ok(())
    }
    fn delete_faces(
        &mut self,
        _: &mut Kernel,
        _: &mut TopologyIntersector<Handle>,
        _: &[Handle],
        keep: Option<&[Handle]>,
    ) -> Result<()> {
        assert!(keep.is_none());
        self.events.push("delete_faces");
        Ok(())
    }
    fn merge_brep(
        &mut self,
        _: &mut Kernel,
        _: Handle,
        _: Handle,
        copy: Option<&mut CopyBrepMap<Handle>>,
        reverse: Option<&mut CopyBrepMap<Handle>>,
        analytics: bool,
    ) -> Result<()> {
        assert!(copy.is_none() && reverse.is_none() && analytics);
        self.events.push("merge");
        Ok(())
    }
    fn do_stitching(
        &mut self,
        _: &mut Kernel,
        _: Handle,
        tolerance: f64,
        vertices: Option<&[Handle]>,
        edges: Option<&[Handle]>,
        settings: StitchSettings,
        min_vertex: Option<&mut f64>,
        min_edge: Option<&mut f64>,
    ) -> Result<StitchAttempt> {
        assert!(
            vertices.is_none() && edges.is_none() && min_vertex.is_none() && min_edge.is_none()
        );
        assert!(settings.do_region_nesting && settings.ignore_problems && settings.validate_result);
        assert!(settings.split_edges_with_vertices && settings.making_manifold_solid);
        assert!(
            !settings.fast_edge_compare
                && !settings.squeeze_small_edges
                && !settings.remove_laminar_slivers
        );
        self.stitch_tolerances.push(tolerance);
        Ok(StitchAttempt {
            status: self.stitch_error.then(failure),
            stitched: 0,
            lamina: self.stitch_lamina,
            max_vertex_gap: 0.0,
            max_edge_gap: 0.0,
        })
    }
    fn remove_topo_edges(&mut self, _: &mut Kernel, _: &mut FilletExecutive<Handle>) -> Result<()> {
        self.events.push("remove_topo");
        Ok(())
    }
    fn orient_trimmed_surfaces(
        &mut self,
        _: &mut Kernel,
        _: Handle,
        outward: bool,
        solid: &mut bool,
        unify: bool,
        faces: Option<&mut Vec<Handle>>,
    ) -> Result<()> {
        assert!(outward && !*solid && !unify && faces.is_none());
        self.events.push("orient");
        Err(failure()) // upstream ignores this ordinary status
    }
    fn calc_corner_vertices(&mut self, _: &mut Kernel, _: Handle) -> Result<()> {
        self.events.push("corner_vertices");
        if self.corner_fail == Some("vertices") {
            Err(failure())
        } else {
            Ok(())
        }
    }
    fn check_degeneracy(&mut self, _: &mut Kernel, _: Handle) -> Result<bool> {
        self.events.push("degeneracy");
        Ok(false)
    }
    fn adjust_topology(&mut self, _: &mut Kernel, _: Handle) -> Result<()> {
        self.events.push("adjust");
        Ok(())
    }
    fn adjust_tangent_corner(&mut self, _: &mut Kernel, _: Handle) -> Result<()> {
        self.events.push("adjust_tangent");
        Ok(())
    }
    fn calc_corner_edges(&mut self, _: &mut Kernel, _: Handle) -> Result<()> {
        self.events.push("corner_edges");
        Ok(())
    }
    fn calc_corner_geom(
        &mut self,
        _: &mut Kernel,
        _: Handle,
        _: &mut FilletExecutive<Handle>,
    ) -> Result<()> {
        self.events.push("corner_geom");
        Ok(())
    }
}

fn state() -> FilletExecutive<Handle> {
    FilletExecutive {
        target: Some(Handle(0)),
        fillet: None,
        corners: vec![],
        solvers: vec![],
        topo_edges: vec![],
        edge_attribute: None,
        initialized: true,
        options: ExecutionOptions {
            piecewise_merge: false,
            topology_insertion: false,
            global_merge: false,
            trimming: false,
            classification: false,
        },
        error: ErrorInfo::default(),
        intersector: TopologyIntersector::new(None, None, 1e-5, 0.5),
        intersection_edges: vec![],
        intersection_vertices: vec![],
    }
}

#[test]
fn uninitialized_changes_only_generic_error_members() {
    let mut executive = state();
    executive.initialized = false;
    executive.error.status = 9;
    let mut access = Access::default();
    let mut k = Kernel::default();
    assert!(matches!(
        executive.do_fillet(&mut k, &mut access),
        Err(SolidError::InvalidInput { .. })
    ));
    assert_eq!(executive.error.code, 1);
    assert_eq!(executive.error.status, 9);
    assert_eq!(
        executive.error.comment,
        "ERROR: fillet parameters not set up."
    );
    assert!(access.events.is_empty());
    assert!(!executive.intersector.doing_fillets);
}

#[test]
fn null_create_status_is_remapped_by_do_filleting_and_boolean_restores() {
    let mut executive = state();
    executive.target = None;
    let mut access = Access::default();
    let mut k = Kernel::default();
    assert!(matches!(
        executive.create_fillet_brep(&mut k, &mut access, None),
        Err(SolidError::Upstream { code: 1006, .. })
    ));
    for prior in [false, true] {
        k.context.set_doing_boolean(prior);
        assert!(matches!(
            executive.do_filleting(&mut k, &mut access, None),
            Err(SolidError::Upstream { code: 1001, .. })
        ));
        assert_eq!(k.context.doing_boolean(), prior);
    }
    assert!(executive.intersector.doing_fillets);
}

#[test]
fn ordinary_smooth_error_is_ignored_and_tail_and_selective_assignment_run() {
    let mut executive = state();
    executive.error.code = 6;
    let mut output = ErrorInfo {
        code: 5,
        comment: "keep comment".into(),
        side_face: Some(Handle(7)),
        ..ErrorInfo::default()
    };
    let mut access = Access {
        process_error: Some(failure()),
        ..Access::default()
    };
    executive
        .create_fillet_brep(&mut Kernel::default(), &mut access, Some(&mut output))
        .unwrap();
    assert_eq!(
        access.events,
        ["smooth", "new", "tolerance", "stitch_faces"]
    );
    assert!(access.fillet_editing);
    assert_eq!(executive.fillet, Some(Handle(1)));
    assert_eq!(output.code, 6);
    assert_eq!(output.comment, "keep comment");
    assert_eq!(output.side_face, Some(Handle(7)));
}

#[test]
fn deferral_is_never_swallowed_or_remapped_and_restores_boolean() {
    let expected = SolidError::unported(
        "SmFilletExecutive::ProcessSmoothEdgeSequences",
        "blend.fillet_executive_p2",
    );
    let mut access = Access {
        process_error: Some(expected.clone()),
        ..Access::default()
    };
    let mut executive = state();
    let mut k = Kernel::default();
    assert_eq!(
        executive.do_filleting(&mut k, &mut access, None),
        Err(expected)
    );
    assert_eq!(access.events, ["smooth"]);
    assert!(!k.context.doing_boolean());
    assert!(executive.intersector.doing_fillets);
    assert!(executive.fillet.is_none());
}

#[test]
fn piecewise_return_restores_both_temporary_flags_on_success_and_error() {
    for error in [
        None,
        Some(failure()),
        Some(SolidError::unported(
            "SmFilletExecutive::DoPiecewiseMerge",
            "blend.fillet_executive_p1",
        )),
    ] {
        let mut executive = state();
        executive.options.piecewise_merge = true;
        let mut access = Access {
            piecewise_error: error.clone(),
            ..Access::default()
        };
        let mut k = Kernel::default();
        assert_eq!(
            executive.do_filleting(&mut k, &mut access, None),
            error.map_or(Ok(()), Err)
        );
        assert!(!access.editing);
        assert!(!k.context.doing_boolean());
        assert_eq!(access.events.last(), Some(&"piecewise"));
        assert!(!access.events.contains(&"edge_map"));
    }
}

#[test]
fn stitching_retries_exactly_five_times_doubles_tolerance_and_ignores_orientation_error() {
    let mut executive = state();
    executive.options.trimming = true;
    let mut access = Access {
        stitch_lamina: 1,
        ..Access::default()
    };
    let mut k = Kernel::default();
    executive.do_filleting(&mut k, &mut access, None).unwrap();
    assert_eq!(access.stitch_tolerances, [0.125, 0.25, 0.5, 1.0, 2.0]);
    assert_eq!(access.events.last(), Some(&"orient"));
    assert!(access.events.contains(&"post"));
    assert!(executive.fillet.is_none());
    assert!(executive.intersector.other_brep().is_none());
    assert!(!access.editing);
}

#[test]
fn final_stitch_failure_sets_bad_merge_and_stops_before_cleanup_postprocess() {
    let mut executive = state();
    executive.options.trimming = true;
    let mut access = Access {
        stitch_error: true,
        ..Access::default()
    };
    let mut k = Kernel::default();
    assert!(matches!(
        executive.do_filleting(&mut k, &mut access, None),
        Err(SolidError::Upstream { code: 1001, .. })
    ));
    assert_eq!(access.stitch_tolerances.len(), 5);
    assert_eq!(executive.error.code, 7);
    assert!(!access.events.contains(&"post"));
    assert!(!access.editing);
    assert!(!k.context.doing_boolean());
}

fn corner() -> CornerView<Handle> {
    CornerView {
        blending: true,
        tangent: true,
        solvers: vec![Handle(3), Handle(4)],
        vertex: Handle(6),
        approx_tolerance: 0.125,
        corner_type: 9,
        name: "corner name".into(),
    }
}

#[test]
fn tangent_corner_constant_radius_uses_both_distances_and_twice_approx_tolerance() {
    for (point, expected) in [
        (0.0, true),
        (0.25, true),
        (0.2500000001, false),
        (5.0, false),
        (9.75, true),
        (10.0, true),
    ] {
        let mut executive = state();
        executive.corners.push(Handle(7));
        let mut access = Access {
            corner_metadata: Some(corner()),
            corner_point: Point3::new(point, 0.0, 0.0),
            ..Access::default()
        };
        executive
            .create_fillet_brep(&mut Kernel::default(), &mut access, None)
            .unwrap();
        assert_eq!(access.constant_results, [expected]);
        assert_eq!(
            access.events,
            [
                "smooth",
                "corner_topology",
                "corner_vertices",
                "degeneracy",
                "adjust",
                "adjust_tangent",
                "corner_vertices",
                "corner_edges",
                "new",
                "tolerance",
                "corner_geom",
                "stitch_faces"
            ]
        );
    }
}

#[test]
fn corner_failure_records_context_and_stops_without_allocating_shell() {
    let mut executive = state();
    executive.corners.push(Handle(7));
    let mut access = Access {
        corner_metadata: Some(corner()),
        corner_fail: Some("vertices"),
        ..Access::default()
    };
    assert!(matches!(
        executive.create_fillet_brep(&mut Kernel::default(), &mut access, None),
        Err(SolidError::Upstream { code: 1001, .. })
    ));
    assert_eq!(executive.error.code, 3);
    assert_eq!(executive.error.status, 2);
    assert_eq!(executive.error.corner_type, 9);
    assert_eq!(executive.error.vertex, Some(Handle(6)));
    assert_eq!(executive.error.surfaces, [None]);
    assert_eq!(executive.error.entity_name, "corner name");
    assert_eq!(
        executive.error.comment,
        "Fillet Corner Error: Unable to Compute Geometry"
    );
    assert!(!access.events.contains(&"new"));
}
