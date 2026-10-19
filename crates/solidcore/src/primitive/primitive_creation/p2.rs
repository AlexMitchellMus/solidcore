//! Partial analytic primitives, tapered extrusions, curve and pipe sweeps,
//! planar-path sweeps, skins and face blends.
//!
//! Upstream: `source/SMLib/src/SmPrimitiveCreation.cpp`,
//! `source/SMLib/inc/SmPrimitiveCreation.h`.
//!
//! | Rust | Upstream |
//! |------|----------|
//! | [`PrimitiveCreation::create_partial_cone`] | six-argument `CreateCone` |
//! | [`PrimitiveCreation::create_partial_sphere`] | four-argument `CreateSphere` |
//! | [`PrimitiveCreation::create_partial_torus`] | five-argument `CreateTorus` |
//! | [`PrimitiveCreation::create_taper_extrude`] | `CreateTaperExtrude` |
//! | [`PrimitiveCreation::create_curve_sweep_from_edges`] | `CreateCurveSweepFromEdges` |
//! | [`PrimitiveCreation::create_curve_sweep_from_faces`] | `CreateCurveSweepFromFaces` |
//! | [`PrimitiveCreation::create_pipe_sweep`] | `CreatePipeSweep` |
//! | [`PrimitiveCreation::create_sweep_along_planar_path`] | eight-argument static `CreateSweepAlongPlanarPath` |
//! | `sweep_profiles_along_path` (private) | `SweepProfilesAlongPath` |
//! | [`PrimitiveCreation::create_skin_primitive`] | nine-argument `CreateSkinPrimitive` |
//! | `skin_primitive_internal` (private) | `sm_CreateSkinPrimitiveInternal` |
//! | [`PrimitiveCreation::create_blend_primitive_between_faces`] | face/edge `CreateBlendPrimitive` |
//! | [`PrimitiveCreationAccess`] | B-rep services these entry points call |
//!
//! Deviations: typed ids and owned curve slots replace pointers; the caller's
//! [`PrimitiveCreationAccess`] resolves B-reps, contexts and services, and its
//! defaults return `Unported` for units that are not ported yet. Every caller
//! tail after such a call is ported. Upstream null dereferences (a missing edge
//! curve, a failed copy) become errors. Ignored statuses keep their outputs and
//! pass only `Unported` on. Debug drawing (`SM_DEBUG_CODE`) is not ported. No
//! recursion; loops are bounded by their input lists.

mod access;
mod blend;
mod curve_sweeps;
mod planar_path;
mod rotational;
mod skin;
mod taper;

pub use super::p1::PrimitiveCreation;
pub use access::PrimitiveCreationAccess;

use crate::curve::curve::{Curve, CurveKind};
use crate::nurbs::b_spline_curve::p1::BSplineCurve;

/// The B-spline an `SmBSplineCurve` cast reads: B-splines and the analytic kinds
/// built on one (`None` for the kinds upstream's cast would misread).
fn as_bspline(curve: &Curve) -> Option<BSplineCurve> {
    match curve.kind() {
        CurveKind::BSpline(c) => Some(c.clone()),
        CurveKind::Line(c) => Some(c.base().clone()),
        CurveKind::Ellipse(c) => Some(c.base().clone()),
        CurveKind::Parabola(c) => Some(c.base().clone()),
        CurveKind::Hyperbola(c) => Some(c.base().clone()),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    //! Entry checks of the p2 primitives on a B-rep with only its infinite region.
    //! Difftest compares these exits too; here a removed check cannot hide behind a
    //! deferred construction.
    use super::*;
    use crate::curve::composite_curve::OffsetMapAttribute;
    use crate::error::{Result, SolidError};
    use crate::geom::attribute::TagData;
    use crate::geom::axis2_placement::Axis2Placement;
    use crate::geom::context::MarkType;
    use crate::geom::object::ObjectRef;
    use crate::kernel::Kernel;
    use crate::math::{Point3, Vec3};
    use crate::nurbs::b_spline_curve::p1::NurbCircleParam;
    use crate::primitive::primitive_creation::{p3::PrimitiveAccess, SweepOptions};
    use crate::sweep::topology_sweep::access::{EdgeRef, FaceRef, SweepAccess, SweepMergeId};
    use crate::sweep::topology_sweep::SweepEntity;
    use crate::topology::brep::p10::TrimType;
    use crate::topology::edge::NotifyEnv;
    use crate::topology::{Brep, BrepId, FaceId, Region, RegionId};

    struct Model {
        k: Kernel,
        b: Brep,
    }
    fn none<T>() -> Result<T> {
        Err(SolidError::invalid_input("p2 test model", "no topology"))
    }
    impl Model {
        fn new() -> Result<Self> {
            let mut b = Brep::new();
            let region = b.add_region_record(Region {
                in_brep: true,
                ..Default::default()
            })?;
            b.set_infinite_region(Some(region));
            Ok(Self {
                k: Kernel::default(),
                b,
            })
        }
    }
    impl SweepAccess for Model {
        fn brep(&self, _: BrepId) -> Result<&Brep> {
            Ok(&self.b)
        }
        fn brep_mut(&mut self, _: BrepId) -> Result<&mut Brep> {
            Ok(&mut self.b)
        }
        fn kernel(&self, _: BrepId) -> Result<&Kernel> {
            Ok(&self.k)
        }
        fn kernel_mut(&mut self, _: BrepId) -> Result<&mut Kernel> {
            Ok(&mut self.k)
        }
        fn with_brep_mut<T>(
            &mut self,
            _: BrepId,
            _: impl FnOnce(
                &mut Kernel,
                &mut Brep,
                &dyn Fn(FaceId) -> ObjectRef,
                &mut NotifyEnv<'_>,
            ) -> Result<T>,
        ) -> Result<T> {
            none()
        }
        fn with_breps_mut<T>(
            &mut self,
            _: BrepId,
            _: BrepId,
            _: impl FnOnce(
                &mut Kernel,
                &mut Brep,
                &mut Brep,
                &dyn Fn(FaceId) -> ObjectRef,
                &mut NotifyEnv<'_>,
            ) -> Result<T>,
        ) -> Result<T> {
            none()
        }
        fn new_brep(&mut self, _: BrepId) -> Result<BrepId> {
            Err(SolidError::unported("SmBrep::SmBrep", "topology.brep_p3"))
        }
        fn delete_brep(&mut self, _: BrepId) {}
        fn delete_merge(&mut self, _: SweepMergeId) {}
        fn tag(&self, _: SweepEntity) -> Result<Option<TagData>> {
            none()
        }
        fn add_tag(&mut self, _: SweepEntity, _: TagData) -> Result<()> {
            none()
        }
        fn curve_tag(&self, _: EdgeRef) -> Result<Option<TagData>> {
            none()
        }
        fn set_editing(&mut self, _: BrepId, editing: bool) -> Result<bool> {
            let previous = self.b.editing_enabled();
            self.b.set_editing_enabled(editing);
            Ok(previous)
        }
        fn lock_mark(&mut self, _: BrepId, _: i32) -> Result<MarkType> {
            none()
        }
        fn unlock_mark(&mut self, _: BrepId, _: MarkType) -> Result<()> {
            none()
        }
        fn is_marked(&self, _: SweepEntity, _: MarkType) -> Result<bool> {
            none()
        }
        fn mark(&mut self, _: SweepEntity, _: MarkType) -> Result<()> {
            none()
        }
        fn collect_faces(&mut self, _: FaceRef, _: MarkType, _: &mut Vec<FaceId>) -> Result<()> {
            none()
        }
    }
    impl PrimitiveAccess for Model {
        fn delete_primitive_faces(&mut self, _: BrepId, _: &[FaceId]) -> Result<()> {
            none()
        }
        fn with_primitive_brep<T>(
            &mut self,
            _: BrepId,
            body: impl FnOnce(&mut Kernel, &mut Brep) -> Result<T>,
        ) -> Result<T> {
            body(&mut self.k, &mut self.b)
        }
        fn axis_inside_face(&mut self, _: FaceRef, _: &Curve) -> Result<bool> {
            none()
        }
    }
    impl PrimitiveCreationAccess for Model {
        fn copy_faces_into(&mut self, _: BrepId, _: &[FaceId], _: BrepId) -> Result<()> {
            none()
        }
        fn attach_offset_map(&mut self, _: FaceRef, _: OffsetMapAttribute) -> Result<()> {
            none()
        }
        fn offset_map(&self, _: FaceRef) -> Result<Option<OffsetMapAttribute>> {
            none()
        }
        fn project_and_trim(
            &mut self,
            _: BrepId,
            _: &BSplineCurve,
            _: Vec3,
            _: Point3,
            _: TrimType,
            _: f64,
        ) -> Result<()> {
            none()
        }
    }

    const B: BrepId = BrepId::from_index(0);
    fn primitive() -> PrimitiveCreation {
        PrimitiveCreation::new(Some(RegionId::from_index(0)), 0, false)
    }
    fn code<T>(result: &Result<T>) -> Option<i32> {
        match result {
            Err(SolidError::Upstream { code, .. }) => Some(*code),
            Err(SolidError::InvalidInput { .. }) => Some(1007),
            _ => None,
        }
    }

    #[test]
    fn partial_cone_checks_height_and_turn_before_edit_mode() -> Result<()> {
        let frame = Axis2Placement::IDENTITY;
        let mut a = Model::new()?;
        let low = primitive().create_partial_cone(&mut a, B, 0.0, 1.0, 1.0, 0.0, 90.0, &frame);
        assert_eq!((code(&low), a.b.editing_enabled()), (Some(1001), false));
        let unturned =
            primitive().create_partial_cone(&mut a, B, 1.0, 1.0, 1.0, 30.0, 30.0, &frame);
        assert_eq!(
            (code(&unturned), a.b.editing_enabled()),
            (Some(1001), false)
        );
        // an axis shorter than the scaled zero fails after edit mode is on, as upstream
        let short =
            primitive().create_partial_cone(&mut a, B, 1.5e-12, 1.0, 1.0, 0.0, 90.0, &frame);
        assert_eq!((code(&short), a.b.editing_enabled()), (Some(1001), true));
        Ok(())
    }

    #[test]
    fn a_skin_of_one_profile_fails_with_cleared_outputs_and_report() -> Result<()> {
        let mut a = Model::new()?;
        let circle = BSplineCurve::create_circle_segment(
            3,
            &Axis2Placement::IDENTITY,
            5.0,
            0.0,
            360.0,
            NurbCircleParam::Quadratic,
        )?;
        let mut curves = vec![Some(Curve::from(circle))];
        let mut starts = vec![FaceId::from_index(9)];
        let mut sides = vec![FaceId::from_index(9)];
        let mut ends = Vec::new();
        let mut transferred = vec![true, true, true];
        let result = primitive().create_skin_primitive(
            &mut a,
            B,
            &[1],
            &mut curves,
            0.0,
            3,
            true,
            &mut starts,
            &mut sides,
            &mut ends,
            &mut transferred,
        );
        assert_eq!(code(&result), Some(1001));
        assert!(starts.is_empty() && sides.is_empty() && ends.is_empty());
        assert_eq!(transferred, vec![false]);
        assert!(curves.first().is_some_and(Option::is_some));
        Ok(())
    }

    #[test]
    fn sweeps_of_nothing_and_a_planar_path_without_paths_fail() -> Result<()> {
        let mut a = Model::new()?;
        let path = BSplineCurve::create_line_segment(3, Vec3::ZERO, Vec3::Z)?;
        let (mut s, mut d, mut e) = (Vec::new(), Vec::new(), Vec::new());
        let edges = primitive().create_curve_sweep_from_edges(
            &mut a,
            B,
            &[],
            &path,
            None,
            None,
            1e-5,
            SweepOptions::default(),
            &mut s,
            &mut d,
            &mut e,
        );
        assert_eq!(code(&edges), Some(1001));
        let faces = primitive().create_curve_sweep_from_faces(
            &mut a,
            B,
            &[],
            &path,
            None,
            None,
            1e-5,
            SweepOptions::default(),
            &mut s,
            &mut d,
            &mut e,
        );
        assert_eq!(code(&faces), Some(1001));
        let planar = PrimitiveCreation::create_sweep_along_planar_path(
            &mut a,
            B,
            &[],
            &[],
            1e-5,
            true,
            true,
            true,
        );
        assert_eq!(code(&planar), Some(1007));
        // without a target region the pipe still resets its lists first
        let mut orphan = PrimitiveCreation::new(None, 0, false);
        let mut stale = vec![FaceId::from_index(1)];
        let pipe = orphan.create_pipe_sweep(
            &mut a, B, 1.0, &path, 0.0, true, false, &mut stale, &mut d, &mut e,
        );
        assert_eq!(code(&pipe), Some(1006));
        assert!(stale.is_empty());
        Ok(())
    }
}
