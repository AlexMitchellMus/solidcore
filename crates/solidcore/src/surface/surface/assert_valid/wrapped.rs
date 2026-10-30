//! The validity checks of the surfaces that wrap another surface (`SmOffsetSurface`,
//! `SmSrfInVolume`) and of the volumes a surface in a volume nests into.
//!
//! Upstream: `source/SMLib/src/SmOffsetSurface.cpp`, `source/SMLib/src/SmSrfInVolume.cpp`,
//! `source/SMLib/src/SmVolume.cpp`, `source/SMLib/src/SmTransform.cpp`,
//! `source/SMLib/src/SmBSplineVolume.cpp`, `source/SMLib/src/SmBendVolume.cpp`,
//! `source/SMLib/src/SmUnbendVolume.cpp`, `source/SMLib/src/SmTwistVolume.cpp`.
//!
//! | Rust                                     | Upstream                              |
//! |------------------------------------------|---------------------------------------|
//! | `offset_surface_assert_valid` (private)  | `SmOffsetSurface::AssertValid`        |
//! | `srf_in_volume_assert_valid` (private)   | `SmSrfInVolume::AssertValid`          |
//! | [`Volume::assert_valid`]                 | `SmVolume::AssertValid` (virtual)     |
//! | `volume_checks` (private)                | `SmVolume::AssertValid`               |
//! | `transform_assert_valid` (private)       | `SmTransform::AssertValid` (an orient map) |
//! | `b_spline_volume_assert_valid` (private) | `SmBSplineVolume::AssertValid`        |
//! | `frame_checks` (private)                 | `SmBendVolume`, `SmUnbendVolume`, `SmTwistVolume::AssertValid` (tests 0..4) |
//!
//! Deviations:
//! - **Contained identities.** Reports about contained objects name keys derived from
//!   the container's key `key`: `8 * key + slot` (wrapping), with slot 1 a surface in a
//!   volume's base surface, 2 its volume, 3 a volume's next map, 4 its orient map and
//!   5 its inverse orient map (`contained_key`). Upstream compares object addresses;
//!   only distinctness matters.
//! - **Context tests** (`SmOffsetSurface` 0, `SmSrfInVolume` 4, 5) compare contexts:
//!   one kernel is the only context, so only their owner-flag part can fail.
//! - **Owner links.** `SmSrfInVolume` 2, 3 read the contained objects' owners, which
//!   upstream sets to this object wherever an owner bit is set (constructors,
//!   `SetOwnerFlag`, the setters, `ReadFromDB`); the port does not model those links,
//!   so the tests are the owner bits. `SmOffsetSurface` 1 is that the offset owns its
//!   base and the base's owner is [`BASE_OWNER`] (a shared base's owner is its caller's).
//! - **Self-intersection** (`SmOffsetSurface` 2): upstream ignores the status of
//!   `GlobalSurfaceIntersect` and counts the curve types it left in the array; the port
//!   counts those of [`Surface::global_surface_intersect_with_partial`] the same way and
//!   propagates only [`SolidError::Unported`].
//! - **`SmBSplineVolume`** test 1 calls `PassesValidityCheck`, which returns `SM_ERR`
//!   (non-zero: passes) in the build without `USE_VOLUMES`; test 0 (`m_pNurb` non-NULL)
//!   always holds (the record is held by value).
//! - **Orient maps** are [`Transform`](crate::geom::transform::Transform) values without maps of their own, so the
//!   `SmVolume` tests 3..8 always hold and an orient map's own check reports nothing.
//! - **Recursion.** A volume's next-map chain is walked with an explicit depth checked
//!   against [`MAX_GEOMETRY_NESTING`]; before a surface in a volume checks its base
//!   surface, its whole owned geometry is checked to nest at most that deep
//!   ([`check_geometry_nesting`]), which bounds the surface recursion.

use super::{a_object, Reporter, LEVEL_GIVEN, UNDEF};
use crate::error::{Result, SolidError};
use crate::geom::assert_array::labels::{
    BEND_VOLUME_TYPE, B_SPLINE_VOLUME_TYPE, OFFSET_SURFACE_TYPE, SRF_IN_VOLUME_TYPE,
    TWIST_VOLUME_TYPE, UNBEND_VOLUME_TYPE, VOLUME_TYPE,
};
use crate::geom::assert_array::{AssertArray, ObjectKey};
use crate::geom::tol::ZONE_TOL_3D;
use crate::kernel::Kernel;
use crate::math::{Vec3, EFF_ZERO};
use crate::surface::offset_surface::BASE_OWNER;
use crate::surface::surface::p2::{check_geometry_nesting, NestedGeometry, MAX_GEOMETRY_NESTING};
use crate::surface::surface::{Surface, SurfaceKind};
use crate::surface::volume::{Volume, VolumeKind};

const OP: &str = "surface::surface::assert_valid::wrapped";

fn invalid(reason: &'static str) -> SolidError {
    SolidError::invalid_input(OP, reason)
}

/// The key of the contained object in `slot` (1..=5, see the module docs) of the
/// object `key`.
const fn contained_key(key: ObjectKey, slot: u64) -> ObjectKey {
    ObjectKey(key.0.wrapping_mul(8).wrapping_add(slot))
}

/// `GetClassString()` of an orient map (an `SmTransform`).
const TRANSFORM_CLASS: &str = "SmTransform";

// ---- SmOffsetSurface --------------------------------------------------------------------

/// `SmOffsetSurface::AssertValid`.
pub(super) fn offset_surface_assert_valid(
    k: &Kernel,
    s: &Surface,
    r: &mut Reporter<'_>,
) -> Result<bool> {
    let SurfaceKind::OffsetSurface(off) = s.kind() else {
        return Err(invalid("not an offset surface"));
    };
    let mut rtn = true;

    // call the base class AssertValid
    rtn &= if r.bases() {
        super::surface_assert_valid(k, s, r)?
    } else {
        true
    };

    /*  0 */
    // m_pSurface shares the same context (one kernel)
    rtn &= r.report(OFFSET_SURFACE_TYPE, 0, 0, true, UNDEF, UNDEF)?;

    /*  1 */
    // Bad m_pSurface Owner
    let owner_ok = off.owns_surface() && off.base_surface().owner() == Some(BASE_OWNER);
    rtn &= r.report(OFFSET_SURFACE_TYPE, 1, 0, owner_ok, UNDEF, UNDEF)?;

    /*  2 */
    // Self-intersecting surface interior (status ignored: the curve types left count)
    let uv_domain = s.natural_uv_domain_with_nurb(k)?;
    let (xsect, status) = s.global_surface_intersect_with_partial(
        k,
        &uv_domain,
        s,
        &uv_domain,
        [false, false],
        None,
        None,
    );
    if let Err(e) = status {
        if e.is_unported() {
            return Err(e);
        }
    }
    // If the surface does not intersect, then pOptCurveTypes is an empty array
    rtn &= r.report(
        OFFSET_SURFACE_TYPE,
        2,
        0,
        xsect.curve_types.is_empty(),
        UNDEF,
        UNDEF,
    )?;

    // all done
    Ok(rtn)
}

// ---- SmSrfInVolume ----------------------------------------------------------------------

/// `SmSrfInVolume::AssertValid` (no base-class call upstream).
pub(super) fn srf_in_volume_assert_valid(
    k: &Kernel,
    s: &Surface,
    r: &mut Reporter<'_>,
) -> Result<bool> {
    let SurfaceKind::SrfInVolume(siv) = s.kind() else {
        return Err(invalid("not a surface in a volume"));
    };
    // bounds the recursion into the base surface (and its own wrappers)
    check_geometry_nesting(NestedGeometry::Surface(s), 0, OP)?;
    let owner_flag = siv.owner_flag();

    /*  0 */
    // test contained surface
    let surface_ok =
        siv.surface()
            .assert_valid(k, r.list.as_deref_mut(), r.level, contained_key(r.key, 1))?;
    let mut rtn = r.report(SRF_IN_VOLUME_TYPE, 0, 0, surface_ok, UNDEF, UNDEF)?;

    /*  1 */
    // test contained volume
    let volume_ok =
        siv.volume()
            .assert_valid(r.list.as_deref_mut(), r.level, contained_key(r.key, 2))?;
    rtn &= r.report(SRF_IN_VOLUME_TYPE, 1, 0, volume_ok, UNDEF, UNDEF)?;

    /*  2 */
    // Bad BaseSurface Owner (the owner is this object whenever bit 1 is set)
    rtn &= r.report(SRF_IN_VOLUME_TYPE, 2, 0, owner_flag & 1 != 0, UNDEF, UNDEF)?;

    /*  3 */
    // Bad BaseVolume Owner (the owner is this object whenever bit 2 is set)
    rtn &= r.report(SRF_IN_VOLUME_TYPE, 3, 0, owner_flag & 2 != 0, UNDEF, UNDEF)?;

    /*  4 */
    // Bad BaseSurface Context (one kernel)
    rtn &= r.report(SRF_IN_VOLUME_TYPE, 4, 0, owner_flag & 1 != 0, UNDEF, UNDEF)?;

    /*  5 */
    // Bad BaseVolume Context (one kernel)
    rtn &= r.report(SRF_IN_VOLUME_TYPE, 5, 0, owner_flag & 2 != 0, UNDEF, UNDEF)?;

    // all done
    Ok(rtn)
}

// ---- volumes ----------------------------------------------------------------------------

impl Volume {
    /// `AssertValid` (virtual): the checks of this volume's class and of its base
    /// classes (unless `level` is [`LEVEL_GIVEN`], which skips the `SmAObject` call only,
    /// as upstream), including the checks of its next map and orient maps, reporting
    /// failures about the volume `key` into `list`. Returns `true` when every test that
    /// ran passed.
    ///
    /// # Errors
    ///
    /// `InvalidInput` for a next-map chain deeper than [`MAX_GEOMETRY_NESTING`]
    /// (the constructors and setters never build one).
    pub fn assert_valid(
        &self,
        list: Option<&mut AssertArray>,
        level: u8,
        key: ObjectKey,
    ) -> Result<bool> {
        volume_assert_valid_at(self, list, level, key, 0)
    }
}

/// [`Volume::assert_valid`] of a volume `depth` maps down a next-map chain.
fn volume_assert_valid_at(
    v: &Volume,
    list: Option<&mut AssertArray>,
    level: u8,
    key: ObjectKey,
    depth: usize,
) -> Result<bool> {
    if depth > MAX_GEOMETRY_NESTING {
        return Err(invalid(
            "volume next-map chain exceeds MAX_GEOMETRY_NESTING",
        ));
    }
    let mut r = Reporter {
        list,
        key,
        level,
        class_string: v.class_string(),
    };
    let r = &mut r;
    match v.kind() {
        // SmVolume, SmTransform (no checks of its own)
        VolumeKind::Base | VolumeKind::Transform(_) => volume_checks(v, r, depth),
        VolumeKind::BSpline(_) => b_spline_volume_assert_valid(v, r, depth),
        VolumeKind::Bend(b) => {
            // check base class
            let mut rtn = volume_checks(v, r, depth)?;
            // Run SmBendVolume checks
            let axes = [v.bend_axis(), v.bend_mid_dir(), v.bend_binormal()];
            rtn &= frame_checks(r, BEND_VOLUME_TYPE, axes)?;
            // BendNeutralDist must be positive definite
            let d = b.bend_neutral_dist();
            let tol = b.singularity_tol();
            rtn &= r.report(BEND_VOLUME_TYPE, 5, 0, d >= tol, EFF_ZERO, d)?;
            // Bend SingularityTol must be positive definite
            rtn &= r.report(BEND_VOLUME_TYPE, 6, 0, tol >= ZONE_TOL_3D, EFF_ZERO, tol)?;
            Ok(rtn)
        }
        VolumeKind::Unbend(u) => {
            // check base class
            let mut rtn = volume_checks(v, r, depth)?;
            // Run SmUnbendVolume checks
            let axes = [v.unbend_axis(), v.unbend_mid_dir(), v.unbend_binormal()];
            rtn &= frame_checks(r, UNBEND_VOLUME_TYPE, axes)?;
            // UnbendNeutralDist must be positive definite
            let d = u.unbend_neutral_dist();
            let tol = u.singularity_tol();
            rtn &= r.report(UNBEND_VOLUME_TYPE, 5, 0, d >= tol, EFF_ZERO, d)?;
            // Unbend SingularityTol must be positive definite
            rtn &= r.report(UNBEND_VOLUME_TYPE, 6, 0, tol >= ZONE_TOL_3D, EFF_ZERO, tol)?;
            Ok(rtn)
        }
        VolumeKind::Twist(t) => {
            // check base class
            let mut rtn = volume_checks(v, r, depth)?;
            // Run SmTwistVolume checks
            let axes = [v.twist_axis(), v.twist_mid_dir(), v.twist_binormal()];
            rtn &= frame_checks(r, TWIST_VOLUME_TYPE, axes)?;
            // Twist PeriodDistance abs value has to be greater than PeriodMinDist
            let p = t.twist_period_dist();
            let ok = p.abs() >= t.period_min_dist();
            rtn &= r.report(TWIST_VOLUME_TYPE, 5, 0, ok, EFF_ZERO, p)?;
            Ok(rtn)
        }
    }
}

/// `SmVolume::AssertValid` of the volume `v` (`r` reports about it).
fn volume_checks(v: &Volume, r: &mut Reporter<'_>, depth: usize) -> Result<bool> {
    // check base class validity
    let mut rtn = true;

    // call the base class AssertValid
    rtn &= if r.level != LEVEL_GIVEN {
        a_object::assert_valid(false)
    } else {
        true
    };

    // check NextMap
    let ok = match v.next_map() {
        None => true,
        Some(next) => volume_assert_valid_at(
            next,
            r.list.as_deref_mut(),
            r.level,
            contained_key(r.key, 3),
            depth + 1,
        )?,
    };
    rtn &= ok;
    r.report(VOLUME_TYPE, 0, 0, rtn, UNDEF, UNDEF)?;

    // check OrientMap
    let ok = match v.orient_map() {
        None => true,
        Some(_) => transform_assert_valid(r.list.as_deref_mut(), r.level, contained_key(r.key, 4))?,
    };
    rtn &= ok;
    r.report(VOLUME_TYPE, 1, 0, rtn, UNDEF, UNDEF)?;

    // check InvOrientMap
    let ok = match v.inv_orient_map() {
        None => true,
        Some(_) => transform_assert_valid(r.list.as_deref_mut(), r.level, contained_key(r.key, 5))?,
    };
    rtn &= ok;
    r.report(VOLUME_TYPE, 2, 0, rtn, UNDEF, UNDEF)?;

    // simple OrientMap check - Orient maps are simple, i.e. they have no Next, Orient, or
    // InvOrient maps (a Transform value holds none)
    for test in 3..=8 {
        r.report(VOLUME_TYPE, test, 0, true, UNDEF, UNDEF)?;
    }

    // all done
    Ok(rtn)
}

/// `SmTransform::AssertValid` of an orient map (a [`Transform`](crate::geom::transform::Transform) without maps: the
/// `SmVolume` checks, which all pass, and none of its own).
fn transform_assert_valid(
    list: Option<&mut AssertArray>,
    level: u8,
    key: ObjectKey,
) -> Result<bool> {
    let mut r = Reporter {
        list,
        key,
        level,
        class_string: TRANSFORM_CLASS,
    };
    let r = &mut r;
    // SmVolume::AssertValid: the base class, no next or orient maps
    let rtn = if r.level != LEVEL_GIVEN {
        a_object::assert_valid(false)
    } else {
        true
    };
    for test in 0..=2 {
        r.report(VOLUME_TYPE, test, 0, rtn, UNDEF, UNDEF)?;
    }
    for test in 3..=8 {
        r.report(VOLUME_TYPE, test, 0, true, UNDEF, UNDEF)?;
    }
    // SmTransform: no checks of its own
    Ok(rtn)
}

/// `SmBSplineVolume::AssertValid`.
fn b_spline_volume_assert_valid(v: &Volume, r: &mut Reporter<'_>, depth: usize) -> Result<bool> {
    // check base class
    let mut rtn = volume_checks(v, r, depth)?;

    // run SmBSplineVolume PassesValidityCheck checks (SM_ERR without USE_VOLUMES: true)
    rtn &= r.report(B_SPLINE_VOLUME_TYPE, 1, 0, true, UNDEF, UNDEF)?;

    // check for nonNULL m_pNurb pointer (held by value)
    rtn &= r.report(B_SPLINE_VOLUME_TYPE, 0, 0, true, UNDEF, UNDEF)?;

    // all done
    Ok(rtn)
}

/// Tests 0..4 of the bend, unbend and twist checks of class `class`: `axes` (axis,
/// mid direction, binormal) are unit vectors forming a right-handed orthonormal frame.
fn frame_checks(r: &mut Reporter<'_>, class: i64, axes: [Vec3; 3]) -> Result<bool> {
    let [axis, mid, binormal] = axes;
    let mut rtn = true;

    // Axis is a unit vector
    let d = (axis.dot(axis) - 1.0).abs();
    rtn &= r.report(class, 0, 0, d < EFF_ZERO, EFF_ZERO, d)?;

    // MidDir is a unit vector
    let d = (mid.dot(mid) - 1.0).abs();
    rtn &= r.report(class, 1, 0, d < EFF_ZERO, EFF_ZERO, d)?;

    // Binormal is a unit vector
    let d = (binormal.dot(binormal) - 1.0).abs();
    rtn &= r.report(class, 2, 0, d < EFF_ZERO, EFF_ZERO, d)?;

    // basis triple product - should be +1
    let triple_product = mid.triple_product(binormal, axis);

    // basis vectors have to be mutually orthogonal
    let d = triple_product.abs() - 1.0;
    rtn &= r.report(class, 3, 0, d < EFF_ZERO, EFF_ZERO, d)?;

    // basis vectors have to be right handed
    let d = triple_product - 1.0;
    rtn &= r.report(class, 4, 0, d < EFF_ZERO, EFF_ZERO, d)?;

    Ok(rtn)
}
