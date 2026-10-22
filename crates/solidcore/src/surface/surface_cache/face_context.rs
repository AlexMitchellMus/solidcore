//! Operation-scoped access to the mutable B-rep and face queries needed by trim caches.
//!
//! Upstream: `source/SMLib/src/SmTrimSrfCache.cpp`, `source/SMLib/src/SmFace.cpp`.
//!
//! | Rust | Upstream |
//! |------|----------|
//! | [`FaceCacheContext`] | `GetFace()->GetBrep()` during `SmTrimSrfCache::BuildTree` |
//! | [`FaceQueries`], [`PreparedFaceQueries`] | `SmFace::PointClassify`, `SmTol::GetZoneTol3d` |
//!
//! Deviations: borrowed operation contexts replace topology pointers and virtual
//! calls; typed entities are mapped by the caller, never decoded from object keys.
//! No topology is retained by a cache; callers provide this borrowed context only while
//! preparing a cache or servicing a face classification.

use core::fmt::{self, Debug};

use crate::error::{EntityRef, Result, SolidError};
use crate::geom::object::ObjectRef;
use crate::intersect::point_class::PointClassification;
use crate::kernel::Kernel;
use crate::math::Point2;
use crate::query::face::p2::{point_classify_prepared, PointClassifyOptions};
use crate::topology::{Brep, FaceId};

/// Face operations needed by a trimmed surface cache.
pub trait FaceQueries<O>: Debug {
    /// Returns the face's 3D zone tolerance.
    fn face_zone_tol_3d(&self, face: O) -> Result<f64>;

    /// Classifies a UV point against a face, returning its class and optional object.
    fn face_point_classify(
        &self,
        k: &Kernel,
        face: O,
        uv: Point2,
        src_tol: f64,
    ) -> Result<(i32, Option<O>)>;

    /// `GlobalPointSolve`'s fast-drop classification uses boundary intersections
    /// on and forced UV classification off, unlike the trim-cache override.
    /// Adapters must explicitly support these flags rather than dropping them.
    fn face_point_classify_fast_drop(
        &self,
        _k: &Kernel,
        _face: O,
        _uv: Point2,
        _src_tol: f64,
    ) -> Result<(i32, Option<O>)> {
        // UNPORTED(query.face_p2): fast-drop caller flags on an opaque adapter.
        Err(SolidError::unported(
            "SmFace::PointClassify fast-drop adapter",
            "query.face_p2",
        ))
    }
}

/// Borrowed mutable B-rep context for preparing one face's trim cache.
pub struct FaceCacheContext<'a, O> {
    /// The mutable topology being prepared.
    pub brep: &'a mut Brep,
    /// Face whose surface is cached.
    pub face: FaceId,
    /// Maps a typed entity to its caller object key.
    pub key_of: &'a dyn Fn(EntityRef) -> ObjectRef,
    /// Resolves an entity to a caller handle.
    pub handle_of: &'a dyn Fn(EntityRef) -> Result<O>,
}

/// Read-only face queries after UV trim preparation has finished.
///
/// This adapter borrows the real topology, not a B-rep copy. Its classifier cannot
/// create trims or re-enter construction of the cache being classified.
pub struct PreparedFaceQueries<'a, O> {
    brep: &'a Brep,
    face: FaceId,
    face_handle: O,
    handle_of: &'a dyn Fn(EntityRef) -> Result<O>,
}

impl<O> Debug for PreparedFaceQueries<'_, O> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("PreparedFaceQueries")
            .field("face", &self.face)
            .finish_non_exhaustive()
    }
}

impl<'a, O: Copy + PartialEq> PreparedFaceQueries<'a, O> {
    /// Borrows one prepared face and the caller's typed-entity handle mapping.
    ///
    /// # Errors
    ///
    /// Those of looking up the face and mapping its handle. Missing trim geometry
    /// is reported by classification rather than synthesized by this adapter.
    pub fn new(
        brep: &'a Brep,
        face: FaceId,
        handle_of: &'a dyn Fn(EntityRef) -> Result<O>,
    ) -> Result<Self> {
        brep.face(face)?;
        Ok(Self {
            brep,
            face,
            face_handle: handle_of(EntityRef::Face(face))?,
            handle_of,
        })
    }

    fn check_face(&self, face: O) -> Result<()> {
        if face == self.face_handle {
            Ok(())
        } else {
            Err(SolidError::invalid_input(
                "surface::surface_cache::face_context::PreparedFaceQueries",
                "query handle does not name the prepared face",
            ))
        }
    }

    fn classify(
        &self,
        k: &Kernel,
        face: O,
        uv: Point2,
        src_tol: f64,
        do_boundary_intersections: bool,
        use_uv_classification: bool,
    ) -> Result<(i32, Option<O>)> {
        self.check_face(face)?;
        let mut pc = PointClassification::new(src_tol);
        point_classify_prepared(
            k,
            self.brep,
            self.face,
            uv,
            &mut pc,
            PointClassifyOptions {
                src_zone_tol_3d: src_tol,
                do_boundary_intersections,
                use_uv_classification,
                ..PointClassifyOptions::default()
            },
        )?;
        Ok((
            pc.point_class().code(),
            pc.object().map(self.handle_of).transpose()?,
        ))
    }
}

impl<O: Copy + PartialEq> FaceQueries<O> for PreparedFaceQueries<'_, O> {
    fn face_zone_tol_3d(&self, face: O) -> Result<f64> {
        self.check_face(face)?;
        crate::topology::face::p2::tolerance(self.brep, self.face)
    }

    fn face_point_classify(
        &self,
        k: &Kernel,
        face: O,
        uv: Point2,
        src_tol: f64,
    ) -> Result<(i32, Option<O>)> {
        // The classification object's tolerance is read by vertex/loop tests;
        // setting only the options' source tolerance is not equivalent.
        self.classify(k, face, uv, src_tol, false, true)
    }

    fn face_point_classify_fast_drop(
        &self,
        k: &Kernel,
        face: O,
        uv: Point2,
        src_tol: f64,
    ) -> Result<(i32, Option<O>)> {
        self.classify(k, face, uv, src_tol, true, false)
    }
}
