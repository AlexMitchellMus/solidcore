//! Live object storage for one fillet executive.
//!
//! Upstream: `source/SMLib/src/SmFilletExecutive.cpp`, `source/SMLib/inc/SmFilletExecutive.h`.
//!
//! | Rust | Upstream |
//! |------|----------|
//! | [`FilletSession`] | `SmFilletExecutive` and its owned objects |
//! | [`SolverObject`] | the concrete `SmFilletSolver` objects |
//!
//! Deviations: one borrowed target and typed owned slots replace native pointers.
//! Slots are never reused. Target curves cannot be replaced while its loading caches
//! are retained: the session exposes only shared access to its borrowed target.
//! Loading creates and retains default full-tree curve caches: only their knot
//! continuities are consumed, calculated at a fixed 0.5 degrees before tessellation
//! (`SmCurveCache.cpp:1508`, `curve_cache.rs:608`). The frozen target makes repeated
//! loads equivalent. An external custom/non-full-tree cache provider is not yet
//! accepted; this is not an arbitrary warmed-cache parity claim.
//! This foundation supports solver construction/loading; it does not substitute a
//! result for the remaining fillet execution stages. No recursion.

use std::collections::BTreeMap;

use super::p1::{ExecutiveObjects, SolverLoadingAccess};
use super::p3::{ErrorInfo, FilletExecutive};
use crate::blend::fillet_corner::p3::FilletCorner;
use crate::blend::fillet_geom::p1::{FilletBrep, FilletCornerId, FilletSolverId};
use crate::blend::fillet_solver::p1::{FilletSolver, FilletSurfaceGenerator};
use crate::blend::fillet_standard_solver::p1::{
    constant_radius_from_edgeuse, ConstantDistanceFs, ConstantRadiusFs, SurfaceSurfaceFs,
    VariableRadiusFs,
};
use crate::curve::curve_cache::CurveCache;
use crate::error::{EntityRef, Result, SolidError};
use crate::geom::attribute::{ids, Attribute, AttributeHandle, AttributeStore, Behavior};
use crate::geom::geom_utility::p2::ContinuityType;
use crate::geom::solution_array::{ObjectKind, SolutionObject};
use crate::intersect::point_class::{ClassObject, ClassObjectKind};
use crate::kernel::Kernel;
use crate::topology::{Brep, EdgeId, EdgeuseId};

const OP: &str = "blend::fillet_executive::session";

fn invalid(reason: &'static str) -> SolidError {
    SolidError::invalid_input(OP, reason)
}

/// Which actual body owns a session topology reference.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum FilletBody {
    /// The caller's borrowed target B-rep.
    Target,
    /// The executive's owned `SmFilletBrep`.
    Pseudo,
}

/// A reference to an actual object in a single executive's object store.
///
/// These identities are session-local. Bare solver/corner IDs and object values from
/// different sessions must not be mixed; they do not encode global provenance.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum FilletObject {
    /// The session's one actual executive.
    Executive,
    /// One of its actual B-reps.
    Body(FilletBody),
    /// Topology in that body's typed arena (different bodies never alias).
    Topology(FilletBody, EntityRef),
    /// The retained curve of an edge; separate from its topology identity.
    Curve(FilletBody, EdgeId),
    /// A concrete solver in the session's solver slots.
    Solver(FilletSolverId),
    /// A concrete corner in the session's corner slots.
    Corner(FilletCornerId),
    /// An actual retained cross-section generator.
    Generator(u32),
    /// An attribute owned by the session's actual attribute store.
    Attribute(AttributeHandle),
}

impl ClassObject for FilletObject {
    fn class_kind(&self) -> ClassObjectKind {
        match self {
            Self::Curve(..) => ClassObjectKind::Curve,
            Self::Topology(_, entity) => match entity {
                EntityRef::Vertex(_) => ClassObjectKind::Vertex,
                EntityRef::Edge(_) => ClassObjectKind::Edge,
                EntityRef::Face(_) => ClassObjectKind::Face,
                EntityRef::Region(_) => ClassObjectKind::Region,
                EntityRef::Loop(_) => ClassObjectKind::Loop,
                EntityRef::Edgeuse(_) => ClassObjectKind::Edgeuse,
                EntityRef::Vertexuse(_) => ClassObjectKind::Vertexuse,
                EntityRef::Faceuse(_) => ClassObjectKind::Faceuse,
                EntityRef::Loopuse(_) => ClassObjectKind::Loopuse,
                EntityRef::Brep(_) | EntityRef::Shell(_) => ClassObjectKind::Other,
            },
            _ => ClassObjectKind::Other,
        }
    }
}

impl SolutionObject for FilletObject {
    fn kind(&self) -> ObjectKind {
        match self {
            Self::Curve(..) => ObjectKind::Curve,
            Self::Topology(_, EntityRef::Vertex(_)) => ObjectKind::Vertex,
            Self::Topology(_, EntityRef::Edge(_)) => ObjectKind::Edge,
            Self::Topology(_, EntityRef::Face(_)) => ObjectKind::Face,
            _ => ObjectKind::Other,
        }
    }
}

/// An owned concrete solver, including its one authoritative mutable base.
#[derive(Debug)]
pub enum SolverObject {
    /// `SmConstantRadiusFS` / its assisted implementation.
    ConstantRadius(Box<ConstantRadiusFs>),
    /// `SmConstantDistanceFS`.
    ConstantDistance(Box<ConstantDistanceFs>),
    /// `SmVariableRadiusFS`, retaining its actual radius law.
    VariableRadius(Box<VariableRadiusFs>),
    /// `SmSurfaceSurfaceFS`.
    SurfaceSurface(Box<SurfaceSurfaceFs>),
}

impl SolverObject {
    /// Borrow the base of this actual concrete solver.
    #[must_use]
    pub fn base(&self) -> &FilletSolver {
        match self {
            Self::ConstantRadius(solver) => &solver.solver,
            Self::ConstantDistance(solver) => &solver.solver,
            Self::VariableRadius(solver) => &solver.curve_based.solver,
            Self::SurfaceSurface(solver) => &solver.base.solver,
        }
    }

    /// Borrow the same base for a native member write.
    pub fn base_mut(&mut self) -> &mut FilletSolver {
        match self {
            Self::ConstantRadius(solver) => &mut solver.solver,
            Self::ConstantDistance(solver) => &mut solver.solver,
            Self::VariableRadius(solver) => &mut solver.curve_based.solver,
            Self::SurfaceSurface(solver) => &mut solver.base.solver,
        }
    }
}

/// Owned objects and the target borrowed by the executive.
#[derive(Debug)]
struct SessionObjects<'target> {
    target: &'target mut Brep,
    pseudo: FilletBrep,
    solvers: Vec<Option<SolverObject>>,
    corners: Vec<FilletCorner>,
    generators: Vec<FilletSurfaceGenerator>,
    attributes: AttributeStore,
    edge_caches: BTreeMap<EdgeId, CurveCache<FilletObject>>,
    target_serial: u64,
}

impl SessionObjects<'_> {
    fn solver(&self, id: FilletSolverId) -> Result<&SolverObject> {
        self.solvers
            .get(id.0 as usize)
            .and_then(Option::as_ref)
            .ok_or_else(|| invalid("solver is not live in this session"))
    }

    fn solver_mut(&mut self, id: FilletSolverId) -> Result<&mut SolverObject> {
        self.solvers
            .get_mut(id.0 as usize)
            .and_then(Option::as_mut)
            .ok_or_else(|| invalid("solver is not live in this session"))
    }
}

/// One actual executive and its live target, pseudo-B-rep, solvers and corners.
#[derive(Debug)]
pub struct FilletSession<'target> {
    executive: FilletExecutive<FilletObject>,
    objects: SessionObjects<'target>,
}

impl<'target> FilletSession<'target> {
    /// Construct the native executive-owned objects, with its ordinary defaults.
    ///
    /// # Errors
    /// Context B-rep construction, attribute allocation or identity exhaustion.
    pub fn new(k: &mut Kernel, target: &'target mut Brep) -> Result<Self> {
        let target_serial = crate::topology::brep::p3::notify_serial(&mut k.context, target)?;
        let mut pseudo = FilletBrep::new();
        // SmFilletBrep derives from SmBrep: use the actual context constructor,
        // including its context tolerances and construction notifications.
        *pseudo.brep_mut() = crate::topology::brep::p3::new_brep(
            k,
            crate::geom::tol::USE_DEFAULT,
            crate::geom::tol::USE_DEFAULT,
        )?;
        let mut attributes = AttributeStore::new();
        let attribute = attributes.insert(Attribute::new(
            ids::ALL_FILLET_EDGES,
            Behavior::StandaloneReference,
        ))?;
        let executive = FilletExecutive::new(
            Some(FilletObject::Body(FilletBody::Target)),
            ExecutiveObjects {
                executive: FilletObject::Executive,
                pseudo_brep: FilletObject::Body(FilletBody::Pseudo),
                edge_attribute: FilletObject::Attribute(attribute),
            },
        );
        Ok(Self {
            executive,
            objects: SessionObjects {
                target,
                pseudo,
                solvers: Vec::new(),
                corners: Vec::new(),
                generators: Vec::new(),
                attributes,
                edge_caches: BTreeMap::new(),
                target_serial,
            },
        })
    }

    /// The actual executive, for shared reads.
    #[must_use]
    pub const fn executive(&self) -> &FilletExecutive<FilletObject> {
        &self.executive
    }

    /// The real target (kept exclusively borrowed until the session is dropped).
    #[must_use]
    pub fn target(&self) -> &Brep {
        self.objects.target
    }

    /// The target entity's actual identity issued by the kernel context.
    /// This is the same canonical identity used by context-owned notifications.
    #[must_use]
    pub fn target_identity(&self, entity: EntityRef) -> crate::geom::object::ObjectRef {
        crate::topology::brep::p3::topology_ref(self.objects.target_serial, entity)
    }

    /// The actual pseudo-B-rep that owns the solver geometries.
    #[must_use]
    pub const fn pseudo_brep(&self) -> &FilletBrep {
        &self.objects.pseudo
    }

    /// The actual attribute store, including the executive's fillet-edge attribute.
    #[must_use]
    pub const fn attributes(&self) -> &AttributeStore {
        &self.objects.attributes
    }

    /// Read a retained concrete solver, without copying its base.
    ///
    /// # Errors
    /// An absent or failed-construction solver slot.
    pub fn solver(&self, id: FilletSolverId) -> Result<&SolverObject> {
        self.objects.solver(id)
    }

    /// Read an actual retained corner.
    ///
    /// # Errors
    /// A corner not issued by this session.
    pub fn corner(&self, id: FilletCornerId) -> Result<&FilletCorner> {
        self.objects
            .corners
            .get(id.0 as usize)
            .ok_or_else(|| invalid("corner is not live in this session"))
    }

    /// Retain a real cross-section generator, as the executive takes ownership.
    ///
    /// # Errors
    /// Generator identity exhaustion.
    pub fn add_surface_generator(
        &mut self,
        generator: FilletSurfaceGenerator,
    ) -> Result<FilletObject> {
        let id = u32::try_from(self.objects.generators.len())
            .map_err(|_| invalid("too many generators"))?;
        self.objects.generators.push(generator);
        let object = FilletObject::Generator(id);
        self.executive.add_fillet_surface_generator(object);
        Ok(object)
    }

    /// Construct `SmConstantRadiusFS` on the actual target sector and pseudo-B-rep.
    ///
    /// # Errors
    /// Solver identity exhaustion or the real native constructor's query failures.
    #[allow(clippy::too_many_arguments)] // the native constructor's parameters
    pub fn new_constant_radius_solver(
        &mut self,
        k: &Kernel,
        edgeuse: EdgeuseId,
        radius: f64,
        approx_tol: f64,
        ang_tol_rad: f64,
        tan_tol_rad: f64,
    ) -> Result<FilletSolverId> {
        let id = FilletSolverId(
            u32::try_from(self.objects.solvers.len()).map_err(|_| invalid("too many solvers"))?,
        );
        // Reserve the identity before construction writes a geometry pointing at it;
        // a later constructor failure must not let another solver reuse that id.
        self.objects.solvers.push(None);
        let solver = constant_radius_from_edgeuse(
            k,
            self.objects.target,
            &mut self.objects.pseudo,
            id,
            approx_tol,
            ang_tol_rad,
            tan_tol_rad,
            radius,
            edgeuse,
            None,
        )?;
        *self
            .objects
            .solvers
            .get_mut(id.0 as usize)
            .ok_or_else(|| invalid("lost solver slot"))? =
            Some(SolverObject::ConstantRadius(Box::new(solver)));
        Ok(id)
    }

    /// Assign a retained generator's descriptor to the actual solver base.
    /// This copies the port's immutable value descriptor; it does not promise native
    /// shared mutable generator-pointer identity.
    ///
    /// # Errors
    /// Unresolved solver or generator identities.
    pub fn set_solver_generator(
        &mut self,
        solver: FilletSolverId,
        generator: FilletObject,
    ) -> Result<()> {
        let FilletObject::Generator(id) = generator else {
            return Err(invalid("object is not a generator"));
        };
        let generator = self
            .objects
            .generators
            .get(id as usize)
            .ok_or_else(|| invalid("generator is not live"))?
            .clone();
        self.objects.solver_mut(solver)?.base_mut().fsg = Some(generator);
        Ok(())
    }

    /// Create and retain the actual corners in native executive order.
    ///
    /// # Errors
    /// A missing pending entry, undefined native null-vertex access, a failed factory
    /// query or exhausted corner identities. Earlier successful corners stay stored.
    pub fn create_corners(
        &mut self,
        k: &Kernel,
        error_info: Option<&mut ErrorInfo<FilletObject>>,
    ) -> Result<()> {
        let SessionObjects {
            target,
            solvers,
            corners,
            ..
        } = &mut self.objects;
        self.executive.create_fillet_corners(
            k,
            |k, executive, vertex, edgeuses| {
                let Some(FilletObject::Topology(FilletBody::Target, EntityRef::Vertex(vertex))) =
                    vertex
                else {
                    return Err(invalid(
                        "native corner factory requires an actual target vertex",
                    ));
                };
                let edgeuses = edgeuses
                    .iter()
                    .map(|object| match object {
                        FilletObject::Topology(FilletBody::Target, EntityRef::Edgeuse(id)) => {
                            Ok(*id)
                        }
                        _ => Err(invalid("corner edge use is not in the target")),
                    })
                    .collect::<Result<Vec<_>>>()?;
                let made = FilletCorner::create_with_solver_lookup(
                    k,
                    target,
                    vertex,
                    &edgeuses,
                    |edgeuse| {
                        let object =
                            FilletObject::Topology(FilletBody::Target, EntityRef::Edgeuse(edgeuse));
                        let Some(FilletObject::Solver(id)) =
                            executive.fillet_solver_of_edgeuse(object)
                        else {
                            return Err(invalid("corner edge use has no loaded solver"));
                        };
                        solvers
                            .get(id.0 as usize)
                            .and_then(Option::as_ref)
                            .map(SolverObject::base)
                            .ok_or_else(|| invalid("corner solver is not live"))
                    },
                )?;
                let Some(corner) = made else {
                    return Ok(None);
                };
                let id = FilletCornerId(
                    u32::try_from(corners.len()).map_err(|_| invalid("too many corners"))?,
                );
                corners.push(corner);
                Ok(Some(FilletObject::Corner(id)))
            },
            error_info,
        )
    }

    /// Load this actual solver through the native executive algorithm.
    ///
    /// # Errors
    /// The real loader's surface, cache, topology or lifecycle errors; prior writes remain.
    pub fn load_solver(&mut self, k: &Kernel, solver: FilletSolverId) -> Result<()> {
        self.executive
            .load_fillet_solver(k, &mut self.objects, FilletObject::Solver(solver))
    }
}

impl SolverLoadingAccess<FilletObject> for SessionObjects<'_> {
    fn load_solver_base_mut(&mut self, solver: FilletObject) -> Result<&mut FilletSolver> {
        let FilletObject::Solver(id) = solver else {
            return Err(invalid("object is not a solver"));
        };
        Ok(self.solver_mut(id)?.base_mut())
    }
    fn load_solver_has_surface(&self, solver: FilletObject, rail: usize) -> Result<bool> {
        let FilletObject::Solver(id) = solver else {
            return Err(invalid("object is not a solver"));
        };
        let base = self.solver(id)?.base();
        if let Some(geom) = base.fg_for_surfaces {
            Ok(self
                .pseudo
                .geom(geom)?
                .offset_surface(rail, Some(base))?
                .is_some())
        } else {
            Ok(base.surface(rail)?.is_some())
        }
    }
    fn load_edge_g1_knots(&mut self, k: &Kernel, edge: FilletObject) -> Result<Vec<f64>> {
        const LOADER: &str = "SmFilletExecutive::LoadFilletSolver";
        let FilletObject::Topology(FilletBody::Target, EntityRef::Edge(id)) = edge else {
            return Err(invalid("solver edge is not in the target"));
        };
        let curve = self
            .target
            .edge(id)?
            .curve
            .as_ref()
            .ok_or_else(|| SolidError::null_pointer(LOADER))?;
        let knots = curve.knots(k, None)?.0;
        if let std::collections::btree_map::Entry::Vacant(entry) = self.edge_caches.entry(id) {
            let cache = crate::curve::curve::p5::make_curve_cache(
                k,
                curve,
                Some(FilletObject::Curve(FilletBody::Target, id)),
            )
            .map_err(|error| {
                if error.is_unported() {
                    error
                } else {
                    SolidError::null_pointer(LOADER)
                }
            })?;
            entry.insert(cache);
        }
        self.edge_caches
            .get(&id)
            .ok_or_else(|| invalid("curve cache not retained"))?
            .continuities()
            .iter()
            .enumerate()
            .filter(|(_, continuity)| **continuity == ContinuityType::G1)
            .map(|(index, _)| {
                knots
                    .get(index)
                    .copied()
                    .ok_or_else(|| invalid("curve knot/cache count mismatch"))
            })
            .collect()
    }
    fn load_solver_edgeuse(
        &self,
        solver: FilletObject,
        rail: usize,
    ) -> Result<Option<FilletObject>> {
        let FilletObject::Solver(id) = solver else {
            return Err(invalid("object is not a solver"));
        };
        Ok(self
            .solver(id)?
            .base()
            .edgeuse(rail)?
            .map(|id| FilletObject::Topology(FilletBody::Target, EntityRef::Edgeuse(id))))
    }
    fn load_edgeuse_edge(&self, edgeuse: FilletObject) -> Result<Option<FilletObject>> {
        let FilletObject::Topology(FilletBody::Target, EntityRef::Edgeuse(id)) = edgeuse else {
            return Err(invalid("edge use is not in the target"));
        };
        Ok(self
            .target
            .edgeuse(id)?
            .edge
            .map(|id| FilletObject::Topology(FilletBody::Target, EntityRef::Edge(id))))
    }
    fn load_edgeuse_vertex(&self, edgeuse: FilletObject) -> Result<Option<FilletObject>> {
        let FilletObject::Topology(FilletBody::Target, EntityRef::Edgeuse(id)) = edgeuse else {
            return Err(invalid("edge use is not in the target"));
        };
        let vertexuse = self
            .target
            .edgeuse(id)?
            .vertexuse
            .ok_or_else(|| invalid("edge use has no vertex use"))?;
        Ok(self
            .target
            .vertexuse(vertexuse)?
            .vertex
            .map(|id| FilletObject::Topology(FilletBody::Target, EntityRef::Vertex(id))))
    }
}
