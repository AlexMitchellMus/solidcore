//! Borrowed offset context for the actual constant-radius rail solver.

use crate::blend::fillet_custom_solver::{
    ConstantRadiusAssistedFs, InitialValues, JacobianLoad, OffsetSurface,
};
use crate::blend::fillet_geom::p1::{FilletGeom, FilletGeomType};
use crate::blend::fillet_standard_solver::p1::constant_radius_load_jacobian;
use crate::blend::fillet_standard_solver::p2::{
    constant_radius_load_initial_values, ConstantRadiusOps,
};
use crate::error::{Result, SolidError};
use crate::geom::extent1d::ExtentNd;
use crate::geom::extent2d::Extent2d;
use crate::geom::matrix::Matrix;
use crate::intersect::surface_intersector::p1::TsectPnt;
use crate::kernel::Kernel;
use crate::surface::surface::Surface;
use crate::topology::{Brep, EdgeuseId};

use super::p1::{FilletSolver, FilletSolverOps};

struct GeometrySolver<'a> {
    actual: &'a mut ConstantRadiusAssistedFs,
    geometry: &'a mut FilletGeom,
}

impl ConstantRadiusOps for GeometrySolver<'_> {
    fn constant_radius(&self) -> &ConstantRadiusAssistedFs {
        self.actual
    }

    fn constant_radius_mut(&mut self) -> &mut ConstantRadiusAssistedFs {
        self.actual
    }

    fn offset_surface(&self, rail: usize) -> Result<&OffsetSurface> {
        self.geometry
            .offset_surface(rail, Some(&self.actual.solver))?
            .ok_or_else(|| SolidError::null_pointer("blend::fillet_solver::surfaces"))
    }

    fn offset_surface_mut(&mut self, rail: usize) -> Result<&mut OffsetSurface> {
        let own = self.geometry.off_surfs.get_mut(rail).ok_or_else(|| {
            SolidError::invalid_input("blend::fillet_solver::surfaces", "rail index past 1")
        })?;
        if let Some(surface) = own {
            return Ok(surface);
        }
        self.actual.surface_mut(rail)
    }

    fn set_offset_by_radius(&mut self, rail: usize) -> Result<()> {
        let distance = self.actual.fillet_radius(rail)? * self.actual.solver.orientation(rail)?;
        ConstantRadiusOps::offset_surface_mut(self, rail)?.set_offset_distance(distance);
        Ok(())
    }

    fn check_tangent_rollover_guess_context(&self) -> Result<()> {
        if self.geometry.fillet_geom_type() == FilletGeomType::TangentRollover {
            // UNPORTED(blend.fillet_standard_solver_p2): SmConstantRadiusFS::LoadInitialValues tangent-rollover vertex-use context
            return Err(SolidError::unported(
                "SmConstantRadiusFS::LoadInitialValues tangent-rollover vertex-use context",
                "blend.fillet_standard_solver_p2",
            ));
        }
        Ok(())
    }

    fn setup_offset_extension(
        &self,
        k: &Kernel,
        surface: &Surface,
        domain: &mut Extent2d,
    ) -> Result<()> {
        FilletSolverOps::setup_offset_extension(self.actual, k, surface, domain)
    }
}

impl FilletSolverOps for GeometrySolver<'_> {
    fn solver(&self) -> &FilletSolver {
        &self.actual.solver
    }
    fn solver_mut(&mut self) -> &mut FilletSolver {
        &mut self.actual.solver
    }

    fn solver_surface(&self, rail: usize) -> Result<Option<&OffsetSurface>> {
        Ok(Some(ConstantRadiusOps::offset_surface(self, rail)?))
    }

    fn solver_surface_mut(&mut self, rail: usize) -> Result<Option<&mut OffsetSurface>> {
        Ok(Some(ConstantRadiusOps::offset_surface_mut(self, rail)?))
    }

    fn jacobian_size(&self) -> Result<usize> {
        FilletSolverOps::jacobian_size(self.actual)
    }

    fn setup_offset_extension(
        &self,
        k: &Kernel,
        surface: &Surface,
        domain: &mut Extent2d,
    ) -> Result<()> {
        FilletSolverOps::setup_offset_extension(self.actual, k, surface, domain)
    }

    fn load_initial_values(
        &mut self,
        k: &Kernel,
        brep: Option<&Brep>,
        do_surf2_calcs: bool,
        rail_index: usize,
        surface1: &mut OffsetSurface,
        surface1_index: usize,
        intervals: &mut ExtentNd,
        periodicities: &mut Vec<bool>,
        guess_t: &mut Vec<f64>,
    ) -> Result<Option<InitialValues>> {
        constant_radius_load_initial_values(
            self,
            k,
            brep,
            do_surf2_calcs,
            rail_index,
            surface1,
            surface1_index,
            intervals,
            periodicities,
            guess_t,
        )
    }

    fn load_jacobian(
        &mut self,
        k: &Kernel,
        num_equations: &mut usize,
        num_parameters: &mut usize,
        load: &mut JacobianLoad<'_>,
        f: &mut [f64],
        jacobian: Option<&mut Matrix>,
    ) -> Result<bool> {
        let mut offsets = [load.surf1_offset, load.surf2_offset];
        let result = constant_radius_load_jacobian(
            self,
            k,
            num_equations,
            num_parameters,
            load.rail_index,
            [load.surface1, load.surface2],
            &mut offsets,
            load.x,
            f,
            jacobian,
        );
        [load.surf1_offset, load.surf2_offset] = offsets;
        result
    }
}

/// Run with one actual solver and its selected geometry; no offset object moves.
#[allow(clippy::too_many_arguments)]
pub(crate) fn rail_edgeuse_intersect(
    actual: &mut ConstantRadiusAssistedFs,
    k: &Kernel,
    geometry: &mut FilletGeom,
    target: &mut Brep,
    rail_index: usize,
    edgeuse: EdgeuseId,
    guess: f64,
    tsect: &mut TsectPnt,
) -> Result<Option<f64>> {
    FilletSolver::rail_edgeuse_intersect(
        &mut GeometrySolver { actual, geometry },
        k,
        target,
        rail_index,
        edgeuse,
        guess,
        tsect,
    )
}
