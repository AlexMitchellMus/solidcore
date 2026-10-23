//! Caller-owned storage for fillet executives and their solver back references.
//!
//! Upstream: `source/SMLib/inc/SmFilletExecutive.h`, `SmFilletSolver.h`.
//!
//! | Rust | Upstream |
//! |------|----------|
//! | [`FilletExecutiveStore`] | the context-owned `SmFilletExecutive` objects |
//! | [`RegisteredExecutiveRef`] | a resolved `SmFilletExecutive*` |
//!
//! Deviations: a store-local index and an identity token replace an object address.
//! The token owns no model data. Copies of an executive lose their registration.

use std::sync::Arc;

use super::p3::FilletExecutive;
use crate::{Result, SolidError};

const OP: &str = "blend::fillet_executive::registration";

/// A solver's reference to an executive retained in a particular store.
/// Only [`FilletExecutiveStore::insert`] issues references.
#[derive(Debug, Clone)]
pub struct RegisteredExecutiveRef {
    store: Arc<()>,
    index: usize,
}

impl PartialEq for RegisteredExecutiveRef {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.store, &other.store) && self.index == other.index
    }
}
impl Eq for RegisteredExecutiveRef {}

/// The executive's own registration. Cloning the executive creates a new object.
#[derive(Debug, Default)]
pub struct ExecutiveRegistration(Option<RegisteredExecutiveRef>);

impl Clone for ExecutiveRegistration {
    fn clone(&self) -> Self {
        Self::default()
    }
}

impl ExecutiveRegistration {
    /// The reference issued when the executive entered its actual store.
    #[must_use]
    pub fn reference(&self) -> Option<RegisteredExecutiveRef> {
        self.0.clone()
    }
}

/// Actual executive values retained by their caller. References from another store
/// are rejected even if both stores allocated the same local index.
#[derive(Debug)]
pub struct FilletExecutiveStore<O> {
    identity: Arc<()>,
    executives: Vec<FilletExecutive<O>>,
}

impl<O> Default for FilletExecutiveStore<O> {
    fn default() -> Self {
        Self {
            identity: Arc::new(()),
            executives: Vec::new(),
        }
    }
}

impl<O> FilletExecutiveStore<O> {
    /// An empty store.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Retain an actual executive and issue its solver back reference.
    pub fn insert(&mut self, mut executive: FilletExecutive<O>) -> RegisteredExecutiveRef {
        let reference = RegisteredExecutiveRef {
            store: Arc::clone(&self.identity),
            index: self.executives.len(),
        };
        executive.registration = ExecutiveRegistration(Some(reference.clone()));
        self.executives.push(executive);
        reference
    }

    /// Resolve a reference to the actual retained executive.
    ///
    /// # Errors
    /// `InvalidInput` for another store's reference or a stale registration.
    pub fn resolve(&self, reference: &RegisteredExecutiveRef) -> Result<&FilletExecutive<O>> {
        self.check(reference)?;
        self.executives
            .get(reference.index)
            .ok_or_else(|| SolidError::invalid_input(OP, "executive reference is not retained"))
    }

    /// Resolve a reference for a mutation of the actual executive.
    ///
    /// # Errors
    /// As [`Self::resolve`].
    pub fn resolve_mut(
        &mut self,
        reference: &RegisteredExecutiveRef,
    ) -> Result<&mut FilletExecutive<O>> {
        self.check(reference)?;
        self.executives
            .get_mut(reference.index)
            .ok_or_else(|| SolidError::invalid_input(OP, "executive reference is not retained"))
    }

    /// Resolve a solver's back pointer to the actual executive it names.
    ///
    /// # Errors
    /// As [`Self::resolve`]; `Unported` for legacy opaque caller state, which never
    /// constitutes a registered object binding.
    pub fn resolve_solver(
        &self,
        solver: &crate::blend::fillet_solver::p1::FilletSolver,
    ) -> Result<Option<&FilletExecutive<O>>> {
        use crate::blend::fillet_solver::p1::FilletExecutiveRef;
        match solver.executive.as_ref() {
            None => Ok(None),
            Some(FilletExecutiveRef::Registered(reference)) => self.resolve(reference).map(Some),
            Some(FilletExecutiveRef::Unresolved) => Err(SolidError::unported(
                "registered fillet executive",
                "blend.fillet_executive_p1",
            )),
        }
    }

    fn check(&self, reference: &RegisteredExecutiveRef) -> Result<()> {
        if !Arc::ptr_eq(&self.identity, &reference.store) {
            return Err(SolidError::invalid_input(
                OP,
                "executive reference belongs to another store",
            ));
        }
        let stored = self
            .executives
            .get(reference.index)
            .and_then(|executive| executive.registration.0.as_ref());
        if stored != Some(reference) {
            return Err(SolidError::invalid_input(
                OP,
                "executive registration is stale",
            ));
        }
        Ok(())
    }
}
