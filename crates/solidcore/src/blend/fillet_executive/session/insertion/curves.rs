//! Actual operation-owned curve copies. Handles never impersonate topology.

use crate::curve::curve::Curve;
use crate::geom::context::Context;
use crate::geom::object::ObjectRef;
use crate::{Result, SolidError};

/// An issued identity and its retained operation-owned Curve slot.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct CurveBinding {
    slot: usize,
    identity: ObjectRef,
}

/// Curves retained for the whole insertion operation, including partial failures.
#[derive(Debug, Default)]
pub struct InsertionCurves {
    curves: Vec<Curve>,
}
impl InsertionCurves {
    /// Move an actual curve copy into operation ownership and issue its identity.
    /// # Errors
    /// Identity-domain issuance failure; the curve is retained on this failure too.
    pub fn insert(&mut self, context: &mut Context, curve: Curve) -> Result<CurveBinding> {
        let slot = self.curves.len();
        self.curves.push(curve);
        let curve = self.curves.get_mut(slot).ok_or_else(|| {
            SolidError::invalid_input(
                "blend::fillet_executive::InsertionCurves",
                "inserted curve slot is absent",
            )
        })?;
        let type_code = curve.type_code();
        let key = context.geometry_key(curve.notify_identity_slot())?;
        Ok(CurveBinding {
            slot,
            identity: ObjectRef { type_code, key },
        })
    }
    /// Borrow the actual retained object, rejecting a foreign/replaced binding.
    #[must_use]
    pub fn get(&self, binding: CurveBinding) -> Option<&Curve> {
        let curve = self.curves.get(binding.slot)?;
        (curve.notify_key() == Some(binding.identity.key)
            && curve.type_code() == binding.identity.type_code)
            .then_some(curve)
    }
    /// Mutably borrow the same retained object after identity validation.
    pub fn get_mut(&mut self, binding: CurveBinding) -> Option<&mut Curve> {
        self.get(binding)?;
        self.curves.get_mut(binding.slot)
    }
    /// Resolve a binding with an explicit failure at a provider capability boundary.
    /// # Errors
    /// A binding from a different store or a replaced object.
    pub fn checked(&self, binding: CurveBinding) -> Result<&Curve> {
        self.get(binding).ok_or_else(|| {
            SolidError::invalid_input(
                "blend::fillet_executive::InsertionCurves",
                "foreign or replaced curve binding",
            )
        })
    }
}
