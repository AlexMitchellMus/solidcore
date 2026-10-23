//! FFI wrappers for unit `geom.grid`.
//!
//! Shim: `harness/cpp/shim/geom_grid.cpp`.
//! Upstream: `source/SMLib/src/SmGrid.cpp` (`SmGrid`, `SmVoxel`, `SmGridElement`,
//! `SmHierarchyElement`).
//!
//! | Rust      | shim                                       |
//! |-----------|--------------------------------------------|
//! | [`run`]   | `sc_geom_grid_ints`, `sc_geom_grid_floats` |
//!
//! Runs a small program against one upstream grid; the encoding is documented
//! in the shim and validated here before any call.
//!
//! | op   | operands                        | meaning                                             |
//! |------|---------------------------------|-----------------------------------------------------|
//! | 0    | `e6 nx ny nz`                   | `SmGrid(bbox, lSize)` (first)                       |
//! | 1    | `type obj e6`                   | `AddElement` (poly B-rep or user data)              |
//! | 2    | `kind surf tol e6`              | `AddHierarchyElement` (base, or a *fake* surface)   |
//! | 3-8  | as the shim                     | queries, 3DDDA, dump                                |
//! | 9    | `h`                             | `ExpandHierarchy` (base or real-surface element)    |
//! | 10   | `s node e6`                     | `AddElement(SM_GE_SURFACE)` of tree node `node`     |
//! | 11   | `0 p q cpts(x y z w)` or `1 o x y scale dom` | declare a surface (variable length)  |
//! | 12   | `s tol e6`                      | `AddHierarchyElement` of a real surface element     |
//!
//! A surface is a clamped single-patch B-spline of degree 1..=3 (finite homogeneous
//! coordinates and Euclidean points `wx / w` of at most [`MAX_SURFACE_COORD`], positive
//! weights) or a plane (independent axes, positive scale and domain). The
//! shim never dereferences the fake surface of op 2, so `expand` of such an element
//! (op 9) and an out-of-range surface index are rejected here.

use crate::validate::bad;
use crate::{call_with_buffer, FfiError};

mod sys {
    extern "C" {
        pub fn sc_geom_grid_ints(
            prog: *const f64,
            n: i64,
            out: *mut i64,
            cap: i64,
            out_len: *mut i64,
        ) -> i32;
        pub fn sc_geom_grid_floats(
            prog: *const f64,
            n: i64,
            out: *mut f64,
            cap: i64,
            out_len: *mut i64,
        ) -> i32;
    }
}

/// Largest object id + 1 (size of the shim's fake object array).
pub const MAX_OBJECTS: i64 = 4096;

/// Largest voxel count per axis accepted.
pub const MAX_AXIS: i64 = 16;

/// Largest absolute coordinate accepted (keeps every box size and voxel
/// address computation finite).
pub const MAX_COORD: f64 = 1e12;

/// Largest number of declared surfaces.
pub const MAX_SURFACES: usize = 4;

/// Largest absolute surface coordinate and weight accepted.
pub const MAX_SURFACE_COORD: f64 = 1e9;

/// Number of operands of op `op`, whose operands start at `rest`; `None` for an
/// unknown op or a malformed variable-length one (op 11).
#[must_use]
pub fn op_size(op: i64, rest: &[f64]) -> Option<usize> {
    Some(match op {
        0 | 2 => 9,
        1 | 10 | 12 => 8,
        3 | 6 => 6,
        4 | 5 => 3,
        7 | 9 => 1,
        8 => 0,
        11 => match rest.first() {
            Some(&1.0) => 16,
            Some(&0.0) => {
                let p = rest.get(1).copied()?;
                let q = rest.get(2).copied()?;
                if !int_in(p, 1, 3) || !int_in(q, 1, 3) {
                    return None;
                }
                3 + 4 * (p as usize + 1) * (q as usize + 1)
            }
            _ => return None,
        },
        _ => return None,
    })
}

fn norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

fn cross(a: &[f64], b: &[f64]) -> [f64; 3] {
    let g = |s: &[f64], i: usize| s.get(i).copied().unwrap_or(f64::NAN);
    [
        g(a, 1) * g(b, 2) - g(a, 2) * g(b, 1),
        g(a, 2) * g(b, 0) - g(a, 0) * g(b, 2),
        g(a, 0) * g(b, 1) - g(a, 1) * g(b, 0),
    ]
}

/// The operands of op 11: a surface the shim can build and the cache can decompose.
fn surface_ok(a: &[f64]) -> bool {
    let surf_coords = |s: &[f64]| s.iter().all(|v| v.is_finite() && v.abs() <= MAX_SURFACE_COORD);
    match a.first().copied() {
        Some(0.0) => {
            let Some(net) = a.get(3..) else { return false };
            let first_w = net.get(3).copied().unwrap_or(f64::NAN);
            // weights: all one (non-rational), or all positive (rational)
            let rational = first_w != 1.0;
            !net.is_empty()
                && net.as_chunks::<4>().0.iter().all(|c| {
                    // the points are homogeneous `(wx, wy, wz, w)`: the Euclidean point
                    // `wx / w` is what the surface evaluates, and must stay small too
                    let w = c.get(3).copied().unwrap_or(f64::NAN);
                    let weight_ok = if rational {
                        w > 0.0 && w <= MAX_SURFACE_COORD
                    } else {
                        w == 1.0
                    };
                    let euclidean: Vec<f64> = c.iter().take(3).map(|x| x / w).collect();
                    weight_ok
                        && surf_coords(c.get(..3).unwrap_or(&[f64::NAN]))
                        && surf_coords(&euclidean)
                })
        }
        Some(1.0) => {
            let (Some(o), Some(x), Some(y), Some(s), Some(d)) = (
                a.get(1..4),
                a.get(4..7),
                a.get(7..10),
                a.get(10..12),
                a.get(12..16),
            ) else {
                return false;
            };
            let c = norm(&cross(x, y));
            surf_coords(o)
                && surf_coords(x)
                && surf_coords(y)
                && surf_coords(s)
                && surf_coords(d)
                && c > 1e-6 * norm(x) * norm(y)
                && norm(x) >= 1e-9
                && norm(y) >= 1e-9
                && s.iter().all(|&v| v > 1e-9)
                && d.first() < d.get(2)
                && d.get(1) < d.get(3)
        }
        _ => false,
    }
}

fn int_in(v: f64, lo: i64, hi: i64) -> bool {
    v.fract() == 0.0 && v >= lo as f64 && v <= hi as f64
}

fn coords_ok(a: &[f64]) -> bool {
    a.iter().all(|v| v.is_finite() && v.abs() <= MAX_COORD)
}

fn validate(prog: &[f64]) -> Result<(), FfiError> {
    let ctor = prog.get(0..10).ok_or_else(|| bad("missing constructor"))?;
    if ctor.first() != Some(&0.0) {
        return Err(bad("constructor must be first"));
    }
    let bbox = ctor.get(1..7).ok_or_else(|| bad("constructor"))?;
    if !coords_ok(bbox) {
        return Err(bad("grid box"));
    }
    for k in 0..3 {
        let (lo, hi) = (bbox.get(k).copied(), bbox.get(k + 3).copied());
        match (lo, hi) {
            // A grid with a zero-size (or inverted) axis makes GetVoxelAddress
            // convert NaN or a negative value to an index: undefined behaviour.
            (Some(l), Some(h)) if h - l > 0.0 && (h - l) * 1e-12 > 0.0 => {}
            _ => return Err(bad("grid box must have positive size")),
        }
    }
    for k in 7..10 {
        if !int_in(ctor.get(k).copied().unwrap_or(0.0), 1, MAX_AXIS) {
            return Err(bad("voxel counts"));
        }
    }
    // per hierarchy element: whether `ExpandHierarchy` may run on it (a base-class
    // element, or a real-surface one; the fake surface of op 2 would be dereferenced)
    let mut expandable: Vec<bool> = Vec::new();
    let mut surfaces = 0usize;
    let mut i = 10usize;
    while i < prog.len() {
        let op = prog
            .get(i)
            .copied()
            .ok_or_else(|| bad("truncated program"))?;
        if !int_in(op, 1, 12) {
            return Err(bad("op"));
        }
        let op = op as i64;
        let rest = prog.get(i + 1..).unwrap_or(&[]);
        let sz = op_size(op, rest).ok_or_else(|| bad("op"))?;
        let a = prog
            .get(i + 1..i + 1 + sz)
            .ok_or_else(|| bad("truncated program"))?;
        let at = |j: usize| a.get(j).copied().unwrap_or(f64::NAN);
        let ok = match op {
            1 => {
                (at(0) == 1.0 || at(0) == 3.0)
                    && int_in(at(1), 0, MAX_OBJECTS - 1)
                    && coords_ok(a.get(2..8).unwrap_or(&[f64::NAN]))
            }
            2 => {
                expandable.push(at(0) == 0.0);
                (at(0) == 0.0 || at(0) == 1.0)
                    && int_in(at(1), 0, MAX_OBJECTS - 1)
                    && coords_ok(a.get(2..9).unwrap_or(&[f64::NAN]))
            }
            3 | 4 | 6 => coords_ok(a),
            5 => a.iter().all(|&v| int_in(v, -1000, 1000)),
            7 => !at(0).is_nan(),
            9 => {
                int_in(at(0), 0, expandable.len() as i64 - 1)
                    && expandable.get(at(0) as usize).copied() == Some(true)
            }
            // a node index outside the surface's tree is answered by the shim (-4)
            10 => {
                int_in(at(0), 0, surfaces as i64 - 1)
                    && at(1).is_finite()
                    && coords_ok(a.get(2..8).unwrap_or(&[f64::NAN]))
            }
            11 => {
                surfaces += 1;
                surfaces <= MAX_SURFACES && surface_ok(a)
            }
            12 => {
                expandable.push(true);
                int_in(at(0), 0, surfaces as i64 - 1)
                    && coords_ok(a.get(1..8).unwrap_or(&[f64::NAN]))
            }
            _ => true,
        };
        if !ok {
            return Err(bad("operands"));
        }
        i += 1 + sz;
    }
    Ok(())
}

/// Runs a grid program; returns its integer and float output streams.
pub fn run(prog: &[f64]) -> Result<(Vec<i64>, Vec<f64>), FfiError> {
    validate(prog)?;
    let n = prog.len() as i64;
    let ints = call_with_buffer(256, |buf: &mut [i64], len| {
        // SAFETY: prog valid for n reads; buf valid for buf.len() writes; len valid.
        unsafe {
            sys::sc_geom_grid_ints(
                prog.as_ptr(),
                n,
                buf.as_mut_ptr(),
                buf.len() as i64,
                &raw mut *len,
            )
        }
    })?;
    let floats = call_with_buffer(256, |buf: &mut [f64], len| {
        // SAFETY: as above.
        unsafe {
            sys::sc_geom_grid_floats(
                prog.as_ptr(),
                n,
                buf.as_mut_ptr(),
                buf.len() as i64,
                &raw mut *len,
            )
        }
    })?;
    Ok((ints, floats))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A grid over the unit cube followed by `ops`.
    fn program(ops: &[f64]) -> Vec<f64> {
        let mut p = vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0];
        p.extend(ops);
        p
    }

    /// Op 11, a bilinear patch from four homogeneous points `(wx, wy, wz, w)`.
    fn bilinear(points: [[f64; 4]; 4]) -> Vec<f64> {
        let mut op = vec![11.0, 0.0, 1.0, 1.0];
        op.extend(points.iter().flatten());
        op
    }

    #[test]
    fn a_well_formed_surface_is_accepted() {
        let unit = [
            [0.0, 0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0, 1.0],
            [1.0, 0.0, 0.0, 1.0],
            [1.0, 1.0, 0.0, 1.0],
        ];
        assert!(validate(&program(&bilinear(unit))).is_ok());
    }

    /// The control points are homogeneous: `wx` and `w` are bounded one by one, but
    /// the point `wx / w` the shim's surface evaluates must stay finite and small.
    #[test]
    fn a_surface_whose_euclidean_point_overflows_is_rejected() {
        let tiny_weight = [
            [1e9, 0.0, 0.0, 1e-300],
            [0.0, 1.0, 0.0, 1.0],
            [1.0, 0.0, 0.0, 1.0],
            [1.0, 1.0, 0.0, 1.0],
        ];
        assert!(validate(&program(&bilinear(tiny_weight))).is_err());
        let large = [
            [1e9, 0.0, 0.0, 1e-3],
            [0.0, 1.0, 0.0, 1.0],
            [1.0, 0.0, 0.0, 1.0],
            [1.0, 1.0, 0.0, 1.0],
        ];
        assert!(validate(&program(&bilinear(large))).is_err());
    }
}
