//! Checked native `primitive.primitive_creation_p2` fixtures: the public
//! `SmPrimitiveCreation` entry points on a fresh native B-rep.
//!
//! Shim: `harness/cpp/shim/primitive_primitive_creation_p2.cpp`. Curves are p1's
//! `[kind, x, y, z, a, b, c]` records (kind 0 a B-spline line, kind 1 a quadratic
//! circular arc centred on xyz, radius a, angles b..c degrees). Frames are
//! `origin, x axis, y axis` (nine values). Layouts (`edit` is the B-rep's initial
//! editing flag):
//!
//! | op | entry point | data |
//! |---:|---|---|
//! | 1 | partial `CreateCone` | `edit h base top start end frame` |
//! | 2 | partial `CreateSphere` | `edit r start end frame` |
//! | 3 | partial `CreateTorus` | `edit major minor start end frame` |
//! | 4 | `CreateSkinPrimitive` (ownership report) | `edit approx degree cap profiles counts... curves...` |
//! | 5 | `CreateCurveSweepFromEdges`, no edges | `edit approx cap translational center perpendicular to_profile path` |
//! | 6 | `CreateCurveSweepFromFaces`, no faces | as op 5 |
//! | 7 | `CreatePipeSweep` | `edit radius approx cap translational path` |
//! | 8 | `CreateTaperExtrude` | `edit height draft caps nx ny nz corner curves...` |
//! | 9 | static `CreateSweepAlongPlanarPath` | `tol cap_start cap_end move profiles paths curves...` |
//!
//! Inputs upstream reads undefined memory for are rejected here: a taper with fewer
//! than two curves (UB-111 for none; one curve can make the sliver test read past
//! the split curve array, the dereference behind UB-112), a curve shorter than
//! `2e-5` (UB-112) or curves its chaining deletes as duplicates and reads again
//! (UB-50, both chaining passes; with the end cap at the offset end, any chain that
//! encloses nothing), a planar path without path curves (`rPathCurves[0]`), skin
//! profile counts past the curves (`r3DCurves[lCount]`), a partial cone or torus that
//! sweeps a profile edge shorter than the zone tolerance (UB-215, use after free), and
//! a sphere, torus or cone a few zone tolerances across, at any turn (HANG-35).
//! Values are finite with magnitude at most `1e6`, frames nondegenerate, lists at
//! most 32 curves.

#![allow(clippy::indexing_slicing)] // fixed fixture layout validated before indexing

use crate::{call_with_buffer, FfiError};

mod sys {
    extern "C" {
        pub fn sc_primitive_primitive_creation_p2_eval(
            op: i32,
            data: *const f64,
            n: i64,
            output: *mut f64,
            cap: i64,
            length: *mut i64,
        ) -> i32;
    }
}

const MAX_CURVES: usize = 32;

fn bad() -> FfiError {
    FfiError::InvalidArgument("primitive p2 layout or geometry".into())
}
fn integer(v: f64, max: f64) -> bool {
    v.is_finite() && v.fract() == 0.0 && (0.0..=max).contains(&v)
}
fn flag(v: f64) -> bool {
    integer(v, 1.0)
}
/// A buildable curve record (as p1's fixtures).
fn curve(v: &[f64]) -> bool {
    let [kind, x, y, z, a, b, c] = v else {
        return false;
    };
    if *kind == 0.0 {
        (a - x).powi(2) + (b - y).powi(2) + (c - z).powi(2) > 0.0
    } else {
        *kind == 1.0
            && *a >= 1e-15
            && (-360.0..=360.0).contains(b)
            && (-360.0..=360.0).contains(c)
            && c > b
            && c - b <= 360.0
    }
}
/// The extent of a curve record: its chord for a line, its radius for an arc.
fn extent(v: &[f64]) -> f64 {
    match v {
        [kind, x, y, z, a, b, c] if *kind == 0.0 => {
            ((a - x).powi(2) + (b - y).powi(2) + (c - z).powi(2)).sqrt()
        }
        [_, _, _, _, a, ..] => *a,
        _ => 0.0,
    }
}
fn curves(data: &[f64]) -> Option<Vec<&[f64]>> {
    let (records, rest) = data.as_chunks::<7>();
    if !rest.is_empty() || records.len() > MAX_CURVES {
        return None;
    }
    let records: Vec<&[f64]> = records.iter().map(<[f64; 7]>::as_slice).collect();
    records.iter().all(|r| curve(r)).then_some(records)
}
fn frame(f: &[f64]) -> bool {
    let [_, _, _, ax, ay, az, bx, by, bz] = f else {
        return false;
    };
    let cross = [ay * bz - az * by, az * bx - ax * bz, ax * by - ay * bx];
    cross.iter().map(|v| v * v).sum::<f64>() > 1e-12
}

/// `SM_EFF_ZERO`.
const EFF_ZERO: f64 = 1e-12;
/// `SM_ZONE_TOL_3D`, the zone tolerance of the fresh native B-rep.
const ZONE_TOL: f64 = 1e-5;
/// From this turn (degrees) up, p3's `CreateRotationalSweep` builds analytic faces
/// instead of sweeping the profile edges; p2's full turns are among them.
const ANALYTIC_TURN: f64 = 359.99999;
/// HANG-35 ends here: every probe from this size up returned in milliseconds, while a
/// sphere of `1e-4` over 359 degrees took 20 s (fuzz allows 20 s per case).
const HANG_MAX: f64 = 1.5e-4;

type V3 = [f64; 3];
fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn along(a: V3, s: f64, b: V3) -> V3 {
    [a[0] + s * b[0], a[1] + s * b[1], a[2] + s * b[2]]
}
fn cross(a: V3, b: V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn norm(a: V3) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}
fn unit(a: V3) -> V3 {
    let n = norm(a);
    [a[0] / n, a[1] / n, a[2] / n]
}
fn max_dim(a: V3) -> f64 {
    a[0].abs().max(a[1].abs()).max(a[2].abs())
}

/// The partial cone's profile lines as upstream builds them (side, base radius,
/// axis, top radius; a radius within the scaled zero is left out), or `None` when
/// the cone fails before sweeping them. The frame is checked nondegenerate first.
fn cone_lines(d: &[f64]) -> Option<Vec<f64>> {
    let (h, base, top, start, f) = (d[1], d[2], d[3], d[4], &d[6..15]);
    if h < EFF_ZERO || (d[5] - start).abs() < EFF_ZERO {
        return None;
    }
    let origin = [f[0], f[1], f[2]];
    let x = unit([f[3], f[4], f[5]]);
    let z = unit(cross(x, [f[6], f[7], f[8]]));
    let (sin, cos) = start.to_radians().sin_cos();
    let x_start = along([0.0; 3], cos, x);
    let x_start = along(x_start, sin, cross(z, x));
    let top_base = along(origin, h, z);
    let start_pt = along(origin, base, x_start);
    let end_pt = along(top_base, top, x_start);
    let tol = EFF_ZERO * (1.0 + max_dim(end_pt) + max_dim(start_pt));
    let (side, axis) = (norm(sub(end_pt, start_pt)), norm(sub(top_base, origin)));
    if side <= tol || axis <= tol {
        return None;
    }
    let lines: Vec<f64> = [
        side,
        norm(sub(start_pt, origin)),
        axis,
        norm(sub(top_base, end_pt)),
    ]
    .into_iter()
    .filter(|l| *l > tol)
    .collect();
    (lines.len() >= 3).then_some(lines)
}

/// UB-215: a turn below [`ANALYTIC_TURN`] sweeps the profile edges, and
/// `SweepEdgeHigh` uses an edge's curve after the face build freed it for being
/// shorter than the zone tolerance. A sphere's line lies on the axis and its arc's
/// face fails the degenerate-surface test first, so spheres never reach it.
fn sweeps_a_short_edge(op: i32, d: &[f64]) -> bool {
    match op {
        1 => {
            d[5] - d[4] < ANALYTIC_TURN
                && cone_lines(d).is_some_and(|l| l.iter().any(|l| *l < ZONE_TOL))
        }
        3 => d[4] - d[3] < ANALYTIC_TURN && d[2] > 0.0 && std::f64::consts::TAU * d[2] < ZONE_TOL,
        _ => false,
    }
}

/// HANG-35: a closed shell a few zone tolerances across takes minutes to place in
/// its region (the ray fire of `FindInfiniteRegion`), whether it is an analytic face
/// or a stitched sweep, at any turn. Cones that fail first (no height or turn, a
/// partial profile that is not built) and zero radii are kept.
fn builds_a_tiny_closed_shell(op: i32, d: &[f64]) -> bool {
    let size = match op {
        1 => {
            let turn = d[5] - d[4];
            let full = (turn - 360.0).abs() < 360.0 * EFF_ZERO;
            if d[1] < EFF_ZERO || turn.abs() < EFF_ZERO || !(full || cone_lines(d).is_some()) {
                return false;
            }
            d[1].max(d[2].abs()).max(d[3].abs())
        }
        2 if d[1] != 0.0 => d[1].abs(),
        3 if d[1] > 0.0 && d[2] > 0.0 => d[1] + d[2],
        _ => return false,
    };
    (ZONE_TOL / 2.0..HANG_MAX).contains(&size)
}

fn dist2(a: V3, b: V3) -> f64 {
    let d = sub(a, b);
    d[0] * d[0] + d[1] * d[1] + d[2] * d[2]
}

/// A curve record's start, end and mid-parameter points (an arc's mid parameter is
/// its mid angle).
fn points(v: &[f64]) -> [V3; 3] {
    match v {
        [kind, x, y, z, a, b, c] if *kind == 0.0 => {
            let (start, end) = ([*x, *y, *z], [*a, *b, *c]);
            [start, end, along(start, 0.5, sub(end, start))]
        }
        [_, x, y, z, r, b, c] => {
            let at = |deg: f64| {
                let (sin, cos) = deg.to_radians().sin_cos();
                [x + r * cos, y + r * sin, *z]
            };
            [at(*b), at(*c), at((b + c) / 2.0)]
        }
        _ => [[0.0; 3]; 3],
    }
}

/// One `SmCurve::FixChainCurvesList` pass as the taper calls it (`dTol = 2 * ZONE_TOL`,
/// deleting duplicates), on the curves' start, end and mid points: the chains in
/// upstream's order, their curves in its order and direction, or `None` when it
/// deletes one (UB-50). A curve is a duplicate when it matches the ends of the chain
/// built so far, not another curve: a square with a repeated side deletes nothing,
/// while a line from a chain's end back to its start, through its last link's mid
/// point, is deleted.
fn chain(curves: Vec<[V3; 3]>) -> Option<Vec<Vec<[V3; 3]>>> {
    let tol2 = (2.0 * ZONE_TOL).powi(2);
    let mut left: Vec<Option<[V3; 3]>> = curves.into_iter().map(Some).collect();
    let mut chains = Vec::new();
    for kk in 0..left.len() {
        let Some(first) = left[kk].take() else {
            continue;
        };
        let mut chained = vec![first];
        let [mut start, mut end, mid] = first;
        let (mut mid_start, mut mid_end) = (mid, mid);
        loop {
            let (mut best, mut next) = (f64::MAX, None);
            for (j, slot) in left.iter().enumerate() {
                let Some([s1, e1, m1]) = *slot else {
                    continue;
                };
                let (s1e, e1e) = (dist2(s1, end), dist2(e1, end));
                let (e1s, s1s) = (dist2(e1, start), dist2(s1, start));
                // the four links in upstream's order; a later, closer one wins
                if s1e < tol2 && s1e < best {
                    if e1s < tol2 && dist2(mid_end, m1) < tol2 {
                        return None;
                    }
                    (best, next) = (s1e, Some((j, 1)));
                }
                if e1e < tol2 && e1e < best {
                    if s1s < tol2 && dist2(mid_end, m1) < tol2 {
                        return None;
                    }
                    (best, next) = (e1e, Some((j, 2)));
                }
                if e1s < tol2 && e1s < best {
                    if s1e < tol2 && dist2(mid_start, m1) < tol2 {
                        return None;
                    }
                    (best, next) = (e1s, Some((j, 3)));
                }
                if s1s < tol2 && s1s < best {
                    if e1e < tol2 && dist2(mid_start, m1) < tol2 {
                        return None;
                    }
                    (best, next) = (s1s, Some((j, 4)));
                }
            }
            let Some((j, link)) = next else {
                break;
            };
            let Some([s1, e1, m1]) = left[j].take() else {
                break;
            };
            // 1 and 2 add at the end, 3 and 4 before the start; 2 and 4 reverse it
            let (s1, e1) = if link % 2 == 0 { (e1, s1) } else { (s1, e1) };
            if link <= 2 {
                chained.push([s1, e1, m1]);
                (end, mid_end) = (e1, m1);
            } else {
                chained.insert(0, [s1, e1, m1]);
                (start, mid_start) = (s1, m1);
            }
        }
        chains.push(chained);
    }
    Some(chains)
}

/// Whether a chain that closes (to the chaining tolerance) is narrower than that
/// tolerance: the vector area of the polygon through its curves' start and mid
/// points against its length. An open chain, such as a repeated side left over from
/// a square, is not tested.
fn encloses_nothing(chain: &[[V3; 3]]) -> bool {
    let (Some(&[start, ..]), Some(&[_, end, _])) = (chain.first(), chain.last()) else {
        return false;
    };
    if dist2(start, end) >= (2.0 * ZONE_TOL).powi(2) {
        return false;
    }
    let points: Vec<V3> = chain.iter().flat_map(|&[s, _, m]| [s, m]).collect();
    let (mut area, mut length) = ([0.0; 3], 0.0);
    for (i, &p) in points.iter().enumerate() {
        let q = points[(i + 1) % points.len()];
        let c = cross(p, q);
        area = [area[0] + c[0], area[1] + c[1], area[2] + c[2]];
        length += norm(sub(q, p));
    }
    0.5 * norm(area) < 2.0 * ZONE_TOL * length
}

/// UB-50 in `CreateTaperExtrude`: it first chains its curves. With a supplied normal
/// it may reverse them all and chain again, depending on their computed normal; that
/// second pass is checked whenever a normal is supplied, so a few inputs it would not
/// reverse are refused too. The end cap at the offset end (`caps` 2 or 3) chains the
/// trimmed offset curves through a copy of their list, so a curve deleted there is
/// freed again at cleanup; which offsets coincide cannot be told without building
/// them, so with that cap a closed chain enclosing nothing (every such probe that
/// crashed) is refused. The offsets' chaining before the caps (`SmPrimitiveCreation.cpp:3884`)
/// is not checked.
fn taper_is_unsafe(caps: f64, normal: &[f64], records: &[&[f64]]) -> bool {
    let Some(chains) = chain(records.iter().map(|r| points(r)).collect()) else {
        return true;
    };
    let supplied = normal.iter().map(|v| v * v).sum::<f64>() >= 0.00001;
    let reversed = || {
        chains
            .iter()
            .flatten()
            .map(|&[s, e, m]| [e, s, m])
            .collect()
    };
    (supplied && chain(reversed()).is_none())
        || (caps >= 2.0 && chains.iter().any(|c| encloses_nothing(c)))
}

/// Validate the layout, flags, counts and the inputs upstream cannot run safely
/// (also used by the Rust runner).
pub fn validate(op: i32, data: &[f64]) -> Result<(), FfiError> {
    if data.is_empty() || !data.iter().all(|v| v.is_finite() && v.abs() <= 1e6) {
        return Err(bad());
    }
    let ok = match op {
        1 => data.len() == 15 && flag(data[0]) && frame(&data[6..15]),
        2 => data.len() == 13 && flag(data[0]) && frame(&data[4..13]),
        3 => data.len() == 14 && flag(data[0]) && frame(&data[5..14]),
        4 => {
            data.len() >= 5
                && flag(data[0])
                && data[1] >= 0.0
                && integer(data[2], 3.0)
                && data[2] >= 1.0
                && flag(data[3])
                && integer(data[4], MAX_CURVES as f64)
                && {
                    let profiles = data[4] as usize;
                    data.len() >= 5 + profiles && {
                        let counts = &data[5..5 + profiles];
                        let total: f64 = counts.iter().sum();
                        counts.iter().all(|c| integer(*c, MAX_CURVES as f64))
                            && curves(&data[5 + profiles..])
                                .is_some_and(|c| total <= c.len() as f64)
                    }
                }
        }
        5 | 6 => {
            data.len() == 14
                && flag(data[0])
                && data[1] >= 0.0
                && data[2..7].iter().all(|v| flag(*v))
                && curve(&data[7..14])
        }
        7 => {
            data.len() == 12
                && flag(data[0])
                && data[2] >= 0.0
                && flag(data[3])
                && flag(data[4])
                && curve(&data[5..12])
        }
        8 => {
            data.len() >= 8
                && flag(data[0])
                && integer(data[3], 3.0)
                && integer(data[7], 2.0)
                && curves(&data[8..]).is_some_and(|c| {
                    c.len() >= 2
                        && c.iter().all(|r| extent(r) >= 2e-5)
                        && !taper_is_unsafe(data[3], &data[4..7], &c)
                })
        }
        9 => {
            data.len() >= 6
                && data[0] > 0.0
                && data[1..4].iter().all(|v| flag(*v))
                && integer(data[4], MAX_CURVES as f64)
                && integer(data[5], MAX_CURVES as f64)
                && data[5] >= 1.0
                && data.len() == 6 + 7 * (data[4] as usize + data[5] as usize)
                && curves(&data[6..]).is_some()
        }
        _ => false,
    };
    // the layout is checked first: the guards index it
    if ok && !sweeps_a_short_edge(op, data) && !builds_a_tiny_closed_shell(op, data) {
        Ok(())
    } else {
        Err(bad())
    }
}

/// Run the native entry point: `status edit doing_boolean` map counts, start, side and
/// end face counts, then the op's extra values.
pub fn eval(op: i32, data: &[f64]) -> Result<Vec<f64>, FfiError> {
    validate(op, data)?;
    call_with_buffer(64, |output, length| {
        // SAFETY: layout, counts, flags and curve records checked above; all
        // input/output buffers remain live for the complete synchronous call.
        unsafe {
            sys::sc_primitive_primitive_creation_p2_eval(
                op,
                data.as_ptr(),
                data.len() as i64,
                output.as_mut_ptr(),
                output.len() as i64,
                &raw mut *length,
            )
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn undefined_native_inputs_are_rejected() {
        // a taper of one curve, or of a curve shorter than 2e-5
        let line = [0., 0., 0., 0., 1., 0., 0.];
        let short = [0., 0., 0., 0., 1e-5, 0., 0.];
        let mut taper = vec![0., 1., 5., 3., 0., 0., 0., 0.];
        taper.extend(line);
        assert!(validate(8, &taper).is_err());
        taper.extend(short);
        assert!(validate(8, &taper).is_err());
        // a planar path without path curves
        let mut path = vec![1e-5, 1., 1., 1., 1., 0.];
        path.extend(line);
        assert!(validate(9, &path).is_err());
        // skin profile counts past the curves
        let mut skin = vec![0., 0., 1., 0., 2., 2., 1.];
        skin.extend(line);
        skin.extend(line);
        assert!(validate(4, &skin).is_err());
    }

    fn taper(caps: f64, normal_z: f64, curves: &[[f64; 7]]) -> Result<(), FfiError> {
        let mut d = vec![0., 1., 5., caps, 0., 0., normal_z, 0.];
        d.extend(curves.iter().flatten());
        validate(8, &d)
    }
    fn line(s: [f64; 3], e: [f64; 3]) -> [f64; 7] {
        [0., s[0], s[1], s[2], e[0], e[1], e[2]]
    }

    #[test]
    fn taper_curves_its_chaining_deletes_are_rejected() {
        // UB-50, each segfaulting natively: two identical lines, a line and its
        // reverse, two identical arcs
        let x = line([0., 0., 0.], [1., 0., 0.]);
        let back = line([1., 0., 0.], [0., 0., 0.]);
        let arc = [1., 0., 0., 0., 1., 0., 180.];
        for caps in [0., 3.] {
            assert!(taper(caps, 1., &[x, x]).is_err());
            assert!(taper(caps, 1., &[x, back]).is_err());
            assert!(taper(caps, 1., &[arc, arc]).is_err());
        }
        // and a line from a chain's end back to its start through its last link's
        // mid point, a repeat of no curve; 2.5e-5 off that mid point it runs
        let tiny = [
            line([0., 0., 0.], [3e-5, 0., 0.]),
            line([3e-5, 0., 0.], [1., 0., 0.]),
            back,
        ];
        for z in [1., -1., 0.] {
            assert!(taper(1., z, &tiny).is_err());
        }
        let short = [
            line([0., 0., 0.], [5e-5, 0., 0.]),
            line([5e-5, 0., 0.], [1., 0., 0.]),
            back,
        ];
        for z in [1., -1., 0.] {
            assert!(taper(1., z, &short).is_ok());
        }
        // a line and its reverse that the first pass chains apart: the reversed second
        // pass deletes one, and segfaults natively under the normal the curves wind
        // against (-Z). Under +Z upstream does not take that pass and runs, but the
        // pass is checked for any supplied normal; a computed normal runs and is kept
        let apart = [
            line([0., 1., 0.], [0., 0., 0.]),
            line([1., 1., 0.], [0., 1., 0.]),
            line([0., 1., 0.], [1., 1., 0.]),
        ];
        assert!(taper(1., -1., &apart).is_err());
        assert!(taper(1., 1., &apart).is_err());
        assert!(taper(1., 0., &apart).is_ok());
        // a link reversed as it joins decides the next duplicate test: these chain
        // apart and run natively (a replay that kept the link's direction deletes)
        let turned = [
            line([0.5, 0., 0.], [0., 1., 0.]),
            line([0.5, 0., 0.], [1., 1., 0.]),
            line([0., 1., 0.], [0.5, 0., 0.]),
        ];
        assert!(taper(1., 0., &turned).is_ok());
        // the end cap at the offset end chains the offsets of a loop that encloses
        // nothing: the short link's loop at side 100 under -Z frees a curve twice
        // natively (caps 2 and 3, not 0 and 1); every closed chain that narrow is
        // refused with that cap, the short link at side 1 under +Z too, which runs.
        // The open chain of the line and its reverse is not tested, and runs
        let far = [
            line([0., 0., 0.], [5e-5, 0., 0.]),
            line([5e-5, 0., 0.], [100., 0., 0.]),
            line([100., 0., 0.], [0., 0., 0.]),
        ];
        for caps in [2., 3.] {
            assert!(taper(caps, -1., &far).is_err());
            assert!(taper(caps, 1., &short).is_err());
        }
        assert!(taper(1., -1., &far).is_ok());
        assert!(taper(3., 0., &apart).is_ok());
        // a square with a repeated or reversed side deletes nothing in either pass,
        // keeps its area and runs natively with both caps, as do two halves of a
        // circle and a half circle's chord
        let [p0, p1, p2, p3] = [[0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.]];
        let square = [line(p0, p1), line(p1, p2), line(p2, p3), line(p3, p0)];
        for side in [
            line(p1, p2),
            line(p2, p1),
            line(p2, p3),
            line(p0, p1),
            line(p3, p0),
        ] {
            for z in [1., -1., 0.] {
                let mut curves = square.to_vec();
                curves.push(side);
                assert!(taper(3., z, &curves).is_ok());
            }
        }
        assert!(taper(3., 1., &[arc, [1., 0., 0., 0., 1., 180., 360.]]).is_ok());
        assert!(taper(3., 1., &[arc, line([-1., 0., 0.], [1., 0., 0.])]).is_ok());
    }

    fn primitive(op: i32, sizes: &[f64], start: f64, end: f64) -> Result<(), FfiError> {
        let mut d = vec![0.];
        d.extend(sizes);
        d.extend([start, end, 0., 0., 0., 1., 0., 0., 0., 1., 0.]);
        validate(op, &d)
    }

    #[test]
    fn short_swept_edges_are_rejected_and_early_exits_kept() {
        // UB-215 at its native edges: a 9e-6 cone line segfaults, 1e-5 does not
        assert!(primitive(1, &[1., 1., 9e-6], 0., 90.).is_err());
        assert!(primitive(1, &[2e-6, 1., 1.], 0., 90.).is_err());
        assert!(primitive(1, &[1., 1., 1e-5], 0., 90.).is_ok());
        // a minor circle of 2 pi 1e-6 segfaults, 2 pi 2e-6 does not
        assert!(primitive(3, &[0., 1e-6], 0., 90.).is_err());
        assert!(primitive(3, &[1., 2e-6], 0., 90.).is_ok());
        // failures before the sweep stay compared: no height, no turn, a zero side,
        // a radius within the scaled zero (not built), a zero minor radius
        assert!(primitive(1, &[0., 1., 1.], 0., 90.).is_ok());
        assert!(primitive(1, &[1., 1., 1e-6], 30., 30.).is_ok());
        assert!(primitive(1, &[1e-12, 1., 1.], 0., 90.).is_ok());
        assert!(primitive(1, &[1., 1., 1e-13], 0., 90.).is_ok());
        assert!(primitive(3, &[1., 0.], 0., 90.).is_ok());
        // analytic faces sweep no edges; spheres never reach the freed curve
        assert!(primitive(1, &[1., 1., 1e-6], 0., 360.).is_ok());
        assert!(primitive(2, &[1e-6], 0., 90.).is_ok());
    }

    #[test]
    fn tiny_closed_shells_are_rejected() {
        // HANG-35: a sphere of 1e-5 takes 87 s, 1e-4 over 359 degrees 20 s; 4e-6 and
        // 1.5e-4 take milliseconds, at any turn
        assert!(primitive(2, &[1e-5], 0., 360.).is_err());
        assert!(primitive(2, &[5e-5], 30., 390.000_000_003_6).is_err());
        assert!(primitive(2, &[1e-4], 0., 359.).is_err());
        assert!(primitive(2, &[1.5e-4], 0., 359.).is_ok());
        assert!(primitive(2, &[4e-6], 0., 360.).is_ok());
        assert!(primitive(1, &[2e-5, 2e-5, 2e-5], 0., 359.).is_err());
        assert!(primitive(1, &[1.5e-4, 1.5e-4, 1.5e-4], 0., 359.).is_ok());
        assert!(primitive(3, &[1e-5, 1e-5], 0., 360.).is_err());
        assert!(primitive(3, &[1., 1e-5], 0., 360.).is_ok());
        // a swept torus is stitched and placed the same way: 3e-6 by 3e-6 over 359
        // degrees takes minutes, 2e-6 by 2e-6 does not
        assert!(primitive(3, &[3e-6, 3e-6], 0., 359.).is_err());
        assert!(primitive(3, &[2e-6, 2e-6], 0., 359.).is_ok());
        // failures before any face is built stay compared
        assert!(primitive(2, &[0.], 0., 360.).is_ok());
        assert!(primitive(3, &[0., 1e-5], 0., 360.).is_ok());
        assert!(primitive(1, &[2e-5, 2e-5, 2e-5], 30., 30.).is_ok());
        assert!(primitive(1, &[1e-13, 2e-5, 2e-5], 0., 360.).is_ok());
    }
}
