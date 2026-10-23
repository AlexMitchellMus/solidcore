//! Where two bodies' surfaces cross.
//!
//! The method is subdivision for completeness and marching for accuracy:
//!
//!   1. Every B-spline is cut into rational Bézier pieces by knot insertion
//!      (The NURBS Book, algorithm A5.1 repeated, as A5.6 does for a curve).
//!   2. Pairs of pieces, one from each body, are halved until they are apart
//!      (bounding boxes, or slabs across a normal cone) or loop-free: their
//!      normal cones hold no parallel pair of normals, so no intersection curve
//!      closes inside the pair, every one crosses transversally, and each has
//!      both ends on the pair's boundary (Sederberg and Meyers's loop test).
//!   3. In a loop-free pair, each boundary curve of each piece is intersected
//!      with the other piece: those points are the seeds. Each arc is marched
//!      from a seed to its exit, a predictor along nA × nB and a Newton
//!      corrector held to a plane across the curve.
//!   4. Arcs are oriented along nA × nB and joined, end to start, into curves.
//!
//! Orientation: with both bodies' normals pointing out, a curve runs along
//! nA × nB. Seen from outside body A, the part of A inside B lies to the left
//! of the curve, and the part outside B to the right (T × nA points out of B).
//!
//! A pair that is still neither apart nor loop-free when both pieces are
//! smaller than `resolution` is a place where the surfaces touch or coincide,
//! or where a piece is degenerate (a collapsed edge at a pole) right at the
//! crossing. That is refused by name, [`Error::Tangent`].
//!
//! Where a patch has an analytic surface (`Patch::analytic`), its implicit
//! equation composed with the other body's pieces (`field`) adds two proofs: a
//! piece apart from the surface leaves the search, and a patch lying on it
//! shares a region with it, whose boundary is the intersection (`coincide`).
//!
//! Ported from Form Studio's `src/kernel/intersect.js` (same author) as of its
//! 2026-10-22 fixes, and grown past it since.

use std::collections::HashMap;

use crate::surface3::Surf;

use crate::algebraic::{Algebraic, Fitter, MAX_DEGREE};
use crate::coincide::{self, Spacing};
use crate::field::Field;
use crate::field::{relate, Relation};
use crate::patch::Patch;
use crate::piece::{
    boxes_apart, cross, curve_at, curve_hull, dot, loop_free, norm, pieces_apart, scale,
    slab_apart, split_curve, sub, unit, Eval, Piece, H, V3,
};
use crate::small::{newton_step, solve_damped, solve_regularized, System};
use crate::trace2d;
use crate::Error;

/// Every option [`intersect_bodies`] takes. Lengths are fractions of the bodies' size.
#[derive(Clone, Copy, Debug)]
pub struct Options {
    /// Newton converges when the two surfaces' points are this close (scaled by the sine of their angle).
    pub tolerance: f64,
    /// Two arc ends this close are one point of a curve.
    pub link: f64,
    /// Pieces smaller than this that are still unresolved are a tangency.
    pub resolution: f64,
    /// The longest marching step.
    pub max_step: f64,
    /// The most a step may turn (radians).
    pub max_turn: f64,
    /// Pairs of pieces visited before the search is refused.
    pub budget: usize,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            tolerance: 1e-11,
            link: 1e-7,
            resolution: 1e-5,
            max_step: 0.02,
            max_turn: 0.06,
            budget: 2_000_000,
        }
    }
}

/// A surface's parameters at a curve point, in its patch's own knot domain.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Param {
    pub patch: usize,
    pub u: f64,
    pub v: f64,
}

/// A point of an intersection curve: the midpoint of the two surfaces' points, and each one's parameters.
#[derive(Clone, Copy, Debug)]
pub struct CurvePoint {
    pub p: V3,
    pub a: Param,
    pub b: Param,
}

/// An intersection curve. A closed curve does not repeat its first point.
#[derive(Clone, Debug)]
pub struct Curve {
    pub closed: bool,
    pub points: Vec<CurvePoint>,
}

/// Counts from one run: `candidates` are the pairs of pieces whose boxes clash,
/// `pairs` every pair (and curve node) visited while refining them.
#[derive(Clone, Copy, Debug, Default)]
pub struct Stats {
    pub pieces_a: usize,
    pub pieces_b: usize,
    pub candidates: usize,
    pub pairs: usize,
    pub leaves: usize,
    pub arcs: usize,
}

#[derive(Clone, Debug)]
pub struct Intersection {
    pub curves: Vec<Curve>,
    /// Isolated points where the surfaces touch.
    pub points: Vec<CurvePoint>,
    /// Points on the curves where branches cross or meet: curves end at a
    /// crossing; two branches tangent to each other pass through their meeting.
    pub vertices: Vec<V3>,
    pub stats: Stats,
}

struct Ctx {
    tolerance: f64,
    link: f64,
    resolution: f64,
    pad: f64,
    max_step: f64,
    max_turn: f64,
    budget: usize,
    visits: usize,
}

impl Ctx {
    fn visit(&mut self) -> Result<(), Error> {
        self.visits += 1;
        if self.visits > self.budget {
            return Err(Error::Tangent {
                at: None,
                detail: format!("the surfaces {EXHAUSTED} {} pairs of pieces; they touch or coincide over a region", self.budget),
            });
        }
        Ok(())
    }
}

/// The curves along which body A's surface crosses body B's. Curves run along
/// nA × nB (see the module head). Two closed bodies give only closed curves;
/// an open curve ends on an open body's boundary.
///
/// A pair refused as tangent is tried again with the implicit forms of the
/// B-spline patches that have one (`algebraic`: a NURBS cylinder, sphere,
/// torus), which put the pair on the field path; fitted only then, so a pair
/// that needs none pays nothing for them.
///
/// Where a patch has no analytic surface, that first try runs on 1/64 of the
/// pair budget: a pair that needs implicit forms shows it by being refused,
/// and need not spend the whole budget first (2 s of a 2.6 s case,
/// ssi-corpus c2n_quartic_contact_parabola). A pair refused only for the
/// short budget, with no forms to try, is tried again on the whole one.
///
/// One implicit form a pair is enough for the field path, so a patch whose
/// every pair already has a field on the other side is not fitted further:
/// fitting is most of a hard pair's time (2.7 s on a patch with no form,
/// ssi-corpus y23_fig1b_two_arcs_and_point). The patches are fitted in step,
/// a degree at a time, so the form of lowest degree is the one found (one
/// patch fitted first to its degree 8 left the other, planar, unfitted,
/// ssi-corpus y23_t6_ex4_tangential_point); within a degree, in the order the
/// field path takes its implicit side (`trace2d::roles`): a patch with a
/// collapsed edge first, then B's (fitted first, a NURBS cylinder was made
/// the implicit side against a NURBS sphere, whose poles were then traced,
/// ssi-corpus c4n_viviani_nurbs).
pub fn intersect_bodies(
    body_a: &[Patch],
    body_b: &[Patch],
    options: &Options,
) -> Result<Intersection, Error> {
    let none = |body: &[Patch]| vec![None; body.len()];
    let bspline = body_a
        .iter()
        .chain(body_b)
        .any(|patch| patch.analytic.is_none());
    let short = Options {
        budget: if bspline {
            (options.budget / 64).max(1)
        } else {
            options.budget
        },
        ..*options
    };
    // Quadric forms first, where neither patch of a pair has an analytic
    // surface: two small fits, and a NURBS paraboloid or cylinder then takes
    // the field path, not subdivision and marching (10 times as fast on two
    // transversal paraboloids, ssi-corpus c1n_parab_parab_circle). Where the
    // other side has one the field path applies already, and the fit is only
    // cost (3.9 -> 10.7 ms, c9n_deg6_circle).
    let quadric = |body: &[Patch], others: &[Patch]| -> Vec<Option<Algebraic>> {
        body.iter()
            .enumerate()
            .map(|(i, p)| {
                if p.analytic.is_some() || others.iter().all(|o| o.analytic.is_some()) {
                    return None;
                }
                let mut fitter = Fitter::new(p, i)?;
                (0..2).find_map(|_| if fitter.open() { fitter.next() } else { None })
            })
            .collect()
    };
    let (pre_a, pre_b) = (quadric(body_a, body_b), quadric(body_b, body_a));
    match intersect_with(body_a, body_b, &short, &pre_a, &pre_b) {
        Err(refused @ Error::Tangent { .. }) => {
            let (mut forms_a, mut forms_b) = (none(body_a), none(body_b));
            let mut fitters: Vec<(bool, usize, Fitter)> = body_a
                .iter()
                .enumerate()
                .map(|(i, p)| (true, i, p))
                .chain(body_b.iter().enumerate().map(|(i, p)| (false, i, p)))
                .filter(|(_, _, p)| p.analytic.is_none())
                .filter_map(|(in_a, i, p)| Fitter::new(p, i).map(|f| (in_a, i, f)))
                .collect();
            fitters.sort_by_key(|(in_a, i, _)| {
                (
                    !(if *in_a { &body_a[*i] } else { &body_b[*i] }).has_collapsed_edge(),
                    *in_a,
                )
            });
            let covered = |body: &[Patch], forms: &[Option<Algebraic>]| {
                body.iter()
                    .zip(forms)
                    .all(|(p, f)| p.analytic.is_some() || f.is_some())
            };
            for _ in 0..MAX_DEGREE {
                for (in_a, i, fitter) in &mut fitters {
                    let (done, pairs_covered) = if *in_a {
                        (forms_a[*i].is_some(), covered(body_b, &forms_b))
                    } else {
                        (forms_b[*i].is_some(), covered(body_a, &forms_a))
                    };
                    if done || pairs_covered || !fitter.open() {
                        continue;
                    }
                    if let Some(form) = fitter.next() {
                        if *in_a {
                            forms_a[*i] = Some(form);
                        } else {
                            forms_b[*i] = Some(form);
                        }
                    }
                }
            }
            if forms_a.iter().chain(&forms_b).all(Option::is_none) {
                if short.budget < options.budget && exhausted(&refused) {
                    return intersect_with(body_a, body_b, options, &pre_a, &pre_b);
                }
                return Err(refused);
            }
            intersect_with(body_a, body_b, options, &forms_a, &forms_b)
        }
        other => other,
    }
}

/// Whether a refusal is the pair budget's (`Ctx::visit`), not a tangency found.
fn exhausted(error: &Error) -> bool {
    matches!(error, Error::Tangent { detail, .. } if detail.contains(EXHAUSTED))
}

/// The words of the pair budget's refusal.
const EXHAUSTED: &str = "could not be separated within";

/// `intersect_bodies` with the patches' proved algebraic forms, where fitted.
fn intersect_with(
    body_a: &[Patch],
    body_b: &[Patch],
    options: &Options,
    forms_a: &[Option<Algebraic>],
    forms_b: &[Option<Algebraic>],
) -> Result<Intersection, Error> {
    let fields_a: Vec<Option<Field>> = body_a
        .iter()
        .zip(forms_a)
        .map(|(p, f)| {
            p.analytic
                .as_ref()
                .map(Field::Analytic)
                .or(f.as_ref().map(Field::Algebraic))
        })
        .collect();
    let fields_b: Vec<Option<Field>> = body_b
        .iter()
        .zip(forms_b)
        .map(|(p, f)| {
            p.analytic
                .as_ref()
                .map(Field::Analytic)
                .or(f.as_ref().map(Field::Algebraic))
        })
        .collect();
    // Each patch's Bézier pieces, made once for every use below (five times,
    // a 31 x 31 span plane's took most of its 39 ms, ssi-corpus
    // c9n_many_spans_plane_sphere).
    let sets_a: Vec<Vec<Piece>> = body_a
        .iter()
        .enumerate()
        .map(|(i, p)| p.bezier_pieces(i))
        .collect::<Result<_, _>>()?;
    let sets_b: Vec<Vec<Piece>> = body_b
        .iter()
        .enumerate()
        .map(|(i, p)| p.bezier_pieces(i))
        .collect::<Result<_, _>>()?;
    let pieces_a: Vec<Piece> = sets_a.iter().flatten().cloned().collect();
    let pieces_b: Vec<Piece> = sets_b.iter().flatten().cloned().collect();
    let mut bounds = [
        f64::INFINITY,
        f64::INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
    ];
    for piece in pieces_a.iter().chain(&pieces_b) {
        let b = piece.bounds();
        for k in 0..3 {
            bounds[k] = bounds[k].min(b[k]);
            bounds[k + 3] = bounds[k + 3].max(b[k + 3]);
        }
    }
    let size = match norm([
        bounds[3] - bounds[0],
        bounds[4] - bounds[1],
        bounds[5] - bounds[2],
    ]) {
        s if s > 0.0 => s,
        _ => 1.0,
    };
    // The pad keeps rounding in the halved control points from parting pieces
    // that touch, so it is a rounding bound (ulps of the largest coordinate),
    // not a multiple of the Newton tolerance: at 10 x tolerance (3.5e-10 on a
    // 3.5 mm body) a loop shallower than that, a paraboloid cut 1e-10 above
    // its apex (ssi-corpus c7_tiny_loop_parab_1e-5), never parted from the
    // plane and was refused as a tangency.
    let magnitude = bounds.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
    // Where a curve grazes a piece's edge the seed there is a double root,
    // found only to the square root of the rounding in the gap (1.4e-8 apart
    // on the two sides of a radius-1 circle inscribed in a patch's square,
    // ssi-corpus c1n_parab_plane_circle), so a link finer than that leaves the
    // arcs unjoined and copies of them undeduplicated. The defaults are above
    // this floor.
    let grazing = 8.0 * (f64::EPSILON * magnitude).sqrt();
    let mut ctx = Ctx {
        tolerance: options.tolerance * size,
        link: (options.link * size).max(grazing),
        resolution: options.resolution * size,
        pad: 1024.0 * f64::EPSILON * magnitude,
        max_step: options.max_step * size,
        max_turn: options.max_turn,
        budget: options.budget,
        visits: 0,
    };

    // Patches lying on one surface share a region, and its boundary is their
    // intersection; the pieces of such a pair are never subdivided, as every
    // pair of them would touch everywhere.
    // Sampled to a 0.01 rad turn at most whatever the marching's limit: a shared
    // boundary's edges are cheap to sample, and `densify` cannot fill them in,
    // their surfaces' normals being parallel.
    let spacing = Spacing {
        max_step: ctx.max_step,
        max_sag: ctx.max_step * ctx.max_turn.min(0.01) / 8.0,
    };
    // Where one patch has an analytic surface, the field path traces the pair
    // on the other's parameters (`trace2d`); a pair it does not resolve is
    // subdivided as before.
    let field = trace2d::Settings {
        resolution: ctx.resolution,
        contact: ctx.tolerance,
        spacing: Spacing {
            max_step: ctx.max_step,
            max_sag: ctx.max_step * ctx.max_turn / 8.0,
        },
    };
    let mut coincident: Vec<(usize, usize)> = Vec::new();
    let mut shared = Vec::new();
    let mut points: Vec<CurvePoint> = Vec::new();
    let mut vertices: Vec<V3> = Vec::new();
    for (ia, pa) in body_a.iter().enumerate() {
        for (ib, pb) in body_b.iter().enumerate() {
            if let Some(arcs) = coincide::boundary(
                pa,
                ia,
                fields_a[ia],
                &sets_a[ia],
                pb,
                ib,
                fields_b[ib],
                &sets_b[ib],
                &spacing,
            )? {
                coincident.push((ia, ib));
                shared.extend(arcs);
            } else if let Some((arcs, touches, crossings)) = trace2d::trace_pair(
                pa,
                ia,
                fields_a[ia],
                &sets_a[ia],
                pb,
                ib,
                fields_b[ib],
                &sets_b[ib],
                &field,
            )? {
                coincident.push((ia, ib));
                shared.extend(arcs);
                vertices.extend(crossings);
                // A touch point on a piece's edge is found in each piece beside it.
                for q in touches {
                    if !points
                        .iter()
                        .any(|other| norm(sub(other.p, q.p)) < ctx.link)
                    {
                        points.push(q);
                    }
                }
            }
        }
    }
    // A piece that the other patch's analytic surface proves apart from it is
    // dropped from the search, whatever its box and normal cone say.
    let analytic_apart = |a: &Piece, b: &Piece| {
        body_b[b.patch]
            .analytic
            .as_ref()
            .is_some_and(|s| relate(a, s) == Relation::Apart)
            || body_a[a.patch]
                .analytic
                .as_ref()
                .is_some_and(|s| relate(b, s) == Relation::Apart)
    };

    // Pairs of pieces, halved until apart or loop-free; depth first, so a
    // tangency is found, and refused, before the search spreads.
    let mut leaves: Vec<(&Piece, &Piece)> = Vec::new();
    let mut candidates = 0;
    for a0 in &pieces_a {
        for b0 in &pieces_b {
            if coincident.contains(&(a0.patch, b0.patch))
                || boxes_apart(a0.bounds(), b0.bounds(), ctx.pad)
            {
                continue;
            }
            candidates += 1;
            let mut stack = vec![(a0, b0)];
            while let Some((a, b)) = stack.pop() {
                if let Err(Error::Tangent { detail, .. }) = ctx.visit() {
                    let at = a.evaluate(0.5, 0.5).point;
                    return Err(Error::Tangent {
                        at: Some(at),
                        detail: format!(
                            "{detail}; the last pair was {} and {}",
                            describe(a),
                            describe(b)
                        ),
                    });
                }
                if pieces_apart(a, b, ctx.pad) {
                    continue;
                }
                if loop_free(a, b) {
                    leaves.push((a, b));
                    continue;
                }
                if analytic_apart(a, b) {
                    continue;
                }
                let split_a = if a.size() >= b.size() {
                    a.size() >= ctx.resolution
                } else {
                    b.size() < ctx.resolution
                };
                let piece = if split_a { a } else { b };
                if piece.size() < ctx.resolution {
                    let at = a.evaluate(0.5, 0.5).point;
                    return Err(Error::Tangent {
                        at: Some(at),
                        detail: format!(
                            "the surfaces touch or coincide near ({:.5}, {:.5}, {:.5})",
                            at[0], at[1], at[2]
                        ),
                    });
                }
                let dir = if piece.extent(0) >= piece.extent(1) {
                    0
                } else {
                    1
                };
                for half in piece.halves(dir) {
                    stack.push(if split_a { (half, b) } else { (a, half) });
                }
            }
        }
    }

    let mut arcs = shared;
    for &(a, b) in &leaves {
        arcs.extend(trace_leaf(
            a,
            b,
            body_a[a.patch].analytic.as_ref(),
            body_b[b.patch].analytic.as_ref(),
            &mut ctx,
        )?);
    }
    let arc_count = arcs.len();
    let curves = join_arcs(arcs, &vertices, &ctx);
    let mut unique_vertices: Vec<V3> = Vec::new();
    for v in vertices {
        if !unique_vertices.iter().any(|&u| norm(sub(u, v)) < ctx.link) {
            unique_vertices.push(v);
        }
    }
    Ok(Intersection {
        curves,
        points,
        vertices: unique_vertices,
        stats: Stats {
            pieces_a: pieces_a.len(),
            pieces_b: pieces_b.len(),
            candidates,
            pairs: ctx.visits,
            leaves: leaves.len(),
            arcs: arc_count,
        },
    })
}

/// A piece for a refusal's detail: its patch, knot rectangle, size and normal cone.
fn describe(piece: &Piece) -> String {
    let cone = match piece.cone() {
        Some(c) => format!(
            "cone {:.3e} rad about ({:.3}, {:.3}, {:.3})",
            c.angle, c.axis[0], c.axis[1], c.axis[2]
        ),
        None => "no normal cone".into(),
    };
    let p = piece.evaluate(0.5, 0.5).point;
    format!(
        "patch {} [{:.6}, {:.6}] x [{:.6}, {:.6}] (size {:.3e}, {cone}, centre ({:.6}, {:.6}, {:.6}))",
        piece.patch, piece.u0, piece.u1, piece.v0, piece.v1, piece.size(), p[0], p[1], p[2]
    )
}

/// The parameter velocity (ds, dt) whose image under [Ss St] is closest to `v`.
/// Damped by 1e-12 of the mean squared derivative, so it stays defined where a
/// derivative vanishes (a pole, where Ss = 0): the velocity there is all along St.
fn velocity(ds: V3, dt: V3, v: V3) -> [f64; 2] {
    let damping = 0.5e-12 * (dot(ds, ds) + dot(dt, dt));
    solve_regularized(&System::columns2(ds, dt, v), damping).map_or([0.0, 0.0], |x| [x[0], x[1]])
}

/// The surface normal Ss × St at local (s, t) of `piece`, evaluated as `e`.
/// Where it vanishes (a collapsed edge: a pole, a cone's apex) it is taken a
/// hair inside the piece, whose direction tends to the surface's own normal.
fn normal(piece: &Piece, s: f64, t: f64, e: &Eval) -> V3 {
    let n = cross(e.ds, e.dt);
    let scale2 = dot(e.ds, e.ds).max(dot(e.dt, e.dt));
    if dot(n, n) > 1e-20 * scale2 * scale2 {
        return n;
    }
    let nudge = |x: f64| x + (0.5 - x) * 2e-6;
    let inside = piece.evaluate(nudge(s), nudge(t));
    cross(inside.ds, inside.dt)
}

// ---------------------------------------------------------------------------
// Curve against surface: the seeds

struct Node<'p> {
    c: Vec<H>,
    t0: f64,
    t1: f64,
    piece: &'p Piece,
}

/// Every isolated point [t, s, u] where the homogeneous curve `curve` (t in
/// [0, 1]) meets the piece `surface`.
///
/// More than 64 distinct roots means the curve lies along the surface: a
/// loop-free pair's boundary curves are low degree and cross the other piece
/// at a handful of isolated points at most. Such a curve is an intersection arc
/// lying on the pair's boundary (a seam, a knot line or a patch edge on the
/// other surface), and its roots are no seeds, so none are returned; the arc
/// is marched from its two ends, which the neighbouring boundary curves meet
/// transversally.
fn curve_surface(curve: &[H], surface: &Piece, ctx: &mut Ctx) -> Result<Vec<[f64; 3]>, Error> {
    let mut roots: Vec<([f64; 3], V3)> = Vec::new();
    let mut stack = vec![Node {
        c: curve.to_vec(),
        t0: 0.0,
        t1: 1.0,
        piece: surface,
    }];
    let residual_at = |x: [f64; 3]| {
        norm(sub(
            curve_at(curve, x[0]).0,
            surface.evaluate(x[1], x[2]).point,
        ))
    };
    // Within tolerance is not converged: where an edge grazes the other surface
    // (a double root, Newton halving its distance each step) the first point
    // within tolerance lies sqrt(tolerance) along the edge from the touch, and a
    // seed there starts its arc that far short of it, so the arcs on the two
    // sides never join (a circle inscribed in a patch's square, ssi-corpus
    // c1n_parab_plane_circle: gaps of 1.5e-5 against a link of 8.5e-7). So
    // Newton runs on until its step stops shrinking, keeping the best point.
    //
    // The best point is a root when the curve's distance off the surface is
    // under the tolerance: the gap divided by the sine of the curve's angle to
    // the surface, so the gap is held to tolerance × sin (above rounding). A
    // near miss, an edge passing 3e-12 above a paraboloid where the tolerance
    // is 3.5e-11 (ssi-corpus c7_tiny_loop_parab_1e-6), is tangent at its closest
    // approach and so never passes; a true grazing root or an edge lying on the
    // surface converges to rounding and does. The sine is floored at 1e-2, so
    // the gap is never held below 1% of the tolerance: exact input is exact only
    // to its floating-point control points, and a line lying on a NURBS
    // cylinder's seam (ssi-corpus c8n_cylN_seam_line) stays 2.4e-14 off it.
    let newton = |t: f64, s: f64, u: f64, tolerance: f64| -> Option<[f64; 3]> {
        let mut x = [t, s, u];
        let mut best = x;
        let mut best_residual = f64::INFINITY;
        for _ in 0..MAX_NEWTON {
            let (c, cd) = curve_at(curve, x[0]);
            let e = surface.evaluate(x[1], x[2]);
            let residual = norm(sub(c, e.point));
            if residual < best_residual {
                best = x;
                best_residual = residual;
            }
            if residual == 0.0 {
                break;
            }
            let mut system = System::new(3, 3);
            let r = sub(e.point, c);
            for k in 0..3 {
                system.j[k] = [cd[k], -e.ds[k], -e.dt[k], 0.0];
                system.r[k] = r[k];
            }
            let Some(step) = solve_damped(&system) else {
                break;
            };
            x = [x[0] + step[0], x[1] + step[1], x[2] + step[2]];
            if x.iter().any(|v| !v.is_finite() || *v < -0.25 || *v > 1.25) {
                break;
            }
            if step[..3].iter().fold(0.0_f64, |m, v| m.max(v.abs())) < 1e-15 {
                let last = residual_at(x);
                if last < best_residual {
                    best = x;
                    best_residual = last;
                }
                break;
            }
        }
        let (c, cd) = curve_at(curve, best[0]);
        let e = surface.evaluate(best[1], best[2]);
        let n = normal(surface, best[1], best[2], &e);
        let sine = dot(cd, n).abs() / (norm(cd) * norm(n));
        let rounding = 64.0 * f64::EPSILON * c.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
        let gap = (tolerance
            * if (0.0..1.0).contains(&sine) {
                sine.max(1e-2)
            } else {
                1.0
            })
        .max(rounding);
        (best_residual < gap).then_some(best)
    };
    while let Some(node) = stack.pop() {
        ctx.visit()?;
        let hull = curve_hull(&node.c);
        let piece = node.piece;
        if boxes_apart(&hull.bounds, piece.bounds(), ctx.pad) {
            continue;
        }
        if let Some(cone) = piece.cone() {
            if slab_apart(&hull.points, piece.points(), cone.axis, ctx.pad) {
                continue;
            }
        }
        let small = hull.size < ctx.resolution && piece.size() < ctx.resolution;
        let mut transversal = false;
        if let (Some(hc), Some(pc)) = (hull.cone, piece.cone()) {
            let theta = dot(hc.axis, pc.axis).clamp(-1.0, 1.0).acos();
            let spread = hc.angle + pc.angle;
            transversal = theta + spread < std::f64::consts::FRAC_PI_2 - 1e-9
                || theta - spread > std::f64::consts::FRAC_PI_2 + 1e-9;
        }
        if transversal || small {
            let (s, u) = surface.local_of(piece, 0.5, 0.5);
            let margin = 1e-9;
            if let Some(root) = newton((node.t0 + node.t1) / 2.0, s, u, ctx.tolerance) {
                if root.iter().all(|&v| v >= -margin && v <= 1.0 + margin) {
                    let clamped = root.map(|v| v.clamp(0.0, 1.0));
                    let point = curve_at(curve, clamped[0]).0;
                    if !roots
                        .iter()
                        .any(|(_, other)| norm(sub(*other, point)) < ctx.link)
                    {
                        roots.push((clamped, point));
                    }
                    if roots.len() > 64 {
                        return Ok(Vec::new());
                    }
                    continue;
                }
            }
            if small {
                continue;
            }
        }
        if hull.size >= piece.size() {
            let (l, r) = split_curve(&node.c);
            let mid = (node.t0 + node.t1) / 2.0;
            stack.push(Node {
                c: l,
                t0: node.t0,
                t1: mid,
                piece,
            });
            stack.push(Node {
                c: r,
                t0: mid,
                t1: node.t1,
                piece,
            });
        } else {
            let dir = if piece.extent(0) >= piece.extent(1) {
                0
            } else {
                1
            };
            for half in piece.halves(dir) {
                stack.push(Node {
                    c: node.c.clone(),
                    t0: node.t0,
                    t1: node.t1,
                    piece: half,
                });
            }
        }
    }
    Ok(roots.into_iter().map(|(x, _)| x).collect())
}

// ---------------------------------------------------------------------------
// One loop-free pair: seeds and arcs

const MAX_NEWTON: usize = 40;
const MARGIN: f64 = 1e-10;

fn outside(q: &[f64; 4]) -> bool {
    q.iter().any(|&x| x < -MARGIN || x > 1.0 + MARGIN)
}

/// Both pieces at q = [sA, tA, sB, tB].
#[derive(Clone, Copy)]
struct State {
    q: [f64; 4],
    a: Eval,
    b: Eval,
    tangent: V3,
}

fn both(a: &Piece, b: &Piece, q: [f64; 4]) -> State {
    let ea = a.evaluate(q[0], q[1]);
    let eb = b.evaluate(q[2], q[3]);
    State {
        q,
        a: ea,
        b: eb,
        tangent: cross(normal(a, q[0], q[1], &ea), normal(b, q[2], q[3], &eb)),
    }
}

struct Solved {
    q: [f64; 4],
    ab: Option<(Eval, Eval)>,
    ok: bool,
}

/// Newton on the pair at q: the two points equal, and either the point on the
/// plane through `target` across `normal` (a marching corrector), or the
/// parameters in `fixed` held where they are (an exit onto a boundary).
fn newton_pair(
    a: &Piece,
    b: &Piece,
    q0: [f64; 4],
    target: Option<(V3, V3)>,
    fixed: &[usize],
    tolerance: f64,
) -> Solved {
    let mut q = q0;
    let mut free = [0usize; 4];
    let mut n_free = 0;
    for k in (0..4).filter(|k| !fixed.contains(k)) {
        free[n_free] = k;
        n_free += 1;
    }
    let free = &free[..n_free];
    for _ in 0..MAX_NEWTON {
        let ea = a.evaluate(q[0], q[1]);
        let eb = b.evaluate(q[2], q[3]);
        let residual = sub(ea.point, eb.point);
        let columns = [ea.ds, ea.dt, scale(eb.ds, -1.0), scale(eb.dt, -1.0)];
        let mut system = System::new(if target.is_some() { 4 } else { 3 }, n_free);
        for row in 0..3 {
            for (i, &k) in free.iter().enumerate() {
                system.j[row][i] = columns[k][row];
            }
            system.r[row] = -residual[row];
        }
        let mut plane_error = 0.0;
        if let Some((point, normal)) = target {
            plane_error = dot(sub(ea.point, point), normal);
            for (i, &k) in free.iter().enumerate() {
                system.j[3][i] = if k < 2 { dot(columns[k], normal) } else { 0.0 };
            }
            system.r[3] = -plane_error;
        }
        // A gap g between the two points leaves the curve point up to g / sin of
        // the surfaces' angle off the curve, so the gap is held to tolerance × sin
        // (above rounding): at the 2e-4 rad crossing of a paraboloid cut 1e-8
        // above its apex (ssi-corpus c7_tiny_loop_parab_1e-4) a gap of 3.5e-11
        // left points 2.4e-7 off the radius-1e-4 loop, and two traces of it too
        // far apart to be taken for one. Where a normal vanishes (a pole) the
        // sine is unknown: 1.
        let na = normal(a, q[0], q[1], &ea);
        let nb = normal(b, q[2], q[3], &eb);
        let sine = norm(cross(na, nb)) / (norm(na) * norm(nb));
        let rounding = 16.0 * f64::EPSILON * ea.point.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
        let gap = (tolerance
            * if (0.0..1.0).contains(&sine) {
                sine
            } else {
                1.0
            })
        .max(rounding);
        if norm(residual) < gap && plane_error.abs() < tolerance {
            return Solved {
                q,
                ab: Some((ea, eb)),
                ok: true,
            };
        }
        let fits = |step: &[f64]| {
            free.iter().enumerate().all(|(i, &k)| {
                (q[k] + step[i]).is_finite() && (-0.5..=1.5).contains(&(q[k] + step[i]))
            })
        };
        let Some(step) = newton_step(&system, fits) else {
            return Solved {
                q,
                ab: None,
                ok: false,
            };
        };
        let mut size: f64 = 0.0;
        for (i, &k) in free.iter().enumerate() {
            q[k] += step[i];
            size = size.max(step[i].abs());
        }
        if q.iter().any(|x| !x.is_finite() || *x < -0.5 || *x > 1.5) {
            return Solved {
                q,
                ab: None,
                ok: false,
            };
        }
        if size < 1e-15 {
            let ea2 = a.evaluate(q[0], q[1]);
            let eb2 = b.evaluate(q[2], q[3]);
            let ok = norm(sub(ea2.point, eb2.point)) < tolerance * 1e3;
            return Solved {
                q,
                ab: Some((ea2, eb2)),
                ok,
            };
        }
    }
    Solved {
        q,
        ab: None,
        ok: false,
    }
}

/// The point a reports: both surfaces' parameters in their patches' knot domains.
fn record(a: &Piece, b: &Piece, q: [f64; 4], ea: &Eval, eb: &Eval) -> CurvePoint {
    CurvePoint {
        p: scale(
            [
                ea.point[0] + eb.point[0],
                ea.point[1] + eb.point[1],
                ea.point[2] + eb.point[2],
            ],
            0.5,
        ),
        a: Param {
            patch: a.patch,
            u: a.u0 + q[0] * (a.u1 - a.u0),
            v: a.v0 + q[1] * (a.v1 - a.v0),
        },
        b: Param {
            patch: b.patch,
            u: b.u0 + q[2] * (b.u1 - b.u0),
            v: b.v0 + q[3] * (b.v1 - b.v0),
        },
    }
}

/// Whether a boundary curve needs its root search against the other piece:
/// not when the other patch's analytic surface proves it lies on that surface
/// (the arc along it is marched from its ends, as for more than 64 roots) or
/// apart from it (no roots). The proof is the composed field's coefficients
/// on the curve, a one-dimensional `field::relate`; the search it replaces
/// found 65 roots one by one before it knew an edge lay on a seam.
fn edge_needs_search(edge: &[H], other: Option<&Surf>) -> bool {
    let Some(surf) = other else { return true };
    let curve = Piece::new(
        0,
        0.0,
        1.0,
        0.0,
        1.0,
        edge.len() - 1,
        0,
        edge.to_vec(),
        [false; 2],
    );
    relate(&curve, surf) == Relation::Unknown
}

/// The arcs of a loop-free pair, each oriented along nA × nB. `surf_a` and
/// `surf_b` are the pieces' patches' analytic surfaces, where known.
fn trace_leaf(
    a: &Piece,
    b: &Piece,
    surf_a: Option<&Surf>,
    surf_b: Option<&Surf>,
    ctx: &mut Ctx,
) -> Result<Vec<Vec<CurvePoint>>, Error> {
    let mut found = Vec::new();
    for e in 0..4 {
        let edge = a.edge(e);
        let roots = if edge_needs_search(&edge, surf_b) {
            curve_surface(&edge, b, ctx)?
        } else {
            Vec::new()
        };
        for x in roots {
            found.push(match e {
                0 => [0.0, x[0], x[1], x[2]],
                1 => [1.0, x[0], x[1], x[2]],
                2 => [x[0], 0.0, x[1], x[2]],
                _ => [x[0], 1.0, x[1], x[2]],
            });
        }
        let edge = b.edge(e);
        let roots = if edge_needs_search(&edge, surf_a) {
            curve_surface(&edge, a, ctx)?
        } else {
            Vec::new()
        };
        for x in roots {
            found.push(match e {
                0 => [x[1], x[2], 0.0, x[0]],
                1 => [x[1], x[2], 1.0, x[0]],
                2 => [x[1], x[2], x[0], 0.0],
                _ => [x[1], x[2], x[0], 1.0],
            });
        }
    }
    // Seeds closer than a link are one seed.
    let mut seeds: Vec<(State, bool)> = Vec::new();
    for q in found {
        let Some(q) = settle(a, b, q) else { continue };
        let state = both(a, b, q);
        if !seeds
            .iter()
            .any(|(seed, _)| norm(sub(seed.a.point, state.a.point)) < ctx.link)
        {
            seeds.push((state, false));
        }
    }
    let mut arcs = Vec::new();
    for i in 0..seeds.len() {
        if seeds[i].1 {
            continue;
        }
        seeds[i].1 = true;
        let start = seeds[i].0;
        let (direction, along_high, along_any) = inward(&start);
        if direction == 0 {
            continue;
        }
        // A seed started along a boundary may only graze it: zero velocity
        // across two boundaries at once is also a curve tangent to both where
        // they meet (a tiny loop touching a line along which both pieces were
        // split, ssi-corpus c7_tiny_loop_parab_1e-6), and the march from there
        // can fail. Such a seed is dropped then; the arc, if it is one, has
        // other seeds, as before boundary arcs were traced at all.
        let mut points = match march(a, b, start, direction as f64, ctx) {
            Ok(points) => points,
            Err(_) if along_any => continue,
            Err(error) => return Err(error),
        };
        let end = points[points.len() - 1].p;
        for other in seeds.iter_mut() {
            if !other.1 && norm(sub(other.0.a.point, end)) < ctx.link {
                other.1 = true;
            }
        }
        if direction < 0 {
            points.reverse();
        }
        // An arc lying on a piece's high boundary (u = u1 or v = v1) where the
        // surface goes on past it is the low boundary of the piece beyond, which
        // traces it: one owner, because the two sides are subdivided apart, end
        // their copies at different points, and the copies are no twins
        // (ssi-corpus c8n_cylN_seam_line, a line on a closed cylinder's seam).
        // It lies there when its seed starts it along that boundary and every
        // point it records stays on it. Neither alone will do: a short arc that
        // enters and leaves by one edge has only its two ends, both on that edge,
        // and a curve that starts tangent to an edge (a circle touching the edge
        // y = 1 of a plane piece at (0, 1, 0), ssi-corpus c8_sphere_plane_poles)
        // has no velocity across it at first and still leaves it at once.
        let owned_elsewhere = (0..4).any(|k| {
            let above = if k < 2 { a.above[k] } else { b.above[k - 2] };
            along_high[k] && above && points.iter().all(|point| on_high(a, b, point, k))
        });
        if !owned_elsewhere
            && points.len() > 1
            && norm(sub(points[0].p, points[points.len() - 1].p)) > ctx.link
        {
            arcs.push(points);
        }
    }
    Ok(arcs)
}

/// A seed on a collapsed edge (a pole, a cone's apex) is one 3D point at every
/// parameter along that edge, but the curve leaves it along one direction
/// only. The seed is moved along the edge to the parameter whose direction off
/// the edge lines up with the curve's tangent; None when no parameter in the
/// piece does, as the curve then does not enter this piece from the point.
/// From the parameter the seed happened to get, the corrector cannot swing
/// round to the right one (ssi-corpus c8_sphere_plane_poles: a pole seed at
/// longitude 45 for a curve leaving along the 90 meridian stalled).
fn settle(a: &Piece, b: &Piece, mut q: [f64; 4]) -> Option<[f64; 4]> {
    let tangent = both(a, b, q).tangent;
    if !(norm(tangent) > 0.0) {
        return Some(q);
    }
    let direction = unit(tangent);
    for (piece, off) in [(a, 0), (b, 2)] {
        let e = piece.evaluate(q[off], q[off + 1]);
        let (ls, lt) = (norm(e.ds), norm(e.dt));
        let on = |x: f64| x <= MARGIN || x >= 1.0 - MARGIN;
        // Which coordinate runs along the collapsed edge, if any.
        let along = if ls <= 1e-9 * lt && on(q[off + 1]) {
            0
        } else if lt <= 1e-9 * ls && on(q[off]) {
            1
        } else {
            continue;
        };
        let alignment = |x: f64| {
            let (s, t) = if along == 0 {
                (x, q[off + 1])
            } else {
                (q[off], x)
            };
            let e = piece.evaluate(s, t);
            let off_edge = if along == 0 { e.dt } else { e.ds };
            dot(unit(off_edge), direction).abs()
        };
        let samples = 64;
        let (mut best, mut score) = (0.0, -1.0);
        for i in 0..=samples {
            let x = i as f64 / samples as f64;
            let value = alignment(x);
            if value > score {
                (best, score) = (x, value);
            }
        }
        // Golden-section refinement around the best sample.
        let (mut lo, mut hi) = (
            (best - 1.0 / samples as f64).max(0.0),
            (best + 1.0 / samples as f64).min(1.0),
        );
        let g = 0.5 * (5.0_f64.sqrt() - 1.0);
        for _ in 0..60 {
            let (x1, x2) = (hi - g * (hi - lo), lo + g * (hi - lo));
            if alignment(x1) >= alignment(x2) {
                hi = x2;
            } else {
                lo = x1;
            }
        }
        let refined = 0.5 * (lo + hi);
        if alignment(refined) > score {
            (best, score) = (refined, alignment(refined));
        }
        if score < 1.0 - 1e-6 {
            return None;
        }
        q[off + along] = best;
    }
    Some(q)
}

/// Whether a curve point lies on the high boundary of local coordinate k of
/// the pair (k = 0, 1 for a's u and v, 2, 3 for b's), read back from its
/// parameters in the patch's knot domain.
fn on_high(a: &Piece, b: &Piece, point: &CurvePoint, k: usize) -> bool {
    let (value, low, high) = match k {
        0 => (point.a.u, a.u0, a.u1),
        1 => (point.a.v, a.v0, a.v1),
        2 => (point.b.u, b.u0, b.u1),
        _ => (point.b.v, b.v0, b.v1),
    };
    (value - low) / (high - low) >= 1.0 - MARGIN
}

/// +1 or -1 when that sign of nA × nB carries a seed into both pieces, 0 when neither does.
///
/// A seed on two or more boundaries (a piece's corner, or an edge of each
/// piece) may also run along a boundary it is on, its velocity across it zero:
/// that is an arc lying on a seam, a knot line or a patch edge, started from
/// one of its ends. A seed on one boundary only must cross into the pair: one
/// with zero velocity there grazes the boundary in the middle of an arc that
/// other seeds trace whole.
///
/// Also returns, per local coordinate, whether the arc starts along that
/// coordinate's high boundary (u = u1 or v = v1 of its piece), and whether it
/// starts along any boundary at all.
fn inward(state: &State) -> (i32, [bool; 4], bool) {
    let q = state.q;
    let va = velocity(state.a.ds, state.a.dt, state.tangent);
    let vb = velocity(state.b.ds, state.b.dt, state.tangent);
    let v = [va[0], va[1], vb[0], vb[1]];
    let speed = v.iter().fold(f64::NEG_INFINITY, |m, x| m.max(x.abs()));
    let still = 1e-9 * speed;
    let on = q
        .iter()
        .filter(|&&x| x <= MARGIN || x >= 1.0 - MARGIN)
        .count();
    for sign in [1.0, -1.0] {
        let mut fine = true;
        let mut along_high = [false; 4];
        let mut along_any = false;
        for k in 0..4 {
            let along = on >= 2 && v[k].abs() <= still;
            if q[k] <= MARGIN {
                if !(sign * v[k] > still || along) {
                    fine = false;
                }
                along_any |= along;
            }
            if q[k] >= 1.0 - MARGIN {
                if !(sign * v[k] < -still || along) {
                    fine = false;
                }
                along_high[k] = along;
                along_any |= along;
            }
        }
        if fine {
            return (sign as i32, along_high, along_any);
        }
    }
    (0, [false; 4], false)
}

/// March from a seed along sign × (nA × nB) until the arc leaves the pair; the points it passes.
fn march(
    a: &Piece,
    b: &Piece,
    start: State,
    sign: f64,
    ctx: &Ctx,
) -> Result<Vec<CurvePoint>, Error> {
    let mut state = start;
    let mut points = vec![record(a, b, state.q, &state.a, &state.b)];
    // A step longer than the pair leaves it at once; the exit is found exactly either way.
    let mut h = ctx.max_step.min(0.5 * a.size().min(b.size()));
    for _ in 0..100_000 {
        let t = scale(unit(state.tangent), sign);
        let va = velocity(state.a.ds, state.a.dt, scale(t, h));
        let vb = velocity(state.b.ds, state.b.dt, scale(t, h));
        let guess = [
            state.q[0] + va[0],
            state.q[1] + va[1],
            state.q[2] + vb[0],
            state.q[3] + vb[1],
        ];
        let target = [0, 1, 2].map(|k| state.a.point[k] + h * t[k]);
        let next = newton_pair(a, b, guess, Some((target, t)), &[], ctx.tolerance);
        let mut accepted = false;
        if next.ok {
            let after = both(a, b, next.q);
            let turn = dot(unit(after.tangent), unit(state.tangent))
                .clamp(-1.0, 1.0)
                .acos();
            let travelled = norm(sub(after.a.point, state.a.point));
            let measured = travelled < 2.0 * h && travelled > 0.25 * h;
            // A step that leaves the pair across a collapsed edge (a pole) is not
            // judged by its turn: past it the surface's extension folds back and
            // its tangent there points the other way, and the exit is solved
            // exactly on the boundary from inside anyway. Anywhere else a large
            // turn is a jump to another branch, and the step is retried shorter.
            // An exit that does not converge is retried from a shorter step,
            // nearer the boundary: a long step can land far enough out that the
            // solve onto the boundary diverges (ssi-corpus
            // c9n_many_spans_plane_sphere, marched at the default 0.06 rad).
            let exiting = measured
                && outside(&next.q)
                && (turn <= ctx.max_turn || crosses_collapsed(a, b, &next.q));
            if exiting {
                if let Ok((q, ea, eb)) = exit_point(a, b, state.q, next.q, ctx) {
                    points.push(record(a, b, q, &ea, &eb));
                    return Ok(points);
                }
            }
            if turn <= ctx.max_turn && measured && !exiting {
                points.push(record(a, b, after.q, &after.a, &after.b));
                state = after;
                accepted = true;
                if turn < ctx.max_turn / 3.0 {
                    h = (h * 1.5).min(ctx.max_step);
                }
            }
        }
        if !accepted {
            h /= 2.0;
            if h < ctx.resolution * 1e-3 {
                return Err(Error::Stalled(state.a.point));
            }
        }
    }
    Err(Error::Unfinished)
}

/// Whether a point outside the pair lies past a collapsed edge of a piece: an
/// edge where the derivative along it vanishes (a pole, a cone's apex).
fn crosses_collapsed(a: &Piece, b: &Piece, beyond: &[f64; 4]) -> bool {
    (0..4).any(|k| {
        if !(beyond[k] < -MARGIN || beyond[k] > 1.0 + MARGIN) {
            return false;
        }
        let (piece, off) = if k < 2 { (a, 0) } else { (b, 2) };
        let e = piece.evaluate(beyond[off].clamp(0.0, 1.0), beyond[off + 1].clamp(0.0, 1.0));
        // Across a v boundary the edge runs along s, and it is collapsed when Ss vanishes.
        let (along, across) = if k % 2 == 1 {
            (e.ds, e.dt)
        } else {
            (e.dt, e.ds)
        };
        norm(along) <= 1e-9 * norm(across)
    })
}

/// Where the arc from q (inside) to beyond (outside) crosses the pair's boundary.
fn exit_point(
    a: &Piece,
    b: &Piece,
    q: [f64; 4],
    beyond: [f64; 4],
    ctx: &Ctx,
) -> Result<([f64; 4], Eval, Eval), Error> {
    let mut guess = beyond;
    let mut fixed: Vec<usize> = Vec::new();
    for _ in 0..4 {
        let mut first = None;
        let mut fraction = f64::INFINITY;
        for k in 0..4 {
            if fixed.contains(&k) {
                continue;
            }
            let bound = if guess[k] < -MARGIN {
                0.0
            } else if guess[k] > 1.0 + MARGIN {
                1.0
            } else {
                continue;
            };
            let f = (bound - q[k]) / (guess[k] - q[k]);
            if f < fraction {
                fraction = f;
                first = Some(k);
            }
        }
        let Some(first) = first else { break };
        let bound = if guess[first] < 0.0 { 0.0 } else { 1.0 };
        let along = fraction.clamp(0.0, 1.0);
        let mut start = [0, 1, 2, 3].map(|k| q[k] + along * (guess[k] - q[k]));
        start[first] = bound;
        fixed.push(first);
        let solved = newton_pair(a, b, start, None, &fixed, ctx.tolerance);
        if !solved.ok {
            return Err(Error::ExitDiverged);
        }
        guess = solved.q;
        if !outside(&guess) {
            let (ea, eb) = solved.ab.expect("a converged solve carries its points");
            return Ok((guess, ea, eb));
        }
    }
    let q2 = guess.map(|x| x.clamp(0.0, 1.0));
    let state = both(a, b, q2);
    Ok((q2, state.a, state.b))
}

/// Arcs joined end to start into curves; an arc traced twice (along a piece
/// boundary, either way) is kept once. No curve runs through a vertex (where branches
/// cross): curves end there, as the edges of the intersection's graph.
fn join_arcs(arcs: Vec<Vec<CurvePoint>>, vertices: &[V3], ctx: &Ctx) -> Vec<Curve> {
    let at_vertex = |p: V3| vertices.iter().any(|&v| norm(sub(v, p)) < ctx.link);
    let cell = ctx.link * 4.0;
    let key_of = |p: V3| [0, 1, 2].map(|k| (p[k] / cell).floor() as i64);
    let mut unique: Vec<Vec<CurvePoint>> = Vec::new();
    for arc in arcs {
        let (first, last, middle) = (arc[0].p, arc[arc.len() - 1].p, arc[arc.len() >> 1].p);
        // Either way along: an edge both patches' fields trace, one from each
        // side, ran opposite ways, and joined end to start with its copy it was
        // a closed loop of twice its length (ssi-corpus
        // y23_t6_ex6_four_boundary_curves).
        let twin = unique.iter().any(|other| {
            let (o_first, o_last) = (other[0].p, other[other.len() - 1].p);
            let ends = (norm(sub(o_first, first)) < ctx.link && norm(sub(o_last, last)) < ctx.link)
                || (norm(sub(o_first, last)) < ctx.link && norm(sub(o_last, first)) < ctx.link);
            ends && norm(sub(other[other.len() >> 1].p, middle)) < ctx.link * 1e3
        });
        if !twin {
            unique.push(arc);
        }
    }
    let mut starts: HashMap<[i64; 3], Vec<usize>> = HashMap::new();
    for (index, arc) in unique.iter().enumerate() {
        let key = key_of(arc[0].p);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    starts
                        .entry([key[0] + dx, key[1] + dy, key[2] + dz])
                        .or_default()
                        .push(index);
                }
            }
        }
    }
    let mut used = vec![false; unique.len()];
    let following = |index: usize, used: &[bool]| -> Option<usize> {
        let end = unique[index][unique[index].len() - 1].p;
        if at_vertex(end) {
            return None;
        }
        starts
            .get(&key_of(end))?
            .iter()
            .copied()
            .find(|&candidate| {
                !used[candidate] && norm(sub(unique[candidate][0].p, end)) < ctx.link
            })
    };
    let mut curves: Vec<Curve> = Vec::new();
    for first in 0..unique.len() {
        if used[first] {
            continue;
        }
        used[first] = true;
        let mut chain = vec![first];
        let mut next = following(first, &used);
        while let Some(n) = next {
            used[n] = true;
            chain.push(n);
            next = following(n, &used);
        }
        let mut points: Vec<CurvePoint> = Vec::new();
        for index in chain {
            let skip = if points.is_empty() { 0 } else { 1 };
            points.extend_from_slice(&unique[index][skip..]);
        }
        let closed = norm(sub(points[0].p, points[points.len() - 1].p)) < ctx.link;
        if closed {
            points.pop();
        }
        curves.push(Curve { closed, points });
    }
    // An open chain may have been started in its middle: join chains whose ends meet.
    loop {
        let mut joined = None;
        'search: for (i, curve) in curves.iter().enumerate() {
            if curve.closed {
                continue;
            }
            let end = curve.points[curve.points.len() - 1].p;
            if at_vertex(end) {
                continue;
            }
            for (j, other) in curves.iter().enumerate() {
                if j != i && !other.closed && norm(sub(other.points[0].p, end)) < ctx.link {
                    joined = Some((i, j));
                    break 'search;
                }
            }
        }
        let Some((i, j)) = joined else { break };
        let other = curves[j].points[1..].to_vec();
        curves[i].points.extend(other);
        let curve = &mut curves[i];
        if norm(sub(
            curve.points[0].p,
            curve.points[curve.points.len() - 1].p,
        )) < ctx.link
        {
            curve.closed = true;
            curve.points.pop();
        }
        curves.remove(j);
    }
    through_vertices(&mut curves, vertices, ctx.link);
    // An open curve within the resolution of a vertex it touches is that
    // vertex's (a tangent stretch's sliver where two loops meet, five points
    // within 1e-4, ssi-corpus y23_fig1e_loops_touching); one of no length is none.
    curves.retain(|c| {
        let first = c.points[0].p;
        let small = c
            .points
            .iter()
            .all(|q| norm(sub(q.p, first)) < ctx.resolution);
        let touches = [first, c.points[c.points.len() - 1].p]
            .iter()
            .any(|&e| at_vertex(e));
        c.points.iter().any(|q| norm(sub(q.p, first)) >= ctx.link)
            && !(small && !c.closed && touches)
    });
    curves
}

/// The curves ending at a vertex joined through it where they run on straight:
/// of the curve ends at one vertex, the two whose directions away from it are
/// most nearly opposite (under -cos 45 degrees) are one branch, joined and
/// reversed as need be, a curve whose two ends pair closing; most opposite
/// first, so the two branches through a crossing each go on as themselves.
/// Split at every crossing, a small loop crossing a curve twice left two short
/// arcs between the same two vertices a lens's width apart, read as one curve
/// lying on another (ssi-corpus y23_t2_ex3_loop_two_curves); the vertices are
/// reported apart (`Intersection::vertices`).
fn through_vertices(curves: &mut Vec<Curve>, vertices: &[V3], link: f64) {
    // The direction away from a curve's end: to its first point some links off.
    let away = |c: &Curve, at_end: bool| -> V3 {
        let n = c.points.len();
        let end = if at_end {
            c.points[n - 1].p
        } else {
            c.points[0].p
        };
        let next = (1..n)
            .map(|k| {
                if at_end {
                    c.points[n - 1 - k].p
                } else {
                    c.points[k].p
                }
            })
            .find(|&p| norm(sub(p, end)) > 16.0 * link)
            .unwrap_or(if at_end {
                c.points[0].p
            } else {
                c.points[n - 1].p
            });
        let d = sub(next, end);
        let l = norm(d);
        if l > 0.0 {
            [d[0] / l, d[1] / l, d[2] / l]
        } else {
            d
        }
    };
    loop {
        let mut best: Option<(f64, usize, bool, usize, bool)> = None;
        for &v in vertices {
            let ends: Vec<(usize, bool)> = curves
                .iter()
                .enumerate()
                .filter(|(_, c)| !c.closed)
                .flat_map(|(i, c)| {
                    let (first, last) = (c.points[0].p, c.points[c.points.len() - 1].p);
                    [
                        (norm(sub(first, v)) < link).then_some((i, false)),
                        (norm(sub(last, v)) < link).then_some((i, true)),
                    ]
                })
                .flatten()
                .collect();
            for (a, &(i, i_end)) in ends.iter().enumerate() {
                for &(j, j_end) in &ends[a + 1..] {
                    let d = dot(away(&curves[i], i_end), away(&curves[j], j_end));
                    if d < -std::f64::consts::FRAC_1_SQRT_2 && best.is_none_or(|b| d < b.0) {
                        best = Some((d, i, i_end, j, j_end));
                    }
                }
            }
        }
        let Some((_, i, i_end, j, j_end)) = best else {
            break;
        };
        if i == j {
            curves[i].closed = true;
            curves[i].points.pop();
            continue;
        }
        // i's paired end last, j's paired end first.
        let mut head = curves[i].points.clone();
        if !i_end {
            head.reverse();
        }
        let mut tail = curves[j].points.clone();
        if j_end {
            tail.reverse();
        }
        head.extend_from_slice(&tail[1..]);
        let (keep, drop) = (i.min(j), i.max(j));
        curves[keep] = Curve {
            closed: false,
            points: head,
        };
        curves.remove(drop);
    }
}
