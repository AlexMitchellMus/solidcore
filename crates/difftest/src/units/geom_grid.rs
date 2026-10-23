//! Runner for unit `geom.grid`: uniform voxel grids and their 3DDDA ray walk.
//!
//! C++: `solidcore_ffi::units::geom_grid` (upstream `source/SMLib/src/SmGrid.cpp`).
//! Rust: `solidcore::geom::grid`.
//!
//! A case is a grid (`SmGrid(bbox, lSize)`) plus a program run against it on both
//! sides: `add` (`AddElement` with poly B-rep or user data), `add_hier`
//! (`AddHierarchyElement` of a base or surface hierarchy element), `find`
//! (`FindVoxelsInBox`), `address` (`GetVoxelAddress`), `voxel_box`
//! (`GetVoxelBBox`), `setup` (`Setup3DDDA`), `step` (`Step3DDDA`), `dump` (every
//! voxel's element and hierarchy lists, `GetType`, `GetElement`, bounding spheres
//! and boxes of the elements, voxel size, voxel tolerance, current mark), `expand`
//! (`ExpandHierarchy`) and `add_surface` (`AddElement` with `SM_GE_SURFACE`).
//!
//! A case may declare real surfaces (`surfaces`, the ray tracer runner's `SurfIn`:
//! a clamped single-patch B-spline or a plane). `add_hier` with `real: s` makes a
//! hierarchy element of surface `s`, `expand` of it runs `ExpandHierarchy` over
//! the surface's cache, and `add_surface { surface, node, bbox }` adds the node
//! `node` of that cache's `GetAllTreeNodes` as a surface element (`-4` when the
//! index is out of range or the node is not an `SM_AD_BEZIER_SURFACE` node). Both
//! sides keep one cache per surface (`smsurf_GetSurfaceCache` returns the surface's
//! one cache, with the polar boxes it has built so far).
//!
//! Outputs: `ints` (statuses, flags, voxel addresses, element type codes and
//! object ids, hierarchy indices, 3DDDA steps; a surface element's id is its
//! surface index) and `floats` (boxes, spheres, ray parameters, the 3DDDA state
//! after every setup and step, and the 51-value leaf record after the sphere and
//! box of a surface element in a dump). Every run has at least [`MIN_GENERATED`]
//! generated cases.

use serde::{Deserialize, Serialize};
use serde_json::Value;

use solidcore::geom::extent2d::Extent2d;
use solidcore::geom::extent3d::Extent3d;
use solidcore::geom::grid::{
    ElementId, ElementSource, Grid, GridElement, GridItem, HierarchyElement, HierarchyId,
    SurfaceGridElement,
};
use solidcore::geom::polar_box::PolarBoxState;
use solidcore::geom::tree::{AuxDataType, NodeId};
use solidcore::geom::vector2d::Vec2;
use solidcore::kernel::Kernel;
use solidcore::math::{Point3, Vec3, BIG_DOUBLE};
use solidcore::surface::plane::Plane;
use solidcore::surface::surface::p2::surface_cache;
use solidcore::surface::surface::Surface;
use solidcore::surface::surface_cache::SurfaceCache;
use solidcore::SolidError;
use solidcore_ffi::units::geom_grid as ffi;

use crate::case::CaseContext;
use crate::runner::default_cases;
use crate::units::intersect_ray_tracer::{rust_surface, SurfIn};
use crate::{Case, HarnessError, Output, Rng, UnitDiff};

/// Generated cases per run, whatever `--seeds` says (it can only raise this).
pub const MIN_GENERATED: usize = 600;

/// `geom.grid` runner.
pub struct Runner;

type P3 = [f64; 3];
type E6 = [f64; 6];

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
enum Kind {
    Poly,
    User,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct GridDesc {
    bbox: E6,
    size: [i64; 3],
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "op", rename_all = "snake_case")]
enum Op {
    Add {
        kind: Kind,
        obj: i64,
        bbox: E6,
    },
    AddHier {
        /// A surface hierarchy element over the fake surface `obj` (it can be dumped, not
        /// expanded) rather than a base-class one.
        surface: bool,
        obj: i64,
        tol: f64,
        bbox: E6,
        /// A surface hierarchy element over the real surface of this index (`surface` and
        /// `obj` are then unused).
        #[serde(default)]
        real: Option<usize>,
    },
    Find {
        bbox: E6,
    },
    Address {
        p: P3,
    },
    VoxelBox {
        a: [i64; 3],
    },
    Setup {
        p: P3,
        v: P3,
    },
    Step {
        t: f64,
    },
    Dump,
    Expand {
        h: i64,
    },
    AddSurface {
        /// The surface whose cache holds the node.
        surface: usize,
        /// Index into the cache tree's `GetAllTreeNodes`.
        node: i64,
        bbox: E6,
    },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct Inputs {
    grid: GridDesc,
    /// The real surfaces of `add_hier { real }` and `add_surface` (none in older cases).
    #[serde(default)]
    surfaces: Vec<SurfIn>,
    ops: Vec<Op>,
}

fn harness(msg: impl Into<String>) -> HarnessError {
    HarnessError::new(msg)
}

fn ext(e: E6) -> Extent3d {
    Extent3d::from_bounds(e[0], e[1], e[2], e[3], e[4], e[5])
}

fn p3(p: P3) -> Point3 {
    Vec3::new(p[0], p[1], p[2])
}

fn put3(o: &mut Vec<f64>, v: Vec3) {
    o.extend([v.x, v.y, v.z]);
}

fn put_ext(o: &mut Vec<f64>, e: &Extent3d) {
    put3(o, e.min());
    put3(o, e.max());
}

/// Upstream status of a Rust error.
fn status_of(e: &SolidError) -> i64 {
    match e {
        SolidError::InvalidInput { .. } => 1007,
        SolidError::Upstream { code, .. } => i64::from(*code),
        _ => 1001,
    }
}

/// The operands of the surface declaration op 11.
fn put_surface_decl(v: &mut Vec<f64>, s: &SurfIn) -> Result<(), HarnessError> {
    v.push(11.0);
    match s {
        SurfIn::Bezier { p, q, net } => {
            if net.len() != p + 1 || net.iter().any(|r| r.len() != q + 1) {
                return Err(harness("surface: control net size"));
            }
            v.extend([0.0, *p as f64, *q as f64]);
            v.extend(net.iter().flatten().flatten());
        }
        SurfIn::Plane {
            o,
            x,
            y,
            scale,
            domain,
        } => {
            v.push(1.0);
            v.extend(o);
            v.extend(x);
            v.extend(y);
            v.extend(scale);
            v.extend(domain);
        }
    }
    Ok(())
}

fn prog(i: &Inputs) -> Result<Vec<f64>, HarnessError> {
    let mut v = vec![0.0];
    v.extend(i.grid.bbox);
    v.extend(i.grid.size.map(|s| s as f64));
    for s in &i.surfaces {
        put_surface_decl(&mut v, s)?;
    }
    for op in &i.ops {
        match op {
            Op::Add { kind, obj, bbox } => {
                v.push(1.0);
                v.push(if *kind == Kind::Poly { 1.0 } else { 3.0 });
                v.push(*obj as f64);
                v.extend(bbox);
            }
            Op::AddHier {
                real: Some(s),
                tol,
                bbox,
                ..
            } => {
                v.extend([12.0, *s as f64, *tol]);
                v.extend(bbox);
            }
            Op::AddHier {
                surface,
                obj,
                tol,
                bbox,
                real: None,
            } => {
                v.extend([2.0, f64::from(u8::from(*surface)), *obj as f64, *tol]);
                v.extend(bbox);
            }
            Op::Find { bbox } => {
                v.push(3.0);
                v.extend(bbox);
            }
            Op::Address { p } => {
                v.push(4.0);
                v.extend(p);
            }
            Op::VoxelBox { a } => {
                v.push(5.0);
                v.extend(a.map(|x| x as f64));
            }
            Op::Setup { p, v: d } => {
                v.push(6.0);
                v.extend(p);
                v.extend(d);
            }
            Op::Step { t } => v.extend([7.0, *t]),
            Op::Dump => v.push(8.0),
            Op::Expand { h } => v.extend([9.0, *h as f64]),
            Op::AddSurface {
                surface,
                node,
                bbox,
            } => {
                v.extend([10.0, *surface as f64, *node as f64]);
                v.extend(bbox);
            }
        }
    }
    Ok(v)
}

/// `(type code, object id)` of a grid element.
fn put_elem(o: &mut Vec<i64>, g: &Grid<i64>, id: ElementId) -> solidcore::Result<()> {
    let e = g.element(id)?;
    o.push(e.type_code());
    o.push(match e.item() {
        GridItem::Polygon(n) => n.0 as i64,
        GridItem::User(obj) => obj,
        // the leaf's surface handle is its index in the case's surface list
        GridItem::Surface(leaf) => leaf.surface,
        GridItem::Hierarchy(_) => -2,
    });
    Ok(())
}

/// The 51-value record of a surface leaf (the shim's `put_leaf`): uv domain (4), box
/// (6), mid point, `du`, `dv` and normal (12), pseudo box basis (9) and intervals (6),
/// polar box state, domain (4) and basis (9) (zeros without a bound).
pub(crate) fn put_leaf<O>(o: &mut Vec<f64>, l: &SurfaceGridElement<O>) {
    let (u0, u1) = (l.uv_domain.min(), l.uv_domain.max());
    o.extend([u0.x, u0.y, u1.x, u1.y]);
    put_ext(o, &l.bbox);
    for v in [l.mid_point, l.mid_du, l.mid_dv, l.mid_normal] {
        put3(o, v);
    }
    for j in 0..3 {
        put3(o, l.pseudo_box.basis_vector(j).unwrap_or(Vec3::ZERO));
    }
    for j in 0..3 {
        let iv = l.pseudo_box.interval(j).unwrap_or_default();
        o.extend([iv.min(), iv.max()]);
    }
    let polar = &l.polar_box;
    let state = match polar.state() {
        PolarBoxState::Uninitialized => 0.0,
        PolarBoxState::PointWithBasis => 1.0,
        PolarBoxState::SinglePoint => 2.0,
        PolarBoxState::Arc => 3.0,
        PolarBoxState::Region => 4.0,
        PolarBoxState::Unbounded => 5.0,
    };
    o.push(state);
    let bounded = (2.0..=4.0).contains(&state);
    let pd = polar.domain();
    o.extend(if bounded {
        [pd.min().x, pd.min().y, pd.max().x, pd.max().y]
    } else {
        [0.0; 4]
    });
    for b in polar.basis() {
        o.extend(if bounded { [b.x, b.y, b.z] } else { [0.0; 3] });
    }
}

fn dump_dda(g: &Grid<i64>, oi: &mut Vec<i64>, of: &mut Vec<f64>) {
    let d = g.dda();
    oi.extend(d.address);
    let next = [d.next_t.x, d.next_t.y, d.next_t.z];
    for (k, n) in next.iter().enumerate() {
        if *n < BIG_DOUBLE {
            oi.push(d.step.get(k).copied().unwrap_or(0));
            oi.push(d.out.get(k).copied().unwrap_or(0));
        }
    }
    put3(of, d.ray_point);
    put3(of, d.ray_vector);
    of.push(d.ray_t);
    put3(of, d.curr_pos);
    put3(of, d.next_t);
    put3(of, d.delta_t);
    put3(of, d.delta_x_vec);
    put3(of, d.delta_y_vec);
    put3(of, d.delta_z_vec);
    for (n, p) in next.iter().zip([d.next_x, d.next_y, d.next_z]) {
        if *n < BIG_DOUBLE {
            put3(of, p);
        }
    }
}

/// Coverage counters of a Rust run (used by the generator test).
#[derive(Debug, Default)]
struct Stats {
    hits: usize,
    multi_hits: usize,
    exits: usize,
    misses: usize,
    errors: usize,
    /// Surface elements in the grids at the end of the runs.
    surface_leaves: usize,
    /// Ray-walk steps that reported a surface element.
    surface_step_hits: usize,
    /// `add_surface` on a Bezier node that succeeded.
    surface_add_ok: usize,
    /// `expand` of a real surface, or `add_surface` of a Bezier node, that did not
    /// return status 0 (the generator keeps its surfaces non-degenerate).
    surface_failures: usize,
}

/// A bad program (the FFI rejects the same inputs).
fn bad_program(reason: &'static str) -> SolidError {
    SolidError::InvalidInput {
        operation: "difftest::geom_grid",
        reason,
    }
}

/// One cache per surface, built on first use (`smsurf_GetSurfaceCache` returns the
/// surface's one cache, whose polar boxes are built lazily). Its handle is the
/// surface's index.
struct Caches {
    surfaces: Vec<Surface>,
    caches: Vec<Option<SurfaceCache<i64>>>,
}

impl Caches {
    fn new(inputs: &[SurfIn]) -> solidcore::Result<Self> {
        let surfaces = inputs
            .iter()
            .map(|s| rust_surface(s).map_err(|_| bad_program("surface")))
            .collect::<solidcore::Result<Vec<_>>>()?;
        let caches = surfaces.iter().map(|_| None).collect();
        Ok(Self { surfaces, caches })
    }

    fn get(&mut self, k: &Kernel, s: usize) -> solidcore::Result<&mut SurfaceCache<i64>> {
        let surface = self.surfaces.get(s).ok_or(bad_program("surface index"))?;
        let slot = self.caches.get_mut(s).ok_or(bad_program("surface index"))?;
        if slot.is_none() {
            *slot = Some(surface_cache::<i64>(k, surface, Some(s as i64))?);
        }
        slot.as_mut().ok_or(bad_program("surface index"))
    }
}

fn run_rust(i: &Inputs, st: &mut Stats) -> solidcore::Result<(Vec<i64>, Vec<f64>)> {
    let k = &Kernel::default();
    let size = i.grid.size.map(|s| s as usize);
    let mut g: Grid<i64> = Grid::new(&ext(i.grid.bbox), size)?;
    let mut caches = Caches::new(&i.surfaces)?;
    // per hierarchy element: the real surface it expands, if any (a base-class element
    // has none and fails before it reads a cache; a fake-surface one cannot expand)
    let mut hier_surface: Vec<(bool, Option<usize>)> = Vec::new();
    let (mut oi, mut of) = (Vec::new(), Vec::new());
    let mut ready = false;
    for op in &i.ops {
        match op {
            Op::Add { kind, obj, bbox } => {
                let src = match kind {
                    Kind::Poly => ElementSource::PolyBrep(NodeId(*obj as usize)),
                    Kind::User => ElementSource::UserData(*obj),
                };
                g.add_element(src, &ext(*bbox))?;
                oi.push(0);
            }
            Op::AddHier {
                surface,
                obj,
                tol,
                bbox,
                real,
            } => {
                let h = match (real, surface) {
                    (Some(s), _) => {
                        if *s >= i.surfaces.len() {
                            return Err(bad_program("surface index"));
                        }
                        HierarchyElement::new_surface(*s as i64, *tol)
                    }
                    (None, true) => HierarchyElement::new_surface(*obj, *tol),
                    (None, false) => HierarchyElement::new(*tol),
                };
                // expandable: a base-class element or a real surface one
                hier_surface.push((real.is_some() || !*surface, *real));
                g.add_hierarchy_element(h, &ext(*bbox))?;
                oi.push(0);
            }
            Op::Find { bbox } => {
                let v = g.find_voxels_in_box(&ext(*bbox));
                oi.push(0);
                oi.push(v.len() as i64);
                for a in v {
                    oi.extend(a.map(|x| x as i64));
                }
            }
            Op::Address { p } => {
                let (a, out) = g.voxel_address(p3(*p));
                oi.extend(a.map(|x| x as i64));
                oi.push(i64::from(out));
            }
            Op::VoxelBox { a } => put_ext(&mut of, &g.voxel_bbox(*a)),
            Op::Setup { p, v } => match g.setup_3ddda(p3(*p), p3(*v)) {
                Ok(miss) => {
                    oi.push(0);
                    oi.push(i64::from(miss));
                    st.misses += usize::from(miss);
                    ready = !miss;
                    if ready {
                        dump_dda(&g, &mut oi, &mut of);
                    }
                }
                Err(e @ SolidError::Unported { .. }) => return Err(e),
                Err(e) => {
                    oi.push(status_of(&e));
                    oi.push(0);
                    st.errors += 1;
                    ready = false;
                }
            },
            Op::Step { t } => {
                if !ready {
                    oi.push(-3);
                    continue;
                }
                match g.step_3ddda(*t) {
                    Ok(h) => {
                        st.hits += usize::from(!h.elements.is_empty());
                        st.multi_hits += usize::from(h.elements.len() > 1);
                        st.exits += usize::from(h.exited_grid);
                        st.surface_step_hits += h
                            .elements
                            .iter()
                            .filter(|&&id| {
                                matches!(g.element(id).map(GridElement::item), Ok(GridItem::Surface(_)))
                            })
                            .count();
                        oi.push(0);
                        oi.push(i64::from(h.exited_grid));
                        oi.push(h.elements.len() as i64);
                        for id in &h.elements {
                            put_elem(&mut oi, &g, *id)?;
                        }
                        of.extend(&h.min_distances);
                        of.extend(&h.max_distances);
                        if h.exited_grid {
                            ready = false;
                        } else {
                            dump_dda(&g, &mut oi, &mut of);
                        }
                    }
                    Err(e @ SolidError::Unported { .. }) => return Err(e),
                    Err(e) => {
                        oi.push(status_of(&e));
                        ready = false;
                    }
                }
            }
            Op::Dump => dump(&g, &mut oi, &mut of)?,
            Op::Expand { h } => {
                let (expandable, real) = usize::try_from(*h)
                    .ok()
                    .and_then(|h| hier_surface.get(h))
                    .copied()
                    .ok_or(bad_program("hierarchy index"))?;
                if !expandable {
                    // the shim would dereference the fake surface
                    return Err(bad_program("expand of a fake-surface element"));
                }
                let id = HierarchyId(*h as usize);
                let result = match real {
                    Some(s) => g.expand_hierarchy(k, id, caches.get(k, s)?),
                    // a base-class element fails before the cache is read
                    None => g.expand_hierarchy(k, id, &mut placeholder_cache(k)?),
                };
                match result {
                    Ok(()) => oi.push(0),
                    Err(e @ SolidError::Unported { .. }) => return Err(e),
                    Err(e) => {
                        st.surface_failures += usize::from(real.is_some());
                        oi.push(status_of(&e));
                    }
                }
            }
            Op::AddSurface {
                surface,
                node,
                bbox,
            } => {
                let cache = caches.get(k, *surface)?;
                let nodes = cache.get_tree(k)?.all_tree_nodes();
                let found = usize::try_from(*node)
                    .ok()
                    .and_then(|n| nodes.get(n))
                    .copied();
                let Some(node) = found else {
                    oi.push(-4);
                    continue;
                };
                let is_bezier = cache.get_tree(k)?.node(node)?.aux_type == AuxDataType::BezierSurface;
                match g.add_surface_node(k, cache, node, &ext(*bbox)) {
                    Ok(_) => {
                        st.surface_add_ok += usize::from(is_bezier);
                        oi.push(0);
                    }
                    Err(e @ SolidError::Unported { .. }) => return Err(e),
                    // the shim answers a node that is not a Bezier patch node itself
                    // (upstream would read garbage); the library guard is this error
                    Err(SolidError::InvalidInput { .. }) if !is_bezier => oi.push(-4),
                    Err(e) => {
                        st.surface_failures += usize::from(is_bezier);
                        oi.push(status_of(&e));
                    }
                }
            }
        }
    }
    st.surface_leaves += g
        .elements()
        .iter()
        .filter(|e| matches!(e.item(), GridItem::Surface(_)))
        .count();
    Ok((oi, of))
}

/// A cache for `expand_hierarchy` on a base-class element, which never reads it.
fn placeholder_cache(k: &Kernel) -> solidcore::Result<SurfaceCache<i64>> {
    let plane = Surface::from(Plane::new(
        Vec3::ZERO,
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec2::new(1.0, 1.0),
        &Extent2d::from_bounds(0.0, 0.0, 1.0, 1.0),
    ));
    surface_cache::<i64>(k, &plane, None)
}

fn dump(g: &Grid<i64>, oi: &mut Vec<i64>, of: &mut Vec<f64>) -> solidcore::Result<()> {
    let [nx, ny, nz] = g.size();
    let mut order: Vec<ElementId> = Vec::new();
    for ix in 0..nx {
        for iy in 0..ny {
            for iz in 0..nz {
                let v = g.voxel([ix, iy, iz])?;
                match v.grid_elements() {
                    None => oi.push(-1),
                    Some(list) => {
                        oi.push(list.len() as i64);
                        for &id in list {
                            put_elem(oi, g, id)?;
                            if !order.contains(&id) {
                                order.push(id);
                            }
                        }
                    }
                }
                match v.hierarchy_elements() {
                    None => oi.push(-1),
                    Some(list) => {
                        oi.push(list.len() as i64);
                        oi.extend(list.iter().map(|h| h.0 as i64));
                    }
                }
            }
        }
    }
    for id in order {
        put_elem(oi, g, id)?;
        let e = g.element(id)?;
        put3(of, e.sphere_center());
        of.push(e.sphere_radius());
        put_ext(of, e.bounding_box());
        if let GridItem::Surface(leaf) = e.item() {
            put_leaf(of, &leaf);
        }
    }
    for h in g.hierarchy_elements() {
        let ge = h.grid_element();
        oi.push(ge.type_code());
        oi.push(match ge.item() {
            GridItem::Hierarchy(Some(id)) => id.0 as i64,
            GridItem::Hierarchy(None) => -1,
            _ => -2,
        });
        put3(of, ge.sphere_center());
        of.push(ge.sphere_radius());
        put_ext(of, ge.bounding_box());
    }
    put3(of, g.voxel_size());
    of.push(g.voxel_tolerance());
    oi.push(g.current_mark() as i64);
    Ok(())
}

// ---- generation -------------------------------------------------------

fn at(a: &[f64], k: usize) -> f64 {
    a.get(k).copied().unwrap_or(0.0)
}

fn set(a: &mut [f64], k: usize, v: f64) {
    if let Some(x) = a.get_mut(k) {
        *x = v;
    }
}

struct Gen {
    lo: P3,
    hi: P3,
    size: [i64; 3],
}

impl Gen {
    fn ext_len(&self, k: usize) -> f64 {
        at(&self.hi, k) - at(&self.lo, k)
    }

    fn diag(&self) -> f64 {
        (0..3)
            .map(|k| self.ext_len(k) * self.ext_len(k))
            .sum::<f64>()
            .sqrt()
    }

    fn center(&self) -> P3 {
        [0, 1, 2].map(|k| 0.5 * (at(&self.lo, k) + at(&self.hi, k)))
    }

    /// A point inside the grid box grown by `grow` (fraction of each side).
    fn point_in(&self, rng: &mut Rng, grow: f64) -> P3 {
        [0, 1, 2].map(|k| {
            let e = self.ext_len(k) * grow;
            rng.uniform(at(&self.lo, k) - e, at(&self.hi, k) + e)
        })
    }

    /// A voxel-lattice point (a corner shared by voxels), possibly on the boundary.
    fn lattice_point(&self, rng: &mut Rng) -> P3 {
        [0, 1, 2].map(|k| {
            let n = self.size.get(k).copied().unwrap_or(1);
            let j = rng.range_i64(0, n);
            at(&self.lo, k) + j as f64 * (self.ext_len(k) / n as f64)
        })
    }

    fn rand_box(&self, rng: &mut Rng) -> E6 {
        let c = self.point_in(rng, 0.2);
        let kind = rng.below(20);
        let mut b = [0.0; 6];
        for k in 0..3 {
            let s = match kind {
                0 => 0.0,
                1 => -rng.uniform(0.0, 0.1) * self.ext_len(k),
                _ => rng.uniform(0.0, 0.4) * self.ext_len(k),
            };
            set(&mut b, k, at(&c, k));
            set(&mut b, k + 3, at(&c, k) + s);
        }
        b
    }
}

fn gen_dir(rng: &mut Rng, g: &Gen, start: P3, boxes: &[E6]) -> P3 {
    if !boxes.is_empty() && rng.chance(0.5) {
        let b = boxes
            .get(rng.below(boxes.len() as u64) as usize)
            .copied()
            .unwrap_or_default();
        let t = [0, 1, 2].map(|k| at(&b, k) + rng.unit() * (at(&b, k + 3) - at(&b, k)));
        let d = [t[0] - start[0], t[1] - start[1], t[2] - start[2]];
        if d.iter().any(|x| *x != 0.0) {
            return d;
        }
    }
    match rng.below(10) {
        0 | 1 => {
            let k = rng.below(3) as usize;
            let mut d = [0.0; 3];
            let sign = if rng.chance(0.5) { 1.0 } else { -1.0 };
            set(&mut d, k, sign * rng.uniform(0.1, 3.0));
            d
        }
        2 => {
            let s = |rng: &mut Rng| if rng.chance(0.5) { 1.0 } else { -1.0 };
            [s(rng), s(rng), s(rng)]
        }
        3 => {
            let k = rng.below(3) as usize;
            let mut d = [
                rng.uniform(-1.0, 1.0),
                rng.uniform(-1.0, 1.0),
                rng.uniform(-1.0, 1.0),
            ];
            set(&mut d, k, 0.0);
            d
        }
        4 if rng.chance(0.3) => [rng.uniform(-1e-7, 1e-7), 0.0, 0.0],
        _ => {
            let t = if rng.chance(0.5) {
                g.point_in(rng, 0.0)
            } else {
                g.lattice_point(rng)
            };
            let d = [t[0] - start[0], t[1] - start[1], t[2] - start[2]];
            if d.iter().all(|x| *x == 0.0) {
                [1.0, 0.5, 0.25]
            } else {
                d
            }
        }
    }
}

/// A random unit vector.
fn unit_dir(rng: &mut Rng) -> P3 {
    for _ in 0..16 {
        let v = [
            rng.uniform(-1.0, 1.0),
            rng.uniform(-1.0, 1.0),
            rng.uniform(-1.0, 1.0),
        ];
        let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        if (0.1..=1.0).contains(&n) {
            return v.map(|c| c / n);
        }
    }
    [1.0, 0.0, 0.0]
}

/// Two independent unit directions.
fn two_dirs(rng: &mut Rng) -> (P3, P3) {
    let x = unit_dir(rng);
    let mut y = unit_dir(rng);
    let dot = x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
    if dot.abs() > 0.8 {
        // nearly parallel: turn `y` to a vector orthogonal to `x` in a coordinate plane
        y = if x[0].abs() < 0.9 {
            [0.0, x[2], -x[1]]
        } else {
            [x[2], 0.0, -x[0]]
        };
        let n = (y[0] * y[0] + y[1] * y[1] + y[2] * y[2]).sqrt();
        y = y.map(|c| c / n);
    }
    (x, y)
}

/// A non-degenerate surface (a clamped single-patch B-spline of degree 1..=3, or a
/// plane) inside the grid box, and the box of its control points or corners grown
/// a little (the hierarchy element's box).
fn gen_grid_surface(rng: &mut Rng, g: &Gen) -> (SurfIn, E6) {
    let c = g.center();
    let side = (0..3)
        .map(|k| g.ext_len(k))
        .fold(f64::INFINITY, f64::min);
    // the half size of the patch: well inside the box
    let a = rng.uniform(0.15, 0.35) * side;
    let c = [0, 1, 2].map(|k| at(&c, k) + rng.uniform(-0.1, 0.1) * g.ext_len(k));
    let (du, dv) = two_dirs(rng);
    let add = |p: P3, d: P3, s: f64| [p[0] + d[0] * s, p[1] + d[1] * s, p[2] + d[2] * s];
    let mut hull: E6 = [f64::INFINITY, f64::INFINITY, f64::INFINITY, -f64::INFINITY, -f64::INFINITY, -f64::INFINITY];
    let mut cover = |p: P3| {
        for k in 0..3 {
            let (lo, hi) = (at(&hull, k), at(&hull, k + 3));
            set(&mut hull, k, lo.min(at(&p, k)));
            set(&mut hull, k + 3, hi.max(at(&p, k)));
        }
    };
    let surf = if rng.chance(0.4) {
        let o = add(add(c, du, -a), dv, -a);
        let s = [2.0 * a, 2.0 * a];
        for (u, v) in [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0)] {
            cover(add(add(o, du, u * s[0]), dv, v * s[1]));
        }
        SurfIn::Plane {
            o,
            x: du,
            y: dv,
            scale: s,
            domain: [0.0, 0.0, 1.0, 1.0],
        }
    } else {
        let p = rng.range_i64(1, 3) as usize;
        let q = rng.range_i64(1, 3) as usize;
        let rational = rng.chance(0.3);
        let net = (0..=p)
            .map(|i| {
                (0..=q)
                    .map(|j| {
                        let s = (i as f64 / p as f64 - 0.5) * 2.0 * a;
                        let t = (j as f64 / q as f64 - 0.5) * 2.0 * a;
                        let n = 0.3 * a;
                        let pt = add(add(c, du, s), dv, t);
                        let pt = [
                            pt[0] + rng.uniform(-n, n),
                            pt[1] + rng.uniform(-n, n),
                            pt[2] + rng.uniform(-n, n),
                        ];
                        cover(pt);
                        let w = if rational { rng.uniform(0.5, 2.0) } else { 1.0 };
                        [pt[0], pt[1], pt[2], w]
                    })
                    .collect()
            })
            .collect();
        SurfIn::Bezier { p, q, net }
    };
    let grow = 0.05 * a;
    let hull = [0, 1, 2, 3, 4, 5].map(|k| at(&hull, k) + if k < 3 { -grow } else { grow });
    (surf, hull)
}

fn gen_case(rng: &mut Rng) -> Inputs {
    let scale = *rng.choose(&[1.0, 10.0, 1000.0, 1e-3]).unwrap_or(&1.0);
    let size = [
        rng.range_i64(1, 8),
        rng.range_i64(1, 8),
        rng.range_i64(1, 8),
    ];
    let (lo, hi) = if rng.chance(0.35) {
        let o = rng.choose(&[0.0, -1.0, 2.0]).copied().unwrap_or(0.0) * scale;
        ([o, o, o], size.map(|n| o + n as f64 * scale))
    } else {
        let c = rng.point3(-scale, scale);
        let e = [0, 1, 2].map(|_| rng.uniform(0.2, 2.0) * scale);
        (c, [c[0] + e[0], c[1] + e[1], c[2] + e[2]])
    };
    let g = Gen { lo, hi, size };
    let mut ops = Vec::new();

    let n_elems = rng.range_i64(0, 16);
    let mut boxes = Vec::new();
    for obj in 0..n_elems {
        let kind = if rng.chance(0.5) {
            Kind::Poly
        } else {
            Kind::User
        };
        let bbox = match boxes.last() {
            // a box overlapping the previous one (several hits in one voxel)
            Some(prev) if rng.chance(0.3) => {
                let p: E6 = *prev;
                let j = [0, 1, 2].map(|k| (at(&p, k + 3) - at(&p, k)) * rng.uniform(-0.5, 0.5));
                [
                    p[0] + j[0],
                    p[1] + j[1],
                    p[2] + j[2],
                    p[3] + j[0],
                    p[4] + j[1],
                    p[5] + j[2],
                ]
            }
            _ => g.rand_box(rng),
        };
        boxes.push(bbox);
        ops.push(Op::Add { kind, obj, bbox });
    }
    let n_hier = rng.range_i64(0, 3);
    for h in 0..n_hier {
        ops.push(Op::AddHier {
            surface: rng.chance(0.5),
            obj: 100 + h,
            tol: rng.uniform(0.0, 0.05) * scale,
            bbox: g.rand_box(rng),
            real: None,
        });
    }
    if n_hier > 0 && rng.chance(0.2) {
        // Only base-class elements: a surface element would defer the case.
        ops.push(Op::Expand { h: -1 });
    }
    // Real surfaces: hierarchy elements over them, expanded before any ray (so the walks
    // meet surface elements), and nodes of their caches added by hand (an index past the
    // tree is `-4` on both sides).
    let mut surfaces = Vec::new();
    // (not at the 1e-3 scale: the surface cache of a random bicubic patch that small
    // grows past 400k tree nodes and does not finish in minutes; see plan/findings.md)
    if rng.chance(0.3) && scale >= 1.0 {
        for s in 0..rng.range_i64(1, 2) {
            let (surf, hull) = gen_grid_surface(rng, &g);
            surfaces.push(surf);
            boxes.push(hull);
            ops.push(Op::AddHier {
                surface: true,
                obj: 0,
                tol: rng.uniform(0.0, 0.02) * scale,
                bbox: hull,
                real: Some(s as usize),
            });
            // the real hierarchy elements follow the base and fake ones
            ops.push(Op::Expand { h: n_hier + s });
        }
        for _ in 0..rng.range_i64(0, 2) {
            let surface = rng.below(surfaces.len() as u64) as usize;
            let node = if rng.chance(0.15) {
                1000
            } else {
                rng.range_i64(0, 4)
            };
            let hull = boxes.last().copied().unwrap_or_default();
            ops.push(Op::AddSurface {
                surface,
                node,
                bbox: hull,
            });
        }
    }
    for _ in 0..rng.range_i64(0, 3) {
        match rng.below(3) {
            0 => ops.push(Op::Find {
                bbox: g.rand_box(rng),
            }),
            1 => {
                let p = if rng.chance(0.3) {
                    g.lattice_point(rng)
                } else {
                    g.point_in(rng, 0.5)
                };
                ops.push(Op::Address { p });
            }
            _ => ops.push(Op::VoxelBox {
                a: size.map(|n| rng.range_i64(-1, n)),
            }),
        }
    }
    if rng.chance(0.3) {
        ops.push(Op::Dump);
    }

    let max_steps = size.iter().sum::<i64>() + 2;
    for _ in 0..rng.range_i64(1, 4) {
        let start = match rng.below(5) {
            0 | 1 => {
                // outside, on a sphere around the grid
                let c = g.center();
                let d = [
                    rng.uniform(-1.0, 1.0),
                    rng.uniform(-1.0, 1.0),
                    rng.uniform(-1.0, 1.0),
                ];
                let r = g.diag() * rng.uniform(0.6, 2.0);
                [c[0] + d[0] * r, c[1] + d[1] * r, c[2] + d[2] * r]
            }
            2 => g.point_in(rng, 0.0),
            3 => g.lattice_point(rng),
            _ => g.point_in(rng, 0.3),
        };
        let v = gen_dir(rng, &g, start, &boxes);
        ops.push(Op::Setup { p: start, v });
        for _ in 0..max_steps {
            let t = if rng.chance(0.8) {
                BIG_DOUBLE
            } else {
                rng.uniform(-0.5, 2.0) * g.diag()
            };
            ops.push(Op::Step { t });
        }
    }
    ops.push(Op::Dump);
    // Resolve the placeholder expand index.
    let mut nh = 0i64;
    let mut base = Vec::new();
    for op in &ops {
        if let Op::AddHier { surface, real, .. } = op {
            if !surface && real.is_none() {
                base.push(nh);
            }
            nh += 1;
        }
    }
    ops.retain_mut(|op| match op {
        Op::Expand { h } if *h < 0 => match base.first() {
            Some(b) => {
                *h = *b;
                true
            }
            None => false,
        },
        _ => true,
    });
    Inputs {
        grid: GridDesc {
            bbox: [lo[0], lo[1], lo[2], hi[0], hi[1], hi[2]],
            size,
        },
        surfaces,
        ops,
    }
}

impl UnitDiff for Runner {
    fn id(&self) -> &'static str {
        "geom.grid"
    }

    fn run_cpp(&self, case: &Case) -> Result<Output, HarnessError> {
        let i: Inputs = case.parse()?;
        let (oi, of) = ffi::run(&prog(&i)?).map_err(|e| harness(e.to_string()))?;
        Ok(Output::new().ints("ints", oi).values("floats", of))
    }

    fn run_rust(&self, case: &Case) -> Result<Output, HarnessError> {
        let i: Inputs = case.parse()?;
        match run_rust(&i, &mut Stats::default()) {
            Ok((oi, of)) => Ok(Output::new().ints("ints", oi).values("floats", of)),
            Err(e) => {
                Ok(Output::unported(&e).unwrap_or_else(|| Output::err(format!("rust_error: {e}"))))
            }
        }
    }

    fn generate(&self, rng: &mut Rng) -> Option<Value> {
        serde_json::to_value(gen_case(rng)).ok()
    }

    fn cases(&self, ctx: &CaseContext) -> Result<Vec<Case>, HarnessError> {
        let ctx = CaseContext {
            seeds: ctx.seeds.max(MIN_GENERATED),
            ..ctx.clone()
        };
        default_cases(self, &ctx)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn generator_covers_walks() {
        let mut st = Stats::default();
        for k in 0..MIN_GENERATED as u64 {
            let mut rng = Rng::new(k);
            let i = gen_case(&mut rng);
            assert!(run_rust(&i, &mut st).is_ok());
        }
        assert!(st.hits > 300, "{st:?}");
        assert!(st.multi_hits > 80, "{st:?}");
        assert!(st.exits > 300, "{st:?}");
        assert!(st.misses > 20, "{st:?}");
        assert!(st.errors > 5, "{st:?}");
    }

    #[test]
    fn generator_covers_surface_elements() {
        let mut st = Stats::default();
        let mut with_surfaces = 0;
        for k in 0..MIN_GENERATED as u64 {
            let mut rng = Rng::new(k);
            let i = gen_case(&mut rng);
            with_surfaces += usize::from(!i.surfaces.is_empty());
            assert!(run_rust(&i, &mut st).is_ok());
        }
        // surface leaves appear in ray-walk hits and in dumps, and no generated
        // `expand` or Bezier `add_surface` fails (the generator keeps its surfaces
        // non-degenerate, as the statuses of `Init` are not propagated upstream)
        assert!(with_surfaces > 100, "{with_surfaces} {st:?}");
        assert!(st.surface_leaves > 100, "{st:?}");
        assert!(st.surface_step_hits > 50, "{st:?}");
        assert!(st.surface_add_ok > 20, "{st:?}");
        assert_eq!(st.surface_failures, 0, "{st:?}");
    }

    fn fixture_inputs(name: &str) -> Inputs {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../fixtures/geom.grid/cases.json");
        let file: Value = serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
        let case = file["cases"]
            .as_array()
            .unwrap()
            .iter()
            .find(|c| c["name"] == name)
            .unwrap();
        serde_json::from_value(case["inputs"].clone()).unwrap()
    }

    /// The two surface fixtures once stayed deferred: both sides must now build the leaves
    /// (status 0) and dump their records, not just agree on an error.
    #[test]
    fn surface_fixtures_dump_leaf_records_on_both_sides() {
        for (name, add_ok) in [
            ("raytracer_surface_hierarchy_expand", 0),
            ("surface_leaf_element", 2),
        ] {
            let i = fixture_inputs(name);
            let mut st = Stats::default();
            let (ri, rf) = run_rust(&i, &mut st).unwrap();
            let (ci, cf) = ffi::run(&prog(&i).unwrap()).unwrap();
            assert!(st.surface_leaves > 0, "{name}: {st:?}");
            assert_eq!(st.surface_failures, 0, "{name}: {st:?}");
            assert_eq!(st.surface_add_ok, add_ok, "{name}: {st:?}");
            assert_eq!(ri, ci, "{name}");
            assert_eq!(rf.len(), cf.len(), "{name}");
            // a leaf record is 51 values; the last dump holds every surface element's
            assert!(rf.len() > st.surface_leaves * 51, "{name}: {}", rf.len());
        }
        // the third `add_surface` of the second fixture names a node past the tree
        let i = fixture_inputs("surface_leaf_element");
        let (ri, _) = run_rust(&i, &mut Stats::default()).unwrap();
        assert_eq!(ri.get(..3), Some(&[0, 0, -4][..]));
    }
}
