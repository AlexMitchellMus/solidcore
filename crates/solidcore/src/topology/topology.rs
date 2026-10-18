//! The topology object base (`SmTopology`): list links, traversal marks, flags, user
//! indices and the notification step that carries them from an object to its split,
//! merged or copied counterparts; plus topology-change tracking through the
//! context's system notify callback (`SmTrackTopologyChanges`).
//!
//! Upstream: `source/SMLib/inc/SmTopology.h`, `source/SMLib/src/SmTopology.cpp`,
//! `source/SMLib/inc/SmTrackTopologyChanges.h`,
//! `source/SMLib/src/SmTrackTopologyChanges.cpp`.
//!
//! | Rust                                         | Upstream                                         |
//! |----------------------------------------------|--------------------------------------------------|
//! | [`TYPE`], [`TYPE_STRING`], [`CLASS_STRING`]  | `SmTopology::GetType`, `GetTypeString`, `GetClassString` |
//! | [`is_kind_of`]                               | `SmTopology::IsKindOf`                           |
//! | [`Topology`], [`Topology::new`]              | `SmTopology` (its data members), `SmTopology()`  |
//! | [`Topology::owner`]                          | `SmTopology::GetOwner`                           |
//! | [`Topology::next`], [`Topology::last`]       | `SmTopology::GetNext`, `GetLast`                 |
//! | [`Topology::set_list_links`]                 | the friend writes of `SmOwningTopology` (`SetListOwner`, `m_pNext`, `m_pLast`) |
//! | [`Topology::is_connected_to`]                | `SmTopology::IsConnectedTo(const SmTopology*)`   |
//! | [`Topology::all_flags`], [`Topology::set_all_flags`] | `SmTopology::GetAllFlags`, `SetAllFlags` |
//! | [`Topology::set_user_index1`], [`Topology::set_user_index2`], [`Topology::set_user_ptr1`] | `SmTopology::SetUserIndex1`, `SetUserIndex2`, `SetUserPtr1` |
//! | [`Topology::user_index1`], [`Topology::user_index2`], [`Topology::user_ptr1`] | `SmTopology::GetUserIndex1`, `GetUserIndex2`, `GetUserPtr1` |
//! | [`Topology::tolerance`]                      | `SmTopology::GetTolerance` (`SM_USE_OLDTOL`)     |
//! | [`Topology::mark`]                           | `SmTopology::Mark`                               |
//! | [`Topology::get_mark`]                       | `SmTopology::GetMark`                            |
//! | [`Topology::is_marked`]                      | `SmTopology::IsMarked`                           |
//! | [`destroy_topology`](crate::topology::proto_topology::p1::destroy_topology) | `SmTopology::~SmTopology` (unit `topology.proto_topology_p1`) |
//! | [`TopologyRecords`]                          | the `SmTopology*` casts of `Notify`'s arguments  |
//! | [`notify`]                                   | `SmTopology::Notify`                             |
//! | [`NotifyEvent`]                              | `SmNotifyEvent` (data only)                      |
//! | [`TopologyChangeLists`]                      | the arrays `SmTopologyChangeCallback` points at  |
//! | [`TopologyChangeLists::reset_arrays`]        | `SmTopologyChangeCallback::ReSetArrays`, `~SmTopologyChangeCallback` |
//! | [`TopologyChangeCallback`]                   | `SmTopologyChangeCallback` (constructor)         |
//! | [`TrackTopologyChanges`]                     | `SmTrackTopologyChanges`                         |
//! | [`TrackTopologyChanges::new`]                | `SmTrackTopologyChanges::SmTrackTopologyChanges` |
//! | [`TrackTopologyChanges::start_tracking`]     | `SmTrackTopologyChanges::StartTracking`          |
//! | [`TrackTopologyChanges::stop_tracking`]      | `SmTrackTopologyChanges::StopTracking`           |
//! | [`TrackTopologyChanges::is_tracking`]        | `SmTrackTopologyChanges::IsTracking`             |
//! | [`TrackTopologyChanges::rm_list`], [`TrackTopologyChanges::add_list`] | `SmTrackTopologyChanges::GetRmList`, `GetAddList` |
//! | [`TrackTopologyChanges::reset_arrays`]       | `SmTrackTopologyChanges::ReSetArrays`            |
//! | [`TrackTopologyChanges::finish`]             | `SmTrackTopologyChanges::~SmTrackTopologyChanges` |
//!
//! Deviations:
//! - **Object model.** There is no class hierarchy: a topology object is a
//!   [`Topology`] record that the derived kinds (faces, edges, vertices, ...) embed,
//!   addressed by the [`ObjectKey`] it is notified under. Its list owner and
//!   neighbours are keys (`None` for null). `Notify` reaches the records of its data
//!   objects through a [`TopologyRecords`] store; a data object that has no record
//!   (upstream would cast a non-topology object, undefined behaviour) is
//!   [`SolidError::InvalidInput`]. The context the object reaches through
//!   `GetContext` is passed explicitly.
//! - `ULONG` marks and flags are `u64` (LP64). `SetUserPtr1` stores an opaque
//!   `void*`; here it is an opaque `u64` (`0` for null).
//! - `IsMarked` with a null context crashes upstream; [`Topology::is_marked`] takes
//!   the context by reference. `Mark` with a null context does nothing, as upstream.
//! - The reference build uses `SM_USE_OLDTOL`, so `GetTolerance` is the constant
//!   `SM_ZONE_TOL_3D` (overridden by the derived kinds in their units).
//! - `IsConnectedTo` is a "pure virtual" stub upstream (an `SE_MSG` compiled out of
//!   the release build, then `FALSE`); the derived kinds override it.
//! - The destructor unlinks the object from its owner's list
//!   (`SmOwningTopology::Remove`); it needs the owner's list data and is
//!   [`destroy_topology`](crate::topology::proto_topology::p1::destroy_topology).
//! - **Change tracking.** Upstream the tracker owns the Rm/Add/Chg/event arrays and
//!   installs a callback object holding pointers to them; `IsTracking` compares the
//!   context's system callback with the address of that object. Here the tracker
//!   and the boxed [`TopologyChangeCallback`] it installs share one
//!   [`TopologyChangeLists`] (`Arc<Mutex<..>>`, the arrays the callback points at),
//!   so the getters return copies; the callback is tagged with a caller-chosen
//!   tracker key and `IsTracking` recognises it by that key. "Was a callback already installed" (`m_pLastNotifyCallback` non-null) is
//!   kept as a flag. The tracker's B-rep only supplies its context, so the tracker
//!   keeps the [`BrepId`] and takes the context as an argument; the destructor is
//!   the explicit [`TrackTopologyChanges::finish`].
//! - `SmTopologyChangeCallback::Execute` (which fills the lists), `AddUnique` and
//!   `Remove` are reached only through the virtual callback dispatch of
//!   `SmObject::Notify`, so the call graph lists them as unused; they are CO-PORTED
//!   here by `topology.brep_p3` (PR #204) for the healer's trackers
//!   (`SmHealData`), as [`TopologyChangeCallback`]'s [`NotifyCallback::execute`]. An
//!   [`ObjectRef`] carries no owner, so `Execute`'s walk from a notified curve or
//!   surface to `GetOwner()` (and from an edge use to `GetEdge()`) reads an owner
//!   table the tracker's user fills ([`TrackTopologyChanges::set_owner`]); an
//!   object missing from it has a NULL owner, so nothing is recorded for it (as for
//!   a detached curve or surface upstream). The walk is bounded at
//!   `MAX_OWNER_STEPS` (16) steps (a cyclic table). Both `Dump`s and the
//!   `SM_DEBUG_CODE` index arrays are compiled out of the reference build and not
//!   ported.
//! - `AssertValid`, `CheckPointerChain`, `Dump`, the obsolete `NewMark*` / `Mark2`
//!   helpers, `UnMark`, `IsUnMarked`, `IsMarkLocked`, the flag bit helpers and the
//!   placement `operator new` / `delete` are not part of this unit.

use std::collections::BTreeMap;
use std::sync::{Arc, Mutex, MutexGuard};

use crate::error::{Result, SolidError};
use crate::geom::a_object;
use crate::geom::attribute::{AttributeStore, ObjectKey};
use crate::geom::context::{Context, MarkType, NotifyCallback};
use crate::geom::object::{NotifyOperation, ObjectRef};
use crate::geom::tol;
use crate::topology::BrepId;

/// `SmTopology_TYPE` (`TOPO_BASE_TYPE + 1`).
pub const TYPE: i64 = 16001;

/// `SmTopology::GetTypeString`.
pub const TYPE_STRING: &str = "SmTopology_TYPE";

/// `SmTopology::GetClassString`.
pub const CLASS_STRING: &str = "SmTopology";

/// `SmTopology::IsKindOf`: [`TYPE`], then the `SmAObject` chain.
#[must_use]
pub const fn is_kind_of(t: i64) -> bool {
    if t == TYPE {
        true
    } else {
        a_object::is_kind_of(t)
    }
}

/// Index of each mark type's value in [`Topology`]'s mark array (`m_lMark`,
/// `m_lMark2`, `m_lMark3`, `m_lMarkIO`, `m_lMarkAssert`); unknown types use the IO
/// mark, as the `default:` branches upstream.
const fn mark_slot(mark_type: MarkType) -> usize {
    match mark_type {
        MarkType::Mark => 0,
        MarkType::Mark2 => 1,
        MarkType::Mark3 => 2,
        MarkType::MarkAssert => 4,
        MarkType::MarkIo | MarkType::NoMark | MarkType::AllMarks => 3,
    }
}

/// The data every topology object carries (`SmTopology`'s members).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct Topology {
    /// `m_pListOwner`.
    list_owner: Option<ObjectKey>,
    /// `m_pNext`.
    next: Option<ObjectKey>,
    /// `m_pLast`.
    last: Option<ObjectKey>,
    /// `m_lMark`, `m_lMark2`, `m_lMark3`, `m_lMarkIO`, `m_lMarkAssert`.
    marks: [u64; 5],
    /// `m_lFlags`.
    flags: u64,
    /// `m_lUserIndex1`.
    user_index1: u64,
    /// `m_lUserIndex2`.
    user_index2: u64,
    /// `m_pUserPtr1` (opaque; `0` is null).
    user_ptr1: u64,
    /// `m_bIsSmallTopology`.
    small: bool,
}

impl Topology {
    /// `SmTopology()`: unlinked, every mark, flag and user value zero.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            list_owner: None,
            next: None,
            last: None,
            marks: [0; 5],
            flags: 0,
            user_index1: 0,
            user_index2: 0,
            user_ptr1: 0,
            small: false,
        }
    }

    /// The owner of the list this object is on (`GetOwner`).
    #[must_use]
    pub const fn owner(&self) -> Option<ObjectKey> {
        self.list_owner
    }

    /// The next object on the owner's list (`GetNext`).
    #[must_use]
    pub const fn next(&self) -> Option<ObjectKey> {
        self.next
    }

    /// The previous object on the owner's list (`GetLast`).
    #[must_use]
    pub const fn last(&self) -> Option<ObjectKey> {
        self.last
    }

    /// Sets the list owner and neighbours (what `SmOwningTopology`, a friend class,
    /// writes when it links or unlinks the object).
    pub fn set_list_links(
        &mut self,
        owner: Option<ObjectKey>,
        next: Option<ObjectKey>,
        last: Option<ObjectKey>,
    ) {
        self.list_owner = owner;
        self.next = next;
        self.last = last;
    }

    /// `IsConnectedTo(const SmTopology*)` of the base class: `false` (the derived
    /// kinds override it).
    #[must_use]
    pub const fn is_connected_to(&self, _target: Option<ObjectKey>) -> bool {
        false
    }

    /// All flags as one bit array (`GetAllFlags`).
    #[must_use]
    pub const fn all_flags(&self) -> u64 {
        self.flags
    }

    /// Replaces all flags (`SetAllFlags`).
    pub fn set_all_flags(&mut self, flags: u64) {
        self.flags = flags;
    }

    /// The first user index (`GetUserIndex1`).
    #[must_use]
    pub const fn user_index1(&self) -> u64 {
        self.user_index1
    }

    /// The second user index (`GetUserIndex2`).
    #[must_use]
    pub const fn user_index2(&self) -> u64 {
        self.user_index2
    }

    /// The user pointer (`GetUserPtr1`; `0` is null).
    #[must_use]
    pub const fn user_ptr1(&self) -> u64 {
        self.user_ptr1
    }

    /// Sets the first user index (`SetUserIndex1`).
    pub fn set_user_index1(&mut self, index: u64) {
        self.user_index1 = index;
    }

    /// Sets the second user index (`SetUserIndex2`).
    pub fn set_user_index2(&mut self, index: u64) {
        self.user_index2 = index;
    }

    /// Sets the user pointer (`SetUserPtr1`; `0` is null).
    pub fn set_user_ptr1(&mut self, ptr: u64) {
        self.user_ptr1 = ptr;
    }

    /// Whether the object is intended to be tiny (`m_bIsSmallTopology`).
    #[must_use]
    pub const fn is_small_topology(&self) -> bool {
        self.small
    }

    /// The object's zone tolerance (`GetTolerance` of the `SM_USE_OLDTOL` build):
    /// [`tol::ZONE_TOL_3D`] for the base class.
    #[must_use]
    pub const fn tolerance(&self) -> f64 {
        tol::ZONE_TOL_3D
    }

    /// Marks the object with the context's current value of `mark_type` (`Mark`).
    /// Nothing happens for [`MarkType::NoMark`] or a null context;
    /// [`MarkType::AllMarks`] sets the IO mark.
    pub fn mark(&mut self, context: Option<&Context>, mark_type: MarkType) {
        if mark_type != MarkType::NoMark {
            let Some(context) = context else {
                // This should not happen
                return;
            };
            let mark = context.current_mark(mark_type);
            if let Some(m) = self.marks.get_mut(mark_slot(mark_type)) {
                *m = mark;
            }
        }
    }

    /// The object's value of `mark_type` (`GetMark`); [`MarkType::NoMark`] and
    /// [`MarkType::AllMarks`] give the IO mark.
    #[must_use]
    pub fn get_mark(&self, mark_type: MarkType) -> u64 {
        self.marks.get(mark_slot(mark_type)).copied().unwrap_or(0)
    }

    /// Whether the object carries the context's current value of `mark_type`
    /// (`IsMarked`).
    #[must_use]
    pub fn is_marked(&self, context: &Context, mark_type: MarkType) -> bool {
        self.get_mark(mark_type) == context.current_mark(mark_type)
    }
}

/// Where [`notify`] finds the [`Topology`] record of an object.
pub trait TopologyRecords {
    /// The record of `key`, if it is a topology object.
    fn topology(&self, key: ObjectKey) -> Option<&Topology>;
    /// Mutable access to the record of `key`.
    fn topology_mut(&mut self, key: ObjectKey) -> Option<&mut Topology>;
}

impl TopologyRecords for BTreeMap<ObjectKey, Topology> {
    fn topology(&self, key: ObjectKey) -> Option<&Topology> {
        self.get(&key)
    }

    fn topology_mut(&mut self, key: ObjectKey) -> Option<&mut Topology> {
        self.get_mut(&key)
    }
}

/// The error of a notified object that has no [`Topology`] record.
const NOT_TOPOLOGY: SolidError = SolidError::InvalidInput {
    operation: "topology::topology::notify",
    reason: "notified object is not a topology object",
};

fn record<R: TopologyRecords + ?Sized>(records: &R, key: ObjectKey) -> Result<Topology> {
    records.topology(key).copied().ok_or(NOT_TOPOLOGY)
}

fn record_mut<R: TopologyRecords + ?Sized>(
    records: &mut R,
    key: ObjectKey,
) -> Result<&mut Topology> {
    records.topology_mut(key).ok_or(NOT_TOPOLOGY)
}

/// Carries the flags (unless `merge`), the three general marks (zeroed instead when
/// `reset_marks`) and the user values of `src` to `dst`: one target of `Notify`.
fn propagate(src: &Topology, dst: &mut Topology, merge: bool, reset_marks: bool) {
    if !merge {
        dst.set_all_flags(src.all_flags());
    }
    if reset_marks {
        // init marks for next mark usage - don't copy marks
        dst.marks[..3].fill(0);
    } else {
        dst.marks[..3].copy_from_slice(&src.marks[..3]);
    }
    dst.user_index1 = src.user_index1;
    dst.user_index2 = src.user_index2;
    dst.user_ptr1 = src.user_ptr1;
}

/// `SmTopology::Notify`: for splits, merges, copies and coincidences, carries the
/// flags, the three general marks and the user values of the source object to the
/// other objects of the operation, then passes the notification to
/// [`a_object::notify`] (attribute handlers and context callbacks).
///
/// | operation      | source (`this`)  | targets            | kind  |
/// |----------------|------------------|--------------------|-------|
/// | `SplitInBrep`  | `data1`          | `data2`, `data3`   |       |
/// | `MergeInBrep`  | `data2`          | `data1`, `data2`   | merge |
/// | `Copy`         | `this`           | `data1`            | copy  |
/// | `Split`        | `this`           | `data1`, `data2`   |       |
/// | `Merge`        | `this`           | `data1`, `data2`   | merge |
/// | `Coincident`   | `data2`          | `data1`            |       |
/// | `Construction` | `this` if `data2`| `data2`            | copy  |
///
/// A merge first ORs the two targets' flags into both and does not otherwise copy
/// flags; a copy resets the target's general marks to zero instead of copying them.
/// A target equal to the source is skipped.
///
/// # Errors
///
/// [`SolidError::InvalidInput`] when an object the operation reads or writes has no
/// record in `records` (upstream casts it to `SmTopology*` regardless); errors of
/// [`a_object::notify`].
pub fn notify<R: TopologyRecords + ?Sized>(
    records: &mut R,
    attributes: &mut AttributeStore,
    context: Option<&mut Context>,
    operation: NotifyOperation,
    this: ObjectRef,
    data1: Option<ObjectRef>,
    data2: Option<ObjectRef>,
    data3: Option<ObjectRef>,
) -> Result<()> {
    notify_with_regions(
        records, attributes, context, operation, this, data1, data2, data3, None,
    )
}

/// Propagate topology marks, then notify attributes with an actual region lookup.
///
/// # Errors
/// Missing topology records, attribute errors or failed region resolution.
pub fn notify_with_regions<R: TopologyRecords + ?Sized>(
    records: &mut R,
    attributes: &mut AttributeStore,
    context: Option<&mut Context>,
    operation: NotifyOperation,
    this: ObjectRef,
    data1: Option<ObjectRef>,
    data2: Option<ObjectRef>,
    data3: Option<ObjectRef>,
    regions: Option<&dyn crate::topology::shell::RegionLookup>,
) -> Result<()> {
    use NotifyOperation as N;
    let key = |o: Option<ObjectRef>| o.map(|r| r.key);
    let mut copy = false;
    let mut merge = false;

    let (p_this, other1, other2) = match operation {
        N::SplitInBrep => (key(data1), key(data2), key(data3)),
        N::MergeInBrep => {
            merge = true;
            (key(data2), key(data1), key(data2))
        }
        N::Copy => {
            copy = true;
            (Some(this.key), key(data1), None)
        }
        N::Split => (Some(this.key), key(data1), key(data2)),
        N::Merge => {
            merge = true;
            (Some(this.key), key(data1), key(data2))
        }
        N::Coincident => (key(data2), key(data1), None),
        N::Construction => {
            if data2.is_some() {
                copy = true;
                (Some(this.key), key(data2), None)
            } else {
                (None, None, None)
            }
        }
        // SM_NO_UNKNOWN: SE_MSG only (compiled out of the release build)
        _ => (None, None, None),
    };

    // for merge
    if let (Some(o1), Some(o2), true) = (other1, other2, merge) {
        let f2 = record(records, o2)?.all_flags();
        let r1 = record_mut(records, o1)?;
        r1.set_all_flags(r1.all_flags() | f2);
        let f1 = r1.all_flags();
        let r2 = record_mut(records, o2)?;
        r2.set_all_flags(f1 | r2.all_flags());
    }

    // propagate the values to pOther1 (a copy resets its marks), then to pOther2
    for (other, reset_marks) in [(other1, copy), (other2, false)] {
        if let (Some(t), Some(o)) = (p_this, other) {
            if o != t {
                let src = record(records, t)?;
                propagate(&src, record_mut(records, o)?, merge, reset_marks);
            }
        }
    }

    // pass the call along to the base class
    a_object::notify_with_regions(
        attributes, context, operation, this, data1, data2, data3, regions,
    )
}

// ---- change tracking ---------------------------------------------------------------

/// One recorded notification (`SmNotifyEvent`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NotifyEvent {
    /// `m_eAction`.
    pub operation: NotifyOperation,
    /// `m_pCaller`.
    pub caller: Option<ObjectRef>,
    /// `m_pData1`.
    pub data1: Option<ObjectRef>,
    /// `m_pData2`.
    pub data2: Option<ObjectRef>,
    /// `m_pData3`.
    pub data3: Option<ObjectRef>,
}

/// The lists a change tracker accumulates (`m_sRmList`, `m_sAddList`, `m_sChgList`,
/// `m_sNotifyList` of `SmTrackTopologyChanges`).
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct TopologyChangeLists {
    /// Objects removed from the topology graph (stale handles).
    pub rm_list: Vec<ObjectKey>,
    /// Objects added and not removed since.
    pub add_list: Vec<ObjectKey>,
    /// Objects whose geometry changed.
    pub chg_list: Vec<ObjectKey>,
    /// The events that changed the lists.
    pub notify_events: Vec<NotifyEvent>,
}

impl TopologyChangeLists {
    /// Empties all four lists (`SmTopologyChangeCallback::ReSetArrays`, also what
    /// `~SmTopologyChangeCallback` does).
    pub fn reset_arrays(&mut self) {
        self.rm_list.clear();
        self.add_list.clear();
        self.chg_list.clear();
        self.notify_events.clear();
    }
}

/// What a tracker and the callback it installs share: the lists (upstream the
/// callback points at the tracker's arrays) and the owners of notified geometry.
#[derive(Debug, Default)]
struct TrackerState {
    lists: TopologyChangeLists,
    /// `GetOwner` of a notified curve or surface and `GetEdge` of a notified edge
    /// use, which an [`ObjectRef`] does not carry (see the module Deviations).
    owners: BTreeMap<ObjectRef, ObjectRef>,
}

type SharedState = Arc<Mutex<TrackerState>>;

/// The shared state, also after a (panic-free, so never expected) poisoning.
fn state(shared: &SharedState) -> MutexGuard<'_, TrackerState> {
    shared
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner)
}

/// `SmFace_TYPE`, `SmEdge_TYPE`, `SmVertex_TYPE`, `SmEdgeuse_TYPE` (no class derives
/// from the first three, so `IsKindOf` is type equality).
const FACE_TYPE: i64 = 16002;
const EDGE_TYPE: i64 = 16016;
const VERTEX_TYPE: i64 = 16021;
const EDGEUSE_TYPE: i64 = 16005;

/// `IsKindOf(SmCurve_TYPE)`: every curve class is in `CURV_BASE_TYPE`'s block.
const fn is_curve_type(t: i64) -> bool {
    14000 <= t && t < 15000
}

/// `IsKindOf(SmSurface_TYPE)`: every surface class is in `SURF_BASE_TYPE`'s block.
const fn is_surface_type(t: i64) -> bool {
    15000 <= t && t < 16000
}

/// More owner steps than any curve / surface / edge use chain has (a bound on a
/// cyclic owner table).
const MAX_OWNER_STEPS: usize = 16;

/// The system notify callback a [`TrackTopologyChanges`] installs
/// (`SmTopologyChangeCallback`), tagged with its tracker's key.
#[derive(Debug, Clone)]
pub struct TopologyChangeCallback {
    /// Key of the tracker that installed this callback.
    tracker: u64,
    /// The tracker's lists and owner table.
    shared: SharedState,
}

impl TopologyChangeCallback {
    /// A callback of tracker `tracker` with lists of its own
    /// (`SmTopologyChangeCallback(...)`).
    #[must_use]
    pub fn new(tracker: u64) -> Self {
        Self {
            tracker,
            shared: SharedState::default(),
        }
    }

    /// Key of the tracker that installed this callback.
    #[must_use]
    pub const fn tracker(&self) -> u64 {
        self.tracker
    }

    /// The lists this callback has built.
    #[must_use]
    pub fn lists(&self) -> TopologyChangeLists {
        state(&self.shared).lists.clone()
    }
}

/// `SmTopologyChangeCallback::AddUnique`: appends `tgt` to `list` unless it is null or
/// already there; `true` when it was added.
fn add_unique(list: &mut Vec<ObjectKey>, tgt: Option<ObjectRef>) -> bool {
    let Some(tgt) = tgt else {
        return false;
    };
    if list.contains(&tgt.key) {
        return false;
    }
    list.push(tgt.key);
    true
}

/// `SmTopologyChangeCallback::Remove`: removes `tgt` from `list` when it is there;
/// `true` when it was removed.
fn remove(list: &mut Vec<ObjectKey>, tgt: Option<ObjectRef>) -> bool {
    let Some(tgt) = tgt else {
        return false;
    };
    match list.iter().position(|&k| k == tgt.key) {
        Some(i) => {
            list.remove(i);
            true
        }
        None => false,
    }
}

impl NotifyCallback for TopologyChangeCallback {
    /// `SmTopologyChangeCallback::Execute`: run once for every `SmObject::Notify`,
    /// builds the Rm, Add and Chg lists (and the event list) of the faces, edges and
    /// vertices a topology edit changes:
    ///
    /// | operation          | tracked object | lists                                        |
    /// |--------------------|----------------|----------------------------------------------|
    /// | `AddToBrep`        | `data1`        | Add `data1`                                  |
    /// | `SplitInBrep`      | `data1`        | Rm `data1`; un-Add, un-Chg `data1`; Add `data2`, `data3` |
    /// | `MergeInBrep`      | `data1`        | Rm `data1`, `data2`; un-Add, un-Chg both; Add `data1` |
    /// | `TrimNoSplitInBrep`| `data1`        | Rm `data1`; un-Chg; Add `data1`              |
    /// | `RmFromBrep`       | `data1`        | Rm `data1`; un-Add, un-Chg                   |
    /// | `ChangeGeometry`   | `caller`       | Chg it                                       |
    /// | `PostEdit`         | `caller` (a curve or surface: its owner) | itself: Rm and Add; an owner: Chg |
    /// | `Split`, `Merge`   | `data3`'s owner| Chg it                                       |
    /// | `Destruction`      | `data1`        | un-Add, un-Chg `data1`                       |
    ///
    /// Only a tracked face, edge or vertex records anything; `RegPropagation` returns
    /// at once. The status is always success (`SmObject::Notify` ignores it).
    fn execute(
        &mut self,
        operation: NotifyOperation,
        caller: ObjectRef,
        data1: Option<ObjectRef>,
        data2: Option<ObjectRef>,
        data3: Option<ObjectRef>,
    ) -> Result<()> {
        use NotifyOperation as N;
        let event = NotifyEvent {
            operation,
            caller: Some(caller),
            data1,
            data2,
            data3,
        };

        if operation == N::RegPropagation {
            return Ok(());
        }

        let mut st = state(&self.shared);
        let mut tracked = data1;
        if matches!(operation, N::ChangeGeometry | N::PostEdit) {
            tracked = Some(caller);
        } else if matches!(operation, N::Split | N::Merge) {
            tracked = data3;
        }

        if (operation == N::PostEdit
            && (is_surface_type(caller.type_code) || is_curve_type(caller.type_code)))
            || matches!(operation, N::Split | N::Merge)
        {
            for _ in 0..MAX_OWNER_STEPS {
                let Some(t) = tracked else { break };
                if is_surface_type(t.type_code)
                    || is_curve_type(t.type_code)
                    || t.type_code == EDGEUSE_TYPE
                {
                    // GetOwner() / GetEdge()
                    tracked = st.owners.get(&t).copied();
                } else {
                    break;
                }
            }
        }

        let Some(t) = tracked else {
            return Ok(());
        };
        if !matches!(t.type_code, FACE_TYPE | EDGE_TYPE | VERTEX_TYPE) {
            return Ok(());
        }
        let l = &mut st.lists;
        match operation {
            N::AddToBrep => {
                l.notify_events.push(event);
                add_unique(&mut l.add_list, data1);
            }
            N::SplitInBrep => {
                l.notify_events.push(event);
                add_unique(&mut l.rm_list, data1);
                remove(&mut l.add_list, data1);
                remove(&mut l.chg_list, data1);
                add_unique(&mut l.add_list, data2);
                add_unique(&mut l.add_list, data3);
            }
            N::MergeInBrep => {
                l.notify_events.push(event);
                add_unique(&mut l.rm_list, data1);
                add_unique(&mut l.rm_list, data2);
                remove(&mut l.add_list, data1);
                remove(&mut l.add_list, data2);
                remove(&mut l.chg_list, data1);
                remove(&mut l.chg_list, data2);
                add_unique(&mut l.add_list, data1);
            }
            N::TrimNoSplitInBrep => {
                l.notify_events.push(event);
                add_unique(&mut l.rm_list, data1);
                remove(&mut l.chg_list, data1);
                add_unique(&mut l.add_list, data1);
            }
            N::RmFromBrep => {
                l.notify_events.push(event);
                add_unique(&mut l.rm_list, data1);
                remove(&mut l.add_list, data1);
                remove(&mut l.chg_list, data1);
            }
            N::ChangeGeometry | N::Split | N::Merge => {
                l.notify_events.push(event);
                add_unique(&mut l.chg_list, Some(t));
            }
            N::PostEdit => {
                l.notify_events.push(event);
                if t == caller {
                    add_unique(&mut l.rm_list, Some(t));
                    add_unique(&mut l.add_list, Some(t));
                } else {
                    add_unique(&mut l.chg_list, Some(t));
                }
            }
            N::Destruction => {
                remove(&mut l.add_list, data1);
                remove(&mut l.chg_list, data1);
            }
            // SM_NO_COINCIDENT, CHANGE_OWNER, CONSTRUCTION, COPY, PRE_EDIT: nothing
            _ => {}
        }
        Ok(())
    }

    fn as_any(&self) -> Option<&dyn core::any::Any> {
        Some(self)
    }
}

/// Tracks changes to a B-rep's topology graph through its context's system notify
/// callback (`SmTrackTopologyChanges`).
#[derive(Debug)]
pub struct TrackTopologyChanges {
    /// `m_pBrep`: the tracked B-rep (`None` for null).
    brep: Option<BrepId>,
    /// Identifies this tracker's callback in a context (upstream: its address).
    key: u64,
    /// The accumulated lists (shared with the installed callback) and owner table.
    shared: SharedState,
    /// `m_pLastNotifyCallback != NULL`: a system callback was already installed when
    /// tracking last started.
    last_notify_callback: bool,
}

impl TrackTopologyChanges {
    /// `SmTrackTopologyChanges(pBrep, bStartNow)`: a tracker of `brep` whose
    /// callback is recognised by `key` (unique among live trackers of a context).
    /// `context` is the B-rep's context. Starts tracking when `start_now`, else just
    /// empties the lists.
    pub fn new(
        brep: Option<BrepId>,
        key: u64,
        start_now: bool,
        context: Option<&mut Context>,
    ) -> Self {
        let mut t = Self {
            brep,
            key,
            shared: SharedState::default(),
            last_notify_callback: false,
        };
        // when asked start tracking now - else just init the Rm and Add arrays.
        if start_now {
            t.start_tracking(context);
        } else {
            t.reset_arrays();
        }
        t
    }

    /// The tracked B-rep.
    #[must_use]
    pub const fn brep(&self) -> Option<BrepId> {
        self.brep
    }

    /// Whether a system callback was already installed when tracking last started
    /// (`m_pLastNotifyCallback != NULL`); [`TrackTopologyChanges::stop_tracking`]
    /// then leaves the context alone.
    #[must_use]
    pub const fn last_notify_callback(&self) -> bool {
        self.last_notify_callback
    }

    /// The key identifying this tracker's callback.
    #[must_use]
    pub const fn key(&self) -> u64 {
        self.key
    }

    /// Empties the lists and, when the B-rep and its context exist and the context
    /// has no system notify callback, installs this tracker's callback
    /// (`StartTracking`). An already installed callback is left in place (upstream
    /// asserts; nested tracking is not supported).
    pub fn start_tracking(&mut self, context: Option<&mut Context>) {
        // init internal arrays
        self.reset_arrays();

        // when m_pBrep and its Context exist
        if let (Some(_), Some(context)) = (self.brep, context) {
            // start TopoGraph Change Tracking
            self.last_notify_callback = context.sys_notify_callback().is_some();
            if !self.last_notify_callback {
                context.set_sys_notify_callback(Some(Box::new(TopologyChangeCallback {
                    tracker: self.key,
                    shared: Arc::clone(&self.shared),
                })));
            }
        }
    }

    /// Whether the context's system notify callback is this tracker's
    /// (`IsTracking`).
    #[must_use]
    pub fn is_tracking(&self, context: Option<&Context>) -> bool {
        let Some(context) = context.filter(|_| self.brep.is_some()) else {
            return false;
        };
        context
            .sys_notify_callback()
            .and_then(NotifyCallback::as_any)
            .and_then(|a| a.downcast_ref::<TopologyChangeCallback>())
            .is_some_and(|c| c.tracker == self.key)
    }

    /// When the B-rep and its context exist and this tracker found no callback
    /// installed when it started, clears the context's system notify callback
    /// (`StopTracking`); otherwise does nothing (upstream asserts).
    pub fn stop_tracking(&mut self, context: Option<&mut Context>) {
        // when m_pBrep and its Context exist
        if let (Some(_), Some(context)) = (self.brep, context) {
            // when LastNotifyCallback was NULL (it was turned on by this object)
            if !self.last_notify_callback {
                // turn off SmObject::Notify() tracking
                context.set_sys_notify_callback(None);
            }
        }
    }

    /// Objects removed from the topology graph since the lists were last reset
    /// (`GetRmList`).
    #[must_use]
    pub fn rm_list(&self) -> Vec<ObjectKey> {
        state(&self.shared).lists.rm_list.clone()
    }

    /// Objects added since the lists were last reset (`GetAddList`).
    #[must_use]
    pub fn add_list(&self) -> Vec<ObjectKey> {
        state(&self.shared).lists.add_list.clone()
    }

    /// All accumulated lists.
    #[must_use]
    pub fn lists(&self) -> TopologyChangeLists {
        state(&self.shared).lists.clone()
    }

    /// Empties the lists (`ReSetArrays`).
    pub fn reset_arrays(&mut self) {
        state(&self.shared).lists.reset_arrays();
    }

    /// Records `owner` as what `GetOwner` of the curve or surface `geometry` (or
    /// `GetEdge` of the edge use `geometry`) returns while this tracker's callback
    /// follows a notified geometry object to its topology (see the module
    /// Deviations). `None` forgets it (a detached object: upstream's NULL owner).
    pub fn set_owner(&mut self, geometry: ObjectRef, owner: Option<ObjectRef>) {
        let mut st = state(&self.shared);
        match owner {
            Some(o) => {
                st.owners.insert(geometry, o);
            }
            None => {
                st.owners.remove(&geometry);
            }
        }
    }

    /// The destructor `~SmTrackTopologyChanges`: stops tracking if this tracker's
    /// callback is installed.
    pub fn finish(mut self, mut context: Option<&mut Context>) {
        if self.is_tracking(context.as_deref()) {
            self.stop_tracking(context.take());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn r(k: u64) -> ObjectRef {
        ObjectRef {
            type_code: TYPE,
            key: ObjectKey(k),
        }
    }

    fn pool(n: u64) -> BTreeMap<ObjectKey, Topology> {
        (0..n).map(|k| (ObjectKey(k), Topology::new())).collect()
    }

    #[test]
    fn types_and_defaults() {
        assert!(is_kind_of(TYPE));
        assert!(is_kind_of(a_object::TYPE));
        assert!(!is_kind_of(a_object::SAG_TYPE));
        let t = Topology::new();
        assert_eq!(t, Topology::default());
        assert_eq!(t.tolerance(), 1.0e-5);
        assert!(!t.is_connected_to(None));
        assert_eq!((t.owner(), t.next(), t.last()), (None, None, None));
    }

    #[test]
    fn marks_follow_context() {
        let ctx = Context::default();
        let mut t = Topology::new();
        assert!(!t.is_marked(&ctx, MarkType::Mark));
        t.mark(Some(&ctx), MarkType::Mark2);
        assert!(t.is_marked(&ctx, MarkType::Mark2));
        assert_eq!(
            t.get_mark(MarkType::Mark2),
            ctx.current_mark(MarkType::Mark2)
        );
        t.mark(Some(&ctx), MarkType::AllMarks);
        assert_eq!(
            t.get_mark(MarkType::MarkIo),
            ctx.current_mark(MarkType::AllMarks)
        );
        let before = t;
        t.mark(None, MarkType::Mark);
        t.mark(Some(&ctx), MarkType::NoMark);
        assert_eq!(t, before);
    }

    #[test]
    fn notify_split_merge_copy() {
        let mut recs = pool(3);
        let mut attrs = AttributeStore::new();
        let src = recs.get_mut(&ObjectKey(0)).unwrap();
        src.set_all_flags(5);
        src.set_user_index1(7);
        src.marks = [1, 2, 3, 4, 5];
        notify(
            &mut recs,
            &mut attrs,
            None,
            NotifyOperation::Split,
            r(0),
            Some(r(1)),
            Some(r(2)),
            None,
        )
        .unwrap();
        for k in [1, 2] {
            let t = recs[&ObjectKey(k)];
            assert_eq!(t.all_flags(), 5);
            assert_eq!(t.user_index1(), 7);
            assert_eq!(t.marks, [1, 2, 3, 0, 0]);
        }
        recs.get_mut(&ObjectKey(2)).unwrap().set_all_flags(8);
        notify(
            &mut recs,
            &mut attrs,
            None,
            NotifyOperation::Merge,
            r(0),
            Some(r(1)),
            Some(r(2)),
            None,
        )
        .unwrap();
        assert_eq!(recs[&ObjectKey(1)].all_flags(), 13);
        assert_eq!(recs[&ObjectKey(2)].all_flags(), 13);
        notify(
            &mut recs,
            &mut attrs,
            None,
            NotifyOperation::Copy,
            r(0),
            Some(r(1)),
            None,
            None,
        )
        .unwrap();
        assert_eq!(recs[&ObjectKey(1)].all_flags(), 5);
        assert_eq!(recs[&ObjectKey(1)].marks, [0, 0, 0, 0, 0]);
        let missing = notify(
            &mut recs,
            &mut attrs,
            None,
            NotifyOperation::Copy,
            r(0),
            Some(r(9)),
            None,
            None,
        );
        assert!(matches!(missing, Err(SolidError::InvalidInput { .. })));
    }

    #[test]
    fn destroy_unlisted_and_linked() {
        use crate::topology::proto_topology::p1::destroy_topology;
        let mut recs: BTreeMap<ObjectKey, Topology> = BTreeMap::new();
        recs.insert(ObjectKey(0), Topology::new());
        assert!(destroy_topology(&mut recs, ObjectKey(0)).is_ok());
        // an owner without list data (never an SmOwningTopology here)
        recs.insert(ObjectKey(1), Topology::new());
        if let Some(r) = recs.get_mut(&ObjectKey(0)) {
            r.set_list_links(Some(ObjectKey(1)), None, None);
        }
        assert!(matches!(
            destroy_topology(&mut recs, ObjectKey(0)),
            Err(SolidError::InvalidInput { .. })
        ));
    }

    #[test]
    fn tracking_installs_and_removes_callback() {
        let mut ctx = Context::default();
        let t = TrackTopologyChanges::new(Some(BrepId::from_index(0)), 7, true, Some(&mut ctx));
        assert!(t.is_tracking(Some(&ctx)));
        let mut other = TrackTopologyChanges::new(Some(BrepId::from_index(0)), 8, false, None);
        other.start_tracking(Some(&mut ctx));
        assert!(!other.is_tracking(Some(&ctx)));
        other.stop_tracking(Some(&mut ctx));
        assert!(t.is_tracking(Some(&ctx)));
        t.finish(Some(&mut ctx));
        assert!(ctx.sys_notify_callback().is_none());
        let none = TrackTopologyChanges::new(None, 1, true, Some(&mut ctx));
        assert!(!none.is_tracking(Some(&ctx)));
        assert!(none.rm_list().is_empty() && none.add_list().is_empty());
    }

    #[test]
    fn tracker_callback_builds_the_lists() {
        use crate::geom::object::notify as obj_notify;
        use NotifyOperation as N;
        let o = |type_code: i64, k: u64| ObjectRef {
            type_code,
            key: ObjectKey(k),
        };
        let (brep, face, edge, v1, v2) = (
            o(16400, 0),
            o(FACE_TYPE, 1),
            o(EDGE_TYPE, 2),
            o(VERTEX_TYPE, 3),
            o(VERTEX_TYPE, 4),
        );
        let region = o(16012, 5);
        let mut ctx = Context::default();
        let mut t = TrackTopologyChanges::new(Some(BrepId::from_index(0)), 3, true, Some(&mut ctx));
        let n = |ctx: &mut Context, op, c, d1, d2, d3| obj_notify(Some(ctx), op, c, d1, d2, d3);
        // SqueezeEdge's sequence: RM edge, MERGE_IN_BREP(survivor, deleted)
        n(&mut ctx, N::PreEdit, brep, Some(brep), None, None);
        n(&mut ctx, N::RmFromBrep, brep, Some(edge), Some(brep), None);
        n(
            &mut ctx,
            N::MergeInBrep,
            brep,
            Some(v1),
            Some(v2),
            Some(brep),
        );
        // regions are not tracked; a destruction only un-adds
        n(
            &mut ctx,
            N::MergeInBrep,
            brep,
            Some(region),
            Some(region),
            Some(brep),
        );
        n(&mut ctx, N::Destruction, v2, Some(v2), None, None);
        assert_eq!(t.rm_list(), vec![edge.key, v1.key, v2.key]);
        assert_eq!(t.add_list(), vec![v1.key]);
        assert_eq!(t.lists().notify_events.len(), 2);
        // a face's post-edit removes and adds it; a surface's goes to its owner's Chg
        let srf = o(15001, 9);
        n(&mut ctx, N::PostEdit, face, Some(face), None, None);
        n(&mut ctx, N::PostEdit, srf, Some(srf), None, None);
        assert!(t.lists().chg_list.is_empty());
        t.set_owner(srf, Some(face));
        n(&mut ctx, N::PostEdit, srf, Some(srf), None, None);
        assert_eq!(t.lists().chg_list, vec![face.key]);
        assert_eq!(t.rm_list(), vec![edge.key, v1.key, v2.key, face.key]);
        // RegPropagation is ignored, ReSetArrays empties the shared lists
        n(&mut ctx, N::RegPropagation, face, Some(face), None, None);
        t.reset_arrays();
        assert!(t.lists().notify_events.is_empty());
        n(&mut ctx, N::AddToBrep, brep, Some(face), Some(brep), None);
        t.stop_tracking(Some(&mut ctx));
        n(&mut ctx, N::AddToBrep, brep, Some(edge), Some(brep), None);
        assert_eq!(t.add_list(), vec![face.key]);
    }
}
