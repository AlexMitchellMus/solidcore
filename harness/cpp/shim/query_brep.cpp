// SPDX-License-Identifier: Apache-2.0
// difftest shim for unit `query.brep`: the SM_API surface evaluation entry points
// SmApiEvaluateSurfacePoint, SmApiEvaluateSurfaceNormal and SmApiEvaluateSurfaceDerivatives
// (source/SM_API/src/SmApiSurfaces.cpp), called through the public API; the decimation
// attribute constructors and MakeCopy overrides (source/SMLib/inc/SmPolyDecimate.h,
// source/SMLib/src/SmPolyDecimate.cpp), whose protected data is reached through derived
// "peek" classes; and SmBrep::CalculateBoundingBox, CalculateTightBoundingBox,
// LockBoundingBox and UnLockBoundingBox (source/SMLib/src/SmBrep.cpp) on a real SmBrep
// built with the shared topology fixture (sc_topology_fixture::apply: unit
// topology.edgeuse_p1's B-rep block and traversal order, which the Rust side imports from
// that unit's export).
//
// The surface reader is a copy of the surface.surface_cache_p3 shim's (another translation
// unit). Layouts are documented in crates/solidcore-ffi/src/units/query_brep.rs, which
// validates every input.
#include "sc_shim.h"
#include "support/topology_fixture.h"

#include "StdAfx.h"
#include <SmApiGeneral.h>
#include <SmApiSurfaces.h>
#include <SmBSplineSurface.h>
#include <SmBrep.h>
#include <SmCacheMgr.h>
#include <SmContext.h>
#include <SmCurveBoundedSurface.h>
#include <SmExtent3d.h>
#include <SmFace.h>
#include <SmExtent2d.h>
#include <SmNurbsSrf.h>
#include <SmPlane.h>
#include <SmPolyDecimate.h>
#include <SmSurface.h>
#include <SmTArray.h>
#include <SmVector2d.h>
#include <SmVector3d.h>
#include <nurbs.h>

#include <limits>
#include <vector>

#include "sc_sm.h"

using sc::put_vec;

namespace sc_query_brep
{
enum
{
    OP_POINT = 0,
    OP_NORMAL = 1,
    OP_DERIVATIVES = 2,
};

struct Rd
{
    const double* in;
    int64_t at;
    double num() { return in[at++]; }
    int64_t idx() { return (int64_t)in[at++]; }
    bool flag() { return in[at++] != 0.0; }
    SmVector3d vec()
    {
        double x = num();
        double y = num();
        double z = num();
        return SmVector3d(x, y, z);
    }
    SmExtent2d domain()
    {
        double a = num();
        double b = num();
        double c = num();
        double d = num();
        return SmExtent2d(a, b, c, d);
    }
};

gw_SURFACE* read_nurb(Rd& rd)
{
    gw_INDEX n = (gw_INDEX)rd.idx();
    gw_INDEX m = (gw_INDEX)rd.idx();
    gw_DEGREE p = (gw_DEGREE)rd.idx();
    gw_DEGREE q = (gw_DEGREE)rd.idx();
    gw_INDEX r = (gw_INDEX)rd.idx();
    gw_INDEX s = (gw_INDEX)rd.idx();
    gw_SURFACE* sur = sm_AllocateNurbSurface(n, m, p, q, r, s);
    for (gw_INDEX i = 0; i <= n; i++)
        for (gw_INDEX j = 0; j <= m; j++)
        {
            sur->net->Pw[i][j].x = rd.num();
            sur->net->Pw[i][j].y = rd.num();
            sur->net->Pw[i][j].z = rd.num();
            sur->net->Pw[i][j].w = rd.num();
        }
    for (gw_INDEX i = 0; i <= r; i++)
        sur->knu->U[i] = rd.num();
    for (gw_INDEX i = 0; i <= s; i++)
        sur->knv->U[i] = rd.num();
    return sur;
}

// `oob form nurb`.
SmBSplineSurface* read_bspline(Rd& rd, const SmContext& ctx)
{
    bool oob = rd.flag();
    SmBSplineSurfaceForm form = (SmBSplineSurfaceForm)rd.idx();
    SmBSplineSurface* s = new (ctx) SmBSplineSurface(read_nurb(rd), FALSE, NULL);
    s->GetOutOfBoundsEnabled() = oob ? TRUE : FALSE;
    s->SetBSplineSurfaceForm(form);
    return s;
}

// A surface `kind ...`: 0 B-spline `oob form nurb`; 1 bounded plane `o x y sx sy D`;
// 2 curve-bounded surface (implicit outer, no boundary curves) over a B-spline block.
SmSurface* read_surface(Rd& rd, const SmContext& ctx)
{
    int64_t kind = rd.idx();
    if (kind == 1)
    {
        SmPoint3d o = rd.vec();
        SmVector3d x = rd.vec();
        SmVector3d y = rd.vec();
        double sx = rd.num();
        double sy = rd.num();
        SmExtent2d d = rd.domain();
        return new (ctx) SmPlane(o, x, y, SmVector2d(sx, sy), d, &ctx);
    }
    SmBSplineSurface* b = read_bspline(rd, ctx);
    if (kind == 0)
        return b;
    // the curve-bounded surface owns (deletes) its basis
    SmTArray<SmCurve*> c3d(ctx), cuv(ctx);
    return new (ctx) SmCurveBoundedSurface(*b, &c3d, &cuv, TRUE);
}

// Peek classes: the decimation attributes' data is protected (SmPolyDecimate is a friend).
struct DecimatePeek : SmDecimateVertexAttr
{
    ULONG& sort() { return m_lSortIndex; }
    SmPolyVertexClass& vclass() { return m_eVertClass; }
    double& distance() { return m_dVertexDistance; }
    double& error() { return m_dMaxFaceError; }
    SmVector3d& normal() { return m_vAveragePlaneNormal; }
};

struct QuadricVertexPeek : SmQuadricVertexAttr
{
    SmVector3d& q(int i) { return m_Q[i]; }
    SmVector3d& b() { return m_b; }
    double& c() { return m_c; }
};

struct QuadricEdgePeek : SmQuadricEdgeAttr
{
    SmPoint3d& point() { return m_vNewPoint; }
    double& error() { return m_dError; }
    ULONG& sort() { return m_lSortIndex; }
};

enum
{
    KIND_DECIMATE_VERTEX = 0,
    KIND_QUADRIC_VERTEX = 1,
    KIND_QUADRIC_EDGE = 2,
};

// The data of a decimation attribute of kind `kind`; the quadric edge sort index only when
// `sort_set` (it is uninitialised until set), else NaN.
void put_data(std::vector<double>& o, int64_t kind, SmAttribute* a, bool sort_set)
{
    if (kind == KIND_DECIMATE_VERTEX)
    {
        DecimatePeek* d = static_cast<DecimatePeek*>(static_cast<SmDecimateVertexAttr*>(a));
        o.push_back((double)d->sort());
        o.push_back((double)d->vclass());
        o.push_back(d->distance());
        o.push_back(d->error());
        put_vec(o, d->normal());
    }
    else if (kind == KIND_QUADRIC_VERTEX)
    {
        QuadricVertexPeek* d = static_cast<QuadricVertexPeek*>(static_cast<SmQuadricVertexAttr*>(a));
        for (int i = 0; i < 3; i++)
            put_vec(o, d->q(i));
        put_vec(o, d->b());
        o.push_back(d->c());
    }
    else
    {
        QuadricEdgePeek* d = static_cast<QuadricEdgePeek*>(static_cast<SmQuadricEdgeAttr*>(a));
        put_vec(o, d->point());
        o.push_back(d->error());
        o.push_back(sort_set ? (double)d->sort() : std::numeric_limits<double>::quiet_NaN());
    }
}

// Sets the data of a decimation attribute from the record (layout in the ffi module).
void set_data(Rd& rd, int64_t kind, SmAttribute* a, bool& sort_set)
{
    if (kind == KIND_DECIMATE_VERTEX)
    {
        DecimatePeek* d = static_cast<DecimatePeek*>(static_cast<SmDecimateVertexAttr*>(a));
        d->sort() = (ULONG)rd.num();
        d->vclass() = (SmPolyVertexClass)rd.idx();
        d->distance() = rd.num();
        d->error() = rd.num();
        d->normal() = rd.vec();
    }
    else if (kind == KIND_QUADRIC_VERTEX)
    {
        QuadricVertexPeek* d = static_cast<QuadricVertexPeek*>(static_cast<SmQuadricVertexAttr*>(a));
        for (int i = 0; i < 3; i++)
            d->q(i) = rd.vec();
        d->b() = rd.vec();
        d->c() = rd.num();
    }
    else
    {
        QuadricEdgePeek* d = static_cast<QuadricEdgePeek*>(static_cast<SmQuadricEdgeAttr*>(a));
        d->point() = rd.vec();
        d->error() = rd.num();
        sort_set = rd.flag();
        double sort = rd.num();
        if (sort_set)
            d->sort() = (ULONG)sort;
    }
}

void put_attribute(std::vector<double>& o, int64_t kind, SmAttribute* a, bool sort_set)
{
    o.push_back((double)a->GetAttributeID());
    o.push_back((double)a->GetBehavior());
    o.push_back((double)a->GetType());
    put_data(o, kind, a, sort_set);
}

} // namespace sc_query_brep

// `kind id has_behavior behavior set [data]`: constructs the attribute (with the default
// behaviour unless `has_behavior`), sets its data when `set`, copies it with MakeCopy and
// returns `id behavior type data` of the original, then of the copy.
SC_EXPORT int32_t sc_query_brep_attribute(const double* in, int64_t n_in, double* out, int64_t cap,
                                          int64_t* out_len)
{
    using namespace sc_query_brep;
    if (!in || n_in < 5)
        return SC_ERR_ARG;
    SC_TRY
    const SmContext& ctx = *SmApiGetOrCreateContext();
    Rd rd{in, 0};
    int64_t kind = rd.idx();
    ULONG id = (ULONG)rd.num();
    bool has_behavior = rd.flag();
    SmAttributeBehaviorType behavior = (SmAttributeBehaviorType)rd.idx();
    bool set = rd.flag();

    SmAttribute* a;
    if (kind == KIND_DECIMATE_VERTEX)
        a = has_behavior ? new (ctx) SmDecimateVertexAttr(id, behavior) : new (ctx) SmDecimateVertexAttr(id);
    else if (kind == KIND_QUADRIC_VERTEX)
        a = has_behavior ? new (ctx) SmQuadricVertexAttr(id, behavior) : new (ctx) SmQuadricVertexAttr(id);
    else if (kind == KIND_QUADRIC_EDGE)
        a = has_behavior ? new (ctx) SmQuadricEdgeAttr(id, behavior) : new (ctx) SmQuadricEdgeAttr(id);
    else
        return SC_ERR_ARG;

    bool sort_set = false;
    if (set)
        set_data(rd, kind, a, sort_set);
    if (rd.at != n_in)
    {
        delete a;
        return SC_ERR_ARG;
    }

    SmAttribute* c = a->MakeCopy(ctx);
    std::vector<double> o;
    put_attribute(o, kind, a, sort_set);
    // every MakeCopy sets or keeps a defined sort index
    put_attribute(o, kind, c, true);
    delete c;
    delete a;
    return sc_copy_out(o.data(), (int64_t)o.size(), out, cap, out_len);
    SC_CATCH
}

// `has_surface [S] u v P Du Dv`: the op's status, then P, Du and Dv as the entry point
// leaves them (the given sentinels where it does not write them).
SC_EXPORT int32_t sc_query_brep_evaluate(int32_t op, const double* in, int64_t n_in, double* out,
                                         int64_t cap, int64_t* out_len)
{
    using namespace sc_query_brep;
    if (op < OP_POINT || op > OP_DERIVATIVES || !in || n_in < 12)
        return SC_ERR_ARG;
    SC_TRY
    const SmContext& ctx = *SmApiGetOrCreateContext();
    Rd rd{in, 0};
    SmSurface* s = rd.flag() ? read_surface(rd, ctx) : NULL;
    double u = rd.num();
    double v = rd.num();
    SmVector2d uv(u, v);
    SmVector3d p = rd.vec();
    SmVector3d du = rd.vec();
    SmVector3d dv = rd.vec();
    if (rd.at != n_in)
    {
        delete s;
        return SC_ERR_ARG;
    }

    SmApiStatus r;
    if (op == OP_POINT)
        r = SmApiEvaluateSurfacePoint(s, uv, p);
    else if (op == OP_NORMAL)
        r = SmApiEvaluateSurfaceNormal(s, uv, p);
    else
        r = SmApiEvaluateSurfaceDerivatives(s, uv, p, du, dv);
    delete s;

    std::vector<double> o;
    o.push_back(r == SM_SUCCESS ? 0.0 : (double)r);
    put_vec(o, p);
    put_vec(o, du);
    put_vec(o, dv);
    return sc_copy_out(o.data(), (int64_t)o.size(), out, cap, out_len);
    SC_CATCH
}

namespace sc_query_brep
{
enum
{
    STEP_BOX = 0,
    STEP_TIGHT = 1,
    STEP_LOCK = 2,
    STEP_UNLOCK = 3,
};

// Appends the status and `min max` of a box seeded with the runner's sentinel box.
template <typename F> void put_box(std::vector<double>& o, F call)
{
    SmExtent3d box;
    box.SetMinMax(SmPoint3d(-17, -19, -23), SmPoint3d(29, 31, 37));
    SmStatus s = call(box);
    o.push_back(s == SM_SUCCESS ? 0.0 : (double)s);
    SmPoint3d lo = box.GetMin(), hi = box.GetMax();
    o.insert(o.end(), {lo.x, lo.y, lo.z, hi.x, hi.y, hi.z});
}

// Runs `n` steps `kind a b c d` on the B-rep (layouts in the ffi module).
int32_t run_steps(const double* args, int64_t n_args, SmBrep& brep, std::vector<double>& o)
{
    if (n_args < 1)
        return SC_ERR_ARG;
    int64_t n = (int64_t)args[0];
    if (n < 0 || n_args != 1 + 5 * n)
        return SC_ERR_ARG;
    for (int64_t i = 0; i < n; i++)
    {
        const double* a = args + 1 + 5 * i;
        int64_t kind = (int64_t)a[0];
        SmBoolean x = a[1] != 0.0 ? TRUE : FALSE;
        SmBoolean y = a[2] != 0.0 ? TRUE : FALSE;
        SmBoolean z = a[3] != 0.0 ? TRUE : FALSE;
        bool defaults = a[4] != 0.0;
        if (kind == STEP_BOX)
            put_box(o, [&](SmExtent3d& b) {
                return defaults ? brep.CalculateBoundingBox(b) : brep.CalculateBoundingBox(b, x, y, z);
            });
        else if (kind == STEP_TIGHT)
            put_box(o, [&](SmExtent3d& b) { return brep.CalculateTightBoundingBox(b); });
        else if (kind == STEP_LOCK)
        {
            SmStatus s = defaults ? brep.LockBoundingBox() : brep.LockBoundingBox(x);
            o.push_back(s == SM_SUCCESS ? 0.0 : (double)s);
        }
        else if (kind == STEP_UNLOCK)
        {
            brep.UnLockBoundingBox();
            o.push_back(0.0);
        }
        else
            return SC_ERR_ARG;
    }
    return 0;
}
} // namespace sc_query_brep

// Builds the B-rep of `path` / `block` and runs the box steps `args` on it; the output is
// each step's values in order.
SC_EXPORT int32_t sc_query_brep_box(const char* path, const double* block, int64_t n_block,
                                    const double* args, int64_t n_args, double* out, int64_t cap,
                                    int64_t* out_len)
{
    if (path == nullptr || block == nullptr || n_block < 1 || args == nullptr || n_args < 1 ||
        out_len == nullptr)
        return SC_ERR_ARG;
    SC_TRY
    std::vector<double> o;
    int32_t rc = sc_topology_fixture::apply(path, block, n_block, [&](const sc_topology_fixture::View& g) -> int32_t {
        // the trimmed surface caches a file read leaves on face surfaces (the port has
        // no object caches)
        for (SmFace* f : g.faces)
        {
            SmSurface* srf = f->GetSurface();
            if (srf == nullptr)
                return SC_ERR_ARG; // SmFace::CalculateBoundingBox dereferences it unchecked
            if (SmCacheMgr::GetObjectCache(SM_OC_TRIMSRF, srf) != nullptr)
                SmCacheMgr::DeleteObjectsCache(SM_OC_TRIMSRF, srf);
        }
        return sc_query_brep::run_steps(args, n_args, g.brep, o);
    });
    if (rc != 0)
        return rc;
    return sc_copy_out(o.data(), (int64_t)o.size(), out, cap, out_len);
    SC_CATCH
}
