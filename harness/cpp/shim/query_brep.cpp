// SPDX-License-Identifier: Apache-2.0
// difftest shim for unit `query.brep`: the SM_API surface evaluation entry points
// SmApiEvaluateSurfacePoint, SmApiEvaluateSurfaceNormal and SmApiEvaluateSurfaceDerivatives
// (source/SM_API/src/SmApiSurfaces.cpp), called through the public API; the decimation
// attribute constructors and MakeCopy overrides (source/SMLib/inc/SmPolyDecimate.h,
// source/SMLib/src/SmPolyDecimate.cpp), whose protected data is reached through derived
// "peek" classes.
//
// The surface reader is a copy of the surface.surface_cache_p3 shim's (another translation
// unit). Layouts are documented in crates/solidcore-ffi/src/units/query_brep.rs, which
// validates every input.
#include "sc_shim.h"

#include "StdAfx.h"
#include <SmApiGeneral.h>
#include <SmApiSurfaces.h>
#include <SmBSplineSurface.h>
#include <SmContext.h>
#include <SmCurveBoundedSurface.h>
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
