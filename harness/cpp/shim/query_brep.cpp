// SPDX-License-Identifier: Apache-2.0
// difftest shim for unit `query.brep`: the SM_API surface evaluation entry points
// SmApiEvaluateSurfacePoint, SmApiEvaluateSurfaceNormal and SmApiEvaluateSurfaceDerivatives
// (source/SM_API/src/SmApiSurfaces.cpp), called through the public API.
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
#include <SmSurface.h>
#include <SmTArray.h>
#include <SmVector2d.h>
#include <SmVector3d.h>
#include <nurbs.h>

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

} // namespace sc_query_brep

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
