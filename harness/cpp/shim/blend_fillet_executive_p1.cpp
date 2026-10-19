// SPDX-License-Identifier: Apache-2.0
// difftest shim for unit `blend.fillet_executive_p1` (the fillet laws). Wraps upstream
// SmFilletLaw / SmLinearFilletLaw / SmBSplineFilletLaw (source/SMLib/src/SmFilletExecutive.cpp
// 5126-6132, source/SMLib/inc/SmFilletExecutive.h 660-899).
//
// Input (flat doubles; validated by crates/solidcore-ffi/src/units/blend_fillet_executive_p1.rs):
//   kind                 0 linear `r0 r1`; 1 B-spline `n p m cpts(x y z w)... knots...`; 2 base
//   has_map e0 e1        the edge map (ignored when has_map == 0 or kind == 2)
//   ncalls               then per call: `t i0 i1 reverse` (parameter, edge interval, orientation)
//   nk kinds[nk]         type codes for IsKindOf
// Output: `0 edge_min edge_max is_linear type_code`, per call `status [v0 v1 v2]`, then
// `has_ext [n p m cpts knots]` (GetExtendedCurve of a B-spline law), then one flag per kind.
#include "sc_shim.h"

#include "StdAfx.h"
#include <SmBSplineCurve.h>
#include <SmContext.h>
#include <SmExtent1d.h>
#include <SmFilletExecutive.h>
#include <SmVector3d.h>
#include <nurbs.h>

#include <vector>

namespace
{
struct Rd
{
    const double* in;
    int64_t at;
    double num() { return in[at++]; }
    int64_t idx() { return (int64_t)in[at++]; }
    bool flag() { return in[at++] != 0.0; }
};

double st(SmStatus s)
{
    return s == SM_SUCCESS ? 0.0 : (double)s;
}

// Reads `n p m cpts knots` into an SmBSplineCurve (SmBSplineCurve(gw_CURVE*, 3, FALSE)).
SmBSplineCurve* read_bcurve(Rd& rd, const SmContext& ctx)
{
    gw_INDEX n = (gw_INDEX)rd.idx();
    gw_DEGREE p = (gw_DEGREE)rd.idx();
    gw_INDEX m = (gw_INDEX)rd.idx();
    gw_CURVE* cur = sm_AllocateNurbCurve(n, p, m);
    for (gw_INDEX i = 0; i <= n; i++)
    {
        cur->pol->Pw[i].x = rd.num();
        cur->pol->Pw[i].y = rd.num();
        cur->pol->Pw[i].z = rd.num();
        cur->pol->Pw[i].w = rd.num();
    }
    for (gw_INDEX i = 0; i <= m; i++)
        cur->knt->U[i] = rd.num();
    return new (ctx) SmBSplineCurve(cur, 3, FALSE);
}

void put_curve(std::vector<double>& o, SmBSplineCurve* c)
{
    gw_CURVE* cur = c ? c->GetOrCreateGwNurbPointer() : NULL;
    if (cur == NULL)
    {
        o.push_back(0.0);
        return;
    }
    o.push_back(1.0);
    o.push_back((double)cur->pol->n);
    o.push_back((double)cur->p);
    o.push_back((double)cur->knt->m);
    for (gw_INDEX i = 0; i <= cur->pol->n; i++)
    {
        o.push_back(cur->pol->Pw[i].x);
        o.push_back(cur->pol->Pw[i].y);
        o.push_back(cur->pol->Pw[i].z);
        o.push_back(cur->pol->Pw[i].w);
    }
    for (gw_INDEX i = 0; i <= cur->knt->m; i++)
        o.push_back(cur->knt->U[i]);
}

} // namespace

SC_EXPORT int32_t sc_blend_fillet_executive_p1_law_eval(const double* in, int64_t n_in, double* out, int64_t cap,
                                                        int64_t* out_len)
{
    (void)n_in;
    SC_TRY
    SmContext ctx;
    std::vector<double> o;
    Rd rd{in, 0};

    int64_t kind = rd.idx();
    SmFilletLaw* law = NULL;
    SmBSplineFilletLaw* bs = NULL;
    if (kind == 0)
    {
        double r0 = rd.num();
        double r1 = rd.num();
        bool has_map = rd.flag();
        double e0 = rd.num();
        double e1 = rd.num();
        SmExtent1d map(e0, e1);
        law = new (ctx) SmLinearFilletLaw(r0, r1, has_map ? &map : NULL);
    }
    else if (kind == 1)
    {
        SmBSplineCurve* lc = read_bcurve(rd, ctx);
        bool has_map = rd.flag();
        double e0 = rd.num();
        double e1 = rd.num();
        bs = new (ctx) SmBSplineFilletLaw(lc, has_map ? SmExtent1d(e0, e1) : SmExtent1d(0.0, 0.0));
        law = bs;
    }
    else
    {
        rd.flag();
        rd.num();
        rd.num();
        law = new (ctx) SmFilletLaw();
    }

    o.push_back(0.0);
    o.push_back(law->GetEdgeMap().GetMin());
    o.push_back(law->GetEdgeMap().GetMax());
    o.push_back(law->IsLinearLaw() ? 1.0 : 0.0);
    o.push_back((double)law->GetType());

    int64_t calls = rd.idx();
    for (int64_t i = 0; i < calls; i++)
    {
        double t = rd.num();
        double i0 = rd.num();
        double i1 = rd.num();
        bool reverse = rd.flag();
        SmExtent1d ivl(i0, i1);
        double v[3] = {0.0, 0.0, 0.0};
        SmStatus s = law->Evaluate(t, ivl, reverse ? TRUE : FALSE, v);
        o.push_back(st(s));
        if (s == SM_SUCCESS)
        {
            o.push_back(v[0]);
            o.push_back(v[1]);
            o.push_back(v[2]);
        }
    }

    put_curve(o, bs ? bs->GetExtendedCurve() : NULL);

    int64_t nk = rd.idx();
    for (int64_t i = 0; i < nk; i++)
        o.push_back(law->IsKindOf((SM_TYPE)rd.idx()) ? 1.0 : 0.0);

    delete law;
    return sc_copy_out(o.data(), (int64_t)o.size(), out, cap, out_len);
    SC_CATCH
}
