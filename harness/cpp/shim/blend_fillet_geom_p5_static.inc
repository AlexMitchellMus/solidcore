// SPDX-License-Identifier: Apache-2.0
// Included at the end of blend_fillet_geom_p1.cpp (after blend_fillet_geom_p4_static.inc), the
// translation unit that compiles upstream SmFilletGeom.cpp: the difftest shim of unit
// `blend.fillet_geom_p5` (SmFilletVertex::CalcFilletIntFillet, SmFilletEdge::CalcCrossSection,
// CalcFilletIntSideFace, CalcEdgeOnExtendedEdge, MakeFilletEdge, SmFilletGeom::MakeRailEdge,
// CreateEndCurve, TrimRailCurves (and TrimCliffRails through it), UpdateSplitFace). Every op
// builds the pseudo B-rep with blend_fillet_geom_p1.cpp's Model and the base solvers of
// blend_fillet_geom_p4_static.inc. Not compiled on its own (blend_fillet_geom_p5.cpp only points
// here).
//
// One entry point, sc_blend_fillet_geom_p5_eval(op, in, n_in, out, cap, out_len). Layouts and
// the refusals (SC_ERR_ARG) are documented in crates/solidcore-ffi/src/units/
// blend_fillet_geom_p5.rs, which validates every input; reads are bounds checked here too.

namespace
{
// Ops (keep in sync with crates/solidcore-ffi/src/units/blend_fillet_geom_p5.rs).
enum : int32_t
{
    P5_INT_FILLET = 0,
    P5_CROSS = 1,
    P5_SIDE = 2,
    P5_EXT_EDGE = 3,
    P5_TRIM = 4,
    P5_SPLIT = 5,
    P5_MAKE_EDGE = 6,
    P5_RAIL_EDGE = 7,
    P5_END = 8,
};

bool p5_is_bspline(const SmCurve* c) { return c != NULL && c->IsKindOf(SmBSplineCurve_TYPE); }

void p5_put_sides(std::vector<double>& o, const Model& M)
{
    for (size_t g = 0; g < M.geoms.size(); ++g)
        o.push_back((double)(M.geoms[g]->*GeomPeek::sides()).GetSize());
}

// CalcCrossSection(ctx, NULL) up to its first unchecked dereference: false when upstream would
// dereference a missing solver or a missing UV curve.
bool p5_cross_runs(const SmContext& ctx, SmFilletEdge* e)
{
    SmFilletVertex* sv = (SmFilletVertex*)e->GetVertex();
    SmFilletVertex* ev = (SmFilletVertex*)e->GetOtherVertex(sv);
    SmFilletEdgeuse* p = (SmFilletEdgeuse*)e->GetPrimaryEdgeuse();
    SmFilletEdgeuse* m = (SmFilletEdgeuse*)p->GetMate();
    SmFilletGeom* g1 = p->GetFilletGeom();
    SmFilletGeom* g2 = m->GetFilletGeom();
    SmSurface* s1 = NULL;
    double tol = 0.0;
    if (g1)
    {
        s1 = g1->GetFilletSurface();
        if (s1 == NULL)
            return true; // NER
        if (g2 != NULL && g2 != g1)
        {
            if (g2->GetFilletSolver() == NULL || g2->GetFilletSurface() == NULL)
                return true; // NER
        }
        if (g1->GetFilletSolver() == NULL)
            return false; // dTol = pFilSolver1->GetThisApproxTol3d()
        tol = g1->GetFilletSolver()->GetThisApproxTol3d();
    }
    else
    {
        SmFace* f = e->GetOriginalFace();
        if (f == NULL || f->GetSurface() == NULL)
            return true; // NER
        s1 = f->GetSurface();
        tol = f->GetTolerance();
    }
    // the same sm_CreateUVCurve call: a success without any curve makes upstream lift NULL
    SmBSplineCurve* c3 = NULL;
    SmBSplineCurve* u1 = NULL;
    SmBSplineCurve* u2 = NULL;
    SmStatus r = sm_CreateUVCurve(ctx, sv, ev, FALSE, s1, 50.0 * tol, c3, u1, u2);
    bool ok = !(r == SM_SUCCESS && c3 == NULL && u1 == NULL);
    delete c3;
    delete u1;
    delete u2;
    return ok;
}

// CalcFilletIntSideFace(ctx, NULL): false when upstream would dereference the NULL corner or
// a missing solver.
bool p5_side_runs(SmFilletEdge* e)
{
    SmFilletVertex* sv = (SmFilletVertex*)e->GetVertex();
    SmFilletVertex* ev = (SmFilletVertex*)e->GetOtherVertex(sv);
    if (!sv->IsProcessed() || !ev->IsProcessed())
        return true;
    SmFilletEdgeuse* p = (SmFilletEdgeuse*)e->GetPrimaryEdgeuse();
    SmFilletGeom* g1 = p->GetFilletGeom();
    if (g1 == NULL)
        return true;
    if (g1->GetFilletSolver() == NULL)
        return false;
    if (!p5_is_bspline(e->GetCurve()) && g1->GetFilletSurface() != NULL)
        return false; // pCorner->IsExtendedSurfacePatch()
    return true;
}

// TrimRailCurves: false when upstream would dereference a missing solver or rail, or reach
// what the harness cannot provide in TrimCliffRails (a cliff vertex to compute, a cliff rail).
bool p5_trim_runs(SmFilletGeom* g)
{
    bool cliff = g->GetFilletGeomType() == SM_FG_CLIFF_ROLLOVER;
    if (cliff && g->GetFilletSolver() == NULL)
        return false; // m_pFilletSolver->GetCreationContext()
    if (g->GetFilletSurface() == NULL)
        return true; // NER
    if (g->GetFilletSolver() == NULL)
        return false;
    for (ULONG ii = 0; ii < 2; ++ii)
        if (g->GetRail(ii) == NULL)
            return false;
    if (cliff)
    {
        const SmTArray<SmFilletEdge*>* lists[2] = {&(g->*GeomPeek::other1()), &(g->*GeomPeek::other2())};
        for (ULONG ii = 0; ii < 2; ++ii)
        {
            std::vector<SmFilletEdge*> all;
            all.push_back(g->GetRail(ii));
            for (ULONG k = 0; k < lists[ii]->GetSize(); ++k)
                all.push_back((*lists[ii])[k]);
            for (size_t k = 0; k < all.size(); ++k)
            {
                SmFilletEdge* r = all[k];
                if (r->GetFilletEdgeType() == SM_FE_CLIFF_RAIL)
                    return false;
                SmFilletVertex* a = (SmFilletVertex*)r->GetVertex();
                SmFilletVertex* b = (SmFilletVertex*)r->GetOtherVertex(a);
                if ((!a->IsProcessed() && a->GetFilletVertexType() == SM_FV_CLIFF_RAIL_X_EDGE) ||
                    (!b->IsProcessed() && b->GetFilletVertexType() == SM_FV_CLIFF_RAIL_X_EDGE))
                    return false;
            }
        }
    }
    return true;
}

void p5_put_rails(std::vector<double>& o, SmFilletGeom* g)
{
    for (ULONG ii = 0; ii < 2; ++ii)
    {
        SmFilletEdge* r = g->GetRail(ii);
        o.push_back(flg(r != NULL));
        if (r == NULL)
            continue;
        put_curve(o, r->GetCurve());
        o.push_back(r->GetInterval().GetMin());
        o.push_back(r->GetInterval().GetMax());
        SmFilletEdgeuse* p = (SmFilletEdgeuse*)r->GetPrimaryEdgeuse();
        SmFilletEdgeuse* m = (SmFilletEdgeuse*)p->GetMate();
        put_curve(o, p->GetUVTrimCurvePointer());
        put_curve(o, m->GetUVTrimCurvePointer());
        put_uv(o, ((SmFilletVertexuse*)p->GetVertexuse())->GetUVPoint());
        put_uv(o, ((SmFilletVertexuse*)m->GetVertexuse())->GetUVPoint());
    }
}

int32_t p5_run(int32_t op, Rd& rd, std::vector<double>& o)
{
    const SmContext& ctx = the_context();
    P4Solvers S;
    Model M(ctx);
    if (!M.build(rd))
        return SC_ERR_ARG;
    if (!p4_read_solvers(rd, ctx, M, S))
        return SC_ERR_ARG;

    switch (op)
    {
    case P5_INT_FILLET:
    {
        SmFilletVertex* v = M.vertex(rd.idx());
        if (rd.bad || v == NULL)
            return SC_ERR_ARG;
        if (v->GetFilletCorner() != NULL)
            return SC_ERR_ARG; // the corner is an opaque tag
        SmStatus r = v->CalcFilletIntFillet();
        o.push_back(st(r));
        p4_put_vertex(o, v);
        return SC_OK;
    }
    case P5_CROSS:
    case P5_SIDE:
    case P5_EXT_EDGE:
    {
        SmFilletEdge* e = M.edge(rd.idx());
        if (rd.bad || e == NULL || e->GetPrimaryEdgeuse() == NULL)
            return SC_ERR_ARG;
        SmStatus r = SM_SUCCESS;
        if (op == P5_CROSS)
        {
            if (!p5_cross_runs(ctx, e))
                return SC_ERR_ARG;
            r = e->CalcCrossSection(ctx, NULL);
            o.push_back(st(r));
            SmFilletEdgeuse* p = (SmFilletEdgeuse*)e->GetPrimaryEdgeuse();
            put_curve(o, e->GetCurve());
            put_curve(o, p->GetUVTrimCurvePointer());
            put_curve(o, p->GetMate()->GetUVTrimCurvePointer());
            p5_put_sides(o, M);
        }
        else if (op == P5_SIDE)
        {
            if (!p5_side_runs(e))
                return SC_ERR_ARG;
            r = e->CalcFilletIntSideFace(ctx, NULL);
            o.push_back(st(r));
            o.push_back((double)e->GetStatus());
            p5_put_sides(o, M);
        }
        else
        {
            r = e->CalcEdgeOnExtendedEdge(ctx, NULL);
            o.push_back(st(r));
        }
        return SC_OK;
    }
    case P5_TRIM:
    {
        SmFilletGeom* g = M.geom(rd.idx());
        if (rd.bad || g == NULL || !p5_trim_runs(g))
            return SC_ERR_ARG;
        SmStatus r = g->TrimRailCurves();
        o.push_back(st(r));
        p5_put_rails(o, g);
        return SC_OK;
    }
    case P5_SPLIT:
    {
        SmFilletGeom* g = M.geom(rd.idx());
        if (rd.bad || g == NULL || g->GetRail(0) == NULL || g->GetRail(1) == NULL)
            return SC_ERR_ARG;
        // SmFilletEdge::UpdateSplitFace -> ShouldSwapFaces dereferences both faces, which
        // this op passes as NULL: only a geometry with no rail edge to update can run (none)
        if (g->GetRail(0) != NULL)
            return SC_ERR_ARG;
        SmStatus r = g->UpdateSplitFace(NULL, NULL);
        o.push_back(st(r));
        return SC_OK;
    }
    case P5_MAKE_EDGE:
    {
        int64_t a = rd.idx();
        int64_t b = rd.idx();
        if (rd.bad)
            return SC_ERR_ARG;
        SmFilletVertex* vs = a < 0 ? NULL : M.vertex(a);
        SmFilletVertex* ve = b < 0 ? NULL : M.vertex(b);
        if ((a >= 0 && vs == NULL) || (b >= 0 && ve == NULL))
            return SC_ERR_ARG;
        SmFilletEdge* ne = NULL;
        SmStatus r = SmFilletEdge::MakeFilletEdge(M.pseudo, vs, ve, ne, NULL);
        o.push_back(st(r));
        o.push_back((double)M.edges().size());
        return SC_OK;
    }
    case P5_RAIL_EDGE:
    {
        SmFilletGeom* g = M.geom(rd.idx());
        int64_t rail = rd.idx();
        if (rd.bad || g == NULL || rail < 0 || rail > 1)
            return SC_ERR_ARG;
        SmFilletEdge* ne = NULL;
        SmStatus r = g->MakeRailEdge(M.pseudo, (ULONG)rail, ne);
        o.push_back(st(r));
        o.push_back((double)M.edges().size());
        return SC_OK;
    }
    case P5_END:
    {
        SmFilletGeom* g = M.geom(rd.idx());
        int64_t end = rd.idx();
        int64_t trim = rd.idx();
        if (rd.bad || g == NULL || end < 0 || end > 1 || trim < 0 || trim > 4)
            return SC_ERR_ARG;
        if (g->GetFilletSurface() != NULL &&
            (g->GetRail(0) == NULL || g->GetRail(1) == NULL || g->GetFilletSolver() == NULL))
            return SC_ERR_ARG; // dereferenced unchecked
        // SM_BT_BLEND runs FindEdgeTangentAtPoint, which casts an edge-classified vertex's
        // classification object to an SmEdge (the model classifies to edge uses)
        if (trim == SM_BT_BLEND)
            return SC_ERR_ARG;
        SmStatus r = g->CreateEndCurve((ULONG)end, (SmBoundaryTrimmingType)trim, M.pseudo);
        o.push_back(st(r));
        o.push_back((double)M.edges().size());
        return SC_OK;
    }
    default:
        return SC_ERR_ARG;
    }
}
} // namespace

SC_EXPORT int32_t sc_blend_fillet_geom_p5_eval(int32_t op, const double* in, int64_t n_in, double* out,
                                               int64_t cap, int64_t* out_len)
{
    if (in == NULL && n_in > 0)
        return SC_ERR_ARG;
    if (out == NULL || out_len == NULL || n_in < 0 || cap < 0)
        return SC_ERR_ARG;
    SC_TRY
    Rd rd(in, n_in);
    // the ops that build surface caches (UB-166 guard in read_surface)
    rd.cache_op = op == P5_CROSS || op == P5_TRIM || op == P5_END || op == P5_SIDE;
    std::vector<double> o;
    int32_t rc = p5_run(op, rd, o);
    if (rc != SC_OK)
        return rc;
    if (rd.at != n_in)
        return SC_ERR_ARG; // trailing input
    return sc_copy_out(o.data(), (int64_t)o.size(), out, cap, out_len);
    SC_CATCH
}
