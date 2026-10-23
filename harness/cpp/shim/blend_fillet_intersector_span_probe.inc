// SPDX-License-Identifier: Apache-2.0
// Unit-local observations for rv_next_clipped_step. Never included in the pinned
// reference. Virtual overrides only copy arguments/results and delegate once;
// the supplemental span reconstruction starts AFTER ComputeNextPoint returns.
#include <SmCrvOnSurf.h>
#include <SmHermiteCurve.h>

namespace
{
void probe_event(std::vector<double>& events, int tag, const std::vector<double>& values)
{
    events.push_back((double)tag);
    events.push_back((double)values.size());
    events.insert(events.end(), values.begin(), values.end());
}

struct ProbeCra : IxCra
{
    using IxCra::IxCra;
    std::vector<double> events;

    SmStatus PointOnPlaneSolve(const SmPoint3d& origin, const SmVector3d& normal,
                               const SmExtent2d& d0, const SmExtent2d& d1,
                               const SmVector2d& uv0, const SmVector2d& uv1,
                               SmBoolean& found, SmTsectPnt& point) override
    {
        std::vector<double> before;
        put_vec(before, origin);
        put_vec(before, normal);
        put_uv(before, uv0);
        put_uv(before, uv1);
        before.push_back(GetSurface(0)->GetOffsetDistance());
        before.push_back(GetSurface(1)->GetOffsetDistance());
        probe_event(events, 1, before);
        SmStatus result = SmConstantRadiusAssistedFS::PointOnPlaneSolve(
            origin, normal, d0, d1, uv0, uv1, found, point);
        std::vector<double> after{st(result), flg(found)};
        put_point(after, point);
        probe_event(events, 2, after);
        return result;
    }
};

struct ProbePeek : Peek
{
    ProbeCra& solver;
    double observed_step = UNSET;
    ProbePeek(const SmSurface& a, const SmExtent2d& da, const SmSurface& b,
              const SmExtent2d& db, ProbeCra& fs) : Peek(a, da, b, db, fs), solver(fs) {}

    SmStatus ComputePointValues(SmPoint2d uv[2], SmTsectPnt& point,
                                SmTsectPnt* previous, double* step) override
    {
        // In the unchanged fixture this is the first call after StepSolve. Its
        // non-virtual implementation has already returned, including rdStep.
        if (previous == NULL || step == NULL)
            throw Bad();
        observed_step = *step;
        std::vector<double> before{observed_step};
        put_point(before, point);
        probe_event(solver.events, 3, before);
        SmStatus result = SmFilletIntersector::ComputePointValues(uv, point, previous, step);
        std::vector<double> after{st(result)};
        put_point(after, point);
        probe_event(solver.events, 4, after);
        return result;
    }
};

// Supplemental replay of SmFilletIntersector.cpp:652-713, 774-839. The
// concrete Hermite/map evaluators are native, but these are NOT its local objects
// or cache state at their original evaluation times. Emit the residual so the
// test can report whether reconstruction reproduces the genuine operation.
void probe_span(std::vector<double>& events, SmSurface& a, SmSurface& b,
                const SmExtent2d& d0, const SmExtent2d& d1,
                const SmTsectPnt& current, const SmTsectPnt& next, double step)
{
    SmHermiteCurve uv0(current.UVPos(0), current.UVDeriv(0)*step,
                       next.UVPos(0), next.UVDeriv(0)*step, 2);
    SmHermiteCurve uv1(current.UVPos(1), current.UVDeriv(1)*step,
                       next.UVPos(1), next.UVDeriv(1)*step, 2);
    uv0.SetContext(NULL);
    uv1.SetContext(NULL);
    SmCrvOnSurf map0(uv0, a, &d0, 0, NULL);
    SmCrvOnSurf map1(uv1, b, &d1, 0, NULL);
    SmPoint3d p0[2], p1[2], q0[2], q1[2];
    if (map0.Evaluate(map0.GetNaturalInterval().GetMin(), 1, TRUE, p0) != SM_SUCCESS ||
        map0.Evaluate(map0.GetNaturalInterval().GetMax(), 1, TRUE, p1) != SM_SUCCESS)
        throw Bad();
    SmHermiteCurve rail0(p0[0], p0[1], p1[0], p1[1]);
    rail0.SetContext(NULL);
    std::vector<double> ends0;
    for (const SmPoint3d& p : {p0[0], p0[1], p1[0], p1[1]}) put_vec(ends0, p);
    probe_event(events, 11, ends0);
    if (map1.Evaluate(map1.GetNaturalInterval().GetMin(), 1, TRUE, q0) != SM_SUCCESS ||
        map1.Evaluate(map1.GetNaturalInterval().GetMax(), 1, TRUE, q1) != SM_SUCCESS)
        throw Bad();
    SmHermiteCurve rail1(q0[0], q0[1], q1[0], q1[1]);
    rail1.SetContext(NULL);
    std::vector<double> ends1;
    for (const SmPoint3d& p : {q0[0], q0[1], q1[0], q1[1]}) put_vec(ends1, p);
    probe_event(events, 12, ends1);
    const double params[] = {0.5, 0.17, 0.83, 0.33, 0.67};
    double maximum = 0.0;
    for (int i = 0; i < 5; ++i)
    {
        SmPoint3d u0, u1, s0, s1, c0, c1;
        const double t = params[i];
        if (uv0.EvaluatePoint(t, u0) != SM_SUCCESS || uv1.EvaluatePoint(t, u1) != SM_SUCCESS)
            throw Bad();
        SmPoint2d v0 = d0.ClampPoint2d(SmPoint2d(u0.x, u0.y));
        SmPoint2d v1 = d1.ClampPoint2d(SmPoint2d(u1.x, u1.y));
        if (a.EvaluatePoint(v0, s0) != SM_SUCCESS || b.EvaluatePoint(v1, s1) != SM_SUCCESS ||
            rail0.EvaluatePoint(t, c0) != SM_SUCCESS || rail1.EvaluatePoint(t, c1) != SM_SUCCESS)
            throw Bad();
        double ds0 = s0.DistanceBetweenSquared(c0);
        double ds1 = s1.DistanceBetweenSquared(c1);
        double distance = smos_Max(ds0, ds1);
        if (distance > maximum) maximum = distance;
        std::vector<double> sample{t};
        put_vec(sample, u0); put_vec(sample, u1);
        put_uv(sample, v0); put_uv(sample, v1);
        put_vec(sample, s0); put_vec(sample, s1);
        put_vec(sample, c0); put_vec(sample, c1);
        put_vec(sample, s0-c0); put_vec(sample, s1-c1);
        sample.insert(sample.end(), {ds0, ds1, distance, maximum});
        probe_event(events, 20+i, sample);
    }
    probe_event(events, 30, {smos_Sqrt(maximum)});
}
} // namespace

// Same validated program as the original entry point. This companion owns fresh
// geometry and still calls run_op -> genuine ComputeNextPoint, not a copied body.
SC_EXPORT int32_t sc_blend_fillet_intersector_probe_next(const double* in, int64_t n_in,
                                                       double* out, int64_t cap, int64_t* out_len)
{
    if (in == NULL || n_in < 2) return SC_ERR_ARG;
    SC_NATIVE_LOCK;
    try
    {
        SmContext ctx;
        Rd rd{in, n_in, 0};
        SmSurface* a = read_surface(rd, ctx);
        bool same = rd.at < rd.n && rd.in[rd.at] == 2.0;
        if (same) ++rd.at;
        SmSurface* b = same ? a : read_surface(rd, ctx);
        SmExtent2d d0 = rd.dom4(), d1 = rd.dom4();
        double r0 = rd.num(), r1 = rd.num();
        bool o0 = rd.flag(), o1 = rd.flag();
        double ftol = rd.num(), fang = rd.num(), ftan = rd.num();
        double ang = rd.num(), tol = rd.num();
        bool clipped = rd.flag();
        double dir = rd.num();
        if (clipped) throw Bad();
        Pts pts;
        int64_t np = rd.idx(0, 16);
        for (int64_t i = 0; i < np; ++i)
        {
            SmTsectPnt* point = new SmTsectPnt();
            pts.v.push_back(point);
            read_point(rd, *point);
        }
        SmBrep* box = NULL;
        if (SmApiCreateBox(SmVector3d(0.0, 0.0, 0.0), 1.0, 1.0, 1.0, box) != SM_SUCCESS || box == NULL)
            throw Bad();
        std::vector<double> result;
        {
            SmTArray<SmEdge*> edges;
            box->GetEdges(edges);
            if (edges.GetSize() == 0 || edges[0]->GetPrimaryEdgeuse() == NULL) throw Bad();
            ProbeCra fs(ctx, ftol, fang, ftan, r0, r1, edges[0]->GetPrimaryEdgeuse());
            SmOffsetSurface* s0 = new (ctx) SmOffsetSurface(r0*(o0 ? 1.0 : -1.0), *a);
            SmOffsetSurface* s1 = new (ctx) SmOffsetSurface(r1*(o1 ? 1.0 : -1.0), *b);
            s0->SetAllowExtension(TRUE); s1->SetAllowExtension(TRUE);
            fs.setup(s0, s1, r0, r1, o0 ? 1.0 : -1.0, o1 ? 1.0 : -1.0);
            ProbePeek si(*a, d0, *b, d1, fs);
            si.set(ang, tol, false, dir);
            if (rd.idx(0, OP_LAST) != OP_NEXT) throw Bad();
            // Peek at the validated point indices, but leave the reader for run_op.
            Rd arguments = rd;
            SmTsectPnt* current = pts.pt(arguments);
            arguments.num();
            SmTsectPnt* next = pts.pt(arguments);
            std::vector<double> original{(double)OP_NEXT};
            run_op(OP_NEXT, rd, pts, si, ctx, original);
            si.put_state(original, fs);
            original.push_back((double)pts.v.size());
            for (SmTsectPnt* point : pts.v) put_point(original, *point);
            result.push_back((double)original.size());
            result.insert(result.end(), original.begin(), original.end());
            if (si.observed_step == UNSET) throw Bad();
            // All extra evaluations occur after the original record is frozen.
            probe_span(fs.events, *a, *b, d0, d1, *current, *next, si.observed_step);
            result.insert(result.end(), fs.events.begin(), fs.events.end());
            for (SmTsectPnt* point : pts.v) point->ReSet();
        }
        delete box;
        if (!same) delete b;
        delete a;
        return sc_copy_out(result.data(), (int64_t)result.size(), out, cap, out_len);
    }
    catch (const Bad&) { return SC_ERR_ARG; }
    catch (...) { return SC_ERR_EXCEPTION; }
}
