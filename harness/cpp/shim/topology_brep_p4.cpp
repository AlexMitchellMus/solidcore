// SPDX-License-Identifier: Apache-2.0
// difftest shim for unit `topology.brep_p4`: SmBrep::DeleteEdge, DeleteVertex,
// ShrinkGeometry, StitchFaces and ~SmBrep (source/SMLib/src/SmBrep.cpp) on real B-reps
// built by the shared fixture (sc_topology_fixture.h).
//
// One entry point, sc_topology_brep_p4_eval(op, path, in, n_in, editing, target, keep,
// tol, flag, out, cap, out_len). `in` is a fixture B-rep block (validated by the Rust
// wrapper through the topology.topology_traverser export, whose `walk(brep, g, true)`
// order names `target` and `keep`). Ops:
//   0 DeleteEdge(edges[target], keep >= 0 ? faces[keep]->GetSurface() : NULL)
//   1 DeleteVertex(vertices[target])
//   2 ShrinkGeometry()
//   3 StitchFaces(tol, .., keep >= 0 ? {faces[keep]} : NULL, flag)
//   4 delete brep (~SmBrep)
// `editing` sets m_bEditingEnabled first. Output: `[status]`, then on success the ten
// record counts of the B-rep walked again (op 3 adds its four outputs first; op 4
// outputs the context's DoingBoolean flag after the destructor). `[-1]` when the
// fixture does not build.
#include "sc_topology_fixture.h"

SC_EXPORT int32_t sc_topology_brep_p4_eval(int32_t op, const char* path, const double* in, int64_t n_in,
                                           int32_t editing, int64_t target, int64_t keep, double tol,
                                           int32_t flag, double* out, int64_t cap, int64_t* out_len)
{
    if (!path || !in || n_in < 1 || !out_len || op < 0 || op > 4)
        return SC_ERR_ARG;
    SC_TRY
    Rd rd{in, 0};
    SmStatus status;
    Built built;
    built.brep = build(rd, path, status);
    std::vector<double> o;
    if (!built.brep)
    {
        o.push_back(-1.0);
        return sc_copy_out(o.data(), (int64_t)o.size(), out, cap, out_len);
    }
    SmBrep* brep = built.brep;
    Graph g;
    walk(brep, g, true);
    brep->m_bEditingEnabled = editing ? TRUE : FALSE;
    SmFace* keep_face = (keep >= 0 && keep < (int64_t)g.faces.size()) ? g.faces[(size_t)keep] : nullptr;
    if (keep >= 0 && !keep_face)
        return SC_ERR_ARG;
    switch (op)
    {
    case 0:
        if (target < 0 || target >= (int64_t)g.edges.size())
            return SC_ERR_ARG;
        status = brep->DeleteEdge(g.edges[(size_t)target], keep_face ? keep_face->GetSurface() : nullptr);
        break;
    case 1:
        if (target < 0 || target >= (int64_t)g.vertices.size())
            return SC_ERR_ARG;
        status = brep->DeleteVertex(g.vertices[(size_t)target]);
        break;
    case 2:
        status = brep->ShrinkGeometry();
        break;
    case 3:
    {
        ULONG stitched = 0, lamina = 0;
        double vgap = 0.0, egap = 0.0;
        SmTArray<SmFace*> faces;
        if (keep_face)
            faces.Add(keep_face);
        status = brep->StitchFaces(tol, stitched, lamina, vgap, egap, keep_face ? &faces : nullptr,
                                   flag ? TRUE : FALSE);
        o.push_back(st(status));
        if (status == SM_SUCCESS)
        {
            o.push_back((double)stitched);
            o.push_back((double)lamina);
            o.push_back(vgap);
            o.push_back(egap);
        }
        break;
    }
    default:
    {
        SmContext* ctx = SM_CONST_CAST(SmContext*, &the_context());
        SmBoolean prev = ctx->GetDoingBoolean();
        delete brep;
        built.brep = nullptr;
        o.push_back(0.0);
        o.push_back(flg(ctx->GetDoingBoolean()));
        ctx->SetDoingBoolean(prev);
        return sc_copy_out(o.data(), (int64_t)o.size(), out, cap, out_len);
    }
    }
    if (op != 3)
        o.push_back(st(status));
    if (status == SM_SUCCESS)
    {
        Graph after;
        walk(brep, after, true);
        for (auto n : {after.regions.size(), after.shells.size(), after.faceuses.size(), after.faces.size(),
                       after.loopuses.size(), after.loops.size(), after.edgeuses.size(), after.edges.size(),
                       after.vertexuses.size(), after.vertices.size()})
            o.push_back((double)n);
    }
    return sc_copy_out(o.data(), (int64_t)o.size(), out, cap, out_len);
    SC_CATCH
}
