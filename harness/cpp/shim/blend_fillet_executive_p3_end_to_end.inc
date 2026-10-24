// SPDX-License-Identifier: Apache-2.0
// Real public fillet API, preserving the selected pre-operation edge and the
// freshly traversed post-operation solid. No executive fixture owns fake objects.
#include "support/primitive_record.h"
#include <SmApiQueries.h>
#include <SmLoop.h>
#include <cmath>
#include <limits>

namespace fillet_real_solid {
inline double status(SmStatus value) { return value == SM_SUCCESS ? 0.0 : double(value); }
inline void section(std::vector<double>& out, const std::vector<double>& data) {
    out.push_back(double(data.size())); out.insert(out.end(), data.begin(), data.end());
}
template<class T> double label(const SmTArray<T*>& list, T* value) {
    if (!value) return -1.0;
    for (ULONG i=0; i<list.GetSize(); ++i) if (list[i] == value) return double(i);
    return -2.0;
}
inline void topology(std::vector<double>& out, SmBrep* b) {
    SmTArray<SmVertex*> vertices; b->GetVertices(vertices);
    SmTArray<SmEdge*> edges; b->GetEdges(edges);
    SmTArray<SmFace*> faces; b->GetFaces(faces);
    out.push_back(double(edges.GetSize()));
    for (ULONG i=0; i<edges.GetSize(); ++i) {
        auto* edge=edges[i];
        out.push_back(edge->GetCurve() ? double(edge->GetCurve()->GetType()) : 0.0);
        SmTArray<SmVertex*> ends; edge->GetVertices(ends);
        out.push_back(double(ends.GetSize()));
        for (ULONG j=0; j<ends.GetSize(); ++j) out.push_back(label(vertices,ends[j]));
        SmTArray<SmEdgeuse*> uses; edge->GetEdgeuses(uses);
        out.push_back(double(uses.GetSize()));
        for (ULONG j=0; j<uses.GetSize(); ++j)
            out.insert(out.end(), {label(faces,uses[j]->GetFace()),double(uses[j]->GetOrientation())});
    }
    out.push_back(double(faces.GetSize()));
    for (ULONG i=0; i<faces.GetSize(); ++i) {
        auto* face=faces[i];
        out.push_back(face->GetSurface() ? double(face->GetSurface()->GetType()) : 0.0);
        SmTArray<SmLoop*> loops; face->GetLoops(loops);
        out.push_back(double(loops.GetSize()));
        for (ULONG j=0; j<loops.GetSize(); ++j) {
            SmTArray<SmEdgeuse*> uses; loops[j]->GetEdgeuses(uses);
            out.push_back(double(uses.GetSize()));
            for (ULONG q=0; q<uses.GetSize(); ++q)
                out.insert(out.end(), {label(edges,uses[q]->GetEdge()),double(uses[q]->GetOrientation())});
        }
    }
}
inline void counts(std::vector<double>& out, SmBrep* b) {
    SmTArray<SmVertex*> vertices; b->GetVertices(vertices);
    SmTArray<SmEdge*> edges; b->GetEdges(edges);
    SmTArray<SmFace*> faces; b->GetFaces(faces);
    out.insert(out.end(), {double(vertices.GetSize()), double(edges.GetSize()), double(faces.GetSize())});
}
}

// Wire: dimensions XYZ, radius, xsect, continuity, thumbweight, edge count,
// then owning-list edge indices. This initial fixture constructor is a box.
SC_EXPORT int32_t sc_blend_fillet_executive_p3_real_solid(
    const double* input, int64_t n, double* out, int64_t cap, int64_t* out_len) {
    if (!input || n < 9) return SC_ERR_ARG;
    for (int64_t i=0; i<n; ++i) if (!std::isfinite(input[i])) return SC_ERR_ARG;
    for (int i=0; i<3; ++i) if (input[i] <= SM_EFF_ZERO) return SC_ERR_ARG;
    // The exclusive power-of-two bound is exact even when ULONG_MAX or
    // long double would round to 2^64. Check before the integer conversion.
    const double ulong_limit=std::ldexp(1.0,std::numeric_limits<ULONG>::digits);
    for (int i : {4,5}) if (input[i] < 0 || std::floor(input[i]) != input[i]
        || input[i] >= ulong_limit) return SC_ERR_ARG;
    if (input[7] < 1 || std::floor(input[7]) != input[7]
        || input[7] != double(n-8)) return SC_ERR_ARG;
    SC_TRY
    SmBrep* raw=nullptr;
    const auto setup=SmApiCreateBox(SmVector3d(0,0,0),input[0],input[1],input[2],raw);
    if (setup != SM_SUCCESS || !raw) return SC_ERR_ARG;
    std::unique_ptr<SmBrep> body(raw);
    SmTArray<SmEdge*> edges, selected; body->GetEdges(edges);
    std::vector<double> before;
    for (int64_t i=8; i<n; ++i) {
        if (input[i] < 0 || input[i] >= edges.GetSize() || std::floor(input[i]) != input[i]) return SC_ERR_ARG;
        auto* edge=edges[ULONG(input[i])]; selected.Add(edge);
        const auto interval=edge->GetInterval();
        before.insert(before.end(), {input[i],interval.GetMin(),interval.GetMax(),double(edge->GetCurve()!=nullptr)});
        if (auto* curve=edge->GetCurve()) {
            before.push_back(double(curve->GetType()));
            for (double fraction : {0.17,0.37,0.83}) {
                SmPoint3d p; const auto s=curve->EvaluatePoint(interval.Evaluate(fraction),p);
                primitive_record::sample(before,s,p);
            }
        }
    }
    const auto result=SmApiFilletEdges(body.get(),selected,input[3],ULONG(input[4]),ULONG(input[5]),input[6]);
    // Selected pointers may no longer be live; all result traversals are fresh.
    std::vector<double> values; fillet_real_solid::section(values,before);
    values.push_back(fillet_real_solid::status(result));
    fillet_real_solid::counts(values,body.get());
    std::vector<double> topology; fillet_real_solid::topology(topology,body.get());
    fillet_real_solid::section(values,topology);
    std::vector<double> geometry; primitive_record::append(geometry,body.get());
    fillet_real_solid::section(values,geometry);
    double area,volume,mass; SmPoint3d centroid; SmVector3d inertia,products;
    const auto mass_status=SmApiBrepComputeMassProperties(body.get(),1e-4,1.0,SmPoint3d(0,0,0),area,volume,mass,centroid,inertia,products);
    values.push_back(fillet_real_solid::status(mass_status));
    if (mass_status == SM_SUCCESS) {
        values.insert(values.end(), {area,volume,mass});
        primitive_record::point(values,centroid);
        primitive_record::point(values,inertia);
        primitive_record::point(values,products);
    }
    return sc_copy_out(values.data(),int64_t(values.size()),out,cap,out_len);
    SC_CATCH
}
