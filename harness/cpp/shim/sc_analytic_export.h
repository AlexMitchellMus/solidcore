// SPDX-License-Identifier: Apache-2.0
// Exact export of analytic face surfaces for the shared B-rep graph exports
// (sc_edgeuse_graph.h, sc_topology_fixture.h): the face surface kinds the Rust import can
// rebuild as their own `SurfaceKind` travel as their upstream database record (binary,
// SM_CURRENT_DATABASE_VERSION), which the importer reads with the ported `ReadFromDB`.
// Only `SmSurfOfRevolution` itself qualifies: its subclasses (`SmCone`, `SmCylinder`,
// `SmSphere`, `SmTorus`) and `SmSurfOfExtrusion` have no `SurfaceKind` in the port (the
// dispatch gap in `Surface::read_from_db`) and keep their B-spline part.
//
// Face surface block (layout in crates/solidcore-ffi/src/units/topology_edgeuse_p1.rs):
//   `2 type version surf_oob gen_oob nbytes byte[nbytes]` then the header's regular
//   surface block (the B-spline part), which the importer reads as the fallback.
// The record is only exported when it is lossless for the C++ object: the surface, its
// generator and its polar curve hold their NURBs (WriteToDB would otherwise create them,
// changing the object the C++ side runs), a C++ ReadFromDB of the bytes succeeds and
// writes back the same bytes, the polar curve survives the round trip, and the polar
// curve's out-of-bounds flag is the reader's default. `m_bOutOfBoundsEnabled` of the
// surface and the generator is not in the record; it is exported beside it.
#ifndef SC_ANALYTIC_EXPORT_H
#define SC_ANALYTIC_EXPORT_H

#include "StdAfx.h"
#include <SmBSplineCurve.h>
#include <SmBSplineSurface.h>
#include <SmConfig.h>
#include <SmContext.h>
#include <SmDatabaseIO.h>
#include <SmEllipse.h>
#include <SmSurfOfRevolution.h>
#include <SmSurface.h>

#include <sstream>
#include <string>
#include <vector>

namespace sc_analytic_export
{
struct SorPeek : SmSurfOfRevolution
{
    static SmPolarConversion SmSurfOfRevolution::*polar() { return &SorPeek::m_vPolarConverter; }
};

// The binary record of `s` (WriteToDB), or false.
inline bool record_of(const SmSurfOfRevolution* s, std::string& bytes)
{
    std::ostringstream oss;
    SmDatabaseIOStream db(SM_BINARY, &oss);
    if (s->WriteToDB(db, SM_CURRENT_DATABASE_VERSION) != SM_SUCCESS)
        return false;
    bytes = oss.str();
    return true;
}

inline bool oob_of(const SmBSplineCurve* c) { return c != nullptr && c->IsOutOfBoundsEnabled() != FALSE; }

// Appends the exact block head (`2 type version surf_oob gen_oob nbytes bytes`) of `s`
// when it can be exported exactly; returns whether it did (the caller then appends the
// regular block either way).
inline bool put_exact_head(std::vector<double>& o, const SmSurface* s, const SmContext& ctx)
{
    if (s == nullptr || s->GetType() != SmSurfOfRevolution_TYPE)
        return false;
    const SmSurfOfRevolution* r = (const SmSurfOfRevolution*)s;
    const SmBSplineCurve* gen = r->GetGenCurve();
    const SmPolarConversion& polar = r->*SorPeek::polar();
    const SmEllipse* pc = polar.GetCurve();
    // WriteToDB creates a missing NURB (GetOrCreateGwNurbPointer): never mutate the object
    if (r->GetGwNurbPointer() == nullptr || gen == nullptr || gen->GetGwNurbPointer() == nullptr ||
        (pc != nullptr && pc->GetGwNurbPointer() == nullptr))
        return false;
    std::string bytes;
    if (!record_of(r, bytes))
        return false;
    // the C++ reader must restore the same object: same record written back, the polar
    // curve kept (a borrowed one is not written) with the reader's out-of-bounds default
    bool ok = false;
    {
        std::istringstream iss(bytes);
        SmDatabaseIOStream db(SM_BINARY, &iss);
        SmSurface* back = nullptr;
        if (SmSurfOfRevolution::ReadFromDB(SmSurfOfRevolution_TYPE, db, ctx, back,
                                           SM_CURRENT_DATABASE_VERSION) == SM_SUCCESS &&
            back != nullptr)
        {
            const SmSurfOfRevolution* rb = (const SmSurfOfRevolution*)back;
            std::string again;
            const SmEllipse* pcb = (rb->*SorPeek::polar()).GetCurve();
            ok = record_of(rb, again) && again == bytes && (pc == nullptr) == (pcb == nullptr) &&
                 (pc == nullptr || oob_of(pc) == oob_of(pcb));
        }
        delete back;
    }
    if (!ok)
        return false;
    o.push_back(2.0);
    o.push_back((double)SmSurfOfRevolution_TYPE);
    o.push_back((double)SM_CURRENT_DATABASE_VERSION);
    o.push_back(r->IsOutOfBoundsEnabled() != FALSE ? 1.0 : 0.0);
    o.push_back(oob_of(gen) ? 1.0 : 0.0);
    o.push_back((double)bytes.size());
    for (char ch : bytes)
        o.push_back((double)(unsigned char)ch);
    return true;
}
} // namespace sc_analytic_export

#endif // SC_ANALYTIC_EXPORT_H
