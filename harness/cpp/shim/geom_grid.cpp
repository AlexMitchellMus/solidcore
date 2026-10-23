// SPDX-License-Identifier: Apache-2.0
// difftest shim for unit `geom.grid`. Wraps upstream SmGrid / SmVoxel / SmGridElement /
// SmHierarchyElement (source/SMLib/src/SmGrid.cpp, source/SMLib/inc/SmGrid.h).
//
// sc_geom_grid_{ints,floats}(prog, n, out, cap, out_len) run the same flat double program
// against one grid and return its integer or float output stream. The encoding is shared
// with crates/solidcore-ffi/src/units/geom_grid.rs, which validates it first (so nothing
// here can index outside the grid or dereference a fake pointer).
//   0  e6 nx ny nz        SmGrid(bbox, lSize)                     (first, exactly once)
//   1  type obj e6        AddElement(fake ptr obj, SM_GE_POLYBREP (1) | SM_GE_USERDATA (3), bbox)
//                         ints: status
//   2  kind surf tol e6   AddHierarchyElement(new SmHierarchyElement(tol) (kind 0) |
//                         new SmSurfaceHierarchyElement(fake surf, tol) (kind 1), bbox)
//                         ints: status  (a fake-surface element can be dumped, never expanded)
//   3  e6                 FindVoxelsInBox                          ints: status, n, n x (ix iy iz)
//   4  p3                 GetVoxelAddress                          ints: ix iy iz outside
//   5  ix iy iz           GetVoxelBBox                             floats: e6
//   6  p3 v3              Setup3DDDA                               ints: status, hits_nothing
//                         [+ dda dump when status 0 and the ray hits the grid]
//   7  t                  Step3DDDA(t)  ints: -3 when the walk is not ready (no hitting setup,
//                         or the last step failed or left the grid); else status, and on
//                         success exited, n, n x (type obj); floats: mins, maxs
//                         [+ dda dump when not exited]
//   8                     dump: per voxel (x-major) n_elems (-1 = null list) + (type obj)...,
//                         n_hier (-1 = null) + hierarchy indices...; per distinct element in
//                         that voxel order: (type obj) ints, centre radius e6 floats (a surface
//                         element then also its 51-value leaf record, below); per
//                         hierarchy element: type, back-ref index ints, centre radius e6
//                         floats; grid: voxel size, voxel tolerance, current mark (int)
//   9  h                  ExpandHierarchy of hierarchy element h (a base-class element or one
//                         made by op 12):  ints: status
//   10 s node e6          AddElement(node, SM_GE_SURFACE, bbox) for node `node` of
//                         smsurf_GetSurfaceCache(surface s)->GetTree()->GetAllTreeNodes():
//                         ints: status, or -4 when the node index is out of range or the node is
//                         not an SM_AD_BEZIER_SURFACE node (upstream casts it and reads garbage)
//   11 kind ...           declares the next surface (before its first use; at most 4):
//                         kind 0: p q cpts((p+1)(q+1) x (x y z w)), a clamped single-patch
//                         SmBSplineSurface of degree p x q; kind 1: o x y scale[2] dom[4], an
//                         SmPlane.  No output.
//   12 s tol e6           AddHierarchyElement(new SmSurfaceHierarchyElement(surface s, tol), bbox)
//                         ints: status
// Surface elements are reported as (SmSurfaceGridElement_TYPE, surface index), their leaf record
// (51 floats, as the ray tracer shim's AddSurfaceToGrid record) follows the sphere and box: uv
// domain (4), box (6), mid point, du, dv, normal (12), pseudo box basis (9) and intervals (6),
// polar box state, domain (4) and basis (9) (zeros when the polar box has no bound).
// dda dump: ints address, then per axis with next_t < SM_BIG_DOUBLE: step out; floats ray
// point, ray vector, ray t, current point, next_t, delta_t, the three delta vectors, then
// per axis with next_t < SM_BIG_DOUBLE: next plane point.
// Objects are fake pointers g_objs + id (stored, never dereferenced); elements are
// reported as (type code, object id).
#include "sc_shim.h"

#include <SmBSplineSurface.h>
#include <SmContext.h>
#include <SmExtent1d.h>
#include <SmExtent2d.h>
#include <SmExtent3d.h>
#include <SmGrid.h>
#include <SmNurbsSrf.h>
#include <SmPlane.h>
#include <SmPolarBox.h>
#include <SmPseudoBox.h>
#include <SmSurface.h>
#include <SmSurfaceCache.h>
#include <SmTArray.h>
#include <SmTree.h>
#include <SmVector2d.h>
#include <SmVector3d.h>
#include <nurbs.h>

#include <map>
#include <set>
#include <vector>

namespace
{
const int64_t kMaxObjects = 4096;
alignas(16) char g_objs[kMaxObjects];

void* obj_ptr(double id) { return reinterpret_cast<void*>(g_objs + (long long)id); }

int64_t obj_id(const void* p) { return p == nullptr ? -1 : (int64_t)(static_cast<const char*>(p) - g_objs); }

SmExtent3d ext(const double* d) { return SmExtent3d(d[0], d[1], d[2], d[3], d[4], d[5]); }

void put3(std::vector<double>& o, const SmVector3d& v)
{
    o.push_back(v.x);
    o.push_back(v.y);
    o.push_back(v.z);
}

void put_ext(std::vector<double>& o, const SmExtent3d& e)
{
    put3(o, e.GetMin());
    put3(o, e.GetMax());
}

int64_t st(SmStatus s) { return s == SM_SUCCESS ? 0 : (int64_t)s; }

const int64_t kMaxSurfaces = 4;

// The operand count of op `op`, whose operands start at `a` (`avail` doubles remain); -1 for an
// unknown or malformed op. Op 11 has a variable length.
int64_t op_size(int op, const double* a, int64_t avail)
{
    switch (op)
    {
    case 0: return 9;
    case 1: return 8;
    case 2: return 9;
    case 3: return 6;
    case 4: return 3;
    case 5: return 3;
    case 6: return 6;
    case 7: return 1;
    case 8: return 0;
    case 9: return 1;
    case 10: return 8;
    case 11:
    {
        if (avail < 1)
            return -1;
        if ((int)a[0] == 1)
            return 16;
        if ((int)a[0] != 0 || avail < 3)
            return -1;
        int64_t p = (int64_t)a[1], q = (int64_t)a[2];
        if (p < 1 || p > 3 || q < 1 || q > 3)
            return -1;
        return 3 + 4 * (p + 1) * (q + 1);
    }
    case 12: return 8;
    default: return -1;
    }
}

// `0 p q cpts` a clamped single-patch B-spline surface (knots 0..0 1..1), or `1 o x y scale dom`
// an SmPlane; the operands of op 11.
SmSurface* make_surface(const double* a, const SmContext& ctx)
{
    if ((int)a[0] == 0)
    {
        gw_DEGREE p = (gw_DEGREE)a[1];
        gw_DEGREE q = (gw_DEGREE)a[2];
        gw_INDEX n = p;
        gw_INDEX m = q;
        gw_SURFACE* sur = sm_AllocateNurbSurface(n, m, p, q, n + p + 1, m + q + 1);
        const double* c = a + 3;
        for (gw_INDEX i = 0; i <= n; i++)
            for (gw_INDEX j = 0; j <= m; j++)
            {
                sur->net->Pw[i][j].x = c[0];
                sur->net->Pw[i][j].y = c[1];
                sur->net->Pw[i][j].z = c[2];
                sur->net->Pw[i][j].w = c[3];
                c += 4;
            }
        for (gw_INDEX i = 0; i <= n + p + 1; i++)
            sur->knu->U[i] = i <= p ? 0.0 : 1.0;
        for (gw_INDEX i = 0; i <= m + q + 1; i++)
            sur->knv->U[i] = i <= q ? 0.0 : 1.0;
        return new (ctx) SmBSplineSurface(sur, FALSE, &ctx);
    }
    return new (ctx) SmPlane(SmPoint3d(a[1], a[2], a[3]), SmVector3d(a[4], a[5], a[6]), SmVector3d(a[7], a[8], a[9]),
                             SmVector2d(a[10], a[11]), SmExtent2d(a[12], a[13], a[14], a[15]), &ctx);
}

// Exposes the protected members of SmGrid read by the dumps.
class ShimGrid : public SmGrid
{
  public:
    ShimGrid(const SmExtent3d& b, ULONG lSize[3]) : SmGrid(b, lSize) {}

    long size(int i) const { return m_lSize[i]; }

    void dump_dda(std::vector<int64_t>& oi, std::vector<double>& of) const
    {
        for (int i = 0; i < 3; i++)
            oi.push_back(m_lAddress[i]);
        const double next[3] = {m_dNextT.x, m_dNextT.y, m_dNextT.z};
        const long step[3] = {m_lStepX, m_lStepY, m_lStepZ};
        const long out[3] = {m_lOutX, m_lOutY, m_lOutZ};
        const SmPoint3d* np[3] = {&m_sNextX, &m_sNextY, &m_sNextZ};
        for (int i = 0; i < 3; i++)
            if (next[i] < SM_BIG_DOUBLE)
            {
                oi.push_back(step[i]);
                oi.push_back(out[i]);
            }
        put3(of, m_sRayPoint);
        put3(of, m_sRayVector);
        of.push_back(m_dRayT);
        put3(of, m_sCurrPoss);
        put3(of, m_dNextT);
        put3(of, m_dDeltaT);
        put3(of, m_sDeltaXVec);
        put3(of, m_sDeltaYVec);
        put3(of, m_sDeltaZVec);
        for (int i = 0; i < 3; i++)
            if (next[i] < SM_BIG_DOUBLE)
                put3(of, *np[i]);
    }

    void dump_grid(std::vector<int64_t>& oi, std::vector<double>& of) const
    {
        put3(of, m_sVoxelSize);
        of.push_back(m_dVoxelTolerance);
        oi.push_back((int64_t)m_lCurrentMark);
    }
};

void put_elem(std::vector<int64_t>& oi, SmGridElement* e, const std::vector<SmSurface*>& surfs)
{
    oi.push_back((int64_t)e->GetType());
    if (e->GetType() == SmPolygonGridElement_TYPE)
        oi.push_back(obj_id(((SmPolygonGridElement*)e)->GetElement()));
    else if (e->GetType() == SmUserGridElement_TYPE)
        oi.push_back(obj_id(((SmUserGridElement*)e)->GetElement()));
    else if (e->GetType() == SmSurfaceGridElement_TYPE)
    {
        // the index of the element's surface (-2: none of ours)
        int64_t idx = -2;
        for (size_t i = 0; i < surfs.size(); i++)
            if (((SmSurfaceGridElement*)e)->GetElement() == surfs[i])
                idx = (int64_t)i;
        oi.push_back(idx);
    }
    else
        oi.push_back(-2);
}

// The 51-value leaf record of a surface grid element (layout in the header comment).
void put_leaf(std::vector<double>& o, const SmSurfaceGridElement& e)
{
    o.push_back(e.m_sUVDomain.GetMin().x);
    o.push_back(e.m_sUVDomain.GetMin().y);
    o.push_back(e.m_sUVDomain.GetMax().x);
    o.push_back(e.m_sUVDomain.GetMax().y);
    put3(o, e.m_sBBox.GetMin());
    put3(o, e.m_sBBox.GetMax());
    put3(o, e.m_sMidPoint);
    put3(o, e.m_sMidDU);
    put3(o, e.m_sMidDV);
    put3(o, e.m_sMidNormal);
    // the pseudo box: basis and intervals
    SmVector3d b[3];
    e.m_sPseudoBox.GetBasis(b[0], b[1], b[2]);
    for (int j = 0; j < 3; j++)
        put3(o, b[j]);
    SmExtent1d iv[3];
    e.m_sPseudoBox.GetIntervals(iv[0], iv[1], iv[2]);
    for (int j = 0; j < 3; j++)
    {
        o.push_back(iv[j].GetMin());
        o.push_back(iv[j].GetMax());
    }
    // the polar box: state, and domain and basis once it has a bound
    int ps = (int)e.m_sPolarBox.GetPolarBoxState();
    o.push_back((double)ps);
    bool bounded = ps >= SM_PS_SINGLE_POINT && ps <= SM_PS_REGION;
    SmExtent2d pd = e.m_sPolarBox.GetDomain();
    o.push_back(bounded ? pd.GetMin().x : 0.0);
    o.push_back(bounded ? pd.GetMin().y : 0.0);
    o.push_back(bounded ? pd.GetMax().x : 0.0);
    o.push_back(bounded ? pd.GetMax().y : 0.0);
    for (ULONG j = 0; j < 3; j++)
    {
        if (bounded)
            put3(o, *e.m_sPolarBox.GetBasis(j));
        else
            put3(o, SmVector3d(0.0, 0.0, 0.0));
    }
}

void put_sphere_box(std::vector<double>& of, const SmGridElement& e)
{
    put3(of, e.GetSphereCenter());
    of.push_back(e.GetSphereRadius());
    put_ext(of, e.GetBoundingBox());
    if (e.GetType() == SmSurfaceGridElement_TYPE)
        put_leaf(of, static_cast<const SmSurfaceGridElement&>(e));
}

int32_t run(const double* p, int64_t n, std::vector<int64_t>& oi, std::vector<double>& of)
{
    if (n < 10 || (int)p[0] != 0)
        return SC_ERR_ARG;
    SmContext sCtx;
    ULONG lSize[3] = {(ULONG)p[7], (ULONG)p[8], (ULONG)p[9]};
    ShimGrid* g = new (sCtx) ShimGrid(ext(p + 1), lSize);
    std::vector<SmHierarchyElement*> hier;
    std::map<const SmHierarchyElement*, int64_t> hidx;
    std::set<const SmHierarchyElement*> expandable; // base-class and real-surface elements
    std::vector<SmSurface*> surfs;
    bool ready = false;

    int64_t i = 10;
    int32_t rc = SC_OK;
    while (i < n)
    {
        int op = (int)p[i];
        int64_t sz = op_size(op, p + i + 1, n - i - 1);
        if (sz < 0 || op == 0 || i + 1 + sz > n)
        {
            rc = SC_ERR_ARG;
            break;
        }
        const double* a = p + i + 1;
        switch (op)
        {
        case 1:
        {
            SmGridElementType t = (int)a[0] == 1 ? SM_GE_POLYBREP : SM_GE_USERDATA;
            oi.push_back(st(g->AddElement(obj_ptr(a[1]), t, ext(a + 2))));
            break;
        }
        case 2:
        {
            SmHierarchyElement* h = (int)a[0] == 0
                                        ? new (sCtx) SmHierarchyElement(a[2])
                                        : new (sCtx) SmSurfaceHierarchyElement((SmSurface*)obj_ptr(a[1]), a[2]);
            hidx[h] = (int64_t)hier.size();
            hier.push_back(h);
            if ((int)a[0] == 0)
                expandable.insert(h);
            oi.push_back(st(g->AddHierarchyElement(h, ext(a + 3))));
            break;
        }
        case 3:
        {
            SmTArray<SmVoxel*> sVoxels;
            oi.push_back(st(g->FindVoxelsInBox(ext(a), sVoxels)));
            oi.push_back((int64_t)sVoxels.GetSize());
            // map each voxel pointer back to its address
            for (ULONG k = 0; k < sVoxels.GetSize(); k++)
            {
                bool found = false;
                for (long ix = 0; ix < g->size(0) && !found; ix++)
                    for (long iy = 0; iy < g->size(1) && !found; iy++)
                        for (long iz = 0; iz < g->size(2) && !found; iz++)
                        {
                            long ad[3] = {ix, iy, iz};
                            if (g->GetVoxel(ad) == sVoxels[k])
                            {
                                oi.push_back(ix);
                                oi.push_back(iy);
                                oi.push_back(iz);
                                found = true;
                            }
                        }
                if (!found)
                {
                    oi.push_back(-1);
                    oi.push_back(-1);
                    oi.push_back(-1);
                }
            }
            break;
        }
        case 4:
        {
            long ad[3] = {0, 0, 0};
            SmBoolean bOut = FALSE;
            SmStatus s = g->GetVoxelAddress(SmPoint3d(a[0], a[1], a[2]), ad, bOut);
            (void)s;
            oi.push_back(ad[0]);
            oi.push_back(ad[1]);
            oi.push_back(ad[2]);
            oi.push_back(bOut ? 1 : 0);
            break;
        }
        case 5:
        {
            long ad[3] = {(long)a[0], (long)a[1], (long)a[2]};
            put_ext(of, g->GetVoxelBBox(ad));
            break;
        }
        case 6:
        {
            SmBoolean bHits = FALSE;
            SmStatus s = g->Setup3DDDA(SmPoint3d(a[0], a[1], a[2]), SmVector3d(a[3], a[4], a[5]), bHits);
            oi.push_back(st(s));
            oi.push_back(bHits ? 1 : 0);
            ready = s == SM_SUCCESS && !bHits;
            if (ready)
                g->dump_dda(oi, of);
            break;
        }
        case 7:
        {
            if (!ready)
            {
                oi.push_back(-3);
                break;
            }
            SmBoolean bExited = FALSE;
            long cur[3] = {0, 0, 0};
            SmTArray<SmGridElement*> sElems;
            SmTArray<double> sMin, sMax;
            SmStatus s = g->Step3DDDA(a[0], bExited, cur, sElems, sMin, sMax);
            oi.push_back(st(s));
            if (s != SM_SUCCESS)
            {
                ready = false;
                break;
            }
            oi.push_back(bExited ? 1 : 0);
            oi.push_back((int64_t)sElems.GetSize());
            for (ULONG k = 0; k < sElems.GetSize(); k++)
                put_elem(oi, sElems[k], surfs);
            for (ULONG k = 0; k < sMin.GetSize(); k++)
                of.push_back(sMin[k]);
            for (ULONG k = 0; k < sMax.GetSize(); k++)
                of.push_back(sMax[k]);
            if (bExited)
                ready = false;
            else
                g->dump_dda(oi, of);
            break;
        }
        case 8:
        {
            std::vector<SmGridElement*> order;
            std::set<SmGridElement*> seen;
            for (long ix = 0; ix < g->size(0); ix++)
                for (long iy = 0; iy < g->size(1); iy++)
                    for (long iz = 0; iz < g->size(2); iz++)
                    {
                        long ad[3] = {ix, iy, iz};
                        SmVoxel* v = g->GetVoxel(ad);
                        SmTArray<SmGridElement*>* pe = v->GetGridElements();
                        if (pe == nullptr)
                            oi.push_back(-1);
                        else
                        {
                            oi.push_back((int64_t)pe->GetSize());
                            for (ULONG k = 0; k < pe->GetSize(); k++)
                            {
                                SmGridElement* e = (*pe)[k];
                                put_elem(oi, e, surfs);
                                if (seen.insert(e).second)
                                    order.push_back(e);
                            }
                        }
                        SmTArray<SmHierarchyElement*>* ph = v->GetHierarchyElements();
                        if (ph == nullptr)
                            oi.push_back(-1);
                        else
                        {
                            oi.push_back((int64_t)ph->GetSize());
                            for (ULONG k = 0; k < ph->GetSize(); k++)
                            {
                                auto it = hidx.find((*ph)[k]);
                                oi.push_back(it == hidx.end() ? -2 : it->second);
                            }
                        }
                    }
            for (SmGridElement* e : order)
            {
                put_elem(oi, e, surfs);
                put_sphere_box(of, *e);
            }
            for (SmHierarchyElement* h : hier)
            {
                SmHierarchyGridElement& ge = h->GetHierarchyGridElement();
                oi.push_back((int64_t)ge.GetType());
                auto it = hidx.find(ge.GetElement());
                oi.push_back(ge.GetElement() == nullptr ? -1 : (it == hidx.end() ? -2 : it->second));
                put_sphere_box(of, ge);
            }
            g->dump_grid(oi, of);
            break;
        }
        case 9:
        {
            // validated by the FFI: a fake-surface element would be dereferenced
            if ((size_t)a[0] >= hier.size() || expandable.count(hier[(size_t)a[0]]) == 0)
            {
                rc = SC_ERR_ARG;
                break;
            }
            // upstream's ExpandHierarchy reads the cache and its tree unchecked (a cache whose
            // tree cannot be built has none): not a program this shim can answer
            SmSurface* pSurf = hier[(size_t)a[0]]->GetSurface();
            if (pSurf != nullptr)
            {
                SmSurfaceCache* pSrfCache = smsurf_GetSurfaceCache(pSurf);
                if (pSrfCache == nullptr || pSrfCache->GetTree() == nullptr)
                {
                    rc = SC_ERR_ARG;
                    break;
                }
            }
            oi.push_back(st(hier[(size_t)a[0]]->ExpandHierarchy(g)));
            break;
        }
        case 10:
        {
            if ((size_t)a[0] >= surfs.size())
            {
                rc = SC_ERR_ARG;
                break;
            }
            // smsurf_GetSurfaceCache returns the surface's one cache
            SmSurfaceCache* pCache = smsurf_GetSurfaceCache(surfs[(size_t)a[0]]);
            if (pCache == nullptr || pCache->GetTree() == nullptr)
            {
                rc = SC_ERR_ARG;
                break;
            }
            SmTArray<SmTreeNode*> sNodes;
            pCache->GetTree()->GetAllTreeNodes(sNodes);
            double idx = a[1];
            if (!(idx >= 0.0 && idx < (double)sNodes.GetSize()) || sNodes[(ULONG)idx]->m_eAuxDataType != SM_AD_BEZIER_SURFACE)
            {
                oi.push_back(-4);
                break;
            }
            oi.push_back(st(g->AddElement(sNodes[(ULONG)idx], SM_GE_SURFACE, ext(a + 2))));
            break;
        }
        case 11:
        {
            if ((int64_t)surfs.size() >= kMaxSurfaces)
            {
                rc = SC_ERR_ARG;
                break;
            }
            surfs.push_back(make_surface(a, sCtx));
            break;
        }
        case 12:
        {
            if ((size_t)a[0] >= surfs.size())
            {
                rc = SC_ERR_ARG;
                break;
            }
            SmHierarchyElement* h = new (sCtx) SmSurfaceHierarchyElement(surfs[(size_t)a[0]], a[1]);
            hidx[h] = (int64_t)hier.size();
            hier.push_back(h);
            expandable.insert(h);
            oi.push_back(st(g->AddHierarchyElement(h, ext(a + 2))));
            break;
        }
        default:
            rc = SC_ERR_ARG;
            break;
        }
        if (rc != SC_OK)
            break;
        i += 1 + sz;
    }
    // the grid first (it owns the hierarchy elements), then the surfaces it refers to
    delete g;
    for (SmSurface* s : surfs)
        delete s;
    return rc;
}
} // namespace

SC_EXPORT int32_t sc_geom_grid_ints(const double* prog, int64_t n, int64_t* out, int64_t cap, int64_t* out_len)
{
    SC_TRY
    std::vector<int64_t> oi;
    std::vector<double> of;
    int32_t s = run(prog, n, oi, of);
    if (s != SC_OK)
        return s;
    return sc_copy_out(oi.data(), (int64_t)oi.size(), out, cap, out_len);
    SC_CATCH
}

SC_EXPORT int32_t sc_geom_grid_floats(const double* prog, int64_t n, double* out, int64_t cap, int64_t* out_len)
{
    SC_TRY
    std::vector<int64_t> oi;
    std::vector<double> of;
    int32_t s = run(prog, n, oi, of);
    if (s != SC_OK)
        return s;
    return sc_copy_out(of.data(), (int64_t)of.size(), out, cap, out_len);
    SC_CATCH
}
