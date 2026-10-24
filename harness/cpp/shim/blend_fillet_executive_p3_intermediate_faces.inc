// Actual native intermediate body allocation/corner/solver face component.
#include <SmFaceuse.h>
#include <SmRegion.h>
namespace fillet_real_solid {
namespace {
struct FaceExecutive : Executive {
    using Executive::Executive;
    SmStatus build_faces() {
        m_pFilletBrep=new(m_crContext) SmBrep();NER(m_pFilletBrep);
        m_pFilletBrep->SetTolerance(m_pTargetBrep->GetTolerance());m_pFilletBrep->m_bEditingEnabled=TRUE;
        for(ULONG i=0;i<m_vFilletCorners.GetSize();++i) {
            auto* corner=m_vFilletCorners[i];
            if(corner->CalcCornerGeom()!=SM_SUCCESS) {
                TCHAR name[SM_TBLOCK_SIZE];corner->GetName(name,SM_TBLOCK_SIZE);
                NoteFilletError(SM_FILERR_VERTEX_PROBLEM,corner->GetCornerType(),SM_FIL_FAILURE,name,corner->GetFilletedVertex(),nullptr,nullptr,nullptr,nullptr,_T("Fillet Corner Error: Corner Geometry Computation Failure"));return SM_ERR;
            }
        }
        for(ULONG i=0;i<m_vFilletSolvers.GetSize();++i) {
            auto* solver=m_vFilletSolvers[i];double tolerance=solver->GetThisApproxTol3d();
            if(m_pFilletBrep->GetTolerance()<tolerance) {
#ifdef SM_USE_NEWTOL
                SM_NEWTOL_LINE m_pFilletBrep->SetZoneTol3d(tolerance);
#else
                SM_OLDTOL_LINE if(m_pFilletBrep->GetTolerance()<tolerance) {m_pFilletBrep->SetTolerance(tolerance);}
#endif
            }
            SmTArray<SmFilletGeom*> geometries;solver->GetFilletGeoms(geometries);
            for(ULONG g=0;g<geometries.GetSize();++g)if(geometries[g]->MakeFaceBrep(m_pFilletBrep)!=SM_SUCCESS) {
                GetFilletErrorInfo()->NoteError(solver->GetEdgeuse(0)->GetEdge(),SM_FIL_FAILURE,SM_FILERR_EDGE_PROBLEM,_T("Fillet Surface Computation Failure: Fillet Face Creation Error"));return SM_ERR;
            }
        }
        return SM_SUCCESS;
    }
};
}
std::pair<std::vector<double>,std::vector<double>> intermediate_faces(SmBrep* body) {
    if(!body)return {{0},{}};
    SmTArray<SmVertex*> vertices;body->GetVertices(vertices);SmTArray<SmEdge*> edges;body->GetEdges(edges);
    SmTArray<SmFace*> faces;body->GetFaces(faces);SmTArray<SmRegion*> regions;body->GetRegions(regions);SmTArray<SmShell*> shells;body->GetShells(shells);
    SmTArray<SmEdgeuse*> uses;for(ULONG i=0;i<edges.GetSize();++i){SmTArray<SmEdgeuse*> own;edges[i]->GetEdgeuses(own);uses.Append(own);}
    std::vector<double> state{1,double(vertices.GetSize()),double(edges.GetSize()),double(faces.GetSize()),double(regions.GetSize()),double(shells.GetSize()),double(body->m_bEditingEnabled)};
    topology(state,body);state.push_back(double(uses.GetSize()));
    for(ULONG i=0;i<uses.GetSize();++i){auto* use=uses[i];auto* vu=use->GetVertexuse();
        state.insert(state.end(),{label(edges,use->GetEdge()),double(use->GetOrientation()),label(vertices,vu?vu->GetVertex():nullptr),label(uses,use->GetMate()),label(uses,use->GetRadial()),label(uses,use->GetCWEdgeuse()),label(uses,use->GetCCWEdgeuse()),label(faces,use->GetFace())});}
    for(ULONG i=0;i<faces.GetSize();++i){SmFaceuse *a=nullptr,*b=nullptr;faces[i]->GetFaceuses(a,b);
        for(auto* use:{a,b}){state.push_back(double(use!=nullptr));if(use)state.insert(state.end(),{double(use->GetOrientation()),label(shells,use->GetShell()),label(regions,use->GetRegion())});}}
    for(ULONG i=0;i<shells.GetSize();++i)state.push_back(label(regions,shells[i]->GetRegion()));
    for(ULONG i=0;i<regions.GetSize();++i)state.push_back(double(regions[i]->IsVoid()));
    std::vector<double> geometry{body->GetTolerance()};primitive_record::append(geometry,body);
    return {std::move(state),std::move(geometry)};
}
}
