// Actual common-mate marking and all-FG original-face selection; no deletion.
namespace fillet_real_solid {
inline void selection_context(std::vector<double>& out,const SmContext* context) {
    for(auto type:{SM_MT_MARK,SM_MT_MARK2,SM_MT_MARK3,SM_MT_MARKIO,SM_MT_MARKASSERT})
        out.insert(out.end(),{double(type),double(context->GetCurrentMark(type)),double(context->IsMarkLocked(type))});
}
template<class T> void selection_mates(std::vector<double>& out,const SmTArray<T*>& originals,const SmTArray<T*>& others,SmTopologyIntersector& ti) {
    std::vector<double> mates{double(originals.GetSize())};
    for(ULONG i=0;i<originals.GetSize();++i){auto* mate=ti.GetOtherMate(originals[i]);
        const auto other=label(others,static_cast<T*>(mate));
        if(mate&&other<0)throw std::runtime_error("selection mate outside actual managed traversal");
        mates.insert(mates.end(),{double(i),double(mate!=nullptr),other});
    }section(out,mates);
}
template<class T> void selection_marks(std::vector<double>& out,const SmTArray<T*>& entities) {
    std::vector<double> marks{double(entities.GetSize())};
    for(ULONG i=0;i<entities.GetSize();++i)
        marks.insert(marks.end(),{double(entities[i]->GetMark(SM_MT_MARK)),double(entities[i]->IsMarked(SM_MT_MARK))});
    section(out,marks);
}
inline std::vector<double> original_selection(FaceExecutive& exec,SmBrep* target,SmStatus& result) {
    auto* context=const_cast<SmContext*>(exec.GetContext());auto* managed=exec.GetFilletBrep();auto& ti=exec.GetTopologyIntersector();
    std::vector<double> out{double(exec.GetDoGlobalMerge())};selection_context(out,context);
    SmTArray<SmVertex*> target_vertices,managed_vertices;target->GetVertices(target_vertices);managed->GetVertices(managed_vertices);
    SmTArray<SmEdge*> target_edges,managed_edges;target->GetEdges(target_edges);managed->GetEdges(managed_edges);
    selection_mates(out,target_vertices,managed_vertices,ti);selection_mates(out,target_edges,managed_edges,ti);
    SmTArray<SmFace*> kept,deleted;
    {
        const auto editing=target->m_bEditingEnabled;
        struct RestoreEdit{SmBrep* target;SmBoolean editing;~RestoreEdit(){target->m_bEditingEnabled=editing;}} restore{target,editing};
        target->m_bEditingEnabled=TRUE;
        SmTArray<SmNewMarkAndLock*> locks;SmNewMarkAndLock fillet_lock(context),original_lock;
        const auto fillet_type=fillet_lock.GetMarkType();auto original_type=fillet_type;locks.Add(&fillet_lock);
        if(context!=target->GetContext()){original_lock.SetContext(const_cast<SmContext*>(target->GetContext()));locks.Add(&original_lock);original_type=original_lock.GetMarkType();}
        for(ULONG i=0;i<target_vertices.GetSize();++i)if(auto* mate=static_cast<SmVertex*>(ti.GetOtherMate(target_vertices[i]))){mate->Mark(fillet_type);target_vertices[i]->Mark(original_type);}
        for(ULONG i=0;i<target_edges.GetSize();++i)if(auto* mate=static_cast<SmEdge*>(ti.GetOtherMate(target_edges[i]))){target_edges[i]->Mark(original_type);mate->Mark(fillet_type);}
        for(ULONG s=0;result==SM_SUCCESS&&s<exec.GetFilletSolvers().GetSize();++s){SmTArray<SmFilletGeom*> geometries;exec.GetFilletSolvers()[s]->GetFilletGeoms(geometries);
            for(ULONG g=0;result==SM_SUCCESS&&g<geometries.GetSize();++g)result=geometries[g]->TrimOriginalSurfaces(kept,deleted,locks);}
    }
    out.push_back(status(result));selection_context(out,context);
    out.insert(out.end(),{double(context->GetDoingBoolean()),double(target->m_bEditingEnabled),double(ti.m_bDoingFillets)});
    SmTArray<SmFace*> target_faces,managed_faces;target->GetFaces(target_faces);managed->GetFaces(managed_faces);
    for(int body=0;body<2;++body){const auto& values=body==0?kept:deleted;const auto& owning=body==0?managed_faces:target_faces;std::vector<double> labels;
        for(ULONG i=0;i<values.GetSize();++i){auto value=label(owning,values[i]);if(value<0)throw std::runtime_error("selection output outside actual body");labels.push_back(value);}section(out,labels);}
    // This ordinary fixture's first available actual lock is SM_MT_MARK.
    // Raw values remain observed, alongside membership and all context lock states.
    selection_marks(out,target_vertices);selection_marks(out,target_edges);selection_marks(out,target_faces);
    selection_marks(out,managed_vertices);selection_marks(out,managed_edges);selection_marks(out,managed_faces);
    return out;
}
}
