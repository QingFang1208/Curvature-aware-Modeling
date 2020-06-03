#pragma once
#include <unordered_map>
#include <algorithm/segMesh.hpp>

typedef Mesh::VertexHandle VH;
typedef Mesh::FaceHandle FH;
typedef Mesh::EdgeHandle EH;
typedef Mesh::HalfedgeHandle HEH;
typedef std::vector<VH> VSet;
typedef std::vector<FH> FSet;
typedef std::vector<HEH> HESet;


class Mesh_Ope
{
public:
	Mesh_Ope(Mesh& ori, OpenMesh::FPropHandleT<int>& seg):mesh_(ori), seg_(seg){}

	~Mesh_Ope();

	Mesh Cut_Seg(int i, std::unordered_map<VH, VH>& old_to_new, std::unordered_map<VH, VH>& new_to_old);
	
	FSet Close_draw_bound(HESet& bound, int& seg_num);

private:

private:
	Mesh& mesh_;
	OpenMesh::FPropHandleT<int>& seg_;


};

