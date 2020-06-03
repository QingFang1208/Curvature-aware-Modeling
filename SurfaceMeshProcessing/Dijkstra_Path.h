#pragma once
#include <queue>
#include <vector>
#include <MeshDefinition.h>

typedef OpenMesh::VertexHandle VH;
typedef OpenMesh::HalfedgeHandle HEH;

struct Node {
	VH v_h;

	int version;

	double dist_;

	bool operator<(const Node& n)const {
		return dist_ > n.dist_;
	}
};

class Dijkstra_Path
{
public:
	Dijkstra_Path(Mesh& mesh) :mesh_(mesh) {}

	~Dijkstra_Path();

	void ComputePath(VH& start, VH& end);

	std::vector<HEH> return_path();

private:
	Mesh& mesh_;

	OpenMesh::VPropHandleT<OpenMesh::HalfedgeHandle> forward_he;
	OpenMesh::VPropHandleT<int> version;
	OpenMesh::VPropHandleT<double> dist;

	std::vector<HEH> Path;
};

