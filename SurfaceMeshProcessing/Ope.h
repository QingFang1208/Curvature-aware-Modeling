#pragma once
#include <MeshDefinition.h>

typedef OpenMesh::VertexHandle VH;

class Ope
{
public:
	Ope(Mesh& input, OpenMesh::EPropHandleT<bool>& is_e, OpenMesh::VPropHandleT<bool>& is_v):
		mesh_(input), is_remesh_e_(is_e), is_remesh_v_(is_v)
	{
		compute_average_edge_length();
	}

	~Ope();

	void Split();

	void Collapse();

	void Equalize_valences();

	void Tangential_relaxation();

private:
	void compute_average_edge_length();

	inline bool is_bound_flip(VH a, VH b, VH c, VH d) {
		if (!mesh_.property(is_remesh_v_, a)
			|| !mesh_.property(is_remesh_v_, b)
			|| !mesh_.property(is_remesh_v_, c)
			|| !mesh_.property(is_remesh_v_, d))
			return true;
		else
		{
			return false;
		}
	}

private:
	double high_, low_;

	Mesh& mesh_;

	OpenMesh::EPropHandleT<bool>& is_remesh_e_;
	OpenMesh::VPropHandleT<bool>& is_remesh_v_;
};

