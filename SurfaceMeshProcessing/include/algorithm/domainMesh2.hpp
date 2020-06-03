#ifndef DOMAINMESH_H
#define DOMAINMESH_H

#include <MeshDefinition.h>
typedef std::vector<Mesh::VertexHandle> VSet;
typedef std::vector<Mesh::EdgeHandle> ESet;
typedef std::vector<Mesh::HalfedgeHandle> HSet;
typedef std::vector<Mesh::FaceHandle> FSet;

template<class iter> class iterRange
{
private:
	iter begin_;
	iter end_;

public:
	iterRange(const iter &bg, const iter &ed) : begin_(bg), end_(ed) {}
	const iter begin() const { return begin_; }
	const iter end() const { return end_; }
	iter begin() { return begin_; }
	iter end() { return end_; }
};
typedef iterRange<VSet::const_iterator> VSetRange;

// connected domian for optmization
class Domain
{
private:
	Mesh &mesh;

	// the index of verts, order num in given part, -1 fixed point, -2 outer part 
	OpenMesh::VPropHandleT<int> vert_idx;
	OpenMesh::VPropHandleT<bool> boundary_vert;
	OpenMesh::FPropHandleT<bool> domain_face;
	VSet vertSet;
	ESet edgeSet;
	HSet halfSet;
	FSet faceSet;
	uint n_var;

public:
	Domain(Mesh &input) : mesh(input) 
	{ 
		mesh.add_property(vert_idx);
		mesh.add_property(domain_face);
		mesh.add_property(boundary_vert);
	}
	~Domain() 
	{ 
		mesh.remove_property(vert_idx); 
		mesh.remove_property(domain_face);
		mesh.remove_property(boundary_vert);
	}

	void build(const VSet &variables, const VSet &anchors)
	{
		clearDomain();
		addVertices(variables);
		addAnchors(anchors);
		buildDomain();
	}
	inline const VSet& vertices() const { return vertSet; }
	inline const ESet& edges() const { return edgeSet; }
	inline const HSet& halfedges() const { return halfSet; }
	inline const FSet& faces() const { return faceSet; }
	inline VSetRange variables() const { return VSetRange(vertSet.begin(), vertSet.begin() + n_var); }
	inline VSetRange anchors() const { return VSetRange(vertSet.begin() + n_var, vertSet.end()); }
	inline int n_variables() const { return n_var; }
	inline int index(const OpenMesh::VertexHandle &v) const { return mesh.property(vert_idx, v); }
	inline int &index(const OpenMesh::VertexHandle &v) { return mesh.property(vert_idx, v); }

	inline bool isVariable(const OpenMesh::VertexHandle &v) const { return index(v) >= 0; }
	inline bool isAnchor(const OpenMesh::VertexHandle &v) const { return index(v) == -1; }
	inline bool inDomain(const OpenMesh::VertexHandle &v) const { return index(v) > -2; }
	inline bool inDomain(const OpenMesh::HalfedgeHandle &h) const 
	{
		return mesh.face_handle(h).is_valid() && inDomain(mesh.face_handle(h));
	}
	inline bool inDomain(const OpenMesh::EdgeHandle &e) const 
	{
		return inDomain(mesh.from_vertex_handle(mesh.halfedge_handle(e, 0))) &&
			inDomain(mesh.to_vertex_handle(mesh.halfedge_handle(e, 0)));
	}
	inline bool inDomain(const OpenMesh::FaceHandle &f) const
	{
		return mesh.property(domain_face, f);
	}
	inline bool onBoundary(const OpenMesh::VertexHandle &v) const
	{
		return mesh.property(boundary_vert, v);
	}
	inline bool onBoundary(const OpenMesh::HalfedgeHandle &h) const
	{
		return !inDomain(h) && inDomain(mesh.opposite_halfedge_handle(h));
	}
	inline bool onBoundary(const OpenMesh::EdgeHandle &e) const
	{
		return onBoundary(mesh.from_vertex_handle(mesh.halfedge_handle(e, 0))) &&
			onBoundary(mesh.to_vertex_handle(mesh.halfedge_handle(e, 0)));
	}

private:
	void clearDomain() { for (auto v : mesh.vertices()) mesh.property(vert_idx, v) = -2; }
	void addVertex(const OpenMesh::VertexHandle &v) { mesh.property(vert_idx, v) = 0; }
	void addAnchor(const OpenMesh::VertexHandle &v) { mesh.property(vert_idx, v) = -1; }
	void addVertices(const VSet &vs) { for (auto v : vs) addVertex(v); }
	void addAnchors(const VSet &vs) { for (auto v : vs) addAnchor(v); }
	inline bool checkOnBoundary(const OpenMesh::VertexHandle &v) const
	{
		if (!inDomain(v)) return false;
		if (mesh.is_boundary(v)) return true;
		for (auto vv : mesh.vv_range(v))
		{
			if (!inDomain(vv)) return true;
		}
		return false;
	}
	void buildDomain()
	{
		VSet achs, vars;
		for (auto v : mesh.vertices())
		{
			switch (index(v))
			{
			case -2:
				break;
			case -1:
				achs.push_back(v);
				break;
			default:
				index(v) = vars.size();
				vars.push_back(v);
				break;
			}
		}
		
		n_var = vars.size();
		vertSet.clear();
		vertSet.insert(vertSet.end(), vars.begin(), vars.end());
		vertSet.insert(vertSet.end(), achs.begin(), achs.end());
		std::cout << "Domain vertices number : " << vertSet.size() << std::endl;

		for (auto v : mesh.vertices())
		{
			mesh.property(boundary_vert, v) = checkOnBoundary(v);
		}

		edgeSet.clear();
		for (auto e : mesh.edges())
		{
			auto he = mesh.halfedge_handle(e, 0);
			bool inDomain = index(mesh.from_vertex_handle(he)) != -2 && index(mesh.to_vertex_handle(he)) != -2;
			if (inDomain)
			{
				edgeSet.push_back(e);
			}
		}

		faceSet.clear();
		for (auto f : mesh.faces())
		{
			bool inDomain = true;
			for (auto fv : mesh.fv_range(f))
			{
				if (index(fv) == -2)
				{
					inDomain = false;
					break;
				}
			}
			mesh.property(domain_face, f) = false;
			if (inDomain)
			{
				faceSet.push_back(f);
				mesh.property(domain_face, f) = true;
			}
		}

		halfSet.clear();
		for (auto he : mesh.halfedges())
		{
			if (inDomain(he))
			{
				halfSet.push_back(he);
			}
		}
	}

};

#endif // !DOMAINMESH_H