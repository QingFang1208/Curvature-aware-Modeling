#ifndef SEGMESH_H
#define SEGMESH_H

#include <MeshDefinition.h>
#include <Eigen/Eigen>
typedef std::vector<Mesh::VertexHandle> VSet;
typedef std::vector<Mesh::EdgeHandle> ESet;
typedef std::vector<Mesh::FaceHandle> FSet;

class segMesh
{
private:
	Mesh &mesh;
	const OpenMesh::FPropHandleT<int> &face_label;
	std::vector<bool> partEdit;

	VSet vert_inners; // verts in editing parts
	std::vector<VSet> vert_contours; // close contours array, counterclockwise
	std::vector<VSet *> vert_outers; // verts in non-editing parts
	Eigen::MatrixXi part_connect; // connection between contours and outers
	std::vector<double> contours_len; // length of contours
	std::vector<bool> contours_select; // true to set gauss curvature zero
	OpenMesh::VPropHandleT<int> dist; // layers shift distances

public:
	segMesh(Mesh &input, const OpenMesh::FPropHandleT<int> &seg, int segnum = -1) : mesh(input), face_label(seg)
	{
		if (segnum < 0)
		{
			for (auto f : mesh.faces())
			{
				if (segnum < mesh.property(seg, f))
					segnum = mesh.property(seg, f);
			}
			segnum = segnum + 1;
		}

		partEdit.resize(segnum);
		mesh.add_property(dist);
	}
	
	~segMesh()
	{
		clearParts();
		mesh.remove_property(dist);
	}

	void setParts(const std::vector<int> &labels) 
	{
		int segnum = -1;
		for (auto f : mesh.faces())
		{
			if (segnum < mesh.property(face_label, f))
				segnum = mesh.property(face_label, f);
		}
		segnum = segnum + 1;
		partEdit.resize(segnum);

		clearParts();
		for (auto label : labels)
		{
			if (label < 0 || label >= partEdit.size()) continue;
			partEdit[label] = true;
		}

		VSet remain;
		processParts(remain);
		processConnect(remain);
	}
	
	inline const VSet& inner() const { return vert_inners; }
	inline const VSet& outer(int i) const { return *(vert_outers[i]); }
	inline int n_outers() const { return vert_outers.size(); }
	inline const VSet& contour(int i) const { return vert_contours[i]; }
	inline int n_contours() const { return vert_contours.size(); }
	inline const std::vector<VSet> &vContours() const { return vert_contours; }
	inline const std::vector<double> & lContours() const { return contours_len; }
	inline const std::vector<bool> & bContours() const { return contours_select; }
	inline std::vector<bool> & bContours() { return contours_select; }
	inline bool isRelated(int i, int j) { return part_connect(i, j) > 0; }

	void extractBroads(std::vector<VSet> &inners, std::vector<VSet> &innerContours,
		std::vector<VSet> &outers, std::vector<VSet> &outerContours, int shift)
	{
		inners.clear();
		outers.clear();
		innerContours.clear();
		outerContours.clear();

		broadContours(shift);
		OpenMesh::VPropHandleT<bool> visited;
		mesh.add_property(visited);

		for (auto v : mesh.vertices())
		{
			mesh.property(visited, v) = false;
		}

		VSet part, partContour;
		for (int i = 0; i < n_contours(); i++)
		{
			for (auto v : contour(i))
			{
				if (broadInnerpart(v, part, partContour, visited, shift))
				{
					inners.push_back(part);
					innerContours.push_back(partContour);
				}
			}
		}

		for (int i = 0; i < n_contours(); i++)
		{
			for (auto v : contour(i))
			{
				mesh.property(visited, v) = false;
			}
		}

		for (int i = 0; i < n_contours(); i++)
		{
			for (auto v : contour(i))
			{
				if (broadOuterpart(v, part, partContour, visited, shift))
				{
					outers.push_back(part);
					outerContours.push_back(partContour);
				}
			}
		}

		mesh.remove_property(visited);
	}

private:
	void clearParts()
	{
		for (auto b : partEdit) b = false;
		vert_inners.clear();
		vert_contours.clear();
		for (auto outer : vert_outers) delete outer;
		vert_outers.clear();
		contours_len.clear();
		contours_select.clear();
	}

	void updateContourLength()
	{
		contours_select.resize(vert_contours.size());
		for (auto cs : contours_select) cs = false;

		contours_len.resize(vert_contours.size());
		for (int i = 0; i < contours_len.size(); i++)
		{
			contours_len[i] = 0.0;
			for (int j = 0; j < vert_contours[i].size(); j++)
			{
				auto v1 = vert_contours[i][j];
				auto v2 = vert_contours[i][(j + 1) % vert_contours[i].size()];
				auto he = mesh.find_halfedge(v1, v2);
				contours_len[i] += mesh.calc_edge_length(he);
			}
		}
	}

	void processParts(VSet &remain)
	{
		remain.clear();
		VSet vert_lines;

		OpenMesh::VPropHandleT<int> indicate;
		OpenMesh::VPropHandleT<bool> visited;
		mesh.add_property(indicate);
		mesh.add_property(visited);

		for (auto v : mesh.vertices())
		{
			bool inEdit = false;
			bool outEdit = false;
			for (auto vf : mesh.vf_range(v))
			{
				int partId = mesh.property(face_label, vf);
				inEdit = inEdit || partEdit[partId];
				outEdit = outEdit || (!partEdit[partId]);
			}
			switch ((int)inEdit - (int)outEdit)
			{
			case 1:
				vert_inners.push_back(v);
				mesh.property(indicate, v) = 1;
				break;
			case 0:
				vert_lines.push_back(v);
				mesh.property(indicate, v) = 0;
				break;
			default:
				remain.push_back(v);
				mesh.property(indicate, v) = -1;
				break;
			}
		}

		for (auto v : vert_lines)
		{
			mesh.property(visited, v) = false;
		}

		VSet stackV;
		VSet contourV;
		for (auto v : vert_lines)
		{
			if (mesh.property(visited, v)) continue;
			contourV.clear();

			mesh.property(visited, v) = true;
			stackV.push_back(v);
			contourV.push_back(v);

			while (!stackV.empty())
			{
				auto sv = stackV.back();
				stackV.pop_back();
				for (auto voh : mesh.voh_range(sv))
				{
					auto vv = mesh.to_vertex_handle(voh);
					if (mesh.property(indicate, vv)!=0 || mesh.property(visited, vv)) continue;

					if (!partEdit[mesh.property(face_label, mesh.face_handle(mesh.opposite_halfedge_handle(voh)))]
						&& partEdit[mesh.property(face_label, mesh.face_handle(voh))])
					{
						mesh.property(visited, vv) = true;
						stackV.push_back(vv);
						contourV.push_back(vv);
						break;
					}
				}
			}

			vert_contours.push_back(contourV);
		}
		updateContourLength();

		mesh.remove_property(indicate);
		mesh.remove_property(visited);
	}

	void processConnect(const VSet &remain)
	{
		Eigen::MatrixXi labelEdge(partEdit.size(), partEdit.size());
		labelEdge.setZero();

		Eigen::VectorXi labelVec(partEdit.size());
		for (auto v : mesh.vertices())
		{
			labelVec.setZero();
			for (auto vf : mesh.vf_range(v))
			{
				labelVec(mesh.property(face_label, vf)) = 1;
			}
			std::vector<int> labels;
			for (int i = 0; i < labelVec.size(); i++)
			{
				if (labelVec(i) == 1) labels.push_back(i);
			}
			for (int i = 0; i < labels.size() - 1; i++)
			{
				for (int j = i + 1; j < labels.size(); j++)
				{
					labelEdge(labels[i], labels[j]) = 1;
					labelEdge(labels[j], labels[i]) = 1;
				}
			}
		}
		
		for (int i = 0; i < partEdit.size(); i++)
		{
			labelVec(i) = -int(partEdit[i]);
		}

		int partNum = 0;
		for (auto contour : vert_contours)
		{
			partNum++;
			std::vector<int> labels;
			for (int i = 0; i < contour.size(); i++)
			{
				auto v1 = contour[i];
				auto v2 = contour[(i + 1) % contour.size()];
				auto he = mesh.find_halfedge(v2, v1);
				int label = mesh.property(face_label, mesh.face_handle(he));
				if (labelVec[label] == 0)
				{
					labels.push_back(label);
					labelVec[label] = partNum;
				}
			}
			//std::cout << labelVec << std::endl;

			for (int i = 0; i < labels.size(); i++)
			{
				for (int j = 0; j < labelEdge.cols(); j++)
				{
					if (labelVec[j] != 0) continue;
					if (labelEdge(labels[i], j) == 1)
					{
						labels.push_back(j);
						labelVec[j] = partNum;
					}
				}
			}
			//std::cout << labelVec << std::endl;
		}

		vert_outers.resize(partNum);
		for (int i = 0; i < vert_outers.size(); i++) vert_outers[i] = new VSet();
		for (auto v : remain)
		{
			auto f = *mesh.vf_begin(v);
			if (!f.is_valid()) continue;
			int partId = labelVec[mesh.property(face_label, f)] - 1;
			if (partId >= 0)
				vert_outers[partId]->push_back(v);
		}

		part_connect.resize(partNum, vert_contours.size());
		part_connect.setZero();
		for (int i = 0; i < vert_contours.size(); i++)
		{
			auto he = mesh.find_halfedge(vert_contours[i][1], vert_contours[i][0]);
			int j = labelVec[mesh.property(face_label, mesh.face_handle(he))] - 1;
			part_connect(j, i) = 1;
		}
	}

	void broadInnerVert(const OpenMesh::VertexHandle &v, int shift)
	{
		std::list<OpenMesh::VertexHandle> broad;
		if (mesh.property(dist, v) <= shift) broad.push_back(v);
		while (!broad.empty())
		{
			auto vt = broad.front();
			broad.pop_front();
			for (auto vv : mesh.vv_range(vt))
			{
				if (mesh.property(dist, vv) > mesh.property(dist, vt) + 1)
				{
					mesh.property(dist, vv) = mesh.property(dist, vt) + 1;
					if (mesh.property(dist, vv) <= shift) broad.push_back(vv);
				}
			}
		}
	}

	void broadOuterVert(const OpenMesh::VertexHandle &v, int shift)
	{
		std::list<OpenMesh::VertexHandle> broad;
		if (mesh.property(dist, v) >= -shift) broad.push_back(v);
		while (!broad.empty())
		{
			auto vt = broad.front();
			broad.pop_front();
			for (auto vv : mesh.vv_range(vt))
			{
				if (mesh.property(dist, vv) < mesh.property(dist, vt) - 1)
				{
					mesh.property(dist, vv) = mesh.property(dist, vt) - 1;
					if (mesh.property(dist, vv) >= -shift) broad.push_back(vv);
				}
			}
		}
	}

	bool broadInnerpart(const OpenMesh::VertexHandle &v, VSet &part, VSet &partcontour, 
		const OpenMesh::VPropHandleT<bool> &visited, int shift)
	{
		if (mesh.property(visited, v)) return false;
		
		part.clear();
		partcontour.clear();
		std::list<OpenMesh::VertexHandle> broad;

		broad.push_back(v);
		partcontour.push_back(v);
		mesh.property(visited, v) = true;

		while (!broad.empty())
		{
			auto vt = broad.front();
			broad.pop_front();

			for (auto vv : mesh.vv_range(vt))
			{
				if (mesh.property(visited, vv) || mesh.property(dist, vv) < 0 || mesh.property(dist, vv) > shift) continue;
				if (mesh.property(dist, vv) == 0 || mesh.property(dist, vv) == shift)
				{
					broad.push_back(vv);
					partcontour.push_back(vv);
					mesh.property(visited, vv) = true;
				}
				else
				{
					broad.push_back(vv);
					part.push_back(vv);
					mesh.property(visited, vv) = true;
				}
			}
		}

		return true;
	}

	bool broadOuterpart(const OpenMesh::VertexHandle &v, VSet &part, VSet &partcontour,
		const OpenMesh::VPropHandleT<bool> &visited, int shift)
	{
		if (mesh.property(visited, v)) return false;

		part.clear();
		partcontour.clear();
		std::list<OpenMesh::VertexHandle> broad;

		broad.push_back(v);
		partcontour.push_back(v);
		mesh.property(visited, v) = true;

		while (!broad.empty())
		{
			auto vt = broad.front();
			broad.pop_front();

			for (auto vv : mesh.vv_range(vt))
			{
				if (mesh.property(visited, vv) || mesh.property(dist, vv) > 0 || mesh.property(dist, vv) < -shift) continue;
				if (mesh.property(dist, vv) == 0 || mesh.property(dist, vv) == -shift)
				{
					broad.push_back(vv);
					partcontour.push_back(vv);
					mesh.property(visited, vv) = true;
				}
				else
				{
					broad.push_back(vv);
					part.push_back(vv);
					mesh.property(visited, vv) = true;
				}
			}
		}

		return true;
	}

	void broadContours(int shift)
	{
		for (auto v : vert_inners)
		{
			mesh.property(dist, v) = 2 * shift;
		}

		for (int i = 0; i < n_outers(); i++)
		{
			for (auto v : outer(i))
			{
				mesh.property(dist, v) = -2 * shift;
			}
		}

		for (int i = 0; i < n_contours(); i++)
		{
			for (auto v : contour(i))
			{
				mesh.property(dist, v) = 0;
			}
		}

		for (int i = 0; i < n_contours(); i++)
		{
			for (auto v : contour(i))
			{
				broadInnerVert(v, shift);
				broadOuterVert(v, shift);
			}
		}
	}

};

#endif // !SEGMESH_H
