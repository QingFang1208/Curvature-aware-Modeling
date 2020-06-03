#ifndef SEGMENTMESH_H
#define SEGMENTMESH_H

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <MeshDefinition.h>
#include <algorithm/Vec3.h>
#include <algorithm/AffineMap.hpp>
#include <algorithm/Deformation.hpp>
typedef std::vector<Mesh::VertexHandle> VSet;
typedef std::vector<Mesh::EdgeHandle> ESet;
typedef std::vector<Mesh::HalfedgeHandle> HSet;
typedef std::vector<Mesh::FaceHandle> FSet;

class SegmentMesh
{
private:
	Mesh &mesh;
	const OpenMesh::FPropHandleT<int> &face_label;
	OpenMesh::VPropHandleT<OpenMesh::VertexHandle> pair_v;
	OpenMesh::HPropHandleT<OpenMesh::HalfedgeHandle> pair_h;
	OpenMesh::VPropHandleT<int> dist;

	std::vector<bool> editing_label;

	Mesh back;
	Mesh edit;
	std::vector<Mesh *> unedit;
	std::vector<HSet> contours; // halfedge handle of edit mesh, the incident face is unvalid
	Eigen::MatrixXi label2label;
	Eigen::MatrixXi contour2label;

public:

	SegmentMesh(Mesh &input, const OpenMesh::FPropHandleT<int> &seg) : mesh(input), face_label(seg)
	{
		mesh.add_property(pair_v);
		mesh.add_property(pair_h);
		mesh.add_property(dist);
		back = mesh;
	}

	~SegmentMesh()
	{
		clear();
		mesh.remove_property(pair_v);
		mesh.remove_property(pair_h);
		mesh.remove_property(dist);
	}

	void update()
	{
		int segnum = -1;

		for (auto f : mesh.faces())
		{
			if (segnum < mesh.property(face_label, f))
				segnum = mesh.property(face_label, f);
		}
		segnum = segnum + 1;

		editing_label.resize(segnum);

		label2label.resize(segnum, segnum);
		label2label.setZero();

		for (auto e : mesh.edges())
		{
			auto f0 = mesh.face_handle(mesh.halfedge_handle(e, 0));
			auto f1 = mesh.face_handle(mesh.halfedge_handle(e, 1));
			if (f0.is_valid() && f1.is_valid())
			{
				int label0 = mesh.property(face_label, f0);
				int label1 = mesh.property(face_label, f1);
				if (label0 != label1)
				{
					label2label(label0, label1) = 1;
					label2label(label1, label0) = 1;
				}
			}
		}

	}

	void select(const std::vector<int> &labels)
	{
		clear();
		if (labels.empty()) return;
		back = mesh;
		updateEditLabel(labels);
		updateEditMesh();
		updateContours();
		updateRemainMesh();

		/*OpenMesh::IO::Options opt_tex = OpenMesh::IO::Options::FaceTexCoord;
		OpenMesh::IO::write_mesh(edit, "./edit.obj", opt_tex);

		for (int i = 0; i < unedit.size(); i++)
		{
			if (unedit[i] == NULL) continue;
			OpenMesh::IO::write_mesh(*unedit[i], "./unedit" + std::to_string(i) + ".obj", opt_tex);
		}*/
	}

	const std::vector<HSet> & hContours() const { return contours; }
	size_t n_contours() const { return contours.size(); }
	size_t n_adjMeshs() const { return unedit.size(); }
	const Mesh &curMesh() const { return edit; }
	Mesh &curMesh() { return edit; }
	const Mesh* adjMesh(size_t idx) const { return unedit[idx]; }
	Mesh* adjMesh(size_t idx) { return unedit[idx]; }
	const std::vector<Mesh *> &adjMeshs() const { return unedit; }
	const HSet& contour(size_t idx) const { return contours[idx]; }

	void assembleMesh(const std::vector<Vec3d> &poss)
	{
		initLocalPos(poss);
		initRemainPos();
	}

	void refreshMesh(const std::vector<int> &labels)
	{
		clear();
		if (labels.empty()) return;
		mesh = back;
		updateEditLabel(labels);
		updateEditMesh();
		updateContours();
		updateRemainMesh();
	}

	void updateMesh()
	{
		for (auto v : edit.vertices())
		{
			mesh.set_point(edit.property(pair_v, v), edit.point(v));
		}

		int maxDist = 5;
		
		std::vector<Mesh *> validMesh(unedit);
		for (int i = 0; i < validMesh.size(); i++)
		{
			if (validMesh[i] == NULL) continue;
			for (int j = i + 1; j < validMesh.size(); j++)
			{
				if (validMesh[j] == validMesh[i])
					validMesh[j] = NULL;
			}
			for (auto v : validMesh[i]->vertices())
			{
				validMesh[i]->property(dist, v) = maxDist;
				mesh.set_point(validMesh[i]->property(pair_v, v), OpenMesh::Vec3d(0, 0, 0));
			}
		}

		for (int i = 0; i < contours.size(); i++)
		{
			if (unedit[i] == NULL) continue;
			std::vector<VH> broad;
			for (int j = 0; j < contours[i].size(); j++)
			{
				auto edit_he = contours[i][j];
				auto unedit_he = mesh.property(pair_h, edit.property(pair_h, edit_he));
				assert(unedit[i]->property(pair_h, unedit_he).idx() == edit.property(pair_h, edit_he).idx());
				auto edit_v = edit.to_vertex_handle(edit_he);
				auto unedit_v = unedit[i]->to_vertex_handle(unedit_he);
				OpenMesh::Vec3d pos1 = edit.point(edit_v);
				OpenMesh::Vec3d pos2 = unedit[i]->point(unedit_v);
				OpenMesh::Vec3d dpos = pos1 - pos2;

				broad.push_back(unedit_v);
				unedit[i]->property(dist, unedit_v) = 0;
				mesh.set_point(unedit[i]->property(pair_v, unedit_v), dpos);
			}

			for(auto v : broad)
			{
				std::list<VH> stack;
				stack.push_back(v);

				while (!stack.empty())
				{
					auto vt = stack.front();
					stack.pop_front();
					for (auto vv : unedit[i]->vv_range(vt))
					{
						if (unedit[i]->property(dist, vv) > unedit[i]->property(dist, vt) + 1)
						{
							unedit[i]->property(dist, vv) = unedit[i]->property(dist, vt) + 1;
							if (unedit[i]->property(dist, vv) < maxDist) stack.push_back(vv);
						}
					}
				}
				
			}
		}

		for (int i = 0; i < validMesh.size(); i++)
		{
			if (validMesh[i] == NULL) continue;
			std::vector<VH> connects;
			for (auto v : validMesh[i]->vertices())
			{
				int tmp = validMesh[i]->property(dist, v);
				if (tmp > 0 && tmp < maxDist) connects.push_back(v);
			}

			std::vector<OpenMesh::Vec3d> Lp(connects.size());
			for (int j = 0; j < 5; j++)
			{
				for (int k = 0; k < connects.size(); k++)
				{
					Lp[k] = OpenMesh::Vec3d(0, 0, 0);
					int valen = validMesh[i]->valence(connects[k]);
					for (auto vv : validMesh[i]->vv_range(connects[k]))
					{
						Lp[k] += mesh.point(validMesh[i]->property(pair_v, vv));
					}
					Lp[k] /= valen;
				}

				for (int k = 0; k < connects.size(); k++)
				{
					mesh.set_point(validMesh[i]->property(pair_v, connects[k]), Lp[k]);
				}
			}
		}

		for (int i = 0; i < validMesh.size(); i++)
		{
			if (validMesh[i] == NULL) continue;
			for (auto v : validMesh[i]->vertices())
			{
				OpenMesh::Vec3d newpos = mesh.point(validMesh[i]->property(pair_v, v)) + validMesh[i]->point(v);
				mesh.set_point(validMesh[i]->property(pair_v, v), newpos);
			}
		}

		mesh.update_face_normals();
		mesh.update_vertex_normals();
	}

private:
	void clear()
	{
		for (auto v : mesh.vertices())
		{
			mesh.property(pair_v, v) = v;
		}

		for (auto he : mesh.halfedges())
		{
			mesh.property(pair_h, he) = he;
		}

		for (auto b : editing_label) b = false;
		edit.clear();
		contours.clear();
		for (int i = 0; i < unedit.size(); i++)
		{
			if (unedit[i] == NULL) continue;
			for (int j = i + 1; j < unedit.size(); j++)
			{
				if (unedit[j] == unedit[i]) unedit[j] = NULL;
			}
			delete unedit[i];
			unedit[i] = NULL;
		}
		unedit.clear();
	}

	void updateEditLabel(const std::vector<int> &editIndex)
	{
		for (auto index : editIndex)
		{
			if (index < 0 || index >= editing_label.size()) continue;
			editing_label[index] = true;
		}
	}

	void updateEditMesh()
	{
		edit = mesh;

		for (auto f : edit.faces())
		{
			if (editing_label[edit.property(face_label, f)]) continue;
			edit.delete_face(f);
		}

		edit.garbage_collection();

		for (auto he : edit.halfedges())
		{
			if (edit.is_boundary(he)) continue;
			auto he2 = edit.property(pair_h, he);
			mesh.property(pair_h, he2) = he;
		}
	}

	void updateContours()
	{
		OpenMesh::HPropHandleT<bool> visited;
		edit.add_property(visited);

		for (auto he : edit.halfedges())
		{
			edit.property(visited, he) = false;
		}

		for (auto he : edit.halfedges())
		{
			if (edit.is_boundary(he) && !edit.property(visited, he))
			{
				HSet stack = HSet(1, he);
				edit.property(visited, he) = true;

				auto hn = edit.next_halfedge_handle(he);
				while (!edit.property(visited, hn))
				{
					stack.push_back(hn);
					edit.property(visited, hn) = true;
					hn = edit.next_halfedge_handle(hn);
				}

				contours.push_back(stack);
			}
		}

		edit.remove_property(visited);

		contour2label.resize(contours.size(), editing_label.size());
		contour2label.setZero();

		for (int i = 0; i < contours.size(); i++)
		{
			for (auto he : contours[i])
			{
				auto adjface = mesh.face_handle(mesh.opposite_halfedge_handle(
					edit.property(pair_h, edit.opposite_halfedge_handle(he))));
				if (adjface.is_valid())
				{
					contour2label(i, mesh.property(face_label, adjface)) = 1;
				}
			}
		}

	}

	void updateRemainMesh()
	{
		std::vector<int> labelparts(editing_label.size(), -2);
		for (int i = 0; i < labelparts.size(); i++)
		{
			if (editing_label[i]) labelparts[i] = -1;
		}

		unedit.resize(contours.size(), NULL);

		for (int i = 0; i < contours.size(); i++)
		{
			std::vector<int> partlabels;
			for (int j = 0; j < labelparts.size(); j++)
			{
				if (contour2label(i, j) == 0) continue;
				switch (labelparts[j])
				{
				case -2:
					partlabels.push_back(j);
					labelparts[j] = i;
					break;
				case -1:
					break;
				default:
					partlabels.push_back(j);
					break;
				}
			}

			if (!partlabels.empty())
			{
				if (labelparts[partlabels[0]] != i)
				{
					unedit[i] = unedit[labelparts[partlabels[0]]];
				}
				else
				{
					for (int j = 0; j < partlabels.size(); j++)
					{
						for (int k = 0; k < labelparts.size(); k++)
						{
							if (label2label(partlabels[j], k) > 0 && labelparts[k] == -2)
							{
								partlabels.push_back(k);
								labelparts[k] = i;
							}
						}
					}

					unedit[i] = new Mesh(mesh);

					for (auto f : unedit[i]->faces())
					{
						if (labelparts[unedit[i]->property(face_label, f)] == i) continue;
						unedit[i]->delete_face(f);
					}

					unedit[i]->garbage_collection();

					for (auto he : unedit[i]->halfedges())
					{
						if (unedit[i]->is_boundary(he)) continue;
						auto he2 = unedit[i]->property(pair_h, he);
						mesh.property(pair_h, he2) = he;
					}
				}

			}
			
		}
	}

	void initLocalPos(const std::vector<Vec3d> &poss)
	{
		// AffineMap transform
		std::vector<Eigen::Vector3d> pts;
		std::vector<double> wts;
		std::vector<Eigen::Vector3d> pt0s;
		std::vector<double> wt0s;

		Eigen::Vector3d p, q, p0, q0, ct(0, 0, 0), ct0(0, 0, 0);
		double pts_w, sum_w1, sum_w2;

		sum_w1 = 0;
		sum_w2 = 0;
		for (int i = 0; i < contours.size(); i++)
		{
			for (int j = 0; j < contours[i].size(); j++)
			{
				p = Eigen::Vector3d(mesh.point(edit.property(pair_v, edit.from_vertex_handle(contours[i][j]))).data());
				q = Eigen::Vector3d(mesh.point(edit.property(pair_v, edit.to_vertex_handle(contours[i][j]))).data());
				pts_w = (p - q).norm();
				wts.push_back(pts_w);
				pts.push_back((p + q) / 2);
				ct += pts_w * (p + q) / 2;
				sum_w1 += pts_w;

				p0 = Eigen::Vector3d(poss[edit.from_vertex_handle(contours[i][j]).idx()].data());
				q0 = Eigen::Vector3d(poss[edit.to_vertex_handle(contours[i][j]).idx()].data());
				pts_w = (p0 - q0).norm();
				wt0s.push_back(pts_w);
				pt0s.push_back((p0 + q0) / 2);
				ct0 += pts_w * (p0 + q0) / 2;
				sum_w2 += pts_w;
			}
		}
		if (sum_w1 > 0 && sum_w2 > 0)
		{
			ct /= sum_w1;
			ct0 /= sum_w2;
		}

		ASAP asap;
		for (int i = 0; i < wts.size(); i++)
		{
			p = pts[i] - ct;
			p0 = pt0s[i] - ct0;
			asap.addCouples(p0, p, wts[i] + wt0s[i]); 
		}
		AffineMap tmp = asap.solve();
		tmp.T = ct - tmp.map(ct0);

		for (auto v : edit.vertices())
		{
			Eigen::Vector3d tpos(poss[v.idx()].data());
			tpos = tmp.map(tpos);
			edit.set_point(v, OpenMesh::Vec3d(tpos[0], tpos[1], tpos[2]));
		}
	}

	void initRemainPos()
	{
		std::vector<Eigen::Vector3d> pts;
		std::vector<double> wts;
		std::vector<Eigen::Vector3d> pt0s;
		std::vector<double> wt0s;
		std::vector<int> contour_ids;

		Eigen::Vector3d p, q, p0, q0, ct(0, 0, 0), ct0(0, 0, 0);
		double pts_w, sum_w1, sum_w2;

		ASAP asap;
		AffineMap tmp;

		std::vector<bool> processed(unedit.size(), false);

		for (int i = 0; i < unedit.size(); i++)
		{
			if (unedit[i] == NULL || processed[i]) continue;

			pts.clear();
			wts.clear();
			pt0s.clear();
			wt0s.clear();

			ct.setZero();
			ct0.setZero();
			sum_w1 = 0;
			sum_w2 = 0;

			contour_ids.clear();
			for (int j = i; j < unedit.size(); j++)
			{
				if (unedit[j] == unedit[i] && processed[j]==false)
				{
					contour_ids.push_back(j);
					processed[j] = true;
				}
			}

			for (auto j : contour_ids)
			{
				for (int k = 0; k < contours[j].size(); k++)
				{
					p = Eigen::Vector3d(edit.point(edit.from_vertex_handle(contours[j][k])).data());
					q = Eigen::Vector3d(edit.point(edit.to_vertex_handle(contours[j][k])).data());
					pts_w = (p - q).norm();
					wts.push_back(pts_w);
					pts.push_back((p + q) / 2);
					ct += pts_w * (p + q) / 2;
					sum_w1 += pts_w;

					p0 = Eigen::Vector3d(mesh.point(edit.property(pair_v, edit.from_vertex_handle(contours[j][k]))).data());
					q0 = Eigen::Vector3d(mesh.point(edit.property(pair_v, edit.to_vertex_handle(contours[j][k]))).data());
					pts_w = (p0 - q0).norm();
					wt0s.push_back(pts_w);
					pt0s.push_back((p0 + q0) / 2);
					ct0 += pts_w * (p0 + q0) / 2;
					sum_w2 += pts_w;
				}
			}
			if (sum_w1 > 0 && sum_w2 > 0)
			{
				ct /= sum_w1;
				ct0 /= sum_w2;
			}

			asap.clear();
			for (int j = 0; j < wts.size(); j++)
			{
				p = pts[j] - ct;
				p0 = pt0s[j] - ct0;
				asap.addCouples(p0, p, wts[j] + wt0s[j]);
			}
			tmp = asap.solve();
			tmp.T = ct - tmp.map(ct0);

			for (auto v : unedit[i]->vertices())
			{
				Eigen::Vector3d tpos(mesh.point(unedit[i]->property(pair_v, v)).data());
				tpos = tmp.map(tpos);
				unedit[i]->set_point(v, OpenMesh::Vec3d(tpos[0], tpos[1], tpos[2]));
			}
		}
	}


	void checkConnect()
	{
		for (auto c : contours)
		{
			for (auto edithe : c)
			{
				//auto uneditmesh.property(pair_h, mesh.opposite_halfedge_handle(edit.property(pair_h, edithe)))
			}
		}
	}
};

#endif // !SEGMENTMESH_H
