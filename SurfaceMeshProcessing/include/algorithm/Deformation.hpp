#ifndef DEFORMATION_H
#define DEFORMATION_H

#include <algorithm/CurvatureDeformation.hpp>

class Spherelize : public CurvatureDeformation
{
private:
	double curv_dilat;
	std::vector<double> contour_lens;

public:
	Spherelize(Mesh &input, const std::vector<std::vector<Mesh::HalfedgeHandle>> &contours) : CurvatureDeformation(input, contours)
	{
		curv_dilat = 2 * M_PI;
		progClip = 1e-3;

		contour_lens.resize(bd_contours.size(), 0);
		for (uint i = 0; i < bd_contours.size(); i++)
		{
			for (uint j = 0; j < bd_contours[i].size(); j++)
			{
				contour_lens[i] += e_Bl[bd_contours[i][j] >> 1];
			}
		}

		std::cout << "Spherelize part Vertices number : " << input.n_vertices() << std::endl;
	}

	void setCurvatureParameters(int slideValue)
	{
		curv_dilat = (slideValue - 0.0) / 10 * M_PI;
		assert(curv_dilat >= 0 && curv_dilat <= 4 * M_PI);
	}

private:
	void calabiAnchors()
	{
		for (auto v : mesh.vertices())
		{
			if (bd_mark[v.idx()]) continue;
			bool bd2 = false;
			for (auto vv : mesh.vv_range(v))
			{
				bd2 = bd2 || bd_mark[vv.idx()];
			}
			if (bd2) continue;
			this->anchorPos(Vectori(1, v.idx()), {});
			break;
		}
	}
	void deformAnchors() {}
	bool checkConditions()
	{
		if (contour_lens.empty() && curv_dilat < 4 * M_PI)
		{
			std::cout << "\nError : For close sphere, the total K_V need to be 4*PI ( SliderPara : 40 )\n";
			return false;
		}
		if (!contour_lens.empty() && curv_dilat >= 4 * M_PI)
		{
			std::cout << "\nWarning : For sphere with boundary, the total K_V need to be less than 4*PI,\
				the value is adjust to 4*PI - 1e-3\n";
			curv_dilat = 4 * M_PI - 1e-3;
		}
		return true;
	}

	void setIntialPosition(const Vectori &achs, const Vector3d &dpos)
	{
		auto v_p0 = mesh.points();
		for (uint i = 0; i < v_p.size(); i++)
		{
			v_p[i].setValue(v_p0[i][0], v_p0[i][1], v_p0[i][2]);
		}

		double *px = v_pos.data();
		double *py = v_pos.data() + var_idx.size();
		double *pz = v_pos.data() + 2 * var_idx.size();
		for (uint i = 0; i < var_idx.size(); i++)
		{
			px[i] = v_p[var_idx[i]][0];
			py[i] = v_p[var_idx[i]][1];
			pz[i] = v_p[var_idx[i]][2];
		}
	}
	void projectPosition()
	{
		/*Vec3d ct(0, 0, 0);
		for (uint i = 0; i < v_p.size(); i++)
		{
			ct += v_p[i];
		}
		ct /= v_p.size();

		double r = 0;
		for (uint i = 0; i < v_p.size(); i++)
		{
			r += (v_p[i] - ct).norm();
		}
		r /= v_p.size();

		for (uint i = 0; i < v_p.size(); i++)
		{
			v_p[i] = ct + (v_p[i] - ct).normalized()*r;
		}*/

		/*Vector3d uvw;
		calcUVWcoordinates(uvw);
		for (uint i = 0; i < v_p.size(); i++)
		{
			v_p[i] = uvw[i];
		}*/
	}
	void updatePosition() {}
	void setTargetGaussCurvature()
	{
		double areaSum = 0.0;
		for (uint i = 0; i < in_idx.size(); i++)
		{
			areaSum += v_s[in_idx[i]];
		}
		double gauss_inner_curv = curv_dilat / areaSum;

		double sqrcircleLength = 0.0;
		for (uint i = 0; i < contour_lens.size(); i++)
		{
			sqrcircleLength += log(1+contour_lens[i])/*pow(contour_lens[i], 1)*/;
		}

		/*for (uint i = 0; i < contour_lens.size(); i++)
		{
			std::cout << log(1+contour_lens[i]) / sqrcircleLength << '\t';
		}*/

		for (uint i = 0; i < in_idx.size(); i++)
		{
			v_Tg[in_idx[i]] = gauss_inner_curv * v_s[in_idx[i]];
		}

		for (uint i = 0; i < bd_idx.size(); i++)
		{
			v_Tg[bd_idx[i]] = 0;
		}

		std::vector<double> contourLen(bd_contours.size(), 0);
		std::vector<double> contourCurv(bd_contours.size(), 0);

		double tmp_el;
		for (uint i = 0; i < bd_contours.size(); i++)
		{
			contourCurv[i] = /*pow(contour_lens[i], 1)*/log(1+contour_lens[i]) / sqrcircleLength * (4 * M_PI - curv_dilat) - 2 * M_PI;
			for (uint j = 0; j < bd_contours[i].size(); j++)
			{
				tmp_el = e_l[bd_contours[i][j] >> 1];
				contourLen[i] += tmp_el;
				v_Tg[hv_idx[bd_contours[i][j]]] += tmp_el;
				v_Tg[hv_idx[bd_contours[i][j] ^ 1]] += tmp_el;
			}

			for (uint j = 0; j < bd_contours[i].size(); j++)
			{
				v_Tg[hv_idx[bd_contours[i][j]]] *= contourCurv[i] / contourLen[i] / 2;
			}
		}
	}
	void setTargetMeanCurvature()
	{
		double areaSum = 0.0;
		for (uint i = 0; i < in_idx.size(); i++)
		{
			areaSum += v_s[in_idx[i]];
		}

		double gauss_inner_curv = curv_dilat / areaSum;
		double mean_inner_curv = sqrt(gauss_inner_curv);

		for (uint i = 0; i < in_idx.size(); i++)
		{
			v_Th[in_idx[i]] = 2 * mean_inner_curv * v_s[in_idx[i]];
		}
	}
	void updateTargetGaussCurvature()
	{
		setTargetGaussCurvature();
	}
	void updateTargetMeanCurvature() {}

	typedef std::pair<Mesh::HalfedgeHandle, double> heapNode;
	struct nodeGreater
	{
		inline bool operator()(const heapNode &x, const heapNode &y) const
		{
			return (x.second > y.second);
		}
	};
	struct nodeLess
	{
		inline bool operator()(const heapNode &x, const heapNode &y) const
		{
			return (x.second < y.second);
		}
	};
	double calcHeError(const Mesh::HalfedgeHandle &he)
	{
		assert(!mesh.is_boundary(he));

		Mesh::HalfedgeHandle he_a = he;
		Mesh::HalfedgeHandle he_b = mesh.next_halfedge_handle(he_a);
		Mesh::HalfedgeHandle he_c = mesh.next_halfedge_handle(he_b);

		const Mesh::VertexHandle &v_B = mesh.from_vertex_handle(he_a);
		const Mesh::VertexHandle &v_C = mesh.to_vertex_handle(he_a);

		return fabs((v_Tg[v_B.idx()] - v_g[v_B.idx()]) * e_l[he_c.idx()/2]) +
			fabs((v_Tg[v_C.idx()] - v_g[v_C.idx()])*e_l[he_b.idx()/2]);
	}
	double calcFootPrint(Vec3d pa, double ra, Vec3d pb, double rb, Vec3d &fp)
	{
		fp = pb - pa;
		double lt = fp.norm();
		if (lt == 0)
		{
			fp = pa;
			return 0;
		}
		double dt = (lt + ra / lt * ra - rb / lt * rb) / 2;
		double dn = 1 - dt / ra * dt / ra;
		if (dn < 0) dn = 0;
		dn = ra * sqrt(dn);
		fp = pa + fp * dt / lt;
		return dn;
	}
	Vec3d calcIntersectLine(Vec3d vn, double dn, Vec3d ve)
	{
		Vec3d n1 = (vn%ve);
		n1.normalize_cond();
		Vec3d n2 = (ve%n1);
		n2.normalize_cond();

		double lambda2 = (1 - vn.sqrnorm() - dn * dn) / 2;
		if ((n2 | vn) != 0) lambda2 /= (n2 | vn);
		double tmp = 1 - lambda2 / dn * lambda2 / dn;
		if (tmp < 0) tmp = 0;
		double lambda1 = dn * sqrt(tmp);

		return vn + n1 * lambda1 + n2 * lambda2;
	}
	void calcUVWcoordinates(std::vector<Vec3d> &uvw)
	{
		uvw.resize(v_p.size());
		Vectorb face_vistor(fh_cot.size(), false);
		Vectord e_Sl(e_l.size(), 0);
		std::vector<heapNode> err_heap;

		double radius = 0.0;
		for (int i = 0; i < v_s.size(); i++)
		{
			radius += v_s[i];
		}
		radius = sqrt(radius / 4 / M_PI);
		for (int i = 0; i < e_Sl.size(); i++) e_Sl[i] = e_l[i] / radius;

		double minCalabi = FLT_MAX;
		double maxCalabi = FLT_MIN;

		for (auto he : mesh.halfedges())
		{
			double tmpCalabi = calcHeError(he);
			if (tmpCalabi < minCalabi)
			{
				minCalabi = tmpCalabi;
			}
			if (tmpCalabi > maxCalabi)
			{
				maxCalabi = tmpCalabi;
			}
		}

		double highratio = 0.7;
		double midCalabi = (1 - highratio)*minCalabi + highratio * maxCalabi;
		double difCalabi = FLT_MAX;
		Mesh::HalfedgeHandle he_(-1);
		for (auto he : mesh.halfedges())
		{
			double tmpCalabi = fabs(calcHeError(he) - midCalabi);
			if (tmpCalabi < difCalabi)
			{
				difCalabi = tmpCalabi;
				he_ = he;
			}
		}

		double l_t = e_Sl[he_.idx() / 2];
		uvw[mesh.from_vertex_handle(he_).idx()].setValue(1, 0, 0);
		uvw[mesh.to_vertex_handle(he_).idx()].setValue(1 - l_t * l_t / 2, l_t * sqrt(4 - l_t * l_t) / 2, 0);

		err_heap.push_back(std::make_pair(he_, calcHeError(he_)));
		std::push_heap(err_heap.begin(), err_heap.end(), nodeLess());
		face_vistor[mesh.face_handle(he_).idx()] = true;

		double fr, tr, dn;
		Vec3d fpos, tpos, fppos;
		OpenMesh::VertexHandle fv_t, tv_t;
		OpenMesh::HalfedgeHandle he_n, he_nn;
		while (!err_heap.empty())
		{
			he_ = err_heap.front().first;
			std::pop_heap(err_heap.begin(), err_heap.end(), nodeLess());
			err_heap.pop_back();

			fv_t = mesh.from_vertex_handle(he_);
			tv_t = mesh.to_vertex_handle(he_);
			fpos = uvw[fv_t.idx()];
			tpos = uvw[tv_t.idx()];

			he_n = mesh.next_halfedge_handle(he_);
			tr = e_Sl[mesh.edge_handle(he_n).idx()];
			he_nn = mesh.next_halfedge_handle(he_n);
			fr = e_Sl[mesh.edge_handle(he_nn).idx()];

			dn = calcFootPrint(fpos, fr, tpos, tr, fppos);
			fppos = calcIntersectLine(fppos, dn, tpos - fpos);
			uvw[mesh.to_vertex_handle(he_n).idx()].setValue(fppos[0], fppos[1], fppos[2]);

			he_ = mesh.opposite_halfedge_handle(he_);
			if (!mesh.is_boundary(he_) && !face_vistor[mesh.face_handle(he_).idx()])
			{
				err_heap.push_back(std::make_pair(he_, calcHeError(he_)));
				std::push_heap(err_heap.begin(), err_heap.end(), nodeLess());
				face_vistor[mesh.face_handle(he_).idx()] = true;
			}

			he_n = mesh.opposite_halfedge_handle(he_n);
			if (!mesh.is_boundary(he_n) && !face_vistor[mesh.face_handle(he_n).idx()])
			{
				err_heap.push_back(std::make_pair(he_n, calcHeError(he_n)));
				std::push_heap(err_heap.begin(), err_heap.end(), nodeLess());
				face_vistor[mesh.face_handle(he_n).idx()] = true;
			}

			he_nn = mesh.opposite_halfedge_handle(he_nn);
			if (!mesh.is_boundary(he_nn) && !face_vistor[mesh.face_handle(he_nn).idx()])
			{
				err_heap.push_back(std::make_pair(he_nn, calcHeError(he_nn)));
				std::push_heap(err_heap.begin(), err_heap.end(), nodeLess());
				face_vistor[mesh.face_handle(he_nn).idx()] = true;
			}
		}

	}
};

class Cylinderize : public CurvatureDeformation
{
private:
	double curv_circle;
	std::vector<bool> contour_selected;

public:
	Cylinderize(Mesh &input, const std::vector<std::vector<Mesh::HalfedgeHandle>> &contours) : CurvatureDeformation(input, contours)
	{
		progClip = 1e-4;
		curv_circle = 0;
		deformH_cap = 200;
		quas_cap = 100;
		contour_selected.resize(contours.size(), false);

		std::cout << "Cylinderize part Vertices number : " << input.n_vertices() << std::endl;
	}

	void setCurvatureParameters(int slideValue)
	{
		curv_circle = (slideValue - 20.0) / 40 * M_PI;
	}

	void selectContours(std::vector<bool> &selected)
	{
		if (contour_selected.size() == selected.size())
			contour_selected = selected;
	}

private:
	void calabiAnchors()
	{
		for (auto v : mesh.vertices())
		{
			if (bd_mark[v.idx()]) continue;
			bool bd2 = false;
			for (auto vv : mesh.vv_range(v))
			{
				bd2 = bd2 || bd_mark[vv.idx()];
			}
			if (bd2) continue;
			this->anchorPos(Vectori(1, v.idx()), {});
			break;
		}
	}
	void deformAnchors() {}
	bool checkConditions()
	{
		int count = 0;
		for (auto b : contour_selected)
			if (b) count++;
		if (count != 2)
		{
			std::cout << "\nError : For Cylinder, the upper and lower base need to be selected\n";
			return false;
		}
		return true;
	}

	void setIntialPosition(const Vectori &achs, const Vector3d &dpos)
	{
		auto v_p0 = mesh.points();
		for (uint i = 0; i < v_p.size(); i++)
		{
			v_p[i].setValue(v_p0[i][0], v_p0[i][1], v_p0[i][2]);
		}

		double *px = v_pos.data();
		double *py = v_pos.data() + var_idx.size();
		double *pz = v_pos.data() + 2 * var_idx.size();
		for (uint i = 0; i < var_idx.size(); i++)
		{
			px[i] = v_p[var_idx[i]][0];
			py[i] = v_p[var_idx[i]][1];
			pz[i] = v_p[var_idx[i]][2];
		}
	}
	void updatePosition() {}
	void projectPosition()
	{
		/*updateCotNormal();

		Eigen::Matrix3d nnT;
		Eigen::Vector3d ct;
		Eigen::Vector3d cv;

		Vec3d nm(0, 0, 0);
		double dm = 0;
		for (uint i = 0; i < in_idx.size(); i++)
		{
			nm += v_cn[in_idx[i]];
			dm += v_p[in_idx[i]] | v_cn[in_idx[i]];
		}
		nm /= in_idx.size();
		dm /= in_idx.size();

		nnT.setZero();
		ct.setZero();
		for (uint i = 0; i < in_idx.size(); i++)
		{
			Vec3d tmpcv = nm - v_cn[in_idx[i]];
			cv = Eigen::Vector3d(tmpcv[0], tmpcv[1], tmpcv[2]);
			nnT += cv * cv.transpose();
			ct += cv * (dm - (v_p[in_idx[i]] | v_cn[in_idx[i]]));
		}

		ct = nnT.inverse()*ct;

		nnT.setZero();
		for (uint i = 0; i < in_idx.size(); i++)
		{
			cv = Eigen::Vector3d(v_cn[in_idx[i]][0], v_cn[in_idx[i]][1], v_cn[in_idx[i]][2]);
			nnT = nnT + cv * cv.transpose();
		}

		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(nnT);
		cv = solver.eigenvectors().col(0);
		cv.normalize();

		Vec3d c_t(ct[0], ct[1], ct[2]), c_v(cv[0], cv[1], cv[2]);

		double r = 0;
		for (uint i = 0; i < v_p.size(); i++)
		{
			Vec3d dp = v_p[i] - c_t;
			r += (dp - c_v * (dp | c_v)).norm();
		}
		r /= v_p.size();

		for (uint i = 0; i < v_p.size(); i++)
		{
			Vec3d dp = v_p[i] - c_t;
			Vec3d nt = c_t + c_v * (dp | c_v);
			v_p[i] = nt + (v_p[i] - nt).normalized()*r;
		}*/
	}
	void setTargetGaussCurvature()
	{
		for (uint i = 0; i < in_idx.size(); i++)
		{
			v_Tg[in_idx[i]] = 0.0;
		}

		for (uint i = 0; i < bd_idx.size(); i++)
		{
			v_Tg[bd_idx[i]] = 0;
		}

		double tmp_el;
		int firstContour = 0;
		for (uint i = 0; i < contour_selected.size(); i++)
		{
			if (contour_selected[i])
			{
				double contourCurv = firstContour == 0 ? curv_circle : -curv_circle;
				double contourLen = 0;
				for (uint j = 0; j < bd_contours[i].size(); j++)
				{
					tmp_el = e_l[bd_contours[i][j] >> 1];
					contourLen += tmp_el;
					v_Tg[hv_idx[bd_contours[i][j]]] += tmp_el;
					v_Tg[hv_idx[bd_contours[i][j] ^ 1]] += tmp_el;
				}
				for (uint j = 0; j < bd_contours[i].size(); j++)
				{
					v_Tg[hv_idx[bd_contours[i][j]]] *= contourCurv / contourLen / 2;
				}
				firstContour++;
			}
			else
			{
				double contourCurv = -2 * M_PI;
				double contourLen = 0;
				for (uint j = 0; j < bd_contours[i].size(); j++)
				{
					tmp_el = e_l[bd_contours[i][j] >> 1];
					contourLen += tmp_el;
					v_Tg[hv_idx[bd_contours[i][j]]] += tmp_el;
					v_Tg[hv_idx[bd_contours[i][j] ^ 1]] += tmp_el;
				}
				for (uint j = 0; j < bd_contours[i].size(); j++)
				{
					v_Tg[hv_idx[bd_contours[i][j]]] *= contourCurv / contourLen / 2;
				}
			}
			
		}
	}
	void setTargetMeanCurvature()
	{
		double contoursLen = 0.0;
		for (uint i = 0; i < contour_selected.size(); i++)
		{
			if (!contour_selected[i]) continue;
			for (uint j = 0; j < bd_contours[i].size(); j++)
			{
				contoursLen += e_l[bd_contours[i][j] >> 1];
			}
		}

		double mean_inner_curv = 2 * M_PI / contoursLen;

		for (uint i = 0; i < in_idx.size(); i++)
		{
			v_Th[in_idx[i]] = 2 * mean_inner_curv * v_s[in_idx[i]];
		}
	}
	void updateTargetGaussCurvature()
	{
		setTargetGaussCurvature();
	}
	void updateTargetMeanCurvature() {}
};

class PlanePolygon : public CurvatureDeformation
{
private:
	int polygon_num;
	int contour_num;

public:
	PlanePolygon(Mesh &input, const std::vector<std::vector<Mesh::HalfedgeHandle>> &contours) : CurvatureDeformation(input, contours)
	{
		polygon_num = 4;
		deformH_cap = 1000;
		lambdaH = 1;
		contour_num = contours.size();

		std::cout << "PlanePolygon part Vertices number : " << input.n_vertices() << std::endl;
	}

	void setCurvatureParameters(int slideValue)
	{
		polygon_num = slideValue / 10 + 3;
	}

private:
	void calabiAnchors()
	{
		for (auto v : mesh.vertices())
		{
			if (bd_mark[v.idx()]) continue;
			bool bd2 = false;
			for (auto vv : mesh.vv_range(v))
			{
				bd2 = bd2 || bd_mark[vv.idx()];
			}
			if (bd2) continue;
			this->anchorPos(Vectori(1, v.idx()), {});
			break;
		}
	}
	void deformAnchors() {}
	bool checkConditions()
	{
		if (contour_num != 1)
		{
			std::cout << "\nError : For PlanePolygon, the contour number must be one\n";
			return false;
		}
		return true;
	}
	void setIntialPosition(const Vectori &achs, const Vector3d &dpos) {}
	void updatePosition()
	{
		std::vector<OpenMesh::Vec2d> uvs;
		calc_uv_coordinates(uvs);
		for (uint i = 0; i < v_p.size(); i++)
		{
			v_p[i].setValue(uvs[i][0], uvs[i][1], 1);
		}

		double *px = v_pos.data();
		double *py = v_pos.data() + var_idx.size();
		double *pz = v_pos.data() + 2 * var_idx.size();
		for (uint i = 0; i < var_idx.size(); i++)
		{
			px[i] = v_p[var_idx[i]][0];
			py[i] = v_p[var_idx[i]][1];
			pz[i] = v_p[var_idx[i]][2];
		}
	}
	void projectPosition() {}
	void setTargetGaussCurvature()
	{
		for (uint i = 0; i < in_idx.size(); i++)
		{
			v_Tg[in_idx[i]] = 0;
		}

		for (uint i = 0; i < bd_idx.size(); i++)
		{
			v_Tg[bd_idx[i]] = 0;
		}

		int corner = (bd_contours[0].size() + 1) / polygon_num;

		for (uint i = 0; i < bd_contours[0].size(); i++)
		{
			if (i%corner == 0)
				v_Tg[hv_idx[bd_contours[0][i]]] = 2 * M_PI / polygon_num;
		}

	}
	void setTargetMeanCurvature()
	{
		for (uint i = 0; i < in_idx.size(); i++)
		{
			v_Th[in_idx[i]] = 0;
		}
	}
	void updateTargetGaussCurvature()
	{
		setTargetGaussCurvature();
	}
	void updateTargetMeanCurvature() {}

	void calc_uv_coordinates(std::vector<OpenMesh::Vec2d> &uvs)
	{
		std::list<uint> stack_hidx;
		Vectorb visited(fh_idx.size(), false);

		uvs.clear();
		uvs.resize(v_p.size(), OpenMesh::Vec2d(0, 0));

		stack_hidx.push_back(fh_idx[0][0]);
		visited[0] = true;
		uvs[hv_idx[fh_idx[0][0]]] = OpenMesh::Vec2d(e_l[fh_idx[0][0] >> 1], 0);

		double fr, tr, dt, dn, l_t;
		OpenMesh::Vec2d fpos, tpos, dpos;
		uint fv_idx, tv_idx, hn_idx, hnn_idx;
		int f_idx;
		
		while (!stack_hidx.empty())
		{
			uint idx_t = stack_hidx.front();
			stack_hidx.pop_front();

			fpos = uvs[hv_idx[idx_t ^ 1]];
			tpos = uvs[hv_idx[idx_t]];
			l_t = e_l[idx_t >> 1];

			hn_idx = mesh.next_halfedge_handle(HH(idx_t)).idx();
			tr = e_l[hn_idx >> 1];
			hnn_idx = mesh.next_halfedge_handle(HH(hn_idx)).idx();
			fr = e_l[hnn_idx >> 1];

			dt = (l_t + fr / l_t * fr - tr / l_t * tr) / 2;
			dn = 1 - dt / fr * dt / fr;
			assert(dn >= 0);
			dn = fr * sqrt(dn);
			dpos = (tpos - fpos).normalize_cond();
			dpos = fpos + dt * dpos + dn * OpenMesh::Vec2d(-dpos[1], dpos[0]);
			uvs[hv_idx[hn_idx]] = dpos;

			f_idx = mesh.face_handle(HH(idx_t ^ 1)).idx();
			if (f_idx >=0 && !visited[f_idx])
			{
				stack_hidx.push_back(idx_t ^ 1);
				visited[f_idx] = true;
			}

			f_idx = mesh.face_handle(HH(hn_idx ^ 1)).idx();
			if (f_idx >= 0 && !visited[f_idx])
			{
				stack_hidx.push_back(hn_idx ^ 1);
				visited[f_idx] = true;
			}
			
			f_idx = mesh.face_handle(HH(hnn_idx ^ 1)).idx();
			if (f_idx >= 0 && !visited[f_idx])
			{
				stack_hidx.push_back(hnn_idx ^ 1);
				visited[f_idx] = true;
			}
		}

	}

};

class PlaneFree : public CurvatureDeformation
{
	int contour_num;

public:
	PlaneFree(Mesh &input, const std::vector<std::vector<Mesh::HalfedgeHandle>> &contours) : CurvatureDeformation(input, contours)
	{
		contour_num = contours.size();
		std::cout << "PlaneFree part Vertices number : " << input.n_vertices() << std::endl;
	}

	void setCurvatureParameters(int slideValue) {}

private:
	void calabiAnchors()
	{
		for (auto v : mesh.vertices())
		{
			if (bd_mark[v.idx()]) continue;
			bool bd2 = false;
			for (auto vv : mesh.vv_range(v))
			{
				bd2 = bd2 || bd_mark[vv.idx()];
			}
			if (bd2) continue;
			this->anchorPos(Vectori(1, v.idx()), {});
			break;
		}
	}
	void deformAnchors()
	{
		this->anchorPos(Vectori(1, 0), {});
	}
	bool checkConditions()
	{
		if (contour_num != 1)
		{
			std::cout << "\nError : For PlanePolygon, the contour number must be one\n";
			return false;
		}
		return true;
	}
	void setIntialPosition(const Vectori &achs, const Vector3d &dpos) {}
	void updatePosition()
	{
		std::vector<OpenMesh::Vec2d> uvs;
		calc_uv_coordinates(uvs);
		for (uint i = 0; i < v_p.size(); i++)
		{
			v_p[i].setValue(uvs[i][0], uvs[i][1], 0);
		}

		double *px = v_pos.data();
		double *py = v_pos.data() + var_idx.size();
		double *pz = v_pos.data() + 2 * var_idx.size();
		for (uint i = 0; i < var_idx.size(); i++)
		{
			px[i] = v_p[var_idx[i]][0];
			py[i] = v_p[var_idx[i]][1];
			pz[i] = v_p[var_idx[i]][2];
		}
	}
	void projectPosition() {}
	void setTargetGaussCurvature()
	{
		for (uint i = 0; i < in_idx.size(); i++)
		{
			v_Tg[in_idx[i]] = 0;
		}

		for (uint i = 0; i < bd_idx.size(); i++)
		{
			v_Tg[bd_idx[i]] = 0;
		}

	}
	void setTargetMeanCurvature()
	{
		for (uint i = 0; i < in_idx.size(); i++)
		{
			v_Th[in_idx[i]] = 0;
		}
	}
	void updateTargetGaussCurvature()
	{
		setTargetGaussCurvature();
	}
	void updateTargetMeanCurvature() {}

	void calc_uv_coordinates(std::vector<OpenMesh::Vec2d> &uvs)
	{
		std::list<uint> stack_hidx;
		Vectorb visited(fh_idx.size(), false);

		uvs.clear();
		uvs.resize(v_p.size(), OpenMesh::Vec2d(0, 0));

		stack_hidx.push_back(fh_idx[0][0]);
		visited[0] = true;
		uvs[hv_idx[fh_idx[0][0]]] = OpenMesh::Vec2d(e_l[fh_idx[0][0] >> 1], 0);

		double fr, tr, dt, dn, l_t;
		OpenMesh::Vec2d fpos, tpos, dpos;
		uint fv_idx, tv_idx, hn_idx, hnn_idx;
		int f_idx;

		while (!stack_hidx.empty())
		{
			uint idx_t = stack_hidx.front();
			stack_hidx.pop_front();

			fpos = uvs[hv_idx[idx_t ^ 1]];
			tpos = uvs[hv_idx[idx_t]];
			l_t = e_l[idx_t >> 1];

			hn_idx = mesh.next_halfedge_handle(HH(idx_t)).idx();
			tr = e_l[hn_idx >> 1];
			hnn_idx = mesh.next_halfedge_handle(HH(hn_idx)).idx();
			fr = e_l[hnn_idx >> 1];

			dt = (l_t + fr / l_t * fr - tr / l_t * tr) / 2;
			dn = 1 - dt / fr * dt / fr;
			assert(dn >= 0);
			dn = fr * sqrt(dn);
			dpos = (tpos - fpos).normalize_cond();
			dpos = fpos + dt * dpos + dn * OpenMesh::Vec2d(-dpos[1], dpos[0]);
			uvs[hv_idx[hn_idx]] = dpos;

			f_idx = mesh.face_handle(HH(idx_t ^ 1)).idx();
			if (f_idx >= 0 && !visited[f_idx])
			{
				stack_hidx.push_back(idx_t ^ 1);
				visited[f_idx] = true;
			}

			f_idx = mesh.face_handle(HH(hn_idx ^ 1)).idx();
			if (f_idx >= 0 && !visited[f_idx])
			{
				stack_hidx.push_back(hn_idx ^ 1);
				visited[f_idx] = true;
			}

			f_idx = mesh.face_handle(HH(hnn_idx ^ 1)).idx();
			if (f_idx >= 0 && !visited[f_idx])
			{
				stack_hidx.push_back(hnn_idx ^ 1);
				visited[f_idx] = true;
			}
		}

	}

};

class SpherePolygon : public CurvatureDeformation
{
private:
	int polygon_num;
	int contour_num;
	double curv_patch;

public:
	SpherePolygon(Mesh &input, const std::vector<std::vector<Mesh::HalfedgeHandle>> &contours) : CurvatureDeformation(input, contours)
	{
		polygon_num = 4;
		contour_num = contours.size();
		curv_patch = M_PI;

		std::cout << "PlanePolygon part Vertices number : " << input.n_vertices() << std::endl;
	}

	void setCurvatureParameters(int slideValue)
	{
		polygon_num = slideValue / 10 + 3;
	}

private:
	void calabiAnchors()
	{
		for (auto v : mesh.vertices())
		{
			if (bd_mark[v.idx()]) continue;
			bool bd2 = false;
			for (auto vv : mesh.vv_range(v))
			{
				bd2 = bd2 || bd_mark[vv.idx()];
			}
			if (bd2) continue;
			this->anchorPos(Vectori(1, v.idx()), {});
			break;
		}
	}
	void deformAnchors() {}
	bool checkConditions()
	{
		if (contour_num != 1)
		{
			std::cout << "\nError : For PlanePolygon, the contour number must be one\n";
			return false;
		}
		return true;
	}
	void setIntialPosition(const Vectori &achs, const Vector3d &dpos)
	{
		auto v_p0 = mesh.points();
		for (uint i = 0; i < v_p.size(); i++)
		{
			v_p[i].setValue(v_p0[i][0], v_p0[i][1], v_p0[i][2]);
		}

		double *px = v_pos.data();
		double *py = v_pos.data() + var_idx.size();
		double *pz = v_pos.data() + 2 * var_idx.size();
		for (uint i = 0; i < var_idx.size(); i++)
		{
			px[i] = v_p[var_idx[i]][0];
			py[i] = v_p[var_idx[i]][1];
			pz[i] = v_p[var_idx[i]][2];
		}
	}
	void updatePosition() {}
	void projectPosition() {}
	void setTargetGaussCurvature()
	{
		double areaSum = 0.0;
		double areaCorner = 0.0;

		for (uint i = 0; i < in_idx.size(); i++)
		{
			areaSum += v_s[in_idx[i]];
		}

		int corner = (bd_contours[0].size() + 1) / polygon_num;
		for (uint i = 0; i < bd_contours[0].size(); i++)
		{
			if (i%corner == 0)
			{
				areaCorner += v_s[hv_idx[bd_contours[0][i]]];
			}
			areaSum += v_s[hv_idx[bd_contours[0][i]]];
		}

		for (uint i = 0; i < in_idx.size(); i++)
		{
			v_Tg[in_idx[i]] = curv_patch / areaSum * v_s[in_idx[i]];
		}
		for (uint i = 0; i < bd_contours[0].size(); i++)
		{
			if (i%corner == 0)
			{
				v_Tg[hv_idx[bd_contours[0][i]]] =
					(2 * M_PI - curv_patch*(1 - areaCorner / areaSum)) / polygon_num;
			}
			else
			{
				v_Tg[hv_idx[bd_contours[0][i]]] =
					curv_patch / areaSum * v_s[hv_idx[bd_contours[0][i]]];
			}
		}

	}
	void setTargetMeanCurvature()
	{
		double areaSum = 0.0;
		for (uint i = 0; i < in_idx.size(); i++)
		{
			areaSum += v_s[in_idx[i]];
		}

		int corner = (bd_contours[0].size() + 1) / polygon_num;
		for (uint i = 0; i < bd_contours[0].size(); i++)
		{
			if (i%corner != 0)
			{
				areaSum += v_s[hv_idx[bd_contours[0][i]]];
			}
		}

		double gauss_inner_curv = curv_patch / areaSum;
		double mean_inner_curv = sqrt(gauss_inner_curv);

		for (uint i = 0; i < in_idx.size(); i++)
		{
			v_Th[in_idx[i]] = 2 * mean_inner_curv * v_s[in_idx[i]];
		}
	}
	void updateTargetGaussCurvature()
	{
		setTargetGaussCurvature();
	}
	void updateTargetMeanCurvature() {}

};

class SmoothControl : public CurvatureDeformation
{
private:
	double smoothness;

public:
	SmoothControl(Mesh &input, const std::vector<std::vector<Mesh::HalfedgeHandle>> &contours) : CurvatureDeformation(input, contours)
	{
		deformH_cap = 30;
		deformG_cap = 0;
		smoothness = 0;
		std::cout << "SmoothControl part Vertices number : " << input.n_vertices() << std::endl;
	}

	void setCurvatureParameters(int slideValue)
	{
		smoothness = (slideValue - 20.0) / 20;
	}

private:
	void calabiAnchors()
	{
		for (auto v : mesh.vertices())
		{
			if (bd_mark[v.idx()]) continue;
			bool bd2 = false;
			for (auto vv : mesh.vv_range(v))
			{
				bd2 = bd2 || bd_mark[vv.idx()];
			}
			if (bd2) continue;
			this->anchorPos(Vectori(1, v.idx()), {});
			break;
		}
	}
	void deformAnchors() {}
	bool checkConditions()
	{
		return true;
	}
	void setIntialPosition(const Vectori &achs, const Vector3d &dpos)
	{
		auto v_p0 = mesh.points();
		for (uint i = 0; i < v_p.size(); i++)
		{
			v_p[i].setValue(v_p0[i][0], v_p0[i][1], v_p0[i][2]);
		}

		double *px = v_pos.data();
		double *py = v_pos.data() + var_idx.size();
		double *pz = v_pos.data() + 2 * var_idx.size();
		for (uint i = 0; i < var_idx.size(); i++)
		{
			px[i] = v_p[var_idx[i]][0];
			py[i] = v_p[var_idx[i]][1];
			pz[i] = v_p[var_idx[i]][2];
		}
	}
	void updatePosition() {}
	void projectPosition() {}
	void setTargetGaussCurvature()
	{
		for (uint i = 0; i < v_Tg.size(); i++)
		{
			v_Tg[i] = v_g[i];
		}

	}
	void setTargetMeanCurvature()
	{
		double minh, curh;
		updateMeanCurvature();

		for (uint i = 0; i < in_idx.size(); i++)
		{
			if (v_g[in_idx[i]] >= 0)
			{
				minh = 2 * sqrt(v_g[in_idx[i]] * v_s[in_idx[i]]);
				curh = v_h[in_idx[i]];
				v_Th[in_idx[i]] = curh + smoothness * (curh - minh);
			}
			else
			{
				minh = 0;
				curh = v_h[in_idx[i]];
				v_Th[in_idx[i]] = curh + smoothness * (curh - minh);
			}
		}
	}
	void updateTargetGaussCurvature() {}
	void updateTargetMeanCurvature() {}

};

class FeatureControl : public CurvatureDeformation
{
private:
	double feature;

public:
	FeatureControl(Mesh &input, const std::vector<std::vector<Mesh::HalfedgeHandle>> &contours) : CurvatureDeformation(input, contours)
	{
		feature = 1;
		std::cout << "FeatureControl part Vertices number : " << input.n_vertices() << std::endl;
	}

	void setCurvatureParameters(int slideValue)
	{
		feature = exp((slideValue - 20.0) / 20 * log(1.2));
	}

private:
	void calabiAnchors()
	{
		for (auto v : mesh.vertices())
		{
			if (bd_mark[v.idx()]) continue;
			bool bd2 = false;
			for (auto vv : mesh.vv_range(v))
			{
				bd2 = bd2 || bd_mark[vv.idx()];
			}
			if (bd2) continue;
			this->anchorPos(Vectori(1, v.idx()), {});
			break;
		}
	}
	void deformAnchors() {}
	bool checkConditions()
	{
		return true;
	}
	void setIntialPosition(const Vectori &achs, const Vector3d &dpos)
	{
		auto v_p0 = mesh.points();
		for (uint i = 0; i < v_p.size(); i++)
		{
			v_p[i].setValue(v_p0[i][0], v_p0[i][1], v_p0[i][2]);
		}

		double *px = v_pos.data();
		double *py = v_pos.data() + var_idx.size();
		double *pz = v_pos.data() + 2 * var_idx.size();
		for (uint i = 0; i < var_idx.size(); i++)
		{
			px[i] = v_p[var_idx[i]][0];
			py[i] = v_p[var_idx[i]][1];
			pz[i] = v_p[var_idx[i]][2];
		}
	}
	void updatePosition() {}
	void projectPosition() {}
	void setTargetGaussCurvature()
	{
		double areaSum = 0.0;
		double gSum = 0;
		for (uint i = 0; i < in_idx.size(); i++)
		{
			gSum += v_g[in_idx[i]];
			areaSum += v_s[in_idx[i]];
		}

		double gMean;
		for (uint i = 0; i < in_idx.size(); i++)
		{
			gMean = gSum / areaSum * v_s[in_idx[i]];
			v_Tg[in_idx[i]] = gMean + feature * (v_g[in_idx[i]] - gMean);
		}

		for (uint i = 0; i < bd_idx.size(); i++)
		{
			v_Tg[bd_idx[i]] = v_g[bd_idx[i]];
		}
		
	}
	void setTargetMeanCurvature()
	{
		updateMeanCurvature();

		for (uint i = 0; i < in_idx.size(); i++)
		{
			v_Th[in_idx[i]] = v_h[in_idx[i]];
		}
	}
	void updateTargetGaussCurvature() {}
	void updateTargetMeanCurvature() {}

};

class DragVertices : public CurvatureDeformation
{
public:
	DragVertices(Mesh &input, const std::vector<std::vector<Mesh::HalfedgeHandle>> &contours) : CurvatureDeformation(input, contours)
	{
		lambdaH = 1;
		deformH_cap = 200;
		deformG_cap = 0;
		std::cout << "DragVertices part Vertices number : " << input.n_vertices() << std::endl;
	}

	void setCurvatureParameters(int slideValue)
	{
		//lambdaH = exp((slideValue - 20.0) / 20 * log(2.0));
	}

private:
	void calabiAnchors() {}
	void deformAnchors() {}
	bool checkConditions()
	{
		if (ach_idx.size() <= 0)
		{
			std::cout << "\nError : For DragVertices, the anchors vertices must be set\n";
			return false;
		}
		return true;
	}
	void setIntialPosition(const Vectori &achs, const Vector3d &dpos)
	{
		auto v_p0 = mesh.points();
		for (uint i = 0; i < v_p.size(); i++)
		{
			v_p[i].setValue(v_p0[i][0], v_p0[i][1], v_p0[i][2]);
		}
		m_achs = achs;
		m_dpos = dpos;
	}
	void updatePosition()
	{
		updateMeanCurvature();
		for (uint i = 0; i < m_achs.size(); i++)
		{
			v_p[m_achs[i]] += m_dpos[i];
		}

		double *px = v_pos.data();
		double *py = v_pos.data() + var_idx.size();
		double *pz = v_pos.data() + 2 * var_idx.size();
		for (uint i = 0; i < var_idx.size(); i++)
		{
			px[i] = v_p[var_idx[i]][0];
			py[i] = v_p[var_idx[i]][1];
			pz[i] = v_p[var_idx[i]][2];
		}
	}
	void projectPosition() {}
	void setTargetGaussCurvature()
	{
		for (uint i = 0; i < v_g.size(); i++)
		{
			v_Tg[i] = v_g[i];
		}
	}
	void setTargetMeanCurvature()
	{
		for (uint i = 0; i < v_h.size(); i++)
		{
			v_Th[i] = v_h[i];
		}
	}
	void updateTargetGaussCurvature() {}
	void updateTargetMeanCurvature() {}

private:
	Vectori m_achs;
	Vector3d m_dpos;

};

class Test : public CurvatureDeformation
{
public:
	Test(Mesh &input, const std::vector<std::vector<Mesh::HalfedgeHandle>> &contours) : CurvatureDeformation(input, contours)
	{
		lambdaH = 500;
		deformH_cap = 20;
		std::cout << "Test part Vertices number : " << input.n_vertices() << std::endl;
	}

	void setCurvatureParameters(int slideValue) {}

private:
	void calabiAnchors()
	{
		this->anchorPos(Vectori(1, v_p.size()/2), {});
	}
	void deformAnchors() {}
	bool checkConditions()
	{
		return true;
	}

	void setIntialPosition(const Vectori &achs, const Vector3d &dpos)
	{
		auto v_p0 = mesh.points();
		for (uint i = 0; i < v_p.size(); i++)
		{
			v_p[i].setValue(v_p0[i][0], v_p0[i][1], v_p0[i][2]);
		}
	}
	void projectPosition() {}
	void updatePosition() 
	{
		auto v_p0 = mesh.points();
		for (uint i = 0; i < v_p.size(); i++)
		{
			v_p[i].setValue(v_p0[i][0], v_p0[i][1], v_p0[i][2]);
			//v_p[i].setValue(v_p0[i][0], v_p0[i][1], sqrt(100 - v_p0[i][0] * v_p0[i][0] - v_p0[i][1] * v_p0[i][1]));
			/*if(v_p0[i][2] > 0)
				v_p[i].setValue(v_p0[i][0], v_p0[i][1], 0.001);
			else
				v_p[i].setValue(v_p0[i][0], v_p0[i][1], -0.001);*/
		}

		double *px = v_pos.data();
		double *py = v_pos.data() + var_idx.size();
		double *pz = v_pos.data() + 2 * var_idx.size();
		for (uint i = 0; i < var_idx.size(); i++)
		{
			px[i] = v_p[var_idx[i]][0];
			py[i] = v_p[var_idx[i]][1];
			pz[i] = v_p[var_idx[i]][2];
		}
	}
	void setTargetGaussCurvature()
	{
		//updateMeanCurvature();
		loadCurvature("./data/complexity/huaping/curvature2.txt");
		/*for (uint i = 0; i < v_g.size(); i++)
		{
			v_Tg[i] = v_g[i];
		}
		writeCurvature("./data/complexity/huaping/curvature2.txt");*/
	}
	void setTargetMeanCurvature()
	{
		/*for (uint i = 0; i < v_g.size(); i++)
		{
			v_Th[i] = v_h[i];
		}*/
	}
	void updateTargetGaussCurvature() {}
	void updateTargetMeanCurvature() {}

	void writeCurvature(std::string filename)
	{
		std::ofstream file(filename);
		if (file)
		{
			for (uint i = 0; i < v_g.size(); i++)
			{
				file << v_g[i] << ' ' << v_h[i] << std::endl;
			}
		}
		file.close();
	}

	void loadCurvature(std::string filename)
	{
		std::ifstream file(filename);
		if (file)
		{
			for (uint i = 0; i < v_g.size(); i++)
			{
				file >> v_Tg[i] >> v_Th[i];
			}
		}
		file.close();
	}
};

#endif // !DEFORMATION_H
