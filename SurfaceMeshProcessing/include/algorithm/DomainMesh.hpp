#ifndef DOMAINMESH_H
#define DOMAINMESH_H

#include <algorithm/Vec3.h>
#include <MeshDefinition.h>
#include <omp.h>

using Vectorb = std::vector<bool>;
using Vectori = std::vector<int>;
using Vectord = std::vector<double>;
using Vector3i = std::vector<Vec3i>;
using Vector3d = std::vector<Vec3d>;
using Vector33d = std::vector<Array3<Vec3d>>;

using VH = OpenMesh::VertexHandle;
using EH = OpenMesh::EdgeHandle;
using FH = OpenMesh::FaceHandle;
using HH = OpenMesh::HalfedgeHandle;

#define COSTHRES 1.0
#define THETAFLIP (M_PI + 5e-2)

// connected domian for optmization
class DomainMesh
{
protected:
	Mesh &mesh;
	Mesh back;

	Vectorb e_b, bd_mark;
	Vector3i fh_idx;
	Vectord v_u, v_g, v_h, v_Bu, v_Tg, v_Th, v_s, e_l, e_Bl, e_theta;
	Vector3d v_p, v_cn, v_n, fh_theta, fh_cot;
	Vector33d fh_d;

	Vectori hv_idx;

public:
	Vectori bd_idx;
	Vectori in_idx;
	std::vector<Vectori> bd_contours;


public:
	DomainMesh(Mesh &input, const std::vector<std::vector<Mesh::HalfedgeHandle>> &contours) : mesh(input), back(mesh)
	{
		size_t vn = mesh.n_vertices();
		size_t en = mesh.n_edges();
		size_t fn = mesh.n_faces();

		e_b = Vectorb(en, false);
		e_l = e_Bl = e_theta = Vectord(en, 0);
		hv_idx = Vectori(2 * en, -1);

		v_u = v_g = v_h = v_s = Vectord(vn, 0);
		v_Bu = v_Tg = v_Th = Vectord(vn, 0);
		v_p = v_cn  = v_n = Vector3d(vn, (Vec3d)0);

		fh_idx = Vector3i(fn, Vec3i(-1));
		fh_theta = fh_cot = Vector3d(fn, (Vec3d)0);
		fh_d = Vector33d(fn, Array3<Vec3d>((Vec3d)0));

		initDomainMesh();

		bd_idx = {};
		in_idx = {};
		in_idx.reserve(v_u.size());
		bd_mark.resize(v_p.size(), false);

		for (int i = 0; i < contours.size(); i++)
		{
			Vectori bd_contour(contours[i].size(), 0);
			for (int j = 0; j < contours[i].size(); j++)
			{
				bd_contour[j] = contours[i][j].idx();
				bd_mark[hv_idx[bd_contour[j]]] = true;
			}
			bd_contours.push_back(bd_contour);
		}

		for (int i = 0; i < v_p.size(); i++)
		{
			if (bd_mark[i]) bd_idx.push_back(i);
			else in_idx.push_back(i);
		}
	}

private:

	inline void calcTheta(double la, double lb, double lc, double &A)
	{
		A = acos((lb / lc + lc / lb - la / lb * la / lc) / 2);
	}

	bool attemptFlip(int idx, std::vector<int> &add)
	{
		auto hcd = mesh.halfedge_handle(idx << 1);
		auto hda = mesh.next_halfedge_handle(hcd);
		auto hac = mesh.next_halfedge_handle(hda);

		auto hdc = mesh.halfedge_handle((idx << 1) ^ 1);
		auto hcb = mesh.next_halfedge_handle(hdc);
		auto hbd = mesh.next_halfedge_handle(hcb);

		auto va = mesh.to_vertex_handle(hda);
		auto vb = mesh.to_vertex_handle(hcb);
		double ua = v_u[va.idx()];
		double ub = v_u[vb.idx()];

		double lab0 = (mesh.point(va) - mesh.point(vb)).norm();
		double lab = exp(ua + ub)*lab0;
		double lbd = e_l[hbd.idx() >> 1];
		double lda = e_l[hda.idx() >> 1];
		double lac = e_l[hac.idx() >> 1];
		double lcb = e_l[hcb.idx() >> 1];

		bool flip_ok = lab + lbd > lda && lab + lda > lbd && lbd + lda > lab
			&& lab + lcb > lac && lab + lac > lcb && lac + lcb > lab;
		if (!flip_ok) return false;

		double abd, bda, dab, cba, bac, acb;
		calcTheta(lab, lbd, lda, bda);
		calcTheta(lbd, lda, lab, dab);
		calcTheta(lda, lab, lbd, abd);
		calcTheta(lcb, lab, lac, bac);
		calcTheta(lab, lac, lcb, acb);
		calcTheta(lac, lcb, lab, cba);
		
		int f1 = mesh.face_handle(hcd).idx();
		int f2 = mesh.face_handle(hdc).idx();

		e_theta[fh_idx[f1][0] >> 1] -= fh_theta[f1][0];
		e_theta[fh_idx[f1][1] >> 1] -= fh_theta[f1][1];
		e_theta[fh_idx[f1][2] >> 1] -= fh_theta[f1][2];
		e_theta[fh_idx[f2][0] >> 1] -= fh_theta[f2][0];
		e_theta[fh_idx[f2][1] >> 1] -= fh_theta[f2][1];
		e_theta[fh_idx[f2][2] >> 1] -= fh_theta[f2][2];

		/*v_g[mesh.to_vertex_handle(HH(fh_idx[f1][0])).idx()] += fh_theta[f1][2];
		v_g[mesh.to_vertex_handle(HH(fh_idx[f1][1])).idx()] += fh_theta[f1][0];
		v_g[mesh.to_vertex_handle(HH(fh_idx[f1][2])).idx()] += fh_theta[f1][1];
		v_g[mesh.to_vertex_handle(HH(fh_idx[f2][0])).idx()] += fh_theta[f2][2];
		v_g[mesh.to_vertex_handle(HH(fh_idx[f2][1])).idx()] += fh_theta[f2][0];
		v_g[mesh.to_vertex_handle(HH(fh_idx[f2][2])).idx()] += fh_theta[f2][1];

		double l[3];
		double varea[3];
		double sarea;
		int maxId;

		l[0] = e_l[fh_idx[f1][0] >> 1];
		l[1] = e_l[fh_idx[f1][1] >> 1];
		l[2] = e_l[fh_idx[f1][2] >> 1];

		maxId = 0;
		if (fh_theta[f1][maxId] < fh_theta[f1][1]) maxId = 1;
		if (fh_theta[f1][maxId] < fh_theta[f1][2]) maxId = 2;

		if (fh_theta[f1][maxId] < M_PI / 2)
		{
			varea[0] = l[0] * l[0] * fh_cot[f1][0] / 8;
			varea[1] = l[1] * l[1] * fh_cot[f1][1] / 8;
			varea[2] = l[2] * l[2] * fh_cot[f1][2] / 8;

			v_s[mesh.to_vertex_handle(HH(fh_idx[f1][1])).idx()] -= varea[1] - varea[2];
			v_s[mesh.to_vertex_handle(HH(fh_idx[f1][2])).idx()] -= varea[2] - varea[0];
			v_s[mesh.to_vertex_handle(HH(fh_idx[f1][0])).idx()] -= varea[0] - varea[1];
		}
		else
		{
			sarea = (l[0] * l[1] * sin(fh_theta[f1][2]) + l[1] * l[2] * sin(fh_theta[f1][0])
				+ l[2] * l[0] * sin(fh_theta[f1][1])) / 12;

			v_s[mesh.to_vertex_handle(HH(fh_idx[f1][1])).idx()] -= sarea;
			v_s[mesh.to_vertex_handle(HH(fh_idx[f1][2])).idx()] -= sarea;
			v_s[mesh.to_vertex_handle(HH(fh_idx[f1][0])).idx()] -= sarea;

			v_s[mesh.to_vertex_handle(HH(fh_idx[f1][maxId])).idx()] += sarea / 2;
			v_s[mesh.to_vertex_handle(HH(fh_idx[f1][maxId] ^ 1)).idx()] += sarea / 2;
		}

		l[0] = e_l[fh_idx[f2][0] >> 1];
		l[1] = e_l[fh_idx[f2][1] >> 1];
		l[2] = e_l[fh_idx[f2][2] >> 1];

		maxId = 0;
		if (fh_theta[f2][maxId] < fh_theta[f2][1]) maxId = 1;
		if (fh_theta[f2][maxId] < fh_theta[f2][2]) maxId = 2;

		if (fh_theta[f2][maxId] < M_PI / 2)
		{
			varea[0] = l[0] * l[0] * fh_cot[f2][0] / 8;
			varea[1] = l[1] * l[1] * fh_cot[f2][1] / 8;
			varea[2] = l[2] * l[2] * fh_cot[f2][2] / 8;

			v_s[mesh.to_vertex_handle(HH(fh_idx[f2][1])).idx()] -= varea[1] - varea[2];
			v_s[mesh.to_vertex_handle(HH(fh_idx[f2][2])).idx()] -= varea[2] - varea[0];
			v_s[mesh.to_vertex_handle(HH(fh_idx[f2][0])).idx()] -= varea[0] - varea[1];
		}
		else
		{
			sarea = (l[0] * l[1] * sin(fh_theta[f2][2]) + l[1] * l[2] * sin(fh_theta[f2][0])
				+ l[2] * l[0] * sin(fh_theta[f2][1])) / 12;

			v_s[mesh.to_vertex_handle(HH(fh_idx[f2][1])).idx()] -= sarea;
			v_s[mesh.to_vertex_handle(HH(fh_idx[f2][2])).idx()] -= sarea;
			v_s[mesh.to_vertex_handle(HH(fh_idx[f2][0])).idx()] -= sarea;

			v_s[mesh.to_vertex_handle(HH(fh_idx[f2][maxId])).idx()] += sarea / 2;
			v_s[mesh.to_vertex_handle(HH(fh_idx[f2][maxId] ^ 1)).idx()] += sarea / 2;
		}*/

		mesh.flip(EH(idx));

		e_Bl[idx] = lab0;
		e_l[idx] = lab;

		hv_idx[hcd.idx()] = va.idx();
		hv_idx[hdc.idx()] = vb.idx();

		fh_idx[f1][0] = hcd.idx();
		fh_idx[f1][1] = hac.idx();
		fh_idx[f1][2] = hcb.idx();
		fh_idx[f2][0] = hdc.idx();
		fh_idx[f2][1] = hbd.idx();
		fh_idx[f2][2] = hda.idx();

		fh_theta[f1][0] = acb;
		fh_theta[f1][1] = cba;
		fh_theta[f1][2] = bac;
		fh_theta[f2][0] = bda;
		fh_theta[f2][1] = dab;
		fh_theta[f2][2] = abd;

		cot(fh_theta[f1], fh_cot[f1]);
		cot(fh_theta[f2], fh_cot[f2]);

		e_theta[fh_idx[f1][0] >> 1] += fh_theta[f1][0];
		e_theta[fh_idx[f1][1] >> 1] += fh_theta[f1][1];
		e_theta[fh_idx[f1][2] >> 1] += fh_theta[f1][2];
		e_theta[fh_idx[f2][0] >> 1] += fh_theta[f2][0];
		e_theta[fh_idx[f2][1] >> 1] += fh_theta[f2][1];
		e_theta[fh_idx[f2][2] >> 1] += fh_theta[f2][2];

		/*v_g[mesh.to_vertex_handle(HH(fh_idx[f1][0])).idx()] -= fh_theta[f1][2];
		v_g[mesh.to_vertex_handle(HH(fh_idx[f1][1])).idx()] -= fh_theta[f1][0];
		v_g[mesh.to_vertex_handle(HH(fh_idx[f1][2])).idx()] -= fh_theta[f1][1];
		v_g[mesh.to_vertex_handle(HH(fh_idx[f2][0])).idx()] -= fh_theta[f2][2];
		v_g[mesh.to_vertex_handle(HH(fh_idx[f2][1])).idx()] -= fh_theta[f2][0];
		v_g[mesh.to_vertex_handle(HH(fh_idx[f2][2])).idx()] -= fh_theta[f2][1];

		l[0] = e_l[fh_idx[f1][0] >> 1];
		l[1] = e_l[fh_idx[f1][1] >> 1];
		l[2] = e_l[fh_idx[f1][2] >> 1];

		maxId = 0;
		if (fh_theta[f1][maxId] < fh_theta[f1][1]) maxId = 1;
		if (fh_theta[f1][maxId] < fh_theta[f1][2]) maxId = 2;

		if (fh_theta[f1][maxId] < M_PI / 2)
		{
			varea[0] = l[0] * l[0] * fh_cot[f1][0] / 8;
			varea[1] = l[1] * l[1] * fh_cot[f1][1] / 8;
			varea[2] = l[2] * l[2] * fh_cot[f1][2] / 8;

			v_s[mesh.to_vertex_handle(HH(fh_idx[f1][1])).idx()] += varea[1] + varea[2];
			v_s[mesh.to_vertex_handle(HH(fh_idx[f1][2])).idx()] += varea[2] + varea[0];
			v_s[mesh.to_vertex_handle(HH(fh_idx[f1][0])).idx()] += varea[0] + varea[1];
		}
		else
		{
			sarea = (l[0] * l[1] * sin(fh_theta[f1][2]) + l[1] * l[2] * sin(fh_theta[f1][0])
				+ l[2] * l[0] * sin(fh_theta[f1][1])) / 12;

			v_s[mesh.to_vertex_handle(HH(fh_idx[f1][1])).idx()] += sarea;
			v_s[mesh.to_vertex_handle(HH(fh_idx[f1][2])).idx()] += sarea;
			v_s[mesh.to_vertex_handle(HH(fh_idx[f1][0])).idx()] += sarea;

			v_s[mesh.to_vertex_handle(HH(fh_idx[f1][maxId])).idx()] -= sarea / 2;
			v_s[mesh.to_vertex_handle(HH(fh_idx[f1][maxId] ^ 1)).idx()] -= sarea / 2;
		}

		l[0] = e_l[fh_idx[f2][0] >> 1];
		l[1] = e_l[fh_idx[f2][1] >> 1];
		l[2] = e_l[fh_idx[f2][2] >> 1];

		maxId = 0;
		if (fh_theta[f2][maxId] < fh_theta[f2][1]) maxId = 1;
		if (fh_theta[f2][maxId] < fh_theta[f2][2]) maxId = 2;

		if (fh_theta[f2][maxId] < M_PI / 2)
		{
			varea[0] = l[0] * l[0] * fh_cot[f2][0] / 8;
			varea[1] = l[1] * l[1] * fh_cot[f2][1] / 8;
			varea[2] = l[2] * l[2] * fh_cot[f2][2] / 8;

			v_s[mesh.to_vertex_handle(HH(fh_idx[f2][1])).idx()] += varea[1] + varea[2];
			v_s[mesh.to_vertex_handle(HH(fh_idx[f2][2])).idx()] += varea[2] + varea[0];
			v_s[mesh.to_vertex_handle(HH(fh_idx[f2][0])).idx()] += varea[0] + varea[1];
		}
		else
		{
			sarea = (l[0] * l[1] * sin(fh_theta[f2][2]) + l[1] * l[2] * sin(fh_theta[f2][0])
				+ l[2] * l[0] * sin(fh_theta[f2][1])) / 12;

			v_s[mesh.to_vertex_handle(HH(fh_idx[f2][1])).idx()] += sarea;
			v_s[mesh.to_vertex_handle(HH(fh_idx[f2][2])).idx()] += sarea;
			v_s[mesh.to_vertex_handle(HH(fh_idx[f2][0])).idx()] += sarea;

			v_s[mesh.to_vertex_handle(HH(fh_idx[f2][maxId])).idx()] -= sarea / 2;
			v_s[mesh.to_vertex_handle(HH(fh_idx[f2][maxId] ^ 1)).idx()] -= sarea / 2;
		}*/

		idx = hbd.idx() >> 1;
		if (!e_b[idx] && e_theta[idx] >= THETAFLIP) { add.push_back(idx); e_b[idx] = true; }
		idx = hda.idx() >> 1;
		if (!e_b[idx] && e_theta[idx] >= THETAFLIP) { add.push_back(idx); e_b[idx] = true; }
		idx = hac.idx() >> 1;
		if (!e_b[idx] && e_theta[idx] >= THETAFLIP) { add.push_back(idx); e_b[idx] = true; }
		idx = hcb.idx() >> 1;
		if (!e_b[idx] && e_theta[idx] >= THETAFLIP) { add.push_back(idx); e_b[idx] = true; }

		return true;
	}

protected:
	void initDomainMesh()
	{
		auto v_p0 = mesh.points();
		for (int i = 0; i < e_l.size(); i++)
		{
			e_Bl[i] = (v_p0[mesh.to_vertex_handle(HH(i << 1)).idx()] -
				v_p0[mesh.to_vertex_handle(HH((i << 1) ^ 1)).idx()]).norm();
		}

		HH tmp;
		for (int i = 0; i < fh_idx.size(); i++)
		{
			tmp = mesh.halfedge_handle(FH(i));
			fh_idx[i][0] = tmp.idx();

			tmp = mesh.next_halfedge_handle(tmp);
			fh_idx[i][1] = tmp.idx();

			tmp = mesh.next_halfedge_handle(tmp);
			fh_idx[i][2] = tmp.idx();
		}

		for (int i = 0; i < hv_idx.size(); i++)
		{
			hv_idx[i] = mesh.to_vertex_handle(HH(i)).idx();
		}
	}

	void edgeFlip()
	{
		std::vector<int> stack;
		stack.reserve(e_b.size());
		updateEdgeTheta();

		for (int i = 0; i < e_b.size(); i++)
		{
			e_b[i] = false;
			if (e_theta[i] < THETAFLIP) continue;
			stack.push_back(i);
			e_b[i] = true;
		}

		OpenMesh::HalfedgeHandle he_t;
		int t;
		while (!stack.empty())
		{
			t = stack.back();
			stack.pop_back();
			e_b[t] = false;

			if (!mesh.is_flip_ok(EH(t))) continue;
			if (mesh.status(EH(t)).feature()) continue;
			attemptFlip(t, stack);
		}

	}

	bool isTriangulable()
	{
		double l[3];
		for (int i = 0; i < fh_theta.size(); i++)
		{
			l[0] = e_l[fh_idx[i][0] >> 1];
			l[1] = e_l[fh_idx[i][1] >> 1];
			l[2] = e_l[fh_idx[i][2] >> 1];

			if (l[0] + l[1] > l[2] && l[0] + l[2] > l[1] && l[1] + l[2] > l[0]) continue;
			return false;
		}

		return true;
	}

	void updateEdgeLength()
	{
		double u0, u1;
		for (int i = 0; i < e_l.size(); i++)
		{
			u0 = v_u[hv_idx[i << 1]];
			u1 = v_u[hv_idx[(i << 1) ^ 1]];
			e_l[i] = exp(u0 + u1)*e_Bl[i];
		}
	}

	void updateHalfEdgeTheta()
	{
		double l[3];
		for (int i = 0; i < fh_theta.size(); i++)
		{
			l[0] = e_l[fh_idx[i][0] >> 1];
			l[1] = e_l[fh_idx[i][1] >> 1];
			l[2] = e_l[fh_idx[i][2] >> 1];

			calcTheta(l[0], l[1], l[2], fh_theta[i][0]);
			calcTheta(l[1], l[0], l[2], fh_theta[i][1]);
			fh_theta[i][2] = M_PI - fh_theta[i][0] - fh_theta[i][1];

			cot(fh_theta[i], fh_cot[i]);
		}
	}

	void updateEdgeTheta()
	{
		for (int i = 0; i < e_theta.size(); i++)
		{
			e_theta[i] = 0;
		}

		for (int i = 0; i < fh_idx.size(); i++)
		{
			e_theta[fh_idx[i][0] >> 1] += fh_theta[i][0];
			e_theta[fh_idx[i][1] >> 1] += fh_theta[i][1];
			e_theta[fh_idx[i][2] >> 1] += fh_theta[i][2];
		}
	}

	void updateGaussCurvature()
	{
		for (int i = 0; i < in_idx.size(); i++)
		{
			v_g[in_idx[i]] = 2 * M_PI;
		}
		for (int i = 0; i < bd_idx.size(); i++)
		{
			v_g[bd_idx[i]] = M_PI;
		}

		for (int i = 0; i < fh_idx.size(); i++)
		{
			v_g[hv_idx[fh_idx[i][0]]] -= fh_theta[i][2];
			v_g[hv_idx[fh_idx[i][1]]] -= fh_theta[i][0];
			v_g[hv_idx[fh_idx[i][2]]] -= fh_theta[i][1];
		}
	}

	void updateCotNormal()
	{
		for (int i = 0; i < v_cn.size(); i++)
		{
			v_cn[i].fill(0);
		}

		for (int i = 0; i < fh_idx.size(); i++)
		{
			sub(v_p[hv_idx[fh_idx[i][0]]], v_p[hv_idx[fh_idx[i][2]]], fh_d[i][0]);
			sub(v_p[hv_idx[fh_idx[i][1]]], v_p[hv_idx[fh_idx[i][0]]], fh_d[i][1]);
			sub(v_p[hv_idx[fh_idx[i][2]]], v_p[hv_idx[fh_idx[i][1]]], fh_d[i][2]);

			fh_d[i][0] *= fh_cot[i][0];
			fh_d[i][1] *= fh_cot[i][1];
			fh_d[i][2] *= fh_cot[i][2];

			subOn(fh_d[i][0], fh_d[i][1], v_cn[hv_idx[fh_idx[i][0]]]);
			subOn(fh_d[i][1], fh_d[i][2], v_cn[hv_idx[fh_idx[i][1]]]);
			subOn(fh_d[i][2], fh_d[i][0], v_cn[hv_idx[fh_idx[i][2]]]);
		}
	}

	void updateNormal()
	{
		for (int i = 0; i < v_cn.size(); i++)
		{
			v_n[i].fill(0);
		}

		for (int i = 0; i < fh_idx.size(); i++)
		{
			sub(v_p[hv_idx[fh_idx[i][0]]], v_p[hv_idx[fh_idx[i][2]]], fh_d[i][0]);
			sub(v_p[hv_idx[fh_idx[i][1]]], v_p[hv_idx[fh_idx[i][0]]], fh_d[i][1]);
			sub(v_p[hv_idx[fh_idx[i][2]]], v_p[hv_idx[fh_idx[i][1]]], fh_d[i][2]);

			v_n[hv_idx[fh_idx[i][0]]] += (fh_d[i][0] % fh_d[i][1]);
			v_n[hv_idx[fh_idx[i][1]]] += (fh_d[i][1] % fh_d[i][2]);
			v_n[hv_idx[fh_idx[i][2]]] += (fh_d[i][2] % fh_d[i][0]);
		}

		for (int i = 0; i < v_cn.size(); i++)
		{
			v_n[i].normalize_cond();
		}
	}

	void updateMeanCurvature()
	{
		updateNormal();
		updateCotNormal();
		for (int i = 0; i < v_h.size(); i++)
		{
			v_h[i] = v_cn[i].norm();
			if ((v_cn[i] | v_n[i]) < 0)
				v_h[i] = -v_h[i];
		}
	}

	void updateMixedArea()
	{
		for (int i = 0; i < v_s.size(); i++)
		{
			v_s[i] = 0;
		}

		Vec3d l, varea;
		double sarea;
		int maxId;

		for (int i = 0; i < fh_idx.size(); i++)
		{
			l[0] = e_l[fh_idx[i][0] >> 1];
			l[1] = e_l[fh_idx[i][1] >> 1];
			l[2] = e_l[fh_idx[i][2] >> 1];

			maxId = 0;
			if (fh_theta[i][maxId] < fh_theta[i][1]) maxId = 1;
			if (fh_theta[i][maxId] < fh_theta[i][2]) maxId = 2;

			if (fh_theta[i][maxId] < M_PI / 2)
			{
				varea[0] = l[0] * l[0] * fh_cot[i][0] / 8;
				varea[1] = l[1] * l[1] * fh_cot[i][1] / 8;
				varea[2] = l[2] * l[2] * fh_cot[i][2] / 8;

				v_s[hv_idx[fh_idx[i][0]]] += varea[0] + varea[1];
				v_s[hv_idx[fh_idx[i][1]]] += varea[1] + varea[2];
				v_s[hv_idx[fh_idx[i][2]]] += varea[2] + varea[0];
			}
			else
			{
				sarea = (l[0] * l[1] * sin(fh_theta[i][2]) + l[1] * l[2] * sin(fh_theta[i][0])
					+ l[2] * l[0] * sin(fh_theta[i][1])) / 24;

				v_s[hv_idx[fh_idx[i][maxId]]] += sarea;
				v_s[hv_idx[fh_idx[i][(maxId - 1) % 3]]] += sarea;
				v_s[hv_idx[fh_idx[i][(maxId + 1) % 3]]] += sarea * 2;
			}
		}
	}
};

#endif // !DOMAINMESH_H