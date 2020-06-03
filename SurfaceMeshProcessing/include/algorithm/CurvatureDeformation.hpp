#ifndef CURVATUREDEFORMATION_H
#define CURVATUREDEFORMATION_H

#ifdef _DEBUG
#define verify(f) assert(f)
#else
#define verify(f) ((void)(f))
#endif

#include <ctime>
#include <fstream>
#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/CholmodSupport>
#include <algorithm/DomainMesh.hpp>
#include <algorithm/AndersonAcceleration.hpp>

class CurvatureDeformation : public DomainMesh
{
protected:
	// parameters in quasi Newton opt
	double shrink = 0.75; // the ratio to shrink step length
	int quas_cap = 300; // max iter number
	int up_cap = 1; // upward terminate condition
	bool increasing = false;
	double curv_err = 1e-6; // the terminate condition
	double progClip = 1e-5;

	// parameters in local shape deformation
	double lambdaH = 1; // mean curvature weight
	int deformH_cap = 60; // With mean curvature max iter number
	int deformG_cap = 20; // No mean curvature max iter number
	int andersonM = 5; // anderson acceleration basis num

	// parameters in ASAP assemble
	int ASAP_cap = 5; // max iter number

protected:
	Vectord v_Dg;
	Eigen::MatrixXd v_pos;
	Vectori var_mark;
	Vectori var_idx;
	Vectori ach_idx;

	std::vector<Eigen::Triplet<double>> dual_lap_trip; // triplets for dual laplacian matrix
	std::vector<Eigen::Triplet<double>> var_lap_trip; // triplets for var pts laplacian matrix
	Eigen::SparseMatrix<double, Eigen::ColMajor> mat_L; // local part of laplace Matrix
	Eigen::SparseMatrix<double, Eigen::ColMajor> mat_C; // local cols part of laplace Matrix
	Eigen::CholmodSimplicialLDLT<Eigen::SparseMatrix<double, Eigen::ColMajor>> cllt;

protected:
	virtual void calabiAnchors() {}
	virtual void deformAnchors() {}
	virtual bool checkConditions() { return false; }
	virtual void setIntialPosition(const Vectori &achs, const Vector3d &dpos) {}
	virtual void updatePosition() {}
	virtual void projectPosition() {}
	virtual void setTargetGaussCurvature() {}
	virtual void setTargetMeanCurvature() {}
	virtual void updateTargetGaussCurvature() {}
	virtual void updateTargetMeanCurvature() {}

public:
	CurvatureDeformation(Mesh &input, const std::vector<std::vector<Mesh::HalfedgeHandle>> &contours) : DomainMesh(input, contours)
	{
		var_mark.resize(v_p.size(), 0);
		anchorPos({}, {});

		/*updateMesh();
		updateMeanCurvature();

		std::ofstream gausscurve;
		gausscurve.open("./meancurve.txt");

		for (int i = 0; i < v_g.size(); i++)
			gausscurve << fabs(v_h[i]) / v_s[i] << std::endl;

		gausscurve.close();*/

	}

	void anchorPos(const Vectori &achs, const Vector3d &dpos)
	{
		var_idx.clear();
		var_idx.reserve(v_u.size());
		ach_idx = achs;

		for (int i = 0; i < var_mark.size(); i++)
			var_mark[i] = 0;
		
		for (int i = 0; i < ach_idx.size(); i++)
			var_mark[ach_idx[i]] = -1;

		for (int i = 0; i < var_mark.size(); i++)
		{
			if (var_mark[i] < 0) continue;
			var_mark[i] = var_idx.size();
			var_idx.push_back(i);
		}

		v_Dg.resize(var_idx.size(), 0);
		v_pos.resize(var_idx.size(), 3);
		setIntialPosition(achs, dpos);
	}

	double run()
	{
		if (!checkConditions()) return 0;

		std::clock_t start, middle, end;
		start = clock();

		if (increasing)
		{
			mesh = back;
			initDomainMesh();
			for (int i = 0; i < v_u.size(); i++)
			{
				v_u[i] = 0;
			}
			increasing = false;
		}

		updateMesh();

		std::cout << "\nCalabi flow ...\n";
		calabiAnchors();
		calabiFlow();
		//mixedQuasi(0);

		middle = clock();
		std::cout << "Calabi flow time : " << (double)(middle - start) << "ms\n";

		updatePosition();

		std::cout << "\nDeformation ...\n";
		deformAnchors();
		deformation();

		projectPosition();

		end = clock();
		std::cout << "Deformation time : " << (double)(end - middle) << "ms\n";
		std::cout << "Total time : " << (double)(end - start) << "ms\n";
		return end - start;
	}

	const Vector3d& deformedPos() const { return v_p; }

	virtual void setCurvatureParameters(int slideValue) {}

private:
	//-------------------------------------------------------
	// ***** ------------- Calabi Flow ------------- ***** //
	double gaussDeviation()
	{
		double l2norm = 0;
		for (int i=0; i<v_Dg.size(); i++)
		{
			v_Dg[i] = v_Tg[var_idx[i]] - v_g[var_idx[i]];
			l2norm += pow(v_Dg[i], 2);
		}
		return l2norm;
	}
	
	void gaussIncrement(double l2norm)
	{
		double scale = 1;
		double mean = progClip * v_Dg.size();
		//std::cout << mean << ' ' << l2norm << std::endl;
		if (l2norm > mean) scale = mean / l2norm;

		for (int i=0; i<v_Dg.size(); i++)
		{
			v_Dg[i] = v_Dg[i] * scale;
		}
	}

	void updateMesh()
	{
		updateEdgeLength();
		updateHalfEdgeTheta();
		//std::cout << "Edge Flipping:";
		edgeFlip();
		//std::cout << "......    ";
		updateMixedArea();
		updateGaussCurvature();
	}

	void assembleLaplaceTriplets()
	{
		dual_lap_trip.clear();
		dual_lap_trip.reserve(9 * fh_idx.size());

		int vid[3];
		for (int i = 0; i < fh_idx.size(); i++)
		{
			vid[0] = hv_idx[fh_idx[i][1]];
			vid[1] = hv_idx[fh_idx[i][2]];
			vid[2] = hv_idx[fh_idx[i][0]];

			dual_lap_trip.push_back(Eigen::Triplet<double>(vid[0], vid[0], fh_cot[i][1] + fh_cot[i][2]));
			dual_lap_trip.push_back(Eigen::Triplet<double>(vid[1], vid[1], fh_cot[i][0] + fh_cot[i][2]));
			dual_lap_trip.push_back(Eigen::Triplet<double>(vid[2], vid[2], fh_cot[i][0] + fh_cot[i][1]));

			dual_lap_trip.push_back(Eigen::Triplet<double>(vid[0], vid[1], -fh_cot[i][2]));
			dual_lap_trip.push_back(Eigen::Triplet<double>(vid[1], vid[0], -fh_cot[i][2]));
			dual_lap_trip.push_back(Eigen::Triplet<double>(vid[1], vid[2], -fh_cot[i][0]));
			dual_lap_trip.push_back(Eigen::Triplet<double>(vid[2], vid[1], -fh_cot[i][0]));
			dual_lap_trip.push_back(Eigen::Triplet<double>(vid[2], vid[0], -fh_cot[i][1]));
			dual_lap_trip.push_back(Eigen::Triplet<double>(vid[0], vid[2], -fh_cot[i][1]));
		}
	}

	void assembleLaplaceMatrix()
	{
		var_lap_trip.clear();
		var_lap_trip.reserve(9 * fh_idx.size());

		int var_row, var_col;
		for (int i = 0; i < dual_lap_trip.size(); i++)
		{
			var_row = var_mark[dual_lap_trip[i].row()];
			var_col = var_mark[dual_lap_trip[i].col()];
			if (var_row < 0 || var_col < 0) continue;
			var_lap_trip.push_back(Eigen::Triplet<double>(var_row, var_col, dual_lap_trip[i].value()));
		}

		mat_L.resize(var_idx.size(), var_idx.size());
		mat_L.setFromTriplets(var_lap_trip.begin(), var_lap_trip.end());
	}

	double calcQuasStep(const double *du)
	{
		double quasStep = 1.0f;
		bool triangulable = false;

		for (int i = 0; i < var_idx.size(); i++)
		{
			v_Bu[var_idx[i]] = v_u[var_idx[i]];
		}

		int count = 0;
		do
		{
			for (int i = 0; i < var_idx.size(); i++)
			{
				v_u[var_idx[i]] = v_Bu[var_idx[i]] + quasStep * du[i];
			}
			updateEdgeLength();
			triangulable = isTriangulable();
			if (!triangulable) quasStep *= shrink;
			count++;
		} while (!triangulable && count < 100);

		for (int i = 0; i < var_idx.size(); i++)
		{
			v_u[var_idx[i]] = v_Bu[var_idx[i]];
		}
		if (count >= 10) return 0;
		return quasStep;
	}

	void calabiFlow()
	{
		int count, upcount;
		double last_err, curr_err;

		double source_areaSum = 0;
		for (int i = 0; i < v_s.size(); i++)
		{
			source_areaSum += v_s[i];
		}

		Eigen::Map<Eigen::VectorXd> du(v_Dg.data(), v_Dg.size());
		setTargetGaussCurvature();
		curr_err = gaussDeviation();

		last_err = curr_err + 1;
		std::cout << "Initial Curvatre Error: " << curr_err << std::endl;
		std::cout << "----------------------------------------------\n";

		clock_t start, end;
		double tick0 = 0, tick1 = 0, tick2 = 0;

		count = 0;
		upcount = 0;
		std::cout << "Quasi Newton phase:\n";
		while (curr_err >= curv_err && count < quas_cap)
		{
			if (curr_err >= last_err) upcount++;
			else upcount = 0;
			increasing = (upcount >= up_cap);
			//if (increasing) break;
			last_err = curr_err;

			start = clock();
			assembleLaplaceTriplets();
			assembleLaplaceMatrix();
			end = clock();
			tick0 += double(end - start);

			start = clock();
			cllt.compute(mat_L);
			end = clock();
			tick1 += double(end - start);

			start = clock();
			gaussIncrement(curr_err);
			du = cllt.solve(du);
			du = du * calcQuasStep(du.data());
			end = clock();
			tick2 += double(end - start);

			for (int i = 0; i < du.size(); i++)
			{
				v_u[var_idx[i]] += du[i];
			}
			
			updateMesh();
			updateTargetGaussCurvature();

			curr_err = gaussDeviation();
			std::cout << "the " << count << "th iter: " << curr_err << std::endl;
			count++;
		}

		//std::cout << "Assemble time : " << tick0 << std::endl;
		//std::cout << "Compute time : " << tick1 << std::endl;
		//std::cout << "Solve time : " << tick2 << std::endl;

		double areaSum = 0;
		for (int i = 0; i < v_s.size(); i++)
		{
			areaSum += v_s[i];
		}

		double minusU = (log(areaSum) - log(source_areaSum)) / 4;
		for (int i = 0; i < v_u.size(); i++)
		{
			v_u[i] -= minusU;
		}

		updateEdgeLength();
		updateMixedArea();
	}

	//-------------------------------------------------------
	// ***** ------------- Comparasion ------------- ***** //

	void linesearch(double curr_err, double &step, Eigen::Map<Eigen::VectorXd> &opgrad)
	{
		double norm = opgrad.norm();
		if (norm == 0) return;
		opgrad /= norm;

		for (int i = 0; i < var_idx.size(); i++)
		{
			v_Bu[var_idx[i]] = v_u[var_idx[i]];
		}

		double new_err, eval;
		bool in_circle;
		int count = 0;
		do
		{
			in_circle = false;
			for (int i = 0; i < var_idx.size(); i++)
			{
				v_u[var_idx[i]] = v_Bu[var_idx[i]] + step * opgrad[i];
			}
			updateEdgeLength();
			if (!isTriangulable())
			{
				step *= 0.5;
				in_circle = true;
			}
			else
			{
				updateHalfEdgeTheta();
				edgeFlip();
				updateMixedArea();
				updateGaussCurvature();
				updateTargetGaussCurvature();
				if (gaussDeviation() >= curr_err - 0.5 * step * norm)
				{
					step *= 0.5;
					in_circle = true;
				}
			}
			count++;
		} while (in_circle && count < 1000);


		for (int i = 0; i < var_idx.size(); i++)
		{
			v_u[var_idx[i]] = v_Bu[var_idx[i]];
		}
	}

	void gradDesc()
	{
		int count;
		double curr_err;

		Eigen::Map<Eigen::VectorXd> du(v_Dg.data(), v_Dg.size());
		setTargetGaussCurvature();
		curr_err = gaussDeviation();
		//std::cout << "Initial Curvatre Error: " << curr_err << std::endl;
		//std::cout << "----------------------------------------------\n";

		clock_t start, end;
		double tick = 0;
		std::vector<double> ticks;
		std::vector<double> errors;
		ticks.push_back(0);
		errors.push_back(curr_err);

		count = 0;
		double step = 1e-2;
		start = clock();
		std::cout << "Grad Desc phase:\n";
		while (count < 100)
		{
			assembleLaplaceTriplets();
			assembleLaplaceMatrix();
			
			du = mat_L * du;
			du.normalize();
			linesearch(curr_err, step, du);

			for (int i = 0; i < du.size(); i++)
			{
				v_u[var_idx[i]] += step *du[i];
			}

			updateEdgeLength();
			updateHalfEdgeTheta();
			edgeFlip();
			updateMixedArea();
			updateGaussCurvature();
			updateTargetGaussCurvature();
			curr_err = gaussDeviation();
			//std::cout << "the " << count << "th iter: " << curr_err << std::endl;
			count++;

			end = clock();
			tick = double(end - start);
			if (count % 3 != 0) continue;
			ticks.push_back(tick);
			errors.push_back(curr_err);
		}

		std::ofstream fout("./data/GradDesc.txt");
		if (fout)
		{
			for (int i = 0; i < errors.size(); i++)
			{
				fout << errors[i] << '\t' << ticks[i] << std::endl;
			}
		}
		fout.close();
	}

	void mixedQuasi(int n)
	{
		int count;
		double curr_err;

		double source_areaSum = 0;
		for (int i = 0; i < v_s.size(); i++)
		{
			source_areaSum += v_s[i];
		}

		Eigen::Map<Eigen::VectorXd> du(v_Dg.data(), v_Dg.size());
		setTargetGaussCurvature();
		curr_err = gaussDeviation();
		std::cout << "Initial Curvatre Error: " << curr_err << std::endl;
		std::cout << "----------------------------------------------\n";

		clock_t start, end;
		double tick = 0;
		std::vector<double> ticks;
		std::vector<double> errors;
		ticks.push_back(0);
		errors.push_back(curr_err);

		count = 0;
		double step = 1e-2;
		start = clock();
		std::cout << "Mixed Quasi phase:\n";
		while (count < n)
		{
			assembleLaplaceTriplets();
			assembleLaplaceMatrix();

			du = mat_L * du;
			du.normalize();
			linesearch(curr_err, step, du);

			for (int i = 0; i < du.size(); i++)
			{
				v_u[var_idx[i]] += step *du[i];
			}

			updateEdgeLength();
			updateHalfEdgeTheta();
			edgeFlip();
			updateMixedArea();
			updateGaussCurvature();
			updateTargetGaussCurvature();

			curr_err = gaussDeviation();
			std::cout << "the " << count << "th iter: " << curr_err << std::endl;
			count++;

			end = clock();
			tick = double(end - start);
			if (count % 3 != 0) continue;
			ticks.push_back(tick);
			errors.push_back(curr_err);
		}

		count = 0;
		while (count < quas_cap)
		{
			assembleLaplaceTriplets();
			assembleLaplaceMatrix();

			cllt.compute(mat_L);
			du = cllt.solve(du);
			du = du * 0.5 * calcQuasStep(du.data());

			for (int i = 0; i < du.size(); i++)
			{
				v_u[var_idx[i]] += du[i];
			}

			updateEdgeLength();
			updateHalfEdgeTheta();
			edgeFlip();
			updateMixedArea();
			updateGaussCurvature();
			updateTargetGaussCurvature();

			curr_err = gaussDeviation();
			std::cout << "the " << count << "th iter: " << curr_err << std::endl;
			count++;

			end = clock();
			tick = double(end - start);
			ticks.push_back(tick);
			errors.push_back(curr_err);
		}

		double areaSum = 0;
		for (int i = 0; i < v_s.size(); i++)
		{
			areaSum += v_s[i];
		}

		double minusU = (log(areaSum) - log(source_areaSum)) / 4;
		for (int i = 0; i < v_u.size(); i++)
		{
			v_u[i] -= minusU;
		}

		updateEdgeLength();
		updateMixedArea();

		/*std::ofstream fout("./data/MixedQuasi" + std::to_string(n) + ".txt");
		if (fout)
		{
			for (int i = 0; i < errors.size(); i++)
			{
				fout << errors[i] << '\t' << ticks[i] << std::endl;
			}
		}
		fout.close();*/
	}

	void progQuasi()
	{
		int count;
		double curr_err;

		Eigen::Map<Eigen::VectorXd> du(v_Dg.data(), v_Dg.size());
		setTargetGaussCurvature();
		curr_err = gaussDeviation();
		std::cout << "Initial Curvatre Error: " << curr_err << std::endl;
		std::cout << "----------------------------------------------\n";

		clock_t start, end;
		double tick = 0;
		std::vector<double> ticks;
		std::vector<double> errors;
		ticks.push_back(0);
		errors.push_back(curr_err);

		count = 0;
		start = clock();
		std::cout << "Prog Quasi phase:\n";
		while (count < 40)
		{
			assembleLaplaceTriplets();
			assembleLaplaceMatrix();

			cllt.compute(mat_L);
			gaussIncrement(curr_err);
			du = cllt.solve(du);
			du = du * calcQuasStep(du.data());

			for (int i = 0; i < du.size(); i++)
			{
				v_u[var_idx[i]] += du[i];
			}

			updateEdgeLength();
			updateHalfEdgeTheta();
			edgeFlip();
			updateMixedArea();
			updateGaussCurvature();
			updateTargetGaussCurvature();

			curr_err = gaussDeviation();
			std::cout << "the " << count << "th iter: " << curr_err << std::endl;
			count++;

			end = clock();
			tick = double(end - start);
			ticks.push_back(tick);
			errors.push_back(curr_err);
		}

		std::ofstream fout("./data/ProgQuasi.txt");
		if (fout)
		{
			for (int i = 0; i < errors.size(); i++)
			{
				fout << errors[i] << '\t' << ticks[i] << std::endl;
			}
		}
		fout.close();
	}

public:
	void runCompare()
	{
		if (ach_idx.size() == 0)
		{
			std::cout << "Need at least one fixed point\n";
			return;
		}


		mesh = back;
		initDomainMesh();
		for (int i = 0; i < v_u.size(); i++)
		{
			v_u[i] = 0;
		}
		updateMesh();

		std::cout << "\nGrad Desc ...\n";
		gradDesc();

		mesh = back;
		initDomainMesh();
		for (int i = 0; i < v_u.size(); i++)
		{
			v_u[i] = 0;
		}
		updateMesh();

		std::cout << "\Mixed Quasi ...\n";
		mixedQuasi(0);

		mesh = back;
		initDomainMesh();
		for (int i = 0; i < v_u.size(); i++)
		{
			v_u[i] = 0;
		}
		updateMesh();

		std::cout << "\Mixed Quasi ...\n";
		mixedQuasi(15);

		mesh = back;
		initDomainMesh();
		for (int i = 0; i < v_u.size(); i++)
		{
			v_u[i] = 0;
		}
		updateMesh();

		std::cout << "\Mixed Quasi ...\n";
		mixedQuasi(30);

		mesh = back;
		initDomainMesh();
		for (int i = 0; i < v_u.size(); i++)
		{
			v_u[i] = 0;
		}
		updateMesh();

		std::cout << "\Prog Quasi ...\n";
		progQuasi();
	}

private:
	//-------------------------------------------------------
	// ***** ------------- Deformation ------------- ***** //
	void assembleLaplaceColMatrix()
	{
		var_lap_trip.clear();
		var_lap_trip.reserve(9 * fh_idx.size());

		int var_col;
		for (int i = 0; i < dual_lap_trip.size(); i++)
		{
			var_col = var_mark[dual_lap_trip[i].col()];
			if (var_col < 0 || bd_mark[dual_lap_trip[i].row()]) continue;
			var_lap_trip.push_back(Eigen::Triplet<double>(dual_lap_trip[i].row(), var_col, dual_lap_trip[i].value()));
		}

		mat_C.resize(var_mark.size(), var_idx.size());
		mat_C.setFromTriplets(var_lap_trip.begin(), var_lap_trip.end());
	}

	inline double massSpring(Vec3d &x, double l)
	{
		double len = x.normalize_cond();
		x *= (len + l) / 2;
		return pow(len - l, 2.0);
	}

	double laplaceSpringNormal(double *Ln)
	{
		double sqr_err = 0.0;
		updateCotNormal();
		updateNormal();

		for (int i = 0; i < bd_idx.size(); i++)
		{
			v_cn[bd_idx[i]].fill(0);
		}

		/*for (int i = 0; i < in_idx.size(); i++)
		{
			if ((v_cn[in_idx[i]] | v_n[in_idx[i]]) < 0)
				v_cn[in_idx[i]] *= (-1);
		}*/

		for (int i = 0; i < in_idx.size(); i++)
		{
			//sqr_err += massSpring(v_cn[in_idx[i]], v_Th[in_idx[i]]);
			if((v_cn[in_idx[i]] | v_n[in_idx[i]]) >= 0)
				sqr_err += massSpring(v_cn[in_idx[i]], v_Th[in_idx[i]]);
			else
				sqr_err += massSpring(v_cn[in_idx[i]], -v_Th[in_idx[i]]);
		}

		int fidx, hidx, nhidx, nnhidx, nvidx, nnvidx;
		for (int i = 0; i < ach_idx.size(); i++)
		{
			for (auto vf : mesh.vf_range(VH(ach_idx[i])))
			{
				fidx = vf.idx();
				for (int j = 0; j < 3; j++)
				{
					if (hv_idx[fh_idx[fidx][j]] == ach_idx[i])
					{
						hidx = j;
						break;
					}
				}

				nhidx = (hidx + 1) % 3;
				nnhidx = (nhidx + 1) % 3;

				nvidx = hv_idx[fh_idx[fidx][nhidx]];
				nnvidx = hv_idx[fh_idx[fidx][nnhidx]];

				multiSub(v_p[ach_idx[i]], fh_cot[fidx][hidx] + fh_cot[fidx][nhidx], v_cn[ach_idx[i]]);
				multiAdd(v_p[ach_idx[i]], fh_cot[fidx][nhidx], v_cn[nvidx]);
				multiAdd(v_p[ach_idx[i]], fh_cot[fidx][hidx], v_cn[nnvidx]);

			}
		}

		for (int i = 0; i < fh_idx.size(); i++)
		{
			sub(v_cn[hv_idx[fh_idx[i][0]]], v_cn[hv_idx[fh_idx[i][2]]], fh_d[i][0]);
			sub(v_cn[hv_idx[fh_idx[i][1]]], v_cn[hv_idx[fh_idx[i][0]]], fh_d[i][1]);
			sub(v_cn[hv_idx[fh_idx[i][2]]], v_cn[hv_idx[fh_idx[i][1]]], fh_d[i][2]);

			fh_d[i][0] *= fh_cot[i][0];
			fh_d[i][1] *= fh_cot[i][1];
			fh_d[i][2] *= fh_cot[i][2];
		}

		for (int i = 0; i < v_cn.size(); i++)
		{
			v_cn[i].fill(0);
		}

		for (int i = 0; i < fh_idx.size(); i++)
		{
			subOn(fh_d[i][0], fh_d[i][1], v_cn[hv_idx[fh_idx[i][0]]]);
			subOn(fh_d[i][1], fh_d[i][2], v_cn[hv_idx[fh_idx[i][1]]]);
			subOn(fh_d[i][2], fh_d[i][0], v_cn[hv_idx[fh_idx[i][2]]]);
		}

		for (int i = 0; i < var_idx.size(); i++)
		{
			Ln[i] = v_cn[var_idx[i]][0];
		}

		Ln = Ln + var_idx.size();
		for (int i = 0; i < var_idx.size(); i++)
		{
			Ln[i] = v_cn[var_idx[i]][1];
		}

		Ln = Ln + var_idx.size();
		for (int i = 0; i < var_idx.size(); i++)
		{
			Ln[i] = v_cn[var_idx[i]][2];
		}

		return sqr_err;
	}

	double edgeSpringDeviation(double *Ed)
	{
		double sqr_err = 0.0;

		for (int i = 0; i < v_cn.size(); i++)
		{
			v_cn[i].fill(0);
		}

		for (int i = 0; i < fh_idx.size(); i++)
		{
			sub(v_p[hv_idx[fh_idx[i][0]]], v_p[hv_idx[fh_idx[i][2]]], fh_d[i][0]);
			sub(v_p[hv_idx[fh_idx[i][1]]], v_p[hv_idx[fh_idx[i][0]]], fh_d[i][1]);
			sub(v_p[hv_idx[fh_idx[i][2]]], v_p[hv_idx[fh_idx[i][1]]], fh_d[i][2]);
			
			sqr_err += massSpring(fh_d[i][0], e_l[fh_idx[i][0] >> 1]);
			sqr_err += massSpring(fh_d[i][1], e_l[fh_idx[i][1] >> 1]);
			sqr_err += massSpring(fh_d[i][2], e_l[fh_idx[i][2] >> 1]);

			fh_d[i][0] *= fh_cot[i][0];
			fh_d[i][1] *= fh_cot[i][1];
			fh_d[i][2] *= fh_cot[i][2];

			subOn(fh_d[i][0], fh_d[i][1], v_cn[hv_idx[fh_idx[i][0]]]);
			subOn(fh_d[i][1], fh_d[i][2], v_cn[hv_idx[fh_idx[i][1]]]);
			subOn(fh_d[i][2], fh_d[i][0], v_cn[hv_idx[fh_idx[i][2]]]);
		}

		int fidx, hidx, nhidx, nnhidx, nvidx, nnvidx;
		for (int i = 0; i < ach_idx.size(); i++)
		{
			for (auto vf : mesh.vf_range(VH(ach_idx[i])))
			{
				fidx = vf.idx();
				for (int j = 0; j < 3; j++)
				{
					if (hv_idx[fh_idx[fidx][j]] == ach_idx[i])
					{
						hidx = j;
						break;
					}
				}

				nhidx = (hidx + 1) % 3;
				nnhidx = (nhidx + 1) % 3;

				nvidx = hv_idx[fh_idx[fidx][nhidx]];
				nnvidx = hv_idx[fh_idx[fidx][nnhidx]];

				multiAdd(v_p[ach_idx[i]], fh_cot[fidx][nhidx], v_cn[nvidx]);
				multiAdd(v_p[ach_idx[i]], fh_cot[fidx][hidx], v_cn[nnvidx]);
			}
		}

		for (int i = 0; i < var_idx.size(); i++)
		{
			Ed[i] = v_cn[var_idx[i]][0];
		}

		Ed = Ed + var_idx.size();
		for (int i = 0; i < var_idx.size(); i++)
		{
			Ed[i] = v_cn[var_idx[i]][1];
		}

		Ed = Ed + var_idx.size();
		for (int i = 0; i < var_idx.size(); i++)
		{
			Ed[i] = v_cn[var_idx[i]][2];
		}

		return sqr_err;
	}

	void deformation()
	{
		if (deformH_cap == 0) return;

		int count, upcount;
		double last_err, curr_err;

		assembleLaplaceTriplets();
		assembleLaplaceMatrix();
		assembleLaplaceColMatrix();

		Eigen::SparseMatrix<double, Eigen::ColMajor> mat_H = lambdaH * mat_C.transpose()*mat_C;
		cllt.compute(mat_H + mat_L);
		assert(cllt.info() == Eigen::Success);
	
		Eigen::MatrixXd Ln(var_idx.size(), 3);
		Eigen::MatrixXd Ed(var_idx.size(), 3);

		AndersonAcceleration aa;
		aa.init(andersonM, var_idx.size() * 3, v_pos.data());

		count = 0;
		upcount = 0;

		setTargetMeanCurvature();
		curr_err = lambdaH * laplaceSpringNormal(Ln.data()) + edgeSpringDeviation(Ed.data());
		last_err = curr_err;
		Ed += lambdaH * Ln;
		std::cout << "Deformation Loss:......    " << "the " << count << "th iter: " << curr_err << std::endl;	
		
		do
		{
			v_pos = cllt.solve(Ed);
			v_pos = Eigen::Map<const Eigen::MatrixXd>(aa.compute(v_pos.data()).data(), var_idx.size(), 3);

			double *px = v_pos.data();
			double *py = v_pos.data() + var_idx.size();
			double *pz = v_pos.data() + 2 * var_idx.size();
			for (int i = 0; i < var_idx.size(); i++)
			{
				v_p[var_idx[i]][0] = px[i];
				v_p[var_idx[i]][1] = py[i];
				v_p[var_idx[i]][2] = pz[i];
			}

			updateTargetMeanCurvature();
			count++;

			aa.replace(v_pos.data());

			curr_err = lambdaH * laplaceSpringNormal(Ln.data()) + edgeSpringDeviation(Ed.data());
			Ed += lambdaH * Ln;

			std::cout << "Deformation Loss:......    " << "the " << count << "th iter: " << curr_err << std::endl;			
		} while (count < deformH_cap);

		if (deformG_cap == 0) return;
		cllt.compute(mat_L);
		assert(cllt.info() == Eigen::Success);

		aa.init(andersonM, var_idx.size() * 3, v_pos.data());
		count = 0;
		setTargetMeanCurvature();
		curr_err = edgeSpringDeviation(Ed.data());
		last_err = curr_err;
		std::cout << "Deformation Loss:......    " << "the " << count << "th iter: " << curr_err << std::endl;
		
		do
		{
			v_pos = cllt.solve(Ed);
			v_pos = Eigen::Map<const Eigen::MatrixXd>(aa.compute(v_pos.data()).data(), var_idx.size(), 3);

			double *px = v_pos.data();
			double *py = v_pos.data() + var_idx.size();
			double *pz = v_pos.data() + 2 * var_idx.size();
			for (int i = 0; i < var_idx.size(); i++)
			{
				v_p[var_idx[i]][0] = px[i];
				v_p[var_idx[i]][1] = py[i];
				v_p[var_idx[i]][2] = pz[i];
			}

			updateTargetMeanCurvature();
			count++;

			aa.replace(v_pos.data());

			curr_err = edgeSpringDeviation(Ed.data());
			std::cout << "Deformation Loss:......    " << "the " << count << "th iter: " << curr_err << std::endl;
		} while (count < deformG_cap);
	}
};

#endif // !CURVATUREDEFORMATION_H