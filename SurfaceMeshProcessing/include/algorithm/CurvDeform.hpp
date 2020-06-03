#ifndef CURVDEFORM_H
#define CURVDEFORM_H

#ifdef _DEBUG
#define verify(f) assert(f)
#else
#define verify(f) ((void)(f))
#endif

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <iostream>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <ctime>
#include <algorithm/AffineMap.hpp>
#include <algorithm/segMesh.hpp>
#include <algorithm/domainMesh.hpp>
#include <algorithm/AndersonAcceleration.hpp>
#define THETATHRES 1 // cos value for smallest angle of mesh
#define THETAFLIP (M_PI + 5e-2) // edge flip threshold for Delaunay triangulation

enum shapeKind
{
	sphere,
	cylinder,
	constantH,
	constantG,
	polygon,
	spherePatch,
	nonIter
};

struct opParas
{
	shapeKind kind = sphere; // the target curvature of optimizing seg
	double curv_dilat = 2.0 * M_PI; // the target curve for spherelize
	double curv_minus = 0.0 * M_PI; // the target mean curve for first circle in cylinderize
	double curv_mean = 0.0;
	double curv_gauss = 0.0 * M_PI;
	double curv_patch = 2.0 * M_PI;
	double curv_smooth = 0.0;

	int polygonNum = 4;

	// parameters in quasi Newton opt
	double shrink = 0.75; // the ratio to shrink step length
	int quas_cap = 30; // max iter number
	double curv_err = 5e-4; // the terminate condition 
	
	// parameters in local shape deformation
	double lambdaH = 0.1; // mean curvature weight
	int deform_cap = 60; // max iter number
	int andersonM = 5; // anderson acceleration basis num

	// parameters in ASAP assemble
	int ASAP_cap = 5; // max iter number
	bool fixBoundary = true; // fixed in global optimization
};

class curvDeform
{
private:
	Mesh &mesh;
	segMesh &seg;
public:
	Domain domain;

private:
	OpenMesh::VPropHandleT<double> para_u; // conformal factor of each vertex
	OpenMesh::EPropHandleT<double> para_l; // current length of each edge
	OpenMesh::HPropHandleT<double> op_theta; // opposite angle value of each halfedge
	OpenMesh::HPropHandleT<double> cot_theta; // cot value of op_theta
	OpenMesh::VPropHandleT<double> para_G; // gauss curvature of each vertex
	OpenMesh::VPropHandleT<double> para_TG; // target gauss curvature of each vertex
	OpenMesh::VPropHandleT<double> para_H; // mean curvature of each vertex
	OpenMesh::VPropHandleT<double> para_TH; // target mean curvature of each vertex
	OpenMesh::VPropHandleT<double> back_u; // backup of para_u in gauss Newton alg.
	OpenMesh::VPropHandleT<Eigen::Vector3d> back_p; // back up position of points on contours
	OpenMesh::VPropHandleT<Eigen::Vector3d> point0; // ASAP iteration begin point
	OpenMesh::VPropHandleT<double> mix_area; // mixed area of each vertex

	OpenMesh::VPropHandleT<double> kmax, kmin;

	OpenMesh::HPropHandleT<Eigen::Vector3d> slack_d; // slack variable in edge keep optimization
	OpenMesh::VPropHandleT<Eigen::Vector3d> slack_n; // slack variable in mean curvature optimization
	OpenMesh::VPropHandleT<AffineMap> affine_M; // AffineMap tranform for ASAP

	OpenMesh::EPropHandleT<bool> flip_indicate;

	std::vector<Eigen::Triplet<double>> dual_lap_trip; // triplets for dual laplacian matrix
	Eigen::SparseMatrix<double, Eigen::ColMajor> global_L; // sparse dual laplacian Matrix
	Eigen::SparseMatrix<double, Eigen::ColMajor> local_L; // local part of laplace Matrix
	Eigen::SparseMatrix<double, Eigen::ColMajor> local_C; // local cols part of laplace Matrix
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double, Eigen::ColMajor>> llt; // simplical cholesky solver

	std::vector<std::function<void(const opParas &)>> set_gauss_func; // functions for set different gauss curvatures
	std::vector<std::function<void(const opParas &)>> set_mean_func; // functions for set different mean curvatures
	std::vector<std::function<void(const opParas &)>> vproject_func; // functions for project contour to shape

private:

public:
	curvDeform(Mesh &input, segMesh &inseg) : mesh(input), seg(inseg), domain(input)
	{
		set_gauss_func.push_back(std::bind(&curvDeform::setSphereGaussCurvature, this, std::placeholders::_1));
		set_gauss_func.push_back(std::bind(&curvDeform::setCylinderGaussCurvature, this, std::placeholders::_1));
		set_gauss_func.push_back(std::bind(&curvDeform::setConstantHGaussCurvature, this, std::placeholders::_1));
		set_gauss_func.push_back(std::bind(&curvDeform::setConstantGGaussCurvature, this, std::placeholders::_1));
		set_gauss_func.push_back(std::bind(&curvDeform::setPolygonGaussCurvature, this, std::placeholders::_1));
		set_gauss_func.push_back(std::bind(&curvDeform::setSpherePatchGaussCurvature, this, std::placeholders::_1));
		set_gauss_func.push_back(std::bind(&curvDeform::setNonIterGaussCurvature, this, std::placeholders::_1));

		set_mean_func.push_back(std::bind(&curvDeform::setSphereMeanCurvature, this, std::placeholders::_1));
		set_mean_func.push_back(std::bind(&curvDeform::setCylinderMeanCurvature, this, std::placeholders::_1));
		set_mean_func.push_back(std::bind(&curvDeform::setConstantHMeanCurvature, this, std::placeholders::_1));
		set_mean_func.push_back(std::bind(&curvDeform::setConstantGMeanCurvature, this, std::placeholders::_1));
		set_mean_func.push_back(std::bind(&curvDeform::setPolygonMeanCurvature, this, std::placeholders::_1));
		set_mean_func.push_back(std::bind(&curvDeform::setSpherePatchMeanCurvature, this, std::placeholders::_1));
		set_mean_func.push_back(std::bind(&curvDeform::setNonIterMeanCurvature, this, std::placeholders::_1));


		vproject_func.push_back(std::bind(&curvDeform::projectSphere, this, std::placeholders::_1));
		vproject_func.push_back(std::bind(&curvDeform::projectCylinder, this, std::placeholders::_1));

		addProperty();
	}

	void run(opParas &op)
	{
		OpenMesh::VPropHandleT<OpenMesh::VertexHandle> pair_v;
		mesh.add_property(pair_v);
		for (auto v : mesh.vertices()) mesh.property(pair_v, v) = v;
		Mesh test = mesh;

		test.delete_vertex(Mesh::VertexHandle(1));
		test.garbage_collection();

		auto test_pts = test.points();
		auto mesh_pts = mesh.points();
		std::cout << test.n_vertices() << std::endl;
		std::cout << test.property(pair_v, OpenMesh::VertexHandle(0)).idx() << std::endl;
		std::cout << test.property(pair_v, OpenMesh::VertexHandle(1)).idx() << std::endl;
		std::cout << test.property(pair_v, OpenMesh::VertexHandle(test.n_vertices() - 1)).idx() << std::endl;

		//std::clock_t start, end;
		//start = clock();

		//VSet vars, achs;
		//initProperty();

		//vars = seg.inner();
		//achs = VSet(1, seg.inner()[0]);
		//for (int i = 0; i < seg.n_contours(); i++)
		//{
		//	vars.insert(vars.end(), seg.contour(i).begin(), seg.contour(i).end());
		//}
		//domain.build(vars, achs);

		//verify(updateDomainTheta(THETATHRES));
		//domainEdgeFlip();
		//updateDomainMixArea();
		//updateDomainGaussCurvature();
		////updateDomainMeanCurvature();

		//std::cout << "Calabi flow ...\n";
		//quasiNewton(op);

		//std::cout << "Local deformation ...\n";
		//op.lambdaH = 0.1;
		//op.deform_cap = 60;
		//shapeOpt(op);

		//std::cout << "Global assemble ...\n";
		//initLocalPos();
		//for (auto c : seg.vContours())
		//{
		//	for (auto v : c)
		//	{
		//		mesh.property(back_p, v) = mesh.point(v);
		//	}
		//}

		//initRemainPos();
		//for (int i = 0; i < seg.n_outers(); i++)
		//{
		//	for (auto v : seg.outer(i))
		//	{
		//		mesh.set_point(v, mesh.property(point0, v));
		//	}
		//}

		//for (auto v : seg.inner())
		//{
		//	mesh.property(point0, v) = mesh.point(v);
		//}
		//
		//std::vector<VSet> inners, innerContours, outers, outerContours;
		//int shiftSize = 10;
		//seg.extractBroads(inners, innerContours, outers, outerContours, shiftSize);

		//for (int i = 0; i < outers.size(); i++)
		//{
		//	domain.build(outers[i], outerContours[i]);
		//	laplacePoint(5);
		//}

		//mesh.update_normals();
		//end = clock();
		//std::cout << "Total time : " << (double)(end - start) << "ms\n";
	}

	~curvDeform()
	{
		removeProperty();
	}

	void estimateCurvature(opParas &op)
	{
		double areaSum = 0.0;
		for (auto v : seg.inner())
		{
			areaSum += mesh.property(mix_area, v);
		}

		double circleArea = 0.0;
		for (auto l : seg.lContours())
		{
			circleArea += l * l / 4 / M_PI;
		}

		op.curv_dilat = areaSum / (areaSum + circleArea) * 4 * M_PI;
	}

	void ARAP(const OpenMesh::VertexHandle &anchor, const Eigen::Vector3d &topos, const opParas &op)
	{
		VSet vars(seg.inner());
		VSet achs;
		for (int i = 0; i < seg.n_contours(); i++)
		{
			if (seg.bContours()[i])
			{
				achs.insert(achs.end(), seg.contour(i).begin(), seg.contour(i).end());
			}
			else
			{
				vars.insert(vars.end(), seg.contour(i).begin(), seg.contour(i).end());
			}
		}
		achs.push_back(anchor);
		domain.build(vars, achs);

		for (auto v : domain.vertices())
		{
			mesh.property(affine_M, v) = AffineMap();
			mesh.property(point0, v) = mesh.point(v);
		}

		for (auto e : domain.edges())
		{
			mesh.property(para_l, e) = mesh.calc_edge_vector(e).norm();
		}

		for (auto f : domain.faces())
		{
			for (auto fh : mesh.fh_range(f))
			{
				verify(calcHeOppositeTheta(fh, mesh.property(op_theta, fh), THETATHRES));
			}
		}

		mesh.set_point(anchor, topos);

		assembleOpt(op, false);
	}

	void constantHSurface(double uniH)
	{
		VSet vars, achs;

		vars = seg.inner();
		for (int i = 0; i < seg.n_contours(); i++)
		{
			if (!seg.bContours()[i]) continue;
			achs.insert(achs.end(), seg.contour(i).begin(), seg.contour(i).end());
		}
		domain.build(vars, achs);

		for (auto v : domain.vertices())
		{
			mesh.property(para_u, v) = 0;
			mesh.property(para_TH, v) = 0;
		}

		int iterNum = 0;
		do {
			updateDomainMesh(THETATHRES);
			assembleLaplaceCol();
			llt.compute(local_C.transpose()*local_C);
			assert(llt.info() == Eigen::Success);

			int count = 0;
			double cur_err = 0.0;
			Eigen::MatrixXd Ln(domain.n_variables(), 3);
			Eigen::MatrixXd pos(domain.n_variables(), 3);

			for (auto v : mesh.vertices())
			{
				mesh.property(point0, v) = mesh.point(v);
			}

			for (auto v : domain.variables())
			{
				pos.row(domain.index(v)) = mesh.point(v);
			}

			AndersonAcceleration aa;
			aa.init(10, domain.n_variables() * 3, pos.data());

			cur_err = updateLslackN(Ln);
			std::cout << "Deform Loss:......    " << "the " << count << "th iter: " << cur_err << std::endl;
			do
			{
				pos = llt.solve(Ln);
				pos = Eigen::Map<const Eigen::MatrixXd>(aa.compute(pos.data()).data(), domain.n_variables(), 3);
				for (auto v : domain.variables())
				{
					mesh.set_point(v, pos.row(domain.index(v)));
				}
				updateDomainMeanCurvature();
				count++;

				aa.replace(pos.data());
				cur_err = updateLslackN(Ln);
				std::cout << "Deform Loss:......    " << "the " << count << "th iter: " << cur_err << std::endl;
			} while (count < 60);
			iterNum++;
		} while (iterNum < 5);

	}

	void planeParalization(opParas &op)
	{
		std::clock_t start, end;
		start = clock();

		VSet vars, achs;
		initProperty();

		vars = seg.inner();
		for (int i = 0; i < seg.n_contours(); i++)
		{
			if (!seg.bContours()[i]) continue;
			achs.insert(achs.end(), seg.contour(i).begin(), seg.contour(i).end());
		}
		domain.build(vars, achs);

		verify(updateDomainTheta(THETATHRES));
		domainEdgeFlip();
		updateDomainMixArea();
		updateDomainGaussCurvature();
		//updateDomainMeanCurvature();

		op.kind = constantG;
		op.curv_gauss = 0;
		op.curv_mean = 0;
		op.deform_cap = 5000;

		std::cout << "Calabi flow ...\n";
		quasiNewton(op);

		vars = seg.inner();
		for (int i = 0; i < seg.n_contours(); i++)
		{
			if (!seg.bContours()[i]) continue;
			vars.insert(vars.end(), seg.contour(i).begin(), seg.contour(i).end());
		}
		achs = VSet(1, seg.inner()[0]);
		domain.build(vars, achs);

		std::cout << "Local deformation ...\n";
		shapeOpt(op);

		mesh.update_normals();
		end = clock();
		std::cout << "Total time : " << (double)(end - start) << "ms\n";
	}

	void spherePatchlization(opParas &op)
	{
		std::clock_t start, end;
		start = clock();

		VSet vars, achs;
		initProperty();

		vars = seg.inner();
		for (int i = 0; i < seg.n_contours(); i++)
		{
			if (!seg.bContours()[i]) continue;
			vars.insert(vars.end(), seg.contour(i).begin(), seg.contour(i).end());
		}
		achs = {};
		achs.push_back(seg.inner()[0]);
		domain.build(vars, achs);

		verify(updateDomainTheta(THETATHRES));
		domainEdgeFlip();
		updateDomainMixArea();
		updateDomainGaussCurvature();
		//updateDomainMeanCurvature();

		op.kind = spherePatch;
		op.polygonNum = 5;
		//op.curv_patch = 3.0 * M_PI;
		op.quas_cap = 20;
		op.deform_cap = 600;

		std::cout << "Calabi flow ...\n";
		quasiNewton(op);

		vars = seg.inner();
		for (int i = 0; i < seg.n_contours(); i++)
		{
			if (!seg.bContours()[i]) continue;
			vars.insert(vars.end(), seg.contour(i).begin(), seg.contour(i).end());
		}
		achs = VSet(1, seg.inner()[0]);
		domain.build(vars, achs);

		std::cout << "Local deformation ...\n";
		shapeOpt(op);

		mesh.update_normals();
		end = clock();
		std::cout << "Total time : " << (double)(end - start) << "ms\n";
	}

	void fastPlanarizetoShape(opParas &op)
	{
		std::clock_t start, end;
		start = clock();

		Mesh output = mesh;

		VSet vars, achs;
		initProperty();

		vars = seg.inner();
		for (int i = 0; i < seg.n_contours(); i++)
		{
			if (!seg.bContours()[i]) continue;
			vars.insert(vars.end(), seg.contour(i).begin(), seg.contour(i).end());
		}
		achs = {};
		achs.push_back(seg.inner()[0]);
		domain.build(vars, achs);

		verify(updateDomainTheta(THETATHRES));
		domainEdgeFlip();
		updateDomainMixArea();
		updateDomainGaussCurvature();
		//updateDomainMeanCurvature();

		op.kind = polygon;
		op.polygonNum = 5;
		op.curv_mean = 0;
		op.deform_cap = 60;

		std::cout << "Calabi flow ...\n";
		quasiNewton(op);

		vars = seg.inner();
		for (int i = 0; i < seg.n_contours(); i++)
		{
			if (!seg.bContours()[i]) continue;
			vars.insert(vars.end(), seg.contour(i).begin(), seg.contour(i).end());
		}
		achs = VSet(1, seg.inner()[0]);
		domain.build(vars, achs);

		std::cout << "Projection ...\n";
		calc_uv_coordinates();

		std::cout << "Local deformation ...\n";
		assembleLaplace();
		llt.compute(local_L);
		assert(llt.info() == Eigen::Success);

		int count = 0;
		double cur_err = 0.0;
		Eigen::MatrixXd Ed(domain.n_variables(), 2);
		Eigen::MatrixXd pos(domain.n_variables(), 2);

		OpenMesh::VPropHandleT<Eigen::Vector2d> uv0;
		OpenMesh::HPropHandleT<Eigen::Vector2d> slack_uv;
		mesh.add_property(uv0);
		mesh.add_property(slack_uv);

		for (auto v : mesh.vertices())
		{
			mesh.property(uv0, v) = mesh.texcoord2D(v);
		}

		for (auto h : mesh.halfedges())
		{
			mesh.property(slack_uv, h) = Eigen::Vector2d(0, 0);
		}

		for (auto v : domain.variables())
		{
			pos.row(domain.index(v)) = mesh.texcoord2D(v);
		}

		AndersonAcceleration aa;
		aa.init(op.andersonM, domain.n_variables() * 2, pos.data());

		cur_err = updateEslackUV(Ed, slack_uv);
		std::cout << "Deform Loss:......    " << "the " << count << "th iter: " << cur_err << std::endl;
		do
		{
			pos = llt.solve(Ed);
			pos = Eigen::Map<const Eigen::MatrixXd>(aa.compute(pos.data()).data(), domain.n_variables(), 2);
			for (auto v : domain.variables())
			{
				mesh.set_texcoord2D(v, pos.row(domain.index(v)));
			}
			count++;

			aa.replace(pos.data());
			cur_err = updateEslackUV(Ed, slack_uv);
			std::cout << "Deform Loss:......    " << "the " << count << "th iter: " << cur_err << std::endl;
		} while (count < op.deform_cap);

		mesh.update_normals();

		mesh.remove_property(uv0);
		mesh.remove_property(slack_uv);

		end = clock();
		std::cout << "Total time : " << (double)(end - start) << "ms\n";

		unifyTexcoords2D();
		for (auto v : domain.vertices())
		{
			output.set_point(v, Eigen::Vector3d(mesh.texcoord2D(v)[0], mesh.texcoord2D(v)[1], 0));
		}
		for (auto f : output.faces())
		{
			bool outDomain = false;
			for (auto fv : mesh.fv_range(f))
			{
				if (domain.index(fv) == -2)
				{
					outDomain = true;
					break;
				}
			}
			if (outDomain)
			{
				output.delete_face(f);
			}
		}

		output.garbage_collection();

		OpenMesh::IO::Options ver_tex = OpenMesh::IO::Options::FaceTexCoord;
		OpenMesh::IO::write_mesh(output, "uv_plane.obj", ver_tex);
	}

	void fastPlanarize(opParas &op)
	{
		std::clock_t start, end;
		start = clock();

		Mesh output = mesh;

		VSet vars, achs;
		initProperty();

		vars = seg.inner();
		for (int i = 0; i < seg.n_contours(); i++)
		{
			if (!seg.bContours()[i]) continue;
			achs.insert(achs.end(), seg.contour(i).begin(), seg.contour(i).end());
		}
		domain.build(vars, achs);

		verify(updateDomainTheta(THETATHRES));
		domainEdgeFlip();
		updateDomainMixArea();
		updateDomainGaussCurvature();
		//updateDomainMeanCurvature();

		op.kind = constantG;
		op.curv_gauss = 0;
		op.curv_mean = 0;
		op.deform_cap = 60;

		std::cout << "Calabi flow ...\n";
		quasiNewton(op);

		vars = seg.inner();
		for (int i = 0; i < seg.n_contours(); i++)
		{
			if (!seg.bContours()[i]) continue;
			vars.insert(vars.end(), seg.contour(i).begin(), seg.contour(i).end());
		}
		achs = VSet(1, seg.inner()[0]);
		domain.build(vars, achs);

		std::cout << "Projection ...\n";
		calc_uv_coordinates();

		std::cout << "Local deformation ...\n";
		assembleLaplace();
		llt.compute(local_L);
		assert(llt.info() == Eigen::Success);

		int count = 0;
		double cur_err = 0.0;
		Eigen::MatrixXd Ed(domain.n_variables(), 2);
		Eigen::MatrixXd pos(domain.n_variables(), 2);

		OpenMesh::VPropHandleT<Eigen::Vector2d> uv0;
		OpenMesh::HPropHandleT<Eigen::Vector2d> slack_uv;
		mesh.add_property(uv0);
		mesh.add_property(slack_uv);

		for (auto v : mesh.vertices())
		{
			mesh.property(uv0, v) = mesh.texcoord2D(v);
		}

		for (auto h : mesh.halfedges())
		{
			mesh.property(slack_uv, h) = Eigen::Vector2d(0, 0);
		}

		for (auto v : domain.variables())
		{
			pos.row(domain.index(v)) = mesh.texcoord2D(v);
		}

		AndersonAcceleration aa;
		aa.init(op.andersonM, domain.n_variables() * 2, pos.data());

		cur_err = updateEslackUV(Ed, slack_uv);
		std::cout << "Deform Loss:......    " << "the " << count << "th iter: " << cur_err << std::endl;
		do
		{
			pos = llt.solve(Ed);
			pos = Eigen::Map<const Eigen::MatrixXd>(aa.compute(pos.data()).data(), domain.n_variables(), 2);
			for (auto v : domain.variables())
			{
				mesh.set_texcoord2D(v, pos.row(domain.index(v)));
			}
			count++;

			aa.replace(pos.data());
			cur_err = updateEslackUV(Ed, slack_uv);
			std::cout << "Deform Loss:......    " << "the " << count << "th iter: " << cur_err << std::endl;
		} while (count < op.deform_cap);

		mesh.update_normals();

		mesh.remove_property(uv0);
		mesh.remove_property(slack_uv);

		end = clock();
		std::cout << "Total time : " << (double)(end - start) << "ms\n";

		unifyTexcoords2D();
		for (auto v : domain.vertices())
		{
			output.set_point(v, Eigen::Vector3d(mesh.texcoord2D(v)[0], mesh.texcoord2D(v)[1], 0));
		}
		for (auto f : output.faces())
		{
			bool outDomain = false;
			for (auto fv : mesh.fv_range(f))
			{
				if (domain.index(fv) == -2)
				{
					outDomain = true;
					break;
				}
			}
			if (outDomain)
			{
				output.delete_face(f);
			}
		}

		output.garbage_collection();

		OpenMesh::IO::Options ver_tex = OpenMesh::IO::Options::FaceTexCoord;
		OpenMesh::IO::write_mesh(output, "uv_plane.obj", ver_tex);
	}

	void smoothbyCurvature(opParas &op)
	{
		std::clock_t start, end;
		start = clock();

		VSet vars, achs;
		initProperty();

		vars = seg.inner();
		for (int i = 0; i < seg.n_contours(); i++)
		{
			if (!seg.bContours()[i]) continue;
			vars.insert(vars.end(), seg.contour(i).begin(), seg.contour(i).end());
		}
		achs = {};
		achs.push_back(seg.inner()[0]);
		domain.build(vars, achs);

		verify(updateDomainTheta(THETATHRES));
		domainEdgeFlip();
		updateDomainMixArea();
		updateDomainGaussCurvature();
		updateDomainMeanCurvature();

		op.kind = nonIter;
		op.lambdaH = 1;

		double vH, vG;
		for (auto v : domain.vertices())
		{
			vH = mesh.property(para_H, v) / mesh.property(mix_area, v);
			vG = mesh.property(para_G, v) / mesh.property(mix_area, v);
			mesh.property(kmax, v) = (vH + sqrt(vH*vH - 4 * vG)) / 2;
			mesh.property(kmin, v) = (vH - sqrt(vH*vH - 4 * vG)) / 2;
		}

		double mean_G = 0;
		for (auto v : seg.inner())
		{
			mean_G += mesh.property(para_G, v);
		}
		mean_G /= seg.inner().size();

		double var_G = 0;
		for (auto v : seg.inner())
		{
			var_G += pow(mesh.property(para_G, v) - mean_G, 2);
		}
		var_G /= seg.inner().size();
		var_G = sqrt(var_G);

		for (auto v : domain.vertices())
		{
			double maxk = mesh.property(kmax, v);
			double mink = mesh.property(kmin, v);
			if (maxk + mink < 0)
			{
				double tmp = mink;
				mink = maxk;
				maxk = tmp;
			}
			double ratio = exp(op.curv_smooth*4);
			mesh.property(para_TH, v) = (maxk * ratio + mink / ratio) / (maxk + mink)*mesh.property(para_H, v);
			mesh.property(para_TG, v) = mesh.property(para_G, v);
		}

		//double sum_G = 0;
		//int num_G = 0;
		//double lambdad = 0.5, lambdal = -0.2;
		//for (auto v : domain.vertices())
		//{
		//	//mesh.property(para_TH, v) = mesh.property(para_H, v);
		//	if (mesh.property(para_G, v) > mean_G + lambdad * var_G || mesh.property(para_G, v) < mean_G - lambdad * var_G)
		//	{
		//		mesh.property(para_TG, v) = (1 + lambdal) * mesh.property(para_G, v);
		//		sum_G += lambdal *mesh.property(para_G, v);
		//	}
		//	else
		//	{
		//		num_G++;
		//	}
		//}

		//for (auto v : domain.vertices())
		//{
		//	if (mesh.property(para_G, v) <= mean_G + lambdad * var_G && mesh.property(para_G, v) >= mean_G - lambdad * var_G)
		//	{
		//		mesh.property(para_TG, v) = mesh.property(para_G, v) - sum_G / num_G;
		//	}
		//}

		//double sum_TG = 0;
		//for (auto v : domain.vertices())
		//{
		//	sum_TG += mesh.property(para_TG, v);
		//}
		//std::cout << "SUM_TG: " << sum_TG << std::endl;

		/*for (auto v : seg.inner())
		{
			mesh.property(para_TG, v) = mesh.property(para_G, v) + op.curv_smooth*(mesh.property(para_G, v) - mean_G);
			mesh.property(para_TH, v) = mesh.property(para_H, v);
		}
		for (auto contour : seg.vContours())
		{
			for (auto v : contour)
			{
				mesh.property(para_TG, v) = mesh.property(para_G, v);
				mesh.property(para_TH, v) = mesh.property(para_H, v);
			}
		}*/

		std::cout << "Calabi flow ...\n";
		quasiNewton(op);

		vars = seg.inner();
		for (int i = 0; i < seg.n_contours(); i++)
		{
			if (!seg.bContours()[i]) continue;
			vars.insert(vars.end(), seg.contour(i).begin(), seg.contour(i).end());
		}
		achs = VSet(1, seg.inner()[0]);
		domain.build(vars, achs);

		std::cout << "Local deformation ...\n";
		shapeOpt(op);

		std::cout << "Global assemble ...\n";
		initLocalPos();
		//vproject_func[op.kind](op);
		for (auto c : seg.vContours())
		{
			for (auto v : c)
			{
				mesh.property(back_p, v) = mesh.point(v);
			}
		}

		initRemainPos();

		for (int i = 0; i < seg.n_outers(); i++)
		{
			for (auto v : seg.outer(i))
			{
				mesh.set_point(v, mesh.property(point0, v));
			}
		}

		for (auto v : seg.inner())
		{
			mesh.property(point0, v) = mesh.point(v);
		}

		std::vector<VSet> inners, innerContours, outers, outerContours;
		int shiftSize = 10;
		seg.extractBroads(inners, innerContours, outers, outerContours, shiftSize);

		for (int i = 0; i < outers.size(); i++)
		{
			domain.build(outers[i], outerContours[i]);
			laplacePoint(5);
		}

		mesh.update_normals();
		end = clock();
		std::cout << "Total time : " << (double)(end - start) << "ms\n";
	}

	void shapeOpt1(VSet move_handle, std::vector<Eigen::Vector3d> move_position, opParas &op)
	{
		op.lambdaH = 1;
		initProperty();
		updateDomainMeanCurvature();
		for (auto v : domain.vertices())
		{
			mesh.property(para_TH, v) = mesh.property(para_H, v);
		}
		assembleLaplace();
		assembleLaplaceCol();
		llt.compute(op.lambdaH*local_C.transpose()*local_C + local_L);
		assert(llt.info() == Eigen::Success);

		int count = 0;
		double cur_err = 0.0;
		Eigen::MatrixXd Ln(domain.n_variables(), 3);
		Eigen::MatrixXd Ed(domain.n_variables(), 3);
		Eigen::MatrixXd pos(domain.n_variables(), 3);

		for (auto v : mesh.vertices())
		{
			mesh.property(point0, v) = mesh.point(v);
		}

		for (size_t i = 0; i < move_handle.size(); i++)
		{
			mesh.set_point(move_handle[i], move_position[i]);
		}

		for (auto v : domain.variables())
		{
			pos.row(domain.index(v)) = mesh.point(v);
		}

		AndersonAcceleration aa;
		aa.init(op.andersonM, domain.n_variables() * 3, pos.data());

		//set_mean_func[op.kind](op);
		cur_err = op.lambdaH * updateLslackN(Ln) + updateEslackD(Ed);
		std::cout << "Deform Loss:......    " << "the " << count << "th iter: " << cur_err << std::endl;
		do
		{
			//std::cout << pos << std::endl;
			pos = llt.solve(op.lambdaH*Ln + Ed);
			//std::cout << Ln << std::endl;
			//std::cout << Ed << std::endl;
			pos = Eigen::Map<const Eigen::MatrixXd>(aa.compute(pos.data()).data(), domain.n_variables(), 3);
			//std::cout << pos << std::endl;
			for (auto v : domain.variables())
			{
				mesh.set_point(v, pos.row(domain.index(v)));
			}
			//updateDomainMeanCurvature();
			count++;

			aa.replace(pos.data());
			cur_err = op.lambdaH * updateLslackN(Ln) + updateEslackD(Ed);
			std::cout << "Deform Loss:......    " << "the " << count << "th iter: " << cur_err << std::endl;
		} while (count < op.deform_cap);
	}


private:
	// ***** --------- Curvature Initialize --------- ***** //

	
	// Input gauss curvature sum should between 0 and 4*M_PI
	void setSphereGaussCurvature(const opParas &op)
	{
		assert(op.curv_dilat >= 0 && op.curv_dilat < 4 * M_PI);
		double areaSum = 0.0;
		for (auto v : seg.inner())
		{
			areaSum += mesh.property(mix_area, v);
		}
		double gauss_inner_curv = op.curv_dilat / areaSum;

		double sqrcircleLength = 0.0;
		for (int i = 0; i < seg.n_contours(); i++)
		{
			sqrcircleLength += pow(seg.lContours()[i], 2);
		}

		for (auto v : seg.inner())
		{
			mesh.property(para_TG, v) = gauss_inner_curv * mesh.property(mix_area, v);
		}

		for (int i = 0; i < seg.vContours().size(); i++)
		{
			double contourCurv = pow(seg.lContours()[i], 2) / sqrcircleLength * (4 * M_PI - op.curv_dilat) - 2 * M_PI;
			double contourLen = 0;
			for (int j = 0; j < seg.vContours()[i].size(); j++)
			{
				auto he = mesh.find_halfedge(seg.vContours()[i][j], seg.vContours()[i][(j + 1) % seg.vContours()[i].size()]);
				contourLen += mesh.property(para_l, mesh.edge_handle(he));
			}
			for (int j = 0; j < seg.vContours()[i].size(); j++)
			{
				double edgeLength = 0;

				auto he = mesh.find_halfedge(seg.vContours()[i][j], seg.vContours()[i][(j + 1) % seg.vContours()[i].size()]);
				edgeLength += mesh.property(para_l, mesh.edge_handle(he));
				he = mesh.find_halfedge(seg.vContours()[i][j],
					seg.vContours()[i][(j - 1 + seg.vContours()[i].size()) % seg.vContours()[i].size()]);
				edgeLength += mesh.property(para_l, mesh.edge_handle(he));

				mesh.property(para_TG, seg.vContours()[i][j]) = contourCurv / contourLen * edgeLength / 2;
			}
		}
	}

	void setCylinderGaussCurvature(const opParas &op)
	{
		double areaSum = 0.0;

		for (auto v : seg.inner())
		{
			areaSum += mesh.property(mix_area, v);
		}

		double gauss_inner_curv = 0;

		for (auto v : seg.inner())
		{
			mesh.property(para_TG, v) = gauss_inner_curv * mesh.property(mix_area, v);
		}

		int firstContour = 0;
		for (int i = 0; i < seg.vContours().size(); i++)
		{
			if (seg.bContours()[i])
			{
				double contourCurv = firstContour == 0 ? op.curv_minus : -op.curv_minus;
				double contourLen = 0;
				for (int j = 0; j < seg.vContours()[i].size(); j++)
				{
					auto he = mesh.find_halfedge(seg.vContours()[i][j], seg.vContours()[i][(j + 1) % seg.vContours()[i].size()]);
					contourLen += mesh.property(para_l, mesh.edge_handle(he));
				}
				for (int j = 0; j < seg.vContours()[i].size(); j++)
				{
					double edgeLength = 0;

					auto he = mesh.find_halfedge(seg.vContours()[i][j], seg.vContours()[i][(j + 1) % seg.vContours()[i].size()]);
					edgeLength += mesh.property(para_l, mesh.edge_handle(he));
					he = mesh.find_halfedge(seg.vContours()[i][j],
						seg.vContours()[i][(j - 1 + seg.vContours()[i].size()) % seg.vContours()[i].size()]);
					edgeLength += mesh.property(para_l, mesh.edge_handle(he));

					mesh.property(para_TG, seg.vContours()[i][j]) = contourCurv / contourLen * edgeLength / 2;
				}
				firstContour++;
			}
			else
			{
				double contourCurv = -2 * M_PI;
				double contourLen = 0;
				for (int j = 0; j < seg.vContours()[i].size(); j++)
				{
					auto he = mesh.find_halfedge(seg.vContours()[i][j], seg.vContours()[i][(j + 1) % seg.vContours()[i].size()]);
					contourLen += mesh.property(para_l, mesh.edge_handle(he));
				}
				for (int j = 0; j < seg.vContours()[i].size(); j++)
				{
					double edgeLength = 0;

					auto he = mesh.find_halfedge(seg.vContours()[i][j], seg.vContours()[i][(j + 1) % seg.vContours()[i].size()]);
					edgeLength += mesh.property(para_l, mesh.edge_handle(he));
					he = mesh.find_halfedge(seg.vContours()[i][j],
						seg.vContours()[i][(j - 1 + seg.vContours()[i].size()) % seg.vContours()[i].size()]);
					edgeLength += mesh.property(para_l, mesh.edge_handle(he));

					mesh.property(para_TG, seg.vContours()[i][j]) = contourCurv / contourLen * edgeLength / 2;
				}
			}
		}
	}

	void setConstantHGaussCurvature(const opParas &op)
	{

	}

	void setConstantGGaussCurvature(const opParas &op)
	{
		for (auto v : seg.inner())
		{
			mesh.property(para_TG, v) = op.curv_gauss*mesh.property(mix_area, v);
		}
	}

	void setPolygonGaussCurvature(const opParas &op)
	{
		assert(seg.n_contours() == 1);

		for (auto v : seg.inner())
		{
			mesh.property(para_TG, v) = 0;
		}

		int corner = (seg.contour(0).size() + 1) / op.polygonNum;
		for (int i = 0; i < seg.contour(0).size(); i++)
		{
			if (i%corner == 0)
			{
				mesh.property(para_TG, seg.contour(0)[i]) = 2 * M_PI / op.polygonNum;
			}
			else
			{
				mesh.property(para_TG, seg.contour(0)[i]) = 0;
			}
		}
	}

	void setSpherePatchGaussCurvature(const opParas &op)
	{
		assert(seg.n_contours() == 1);

		double areaSum = 0.0;
		double areaCorner = 0.0;

		for (auto v : seg.inner())
		{
			areaSum += mesh.property(mix_area, v);
		}
		int corner = (seg.contour(0).size() + 1) / op.polygonNum;
		for (int i = 0; i < seg.contour(0).size(); i++)
		{
			if (i%corner == 0)
			{
				areaCorner += mesh.property(mix_area, seg.contour(0)[i]);
			}
			areaSum += mesh.property(mix_area, seg.contour(0)[i]);
		}

		for (auto v : seg.inner())
		{
			mesh.property(para_TG, v) = op.curv_patch / areaSum * mesh.property(mix_area, v);
		}
		for (int i = 0; i < seg.contour(0).size(); i++)
		{
			if (i%corner == 0)
			{
				//mesh.property(para_TG, seg.contour(0)[i]) = 
				//	M_PI - (op.curv_patch + (op.polygonNum - 2)*M_PI) / op.polygonNum;
				mesh.property(para_TG, seg.contour(0)[i]) =
					(2 * M_PI - op.curv_patch*(1 - areaCorner / areaSum)) / op.polygonNum;
			}
			else
			{
				mesh.property(para_TG, seg.contour(0)[i]) = 
					op.curv_patch / areaSum * mesh.property(mix_area, seg.contour(0)[i]);
			}
		}
	}

	void setNonIterGaussCurvature(const opParas &op)
	{
		/*double sum_TG = 0;
		for (auto v : domain.vertices())
		{
			mesh.property(para_TG, v) = mesh.property(kmax, v)*mesh.property(kmin, v) * mesh.property(mix_area, v);
			mesh.property(para_TH, v) = (mesh.property(kmax, v) + mesh.property(kmin, v)) * mesh.property(mix_area, v);
			sum_TG += mesh.property(para_TG, v);
		}
		std::cout << sum_TG << std::endl;

		for (auto v : domain.vertices())
		{
			mesh.property(para_TG, v) /= sum_TG * 2 * M_PI;
			mesh.property(para_TH, v) /= sum_TG * 2 * M_PI;
		}*/
	}

	void setSphereMeanCurvature(const opParas &op)
	{
		assert(op.curv_dilat >= 0 && op.curv_dilat < 4 * M_PI);
		double areaSum = 0.0;
		for (auto v : seg.inner())
		{
			areaSum += mesh.property(mix_area, v);
		}

		double gauss_inner_curv = op.curv_dilat / areaSum;
		double mean_inner_curv = sqrt(gauss_inner_curv);
		double mean_edges_curv = 0;

		for (auto v : seg.inner())
		{
			mesh.property(para_TH, v) = 2 * mean_inner_curv * mesh.property(mix_area, v);
		}

		for (int i = 0; i < seg.vContours().size(); i++)
		{
			for (int j = 0; j < seg.vContours()[i].size(); j++)
			{
				mesh.property(para_TH, seg.vContours()[i][j]) = mean_edges_curv;
			}
		}
	}
	
	void setCylinderMeanCurvature(const opParas &op)
	{
		double contoursLen = 0.0;
		for (int i=0; i<seg.vContours().size(); i++)
		{
			if (!seg.bContours()[i]) continue;
			for (int j = 0; j < seg.vContours()[i].size(); j++)
			{
				auto he = mesh.find_halfedge(seg.vContours()[i][j], seg.vContours()[i][(j + 1) % seg.vContours()[i].size()]);
				contoursLen += mesh.property(para_l, mesh.edge_handle(he));
			}
		}

		double mean_inner_curv = 2 * M_PI / contoursLen;

		for (auto v : seg.inner())
		{
			mesh.property(para_TH, v) = 2 * mean_inner_curv * mesh.property(mix_area, v);
		}

		for (int i = 0; i < seg.vContours().size(); i++)
		{
			if (seg.bContours()[i])
			{
				for (auto v : seg.vContours()[i])
				{
					mesh.property(para_TH, v) = 0;
				}
			}
			else
			{
				for (int j = 0; j < seg.vContours()[i].size(); j++)
				{
					mesh.property(para_TH, seg.vContours()[i][j]) = 0;
				}
			}
		}
	}

	void setConstantHMeanCurvature(const opParas &op)
	{

	}

	void setConstantGMeanCurvature(const opParas &op)
	{
		for (auto v : domain.variables())
		{
			mesh.property(para_TH, v) = op.curv_mean*mesh.property(mix_area, v);
		}
	}

	void setPolygonMeanCurvature(const opParas &op)
	{
		for (auto v : seg.inner())
		{
			mesh.property(para_TH, v) = 0;
		}
	}

	void setSpherePatchMeanCurvature(const opParas &op)
	{
		double areaSum = 0.0;
		for (auto v : seg.inner())
		{
			areaSum += mesh.property(mix_area, v);
		}
		int corner = (seg.contour(0).size() + 1) / op.polygonNum;
		for (int i = 0; i < seg.contour(0).size(); i++)
		{
			if (i%corner != 0)
			{
				areaSum += mesh.property(mix_area, seg.contour(0)[i]);
			}
		}

		double gauss_inner_curv = op.curv_patch / areaSum;
		double mean_inner_curv = sqrt(gauss_inner_curv);

		for (auto v : seg.inner())
		{
			mesh.property(para_TH, v) = 2 * mean_inner_curv * mesh.property(mix_area, v);
		}
	}

	void setNonIterMeanCurvature(const opParas &op)
	{
	}

	void projectSphere(const opParas &op)
	{
		mesh.update_normals();

		Eigen::Matrix3d nnT;
		Eigen::Vector3d ct;
		Eigen::Vector3d cv;

		Eigen::Vector3d nm(0, 0, 0);
		double dm = 0;
		for (auto v : seg.inner())
		{
			nm += mesh.normal(v);
			dm += mesh.point(v).dot(mesh.normal(v));
		}
		nm /= seg.inner().size();
		dm /= seg.inner().size();

		nnT.setZero();
		ct.setZero();
		for (auto v : seg.inner())
		{
			cv = nm - mesh.normal(v);
			nnT += cv * cv.transpose();
			ct += cv * (dm - mesh.point(v).dot(mesh.normal(v)));
		}

		ct = nnT.inverse()*ct;

		for (auto c : seg.vContours())
		{
			Eigen::Vector3d cct0(0, 0, 0);
			Eigen::Vector3d cct(0, 0, 0);
			for (auto v : c)
			{
				cct += mesh.point(v);
				cct0 += mesh.property(point0, v);
			}
			cct0 /= c.size();
			cct /= c.size();

			double r = 0;
			for (auto v : c)
			{
				r += (mesh.point(v) - cct).norm();
			}
			r /= c.size();

			Eigen::Vector3d cv = cct - ct;
			cv.normalize();
			for (auto v : c)
			{
				Eigen::Vector3d ntmp = mesh.property(point0, v) - cct0;
				Eigen::Vector3d vtmp = ntmp.dot(cv)*cv;
				ntmp -= vtmp;
				vtmp = cct + vtmp + ntmp / ntmp.norm()*r - ct;
				ntmp = mesh.point(v) - ct;
				mesh.property(back_p, v) = ct + vtmp / vtmp.norm() * ntmp.norm();
			}
		}
	}

	void projectCylinder(const opParas &op)
	{
		mesh.update_normals();
		
		Eigen::Matrix3d nnT;
		Eigen::Vector3d ct;
		Eigen::Vector3d cv;

		Eigen::Vector3d nm(0, 0, 0);
		double dm = 0;
		for (auto v : seg.inner())
		{
			nm += mesh.normal(v);
			dm += mesh.point(v).dot(mesh.normal(v));
		}
		nm /= seg.inner().size();
		dm /= seg.inner().size();

		nnT.setZero();
		ct.setZero();
		for (auto v : seg.inner())
		{
			cv = nm - mesh.normal(v);
			nnT += cv * cv.transpose();
			ct += cv * (dm - mesh.point(v).dot(mesh.normal(v)));
		}

		ct = nnT.inverse()*ct;

		nnT.setZero();
		for (auto v : seg.inner())
		{
			cv = mesh.normal(v);
			nnT = nnT + cv * cv.transpose();
		}

		Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(nnT);
		cv = solver.eigenvectors().col(0);
		cv.normalize();

		for (auto c : seg.vContours())
		{
			Eigen::Vector3d cct0(0, 0, 0);
			Eigen::Vector3d cct(0, 0, 0);
			for (auto v : c)
			{
				cct += mesh.point(v);
				cct0 += mesh.property(point0, v);
			}
			cct0 /= c.size();
			cct /= c.size();

			double r = 0;
			for (auto v : c)
			{
				r += (mesh.point(v) - cct).norm();
			}
			r /= c.size();

			Eigen::Vector3d shift = ct + (cct0 - ct).dot(cv)*cv;

			for (auto v : c)
			{
				Eigen::Vector3d ntmp = mesh.property(point0, v) - cct0;
				Eigen::Vector3d vtmp = ntmp.dot(cv)*cv;
				ntmp -= vtmp;
				mesh.property(back_p, v) = shift + vtmp + ntmp / ntmp.norm()*r;
			}
		}
	}

	//-------------------------------------------------------
	// ***** ---------- Local Calabi Flow ---------- ***** //
	void curve_progressive(Eigen::VectorXd &differ, double maxErr)
	{
		assert(differ.size() == domain.n_variables());

		double clip = 0.3;
		double scale;
		if (maxErr > clip) scale = clip / maxErr;
		else scale = 1;

		for (auto v : domain.variables())
		{
			differ[domain.index(v)] = (mesh.property(para_TG, v) - mesh.property(para_G, v)) *scale;
		}
	}

	void curve_progressive(Eigen::VectorXd &differ, int step = 4)
	{
		assert(differ.size() == domain.n_variables());

		for (auto v : domain.variables())
		{
			differ[domain.index(v)] = (mesh.property(para_TG, v) - mesh.property(para_G, v)) / step;
		}
	}
	void curve_deviation(Eigen::VectorXd &deviat)
	{
		assert(deviat.size() == domain.n_variables());

		for (auto v : domain.variables())
		{
			deviat[domain.index(v)] = mesh.property(para_TG, v) - mesh.property(para_G, v);
		}
	}
	Eigen::VectorXd curve_progressive()
	{
		Eigen::VectorXd rt(domain.n_variables());
		curve_progressive(rt);
		return rt;
	}
	Eigen::VectorXd curve_deviation()
	{
		Eigen::VectorXd rt(domain.n_variables());
		curve_deviation(rt);
		return rt;
	}

	// request update op_theta before
	void domainEdgeFlip()
	{
		std::vector<Mesh::EdgeHandle> delaunay_stack;
		delaunay_stack.reserve(domain.edges().size());

		bool indicate = false;
		int count = 0;
		do
		{
			for (auto e : domain.edges())
			{
				if (domain.onBoundary(e)) continue;
				delaunay_stack.push_back(e);
				mesh.property(flip_indicate, e) = true;
			}
			indicate = false;

			OpenMesh::HalfedgeHandle he_t;
			OpenMesh::EdgeHandle e_t, e_n;
			while (!delaunay_stack.empty())
			{
				e_t = delaunay_stack.back();
				delaunay_stack.pop_back();
				mesh.property(flip_indicate, e_t) = false;

				if (mesh.property(op_theta, mesh.halfedge_handle(e_t, 0))
					+ mesh.property(op_theta, mesh.halfedge_handle(e_t, 1)) > THETAFLIP)
				{
					if (!mesh.is_flip_ok(e_t)) continue;
					if (mesh.status(e_t).feature()) continue;				
					indicate = attemptFlip(e_t, delaunay_stack);
				}
			}
			count++;
		} while (indicate && count < 5);
	}

	void updateDomainEdgeLength()
	{
		for (auto e : domain.edges())
		{
			mesh.property(para_l, e) = calcELength(e);
		}
	}

	// request update edge length before
	bool updateDomainTheta(double thetaThres)
	{
		bool inQuality = true;
		for (auto he : domain.halfedges())
		{
			inQuality = calcHeOppositeTheta(he, mesh.property(op_theta, he), thetaThres) && inQuality;
			mesh.property(cot_theta, he) = 1.0 / tan(mesh.property(op_theta, he));
		}
		return inQuality;
	}

	// request update op_theta and mix_area before
	void updateDomainGaussCurvature()
	{
		double gauss_curv;
		OpenMesh::HalfedgeHandle he_t;

		for (auto v : domain.vertices())
		{
			mesh.property(para_G, v) = domain.onBoundary(v) ? M_PI : 2 * M_PI;
		}

		for (auto he : domain.halfedges())
		{
			mesh.property(para_G, mesh.from_vertex_handle(he)) -= mesh.property(op_theta, mesh.next_halfedge_handle(he));
		}
	}

	void updateDomainMeanCurvature()
	{
		OpenMesh::VertexHandle fromV, toV;
		double cot_w;

		for (auto v : domain.vertices())
		{
			mesh.property(slack_n, v).setZero();
		}

		for (auto he : domain.halfedges())
		{
			cot_w = mesh.property(cot_theta, he);
			fromV = mesh.from_vertex_handle(he);
			toV = mesh.to_vertex_handle(he);
			mesh.property(slack_n, fromV) += cot_w * (mesh.point(fromV) - mesh.point(toV));
			mesh.property(slack_n, toV) += cot_w * (mesh.point(toV) - mesh.point(fromV));
		}

		for (auto v : domain.vertices())
		{
			mesh.property(para_H, v) = mesh.property(slack_n, v).norm();
		}
	}

	// request update op_theta before
	void updateDomainMixArea()
	{
		OpenMesh::HalfedgeHandle he[3];
		double helen[3];
		double theta[3];
		double varea[3];
		int maxId;
		double sarea;
		
		for (auto v : domain.vertices())
		{
			mesh.property(mix_area, v) = 0;
		}

		for (auto f : domain.faces())
		{
			he[0] = mesh.halfedge_handle(f);
			helen[0] = mesh.property(para_l, mesh.edge_handle(he[0]));
			theta[0] = mesh.property(op_theta, he[0]);
			varea[0] = helen[0] * helen[0] * mesh.property(cot_theta, he[0]) / 8;

			he[1] = mesh.next_halfedge_handle(he[0]);
			helen[1] = mesh.property(para_l, mesh.edge_handle(he[1]));
			theta[1] = mesh.property(op_theta, he[1]);
			varea[1] = helen[1] * helen[1] * mesh.property(cot_theta, he[1]) / 8;

			he[2] = mesh.next_halfedge_handle(he[1]);
			helen[2] = mesh.property(para_l, mesh.edge_handle(he[2]));
			theta[2] = mesh.property(op_theta, he[2]);
			varea[2] = helen[2] * helen[2] * mesh.property(cot_theta, he[2]) / 8;

			maxId = 0;
			if (helen[maxId] < helen[1]) maxId = 1;
			if (helen[maxId] < helen[2]) maxId = 2;

			if (theta[maxId] < M_PI / 2)
			{
				mesh.property(mix_area, mesh.to_vertex_handle(he[0])) += varea[0] + varea[1];
				mesh.property(mix_area, mesh.to_vertex_handle(he[1])) += varea[1] + varea[2];
				mesh.property(mix_area, mesh.to_vertex_handle(he[2])) += varea[2] + varea[0];
			}
			else
			{
				sarea = (helen[0] * helen[1] * sin(theta[2]) + helen[1] * helen[2] * sin(theta[0])
					+ helen[2] * helen[0] * sin(theta[1])) / 12;

				mesh.property(mix_area, mesh.to_vertex_handle(he[0])) += sarea;
				mesh.property(mix_area, mesh.to_vertex_handle(he[1])) += sarea;
				mesh.property(mix_area, mesh.to_vertex_handle(he[2])) += sarea;

				mesh.property(mix_area, mesh.from_vertex_handle(he[maxId])) -= sarea / 2;
				mesh.property(mix_area, mesh.to_vertex_handle(he[maxId])) -= sarea / 2;
			}
		}
	}

	void updateDomainMesh(double thetaThres)
	{
		std::cout << "Edge Flipping:";
		updateDomainEdgeLength();
		verify(updateDomainTheta(thetaThres));
		domainEdgeFlip();
		std::cout << "......    ";
		updateDomainMixArea();
		updateDomainGaussCurvature();
	}

	void assembleLaplace()
	{
		int fromIndex, toIndex;
		OpenMesh::HalfedgeHandle he_t;
		double l_ii, l_ij;

		dual_lap_trip.clear();
		for (auto v : domain.variables())
		{
			l_ii = 0.0;
			toIndex = domain.index(v);
			assert(toIndex >= 0 && toIndex < domain.n_variables());
			for (auto vih : mesh.vih_range(v))
			{
				fromIndex = domain.index(mesh.from_vertex_handle(vih));
				if (fromIndex < -1) continue;

				l_ij = 0;

				if (domain.inDomain(vih))
					l_ij -= mesh.property(cot_theta, vih);
				he_t = mesh.opposite_halfedge_handle(vih);
				if (domain.inDomain(he_t))
					l_ij -= mesh.property(cot_theta, he_t);

				l_ii -= l_ij;

				if (fromIndex < 0) continue;
				assert(fromIndex >= 0 && fromIndex < domain.n_variables());
				dual_lap_trip.push_back(Eigen::Triplet<double>(
					toIndex, fromIndex, l_ij));
			}
			dual_lap_trip.push_back(Eigen::Triplet<double>(
				toIndex, toIndex, l_ii));
		}
		local_L.resize(domain.n_variables(), domain.n_variables());
		local_L.setFromTriplets(dual_lap_trip.begin(), dual_lap_trip.end());
	}

	void assembleLaplaceCol()
	{
		int fromIndex, toIndex;
		OpenMesh::HalfedgeHandle he_t;
		double l_ii, l_ij;

		dual_lap_trip.clear();
		for (int i = 0; i < domain.vertices().size(); i++)
		{
			l_ii = 0.0;
			if (domain.onBoundary(domain.vertices()[i])) continue;
			toIndex = domain.index(domain.vertices()[i]);
			for (auto vih : mesh.vih_range(domain.vertices()[i]))
			{
				fromIndex = domain.index(mesh.from_vertex_handle(vih));
				if (fromIndex < -1) continue;

				l_ij = 0;

				if (domain.inDomain(vih))
					l_ij -= mesh.property(cot_theta, vih);
				he_t = mesh.opposite_halfedge_handle(vih);
				if (domain.inDomain(he_t))
					l_ij -= mesh.property(cot_theta, he_t);

				l_ii -= l_ij;

				if (fromIndex < 0) continue;
				assert(fromIndex >= 0 && fromIndex < domain.n_variables());
				dual_lap_trip.push_back(Eigen::Triplet<double>(
					i, fromIndex, l_ij));
			}
			if (toIndex < 0) continue;
			dual_lap_trip.push_back(Eigen::Triplet<double>(
				i, toIndex, l_ii));
		}
		local_C.resize(domain.vertices().size(), domain.n_variables());
		local_C.setFromTriplets(dual_lap_trip.begin(), dual_lap_trip.end());
	}

	double calcQuasStep(const Eigen::VectorXd &du, double thetaThres, double shrink)
	{
		double quasStep = 1.0f;
		bool triangulable = false;

		for (auto v : domain.variables())
		{
			mesh.property(back_u, v) = mesh.property(para_u, v);
		}

		int count = 0;
		do
		{
			for (auto v : domain.variables())
			{
				mesh.property(para_u, v) = mesh.property(back_u, v)
					+ quasStep * du(domain.index(v));
			}
			updateDomainEdgeLength();
			triangulable = updateDomainTheta(thetaThres);
			if (!triangulable) quasStep *= shrink;
			count++;
		} while (!triangulable && count < 10);

		for (auto v : domain.variables())
		{
			mesh.property(para_u, v) = mesh.property(back_u, v);
		}
		if (count >= 10) return 0;
		return quasStep;
	}

	void quasiNewton(const opParas &op)
	{
		int count, upcount;
		double last_err, curr_err;
		Eigen::VectorXd du(domain.n_variables());

		double source_areaSum = 0;
		for (auto v : domain.vertices())
		{
			source_areaSum += mesh.property(mix_area, v);
		}

		set_gauss_func[op.kind](op);
		curr_err = curve_deviation().lpNorm<Eigen::Infinity>();
		last_err = curr_err;
		std::cout << "Initial Curvatre Error: " << curr_err << std::endl;
		std::cout << "----------------------------------------------\n";

		count = 0;
		upcount = 0;
		std::cout << "Quasi Newton phase:\n";
		while (curr_err >= op.curv_err && count < op.quas_cap && upcount < 2)
		{
			if (curr_err > last_err) upcount++;
			else upcount = 0;
			last_err = curr_err;

			assembleLaplace();
			llt.compute(local_L);
			assert(llt.info() == Eigen::Success);

			curve_progressive(du, curr_err);
			//curve_progressive(du, (op.quas_cap - count) / 5 + 1);
			du = llt.solve(du);
			du = du * calcQuasStep(du, THETATHRES, op.shrink);

			for (auto v : domain.variables())
			{
				mesh.property(para_u, v) += du(domain.index(v));
			}
			updateDomainMesh(THETATHRES);

			set_gauss_func[op.kind](op);
			curr_err = curve_deviation().lpNorm<Eigen::Infinity>();
			std::cout << "the " << count << "th iter: " << curr_err << std::endl;
			count++;
		}

		double areaSum = 0;
		for (auto v : domain.vertices())
		{
			areaSum += mesh.property(mix_area, v);
		}

		double minusU = (log(areaSum) - log(source_areaSum)) / 4;
		for (auto v : domain.vertices())
		{
			mesh.property(para_u, v) -= minusU;
		}
		updateDomainMesh(THETATHRES);

		set_gauss_func[op.kind](op);
		curr_err = curve_deviation().lpNorm<Eigen::Infinity>();
		std::cout << "the " << count << "th iter: " << curr_err << std::endl;
	}

	void shapeOpt(const opParas &op)
	{
		assembleLaplace();
		assembleLaplaceCol();
		llt.compute(op.lambdaH*local_C.transpose()*local_C + local_L);
		assert(llt.info() == Eigen::Success);

		int count = 0;
		double cur_err = 0.0;
		Eigen::MatrixXd Ln(domain.n_variables(), 3);
		Eigen::MatrixXd Ed(domain.n_variables(), 3);
		Eigen::MatrixXd pos(domain.n_variables(), 3);

		for (auto v : mesh.vertices())
		{
			mesh.property(point0, v) = mesh.point(v);
		}

		for (auto v : domain.variables())
		{
			pos.row(domain.index(v)) = mesh.point(v);
		}

		AndersonAcceleration aa;
		aa.init(op.andersonM, domain.n_variables() * 3, pos.data());

		//mesh.request_face_normals();
		//mesh.update_normals();

		set_mean_func[op.kind](op);
		cur_err = op.lambdaH * updateLslackN(Ln) + updateEslackD(Ed);
		std::cout << "Deform Loss:......    " << "the " << count << "th iter: " << cur_err << std::endl;
		do
		{
			pos = llt.solve(op.lambdaH*Ln + Ed);
			pos = Eigen::Map<const Eigen::MatrixXd>(aa.compute(pos.data()).data(), domain.n_variables(), 3);
			for (auto v : domain.variables())
			{
				mesh.set_point(v, pos.row(domain.index(v)));
			}
			//updateDomainMeanCurvature();
			//mesh.update_normals();
			count++;

			aa.replace(pos.data());
			cur_err = op.lambdaH * updateLslackN(Ln) + updateEslackD(Ed);
			std::cout << "Deform Loss:......    " << "the " << count << "th iter: " << cur_err << std::endl;
		} while (count < op.deform_cap);
	}

	//-------------------------------------------------------
	// ***** ---------- global position opt ---------- ***** //

	void initLocalPos()
	{
		// AffineMap transform
		Eigen::Vector3d p, p0, ct(0, 0, 0), ct0(0, 0, 0);
		double pts_w, sum_w1, sum_w2;

		sum_w1 = 0;
		sum_w2 = 0;
		for (int i = 0; i < seg.n_contours(); i++)
		{
			if (!seg.bContours()[i]) continue;
			for (int j=0; j < seg.contour(i).size(); j++)
			{
				pts_w = (mesh.property(point0, seg.contour(i)[j]) -
					mesh.property(point0, seg.contour(i)[(j + 1) % seg.contour(i).size()])).norm();
				pts_w += (mesh.property(point0, seg.contour(i)[j]) -
					mesh.property(point0, seg.contour(i)[(j + seg.contour(i).size() - 1) % seg.contour(i).size()])).norm();
				ct += pts_w * mesh.property(point0, seg.contour(i)[j]);
				sum_w1 += pts_w;

				pts_w = (mesh.point(seg.contour(i)[j]) -
					mesh.point(seg.contour(i)[(j + 1) % seg.contour(i).size()])).norm();
				pts_w += (mesh.point(seg.contour(i)[j]) -
					mesh.point(seg.contour(i)[(j + seg.contour(i).size() - 1) % seg.contour(i).size()])).norm();
				ct0 += pts_w * mesh.point(seg.contour(i)[j]);
				sum_w2 += pts_w;
			}
		}
		if (sum_w1 > 0 && sum_w2 > 0)
		{
			ct /= sum_w1;
			ct0 /= sum_w2;
		}
		
		ASAP asap;
		for (int i = 0; i < seg.n_contours(); i++)
		{
			if (!seg.bContours()[i]) continue;
			for (auto v : seg.contour(i))
			{
				p = mesh.property(point0, v) - ct;
				p0 = mesh.point(v) - ct0;
				asap.addCouples(p0, p, 1); // may have problem
			}
		}
		AffineMap tmp = asap.solve();

		tmp.T = ct - tmp.map(ct0);

		for (auto v : domain.vertices())
		{
			mesh.set_point(v, tmp.map(mesh.point(v)));
		}
	}

	void initRemainPos()
	{
		Eigen::Vector3d p, p0, ct(0, 0, 0), ct0(0, 0, 0);
		double pts_w, sum_w1, sum_w2;
		ASAP asap;
		AffineMap tmp;

		for (int i = 0; i < seg.n_outers(); i++)
		{
			ct.setZero();
			ct0.setZero();
			sum_w1 = 0;
			sum_w2 = 0;

			for (int j = 0; j < seg.n_contours(); j++)
			{
				if (!seg.isRelated(i, j)) continue;
				for (int k = 0; k < seg.contour(j).size(); k++)
				{
					pts_w = (mesh.property(point0, seg.contour(j)[k]) -
						mesh.property(point0, seg.contour(j)[(k + 1) % seg.contour(j).size()])).norm();
					pts_w += (mesh.property(point0, seg.contour(j)[k]) -
						mesh.property(point0, seg.contour(j)[(k + seg.contour(j).size() - 1) % seg.contour(j).size()])).norm();
					ct0 += pts_w * mesh.property(point0, seg.contour(j)[k]);
					sum_w1 += pts_w;

					pts_w = (mesh.property(back_p, seg.contour(j)[k]) -
						mesh.property(back_p, seg.contour(j)[(k + 1) % seg.contour(j).size()])).norm();
					pts_w += (mesh.property(back_p, seg.contour(j)[k]) -
						mesh.property(back_p, seg.contour(j)[(k + seg.contour(j).size() - 1) % seg.contour(j).size()])).norm();
					ct += pts_w * mesh.property(back_p, seg.contour(j)[k]);
					sum_w2 += pts_w;
				}
			}
			if (sum_w1 > 0 && sum_w2 > 0)
			{
				ct0 /= sum_w1;
				ct /= sum_w2;
			}

			asap.clear();
			for (int j = 0; j < seg.n_contours(); j++)
			{
				if (!seg.isRelated(i, j)) continue;
				for (auto v : seg.contour(j))
				{
					p0 = mesh.property(point0, v) - ct0;
					p = mesh.property(back_p, v) - ct;
					asap.addCouples(p0, p, 1); // may have problem
				}
			}
			tmp = asap.solve();
			tmp.T = ct - tmp.map(ct0);

			for (auto v : seg.outer(i))
			{
				mesh.property(point0, v) = tmp.map(mesh.property(point0, v));
			}
		}
	}

	void initASAPS()
	{
		for (auto v : domain.vertices())
		{
			if (!domain.onBoundary(v)) continue;
			double l = 0, l0 = 0;
			for (auto vih : mesh.vih_range(v))
			{
				if (domain.onBoundary(vih))
				{
					l += (mesh.point(mesh.from_vertex_handle(vih)) - mesh.point(v)).norm();
					l0 += (mesh.property(point0, mesh.from_vertex_handle(vih)) -
						mesh.property(point0, v)).norm();
				}
			}
			mesh.property(affine_M, v).S = l / l0;
		}
	}

	void laplacePoint()
	{
		int fromIndex, toIndex;
		OpenMesh::HalfedgeHandle he_t;
		double l_ii, l_ij;

		dual_lap_trip.clear();
		for (auto v : domain.variables())
		{
			l_ii = 0.0;
			toIndex = domain.index(v);
			assert(toIndex >= 0 && toIndex < domain.n_variables());
			for (auto vih : mesh.vih_range(v))
			{
				fromIndex = domain.index(mesh.from_vertex_handle(vih));
				if (fromIndex < -1) continue;

				l_ij = 0;

				if (domain.inDomain(vih))
					l_ij -= 1;
				he_t = mesh.opposite_halfedge_handle(vih);
				if (domain.inDomain(he_t))
					l_ij -= 1;

				l_ii -= l_ij;

				if (fromIndex < 0) continue;
				assert(fromIndex >= 0 && fromIndex < domain.n_variables());
				dual_lap_trip.push_back(Eigen::Triplet<double>(
					toIndex, fromIndex, l_ij));
			}
			dual_lap_trip.push_back(Eigen::Triplet<double>(
				toIndex, toIndex, l_ii));
		}
		local_L.resize(domain.n_variables(), domain.n_variables());
		local_L.setFromTriplets(dual_lap_trip.begin(), dual_lap_trip.end());
	}

	void laplacePoint(int cap)
	{
		OpenMesh::VPropHandleT<Eigen::Vector3d> Lp;
		mesh.add_property(Lp);
		OpenMesh::VertexHandle fromV;

		for (int i = 0; i < cap; i++)
		{
			for (auto v : domain.vertices())
			{
				if (domain.onBoundary(v)) continue;
				int valen = mesh.valence(v);
				Eigen::Vector3d dp(0, 0, 0);
				Eigen::Vector3d dn(0, 0, 0);
				mesh.property(Lp, v).setZero();
				for (auto vih : mesh.vih_range(v))
				{
					fromV = mesh.from_vertex_handle(vih);
					dp += (mesh.point(fromV) - mesh.point(v)) / valen;

					double w = 0;
					w += mesh.property(cot_theta, vih);
					auto he_t = mesh.opposite_halfedge_handle(vih);
					w += mesh.property(cot_theta, he_t);
					dn += w * (mesh.point(fromV) - mesh.point(v));
				}

				dn.normalize();
				mesh.property(Lp, v) = dp; // -(dp | dn) * dn;
			}

			for (auto v : domain.vertices())
			{
				if (domain.onBoundary(v)) continue;
				mesh.set_point(v, mesh.point(v) + mesh.property(Lp, v));
			}

			/*for (auto e : domain.edges())
			{
				mesh.property(para_l, e) = mesh.calc_edge_length(e);
			}
			updateDomainEdgeLength();
			updateDomainTheta(THETATHRES);*/
		}

		mesh.remove_property(Lp);
	}

	void laplaceASAPS(int cap)
	{
		OpenMesh::VPropHandleT<double> affine_S;
		mesh.add_property(affine_S);
		OpenMesh::VertexHandle fromV;

		for (int i = 0; i < cap; i++)
		{
			for (auto v : domain.vertices())
			{
				if (domain.onBoundary(v)) continue;
				int valen = mesh.valence(v);
				mesh.property(affine_S, v) = 0;
				for (auto vih : mesh.vih_range(v))
				{
					fromV = mesh.from_vertex_handle(vih);
					mesh.property(affine_S, v) += mesh.property(affine_M, fromV).S;
				}

				mesh.property(affine_S, v) /= valen;
			}

			for (auto v : domain.vertices())
			{
				if (domain.onBoundary(v)) continue;
				mesh.property(affine_M, v).S = mesh.property(affine_S, v);
			}
		}

		mesh.remove_property(affine_S);
	}

	void assembleOpt(const opParas &op, bool allowscale)
	{
		assembleLaplace();
		llt.compute(local_L);
		assert(llt.info() == Eigen::Success);

		Eigen::MatrixXd b(domain.n_variables(), 3);
		Eigen::MatrixXd pos(domain.n_variables(), 3);


		if (!allowscale)
		{
			initASAPS();
			laplaceASAPS(5);
		}

		/*for (auto v : domain.vertices())
		{
			if(!domain.onBoundary(v))
				std::cout << mesh.property(AffineMap, v).S << std::endl;
		}*/


		for (auto v : domain.variables())
		{
			pos.row(domain.index(v)) = mesh.point(v);
		}

		//AndersonAcceleration aa;
		//aa.init(8, var_vert.size() * 3, pos.data());

		int count = 0;
		do
		{
			updateASAPS(allowscale);
			updateASAPb(b);
			pos = llt.solve(b);
			//pos = Eigen::Map<const Eigen::MatrixXd>(aa.compute(pos.data()).data(), var_vert.size(), 3);
			for (auto v : domain.variables())
			{
				mesh.set_point(v, pos.row(domain.index(v)));
			}
			count++;
			//aa.replace(pos.data());
		} while (count < op.ASAP_cap);
	}

	void initProperty()
	{
		for (auto v : mesh.vertices())
		{
			mesh.property(para_u, v) = 0;
			mesh.property(slack_n, v).setZero();
			mesh.property(affine_M, v) = AffineMap();
			mesh.property(para_TH, v) = 0;
		}

		for (auto e : mesh.edges())
		{
			mesh.property(para_l, e) = mesh.calc_edge_vector(e).norm();
		}

		for (auto he : mesh.halfedges())
		{
			mesh.property(slack_d, he).setZero();
			if (mesh.is_boundary(he)) continue;
			verify(calcHeOppositeTheta(he, mesh.property(op_theta, he), THETATHRES));
			mesh.property(cot_theta, he) = 1.0 / mesh.property(op_theta, he);
		}
	}

	void addProperty()
	{
		mesh.add_property(para_u);
		mesh.add_property(para_l);
		mesh.add_property(op_theta);
		mesh.add_property(cot_theta);
		mesh.add_property(para_G);
		mesh.add_property(para_TG);
		mesh.add_property(para_H);
		mesh.add_property(para_TH);
		mesh.add_property(back_u);
		mesh.add_property(back_p);
		mesh.add_property(point0);
		mesh.add_property(mix_area);
		mesh.add_property(slack_d);
		mesh.add_property(slack_n);
		mesh.add_property(affine_M);
		mesh.add_property(flip_indicate);

		mesh.add_property(kmax);
		mesh.add_property(kmin);
	}

	void removeProperty()
	{
		mesh.remove_property(para_u);
		mesh.remove_property(para_l);
		mesh.remove_property(op_theta);
		mesh.remove_property(cot_theta);
		mesh.remove_property(para_G);
		mesh.remove_property(para_TG);
		mesh.remove_property(para_H);
		mesh.remove_property(para_TH);
		mesh.remove_property(back_u);
		mesh.remove_property(back_p);
		mesh.remove_property(point0);
		mesh.remove_property(mix_area);
		mesh.remove_property(slack_d);
		mesh.remove_property(slack_n);
		mesh.remove_property(affine_M);
		mesh.remove_property(flip_indicate);

		mesh.remove_property(kmax);
		mesh.remove_property(kmin);
	}

	//-------------------------------------------------------

	double updateLslackN(Eigen::MatrixXd &Ln)
	{
		Ln.setZero();
		double ms_err = 0.0;
		OpenMesh::VertexHandle fromV, toV;
		double cot_w;

		for (auto v : domain.vertices())
		{
			mesh.property(slack_n, v).setZero();
		}

		for (auto he : domain.halfedges())
		{
			cot_w = mesh.property(cot_theta, he);
			fromV = mesh.from_vertex_handle(he);
			toV = mesh.to_vertex_handle(he);
			mesh.property(slack_n, fromV) += cot_w * (mesh.point(fromV) - mesh.point(toV));
			mesh.property(slack_n, toV) += cot_w * (mesh.point(toV) - mesh.point(fromV));
		}

		for (auto v : domain.vertices())
		{
			if (domain.onBoundary(v))
			{
				mesh.property(slack_n, v).setZero();
				continue;
			}

			double tmp = massSpring(mesh.property(slack_n, v), mesh.property(slack_n, v), mesh.property(para_TH, v));
			ms_err = (tmp > ms_err) ? tmp : ms_err;
		}

		for (auto v : domain.anchors())
		{
			for (auto vih : mesh.vih_range(v))
			{
				if (!domain.inDomain(vih)) continue;
				cot_w = mesh.property(cot_theta, vih);
				fromV = mesh.from_vertex_handle(vih);
				mesh.property(slack_n, fromV) += cot_w * mesh.point(v);
				mesh.property(slack_n, v) -= cot_w * mesh.point(v);
			}

			for (auto voh : mesh.voh_range(v))
			{
				if (!domain.inDomain(voh)) continue;
				cot_w = mesh.property(cot_theta, voh);
				toV = mesh.to_vertex_handle(voh);
				mesh.property(slack_n, toV) += cot_w * mesh.point(v);
				mesh.property(slack_n, v) -= cot_w * mesh.point(v);
			}
		}

		for (auto he : domain.halfedges())
		{
			fromV = mesh.from_vertex_handle(he);
			toV = mesh.to_vertex_handle(he);
			mesh.property(slack_d, he) = mesh.property(slack_n, toV) - mesh.property(slack_n, fromV);
		}

		for (auto v : domain.vertices())
		{
			mesh.property(slack_n, v).setZero();
		}

		for (auto he : domain.halfedges())
		{
			cot_w = mesh.property(cot_theta, he);
			fromV = mesh.from_vertex_handle(he);
			toV = mesh.to_vertex_handle(he);
			mesh.property(slack_n, fromV) -= cot_w * mesh.property(slack_d, he);
			mesh.property(slack_n, toV) += cot_w * mesh.property(slack_d, he);
		}

		for (auto v : domain.variables())
		{
			Ln.row(domain.index(v)) = mesh.property(slack_n, v);
		}

		return ms_err;
	}

	double updateEslackUV(Eigen::MatrixXd &Ed, OpenMesh::HPropHandleT<Eigen::Vector2d> &slack_uv)
	{
		Ed.setZero();
		OpenMesh::HalfedgeHandle he_t;
		OpenMesh::VertexHandle fromV;

		double ms_err = 0.0;

		for (auto e : domain.edges())
		{
			for (int i = 0; i < 2; i++)
			{
				he_t = mesh.halfedge_handle(e, i);
				Eigen::Vector2d Ep = mesh.texcoord2D(mesh.to_vertex_handle(he_t)) 
					- mesh.texcoord2D(mesh.from_vertex_handle(he_t));
				double tmp = massSpring(Ep, mesh.property(slack_uv, he_t), mesh.property(para_l, e));
				ms_err = (tmp > ms_err) ? tmp : ms_err;
			}
		}

		for (auto v : domain.variables())
		{
			Eigen::Vector2d tmp(0, 0);
			for (auto vih : mesh.vih_range(v))
			{
				fromV = mesh.from_vertex_handle(vih);
				double cot_w = 0;

				if (domain.inDomain(vih))
					cot_w += mesh.property(cot_theta, vih);
				he_t = mesh.opposite_halfedge_handle(vih);
				if (domain.inDomain(he_t))
					cot_w += mesh.property(cot_theta, he_t);

				tmp += cot_w * mesh.property(slack_uv, vih);
				if (domain.isAnchor(fromV))
					tmp += cot_w * mesh.texcoord2D(fromV);
			}
			Ed.row(domain.index(v)) = tmp;
		}

		return ms_err;
	}

	double updateEslackD(Eigen::MatrixXd &Ed)
	{
		Ed.setZero();
		OpenMesh::VertexHandle fromV, toV;
		double cot_w;

		double ms_err = 0.0;

		for (auto he : domain.halfedges())
		{
			mesh.property(slack_d, he) = mesh.point(mesh.to_vertex_handle(he)) - mesh.point(mesh.from_vertex_handle(he));
			double tmp = massSpring(mesh.property(slack_d, he), mesh.property(slack_d, he),
				mesh.property(para_l, mesh.edge_handle(he)));
			ms_err = (tmp > ms_err) ? tmp : ms_err;
		}

		for (auto v : domain.vertices())
		{
			mesh.property(slack_n, v).setZero();
		}

		for (auto he : domain.halfedges())
		{
			cot_w = mesh.property(cot_theta, he);
			fromV = mesh.from_vertex_handle(he);
			toV = mesh.to_vertex_handle(he);
			mesh.property(slack_n, fromV) -= cot_w * mesh.property(slack_d, he);
			mesh.property(slack_n, toV) += cot_w * mesh.property(slack_d, he);
		}

		for (auto v : domain.anchors())
		{
			for (auto vih : mesh.vih_range(v))
			{
				if (!domain.inDomain(vih)) continue;
				cot_w = mesh.property(cot_theta, vih);
				fromV = mesh.from_vertex_handle(vih);
				mesh.property(slack_n, fromV) += cot_w * mesh.point(v);
			}

			for (auto voh : mesh.voh_range(v))
			{
				if (!domain.inDomain(voh)) continue;
				cot_w = mesh.property(cot_theta, voh);
				toV = mesh.to_vertex_handle(voh);
				mesh.property(slack_n, toV) += cot_w * mesh.point(v);
			}
		}

		for (auto v : domain.variables())
		{
			Ed.row(domain.index(v)) = mesh.property(slack_n, v);
		}

		return ms_err;
	}

	void updateASAPb(Eigen::MatrixXd &b)
	{
		b.setZero();
		double l_ij;
		OpenMesh::HalfedgeHandle he_t;
		OpenMesh::VertexHandle fromV;

		for (auto v : domain.variables())
		{
			Eigen::Vector3d tmp(0, 0, 0);
			Eigen::Vector3d sum(0, 0, 0);
			for (auto vih : mesh.vih_range(v))
			{
				l_ij = 0;
				fromV = mesh.from_vertex_handle(vih);
				tmp = mesh.property(point0, v) - mesh.property(point0, fromV);
	
				if (domain.inDomain(vih))
					l_ij += mesh.property(cot_theta, vih);
				he_t = mesh.opposite_halfedge_handle(vih);
				if (domain.inDomain(he_t))
					l_ij += mesh.property(cot_theta, he_t);

				const AffineMap &toaf = mesh.property(affine_M, v);
				const AffineMap &fromaf = mesh.property(affine_M, fromV);
				sum += l_ij / 2 * (toaf.S*toaf.R + fromaf.S * fromaf.R)*tmp;

				if (domain.isAnchor(fromV))
					sum += l_ij * mesh.point(fromV);
			}
			b.row(domain.index(v)) = sum;
		}
	}

	void updateASAPS(bool allowscale)
	{
		Eigen::Vector3d p, p0;
		double l_ij = 0;
		OpenMesh::HalfedgeHandle he_t;
		OpenMesh::VertexHandle fromV;

		for (auto v : domain.vertices())
		{
			ASAP asap(10);
			double s = allowscale ? 1 : mesh.property(affine_M, v).S;

			for (auto vih : mesh.vih_range(v))
			{
				fromV = mesh.from_vertex_handle(vih);
				p = mesh.point(v) - mesh.point(fromV);
				p0 = mesh.property(point0, v) - mesh.property(point0, fromV);

				l_ij = 0;

				if (domain.inDomain(vih))
					l_ij += 1;
				he_t = mesh.opposite_halfedge_handle(vih);
				if (domain.inDomain(he_t))
					l_ij += 1;
				
				asap.addCouples(p0, p, l_ij);
			}

			mesh.property(affine_M, v) = asap.solve(s);
		}
	}

	//-------------------------------------------------------

	/**** ----- inline function for meta operation ----- ****/

	inline Eigen::Vector3d calcLP(const Mesh::VertexHandle &v)
	{
		double l_ii = 0.0;
		double l_ij = 0.0;
		OpenMesh::HalfedgeHandle he_t;
		Eigen::Vector3d rt(0, 0, 0);
		OpenMesh::VertexHandle fromV;

		for (auto vih : mesh.vih_range(v))
		{
			fromV = mesh.from_vertex_handle(vih);
			l_ij = 0;

			if (domain.inDomain(vih))
				l_ij -= mesh.property(cot_theta, vih);
			he_t = mesh.opposite_halfedge_handle(vih);
			if (domain.inDomain(he_t))
				l_ij -= mesh.property(cot_theta, he_t);

			l_ii -= l_ij;
			rt += l_ij * mesh.point(fromV);
		}
		return rt + l_ii * mesh.point(v);
	}

	inline Eigen::Vector3d calcEp(const Mesh::HalfedgeHandle &h)
	{
		return mesh.point(mesh.to_vertex_handle(h)) - mesh.point(mesh.from_vertex_handle(h));
	}

	inline Eigen::Vector3d massSpring(Eigen::Vector3d v, double l)
	{
		return (v.norm() + l) / 2 * v.normalized();
	}

	inline Eigen::Vector2d massSpring(Eigen::Vector2d v, double l)
	{
		return (v.norm() + l) / 2 * v.normalized();
	}

	inline double massSpring(Eigen::Vector3d v, Eigen::Vector3d &x, double l)
	{
		x = (v.norm() + l) / 2 * v.normalized();
		return fabs(x.norm() - l);
	}

	inline double massSpring(Eigen::Vector2d v, Eigen::Vector2d &x, double l)
	{
		x = (v.norm() + l) / 2 * v.normalized();
		return fabs(x.norm() - l);
	}

	inline double calcVertexMixedArea(const Mesh::VertexHandle &v)
	{
		OpenMesh::HalfedgeHandle he_a, he_b;
		double mixArea = 0.0;
		double alpha, beta, gama;

		for (auto vih : mesh.vih_range(v))
		{
			if (!domain.inDomain(vih)) continue;
			he_a = vih;
			he_b = mesh.next_halfedge_handle(he_a);
			alpha = mesh.property(op_theta, he_a);
			beta = mesh.property(op_theta, he_b);
			gama = M_PI - alpha - beta;
			if (alpha < M_PI / 2 && beta < M_PI / 2 && gama < M_PI / 2)
			{
				mixArea += pow(mesh.property(para_l, mesh.edge_handle(he_a)), 2.0)*tan(M_PI / 2 - alpha) / 8;
				mixArea += pow(mesh.property(para_l, mesh.edge_handle(he_b)), 2.0)*tan(M_PI / 2 - beta) / 8;
			}
			else
			{
				if (gama > M_PI / 2)
				{
					mixArea += mesh.property(para_l, mesh.edge_handle(he_a))*
						mesh.property(para_l, mesh.edge_handle(he_b))*sin(gama) / 4;
				}
				else
				{
					mixArea += mesh.property(para_l, mesh.edge_handle(he_a))*
						mesh.property(para_l, mesh.edge_handle(he_b))*sin(gama) / 8;
				}
			}
		}
		return mixArea;
	}

	inline double calcELength(const Mesh::EdgeHandle &e_input)
	{
		OpenMesh::HalfedgeHandle he_t = mesh.halfedge_handle(e_input, 0);
		double ui = mesh.property(para_u, mesh.from_vertex_handle(he_t));
		double uj = mesh.property(para_u, mesh.to_vertex_handle(he_t));
		return exp(ui + uj)*mesh.calc_edge_vector(e_input).norm();
	}

	inline bool calcTheta(double la, double lb, double lc, double &A, double thres)
	{
		double cos_value = (lb / lc + lc / lb - la / lb * la / lc) / 2;
		bool inQuality = cos_value <= thres && cos_value >= -thres;
		if (cos_value > thres)
		{
			cos_value = thres;
		}
		else if (cos_value < -thres)
		{
			cos_value = -thres;
		}
		A = acos(cos_value);
		return inQuality;
	}

	inline bool calcHeOppositeTheta(const Mesh::HalfedgeHandle &he_input, double &theta, double thres)
	{
		OpenMesh::HalfedgeHandle he_a = he_input;
		OpenMesh::HalfedgeHandle he_b = mesh.next_halfedge_handle(he_a);
		OpenMesh::HalfedgeHandle he_c = mesh.next_halfedge_handle(he_b);

		double la = mesh.property(para_l, mesh.edge_handle(he_a));
		double lb = mesh.property(para_l, mesh.edge_handle(he_b));
		double lc = mesh.property(para_l, mesh.edge_handle(he_c));

		return calcTheta(la, lb, lc, theta, thres);
	}

	inline bool attemptFlip(const OpenMesh::EdgeHandle &e, ESet &add)
	{
		auto hcd = mesh.halfedge_handle(e, 0);
		auto hda = mesh.next_halfedge_handle(hcd);
		auto hac = mesh.next_halfedge_handle(hda);

		auto hdc = mesh.halfedge_handle(e, 1);
		auto hcb = mesh.next_halfedge_handle(hdc);
		auto hbd = mesh.next_halfedge_handle(hcb);

		auto va = mesh.opposite_vh(hcd);
		auto vb = mesh.opposite_vh(hdc);
		double ua = mesh.property(para_u, va);
		double ub = mesh.property(para_u, vb);

		double lab = exp(ua + ub)*(mesh.point(va) - mesh.point(vb)).norm();
		double lbd = mesh.property(para_l, mesh.edge_handle(hbd));
		double lda = mesh.property(para_l, mesh.edge_handle(hda));
		double lac = mesh.property(para_l, mesh.edge_handle(hac));
		double lcb = mesh.property(para_l, mesh.edge_handle(hcb));

		bool flip_ok = true;
		double abd, bda, dab, cba, bac, acb;
		
		flip_ok = flip_ok && calcTheta(lab, lbd, lda, bda, THETATHRES);
		flip_ok = flip_ok && calcTheta(lbd, lda, lab, dab, THETATHRES);
		flip_ok = flip_ok && calcTheta(lda, lab, lbd, abd, THETATHRES);

		flip_ok = flip_ok && calcTheta(lcb, lab, lac, bac, THETATHRES);
		flip_ok = flip_ok && calcTheta(lab, lac, lcb, acb, THETATHRES);
		flip_ok = flip_ok && calcTheta(lac, lcb, lab, cba, THETATHRES);

		if (!flip_ok) return false;

		mesh.flip(e);
		mesh.property(para_l, e) = lab;

		mesh.property(op_theta, hcd) = acb;
		mesh.property(op_theta, hdc) = bda;
		mesh.property(op_theta, hbd) = dab;
		mesh.property(op_theta, hda) = abd;
		mesh.property(op_theta, hac) = cba;
		mesh.property(op_theta, hcb) = bac;

		OpenMesh::EdgeHandle tmp;
		tmp = mesh.edge_handle(hbd);
		if (!mesh.property(flip_indicate, tmp) && !domain.onBoundary(tmp))
			add.push_back(tmp);
		tmp = mesh.edge_handle(hda);
		if (!mesh.property(flip_indicate, tmp) && !domain.onBoundary(tmp))
			add.push_back(tmp);
		tmp = mesh.edge_handle(hac);
		if (!mesh.property(flip_indicate, tmp) && !domain.onBoundary(tmp))
			add.push_back(tmp);
		tmp = mesh.edge_handle(hcb);
		if (!mesh.property(flip_indicate, tmp) && !domain.onBoundary(tmp))
			add.push_back(tmp);

		return true;
	}

	//--------------------------------------------------------------------
	void calc_uv_coordinates()
	{
		mesh.request_vertex_texcoords2D();

		OpenMesh::FPropHandleT<bool> face_vistor;
		std::list<Mesh::HalfedgeHandle> stack_he;
		mesh.add_property(face_vistor);
		for (auto f : mesh.faces())
		{
			mesh.property(face_vistor, f) = true;
		}
		for (auto f : domain.faces())
		{
			mesh.property(face_vistor, f) = false;
		}

		auto he_t = mesh.halfedge_handle(domain.faces()[0]);
		double l_t = mesh.property(para_l, mesh.edge_handle(he_t));
		mesh.set_texcoord2D(mesh.from_vertex_handle(he_t), Eigen::Vector2d(0, 0));
		mesh.set_texcoord2D(mesh.to_vertex_handle(he_t), Eigen::Vector2d(l_t, 0));

		stack_he.push_back(he_t);
		mesh.property(face_vistor, mesh.face_handle(he_t)) = true;

		double fr, tr, dt, dn;
		Eigen::Vector2d fpos, tpos, dpos;
		OpenMesh::VertexHandle fv_t, tv_t;
		OpenMesh::HalfedgeHandle he_n, he_nn;
		while (!stack_he.empty())
		{
			he_t = stack_he.front();
			stack_he.pop_front();

			fv_t = mesh.from_vertex_handle(he_t);
			tv_t = mesh.to_vertex_handle(he_t);
			fpos = mesh.texcoord2D(fv_t);
			tpos = mesh.texcoord2D(tv_t);
			l_t = mesh.property(para_l, mesh.edge_handle(he_t));

			he_n = mesh.next_halfedge_handle(he_t);
			tr = mesh.property(para_l, mesh.edge_handle(he_n));
			he_nn = mesh.next_halfedge_handle(he_n);
			fr = mesh.property(para_l, mesh.edge_handle(he_nn));

			dt = (l_t + fr / l_t * fr - tr / l_t * tr) / 2;
			dn = 1 - dt / fr * dt / fr;
			assert(dn >= 0);
			dn = fr * sqrt(dn);
			dpos = (tpos - fpos).normalized();
			dpos = fpos + dt * dpos + dn * Eigen::Vector2d(-dpos[1], dpos[0]);
			mesh.set_texcoord2D(mesh.to_vertex_handle(he_n), dpos);

			he_t = mesh.opposite_halfedge_handle(he_t);
			if (mesh.face_handle(he_t).is_valid() &&
				!mesh.property(face_vistor, mesh.face_handle(he_t)))
			{
				stack_he.push_back(he_t);
				mesh.property(face_vistor, mesh.face_handle(he_t)) = true;
			}
			he_n = mesh.opposite_halfedge_handle(he_n);
			if (mesh.face_handle(he_n).is_valid() &&
				!mesh.property(face_vistor, mesh.face_handle(he_n)))
			{
				stack_he.push_back(he_n);
				mesh.property(face_vistor, mesh.face_handle(he_n)) = true;
			}
			he_nn = mesh.opposite_halfedge_handle(he_nn);
			if (mesh.face_handle(he_nn).is_valid() &&
				!mesh.property(face_vistor, mesh.face_handle(he_nn)))
			{
				stack_he.push_back(he_nn);
				mesh.property(face_vistor, mesh.face_handle(he_nn)) = true;
			}
		}

		mesh.remove_property(face_vistor);
		//unifyTexcoords2D();
	}

	void unifyTexcoords2D()
	{
		if (!mesh.has_vertex_texcoords2D())
		{
			std::cout << "No texcoords2D\n";
			exit(-1);
		}

		Eigen::Vector2d uvmin(DBL_MAX, DBL_MAX), uvmax(DBL_MIN, DBL_MIN), uvtmp;
		for (auto v : domain.vertices())
		{
			uvtmp = mesh.texcoord2D(v);
			uvmin = uvmin.cwiseMin(uvtmp);
			uvmax = uvmax.cwiseMax(uvtmp);
		}

		double scale = fmax(uvmax[0] - uvmin[0], uvmax[1] - uvmin[1]);
		for (auto v : domain.vertices())
		{
			uvtmp = mesh.texcoord2D(v);
			uvtmp = (uvtmp - uvmin) / scale;
			mesh.set_texcoord2D(v, uvtmp);
		}
	}
};

#endif // !CURVDEFORM_H