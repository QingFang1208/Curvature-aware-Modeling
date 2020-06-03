#ifndef AFFINEMAP_H
#define AFFINEMAP_H

#include <Eigen/Eigen>
#include <Eigen/Dense>

struct AffineMap
{
	Eigen::Matrix3d R;
	Eigen::Vector3d T;
	double S;
	AffineMap() { R.setIdentity(); T.setZero(); S = 1; }
	AffineMap(const Eigen::Matrix3d &IR, const double &IS) : R(IR), S(IS) { T.setZero(); }
	Eigen::Vector3d map(const Eigen::Vector3d &p) { return R * p*S + T; }
};

class ASAP
{
private:
	Eigen::Matrix3d cov;
	std::vector<Eigen::Vector3d> src;
	std::vector<Eigen::Vector3d> tar;
	std::vector<double> w;

public:
	ASAP(int cap = 10) { src.reserve(cap); tar.reserve(cap); w.reserve(cap); }
	void addCouples(const Eigen::Vector3d &source, const Eigen::Vector3d &target, double weight = 1)
	{
		src.push_back(source);
		tar.push_back(target);
		w.push_back(weight);
	}

	void clear()
	{
		src.clear();
		tar.clear();
		w.clear();
	}

	AffineMap solve(double s = -1)
	{
		if (src.empty()) return AffineMap();

		cov.setZero();
		for (int i = 0; i < src.size(); i++)
		{
			cov += w[i] * src[i] * tar[i].transpose();
		}
		
		Eigen::JacobiSVD<Eigen::MatrixXd> svd(cov, Eigen::ComputeThinU | Eigen::ComputeThinV);
		Eigen::MatrixXd U = svd.matrixU();
		Eigen::MatrixXd V = svd.matrixV();
		Eigen::MatrixXd R = V * U.transpose();
		
		if (R.determinant() < 0)
		{
			U.col(2) = -U.col(2);
			R = V * U.transpose();
			assert(R.determinant() > 0);
		}

		if (s > 0)
			return AffineMap(R, s);

		double nume = 0, deno = 0;
		for (int i = 0; i < src.size(); i++)
		{
			nume += w[i] * tar[i].dot(R*src[i]);
			deno += w[i] * src[i].dot(src[i]);
		}

		return AffineMap(R, nume / deno);
	}
};

#endif // !AFFINEMAP_H