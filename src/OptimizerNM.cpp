#include "OptimizerNM.h"
#include "Objective.h"
#include <iostream>

using namespace std;
using namespace Eigen;

OptimizerNM::OptimizerNM() :
	tol(1e-6),
	iterMax(100),
	iter(0)
{
	
}

OptimizerNM::OptimizerNM(int n) :
	tol(1e-3),
	iterMax(10 * n),
	iter(0)
{
}

OptimizerNM::~OptimizerNM()
{
	
}

VectorXd OptimizerNM::optimize(const shared_ptr<Objective> objective, const VectorXd &xInit)
{
	int n = xInit.rows();
	VectorXd x = xInit;
	VectorXd g(n);
	MatrixXd H(n, n);
	iter = 0;
	for (int i = 0; i < iterMax; i++) {
		double f = objective->evalObjective(x, g, H);
		VectorXd dx = -1 * H.ldlt().solve(g);
		x += dx;
		if (dx.norm() < tol) {
			iter = i + 1;
			break;
		}
	}
	return x;
}
