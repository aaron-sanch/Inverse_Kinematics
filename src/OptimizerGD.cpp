#include "OptimizerGD.h"
#include "Objective.h"
#include <iostream>

using namespace std;
using namespace Eigen;

OptimizerGD::OptimizerGD() :
	alpha(1.0),
	tol(1e-6),
	iterMax(100),
	iter(0)
{
	
}

OptimizerGD::~OptimizerGD()
{
	
}

VectorXd OptimizerGD::optimize(const shared_ptr<Objective> objective, const VectorXd &xInit)
{
	int n = xInit.rows();
	VectorXd x = xInit;
	VectorXd g(n);
	iter = 0;
	for (int i = 0; i < iterMax; i++) {
		double f = objective->evalObjective(x, g);
		VectorXd dx = -1 * alpha * g;
		x += dx;
		if (dx.norm() < tol) {
			iter = i + 1;
			break;
		}
	}
	return x;
}
