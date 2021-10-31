#include "OptimizerGDLS.h"
#include "Objective.h"
#include <iostream>

using namespace std;
using namespace Eigen;

OptimizerGDLS::OptimizerGDLS() :
	alphaInit(1e-4),
	gamma(0.5),
	tol(1e-3),
	iterMax(100),
	iter(0)
{
	
}

OptimizerGDLS::OptimizerGDLS(int n) :
	alphaInit(1e-4),
	gamma(0.5),
	tol(1e-3),
	iterMax(10 * n),
	iter(0)
{

}

OptimizerGDLS::~OptimizerGDLS()
{
	
}

VectorXd OptimizerGDLS::optimize(const shared_ptr<Objective> objective, const VectorXd &xInit)
{
	int n = xInit.rows();
	if (n <= 2) {
		alphaInit = 1;
	}
	VectorXd x = xInit;
	VectorXd g(n);
	iter = 0;
	for (int i = 0; i < iterMax; i++) {
		double f = objective->evalObjective(x, g);
		double alpha = alphaInit;
		VectorXd dx(n);
		for (int j = 0; j < 20; j++) {
			dx = -1 * alpha * g;
			double fnew = objective->evalObjective(x + dx);
			if (fnew < f) {
				break;
			}
			alpha *= gamma;
		}
		x += dx;
		if (dx.norm() < tol) {
			iter = i + 1;
			break;
		}
	}
	return x;
}
