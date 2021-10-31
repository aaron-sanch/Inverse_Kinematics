#include "ObjectiveLab.h"
#include <cmath>

using namespace Eigen;

ObjectiveLab::ObjectiveLab()
{
	
}

ObjectiveLab::~ObjectiveLab()
{
	
}

double ObjectiveLab::evalObjective(const VectorXd &x) const
{
	return (0.5 * sin(x.x())) + (0.5 * sin(x.y())) + (1 / 20.0 * x.x() * x.x()) + (1 / 20.0 * x.x() * x.y()) + (1 / 20.0 * x.y() * x.y());
}

double ObjectiveLab::evalObjective(const VectorXd &x, VectorXd &g) const
{
	g << ((0.5 * cos(x.x())) + (1 / 10.0 * x.x()) + (1 / 20.0 * x.y())),
		((0.5 * cos(x.y())) + (1 / 10.0 * x.y()) + (1 / 20.0 * x.x()));
	return (0.5 * sin(x.x())) + (0.5 * sin(x.y())) + (1 / 20.0 * x.x() * x.x()) + (1 / 20.0 * x.x() * x.y()) + (1 / 20.0 * x.y() * x.y());
}

double ObjectiveLab::evalObjective(const VectorXd &x, VectorXd &g, MatrixXd &H) const
{
	g << ((0.5 * cos(x.x())) + (1 / 10.0 * x.x()) + (1 / 20.0 * x.y())),
		((0.5 * cos(x.y())) + (1 / 10.0 * x.y()) + (1 / 20.0 * x.x()));
	H << ((-0.5 * sin(x.x())) + 1 / 10.0), (1 / 20.0),
		(1 / 20.0), ((-0.5 * sin(x.y())) + 1 / 10.0);
	return (0.5 * sin(x.x())) + (0.5 * sin(x.y())) + (1 / 20.0 * x.x() * x.x()) + (1 / 20.0 * x.x() * x.y()) + (1 / 20.0 * x.y() * x.y());
}
