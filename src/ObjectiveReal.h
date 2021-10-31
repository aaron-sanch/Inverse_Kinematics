#pragma once
#ifndef OBJECTIVE_REAL_H
#define OBJECTIVE_REAL_H

#include "Objective.h"

class ObjectiveReal : public Objective
{
private:
	Eigen::Vector3d pTar;
	int n;
	int wtar = 1000;
	int wreg = 1;

public:
	ObjectiveReal();
	ObjectiveReal(int nlinks, Eigen::Vector3d pTarget);

	Eigen::VectorXd calculate_p(const int curr_iter, const int size, const Eigen::VectorXd& theta) const;
	Eigen::MatrixXd calculate_p_prime(const int size, const Eigen::VectorXd& theta) const;
	Eigen::MatrixXd calculate_p_prime_2(const int size, const Eigen::VectorXd& theta) const;
	virtual ~ObjectiveReal();
	virtual double evalObjective(const Eigen::VectorXd& theta) const;
	virtual double evalObjective(const Eigen::VectorXd& theta, Eigen::VectorXd& g) const;
	virtual double evalObjective(const Eigen::VectorXd& theta, Eigen::VectorXd& g, Eigen::MatrixXd& H) const;
};

#endif
