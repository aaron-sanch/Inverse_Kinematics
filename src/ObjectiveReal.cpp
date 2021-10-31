#include "ObjectiveReal.h"
#include <cmath>
#include <iostream>

using namespace Eigen;
using namespace std;


ObjectiveReal::ObjectiveReal(): 
	n(0),
	pTar(0.0, 0.0, 0.0)
{
}

ObjectiveReal::ObjectiveReal(int nlink, Eigen::Vector3d pTarget):
	n(nlink),
	pTar(pTarget)
{
}


Eigen::VectorXd ObjectiveReal::calculate_p(const int curr_iter, const int size, const Eigen::VectorXd& theta) const
{
	VectorXd p(3);
	p << 1.0, 0.0, 1.0;
	VectorXd end_effector(3);
	end_effector << 1.0, 0.0, 1.0;
	Matrix3d R;
	R << cos(theta(curr_iter)), -1 * sin(theta(curr_iter)), 0,
		sin(theta(curr_iter)), cos(theta(curr_iter)), 0,
		0, 0, 1;
	Matrix3d T;
	T.setIdentity();
	if (curr_iter != 0) {
		T << 1.0, 0, 1.0,
			0.0, 1.0, 0.0,
			0.0, 0.0, 1.0;
	}
	if (curr_iter == (size - 1)) {
		p = (T * R * end_effector);
		return p;
	}
	VectorXd temp = calculate_p(curr_iter + 1, size, theta);
	VectorXd curr(3);
	curr = T * R * temp;
	p(0) = curr(0);
	p(1) = curr(1);
	p(2) = curr(2);
	return p;
}

Eigen::MatrixXd ObjectiveReal::calculate_p_prime(const int size, const Eigen::VectorXd& theta) const
{
	Vector3d p;
	Vector3d end_effector;
	end_effector << 1.0, 0.0, 1.0;
	MatrixXd p_prime(3, size);
	p_prime.setIdentity();
	for (int i = 0; i < size; i++) {
		Vector3d p_curr;
		for (int j = size - 1; j >= 0; j--) {
			MatrixXd R(3, 3);
			R << cos(theta(j)), -1 * sin(theta(j)), 0,
				sin(theta(j)), cos(theta(j)), 0,
				0, 0, 1;
			MatrixXd T(3, 3);
			T.setIdentity();
			if (j != 0) {
				T << 1.0, 0, 1.0,
					0.0, 1.0, 0.0,
					0.0, 0.0, 1.0;
			}
			if (i == j) {
				R << -sin(theta(j)), -cos(theta(j)), 0,
					cos(theta(j)), -sin(theta(j)), 0,
					0, 0, 0;
			}
			if (j == (size - 1)) {
				p_curr = T * R * end_effector;
			}
			else {
				p_curr = T * R * p_curr;
			}
			
		}
		p_prime.col(i) = p_curr;
	}
	return p_prime;
}

Eigen::MatrixXd ObjectiveReal::calculate_p_prime_2(const int size, const Eigen::VectorXd& theta) const
{
	Vector3d p;
	Vector3d end_effector;
	end_effector << 1.0, 0.0, 1.0;
	MatrixXd p_prime_2(3 * size, size);
	p_prime_2.setIdentity();
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			Vector3d p_curr;
			for (int l = size -1; l >= 0; l--) {
				MatrixXd R(3, 3);
				R << cos(theta(l)), -1 * sin(theta(l)), 0,
					sin(theta(l)), cos(theta(l)), 0,
					0, 0, 1;
				MatrixXd T(3, 3);
				T.setIdentity();
				if (l != 0) {
					T << 1.0, 0, 1.0,
						0.0, 1.0, 0.0,
						0.0, 0.0, 1.0;
				}
				if (i == j && j == l) {
					R << -cos(theta(l)), sin(theta(l)), 0,
						-sin(theta(l)), -cos(theta(l)), 0,
						0, 0, 0;
				}
				else if (l == j || l == i) {
					R << -sin(theta(l)), -cos(theta(l)), 0,
						cos(theta(l)), -sin(theta(l)), 0,
						0, 0, 0;
				}
				if (l == (size - 1)) {
					p_curr = T * R * end_effector;
				}
				else {
					p_curr = T * R * p_curr;
				}
			}
			p_prime_2.col(j)(i*3) = p_curr(0);
			p_prime_2.col(j)(i*3 + 1) = p_curr(1);
			p_prime_2.col(j)(i*3 + 2) = p_curr(2);
		}
	}
	return p_prime_2;
}


ObjectiveReal::~ObjectiveReal()
{
}


double ObjectiveReal::evalObjective(const Eigen::VectorXd& theta) const
{
	VectorXd deltaP(3);
	VectorXd p = calculate_p(0, theta.rows(),theta);
	deltaP = p - pTar;
	double pVal = deltaP.transpose() * deltaP;
	double thetaVal = theta.transpose() * theta;
	return ((1 / 2.0 * wtar * pVal) + (1 / 2.0 * wreg * thetaVal));
}

double ObjectiveReal::evalObjective(const Eigen::VectorXd& theta, Eigen::VectorXd& g) const
{
	VectorXd deltaP(3);
	VectorXd p = calculate_p(0, theta.rows(), theta);
	MatrixXd p_prime = calculate_p_prime(theta.rows(), theta);
	deltaP = p - pTar;
	g = (wtar * deltaP.transpose() * p_prime).transpose() + (wreg * theta);
	double pVal = deltaP.transpose() * deltaP;
	double thetaVal = theta.transpose() * theta;
	return ((1 / 2.0 * wtar * pVal) + (1 / 2.0 * wreg * thetaVal));
}

double ObjectiveReal::evalObjective(const Eigen::VectorXd& theta, Eigen::VectorXd& g, Eigen::MatrixXd& H) const
{
	VectorXd deltaP(3);
	VectorXd p = calculate_p(0, theta.rows(), theta);
	MatrixXd p_prime = calculate_p_prime(theta.rows(), theta);
	MatrixXd p_prime_2 = calculate_p_prime_2(theta.rows(), theta);
	deltaP = p - pTar;
	g = (wtar * deltaP.transpose() * p_prime).transpose() + (wreg * theta);
	double pVal = deltaP.transpose() * deltaP;
	double thetaVal = theta.transpose() * theta;
	MatrixXd p_h_1 = p_prime.transpose() * p_prime_2;
	MatrixXd p_h_2 = deltaP.transpose() * p_prime_2;
	MatrixXd identity(n, n);
	identity.setIdentity();
	MatrixXd p_for_H(n, n);
	//weird circle dot multiplication
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			VectorXd to_mult(3);
			to_mult(0) = p_prime_2(i * 3 , j);
			to_mult(1) = p_prime_2(i * 3 + 1, j);
			to_mult(2) = p_prime_2(i * 3 + 2, j);
			double sum_p_2 = (deltaP.transpose() * to_mult);
			p_for_H(i, j) = sum_p_2;
		}
	}
	H = (wtar * ((p_prime.transpose() * p_prime) + p_for_H)) + (wreg * identity);
 	return ((1 / 2.0 * wtar * pVal) + (1 / 2.0 * wreg * thetaVal));
}
