#include <assert.h>
#include <iostream>

#define GLEW_STATIC
#include <GL/glew.h>

#include <glm/gtc/type_ptr.hpp>
#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>

#include "Link.h"
#include "Shape.h"
#include "MatrixStack.h"
#include "Program.h"
#include "OptimizerGDLS.h"
#include "OptimizerNM.h"
#include "ObjectiveReal.h"

using namespace std;
using namespace Eigen;

std::vector<std::shared_ptr<Link>> Link::links;
Eigen::VectorXd Link::thetas;

Link::Link() : 
	parent(NULL),
	child(NULL),
	depth(0),
	position(0.0, 0.0),
	angle(0.0),
	meshMat(Matrix4d::Identity()),
	end_effector(1.0, 0.0, 1.0),
	wtar(1000),
	wreg(1),
	NLINKS(1)
{
}

Link::Link(int n) :
	parent(NULL),
	child(NULL),
	depth(0),
	position(0.0, 0.0),
	angle(0.0),
	meshMat(Matrix4d::Identity()),
	end_effector(1.0, 0.0, 1.0),
	wtar(1000),
	wreg(1),
	NLINKS(n)
{
	
}

Link::~Link()
{
	
}

void Link::addChild(shared_ptr<Link> child)
{
	child->parent = shared_from_this();
	child->depth = depth + 1;
	this->child = child;
}

void Link::draw(const shared_ptr<Program> prog, shared_ptr<MatrixStack> MV, const shared_ptr<Shape> shape) const
{
	assert(prog);
	assert(MV);
	assert(shape);
	
	MV->pushMatrix();
	Matrix4d jointMat;
	jointMat << cos(angle), -1 * sin(angle), 0, position.x(),
				sin(angle), cos(angle), 0, position.y(),
				0, 0, 1, 0,
				0, 0, 0, 1;
	MV->multMatrix(jointMat);
		MV->pushMatrix();
		MV->multMatrix(meshMat);
		glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
		shape->draw();
		MV->popMatrix();
	if (child != nullptr) {
		child->draw(prog, MV, shape);
	}
	
	MV->popMatrix();
}


void Link::set_links(std::vector<std::shared_ptr<Link>> l)
{
	links = l;
}

void Link::set_theta(Eigen::VectorXd t) {
	thetas = t;
}

void Link::set_thetas()
{
	// accounts for wrapping
	for (int i = 0; i < links.size(); i++) {
		while (thetas(i) > M_PI) {
			thetas(i) -= (2 * M_PI);
		}
		while (thetas(i) < -M_PI) {
			thetas(i) += (2 * M_PI);
		}
		links.at(i)->setAngle(thetas(i));
	}
}

void Link::calculate_thetas(Eigen::Vector3d& targetPos)
{
	int NLINKS = links.at(0)->getn();
	OptimizerGDLS op(NLINKS); // Need to move NLINKS here
	OptimizerNM opNM(NLINKS);
	Vector3d target;
	target << 0.0, 1.0, 1.0;
	auto objective = make_shared<ObjectiveReal>(NLINKS, targetPos);
	VectorXd xInit(NLINKS);
	VectorXd g(NLINKS);
	MatrixXd H(NLINKS, NLINKS);
	xInit = thetas;
	double f = objective->evalObjective(xInit, g, H);

	VectorXd theta_gdls = op.optimize(objective, xInit);
	VectorXd theta_nm = opNM.optimize(objective, theta_gdls);

	double f_nm = objective->evalObjective(theta_nm);
	if (f_nm < f) {
		for (int i = 0; i < thetas.rows(); i++) {
			thetas(i) = theta_nm(i);
		}
	}
	else {
		for (int i = 0; i < thetas.rows(); i++) {
			thetas(i) = theta_gdls(i);
		}
	}

	for (int i = 0; i < thetas.rows(); i++) {
		thetas(i) = theta_gdls(i);
	}
	set_thetas();
}
