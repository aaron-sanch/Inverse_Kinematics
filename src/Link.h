#pragma once
#ifndef Link_H
#define Link_H

#include <memory>

#include <Eigen/Dense>

class Shape;
class MatrixStack;
class Program;

class Link : public std::enable_shared_from_this<Link>
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	
	Link();
	Link(int n);
	virtual ~Link();

	void addChild(std::shared_ptr<Link> child);
	std::shared_ptr<Link> getChild() const { return child; }
	std::shared_ptr<Link> getParent() const { return parent; }
	int getDepth() const { return depth; }

	void setPosition(double x, double y) { position << x, y; }
	const Eigen::Vector2d &getPosition() const { return position; }

	void setAngle(double a) { angle = a; }
	double getAngle() const { return angle; }

	void setMeshMatrix(const Eigen::Matrix4d &M) { meshMat = M; }

	void draw(const std::shared_ptr<Program> prog, std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Shape> shape) const;

	int getn() { return NLINKS; }

	static void set_links(std::vector<std::shared_ptr<Link>> l);
	static void set_theta(Eigen::VectorXd t);
	static void set_thetas();
	static void calculate_thetas(Eigen::Vector3d& targetPos);
	
private:
	std::shared_ptr<Link> parent;
	std::shared_ptr<Link> child;
	int depth;
	Eigen::Vector2d position;
	double angle;
	Eigen::Matrix4d meshMat;
	Eigen::Vector3d end_effector;
	double wtar;
	double wreg;
	int NLINKS;
	
	static Eigen::VectorXd thetas;
	static std::vector<std::shared_ptr<Link>> links;
};

#endif
