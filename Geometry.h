/*---
  Flow*: A Verification Tool for Cyber-Physical Systems.
  Authors: Xin Chen, Sriram Sankaranarayanan, and Erika Abraham.
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.
  
  The code is released as is under the GNU General Public License (GPL).
---*/

#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include "Constraints.h"

namespace flowstar
{

class Polyhedron
{
public:
	std::vector<LinearConstraint> constraints;

public:
	Polyhedron();
	Polyhedron(const std::vector<LinearConstraint> & cs);
	Polyhedron(const Polyhedron & P);
//	Polyhedron(Matrix<double> & A, Matrix<double> & b);
	Polyhedron(const std::vector<std::vector<Real> > & A, const std::vector<Real> & B);
	~Polyhedron();

	Polyhedron(const std::vector<std::vector<Real> > & template_matrix, const TaylorModelVec<Real> & tmv, const std::vector<Interval> & domain);

	double rho(const std::vector<Real> & l) const;
//	double rho(Matrix<Real> & l) const;
	void tightenConstraints();
	bool empty() const;
	void get(std::vector<std::vector<Real> > & A, std::vector<Real> & B) const;

	void output(std::ostream & os, const Variables & stateVars) const;

	Polyhedron & operator = (const Polyhedron & P);
};





class Parallelotope									// H-Representation
{
public:
	Matrix<Real> paraTemplate;					// only half of the facet normals are kept
	Matrix<Real> b;

	Parallelotope(const Matrix<Real> & A, const Matrix<Real> & B);
	Parallelotope(const Parallelotope & P);
	~Parallelotope();

	void center(std::vector<Real> & c);				// center of the parallelotope.
	void output(std::ostream & os);

	// converse a parallelotope to a Taylor model, the domain is normalized
	void toTaylorModel(TaylorModelVec<Real> & result);

	Parallelotope & operator = (const Parallelotope & P);
};





class Zonotope
{
public:
	Matrix<Real> center;
	std::list<Matrix<Real> > generators;

public:
	Zonotope();
	Zonotope(const Matrix<Real> & c, const std::list<Matrix<Real> > & G);
	Zonotope(const Zonotope & zonotope);
	Zonotope(const std::vector<Interval> & box);
	Zonotope(Matrix<Interval> & box);
	Zonotope(const unsigned int d);
	~Zonotope();

	bool isEmpty() const;
//	bool isIntersected(Matrix<Interval> & box) const;

	unsigned int numOfGen() const;

	void simplify();

	void intEval(Matrix<Interval> & range);

//	void output(FILE *fp) const;


	Zonotope & operator += (const Zonotope & Z);

	friend Zonotope operator + (const Zonotope & Z1, const Zonotope & Z2);
	friend Zonotope operator * (const Matrix<Real> & A, const Zonotope & Z);

/*
	void linearTrans(Zonotope & result, const iMatrix & map) const;
	void linearTrans_assign(const iMatrix & map);

	void MinSum(Zonotope & result, const Zonotope & zonotope) const;
	void MinSum_assign(const Zonotope & zonotope);
*/


	// translate the zonotope to a Taylor model with a zero remainder
	// the first variable is NOT reserved by t
	void toPolynomial(std::vector<Polynomial<Real> > & result);

	Zonotope & operator = (const Zonotope & zonotope);
	Zonotope & operator = (const std::vector<Interval> & box);
	Zonotope & operator = (Matrix<Interval> & box);


/*
	bool belongsto(const std::vector<double> & x);
	int contract(const LinearConstraint & constraint);

	void to2DBox(Matrix<Interval> & box, const int x, const int y);
	void intervalRange(Interval & range, const int x);

	void plot(FILE *fp, const int x, const int y);	// only for 2D zonotopes
*/
};


}

#endif /* GEOMETRY_H_ */
