/*---
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.
  The code is released as is under the GNU General Public License (GPL).
---*/

#ifndef CONSTRAINTS_H_
#define CONSTRAINTS_H_

#include "TaylorModel.h"
#include "expression.h"

namespace flowstar
{

// Define the classes of linear and nonlinear constraints

class LinearConstraint		// A x <= B
{
public:
	std::vector<Real> A;
	Real B;

public:
	LinearConstraint();
	LinearConstraint(const std::vector<Real> & a, const Real & b);
	LinearConstraint(const LinearConstraint & lc);
	~LinearConstraint();

	void output(std::ostream & os, const Variables & stateVars) const;

	LinearConstraint & operator = (const LinearConstraint & lc);
};

class PolynomialConstraint	// p(x) <= b
{
public:
	Polynomial<Real> p;
	HornerForm<Real> hf;		// a HornerForm of p
	Real B;

public:
	PolynomialConstraint();
	PolynomialConstraint(const Polynomial<Real> & polynomial, const Real & b);
	PolynomialConstraint(const PolynomialConstraint & pc);
//	PolynomialConstraint(const std::string & strPolynomial, const Variables & vars);
	~PolynomialConstraint();

	void output(std::ostream & os, const Variables & stateVars) const;

	PolynomialConstraint & operator = (const PolynomialConstraint & pc);
};

class Constraint
{
public:
	Expression_AST<Real> expression;
	Real bound;

public:
	Constraint();
	Constraint(const Expression_AST<Real> & exp, const Real & b);
	Constraint(const std::string & strExpression);
	Constraint(const Constraint & constraint);
	~Constraint();

	void output(std::ostream & os, const Variables & stateVars) const;

	Constraint & operator = (const Constraint & constraint);
};


}

#endif /* CONSTRAINTS_H_ */
