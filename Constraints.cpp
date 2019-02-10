/*---
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.
  The code is released as is under the GNU General Public License (GPL).
---*/

#include "Constraints.h"

using namespace flowstar;

LinearConstraint::LinearConstraint()
{
}

LinearConstraint::LinearConstraint(const std::vector<Real> & a, const Real & b)
{
	A = a;
	B = b;
}

LinearConstraint::LinearConstraint(const LinearConstraint & lc)
{
	A = lc.A;
	B = lc.B;
}

LinearConstraint::~LinearConstraint()
{
}

void LinearConstraint::output(std::ostream & os, const Variables & stateVars) const
{
	int d = A.size();
	for(int i=0; i<d-1; ++i)
	{
		os << "(" << A[i] << "*" << stateVars.varNames[i] << ") + ";
	}

	os << "(" << A[d-1] << "*" << stateVars.varNames[d-1] << " <= " << B << std::endl;
}

LinearConstraint & LinearConstraint::operator = (const LinearConstraint & lc)
{
	if(this == &lc)
		return *this;

	A = lc.A;
	B = lc.B;

	return *this;
}









PolynomialConstraint::PolynomialConstraint()
{
}

PolynomialConstraint::PolynomialConstraint(const Polynomial<Real> & polynomial, const Real & b)
{
	p = polynomial;
	p.toHornerForm(hf);
	B = b;
}

PolynomialConstraint::PolynomialConstraint(const PolynomialConstraint & pc)
{
	p = pc.p;
	hf = pc.hf;
	B = pc.B;
}

PolynomialConstraint::~PolynomialConstraint()
{
}

void PolynomialConstraint::output(std::ostream & os, const Variables & stateVars) const
{
	p.output_constraint(os, stateVars);
	os << " <= " << B << std::endl;
}

PolynomialConstraint & PolynomialConstraint::operator = (const PolynomialConstraint & pc)
{
	if(this == &pc)
		return *this;

	p = pc.p;
	hf = pc.hf;
	B = pc.B;

	return *this;
}











Constraint::Constraint()
{
}

Constraint::Constraint(const Expression_AST<Real> & exp, const Real & b)
{
	expression	= exp;
	bound		= b;
}

Constraint::Constraint(const std::string & strExpression)
{
	expression_ast_setting.clear();

	std::string prefix(str_prefix_expression_ast);
	std::string suffix(str_suffix);

	expression_ast_setting.strExpression = prefix + strExpression + suffix;

	parseExpression();

	expression_ast_setting.result.toReal(expression);
}

Constraint::Constraint(const Constraint & constraint)
{
	expression	= constraint.expression;
	bound		= constraint.bound;
}

Constraint::~Constraint()
{
}

void Constraint::output(std::ostream & os, const Variables & stateVars) const
{
	expression.output(os, stateVars);
	os << " <= " << bound;
}

Constraint & Constraint::operator = (const Constraint & constraint)
{
	if(this == &constraint)
		return *this;

	expression	= constraint.expression;
	bound		= constraint.bound;

	return *this;
}


