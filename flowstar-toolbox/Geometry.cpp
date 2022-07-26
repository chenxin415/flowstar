/*---
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.
  The code is released as is under the GNU General Public License (GPL).
---*/

#include "Geometry.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

using namespace flowstar;

// class Polyhedron

Polyhedron::Polyhedron()
{
}

Polyhedron::Polyhedron(const std::vector<LinearConstraint> & cs)
{
	constraints = cs;
}

Polyhedron::Polyhedron(const Polyhedron & P):constraints(P.constraints)
{
}
/*
Polyhedron::Polyhedron(Matrix<double> & A, Matrix<double> & b)
{
	int rows = A.rows();
	int cols = A.cols();

	for(int i=0; i<rows; ++i)
	{
		std::vector<Interval> row;

		for(int j=0; j<cols; ++j)
		{
			row.push_back(A[i][j]);
		}

		LinearConstraint lc(row, b[i][0]);
		constraints.push_back(lc);
	}
}
*/
Polyhedron::Polyhedron(const std::vector<std::vector<Real> > & A, const std::vector<Real> & B)
{
	for(int i=0; i<A.size(); ++i)
	{
		LinearConstraint lc(A[i], B[i]);
		constraints.push_back(lc);
	}
}

Polyhedron::~Polyhedron()
{
	constraints.clear();
}

Polyhedron::Polyhedron(const std::vector<std::vector<Real> > & template_matrix, const TaylorModelVec<Real> & tmv, const std::vector<Interval> & domain)
{
	for(unsigned int i=0; i<template_matrix.size(); ++i)
	{
		LinearConstraint lc(template_matrix[i], tmv.rho(template_matrix[i], domain));
		constraints.push_back(lc);
	}
}

double Polyhedron::rho(const std::vector<Real> & l) const
{
	int d = l.size();
	int n = constraints.size();
	int size = n*d;

	int *rowInd = new int[ 1 + size ];
	int *colInd = new int[ 1 + size ];
	double *coes = new double [ 1 + size ];

	glp_term_out(GLP_OFF);

	glp_prob *lp;
	lp = glp_create_prob();
	glp_set_obj_dir(lp, GLP_MAX);

	glp_add_rows(lp, n);

	for(int i=1; i<=n; ++i)
	{
		glp_set_row_bnds(lp, i, GLP_UP, 0.0, constraints[i-1].B.toDouble());
	}

	glp_add_cols(lp, d);
	for(int i=1; i<=d; ++i)
	{
		glp_set_col_bnds(lp, i, GLP_FR, 0.0, 0.0);
		glp_set_obj_coef(lp, i, l[i-1].toDouble());
	}

	for(int i=1; i<=n; ++i)
	{
		for(int j=1; j<=d; ++j)
		{
			int pos = j + (i-1)*d;
			rowInd[pos] = i;
			colInd[pos] = j;
			coes[pos] = constraints[i-1].A[j-1].toDouble();
		}
	}

	glp_load_matrix(lp, size, rowInd, colInd, coes);
	glp_simplex(lp, NULL);
	double result = glp_get_obj_val(lp);
	int status = glp_get_status(lp);

	if(status == GLP_INFEAS || status == GLP_NOFEAS)
	{
		result = INVALID;
	}
	else if(status == GLP_UNBND)
	{
		result = UNBOUNDED;
	}

	glp_delete_prob(lp);
	delete[] rowInd;
	delete[] colInd;
	delete[] coes;

	return result;
}
/*
Interval Polyhedron::rho(Matrix<Real> & l) const
{
	int d = l.cols();
	int n = constraints.size();
	int size = n*d;

	int *rowInd = new int[ 1 + size ];
	int *colInd = new int[ 1 + size ];
	double *coes = new double [ 1 + size ];

	glp_term_out(GLP_OFF);

	glp_prob *lp;
	lp = glp_create_prob();
	glp_set_obj_dir(lp, GLP_MAX);

	glp_add_rows(lp, n);

	for(int i=1; i<=n; ++i)
	{
		glp_set_row_bnds(lp, i, GLP_UP, 0.0, constraints[i-1].B.midpoint());
	}

	glp_add_cols(lp, d);
	for(int i=1; i<=d; ++i)
	{
		glp_set_col_bnds(lp, i, GLP_FR, 0.0, 0.0);

		Interval intTemp = l[0][i-1];
		glp_set_obj_coef(lp, i, intTemp.midpoint());
	}

	for(int i=1; i<=n; ++i)
	{
		for(int j=1; j<=d; ++j)
		{
			int pos = j + (i-1)*d;
			rowInd[pos] = i;
			colInd[pos] = j;
			coes[pos] = constraints[i-1].A[j-1].midpoint();
		}
	}

	glp_load_matrix(lp, size, rowInd, colInd, coes);
	glp_simplex(lp, NULL);
	double result = glp_get_obj_val(lp);
	int status = glp_get_status(lp);

	if(status == GLP_INFEAS || status == GLP_NOFEAS)
	{
		result = INVALID;
	}
	else if(status == GLP_UNBND)
	{
		result = UNBOUNDED;
	}

	glp_delete_prob(lp);
	delete[] rowInd;
	delete[] colInd;
	delete[] coes;

	Interval intTemp(result);
	return intTemp;
}
*/
void Polyhedron::tightenConstraints()
{
	for(int i=0; i<constraints.size(); ++i)
	{
		Real tmp = rho(constraints[i].A);

		if(tmp < constraints[i].B)
		{
			constraints[i].B = tmp;
		}
	}
}

bool Polyhedron::empty() const
{
	if(constraints.size() == 0)
	{
		return false;
	}
	else
	{
		int d = constraints.begin()->A.size();

		std::vector<Real> l(d);
		l[0] = 1;

		double result = this->rho(l);

		if(result <= INVALID + THRESHOLD_HIGH)
			return true;
		else
			return false;
	}
}

void Polyhedron::get(std::vector<std::vector<Real> > & A, std::vector<Real> & B) const
{
	A.clear();
	B.clear();

	for(int i=0; i<constraints.size(); ++i)
	{
		A.push_back(constraints[i].A);
		B.push_back(constraints[i].B);
	}
}

void Polyhedron::output(std::ostream & os, const Variables & stateVars) const
{
	for(int i=0; i<constraints.size(); ++i)
	{
		constraints[i].output(os, stateVars);
	}

	os << std::endl;
}

Polyhedron & Polyhedron::operator = (const Polyhedron & P)
{
	if(this == &P)
		return *this;

	constraints = P.constraints;
	return *this;
}































// class Parallelotope

Parallelotope::Parallelotope(const Matrix<Real> & A, const Matrix<Real> & B)
{
	paraTemplate = A;
	b = B;
}

Parallelotope::Parallelotope(const Parallelotope & P)
{
	paraTemplate = P.paraTemplate;
	b = P.b;
}

Parallelotope::~Parallelotope()
{
}

void Parallelotope::center(std::vector<Real> & c)
{
	int d = paraTemplate.cols();

	gsl_vector *r = gsl_vector_alloc(d);
	for(int i=0; i<d; ++i)
		gsl_vector_set( r, i, ((b[i][0] - b[i+d][0])/2).toDouble());

	// We use GSL to solve the linear equations B x = r.

	gsl_matrix *B = gsl_matrix_alloc(d, d);

	for(int i=0; i<d; ++i)
	{
		for(int j=0; j<d; ++j)
		{
			gsl_matrix_set(B, i, j, paraTemplate[i][j].toDouble());
		}
	}

	gsl_vector *x = gsl_vector_alloc(d);

	gsl_linalg_HH_solve(B, r, x);

	for(int i=0; i<d; ++i)
	{
		c[i] = gsl_vector_get(x,i);
	}

	gsl_vector_free(r);
	gsl_matrix_free(B);
	gsl_vector_free(x);
}

void Parallelotope::output(std::ostream & os)
{
	int rows = paraTemplate.rows();
	int cols = rows;
	int rangeDim = rows;

	for(int i=0; i<rows; ++i)
	{
		os << "[ ";
		for(int j=0; j<cols-1; ++j)
		{
			os << paraTemplate[i][j] << ",\t";
		}

		os << paraTemplate[i][cols-1] << " ]\t<=\t" << b[i][0] << std::endl;
	}

	for(int i=0; i<rows; ++i)
	{
		os << "[ ";
		for(int j=0; j<cols-1; ++j)
		{
			os << -paraTemplate[i][j] << ",\t";
		}

		os << -paraTemplate[i][cols-1] << " ]\t<=\t" << b[i+rangeDim][0] << std::endl;
	}
}

void Parallelotope::toTaylorModel(TaylorModelVec<Real> & result)
{
	int rangeDim = paraTemplate.rows();
	int domainDim = rangeDim + 1;

	// 1: we converse the center point to a Taylor model

	std::vector<Real> colVecCenter(rangeDim);
	center(colVecCenter);

	TaylorModelVec<Real> tmvCenter(colVecCenter, domainDim);

	// 2: we center the parallelotope at 0
	std::vector<Real> colVecDiff = paraTemplate * colVecCenter;

	// since a parallelotope is symmetric, we only need to consider half of the intercepts
	std::vector<double> new_b(rangeDim, 1);

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		new_b[i] = (b[i][0] - colVecDiff[i]).toDouble();
	}

	// 3: compute the generators.
	Matrix<double> generators(rangeDim, rangeDim);
	std::vector<int> zeroRows;	// the row indices for zero intercepts

	for(int i=0; i<rangeDim; ++i)
	{
		if(new_b[i] <= THRESHOLD_LOW && new_b[i] >= -THRESHOLD_LOW)	// zero
		{
			zeroRows.push_back(i);

			for(int j=0; j<rangeDim; ++j)
			{
				generators[i][j] = paraTemplate[i][j].toDouble();
			}
		}
		else
		{
			for(int j=0; j<rangeDim; ++j)
			{
				generators[i][j] = (paraTemplate[i][j] / new_b[i]).toDouble();
			}
		}
	}

	generators.inverse_assign();

	Matrix<Real> tmv_coefficients(rangeDim, domainDim);

	for(int j=0, k=0; j<rangeDim; ++j)
	{
		if(k < zeroRows.size() && j == zeroRows[k])	// neglect the zero length generators
		{
			++k;
		}
		else
		{
			for(int i=0; i<rangeDim; ++i)
			{
				tmv_coefficients[i][j+1] = generators[i][j];
			}
		}
	}

	TaylorModelVec<Real> tmvParallelotope(tmv_coefficients);
	result = tmvParallelotope + tmvCenter;
}

Parallelotope & Parallelotope::operator = (const Parallelotope & P)
{
	if(this == &P)
		return *this;

	paraTemplate = P.paraTemplate;
	b = P.b;
	return *this;
}









// class zonotope

Zonotope::Zonotope()
{
}

Zonotope::Zonotope(const Matrix<Real> & c, const std::list<Matrix<Real> > & G)
{
	center = c;
	generators = G;
}

Zonotope::Zonotope(const std::vector<Real> & c, const std::vector<std::vector<Real> > & G)
{
	if(c.size() > 0)
	{
		Matrix<Real> rm_c(c.size(), 1);

		for(int i=0; i<c.size(); ++i)
		{
			rm_c[i][0] = c[i];
		}

		center = rm_c;

		if(G.size() > 0)
		{
			if(G[0].size() != c.size())
			{
				printf("Creating a Zonotope: Generators should have the same dimension as the center.\n");
			}
			else
			{
				for(int i=0; i<G.size(); ++i)
				{
					Matrix<Real> rm_g(c.size(), 1);

					for(int j=0; j<G[i].size(); ++j)
					{
						rm_g[j][0] = G[i][j];
					}

					generators.push_back(rm_g);
				}
			}
		}
	}
}

Zonotope::Zonotope(const Zonotope & zonotope)
{
	center = zonotope.center;
	generators = zonotope.generators;
}

Zonotope::Zonotope(const std::vector<Interval> & box)
{
	int d = box.size();
	Matrix<Real> rm_center(d,1);
	center = rm_center;

	for(int i=0; i<d; ++i)
	{
		Matrix<Real> rm_generator(d,1);
		box[i].toCenterForm(center[i][0], rm_generator[i][0]);

		generators.push_back(rm_generator);
	}
}

Zonotope::Zonotope(Matrix<Interval> & box)
{
	int d = box.rows();
	Matrix<Real> rm_center(d,1);
	center = rm_center;

	for(int i=0; i<d; ++i)
	{
		Matrix<Real> rm_generator(d,1);
		box[i][0].toCenterForm(center[i][0], rm_generator[i][0]);

		generators.push_back(rm_generator);
	}
}

Zonotope::Zonotope(const unsigned int d)
{
	Matrix<Real> rm_zero(d, 1);
	center = rm_zero;
}

Zonotope::~Zonotope()
{
	generators.clear();
}

bool Zonotope::isEmpty() const
{
	if(center.rows() == 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

unsigned int Zonotope::numOfGen() const
{
	return generators.size();
}

void Zonotope::simplify()
{
	if(center.cols() < 1)
		return;

	std::list<Matrix<Real> > result1;

	int n = center.rows();
	int m = 10*n;		// we will select m generators

	if(generators.size() < m)
		return;

/*
	std::cout << "Before elimination:" << std::endl;
	std::list<Matrix<Real> >::iterator iter_test = generators.begin();
	for(; iter_test != generators.end(); ++iter_test)
	{
		for(int i=0; i<n; ++i)
		{
			std::cout << (*iter_test)[i][0] << ",\t";
		}

		std::cout << std::endl;
	}
*/


	for(unsigned int i=0; i<m; ++i)
	{
		std::list<Matrix<Real> >::iterator iter = generators.begin();
		std::list<Matrix<Real> >::iterator iter_max = iter;
		++iter;

		Real max = iter_max->norm(1) - iter_max->norm(0);

		for(; iter != generators.end(); ++iter)
		{
			Real norm = iter->norm(1) - iter->norm(0);

			if(max < norm)
			{
				max = norm;
				iter_max = iter;
			}
		}

		result1.push_back(*iter_max);
		generators.erase(iter_max);
	}



/*
	std::cout << "After elimination:" << std::endl;
	iter_test = generators.begin();
	for(; iter_test != generators.end(); ++iter_test)
	{
		for(int i=0; i<n; ++i)
		{
			std::cout << (*iter_test)[i][0] << ",\t";
		}

		std::cout << std::endl;
	}
*/

	std::list<Matrix<Real> > result2;
	Matrix<Real> vecZero(n, 1);

	for(unsigned int i=0; i<n; ++i)
	{
		result2.push_back(vecZero);
	}

	std::list<Matrix<Real> >::iterator iter = generators.begin();
	for(; iter != generators.end(); ++iter)
	{
		std::list<Matrix<Real> >::iterator iter2 = result2.begin();
		for(int j=0; iter2 != result2.end(); ++iter2, ++j)
		{
			(*iter2)[j][0] += (*iter)[j][0].abs();
		}
	}

	result1.splice(result1.end(), result2);
	generators = result1;

}

void Zonotope::intEval(Matrix<Interval> & range)
{
	int d = center.rows();
	range = center;

	std::list<Matrix<Real> >::iterator iter = generators.begin();
	for(; iter != generators.end(); ++iter)
	{
		for(int j=0; j<d; ++j)
		{
			double r = (*iter)[j][0].mag();
			range[j][0].bloat(r);
		}
	}
}

Zonotope & Zonotope::operator += (const Zonotope & Z)
{
	if(Z.isEmpty())
		return *this;
	else if(this->isEmpty())
	{
		*this = Z;
		return *this;
	}
	else
	{
		center += Z.center;

		std::list<Matrix<Real> >::const_iterator iter = Z.generators.begin();

		for(; iter != Z.generators.end(); ++iter)
		{
			generators.push_back(*iter);
		}

		return *this;
	}
}

Zonotope & Zonotope::operator = (const Zonotope & zonotope)
{
	if(this == &zonotope)
		return *this;

	center = zonotope.center;
	generators = zonotope.generators;

	return *this;
}

Zonotope & Zonotope::operator = (const std::vector<Interval> & box)
{
	int d = box.size();
	Matrix<Real> rm_center(d,1);
	center = rm_center;

	generators.clear();

	for(int i=0; i<d; ++i)
	{
		Matrix<Real> rm_generator(d,1);
		box[i].toCenterForm(center[i][0], rm_generator[i][0]);

		generators.push_back(rm_generator);
	}

	return *this;
}

Zonotope & Zonotope::operator = (Matrix<Interval> & box)
{
	int d = box.rows();
	Matrix<Real> rm_center(d,1);
	center = rm_center;

	generators.clear();

	for(int i=0; i<d; ++i)
	{
		Matrix<Real> rm_generator(d,1);
		box[i][0].toCenterForm(center[i][0], rm_generator[i][0]);

		generators.push_back(rm_generator);
	}

	return *this;
}

void Zonotope::toPolynomial(std::vector<Polynomial<Real> > & result)
{
	int rangeDim = center.rows();
	int domainDim = generators.size();

	Interval intZero;

	result.clear();

	for(int i=0; i<rangeDim; ++i)
	{
		std::vector<Real> coefficients;

		std::list<Matrix<Real> >::iterator iter = generators.begin();

		for(; iter != generators.end(); ++iter)
		{
			coefficients.push_back((*iter)[i][0]);
		}

		Polynomial<Real> polyTmp(coefficients);
		Polynomial<Real> constant(center[i][0], domainDim);

		result.push_back(polyTmp + constant);
	}
}



namespace flowstar
{

Zonotope operator + (const Zonotope & Z1, const Zonotope & Z2)
{
	Zonotope result = Z1;

	result += Z2;

	return result;
}

Zonotope operator * (const Matrix<Real> & A, const Zonotope & Z)
{
	if(Z.isEmpty())
		return Z;

	Zonotope result;

	result.center = A * Z.center;

	std::list<Matrix<Real> >::const_iterator iter = Z.generators.begin();
	for(; iter != Z.generators.end(); ++iter)
	{
		result.generators.push_back(A * (*iter));
	}

	return result;
}

}
