/*---
  Flow*: A Verification Tool for Cyber-Physical Systems.
  Authors: Xin Chen, Sriram Sankaranarayanan, and Erika Abraham.
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
/*
Zonotope::Zonotope(iMatrix2 & box)
{
	int d = box.rows();
	iMatrix2 im2_center(d,1);
	center = im2_center;

	for(int i=0; i<d; ++i)
	{
		center.center[i][0] = box.center[i][0];

		iMatrix2 im2_generator(d,1);
		im2_generator.center[i][0] = box.radius[i][0];

		generators.push_back(im2_generator);
	}
}
*/
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
	std::list<Matrix<Real> > result;

	int d = center.rows();
	Matrix<Real> vecZero(d, 1);

	for(int i=0; i<d; ++i)
	{
		result.push_back(vecZero);
	}

	std::list<Matrix<Real> >::iterator iter = generators.begin();
	for(; iter != generators.end(); ++iter)
	{
		std::list<Matrix<Real> >::iterator iter2 = result.begin();
		for(int j=0; iter2 != result.end(); ++iter2, ++j)
		{
			(*iter2)[j][0] += (*iter)[j][0];
		}
	}

	generators = result;
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
/*
void Zonotope::output(FILE *fp) const
{
	fprintf(fp, "Center:\n");

	center.output(fp);
	fprintf(fp, "\nGenerators:\n");

	std::list<iMatrix>::const_iterator iter = generators.begin();
	for(; iter != generators.end(); ++iter)
	{
		(*iter).output(fp);
		fprintf(fp,"\n");
	}
}
*/


Zonotope & Zonotope::operator += (const Zonotope & Z)
{
	center += Z.center;

	std::list<Matrix<Real> >::const_iterator iter = Z.generators.begin();

	for(; iter != Z.generators.end(); ++iter)
	{
		generators.push_back(*iter);
	}

	return *this;
}










/*
void Zonotope::linearTrans(Zonotope & result, const iMatrix & map) const
{
	result.center = map * center;

	result.generators.clear();

	std::list<iMatrix>::const_iterator iter = generators.begin();
	for(; iter != generators.end(); ++iter)
	{
		iMatrix im_tmp = map * (*iter);
		result.generators.push_back(im_tmp);
	}
}

void Zonotope::linearTrans_assign(const iMatrix & map)
{
	center = map * center;

	std::list<iMatrix>::iterator iter = generators.begin();
	for(; iter != generators.end(); ++iter)
	{
		iMatrix im_tmp = map * (*iter);
		(*iter) = im_tmp;
	}
}

void Zonotope::MinSum(Zonotope & result, const Zonotope & zonotope) const
{
	result.center = center + zonotope.center;
	result.generators = generators;

	std::list<iMatrix>::const_iterator iter = zonotope.generators.begin();

	for(; iter != zonotope.generators.end(); ++iter)
	{
		result.generators.push_back(*iter);
	}
}

void Zonotope::MinSum_assign(const Zonotope & zonotope)
{
	center += zonotope.center;

	std::list<iMatrix>::const_iterator iter = zonotope.generators.begin();

	for(; iter != zonotope.generators.end(); ++iter)
	{
		generators.push_back(*iter);
	}
}
*/



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

/*
Zonotope & Zonotope::operator = (iMatrix2 & box)
{
	int d = box.rows();

	iMatrix2 im2_center(d,1);

	for(int i=0; i<d; ++i)
	{
		im2_center.center[i][0] = box.center[i][0];

		iMatrix2 im2_generator(d,1);
		im2_generator.center[i][0] = box.radius[i][0];
		generators.push_back(im2_generator);
	}

	center = im2_center;

	return *this;
}
*/

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

/*
bool Zonotope::belongsto(const std::vector<double> & x)
{
	int d = generators.size();
	int n = x.size();
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
		double tmp = x[i-1] - center[i-1][0].midpoint();
		glp_set_row_bnds(lp, i, GLP_FX, tmp, 0.0);
	}

	glp_add_cols(lp, d);
	for(int i=1; i<=d; ++i)
	{
		glp_set_col_bnds(lp, i, GLP_DB, -1.0, 1.0);
		glp_set_obj_coef(lp, i, 1.0);
	}

	for(int i=1; i<=n; ++i)
	{
		std::list<iMatrix>::iterator iter = generators.begin();
		for(int j=1; j<=d; ++j)
		{
			int pos = j + (i-1)*d;
			rowInd[pos] = i;
			colInd[pos] = j;
			coes[pos] = (*iter)[i-1][0].midpoint();
			++iter;
		}
	}

	glp_load_matrix(lp, size, rowInd, colInd, coes);
	glp_simplex(lp, NULL);

	int status = glp_get_status(lp);

	glp_delete_prob(lp);
	delete[] rowInd;
	delete[] colInd;
	delete[] coes;

	if(status == GLP_OPT || status == GLP_FEAS)
	{
		return true;
	}
	else
	{
		return false;
	}
}

int Zonotope::contract(const LinearConstraint & constraint)
{
	int rangeDim = center.rows();
	int domainDim = generators.size();

	std::vector<Polynomial> polynomials;
	toPolynomial(polynomials);

	std::vector<Interval> domain;
	Interval intUnit(-1,1);
	for(int i=0; i<domainDim; ++i)
	{
		domain.push_back(intUnit);
	}

	Polynomial obj;
	for(int i=0; i<constraint.A.size(); ++i)
	{
		Polynomial polyTemp = polynomials[i];
		polyTemp.mul_assign(constraint.A[i]);
		obj += polyTemp;
	}

	HornerForm obj_hf;
	obj.toHornerForm(obj_hf);


	bool bvalid = true;
	bool bcontinue = true;
	Interval W;
	Interval intZero;

	for(; bcontinue; )
	{
		std::vector<Interval> oldDomain = domain;

		// contract the domain
		for(int i=0; i<domainDim; ++i)
		{
			Interval newInt = domain[i];

			newInt.width(W);

			// search an approximation for the lower bound
			for(; W.greaterThan(DC_THRESHOLD_SEARCH);)
			{
				Interval intLeft;
				Interval intRight;
				newInt.split(intLeft, intRight);

				std::vector<Interval> newDomain = domain;
				newDomain[i] = intLeft;

				Interval intTemp;
				obj_hf.intEval(intTemp, newDomain);

				if(intTemp.greaterThan(constraint.B))
				{
					// no intersection on the left half
					newInt = intRight;
					newInt.width(W);
					break;
				}
				else if(intTemp.smallereq(constraint.B))
				{
					// do not need to apply domain contraction w.r.t. the current constraint
					newInt = intLeft;
					newInt.width(W);
				}
				else
				{
					// refine the interval
					newInt = intLeft;
					newInt.width(W);
				}
			}


			// set the lower bound
			Interval Inf;
			newInt.inf(Inf);
			domain[i].setInf(Inf);

			newInt = domain[i];

			newInt.width(W);

			// search an approximation for the upper bound
			for(; W.greaterThan(DC_THRESHOLD_SEARCH);)
			{
				Interval intLeft;
				Interval intRight;
				newInt.split(intLeft, intRight);

				std::vector<Interval> newDomain = domain;
				newDomain[i] = intRight;

				Interval intTemp;
				obj_hf.intEval(intTemp, newDomain);

				if(intTemp.greaterThan(constraint.B))
				{
					// no intersection on the right half
					newInt = intLeft;
					newInt.width(W);
					break;
				}
				else if(intTemp.smallereq(constraint.B))
				{
					// do not need to apply domain contraction w.r.t. the current constraint
					newInt = intRight;
					newInt.width(W);
				}
				else
				{
					// refine the interval
					newInt = intRight;
					newInt.width(W);
				}
			}

			Interval Sup;
			newInt.sup(Sup);
			domain[i].setSup(Sup);	// set the upper bound

			if(!domain[i].valid())
			{
				bvalid = false;
				break;
			}
		}

		if(!bvalid)
		{
			break;
		}

		bcontinue = false;
		for(int i=0; i<domainDim; ++i)
		{
			if(oldDomain[i].widthRatio(domain[i]) <= DC_THRESHOLD_IMPROV)
			{
				bcontinue = true;
				break;
			}
		}
	}

	if(!bvalid)
	{
		return -1;
	}
	else
	{
		// compute the contracted zonotope

		for(int i=0; i<rangeDim; ++i)
		{
			std::list<iMatrix>::iterator iter = generators.begin();
			for(int j=0; j<domain.size(); ++j)
			{
				Interval M;
				domain[j].remove_midpoint(M);

				center[i][0] += (*iter)[i][0] * M;

				(*iter)[i][0] *= domain[j].sup();

				++iter;
			}
		}

		return 1;
	}
}



void Zonotope::to2DBox(iMatrix & box, const int x, const int y)
{
	std::list<iMatrix>::iterator iter = generators.begin();

	double radius_x = 0, radius_y = 0;

	for(; iter != generators.end(); ++iter)
	{
		radius_x += (*iter)[x][0].mag();
		radius_y += (*iter)[y][0].mag();
	}

	double center_x = center[x][0].midpoint();
	double center_y = center[y][0].midpoint();

	Interval int_x(center_x - radius_x, center_x + radius_x);
	box[0][0] = int_x;

	Interval int_y(center_y - radius_y, center_y + radius_y);
	box[1][0] = int_y;
}

void Zonotope::intervalRange(Interval & range, const int x)
{
	double c = center[x][0].midpoint();
	double r = 0;

	std::list<iMatrix>::iterator iter = generators.begin();
	for(; iter != generators.end(); ++iter)
	{
		r += (*iter)[x][0].mag();
	}

	range.setSup(c + r);
	range.setInf(c - r);
}

void Zonotope::plot(FILE *fp, const int x, const int y)
{



	iMatrix box(2,1);
	to2DBox(box, x, y);

	fprintf(fp,"plot( [%e,%e,%e,%e,%e] , [%e,%e,%e,%e,%e] , 'color' , 'b');\nhold on;\nclear;\n",
			box[0][0].inf(), box[0][0].sup(), box[0][0].sup(), box[0][0].inf(), box[0][0].inf(), box[1][0].inf(),
			box[1][0].inf(), box[1][0].sup(), box[1][0].sup(), box[1][0].inf());

}
*/

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
