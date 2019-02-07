/*---
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.
  The code is released as is under the GNU General Public License (GPL).
---*/

#ifndef TAYLORMODEL_H_
#define TAYLORMODEL_H_

#include "Polynomial.h"
#include "settings.h"


namespace flowstar
{


inline void exp_taylor_remainder(Interval & result, const Interval & tmRange, const unsigned int order, const Global_Computation_Setting & setting)
{
	Interval intProd = tmRange.pow(order);

	Interval J(0,1);
	J *= tmRange;
	J.exp_assign();

	result = setting.factorial_rec[order] * intProd * J;
}

inline void rec_taylor_remainder(Interval & result, const Interval & tmRange, const unsigned int order, const Global_Computation_Setting & setting)
{
	Interval J(0,1);
	J *= tmRange;
	J += 1;
	J.rec_assign();

	Interval intProd = J;
	intProd *= tmRange;
	intProd *= -1;

	result = intProd.pow(order);
	result *= J;
}

inline void sin_taylor_remainder(Interval & result, const Interval & C, const Interval & tmRange, const unsigned int order, const Global_Computation_Setting & setting)
{
	Interval intProd = tmRange.pow(order);

	Interval J(0,1);
	J *= tmRange;
	J += C;

	int k = order % 4;

	switch(k)
	{
	case 0:
		J.sin_assign();
		break;
	case 1:
		J.cos_assign();
		break;
	case 2:
		J.sin_assign();
		J.inv_assign();
		break;
	case 3:
		J.cos_assign();
		J.inv_assign();
		break;
	}

	result = setting.factorial_rec[order] * intProd * J;
}

inline void cos_taylor_remainder(Interval & result, const Interval & C, const Interval & tmRange, const unsigned int order, const Global_Computation_Setting & setting)
{
	Interval intProd = tmRange.pow(order);

	Interval J(0,1);
	J *= tmRange;
	J += C;

	int k = order % 4;

	switch(k)
	{
	case 0:
		J.cos_assign();
		break;
	case 1:
		J.sin_assign();
		J.inv_assign();
		break;
	case 2:
		J.cos_assign();
		J.inv_assign();
		break;
	case 3:
		J.sin_assign();
		break;
	}

	result = setting.factorial_rec[order] * intProd * J;
}

inline void log_taylor_remainder(Interval & result, const Interval & tmRange, const int order)
{
	Interval J(0,1);
	J *= tmRange;
	J += 1;
	J.rec_assign();

	Interval I = tmRange;
	I *= J;

	result = I.pow(order);

	result /= order;

	if((order+1)%2 == 1)		// order+1 is odd
	{
		result *= -1;
	}
}

inline void sqrt_taylor_remainder(Interval & result, const Interval & tmRange, const int order, const Global_Computation_Setting & setting)
{
	Interval I(0,1);
	I *= tmRange;
	I += 1;
	I.rec_assign();

	Interval intTemp;
	I.sqrt(intTemp);

	I *= tmRange;
	I /= 2;

	Interval intProd = I.pow(order-1);

	intProd /= intTemp;
	intProd *= tmRange;
	intProd /= 2;

	result = setting.double_factorial[2*order-3] * setting.factorial_rec[order] * intProd;

	if(order % 2 == 0)
	{
		result *= -1;
	}
}








template <class DATA_TYPE>
class HornerForm;

template <class DATA_TYPE>
class TaylorModelVec;

template <class DATA_TYPE>
class Expression_AST;



template <class DATA_TYPE>
class TaylorModel			// Taylor models: R^n -> R. We use t to denote the time variable and x to denote the state variable.
{
public:
	Polynomial<DATA_TYPE> expansion;	// Taylor expansion
	Interval remainder;					// remainder interval

public:
	TaylorModel();
	TaylorModel(const DATA_TYPE & c, const unsigned int numVars);				// constant
	TaylorModel(const Polynomial<DATA_TYPE> & polyExp);
	TaylorModel(const Polynomial<DATA_TYPE> & polyExp, const Interval & I);		// Taylor model (P,I)

	TaylorModel(Matrix<DATA_TYPE> & coefficients);
	TaylorModel(Matrix<DATA_TYPE> & coefficients, const Interval & I);
	TaylorModel(const std::vector<DATA_TYPE> & coefficients);
	TaylorModel(const std::vector<DATA_TYPE> & coefficients, const Interval & I);

	TaylorModel(const DATA_TYPE *pCoefficients, const unsigned int numVars);
	TaylorModel(const DATA_TYPE *pCoefficients, const unsigned int numVars, const Interval & I);

	TaylorModel(const UnivariateTaylorModel<DATA_TYPE> & utm, const unsigned int numVars, const bool dummy);
	TaylorModel(const TaylorModel<DATA_TYPE> & tm);
	~TaylorModel();

//	TaylorModel(const std::string & strPolynomial, const Variables & vars);
//	TaylorModel(const std::string & strPolynomial, const Interval & rem, const Variables & vars);

	void clear();
	void output(std::ostream & os, const Variables & vars) const;

	void constant(DATA_TYPE & c) const;									// Return the constant part of the expansion.

	template <class DATA_TYPE2>
	void intEval(Interval & result, const std::vector<DATA_TYPE2> & domain) const;

	template <class DATA_TYPE2>
	void intEvalNormal(Interval & result, const std::vector<DATA_TYPE2> & step_exp_table) const;

	template <class DATA_TYPE2>
	void ctrunc(const std::vector<DATA_TYPE2> & domain, const unsigned int order);

	void nctrunc(const unsigned int order);

	template <class DATA_TYPE2>
	void ctrunc_normal(const std::vector<DATA_TYPE2> & step_exp_table, const unsigned int order);

	TaylorModel<DATA_TYPE> & operator = (const TaylorModel<DATA_TYPE> & tm);
	TaylorModel<DATA_TYPE> & operator = (const Polynomial<DATA_TYPE> & p);

	TaylorModel<DATA_TYPE> & operator += (const TaylorModel<DATA_TYPE> & tm);
	TaylorModel<DATA_TYPE> & operator += (const Polynomial<DATA_TYPE> & p);
	TaylorModel<DATA_TYPE> & operator -= (const TaylorModel<DATA_TYPE> & tm);
	TaylorModel<DATA_TYPE> & operator -= (const Polynomial<DATA_TYPE> & p);

	TaylorModel<DATA_TYPE> operator + (const TaylorModel<DATA_TYPE> & tm) const;
	TaylorModel<DATA_TYPE> operator + (const Polynomial<DATA_TYPE> & p) const;
	TaylorModel<DATA_TYPE> operator - (const TaylorModel<DATA_TYPE> & tm) const;
	TaylorModel<DATA_TYPE> operator - (const Polynomial<DATA_TYPE> & p) const;

	TaylorModel<DATA_TYPE> & operator *= (const DATA_TYPE & c);
	TaylorModel<DATA_TYPE> & operator /= (const DATA_TYPE & c);
	TaylorModel<DATA_TYPE> operator * (const DATA_TYPE & c) const;
	TaylorModel<DATA_TYPE> operator / (const DATA_TYPE & c) const;


	template <class DATA_TYPE2>
	void mul_ctrunc(TaylorModel<DATA_TYPE> & result, const TaylorModel<DATA_TYPE> & tm, const std::vector<DATA_TYPE2> & domain, const unsigned int order, const Interval & cutoff_threshold) const;

	template <class DATA_TYPE2>
	void mul_ctrunc_assign(const TaylorModel<DATA_TYPE> & tm, const std::vector<DATA_TYPE2> & domain, const unsigned int order, const Interval & cutoff_threshold);

	template <class DATA_TYPE2>
	void mul_ctrunc_normal(TaylorModel<DATA_TYPE> & result, const TaylorModel<DATA_TYPE> & tm, const std::vector<DATA_TYPE2> & step_exp_table, const unsigned int order, const Interval & cutoff_threshold) const;

	template <class DATA_TYPE2>
	void mul_ctrunc_normal_assign(const TaylorModel<DATA_TYPE> & tm, const std::vector<DATA_TYPE2> & step_exp_table, const unsigned int order, const Interval & cutoff_threshold);

	void mul_no_remainder(TaylorModel<DATA_TYPE> & result, const TaylorModel<DATA_TYPE> & tm, const unsigned int order, const Interval & cutoff_threshold) const;
//	void mul_no_remainder_no_cutoff(TaylorModel & result, const TaylorModel & tm, const int order) const;


//	void mul_ctrunc_assign(const TaylorModel & tm, const std::vector<Interval> & domain, const int order, const Interval & cutoff_threshold);
//	void mul_ctrunc_normal_assign(const TaylorModel & tm, const std::vector<Interval> & step_exp_table, const int order, const Interval & cutoff_threshold);
	void mul_no_remainder_assign(const TaylorModel<DATA_TYPE> & tm, const unsigned int order, const Interval & cutoff_threshold);
//	void mul_no_remainder_no_cutoff_assign(const TaylorModel & tm, const int order);


//	void mul_insert(TaylorModel & result, const TaylorModel & tm, const Interval & tmPolyRange, const std::vector<Interval> & domain, const Interval & cutoff_threshold) const;
//	void mul_insert_normal(TaylorModel & result, const TaylorModel & tm, const Interval & tmPolyRange, const std::vector<Interval> & step_exp_table, const Interval & cutoff_threshold) const;
	template <class DATA_TYPE2>
	void mul_insert_ctrunc(TaylorModel<DATA_TYPE> & result, const TaylorModel<DATA_TYPE> & tm, const DATA_TYPE2 & tmPolyRange, const std::vector<DATA_TYPE2> & domain, const unsigned int order, const Interval & cutoff_threshold) const;

	template <class DATA_TYPE2, class DATA_TYPE3>
	void mul_insert_ctrunc_normal(TaylorModel<DATA_TYPE> & result, const TaylorModel<DATA_TYPE> & tm, const DATA_TYPE2 & tmPolyRange, const std::vector<DATA_TYPE3> & step_exp_table, const unsigned int order, const Interval & cutoff_threshold) const;

//	void mul_insert_ctrunc_normal_no_cutoff(TaylorModel & result, const TaylorModel & tm, const Interval & tmPolyRange, const std::vector<Interval> & step_exp_table, const int order) const;
	template <class DATA_TYPE2, class DATA_TYPE3>
	void mul_insert_ctrunc_normal(TaylorModel<DATA_TYPE> & result, Interval & tm1, Interval & intTrunc, const TaylorModel<DATA_TYPE> & tm, const DATA_TYPE2 & tmPolyRange, const std::vector<DATA_TYPE3> & step_exp_table, const unsigned int order, const Interval & cutoff_threshold) const;
//	void mul_insert_ctrunc_normal_no_cutoff(TaylorModel & result, Interval & tm1, Interval & intTrunc, const TaylorModel & tm, const Interval & tmPolyRange, const std::vector<Interval> & step_exp_table, const int order) const;


//	void mul_insert_assign(const TaylorModel & tm, const Interval & tmPolyRange, const std::vector<Interval> & domain, const Interval & cutoff_threshold);
//	void mul_insert_normal_assign(const TaylorModel & tm, const Interval & tmPolyRange, const std::vector<Interval> & step_exp_table, const Interval & cutoff_threshold);

	template <class DATA_TYPE2>
	void mul_insert_ctrunc_assign(const TaylorModel<DATA_TYPE> & tm, const DATA_TYPE2 & tmPolyRange, const std::vector<DATA_TYPE2> & domain, const unsigned int order, const Interval & cutoff_threshold);

	template <class DATA_TYPE2, class DATA_TYPE3>
	void mul_insert_ctrunc_normal_assign(const TaylorModel<DATA_TYPE> & tm, const DATA_TYPE2 & tmPolyRange, const std::vector<DATA_TYPE3> & step_exp_table, const unsigned int order, const Interval & cutoff_threshold);
//	void mul_insert_ctrunc_normal_no_cutoff_assign(const TaylorModel & tm, const Interval & tmPolyRange, const std::vector<Interval> & step_exp_table, const int order);

	template <class DATA_TYPE2, class DATA_TYPE3>
	void mul_insert_ctrunc_normal_assign(Interval & tm1, Interval & intTrunc, const TaylorModel<DATA_TYPE> & tm, const DATA_TYPE2 & tmPolyRange, const std::vector<DATA_TYPE3> & step_exp_table, const unsigned int order, const Interval & cutoff_threshold);
//	void mul_insert_ctrunc_normal_no_cutoff_assign(Interval & tm1, Interval & intTrunc, const TaylorModel & tm, const Interval & tmPolyRange, const std::vector<Interval> & step_exp_table, const int order);


	void derivative(TaylorModel<DATA_TYPE> & result, const unsigned int varIndex) const;		// derivative with respect to a variable

	// Lie derivative, the vector field is given by f
	void LieDerivative(TaylorModel<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & f, const unsigned int order, const Interval & cutoff_threshold) const;

	void integral_time(TaylorModel<DATA_TYPE> & result, const Interval & I) const;				// Integral with respect to t
	void integral_time(TaylorModel<DATA_TYPE> & result) const;

	void linearCoefficients(Matrix<DATA_TYPE> & coefficients, const unsigned int row) const;
	void linearCoefficients(std::vector<DATA_TYPE> & coefficients) const;


	void toHornerForm(HornerForm<DATA_TYPE> & hf, Interval & I) const;


//	void insert(TaylorModel & result, const TaylorModelVec & vars, const std::vector<Interval> & varsPolyRange, const std::vector<Interval> & domain, const Interval & cutoff_threshold) const;
//	void insert_normal(TaylorModel & result, const TaylorModelVec & vars, const std::vector<Interval> & varsPolyRange, const std::vector<Interval> & step_exp_table, const int numVars, const Interval & cutoff_threshold) const;

	template <class DATA_TYPE2>
	void insert_ctrunc(TaylorModel<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & vars, const std::vector<DATA_TYPE2> & varsPolyRange, const std::vector<DATA_TYPE2> & domain, const unsigned int order, const Interval & cutoff_threshold) const;

	void insert_no_remainder(TaylorModel<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & vars, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold) const;
//	void insert_no_remainder_no_cutoff(TaylorModel & result, const TaylorModelVec & vars, const int numVars, const int order) const;

	template <class DATA_TYPE2, class DATA_TYPE3>
	void insert_ctrunc_normal(TaylorModel<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & vars, const std::vector<DATA_TYPE2> & varsPolyRange, const std::vector<DATA_TYPE3> & step_exp_table, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold) const;

//	void insert_ctrunc_normal_no_cutoff(TaylorModel & result, const TaylorModelVec & vars, const std::vector<Interval> & varsPolyRange, const std::vector<Interval> & step_exp_table, const int numVars, const int order) const;

	template <class DATA_TYPE2>
	void evaluate_time(TaylorModel<DATA_TYPE> & result, const std::vector<DATA_TYPE2> & step_exp_table) const;			// evaluate the Taylor model at time t

	void mul(TaylorModel<DATA_TYPE> & result, const unsigned int varIndex, const unsigned int degree) const;		// multiplied by a term x^d
	void mul_assign(const unsigned int varIndex, const unsigned int degree);

	void mul_assign(const unsigned int varIndex, const Interval & range);

	void rmConstant();
	void decompose(TaylorModel<DATA_TYPE> & linear, TaylorModel<DATA_TYPE> & other) const;

	template <class DATA_TYPE2>
	void cutoff_normal(const std::vector<DATA_TYPE2> & step_exp_table, const Interval & cutoff_threshold);

	template <class DATA_TYPE2>
	void cutoff(const std::vector<DATA_TYPE2> & domain, const Interval & cutoff_threshold);

	void cutoff(const Interval & cutoff_threshold);

	unsigned int degree() const;
	bool isZero() const;

//	void center_nc();

	void rmZeroTerms(const std::vector<unsigned int> & indices);

	void normalize(std::vector<Interval> & domain, const Interval & cutoff_threshold);

	template <class DATA_TYPE2>
	void polyRange(Interval & result, const std::vector<DATA_TYPE2> & domain) const;

	template <class DATA_TYPE2>
	void polyRangeNormal(Interval & result, const std::vector<DATA_TYPE2> & step_exp_table) const;

	Interval getRemainder() const;
	void getExpansion(Polynomial<DATA_TYPE> & p) const;

	void exp_taylor(TaylorModel<DATA_TYPE> & result, std::list<Interval> & ranges, const std::vector<Interval> & step_exp_table, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold, const Global_Computation_Setting & setting) const;
	void rec_taylor(TaylorModel<DATA_TYPE> & result, std::list<Interval> & ranges, const std::vector<Interval> & step_exp_table, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold, const Global_Computation_Setting & setting) const;
	void sin_taylor(TaylorModel<DATA_TYPE> & result, std::list<Interval> & ranges, const std::vector<Interval> & step_exp_table, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold, const Global_Computation_Setting & setting) const;
	void cos_taylor(TaylorModel<DATA_TYPE> & result, std::list<Interval> & ranges, const std::vector<Interval> & step_exp_table, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold, const Global_Computation_Setting & setting) const;
	void log_taylor(TaylorModel<DATA_TYPE> & result, std::list<Interval> & ranges, const std::vector<Interval> & step_exp_table, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold) const;
	void sqrt_taylor(TaylorModel<DATA_TYPE> & result, std::list<Interval> & ranges, const std::vector<Interval> & step_exp_table, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold, const Global_Computation_Setting & setting) const;

	void exp_taylor(TaylorModel<DATA_TYPE> & result, const std::vector<Interval> & domain, const unsigned int order, const Interval & cutoff_threshold, const Global_Computation_Setting & setting) const;
	void rec_taylor(TaylorModel<DATA_TYPE> & result, const std::vector<Interval> & domain, const unsigned int order, const Interval & cutoff_threshold, const Global_Computation_Setting & setting) const;
	void sin_taylor(TaylorModel<DATA_TYPE> & result, const std::vector<Interval> & domain, const unsigned int order, const Interval & cutoff_threshold, const Global_Computation_Setting & setting) const;
	void cos_taylor(TaylorModel<DATA_TYPE> & result, const std::vector<Interval> & domain, const unsigned int order, const Interval & cutoff_threshold, const Global_Computation_Setting & setting) const;
	void log_taylor(TaylorModel<DATA_TYPE> & result, const std::vector<Interval> & domain, const unsigned int order, const Interval & cutoff_threshold) const;
	void sqrt_taylor(TaylorModel<DATA_TYPE> & result, const std::vector<Interval> & domain, const unsigned int order, const Interval & cutoff_threshold, const Global_Computation_Setting & setting) const;




/*
	void exp_taylor(TaylorModel & result, const std::vector<Interval> & step_exp_table, const int numVars, const int order, const Interval & cutoff_threshold) const;
	void rec_taylor(TaylorModel & result, const std::vector<Interval> & step_exp_table, const int numVars, const int order, const Interval & cutoff_threshold) const;
	void sin_taylor(TaylorModel & result, const std::vector<Interval> & step_exp_table, const int numVars, const int order, const Interval & cutoff_threshold) const;
	void cos_taylor(TaylorModel & result, const std::vector<Interval> & step_exp_table, const int numVars, const int order, const Interval & cutoff_threshold) const;
	void log_taylor(TaylorModel & result, const std::vector<Interval> & step_exp_table, const int numVars, const int order, const Interval & cutoff_threshold) const;
	void sqrt_taylor(TaylorModel & result, const std::vector<Interval> & step_exp_table, const int numVars, const int order, const Interval & cutoff_threshold) const;
*/
/*
	// ================== API ==================
	void exp_taylor(TaylorModel & result, const int order, const Taylor_Model_Computation_Setting & setting) const;
	void rec_taylor(TaylorModel & result, const int order, const Taylor_Model_Computation_Setting & setting) const;
	void sin_taylor(TaylorModel & result, const int order, const Taylor_Model_Computation_Setting & setting) const;
	void cos_taylor(TaylorModel & result, const int order, const Taylor_Model_Computation_Setting & setting) const;
	void log_taylor(TaylorModel & result, const int order, const Taylor_Model_Computation_Setting & setting) const;
	void sqrt_taylor(TaylorModel & result, const int order, const Taylor_Model_Computation_Setting & setting) const;

	void exp_taylor(TaylorModel & result, const int order, const std::vector<Interval> & domain, const Interval & cutoff_threshold) const;
	void rec_taylor(TaylorModel & result, const int order, const std::vector<Interval> & domain, const Interval & cutoff_threshold) const;
	void sin_taylor(TaylorModel & result, const int order, const std::vector<Interval> & domain, const Interval & cutoff_threshold) const;
	void cos_taylor(TaylorModel & result, const int order, const std::vector<Interval> & domain, const Interval & cutoff_threshold) const;
	void log_taylor(TaylorModel & result, const int order, const std::vector<Interval> & domain, const Interval & cutoff_threshold) const;
	void sqrt_taylor(TaylorModel & result, const int order, const std::vector<Interval> & domain, const Interval & cutoff_threshold) const;


	void add(TaylorModel & result, const TaylorModel & tm) const;			// addition
	void sub(TaylorModel & result, const TaylorModel & tm) const;			// subtraction
	void mul(TaylorModel & result, const TaylorModel & tm, const int order, const Taylor_Model_Computation_Setting & setting) const;
	void intEval(Interval & result, const Taylor_Model_Computation_Setting & setting) const;
	// =========================================
*/


/*
	void extend(const int num);
	void extend();

	void substitute(TaylorModel & result, const std::vector<int> & varIDs, const std::vector<Interval> & intVals) const;
	void substitute_with_precond(const std::vector<bool> & substitution, const std::vector<Interval> & step_exp_table);
	void substitute_with_precond_no_remainder(const std::vector<bool> & substitution);
*/


	template <class DATA_TYPE2>
	friend class HornerForm;

	template <class DATA_TYPE2>
	friend class Polynomial;

	template <class DATA_TYPE2>
	friend class TaylorModelVec;

	friend class Flowpipe;
//	friend class ContinuousSystem;
//	friend class ContinuousReachability;
//	friend class HybridSystem;
//	friend class HybridReachability;
};



template <class DATA_TYPE>
TaylorModel<DATA_TYPE>::TaylorModel()
{
}

template <class DATA_TYPE>
TaylorModel<DATA_TYPE>::TaylorModel(const DATA_TYPE & c, const unsigned int numVars)
{
	Polynomial<DATA_TYPE> p(c, numVars);
	expansion = p;
}

template <class DATA_TYPE>
TaylorModel<DATA_TYPE>::TaylorModel(const Polynomial<DATA_TYPE> & polyExp)
{
	expansion = polyExp;
}

template <class DATA_TYPE>
TaylorModel<DATA_TYPE>::TaylorModel(const Polynomial<DATA_TYPE> & polyExp, const Interval & I)
{
	expansion = polyExp;
	remainder = I;
}

template <class DATA_TYPE>
TaylorModel<DATA_TYPE>::TaylorModel(Matrix<DATA_TYPE> & coefficients)
{
	Polynomial<DATA_TYPE> p(coefficients);
	expansion = p;
}

template <class DATA_TYPE>
TaylorModel<DATA_TYPE>::TaylorModel(Matrix<DATA_TYPE> & coefficients, const Interval & I)
{
	Polynomial<DATA_TYPE> p(coefficients);
	expansion = p;
	remainder = I;
}

template <class DATA_TYPE>
TaylorModel<DATA_TYPE>::TaylorModel(const std::vector<DATA_TYPE> & coefficients)
{
	Polynomial<DATA_TYPE> p(coefficients);
	expansion = p;
}

template <class DATA_TYPE>
TaylorModel<DATA_TYPE>::TaylorModel(const std::vector<DATA_TYPE> & coefficients, const Interval & I)
{
	Polynomial<DATA_TYPE> p(coefficients);
	expansion = p;
	remainder = I;
}

template <class DATA_TYPE>
TaylorModel<DATA_TYPE>::TaylorModel(const DATA_TYPE *pCoefficients, const unsigned int numVars)
{
	Polynomial<DATA_TYPE> p(pCoefficients, numVars);
	expansion = p;
}

template <class DATA_TYPE>
TaylorModel<DATA_TYPE>::TaylorModel(const DATA_TYPE *pCoefficients, const unsigned int numVars, const Interval & I)
{
	Polynomial<DATA_TYPE> p(pCoefficients, numVars);
	expansion = p;
	remainder = I;
}

template <class DATA_TYPE>
TaylorModel<DATA_TYPE>::TaylorModel(const UnivariateTaylorModel<DATA_TYPE> & utm, const unsigned int numVars, const bool dummy)
{
	std::vector<unsigned int> degrees(numVars);

	for(unsigned int i=0; i<utm.expansion.coefficients.size(); ++i)
	{
		if(utm.expansion.coefficients[i] != 0)
		{
			Term<DATA_TYPE> term(utm.expansion.coefficients[i], degrees);
			term.degrees[0] = i;
			term.d = i;

			expansion.terms.push_back(term);
		}
	}

	remainder = utm.remainder;
}

template <class DATA_TYPE>
TaylorModel<DATA_TYPE>::TaylorModel(const TaylorModel<DATA_TYPE> & tm)
{
	expansion = tm.expansion;
	remainder = tm.remainder;
}

template <class DATA_TYPE>
TaylorModel<DATA_TYPE>::~TaylorModel()
{
}

template <class DATA_TYPE>
void TaylorModel<DATA_TYPE>::clear()
{
	expansion.clear();
	remainder = 0;
}

template <class DATA_TYPE>
void TaylorModel<DATA_TYPE>::output(std::ostream & os, const Variables & vars) const
{
	expansion.output(os, vars);
	os << " + " << remainder;
}

template <class DATA_TYPE>
void TaylorModel<DATA_TYPE>::constant(DATA_TYPE & c) const
{
	expansion.constant(c);
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModel<DATA_TYPE>::intEval(Interval & result, const std::vector<DATA_TYPE2> & domain) const
{
	expansion.evaluate(result, domain);
	result += remainder;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModel<DATA_TYPE>::intEvalNormal(Interval & result, const std::vector<DATA_TYPE2> & step_exp_table) const
{
	expansion.intEvalNormal(result, step_exp_table);
	result += remainder;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModel<DATA_TYPE>::ctrunc(const std::vector<DATA_TYPE2> & domain, const unsigned int order)
{
	Interval I;
	expansion.ctrunc(I, domain, order);
	remainder += I;
}

template <class DATA_TYPE>
void TaylorModel<DATA_TYPE>::nctrunc(const unsigned int order)
{
	expansion.nctrunc(order);
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModel<DATA_TYPE>::ctrunc_normal(const std::vector<DATA_TYPE2> & step_exp_table, const unsigned int order)
{
	Interval I;
	expansion.ctrunc_normal(I, step_exp_table, order);
	remainder += I;
}

template <class DATA_TYPE>
TaylorModel<DATA_TYPE> & TaylorModel<DATA_TYPE>::operator = (const TaylorModel<DATA_TYPE> & tm)
{
	if(this == &tm)
		return *this;

	expansion = tm.expansion;
	remainder = tm.remainder;

	return *this;
}

template <class DATA_TYPE>
TaylorModel<DATA_TYPE> & TaylorModel<DATA_TYPE>::operator = (const Polynomial<DATA_TYPE> & p)
{
	expansion = p;
	remainder = 0;

	return *this;
}

template <class DATA_TYPE>
TaylorModel<DATA_TYPE> & TaylorModel<DATA_TYPE>::operator += (const TaylorModel<DATA_TYPE> & tm)
{
	expansion += tm.expansion;
	remainder += tm.remainder;

	return *this;
}

template <class DATA_TYPE>
TaylorModel<DATA_TYPE> & TaylorModel<DATA_TYPE>::operator += (const Polynomial<DATA_TYPE> & p)
{
	expansion += p;

	return *this;
}

template <class DATA_TYPE>
TaylorModel<DATA_TYPE> & TaylorModel<DATA_TYPE>::operator -= (const TaylorModel<DATA_TYPE> & tm)
{
	expansion -= tm.expansion;
	remainder -= tm.remainder;

	return *this;
}

template <class DATA_TYPE>
TaylorModel<DATA_TYPE> & TaylorModel<DATA_TYPE>::operator -= (const Polynomial<DATA_TYPE> & p)
{
	expansion -= p;

	return *this;
}

template <class DATA_TYPE>
TaylorModel<DATA_TYPE> TaylorModel<DATA_TYPE>::operator + (const TaylorModel<DATA_TYPE> & tm) const
{
	TaylorModel<DATA_TYPE> result = *this;
	result += tm;

	return result;
}

template <class DATA_TYPE>
TaylorModel<DATA_TYPE> TaylorModel<DATA_TYPE>::operator + (const Polynomial<DATA_TYPE> & p) const
{
	TaylorModel<DATA_TYPE> result = *this;
	result += p;

	return result;
}

template <class DATA_TYPE>
TaylorModel<DATA_TYPE> TaylorModel<DATA_TYPE>::operator - (const TaylorModel<DATA_TYPE> & tm) const
{
	TaylorModel<DATA_TYPE> result = *this;
	result -= tm;

	return result;
}

template <class DATA_TYPE>
TaylorModel<DATA_TYPE> TaylorModel<DATA_TYPE>::operator - (const Polynomial<DATA_TYPE> & p) const
{
	TaylorModel<DATA_TYPE> result = *this;
	result -= p;

	return result;
}

template <class DATA_TYPE>
TaylorModel<DATA_TYPE> & TaylorModel<DATA_TYPE>::operator *= (const DATA_TYPE & c)
{
	expansion *= c;
	remainder *= c;

	return *this;
}

template <class DATA_TYPE>
TaylorModel<DATA_TYPE> & TaylorModel<DATA_TYPE>::operator /= (const DATA_TYPE & c)
{
	expansion /= c;
	remainder /= c;

	return *this;
}

template <class DATA_TYPE>
TaylorModel<DATA_TYPE> TaylorModel<DATA_TYPE>::operator * (const DATA_TYPE & c) const
{
	TaylorModel<DATA_TYPE> result = *this;
	result *= c;

	return result;
}

template <class DATA_TYPE>
TaylorModel<DATA_TYPE> TaylorModel<DATA_TYPE>::operator / (const DATA_TYPE & c) const
{
	TaylorModel<DATA_TYPE> result = *this;
	result /= c;

	return result;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModel<DATA_TYPE>::mul_ctrunc(TaylorModel<DATA_TYPE> & result, const TaylorModel<DATA_TYPE> & tm, const std::vector<DATA_TYPE2> & domain, const unsigned int order, const Interval & cutoff_threshold) const
{
	Polynomial<DATA_TYPE> P1xP2;
	Interval P1xI2, P2xI1, I1xI2;

	P1xP2 = expansion * tm.expansion;

	if(tm.remainder != 0)
	{
		expansion.evaluate(P1xI2, domain);
		P1xI2 *= tm.remainder;
	}

	if(remainder != 0)
	{
		tm.expansion.evaluate(P2xI1, domain);
		P2xI1 *= remainder;
	}

	I1xI2 = remainder * tm.remainder;

	result.expansion = P1xP2;

	result.remainder = I1xI2;
	result.remainder += P2xI1;
	result.remainder += P1xI2;

	result.ctrunc(domain, order);
	result.cutoff(domain, cutoff_threshold);
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModel<DATA_TYPE>::mul_ctrunc_assign(const TaylorModel<DATA_TYPE> & tm, const std::vector<DATA_TYPE2> & domain, const unsigned int order, const Interval & cutoff_threshold)
{
	Polynomial<DATA_TYPE> P1xP2;
	Interval P1xI2, P2xI1, I1xI2;

	P1xP2 = expansion * tm.expansion;

	if(tm.remainder != 0)
	{
		expansion.evaluate(P1xI2, domain);
		P1xI2 *= tm.remainder;
	}

	if(remainder != 0)
	{
		tm.expansion.evaluate(P2xI1, domain);
		P2xI1 *= remainder;
	}

	I1xI2 = remainder * tm.remainder;

	expansion = P1xP2;
	remainder = I1xI2 + P2xI1 + P1xI2;

	ctrunc(domain, order);
	cutoff(domain, cutoff_threshold);
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModel<DATA_TYPE>::mul_ctrunc_normal(TaylorModel<DATA_TYPE> & result, const TaylorModel<DATA_TYPE> & tm, const std::vector<DATA_TYPE2> & step_exp_table, const unsigned int order, const Interval & cutoff_threshold) const
{
	Polynomial<DATA_TYPE> P1xP2;
	Interval P1xI2, P2xI1, I1xI2;

	P1xP2 = expansion * tm.expansion;

	if(tm.remainder != 0)
	{
		expansion.evaluate_normal(P1xI2, step_exp_table);
		P1xI2 *= tm.remainder;
	}

	if(remainder != 0)
	{
		tm.expansion.evaluate_normal(P2xI1, step_exp_table);
		P2xI1 *= remainder;
	}

	I1xI2 = remainder * tm.remainder;

	result.expansion = P1xP2;

	result.remainder = I1xI2;
	result.remainder += P2xI1;
	result.remainder += P1xI2;

	result.ctrunc_normal(step_exp_table, order);
	result.cutoff_normal(step_exp_table, cutoff_threshold);
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModel<DATA_TYPE>::mul_ctrunc_normal_assign(const TaylorModel<DATA_TYPE> & tm, const std::vector<DATA_TYPE2> & step_exp_table, const unsigned int order, const Interval & cutoff_threshold)
{
	Polynomial<DATA_TYPE> P1xP2;
	Interval P1xI2, P2xI1, I1xI2;

	P1xP2 = expansion * tm.expansion;

	if(tm.remainder != 0)
	{
		expansion.evaluate_normal(P1xI2, step_exp_table);
		P1xI2 *= tm.remainder;
	}

	if(remainder != 0)
	{
		tm.expansion.evaluate_normal(P2xI1, step_exp_table);
		P2xI1 *= remainder;
	}

	I1xI2 = remainder * tm.remainder;

	expansion = P1xP2;
	remainder = I1xI2 + P2xI1 + P1xI2;

	ctrunc_normal(step_exp_table, order);
	cutoff_normal(step_exp_table, cutoff_threshold);
}


template <class DATA_TYPE>
void TaylorModel<DATA_TYPE>::mul_no_remainder(TaylorModel<DATA_TYPE> & result, const TaylorModel<DATA_TYPE> & tm, const unsigned int order, const Interval & cutoff_threshold) const
{
	result.expansion = expansion * tm.expansion;
	result.expansion.nctrunc(order);
	result.expansion.cutoff(cutoff_threshold);
	result.remainder = 0;
}

template <class DATA_TYPE>
void TaylorModel<DATA_TYPE>::mul_no_remainder_assign(const TaylorModel<DATA_TYPE> & tm, const unsigned int order, const Interval & cutoff_threshold)
{
	expansion *= tm.expansion;
	expansion.nctrunc(order);
	expansion.cutoff(cutoff_threshold);
	remainder = 0;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModel<DATA_TYPE>::mul_insert_ctrunc(TaylorModel<DATA_TYPE> & result, const TaylorModel<DATA_TYPE> & tm, const DATA_TYPE2 & tmPolyRange, const std::vector<DATA_TYPE2> & domain, const unsigned int order, const Interval & cutoff_threshold) const
{
	Polynomial<DATA_TYPE> P1xP2;
	Interval P1xI2, P2xI1, I1xI2;

	P1xP2 = expansion * tm.expansion;

	if(tm.remainder != 0)
	{
		expansion.evaluate(P1xI2, domain);
		P1xI2 *= tm.remainder;
	}

	if(remainder != 0)
	{
		P2xI1 = tmPolyRange * remainder;
	}

	I1xI2 = remainder * tm.remainder;

	result.expansion = P1xP2;

	result.remainder = I1xI2;
	result.remainder += P2xI1;
	result.remainder += P1xI2;

	result.ctrunc(domain, order);
	result.cutoff(domain, cutoff_threshold);
}

template <class DATA_TYPE>
template <class DATA_TYPE2, class DATA_TYPE3>
void TaylorModel<DATA_TYPE>::mul_insert_ctrunc_normal(TaylorModel<DATA_TYPE> & result, const TaylorModel<DATA_TYPE> & tm, const DATA_TYPE2 & tmPolyRange, const std::vector<DATA_TYPE3> & step_exp_table, const unsigned int order, const Interval & cutoff_threshold) const
{
	Polynomial<DATA_TYPE> P1xP2;
	Interval P1xI2, P2xI1, I1xI2;

	P1xP2 = expansion * tm.expansion;

	if(tm.remainder != 0)
	{
		expansion.intEvalNormal(P1xI2, step_exp_table);
		P1xI2 *= tm.remainder;
	}

	if(remainder != 0)
	{
		P2xI1 = tmPolyRange * remainder;
	}

	I1xI2 = remainder * tm.remainder;

	result.expansion = P1xP2;

	result.remainder = I1xI2;
	result.remainder += P2xI1;
	result.remainder += P1xI2;

	result.ctrunc_normal(step_exp_table, order);
	result.cutoff_normal(step_exp_table, cutoff_threshold);
}

template <class DATA_TYPE>
template <class DATA_TYPE2, class DATA_TYPE3>
void TaylorModel<DATA_TYPE>::mul_insert_ctrunc_normal(TaylorModel<DATA_TYPE> & result, Interval & tm1, Interval & intTrunc, const TaylorModel<DATA_TYPE> & tm, const DATA_TYPE2 & tmPolyRange, const std::vector<DATA_TYPE3> & step_exp_table, const unsigned int order, const Interval & cutoff_threshold) const
{
	Polynomial<DATA_TYPE> P1xP2;
	Interval P1xI2, P2xI1, I1xI2;

	P1xP2 = expansion * tm.expansion;

	tm1 = 0;
	intTrunc = 0;

	if(tm.remainder != 0)
	{
		expansion.intEvalNormal(P1xI2, step_exp_table);
		tm1 = P1xI2;
		P1xI2 *= tm.remainder;
	}

	if(remainder != 0)
	{
		P2xI1 = tmPolyRange * remainder;
	}

	I1xI2 = remainder * tm.remainder;

	result.expansion = P1xP2;

	result.remainder = I1xI2;
	result.remainder += P2xI1;
	result.remainder += P1xI2;

	result.expansion.ctrunc_normal(intTrunc, step_exp_table, order);

	Interval intRound;
	result.expansion.cutoff_normal(intRound, step_exp_table, cutoff_threshold);

	intTrunc += intRound;

	result.remainder += intTrunc;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModel<DATA_TYPE>::mul_insert_ctrunc_assign(const TaylorModel<DATA_TYPE> & tm, const DATA_TYPE2 & tmPolyRange, const std::vector<DATA_TYPE2> & domain, const unsigned int order, const Interval & cutoff_threshold)
{
	TaylorModel<DATA_TYPE> result;
	mul_insert_ctrunc(result, tm, tmPolyRange, domain, order, cutoff_threshold);
	*this = result;
}

template <class DATA_TYPE>
template <class DATA_TYPE2, class DATA_TYPE3>
void TaylorModel<DATA_TYPE>::mul_insert_ctrunc_normal_assign(const TaylorModel<DATA_TYPE> & tm, const DATA_TYPE2 & tmPolyRange, const std::vector<DATA_TYPE3> & step_exp_table, const unsigned int order, const Interval & cutoff_threshold)
{
	TaylorModel<DATA_TYPE> result;
	mul_insert_ctrunc_normal(result, tm, tmPolyRange, step_exp_table, order, cutoff_threshold);
	*this = result;
}

template <class DATA_TYPE>
template <class DATA_TYPE2, class DATA_TYPE3>
void TaylorModel<DATA_TYPE>::mul_insert_ctrunc_normal_assign(Interval & tm1, Interval & intTrunc, const TaylorModel<DATA_TYPE> & tm, const DATA_TYPE2 & tmPolyRange, const std::vector<DATA_TYPE3> & step_exp_table, const unsigned int order, const Interval & cutoff_threshold)
{
	TaylorModel<DATA_TYPE> result;
	mul_insert_ctrunc_normal(result, tm1, intTrunc, tm, tmPolyRange, step_exp_table, order, cutoff_threshold);
	*this = result;
}

template <class DATA_TYPE>
void TaylorModel<DATA_TYPE>::derivative(TaylorModel<DATA_TYPE> & result, const unsigned int varIndex) const
{
	expansion.derivative(result.expansion, varIndex);
	remainder = 0;
}

template <class DATA_TYPE>
void TaylorModel<DATA_TYPE>::LieDerivative(TaylorModel<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & f, const unsigned int order, const Interval & cutoff_threshold) const
{
	expansion.LieDerivative(result.expansion, f);
	expansion.nctrunc(order);
	expansion.cutoff(cutoff_threshold);
	remainder = 0;
}

template <class DATA_TYPE>
void TaylorModel<DATA_TYPE>::integral_time(TaylorModel<DATA_TYPE> & result, const Interval & I) const
{
	result.expansion = expansion;
	result.expansion.integral_time();
	result.remainder = remainder * I;
}

template <class DATA_TYPE>
void TaylorModel<DATA_TYPE>::integral_time(TaylorModel<DATA_TYPE> & result) const
{
	result.expansion = expansion;
	result.expansion.integral_time();
//	result.remainder = 0;
}

template <class DATA_TYPE>
void TaylorModel<DATA_TYPE>::linearCoefficients(Matrix<DATA_TYPE> & coefficients, const unsigned int row) const
{
	expansion.linearCoefficients(coefficients, row);
}

template <class DATA_TYPE>
void TaylorModel<DATA_TYPE>::linearCoefficients(std::vector<DATA_TYPE> & coefficients) const
{
	expansion.linearCoefficients(coefficients);
}

template <class DATA_TYPE>
void TaylorModel<DATA_TYPE>::toHornerForm(HornerForm<DATA_TYPE> & hf, Interval & I) const
{
	expansion.toHornerForm(hf);
	I = remainder;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModel<DATA_TYPE>::insert_ctrunc(TaylorModel<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & vars, const std::vector<DATA_TYPE2> & varsPolyRange, const std::vector<DATA_TYPE2> & domain, const unsigned int order, const Interval & cutoff_threshold) const
{
	if(vars.tms.size() == 0)
	{
		result = *this;

		typename std::list<Term<DATA_TYPE> >::iterator iter;

		for(iter = result.expansion.terms.begin(); iter != result.expansion.terms.end(); )
		{
			if( ((iter->d) - (iter->degrees[0])) > 0 )
			{
				iter = result.expansion.terms.erase(iter);
			}
			else
			{
				++iter;
			}
		}
	}
	else
	{
		HornerForm<DATA_TYPE> hf;
		expansion.toHornerForm(hf);

		hf.insert_ctrunc(result, vars, varsPolyRange, domain, order, cutoff_threshold);
		result.remainder += remainder;
	}
}

template <class DATA_TYPE>
void TaylorModel<DATA_TYPE>::insert_no_remainder(TaylorModel<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & vars, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold) const
{
	if(vars.tms.size() == 0)
	{
		result = *this;

		typename std::list<Term<DATA_TYPE> >::iterator iter;

		for(iter = result.expansion.terms.begin(); iter != result.expansion.terms.end();)
		{
			if( ((iter->d) - (iter->degrees[0])) > 0 )
			{
				iter = result.expansion.terms.erase(iter);
			}
			else
			{
				++iter;
			}
		}
	}
	else
	{
		HornerForm<DATA_TYPE> hf;
		expansion.toHornerForm(hf);

		hf.insert_no_remainder(result, vars, numVars, order, cutoff_threshold);
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2, class DATA_TYPE3>
void TaylorModel<DATA_TYPE>::insert_ctrunc_normal(TaylorModel<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & vars, const std::vector<DATA_TYPE2> & varsPolyRange, const std::vector<DATA_TYPE3> & step_exp_table, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold) const
{
	if(vars.tms.size() == 0)
	{
		result = *this;

		typename std::list<Term<DATA_TYPE> >::iterator iter;

		for(iter = result.expansion.terms.begin(); iter != result.expansion.terms.end(); )
		{
			if( ((iter->d) - (iter->degrees[0])) > 0 )
			{
				iter = result.expansion.terms.erase(iter);
			}
			else
			{
				++iter;
			}
		}
	}
	else
	{
		HornerForm<DATA_TYPE> hf;
		expansion.toHornerForm(hf);

		hf.insert_ctrunc_normal(result, vars, varsPolyRange, step_exp_table, numVars, order, cutoff_threshold);
		result.remainder += remainder;
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModel<DATA_TYPE>::evaluate_time(TaylorModel<DATA_TYPE> & result, const std::vector<DATA_TYPE2> & step_exp_table) const
{
	expansion.evaluate_time(result.expansion, step_exp_table);
	result.remainder = remainder;
}

template <class DATA_TYPE>
void TaylorModel<DATA_TYPE>::mul(TaylorModel<DATA_TYPE> & result, const unsigned int varIndex, const unsigned int degree) const
{
	expansion.mul(result.expansion, varIndex, degree);
	result.remainder = remainder;
}

template <class DATA_TYPE>
void TaylorModel<DATA_TYPE>::mul_assign(const unsigned int varIndex, const unsigned int degree)
{
	expansion.mul_assign(varIndex, degree);
}

template <class DATA_TYPE>
void TaylorModel<DATA_TYPE>::mul_assign(const unsigned int varIndex, const Interval & range)
{
	expansion.mul_assign(varIndex, 1);
	remainder *= range;
}

template <class DATA_TYPE>
void TaylorModel<DATA_TYPE>::rmConstant()
{
	expansion.rmConstant();
}

template <class DATA_TYPE>
void TaylorModel<DATA_TYPE>::decompose(TaylorModel<DATA_TYPE> & linear, TaylorModel<DATA_TYPE> & other) const
{
	expansion.decompose(linear.expansion, other.expansion);
	linear.remainder = 0;
	other.remainder = remainder;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModel<DATA_TYPE>::cutoff_normal(const std::vector<DATA_TYPE2> & step_exp_table, const Interval & cutoff_threshold)
{
	Interval I;
	expansion.cutoff_normal(I, step_exp_table, cutoff_threshold);
	remainder += I;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModel<DATA_TYPE>::cutoff(const std::vector<DATA_TYPE2> & domain, const Interval & cutoff_threshold)
{
	Interval I;
	expansion.cutoff(I, domain, cutoff_threshold);
	remainder += I;
}

template <class DATA_TYPE>
void TaylorModel<DATA_TYPE>::cutoff(const Interval & cutoff_threshold)
{
	expansion.cutoff(cutoff_threshold);
}

template <class DATA_TYPE>
unsigned int TaylorModel<DATA_TYPE>::degree() const
{
	return expansion.degree();
}

template <class DATA_TYPE>
bool TaylorModel<DATA_TYPE>::isZero() const
{
	if(expansion.terms.size() == 0)
	{
		return true;
	}
	else if(expansion.isZero())
	{
		if(remainder == 0)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	else
	{
		return false;
	}
}

template <class DATA_TYPE>
void TaylorModel<DATA_TYPE>::rmZeroTerms(const std::vector<unsigned int> & indices)
{
	expansion.rmZeroTerms(indices);
}

template <class DATA_TYPE>
void TaylorModel<DATA_TYPE>::normalize(std::vector<Interval> & domain, const Interval & cutoff_threshold)
{
	unsigned int domainDim = domain.size();
	unsigned int rangeDim = domainDim - 1;

	// compute the center of the original domain and make it origin-centered
	std::vector<Real> center(domainDim);
	for(unsigned int i=1; i<domainDim; ++i)		// we omit the time dimension
	{
		Real c;
		domain[i].remove_midpoint(c);
		center[i] = c;
	}

	// compute the scalars
	std::vector<std::vector<DATA_TYPE> > coefficients(rangeDim, std::vector<DATA_TYPE>(domainDim));

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		DATA_TYPE c;
		domain[i].mag(c);
		coefficients[i][i+1] = c;
	}

	TaylorModelVec<DATA_TYPE> newVars(coefficients);
	for(unsigned int i=0; i<rangeDim; ++i)
	{
		TaylorModel<DATA_TYPE> tmTemp(center[i], domainDim);
		newVars.tms[i] += tmTemp;
	}

	Interval intUnit(-1,1);
	for(unsigned int i=1; i<domainDim; ++i)
	{
		domain[i] = intUnit;
	}

	TaylorModel<DATA_TYPE> tmTmp;
	insert_no_remainder(tmTmp, newVars, domainDim, degree(), cutoff_threshold);
	expansion = tmTmp.expansion;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModel<DATA_TYPE>::polyRange(Interval & result, const std::vector<DATA_TYPE2> & domain) const
{
	expansion.evaluate(result, domain);
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModel<DATA_TYPE>::polyRangeNormal(Interval & result, const std::vector<DATA_TYPE2> & step_exp_table) const
{
	expansion.intEvalNormal(result, step_exp_table);
}

template <class DATA_TYPE>
Interval TaylorModel<DATA_TYPE>::getRemainder() const
{
	return remainder;
}

template <class DATA_TYPE>
void TaylorModel<DATA_TYPE>::getExpansion(Polynomial<DATA_TYPE> & p) const
{
	p = expansion;
}

template <class DATA_TYPE>
void TaylorModel<DATA_TYPE>::exp_taylor(TaylorModel<DATA_TYPE> & result, std::list<Interval> & ranges, const std::vector<Interval> & step_exp_table, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold, const Global_Computation_Setting & setting) const
{
	DATA_TYPE const_part;

	TaylorModel<DATA_TYPE> tmF = *this;

	// remove the center point of tm
	tmF.constant(const_part);
	tmF.rmConstant();			// F = tm - c

	const_part.exp_assign();	// exp(c)

	if(tmF.isZero())			// tm = c
	{
		TaylorModel<DATA_TYPE> tmExp(const_part, numVars);
		result = tmExp;

		Interval invalid(1,-1);
		ranges.push_back(invalid);

		return;
	}

	ranges.push_back(const_part);			// keep the unchanged part

	Polynomial<DATA_TYPE> polyOne(1, numVars);

	// to compute the expression 1 + F + (1/2!)F^2 + ... + (1/k!)F^k,
	// we evaluate its Horner form (...((1/(k-1))((1/k)*F+1)*F + 1) ... + 1)

	result.expansion = polyOne;
	result.remainder = 0;

	Interval tmFPolyRange;
	tmF.polyRangeNormal(tmFPolyRange, step_exp_table);

	for(int i=order; i>0; --i)
	{
		result /= i;

		Interval tm1Poly, intTrunc;
		result.mul_insert_ctrunc_normal_assign(tm1Poly, intTrunc, tmF, tmFPolyRange, step_exp_table, order, cutoff_threshold);

		ranges.push_back(tm1Poly);			// keep the unchanged part
		ranges.push_back(tmFPolyRange);		// keep the unchanged part
		ranges.push_back(intTrunc);			// keep the unchanged part

		result.expansion = polyOne;
	}

	result *= const_part;

	Interval intCutoff;
	result.expansion.cutoff_normal(intCutoff, step_exp_table, cutoff_threshold);
	ranges.push_back(intCutoff);			// keep the unchanged part
	result.remainder += intCutoff;

	Interval rem, tmRange;
	ranges.push_back(tmFPolyRange);			// keep the unchanged part
	tmRange = tmFPolyRange + tmF.remainder;
	exp_taylor_remainder(rem, tmRange, order+1, setting);

	result.remainder += const_part * rem;
}

template <class DATA_TYPE>
void TaylorModel<DATA_TYPE>::rec_taylor(TaylorModel<DATA_TYPE> & result, std::list<Interval> & ranges, const std::vector<Interval> & step_exp_table, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold, const Global_Computation_Setting & setting) const
{
	DATA_TYPE const_part;

	TaylorModel<DATA_TYPE> tmF = *this;

	// remove the center point of tm
	tmF.constant(const_part);
	tmF.rmConstant();			// F = tm - c

	const_part.rec_assign();	// 1/c

	if(tmF.isZero())			// tm = c
	{
		TaylorModel<DATA_TYPE> tmRec(const_part, numVars);
		result = tmRec;

		Interval invalid(1,-1);
		ranges.push_back(invalid);

		return;
	}

	Polynomial<DATA_TYPE> polyOne(1, numVars);
	TaylorModel<DATA_TYPE> tmF_c = tmF * const_part;

	ranges.push_back(const_part);			// keep the unchanged part


	// to compute the expression 1 - F/c + (F/c)^2 - ... + (-1)^k (F/c)^k,
	// we evaluate its Horner form (-1)*(...((-1)*(-F/c + 1)*F/c + 1)...) + 1

	result.expansion = polyOne;
	result.remainder = 0;

	Interval tmF_cPolyRange;
	tmF_c.polyRangeNormal(tmF_cPolyRange, step_exp_table);

	for(int i=order; i>0; --i)
	{
		result *= -1;

		Interval tm1Poly, intTrunc;
		result.mul_insert_ctrunc_normal_assign(tm1Poly, intTrunc, tmF_c, tmF_cPolyRange, step_exp_table, order, cutoff_threshold);

		ranges.push_back(tm1Poly);			// keep the unchanged part
		ranges.push_back(tmF_cPolyRange);	// keep the unchanged part
		ranges.push_back(intTrunc);			// keep the unchanged part

		result.expansion += polyOne;
	}

	result *= const_part;

	Interval intCutoff;
	result.expansion.cutoff_normal(intCutoff, step_exp_table, cutoff_threshold);
	ranges.push_back(intCutoff);			// keep the unchanged part
	result.remainder += intCutoff;

	Interval rem, tmF_cRange;
	ranges.push_back(tmF_cPolyRange);		// keep the unchanged part
	tmF_cRange = tmF_cPolyRange + tmF_c.remainder;

	rec_taylor_remainder(rem, tmF_cRange, order+1, setting);

	result.remainder += rem * const_part;
}

template <class DATA_TYPE>
void TaylorModel<DATA_TYPE>::sin_taylor(TaylorModel<DATA_TYPE> & result, std::list<Interval> & ranges, const std::vector<Interval> & step_exp_table, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold, const Global_Computation_Setting & setting) const
{
	DATA_TYPE const_part;

	TaylorModel<DATA_TYPE> tmF = *this;

	// remove the center point of tm
	tmF.constant(const_part);
	tmF.rmConstant();			// F = tm - c

	if(tmF.isZero())			// tm = c
	{
		const_part.sin_assign();
		TaylorModel<DATA_TYPE> tmSin(const_part, numVars);
		result = tmSin;

		Interval invalid(1,-1);
		ranges.push_back(invalid);

		return;
	}

	ranges.push_back(const_part);				// keep the unchanged part

	DATA_TYPE sinc, cosc, msinc, mcosc;
	const_part.sin(sinc);
	const_part.cos(cosc);

	msinc = -sinc;
	mcosc = -cosc;

	TaylorModel<DATA_TYPE> tmTmp(sinc, numVars);
	result = tmTmp;

	Interval tmFPolyRange;
	tmF.polyRangeNormal(tmFPolyRange, step_exp_table);

	int k = 1;

	TaylorModel<DATA_TYPE> tmPowerTmF(1, numVars);

	for(int i=1; i<=order; ++i, ++k)
	{
		k %= 4;

		switch(k)
		{
		case 0:
		{
			Interval tm1Poly, intTrunc;
			tmPowerTmF.mul_insert_ctrunc_normal_assign(tm1Poly, intTrunc, tmF, tmFPolyRange, step_exp_table, order, cutoff_threshold);

			ranges.push_back(tm1Poly);			// keep the unchanged part
			ranges.push_back(tmFPolyRange);		// keep the unchanged part
			ranges.push_back(intTrunc);			// keep the unchanged part

			Real tmp = sinc / i;
			ranges.push_back(tmp);				// keep the unchanged part

			tmPowerTmF *= tmp;
			result += tmPowerTmF;

			break;
		}
		case 1:
		{
			Interval tm1Poly, intTrunc;
			tmPowerTmF.mul_insert_ctrunc_normal_assign(tm1Poly, intTrunc, tmF, tmFPolyRange, step_exp_table, order, cutoff_threshold);

			ranges.push_back(tm1Poly);			// keep the unchanged part
			ranges.push_back(tmFPolyRange);		// keep the unchanged part
			ranges.push_back(intTrunc);			// keep the unchanged part

			Real tmp = cosc / i;
			ranges.push_back(tmp);				// keep the unchanged part

			tmPowerTmF *= tmp;
			result += tmPowerTmF;

			break;
		}
		case 2:
		{
			Interval tm1Poly, intTrunc;
			tmPowerTmF.mul_insert_ctrunc_normal_assign(tm1Poly, intTrunc, tmF, tmFPolyRange, step_exp_table, order, cutoff_threshold);

			ranges.push_back(tm1Poly);			// keep the unchanged part
			ranges.push_back(tmFPolyRange);		// keep the unchanged part
			ranges.push_back(intTrunc);			// keep the unchanged part

			Real tmp = msinc / i;
			ranges.push_back(tmp);				// keep the unchanged part

			tmPowerTmF *= tmp;
			result += tmPowerTmF;

			break;
		}
		case 3:
		{
			Interval tm1Poly, intTrunc;
			tmPowerTmF.mul_insert_ctrunc_normal_assign(tm1Poly, intTrunc, tmF, tmFPolyRange, step_exp_table, order, cutoff_threshold);

			ranges.push_back(tm1Poly);			// keep the unchanged part
			ranges.push_back(tmFPolyRange);		// keep the unchanged part
			ranges.push_back(intTrunc);			// keep the unchanged part

			Real tmp = mcosc / i;
			ranges.push_back(tmp);				// keep the unchanged part

			tmPowerTmF *= tmp;
			result += tmPowerTmF;

			break;
		}
		}
	}

	Interval intCutoff;
	result.expansion.cutoff_normal(intCutoff, step_exp_table, cutoff_threshold);
	ranges.push_back(intCutoff);				// keep the unchanged part
	result.remainder += intCutoff;

	// evaluate the remainder
	Interval tmRange, rem;
	tmRange = tmFPolyRange + tmF.remainder;

	ranges.push_back(tmFPolyRange);				// keep the unchanged part
	sin_taylor_remainder(rem, const_part, tmRange, order+1, setting);

	result.remainder += rem;
}

template <class DATA_TYPE>
void TaylorModel<DATA_TYPE>::cos_taylor(TaylorModel<DATA_TYPE> & result, std::list<Interval> & ranges, const std::vector<Interval> & step_exp_table, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold, const Global_Computation_Setting & setting) const
{
	DATA_TYPE const_part;

	TaylorModel<DATA_TYPE> tmF = *this;

	// remove the center point of tm
	tmF.constant(const_part);
	tmF.rmConstant();			// F = tm - c

	if(tmF.isZero())			// tm = c
	{
		const_part.cos_assign();
		TaylorModel<DATA_TYPE> tmCos(const_part, numVars);
		result = tmCos;

		Interval invalid(1,-1);
		ranges.push_back(invalid);

		return;
	}

	ranges.push_back(const_part);				// keep the unchanged part

	DATA_TYPE sinc, cosc, msinc, mcosc;
	const_part.sin(sinc);
	const_part.cos(cosc);

	msinc = -sinc;
	mcosc = -cosc;


	TaylorModel<DATA_TYPE> tmTmp(cosc, numVars);
	result = tmTmp;

	Interval tmFPolyRange;
	tmF.polyRangeNormal(tmFPolyRange, step_exp_table);

	int k = 1;

	TaylorModel tmPowerTmF(1, numVars);

	for(int i=1; i<=order; ++i, ++k)
	{
		k %= 4;

		switch(k)
		{
		case 0:
		{
			Interval tm1Poly, intTrunc;
			tmPowerTmF.mul_insert_ctrunc_normal_assign(tm1Poly, intTrunc, tmF, tmFPolyRange, step_exp_table, order, cutoff_threshold);

			ranges.push_back(tm1Poly);			// keep the unchanged part
			ranges.push_back(tmFPolyRange);		// keep the unchanged part
			ranges.push_back(intTrunc);			// keep the unchanged part

			Real tmp = cosc / i;
			ranges.push_back(tmp);				// keep the unchanged part

			tmPowerTmF *= tmp;
			result += tmPowerTmF;

			break;
		}
		case 1:
		{
			Interval tm1Poly, intTrunc;
			tmPowerTmF.mul_insert_ctrunc_normal_assign(tm1Poly, intTrunc, tmF, tmFPolyRange, step_exp_table, order, cutoff_threshold);

			ranges.push_back(tm1Poly);			// keep the unchanged part
			ranges.push_back(tmFPolyRange);		// keep the unchanged part
			ranges.push_back(intTrunc);			// keep the unchanged part

			Real tmp = msinc / i;
			ranges.push_back(tmp);				// keep the unchanged part

			tmPowerTmF *= tmp;
			result += tmPowerTmF;

			break;
		}
		case 2:
		{
			Interval tm1Poly, intTrunc;
			tmPowerTmF.mul_insert_ctrunc_normal_assign(tm1Poly, intTrunc, tmF, tmFPolyRange, step_exp_table, order, cutoff_threshold);

			ranges.push_back(tm1Poly);			// keep the unchanged part
			ranges.push_back(tmFPolyRange);		// keep the unchanged part
			ranges.push_back(intTrunc);			// keep the unchanged part

			Real tmp = mcosc / i;
			ranges.push_back(tmp);				// keep the unchanged part

			tmPowerTmF *= tmp;
			result += tmPowerTmF;

			break;
		}
		case 3:
		{
			Interval tm1Poly, intTrunc;
			tmPowerTmF.mul_insert_ctrunc_normal_assign(tm1Poly, intTrunc, tmF, tmFPolyRange, step_exp_table, order, cutoff_threshold);

			ranges.push_back(tm1Poly);			// keep the unchanged part
			ranges.push_back(tmFPolyRange);		// keep the unchanged part
			ranges.push_back(intTrunc);			// keep the unchanged part

			Real tmp = sinc / i;
			ranges.push_back(tmp);				// keep the unchanged part

			tmPowerTmF *= tmp;
			result += tmPowerTmF;

			break;
		}
		}
	}

	Interval intCutoff;
	result.expansion.cutoff_normal(intCutoff, step_exp_table, cutoff_threshold);
	ranges.push_back(intCutoff);				// keep the unchanged part
	result.remainder += intCutoff;

	// evaluate the remainder
	Interval tmRange, rem;
	tmRange = tmFPolyRange + tmF.remainder;

	ranges.push_back(tmFPolyRange);				// keep the unchanged part
	cos_taylor_remainder(rem, const_part, tmRange, order+1, setting);

	result.remainder += rem;
}

template <class DATA_TYPE>
void TaylorModel<DATA_TYPE>::log_taylor(TaylorModel<DATA_TYPE> & result, std::list<Interval> & ranges, const std::vector<Interval> & step_exp_table, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold) const
{
	DATA_TYPE const_part;

	TaylorModel<DATA_TYPE> tmF = *this;

	// remove the center point of tm
	tmF.constant(const_part);
	tmF.rmConstant();			// F = tm - c

	DATA_TYPE C = const_part;
	ranges.push_back(const_part);			// keep the unchanged part

	const_part.log_assign();	// log(c)

	if(tmF.isZero())			// tm = c
	{
		TaylorModel<DATA_TYPE> tmLog(const_part, numVars);
		result = tmLog;

		Interval invalid(1,-1);
		ranges.push_back(invalid);

		return;
	}

	TaylorModel<DATA_TYPE> tmF_c = tmF / C;
	result = tmF_c / order;


	Interval tmF_cPolyRange;
	tmF_c.polyRangeNormal(tmF_cPolyRange, step_exp_table);

	for(int i=order-1; i>=1; --i)
	{
		Real tmp = 1;
		tmp /= i;

		TaylorModel<DATA_TYPE> tmTmp(tmp, numVars);
		result -= tmTmp;
		result *= -1;

		Interval tm1Poly, intTrunc;
		result.mul_insert_ctrunc_normal_assign(tm1Poly, intTrunc, tmF_c, tmF_cPolyRange, step_exp_table, order, cutoff_threshold);

		ranges.push_back(tm1Poly);			// keep the unchanged part
		ranges.push_back(tmF_cPolyRange);	// keep the unchanged part
		ranges.push_back(intTrunc);			// keep the unchanged part
	}

	TaylorModel const_part_tm(const_part, numVars);
	result += const_part_tm;

	Interval intCutoff;
	result.expansion.cutoff_normal(intCutoff, step_exp_table, cutoff_threshold);
	ranges.push_back(intCutoff);			// keep the unchanged part
	result.remainder += intCutoff;

	Interval rem, tmF_cRange;
	ranges.push_back(tmF_cPolyRange);		// keep the unchanged part
	tmF_cRange = tmF_cPolyRange + tmF_c.remainder;

	log_taylor_remainder(rem, tmF_cRange, order+1);

	result.remainder += rem;
}

template <class DATA_TYPE>
void TaylorModel<DATA_TYPE>::sqrt_taylor(TaylorModel<DATA_TYPE> & result, std::list<Interval> & ranges, const std::vector<Interval> & step_exp_table, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold, const Global_Computation_Setting & setting) const
{
	DATA_TYPE const_part;

	TaylorModel<DATA_TYPE> tmF = *this;

	// remove the center point of tm
	tmF.constant(const_part);
	tmF.rmConstant();			// F = tm - c

	DATA_TYPE C = const_part;
	ranges.push_back(const_part);			// keep the unchanged part

	const_part.sqrt_assign();	// sqrt(c)

	if(tmF.isZero())			// tm = c
	{
		TaylorModel<DATA_TYPE> tmSqrt(const_part, numVars);
		result = tmSqrt;

		Interval invalid(1,-1);
		ranges.push_back(invalid);

		return;
	}

	TaylorModel<DATA_TYPE> tmF_2c = tmF / (2*C);
	Polynomial<DATA_TYPE> polyOne(1, numVars);

	result = tmF_2c;

	Interval tmF_2cPolyRange;
	tmF_2c.polyRangeNormal(tmF_2cPolyRange, step_exp_table);

	for(int i=order, j=2*order-3; i>=2; --i, j-=2)
	{
		DATA_TYPE tmp = j;
		tmp /= -i;
		result *= tmp;
		result.expansion += polyOne;

		Interval tm1Poly, intTrunc;
		result.mul_insert_ctrunc_normal_assign(tm1Poly, intTrunc, tmF_2c, tmF_2cPolyRange, step_exp_table, order, cutoff_threshold);

		ranges.push_back(tm1Poly);			// keep the unchanged part
		ranges.push_back(tmF_2cPolyRange);	// keep the unchanged part
		ranges.push_back(intTrunc);			// keep the unchanged part
	}

	result.expansion += polyOne;

	result *= const_part;

	Interval intCutoff;
	result.expansion.cutoff_normal(intCutoff, step_exp_table, cutoff_threshold);
	ranges.push_back(intCutoff);			// keep the unchanged part
	result.remainder += intCutoff;

	Interval rem, tmF_cRange = tmF_2cPolyRange;
	tmF_cRange *= 2;

	ranges.push_back(tmF_cRange);			// keep the unchanged part

	Interval tmF_c_remainder = tmF_2c.remainder;
	tmF_c_remainder *= 2;
	tmF_cRange += tmF_c_remainder;

	sqrt_taylor_remainder(rem, tmF_cRange, order+1, setting);

	result.remainder += rem * const_part;
}

template <class DATA_TYPE>
void TaylorModel<DATA_TYPE>::exp_taylor(TaylorModel<DATA_TYPE> & result, const std::vector<Interval> & domain, const unsigned int order, const Interval & cutoff_threshold, const Global_Computation_Setting & setting) const
{
	DATA_TYPE const_part;

	TaylorModel<DATA_TYPE> tmF = *this;

	// remove the center point of tm
	tmF.constant(const_part);
	tmF.rmConstant();			// F = tm - c

	const_part.exp_assign();	// exp(c)

	unsigned int numVars = domain.size();

	if(tmF.isZero())			// tm = c
	{
		TaylorModel<DATA_TYPE> tmExp(const_part, numVars);
		result = tmExp;

		return;
	}

	Polynomial<DATA_TYPE> polyOne(1, numVars);

	// to compute the expression 1 + F + (1/2!)F^2 + ... + (1/k!)F^k,
	// we evaluate its Horner form (...((1/(k-1))((1/k)*F+1)*F + 1) ... + 1)

	result.expansion = polyOne;
	result.remainder = 0;

	Interval tmFPolyRange;
	tmF.polyRange(tmFPolyRange, domain);

	for(int i=order; i>0; --i)
	{
		result /= i;

		result.mul_insert_ctrunc_assign(tmF, tmFPolyRange, domain, order, cutoff_threshold);

		result.expansion = polyOne;
	}

	result *= const_part;

	result.cutoff(domain, cutoff_threshold);

	Interval rem, tmRange;
	tmRange = tmFPolyRange + tmF.remainder;
	exp_taylor_remainder(rem, tmRange, order+1, setting);

	result.remainder += const_part * rem;
}

template <class DATA_TYPE>
void TaylorModel<DATA_TYPE>::rec_taylor(TaylorModel<DATA_TYPE> & result, const std::vector<Interval> & domain, const unsigned int order, const Interval & cutoff_threshold, const Global_Computation_Setting & setting) const
{
	DATA_TYPE const_part;

	TaylorModel<DATA_TYPE> tmF = *this;

	// remove the center point of tm
	tmF.constant(const_part);
	tmF.rmConstant();			// F = tm - c

	unsigned int numVars = domain.size();

	const_part.rec_assign();	// 1/c

	if(tmF.isZero())			// tm = c
	{
		TaylorModel<DATA_TYPE> tmRec(const_part, numVars);
		result = tmRec;

		return;
	}

	Polynomial<DATA_TYPE> polyOne(1, numVars);
	TaylorModel<DATA_TYPE> tmF_c = tmF * const_part;


	// to compute the expression 1 - F/c + (F/c)^2 - ... + (-1)^k (F/c)^k,
	// we evaluate its Horner form (-1)*(...((-1)*(-F/c + 1)*F/c + 1)...) + 1

	result.expansion = polyOne;
	result.remainder = 0;

	Interval tmF_cPolyRange;
	tmF_c.polyRange(tmF_cPolyRange, domain);

	for(int i=order; i>0; --i)
	{
		result *= -1;

		result.mul_insert_ctrunc_assign(tmF_c, tmF_cPolyRange, domain, order, cutoff_threshold);

		result.expansion += polyOne;
	}

	result *= const_part;

	result.cutoff(domain, cutoff_threshold);

	Interval rem, tmF_cRange;
	tmF_cRange = tmF_cPolyRange + tmF_c.remainder;

	rec_taylor_remainder(rem, tmF_cRange, order+1, setting);

	result.remainder += rem * const_part;
}

template <class DATA_TYPE>
void TaylorModel<DATA_TYPE>::sin_taylor(TaylorModel<DATA_TYPE> & result, const std::vector<Interval> & domain, const unsigned int order, const Interval & cutoff_threshold, const Global_Computation_Setting & setting) const
{
	DATA_TYPE const_part;

	TaylorModel<DATA_TYPE> tmF = *this;

	// remove the center point of tm
	tmF.constant(const_part);
	tmF.rmConstant();			// F = tm - c

	unsigned int numVars = domain.size();

	if(tmF.isZero())			// tm = c
	{
		const_part.sin_assign();
		TaylorModel<DATA_TYPE> tmSin(const_part, numVars);
		result = tmSin;

		return;
	}

	DATA_TYPE sinc, cosc, msinc, mcosc;
	const_part.sin(sinc);
	const_part.cos(cosc);

	msinc = -sinc;
	mcosc = -cosc;

	TaylorModel<DATA_TYPE> tmTmp(sinc, numVars);
	result = tmTmp;

	Interval tmFPolyRange;
	tmF.polyRange(tmFPolyRange, domain);

	int k = 1;

	TaylorModel<DATA_TYPE> tmPowerTmF(1, numVars);

	for(int i=1; i<=order; ++i, ++k)
	{
		k %= 4;

		switch(k)
		{
		case 0:
		{
			tmPowerTmF.mul_insert_ctrunc_assign(tmF, tmFPolyRange, domain, order, cutoff_threshold);

			Real tmp = sinc / i;

			tmPowerTmF *= tmp;
			result += tmPowerTmF;

			break;
		}
		case 1:
		{
			tmPowerTmF.mul_insert_ctrunc_assign(tmF, tmFPolyRange, domain, order, cutoff_threshold);

			Real tmp = cosc / i;

			tmPowerTmF *= tmp;
			result += tmPowerTmF;

			break;
		}
		case 2:
		{
			tmPowerTmF.mul_insert_ctrunc_assign(tmF, tmFPolyRange, domain, order, cutoff_threshold);

			Real tmp = msinc / i;

			tmPowerTmF *= tmp;
			result += tmPowerTmF;

			break;
		}
		case 3:
		{
			tmPowerTmF.mul_insert_ctrunc_assign(tmF, tmFPolyRange, domain, order, cutoff_threshold);

			Real tmp = mcosc / i;

			tmPowerTmF *= tmp;
			result += tmPowerTmF;

			break;
		}
		}
	}

	result.cutoff(domain, cutoff_threshold);

	// evaluate the remainder
	Interval tmRange, rem;
	tmRange = tmFPolyRange + tmF.remainder;

	sin_taylor_remainder(rem, const_part, tmRange, order+1, setting);

	result.remainder += rem;
}

template <class DATA_TYPE>
void TaylorModel<DATA_TYPE>::cos_taylor(TaylorModel<DATA_TYPE> & result, const std::vector<Interval> & domain, const unsigned int order, const Interval & cutoff_threshold, const Global_Computation_Setting & setting) const
{
	DATA_TYPE const_part;

	TaylorModel<DATA_TYPE> tmF = *this;

	// remove the center point of tm
	tmF.constant(const_part);
	tmF.rmConstant();			// F = tm - c

	unsigned int numVars = domain.size();

	if(tmF.isZero())			// tm = c
	{
		const_part.cos_assign();
		TaylorModel<DATA_TYPE> tmCos(const_part, numVars);
		result = tmCos;

		return;
	}

	DATA_TYPE sinc, cosc, msinc, mcosc;
	const_part.sin(sinc);
	const_part.cos(cosc);

	msinc = -sinc;
	mcosc = -cosc;


	TaylorModel<DATA_TYPE> tmTmp(cosc, numVars);
	result = tmTmp;

	Interval tmFPolyRange;
	tmF.polyRange(tmFPolyRange, domain);

	int k = 1;

	TaylorModel tmPowerTmF(1, numVars);

	for(int i=1; i<=order; ++i, ++k)
	{
		k %= 4;

		switch(k)
		{
		case 0:
		{
			tmPowerTmF.mul_insert_ctrunc_assign(tmF, tmFPolyRange, domain, order, cutoff_threshold);

			Real tmp = cosc / i;

			tmPowerTmF *= tmp;
			result += tmPowerTmF;

			break;
		}
		case 1:
		{
			tmPowerTmF.mul_insert_ctrunc_normal_assign(tmF, tmFPolyRange, domain, order, cutoff_threshold);

			Real tmp = msinc / i;

			tmPowerTmF *= tmp;
			result += tmPowerTmF;

			break;
		}
		case 2:
		{
			tmPowerTmF.mul_insert_ctrunc_normal_assign(tmF, tmFPolyRange, domain, order, cutoff_threshold);

			Real tmp = mcosc / i;

			tmPowerTmF *= tmp;
			result += tmPowerTmF;

			break;
		}
		case 3:
		{
			tmPowerTmF.mul_insert_ctrunc_normal_assign(tmF, tmFPolyRange, domain, order, cutoff_threshold);

			Real tmp = sinc / i;

			tmPowerTmF *= tmp;
			result += tmPowerTmF;

			break;
		}
		}
	}

	result.cutoff(domain, cutoff_threshold);

	// evaluate the remainder
	Interval tmRange, rem;
	tmRange = tmFPolyRange + tmF.remainder;

	cos_taylor_remainder(rem, const_part, tmRange, order+1, setting);

	result.remainder += rem;
}

template <class DATA_TYPE>
void TaylorModel<DATA_TYPE>::log_taylor(TaylorModel<DATA_TYPE> & result, const std::vector<Interval> & domain, const unsigned int order, const Interval & cutoff_threshold) const
{
	DATA_TYPE const_part;

	TaylorModel<DATA_TYPE> tmF = *this;

	// remove the center point of tm
	tmF.constant(const_part);
	tmF.rmConstant();			// F = tm - c

	unsigned int numVars = domain.size();

	DATA_TYPE C = const_part;

	const_part.log_assign();	// log(c)

	if(tmF.isZero())			// tm = c
	{
		TaylorModel<DATA_TYPE> tmLog(const_part, numVars);
		result = tmLog;

		return;
	}

	TaylorModel<DATA_TYPE> tmF_c = tmF / C;
	result = tmF_c / order;

	Interval tmF_cPolyRange;
	tmF_c.polyRange(tmF_cPolyRange, domain);

	for(int i=order-1; i>=1; --i)
	{
		Real tmp = 1;
		tmp /= i;

		TaylorModel<DATA_TYPE> tmTmp(tmp, numVars);
		result -= tmTmp;
		result *= -1;

		result.mul_insert_ctrunc_assign(tmF_c, tmF_cPolyRange, domain, order, cutoff_threshold);
	}

	TaylorModel const_part_tm(const_part, numVars);
	result += const_part_tm;

	result.cutoff(domain, cutoff_threshold);

	Interval rem, tmF_cRange;
	tmF_cRange = tmF_cPolyRange + tmF_c.remainder;

	log_taylor_remainder(rem, tmF_cRange, order+1);

	result.remainder += rem;
}

template <class DATA_TYPE>
void TaylorModel<DATA_TYPE>::sqrt_taylor(TaylorModel<DATA_TYPE> & result, const std::vector<Interval> & domain, const unsigned int order, const Interval & cutoff_threshold, const Global_Computation_Setting & setting) const
{
	DATA_TYPE const_part;

	TaylorModel<DATA_TYPE> tmF = *this;

	// remove the center point of tm
	tmF.constant(const_part);
	tmF.rmConstant();			// F = tm - c

	unsigned int numVars = domain.size();

	DATA_TYPE C = const_part;

	const_part.sqrt_assign();	// sqrt(c)

	if(tmF.isZero())			// tm = c
	{
		TaylorModel<DATA_TYPE> tmSqrt(const_part, numVars);
		result = tmSqrt;

		return;
	}

	TaylorModel<DATA_TYPE> tmF_2c = tmF / (2*C);
	Polynomial<DATA_TYPE> polyOne(1, numVars);

	result = tmF_2c;

	Interval tmF_2cPolyRange;
	tmF_2c.polyRange(tmF_2cPolyRange, domain);

	for(int i=order, j=2*order-3; i>=2; --i, j-=2)
	{
		DATA_TYPE tmp = j;
		tmp /= -i;
		result *= tmp;
		result.expansion += polyOne;

		result.mul_insert_ctrunc_assign(tmF_2c, tmF_2cPolyRange, domain, order, cutoff_threshold);
	}

	result.expansion += polyOne;

	result *= const_part;

	result.cutoff(domain, cutoff_threshold);

	Interval rem, tmF_cRange = tmF_2cPolyRange;
	tmF_cRange *= 2;

	Interval tmF_c_remainder = tmF_2c.remainder;
	tmF_c_remainder *= 2;
	tmF_cRange += tmF_c_remainder;

	sqrt_taylor_remainder(rem, tmF_cRange, order+1, setting);

	result.remainder += rem * const_part;
}















template <class DATA_TYPE>
class TaylorModelVec			// Taylor models: R^n -> R^m
{
public:
	std::vector<TaylorModel<DATA_TYPE> > tms;

public:
	TaylorModelVec();
	TaylorModelVec(const std::vector<TaylorModel<DATA_TYPE> > & tms_input);
	TaylorModelVec(const std::vector<DATA_TYPE> & constants, const unsigned int numVars);

	TaylorModelVec(Matrix<DATA_TYPE> & coefficients);
	TaylorModelVec(const std::vector<DATA_TYPE> & coefficients);
	TaylorModelVec(Matrix<DATA_TYPE> & coefficients, const std::vector<Interval> & remainders);

	TaylorModelVec(const std::vector<std::vector<DATA_TYPE> > & coefficients);
	TaylorModelVec(const std::vector<std::vector<DATA_TYPE> > & coefficients, const std::vector<Interval> & remainders);

	TaylorModelVec(const std::vector<Interval> & box, std::vector<Interval> & domain);
	TaylorModelVec(Matrix<Interval> & box, std::vector<Interval> & domain);

	TaylorModelVec(const unsigned int dim);
	TaylorModelVec(const TaylorModelVec<DATA_TYPE> & tmv);
	~TaylorModelVec();

	void clear();

//	void dump_interval(FILE *fp, const std::vector<std::string> & stateVarNames, const std::vector<std::string> & tmVarNames) const;
//	void dump_constant(FILE *fp, const std::vector<std::string> & stateVarNames, const std::vector<std::string> & tmVarNames) const;

	void output(std::ostream & os, const Variables & stateVars, const Variables & tmVars) const;
	void constant(std::vector<DATA_TYPE> & c) const;

	template <class DATA_TYPE2>
	void intEval(std::vector<Interval> & result, const std::vector<DATA_TYPE2> & domain) const;

	template <class DATA_TYPE2>
	void intEval(std::vector<Interval> & result, const std::vector<DATA_TYPE2> & domain, const std::vector<unsigned int> & varIDs) const;

	template <class DATA_TYPE2>
	void intEvalNormal(std::vector<Interval> & result, const std::vector<DATA_TYPE2> & step_exp_table) const;

	template <class DATA_TYPE2>
	void intEvalNormal(std::vector<Interval> & result, const std::vector<DATA_TYPE2> & step_exp_table, const std::vector<unsigned int> & varIDs) const;

	template <class DATA_TYPE2>
	void ctrunc(const std::vector<DATA_TYPE2> & domain, const unsigned int order);

	void nctrunc(const unsigned int order);

	template <class DATA_TYPE2>
	void ctrunc_normal(const std::vector<DATA_TYPE2> & step_exp_table, const unsigned int order);

	template <class DATA_TYPE2>
	void ctrunc(const std::vector<DATA_TYPE2> & domain, const std::vector<unsigned int> & orders);

	void nctrunc(const std::vector<unsigned int> & orders);

	template <class DATA_TYPE2>
	void ctrunc_normal(const std::vector<DATA_TYPE2> & step_exp_table, const std::vector<unsigned int> & orders);


	TaylorModelVec<DATA_TYPE> & operator = (const TaylorModelVec<DATA_TYPE> & tmv);
	TaylorModelVec<DATA_TYPE> & operator = (const std::vector<Polynomial<DATA_TYPE> > & pv);

	TaylorModelVec<DATA_TYPE> & operator += (const TaylorModelVec<DATA_TYPE> & tmv);
	TaylorModelVec<DATA_TYPE> & operator -= (const TaylorModelVec<DATA_TYPE> & tmv);

	TaylorModelVec<DATA_TYPE> operator + (const TaylorModelVec<DATA_TYPE> & tmv) const;
	TaylorModelVec<DATA_TYPE> operator - (const TaylorModelVec<DATA_TYPE> & tmv) const;

	TaylorModelVec<DATA_TYPE> & operator *= (const DATA_TYPE & c);
	TaylorModelVec<DATA_TYPE> & operator /= (const DATA_TYPE & c);
	TaylorModelVec<DATA_TYPE> operator * (const DATA_TYPE & c) const;
	TaylorModelVec<DATA_TYPE> operator / (const DATA_TYPE & c) const;


	void derivative(TaylorModelVec<DATA_TYPE> & result, const unsigned int varIndex) const;

	void LieDerivative(TaylorModelVec<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & f, const unsigned int order, const Interval & cutoff_threshold) const;
	void LieDerivative(TaylorModelVec<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & f, const std::vector<unsigned int> & orders, const Interval & cutoff_threshold) const;

	void integral_time(TaylorModelVec<DATA_TYPE> & result, const Interval & I) const;
	void integral_time(TaylorModelVec<DATA_TYPE> & result) const;

	void linearCoefficients(Matrix<DATA_TYPE> & coefficients) const;
	void linearCoefficients(std::vector<std::vector<DATA_TYPE> > & coefficients) const;

	void rmZeroTerms(const std::vector<unsigned int> & indices);



//	void insert(TaylorModelVec<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & vars, const std::vector<DATA_TYPE2> & varsPolyRange, const std::vector<DATA_TYPE2> & domain, const Interval & cutoff_threshold) const;
//	void insert_normal(TaylorModelVec & result, const TaylorModelVec & vars, const std::vector<Interval> & varsPolyRange, const std::vector<Interval> & step_exp_table, const int numVars, const Interval & cutoff_threshold) const;

	template <class DATA_TYPE2>
	void insert_ctrunc(TaylorModelVec<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & vars, const std::vector<DATA_TYPE2> & varsPolyRange, const std::vector<DATA_TYPE2> & domain, const unsigned int order, const Interval & cutoff_threshold) const;

	void insert_no_remainder(TaylorModelVec<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & vars, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold) const;

	template <class DATA_TYPE2, class DATA_TYPE3>
	void insert_ctrunc_normal(TaylorModelVec<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & vars, const std::vector<DATA_TYPE2> & varsPolyRange, const std::vector<DATA_TYPE3> & step_exp_table, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold) const;

//	void insert_ctrunc_normal_no_cutoff(TaylorModelVec & result, const TaylorModelVec & vars, const std::vector<Interval> & varsPolyRange, const std::vector<Interval> & step_exp_table, const int numVars, const int order) const;

	template <class DATA_TYPE2>
	void insert_ctrunc(TaylorModelVec<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & vars, const std::vector<DATA_TYPE2> & varsPolyRange, const std::vector<DATA_TYPE2> & domain, const std::vector<unsigned int> & orders, const Interval & cutoff_threshold) const;

	void insert_no_remainder(TaylorModelVec<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & vars, const unsigned int numVars, const std::vector<unsigned int> & orders, const Interval & cutoff_threshold) const;

	template <class DATA_TYPE2, class DATA_TYPE3>
	void insert_ctrunc_normal(TaylorModelVec<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & vars, const std::vector<DATA_TYPE2> & varsPolyRange, const std::vector<DATA_TYPE3> & step_exp_table, const unsigned int numVars, const std::vector<unsigned int> & orders, const Interval & cutoff_threshold) const;

//	void insert_no_remainder_no_cutoff(TaylorModelVec & result, const TaylorModelVec & vars, const int numVars, const int order) const;

	template <class DATA_TYPE2>
	void evaluate_time(TaylorModelVec<DATA_TYPE> & result, const std::vector<DATA_TYPE2> & step_exp_table) const;

//	void mul(TaylorModelVec & result, const int varIndex, const int degree) const;
//	void mul_assign(const int varIndex, const int degree);

//	void linearTrans(TaylorModelVec & result, Matrix<double> & A) const;		// linear transformation
//	void linearTrans(TaylorModelVec & result, rMatrix & A) const;
//	void linearTrans_assign(Matrix<double> & A);
//	void linearTrans_assign(rMatrix & A);

	void scale(TaylorModelVec<DATA_TYPE> & result, const std::vector<DATA_TYPE> & S);
	void scale_assign(const std::vector<DATA_TYPE> & S);

	void rmConstant();
	void decompose(TaylorModelVec<DATA_TYPE> & linear, TaylorModelVec<DATA_TYPE> & other) const;

	template <class DATA_TYPE2>
	void cutoff_normal(const std::vector<DATA_TYPE2> & step_exp_table, const Interval & cutoff_threshold);

	template <class DATA_TYPE2>
	void cutoff(const std::vector<DATA_TYPE2> & domain, const Interval & cutoff_threshold);

	void cutoff(const Interval & cutoff_threshold);

	double rho(const std::vector<Real> & l, const std::vector<Interval> & domain) const;
	double rho_normal(const std::vector<Real> & l, const std::vector<Interval> & step_exp_table) const;

//	template <class DATA_TYPE2>
//	void cutoff_normal(Matrix<Interval> & M, const std::vector<DATA_TYPE2> & step_exp_table, const Interval & cutoff_threshold);

//	void center_nc();
	void Expansion(std::vector<Polynomial<DATA_TYPE> > & polys) const;
	void Remainder(Matrix<Interval> & rem) const;

//	void extend(const int num);
//	void extend();


	template <class DATA_TYPE2>
	void Picard_no_remainder(TaylorModelVec<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & x0, const std::vector<HornerForm<DATA_TYPE2> > & ode, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold) const;

	template <class DATA_TYPE2>
	void Picard_no_remainder_assign(const TaylorModelVec<DATA_TYPE> & x0, const std::vector<HornerForm<DATA_TYPE2> > & ode, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold);

	template <class DATA_TYPE2>
	void Picard_no_remainder(TaylorModelVec<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & x0, const std::vector<HornerForm<DATA_TYPE2> > & ode, const unsigned int numVars, const std::vector<unsigned int> & orders, const Interval & cutoff_threshold) const;

	template <class DATA_TYPE2>
	void Picard_no_remainder_assign(const TaylorModelVec<DATA_TYPE> & x0, const std::vector<HornerForm<DATA_TYPE2> > & ode, const unsigned int numVars, const std::vector<unsigned int> & orders, const Interval & cutoff_threshold);


	template <class DATA_TYPE2, class DATA_TYPE3>
	void Picard_ctrunc_normal(TaylorModelVec<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & x0, const std::vector<HornerForm<DATA_TYPE2> > & ode, const std::vector<DATA_TYPE3> & varsPolyRange, const std::vector<DATA_TYPE3> & step_exp_table, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold, std::list<Interval> & intermediate_ranges) const;

	template <class DATA_TYPE2>
	void Picard_ctrunc_normal_remainder(std::vector<Interval> & result, const std::vector<HornerForm<DATA_TYPE2> > & ode, const Interval & timeStep, std::list<Interval> & intermediate_ranges) const;

	template <class DATA_TYPE2, class DATA_TYPE3>
	void Picard_ctrunc_normal(TaylorModelVec<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & x0, const std::vector<HornerForm<DATA_TYPE2> > & ode, const std::vector<DATA_TYPE3> & varsPolyRange, const std::vector<DATA_TYPE3> & step_exp_table, const unsigned int numVars, const std::vector<unsigned int> & orders, const Interval & cutoff_threshold, std::list<Interval> & intermediate_ranges) const;


	template <class DATA_TYPE2>
	void Picard_no_remainder(TaylorModelVec<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & x0, const std::vector<Expression_AST<DATA_TYPE2> > & ode, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold) const;

	template <class DATA_TYPE2>
	void Picard_no_remainder_assign(const TaylorModelVec<DATA_TYPE> & x0, const std::vector<Expression_AST<DATA_TYPE2> > & ode, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold);

	template <class DATA_TYPE2, class DATA_TYPE3>
	void Picard_ctrunc_normal(TaylorModelVec<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & x0, const std::vector<Expression_AST<DATA_TYPE2> > & ode, const std::vector<DATA_TYPE3> & step_exp_table, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold, std::list<Interval> & intermediate_ranges, const Global_Computation_Setting & setting) const;

	template <class DATA_TYPE2>
	void Picard_ctrunc_normal_remainder(std::vector<Interval> & result, const std::vector<Expression_AST<DATA_TYPE2> > & ode, const Interval & timeStep, const unsigned int order, std::list<Interval> & intermediate_ranges, const Global_Computation_Setting & setting) const;





//	void Picard_no_remainder(TaylorModelVec & result, const TaylorModelVec & x0, const std::vector<HornerForm> & ode, const std::vector<std::vector<bool> > & substitution, const int numVars, const int order, const Interval & cutoff_threshold) const;
//	void Picard_no_remainder_no_cutoff(TaylorModelVec & result, const TaylorModelVec & x0, const std::vector<HornerForm> & ode, const int numVars, const int order) const;
//	void Picard_no_remainder_no_cutoff(TaylorModelVec & result, const TaylorModelVec & x0, const std::vector<HornerForm> & ode, const std::vector<TaylorModelVec> & sub_tmvs, const int numVars, const int order) const;


//	void Picard_no_remainder_assign(const TaylorModelVec & x0, const std::vector<HornerForm> & ode, const std::vector<std::vector<bool> > & substitution, const int numVars, const int order, const Interval & cutoff_threshold);

//	void Picard_no_remainder_no_cutoff_assign(const TaylorModelVec & x0, const std::vector<HornerForm> & ode, const int numVars, const int order);
//	void Picard_no_remainder_no_cutoff_assign(const TaylorModelVec & x0, const std::vector<HornerForm> & ode, const std::vector<TaylorModelVec> & sub_tmvs, const int numVars, const int order);

//	void Picard_ctrunc_normal(TaylorModelVec & result, const TaylorModelVec & x0, const std::vector<Interval> & polyRange, const std::vector<HornerForm> & ode, const std::vector<Interval> & step_exp_table, const int numVars, const int order, const Interval & cutoff_threshold) const;
//	void Picard_ctrunc_normal_assign(const TaylorModelVec & x0, const std::vector<Interval> & polyRange, const std::vector<HornerForm> & ode, const std::vector<Interval> & step_exp_table, const int numVars, const int order, const Interval & cutoff_threshold);

//	void Picard_ctrunc_normal(TaylorModelVec & result, std::vector<RangeTree *> & trees, const TaylorModelVec & x0, const std::vector<Interval> & polyRange, const std::vector<HornerForm> & ode, const std::vector<Interval> & step_exp_table, const int numVars, const int order, const Interval & cutoff_threshold, const std::vector<bool> & constant) const;
//	void Picard_ctrunc_normal(TaylorModelVec & result, std::vector<RangeTree *> & trees, const TaylorModelVec & x0, const std::vector<Interval> & polyRange, const std::vector<HornerForm> & ode, const std::vector<Interval> & step_exp_table, const int numVars, const std::vector<int> & orders, const Interval & cutoff_threshold, const std::vector<bool> & constant) const;

//	void Picard_ctrunc_normal(TaylorModelVec & result, std::vector<RangeTree *> & trees, std::vector<std::vector<Interval> > & trunc_parts, const TaylorModelVec & x0, const std::vector<Interval> & polyRange, const std::vector<HornerForm> & ode, const std::vector<std::vector<bool> > & substitution, const std::vector<Interval> & step_exp_table, const int numVars, const int order, const Interval & cutoff_threshold) const;
//	void Picard_only_remainder(std::vector<Interval> & result, std::vector<RangeTree *> & trees, std::vector<std::vector<Interval> > & trunc_parts, const TaylorModelVec & x0, const std::vector<HornerForm> & ode, const Interval & timeStep) const;

//	void Picard_ctrunc_normal_no_cutoff(TaylorModelVec & result, std::vector<RangeTree *> & trees, const TaylorModelVec & x0, const std::vector<Interval> & polyRange, const std::vector<HornerForm> & ode, const std::vector<Interval> & step_exp_table, const int numVars, const int order) const;
//	void Picard_ctrunc_normal_no_cutoff(TaylorModelVec & result, std::vector<RangeTree *> & trees, const TaylorModelVec & x0, const std::vector<Interval> & polyRange, const std::vector<HornerForm> & ode, const std::vector<TaylorModelVec> & sub_tmv, const std::vector<Interval> & step_exp_table, const int numVars, const int order) const;

//	void Picard_only_remainder(std::vector<Interval> & result, std::vector<RangeTree *> & trees, const TaylorModelVec & x0, const std::vector<HornerForm> & ode, const Interval & timeStep, const std::vector<bool> & constant) const;

//	void Picard_no_remainder(TaylorModelVec & result, const TaylorModelVec & x0, const std::vector<HornerForm> & ode, const int numVars, const std::vector<int> & orders, const std::vector<bool> & bIncreased, const Interval & cutoff_threshold) const;
//	void Picard_no_remainder_assign(const TaylorModelVec & x0, const std::vector<HornerForm> & ode, const int numVars, const std::vector<int> & orders, const std::vector<bool> & bIncreased, const Interval & cutoff_threshold);
//	void Picard_ctrunc_normal(TaylorModelVec & result, const TaylorModelVec & x0, const std::vector<Interval> & polyRange, const std::vector<HornerForm> & ode, const std::vector<Interval> & step_exp_table, const int numVars, const std::vector<int> & orders, const Interval & cutoff_threshold) const;
//	void Picard_ctrunc_normal_assign(const TaylorModelVec & x0, const std::vector<Interval> & polyRange, const std::vector<HornerForm> & ode, const std::vector<Interval> & step_exp_table, const int numVars, const std::vector<int> & orders, const Interval & cutoff_threshold);
/*
	// using Taylor approximation
	void Picard_non_polynomial_taylor_no_remainder(TaylorModelVec & result, const TaylorModelVec & x0, const std::vector<std::string> & strOde, const int order, const Interval & cutoff_threshold) const;
	void Picard_non_polynomial_taylor_no_remainder_assign(const TaylorModelVec & x0, const std::vector<std::string> & strOde, const int order, const Interval & cutoff_threshold);

	void Picard_non_polynomial_taylor_no_remainder(TaylorModelVec & result, const TaylorModelVec & x0, const std::vector<std::string> & strOde, const std::vector<int> & orders, const std::vector<bool> & bIncreased, const Interval & cutoff_threshold) const;
	void Picard_non_polynomial_taylor_no_remainder_assign(const TaylorModelVec & x0, const std::vector<std::string> & strOde, const std::vector<int> & orders, const std::vector<bool> & bIncreased, const Interval & cutoff_threshold);

	void Picard_non_polynomial_taylor_ctrunc_normal(TaylorModelVec & result, const TaylorModelVec & x0, const std::vector<std::string> & strOde, const std::vector<Interval> & step_exp_table, const int order, const Interval & cutoff_threshold, const std::vector<bool> & constant, const std::vector<Interval> & constant_part) const;
//	void Picard_non_polynomial_taylor_ctrunc_normal_assign(const TaylorModelVec & x0, const std::vector<std::string> & strOde, const std::vector<Interval> & step_exp_table, const int order, const Interval & cutoff_threshold);

	void Picard_non_polynomial_taylor_ctrunc_normal(TaylorModelVec & result, const TaylorModelVec & x0, const std::vector<std::string> & strOde, const std::vector<Interval> & step_exp_table, const std::vector<int> & orders, const Interval & cutoff_threshold, const std::vector<bool> & constant, const std::vector<Interval> & constant_part) const;
//	void Picard_non_polynomial_taylor_ctrunc_normal_assign(const TaylorModelVec & x0, const std::vector<std::string> & strOde, const std::vector<Interval> & step_exp_table, const std::vector<int> & orders, const Interval & cutoff_threshold);

	void Picard_non_polynomial_taylor_only_remainder(std::vector<Interval> & result, const TaylorModelVec & x0, const std::vector<std::string> & strOde, const Interval & timeStep, const int order, const std::vector<bool> & constant) const;
	void Picard_non_polynomial_taylor_only_remainder(std::vector<Interval> & result, const TaylorModelVec & x0, const std::vector<std::string> & strOde, const Interval & timeStep, const std::vector<int> & orders, const std::vector<bool> & constant) const;
*/

/*
	// ================ Picard operation for the new expression data structure ================

	void Picard_trunc_no_remainder(TaylorModelVec & result, const TaylorModelVec & x0, const std::vector<Expression_AST> & ode, const int numVars, const int order, const Interval & cutoff_threshold, const std::vector<bool> & constant, const std::vector<Interval> & constant_part) const;
	void Picard_trunc_no_remainder_assign(const TaylorModelVec & x0, const std::vector<Expression_AST> & ode, const int numVars, const int order, const Interval & cutoff_threshold, const std::vector<bool> & constant, const std::vector<Interval> & constant_part);
	void Picard_ctrunc_normal(TaylorModelVec & result, const TaylorModelVec & x0, const std::vector<Expression_AST> & ode, const std::vector<Interval> & step_exp_table, const int order, const int numVars, const Interval & cutoff_threshold, const std::vector<bool> & constant, const std::vector<Interval> & constant_part, std::list<Interval> & intermediate_ranges) const;
	void Picard_ctrunc_normal_remainder(std::vector<Interval> & result, const TaylorModelVec & x0, const std::vector<Expression_AST> & ode, const Interval & timeStep, const int order, const std::vector<bool> & constant, std::list<Interval> & intermediate_ranges) const;

	// ========================================================================================
*/

	void normalize(std::vector<Interval> & domain, const Interval & cutoff_threshold);

	template <class DATA_TYPE2>
	void polyRange(std::vector<Interval> & result, const std::vector<DATA_TYPE2> & domain) const;

	template <class DATA_TYPE2>
	void polyRangeNormal(std::vector<Interval> & result, const std::vector<DATA_TYPE2> & step_exp_table) const;

	void get_samples(Matrix<DATA_TYPE> & samples) const;

	friend TaylorModelVec<DATA_TYPE> operator * (const Matrix<DATA_TYPE> & A, const TaylorModelVec<DATA_TYPE> & tmv);
};


template <class DATA_TYPE>
TaylorModelVec<DATA_TYPE>::TaylorModelVec()
{
}

template <class DATA_TYPE>
TaylorModelVec<DATA_TYPE>::TaylorModelVec(const std::vector<TaylorModel<DATA_TYPE> > & tms_input)
{
	tms = tms_input;
}

template <class DATA_TYPE>
TaylorModelVec<DATA_TYPE>::TaylorModelVec(const std::vector<DATA_TYPE> & constants, const unsigned int numVars)
{
	for(unsigned int i=0; i<constants.size(); ++i)
	{
		TaylorModel<DATA_TYPE> tmTemp(constants[i], numVars);
		tms.push_back(tmTemp);
	}
}

template <class DATA_TYPE>
TaylorModelVec<DATA_TYPE>::TaylorModelVec(Matrix<DATA_TYPE> & coefficients)
{
	unsigned int rows = coefficients.rows();
	unsigned int cols = coefficients.cols();

	for(unsigned int i=0; i<rows; ++i)
	{
		DATA_TYPE *p = coefficients.getRowVecRef(i);
		TaylorModel<DATA_TYPE> tmTmp(p, cols);
		tms.push_back(tmTmp);
	}
}

template <class DATA_TYPE>
TaylorModelVec<DATA_TYPE>::TaylorModelVec(const std::vector<DATA_TYPE> & coefficients)
{
	unsigned int numVars = coefficients.size() + 1;
	for(int i=0; i<coefficients.size(); ++i)
	{
		TaylorModel<DATA_TYPE> tmTmp(coefficients[i], numVars);
		tmTmp.expansion.mul_assign(i+1, 1);
		tms.push_back(tmTmp);
	}
}

template <class DATA_TYPE>
TaylorModelVec<DATA_TYPE>::TaylorModelVec(Matrix<DATA_TYPE> & coefficients, const std::vector<Interval> & remainders)
{
	unsigned int rows = coefficients.rows();
	unsigned int cols = coefficients.cols();

	for(unsigned int i=0; i<rows; ++i)
	{
		DATA_TYPE *p = coefficients.getRowVecRef(i);
		TaylorModel<DATA_TYPE> tmTmp(p, cols, remainders[i]);
		tms.push_back(tmTmp);
	}
}

template <class DATA_TYPE>
TaylorModelVec<DATA_TYPE>::TaylorModelVec(const std::vector<std::vector<DATA_TYPE> > & coefficients)
{
	for(unsigned int i=0; i<coefficients.size(); ++i)
	{
		TaylorModel<DATA_TYPE> tmTmp(coefficients[i]);
		tms.push_back(tmTmp);
	}
}

template <class DATA_TYPE>
TaylorModelVec<DATA_TYPE>::TaylorModelVec(const std::vector<std::vector<DATA_TYPE> > & coefficients, const std::vector<Interval> & remainders)
{
	for(unsigned int i=0; i<coefficients.size(); ++i)
	{
		TaylorModel<DATA_TYPE> tmTemp(coefficients[i], remainders[i]);
		tms.push_back(tmTemp);
	}
}

template <class DATA_TYPE>
TaylorModelVec<DATA_TYPE>::TaylorModelVec(const std::vector<Interval> & box, std::vector<Interval> & domain)
{
	unsigned int rangeDim = box.size();
	unsigned int domainDim = rangeDim + 1;
	domain.resize(domainDim, 0);
	domain[0] = 0;

	std::vector<DATA_TYPE> center(rangeDim);
	std::vector<DATA_TYPE> scalars(rangeDim);
	for(unsigned int i=0; i<rangeDim; ++i)
	{
		DATA_TYPE c, r;
		box[i].toCenterForm(c, r);
		center[i] = c;
		scalars[i] = r;
	}

	std::vector<std::vector<DATA_TYPE> > coefficients(rangeDim, std::vector<DATA_TYPE>(domainDim));

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		coefficients[i][i+1] = scalars[i];
	}

	TaylorModelVec<DATA_TYPE> newVars(coefficients);
	tms = newVars.tms;

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		TaylorModel<DATA_TYPE> tmTmp(center[i], domainDim);
		tms[i] += tmTmp;
	}

	Interval intUnit(-1,1);
	for(unsigned int i=1; i<domainDim; ++i)
	{
		domain[i] = intUnit;
	}
}

template <class DATA_TYPE>
TaylorModelVec<DATA_TYPE>::TaylorModelVec(Matrix<Interval> & box, std::vector<Interval> & domain)
{
	unsigned int rangeDim = box.cols();
	unsigned int domainDim = rangeDim + 1;
	domain.resize(domainDim, 0);
	domain[0] = 0;

	// box should always be a column vector
	Interval *p = box.getRowVecRef(0);

	std::vector<DATA_TYPE> center(rangeDim);
	std::vector<DATA_TYPE> scalars(rangeDim);
	for(unsigned int i=0; i<rangeDim; ++i)
	{
		DATA_TYPE c, r;
		(p + i)->toCenterForm(c, r);
		center[i] = c;
		scalars[i] = r;
	}

	std::vector<std::vector<DATA_TYPE> > coefficients(rangeDim, std::vector<DATA_TYPE>(domainDim));

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		coefficients[i][i+1] = scalars[i];
	}

	TaylorModelVec<DATA_TYPE> newVars(coefficients);
	for(unsigned int i=0; i<rangeDim; ++i)
	{
		TaylorModel<DATA_TYPE> tmTmp(center[i], domainDim);
		newVars.tms[i] += tmTmp;
	}

	Interval intUnit(-1,1);
	for(unsigned int i=1; i<domainDim; ++i)
	{
		domain[i] = intUnit;
	}
}

template <class DATA_TYPE>
TaylorModelVec<DATA_TYPE>::TaylorModelVec(const unsigned int dim)
{
	unsigned int domainDim = dim + 1;
	for(unsigned int i=1; i<domainDim; ++i)
	{
		TaylorModel<DATA_TYPE> tmTmp(1, domainDim);
		tmTmp.mul_assign(i, 1);
		tms.push_back(tmTmp);
	}
}

template <class DATA_TYPE>
TaylorModelVec<DATA_TYPE>::TaylorModelVec(const TaylorModelVec<DATA_TYPE> & tmv)
{
	tms = tmv.tms;
}

template <class DATA_TYPE>
TaylorModelVec<DATA_TYPE>::~TaylorModelVec()
{
}

template <class DATA_TYPE>
void TaylorModelVec<DATA_TYPE>::clear()
{
	tms.clear();
}

template <class DATA_TYPE>
void TaylorModelVec<DATA_TYPE>::output(std::ostream & os, const Variables & stateVars, const Variables & tmVars) const
{
	for(unsigned int i=0; i<tms.size(); ++i)
	{
		os << stateVars.varNames[i] << " = ";
		tms[i].output(os, tmVars);
		os << "\n\n";
	}
}

template <class DATA_TYPE>
void TaylorModelVec<DATA_TYPE>::constant(std::vector<DATA_TYPE> & c) const
{
	c.clear();

	for(unsigned int i=0; i<tms.size(); ++i)
	{
		DATA_TYPE tmp;
		tms[i].constant(tmp);
		c.push_back(tmp);
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModelVec<DATA_TYPE>::intEval(std::vector<Interval> & result, const std::vector<DATA_TYPE2> & domain) const
{
	result.clear();
	for(unsigned int i=0; i<tms.size(); ++i)
	{
		Interval I;
		tms[i].intEval(I, domain);
		result.push_back(I);
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModelVec<DATA_TYPE>::intEval(std::vector<Interval> & result, const std::vector<DATA_TYPE2> & domain, const std::vector<unsigned int> & varIDs) const
{
	result.clear();
	for(unsigned int i=0; i<varIDs.size(); ++i)
	{
		Interval I;
		tms[varIDs[i]].intEval(I, domain);
		result.push_back(I);
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModelVec<DATA_TYPE>::intEvalNormal(std::vector<Interval> & result, const std::vector<DATA_TYPE2> & step_exp_table) const
{
	result.clear();
	for(unsigned int i=0; i<tms.size(); ++i)
	{
		Interval I;
		tms[i].intEvalNormal(I, step_exp_table);
		result.push_back(I);
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModelVec<DATA_TYPE>::intEvalNormal(std::vector<Interval> & result, const std::vector<DATA_TYPE2> & step_exp_table, const std::vector<unsigned int> & varIDs) const
{
	result.clear();
	for(unsigned int i=0; i<varIDs.size(); ++i)
	{
		Interval I;
		tms[varIDs[i]].intEvalNormal(I, step_exp_table);
		result.push_back(I);
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModelVec<DATA_TYPE>::ctrunc(const std::vector<DATA_TYPE2> & domain, const unsigned int order)
{
	for(unsigned int i=0; i<tms.size(); ++i)
	{
		tms[i].ctrunc(domain, order);
	}
}

template <class DATA_TYPE>
void TaylorModelVec<DATA_TYPE>::nctrunc(const unsigned int order)
{
	for(unsigned int i=0; i<tms.size(); ++i)
	{
		tms[i].nctrunc(order);
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModelVec<DATA_TYPE>::ctrunc_normal(const std::vector<DATA_TYPE2> & step_exp_table, const unsigned int order)
{
	for(unsigned int i=0; i<tms.size(); ++i)
	{
		tms[i].ctrunc_normal(step_exp_table, order);
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModelVec<DATA_TYPE>::ctrunc(const std::vector<DATA_TYPE2> & domain, const std::vector<unsigned int> & orders)
{
	for(unsigned int i=0; i<tms.size(); ++i)
	{
		tms[i].ctrunc(domain, orders[i]);
	}
}

template <class DATA_TYPE>
void TaylorModelVec<DATA_TYPE>::nctrunc(const std::vector<unsigned int> & orders)
{
	for(unsigned int i=0; i<tms.size(); ++i)
	{
		tms[i].nctrunc(orders[i]);
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModelVec<DATA_TYPE>::ctrunc_normal(const std::vector<DATA_TYPE2> & step_exp_table, const std::vector<unsigned int> & orders)
{
	for(unsigned int i=0; i<tms.size(); ++i)
	{
		tms[i].ctrunc_normal(step_exp_table, orders[i]);
	}
}

template <class DATA_TYPE>
TaylorModelVec<DATA_TYPE> & TaylorModelVec<DATA_TYPE>::operator = (const TaylorModelVec<DATA_TYPE> & tmv)
{
	if(this == &tmv)
		return *this;

	tms = tmv.tms;
	return *this;
}

template <class DATA_TYPE>
TaylorModelVec<DATA_TYPE> & TaylorModelVec<DATA_TYPE>::operator = (const std::vector<Polynomial<DATA_TYPE> > & pv)
{
	if(pv.size() != tms.size())
	{
		printf("Dimensions do not match.\n");
		return *this;
	}
	else
	{
		for(unsigned int i=0; i<pv.size(); ++i)
		{
			tms[i].expansion = pv[i];
		}

		return *this;
	}
}

template <class DATA_TYPE>
TaylorModelVec<DATA_TYPE> & TaylorModelVec<DATA_TYPE>::operator += (const TaylorModelVec<DATA_TYPE> & tmv)
{
	if(tms.size() != tmv.tms.size())
	{
		printf("Dimensions do not match.\n");
		return *this;
	}
	else
	{
		for(unsigned int i=0; i<tms.size(); ++i)
		{
			tms[i] += tmv.tms[i];
		}

		return *this;
	}
}

template <class DATA_TYPE>
TaylorModelVec<DATA_TYPE> & TaylorModelVec<DATA_TYPE>::operator -= (const TaylorModelVec<DATA_TYPE> & tmv)
{
	if(tms.size() != tmv.tms.size())
	{
		printf("Dimensions do not match.\n");
		return *this;
	}
	else
	{
		for(unsigned int i=0; i<tms.size(); ++i)
		{
			tms[i] -= tmv.tms[i];
		}

		return *this;
	}
}

template <class DATA_TYPE>
TaylorModelVec<DATA_TYPE> TaylorModelVec<DATA_TYPE>::operator + (const TaylorModelVec<DATA_TYPE> & tmv) const
{
	TaylorModelVec<DATA_TYPE> result = *this;
	result += tmv;
	return result;
}

template <class DATA_TYPE>
TaylorModelVec<DATA_TYPE> TaylorModelVec<DATA_TYPE>::operator - (const TaylorModelVec<DATA_TYPE> & tmv) const
{
	TaylorModelVec<DATA_TYPE> result = *this;
	result -= tmv;
	return result;
}

template <class DATA_TYPE>
TaylorModelVec<DATA_TYPE> & TaylorModelVec<DATA_TYPE>::operator *= (const DATA_TYPE & c)
{
	if(c == 0)
	{
		this->clear();
		return *this;
	}
	else
	{
		for(unsigned int i=0; i<tms.size(); ++i)
		{
			tms[i] *= c;
		}

		return *this;
	}
}

template <class DATA_TYPE>
TaylorModelVec<DATA_TYPE> & TaylorModelVec<DATA_TYPE>::operator /= (const DATA_TYPE & c)
{
	for(unsigned int i=0; i<tms.size(); ++i)
	{
		tms[i] /= c;
	}

	return *this;
}

template <class DATA_TYPE>
TaylorModelVec<DATA_TYPE> TaylorModelVec<DATA_TYPE>::operator * (const DATA_TYPE & c) const
{
	TaylorModelVec<DATA_TYPE> result = *this;
	result *= c;
	return result;
}

template <class DATA_TYPE>
TaylorModelVec<DATA_TYPE> TaylorModelVec<DATA_TYPE>::operator / (const DATA_TYPE & c) const
{
	TaylorModelVec<DATA_TYPE> result = *this;
	result /= c;
	return result;
}

template <class DATA_TYPE>
void TaylorModelVec<DATA_TYPE>::derivative(TaylorModelVec<DATA_TYPE> & result, const unsigned int varIndex) const
{
	result.clear();
	for(unsigned int i=0; i<tms.size(); ++i)
	{
		TaylorModel<DATA_TYPE> tmTmp;
		tms[i].derivative(tmTmp, varIndex);
		result.tms.push_back(tmTmp);
	}
}

template <class DATA_TYPE>
void TaylorModelVec<DATA_TYPE>::LieDerivative(TaylorModelVec<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & f, const unsigned int order, const Interval & cutoff_threshold) const
{
	result.clear();
	for(unsigned int i=0; i<tms.size(); ++i)
	{
		TaylorModel<DATA_TYPE> tmTmp;
		tms[i].LieDerivative(tmTmp, f, order, cutoff_threshold);
		result.tms.push_back(tmTmp);
	}
}

template <class DATA_TYPE>
void TaylorModelVec<DATA_TYPE>::LieDerivative(TaylorModelVec<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & f, const std::vector<unsigned int> & orders, const Interval & cutoff_threshold) const
{
	result.clear();
	for(unsigned int i=0; i<tms.size(); ++i)
	{
		TaylorModel<DATA_TYPE> tmTmp;
		tms[i].LieDerivative(tmTmp, f, orders[i], cutoff_threshold);
		result.tms.push_back(tmTmp);
	}
}

template <class DATA_TYPE>
void TaylorModelVec<DATA_TYPE>::integral_time(TaylorModelVec<DATA_TYPE> & result, const Interval & I) const
{
	result.clear();
	for(unsigned int i=0; i<tms.size(); ++i)
	{
		TaylorModel<DATA_TYPE> tmTmp;
		tms[i].integral_time(tmTmp, I);
		result.tms.push_back(tmTmp);
	}
}

template <class DATA_TYPE>
void TaylorModelVec<DATA_TYPE>::integral_time(TaylorModelVec<DATA_TYPE> & result) const
{
	result.clear();
	for(unsigned int i=0; i<tms.size(); ++i)
	{
		TaylorModel<DATA_TYPE> tmTmp;
		tms[i].integral_time(tmTmp);
		result.tms.push_back(tmTmp);
	}
}

template <class DATA_TYPE>
void TaylorModelVec<DATA_TYPE>::linearCoefficients(Matrix<DATA_TYPE> & coefficients) const
{
	for(unsigned int i=0; i<tms.size(); ++i)
	{
		tms[i].linearCoefficients(coefficients, i);
	}
}

template <class DATA_TYPE>
void TaylorModelVec<DATA_TYPE>::linearCoefficients(std::vector<std::vector<DATA_TYPE> > & coefficients) const
{
	for(unsigned int i=0; i<tms.size(); ++i)
	{
		tms[i].linearCoefficients(coefficients[i]);
	}
}

template <class DATA_TYPE>
void TaylorModelVec<DATA_TYPE>::rmZeroTerms(const std::vector<unsigned int> & indices)
{
	if(indices.size() != 0)
	{
		for(unsigned int i=0; i<tms.size(); ++i)
		{
			tms[i].rmZeroTerms(indices);
		}
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModelVec<DATA_TYPE>::insert_ctrunc(TaylorModelVec<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & vars, const std::vector<DATA_TYPE2> & varsPolyRange, const std::vector<DATA_TYPE2> & domain, const unsigned int order, const Interval & cutoff_threshold) const
{
	result.clear();

	for(unsigned int i=0; i<tms.size(); ++i)
	{
		TaylorModel<DATA_TYPE> tmTmp;
		tms[i].insert_ctrunc(tmTmp, vars, varsPolyRange, domain, order, cutoff_threshold);
		result.tms.push_back(tmTmp);
	}
}

template <class DATA_TYPE>
void TaylorModelVec<DATA_TYPE>::insert_no_remainder(TaylorModelVec<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & vars, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold) const
{
	result.clear();

	for(unsigned int i=0; i<tms.size(); ++i)
	{
		TaylorModel<DATA_TYPE> tmTmp;
		tms[i].insert_no_remainder(tmTmp, vars, numVars, order, cutoff_threshold);
		result.tms.push_back(tmTmp);
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2, class DATA_TYPE3>
void TaylorModelVec<DATA_TYPE>::insert_ctrunc_normal(TaylorModelVec<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & vars, const std::vector<DATA_TYPE2> & varsPolyRange, const std::vector<DATA_TYPE3> & step_exp_table, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold) const
{
	result.clear();

	for(unsigned int i=0; i<tms.size(); ++i)
	{
		TaylorModel<DATA_TYPE> tmTmp;
		tms[i].insert_ctrunc_normal(tmTmp, vars, varsPolyRange, step_exp_table, numVars, order, cutoff_threshold);
		result.tms.push_back(tmTmp);
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModelVec<DATA_TYPE>::insert_ctrunc(TaylorModelVec<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & vars, const std::vector<DATA_TYPE2> & varsPolyRange, const std::vector<DATA_TYPE2> & domain, const std::vector<unsigned int> & orders, const Interval & cutoff_threshold) const
{
	result.clear();

	for(unsigned int i=0; i<tms.size(); ++i)
	{
		TaylorModel<DATA_TYPE> tmTmp;
		tms[i].insert_ctrunc(tmTmp, vars, varsPolyRange, domain, orders[i], cutoff_threshold);
		result.tms.push_back(tmTmp);
	}
}

template <class DATA_TYPE>
void TaylorModelVec<DATA_TYPE>::insert_no_remainder(TaylorModelVec<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & vars, const unsigned int numVars, const std::vector<unsigned int> & orders, const Interval & cutoff_threshold) const
{
	result.clear();

	for(unsigned int i=0; i<tms.size(); ++i)
	{
		TaylorModel<DATA_TYPE> tmTmp;
		tms[i].insert_no_remainder(tmTmp, vars, numVars, orders[i], cutoff_threshold);
		result.tms.push_back(tmTmp);
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2, class DATA_TYPE3>
void TaylorModelVec<DATA_TYPE>::insert_ctrunc_normal(TaylorModelVec<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & vars, const std::vector<DATA_TYPE2> & varsPolyRange, const std::vector<DATA_TYPE3> & step_exp_table, const unsigned int numVars, const std::vector<unsigned int> & orders, const Interval & cutoff_threshold) const
{
	result.clear();

	for(unsigned int i=0; i<tms.size(); ++i)
	{
		TaylorModel<DATA_TYPE> tmTmp;
		tms[i].insert_ctrunc_normal(tmTmp, vars, varsPolyRange, step_exp_table, numVars, orders[i], cutoff_threshold);
		result.tms.push_back(tmTmp);
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModelVec<DATA_TYPE>::evaluate_time(TaylorModelVec<DATA_TYPE> & result, const std::vector<DATA_TYPE2> & step_exp_table) const
{
	result.clear();

	for(unsigned int i=0; i<tms.size(); ++i)
	{
		TaylorModel<DATA_TYPE> tmTmp;
		tms[i].evaluate_time(tmTmp, step_exp_table);
		result.tms.push_back(tmTmp);
	}
}

template <class DATA_TYPE>
void TaylorModelVec<DATA_TYPE>::scale(TaylorModelVec<DATA_TYPE> & result, const std::vector<DATA_TYPE> & S)
{
	result = *this;
	result.scale_assign(S);
}

template <class DATA_TYPE>
void TaylorModelVec<DATA_TYPE>::scale_assign(const std::vector<DATA_TYPE> & S)
{
	for(unsigned int i=0; i<tms.size(); ++i)
	{
		tms[i] *= S[i];
	}
}

template <class DATA_TYPE>
void TaylorModelVec<DATA_TYPE>::rmConstant()
{
	for(unsigned int i=0; i<tms.size(); ++i)
	{
		tms[i].rmConstant();
	}
}

template <class DATA_TYPE>
void TaylorModelVec<DATA_TYPE>::decompose(TaylorModelVec<DATA_TYPE> & linear, TaylorModelVec<DATA_TYPE> & other) const
{
	linear.tms.clear();
	other.tms.clear();

	for(unsigned int i=0; i<tms.size(); ++i)
	{
		TaylorModel<DATA_TYPE> tm_linear, tm_other;
		tms[i].decompose(tm_linear, tm_other);

		linear.tms.push_back(tm_linear);
		other.tms.push_back(tm_other);
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModelVec<DATA_TYPE>::cutoff_normal(const std::vector<DATA_TYPE2> & step_exp_table, const Interval & cutoff_threshold)
{
	for(unsigned int i=0; i<tms.size(); ++i)
	{
		tms[i].cutoff_normal(step_exp_table, cutoff_threshold);
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModelVec<DATA_TYPE>::cutoff(const std::vector<DATA_TYPE2> & domain, const Interval & cutoff_threshold)
{
	for(unsigned int i=0; i<tms.size(); ++i)
	{
		tms[i].cutoff(domain, cutoff_threshold);
	}
}

template <class DATA_TYPE>
void TaylorModelVec<DATA_TYPE>::cutoff(const Interval & cutoff_threshold)
{
	for(unsigned int i=0; i<tms.size(); ++i)
	{
		tms[i].cutoff(cutoff_threshold);
	}
}
/*
template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModelVec<DATA_TYPE>::cutoff_normal(Matrix<Interval> & M, const std::vector<DATA_TYPE2> & step_exp_table, const Interval & cutoff_threshold)
{
	int d = tms.size();
	Matrix<Interval> tmp(d, 1);

	for(int i=0; i<tms.size(); ++i)
	{
		Interval I;
		tms[i].expansion.cutoff_normal(I, step_exp_table, cutoff_threshold);
		tmp[i][0] = I;
		tms[i].remainder += I;
	}

	M = tmp;
}
*/
template <class DATA_TYPE>
void TaylorModelVec<DATA_TYPE>::Expansion(std::vector<Polynomial<DATA_TYPE> > & polys) const
{
	polys.clear();

	for(unsigned int i=0; i<tms.size(); ++i)
	{
		polys.push_back(tms[i].expansion);
	}
}

template <class DATA_TYPE>
void TaylorModelVec<DATA_TYPE>::Remainder(Matrix<Interval> & rem) const
{
	Interval *p = rem.getRowVecRef(0);
	for(int i=0; i<tms.size(); ++i)
	{
		*(p + i) = tms[i].remainder;
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModelVec<DATA_TYPE>::Picard_no_remainder(TaylorModelVec<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & x0, const std::vector<HornerForm<DATA_TYPE2> > & ode, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold) const
{
	TaylorModelVec<DATA_TYPE> tmvTmp;

	unsigned int k = order > 1 ? (order - 1) : 1;

	for(unsigned int i=0; i<ode.size(); ++i)
	{
		TaylorModel<DATA_TYPE> tmTmp;
		ode[i].insert_no_remainder(tmTmp, *this, numVars, k, cutoff_threshold);
		tmvTmp.tms.push_back(tmTmp);
	}

	TaylorModelVec<DATA_TYPE> tmvTmp2;
	tmvTmp.integral_time(tmvTmp2);

	result = x0 + tmvTmp2;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModelVec<DATA_TYPE>::Picard_no_remainder_assign(const TaylorModelVec<DATA_TYPE> & x0, const std::vector<HornerForm<DATA_TYPE2> > & ode, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold)
{
	TaylorModelVec<DATA_TYPE> result;
	Picard_no_remainder(result, x0, ode, numVars, order, cutoff_threshold);
	*this = result;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModelVec<DATA_TYPE>::Picard_no_remainder(TaylorModelVec<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & x0, const std::vector<HornerForm<DATA_TYPE2> > & ode, const unsigned int numVars, const std::vector<unsigned int> & orders, const Interval & cutoff_threshold) const
{
	TaylorModelVec<DATA_TYPE> tmvTmp;

	for(unsigned int i=0; i<ode.size(); ++i)
	{
		unsigned int k = orders[i] > 1 ? (orders[i] - 1) : 1;

		TaylorModel<DATA_TYPE> tmTmp;
		ode[i].insert_no_remainder(tmTmp, *this, numVars, k, cutoff_threshold);
		tmvTmp.tms.push_back(tmTmp);
	}

	TaylorModelVec<DATA_TYPE> tmvTmp2;
	tmvTmp.integral_time(tmvTmp2);

	result = x0 + tmvTmp2;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModelVec<DATA_TYPE>::Picard_no_remainder_assign(const TaylorModelVec<DATA_TYPE> & x0, const std::vector<HornerForm<DATA_TYPE2> > & ode, const unsigned int numVars, const std::vector<unsigned int> & orders, const Interval & cutoff_threshold)
{
	TaylorModelVec<DATA_TYPE> result;
	Picard_no_remainder(result, x0, ode, numVars, orders, cutoff_threshold);
	*this = result;
}

template <class DATA_TYPE>
template <class DATA_TYPE2, class DATA_TYPE3>
void TaylorModelVec<DATA_TYPE>::Picard_ctrunc_normal(TaylorModelVec<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & x0, const std::vector<HornerForm<DATA_TYPE2> > & ode, const std::vector<DATA_TYPE3> & varsPolyRange, const std::vector<DATA_TYPE3> & step_exp_table, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold, std::list<Interval> & intermediate_ranges) const
{
	TaylorModelVec<DATA_TYPE> tmvTmp;

	unsigned int k = order - 1;

	for(unsigned int i=0; i<ode.size(); ++i)
	{
		TaylorModel<DATA_TYPE> tmTmp;
		ode[i].insert_ctrunc_normal(tmTmp, intermediate_ranges, *this, varsPolyRange, step_exp_table, numVars, k, cutoff_threshold);
		tmvTmp.tms.push_back(tmTmp);
	}

	TaylorModelVec<DATA_TYPE> tmvTmp2;
	tmvTmp.integral_time(tmvTmp2, step_exp_table[1]);
	result = x0 + tmvTmp2;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModelVec<DATA_TYPE>::Picard_ctrunc_normal_remainder(std::vector<Interval> & result, const std::vector<HornerForm<DATA_TYPE2> > & ode, const Interval & timeStep, std::list<Interval> & intermediate_ranges) const
{
	std::list<Interval>::iterator iter = intermediate_ranges.begin();

	result.clear();

	for(int i=0; i<ode.size(); ++i)
	{
		Interval intTmp;
		ode[i].insert_only_remainder(intTmp, iter, *this, timeStep);

		intTmp *= timeStep;
		result.push_back(intTmp);
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2, class DATA_TYPE3>
void TaylorModelVec<DATA_TYPE>::Picard_ctrunc_normal(TaylorModelVec<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & x0, const std::vector<HornerForm<DATA_TYPE2> > & ode, const std::vector<DATA_TYPE3> & varsPolyRange, const std::vector<DATA_TYPE3> & step_exp_table, const unsigned int numVars, const std::vector<unsigned int> & orders, const Interval & cutoff_threshold, std::list<Interval> & intermediate_ranges) const
{
	TaylorModelVec<DATA_TYPE> tmvTmp;

	for(unsigned int i=0; i<ode.size(); ++i)
	{
		unsigned int k = orders[i] - 1;

		TaylorModel<DATA_TYPE> tmTmp;
		ode[i].insert_ctrunc_normal(tmTmp, intermediate_ranges, *this, varsPolyRange, step_exp_table, numVars, k, cutoff_threshold);
		tmvTmp.tms.push_back(tmTmp);
	}

	TaylorModelVec<DATA_TYPE> tmvTmp2;
	tmvTmp.integral_time(tmvTmp2, step_exp_table[1]);

	result = x0 + tmvTmp2;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModelVec<DATA_TYPE>::Picard_no_remainder(TaylorModelVec<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & x0, const std::vector<Expression_AST<DATA_TYPE2> > & ode, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold) const
{
	TaylorModelVec<DATA_TYPE> tmvTmp;

	unsigned int k = order - 1;

	for(int i=0; i<ode.size(); ++i)
	{
		TaylorModel<DATA_TYPE> tmTmp;
		ode[i].evaluate_no_remainder(tmTmp, this->tms, k, cutoff_threshold, numVars);
		tmvTmp.tms.push_back(tmTmp);
	}

	TaylorModelVec tmvTmp2;
	tmvTmp.integral_time(tmvTmp2);

	result = x0 + tmvTmp2;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModelVec<DATA_TYPE>::Picard_no_remainder_assign(const TaylorModelVec<DATA_TYPE> & x0, const std::vector<Expression_AST<DATA_TYPE2> > & ode, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold)
{
	TaylorModelVec<DATA_TYPE> result;
	Picard_no_remainder(result, x0, ode, numVars, order, cutoff_threshold);
	*this = result;
}

template <class DATA_TYPE>
template <class DATA_TYPE2, class DATA_TYPE3>
void TaylorModelVec<DATA_TYPE>::Picard_ctrunc_normal(TaylorModelVec<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & x0, const std::vector<Expression_AST<DATA_TYPE2> > & ode, const std::vector<DATA_TYPE3> & step_exp_table, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold, std::list<Interval> & intermediate_ranges, const Global_Computation_Setting & setting) const
{
	TaylorModelVec<DATA_TYPE> tmvTmp;

	unsigned int k = order - 1;

	for(int i=0; i<ode.size(); ++i)
	{
		TaylorModel<DATA_TYPE> tmTmp;
		ode[i].evaluate(tmTmp, this->tms, k, step_exp_table, cutoff_threshold, numVars, intermediate_ranges, setting);
		tmvTmp.tms.push_back(tmTmp);
	}

	TaylorModelVec tmvTmp2;
	tmvTmp.integral_time(tmvTmp2, step_exp_table[1]);

	result = x0 + tmvTmp2;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModelVec<DATA_TYPE>::Picard_ctrunc_normal_remainder(std::vector<Interval> & result, const std::vector<Expression_AST<DATA_TYPE2> > & ode, const Interval & timeStep, const unsigned int order, std::list<Interval> & intermediate_ranges, const Global_Computation_Setting & setting) const
{
	std::list<Interval>::iterator iter = intermediate_ranges.begin();

	result.clear();

	unsigned int k = order - 1;

	for(int i=0; i<ode.size(); ++i)
	{
		Interval intTmp;
		ode[i].evaluate_remainder(intTmp, this->tms, k, iter, setting);

		intTmp *= timeStep;
		result.push_back(intTmp);
	}
}


template <class DATA_TYPE>
void TaylorModelVec<DATA_TYPE>::normalize(std::vector<Interval> & domain, const Interval & cutoff_threshold)
{
	unsigned int domainDim = domain.size();
	unsigned int rangeDim = domainDim - 1;

	// compute the center of the original domain and make it origin-centered
	std::vector<Real> center(domainDim);
	for(unsigned int i=1; i<domainDim; ++i)		// we omit the time dimension
	{
		Real c;
		domain[i].remove_midpoint(c);
		center[i] = c;
	}

	// compute the scalars
	std::vector<std::vector<DATA_TYPE> > coefficients(rangeDim, std::vector<DATA_TYPE>(domainDim));

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		DATA_TYPE c;
		domain[i].mag(c);
		coefficients[i][i+1] = c;
	}

	TaylorModelVec<DATA_TYPE> newVars(coefficients);
	for(unsigned int i=0; i<rangeDim; ++i)
	{
		TaylorModel<DATA_TYPE> tmTemp(center[i], domainDim);
		newVars.tms[i] += tmTemp;
	}

	Interval intUnit(-1,1);
	for(unsigned int i=1; i<domainDim; ++i)
	{
		domain[i] = intUnit;
	}

	for(unsigned int i=0; i<tms.size(); ++i)
	{
		TaylorModel<DATA_TYPE> tmTmp;
		tms[i].insert_no_remainder(tmTmp, newVars, domainDim, tms[i].degree(), cutoff_threshold);
		tms[i].expansion = tmTmp.expansion;
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModelVec<DATA_TYPE>::polyRange(std::vector<Interval> & result, const std::vector<DATA_TYPE2> & domain) const
{
	result.clear();

	for(unsigned int i=0; i<tms.size(); ++i)
	{
		Interval intTemp;
		tms[i].polyRange(intTemp, domain);
		result.push_back(intTemp);
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void TaylorModelVec<DATA_TYPE>::polyRangeNormal(std::vector<Interval> & result, const std::vector<DATA_TYPE2> & step_exp_table) const
{
	result.clear();

	for(unsigned int i=0; i<tms.size(); ++i)
	{
		Interval intTemp;
		tms[i].polyRangeNormal(intTemp, step_exp_table);
		result.push_back(intTemp);
	}
}

template <class DATA_TYPE>
void TaylorModelVec<DATA_TYPE>::get_samples(Matrix<DATA_TYPE> & samples) const
{
	int rangeDim = tms.size();
	int col = 1;

	// samples should be given as a (1 + 2*rangeDim) \times rangeDim real matrix

	// center point
	for(int i=0; i<rangeDim; ++i)
	{
		DATA_TYPE tmp;
		tms[i].constant(tmp);
		samples[i][0] = tmp;
	}

	std::vector<HornerForm<DATA_TYPE> > hfs;
	for(int i=0; i<rangeDim; ++i)
	{
		HornerForm<DATA_TYPE> hf;
		tms[i].expansion.toHornerForm(hf);
		hfs.push_back(hf);
	}

	// prepare the domain
	std::vector<Interval> domain(rangeDim + 1);

	// center point of every facet of the domain box
	int last_pos = -1;
	for(int i=1; i<=rangeDim; ++i)
	{
		if(last_pos > 0)
		{
			domain[last_pos] = 0;
		}

		domain[i] = 1;
		last_pos = i;

		// 0 ... 1 ... 0
		for(int j=0; j<rangeDim; ++j)
		{
			Interval tmp;
			DATA_TYPE sample_j;
			hfs[j].evaluate(tmp, domain);
			tmp.sup(sample_j);

			samples[j][col] = sample_j;
		}

		++col;

		// 0 ... -1 ... 0
		domain[i] = -1;
		for(int j=0; j<rangeDim; ++j)
		{
			Interval tmp;
			DATA_TYPE sample_j;
			hfs[j].evaluate(tmp, domain);
			tmp.midpoint(sample_j);

			samples[j][col] = sample_j;
		}

		++col;
	}
}

template <class DATA_TYPE>
double TaylorModelVec<DATA_TYPE>::rho(const std::vector<Real> & l, const std::vector<Interval> & domain) const
{
	unsigned int d = l.size();
	TaylorModel<DATA_TYPE> tmp1;

	for(int i=0; i<d; ++i)
	{
		TaylorModel<DATA_TYPE> tmp2 = tms[i] * l[i];
		tmp1 += tmp2;
	}

	Interval range;
	tmp1.intEval(range, domain);

	return range.sup();
}

template <class DATA_TYPE>
double TaylorModelVec<DATA_TYPE>::rho_normal(const std::vector<Real> & l, const std::vector<Interval> & step_exp_table) const
{
	unsigned int d = l.size();
	TaylorModel<DATA_TYPE> tmp1;

	for(int i=0; i<d; ++i)
	{
		TaylorModel<DATA_TYPE> tmp2 = tms[i] * l[i];
		tmp1 += tmp2;
	}

	Interval range;
	tmp1.intEvalNormal(range, step_exp_table);

	return range.sup();
}


template <class DATA_TYPE>
TaylorModelVec<DATA_TYPE> operator * (const Matrix<DATA_TYPE> & A, const TaylorModelVec<DATA_TYPE> & tmv)
{
	if(A.size2 != tmv.tms.size())
	{
		printf("Matrix multiplication: Dimensions do not match.\n");
		exit(1);
	}

	TaylorModelVec<DATA_TYPE> result;

	for(unsigned int i=0, pos=0; i<A.size1; ++i, pos+=A.size2)
	{
		TaylorModel<DATA_TYPE> tmTmp;

		for(unsigned int j=0; j<A.size2; ++j)
		{
			tmTmp += tmv.tms[j] * A.data[pos + j];
		}

		result.tms.push_back(tmTmp);
	}

	return result;
}






















// c + (...)*x1 + (...)*x2 + ... + (...)*xn
template <class DATA_TYPE>
class HornerForm
{
protected:
	DATA_TYPE constant_part;
	std::vector<HornerForm<DATA_TYPE> > nonconstant_part;

public:
	HornerForm();
	HornerForm(const DATA_TYPE & c);
	HornerForm(const HornerForm<DATA_TYPE> & hf);
	~HornerForm();

	void clear();

	template <class DATA_TYPE2, class DATA_TYPE3>
	void evaluate(DATA_TYPE2 & result, const std::vector<DATA_TYPE3> & domain) const;

//	template <class DATA_TYPE3, class DATA_TYPE4>
//	void insert(TaylorModel & result, const TaylorModelVec & vars, const std::vector<DATA_TYPE3> & varsPolyRange, const std::vector<DATA_TYPE4> & domain, const Interval & cutoff_threshold) const;

//	template <class DATA_TYPE3, class DATA_TYPE4>
//	void insert_normal(TaylorModel & result, const TaylorModelVec & vars, const std::vector<DATA_TYPE3> & varsPolyRange, const std::vector<DATA_TYPE4> & step_exp_table, const unsigned int numVars, const Interval & cutoff_threshold) const;


	template <class DATA_TYPE2, class DATA_TYPE3>
	void insert_ctrunc_normal(TaylorModel<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & vars, const std::vector<DATA_TYPE2> & varsPolyRange, const std::vector<DATA_TYPE3> & step_exp_table, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold) const;

	template <class DATA_TYPE2, class DATA_TYPE3, class DATA_TYPE4>
	void insert_ctrunc_normal(TaylorModel<DATA_TYPE2> & result, std::list<Interval> & intermediate_ranges, const TaylorModelVec<DATA_TYPE2> & vars, const std::vector<DATA_TYPE3> & varsPolyRange, const std::vector<DATA_TYPE4> & step_exp_table, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold) const;

	template <class DATA_TYPE2>
	void insert_ctrunc(TaylorModel<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & vars, const std::vector<DATA_TYPE2> & varsPolyRange, const std::vector<DATA_TYPE2> & domain, const unsigned int order, const Interval & cutoff_threshold) const;

	template <class DATA_TYPE2>
	void insert_no_remainder(TaylorModel<DATA_TYPE2> & result, const TaylorModelVec<DATA_TYPE2> & vars, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold) const;

	template <class DATA_TYPE2, class DATA_TYPE3>
	void insert_only_remainder(Interval & result, std::list<Interval>::iterator & range_iter, const TaylorModelVec<DATA_TYPE2> & vars, const DATA_TYPE3 & timeStep) const;

	HornerForm<DATA_TYPE> & operator = (const HornerForm<DATA_TYPE> & hf);

	void output(std::ostream & os, const Variables & vars) const;

	template <class DATA_TYPE2>
	friend class HornerForm;

	template <class DATA_TYPE2>
	friend class Polynomial;

	template <class DATA_TYPE2>
	friend class TaylorModel;

	template <class DATA_TYPE2>
	friend class TaylorModelVec;
};


template <class DATA_TYPE>
HornerForm<DATA_TYPE>::HornerForm()
{
}

template <class DATA_TYPE>
HornerForm<DATA_TYPE>::HornerForm(const DATA_TYPE & c)
{
	constant_part = c;
}

template <class DATA_TYPE>
HornerForm<DATA_TYPE>::HornerForm(const HornerForm<DATA_TYPE> & hf)
{
	constant_part = hf.constant_part;
	nonconstant_part = hf.nonconstant_part;
}

template <class DATA_TYPE>
HornerForm<DATA_TYPE>::~HornerForm()
{
}

template <class DATA_TYPE>
void HornerForm<DATA_TYPE>::clear()
{
	constant_part = 0;
	nonconstant_part.clear();
}

template <class DATA_TYPE>
template <class DATA_TYPE2, class DATA_TYPE3>
void HornerForm<DATA_TYPE>::evaluate(DATA_TYPE2 & result, const std::vector<DATA_TYPE3> & domain) const
{
	result = constant_part;

	for(unsigned int i=0; i<nonconstant_part.size(); ++i)
	{
		DATA_TYPE2 tmp;
		nonconstant_part[i].evaluate(tmp, domain);
		tmp *= domain[i];
		result += tmp;
	}
}

/*
template <class DATA_TYPE>
template <class DATA_TYPE3, class DATA_TYPE4>
void HornerForm<DATA_TYPE>::insert(TaylorModel & result, const TaylorModelVec & vars, const std::vector<DATA_TYPE3> & varsPolyRange, const std::vector<DATA_TYPE4> & domain, const Interval & cutoff_threshold) const
{
	unsigned int numVars = domain.size();

	result.clear();

	if(constant_part != 0)
	{
		TaylorModel tmp(constant_part, numVars);
		result = tmp;
	}

	if(nonconstant_part.size() > 0)						// the first variable is t
	{
		TaylorModel tmp;
		nonconstant_part[0].insert(tmp, vars, varsPolyRange, domain, cutoff_threshold);

		tmp.expansion.mul_assign(0,1);					// multiplied by t
		tmp.remainder *= domain[0];
		result.add_assign(tmp);

		for(int i=1; i<nonconstant_part.size(); ++i)
		{
			nonconstant_part[i].insert(tmp, vars, varsPolyRange, domain, cutoff_threshold);	// recursive call
			tmp.mul_insert_assign(vars.tms[i-1], varsPolyRange[i-1], domain, cutoff_threshold);
			result.add_assign(tmp);
		}
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE3, class DATA_TYPE4>
void HornerForm<DATA_TYPE>::insert_normal(TaylorModel & result, const TaylorModelVec & vars, const std::vector<DATA_TYPE3> & varsPolyRange, const std::vector<DATA_TYPE4> & step_exp_table, const unsigned int numVars, const Interval & cutoff_threshold) const
{
	result.clear();

	if(constant_part != 0)
	{
		TaylorModel tmp(constant_part, numVars);
		result = tmp;
	}

	if(nonconstant_part.size() > 0)						// the first variable is t
	{
		TaylorModel tmp;
		nonconstant_part[0].insert_normal(tmp, vars, varsPolyRange, step_exp_table, numVars, cutoff_threshold);

		tmp.expansion.mul_assign(0,1);					// multiplied by t
		tmp.remainder *= step_exp_table[1];
		result.add_assign(tmp);

		for(int i=1; i<nonconstant_part.size(); ++i)
		{
			nonconstant_part[i].insert_normal(tmp, vars, varsPolyRange, step_exp_table, numVars, cutoff_threshold);	// recursive call
			tmp.mul_insert_normal_assign(vars.tms[i-1], varsPolyRange[i-1], step_exp_table, cutoff_threshold);
			result.add_assign(tmp);
		}
	}
}
*/

template <class DATA_TYPE>
template <class DATA_TYPE2, class DATA_TYPE3>
void HornerForm<DATA_TYPE>::insert_ctrunc_normal(TaylorModel<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & vars, const std::vector<DATA_TYPE2> & varsPolyRange, const std::vector<DATA_TYPE3> & step_exp_table, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold) const
{
	result.clear();

	if(constant_part != 0)
	{
		TaylorModel<DATA_TYPE> tmp(constant_part, numVars);
		result = tmp;
	}

	if(nonconstant_part.size() > 0)				// the first variable is t
	{
		TaylorModel<DATA_TYPE> tmp;
		nonconstant_part[0].insert_ctrunc_normal(tmp, vars, varsPolyRange, step_exp_table, numVars, order, cutoff_threshold);

		tmp.expansion.mul_assign(0,1);			// multiplied by t
		tmp.remainder *= step_exp_table[1];

		tmp.ctrunc_normal(step_exp_table, order);
		result += tmp;

		for(int i=1; i<nonconstant_part.size(); ++i)
		{
			nonconstant_part[i].insert_ctrunc_normal(tmp, vars, varsPolyRange, step_exp_table, numVars, order, cutoff_threshold);	// recursive call

			tmp.mul_insert_ctrunc_normal_assign(vars.tms[i-1], varsPolyRange[i-1], step_exp_table, order, cutoff_threshold);

			result += tmp;
		}
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2, class DATA_TYPE3, class DATA_TYPE4>
void HornerForm<DATA_TYPE>::insert_ctrunc_normal(TaylorModel<DATA_TYPE2> & result, std::list<Interval> & intermediate_ranges, const TaylorModelVec<DATA_TYPE2> & vars, const std::vector<DATA_TYPE3> & varsPolyRange, const std::vector<DATA_TYPE4> & step_exp_table, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold) const
{
	result.clear();

	if(constant_part != 0)
	{
		TaylorModel<DATA_TYPE2> tmp(constant_part, numVars);
		result = tmp;
	}

	if(nonconstant_part.size() > 0)					// the first variable is t
	{
		TaylorModel<DATA_TYPE2> tmp;

		nonconstant_part[0].insert_ctrunc_normal(tmp, intermediate_ranges, vars, varsPolyRange, step_exp_table, numVars, order, cutoff_threshold);

		tmp.expansion.mul_assign(0,1);				// multiplied by t
		tmp.remainder *= step_exp_table[1];

		Interval intTrunc;
		tmp.expansion.ctrunc_normal(intTrunc, step_exp_table, order);
		tmp.remainder += intTrunc;

		intermediate_ranges.push_back(intTrunc);

		result += tmp;

		for(int i=1; i<nonconstant_part.size(); ++i)
		{
			TaylorModel<DATA_TYPE2> tmp;

			nonconstant_part[i].insert_ctrunc_normal(tmp, intermediate_ranges, vars, varsPolyRange, step_exp_table, numVars, order, cutoff_threshold);	// recursive call

			Interval tm1, intTrunc2;
			tmp.mul_insert_ctrunc_normal_assign(tm1, intTrunc2, vars.tms[i-1], varsPolyRange[i-1], step_exp_table, order, cutoff_threshold); 	// here coefficient_range = tm1

			intermediate_ranges.push_back(tm1);
			intermediate_ranges.push_back(varsPolyRange[i-1]);
			intermediate_ranges.push_back(intTrunc2);

			result += tmp;
		}
	}
}

template <>
template <>
inline void HornerForm<Interval>::insert_ctrunc_normal<Real, Interval, Interval>(TaylorModel<Real> & result, std::list<Interval> & intermediate_ranges, const TaylorModelVec<Real> & vars, const std::vector<Interval> & varsPolyRange, const std::vector<Interval> & step_exp_table, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold) const
{
	result.clear();

	if(constant_part != 0)
	{
		Real c;
		Interval I = constant_part;
		I.remove_midpoint(c);
		TaylorModel<Real> tmp(c, numVars);
		result = tmp;
		intermediate_ranges.push_back(I);
	}
	else
	{
		Interval invalidInt(1,-1);
		intermediate_ranges.push_back(invalidInt);
	}

	if(nonconstant_part.size() > 0)					// the first variable is t
	{
		TaylorModel<Real> tmp;

		nonconstant_part[0].insert_ctrunc_normal(tmp, intermediate_ranges, vars, varsPolyRange, step_exp_table, numVars, order, cutoff_threshold);

		tmp.expansion.mul_assign(0,1);				// multiplied by t
		tmp.remainder *= step_exp_table[1];

		Interval intTrunc;
		tmp.expansion.ctrunc_normal(intTrunc, step_exp_table, order);
		tmp.remainder += intTrunc;

		intermediate_ranges.push_back(intTrunc);

		result += tmp;

		for(int i=1; i<nonconstant_part.size(); ++i)
		{
			TaylorModel<Real> tmp;

			nonconstant_part[i].insert_ctrunc_normal(tmp, intermediate_ranges, vars, varsPolyRange, step_exp_table, numVars, order, cutoff_threshold);	// recursive call

			Interval tm1, intTrunc2;
			tmp.mul_insert_ctrunc_normal_assign(tm1, intTrunc2, vars.tms[i-1], varsPolyRange[i-1], step_exp_table, order, cutoff_threshold); 	// here coefficient_range = tm1

			intermediate_ranges.push_back(tm1);
			intermediate_ranges.push_back(varsPolyRange[i-1]);
			intermediate_ranges.push_back(intTrunc2);

			result += tmp;
		}
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void HornerForm<DATA_TYPE>::insert_ctrunc(TaylorModel<DATA_TYPE> & result, const TaylorModelVec<DATA_TYPE> & vars, const std::vector<DATA_TYPE2> & varsPolyRange, const std::vector<DATA_TYPE2> & domain, const unsigned int order, const Interval & cutoff_threshold) const
{
	unsigned int numVars = domain.size();

	result.clear();

	if(constant_part != 0)
	{
		TaylorModel<DATA_TYPE> tmp(constant_part, numVars);
		result = tmp;
	}

	if(nonconstant_part.size() > 0)						// the first variable is t
	{
		TaylorModel<DATA_TYPE> tmp;
		nonconstant_part[0].insert_ctrunc(tmp, vars, varsPolyRange, domain, order, cutoff_threshold);

		tmp.expansion.mul_assign(0,1);					// multiplied by t
		tmp.remainder *= domain[0];

		tmp.ctrunc(domain, order);
		result += tmp;

		for(int i=1; i<nonconstant_part.size(); ++i)
		{
			nonconstant_part[i].insert_ctrunc(tmp, vars, varsPolyRange, domain, order, cutoff_threshold);
			tmp.mul_insert_ctrunc_assign(vars.tms[i-1], varsPolyRange[i-1], domain, order, cutoff_threshold);
			result += tmp;
		}
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void HornerForm<DATA_TYPE>::insert_no_remainder(TaylorModel<DATA_TYPE2> & result, const TaylorModelVec<DATA_TYPE2> & vars, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold) const
{
	result.clear();

	if(constant_part != 0)
	{
		TaylorModel<DATA_TYPE2> tmp(constant_part, numVars);
		result = tmp;
	}

	if(nonconstant_part.size() > 0)						// the first variable is t
	{
		TaylorModel<DATA_TYPE2> tmp;
		nonconstant_part[0].insert_no_remainder(tmp, vars, numVars, order, cutoff_threshold);

		tmp.expansion.mul_assign(0,1);					// multiplied by t
		tmp.nctrunc(order);
		result.expansion += tmp.expansion;

		for(int i=1; i<nonconstant_part.size(); ++i)
		{
			nonconstant_part[i].insert_no_remainder(tmp, vars, numVars, order, cutoff_threshold);
			tmp.mul_no_remainder_assign(vars.tms[i-1], order, cutoff_threshold);
			result.expansion += tmp.expansion;
		}
	}
}

template <>
template <>
inline void HornerForm<Interval>::insert_no_remainder<Real>(TaylorModel<Real> & result, const TaylorModelVec<Real> & vars, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold) const
{
	result.clear();

	if(constant_part != 0)
	{
		TaylorModel<Real> tmp(constant_part.toReal(), numVars);
		result = tmp;
	}

	if(nonconstant_part.size() > 0)						// the first variable is t
	{
		TaylorModel<Real> tmp;
		nonconstant_part[0].insert_no_remainder(tmp, vars, numVars, order, cutoff_threshold);

		tmp.expansion.mul_assign(0,1);					// multiplied by t
		tmp.nctrunc(order);
		result.expansion += tmp.expansion;

		for(int i=1; i<nonconstant_part.size(); ++i)
		{
			nonconstant_part[i].insert_no_remainder(tmp, vars, numVars, order, cutoff_threshold);
			tmp.mul_no_remainder_assign(vars.tms[i-1], order, cutoff_threshold);
			result.expansion += tmp.expansion;
		}
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2, class DATA_TYPE3>
void HornerForm<DATA_TYPE>::insert_only_remainder(Interval & result, std::list<Interval>::iterator & range_iter, const TaylorModelVec<DATA_TYPE2> & vars, const DATA_TYPE3 & timeStep) const
{
	result = 0;

	if(nonconstant_part.size() > 0)						// the first variable is t
	{
		Interval tmp1;
		nonconstant_part[0].insert_only_remainder(tmp1, range_iter, vars, timeStep);
		tmp1 *= timeStep;

		tmp1 += (*range_iter);
		result += tmp1;

		++range_iter;

		for(int i=1; i<nonconstant_part.size(); ++i)
		{
			Interval tmp2;
			nonconstant_part[i].insert_only_remainder(tmp2, range_iter, vars, timeStep);

			Interval newRemainder = (*range_iter) * vars.tms[i-1].remainder;
			++range_iter;
			newRemainder += (*range_iter) * tmp2;
			newRemainder += vars.tms[i-1].remainder * tmp2;
			++range_iter;
			newRemainder += (*range_iter);

			result += newRemainder;
			++range_iter;
		}
	}
}

template <>
template <>
inline void HornerForm<Interval>::insert_only_remainder<Real, Interval>(Interval & result, std::list<Interval>::iterator & range_iter, const TaylorModelVec<Real> & vars, const Interval & timeStep) const
{
	if(range_iter->valid())
	{
		result = *range_iter;
	}
	else
	{
		result = 0;
	}

	++range_iter;

	if(nonconstant_part.size() > 0)						// the first variable is t
	{
		Interval tmp1;
		nonconstant_part[0].insert_only_remainder(tmp1, range_iter, vars, timeStep);
		tmp1 *= timeStep;

		tmp1 += (*range_iter);
		result += tmp1;

		++range_iter;

		for(int i=1; i<nonconstant_part.size(); ++i)
		{
			Interval tmp2;
			nonconstant_part[i].insert_only_remainder(tmp2, range_iter, vars, timeStep);

			Interval newRemainder = (*range_iter) * vars.tms[i-1].remainder;
			++range_iter;
			newRemainder += (*range_iter) * tmp2;
			newRemainder += vars.tms[i-1].remainder * tmp2;
			++range_iter;
			newRemainder += (*range_iter);

			result += newRemainder;
			++range_iter;
		}
	}
}

template <class DATA_TYPE>
HornerForm<DATA_TYPE> & HornerForm<DATA_TYPE>::operator = (const HornerForm<DATA_TYPE> & hf)
{
	if(this == &hf)
		return *this;

	constant_part = hf.constant_part;
	nonconstant_part = hf.nonconstant_part;

	return *this;
}

template <class DATA_TYPE>
void HornerForm<DATA_TYPE>::output(std::ostream & os, const Variables & vars) const
{
	unsigned int numVars = nonconstant_part.size();

	bool bPlus = false;

	os << " ( ";
	if(constant_part != 0)
	{
		bPlus = true;
		os << constant_part;
	}

	if(numVars == 0)
	{
		os << " ) ";
		return;
	}

	for(int i=0; i<numVars; ++i)
	{
		if(nonconstant_part[i].nonconstant_part.size() != 0 || nonconstant_part[i].constant_part != 0)
		{
			if(bPlus)
				os << " + ";
			else
				bPlus = true;

			nonconstant_part[i].output(os, vars);
			os << " * " << vars.varNames[i];
		}
	}

	os << " ) ";
}


}



#endif /* TAYLORMODEL_H_ */
