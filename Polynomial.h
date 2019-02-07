/*---
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.
  The code is released as is under the GNU General Public License (GPL).
---*/

#ifndef POLYNOMIAL_H_
#define POLYNOMIAL_H_

#include "Term.h"
#include "Matrix.h"

void parseMultivariatePolynomial();

using namespace flowstar;
extern Variables stateVars;
extern Variables tmVars;

namespace flowstar
{


template <class DATA_TYPE>
class Matrix;

template <class DATA_TYPE>
class HornerForm;

template <class DATA_TYPE>
class Polynomial														// polynomials in monomial form
{
protected:
	std::list<Term<DATA_TYPE> > terms;

public:
	Polynomial();														// empty polynomial
	Polynomial(const DATA_TYPE & c, const unsigned int numVars);
	Polynomial(Matrix<DATA_TYPE> & coefficients);						// linear polynomial with the given coefficients, the input matrix is a row vector
	Polynomial(const std::vector<DATA_TYPE> & coefficients);			// linear polynomial with the given coefficients, the input matrix is a row vector
	Polynomial(const DATA_TYPE *pCoefficients, const unsigned int numVars);
	Polynomial(const Term<DATA_TYPE> & term);
	Polynomial(const std::list<Term<DATA_TYPE> > & term_list);
	Polynomial(const unsigned int varID, const unsigned int degree, const unsigned int numVars);

//	Polynomial(const UnivariatePolynomial & up, const int numVars);

	Polynomial(const Polynomial<DATA_TYPE> & polynomial);

	~Polynomial();

	void toHornerForm(HornerForm<DATA_TYPE> & hf) const;

	void reorder();														// sort the terms.
	void clear();

	void toReal(Polynomial<Real> & realPoly);
	Polynomial(const std::string & strPolynomial);

	void toString(std::string & result, const Variables & vars) const;

	void output(std::ostream & os, const Variables & vars) const;
	void output_constraint(std::ostream & os, const Variables & vars) const;

	void constant(DATA_TYPE & c) const;									// constant part of the polynomial

//	void toHornerForm(Expression_AST<DATA_TYPE> & hf) const;

	template <class DATA_TYPE2, class DATA_TYPE3>
	void evaluate(DATA_TYPE2 & result, const std::vector<DATA_TYPE3> & domain) const;

	template <class DATA_TYPE2>
	void evaluate_time(Polynomial<DATA_TYPE> & result, const std::vector<DATA_TYPE2> & step_exp_table) const;

	template <class DATA_TYPE2>
	void intEvalNormal(Interval & result, const std::vector<DATA_TYPE2> & step_exp_table) const;

	void pow(Polynomial<DATA_TYPE> & result, const unsigned int degree) const;
	void pow_assign(const unsigned int degree);
	void pow(Polynomial<DATA_TYPE> & result, const unsigned int degree, const unsigned int order) const;
	void pow_assign(const unsigned int degree, const unsigned int order);

//	void center();

	void mul_assign(const unsigned int varIndex, const unsigned int degree);		// multiplied by a term x^d
	void mul(Polynomial<DATA_TYPE> result, const unsigned int varIndex, const unsigned int degree) const;

	Polynomial<DATA_TYPE> & operator = (const Polynomial<DATA_TYPE> & polynomial);
	Polynomial<DATA_TYPE> & operator = (const Term<DATA_TYPE> & term);

	Polynomial<DATA_TYPE> & operator += (const Polynomial<DATA_TYPE> & polynomial);
	Polynomial<DATA_TYPE> & operator += (const Term<DATA_TYPE> & term);
	Polynomial<DATA_TYPE> & operator -= (const Polynomial<DATA_TYPE> & polynomial);
	Polynomial<DATA_TYPE> & operator -= (const Term<DATA_TYPE> & term);

	Polynomial<DATA_TYPE> & operator *= (const Polynomial<DATA_TYPE> & polynomial);
	Polynomial<DATA_TYPE> & operator *= (const Term<DATA_TYPE> & term);
	Polynomial<DATA_TYPE> & operator *= (const DATA_TYPE & c);

	Polynomial<DATA_TYPE> & operator /= (const DATA_TYPE & c);

	Polynomial<DATA_TYPE> operator + (const Polynomial<DATA_TYPE> & polynomial) const;
	Polynomial<DATA_TYPE> operator + (const Term<DATA_TYPE> & term) const;
	Polynomial<DATA_TYPE> operator - (const Polynomial<DATA_TYPE> & polynomial) const;
	Polynomial<DATA_TYPE> operator - (const Term<DATA_TYPE> & term) const;

	Polynomial<DATA_TYPE> operator * (const Polynomial<DATA_TYPE> & polynomial) const;
	Polynomial<DATA_TYPE> operator * (const Term<DATA_TYPE> & term) const;
	Polynomial<DATA_TYPE> operator * (const DATA_TYPE & c) const;

	Polynomial<DATA_TYPE> operator / (const DATA_TYPE & c) const;

	template <class DATA_TYPE2, class DATA_TYPE3>
	void ctrunc(DATA_TYPE2 & remainder, const std::vector<DATA_TYPE3> & domain, const unsigned int order);

	template <class DATA_TYPE2>
	void ctrunc_normal(Interval & remainder, const std::vector<DATA_TYPE2> & step_exp_table, const unsigned int order);

	void nctrunc(const unsigned int order);

	void linearCoefficients(Matrix<DATA_TYPE> & coefficients, const unsigned int row) const;
	void linearCoefficients(std::vector<DATA_TYPE> & coefficients) const;

	void rmConstant();							// remove the constant part

	void decompose(Polynomial<DATA_TYPE> & linear, Polynomial<DATA_TYPE> & other) const;

	unsigned int degree() const;				// degree of the polynomial
//	unsigned int degree_wo_t() const;			// degree of the polynomial without the time variable

//	bool isLinear_wo_t() const;
	bool isZero() const;

	void rmZeroTerms(const std::vector<unsigned int> & indices);

	void integral_time();

	template <class DATA_TYPE2>
	void cutoff_normal(Interval & intRem, const std::vector<DATA_TYPE2> & step_exp_table, const Interval & cutoff_threshold);

	template <class DATA_TYPE2>
	void cutoff(Interval & intRem, const std::vector<DATA_TYPE2> & domain, const Interval & cutoff_threshold);

	void cutoff(const Interval & cutoff_threshold);

	void derivative(Polynomial<DATA_TYPE> & result, const unsigned int varIndex) const;						// derivative with respect to a variable
	void LieDerivative(Polynomial<DATA_TYPE> & result, const std::vector<Polynomial<DATA_TYPE> > & f) const;	// Lie derivative without truncation

//	void insert(TaylorModel & result, const TaylorModelVec & vars, const std::vector<Interval> & varsPolyRange, const std::vector<Interval> & domain, const Interval & cutoff_threshold) const;
//	void insert_normal(TaylorModel & result, const TaylorModelVec & vars, const std::vector<Interval> & varsPolyRange, const std::vector<Interval> & step_exp_table, const int numVars, const Interval & cutoff_threshold) const;

//	void sub(Polynomial & result, const Polynomial & P, const int order) const;		// compute the subtraction of the monomials with some order

	void exp_taylor(Polynomial<DATA_TYPE> & result, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold) const;
	void rec_taylor(Polynomial<DATA_TYPE> & result, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold) const;
	void sin_taylor(Polynomial<DATA_TYPE> & result, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold) const;
	void cos_taylor(Polynomial<DATA_TYPE> & result, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold) const;
	void log_taylor(Polynomial<DATA_TYPE> & result, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold) const;
	void sqrt_taylor(Polynomial<DATA_TYPE> & result, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold) const;

/*
	void substitute(const int varID, const Interval & intVal);									// substitute a variable by an interval
	void substitute(const std::vector<int> & varIDs, const std::vector<Interval> & intVals);	// substitute a set of variables by intervals

	void substitute_with_precond(Interval & intRem, const std::vector<bool> & substitution, const std::vector<Interval> & step_exp_table);
	void substitute_with_precond_no_remainder(const std::vector<bool> & substitution);

	void simplification_in_decomposition(const std::vector<bool> & substitution);

	void substitute(Polynomial & result, const int varID, const Interval & intVal) const;
	void substitute(Polynomial & result, const std::vector<int> & varIDs, const std::vector<Interval> & intVals) const;
*/


//	void extend(const int num);		// current dim -> dim + num
//	void extend();					// current dim -> dim + 1


	template <class DATA_TYPE2>
	friend class Polynomial;

	template <class DATA_TYPE2>
	friend class TaylorModel;

	template <class DATA_TYPE2>
	friend class TaylorModelVec;

//	friend class Flowpipe;
//	friend class ContinuousSystem;
};


template <class DATA_TYPE>
Polynomial<DATA_TYPE>::Polynomial()
{
}

template <class DATA_TYPE>
Polynomial<DATA_TYPE>::Polynomial(const DATA_TYPE & c, const unsigned int numVars)
{
	if(c != 0)
	{
		Term<DATA_TYPE> term(c, numVars);
		terms.push_back(term);
	}
}

template <class DATA_TYPE>
Polynomial<DATA_TYPE>::Polynomial(Matrix<DATA_TYPE> & coefficients)
{
	int numVars = coefficients.cols();

	for(int i=numVars-1; i>=0; --i)
	{
		DATA_TYPE tmp = coefficients[0][i];
		if(tmp != 0)
		{
			Term<DATA_TYPE> term(tmp, numVars);
			term.degrees[i] = 1;
			term.d = 1;
			terms.push_back(term);
		}
	}
}

template <class DATA_TYPE>
Polynomial<DATA_TYPE>::Polynomial(const std::vector<DATA_TYPE> & coefficients)
{
	int numVars = coefficients.size();

	for(int i=numVars-1; i>=0; --i)
	{
		DATA_TYPE tmp = coefficients[i];
		if(tmp != 0)
		{
			Term<DATA_TYPE> term(tmp, numVars);
			term.degrees[i] = 1;
			term.d = 1;
			terms.push_back(term);
		}
	}
}

template <class DATA_TYPE>
Polynomial<DATA_TYPE>::Polynomial(const DATA_TYPE *pCoefficients, const unsigned int numVars)
{
	for(int i=numVars-1; i>=0; --i)
	{
		DATA_TYPE tmp = *(pCoefficients + i);
		if(tmp != 0)
		{
			Term<DATA_TYPE> term(tmp, numVars);
			term.degrees[i] = 1;
			term.d = 1;
			terms.push_back(term);
		}
	}
}

template <class DATA_TYPE>
Polynomial<DATA_TYPE>::Polynomial(const Term<DATA_TYPE> & term)
{
	terms.push_back(term);
}

template <class DATA_TYPE>
Polynomial<DATA_TYPE>::Polynomial(const std::list<Term<DATA_TYPE> > & term_list)
{
	terms = term_list;
}

template <class DATA_TYPE>
Polynomial<DATA_TYPE>::Polynomial(const unsigned int varID, const unsigned int degree, const unsigned int numVars)
{
	Term<DATA_TYPE> term(1, numVars);
	term.degrees[varID] = degree;
	term.d = degree;
	terms.push_back(term);
}

template <class DATA_TYPE>
Polynomial<DATA_TYPE>::Polynomial(const Polynomial<DATA_TYPE> & polynomial)
{
	terms = polynomial.terms;
}

template <class DATA_TYPE>
Polynomial<DATA_TYPE>::~Polynomial()
{
}

template <class DATA_TYPE>
void Polynomial<DATA_TYPE>::toHornerForm(HornerForm<DATA_TYPE> & hf) const
{
	hf.clear();

	if(terms.size() == 0)
		return;

	unsigned int numVars = (terms.begin())->degrees.size();

	std::list<Term<DATA_TYPE> > term_list = terms;
	typename std::list<Term<DATA_TYPE> >::iterator iter = term_list.begin();

	if(iter->d == 0)
	{
		hf.constant_part = iter->coefficient;
		iter = term_list.erase(iter);

		if(term_list.size() == 0)
			return;
	}

	std::vector<std::list<Term<DATA_TYPE> > > term_list_table(numVars);

	for(iter = term_list.begin(); iter != term_list.end(); ++iter)
	{
		unsigned int index = 0;

		for(unsigned int i=0; i<numVars; ++i)
		{
			if(iter->degrees[i] > 0)
			{
				index = i;
				break;
			}
		}

		iter->degrees[index] -= 1;
		iter->d -= 1;
		term_list_table[index].push_back(*iter);
	}

	for(unsigned int i=0; i<numVars; ++i)
	{
		Polynomial<DATA_TYPE> polyTmp(term_list_table[i]);
		HornerForm<DATA_TYPE> hfTmp;
		polyTmp.toHornerForm(hfTmp);
		hf.nonconstant_part.push_back(hfTmp);
	}
}

template <class DATA_TYPE>
void Polynomial<DATA_TYPE>::reorder()
{
	terms.sort();
}

template <class DATA_TYPE>
void Polynomial<DATA_TYPE>::clear()
{
	terms.clear();
}

template <>
inline void Polynomial<Interval>::toReal(Polynomial<Real> & realPoly)
{
	realPoly.terms.clear();

	std::list<Term<Interval> >::const_iterator iter;
	for(iter = terms.begin(); iter != terms.end(); ++iter)
	{
		Term<Real> term;
		term.coefficient = iter->coefficient.toReal();
		term.degrees = iter->degrees;
		term.d = iter->d;

		realPoly.terms.push_back(term);
	}
}

template <>
inline Polynomial<Interval>::Polynomial(const std::string & strPolynomial)
{
	multivariate_polynomial_setting.clear();

	std::string prefix(str_prefix_multivariate_polynomial);
	std::string suffix(str_suffix);

	multivariate_polynomial_setting.strPolynomial = prefix + strPolynomial + suffix;

	parseMultivariatePolynomial();

	*this = multivariate_polynomial_setting.result;
}

template <>
inline Polynomial<Real>::Polynomial(const std::string & strPolynomial)
{
	multivariate_polynomial_setting.clear();

	std::string prefix(str_prefix_multivariate_polynomial);
	std::string suffix(str_suffix);

	multivariate_polynomial_setting.strPolynomial = prefix + strPolynomial + suffix;

	parseMultivariatePolynomial();

	multivariate_polynomial_setting.result.toReal(*this);
}

template <class DATA_TYPE>
void Polynomial<DATA_TYPE>::toString(std::string & result, const Variables & vars) const
{
	if(terms.size() == 0)
	{
		result = "0";
		return;
	}
	else
	{
		result.clear();
		typename std::list<Term<DATA_TYPE> >::const_iterator iter = terms.begin(), iter_last = terms.end();

		--iter_last;

		for(; iter != iter_last; ++iter)
		{
			std::string tmp;
			iter->toString(tmp, vars.varNames);
			result += tmp + " + ";
		}

		std::string tmp;
		iter_last->toString(tmp, vars.varNames);
		result += tmp;
	}
}

template <class DATA_TYPE>
void Polynomial<DATA_TYPE>::output(std::ostream & os, const Variables & vars) const
{
	if(terms.size() == 0)
	{
		os << 0;
		return;
	}

	typename std::list<Term<DATA_TYPE> >::const_iterator iter = terms.begin(), iter_last = terms.end();

	--iter_last;

	for(; iter != iter_last; ++iter)
	{
		iter->output(os, vars);
		os << " + ";
	}

	iter_last->output(os, vars);
}

template <class DATA_TYPE>
void Polynomial<DATA_TYPE>::output_constraint(std::ostream & os, const Variables & vars) const
{
	if(terms.size() == 0)
	{
		os << 0;
		return;
	}

	typename std::list<Term<DATA_TYPE> >::const_iterator iter = terms.begin(), iter_last = terms.end();

	--iter_last;

	for(; iter != iter_last; ++iter)
	{
		iter->output_constraint(os, vars);
		os << " + ";
	}

	iter_last->output_constraint(os, vars);
}

template <class DATA_TYPE>
void Polynomial<DATA_TYPE>::constant(DATA_TYPE & c) const
{
	if(terms.size() > 0 && (terms.begin())->d == 0)
	{
		c = (terms.begin())->coefficient;
	}
	else
	{
		c = 0;
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2, class DATA_TYPE3>
void Polynomial<DATA_TYPE>::evaluate(DATA_TYPE2 & result, const std::vector<DATA_TYPE3> & domain) const
{
	HornerForm<DATA_TYPE> hf;
	toHornerForm(hf);
	hf.evaluate(result, domain);
}


template <class DATA_TYPE>
template <class DATA_TYPE2>
void Polynomial<DATA_TYPE>::evaluate_time(Polynomial<DATA_TYPE> & result, const std::vector<DATA_TYPE2> & step_exp_table) const
{
	result.clear();

	if(terms.size() == 0)
		return;

	typename std::list<Term<DATA_TYPE> >::const_iterator iter;

	if(step_exp_table[1] == 0 || step_exp_table.size() == 0)		// t = 0
	{
		for(iter = terms.begin(); iter != terms.end(); ++iter)
		{
			if(iter->degrees[0] == 0)
			{
				result += *iter;
			}
		}
	}
	else
	{
		for(iter = terms.begin(); iter != terms.end(); ++iter)
		{
			Term<DATA_TYPE> term = *iter;
			unsigned int tmp = term.degrees[0];

			if(tmp > 0)
			{
				term.coefficient *= step_exp_table[tmp];
				term.d -= tmp;
				term.degrees[0] = 0;
			}

			result += term;
		}
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void Polynomial<DATA_TYPE>::intEvalNormal(Interval & result, const std::vector<DATA_TYPE2> & step_exp_table) const
{
	result = 0;

	typename std::list<Term<DATA_TYPE> >::const_iterator iter;

	for(iter = terms.begin(); iter != terms.end(); ++iter)
	{
		Interval intTemp;
		iter->intEvalNormal(intTemp, step_exp_table);

		result += intTemp;
	}
}

template <class DATA_TYPE>
void Polynomial<DATA_TYPE>::pow(Polynomial<DATA_TYPE> & result, const unsigned int degree) const
{
	Polynomial<DATA_TYPE> tmp = *this;
	result = *this;

	for(int d = degree - 1; d > 0;)
	{
		if(d & 1)
		{
			result *= tmp;
		}

		d >>= 1;

		if(d > 0)
		{
			tmp *= tmp;
		}
	}
}

template <class DATA_TYPE>
void Polynomial<DATA_TYPE>::pow_assign(const unsigned int degree)
{
	Polynomial<DATA_TYPE> tmp = *this;
	Polynomial<DATA_TYPE> result = *this;

	for(int d = degree - 1; d > 0;)
	{
		if(d & 1)
		{
			result *= tmp;
		}

		d >>= 1;

		if(d > 0)
		{
			tmp *= tmp;
		}
	}

	*this = result;
}

template <class DATA_TYPE>
void Polynomial<DATA_TYPE>::pow(Polynomial<DATA_TYPE> & result, const unsigned int degree, const unsigned int order) const
{
	Polynomial<DATA_TYPE> p = *this;
	p.nctrunc(order);

	Polynomial<DATA_TYPE> tmp = p;
	result = p;

	for(int d = degree - 1; d > 0;)
	{
		if(d & 1)
		{
			result *= tmp;
			result.nctrunc(order);
		}

		d >>= 1;

		if(d > 0)
		{
			tmp *= tmp;
			tmp.nctrunc(order);
		}
	}
}

template <class DATA_TYPE>
void Polynomial<DATA_TYPE>::pow_assign(const unsigned int degree, const unsigned int order)
{
	Polynomial<DATA_TYPE> p = *this;
	p.nctrunc(order);

	Polynomial<DATA_TYPE> tmp = p;
	Polynomial<DATA_TYPE> result = p;

	for(int d = degree - 1; d > 0;)
	{
		if(d & 1)
		{
			result *= tmp;
			result.nctrunc(order);
		}

		d >>= 1;

		if(d > 0)
		{
			tmp *= tmp;
			tmp.nctrunc(order);
		}
	}

	*this = result;
}

template <class DATA_TYPE>
void Polynomial<DATA_TYPE>::mul_assign(const unsigned int varIndex, const unsigned int degree)
{
	typename std::list<Term<DATA_TYPE> >::iterator iter;

	for(iter = terms.begin(); iter != terms.end(); ++iter)
	{
		iter->degrees[varIndex] += degree;
		iter->d += degree;
	}
}

template <class DATA_TYPE>
void Polynomial<DATA_TYPE>::mul(Polynomial<DATA_TYPE> result, const unsigned int varIndex, const unsigned int degree) const
{
	result = *this;
	result.mul_assign(varIndex, degree);
}

template <class DATA_TYPE>
Polynomial<DATA_TYPE> & Polynomial<DATA_TYPE>::operator = (const Polynomial<DATA_TYPE> & polynomial)
{
	if(this == &polynomial)
		return *this;

	terms = polynomial.terms;
	return *this;
}

template <class DATA_TYPE>
Polynomial<DATA_TYPE> & Polynomial<DATA_TYPE>::operator = (const Term<DATA_TYPE> & term)
{
	terms.clear();
	terms.push_back(term);

	return *this;
}

template <class DATA_TYPE>
Polynomial<DATA_TYPE> & Polynomial<DATA_TYPE>::operator += (const Polynomial<DATA_TYPE> & polynomial)
{
	Polynomial<DATA_TYPE> result;

	typename std::list<Term<DATA_TYPE> >::const_iterator iterA;			// polynomial A
	typename std::list<Term<DATA_TYPE> >::const_iterator iterB;			// polynomial B

	for(iterA = terms.begin(), iterB = polynomial.terms.begin(); ; )
	{
		if(iterA == terms.end() || iterB == polynomial.terms.end())
			break;

		if((*iterA) < (*iterB))
	    {
			result.terms.push_back(*iterA);
			++iterA;
	    }
		else if((*iterB) < (*iterA))
	    {
			result.terms.push_back(*iterB);
			++iterB;
	    }
		else
		{
			DATA_TYPE tmp;
			tmp = iterA->coefficient + iterB->coefficient;

			if(tmp != 0)
			{
				Term<DATA_TYPE> term(*iterA);
				term.coefficient = tmp;
				result.terms.push_back(term);
			}

			++iterA;
			++iterB;
		}
	}

	if(iterA == terms.end() && iterB != polynomial.terms.end())
	{
		for(; iterB != polynomial.terms.end(); ++iterB)
			result.terms.push_back(*iterB);
	}
	else if(iterA != terms.end() && iterB == polynomial.terms.end())
	{
		for(; iterA != terms.end(); ++iterA)
			result.terms.push_back(*iterA);
	}

	*this = result;
	return *this;
}

template <class DATA_TYPE>
Polynomial<DATA_TYPE> & Polynomial<DATA_TYPE>::operator += (const Term<DATA_TYPE> & term)
{
	if(term.coefficient == 0)
	{
		return *this;
	}

	typename std::list<Term<DATA_TYPE> >::iterator iter;
	bool bAdded = false;

	for(iter = terms.begin(); iter != terms.end(); ++iter)
	{
		if(term < *iter)
		{
			terms.insert(iter, term);
			bAdded = true;
			break;
		}
		else if(term == *iter)
		{
			(*iter) += term;

			if(iter->coefficient == 0)
			{
				iter = terms.erase(iter);
			}

			bAdded = true;
			break;
		}
	}

	if(!bAdded)
	{
		terms.push_back(term);
	}

	return *this;
}

template <class DATA_TYPE>
Polynomial<DATA_TYPE> & Polynomial<DATA_TYPE>::operator -= (const Polynomial<DATA_TYPE> & polynomial)
{
	Polynomial<DATA_TYPE> result;

	typename std::list<Term<DATA_TYPE> >::const_iterator iterA;			// polynomial A
	typename std::list<Term<DATA_TYPE> >::const_iterator iterB;			// polynomial B

	for(iterA = terms.begin(), iterB = polynomial.terms.begin(); ; )
	{
		if(iterA == terms.end() || iterB == polynomial.terms.end())
			break;

		if((*iterA) < (*iterB))
	    {
			result.terms.push_back(*iterA);
			++iterA;
	    }
		else if((*iterB) < (*iterA))
	    {
			result.terms.push_back(-(*iterB));
			++iterB;
	    }
		else
		{
			DATA_TYPE tmp;
			tmp = iterA->coefficient - iterB->coefficient;

			if(tmp != 0)
			{
				Term<DATA_TYPE> term(*iterA);
				term.coefficient = tmp;
				result.terms.push_back(term);
			}

			++iterA;
			++iterB;
		}
	}

	if(iterA == terms.end() && iterB != polynomial.terms.end())
	{
		for(; iterB != polynomial.terms.end(); ++iterB)
			result.terms.push_back(-(*iterB));
	}
	else if(iterA != terms.end() && iterB == polynomial.terms.end())
	{
		for(; iterA != terms.end(); ++iterA)
			result.terms.push_back(*iterA);
	}

	*this = result;
	return *this;
}

template <class DATA_TYPE>
Polynomial<DATA_TYPE> & Polynomial<DATA_TYPE>::operator -= (const Term<DATA_TYPE> & term)
{
	if(term.coefficient == 0)
	{
		return *this;
	}

	typename std::list<Term<DATA_TYPE> >::iterator iter;
	bool bAdded = false;

	for(iter = terms.begin(); iter != terms.end(); ++iter)
	{
		if(term < *iter)
		{
			terms.insert(iter, -term);
			bAdded = true;
			break;
		}
		else if(term == *iter)
		{
			(*iter) -= term;

			if(iter->coefficient == 0)
			{
				iter = terms.erase(iter);
			}

			bAdded = true;
			break;
		}
	}

	if(!bAdded)
	{
		terms.push_back(-term);
	}

	return *this;
}

template <class DATA_TYPE>
Polynomial<DATA_TYPE> & Polynomial<DATA_TYPE>::operator *= (const Polynomial<DATA_TYPE> & polynomial)
{
	Polynomial<DATA_TYPE> result;

	if((terms.size() == 0) || (polynomial.terms.size() == 0))
	{
		this->clear();
		return *this;
	}

	typename std::list<Term<DATA_TYPE> >::const_iterator iterB;

	for(iterB = polynomial.terms.begin(); iterB != polynomial.terms.end(); ++iterB)
	{
		Polynomial<DATA_TYPE> tmp = *this;
		tmp *= *iterB;
		result += tmp;
	}

	*this = result;
	return *this;
}

template <class DATA_TYPE>
Polynomial<DATA_TYPE> & Polynomial<DATA_TYPE>::operator *= (const Term<DATA_TYPE> & term)
{
	if(term.coefficient == 0)
	{
		clear();
	}
	else
	{
		typename std::list<Term<DATA_TYPE> >::iterator iter;
		for(iter = terms.begin(); iter != terms.end(); ++iter)
		{
			(*iter) *= term;
		}
	}

	return *this;
}

template <class DATA_TYPE>
Polynomial<DATA_TYPE> & Polynomial<DATA_TYPE>::operator *= (const DATA_TYPE & c)
{
	if(c == 0)
	{
		clear();
	}
	else
	{
		typename std::list<Term<DATA_TYPE> >::iterator iter;
		for(iter = terms.begin(); iter != terms.end(); ++iter)
		{
			iter->coefficient *= c;
		}
	}

	return *this;
}

template <class DATA_TYPE>
Polynomial<DATA_TYPE> & Polynomial<DATA_TYPE>::operator /= (const DATA_TYPE & c)
{
	typename std::list<Term<DATA_TYPE> >::iterator iter;
	for(iter = terms.begin(); iter != terms.end(); ++iter)
	{
		iter->coefficient /= c;
	}

	return *this;
}

template <class DATA_TYPE>
Polynomial<DATA_TYPE> Polynomial<DATA_TYPE>::operator + (const Polynomial<DATA_TYPE> & polynomial) const
{
	Polynomial<DATA_TYPE> result;

	typename std::list<Term<DATA_TYPE> >::const_iterator iterA;			// polynomial A
	typename std::list<Term<DATA_TYPE> >::const_iterator iterB;			// polynomial B

	for(iterA = terms.begin(), iterB = polynomial.terms.begin(); ; )
	{
		if(iterA == terms.end() || iterB == polynomial.terms.end())
			break;

		if((*iterA) < (*iterB))
	    {
			result.terms.push_back(*iterA);
			++iterA;
	    }
		else if((*iterB) < (*iterA))
	    {
			result.terms.push_back(*iterB);
			++iterB;
	    }
		else
		{
			DATA_TYPE tmp;
			tmp = iterA->coefficient + iterB->coefficient;

			if(tmp != 0)
			{
				Term<DATA_TYPE> term(*iterA);
				term.coefficient = tmp;
				result.terms.push_back(term);
			}

			++iterA;
			++iterB;
		}
	}

	if(iterA == terms.end() && iterB != polynomial.terms.end())
	{
		for(; iterB != polynomial.terms.end(); ++iterB)
			result.terms.push_back(*iterB);
	}
	else if(iterA != terms.end() && iterB == polynomial.terms.end())
	{
		for(; iterA != terms.end(); ++iterA)
			result.terms.push_back(*iterA);
	}

	return result;
}

template <class DATA_TYPE>
Polynomial<DATA_TYPE> Polynomial<DATA_TYPE>::operator + (const Term<DATA_TYPE> & term) const
{
	Polynomial<DATA_TYPE> result = *this;
	result += term;
	return result;
}

template <class DATA_TYPE>
Polynomial<DATA_TYPE> Polynomial<DATA_TYPE>::operator - (const Polynomial<DATA_TYPE> & polynomial) const
{
	Polynomial<DATA_TYPE> result;

	typename std::list<Term<DATA_TYPE> >::const_iterator iterA;			// polynomial A
	typename std::list<Term<DATA_TYPE> >::const_iterator iterB;			// polynomial B

	for(iterA = terms.begin(), iterB = polynomial.terms.begin(); ; )
	{
		if(iterA == terms.end() || iterB == polynomial.terms.end())
			break;

		if((*iterA) < (*iterB))
	    {
			result.terms.push_back(*iterA);
			++iterA;
	    }
		else if((*iterB) < (*iterA))
	    {
			result.terms.push_back(-(*iterB));
			++iterB;
	    }
		else
		{
			DATA_TYPE tmp;
			tmp = iterA->coefficient - iterB->coefficient;

			if(tmp != 0)
			{
				Term<DATA_TYPE> term(*iterA);
				term.coefficient = tmp;
				result.terms.push_back(term);
			}

			++iterA;
			++iterB;
		}
	}

	if(iterA == terms.end() && iterB != polynomial.terms.end())
	{
		for(; iterB != polynomial.terms.end(); ++iterB)
			result.terms.push_back(-(*iterB));
	}
	else if(iterA != terms.end() && iterB == polynomial.terms.end())
	{
		for(; iterA != terms.end(); ++iterA)
			result.terms.push_back(*iterA);
	}

	return result;
}

template <class DATA_TYPE>
Polynomial<DATA_TYPE> Polynomial<DATA_TYPE>::operator - (const Term<DATA_TYPE> & term) const
{
	Polynomial<DATA_TYPE> result = *this;
	result -= term;
	return result;
}

template <class DATA_TYPE>
Polynomial<DATA_TYPE> Polynomial<DATA_TYPE>::operator * (const Polynomial<DATA_TYPE> & polynomial) const
{
	Polynomial<DATA_TYPE> result;

	if((terms.size() == 0) || (polynomial.terms.size() == 0))
	{
		return result;
	}

	typename std::list<Term<DATA_TYPE> >::const_iterator iterB;

	for(iterB = polynomial.terms.begin(); iterB != polynomial.terms.end(); ++iterB)
	{
		Polynomial<DATA_TYPE> tmp = *this;
		tmp *= *iterB;
		result += tmp;
	}

	return result;
}

template <class DATA_TYPE>
Polynomial<DATA_TYPE> Polynomial<DATA_TYPE>::operator * (const Term<DATA_TYPE> & term) const
{
	Polynomial<DATA_TYPE> result;
	if(term.coefficient == 0)
	{
		return result;
	}
	else
	{
		result = *this;
		result *= term;
		return result;
	}
}

template <class DATA_TYPE>
Polynomial<DATA_TYPE> Polynomial<DATA_TYPE>::operator * (const DATA_TYPE & c) const
{
	Polynomial<DATA_TYPE> result;
	if(c == 0)
	{
		return result;
	}
	else
	{
		result = *this;
		result *= c;
		return result;
	}
}

template <class DATA_TYPE>
Polynomial<DATA_TYPE> Polynomial<DATA_TYPE>::operator / (const DATA_TYPE & c) const
{
	Polynomial<DATA_TYPE> result = *this;
	result /= c;
	return result;
}

template <class DATA_TYPE>
template <class DATA_TYPE2, class DATA_TYPE3>
void Polynomial<DATA_TYPE>::ctrunc(DATA_TYPE2 & remainder, const std::vector<DATA_TYPE3> & domain, const unsigned int order)
{
	Polynomial<DATA_TYPE> polyTmp;
	Term<DATA_TYPE> termTmp;

	for(; terms.size() > 0;)
	{
		termTmp = terms.back();

		if(termTmp.d > order)
		{
			polyTmp.terms.insert(polyTmp.terms.begin(), termTmp);
			terms.pop_back();
		}
		else
		{
			break;
		}
	}

	polyTmp.evaluate(remainder, domain);
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void Polynomial<DATA_TYPE>::ctrunc_normal(Interval & remainder, const std::vector<DATA_TYPE2> & step_exp_table, const unsigned int order)
{
	Polynomial<DATA_TYPE> polyTmp;
	Term<DATA_TYPE> termTmp;

	for(; terms.size() > 0;)
	{
		termTmp = terms.back();

		if(termTmp.d > order)
		{
			polyTmp.terms.insert(polyTmp.terms.begin(), termTmp);
			terms.pop_back();
		}
		else
		{
			break;
		}
	}

	polyTmp.intEvalNormal(remainder, step_exp_table);
}

template <class DATA_TYPE>
void Polynomial<DATA_TYPE>::nctrunc(const unsigned int order)
{
	for(; terms.size() > 0;)
	{
		if(terms.back().d > order)
		{
			terms.pop_back();
		}
		else
		{
			break;
		}
	}
}

template <class DATA_TYPE>
void Polynomial<DATA_TYPE>::linearCoefficients(Matrix<DATA_TYPE> & coefficients, const unsigned int row) const
{
	typename std::list<Term<DATA_TYPE> >::const_iterator iter;

	for(iter = terms.begin(); iter != terms.end(); ++iter)
	{
		unsigned int i;

		if(iter->d > 1)
			break;

		if(iter->isLinear(i))
		{
			if(i != 0)		// variable t is not considered
			{
				coefficients[row][i-1] = iter->coefficient;
			}
		}
	}
}

template <class DATA_TYPE>
void Polynomial<DATA_TYPE>::linearCoefficients(std::vector<DATA_TYPE> & coefficients) const
{
	typename std::list<Term<DATA_TYPE> >::const_iterator iter;

	for(iter = terms.begin(); iter != terms.end(); ++iter)
	{
		unsigned int i;

		if(iter->d > 1)
			break;

		if(iter->isLinear(i))
		{
			if(i != 0)		// variable t is not considered
			{
				coefficients[i-1] = iter->coefficient;
			}
		}
	}
}

template <class DATA_TYPE>
void Polynomial<DATA_TYPE>::rmConstant()
{
	if(terms.size() > 0 && terms.front().d == 0)
	{
		terms.erase( terms.begin() );
	}
}

template <class DATA_TYPE>
void Polynomial<DATA_TYPE>::decompose(Polynomial<DATA_TYPE> & linear, Polynomial<DATA_TYPE> & other) const
{
	typename std::list<Term<DATA_TYPE> >::const_iterator iter;

	linear.terms.clear();
	other.terms.clear();

	for(iter=terms.begin(); iter!=terms.end(); ++iter)
	{
		if(iter->d != 1)
		{
			other.terms.push_back(*iter);
		}
		else
		{
			linear.terms.push_back(*iter);
		}
	}
}

template <class DATA_TYPE>
unsigned int Polynomial<DATA_TYPE>::degree() const
{
	if(terms.size() > 0)
	{
		return terms.back().d;
	}
	else
	{
		return 0;
	}
}

template <class DATA_TYPE>
bool Polynomial<DATA_TYPE>::isZero() const
{
	if(terms.size() == 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

template <class DATA_TYPE>
void Polynomial<DATA_TYPE>::rmZeroTerms(const std::vector<unsigned int> & indices)
{
	if(indices.size() == 0)
	{
		return;
	}

	typename std::list<Term<DATA_TYPE> >::iterator iter = terms.begin();

	for(; iter != terms.end();)
	{
		bool bDeleted = false;

		for(unsigned int i=0; i<indices.size(); ++i)
		{
			if(iter->degrees[indices[i]] > 0)
			{
				iter = terms.erase(iter);
				bDeleted = true;
				break;
			}
		}

		if(bDeleted == false)
		{
			++iter;
		}
	}
}

template <class DATA_TYPE>
void Polynomial<DATA_TYPE>::integral_time()
{
	typename std::list<Term<DATA_TYPE> >::iterator iter = terms.begin();

	for(; iter != terms.end(); ++iter)
	{
		if(iter->degrees[0] > 0)
		{
			iter->degrees[0] += 1;
			iter->d += 1;
			iter->coefficient /= iter->degrees[0];
		}
		else
		{
			iter->degrees[0] += 1;
			iter->d += 1;
		}
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void Polynomial<DATA_TYPE>::cutoff_normal(Interval & intRem, const std::vector<DATA_TYPE2> & step_exp_table, const Interval & cutoff_threshold)
{
	Polynomial<DATA_TYPE> polyTmp;

	typename std::list<Term<DATA_TYPE> >::iterator iter;
	for(iter = terms.begin(); iter != terms.end(); )
	{
		if(iter->coefficient.belongsTo(cutoff_threshold))
		{
			polyTmp.terms.push_back(*iter);
			iter = terms.erase(iter);
		}
		else
		{
			++iter;
		}
	}

	polyTmp.intEvalNormal(intRem, step_exp_table);
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void Polynomial<DATA_TYPE>::cutoff(Interval & intRem, const std::vector<DATA_TYPE2> & domain, const Interval & cutoff_threshold)
{
	Polynomial<DATA_TYPE> polyTmp;

	typename std::list<Term<DATA_TYPE> >::iterator iter;
	for(iter = terms.begin(); iter != terms.end(); )
	{
		if(iter->coefficient.belongsTo(cutoff_threshold))
		{
			polyTmp.terms.push_back(*iter);
			iter = terms.erase(iter);
		}
		else
		{
			++iter;
		}
	}

	polyTmp.evaluate(intRem, domain);
}

template <class DATA_TYPE>
void Polynomial<DATA_TYPE>::cutoff(const Interval & cutoff_threshold)
{
	typename std::list<Term<DATA_TYPE> >::iterator iter;
	for(iter = terms.begin(); iter != terms.end(); )
	{
		if(iter->coefficient.belongsTo(cutoff_threshold))
		{
			iter = terms.erase(iter);
		}
		else
		{
			++iter;
		}
	}
}

template <class DATA_TYPE>
void Polynomial<DATA_TYPE>::derivative(Polynomial<DATA_TYPE> & result, const unsigned int varIndex) const
{
	result = *this;

	typename std::list<Term<DATA_TYPE> >::iterator iter;

	for(iter = result.terms.begin(); iter != result.terms.end(); )
	{
		if(iter->degrees[varIndex] > 0)
		{
			iter->coefficient *= iter->degrees[varIndex];
			iter->degrees[varIndex] -= 1;
			iter->d -= 1;
			++iter;
		}
		else
		{
			iter = result.terms.erase(iter);
		}
	}
}

template <class DATA_TYPE>
void Polynomial<DATA_TYPE>::LieDerivative(Polynomial<DATA_TYPE> & result, const std::vector<Polynomial<DATA_TYPE> > & f) const
{
	derivative(result, 0);

	unsigned int n = f.size();

	for(unsigned int i=0; i<n; ++i)
	{
		Polynomial<DATA_TYPE> p;
		derivative(p, i+1);
		p *= f[i];
		result += p;
	}
}

template <class DATA_TYPE>
void Polynomial<DATA_TYPE>::exp_taylor(Polynomial<DATA_TYPE> & result, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold) const
{
	DATA_TYPE const_part;

	Polynomial<DATA_TYPE> F = *this;

	// remove the center point of tm
	F.constant(const_part);
	F.rmConstant();				// F = tm - c

	const_part.exp_assign();	// exp(c)

	if(F.isZero())				// tm = c
	{
		Polynomial<DATA_TYPE> polyExp(const_part, numVars);
		result = polyExp;
		return;
	}

	Polynomial<DATA_TYPE> polyOne(1, numVars);

	// to compute the expression 1 + F + (1/2!)F^2 + ... + (1/k!)F^k,
	// we evaluate its Horner form (...((1/(k-1))((1/k)*F+1)*F + 1) ... + 1)

	result = polyOne;

	for(int i=order; i>0; --i)
	{
		result /= i;
		result *= F;

		result.nctrunc(order);
		result.cutoff(cutoff_threshold);

		result += polyOne;
	}

	result *= const_part;
}

template <class DATA_TYPE>
void Polynomial<DATA_TYPE>::rec_taylor(Polynomial<DATA_TYPE> & result, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold) const
{
	DATA_TYPE const_part;

	Polynomial F = *this;

	// remove the center point of tm
	F.constant(const_part);
	F.rmConstant();				// F = tm - c

	const_part.rec_assign();	// 1/c

	if(F.isZero())				// tm = c
	{
		Polynomial<DATA_TYPE> polyExp(const_part, numVars);
		result = polyExp;
		return;
	}

	Polynomial<DATA_TYPE> polyOne(1, numVars);
	Polynomial<DATA_TYPE> F_c = F * const_part;

	// to compute the expression 1 - F/c + (F/c)^2 - ... + (-1)^k (F/c)^k,
	// we evaluate its Horner form (-1)*(...((-1)*(-F/c + 1)*F/c + 1)...) + 1

	result = polyOne;

	for(int i=order; i>0; --i)
	{
		result *= -1;
		result *= F_c;

		result.nctrunc(order);
		result.cutoff(cutoff_threshold);

		result += polyOne;
	}

	result *= const_part;
}

template <class DATA_TYPE>
void Polynomial<DATA_TYPE>::sin_taylor(Polynomial<DATA_TYPE> & result, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold) const
{
	DATA_TYPE const_part;

	Polynomial<DATA_TYPE> F = *this;

	// remove the center point of tm
	F.constant(const_part);
	F.rmConstant();			// F = tm - c

	if(F.isZero())			// tm = c
	{
		const_part.sin_assign();
		Polynomial<DATA_TYPE> polyExp(const_part, numVars);
		result = polyExp;
		return;
	}

	DATA_TYPE sinc, cosc, msinc, mcosc;
	const_part.sin(sinc);
	const_part.cos(cosc);

	msinc = -sinc;
	mcosc = -cosc;

	Polynomial<DATA_TYPE> polyTmp(sinc, numVars);
	result = polyTmp;

	int k = 1;

	Polynomial<DATA_TYPE> polyPowerF(1, numVars);

	for(unsigned int i=1; i<=order; ++i, ++k)
	{
		k %= 4;

		switch(k)
		{
		case 0:
		{
			polyPowerF *= F;
			polyPowerF.nctrunc(order);
			polyPowerF.cutoff(cutoff_threshold);
			polyPowerF *= (sinc / i);

			result += polyPowerF;

			break;
		}
		case 1:
		{
			polyPowerF *= F;
			polyPowerF.nctrunc(order);
			polyPowerF.cutoff(cutoff_threshold);
			polyPowerF *= (cosc / i);

			result += polyPowerF;

			break;
		}
		case 2:
		{
			polyPowerF *= F;
			polyPowerF.nctrunc(order);
			polyPowerF.cutoff(cutoff_threshold);
			polyPowerF *= (msinc / i);

			result += polyPowerF;

			break;
		}
		case 3:
		{
			polyPowerF *= F;
			polyPowerF.nctrunc(order);
			polyPowerF.cutoff(cutoff_threshold);
			polyPowerF *= (mcosc / i);

			result += polyPowerF;

			break;
		}
		}
	}

	result.cutoff(cutoff_threshold);
}

template <class DATA_TYPE>
void Polynomial<DATA_TYPE>::cos_taylor(Polynomial<DATA_TYPE> & result, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold) const
{
	DATA_TYPE const_part;

	Polynomial<DATA_TYPE> F = *this;

	// remove the center point of tm
	F.constant(const_part);
	F.rmConstant();			// F = tm - c

	if(F.isZero())			// tm = c
	{
		const_part.cos_assign();
		Polynomial<DATA_TYPE> polyExp(const_part, numVars);
		result = polyExp;
		return;
	}

	DATA_TYPE sinc, cosc, msinc, mcosc;
	const_part.sin(sinc);
	const_part.cos(cosc);

	msinc = -sinc;
	mcosc = -cosc;

	Polynomial<DATA_TYPE> polyTemp(cosc, numVars);
	result = polyTemp;

	int k = 1;

	Polynomial polyPowerF(1, numVars);

	for(unsigned int i=1; i<=order; ++i, ++k)
	{
		k %= 4;

		switch(k)
		{
		case 0:
		{
			polyPowerF *= F;
			polyPowerF.nctrunc(order);
			polyPowerF.cutoff(cutoff_threshold);
			polyPowerF *= (cosc / i);

			result += polyPowerF;

			break;
		}
		case 1:
		{
			polyPowerF *= F;
			polyPowerF.nctrunc(order);
			polyPowerF.cutoff(cutoff_threshold);
			polyPowerF *= (msinc / i);

			result += polyPowerF;

			break;
		}
		case 2:
		{
			polyPowerF *= F;
			polyPowerF.nctrunc(order);
			polyPowerF.cutoff(cutoff_threshold);
			polyPowerF *= (mcosc / i);

			result += polyPowerF;

			break;
		}
		case 3:
		{
			polyPowerF *= F;
			polyPowerF.nctrunc(order);
			polyPowerF.cutoff(cutoff_threshold);
			polyPowerF *= (sinc / i);

			result += polyPowerF;

			break;
		}
		}
	}

	result.cutoff(cutoff_threshold);
}

template <class DATA_TYPE>
void Polynomial<DATA_TYPE>::log_taylor(Polynomial<DATA_TYPE> & result, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold) const
{
	DATA_TYPE const_part;

	Polynomial<DATA_TYPE> F = *this;

	// remove the center point of tm
	F.constant(const_part);
	F.rmConstant();				// F = tm - c

	DATA_TYPE C = const_part;

	const_part.log_assign();	// log(c)

	if(F.isZero())				// tm = c
	{
		Polynomial<DATA_TYPE> polyLog(const_part, numVars);
		result = polyLog;

		return;
	}

	Polynomial<DATA_TYPE> F_c = F / C;
	result = F_c / order;

	for(int i=order-1; i>=1; --i)
	{
		DATA_TYPE tmp = 1;
		tmp /= i;
		Polynomial<DATA_TYPE> polyTmp(tmp, numVars);

		result -= polyTmp;
		result *= -1;

		result *= F_c;
		result.nctrunc(order);
		result.cutoff(cutoff_threshold);
	}

	Polynomial<DATA_TYPE> const_part_poly(const_part, numVars);
	result += const_part_poly;
}

template <class DATA_TYPE>
void Polynomial<DATA_TYPE>::sqrt_taylor(Polynomial<DATA_TYPE> & result, const unsigned int numVars, const unsigned int order, const Interval & cutoff_threshold) const
{
	DATA_TYPE const_part;

	Polynomial<DATA_TYPE> F = *this;

	// remove the center point of tm
	F.constant(const_part);
	F.rmConstant();				// F = tm - c

	DATA_TYPE C = const_part;
	const_part.sqrt_assign();	// sqrt(c)

	if(F.isZero())				// tm = c
	{
		Polynomial<DATA_TYPE> polySqrt(const_part, numVars);
		result = polySqrt;

		return;
	}

	Polynomial<DATA_TYPE> F_2c = F / (2 * C);	// F/2c
	Polynomial<DATA_TYPE> polyOne(1, numVars);

	result = F_2c;

	for(int i=order, j=2*order-3; i>=2; --i, j-=2)
	{
		DATA_TYPE tmp = j;
		tmp /= -i;
		result *= tmp;

		result += polyOne;
		result *= F_2c;
		result.nctrunc(order);
		result.cutoff(cutoff_threshold);
	}

	result += polyOne;
	result *= const_part;
}


}

#endif /* POLYNOMIAL_H_ */
