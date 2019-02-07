/*---
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.
  The code is released as is under the GNU General Public License (GPL).
---*/

#ifndef TERM_H_
#define TERM_H_

#include "Interval.h"
#include "Variables.h"

namespace flowstar
{

template <class DATA_TYPE>
class Term
{
protected:
	DATA_TYPE coefficient;					// the coefficient of the term
	std::vector<unsigned int> degrees;		// the degrees of the variables, e.g., [2,0,4] is the notation for x1^2 x3^4
	unsigned int d;			        		// the degree of the term, it is the sum of the values in degrees.

public:
	Term();		// empty term.
	Term(const DATA_TYPE & c, const std::vector<unsigned int> & degs);
	Term(const Term<DATA_TYPE> & term);
	Term(const DATA_TYPE & c);				// a constant
	Term(const DATA_TYPE & c, const unsigned int numVars);
	~Term();

	unsigned int degree() const;		// degree of the term
	unsigned int dimension() const;		// dimension of the term

	void evaluate(DATA_TYPE & result, const std::vector<DATA_TYPE> & domain) const;	// evaluation of the term

	// interval evaluation of the term, we assume that the domain is normalized to [0,s] x [-1,1]^(d-1)
	template <class DATA_TYPE2>
	void intEvalNormal(Interval & result, const std::vector<DATA_TYPE2> & step_exp_table) const;

	Term<DATA_TYPE> & operator = (const Term<DATA_TYPE> & term);

	Term<DATA_TYPE> & operator += (const Term<DATA_TYPE> & term);			// we assume the two terms can be added up
	Term<DATA_TYPE> & operator -= (const Term<DATA_TYPE> & term);
	Term<DATA_TYPE> & operator *= (const Term<DATA_TYPE> & term);

	Term<DATA_TYPE> operator - () const;
	Term<DATA_TYPE> operator + (const Term<DATA_TYPE> & term) const;
	Term<DATA_TYPE> operator - (const Term<DATA_TYPE> & term) const;
	Term<DATA_TYPE> operator * (const Term<DATA_TYPE> & term) const;

	bool isLinear(unsigned int & index) const;							// Check if the degree of the term is 1. If so then return the index of the variable of degree 1.

	void toString(std::string & result, const Variables & vars) const;
	void output(std::ostream & os, const Variables & vars) const;
	void output_constraint(std::ostream & os, const Variables & vars) const;

	bool operator < (const Term<DATA_TYPE> & term) const;					// Define a partial order over the terms
	bool operator == (const Term<DATA_TYPE> & term) const;

	bool center();

	int cutoff(Term<DATA_TYPE> & remainder, const Interval & cutoff_threshold);
	int cutoff(const Interval & cutoff_threshold);

/*
	void substitute(const int varID, const DATA_TYPE & value);											// substitute a variable by an Interval
	void substitute(const std::vector<unsigned int> & varIDs, const std::vector<DATA_TYPE> & values);	// substitute a set of variables by intervals

	bool substitute_with_precond(const std::vector<bool> & substitution);

	void substitute(Term & result, const int varID, const Interval & intVal) const;
	void substitute(Term & result, const std::vector<int> & varIDs, const std::vector<Interval> & intVals) const;
*/

	void extend(const unsigned int num);
	void extend();

	template <class DATA_TYPE2>
	friend class Polynomial;

	template <class DATA_TYPE2>
	friend class TaylorModel;


//	friend class TaylorModelVec;
//	friend class UnivariatePolynomial;
//	friend class Matrix;
//	friend class Flowpipe;
};

template <class DATA_TYPE>
Term<DATA_TYPE>::Term()
{
	d = 0;
}

template <class DATA_TYPE>
Term<DATA_TYPE>::Term(const DATA_TYPE & c, const std::vector<unsigned int> & degs) : coefficient(c), degrees(degs)
{
	for(int i=0; i<degs.size(); ++i)
	{
		d += degs[i];
	}
}

template <class DATA_TYPE>
Term<DATA_TYPE>::Term(const Term & term) : coefficient(term.coefficient), degrees(term.degrees), d(term.d)
{
}

template <class DATA_TYPE>
Term<DATA_TYPE>::Term(const DATA_TYPE & c) : coefficient(c), d(0)
{
}

template <class DATA_TYPE>
Term<DATA_TYPE>::Term(const DATA_TYPE & c, const unsigned int numVars)
{
	coefficient = c;
	degrees.resize(numVars, 0);
	d = 0;
}

template <class DATA_TYPE>
Term<DATA_TYPE>::~Term()
{
}

template <class DATA_TYPE>
unsigned int Term<DATA_TYPE>::degree() const
{
	return d;
}

template <class DATA_TYPE>
unsigned int Term<DATA_TYPE>::dimension() const
{
	return degrees.size();
}

template <class DATA_TYPE>
void Term<DATA_TYPE>::evaluate(DATA_TYPE & result, const std::vector<DATA_TYPE> & domain) const
{
	result = coefficient;

	for(int i=0; i<degrees.size(); ++i)
	{
		DATA_TYPE tmp1 = domain[i];
		DATA_TYPE tmp2 = tmp1;

		for(int d = degrees[i] - 1; d > 0;)
		{
			if(d & 1)
			{
				tmp2 *= tmp1;
			}

			d >>= 1;

			if(d > 0)
			{
				tmp1 *= tmp1;
			}
		}

		result *= tmp2;
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void Term<DATA_TYPE>::intEvalNormal(Interval & result, const std::vector<DATA_TYPE2> & step_exp_table) const
{
	result = 0;

	if(degrees.size() == 0)
		return;

	result = coefficient;
	result *= step_exp_table[degrees[0]];

	Interval evenInt(0,1), oddInt(-1,1);
	Interval intFactor(1);
	bool bSet = false;

	for(unsigned int i=1; i<degrees.size(); ++i)
	{
		if(degrees[i] == 0)			// degree is zero
		{
			continue;
		}
		else if(degrees[i]%2 == 0)	// degree is an even number
		{
			if(!bSet)
			{
				intFactor = evenInt;
				bSet = true;
			}
		}
		else						// degree is an odd number
		{
			intFactor = oddInt;
			break;
		}
	}

	result *= intFactor;
}

template <class DATA_TYPE>
Term<DATA_TYPE> & Term<DATA_TYPE>::operator = (const Term<DATA_TYPE> & term)
{
	if(this == &term)
		return *this;

	coefficient = term.coefficient;
	degrees = term.degrees;
	d = term.d;

	return *this;
}

template <class DATA_TYPE>
Term<DATA_TYPE> & Term<DATA_TYPE>::operator += (const Term<DATA_TYPE> & term)
{
	coefficient += term.coefficient;
	return *this;
}

template <class DATA_TYPE>
Term<DATA_TYPE> & Term<DATA_TYPE>::operator -= (const Term<DATA_TYPE> & term)
{
	coefficient -= term.coefficient;
	return *this;
}

template <class DATA_TYPE>
Term<DATA_TYPE> & Term<DATA_TYPE>::operator *= (const Term<DATA_TYPE> & term)
{
	coefficient *= term.coefficient;

	for(unsigned int i=0; i<degrees.size(); ++i)
	{
		degrees[i] += term.degrees[i];
	}

	d += term.d;
	return *this;
}

template <class DATA_TYPE>
Term<DATA_TYPE> Term<DATA_TYPE>::operator - () const
{
	Term<DATA_TYPE> result;
	result.coefficient = -(this->coefficient);
	result.degrees = this->degrees;
	result.d = this->d;

	return result;
}

template <class DATA_TYPE>
Term<DATA_TYPE> Term<DATA_TYPE>::operator + (const Term<DATA_TYPE> & term) const
{
	Term<DATA_TYPE> result = *this;
	result += term;
	return result;
}

template <class DATA_TYPE>
Term<DATA_TYPE> Term<DATA_TYPE>::operator - (const Term<DATA_TYPE> & term) const
{
	Term<DATA_TYPE> result = *this;
	result -= term;
	return result;
}

template <class DATA_TYPE>
Term<DATA_TYPE> Term<DATA_TYPE>::operator * (const Term<DATA_TYPE> & term) const
{
	Term<DATA_TYPE> result = *this;
	result *= term;
	return result;
}

template <class DATA_TYPE>
bool Term<DATA_TYPE>::isLinear(unsigned int & index) const
{
	if(d == 1)
	{
		for(int i=0; i<degrees.size(); ++i)
		{
			if(degrees[i] == 1)
			{
				index = i;
				return true;
			}
		}
	}

	return false;
}

template <class DATA_TYPE>
void Term<DATA_TYPE>::toString(std::string & result, const Variables & vars) const
{
	result = '(' + coefficient.toString();

	for(int i=0; i<degrees.size(); i++)
	{
		if(degrees[i] != 0)
		{
			if(degrees[i] == 1)
			{
				result += " * " + vars.varNames[i];
			}
			else
			{
				result += " * " + vars.varNames[i] + "^" + std::to_string(degrees[i]);
			}
		}
	}

	result += ')';
}

template <class DATA_TYPE>
void Term<DATA_TYPE>::output(std::ostream & os, const Variables & vars) const
{
	std::string str;
	toString(str, vars);
	os << str;
}

template <class DATA_TYPE>
void Term<DATA_TYPE>::output_constraint(std::ostream & os, const Variables & vars) const
{
	os << "(" << coefficient << " * ";

	for(int i=1; i<degrees.size(); i++)
	{
		if(degrees[i] != 0)
		{
			if(degrees[i] == 1)
			{
				os << vars.varNames[i-1];
			}
			else
			{
				os << vars.varNames[i-1] << "^" << degrees[i];
			}
		}
	}

	os << ")";
}

template <class DATA_TYPE>
bool Term<DATA_TYPE>::operator == (const Term<DATA_TYPE> & term) const
{
	if (d == term.d)
	{
		for(int i=0; i<degrees.size(); i++)
		{
			if(degrees[i] != term.degrees[i])
				return false;
		}

		return true;	// The two terms are identical without considering the coefficients.
	}
	else
		return false;
}

template <class DATA_TYPE>
bool Term<DATA_TYPE>::operator < (const Term<DATA_TYPE> & term) const
{
	if(d < term.d)
	{
		return true;
	}
	else if(d > term.d)
	{
		return false;
	}
	else
	{
		for(int i=0; i<degrees.size(); ++i)
		{
			if(degrees[i] < term.degrees[i])
				return true;
			else if(degrees[i] > term.degrees[i])
				return false;
		}
	}

	return false;	// =
}

template <class DATA_TYPE>
int Term<DATA_TYPE>::cutoff(Term<DATA_TYPE> & remainder, const Interval & cutoff_threshold)
{
	if(coefficient.subseteq(cutoff_threshold))
	{
		return 2;
	}
	else
	{
		return 0;
	}
}

template <>
inline int Term<Interval>::cutoff(Term<Interval> & remainder, const Interval & cutoff_threshold)
{
	if(coefficient.width() >= MAX_WIDTH)
	{
		Interval midpoint;
		remainder = *this;
		remainder.coefficient.remove_midpoint(midpoint);
		coefficient = midpoint;
		return 1;
	}
	else if(coefficient.subseteq(cutoff_threshold))
	{
		return 2;
	}
	else
	{
		return 0;
	}
}

template <class DATA_TYPE>
int Term<DATA_TYPE>::cutoff(const Interval & cutoff_threshold)
{
	if(coefficient.subseteq(cutoff_threshold))
	{
		return 2;
	}
	else
	{
		return 0;
	}
}

template <>
inline int Term<Interval>::cutoff(const Interval & cutoff_threshold)
{
	if(coefficient.width() >= MAX_WIDTH)
	{
		coefficient = coefficient.midpoint();
		return 1;
	}
	else if(coefficient.subseteq(cutoff_threshold))
	{
		return 2;
	}
	else
	{
		return 0;
	}
}

template <class DATA_TYPE>
bool Term<DATA_TYPE>::center()
{
	return true;
}

template <>
inline bool Term<Interval>::center()
{
	Interval midpoint, intZero;

	coefficient.midpoint(midpoint);
	if(midpoint.subseteq(intZero))
	{
		return false;
	}
	else
	{
		coefficient = midpoint;
		return true;
	}
}

/*
void Term::substitute(const int varID, const Interval & intVal)
{
	int deg = degrees[varID];

	if(deg > 0)
	{
		d -= deg;
		degrees[varID] = 0;

		Interval pow = intVal.pow(deg);
		coefficient *= pow;
	}
}

void Term::substitute(const std::vector<int> & varIDs, const std::vector<Interval> & intVals)
{
	for(int i=0; i<varIDs.size(); ++i)
	{
		int deg = degrees[varIDs[i]];

		if(deg > 0)
		{
			d -= deg;
			degrees[varIDs[i]] = 0;

			Interval pow = intVals[i].pow(deg);
			coefficient *= pow;
		}
	}
}

bool Term::substitute_with_precond(const std::vector<bool> & substitution)
{
	int multi = 0;
	bool breformulate = false;

	for(int i=1; i<substitution.size(); ++i)
	{
		if(substitution[i])
		{
			int deg = degrees[i];

			if(deg > 1)
			{
				return 2;
			}

			if(deg >= 1)
			{
				if(multi <= 1)
				{
					++multi;
				}
				else
				{
					return true;
				}
			}
		}
	}

	return false;
}

void Term::substitute(Term & result, const int varID, const Interval & intVal) const
{
	result = *this;
	result.substitute(varID, intVal);
}

void Term::substitute(Term & result, const std::vector<int> & varIDs, const std::vector<Interval> & intVals) const
{
	result = *this;
	result.substitute(varIDs, intVals);
}
*/

template <class DATA_TYPE>
void Term<DATA_TYPE>::extend(const unsigned int num)
{
	for(unsigned int i=0; i<num; ++i)
	{
		degrees.push_back(0);
	}
}

template <class DATA_TYPE>
void Term<DATA_TYPE>::extend()
{
	degrees.insert(degrees.begin(), 0);
}

}

#endif /* TERM_H_ */
