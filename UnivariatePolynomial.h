/*---
  Flow*: A Verification Tool for Cyber-Physical Systems.
  Authors: Xin Chen, Sriram Sankaranarayanan, and Erika Abraham.
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.

  The code is released as is under the GNU General Public License (GPL).
---*/

#ifndef UNIVARIATEPOLYNOMIAL_H_
#define UNIVARIATEPOLYNOMIAL_H_

#include "Interval.h"
#include "settings.h"

//void parseUnivariatePolynomial(const std::string & strPolynomial);

namespace flowstar
{

template <class DATA_TYPE>
class UnivariatePolynomial;

template <class DATA_TYPE>
class UnivariateTaylorModel;

extern UnivariatePolynomial<Real> up_parseresult;

template <class DATA_TYPE>
std::ostream & operator << (std::ostream & os, const UnivariatePolynomial<DATA_TYPE> & up);


template <class DATA_TYPE>
class UnivariatePolynomial
{
protected:
	std::vector<DATA_TYPE> coefficients;

public:
	UnivariatePolynomial();
	UnivariatePolynomial(const std::vector<DATA_TYPE> & coeffs);
	UnivariatePolynomial(const UnivariatePolynomial<DATA_TYPE> & up);
	UnivariatePolynomial(const DATA_TYPE & c);
	UnivariatePolynomial(const double c);
	UnivariatePolynomial(const DATA_TYPE & c, const unsigned int degree);

	~UnivariatePolynomial();

//	bool set(const std::string & strPolynomial);
	void clear();
	void set2zero();
	unsigned int degree() const;
	bool isZero() const;

	template <class DATA_TYPE2, class DATA_TYPE3>
	void round(DATA_TYPE2 & remainder, const DATA_TYPE3 & val);

	template <class DATA_TYPE2, class DATA_TYPE3>
	void evaluate(DATA_TYPE2 & result, const std::vector<DATA_TYPE3> & val_exp_table) const;

	template <class DATA_TYPE2, class DATA_TYPE3>
	void evaluate(DATA_TYPE2 & result, const DATA_TYPE3 & val) const;

	void integral();
	void times_x(const unsigned int order);

	void pow(UnivariatePolynomial<DATA_TYPE> & result, const unsigned int exponent) const;

	template <class DATA_TYPE2, class DATA_TYPE3>
	void ctrunc(DATA_TYPE2 & remainder, const unsigned int order, const std::vector<DATA_TYPE3> & val_exp_table);

	template <class DATA_TYPE2, class DATA_TYPE3>
	void ctrunc(DATA_TYPE2 & remainder, const unsigned int order, const DATA_TYPE3 & val);

	template <class DATA_TYPE2, class DATA_TYPE3>
	void ctrunc(DATA_TYPE3 & remainder1, DATA_TYPE3 & remainder2, const unsigned int order, const std::vector<DATA_TYPE3> & val1_exp_table, const std::vector<DATA_TYPE3> & val2_exp_table);

//	template <class DATA_TYPE2, class DATA_TYPE3>
//	void ctrunc(DATA_TYPE2 & remainder1, DATA_TYPE3 & remainder2, const unsigned int order, const DATA_TYPE3 & val1, const DATA_TYPE3 & val2);

	void nctrunc(const unsigned int order);

//	void substitute(UnivariatePolynomial<DATA_TYPE> & result, const std::vector<UnivariatePolynomial<DATA_TYPE> > & x_exp_table) const;
	void substitute(UnivariatePolynomial<DATA_TYPE> & result, const UnivariatePolynomial<DATA_TYPE> & x) const;


	UnivariatePolynomial<DATA_TYPE> & operator = (const UnivariatePolynomial<DATA_TYPE> & up);
	UnivariatePolynomial<DATA_TYPE> & operator = (const DATA_TYPE & c);

	template <class DATA_TYPE2>
	UnivariatePolynomial<DATA_TYPE> & operator += (const UnivariatePolynomial<DATA_TYPE2> & up);

	template <class DATA_TYPE2>
	UnivariatePolynomial<DATA_TYPE> & operator -= (const UnivariatePolynomial<DATA_TYPE2> & up);

	template <class DATA_TYPE2>
	UnivariatePolynomial<DATA_TYPE> & operator *= (const UnivariatePolynomial<DATA_TYPE2> & up);

	template <class DATA_TYPE2>
	UnivariatePolynomial<DATA_TYPE> & operator += (const DATA_TYPE2 & c);

	template <class DATA_TYPE2>
	UnivariatePolynomial<DATA_TYPE> & operator -= (const DATA_TYPE2 & c);

	template <class DATA_TYPE2>
	UnivariatePolynomial<DATA_TYPE> & operator *= (const DATA_TYPE2 & c);

	template <class DATA_TYPE2>
	UnivariatePolynomial<DATA_TYPE> & operator /= (const DATA_TYPE2 & c);


	template <class DATA_TYPE2>
	UnivariatePolynomial<DATA_TYPE> operator + (const UnivariatePolynomial<DATA_TYPE2> & up) const;

	template <class DATA_TYPE2>
	UnivariatePolynomial<DATA_TYPE> operator - (const UnivariatePolynomial<DATA_TYPE2> & up) const;

	template <class DATA_TYPE2>
	UnivariatePolynomial<DATA_TYPE> operator * (const UnivariatePolynomial<DATA_TYPE2> & up) const;

	template <class DATA_TYPE2>
	UnivariatePolynomial<DATA_TYPE> operator + (const DATA_TYPE2 & c) const;

	template <class DATA_TYPE2>
	UnivariatePolynomial<DATA_TYPE> operator - (const DATA_TYPE2 & c) const;

	template <class DATA_TYPE2>
	UnivariatePolynomial<DATA_TYPE> operator * (const DATA_TYPE2 & c) const;

	template <class DATA_TYPE2>
	UnivariatePolynomial<DATA_TYPE> operator / (const DATA_TYPE2 & c) const;


	template <class DATA_TYPE2>
	UnivariateTaylorModel<DATA_TYPE2> operator + (const UnivariateTaylorModel<DATA_TYPE2> & utm) const;

	template <class DATA_TYPE2>
	UnivariateTaylorModel<DATA_TYPE2> operator - (const UnivariateTaylorModel<DATA_TYPE2> & utm) const;

	template <class DATA_TYPE2>
	UnivariateTaylorModel<DATA_TYPE2> operator * (const UnivariateTaylorModel<DATA_TYPE2> & utm) const;


	template <class DATA_TYPE2>
	operator UnivariatePolynomial<DATA_TYPE2> () const;


	friend std::ostream & operator << <DATA_TYPE> (std::ostream & os, const UnivariatePolynomial<DATA_TYPE> & up);

	template <class DATA_TYPE2>
	friend class UnivariatePolynomial;

	template <class DATA_TYPE2>
	friend class UnivariateTaylorModel;

	template <class DATA_TYPE2>
	friend class TaylorModel;
};


template <class DATA_TYPE>
UnivariatePolynomial<DATA_TYPE>::UnivariatePolynomial()
{
}

template <class DATA_TYPE>
UnivariatePolynomial<DATA_TYPE>::UnivariatePolynomial(const std::vector<DATA_TYPE> & coeffs)
{
	coefficients = coeffs;
}

template <class DATA_TYPE>
UnivariatePolynomial<DATA_TYPE>::UnivariatePolynomial(const UnivariatePolynomial<DATA_TYPE> & up)
{
	coefficients = up.coefficients;
}

template <class DATA_TYPE>
UnivariatePolynomial<DATA_TYPE>::UnivariatePolynomial(const DATA_TYPE & c)
{
	coefficients.push_back(c);
}

template <class DATA_TYPE>
UnivariatePolynomial<DATA_TYPE>::UnivariatePolynomial(const double c)
{
	coefficients.push_back(c);
}

template <class DATA_TYPE>
UnivariatePolynomial<DATA_TYPE>::UnivariatePolynomial(const DATA_TYPE & c, const unsigned int degree)
{
	for(unsigned int i=0; i<degree; ++i)
	{
		coefficients.push_back(0);
	}

	coefficients.push_back((DATA_TYPE)c);
}

template <class DATA_TYPE>
UnivariatePolynomial<DATA_TYPE>::~UnivariatePolynomial()
{
}
/*
template <class DATA_TYPE>
bool UnivariatePolynomial<DATA_TYPE>::set(const std::string & strPolynomial)
{
	std::string prefix(str_prefix_univariate_polynomial);
	std::string suffix(str_suffix);

	std::string input = prefix + strPolynomial + suffix;

	parseUnivariatePolynomial(input);		// call the parser

	*this = (UnivariatePolynomial<DATA_TYPE>)(up_parseresult);

	return true;
}
*/
template <class DATA_TYPE>
void UnivariatePolynomial<DATA_TYPE>::clear()
{
	coefficients.clear();
}

template <class DATA_TYPE>
void UnivariatePolynomial<DATA_TYPE>::set2zero()
{
	coefficients.clear();
	coefficients.push_back(0);
}

template <class DATA_TYPE>
unsigned int UnivariatePolynomial<DATA_TYPE>::degree() const
{
	if(coefficients.size() == 0)
	{
		return 0;
	}
	else
	{
		return coefficients.size() - 1;
	}
}

template <class DATA_TYPE>
bool UnivariatePolynomial<DATA_TYPE>::isZero() const
{
	unsigned int n = coefficients.size();

	if(n == 0)
	{
		return true;
	}
	else
	{
		bool bZero = true;
		for(unsigned int i=0; i<n; ++i)
		{
			if(!(coefficients[i] == 0))
			{
				bZero = false;
				break;
			}
		}

		return bZero;
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2, class DATA_TYPE3>
void UnivariatePolynomial<DATA_TYPE>::round(DATA_TYPE2 & remainder, const DATA_TYPE3 & val)
{
	UnivariatePolynomial<DATA_TYPE> upTemp;

	for(unsigned int i=0; i<coefficients.size(); ++i)
	{
		DATA_TYPE I;
		coefficients[i].round(I);
		upTemp.coefficients.push_back(I);
	}

	upTemp.evaluate(remainder, val);
}

template <class DATA_TYPE>
template <class DATA_TYPE2, class DATA_TYPE3>
void UnivariatePolynomial<DATA_TYPE>::evaluate(DATA_TYPE2 & result, const std::vector<DATA_TYPE3> & val_exp_table) const
{
	if(coefficients.size() == 0)
	{
		result = 0;
	}
	else
	{
		result = coefficients[0];

		for(unsigned int i=1; i<coefficients.size(); ++i)
		{
			result += coefficients[i] * val_exp_table[i];
		}
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2, class DATA_TYPE3>
void UnivariatePolynomial<DATA_TYPE>::evaluate(DATA_TYPE2 & result, const DATA_TYPE3 & val) const
{
	if(coefficients.size() == 0)
	{
		result = 0;
	}
	else
	{
		result = coefficients[coefficients.size()-1];

		for(int i=coefficients.size()-2; i>=0; --i)
		{
			result = result * val + coefficients[i];
		}
	}
}

template <class DATA_TYPE>
void UnivariatePolynomial<DATA_TYPE>::integral()
{
	if(this->isZero())
	{
		return;
	}
	else if(coefficients.size() == 1)
	{
		coefficients.push_back(coefficients[0]);
		coefficients[0] = 0;
	}
	else
	{
		coefficients.push_back(0);

		for(unsigned int i=coefficients.size()-1; i>=1; --i)
		{
			coefficients[i] = coefficients[i-1] / i;
		}

		coefficients[0] = 0;
	}
}

template <class DATA_TYPE>
void UnivariatePolynomial<DATA_TYPE>::times_x(const unsigned int order)
{
	if(order != 0 && !this->isZero())
	{
		std::vector<DATA_TYPE> vecTemp;

		for(unsigned int i=0; i<order; ++i)
		{
			vecTemp.push_back(0);
		}

		for(unsigned int i=0; i<coefficients.size(); ++i)
		{
			vecTemp.push_back(coefficients[i]);
		}

		coefficients = vecTemp;
	}
}

template <class DATA_TYPE>
void UnivariatePolynomial<DATA_TYPE>::pow(UnivariatePolynomial<DATA_TYPE> & result, const unsigned int exponent) const
{
	if(exponent == 0)
	{
		UnivariatePolynomial<DATA_TYPE> upTemp = 1;
		result = upTemp;
	}
	else
	{
		UnivariatePolynomial<DATA_TYPE> upTemp = *this;
		result = *this;

		for(unsigned int d = exponent - 1; d > 0;)
		{
			if(d & 1)
			{
				result *= upTemp;
			}

			d >>= 1;

			if(d > 0)
			{
				upTemp *= upTemp;
			}
		}
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2, class DATA_TYPE3>
void UnivariatePolynomial<DATA_TYPE>::ctrunc(DATA_TYPE2 & remainder, const unsigned int order, const std::vector<DATA_TYPE3> & val_exp_table)
{
	UnivariatePolynomial<DATA_TYPE> trunc_part;

	for(unsigned int i=order+1; i<coefficients.size(); ++i)
	{
		trunc_part.coefficients.push_back(coefficients[i]);
	}

	for(int i=coefficients.size()-1; i>order; --i)
	{
		coefficients.pop_back();
	}

	trunc_part.evaluate(remainder, val_exp_table);

	remainder *= val_exp_table[order+1];
}

template <class DATA_TYPE>
template <class DATA_TYPE2, class DATA_TYPE3>
void UnivariatePolynomial<DATA_TYPE>::ctrunc(DATA_TYPE2 & remainder, const unsigned int order, const DATA_TYPE3 & val)
{
	UnivariatePolynomial<DATA_TYPE> trunc_part;

	for(unsigned int i=order+1; i<coefficients.size(); ++i)
	{
		trunc_part.coefficients.push_back(coefficients[i]);
	}

	for(int i=coefficients.size()-1; i>order; --i)
	{
		coefficients.pop_back();
	}

	trunc_part.evaluate(remainder, val);

	DATA_TYPE3 tmp = val;
	tmp.pow_assign(order + 1);

	remainder *= tmp;
}

template <class DATA_TYPE>
template <class DATA_TYPE2, class DATA_TYPE3>
void UnivariatePolynomial<DATA_TYPE>::ctrunc(DATA_TYPE3 & remainder1, DATA_TYPE3 & remainder2, const unsigned int order, const std::vector<DATA_TYPE3> & val1_exp_table, const std::vector<DATA_TYPE3> & val2_exp_table)
{
	UnivariatePolynomial<DATA_TYPE> trunc_part;

	for(unsigned int i=order+1; i<coefficients.size(); ++i)
	{
		trunc_part.coefficients.push_back(coefficients[i]);
	}

	for(int i=coefficients.size()-1; i>order; --i)
	{
		coefficients.pop_back();
	}

	trunc_part.evaluate(remainder1, val1_exp_table);
	trunc_part.evaluate(remainder2, val2_exp_table);

	remainder1 *= val1_exp_table[order+1];
	remainder2 *= val2_exp_table[order+1];
}

template <class DATA_TYPE>
void UnivariatePolynomial<DATA_TYPE>::nctrunc(const unsigned int order)
{
	for(int i=coefficients.size()-1; i>=order; --i)
	{
		coefficients.pop_back();
	}
}

template <class DATA_TYPE>
void UnivariatePolynomial<DATA_TYPE>::substitute(UnivariatePolynomial<DATA_TYPE> & result, const UnivariatePolynomial<DATA_TYPE> & x) const
{
	result.coefficients.clear();

	if(coefficients.size() > 0)
	{
		result.coefficients.push_back(coefficients[coefficients.size()-1]);

		for(int i=coefficients.size()-2; i>=0; --i)
		{
			result = result * x + coefficients[i];
		}
	}
}

template <class DATA_TYPE>
UnivariatePolynomial<DATA_TYPE> & UnivariatePolynomial<DATA_TYPE>::operator = (const UnivariatePolynomial<DATA_TYPE> & up)
{
	if(this == &up)
		return *this;

	coefficients = up.coefficients;

	return *this;
}

template <class DATA_TYPE>
UnivariatePolynomial<DATA_TYPE> & UnivariatePolynomial<DATA_TYPE>::operator = (const DATA_TYPE & c)
{
	coefficients.clear();
	coefficients.push_back(c);

	return *this;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariatePolynomial<DATA_TYPE> & UnivariatePolynomial<DATA_TYPE>::operator += (const UnivariatePolynomial<DATA_TYPE2> & up)
{
	unsigned int n = coefficients.size();
	unsigned int m = up.coefficients.size();

	if(n <= m)
	{
		for(unsigned int i=0; i<n; ++i)
		{
			coefficients[i] += up.coefficients[i];
		}

		for(unsigned int i=n; i<m; ++i)
		{
			coefficients.push_back(up.coefficients[i]);
		}
	}
	else
	{
		for(unsigned int i=0; i<m; ++i)
		{
			coefficients[i] += up.coefficients[i];
		}
	}

	return *this;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariatePolynomial<DATA_TYPE> & UnivariatePolynomial<DATA_TYPE>::operator -= (const UnivariatePolynomial<DATA_TYPE2> & up)
{
	unsigned int n = coefficients.size();
	unsigned int m = up.coefficients.size();

	if(n <= m)
	{
		for(unsigned int i=0; i<n; ++i)
		{
			coefficients[i] -= up.coefficients[i];
		}

		for(unsigned int i=n; i<m; ++i)
		{
			coefficients.push_back(-up.coefficients[i]);
		}
	}
	else
	{
		for(unsigned int i=0; i<m; ++i)
		{
			coefficients[i] -= up.coefficients[i];
		}
	}

	return *this;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariatePolynomial<DATA_TYPE> & UnivariatePolynomial<DATA_TYPE>::operator *= (const UnivariatePolynomial<DATA_TYPE2> & up)
{
	if(this->isZero())
	{
		return *this;
	}

	if(up.isZero())
	{
		this->set2zero();
		return *this;
	}

	unsigned int n = coefficients.size();
	unsigned int m = up.coefficients.size();

	std::vector<DATA_TYPE> result;

	for(unsigned int i=0; i<n; ++i)
	{
		for(unsigned int j=0; j<m; ++j)
		{
			unsigned int degree = i + j;
			if(result.size() <= degree)
			{
				result.push_back(coefficients[i] * up.coefficients[j]);
			}
			else
			{
				result[degree] += coefficients[i] * up.coefficients[j];
			}
		}
	}

	coefficients = result;

	return *this;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariatePolynomial<DATA_TYPE> & UnivariatePolynomial<DATA_TYPE>::operator += (const DATA_TYPE2 & c)
{
	coefficients[0] += c;
	return *this;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariatePolynomial<DATA_TYPE> & UnivariatePolynomial<DATA_TYPE>::operator -= (const DATA_TYPE2 & c)
{
	coefficients[0] -= c;
	return *this;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariatePolynomial<DATA_TYPE> & UnivariatePolynomial<DATA_TYPE>::operator *= (const DATA_TYPE2 & c)
{
	if(this->isZero())
	{
		return *this;
	}

	if(c == 0)
	{
		this->set2zero();
		return *this;
	}

	for(unsigned int i=0; i<coefficients.size(); ++i)
	{
		coefficients[i] *= c;
	}

	return *this;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariatePolynomial<DATA_TYPE> & UnivariatePolynomial<DATA_TYPE>::operator /= (const DATA_TYPE2 & c)
{
	for(unsigned int i=0; i<coefficients.size(); ++i)
	{
		coefficients[i] /= c;
	}

	return *this;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariatePolynomial<DATA_TYPE> UnivariatePolynomial<DATA_TYPE>::operator + (const UnivariatePolynomial<DATA_TYPE2> & up) const
{
	unsigned int n = coefficients.size();
	unsigned int m = up.coefficients.size();

	UnivariatePolynomial<DATA_TYPE> result;

	if(n <= m)
	{
		for(unsigned int i=0; i<n; ++i)
		{
			result.coefficients.push_back(coefficients[i] + up.coefficients[i]);
		}

		for(unsigned int i=n; i<m; ++i)
		{
			result.coefficients.push_back(up.coefficients[i]);
		}
	}
	else
	{
		for(unsigned int i=0; i<m; ++i)
		{
			result.coefficients.push_back(coefficients[i] + up.coefficients[i]);
		}

		for(unsigned int i=m; i<n; ++i)
		{
			result.coefficients.push_back(coefficients[i]);
		}
	}

	return result;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariatePolynomial<DATA_TYPE> UnivariatePolynomial<DATA_TYPE>::operator - (const UnivariatePolynomial<DATA_TYPE2> & up) const
{
	unsigned int n = coefficients.size();
	unsigned int m = up.coefficients.size();

	UnivariatePolynomial<DATA_TYPE> result;

	if(n <= m)
	{
		for(unsigned int i=0; i<n; ++i)
		{
			result.coefficients.push_back(coefficients[i] - up.coefficients[i]);
		}

		for(unsigned int i=n; i<m; ++i)
		{
			result.coefficients.push_back(-up.coefficients[i]);
		}
	}
	else
	{
		for(unsigned int i=0; i<m; ++i)
		{
			result.coefficients.push_back(coefficients[i] + up.coefficients[i]);
		}

		for(unsigned int i=m; i<n; ++i)
		{
			result.coefficients.push_back(coefficients[i]);
		}
	}

	return result;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariatePolynomial<DATA_TYPE> UnivariatePolynomial<DATA_TYPE>::operator * (const UnivariatePolynomial<DATA_TYPE2> & up) const
{
	UnivariatePolynomial<DATA_TYPE> result;

	if(this->isZero() || up.isZero())
	{
		return result;
	}

	unsigned int n = coefficients.size();
	unsigned int m = up.coefficients.size();

	for(unsigned int i=0; i<n; ++i)
	{
		for(unsigned int j=0; j<m; ++j)
		{
			unsigned int degree = i + j;
			if(result.coefficients.size() <= degree)
			{
				result.coefficients.push_back(coefficients[i] * up.coefficients[j]);
			}
			else
			{
				result.coefficients[degree] += coefficients[i] * up.coefficients[j];
			}
		}
	}

	return result;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariatePolynomial<DATA_TYPE> UnivariatePolynomial<DATA_TYPE>::operator + (const DATA_TYPE2 & c) const
{
	UnivariatePolynomial<DATA_TYPE> result = *this;
	result.coefficients[0] += c;

	return result;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariatePolynomial<DATA_TYPE> UnivariatePolynomial<DATA_TYPE>::operator - (const DATA_TYPE2 & c) const
{
	UnivariatePolynomial<DATA_TYPE> result = *this;
	result.coefficients[0] -= c;

	return result;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariatePolynomial<DATA_TYPE> UnivariatePolynomial<DATA_TYPE>::operator * (const DATA_TYPE2 & c) const
{
	UnivariatePolynomial<DATA_TYPE> result;

	if(this->isZero() || c == 0)
	{
		return result;
	}

	for(unsigned int i=0; i<coefficients.size(); ++i)
	{
		result.coefficients.push_back(coefficients[i] * c);
	}

	return result;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariatePolynomial<DATA_TYPE> UnivariatePolynomial<DATA_TYPE>::operator / (const DATA_TYPE2 & c) const
{
	UnivariatePolynomial<DATA_TYPE> result;

	for(unsigned int i=0; i<coefficients.size(); ++i)
	{
		result.coefficients.push_back(coefficients[i] / c);
	}

	return result;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariateTaylorModel<DATA_TYPE2> UnivariatePolynomial<DATA_TYPE>::operator + (const UnivariateTaylorModel<DATA_TYPE2> & utm) const
{
	UnivariateTaylorModel<DATA_TYPE2> result = utm;
	result += *this;

	return result;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariateTaylorModel<DATA_TYPE2> UnivariatePolynomial<DATA_TYPE>::operator - (const UnivariateTaylorModel<DATA_TYPE2> & utm) const
{
	UnivariateTaylorModel<DATA_TYPE2> result = utm;
	result -= *this;

	return result;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariateTaylorModel<DATA_TYPE2> UnivariatePolynomial<DATA_TYPE>::operator * (const UnivariateTaylorModel<DATA_TYPE2> & utm) const
{
	UnivariateTaylorModel<DATA_TYPE2> result = utm;
	result *= *this;

	return result;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariatePolynomial<DATA_TYPE>::operator UnivariatePolynomial<DATA_TYPE2> () const
{
	UnivariatePolynomial<DATA_TYPE2> result;

	unsigned int n = coefficients.size();

	for(unsigned int i=0; i<n; ++i)
	{
		DATA_TYPE2 tmp(coefficients[i]);
		result.coefficients.push_back(tmp);
	}

	return result;
}


template <class DATA_TYPE>
std::ostream & operator << (std::ostream & os, const UnivariatePolynomial<DATA_TYPE> & up)
{
	if(up.coefficients.size() == 0)
	{
		std::cout << "0";
		return os;
	}

	bool bfirst = true;

	if(!(up.coefficients[0] == 0))
	{
		up.coefficients[0].dump(stdout);
		bfirst = false;
	}

	for(unsigned int i=1; i<up.coefficients.size(); ++i)
	{
		if(!(up.coefficients[i] == 0))
		{
			if(!bfirst)
			{
				std::cout << " + ";
			}
			else
			{
				bfirst = false;
			}

			std::cout << " ( " << up.coefficients[i] << " ) " << "*t";

			if(i > 1)
			{
				std::cout << '^' << i;
			}
		}
	}

	if(bfirst)
	{
		std::cout << "0";
	}

	return os;
}

}

#endif /* UNIVARIATEPOLYNOMIAL_H_ */
