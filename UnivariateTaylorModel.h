/*---
  Flow*: A Verification Tool for Cyber-Physical Systems.
  Authors: Xin Chen, Sriram Sankaranarayanan, and Erika Abraham.
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.

  The code is released as is under the GNU General Public License (GPL).
---*/

#ifndef UNIVARIATETAYLORMODEL_H_
#define UNIVARIATETAYLORMODEL_H_

#include "UnivariatePolynomial.h"


namespace flowstar
{

template <class DATA_TYPE>
class UnivariateTaylorModel;

template <class DATA_TYPE>
std::ostream & operator << (std::ostream & os, const UnivariateTaylorModel<DATA_TYPE> & utm);


template <class DATA_TYPE>
class UnivariateTaylorModel
{
protected:
	UnivariatePolynomial<DATA_TYPE> expansion;
	Interval remainder;

public:
	UnivariateTaylorModel();
	UnivariateTaylorModel(const DATA_TYPE & c);
	UnivariateTaylorModel(const double c);
	UnivariateTaylorModel(const Interval & I);
	UnivariateTaylorModel(const UnivariatePolynomial<DATA_TYPE> & up);
	UnivariateTaylorModel(const UnivariatePolynomial<DATA_TYPE> & up, const Interval & I);
	UnivariateTaylorModel(const UnivariateTaylorModel<DATA_TYPE> & utm);
	~UnivariateTaylorModel();

	void evaluate(Interval & result) const;

	template <class DATA_TYPE2>
	void evaluate(Interval & result, const std::vector<DATA_TYPE2> & val_exp_table) const;

	template <class DATA_TYPE2>
	void evaluate(Real & result, const std::vector<DATA_TYPE2> & val_exp_table) const;

	template <class DATA_TYPE2>
	void evaluate(UnivariateTaylorModel<DATA_TYPE> & result, const std::vector<DATA_TYPE2> & val_exp_table) const;

	template <class DATA_TYPE2>
	void ctrunc(const unsigned int order, const std::vector<DATA_TYPE2> & val_exp_table);

	void ctrunc(const unsigned int order);
	void nctrunc(const unsigned int order);

	template <class DATA_TYPE2>
	void integral(const DATA_TYPE2 & val);

	void setRemainder(const Interval & I);
	void getRemainder(Interval & I) const;
	void addRemainder(const Interval & I);
	double remainderSize() const;


	UnivariateTaylorModel<DATA_TYPE> & operator = (const UnivariateTaylorModel<DATA_TYPE> & utm);
	UnivariateTaylorModel<DATA_TYPE> & operator = (const UnivariatePolynomial<DATA_TYPE> & up);
	UnivariateTaylorModel<DATA_TYPE> & operator = (const DATA_TYPE & c);

	template <class DATA_TYPE2>
	UnivariateTaylorModel<DATA_TYPE> & operator += (const UnivariateTaylorModel<DATA_TYPE2> & utm);

	template <class DATA_TYPE2>
	UnivariateTaylorModel<DATA_TYPE> & operator -= (const UnivariateTaylorModel<DATA_TYPE2> & utm);

	template <class DATA_TYPE2>
	UnivariateTaylorModel<DATA_TYPE> & operator *= (const UnivariateTaylorModel<DATA_TYPE2> & utm);

	template <class DATA_TYPE2>
	UnivariateTaylorModel<DATA_TYPE> & operator += (const UnivariatePolynomial<DATA_TYPE2> & up);

	template <class DATA_TYPE2>
	UnivariateTaylorModel<DATA_TYPE> & operator -= (const UnivariatePolynomial<DATA_TYPE2> & up);

	template <class DATA_TYPE2>
	UnivariateTaylorModel<DATA_TYPE> & operator *= (const UnivariatePolynomial<DATA_TYPE2> & up);


	template <class DATA_TYPE2>
	UnivariateTaylorModel<DATA_TYPE> & operator += (const DATA_TYPE2 & c);

	template <class DATA_TYPE2>
	UnivariateTaylorModel<DATA_TYPE> & operator -= (const DATA_TYPE2 & c);

	template <class DATA_TYPE2>
	UnivariateTaylorModel<DATA_TYPE> & operator *= (const DATA_TYPE2 & c);

	template <class DATA_TYPE2>
	UnivariateTaylorModel<DATA_TYPE> & operator /= (const DATA_TYPE2 & c);


	template <class DATA_TYPE2>
	UnivariateTaylorModel<DATA_TYPE> operator + (const UnivariateTaylorModel<DATA_TYPE2> & utm) const;

	template <class DATA_TYPE2>
	UnivariateTaylorModel<DATA_TYPE> operator - (const UnivariateTaylorModel<DATA_TYPE2> & utm) const;

	template <class DATA_TYPE2>
	UnivariateTaylorModel<DATA_TYPE> operator * (const UnivariateTaylorModel<DATA_TYPE2> & utm) const;

	template <class DATA_TYPE2>
	UnivariateTaylorModel<DATA_TYPE> operator + (const UnivariatePolynomial<DATA_TYPE2> & up) const;

	template <class DATA_TYPE2>
	UnivariateTaylorModel<DATA_TYPE> operator - (const UnivariatePolynomial<DATA_TYPE2> & up) const;

	template <class DATA_TYPE2>
	UnivariateTaylorModel<DATA_TYPE> operator * (const UnivariatePolynomial<DATA_TYPE2> & up) const;


	template <class DATA_TYPE2>
	UnivariateTaylorModel<DATA_TYPE> operator + (const DATA_TYPE2 & c) const;

	template <class DATA_TYPE2>
	UnivariateTaylorModel<DATA_TYPE> operator - (const DATA_TYPE2 & c) const;

	template <class DATA_TYPE2>
	UnivariateTaylorModel<DATA_TYPE> operator * (const DATA_TYPE2 & c) const;

	template <class DATA_TYPE2>
	UnivariateTaylorModel<DATA_TYPE> operator / (const DATA_TYPE2 & c) const;

	bool operator != (const int n) const;


	template <class DATA_TYPE2>
	operator UnivariateTaylorModel<DATA_TYPE2> () const;

	friend std::ostream & operator << <DATA_TYPE> (std::ostream & os, const UnivariateTaylorModel<DATA_TYPE> & utm);


	template <class DATA_TYPE2>
	friend class UnivariateTaylorModel;

	template <class DATA_TYPE2>
	friend class TaylorModel;

	friend class Real;
//	friend class LinearFlowpipe;
};





template <class DATA_TYPE>
UnivariateTaylorModel<DATA_TYPE>::UnivariateTaylorModel()
{
}

template <class DATA_TYPE>
UnivariateTaylorModel<DATA_TYPE>::UnivariateTaylorModel(const DATA_TYPE & c)
{
	UnivariatePolynomial<DATA_TYPE> tmp(c);
	expansion = tmp;
}

template <class DATA_TYPE>
UnivariateTaylorModel<DATA_TYPE>::UnivariateTaylorModel(const double c)
{
	UnivariatePolynomial<DATA_TYPE> tmp(c);
	expansion = tmp;
}

template <class DATA_TYPE>
UnivariateTaylorModel<DATA_TYPE>::UnivariateTaylorModel(const Interval & I)
{
	Interval J = I;
	Real c;

	J.remove_midpoint(c);

	UnivariatePolynomial<DATA_TYPE> tmp(c);
	expansion = tmp;
	remainder = J;
}

template <class DATA_TYPE>
UnivariateTaylorModel<DATA_TYPE>::UnivariateTaylorModel(const UnivariatePolynomial<DATA_TYPE> & up)
{
	expansion = up;
}

template <class DATA_TYPE>
UnivariateTaylorModel<DATA_TYPE>::UnivariateTaylorModel(const UnivariatePolynomial<DATA_TYPE> & up, const Interval & I)
{
	expansion = up;
	remainder = I;
}

template <class DATA_TYPE>
UnivariateTaylorModel<DATA_TYPE>::UnivariateTaylorModel(const UnivariateTaylorModel<DATA_TYPE> & utm)
{
	expansion = utm.expansion;
	remainder = utm.remainder;
}

template <class DATA_TYPE>
UnivariateTaylorModel<DATA_TYPE>::~UnivariateTaylorModel()
{
}

template <class DATA_TYPE>
void UnivariateTaylorModel<DATA_TYPE>::evaluate(Interval & result) const
{
	expansion.evaluate(result, interval_utm_setting.val_exp_table);
	result += remainder;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void UnivariateTaylorModel<DATA_TYPE>::evaluate(Interval & result, const std::vector<DATA_TYPE2> & val_exp_table) const
{
	expansion.evaluate(result, val_exp_table);
	result += remainder;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void UnivariateTaylorModel<DATA_TYPE>::evaluate(Real & result, const std::vector<DATA_TYPE2> & val_exp_table) const
{
	if(expansion.coefficients.size() == 0)
	{
		result = 0;
	}
	else
	{
		result = expansion.coefficients[0];

		for(unsigned int i=1; i<expansion.coefficients.size(); ++i)
		{
			result += expansion.coefficients[i] * (val_exp_table[i].toReal());
		}
	}
}

template <>
template <>
inline void UnivariateTaylorModel<Real>::evaluate<Real>(Real & result, const std::vector<Real> & val_exp_table) const
{
	if(expansion.coefficients.size() == 0)
	{
		result = 0;
	}
	else
	{
		result = expansion.coefficients[0];

		for(unsigned int i=1; i<expansion.coefficients.size(); ++i)
		{
			result += expansion.coefficients[i] * val_exp_table[i];
		}
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void UnivariateTaylorModel<DATA_TYPE>::evaluate(UnivariateTaylorModel<DATA_TYPE> & result, const std::vector<DATA_TYPE2> & val_exp_table) const
{
	Interval I;
	expansion.evaluate(I, val_exp_table);

	Real c;
	I.remove_midpoint(c);

	result.expansion = c;
	result.remainder = remainder + I;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void UnivariateTaylorModel<DATA_TYPE>::ctrunc(const unsigned int order, const std::vector<DATA_TYPE2> & val_exp_table)
{
	Interval tmp;
	expansion.ctrunc(tmp, order, val_exp_table);
	remainder += tmp;
}

template <class DATA_TYPE>
void UnivariateTaylorModel<DATA_TYPE>::ctrunc(const unsigned int order)
{
	Interval tmp;
	expansion.ctrunc(tmp, order, interval_utm_setting.val_exp_table);
	remainder += tmp;
}

template <class DATA_TYPE>
void UnivariateTaylorModel<DATA_TYPE>::nctrunc(const unsigned int order)
{
	expansion.nctrunc(order);
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void UnivariateTaylorModel<DATA_TYPE>::integral(const DATA_TYPE2 & val)
{
	expansion.integral();
	remainder *= val;
}

template <class DATA_TYPE>
void UnivariateTaylorModel<DATA_TYPE>::setRemainder(const Interval & I)
{
	remainder = I;
}

template <class DATA_TYPE>
void UnivariateTaylorModel<DATA_TYPE>::getRemainder(Interval & I) const
{
	I = remainder;
}

template <class DATA_TYPE>
void UnivariateTaylorModel<DATA_TYPE>::addRemainder(const Interval & I)
{
	remainder += I;
}

template <class DATA_TYPE>
double UnivariateTaylorModel<DATA_TYPE>::remainderSize() const
{
	return remainder.mag();
}

template <class DATA_TYPE>
UnivariateTaylorModel<DATA_TYPE> & UnivariateTaylorModel<DATA_TYPE>::operator = (const UnivariateTaylorModel<DATA_TYPE> & utm)
{
	if(this == &utm)
		return *this;

	expansion = utm.expansion;
	remainder = utm.remainder;

	return *this;
}

template <class DATA_TYPE>
UnivariateTaylorModel<DATA_TYPE> & UnivariateTaylorModel<DATA_TYPE>::operator = (const UnivariatePolynomial<DATA_TYPE> & up)
{
	expansion = up;
	remainder = 0;
}

template <class DATA_TYPE>
UnivariateTaylorModel<DATA_TYPE> & UnivariateTaylorModel<DATA_TYPE>::operator = (const DATA_TYPE & c)
{
	expansion = c;
	remainder = 0;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariateTaylorModel<DATA_TYPE> & UnivariateTaylorModel<DATA_TYPE>::operator += (const UnivariateTaylorModel<DATA_TYPE2> & utm)
{
	expansion += utm.expansion;
	remainder += utm.remainder;

	return *this;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariateTaylorModel<DATA_TYPE> & UnivariateTaylorModel<DATA_TYPE>::operator -= (const UnivariateTaylorModel<DATA_TYPE2> & utm)
{
	expansion -= utm.expansion;
	remainder -= utm.remainder;

	return *this;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariateTaylorModel<DATA_TYPE> & UnivariateTaylorModel<DATA_TYPE>::operator *= (const UnivariateTaylorModel<DATA_TYPE2> & utm)
{
	Interval tmp1;
	if(!(remainder == 0))
	{
		utm.evaluate(tmp1, interval_utm_setting.val_exp_table);
		tmp1 *= remainder;
	}

	Interval tmp2;
	if(!(utm.remainder == 0))
	{
		expansion.evaluate(tmp2, interval_utm_setting.val_exp_table);
		tmp2 *= utm.remainder;
	}

	Interval tmp3;
	expansion *= utm.expansion;
	expansion.ctrunc(tmp3, interval_utm_setting.order, interval_utm_setting.val_exp_table);

	remainder = tmp1 + tmp2 + tmp3;

	return *this;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariateTaylorModel<DATA_TYPE> & UnivariateTaylorModel<DATA_TYPE>::operator += (const UnivariatePolynomial<DATA_TYPE2> & up)
{
	expansion += up;

	return *this;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariateTaylorModel<DATA_TYPE> & UnivariateTaylorModel<DATA_TYPE>::operator -= (const UnivariatePolynomial<DATA_TYPE2> & up)
{
	expansion -= up;

	return *this;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariateTaylorModel<DATA_TYPE> & UnivariateTaylorModel<DATA_TYPE>::operator *= (const UnivariatePolynomial<DATA_TYPE2> & up)
{
	Interval tmp;

	if(!(remainder == 0))
	{
		up.evaluate(tmp, interval_utm_setting.val_exp_table);
		remainder += tmp * remainder;
	}

	expansion *= up;
	expansion.ctrunc(tmp, interval_utm_setting.order, interval_utm_setting.val_exp_table);
	remainder += tmp;

	return *this;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariateTaylorModel<DATA_TYPE> & UnivariateTaylorModel<DATA_TYPE>::operator += (const DATA_TYPE2 & c)
{
	expansion += c;
	return *this;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariateTaylorModel<DATA_TYPE> & UnivariateTaylorModel<DATA_TYPE>::operator -= (const DATA_TYPE2 & c)
{
	expansion -= c;
	return *this;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariateTaylorModel<DATA_TYPE> & UnivariateTaylorModel<DATA_TYPE>::operator *= (const DATA_TYPE2 & c)
{
	expansion *= c;
	remainder *= c;
	return *this;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariateTaylorModel<DATA_TYPE> & UnivariateTaylorModel<DATA_TYPE>::operator /= (const DATA_TYPE2 & c)
{
	expansion /= c;
	remainder /= c;
	return *this;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariateTaylorModel<DATA_TYPE> UnivariateTaylorModel<DATA_TYPE>::operator + (const UnivariateTaylorModel<DATA_TYPE2> & utm) const
{
	UnivariateTaylorModel<DATA_TYPE> result;

	result.expansion = expansion + utm.expansion;
	result.remainder = remainder + utm.remainder;

	return result;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariateTaylorModel<DATA_TYPE> UnivariateTaylorModel<DATA_TYPE>::operator - (const UnivariateTaylorModel<DATA_TYPE2> & utm) const
{
	UnivariateTaylorModel<DATA_TYPE> result;

	result.expansion = expansion - utm.expansion;
	result.remainder = remainder - utm.remainder;

	return result;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariateTaylorModel<DATA_TYPE> UnivariateTaylorModel<DATA_TYPE>::operator * (const UnivariateTaylorModel<DATA_TYPE2> & utm) const
{
	UnivariateTaylorModel<DATA_TYPE> result = *this;
	result *= utm;
	return result;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariateTaylorModel<DATA_TYPE> UnivariateTaylorModel<DATA_TYPE>::operator + (const UnivariatePolynomial<DATA_TYPE2> & up) const
{
	UnivariateTaylorModel<DATA_TYPE> result = *this;
	result += up;
	return result;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariateTaylorModel<DATA_TYPE> UnivariateTaylorModel<DATA_TYPE>::operator - (const UnivariatePolynomial<DATA_TYPE2> & up) const
{
	UnivariateTaylorModel<DATA_TYPE> result = *this;
	result -= up;
	return result;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariateTaylorModel<DATA_TYPE> UnivariateTaylorModel<DATA_TYPE>::operator * (const UnivariatePolynomial<DATA_TYPE2> & up) const
{
	UnivariateTaylorModel<DATA_TYPE> result = *this;
	result *= up;
	return result;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariateTaylorModel<DATA_TYPE> UnivariateTaylorModel<DATA_TYPE>::operator + (const DATA_TYPE2 & c) const
{
	UnivariateTaylorModel<DATA_TYPE> result = *this;
	result += c;
	return result;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariateTaylorModel<DATA_TYPE> UnivariateTaylorModel<DATA_TYPE>::operator - (const DATA_TYPE2 & c) const
{
	UnivariateTaylorModel<DATA_TYPE> result = *this;
	result -= c;
	return result;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariateTaylorModel<DATA_TYPE> UnivariateTaylorModel<DATA_TYPE>::operator * (const DATA_TYPE2 & c) const
{
	UnivariateTaylorModel<DATA_TYPE> result = *this;
	result *= c;
	return result;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariateTaylorModel<DATA_TYPE> UnivariateTaylorModel<DATA_TYPE>::operator / (const DATA_TYPE2 & c) const
{
	UnivariateTaylorModel<DATA_TYPE> result = *this;
	result /= c;
	return result;
}

template <class DATA_TYPE>
bool UnivariateTaylorModel<DATA_TYPE>::operator != (const int n) const
{
	if(remainder == 0)
	{
		if(expansion.coefficients.size() == 0 && n == 0)
		{
			return false;
		}
		else if(expansion.coefficients.size() == 1 && expansion.coefficients[0] == n)
		{
			return false;
		}
	}

	return true;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
UnivariateTaylorModel<DATA_TYPE>::operator UnivariateTaylorModel<DATA_TYPE2> () const
{
	UnivariateTaylorModel<DATA_TYPE2> result;

	result.expansion = (UnivariatePolynomial<DATA_TYPE2>)(expansion);
	result.remainder = remainder;

	return result;
}

template <class DATA_TYPE>
std::ostream & operator << (std::ostream & os, const UnivariateTaylorModel<DATA_TYPE> & utm)
{
	os << utm.expansion << " + " << utm.remainder;

	return os;
}

}



#endif /* UNIVARIATETAYLORMODEL_H_ */
