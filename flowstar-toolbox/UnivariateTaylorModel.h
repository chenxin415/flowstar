/*---
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
class TaylorModel;

template <class DATA_TYPE>
std::ostream & operator << (std::ostream & os, const UnivariateTaylorModel<DATA_TYPE> & utm);


template <class DATA_TYPE>
class UnivariateTaylorModel
{
public:
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

	void toMultivariateTaylorModel(TaylorModel<DATA_TYPE> & result, const unsigned int numVars) const;

	template <class DATA_TYPE2>
	void evaluate(Interval & result, const DATA_TYPE2 & val) const;

	template <class DATA_TYPE2>
	void evaluate(UnivariateTaylorModel<DATA_TYPE> & result, const DATA_TYPE2 & val) const;

	template <class DATA_TYPE2>
	void evaluate_assign(const DATA_TYPE2 & val);

	void ctrunc(const unsigned int order);
	void nctrunc(const unsigned int order);

	template <class DATA_TYPE2>
	void integral(const DATA_TYPE2 & val);

	void intIntegral(const Interval & sup);
	void invIntegral(const Interval & sup);

	void setRemainder(const Interval & I);
	void getRemainder(Interval & I) const;
	void addRemainder(const Interval & I);
	double remainderSize() const;

	void constant(DATA_TYPE & c) const;

	void clear();

	UnivariateTaylorModel<DATA_TYPE> & operator = (const UnivariateTaylorModel<DATA_TYPE> & utm);
	UnivariateTaylorModel<DATA_TYPE> & operator = (const UnivariatePolynomial<DATA_TYPE> & up);
	UnivariateTaylorModel<DATA_TYPE> & operator = (const DATA_TYPE & c);
	UnivariateTaylorModel<DATA_TYPE> & operator = (const double c);

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

	void pow(UnivariateTaylorModel<DATA_TYPE> & result, const unsigned int n) const;

	void exp_taylor(UnivariateTaylorModel<DATA_TYPE> & result, const Interval & t, const unsigned int order, const Global_Setting & setting) const;
	void rec_taylor(UnivariateTaylorModel<DATA_TYPE> & result, const Interval & t, const unsigned int order, const Global_Setting & setting) const;
	void sin_taylor(UnivariateTaylorModel<DATA_TYPE> & result, const Interval & t, const unsigned int order, const Global_Setting & setting) const;
	void cos_taylor(UnivariateTaylorModel<DATA_TYPE> & result, const Interval & t, const unsigned int order, const Global_Setting & setting) const;
	void log_taylor(UnivariateTaylorModel<DATA_TYPE> & result, const Interval & t, const unsigned int order, const Global_Setting & setting) const;
	void sqrt_taylor(UnivariateTaylorModel<DATA_TYPE> & result, const Interval & t, const unsigned int order, const Global_Setting & setting) const;

	void sigmoid_taylor(UnivariateTaylorModel<DATA_TYPE> & result, const Interval & t, const unsigned int order, const Global_Setting & setting) const;
	void tanh_taylor(UnivariateTaylorModel<DATA_TYPE> & result, const Interval & t, const unsigned int order, const Global_Setting & setting) const;


	template <class DATA_TYPE2>
	operator UnivariateTaylorModel<DATA_TYPE2> () const;

	friend std::ostream & operator << <DATA_TYPE> (std::ostream & os, const UnivariateTaylorModel<DATA_TYPE> & utm);


	template <class DATA_TYPE2>
	friend class UnivariateTaylorModel;

	template <class DATA_TYPE2>
	friend class TaylorModel;

	friend class Real;
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
void UnivariateTaylorModel<DATA_TYPE>::toMultivariateTaylorModel(TaylorModel<DATA_TYPE> & result, const unsigned int numVars) const
{
	expansion.toMultivariatePolynomial(result.expansion, numVars);
	result.remainder = remainder;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void UnivariateTaylorModel<DATA_TYPE>::evaluate(Interval & result, const DATA_TYPE2 & val) const
{
	expansion.evaluate(result, val);
	result += remainder;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void UnivariateTaylorModel<DATA_TYPE>::evaluate(UnivariateTaylorModel<DATA_TYPE> & result, const DATA_TYPE2 & val) const
{
	Interval I;
	expansion.evaluate(I, val);

	Real c;
	I.remove_midpoint(c);

	result.expansion = c;
	result.remainder = remainder + I;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void UnivariateTaylorModel<DATA_TYPE>::evaluate_assign(const DATA_TYPE2 & val)
{
	DATA_TYPE2 c;
	expansion.evaluate(c, val);

	expansion = c;
}

template <class DATA_TYPE>
void UnivariateTaylorModel<DATA_TYPE>::ctrunc(const unsigned int order)
{
	Interval tmp;
	expansion.ctrunc(tmp, order, interval_utm_setting.val);
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
void UnivariateTaylorModel<DATA_TYPE>::intIntegral(const Interval & sup)
{
	UnivariateTaylorModel<DATA_TYPE> utm_sup(sup);

	expansion.integral();
	remainder *= sup;

	UnivariateTaylorModel<DATA_TYPE> utm_exp_sup;
	expansion.evaluate(utm_exp_sup, utm_sup);

	*this = utm_exp_sup - (*this);
}

template <class DATA_TYPE>
void UnivariateTaylorModel<DATA_TYPE>::invIntegral(const Interval & sup)
{
	Real t_end;
	sup.sup(t_end);

	expansion.integral();
	remainder *= sup;

	Real r;
	expansion.evaluate(r, t_end);
	*this *= -1;

	this->expansion.coefficients[0] += r;
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
void UnivariateTaylorModel<DATA_TYPE>::constant(DATA_TYPE & c) const
{
	if(expansion.coefficients.size() > 0)
		c = expansion.coefficients[0];
}

template <class DATA_TYPE>
void UnivariateTaylorModel<DATA_TYPE>::clear()
{
	expansion.clear();
	remainder = 0;
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

	return *this;
}

template <class DATA_TYPE>
UnivariateTaylorModel<DATA_TYPE> & UnivariateTaylorModel<DATA_TYPE>::operator = (const DATA_TYPE & c)
{
	expansion = c;
	remainder = 0;

	return *this;
}

template <class DATA_TYPE>
UnivariateTaylorModel<DATA_TYPE> & UnivariateTaylorModel<DATA_TYPE>::operator = (const double c)
{
	expansion = c;
	remainder = 0;

	return *this;
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
		utm.evaluate(tmp1, interval_utm_setting.val);
		tmp1 *= remainder;
	}

	Interval tmp2;
	if(!(utm.remainder == 0))
	{
		expansion.evaluate(tmp2, interval_utm_setting.val);
		tmp2 *= utm.remainder;
	}

	Interval tmp3;
	expansion *= utm.expansion;
	expansion.ctrunc(tmp3, interval_utm_setting.order, interval_utm_setting.val);

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
		up.evaluate(tmp, interval_utm_setting.val);
		remainder += tmp * remainder;
	}

	expansion *= up;
	expansion.ctrunc(tmp, interval_utm_setting.order, interval_utm_setting.val);
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
void UnivariateTaylorModel<DATA_TYPE>::pow(UnivariateTaylorModel<DATA_TYPE> & result, const unsigned int n) const
{
	UnivariateTaylorModel<DATA_TYPE> utm = *this;

	UnivariateTaylorModel<DATA_TYPE> tmp = utm;
	result = utm;

	for(int d = n - 1; d > 0;)
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
void UnivariateTaylorModel<DATA_TYPE>::exp_taylor(UnivariateTaylorModel<DATA_TYPE> & result, const Interval & t, const unsigned int order, const Global_Setting & setting) const
{
	DATA_TYPE const_part;

	UnivariateTaylorModel<DATA_TYPE> utmF = *this;

	Interval tmRange;
	utmF.evaluate(tmRange, t);

	tmRange.remove_midpoint(const_part);

	utmF -= const_part;

	const_part.exp_assign();	// exp(c)


	if(tmRange.isZero())			// tm = c
	{
		UnivariateTaylorModel<DATA_TYPE> utmExp(const_part);
		result = utmExp;

		return;
	}

	UnivariatePolynomial<DATA_TYPE> polyOne(1);

	// to compute the expression 1 + F + (1/2!)F^2 + ... + (1/k!)F^k,
	// we evaluate its Horner form (...((1/(k-1))((1/k)*F+1)*F + 1) ... + 1)

	result.expansion = polyOne;
	result.remainder = 0;

	for(int i=order; i>0; --i)
	{
		result /= i;
		result *= utmF;
		result.expansion += polyOne;
	}

	result *= const_part;

	Interval rem;
	exp_taylor_remainder(rem, tmRange, order+1, setting);

	result.remainder += const_part * rem;
}

template <class DATA_TYPE>
void UnivariateTaylorModel<DATA_TYPE>::rec_taylor(UnivariateTaylorModel<DATA_TYPE> & result, const Interval & t, const unsigned int order, const Global_Setting & setting) const
{
	DATA_TYPE const_part;

	UnivariateTaylorModel<DATA_TYPE> utmF = *this;

	Interval tmRange;
	utmF.evaluate(tmRange, t);

	tmRange.remove_midpoint(const_part);

	Interval c_f = const_part;

	utmF -= const_part;

	const_part.rec_assign();	// 1/c

	if(tmRange.isZero())			// tm = c
	{
		UnivariateTaylorModel<DATA_TYPE> utmRec(const_part);
		result = utmRec;

		return;
	}

	UnivariatePolynomial<DATA_TYPE> polyOne(1);
	UnivariateTaylorModel<DATA_TYPE> utmF_c = utmF * const_part;


	// to compute the expression 1 - F/c + (F/c)^2 - ... + (-1)^k (F/c)^k,
	// we evaluate its Horner form (-1)*(...((-1)*(-F/c + 1)*F/c + 1)...) + 1

	result.expansion = polyOne;
	result.remainder = 0;


	for(int i=order; i>0; --i)
	{
		result *= -1;
		result *= utmF_c;
		result.expansion += polyOne;
	}

	result *= const_part;

	Interval rem;

	rec_taylor_remainder(rem, c_f, tmRange, order+1, setting);

	result.remainder += rem * const_part;
}

template <class DATA_TYPE>
void UnivariateTaylorModel<DATA_TYPE>::sin_taylor(UnivariateTaylorModel<DATA_TYPE> & result, const Interval & t, const unsigned int order, const Global_Setting & setting) const
{
	DATA_TYPE const_part;

	UnivariateTaylorModel<DATA_TYPE> tmF = *this;

	Interval tmRange;
	tmF.evaluate(tmRange, t);

	tmRange.remove_midpoint(const_part);

	tmF -= const_part;

	if(tmRange.isZero())
	{
		const_part.sin_assign();
		result.clear();
		result.expansion.coefficients.push_back(const_part);

		return;
	}

	DATA_TYPE sinc, cosc, msinc, mcosc;
	sinc = const_part.sin();
	cosc = const_part.cos();

	msinc = -sinc;
	mcosc = -cosc;

	UnivariateTaylorModel<DATA_TYPE> utmTmp(sinc);
	result = utmTmp;


	int k = 1;

	UnivariateTaylorModel<DATA_TYPE> utmPowerTmF(1);

	for(int i=1; i<=order; ++i, ++k)
	{
		k %= 4;

		utmPowerTmF *= tmF;

		UnivariateTaylorModel<DATA_TYPE> utmTmp2 = utmPowerTmF;

		switch(k)
		{
		case 0:
		{
			utmTmp2 *= sinc / i;
			break;
		}
		case 1:
		{
			utmTmp2 *= cosc / i;
			break;
		}
		case 2:
		{
			utmTmp2*= msinc / i;
			break;
		}
		case 3:
		{
			utmTmp2 *= mcosc / i;
			break;
		}
		}

		result += utmTmp2;
	}

	Interval rem;

	sin_taylor_remainder(rem, const_part, tmRange, order+1, setting);

	result.remainder += rem;
}

template <class DATA_TYPE>
void UnivariateTaylorModel<DATA_TYPE>::cos_taylor(UnivariateTaylorModel<DATA_TYPE> & result, const Interval & t, const unsigned int order, const Global_Setting & setting) const
{
	DATA_TYPE const_part;

	UnivariateTaylorModel<DATA_TYPE> utmF = *this;

	Interval tmRange;
	utmF.evaluate(tmRange, t);

	tmRange.remove_midpoint(const_part);

	utmF -= const_part;

	if(tmRange.isZero())
	{
		const_part.cos_assign();
		result.clear();
		result.expansion.coefficients.push_back(const_part);

		return;
	}

	DATA_TYPE sinc, cosc, msinc, mcosc;
	sinc = const_part.sin();
	cosc = const_part.cos();

	msinc = -sinc;
	mcosc = -cosc;


	UnivariateTaylorModel<DATA_TYPE> utmTmp(cosc);
	result = utmTmp;


	int k = 1;

	UnivariateTaylorModel<DATA_TYPE> utmPowerTmF(1);

	for(int i=1; i<=order; ++i, ++k)
	{
		k %= 4;

		utmPowerTmF *= utmF;

		UnivariateTaylorModel<DATA_TYPE> utmTmp2 = utmPowerTmF;

		switch(k)
		{
		case 0:
		{
			utmTmp2 *= cosc / i;
			break;
		}
		case 1:
		{
			utmTmp2 *= msinc / i;
			break;
		}
		case 2:
		{
			utmTmp2 *= mcosc / i;
			break;
		}
		case 3:
		{
			utmTmp2 *= sinc / i;
			break;
		}
		}

		result += utmTmp2;
	}

	Interval rem;

	cos_taylor_remainder(rem, const_part, tmRange, order+1, setting);

	result.remainder += rem;
}

template <class DATA_TYPE>
void UnivariateTaylorModel<DATA_TYPE>::log_taylor(UnivariateTaylorModel<DATA_TYPE> & result, const Interval & t, const unsigned int order, const Global_Setting & setting) const
{
	DATA_TYPE const_part;

	UnivariateTaylorModel<DATA_TYPE> utmF = *this;

	Interval tmRange;
	utmF.evaluate(tmRange, t);

	tmRange.remove_midpoint(const_part);

	utmF -= const_part;


//	DATA_TYPE C = const_part;

	const_part.log_assign();	// log(c)

	if(tmRange.isZero())			// tm = c
	{
		UnivariateTaylorModel<DATA_TYPE> utmLog(const_part);
		result = utmLog;

		return;
	}

	UnivariateTaylorModel<DATA_TYPE> utmF_c = utmF / const_part;
	result = utmF_c / order;


	for(int i=order-1; i>=1; --i)
	{
		result -= (1/i);
		result *= -1;

		result *= utmF_c;
	}

	UnivariateTaylorModel<DATA_TYPE> const_part_utm(const_part);
	result += const_part_utm;

	Interval rem;

	log_taylor_remainder(rem, tmRange / const_part, order+1);

	result.remainder += rem;
}

template <class DATA_TYPE>
void UnivariateTaylorModel<DATA_TYPE>::sqrt_taylor(UnivariateTaylorModel<DATA_TYPE> & result, const Interval & t, const unsigned int order, const Global_Setting & setting) const
{
	DATA_TYPE const_part;

	UnivariateTaylorModel<DATA_TYPE> utmF = *this;

	Interval tmRange;
	utmF.evaluate(tmRange, t);

	tmRange.remove_midpoint(const_part);

	utmF -= const_part;

	const_part.sqrt_assign();	// sqrt(c)

	if(tmRange.isZero())			// tm = c
	{
		UnivariateTaylorModel<DATA_TYPE> utmSqrt(const_part);
		result = utmSqrt;

		return;
	}

	UnivariateTaylorModel<DATA_TYPE> utmF_2c = utmF / (2*const_part);
	UnivariatePolynomial<DATA_TYPE> polyOne(1);

	result = utmF_2c;

	for(int i=order, j=2*order-3; i>=2; --i, j-=2)
	{
		result *= -(j / i);
		result.expansion += polyOne;
		result *= utmF_2c;
	}

	result.expansion += polyOne;

	result *= const_part;

	Interval rem;

	sqrt_taylor_remainder(rem, tmRange / const_part, order+1, setting);

	result.remainder += rem * const_part;
}

template <class DATA_TYPE>
void UnivariateTaylorModel<DATA_TYPE>::sigmoid_taylor(UnivariateTaylorModel<DATA_TYPE> & result, const Interval & t, const unsigned int order, const Global_Setting & setting) const
{
	UnivariateTaylorModel<DATA_TYPE> utm_x = *this;
	utm_x *= -1;

	UnivariateTaylorModel<DATA_TYPE> utm_tmp1;
	utm_x.exp_taylor(utm_tmp1, t, order, setting);

	utm_tmp1 += 1;

	utm_tmp1.rec_taylor(result, t, order, setting);
}

template <class DATA_TYPE>
void UnivariateTaylorModel<DATA_TYPE>::tanh_taylor(UnivariateTaylorModel<DATA_TYPE> & result, const Interval & t, const unsigned int order, const Global_Setting & setting) const
{
	UnivariateTaylorModel<DATA_TYPE> utm_2x = *this;
	utm_2x *= 2;

	UnivariateTaylorModel<DATA_TYPE> utm_exp_2x;
	utm_2x.exp_taylor(utm_exp_2x, t, order, setting);

	UnivariateTaylorModel<DATA_TYPE> utm_nom = utm_exp_2x - 1;
	UnivariateTaylorModel<DATA_TYPE> utm_dom = utm_exp_2x + 1;

	UnivariateTaylorModel<DATA_TYPE> utm_rec_dom;
	utm_dom.rec_taylor(utm_rec_dom, t, order, setting);

	result = utm_nom * utm_rec_dom;
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
