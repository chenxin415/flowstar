/*---
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.
  The code is released as is under the GNU General Public License (GPL).
---*/

#ifndef SETTINGS_H_
#define SETTINGS_H_

#include "Interval.h"
#include "Variables.h"


namespace flowstar
{

template <class DATA_TYPE>
class Polynomial;

template <class DATA_TYPE>
class Expression;


inline void compute_factorial_rec(std::vector<Real> & factorial_rec, const unsigned int order)
{
	factorial_rec.push_back(1);
	Real r(1);

	for(unsigned int i=1; i<=order; ++i)
	{
		r /= i;
		factorial_rec.push_back(r);
	}
}

inline void compute_power_4(std::vector<Real> & power_4, const unsigned int order)
{
	power_4.push_back(1);
	Real r(1);

	for(unsigned int i=1; i<=order; ++i)
	{
		r *= 4;
		power_4.push_back(r);
	}
}

inline void compute_double_factorial(std::vector<Real> & double_factorial, const unsigned int order)
{
	double_factorial.push_back(1);
	double_factorial.push_back(1);

	Real odd(1), even(1);

	for(unsigned int i=2; i<=order; ++i)
	{
		if(i%2 == 0)
		{
			even *= i;
			double_factorial.push_back(even);
		}
		else
		{
			odd *= i;
			double_factorial.push_back(odd);
		}
	}
}




template <class DATA_TYPE>
class UTM_Setting
{
public:
	unsigned int order;
	DATA_TYPE val;

public:
	UTM_Setting();
	UTM_Setting(const UTM_Setting<DATA_TYPE> & setting);
	~UTM_Setting();

	UTM_Setting<DATA_TYPE> & operator = (const UTM_Setting<DATA_TYPE> & setting);
};

template <class DATA_TYPE>
UTM_Setting<DATA_TYPE>::UTM_Setting()
{
	order = 0;
}

template <class DATA_TYPE>
UTM_Setting<DATA_TYPE>::UTM_Setting(const UTM_Setting<DATA_TYPE> & setting)
{
	order = setting.order;
	val = setting.val;
}

template <class DATA_TYPE>
UTM_Setting<DATA_TYPE>::~UTM_Setting()
{
}

template <class DATA_TYPE>
UTM_Setting<DATA_TYPE> & UTM_Setting<DATA_TYPE>::operator = (const UTM_Setting<DATA_TYPE> & setting)
{
	if(this == &setting)
		return *this;

	order = setting.order;
	val = setting.val;

	return *this;
}



class Taylor_Model_Setting
{
public:
	Variables variables;		// state variables
	Parameters parameters;		// constant parameters
	Interval cutoff_threshold;
	std::vector<Interval> domain;
	std::vector<Interval> remainder_estimation;
	unsigned int order;
	std::vector<Interval> step_exp_table;
	std::vector<Real> step_end_exp_table;

	double step_min;
	double step_max;
	unsigned int order_min;
	unsigned int order_max;

public:
	Taylor_Model_Setting();
	Taylor_Model_Setting(const Variables & vars);
	Taylor_Model_Setting(const Variables & vars, const Parameters & pars);
	Taylor_Model_Setting(const Taylor_Model_Setting & setting);
	~Taylor_Model_Setting();

	Taylor_Model_Setting & operator = (const Taylor_Model_Setting & setting);

	void initializeAdaptiveSettings(const double delta_min, const double delta_max, const unsigned int k_min, const unsigned int k_max);

	bool setCutoff(const Interval & cutoff);
	bool setStepsize(const double step, const unsigned int k);

	bool resetOrder(const double step, const unsigned int k);
	bool resetOrder(const unsigned int k);

	void setRemainderEstimation(const std::vector<Interval> & intVec);
	void setDomain(const std::vector<Interval> & intVec);

	void clear();

	bool bAdaptiveStepSize() const;
	bool bAdaptiveOrder() const;

	friend class Flowpipe;
	friend class Continuous_Reachability;
};



class Global_Setting
{
public:
	std::vector<Real> factorial_rec;
	std::vector<Real> power_4;
	std::vector<Real> double_factorial;

public:
	Global_Setting();
	Global_Setting(const Global_Setting & setting);
	~Global_Setting();

	bool resetOrder(const unsigned int k);

	Global_Setting & operator = (const Global_Setting & setting);

	void prepareForReachability(const unsigned int maxOrder);
};



template <class DATA_TYPE>
class Multivariate_Polynomial_Setting
{
public:
	std::string strPolynomial;
	Polynomial<DATA_TYPE> result;
	bool bDeterministic;
	Variables *pVars;

public:
	Multivariate_Polynomial_Setting();
	Multivariate_Polynomial_Setting(const Multivariate_Polynomial_Setting & setting);
	~Multivariate_Polynomial_Setting();

	Multivariate_Polynomial_Setting & operator = (const Multivariate_Polynomial_Setting & setting);

	void clear();
};

template <class DATA_TYPE>
Multivariate_Polynomial_Setting<DATA_TYPE>::Multivariate_Polynomial_Setting()
{
	bDeterministic = true;
	pVars = NULL;
}

template <class DATA_TYPE>
Multivariate_Polynomial_Setting<DATA_TYPE>::Multivariate_Polynomial_Setting(const Multivariate_Polynomial_Setting<DATA_TYPE> & setting)
{
	strPolynomial = setting.strPolynomial;
	result = setting.result;
	bDeterministic = setting.bDeterministic;
	pVars = setting.pVars;
}

template <class DATA_TYPE>
Multivariate_Polynomial_Setting<DATA_TYPE>::~Multivariate_Polynomial_Setting()
{
}

template <class DATA_TYPE>
Multivariate_Polynomial_Setting<DATA_TYPE> & Multivariate_Polynomial_Setting<DATA_TYPE>::operator = (const Multivariate_Polynomial_Setting<DATA_TYPE> & setting)
{
	if(this == &setting)
		return *this;

	strPolynomial = setting.strPolynomial;
	result = setting.result;
	bDeterministic = setting.bDeterministic;
	pVars = setting.pVars;

	return *this;
}

template <class DATA_TYPE>
void Multivariate_Polynomial_Setting<DATA_TYPE>::clear()
{
	result.clear();
	bDeterministic = true;
	pVars = NULL;
}





template <class DATA_TYPE>
class Expression_Setting
{
public:
	std::string strExpression;
	Expression<DATA_TYPE> result;
	bool bDeterministic;
	Variables *pVars;

public:
	Expression_Setting();
	Expression_Setting(const Expression_Setting & setting);
	~Expression_Setting();

	Expression_Setting & operator = (const Expression_Setting & setting);

	void clear();
};

template <class DATA_TYPE>
Expression_Setting<DATA_TYPE>::Expression_Setting()
{
	bDeterministic = true;
	pVars = NULL;
}

template <class DATA_TYPE>
Expression_Setting<DATA_TYPE>::Expression_Setting(const Expression_Setting<DATA_TYPE> & setting)
{
	strExpression = setting.strExpression;
	result = setting.result;
	bDeterministic = setting.bDeterministic;
	pVars = setting.pVars;
}

template <class DATA_TYPE>
Expression_Setting<DATA_TYPE>::~Expression_Setting()
{
}

template <class DATA_TYPE>
Expression_Setting<DATA_TYPE> & Expression_Setting<DATA_TYPE>::operator = (const Expression_Setting<DATA_TYPE> & setting)
{
	if(this == &setting)
		return *this;

	strExpression = setting.strExpression;
	result = setting.result;
	bDeterministic = setting.bDeterministic;
	pVars = setting.pVars;

	return *this;
}

template <class DATA_TYPE>
void Expression_Setting<DATA_TYPE>::clear()
{
	result.clear();
	bDeterministic = true;
	pVars = NULL;
}


extern UTM_Setting<Interval> interval_utm_setting;
extern Multivariate_Polynomial_Setting<Interval> multivariate_polynomial_setting;
extern Expression_Setting<Interval> expression_setting;


inline void exp_taylor_remainder(Interval & result, const Interval & tmRange, const unsigned int order, const Global_Setting & setting)
{
	Interval intProd = tmRange.pow(order);

	Interval J = tmRange;
	J.exp_assign();

	result = setting.factorial_rec[order] * intProd * J;
}

inline void rec_taylor_remainder(Interval & result, const Interval & const_part, const Interval & tmRange, const unsigned int order, const Global_Setting & setting)
{
	Interval original_range = const_part + tmRange;

	Interval dom = original_range;
	Interval nom = tmRange;

	result = nom / dom;
	result.pow_assign(order);

	result /= dom;
/*
	Interval dom = original_range.pow(order + 1);

	Interval nom = tmRange.pow(order);

	result = nom / dom;
*/
	if(order % 2 == 1)
	{
		result *= -1;
	}
}

inline void sin_taylor_remainder(Interval & result, const Interval & C, const Interval & tmRange, const unsigned int order, const Global_Setting & setting)
{
	Interval intProd = tmRange.pow(order);

	Interval J = tmRange;
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

inline void cos_taylor_remainder(Interval & result, const Interval & C, const Interval & tmRange, const unsigned int order, const Global_Setting & setting)
{
	Interval intProd = tmRange.pow(order);

	Interval J = tmRange;
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
	result = tmRange;
	result += 1;
	result.rec_assign();

	result *= tmRange;

	result.pow_assign(order);

	result /= order;

	if((order+1)%2 == 1)		// order+1 is odd
	{
		result *= -1;
	}
}

inline void sqrt_taylor_remainder(Interval & result, const Interval & tmRange, const int order, const Global_Setting & setting)
{
	Interval I = tmRange;
	I += 1;
	I.rec_assign();
	I.sqrt_assign();

	I *= tmRange;
	I /= 2;

	result = I.pow(order);

	result *= setting.double_factorial[2*order-3] * setting.factorial_rec[order];

	if(order % 2 == 0)
	{
		result *= -1;
	}
}



}


#endif /* SETTINGS_H_ */
