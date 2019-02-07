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
class Expression_AST;


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
	std::vector<DATA_TYPE> val_exp_table;

public:
	UTM_Setting();
	UTM_Setting(const UTM_Setting<DATA_TYPE> & setting);
	~UTM_Setting();

	bool setValue(const DATA_TYPE & value, const unsigned int k);
	bool setOrder(const unsigned int k);
	bool resetOrder(const DATA_TYPE & value, const unsigned int k);

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
	val_exp_table = setting.val_exp_table;
}

template <class DATA_TYPE>
UTM_Setting<DATA_TYPE>::~UTM_Setting()
{
}

template <class DATA_TYPE>
bool UTM_Setting<DATA_TYPE>::setValue(const DATA_TYPE & value, const unsigned int k)
{
	if(k == 0)
	{
		std::cout << "The order should be a positive integer." <<std::endl;
		return false;
	}
	else
	{
		val_exp_table.clear();
		order = k;

		val_exp_table.push_back(1);
		val_exp_table.push_back(value);

		DATA_TYPE tmp = value;
		unsigned int nec_order = 2*k + 1;

		for(unsigned int i=2; i<=nec_order; ++i)
		{
			tmp *= value;
			val_exp_table.push_back(tmp);
		}

		return true;
	}
}

template <class DATA_TYPE>
bool UTM_Setting<DATA_TYPE>::resetOrder(const DATA_TYPE & value, const unsigned int k)
{
	if(k == 0)
	{
		std::cout << "The order should be a positive integer." <<std::endl;
		return false;
	}
	else if(k < order)
	{
		order = k;
		return true;
	}
	else if(order == 0)
	{
		setValue(value, k);
		return true;
	}
	else
	{
		unsigned int diff = k - order;
		DATA_TYPE tmp;

		tmp = val_exp_table.back();

		for(unsigned int i=0; i<diff; ++i)
		{
			tmp *= value;
			val_exp_table.push_back(tmp);
		}

		order = k;

		return true;
	}
}

template <class DATA_TYPE>
bool UTM_Setting<DATA_TYPE>::setOrder(const unsigned int k)
{
	if(k == 0)
	{
		std::cout << "The order should be a positive integer." <<std::endl;
		return false;
	}
	else
	{
		order = k;
		return true;
	}
}

template <class DATA_TYPE>
UTM_Setting<DATA_TYPE> & UTM_Setting<DATA_TYPE>::operator = (const UTM_Setting<DATA_TYPE> & setting)
{
	if(this == &setting)
		return *this;

	order = setting.order;
	val_exp_table = setting.val_exp_table;

	return *this;
}



class Taylor_Model_Computation_Setting
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

	unsigned int queue_size;

public:
	Taylor_Model_Computation_Setting();
	Taylor_Model_Computation_Setting(const Variables & vars);
	Taylor_Model_Computation_Setting(const Variables & vars, const Parameters & pars);
	Taylor_Model_Computation_Setting(const Taylor_Model_Computation_Setting & setting);
	~Taylor_Model_Computation_Setting();

	Taylor_Model_Computation_Setting & operator = (const Taylor_Model_Computation_Setting & setting);

	void initializeAdaptiveSettings(const double delta_min, const double delta_max, const unsigned int k_min, const unsigned int k_max);

	bool setCutoff(const Interval & cutoff);
	bool setStepsize(const double step, const unsigned int k);

	bool resetOrder(const double step, const unsigned int k);
	bool resetOrder(const unsigned int k);

	void setRemainderEstimation(const std::vector<Interval> & intVec);
	void setDomain(const std::vector<Interval> & intVec);
	void setQueueSize(const unsigned int m);

	void clear();

	bool bAdaptiveStepSize() const;
	bool bAdaptiveOrder() const;

	friend class Flowpipe;
	friend class Continuous_Reachability;
};



class Global_Computation_Setting
{
public:
	std::vector<Real> factorial_rec;
	std::vector<Real> power_4;
	std::vector<Real> double_factorial;

public:
	Global_Computation_Setting();
	Global_Computation_Setting(const Global_Computation_Setting & setting);
	~Global_Computation_Setting();

	bool resetOrder(const unsigned int k);

	Global_Computation_Setting & operator = (const Global_Computation_Setting & setting);

	void prepareForReachability(const unsigned int maxOrder);
};



template <class DATA_TYPE>
class Multivariate_Polynomial_Setting
{
public:
	std::string strPolynomial;
	Polynomial<DATA_TYPE> result;
	bool bDeterministic;

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
}

template <class DATA_TYPE>
Multivariate_Polynomial_Setting<DATA_TYPE>::Multivariate_Polynomial_Setting(const Multivariate_Polynomial_Setting<DATA_TYPE> & setting)
{
	strPolynomial = setting.strPolynomial;
	result = setting.result;
	bDeterministic = setting.bDeterministic;
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

	return *this;
}

template <class DATA_TYPE>
void Multivariate_Polynomial_Setting<DATA_TYPE>::clear()
{
	result.clear();
	bDeterministic = true;
}





template <class DATA_TYPE>
class Expression_AST_Setting
{
public:
	std::string strExpression;
	Expression_AST<DATA_TYPE> result;
	bool bDeterministic;

public:
	Expression_AST_Setting();
	Expression_AST_Setting(const Expression_AST_Setting & setting);
	~Expression_AST_Setting();

	Expression_AST_Setting & operator = (const Expression_AST_Setting & setting);

	void clear();
};

template <class DATA_TYPE>
Expression_AST_Setting<DATA_TYPE>::Expression_AST_Setting()
{
	bDeterministic = true;
}

template <class DATA_TYPE>
Expression_AST_Setting<DATA_TYPE>::Expression_AST_Setting(const Expression_AST_Setting<DATA_TYPE> & setting)
{
	strExpression = setting.strExpression;
	result = setting.result;
	bDeterministic = setting.bDeterministic;
}

template <class DATA_TYPE>
Expression_AST_Setting<DATA_TYPE>::~Expression_AST_Setting()
{
}

template <class DATA_TYPE>
Expression_AST_Setting<DATA_TYPE> & Expression_AST_Setting<DATA_TYPE>::operator = (const Expression_AST_Setting<DATA_TYPE> & setting)
{
	if(this == &setting)
		return *this;

	strExpression = setting.strExpression;
	result = setting.result;
	bDeterministic = setting.bDeterministic;

	return *this;
}

template <class DATA_TYPE>
void Expression_AST_Setting<DATA_TYPE>::clear()
{
	result.clear();
	bDeterministic = true;
}


extern UTM_Setting<Interval> interval_utm_setting;
extern Multivariate_Polynomial_Setting<Interval> multivariate_polynomial_setting;
extern Expression_AST_Setting<Interval> expression_ast_setting;

}


#endif /* SETTINGS_H_ */
