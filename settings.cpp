/*---
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.
  The code is released as is under the GNU General Public License (GPL).
---*/

#include "settings.h"
#include "expression.h"

using namespace flowstar;

namespace flowstar
{

UTM_Setting<Interval> interval_utm_setting;
Multivariate_Polynomial_Setting<Interval> multivariate_polynomial_setting;
Expression_AST_Setting<Interval> expression_ast_setting;


Taylor_Model_Computation_Setting::Taylor_Model_Computation_Setting()
{
	order = 0;

	step_min = 0;
	step_max = 0;
	order_min = 0;
	order_max = 0;
	queue_size = 0;
}

Taylor_Model_Computation_Setting::Taylor_Model_Computation_Setting(const Variables & vars)
{
	variables = vars;
	order = 0;

	step_min = 0;
	step_max = 0;
	order_min = 0;
	order_max = 0;
	queue_size = 0;
}

Taylor_Model_Computation_Setting::Taylor_Model_Computation_Setting(const Variables & vars, const Parameters & pars)
{
	variables = vars;
	parameters = pars;
	order = 0;

	step_min = 0;
	step_max = 0;
	order_min = 0;
	order_max = 0;
	queue_size = 0;
}

Taylor_Model_Computation_Setting::Taylor_Model_Computation_Setting(const Taylor_Model_Computation_Setting & setting)
{
	variables			= setting.variables;
	parameters			= setting.parameters;
	order				= setting.order;
	remainder_estimation= setting.remainder_estimation;
	domain				= setting.domain;
	step_exp_table		= setting.step_exp_table;
	step_end_exp_table	= setting.step_end_exp_table;
	cutoff_threshold	= setting.cutoff_threshold;

	step_min			= setting.step_min;
	step_max			= setting.step_max;
	order_min			= setting.order_min;
	order_max			= setting.order_max;

	queue_size			= setting.queue_size;
}

Taylor_Model_Computation_Setting::~Taylor_Model_Computation_Setting()
{
}

Taylor_Model_Computation_Setting & Taylor_Model_Computation_Setting::operator = (const Taylor_Model_Computation_Setting & setting)
{
	if(this == &setting)
		return *this;

	variables			= setting.variables;
	parameters			= setting.parameters;
	order				= setting.order;
	remainder_estimation= setting.remainder_estimation;
	domain				= setting.domain;
	step_exp_table		= setting.step_exp_table;
	step_end_exp_table	= setting.step_end_exp_table;
	cutoff_threshold	= setting.cutoff_threshold;

	step_min			= setting.step_min;
	step_max			= setting.step_max;
	order_min			= setting.order_min;
	order_max			= setting.order_max;

	queue_size			= setting.queue_size;

	return *this;
}

void Taylor_Model_Computation_Setting::initializeAdaptiveSettings(const double delta_min, const double delta_max, const unsigned int k_min, const unsigned int k_max)
{
	step_min = delta_min;
	step_max = delta_max;
	order_min = k_min;
	order_max = k_max;
}

bool Taylor_Model_Computation_Setting::setCutoff(const Interval & cutoff)
{
	if(cutoff.valid())
	{
		cutoff_threshold = cutoff;
		return true;
	}
	else
	{
		return false;
	}
}

bool Taylor_Model_Computation_Setting::setStepsize(const double step, const unsigned int k)
{
	if(step <= 0 || k < 2)
	{
		std::cout << "The stepsize should be a positive number and the order should be an integer > 1." <<std::endl;
		return false;
	}

	step_exp_table.clear();
	step_end_exp_table.clear();

	Interval intStep(0, step);
	Interval tmp = intStep;

	step_exp_table.push_back(1);
	step_exp_table.push_back(intStep);

	step_end_exp_table.push_back(1);
	step_end_exp_table.push_back(step);

	order = k;
	unsigned int n = 2*(k + 1) + 1;

	for(unsigned int i=2; i<=n; ++i)
	{
		tmp *= intStep;
		step_exp_table.push_back(tmp);

		Real c;
		tmp.sup(c);
		step_end_exp_table.push_back(c);
	}

	return true;
}

bool Taylor_Model_Computation_Setting::resetOrder(const double step, const unsigned int k)
{
	if(k < 2)
	{
		std::cout << "The order should be an integer > 1." <<std::endl;
		return false;
	}
	else if(k < order)
	{
		order = k;
		return true;
	}
	else if(order == 0)
	{
		setStepsize(step, k);
		return true;
	}
	else
	{
		unsigned int n = 2*(k + 1) + 1;
		unsigned int currentSize = step_exp_table.size();
		if(currentSize < n)
		{
			Interval step = step_exp_table[1];
			Interval tmp = step_exp_table.back();

			for(unsigned int i = currentSize-1; i<=k; ++i)
			{
				tmp *= step;
				step_exp_table.push_back(tmp);
				step_end_exp_table.push_back(tmp.sup());
			}
		}

		order = k;

		return true;
	}
}

bool Taylor_Model_Computation_Setting::resetOrder(const unsigned int k)
{
	if(k < 2)
	{
		return false;
	}

	if(k <= order)
	{
		order = k;
		return true;
	}
	else
	{
		unsigned int n = 2*(k + 1) + 1;
		unsigned int currentSize = step_exp_table.size();
		if(currentSize < n)
		{
			Interval step = step_exp_table[1];
			Interval tmp = step_exp_table.back();

			for(unsigned int i = currentSize-1; i<=k; ++i)
			{
				tmp *= step;
				step_exp_table.push_back(tmp);
				step_end_exp_table.push_back(tmp.sup());
			}
		}

		order = k;

		return true;
	}
}

void Taylor_Model_Computation_Setting::setRemainderEstimation(const std::vector<Interval> & intVec)
{
	remainder_estimation = intVec;
}

void Taylor_Model_Computation_Setting::setDomain(const std::vector<Interval> & intVec)
{
	domain = intVec;
}

void Taylor_Model_Computation_Setting::setQueueSize(const unsigned int m)
{
	queue_size = m;
}

void Taylor_Model_Computation_Setting::clear()
{
	order = 0;
	step_exp_table.clear();
	step_end_exp_table.clear();
	remainder_estimation.clear();
	domain.clear();
}

bool Taylor_Model_Computation_Setting::bAdaptiveStepSize() const
{
	if(step_min > 0)
		return true;
	else
		return false;
}

bool Taylor_Model_Computation_Setting::bAdaptiveOrder() const
{
	if(order_max > 0)
		return true;
	else
		return false;
}





Global_Computation_Setting::Global_Computation_Setting()
{
}

Global_Computation_Setting::Global_Computation_Setting(const Global_Computation_Setting & setting)
{
	factorial_rec = setting.factorial_rec;
	power_4 = setting.power_4;
	double_factorial = setting.double_factorial;
}

Global_Computation_Setting::~Global_Computation_Setting()
{
}

bool Global_Computation_Setting::resetOrder(const unsigned int k)
{
	unsigned int currentOrder = factorial_rec.size() - 1;

	if(k < 2)
	{
		return false;
	}
	else if(k <= currentOrder)
	{
		return true;
	}

	Real r = factorial_rec.back();

	for(unsigned int i=currentOrder; i<=k; ++i)
	{
		r /= i;
		factorial_rec.push_back(r);
	}

	r = power_4.back();

	for(unsigned int i=currentOrder; i<=k; ++i)
	{
		r *= 4;
		power_4.push_back(r);
	}

	Real odd;
	Real even;

	if(currentOrder % 2 == 0)
	{
		odd = double_factorial[currentOrder-1];
		even = double_factorial[currentOrder];
	}
	else
	{
		odd = double_factorial[currentOrder];
		even = double_factorial[currentOrder-1];
	}

	for(unsigned int i=currentOrder; i<=k; ++i)
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

	return true;
}

Global_Computation_Setting & Global_Computation_Setting::operator = (const Global_Computation_Setting & setting)
{
	if(this == &setting)
		return *this;

	factorial_rec = setting.factorial_rec;
	power_4 = setting.power_4;
	double_factorial = setting.double_factorial;

	return *this;
}

void Global_Computation_Setting::prepareForReachability(const unsigned int maxOrder)
{
	compute_factorial_rec(factorial_rec, maxOrder+5);
	compute_power_4(power_4, maxOrder+5);
	compute_double_factorial(double_factorial, 2*maxOrder+10);
}



}




