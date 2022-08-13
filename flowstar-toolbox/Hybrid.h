/*---
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.
  The code is released as is under the GNU General Public License (GPL).
---*/


#ifndef HYBRID_H_
#define HYBRID_H_

#include "Continuous.h"



/*
 * class for state feedback control systems
 * The dynamics is defined by an ODE x' = f(x,u) such that the control variable(s) u
 * is updated every delta time (including t = 0) by the feedback law u = g(x).
 */

template <class DATA_TYPE>
class Feedback
{
protected:
	ODE<DATA_TYPE> ode;
	std::vector<unsigned int> control_var_IDs;
	double delta;

	std::vector<Expression<DATA_TYPE> > exprs;

public:
	Feedback(const Variables & vars, const double control_stepsize, const std::vector<std::string> & str_ode, const std::vector<std::string> & str_feedback_law);
	~Feedback() {} ;

	bool setFeedbackLaw(const std::vector<std::string> & str_feedback_law);

	void reach(Result_of_Reachability & result, const Flowpipe & initialSet, const double T, Computational_Setting & setting, const std::vector<Constraint> & safeSet) const;
	void reach(Result_of_Reachability & result, const Flowpipe & initialSet, const double T, Computational_Setting & setting, const std::vector<Constraint> & safeSet, Symbolic_Remainder & symbolic_remainder) const;
};



template <class DATA_TYPE>
Feedback<DATA_TYPE>::Feedback(const Variables & vars, const double control_stepsize, const std::vector<std::string> & str_ode, const std::vector<std::string> & str_feedback_law)
{
	ode.stateVars = vars;

	Expression<DATA_TYPE> zero(DATA_TYPE(0));
	ode.expressions.resize(vars.size(), zero);

	int numOfcontrolVars = str_feedback_law.size();

	for(int i=str_ode.size(); i<vars.size(); ++i)
	{
		control_var_IDs.push_back(i);
	}

	for(int i=0; i<str_ode.size(); ++i)
	{
		Expression<DATA_TYPE> expr(str_ode[i], ode.stateVars);
		ode.expressions[i] = expr;
	}

	if(str_feedback_law.size() != control_var_IDs.size())
	{
		printf("Setting Feedback Law: Variable numbers do not match.\n");
	}
	else
	{
		for(int i=0; i<str_feedback_law.size(); ++i)
		{
			Expression<DATA_TYPE> expr(str_feedback_law[i], ode.stateVars);
			exprs.push_back(expr);
		}
	}

	delta = control_stepsize;
}

template <class DATA_TYPE>
bool Feedback<DATA_TYPE>::setFeedbackLaw(const std::vector<std::string> & str_feedback_law)
{
	if(str_feedback_law.size() != control_var_IDs.size())
	{
		printf("Setting Feedback Law: Variable numbers do not match.\n");
		return false;
	}
	else
	{
		for(int i=0; i<str_feedback_law.size(); ++i)
		{
			Expression<DATA_TYPE> expr(str_feedback_law[i], ode.stateVars);
			exprs.push_back(expr);
		}

		return true;
	}
}

template <class DATA_TYPE>
void Feedback<DATA_TYPE>::reach(Result_of_Reachability & result, const Flowpipe & initialSet, const double T, Computational_Setting & setting, const std::vector<Constraint> & safeSet) const
{
	Flowpipe step_initial_set = initialSet;
	double remaining_time = T;
	int k = 1;

	for(double t = THRESHOLD_HIGH; t<T; t+=delta)
	{
		// computing the control input range
		std::vector<TaylorModel<Real> > tmv_state;
		int num_of_state_vars = step_initial_set.tmvPre.tms.size() - control_var_IDs.size();

		for(int i=0; i<num_of_state_vars; ++i)
		{
			tmv_state.push_back(step_initial_set.tmvPre.tms[i]);
		}

		for(int i=0; i<exprs.size(); ++i)
		{
			TaylorModel<Real> tmTemp;
			exprs[i].evaluate(tmTemp, tmv_state, setting.tm_setting.order, step_initial_set.domain, setting.tm_setting.cutoff_threshold, setting.g_setting);

			step_initial_set.tmvPre.tms[control_var_IDs[i]] = tmTemp;
		}

		remaining_time = T - t;

		// computing the flowpipes in one control step
		if(remaining_time >= delta)
		{
			ode.reach(result, step_initial_set, delta, setting, safeSet);

			if(result.status > 3)
			{
				printf("Feedback: Cannot complete the reachable set computation.\n");
				return;
			}
			else if(result.status == COMPLETED_UNSAFE)
			{
				printf("Feedback: The system is unsafe.\n");
				return;
			}
		}
		else
		{
			ode.reach(result, step_initial_set, remaining_time, setting, safeSet);

			if(result.status > 3)
			{
				printf("Feedback: Cannot complete the reachable set computation.\n");
				return;
			}
			else if(result.status == COMPLETED_UNSAFE)
			{
				printf("Feedback: The system is unsafe.\n");
				return;
			}
		}

		step_initial_set = result.fp_end_of_time;

		printf("Step %d: Done.\n", k++);
	}
}

template <class DATA_TYPE>
void Feedback<DATA_TYPE>::reach(Result_of_Reachability & result, const Flowpipe & initialSet, const double T, Computational_Setting & setting, const std::vector<Constraint> & safeSet, Symbolic_Remainder & symbolic_remainder) const
{
	Flowpipe step_initial_set = initialSet;
	double remaining_time = T;
	int k = 1;

	for(double t = THRESHOLD_HIGH; t<T; t+=delta)
	{
		// computing the control input range
		std::vector<TaylorModel<Real> > tmv_state;
		int num_of_state_vars = step_initial_set.tmvPre.tms.size() - control_var_IDs.size();

		for(int i=0; i<num_of_state_vars; ++i)
		{
			tmv_state.push_back(step_initial_set.tmvPre.tms[i]);
		}

		for(int i=0; i<exprs.size(); ++i)
		{
			TaylorModel<Real> tmTemp;
			exprs[i].evaluate(tmTemp, tmv_state, setting.tm_setting.order, step_initial_set.domain, setting.tm_setting.cutoff_threshold, setting.g_setting);

			step_initial_set.tmvPre.tms[control_var_IDs[i]] = tmTemp;
		}

		remaining_time = T - t;

		// computing the flowpipes in one control step
		if(remaining_time >= delta)
		{
			ode.reach(result, step_initial_set, delta, setting, safeSet, symbolic_remainder);

			if(result.status > 3)
			{
				printf("Feedback: Cannot complete the reachable set computation.\n");
				return;
			}
			else if(result.status == COMPLETED_UNSAFE)
			{
				printf("Feedback: The system is unsafe.\n");
				return;
			}
		}
		else
		{
			ode.reach(result, step_initial_set, remaining_time, setting, safeSet, symbolic_remainder);

			if(result.status > 3)
			{
				printf("Feedback: Cannot complete the reachable set computation.\n");
				return;
			}
			else if(result.status == COMPLETED_UNSAFE)
			{
				printf("Feedback: The system is unsafe.\n");
				return;
			}
		}

		step_initial_set = result.fp_end_of_time;

		printf("Step %d: Done.\n", k++);
	}
}







class Initial_Set_Configuration
{
public:
	Flowpipe initialSet;
	double remaining_time;
	int previous_mode;

public:
	Initial_Set_Configuration();
	Initial_Set_Configuration(const Flowpipe & fp, const double t, const int mode);
	~Initial_Set_Configuration();

	Initial_Set_Configuration & operator = (const Initial_Set_Configuration & isc);
};




class Queue_of_Initial_Set
{
protected:
	std::list<Initial_Set_Configuration> iscs;

public:
	Queue_of_Initial_Set();
	~Queue_of_Initial_Set();

	bool isEmpty() const;
	void enqueue(const Initial_Set_Configuration & isc);
	void dequeue(Initial_Set_Configuration & isc);
};



#endif /* HYBRID_H_ */
