/*---
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.
  The code is released as is under the GNU General Public License (GPL).
---*/

#ifndef CONTINUOUS_H_
#define CONTINUOUS_H_

#include "Geometry.h"

namespace flowstar
{

class Flowpipe;

class Symbolic_Remainder
{
public:
	std::vector<Matrix<Interval> > J;
	std::vector<Matrix<Real> > Phi_L;
	std::vector<Real> scalars;
	std::vector<Polynomial<Real> > polynomial_of_initial_set;
	unsigned int max_size;

public:
	Symbolic_Remainder();
	Symbolic_Remainder(const Flowpipe & initial_set, const unsigned int s);
	Symbolic_Remainder(const Symbolic_Remainder & symbolic_remainder);
	~Symbolic_Remainder();

	void reset(const Flowpipe & initial_set);

	Symbolic_Remainder & operator = (const Symbolic_Remainder & symbolic_remainder);
};



class Computational_Setting
{
public:
	Taylor_Model_Setting tm_setting;
	Global_Setting g_setting;
	Symbolic_Remainder symbolic_remainder;
	double time;
	bool bPrint;

public:
	Computational_Setting();
	Computational_Setting(const Computational_Setting & setting);
	~Computational_Setting();

	void clear();

	bool setTime(const double t);
	bool setFixedStepsize(const double step, const unsigned int order);
	bool setFixedStepsize(const double step, const unsigned int order_min, const unsigned int order_max);
	bool setAdaptiveStepsize(const double step_min, const double step_max, const unsigned int order);

	bool setCutoffThreshold(const double threshold);
	void setRemainderEstimation(const std::vector<Interval> & estimation);

	bool resetOrder(const unsigned int order);
	bool resetOrder(const unsigned int order_min, const unsigned int order_max);

	void printOn();
	void printOff();

	void prepare();

	Computational_Setting & operator = (const Computational_Setting & setting);
};




class Flowpipe					// A flowpipe is represented by a composition of two Taylor models. The left Taylor model is the preconditioning part.
{
public:
	TaylorModelVec<Real> tmvPre;
	TaylorModelVec<Real> tmv;

	std::vector<Interval> domain;	// domain of TMV_right, the first variable is t

public:
	Flowpipe();
	Flowpipe(const std::vector<Interval> & box);
	Flowpipe(const TaylorModelVec<Real> & tmv_flowpipe, const std::vector<Interval> & flowpipe_domain);
	Flowpipe(const TaylorModelVec<Real> & tmv_flowpipe, const std::vector<Interval> & flowpipe_domain, const Interval & cutoff_threshold);
	Flowpipe(const Flowpipe & flowpipe);
	~Flowpipe();

	void clear();

	void compose(TaylorModelVec<Real> & result, const unsigned int order, const Interval & cutoff_threshold) const;
	void compose(TaylorModelVec<Real> & result, const std::vector<unsigned int> & orders, const Interval & cutoff_threshold) const;
	void compose(TaylorModelVec<Real> & result, const std::vector<unsigned int> & outputAxes, const unsigned int order, const Interval & cutoff_threshold) const;

	void compose_normal(TaylorModelVec<Real> & result, const std::vector<Interval> & step_exp_table, const unsigned int order, const Interval & cutoff_threshold) const;
	void compose_normal(TaylorModelVec<Real> & result, const std::vector<Interval> & step_exp_table, const std::vector<unsigned int> & orders, const Interval & cutoff_threshold) const;
	void compose_normal(TaylorModelVec<Real> & result, const std::vector<unsigned int> & outputAxes, const unsigned int order, const std::vector<Interval> & step_exp_table, const Interval & cutoff_threshold) const;

	void intEval(std::vector<Interval> & result, const unsigned int order, const Interval & cutoff_threshold) const;
	void intEvalNormal(std::vector<Interval> & result, const std::vector<Interval> & step_exp_table, const unsigned int order, const Interval & cutoff_threshold) const;

	void normalize(const Interval & cutoff_threshold);

	int safetyChecking(const std::vector<Constraint> & unsafeSet, const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting) const;

	bool isInTarget(const std::vector<Constraint> & targetSet, const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting) const;
	bool isInTarget(const std::vector<Constraint> & targetSet, const Computational_Setting & setting) const;

	Flowpipe & operator = (const Flowpipe & flowpipe);


	// interval remainders
	// fixed step sizes and orders
	int advance_deterministic(Flowpipe & result, const std::vector<Expression<Real> > & ode, const Taylor_Model_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Setting & g_setting) const;

	int advance_nondeterministic(Flowpipe & result, const std::vector<Expression<Interval> > & ode, const Taylor_Model_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Setting & g_setting) const;


	// adaptive step sizes and fixed orders
	int advance_deterministic_adaptive_stepsize(Flowpipe & result, const std::vector<Expression<Real> > & ode, const double new_stepsize, Taylor_Model_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Setting & g_setting) const;

	int advance_nondeterministic_adaptive_stepsize(Flowpipe & result, const std::vector<Expression<Interval> > & ode, const double new_stepsize, Taylor_Model_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Setting & g_setting) const;


	// fixed step sizes and adaptive orders
	int advance_deterministic_adaptive_order(Flowpipe & result, const std::vector<Expression<Real> > & ode, Taylor_Model_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Setting & g_setting) const;

	int advance_nondeterministic_adaptive_order(Flowpipe & result, const std::vector<Expression<Interval> > & ode, Taylor_Model_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Setting & g_setting) const;


	// symbolic remainders
	// fixed step sizes and orders
	int advance_deterministic(Flowpipe & result, const std::vector<Expression<Real> > & ode, const Taylor_Model_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Setting & g_setting, Symbolic_Remainder & symbolic_remainder) const;

	int advance_nondeterministic(Flowpipe & result, const std::vector<Expression<Interval> > & ode, const Taylor_Model_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Setting & g_setting, Symbolic_Remainder & symbolic_remainder) const;


	// adaptive step sizes and fixed orders
	int advance_deterministic_adaptive_stepsize(Flowpipe & result, const std::vector<Expression<Real> > & ode, const double new_stepsize, Taylor_Model_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Setting & g_setting, Symbolic_Remainder & symbolic_remainder) const;

	int advance_nondeterministic_adaptive_stepsize(Flowpipe & result, const std::vector<Expression<Interval> > & ode, const double new_stepsize, Taylor_Model_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Setting & g_setting, Symbolic_Remainder & symbolic_remainder) const;


	// fixed step sizes and adaptive orders
	int advance_deterministic_adaptive_order(Flowpipe & result, const std::vector<Expression<Real> > & ode, Taylor_Model_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Setting & g_setting, Symbolic_Remainder & symbolic_remainder) const;

	int advance_nondeterministic_adaptive_order(Flowpipe & result, const std::vector<Expression<Interval> > & ode, Taylor_Model_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Setting & g_setting, Symbolic_Remainder & symbolic_remainder) const;
};




class LinearFlowpipe
{
public:
	Matrix<UnivariateTaylorModel<Real> > Phi;
	Matrix<UnivariateTaylorModel<Real> > Psi;

	Zonotope tv_remainder;

	Matrix<Interval> remainder_constraints;

public:
	LinearFlowpipe();
	LinearFlowpipe(const LinearFlowpipe & flowpipe);
	~LinearFlowpipe();

	int safetyChecking(const std::vector<Constraint> & unsafeSet, const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting,
			const TaylorModelVec<Real> & tmv_of_X0, const std::vector<Interval> & polyRangeX0, const std::vector<Interval> & range_of_X0, const std::vector<Interval> & domain);

	void evaluate(TaylorModelVec<Real> & result, const TaylorModelVec<Real> & tmv_of_X0, const std::vector<Interval> & polyRangeX0, const std::vector<Interval> & range_of_X0, const std::vector<Interval> & domain, const Taylor_Model_Setting & tm_setting);
	void evaluate(TaylorModelVec<Real> & result, const TaylorModelVec<Real> & initialSet, const std::vector<Interval> & domain, const unsigned int order, const Interval & cutoff_threshold);
	void evaluate(TaylorModelVec<Real> & result, const std::vector<unsigned int> & outputAxes, const TaylorModelVec<Real> & tmv_of_X0, const std::vector<Interval> & polyRangeX0, const std::vector<Interval> & range_of_X0, const std::vector<Interval> & domain, const Taylor_Model_Setting & tm_setting);

	LinearFlowpipe & operator = (const LinearFlowpipe & flowpipe);

	friend class Flowpipe;
	friend class Linear_Time_Invariant_Dynamics;
	friend class Linear_Time_Varying_Dynamics;
};





class Result_of_Reachability
{
public:
	int status;
	unsigned long num_of_flowpipes;

	Flowpipe fp_end_of_time;
	TaylorModelVec<Real> tmv_end_of_time;

	LinearFlowpipe lfp_end_of_time;		// the last linear flowpipe
	Flowpipe contracted_initial_set;	// the contracted initial set for the last linear flowpipe if there is an invariant

	std::list<LinearFlowpipe> linear_flowpipes;
	std::list<Flowpipe> nonlinear_flowpipes;

	std::list<TaylorModelVec<Real> > tmv_flowpipes;
	std::list<std::vector<Interval> > tmv_flowpipes_domains;

	std::list<unsigned int> orders_of_flowpipes;
	std::list<int> safety_of_flowpipes;
	std::list<bool> contraction_of_flowpipes;

public:
	Result_of_Reachability();
	Result_of_Reachability(const Result_of_Reachability & result);
	~Result_of_Reachability();

	void clear();
	void merge(const Result_of_Reachability & result);

	void transformToTaylorModels(const Taylor_Model_Setting & tm_setting, const bool bPrint);
	void transformToTaylorModels(const Computational_Setting & c_setting);

	void computeBoxOverapproximations(std::list<std::vector<Interval> > & boxes, const Taylor_Model_Setting & tm_setting, const bool bPrint);
	void computeBoxOverapproximations(std::list<std::vector<Interval> > & boxes, const Computational_Setting & c_setting);

	void computeDiscreteBoxOverapproximations(std::list<std::vector<Interval> > & boxes, const Taylor_Model_Setting & tm_setting, const bool bPrint);
	void computeDiscreteBoxOverapproximations(std::list<std::vector<Interval> > & boxes, const Computational_Setting & c_setting);

	void transformToTaylorModels(const Taylor_Model_Setting & tm_setting, const bool bPrint, const Flowpipe & initialSet);
	void transformToTaylorModels(const Computational_Setting & c_setting, const Flowpipe & initialSet);


	Result_of_Reachability & operator = (const Result_of_Reachability & result);
};








// base class for continuous dynamics
class Dynamics
{
public:
	Dynamics()
	{
	}

	virtual ~Dynamics()
	{
	}

	friend class Flowpipe;
};









class Linear_Time_Invariant_Dynamics : public Dynamics
{
protected:
	Matrix<Real> 							rm_dyn_A;
	Matrix<UnivariateTaylorModel<Real> > 	utm_dyn_B;
	Matrix<Real> 							rm_dyn_C;
	Matrix<bool> 							connectivity;

public:
	Linear_Time_Invariant_Dynamics(const Matrix<Real> & A, const Matrix<UnivariateTaylorModel<Real> > & B);
	Linear_Time_Invariant_Dynamics(const Matrix<Real> & A, const Matrix<UnivariateTaylorModel<Real> > & B, const Matrix<Real> & C);
	Linear_Time_Invariant_Dynamics(const Linear_Time_Invariant_Dynamics & dynamics);
	virtual ~Linear_Time_Invariant_Dynamics();

	Linear_Time_Invariant_Dynamics & operator = (const Linear_Time_Invariant_Dynamics & dynamics);


	void evaluate(Interval & result, const unsigned int varID, const TaylorModelVec<Real> & tmv_range, const std::vector<Interval> & domain);
	void evaluate(Interval & result, const std::vector<Expression<Real> > & coeff_of_Lie_deriv, const TaylorModelVec<Real> & tmv_range, const std::vector<Interval> & domain, const unsigned int order, const Computational_Setting & setting);


	int reach(std::list<LinearFlowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const int zono_order, const Flowpipe & initialSet, const Taylor_Model_Setting & tm_setting,
			const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput);

	int reach_one_step(TaylorModelVec<Real> & result, const double stepsize, const TaylorModelVec<Real> & initialSet,
			std::vector<Interval> & domain, const unsigned int order, const Interval & cutoff_threshold);

	int reach(std::list<LinearFlowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const int zono_order, LinearFlowpipe & l_initial_set, const Flowpipe & g_initial_set,
			const Taylor_Model_Setting & tm_setting,
			const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput);


	int reach_inv(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<std::vector<Interval> > & domains, LinearFlowpipe & last_flowpipe, Flowpipe & contracted_initialSet, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const double time, const int zono_order, const Flowpipe & initialSet, const std::vector<Constraint> & invariant,
			const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking, const bool bPlot, const bool bTMOutput);

	int reach_inv(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<std::vector<Interval> > & domains, LinearFlowpipe & last_flowpipe, Flowpipe & contracted_initialSet, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const double time, const int zono_order, LinearFlowpipe & l_initial_set, const Flowpipe & g_initial_set,
			const std::vector<Constraint> & invariant, const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput);


	void reach(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const std::vector<Constraint> & unsafeSet, const int zono_order = -1);
	void reach(Result_of_Reachability & result, Computational_Setting & setting, LinearFlowpipe & l_initial_set, const Flowpipe & g_initial_set, const std::vector<Constraint> & unsafeSet, const int zono_order = -1);

	void reach_inv(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, Flowpipe & contracted_initialSet, const std::vector<Constraint> & invariant, const std::vector<Constraint> & unsafeSet, const int zono_order = -1);
	void reach_inv(Result_of_Reachability & result, Computational_Setting & setting, LinearFlowpipe & l_initial_set, const Flowpipe & g_initial_set, Flowpipe & contracted_initialSet, const std::vector<Constraint> & invariant, const std::vector<Constraint> & unsafeSet, const int zono_order = -1);
};



class Linear_Time_Varying_Dynamics : public Dynamics
{
protected:
	Matrix<UnivariatePolynomial<Real> >		upm_dyn_A;
	Matrix<UnivariatePolynomial<Real> >		upm_dyn_B;
	Matrix<UnivariatePolynomial<Real> >		upm_dyn_tv;
	Matrix<bool> 							connectivity;
	bool									bAuto;

	Matrix<Interval>	 					uncertain_range;

public:
	Linear_Time_Varying_Dynamics(const Matrix<UnivariatePolynomial<Real> > & A, const Matrix<UnivariatePolynomial<Real> > & B, const Matrix<UnivariatePolynomial<Real> > & C);
	Linear_Time_Varying_Dynamics(const Linear_Time_Varying_Dynamics & dynamics);
	virtual ~Linear_Time_Varying_Dynamics();

	Linear_Time_Varying_Dynamics & operator = (const Linear_Time_Varying_Dynamics & dynamics);

	void evaluate(Interval & result, const unsigned int varID, const TaylorModelVec<Real> & tmv_range, const std::vector<Interval> & domain, const Real & t_lb, const unsigned int order, const Interval & cutoff_threshold);
//	void evaluate(Interval & result, const std::vector<Expression<Real> > & coeff_of_Lie_deriv, const TaylorModelVec<Real> & tmv_range, const std::vector<Interval> & domain, const Real & t_lb, const unsigned int order, const Computational_Setting & setting);

	int reach(std::list<LinearFlowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, const int zono_order, const Taylor_Model_Setting & tm_setting,
			const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput);

	int reach_inv(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, const int zono_order, const std::vector<Constraint> & invariant,
			const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking, const bool bPlot, const bool bTMOutput);

	void reach(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const std::vector<Constraint> & unsafeSet, const int zono_order = -1);
	void reach_inv(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const std::vector<Constraint> & invariant, const std::vector<Constraint> & unsafeSet, const int zono_order = -1);
};


class Deterministic_Continuous_Dynamics : public Dynamics
{
protected:
	std::vector<Expression<Real> >	expressions;

public:
	Deterministic_Continuous_Dynamics(const std::vector<Expression<Real> > & dynamics);
	Deterministic_Continuous_Dynamics(const Deterministic_Continuous_Dynamics & dynamics);
	virtual ~Deterministic_Continuous_Dynamics();

	Deterministic_Continuous_Dynamics & operator = (const Deterministic_Continuous_Dynamics & dynamics);

	void evaluate(Interval & result, const unsigned int varID, const TaylorModelVec<Real> & tmv_range, const std::vector<Interval> & domain, const unsigned int order, const Computational_Setting & setting);

	int reach(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, const Taylor_Model_Setting & tm_setting,
			const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;

	int reach_adaptive_stepsize(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, Taylor_Model_Setting & tm_setting,
			const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;

	int reach_adaptive_order(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, Taylor_Model_Setting & tm_setting,
			const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;

	int reach_symbolic_remainder(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, Symbolic_Remainder & symbolic_remainder, const Taylor_Model_Setting & tm_setting,
			const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;

	int reach_symbolic_remainder_adaptive_stepsize(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, Symbolic_Remainder & symbolic_remainder, Taylor_Model_Setting & tm_setting,
			const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;

	int reach_symbolic_remainder_adaptive_order(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, Symbolic_Remainder & symbolic_remainder, Taylor_Model_Setting & tm_setting,
			const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;


	void reach(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const std::vector<Constraint> & unsafeSet) const;
	void reach_sr(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const std::vector<Constraint> & unsafeSet, Symbolic_Remainder & symbolic_remainder) const;


	// reachability computation in a given invariant set
	int reach_inv(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<std::vector<Interval> > & domains, Flowpipe & last_flowpipe, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, const std::vector<Constraint> & invariant,
			const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking, const bool bPlot, const bool bTMOutput) const;

	int reach_inv_adaptive_stepsize(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<std::vector<Interval> > & domains, Flowpipe & last_flowpipe, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, const std::vector<Constraint> & invariant,
			Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking, const bool bPlot, const bool bTMOutput) const;

	int reach_inv_adaptive_order(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<std::vector<Interval> > & domains, Flowpipe & last_flowpipe, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, const std::vector<Constraint> & invariant,
			Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking, const bool bPlot, const bool bTMOutput) const;


	int reach_inv_symbolic_remainder(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<std::vector<Interval> > & domains, Flowpipe & last_flowpipe, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, const std::vector<Constraint> & invariant, Symbolic_Remainder & symbolic_remainder,
			const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking, const bool bPlot, const bool bTMOutput) const;

	int reach_inv_symbolic_remainder_adaptive_stepsize(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<std::vector<Interval> > & domains, Flowpipe & last_flowpipe, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, const std::vector<Constraint> & invariant, Symbolic_Remainder & symbolic_remainder,
			Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking, const bool bPlot, const bool bTMOutput) const;

	int reach_inv_symbolic_remainder_adaptive_order(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<std::vector<Interval> > & domains, Flowpipe & last_flowpipe, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, const std::vector<Constraint> & invariant, Symbolic_Remainder & symbolic_remainder,
			Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking, const bool bPlot, const bool bTMOutput) const;


	void reach_inv(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const std::vector<Constraint> & invariant, const std::vector<Constraint> & unsafeSet) const;
	void reach_inv_sr(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const std::vector<Constraint> & invariant, const std::vector<Constraint> & unsafeSet, Symbolic_Remainder & symbolic_remainder) const;
};



class Nondeterministic_Continuous_Dynamics : public Dynamics
{
protected:
	std::vector<Expression<Interval> >	expressions;

public:
	Nondeterministic_Continuous_Dynamics(const std::vector<Expression<Interval> > & dynamics);
	Nondeterministic_Continuous_Dynamics(const Nondeterministic_Continuous_Dynamics & dynamics);
	virtual ~Nondeterministic_Continuous_Dynamics();

	Nondeterministic_Continuous_Dynamics & operator = (const Nondeterministic_Continuous_Dynamics & dynamics);

	void evaluate(Interval & result, const unsigned int varID, const TaylorModelVec<Real> & tmv_range, const std::vector<Interval> & domain, const unsigned int order, const Computational_Setting & setting);

	int reach(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, const Taylor_Model_Setting & tm_setting,
			const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;

	int reach_adaptive_stepsize(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, Taylor_Model_Setting & tm_setting,
			const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;

	int reach_adaptive_order(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, Taylor_Model_Setting & tm_setting,
			const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;

	int reach_symbolic_remainder(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, const Taylor_Model_Setting & tm_setting, Symbolic_Remainder & symbolic_remainder,
			const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;

	int reach_symbolic_remainder_adaptive_stepsize(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, Taylor_Model_Setting & tm_setting, Symbolic_Remainder & symbolic_remainder,
			const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;

	int reach_symbolic_remainder_adaptive_order(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, Taylor_Model_Setting & tm_setting, Symbolic_Remainder & symbolic_remainder,
			const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;

	void reach(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const std::vector<Constraint> & unsafeSet) const;
	void reach_sr(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const std::vector<Constraint> & unsafeSet, Symbolic_Remainder & symbolic_remainder) const;


	// reachability computation in a given invariant set
	int reach_inv(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, const std::vector<Constraint> & invariant,
			const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking, const bool bPlot, const bool bTMOutput) const;

	int reach_inv_adaptive_stepsize(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, const std::vector<Constraint> & invariant,
			Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking, const bool bPlot, const bool bTMOutput) const;

	int reach_inv_adaptive_order(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, const std::vector<Constraint> & invariant,
			Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking, const bool bPlot, const bool bTMOutput) const;


	int reach_inv_symbolic_remainder(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, const std::vector<Constraint> & invariant, Symbolic_Remainder & symbolic_remainder,
			const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking, const bool bPlot, const bool bTMOutput) const;

	int reach_inv_symbolic_remainder_adaptive_stepsize(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, const std::vector<Constraint> & invariant, Symbolic_Remainder & symbolic_remainder,
			Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking, const bool bPlot, const bool bTMOutput) const;

	int reach_inv_symbolic_remainder_adaptive_order(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, const std::vector<Constraint> & invariant, Symbolic_Remainder & symbolic_remainder,
			Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking, const bool bPlot, const bool bTMOutput) const;



	void reach_inv(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const std::vector<Constraint> & invariant, const std::vector<Constraint> & unsafeSet) const;
	void reach_inv_sr(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const std::vector<Constraint> & invariant, const std::vector<Constraint> & unsafeSet, Symbolic_Remainder & symbolic_remainder) const;
};




class Plot_Setting
{
public:
	Variables variables;
	std::vector<unsigned int> outputDims;
	unsigned int type_of_file;
	unsigned int type_of_object;
	unsigned int num_of_pieces;
	bool bProjected;
	bool bPrint;
	bool bDiscrete;

public:
	Plot_Setting();
	Plot_Setting(const Variables & vars);
	Plot_Setting(const Plot_Setting & setting);
	~Plot_Setting();

	Plot_Setting & operator = (const Plot_Setting & setting);

	void setOutputDims(const std::string & x, const std::string & y);
	void setFileType(const unsigned int type);
	void setObjectType(const unsigned int type);
	void setNumOfPieces(const unsigned int n);

	void printOn();
	void printOff();

	void discreteOutput();
	void continuousOutput();

	void plot_2D(const std::string & path, const std::string & fileName, const Result_of_Reachability & result) const;

	void plot_2D_MATLAB(const std::string & path, const std::string & fileName, const Result_of_Reachability & result) const;
	void plot_2D_interval_MATLAB(const std::string & path, const std::string & fileName, const Result_of_Reachability & result) const;
	void plot_2D_octagon_MATLAB(const std::string & path, const std::string & fileName, const Result_of_Reachability & result) const;
	void plot_2D_grids_MATLAB(const std::string & path, const std::string & fileName, const unsigned int num, const Result_of_Reachability & result) const;

	void plot_2D_GNUPLOT(const std::string & path, const std::string & fileName, const Result_of_Reachability & result) const;
	void plot_2D_interval_GNUPLOT(const std::string & path, const std::string & fileName, const Result_of_Reachability & result) const;
	void plot_2D_octagon_GNUPLOT(const std::string & path, const std::string & fileName, const Result_of_Reachability & result) const;
	void plot_2D_grids_GNUPLOT(const std::string & path, const std::string & fileName, const unsigned int num, const Result_of_Reachability & result) const;
};


// pool for candidate vectors
class CandidateVectorPool
{
public:
	std::vector<std::vector<double> > vectors;
	Matrix<double> weightMat;

public:
	CandidateVectorPool();
	CandidateVectorPool(const std::vector<std::vector<double> > & candidates);
	CandidateVectorPool(const CandidateVectorPool & pool);
	~CandidateVectorPool();

	unsigned int size() const;

	CandidateVectorPool & operator = (const CandidateVectorPool & pool);
};


// for template selection
class FactorTab
{
public:
	int index;
	double factor;
	double intercept;
public:
	FactorTab();
	FactorTab(const int i, const double & f, const double & interc);
	~FactorTab();

	friend bool compareFactor(const FactorTab & a, const FactorTab & b);
	friend bool compareIntercept(const FactorTab & a, const FactorTab & b);
};



/* data structure to keep the information of all flowpipe/guard intersections
 * consecutive flowpipes are kept by groups (vectors)
 */

class Intersected_Flowpipes
{
public:
	std::vector<std::vector<TaylorModelVec<Real> > > flowpipes;		// Taylor model vectors
	std::vector<std::vector<std::vector<Interval> > > domains;		// contrated domains for the Taylor model vectors
	std::vector<Real> start_t;										// time of the first nonempty intersection
	std::vector<Real> durations;									// time duration of the consecutive flowpipes

public:
	Intersected_Flowpipes();
	Intersected_Flowpipes(const Intersected_Flowpipes & intersections);
	~Intersected_Flowpipes();

	Intersected_Flowpipes & operator = (const Intersected_Flowpipes & intersections);
	void clear();

	unsigned int numOfGroups() const;

	void interval_aggregation(Flowpipe & result) const;
	void interval_aggregation(std::vector<Flowpipe> & result) const;

	// all the vectors in the pool should be normalized to length 1
	void parallelotope_aggregation(Flowpipe & result, CandidateVectorPool & pool) const;
	void parallelotope_aggregation(std::vector<Flowpipe> & result, CandidateVectorPool & pool) const;
};


int safetyChecking(const TaylorModelVec<Real> & tmv, const std::vector<Interval> & domain, const std::vector<Constraint> & unsafeSet, const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting);

void gridBox(std::list<std::vector<Interval> > & grids, const std::vector<Interval> & box, const unsigned int num);

int remainder_contraction_int(const std::vector<Interval> & polyRange, std::vector<Interval> & remainders, const std::vector<Constraint> & constraints);
int domain_contraction_int(const TaylorModelVec<Real> & tmv_flowpipe, std::vector<Interval> & domain, const std::vector<Constraint> & constraints, const unsigned int order, const Interval & cutoff_threshold, const Global_Setting & g_setting);

unsigned int findProperOrder(Real & error, const Real & max, const Real & min, const Real & tolerance, const unsigned int start_order);

void check_connectivities(Matrix<bool> & result, Matrix<bool> & adjMatrix);

void compute_one_step_trans(Matrix<UnivariateTaylorModel<Real> > & utm_Phi_t, Matrix<Real> & rm_Phi_t, Matrix<UnivariateTaylorModel<Real> > & utm_Psi_t, Matrix<Real> & rm_Psi_t,
		Matrix<Interval> & tv_part, const Matrix<UnivariatePolynomial<Real> > & A_t, const Matrix<UnivariatePolynomial<Real> > & B_t, const Matrix<UnivariatePolynomial<Real> > & tv_t,
		Matrix<bool> & connectivity, const bool bAuto, const UnivariatePolynomial<Real> & up_t, const unsigned int order, std::vector<Real> & step_end_exp_table);


// flowpipe/guard intersections
void intersect_a_guard(Intersected_Flowpipes & result, const Result_of_Reachability & flowpipes, const std::vector<Constraint> & guard, const bool boundary_of_invariant, const Computational_Setting & setting);


void compute_weight_matrix(Matrix<double> & weightMat, const std::vector<std::vector<double> > & candidate_vectors);
bool check_validity(Matrix<double> & matTemplate, const std::vector<double> & vec, const int rank);
bool select_a_vector(FactorTab & lst_selected, std::list<FactorTab> & lst_unselected, Matrix<double> & matTemplate, const std::vector<std::vector<double> > & candidate_vectors, int & rank);



void eliminate_t(TaylorModelVec<Real> & tmv_flowpipe, std::vector<Interval> & domain);

void create_initial_set(Flowpipe & initial_set, TaylorModelVec<Real> & tmv_flowpipe, std::vector<Interval> & domain, const unsigned int order, const Interval & cutoff_threshold);

void merge_consecutive_flowpipes(Flowpipe & result, std::vector<TaylorModelVec<Real> > & flowpipes, std::vector<std::vector<Interval> > & flowpipe_domains, Linear_Time_Invariant_Dynamics & LTI_dynamics,
		const std::vector<Constraint> & invariant, const unsigned int order, const Computational_Setting & setting);





// only for testing
void test_domain_contraction(Result_of_Reachability & contraction_result, Result_of_Reachability & reachability_result, const std::vector<Constraint> & constraints, const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting);

}

#endif /* CONTINUOUS_H_ */
