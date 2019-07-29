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

public:
	Symbolic_Remainder();
	Symbolic_Remainder(const Flowpipe & initial_set);
	Symbolic_Remainder(const Symbolic_Remainder & symbolic_remainder);
	~Symbolic_Remainder();

	void reset(const Flowpipe & initial_set);

	Symbolic_Remainder & operator = (const Symbolic_Remainder & symbolic_remainder);
};



class Computational_Setting
{
public:
	Taylor_Model_Computation_Setting tm_setting;
	Global_Computation_Setting g_setting;
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
	void setQueueSize(const unsigned int m);

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

	int safetyChecking(const std::vector<Constraint> & unsafeSet, const Taylor_Model_Computation_Setting & tm_setting, const Global_Computation_Setting & g_setting) const;

	Flowpipe & operator = (const Flowpipe & flowpipe);


	// interval remainders
	// fixed step sizes and orders
	int advance_deterministic(Flowpipe & result, const std::vector<Expression_AST<Real> > & ode, const Taylor_Model_Computation_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Computation_Setting & g_setting) const;

	int advance_nondeterministic(Flowpipe & result, const std::vector<Expression_AST<Interval> > & ode, const Taylor_Model_Computation_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Computation_Setting & g_setting) const;


	// adaptive step sizes and fixed orders
	int advance_deterministic_adaptive_stepsize(Flowpipe & result, const std::vector<Expression_AST<Real> > & ode, Taylor_Model_Computation_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Computation_Setting & g_setting) const;

	int advance_nondeterministic_adaptive_stepsize(Flowpipe & result, const std::vector<Expression_AST<Interval> > & ode, Taylor_Model_Computation_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Computation_Setting & g_setting) const;


	// fixed step sizes and adaptive orders
	int advance_deterministic_adaptive_order(Flowpipe & result, const std::vector<Expression_AST<Real> > & ode, Taylor_Model_Computation_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Computation_Setting & g_setting) const;

	int advance_nondeterministic_adaptive_order(Flowpipe & result, const std::vector<Expression_AST<Interval> > & ode, Taylor_Model_Computation_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Computation_Setting & g_setting) const;


	// symbolic remainders
	// fixed step sizes and orders
	int advance_deterministic(Flowpipe & result, const std::vector<Expression_AST<Real> > & ode, const Taylor_Model_Computation_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Computation_Setting & g_setting, Symbolic_Remainder & symbolic_remainder) const;

	int advance_nondeterministic(Flowpipe & result, const std::vector<Expression_AST<Interval> > & ode, const Taylor_Model_Computation_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Computation_Setting & g_setting, Symbolic_Remainder & symbolic_remainder) const;


	// adaptive step sizes and fixed orders
	int advance_deterministic_adaptive_stepsize(Flowpipe & result, const std::vector<Expression_AST<Real> > & ode, Taylor_Model_Computation_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Computation_Setting & g_setting, Symbolic_Remainder & symbolic_remainder) const;

	int advance_nondeterministic_adaptive_stepsize(Flowpipe & result, const std::vector<Expression_AST<Interval> > & ode, Taylor_Model_Computation_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Computation_Setting & g_setting, Symbolic_Remainder & symbolic_remainder) const;


	// fixed step sizes and adaptive orders
	int advance_deterministic_adaptive_order(Flowpipe & result, const std::vector<Expression_AST<Real> > & ode, Taylor_Model_Computation_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Computation_Setting & g_setting, Symbolic_Remainder & symbolic_remainder) const;

	int advance_nondeterministic_adaptive_order(Flowpipe & result, const std::vector<Expression_AST<Interval> > & ode, Taylor_Model_Computation_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Computation_Setting & g_setting, Symbolic_Remainder & symbolic_remainder) const;
};




class LinearFlowpipe
{
protected:
	Matrix<UnivariateTaylorModel<Real> > Phi;
	Matrix<UnivariateTaylorModel<Real> > Psi;

	Zonotope tv_remainder;

public:
	LinearFlowpipe();
	LinearFlowpipe(const LinearFlowpipe & flowpipe);
	~LinearFlowpipe();

	int safetyChecking(const std::vector<Constraint> & unsafeSet, const Taylor_Model_Computation_Setting & tm_setting, const Global_Computation_Setting & g_setting,
			const TaylorModelVec<Real> & tmv_of_X0, const std::vector<Interval> & polyRangeX0, const std::vector<Interval> & range_of_X0, const std::vector<Interval> & domain);

	void evaluate(TaylorModelVec<Real> & result, const TaylorModelVec<Real> & tmv_of_X0, const std::vector<Interval> & polyRangeX0, const std::vector<Interval> & range_of_X0, const std::vector<Interval> & domain, const Taylor_Model_Computation_Setting & tm_setting);
	void evaluate(TaylorModelVec<Real> & result, const std::vector<unsigned int> & outputAxes, const TaylorModelVec<Real> & tmv_of_X0, const std::vector<Interval> & polyRangeX0, const std::vector<Interval> & range_of_X0, const std::vector<Interval> & domain, const Taylor_Model_Computation_Setting & tm_setting);

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

	std::list<Flowpipe> nonlinear_flowpipes;
	std::list<TaylorModelVec<Real> > tmv_flowpipes;
	std::list<unsigned int> orders_of_flowpipes;
	std::list<int> safety_of_flowpipes;

public:
	Result_of_Reachability();
	Result_of_Reachability(const Result_of_Reachability & result);
	~Result_of_Reachability();

	void clear();
	void transformToTaylorModels(const Taylor_Model_Computation_Setting & tm_setting, const bool bPrint);
	void transformToTaylorModels(const Computational_Setting & c_setting);

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

	virtual int reach_LTI(std::list<LinearFlowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, const Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) = 0;

	virtual int reach_LTV(std::list<LinearFlowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, const Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) = 0;

	virtual int reach(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, const Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const = 0;

	virtual int reach_adaptive_stepsize(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const = 0;

	virtual int reach_adaptive_order(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const = 0;

	virtual int reach_symbolic_remainder(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, const Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const = 0;

	virtual int reach_symbolic_remainder_adaptive_stepsize(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const = 0;

	virtual int reach_symbolic_remainder_adaptive_order(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const = 0;

	friend class Flowpipe;
};









class Linear_Time_Invariant_Dynamics : public Dynamics
{
protected:
	Matrix<Real> 							rm_dyn_A;
	Matrix<UnivariateTaylorModel<Real> > 	utm_dyn_B;
	Matrix<bool> 							connectivity;
	bool									bAuto;

public:
	Linear_Time_Invariant_Dynamics(const Matrix<Real> & A, const Matrix<UnivariateTaylorModel<Real> > & B);
	Linear_Time_Invariant_Dynamics(const Linear_Time_Invariant_Dynamics & dynamics);
	virtual ~Linear_Time_Invariant_Dynamics();

	Linear_Time_Invariant_Dynamics & operator = (const Linear_Time_Invariant_Dynamics & dynamics);

	virtual int reach_LTI(std::list<LinearFlowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, const Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput);

	virtual int reach_LTV(std::list<LinearFlowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, const Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput);

	virtual int reach(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, const Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;

	virtual int reach_adaptive_stepsize(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;

	virtual int reach_adaptive_order(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;

	virtual int reach_symbolic_remainder(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, const Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;

	virtual int reach_symbolic_remainder_adaptive_stepsize(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;

	virtual int reach_symbolic_remainder_adaptive_order(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;
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

	virtual int reach_LTI(std::list<LinearFlowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, const Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput);

	virtual int reach_LTV(std::list<LinearFlowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, const Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput);

	virtual int reach(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, const Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;

	virtual int reach_adaptive_stepsize(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;

	virtual int reach_adaptive_order(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;

	virtual int reach_symbolic_remainder(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, const Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;

	virtual int reach_symbolic_remainder_adaptive_stepsize(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;

	virtual int reach_symbolic_remainder_adaptive_order(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;
};


class Deterministic_Continuous_Dynamics : public Dynamics
{
protected:
	std::vector<Expression_AST<Real> >	expressions;

public:
	Deterministic_Continuous_Dynamics(const std::vector<Expression_AST<Real> > & dynamics);
	Deterministic_Continuous_Dynamics(const Deterministic_Continuous_Dynamics & dynamics);
	virtual ~Deterministic_Continuous_Dynamics();

	Deterministic_Continuous_Dynamics & operator = (const Deterministic_Continuous_Dynamics & dynamics);

	virtual int reach_LTI(std::list<LinearFlowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, const Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput);

	virtual int reach_LTV(std::list<LinearFlowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, const Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput);

	virtual int reach(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, const Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;

	virtual int reach_adaptive_stepsize(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;

	virtual int reach_adaptive_order(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;

	virtual int reach_symbolic_remainder(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, const Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;

	virtual int reach_symbolic_remainder_adaptive_stepsize(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;

	virtual int reach_symbolic_remainder_adaptive_order(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;

	void reach(Result_of_Reachability & result, Computational_Setting & setting, const std::vector<Flowpipe> & initialSets, const std::vector<Constraint> & unsafeSet) const;

	void reach(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const std::vector<Constraint> & unsafeSet) const;
};



class Nondeterministic_Continuous_Dynamics : public Dynamics
{
protected:
	std::vector<Expression_AST<Interval> >	expressions;

public:
	Nondeterministic_Continuous_Dynamics(const std::vector<Expression_AST<Interval> > & dynamics);
	Nondeterministic_Continuous_Dynamics(const Nondeterministic_Continuous_Dynamics & dynamics);
	virtual ~Nondeterministic_Continuous_Dynamics();

	Nondeterministic_Continuous_Dynamics & operator = (const Nondeterministic_Continuous_Dynamics & dynamics);

	virtual int reach_LTI(std::list<LinearFlowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, const Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput);

	virtual int reach_LTV(std::list<LinearFlowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, const Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput);

	virtual int reach(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, const Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;

	virtual int reach_adaptive_stepsize(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;

	virtual int reach_adaptive_order(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;

	virtual int reach_symbolic_remainder(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, const Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;

	virtual int reach_symbolic_remainder_adaptive_stepsize(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;

	virtual int reach_symbolic_remainder_adaptive_order(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, Taylor_Model_Computation_Setting & tm_setting,
			const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;

	void reach(Result_of_Reachability & result, Computational_Setting & setting, const std::vector<Flowpipe> & initialSets, const std::vector<Constraint> & unsafeSet) const;

	void reach(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const std::vector<Constraint> & unsafeSet) const;
};




class Plot_Setting
{
public:
	std::vector<unsigned int> outputDims;
	unsigned int type_of_file;
	unsigned int type_of_object;
	unsigned int num_of_pieces;
	bool bProjected;
	bool bPrint;

public:
	Plot_Setting();
	Plot_Setting(const Plot_Setting & setting);
	~Plot_Setting();

	Plot_Setting & operator = (const Plot_Setting & setting);

	void setOutputDims(const unsigned int x, const unsigned int y);
	void setOutputDims(const std::vector<unsigned int> & dims);
	void setFileType(const unsigned int type);
	void setObjectType(const unsigned int type);
	void setNumOfPieces(const unsigned int n);

	void printOn();
	void printOff();

	void plot_2D(const std::string & fileName, const Result_of_Reachability & result) const;

	void plot_2D_MATLAB(const std::string & fileName, const Result_of_Reachability & result) const;
	void plot_2D_interval_MATLAB(const std::string & fileName, const Result_of_Reachability & result) const;
	void plot_2D_octagon_MATLAB(const std::string & fileName, const Result_of_Reachability & result) const;
	void plot_2D_grids_MATLAB(const std::string & fileName, const unsigned int num, const Result_of_Reachability & result) const;

	void plot_2D_GNUPLOT(const std::string & fileName, const Result_of_Reachability & result) const;
	void plot_2D_interval_GNUPLOT(const std::string & fileName, const Result_of_Reachability & result) const;
	void plot_2D_octagon_GNUPLOT(const std::string & fileName, const Result_of_Reachability & result) const;
	void plot_2D_grids_GNUPLOT(const std::string & fileName, const unsigned int num, const Result_of_Reachability & result) const;
};






class Continuous_Reachability_Problem_Description
{
public:
	Variables stateVars;						// state variables
	Variables tmVars;							// Taylor model variables
	double step_max;							// stepsize
	double step_min;							// minimum stepsize when adptive stepsize is in use
	unsigned int order_min;						// Taylor model order
	unsigned int order_max;						// highest order when adptive order is in use

	double time;								// time horizon

	std::vector<Interval> remainder_estimation;	// remainder estimation for varying time step
	Interval cutoff_threshold;					// cutoff_threshold

	bool bSafetyChecking;
	bool bPlot;
	bool bTMOutput;
	bool bPrint;

	bool bSymbolicRemainder;
	unsigned int queue_size;

	std::vector<Constraint>		unsafeSet;

	int type_of_dynamics;

	std::vector<Expression_AST<Real> >		deterministic_dynamics;
	std::vector<Expression_AST<Interval> >	nondeterministic_dynamics;

	Matrix<Real> 							rm_dyn_A;
	Matrix<UnivariateTaylorModel<Real> > 	utm_dyn_B;

	Matrix<UnivariatePolynomial<Real> >		upm_dyn_A;
	Matrix<UnivariatePolynomial<Real> >		upm_dyn_B;
	Matrix<UnivariatePolynomial<Real> >		upm_dyn_tv;
	Matrix<Interval> 						im_uncertain_range;

	std::vector<Flowpipe> 					initialSets;

	Plot_Setting							plot_setting;
	bool									bDeterministic;
	std::string								fileName;

public:
	Continuous_Reachability_Problem_Description();
	Continuous_Reachability_Problem_Description(const Continuous_Reachability_Problem_Description & description);
	~Continuous_Reachability_Problem_Description();

	Continuous_Reachability_Problem_Description & operator = (const Continuous_Reachability_Problem_Description & description);

	void setStateVars(const Variables & vars);
	void setTMVars(const Variables & vars);

	bool setTimeHorizon(const double t);

	bool setFixedStepsize(const double delta);
	bool setAdaptiveStepsize(const double delta_min, const double delta_max);

	bool setFixedOrder(const unsigned int k);
	bool setAdaptiveOrder(const unsigned int k_min, const unsigned int k_max);

	bool setRemainderEstimation(const std::vector<Interval> & intVec);
	bool setCutoff(const Interval & cutoff);
	bool setPrecision(const unsigned int prec);

	void printOn();
	void printOff();

	void safetyCheckingOn();
	void safetyCheckingOff();

	void setOutputDims(const unsigned int x, const unsigned int y);
	void setFileType(const unsigned int type);
	void setObjectType(const unsigned int type);
	void setNumOfPieces(const unsigned int n);

	void plotOn();
	void plotOff();

	void tmOutputOn();
	void tmOutputOff();

	void setUnsafe(const std::vector<Constraint> & unsafe_constraints);
	void setInitialSets(const std::vector<Flowpipe> & flowpipes);

	void setFileName(const std::string & str);
};





class Continuous_Reachability
{
public:
	int									type_of_dynamics;

	Dynamics							*pDynamics;
	Taylor_Model_Computation_Setting	*p_tm_setting;
	Global_Computation_Setting			*p_g_setting;
	Plot_Setting						*p_p_setting;

	std::list<LinearFlowpipe> 			linear_flowpipes;

	Result_of_Reachability				result_of_reachability;

	double								time;
	bool								bSafetyChecking;
	bool								bPlot;
	bool								bTMOutput;
	bool								bPrint;
	bool								bSymbolicRemainder;

	std::vector<Constraint>				unsafeSet;
	std::vector<Flowpipe>				initialSets;

	std::string							fileName;

public:
	Continuous_Reachability();
	Continuous_Reachability(const Continuous_Reachability_Problem_Description & problem_description);
	~Continuous_Reachability();

	void setup(const Continuous_Reachability_Problem_Description & problem_description);

	unsigned long run();

	int safetyChecking();

	void prepareForPlotting();
	void prepareForTMOutput();

	void plot_2D() const;
	void tmOutput(std::ostream & os) const;
};


int safetyChecking(const TaylorModelVec<Real> & tmv, const std::vector<Interval> & domain, const std::vector<Constraint> & unsafeSet, const Taylor_Model_Computation_Setting & tm_setting, const Global_Computation_Setting & g_setting);

void gridBox(std::list<std::vector<Interval> > & grids, const std::vector<Interval> & box, const unsigned int num);

int contract_remainder(const std::vector<Interval> & polyRange, std::vector<Interval> & remainders, const std::vector<Constraint> & constraints);

unsigned int findProperOrder(Real & error, const Real & max, const Real & min, const Real & tolerance, const unsigned int start_order);

void check_connectivities(Matrix<bool> & result, Matrix<bool> & adjMatrix);

void compute_one_step_trans(Matrix<UnivariateTaylorModel<Real> > & utm_Phi_t, Matrix<Real> & rm_Phi_t, Matrix<UnivariateTaylorModel<Real> > & utm_Psi_t, Matrix<Real> & rm_Psi_t,
		Matrix<Interval> & tv_part, const Matrix<UnivariatePolynomial<Real> > & A_t, const Matrix<UnivariatePolynomial<Real> > & B_t, const Matrix<UnivariatePolynomial<Real> > & tv_t,
		Matrix<bool> & connectivity, const bool bAuto, const UnivariatePolynomial<Real> & up_t, const unsigned int order, std::vector<Real> & step_end_exp_table);

}

#endif /* CONTINUOUS_H_ */
