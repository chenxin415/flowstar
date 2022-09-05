/*---
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.
  The code is released as is under the GNU General Public License (GPL).
---*/

#ifndef CONTINUOUS_H_
#define CONTINUOUS_H_

#include "Geometry.h"


namespace flowstar
{

class TaylorModelFlowpipe;
class Flowpipe;
class Result_of_Reachability;

class Symbolic_Remainder
{
public:
	std::vector<Matrix<Interval> > J;
	std::vector<Matrix<Real> > Phi_L;
	std::vector<Real> scalars;
	unsigned int max_size;

public:
	Symbolic_Remainder();
	Symbolic_Remainder(const Flowpipe & initialSet, const unsigned int s);
	Symbolic_Remainder(const Symbolic_Remainder & symbolic_remainder);
	~Symbolic_Remainder();

	void reset(const unsigned int dim);

	Symbolic_Remainder & operator = (const Symbolic_Remainder & symbolic_remainder);
};



class Computational_Setting
{
public:
	Taylor_Model_Setting tm_setting;
	Global_Setting g_setting;
	bool bPrint;

	unsigned int max_order;

public:
	Computational_Setting(const Variables & vars);
	Computational_Setting(const Computational_Setting & setting);
	~Computational_Setting();

	void clear();

	bool setFixedStepsize(const double step, const unsigned int order);
	bool setFixedStepsize(const double step, const unsigned int order_min, const unsigned int order_max);
	bool setAdaptiveStepsize(const double step_min, const double step_max, const unsigned int order);

	bool setCutoffThreshold(const double threshold);
	void setRemainderEstimation(const std::vector<Interval> & estimation);

	// no need to call this function when the expected max order is same as the max flowpipe order
	bool setMaxOrder(const unsigned int order);

	void printOn();
	void printOff();

	Computational_Setting & operator = (const Computational_Setting & setting);
};




// A flowpipe is represented by a composition of two Taylor models. The left Taylor model is the preconditioning part.
class Flowpipe
{
public:
	TaylorModelVec<Real> tmvPre;
	TaylorModelVec<Real> tmv;
	std::vector<Interval> domain;	// domain of TMV_right, the first variable is t
	int safety;
	bool bConstrained;				// whether this flowpipe is contracted by an invariant

public:
	Flowpipe();
	Flowpipe(const std::vector<Interval> & box);
	Flowpipe(Zonotope & zonotope);
	Flowpipe(const TaylorModelFlowpipe & fp);
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

	int safetyChecking(const std::vector<Constraint> & safeSet, const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting) const;
	int unsafetyChecking(const std::vector<Constraint> & unsafeSet, const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting) const;

	bool isInTarget(const std::vector<Constraint> & targetSet, const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting) const;
	bool isInTarget(const std::vector<Constraint> & targetSet, const Computational_Setting & setting) const;

	Flowpipe & operator = (const Flowpipe & flowpipe);


	// interval remainders
	// fixed step sizes and orders
	int advance(Flowpipe & result, const std::vector<Expression<Real> > & ode, const Taylor_Model_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Setting & g_setting) const;

	int advance(Flowpipe & result, const std::vector<Expression<Interval> > & ode, const Taylor_Model_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Setting & g_setting) const;


	// adaptive step sizes and fixed orders
	int advance_adaptive_stepsize(Flowpipe & result, const std::vector<Expression<Real> > & ode, const double new_stepsize, Taylor_Model_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Setting & g_setting) const;

	int advance_adaptive_stepsize(Flowpipe & result, const std::vector<Expression<Interval> > & ode, const double new_stepsize, Taylor_Model_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Setting & g_setting) const;


	// fixed step sizes and adaptive orders
	int advance_adaptive_order(Flowpipe & result, const std::vector<Expression<Real> > & ode, Taylor_Model_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Setting & g_setting) const;

	int advance_adaptive_order(Flowpipe & result, const std::vector<Expression<Interval> > & ode, Taylor_Model_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Setting & g_setting) const;


	// symbolic remainders
	// fixed step sizes and orders
	int advance(Flowpipe & result, const std::vector<Expression<Real> > & ode, const Taylor_Model_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Setting & g_setting, Symbolic_Remainder & symbolic_remainder) const;

	int advance(Flowpipe & result, const std::vector<Expression<Interval> > & ode, const Taylor_Model_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Setting & g_setting, Symbolic_Remainder & symbolic_remainder) const;


	// adaptive step sizes and fixed orders
	int advance_adaptive_stepsize(Flowpipe & result, const std::vector<Expression<Real> > & ode, const double new_stepsize, Taylor_Model_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Setting & g_setting, Symbolic_Remainder & symbolic_remainder) const;

	int advance_adaptive_stepsize(Flowpipe & result, const std::vector<Expression<Interval> > & ode, const double new_stepsize, Taylor_Model_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Setting & g_setting, Symbolic_Remainder & symbolic_remainder) const;


	// fixed step sizes and adaptive orders
	int advance_adaptive_order(Flowpipe & result, const std::vector<Expression<Real> > & ode, Taylor_Model_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Setting & g_setting, Symbolic_Remainder & symbolic_remainder) const;

	int advance_adaptive_order(Flowpipe & result, const std::vector<Expression<Interval> > & ode, Taylor_Model_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Setting & g_setting, Symbolic_Remainder & symbolic_remainder) const;
};






class TaylorModelFlowpipe
{
public:
	TaylorModelVec<Real> tmv_flowpipe;	// Taylor model (vector) for the flowmap function
	std::vector<Interval> domain;		// domain of the Taylor model
	int safety;							// safety of the flowpipe if a safe set is specified
	bool bConstrained;					// whether this flowpipe is contracted by an invariant

public:
	TaylorModelFlowpipe();
	TaylorModelFlowpipe(const TaylorModelVec<Real> & tmv, const std::vector<Interval> & domain);
	TaylorModelFlowpipe(const TaylorModelFlowpipe & flowpipe);
	~TaylorModelFlowpipe();

	TaylorModelFlowpipe & operator = (const TaylorModelFlowpipe & flowpipe);
};





// pool for keeping Taylor model flowpipes
class TaylorModelFlowpipes
{
public:
	std::list<TaylorModelFlowpipe> tmv_flowpipes;

public:
	TaylorModelFlowpipes();
	TaylorModelFlowpipes(const TaylorModelFlowpipes & flowpipes);
	~TaylorModelFlowpipes();

	TaylorModelFlowpipes & operator = (const TaylorModelFlowpipes & flowpipes);

	void clear();
	unsigned int size() const;
	void merge(const TaylorModelFlowpipes & flowpipes);
};






class LinearFlowmap
{
public:
	Matrix<UnivariateTaylorModel<Real> > Phi;	// state-transition matrix
	Matrix<UnivariateTaylorModel<Real> > Psi;	// transition matrix for constants
	Matrix<UnivariateTaylorModel<Real> > Omega;	// transition matrix for constant parameters

	Zonotope tv_remainder;	// zonotope enclosure for the range-bounded uncertainties

	Matrix<Interval> interval_remainder;	// interval range of the range-bounded uncertainties

public:
	LinearFlowmap();
	LinearFlowmap(const LinearFlowmap & flowmap);
	~LinearFlowmap();

	int safetyChecking(const std::vector<Constraint> & safeSet, const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting,
			const TaylorModelVec<Real> & tmv_of_X0, const std::vector<Interval> & polyRangeX0, const std::vector<Interval> & range_of_X0, const std::vector<Interval> & domain);

	void evaluate(TaylorModelVec<Real> & result, const TaylorModelVec<Real> & tmv_of_X0, const std::vector<Interval> & polyRangeX0, const std::vector<Interval> & range_of_X0, const std::vector<Interval> & domain, const Taylor_Model_Setting & tm_setting);

	LinearFlowmap & operator = (const LinearFlowmap & flowmap);

	friend class FlowmapAbstraction;
};


// class for the flowmap abstraction of a LTI or LTV ODE
class FlowmapAbstraction
{
public:
	std::vector<LinearFlowmap> flowmaps;
	double stepsize;

public:
	FlowmapAbstraction();
	FlowmapAbstraction(const FlowmapAbstraction & abstraction);
	~FlowmapAbstraction();

	FlowmapAbstraction & operator = (const FlowmapAbstraction & abstraction);

	void clear();

	void compose(FlowmapAbstraction & result, const FlowmapAbstraction & abstraction, const Computational_Setting & setting, const int zono_order = -1);

	int reach(TaylorModelFlowpipes & result, const Flowpipe & initialSet, const Computational_Setting & setting, const std::vector<Constraint> & safeSet);

	// reachability computation in an invariant
	// the initial set might be contracted by the invariant
	int reach_inv(TaylorModelFlowpipes & result, Flowpipe & initialSet, const std::vector<Constraint> & invariant, const Computational_Setting & setting, const std::vector<Constraint> & safeSet);

	int safetychecking(std::vector<int> & safety, const Flowpipe & initialSet, const std::vector<Constraint> & safeSet, const Computational_Setting & setting);
};








class Result_of_Reachability
{
public:
	int status;

	Flowpipe fp_end_of_time;
	TaylorModelFlowpipe tmv_fp_end_of_time;

	std::list<Flowpipe> flowpipes;

	TaylorModelFlowpipes tmv_flowpipes;


public:
	Result_of_Reachability();
	Result_of_Reachability(const Result_of_Reachability & result);
	~Result_of_Reachability();

	void clear();
	void merge(const Result_of_Reachability & result);

	Flowpipe & safetyChecking(const std::vector<Constraint> & safeSet, const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting);
	Flowpipe & unsafetyChecking(const std::vector<Constraint> & unsafeSet, const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting);

	// transforming all of the flowpipes to Taylor model flowpipes and appending them to tmv_flowpipes
	void transformToTaylorModels(const Computational_Setting & setting);

	bool isSafe() const;
	bool isUnsafe() const;
	bool isCompleted() const;

	Result_of_Reachability & operator = (const Result_of_Reachability & result);

protected:
	void transformToTaylorModels(const Taylor_Model_Setting & tm_setting, const bool bPrint);
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







/*
 * class of linear time-invariant ODEs: x' = A x + B + C p + D v where
 * x are the state variables, p are the time-invariant parameters
 * B is a constant matrix
 * D is the constant matrix for the time-varying uncertainties in [-1,1]
 */


class LTI_ODE : public Dynamics
{
protected:
	Matrix<Real> 							rm_dyn_A;
	Matrix<Real> 							rm_dyn_B;
	Matrix<Real> 							rm_dyn_C;
	Matrix<Real> 							rm_dyn_D;

	Matrix<bool> 							connectivity;

public:
	LTI_ODE(const Matrix<Real> & A);
	LTI_ODE(const Matrix<Real> & A, const Matrix<Real> & B);
	LTI_ODE(const Matrix<Real> & A, const Matrix<Real> & B, const Matrix<Real> & C);
	LTI_ODE(const Matrix<Real> & A, const Matrix<Real> & B, const Matrix<Real> & C, const Matrix<Real> & D);
	LTI_ODE(const LTI_ODE & lti_ode);
	virtual ~LTI_ODE() {} ;

	LTI_ODE & operator = (const LTI_ODE & lti_ode);


//	void evaluate(Interval & result, const unsigned int varID, const TaylorModelVec<Real> & tmv_range, const std::vector<Interval> & domain);
//	void evaluate(Interval & result, const std::vector<Expression<Real> > & coeff_of_Lie_deriv, const TaylorModelVec<Real> & tmv_range, const std::vector<Interval> & domain, const unsigned int order, const Computational_Setting & setting);


	void abstract(FlowmapAbstraction & abstraction, const int N, const int zono_order, const Taylor_Model_Setting & tm_setting,
			const Global_Setting & g_setting);

	void abstract(FlowmapAbstraction & abstraction, const int N, Computational_Setting & setting, const int zono_order = -1);

	int compute_one_flowpipe(TaylorModelFlowpipe & result, const Real & stepsize, const TaylorModelFlowpipe & initialSet, const unsigned int order, const Taylor_Model_Setting & tm_setting);
};




/*
 * class of linear time-invariant ODEs: x' = A(t) x + B(t) + C(t) p + D(t) v where
 * x are the state variables, p are the time-invariant parameters
 * B is a univariate matrix
 * D is the univariate matrix for the time-varying uncertainties in [-1,1]
 */

class LTV_ODE : public Dynamics
{
public:
	Matrix<Expression<Real> >		expr_dyn_A;
	Matrix<Expression<Real> >		expr_dyn_B;
	Matrix<Expression<Real> >		expr_dyn_C;
	Matrix<Expression<Real> >		expr_dyn_D;

	Matrix<bool> 					connectivity;

public:
	LTV_ODE(const Matrix<Expression<Real> > & A);
	LTV_ODE(const Matrix<Expression<Real> > & A, const Matrix<Expression<Real> > & B);
	LTV_ODE(const Matrix<Expression<Real> > & A, const Matrix<Expression<Real> > & B, const Matrix<Expression<Real> > & C);
	LTV_ODE(const Matrix<Expression<Real> > & A, const Matrix<Expression<Real> > & B, const Matrix<Expression<Real> > & C, const Matrix<Expression<Real> > & D);
	LTV_ODE(const LTV_ODE & ltv_ode);
	virtual ~LTV_ODE();

	LTV_ODE & operator = (const LTV_ODE & ltv_ode);

	void abstract(FlowmapAbstraction & abstraction, const double t0, const int N, const int zono_order, const Taylor_Model_Setting & tm_setting,
			const Global_Setting & g_setting);

	void abstract(FlowmapAbstraction & abstraction, const double t0, const int N, Computational_Setting & setting, const int zono_order = -1);

/*
	void evaluate(Interval & result, const unsigned int varID, const TaylorModelVec<Real> & tmv_range, const std::vector<Interval> & domain, const Real & t_lb, const unsigned int order, const Interval & cutoff_threshold);
//	void evaluate(Interval & result, const std::vector<Expression<Real> > & coeff_of_Lie_deriv, const TaylorModelVec<Real> & tmv_range, const std::vector<Interval> & domain, const Real & t_lb, const unsigned int order, const Computational_Setting & setting);

	int reach(std::list<LinearFlowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, const int zono_order, const Taylor_Model_Setting & tm_setting,
			const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput);

	int reach_inv(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, const int zono_order, const std::vector<Constraint> & invariant,
			const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking, const bool bPlot, const bool bTMOutput);

	void reach(Result_of_Reachability & result, const Flowpipe & initialSet, const double T, Computational_Setting & setting, const std::vector<Constraint> & unsafeSet, const int zono_order = -1);
	void reach_inv(Result_of_Reachability & result, const Flowpipe & initialSet, const double T, const std::vector<Constraint> & invariant, Computational_Setting & setting, const std::vector<Constraint> & unsafeSet, const int zono_order = -1);
*/
};







int safetyChecking(const TaylorModelVec<Real> & tmv, const std::vector<Interval> & domain, const std::vector<Constraint> & safeSet, const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting);
int unsafetyChecking(const TaylorModelVec<Real> & tmv, const std::vector<Interval> & domain, const std::vector<Constraint> & unsafeSet, const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting);

int remainder_contraction_int(const std::vector<Interval> & polyRange, std::vector<Interval> & remainders, const std::vector<Constraint> & constraints);
int domain_contraction_int(const TaylorModelVec<Real> & tmv_flowpipe, std::vector<Interval> & domain, const std::vector<Constraint> & constraints, const unsigned int order, const Interval & cutoff_threshold, const Global_Setting & g_setting);




template <class DATA_TYPE>
class ODE : public Dynamics
{
public:
	Variables stateVars;
	std::vector<Expression<DATA_TYPE> >	expressions;

public:
	ODE() {} ;
	ODE(const Variables & vars);
	ODE(const std::vector<std::string> & str_ode, Variables & vars);
	ODE(const ODE<DATA_TYPE> & ode);
	virtual ~ODE() {} ;

	ODE<DATA_TYPE> & operator = (const ODE<DATA_TYPE> & ode);

	bool defDeriv(const std::string & varName, const Expression<DATA_TYPE> & deriv);

	void evaluate(Interval & result, const unsigned int varID, const TaylorModelVec<Real> & tmv_range, const std::vector<Interval> & domain, const unsigned int order, const Computational_Setting & setting);

protected:
	int reach(std::list<Flowpipe> & flowpipes, const double time, const Flowpipe & initialSet, const Taylor_Model_Setting & tm_setting,
			const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & safeSet, const bool bSafetyChecking) const;

	int reach_adaptive_stepsize(std::list<Flowpipe> & flowpipes, const double time, const Flowpipe & initialSet, Taylor_Model_Setting & tm_setting,
			const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & safeSet, const bool bSafetyChecking) const;

	int reach_adaptive_order(std::list<Flowpipe> & flowpipes, const double time, const Flowpipe & initialSet, Taylor_Model_Setting & tm_setting,
			const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & safeSet, const bool bSafetyChecking) const;

	int reach_symbolic_remainder(std::list<Flowpipe> & flowpipes, const double time, const Flowpipe & initialSet, Symbolic_Remainder & symbolic_remainder, const Taylor_Model_Setting & tm_setting,
			const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & safeSet, const bool bSafetyChecking) const;

	int reach_symbolic_remainder_adaptive_stepsize(std::list<Flowpipe> & flowpipes, const double time, const Flowpipe & initialSet, Symbolic_Remainder & symbolic_remainder, Taylor_Model_Setting & tm_setting,
			const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & safeSet, const bool bSafetyChecking) const;

	int reach_symbolic_remainder_adaptive_order(std::list<Flowpipe> & flowpipes, const double time, const Flowpipe & initialSet, Symbolic_Remainder & symbolic_remainder, Taylor_Model_Setting & tm_setting,
			const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & safeSet, const bool bSafetyChecking) const;


	// reachability computation in a given invariant set
	int reach_inv(TaylorModelFlowpipes & tmv_flowpipes, Flowpipe & last_flowpipe, const double time, const Flowpipe & initialSet, const std::vector<Constraint> & invariant,
			const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & safeSet, const bool bSafetyChecking) const;

	int reach_inv_adaptive_stepsize(TaylorModelFlowpipes & tmv_flowpipes, Flowpipe & last_flowpipe, const double time, const Flowpipe & initialSet, const std::vector<Constraint> & invariant,
			Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & safeSet, const bool bSafetyChecking) const;

	int reach_inv_adaptive_order(TaylorModelFlowpipes & tmv_flowpipes, Flowpipe & last_flowpipe, const double time, const Flowpipe & initialSet, const std::vector<Constraint> & invariant,
			Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & safeSet, const bool bSafetyChecking) const;


	int reach_inv_symbolic_remainder(TaylorModelFlowpipes & tmv_flowpipes, Flowpipe & last_flowpipe, const double time, const Flowpipe & initialSet, const std::vector<Constraint> & invariant, Symbolic_Remainder & symbolic_remainder,
			const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & safeSet, const bool bSafetyChecking) const;

	int reach_inv_symbolic_remainder_adaptive_stepsize(TaylorModelFlowpipes & tmv_flowpipes, Flowpipe & last_flowpipe, const double time, const Flowpipe & initialSet, const std::vector<Constraint> & invariant, Symbolic_Remainder & symbolic_remainder,
			Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & safeSet, const bool bSafetyChecking) const;

	int reach_inv_symbolic_remainder_adaptive_order(TaylorModelFlowpipes & tmv_flowpipes, Flowpipe & last_flowpipe, const double time, const Flowpipe & initialSet, const std::vector<Constraint> & invariant, Symbolic_Remainder & symbolic_remainder,
			Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & safeSet, const bool bSafetyChecking) const;


public:
	void reach(Result_of_Reachability & result, const Flowpipe & initialSet, const double T, Computational_Setting & setting, const std::vector<Constraint> & safeSet) const;
	void reach(Result_of_Reachability & result, const Flowpipe & initialSet, const double T, Computational_Setting & setting, const std::vector<Constraint> & safeSet, Symbolic_Remainder & symbolic_remainder) const;

	void reach_inv(Result_of_Reachability & result, const Flowpipe & initialSet, const double T, const std::vector<Constraint> & invariant, Computational_Setting & setting, const std::vector<Constraint> & safeSet) const;
	void reach_inv(Result_of_Reachability & result, const Flowpipe & initialSet, const double T, const std::vector<Constraint> & invariant, Computational_Setting & setting, const std::vector<Constraint> & safeSet, Symbolic_Remainder & symbolic_remainder) const;

};



template <class DATA_TYPE>
ODE<DATA_TYPE>::ODE(const Variables & vars)
{
	stateVars = vars;

	Expression<DATA_TYPE> zero(DATA_TYPE(0));
	expressions.resize(vars.size(), zero);
}

template <class DATA_TYPE>
ODE<DATA_TYPE>::ODE(const std::vector<std::string> & str_ode, Variables & vars)
{
	if(str_ode.size() > vars.size())
	{
		printf("ODE: There are more derivatives than the state variables.\n");

		Expression<DATA_TYPE> zero(DATA_TYPE(0));
		expressions.resize(vars.size(), zero);
	}
	else
	{
		Expression<DATA_TYPE> zero(DATA_TYPE(0));
		expressions.resize(vars.size(), zero);

		for(int i=0; i<str_ode.size(); ++i)
		{
			Expression<DATA_TYPE> deriv(str_ode[i], vars);
			expressions[i] = deriv;
		}
	}
}

template <class DATA_TYPE>
ODE<DATA_TYPE>::ODE(const ODE<DATA_TYPE> & ode)
{
	stateVars	= ode.stateVars;
	expressions	= ode.expressions;
}

template <class DATA_TYPE>
ODE<DATA_TYPE> & ODE<DATA_TYPE>::operator = (const ODE<DATA_TYPE> & ode)
{
	if(this == &ode)
		return *this;

	stateVars	= ode.stateVars;
	expressions	= ode.expressions;

	return *this;
}

template <class DATA_TYPE>
bool ODE<DATA_TYPE>::defDeriv(const std::string & varName, const Expression<DATA_TYPE> & deriv)
{
	int varID = stateVars.getIDForVar(varName);

	if(varID < 0)
	{
		printf("State variable %s is not declared.\n", varName.c_str());
		return false;
	}

	expressions[varID] = deriv;

	return true;
}

template <class DATA_TYPE>
void ODE<DATA_TYPE>::evaluate(Interval & result, const unsigned int varID, const TaylorModelVec<Real> & tmv_range, const std::vector<Interval> & domain, const unsigned int order, const Computational_Setting & setting)
{
	TaylorModel<Real> tm;
	expressions[varID].evaluate(tm, tmv_range.tms, order, domain, setting.tm_setting.cutoff_threshold, setting.g_setting);

	tm.intEval(result, domain);
}

template <class DATA_TYPE>
int ODE<DATA_TYPE>::reach(std::list<Flowpipe> & flowpipes, const double time, const Flowpipe & initialSet, const Taylor_Model_Setting & tm_setting,
		const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & safeSet, const bool bSafetyChecking) const
{
	std::vector<Constraint> dummy_invariant;

	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = currentFlowpipe.advance(newFlowpipe, expressions, tm_setting, dummy_invariant, g_setting);

		if(res == 1)
		{
			double remaining_time = time - t;

			if(remaining_time < step)
			{
				step = remaining_time;
				newFlowpipe.domain[0].setSup(remaining_time);
			}

			if(bSafetyChecking)
			{
				int safety = newFlowpipe.safetyChecking(safeSet, tm_setting, g_setting);

				newFlowpipe.safety = safety;
				flowpipes.push_back(newFlowpipe);

				if(safety == UNSAFE)
				{
					return COMPLETED_UNSAFE;
				}
				else if(safety == UNKNOWN && checking_result == COMPLETED_SAFE)
				{
					checking_result = COMPLETED_UNKNOWN;
				}
			}
			else
			{
				newFlowpipe.safety = SAFE;
				flowpipes.push_back(newFlowpipe);
			}

			currentFlowpipe = newFlowpipe;

			t += step;

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("order = %d\n", tm_setting.order);
			}
		}
		else
		{
			switch(checking_result)
			{
			case COMPLETED_SAFE:
				return UNCOMPLETED_SAFE;
			case COMPLETED_UNSAFE:
				return UNCOMPLETED_UNSAFE;
			case COMPLETED_UNKNOWN:
				return UNCOMPLETED_UNKNOWN;
			}
		}
	}

	return checking_result;
}

template <class DATA_TYPE>
int ODE<DATA_TYPE>::reach_adaptive_stepsize(std::list<Flowpipe> & flowpipes, const double time, const Flowpipe & initialSet, Taylor_Model_Setting & tm_setting,
		const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & safeSet, const bool bSafetyChecking) const
{
	std::vector<Constraint> dummy_invariant;

	int checking_result = COMPLETED_SAFE;

	Flowpipe newFlowpipe, currentFlowpipe = initialSet;
	double new_stepsize = -1;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = currentFlowpipe.advance_adaptive_stepsize(newFlowpipe, expressions, new_stepsize, tm_setting, dummy_invariant, g_setting);

		if(res == 1)
		{
			double current_stepsize = tm_setting.step_exp_table[1].sup();

			double remaining_time = time - t;

			if(remaining_time < current_stepsize)
			{
				current_stepsize = remaining_time;
				newFlowpipe.domain[0].setSup(remaining_time);
			}

			if(bSafetyChecking)
			{
				int safety = newFlowpipe.safetyChecking(safeSet, tm_setting, g_setting);

				newFlowpipe.safety = safety;
				flowpipes.push_back(newFlowpipe);

				if(safety == UNSAFE)
				{
					return COMPLETED_UNSAFE;
				}
				else if(safety == UNKNOWN && checking_result == COMPLETED_SAFE)
				{
					checking_result = COMPLETED_UNKNOWN;
				}
			}
			else
			{
				newFlowpipe.safety = SAFE;
				flowpipes.push_back(newFlowpipe);
			}

			currentFlowpipe = newFlowpipe;


			t += current_stepsize;

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", current_stepsize);
				printf("order = %d\n", tm_setting.order);
			}

			new_stepsize = current_stepsize * LAMBDA_UP;
			if(new_stepsize > tm_setting.step_max - THRESHOLD_HIGH)
			{
				new_stepsize = -1;
			}
		}
		else
		{
			switch(checking_result)
			{
			case COMPLETED_SAFE:
				return UNCOMPLETED_SAFE;
			case COMPLETED_UNSAFE:
				return UNCOMPLETED_UNSAFE;
			case COMPLETED_UNKNOWN:
				return UNCOMPLETED_UNKNOWN;
			}
		}
	}

	return checking_result;
}

template <class DATA_TYPE>
int ODE<DATA_TYPE>::reach_adaptive_order(std::list<Flowpipe> & flowpipes, const double time, const Flowpipe & initialSet, Taylor_Model_Setting & tm_setting,
		const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & safeSet, const bool bSafetyChecking) const
{
	std::vector<Constraint> dummy_invariant;

	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = currentFlowpipe.advance_adaptive_order(newFlowpipe, expressions, tm_setting, dummy_invariant, g_setting);

		if(res == 1)
		{
			double remaining_time = time - t;

			if(remaining_time < step)
			{
				step = remaining_time;
				newFlowpipe.domain[0].setSup(remaining_time);
			}

			if(bSafetyChecking)
			{
				int safety = newFlowpipe.safetyChecking(safeSet, tm_setting, g_setting);

				newFlowpipe.safety = safety;
				flowpipes.push_back(newFlowpipe);

				if(safety == UNSAFE)
				{
					return COMPLETED_UNSAFE;
				}
				else if(safety == UNKNOWN && checking_result == COMPLETED_SAFE)
				{
					checking_result = COMPLETED_UNKNOWN;
				}
			}
			else
			{
				newFlowpipe.safety = SAFE;
				flowpipes.push_back(newFlowpipe);
			}

			currentFlowpipe = newFlowpipe;

			t += step;

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("order = %d\n", tm_setting.order);
			}

			if(tm_setting.order > tm_setting.order_min)
			{
				--tm_setting.order;
			}
		}
		else
		{
			switch(checking_result)
			{
			case COMPLETED_SAFE:
				return UNCOMPLETED_SAFE;
			case COMPLETED_UNSAFE:
				return UNCOMPLETED_UNSAFE;
			case COMPLETED_UNKNOWN:
				return UNCOMPLETED_UNKNOWN;
			}
		}
	}


	return checking_result;
}

template <class DATA_TYPE>
int ODE<DATA_TYPE>::reach_symbolic_remainder(std::list<Flowpipe> & flowpipes, const double time, const Flowpipe & initialSet, Symbolic_Remainder & symbolic_remainder, const Taylor_Model_Setting & tm_setting,
		const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & safeSet, const bool bSafetyChecking) const
{
	std::vector<Constraint> dummy_invariant;

	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = currentFlowpipe.advance(newFlowpipe, expressions, tm_setting, dummy_invariant, g_setting, symbolic_remainder);

		if(res == 1)
		{
			double remaining_time = time - t;

			if(remaining_time < step)
			{
				step = remaining_time;
				newFlowpipe.domain[0].setSup(remaining_time);
			}

			if(bSafetyChecking)
			{
				int safety = newFlowpipe.safetyChecking(safeSet, tm_setting, g_setting);

				newFlowpipe.safety = safety;
				flowpipes.push_back(newFlowpipe);

				if(safety == UNSAFE)
				{
					return COMPLETED_UNSAFE;
				}
				else if(safety == UNKNOWN && checking_result == COMPLETED_SAFE)
				{
					checking_result = COMPLETED_UNKNOWN;
				}
			}
			else
			{
				newFlowpipe.safety = SAFE;
				flowpipes.push_back(newFlowpipe);
			}

			currentFlowpipe = newFlowpipe;


			t += step;

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("order = %d\n", tm_setting.order);
			}

			if(symbolic_remainder.J.size() >= symbolic_remainder.max_size)
			{
				symbolic_remainder.reset(currentFlowpipe.tmvPre.tms.size());
			}
		}
		else
		{
			switch(checking_result)
			{
			case COMPLETED_SAFE:
				return UNCOMPLETED_SAFE;
			case COMPLETED_UNSAFE:
				return UNCOMPLETED_UNSAFE;
			case COMPLETED_UNKNOWN:
				return UNCOMPLETED_UNKNOWN;
			}
		}
	}

	return checking_result;
}

template <class DATA_TYPE>
int ODE<DATA_TYPE>::reach_symbolic_remainder_adaptive_stepsize(std::list<Flowpipe> & flowpipes, const double time, const Flowpipe & initialSet, Symbolic_Remainder & symbolic_remainder, Taylor_Model_Setting & tm_setting,
		const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & safeSet, const bool bSafetyChecking) const
{
	std::vector<Constraint> dummy_invariant;

	int checking_result = COMPLETED_SAFE;

	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	tm_setting.setStepsize(tm_setting.step_max, tm_setting.order);
	double new_stepsize = -1;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = currentFlowpipe.advance_adaptive_stepsize(newFlowpipe, expressions, new_stepsize, tm_setting, dummy_invariant, g_setting, symbolic_remainder);

		if(res == 1)
		{
			double current_stepsize = tm_setting.step_exp_table[1].sup();

			double remaining_time = time - t;

			if(remaining_time < current_stepsize)
			{
				current_stepsize = remaining_time;
				newFlowpipe.domain[0].setSup(remaining_time);
			}

			if(bSafetyChecking)
			{
				int safety = newFlowpipe.safetyChecking(safeSet, tm_setting, g_setting);

				newFlowpipe.safety = safety;
				flowpipes.push_back(newFlowpipe);

				if(safety == UNSAFE)
				{
					return COMPLETED_UNSAFE;
				}
				else if(safety == UNKNOWN && checking_result == COMPLETED_SAFE)
				{
					checking_result = COMPLETED_UNKNOWN;
				}
			}
			else
			{
				newFlowpipe.safety = SAFE;
				flowpipes.push_back(newFlowpipe);
			}

			currentFlowpipe = newFlowpipe;

			if(symbolic_remainder.J.size() >= symbolic_remainder.max_size)
			{
				symbolic_remainder.reset(currentFlowpipe.tmvPre.tms.size());
			}

			t += current_stepsize;

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", current_stepsize);
				printf("order = %d\n", tm_setting.order);
			}

			new_stepsize = current_stepsize * LAMBDA_UP;
			if(new_stepsize > tm_setting.step_max - THRESHOLD_HIGH)
			{
				new_stepsize = -1;
			}
		}
		else
		{
			switch(checking_result)
			{
			case COMPLETED_SAFE:
				return UNCOMPLETED_SAFE;
			case COMPLETED_UNSAFE:
				return UNCOMPLETED_UNSAFE;
			case COMPLETED_UNKNOWN:
				return UNCOMPLETED_UNKNOWN;
			}
		}
	}

	return checking_result;
}

template <class DATA_TYPE>
int ODE<DATA_TYPE>::reach_symbolic_remainder_adaptive_order(std::list<Flowpipe> & flowpipes, const double time, const Flowpipe & initialSet, Symbolic_Remainder & symbolic_remainder, Taylor_Model_Setting & tm_setting,
		const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & safeSet, const bool bSafetyChecking) const
{
	std::vector<Constraint> dummy_invariant;

	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = currentFlowpipe.advance_adaptive_order(newFlowpipe, expressions, tm_setting, dummy_invariant, g_setting, symbolic_remainder);

		if(res == 1)
		{
			double remaining_time = time - t;

			if(remaining_time < step)
			{
				step = remaining_time;
				newFlowpipe.domain[0].setSup(remaining_time);
			}

			if(bSafetyChecking)
			{
				int safety = newFlowpipe.safetyChecking(safeSet, tm_setting, g_setting);

				newFlowpipe.safety = safety;
				flowpipes.push_back(newFlowpipe);

				if(safety == UNSAFE)
				{
					return COMPLETED_UNSAFE;
				}
				else if(safety == UNKNOWN && checking_result == COMPLETED_SAFE)
				{
					checking_result = COMPLETED_UNKNOWN;
				}
			}
			else
			{
				newFlowpipe.safety = SAFE;
				flowpipes.push_back(newFlowpipe);
			}

			currentFlowpipe = newFlowpipe;

			t += step;

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("order = %d\n", tm_setting.order);
			}

			if(symbolic_remainder.J.size() >= symbolic_remainder.max_size)
			{
				symbolic_remainder.reset(currentFlowpipe.tmvPre.tms.size());
			}

			if(tm_setting.order > tm_setting.order_min)
			{
				--tm_setting.order;
			}
		}
		else
		{
			switch(checking_result)
			{
			case COMPLETED_SAFE:
				return UNCOMPLETED_SAFE;
			case COMPLETED_UNSAFE:
				return UNCOMPLETED_UNSAFE;
			case COMPLETED_UNKNOWN:
				return UNCOMPLETED_UNKNOWN;
			}
		}
	}

	return checking_result;
}

template <class DATA_TYPE>
void ODE<DATA_TYPE>::reach(Result_of_Reachability & result, const Flowpipe & initialSet, const double T, Computational_Setting & setting, const std::vector<Constraint> & safeSet) const
{
	if(T <= 0)
	{
		printf("Time horizon should be positive.\n");
		return;
	}

	bool bSafetyChecking = false;

	if(safeSet.size() > 0)
	{
		bSafetyChecking = true;
	}

	if(setting.tm_setting.step_min > 0)
	{
		// adaptive stepsizes
		result.status = reach_adaptive_stepsize(result.flowpipes, T, initialSet, setting.tm_setting, setting.g_setting, setting.bPrint, safeSet, bSafetyChecking);
	}
	else if(setting.tm_setting.order_max > 0)
	{
		// adaptive orders
		result.status = reach_adaptive_order(result.flowpipes, T, initialSet, setting.tm_setting, setting.g_setting, setting.bPrint, safeSet, bSafetyChecking);
	}
	else
	{
		// fixed stepsizes and orders
		result.status = reach(result.flowpipes, T, initialSet, setting.tm_setting, setting.g_setting, setting.bPrint, safeSet, bSafetyChecking);
	}


	if(result.flowpipes.size() > 0)
	{
		Flowpipe fpTmp = result.flowpipes.back();
		result.fp_end_of_time = fpTmp;

		Real t;
		fpTmp.domain[0].sup(t);

		std::vector<Real> realVec;
		realVec.push_back(1);
		realVec.push_back(t);

		Real tmp = t;

		for(unsigned int i=2; i<=setting.tm_setting.step_end_exp_table.size(); ++i)
		{
			tmp *= t;
			realVec.push_back(tmp);
		}


		fpTmp.tmvPre.evaluate_time(result.fp_end_of_time.tmvPre, realVec);
	}
}

template <class DATA_TYPE>
void ODE<DATA_TYPE>::reach(Result_of_Reachability & result, const Flowpipe & initialSet, const double T, Computational_Setting & setting, const std::vector<Constraint> & safeSet, Symbolic_Remainder & symbolic_remainder) const
{
	if(T <= 0)
	{
		printf("Time horizon should be positive.\n");
		return;
	}

	bool bSafetyChecking = false;

	if(safeSet.size() > 0)
	{
		bSafetyChecking = true;
	}

	if(setting.tm_setting.step_min > 0)
	{
		// adaptive stepsizes
		result.status = reach_symbolic_remainder_adaptive_stepsize(result.flowpipes, T, initialSet, symbolic_remainder, setting.tm_setting, setting.g_setting, setting.bPrint, safeSet, bSafetyChecking);
	}
	else if(setting.tm_setting.order_max > 0)
	{
		// adaptive orders
		result.status = reach_symbolic_remainder_adaptive_order(result.flowpipes, T, initialSet, symbolic_remainder, setting.tm_setting, setting.g_setting, setting.bPrint, safeSet, bSafetyChecking);
	}
	else
	{
		// fixed stepsizes and orders
		result.status = reach_symbolic_remainder(result.flowpipes, T, initialSet, symbolic_remainder, setting.tm_setting, setting.g_setting, setting.bPrint, safeSet, bSafetyChecking);
	}

	if(result.flowpipes.size() > 0)
	{
		Flowpipe fpTmp = result.flowpipes.back();
		result.fp_end_of_time = fpTmp;

		Real t;
		fpTmp.domain[0].sup(t);

		std::vector<Real> realVec;
		realVec.push_back(1);
		realVec.push_back(t);

		Real tmp = t;

		for(unsigned int i=2; i<=setting.tm_setting.step_end_exp_table.size(); ++i)
		{
			tmp *= t;
			realVec.push_back(tmp);
		}


		fpTmp.tmvPre.evaluate_time(result.fp_end_of_time.tmvPre, realVec);
	}
}

template <class DATA_TYPE>
int ODE<DATA_TYPE>::reach_inv(TaylorModelFlowpipes & flowpipes, Flowpipe & last_flowpipe, const double time, const Flowpipe & initialSet, const std::vector<Constraint> & invariant,
		const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & safeSet, const bool bSafetyChecking) const
{
	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	Flowpipe new_flowpipe, current_flowpipe = initialSet;
	bool bContracted = initialSet.bConstrained;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = current_flowpipe.advance(new_flowpipe, expressions, tm_setting, invariant, g_setting);

		bool step_changed = false;
		double remaining_time = time - t;

		if(remaining_time < step)
		{
			step = remaining_time;
			new_flowpipe.domain[0].setSup(remaining_time);
			step_changed = true;
		}

		if(res == 1)
		{
			// constrain the flowpipe by the invariant
			TaylorModelVec<Real> tmv_flowpipe;

			if(step_changed)
			{
				new_flowpipe.compose(tmv_flowpipe, tm_setting.order, tm_setting.cutoff_threshold);
			}
			else
			{
				new_flowpipe.compose_normal(tmv_flowpipe, tm_setting.step_exp_table, tm_setting.order, tm_setting.cutoff_threshold);
			}

			std::vector<Interval> contracted_domain = new_flowpipe.domain;
			int type = domain_contraction_int(tmv_flowpipe, contracted_domain, invariant, tm_setting.order, tm_setting.cutoff_threshold, g_setting);

			switch(type)
			{
			case UNSAT:		// the intersection is empty
				last_flowpipe = current_flowpipe;
				return checking_result;
			case SAT:		// the domain is not contracted
			{
				TaylorModelFlowpipe flowpipe;
				flowpipe.tmv_flowpipe = tmv_flowpipe;
				flowpipe.domain = contracted_domain;
				flowpipe.bConstrained = bContracted;

				if(bSafetyChecking)
				{
					int safety = safetyChecking(tmv_flowpipe, contracted_domain, safeSet, tm_setting, g_setting);

					if(safety == SAFE)
					{
						flowpipe.safety = SAFE;
						flowpipes.tmv_flowpipes.push_back(flowpipe);
					}
					else if(safety == UNSAFE)
					{
						if(!bContracted)
						{
							flowpipe.safety = UNSAFE;
							flowpipes.tmv_flowpipes.push_back(flowpipe);
							return COMPLETED_UNSAFE;
						}
						else
						{
							flowpipe.safety = UNKNOWN;
							flowpipes.tmv_flowpipes.push_back(flowpipe);
							return COMPLETED_UNKNOWN;
						}
					}
					else
					{
						flowpipe.safety = UNKNOWN;
						flowpipes.tmv_flowpipes.push_back(flowpipe);

						if(checking_result == COMPLETED_SAFE)
						{
							checking_result = COMPLETED_UNKNOWN;
						}
					}
				}
				else
				{
					flowpipe.safety = SAFE;
					flowpipes.tmv_flowpipes.push_back(flowpipe);
				}

				current_flowpipe = new_flowpipe;
				break;
			}
			case CONTRACTED: 	// the domain is contracted but the time interval is not
			{
				bContracted = true;

				TaylorModelFlowpipe flowpipe;
				flowpipe.tmv_flowpipe = tmv_flowpipe;
				flowpipe.domain = contracted_domain;
				flowpipe.bConstrained = true;

				new_flowpipe.domain = contracted_domain;
				new_flowpipe.normalize(tm_setting.cutoff_threshold);

				if(bSafetyChecking)
				{
					int safety = safetyChecking(tmv_flowpipe, contracted_domain, safeSet, tm_setting, g_setting);

					if(safety == SAFE)
					{
						flowpipe.safety = SAFE;
						flowpipes.tmv_flowpipes.push_back(flowpipe);
					}
					else if(safety == UNSAFE)
					{
						flowpipe.safety = UNKNOWN;
						flowpipes.tmv_flowpipes.push_back(flowpipe);
						return COMPLETED_UNKNOWN;
					}
					else
					{
						flowpipe.safety = UNKNOWN;
						flowpipes.tmv_flowpipes.push_back(flowpipe);

						if(checking_result == COMPLETED_SAFE)
						{
							checking_result = COMPLETED_UNKNOWN;
						}
					}
				}
				else
				{
					flowpipe.safety = SAFE;
					flowpipes.tmv_flowpipes.push_back(flowpipe);
				}

				current_flowpipe = new_flowpipe;

				break;
			}
			case TIME_RANGE_CONTRACTED: 	// time interval is contracted
			{
				Real zero(0);

				if(contracted_domain[0].greaterThan(zero))
				{
					return checking_result;
				}
				else
				{
					bContracted = true;

					TaylorModelFlowpipe flowpipe;
					flowpipe.tmv_flowpipe = tmv_flowpipe;
					flowpipe.domain = contracted_domain;
					flowpipe.bConstrained = true;

					new_flowpipe.domain = contracted_domain;
					new_flowpipe.normalize(tm_setting.cutoff_threshold);

					if(bSafetyChecking)
					{
						int safety = safetyChecking(tmv_flowpipe, contracted_domain, safeSet, tm_setting, g_setting);

						if(safety == SAFE)
						{
							flowpipe.safety = SAFE;
							flowpipes.tmv_flowpipes.push_back(flowpipe);
						}
						else
						{
							flowpipe.safety = UNKNOWN;
							flowpipes.tmv_flowpipes.push_back(flowpipe);

							if(checking_result == COMPLETED_SAFE)
							{
								checking_result = COMPLETED_UNKNOWN;
							}
						}
					}
					else
					{
						flowpipe.safety = SAFE;
						flowpipes.tmv_flowpipes.push_back(flowpipe);
					}

					t += contracted_domain[0].sup();

					if(bPrint)
					{
						printf("time = %f,\t", t);
						printf("step = %f,\t", contracted_domain[0].sup());
						printf("order = %d\n", tm_setting.order);
					}

					last_flowpipe = new_flowpipe;

					return checking_result;
				}
			}
			}


			t += contracted_domain[0].sup();

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("order = %d\n", tm_setting.order);
			}
		}
		else
		{
			switch(checking_result)
			{
			case COMPLETED_SAFE:
				return UNCOMPLETED_SAFE;
			case COMPLETED_UNSAFE:
				return UNCOMPLETED_UNSAFE;
			case COMPLETED_UNKNOWN:
				return UNCOMPLETED_UNKNOWN;
			}
		}
	}

	last_flowpipe = current_flowpipe;

	return checking_result;
}

template <class DATA_TYPE>
int ODE<DATA_TYPE>::reach_inv_adaptive_stepsize(TaylorModelFlowpipes & flowpipes, Flowpipe & last_flowpipe, const double time, const Flowpipe & initialSet, const std::vector<Constraint> & invariant,
		Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & safeSet, const bool bSafetyChecking) const
{
	int checking_result = COMPLETED_SAFE;

	Flowpipe new_flowpipe, current_flowpipe = initialSet;
	bool bContracted = initialSet.bConstrained;
	double new_stepsize = -1;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = current_flowpipe.advance_adaptive_stepsize(new_flowpipe, expressions, new_stepsize, tm_setting, invariant, g_setting);

		double current_stepsize = tm_setting.step_exp_table[1].sup();

		bool step_changed = false;
		double remaining_time = time - t;

		if(remaining_time < current_stepsize)
		{
			current_stepsize = remaining_time;
			new_flowpipe.domain[0].setSup(remaining_time);
			step_changed = true;
		}

		if(res == 1)
		{
			// constrain the flowpipe by the invariant
			TaylorModelVec<Real> tmv_flowpipe;

			if(step_changed)
			{
				new_flowpipe.compose(tmv_flowpipe, tm_setting.order, tm_setting.cutoff_threshold);
			}
			else
			{
				new_flowpipe.compose_normal(tmv_flowpipe, tm_setting.step_exp_table, tm_setting.order, tm_setting.cutoff_threshold);
			}

			std::vector<Interval> contracted_domain = new_flowpipe.domain;
			int type = domain_contraction_int(tmv_flowpipe, contracted_domain, invariant, tm_setting.order, tm_setting.cutoff_threshold, g_setting);

			switch(type)
			{
			case UNSAT:		// the intersection is empty
				last_flowpipe = current_flowpipe;
				return checking_result;
			case SAT:		// the domain is not contracted
			{
				TaylorModelFlowpipe flowpipe;
				flowpipe.tmv_flowpipe = tmv_flowpipe;
				flowpipe.domain = contracted_domain;
				flowpipe.bConstrained = bContracted;

				if(bSafetyChecking)
				{
					int safety = safetyChecking(tmv_flowpipe, contracted_domain, safeSet, tm_setting, g_setting);

					if(safety == SAFE)
					{
						flowpipe.safety = SAFE;
						flowpipes.tmv_flowpipes.push_back(flowpipe);
					}
					else if(safety == UNSAFE)
					{
						if(!bContracted)
						{
							flowpipe.safety = UNSAFE;
							flowpipes.tmv_flowpipes.push_back(flowpipe);
							return COMPLETED_UNSAFE;
						}
						else
						{
							flowpipe.safety = UNKNOWN;
							flowpipes.tmv_flowpipes.push_back(flowpipe);
							return COMPLETED_UNKNOWN;
						}
					}
					else
					{
						flowpipe.safety = UNKNOWN;
						flowpipes.tmv_flowpipes.push_back(flowpipe);

						if(checking_result == COMPLETED_SAFE)
						{
							checking_result = COMPLETED_UNKNOWN;
						}
					}
				}
				else
				{
					flowpipe.safety = SAFE;
					flowpipes.tmv_flowpipes.push_back(flowpipe);
				}

				current_flowpipe = new_flowpipe;
				break;
			}
			case CONTRACTED: 	// the domain is contracted but the time interval is not
			{
				bContracted = true;

				TaylorModelFlowpipe flowpipe;
				flowpipe.tmv_flowpipe = tmv_flowpipe;
				flowpipe.domain = contracted_domain;
				flowpipe.bConstrained = true;

				new_flowpipe.domain = contracted_domain;
				new_flowpipe.normalize(tm_setting.cutoff_threshold);

				if(bSafetyChecking)
				{
					int safety = safetyChecking(tmv_flowpipe, contracted_domain, safeSet, tm_setting, g_setting);

					if(safety == SAFE)
					{
						flowpipe.safety = SAFE;
						flowpipes.tmv_flowpipes.push_back(flowpipe);
					}
					else if(safety == UNSAFE)
					{
						flowpipe.safety = UNKNOWN;
						flowpipes.tmv_flowpipes.push_back(flowpipe);
						return COMPLETED_UNKNOWN;
					}
					else
					{
						flowpipe.safety = UNKNOWN;
						flowpipes.tmv_flowpipes.push_back(flowpipe);

						if(checking_result == COMPLETED_SAFE)
						{
							checking_result = COMPLETED_UNKNOWN;
						}
					}
				}
				else
				{
					flowpipe.safety = SAFE;
					flowpipes.tmv_flowpipes.push_back(flowpipe);
				}

				current_flowpipe = new_flowpipe;

				break;
			}
			case TIME_RANGE_CONTRACTED: 	// time interval is contracted
			{
				Real zero(0);

				if(contracted_domain[0].greaterThan(zero))
				{
					return checking_result;
				}
				else
				{
					bContracted = true;

					TaylorModelFlowpipe flowpipe;
					flowpipe.tmv_flowpipe = tmv_flowpipe;
					flowpipe.domain = contracted_domain;
					flowpipe.bConstrained = true;

					new_flowpipe.domain = contracted_domain;
					new_flowpipe.normalize(tm_setting.cutoff_threshold);

					if(bSafetyChecking)
					{
						int safety = safetyChecking(tmv_flowpipe, contracted_domain, safeSet, tm_setting, g_setting);

						if(safety == SAFE)
						{
							flowpipe.safety = SAFE;
							flowpipes.tmv_flowpipes.push_back(flowpipe);
						}
						else
						{
							flowpipe.safety = UNKNOWN;
							flowpipes.tmv_flowpipes.push_back(flowpipe);

							if(checking_result == COMPLETED_SAFE)
							{
								checking_result = COMPLETED_UNKNOWN;
							}
						}
					}
					else
					{
						flowpipe.safety = SAFE;
						flowpipes.tmv_flowpipes.push_back(flowpipe);
					}

					t += contracted_domain[0].sup();

					if(bPrint)
					{
						printf("time = %f,\t", t);
						printf("step = %f,\t", contracted_domain[0].sup());
						printf("order = %d\n", tm_setting.order);
					}

					last_flowpipe = new_flowpipe;

					return checking_result;
				}
			}
			}

			t += current_stepsize;

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", current_stepsize);
				printf("order = %d\n", tm_setting.order);
			}

			new_stepsize = current_stepsize * LAMBDA_UP;
			if(new_stepsize > tm_setting.step_max - THRESHOLD_HIGH)
			{
				new_stepsize = -1;
			}
		}
		else
		{
			switch(checking_result)
			{
			case COMPLETED_SAFE:
				return UNCOMPLETED_SAFE;
			case COMPLETED_UNSAFE:
				return UNCOMPLETED_UNSAFE;
			case COMPLETED_UNKNOWN:
				return UNCOMPLETED_UNKNOWN;
			}
		}
	}

	last_flowpipe = current_flowpipe;

	return checking_result;
}

template <class DATA_TYPE>
int ODE<DATA_TYPE>::reach_inv_adaptive_order(TaylorModelFlowpipes & flowpipes, Flowpipe & last_flowpipe, const double time, const Flowpipe & initialSet, const std::vector<Constraint> & invariant,
		Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & safeSet, const bool bSafetyChecking) const
{
	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	Flowpipe new_flowpipe, current_flowpipe = initialSet;
	bool bContracted = initialSet.bConstrained;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = current_flowpipe.advance_adaptive_order(new_flowpipe, expressions, tm_setting, invariant, g_setting);

		bool step_changed = false;
		double remaining_time = time - t;

		if(remaining_time < step)
		{
			step = remaining_time;
			new_flowpipe.domain[0].setSup(remaining_time);
			step_changed = true;
		}

		if(res == 1)
		{
			// constrain the flowpipe by the invariant
			TaylorModelVec<Real> tmv_flowpipe;

			if(step_changed)
			{
				new_flowpipe.compose(tmv_flowpipe, tm_setting.order, tm_setting.cutoff_threshold);
			}
			else
			{
				new_flowpipe.compose_normal(tmv_flowpipe, tm_setting.step_exp_table, tm_setting.order, tm_setting.cutoff_threshold);
			}

			std::vector<Interval> contracted_domain = new_flowpipe.domain;
			int type = domain_contraction_int(tmv_flowpipe, contracted_domain, invariant, tm_setting.order, tm_setting.cutoff_threshold, g_setting);

			switch(type)
			{
			case UNSAT:		// the intersection is empty
				last_flowpipe = current_flowpipe;
				return checking_result;
			case SAT:		// the domain is not contracted
			{
				TaylorModelFlowpipe flowpipe;
				flowpipe.tmv_flowpipe = tmv_flowpipe;
				flowpipe.domain = contracted_domain;
				flowpipe.bConstrained = bContracted;

				if(bSafetyChecking)
				{
					int safety = safetyChecking(tmv_flowpipe, contracted_domain, safeSet, tm_setting, g_setting);

					if(safety == SAFE)
					{
						flowpipe.safety = SAFE;
						flowpipes.tmv_flowpipes.push_back(flowpipe);
					}
					else if(safety == UNSAFE)
					{
						if(!bContracted)
						{
							flowpipe.safety = UNSAFE;
							flowpipes.tmv_flowpipes.push_back(flowpipe);
							return COMPLETED_UNSAFE;
						}
						else
						{
							flowpipe.safety = UNKNOWN;
							flowpipes.tmv_flowpipes.push_back(flowpipe);
							return COMPLETED_UNKNOWN;
						}
					}
					else
					{
						flowpipe.safety = UNKNOWN;
						flowpipes.tmv_flowpipes.push_back(flowpipe);

						if(checking_result == COMPLETED_SAFE)
						{
							checking_result = COMPLETED_UNKNOWN;
						}
					}
				}
				else
				{
					flowpipe.safety = SAFE;
					flowpipes.tmv_flowpipes.push_back(flowpipe);
				}

				current_flowpipe = new_flowpipe;
				break;
			}
			case CONTRACTED: 	// the domain is contracted but the time interval is not
			{
				bContracted = true;

				TaylorModelFlowpipe flowpipe;
				flowpipe.tmv_flowpipe = tmv_flowpipe;
				flowpipe.domain = contracted_domain;
				flowpipe.bConstrained = true;

				new_flowpipe.domain = contracted_domain;
				new_flowpipe.normalize(tm_setting.cutoff_threshold);

				if(bSafetyChecking)
				{
					int safety = safetyChecking(tmv_flowpipe, contracted_domain, safeSet, tm_setting, g_setting);

					if(safety == SAFE)
					{
						flowpipe.safety = SAFE;
						flowpipes.tmv_flowpipes.push_back(flowpipe);
					}
					else if(safety == UNSAFE)
					{
						flowpipe.safety = UNKNOWN;
						flowpipes.tmv_flowpipes.push_back(flowpipe);
						return COMPLETED_UNKNOWN;
					}
					else
					{
						flowpipe.safety = UNKNOWN;
						flowpipes.tmv_flowpipes.push_back(flowpipe);

						if(checking_result == COMPLETED_SAFE)
						{
							checking_result = COMPLETED_UNKNOWN;
						}
					}
				}
				else
				{
					flowpipe.safety = SAFE;
					flowpipes.tmv_flowpipes.push_back(flowpipe);
				}

				current_flowpipe = new_flowpipe;

				break;
			}
			case TIME_RANGE_CONTRACTED: 	// time interval is contracted
			{
				Real zero(0);

				if(contracted_domain[0].greaterThan(zero))
				{
					return checking_result;
				}
				else
				{
					bContracted = true;

					TaylorModelFlowpipe flowpipe;
					flowpipe.tmv_flowpipe = tmv_flowpipe;
					flowpipe.domain = contracted_domain;
					flowpipe.bConstrained = true;

					new_flowpipe.domain = contracted_domain;
					new_flowpipe.normalize(tm_setting.cutoff_threshold);

					if(bSafetyChecking)
					{
						int safety = safetyChecking(tmv_flowpipe, contracted_domain, safeSet, tm_setting, g_setting);

						if(safety == SAFE)
						{
							flowpipe.safety = SAFE;
							flowpipes.tmv_flowpipes.push_back(flowpipe);
						}
						else
						{
							flowpipe.safety = UNKNOWN;
							flowpipes.tmv_flowpipes.push_back(flowpipe);

							if(checking_result == COMPLETED_SAFE)
							{
								checking_result = COMPLETED_UNKNOWN;
							}
						}
					}
					else
					{
						flowpipe.safety = SAFE;
						flowpipes.tmv_flowpipes.push_back(flowpipe);
					}

					t += contracted_domain[0].sup();

					if(bPrint)
					{
						printf("time = %f,\t", t);
						printf("step = %f,\t", contracted_domain[0].sup());
						printf("order = %d\n", tm_setting.order);
					}

					last_flowpipe = new_flowpipe;

					return checking_result;
				}
			}
			}


			t += contracted_domain[0].sup();

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("order = %d\n", tm_setting.order);
			}

			if(tm_setting.order > tm_setting.order_min)
			{
				--tm_setting.order;
			}
		}
		else
		{
			switch(checking_result)
			{
			case COMPLETED_SAFE:
				return UNCOMPLETED_SAFE;
			case COMPLETED_UNSAFE:
				return UNCOMPLETED_UNSAFE;
			case COMPLETED_UNKNOWN:
				return UNCOMPLETED_UNKNOWN;
			}
		}
	}

	last_flowpipe = current_flowpipe;

	return checking_result;
}

template <class DATA_TYPE>
int ODE<DATA_TYPE>::reach_inv_symbolic_remainder(TaylorModelFlowpipes & flowpipes, Flowpipe & last_flowpipe, const double time, const Flowpipe & initialSet, const std::vector<Constraint> & invariant, Symbolic_Remainder & symbolic_remainder,
		const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & safeSet, const bool bSafetyChecking) const
{
	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	Flowpipe new_flowpipe, current_flowpipe = initialSet;
	bool bContracted = initialSet.bConstrained;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = current_flowpipe.advance(new_flowpipe, expressions, tm_setting, invariant, g_setting, symbolic_remainder);

		bool step_changed = false;
		double remaining_time = time - t;

		if(remaining_time < step)
		{
			step = remaining_time;
			new_flowpipe.domain[0].setSup(remaining_time);
			step_changed = true;
		}

		if(res == 1)
		{
			// constrain the flowpipe by the invariant
			TaylorModelVec<Real> tmv_flowpipe;

			if(step_changed)
			{
				new_flowpipe.compose(tmv_flowpipe, tm_setting.order, tm_setting.cutoff_threshold);
			}
			else
			{
				new_flowpipe.compose_normal(tmv_flowpipe, tm_setting.step_exp_table, tm_setting.order, tm_setting.cutoff_threshold);
			}

			std::vector<Interval> contracted_domain = new_flowpipe.domain;
			int type = domain_contraction_int(tmv_flowpipe, contracted_domain, invariant, tm_setting.order, tm_setting.cutoff_threshold, g_setting);

			switch(type)
			{
			case UNSAT:		// the intersection is empty
				last_flowpipe = current_flowpipe;
				return checking_result;
			case SAT:		// the domain is not contracted
			{
				TaylorModelFlowpipe flowpipe;
				flowpipe.tmv_flowpipe = tmv_flowpipe;
				flowpipe.domain = contracted_domain;
				flowpipe.bConstrained = bContracted;

				if(bSafetyChecking)
				{
					int safety = safetyChecking(tmv_flowpipe, contracted_domain, safeSet, tm_setting, g_setting);

					if(safety == SAFE)
					{
						flowpipe.safety = SAFE;
						flowpipes.tmv_flowpipes.push_back(flowpipe);
					}
					else if(safety == UNSAFE)
					{
						if(!bContracted)
						{
							flowpipe.safety = UNSAFE;
							flowpipes.tmv_flowpipes.push_back(flowpipe);
							return COMPLETED_UNSAFE;
						}
						else
						{
							flowpipe.safety = UNKNOWN;
							flowpipes.tmv_flowpipes.push_back(flowpipe);
							return COMPLETED_UNKNOWN;
						}
					}
					else
					{
						flowpipe.safety = UNKNOWN;
						flowpipes.tmv_flowpipes.push_back(flowpipe);

						if(checking_result == COMPLETED_SAFE)
						{
							checking_result = COMPLETED_UNKNOWN;
						}
					}
				}
				else
				{
					flowpipe.safety = SAFE;
					flowpipes.tmv_flowpipes.push_back(flowpipe);
				}

				current_flowpipe = new_flowpipe;
				break;
			}
			case CONTRACTED: 	// the domain is contracted but the time interval is not
			{
				bContracted = true;

				TaylorModelFlowpipe flowpipe;
				flowpipe.tmv_flowpipe = tmv_flowpipe;
				flowpipe.domain = contracted_domain;
				flowpipe.bConstrained = true;

				new_flowpipe.domain = contracted_domain;
				new_flowpipe.normalize(tm_setting.cutoff_threshold);

				symbolic_remainder.reset(new_flowpipe.tmvPre.tms.size());

				if(bSafetyChecking)
				{
					int safety = safetyChecking(tmv_flowpipe, contracted_domain, safeSet, tm_setting, g_setting);

					if(safety == SAFE)
					{
						flowpipe.safety = SAFE;
						flowpipes.tmv_flowpipes.push_back(flowpipe);
					}
					else if(safety == UNSAFE)
					{
						flowpipe.safety = UNKNOWN;
						flowpipes.tmv_flowpipes.push_back(flowpipe);
						return COMPLETED_UNKNOWN;
					}
					else
					{
						flowpipe.safety = UNKNOWN;
						flowpipes.tmv_flowpipes.push_back(flowpipe);

						if(checking_result == COMPLETED_SAFE)
						{
							checking_result = COMPLETED_UNKNOWN;
						}
					}
				}
				else
				{
					flowpipe.safety = SAFE;
					flowpipes.tmv_flowpipes.push_back(flowpipe);
				}

				current_flowpipe = new_flowpipe;

				break;
			}
			case TIME_RANGE_CONTRACTED: 	// time interval is contracted
			{
				Real zero(0);

				if(contracted_domain[0].greaterThan(zero))
				{
					return checking_result;
				}
				else
				{
					bContracted = true;

					TaylorModelFlowpipe flowpipe;
					flowpipe.tmv_flowpipe = tmv_flowpipe;
					flowpipe.domain = contracted_domain;
					flowpipe.bConstrained = true;

					new_flowpipe.domain = contracted_domain;
					new_flowpipe.normalize(tm_setting.cutoff_threshold);

					symbolic_remainder.reset(new_flowpipe.tmvPre.tms.size());

					if(bSafetyChecking)
					{
						int safety = safetyChecking(tmv_flowpipe, contracted_domain, safeSet, tm_setting, g_setting);

						if(safety == SAFE)
						{
							flowpipe.safety = SAFE;
							flowpipes.tmv_flowpipes.push_back(flowpipe);
						}
						else
						{
							flowpipe.safety = UNKNOWN;
							flowpipes.tmv_flowpipes.push_back(flowpipe);

							if(checking_result == COMPLETED_SAFE)
							{
								checking_result = COMPLETED_UNKNOWN;
							}
						}
					}
					else
					{
						flowpipe.safety = SAFE;
						flowpipes.tmv_flowpipes.push_back(flowpipe);
					}

					t += contracted_domain[0].sup();

					if(bPrint)
					{
						printf("time = %f,\t", t);
						printf("step = %f,\t", contracted_domain[0].sup());
						printf("order = %d\n", tm_setting.order);
					}

					last_flowpipe = new_flowpipe;

					return checking_result;
				}
			}
			}


			t += contracted_domain[0].sup();

			if(symbolic_remainder.J.size() >= symbolic_remainder.max_size)
			{
				symbolic_remainder.reset(current_flowpipe.tmvPre.tms.size());
			}

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("order = %d\n", tm_setting.order);
			}
		}
		else
		{
			switch(checking_result)
			{
			case COMPLETED_SAFE:
				return UNCOMPLETED_SAFE;
			case COMPLETED_UNSAFE:
				return UNCOMPLETED_UNSAFE;
			case COMPLETED_UNKNOWN:
				return UNCOMPLETED_UNKNOWN;
			}
		}
	}

	last_flowpipe = current_flowpipe;

	return checking_result;
}

template <class DATA_TYPE>
int ODE<DATA_TYPE>::reach_inv_symbolic_remainder_adaptive_stepsize(TaylorModelFlowpipes & flowpipes, Flowpipe & last_flowpipe, const double time, const Flowpipe & initialSet, const std::vector<Constraint> & invariant, Symbolic_Remainder & symbolic_remainder,
		Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & safeSet, const bool bSafetyChecking) const
{
	int checking_result = COMPLETED_SAFE;

	Flowpipe new_flowpipe, current_flowpipe = initialSet;
	bool bContracted = initialSet.bConstrained;
	double new_stepsize = -1;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = current_flowpipe.advance_adaptive_stepsize(new_flowpipe, expressions, new_stepsize, tm_setting, invariant, g_setting, symbolic_remainder);

		double current_stepsize = tm_setting.step_exp_table[1].sup();

		bool step_changed = false;
		double remaining_time = time - t;

		if(remaining_time < current_stepsize)
		{
			current_stepsize = remaining_time;
			new_flowpipe.domain[0].setSup(remaining_time);
			step_changed = true;
		}

		if(res == 1)
		{
			// constrain the flowpipe by the invariant
			TaylorModelVec<Real> tmv_flowpipe;

			if(step_changed)
			{
				new_flowpipe.compose(tmv_flowpipe, tm_setting.order, tm_setting.cutoff_threshold);
			}
			else
			{
				new_flowpipe.compose_normal(tmv_flowpipe, tm_setting.step_exp_table, tm_setting.order, tm_setting.cutoff_threshold);
			}

			std::vector<Interval> contracted_domain = new_flowpipe.domain;
			int type = domain_contraction_int(tmv_flowpipe, contracted_domain, invariant, tm_setting.order, tm_setting.cutoff_threshold, g_setting);

			switch(type)
			{
			case UNSAT:		// the intersection is empty
				last_flowpipe = current_flowpipe;
				return checking_result;
			case SAT:		// the domain is not contracted
			{
				TaylorModelFlowpipe flowpipe;
				flowpipe.tmv_flowpipe = tmv_flowpipe;
				flowpipe.domain = contracted_domain;
				flowpipe.bConstrained = bContracted;

				if(bSafetyChecking)
				{
					int safety = safetyChecking(tmv_flowpipe, contracted_domain, safeSet, tm_setting, g_setting);

					if(safety == SAFE)
					{
						flowpipe.safety = SAFE;
						flowpipes.tmv_flowpipes.push_back(flowpipe);
					}
					else if(safety == UNSAFE)
					{
						if(!bContracted)
						{
							flowpipe.safety = UNSAFE;
							flowpipes.tmv_flowpipes.push_back(flowpipe);
							return COMPLETED_UNSAFE;
						}
						else
						{
							flowpipe.safety = UNKNOWN;
							flowpipes.tmv_flowpipes.push_back(flowpipe);
							return COMPLETED_UNKNOWN;
						}
					}
					else
					{
						flowpipe.safety = UNKNOWN;
						flowpipes.tmv_flowpipes.push_back(flowpipe);

						if(checking_result == COMPLETED_SAFE)
						{
							checking_result = COMPLETED_UNKNOWN;
						}
					}
				}
				else
				{
					flowpipe.safety = SAFE;
					flowpipes.tmv_flowpipes.push_back(flowpipe);
				}

				current_flowpipe = new_flowpipe;
				break;
			}
			case CONTRACTED: 	// the domain is contracted but the time interval is not
			{
				bContracted = true;

				TaylorModelFlowpipe flowpipe;
				flowpipe.tmv_flowpipe = tmv_flowpipe;
				flowpipe.domain = contracted_domain;
				flowpipe.bConstrained = true;

				new_flowpipe.domain = contracted_domain;
				new_flowpipe.normalize(tm_setting.cutoff_threshold);

				symbolic_remainder.reset(new_flowpipe.tmvPre.tms.size());

				if(bSafetyChecking)
				{
					int safety = safetyChecking(tmv_flowpipe, contracted_domain, safeSet, tm_setting, g_setting);

					if(safety == SAFE)
					{
						flowpipe.safety = SAFE;
						flowpipes.tmv_flowpipes.push_back(flowpipe);
					}
					else if(safety == UNSAFE)
					{
						flowpipe.safety = UNKNOWN;
						flowpipes.tmv_flowpipes.push_back(flowpipe);
						return COMPLETED_UNKNOWN;
					}
					else
					{
						flowpipe.safety = UNKNOWN;
						flowpipes.tmv_flowpipes.push_back(flowpipe);

						if(checking_result == COMPLETED_SAFE)
						{
							checking_result = COMPLETED_UNKNOWN;
						}
					}
				}
				else
				{
					flowpipe.safety = SAFE;
					flowpipes.tmv_flowpipes.push_back(flowpipe);
				}

				current_flowpipe = new_flowpipe;

				break;
			}
			case TIME_RANGE_CONTRACTED: 	// time interval is contracted
			{
				Real zero(0);

				if(contracted_domain[0].greaterThan(zero))
				{
					return checking_result;
				}
				else
				{
					bContracted = true;

					TaylorModelFlowpipe flowpipe;
					flowpipe.tmv_flowpipe = tmv_flowpipe;
					flowpipe.domain = contracted_domain;
					flowpipe.bConstrained = true;

					new_flowpipe.domain = contracted_domain;
					new_flowpipe.normalize(tm_setting.cutoff_threshold);

					symbolic_remainder.reset(new_flowpipe.tmvPre.tms.size());

					if(bSafetyChecking)
					{
						int safety = safetyChecking(tmv_flowpipe, contracted_domain, safeSet, tm_setting, g_setting);

						if(safety == SAFE)
						{
							flowpipe.safety = SAFE;
							flowpipes.tmv_flowpipes.push_back(flowpipe);
						}
						else
						{
							flowpipe.safety = UNKNOWN;
							flowpipes.tmv_flowpipes.push_back(flowpipe);

							if(checking_result == COMPLETED_SAFE)
							{
								checking_result = COMPLETED_UNKNOWN;
							}
						}
					}
					else
					{
						flowpipe.safety = SAFE;
						flowpipes.tmv_flowpipes.push_back(flowpipe);
					}

					t += contracted_domain[0].sup();

					if(bPrint)
					{
						printf("time = %f,\t", t);
						printf("step = %f,\t", contracted_domain[0].sup());
						printf("order = %d\n", tm_setting.order);
					}

					last_flowpipe = new_flowpipe;

					return checking_result;
				}
			}
			}

			t += current_stepsize;

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", current_stepsize);
				printf("order = %d\n", tm_setting.order);
			}

			new_stepsize = current_stepsize * LAMBDA_UP;
			if(new_stepsize > tm_setting.step_max - THRESHOLD_HIGH)
			{
				new_stepsize = -1;
			}

			if(symbolic_remainder.J.size() >= symbolic_remainder.max_size)
			{
				symbolic_remainder.reset(current_flowpipe.tmvPre.tms.size());
			}
		}
		else
		{
			switch(checking_result)
			{
			case COMPLETED_SAFE:
				return UNCOMPLETED_SAFE;
			case COMPLETED_UNSAFE:
				return UNCOMPLETED_UNSAFE;
			case COMPLETED_UNKNOWN:
				return UNCOMPLETED_UNKNOWN;
			}
		}
	}

	last_flowpipe = current_flowpipe;

	return checking_result;
}

template <class DATA_TYPE>
int ODE<DATA_TYPE>::reach_inv_symbolic_remainder_adaptive_order(TaylorModelFlowpipes & flowpipes, Flowpipe & last_flowpipe, const double time, const Flowpipe & initialSet, const std::vector<Constraint> & invariant, Symbolic_Remainder & symbolic_remainder,
		Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & safeSet, const bool bSafetyChecking) const
{
	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	Flowpipe new_flowpipe, current_flowpipe = initialSet;
	bool bContracted = initialSet.bConstrained;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = current_flowpipe.advance_adaptive_order(new_flowpipe, expressions, tm_setting, invariant, g_setting, symbolic_remainder);

		bool step_changed = false;
		double remaining_time = time - t;

		if(remaining_time < step)
		{
			step = remaining_time;
			new_flowpipe.domain[0].setSup(remaining_time);
			step_changed = true;
		}

		if(res == 1)
		{
			// constrain the flowpipe by the invariant
			TaylorModelVec<Real> tmv_flowpipe;

			if(step_changed)
			{
				new_flowpipe.compose(tmv_flowpipe, tm_setting.order, tm_setting.cutoff_threshold);
			}
			else
			{
				new_flowpipe.compose_normal(tmv_flowpipe, tm_setting.step_exp_table, tm_setting.order, tm_setting.cutoff_threshold);
			}

			std::vector<Interval> contracted_domain = new_flowpipe.domain;
			int type = domain_contraction_int(tmv_flowpipe, contracted_domain, invariant, tm_setting.order, tm_setting.cutoff_threshold, g_setting);

			switch(type)
			{
			case UNSAT:		// the intersection is empty
				last_flowpipe = current_flowpipe;
				return checking_result;
			case SAT:		// the domain is not contracted
			{
				TaylorModelFlowpipe flowpipe;
				flowpipe.tmv_flowpipe = tmv_flowpipe;
				flowpipe.domain = contracted_domain;
				flowpipe.bConstrained = bContracted;

				if(bSafetyChecking)
				{
					int safety = safetyChecking(tmv_flowpipe, contracted_domain, safeSet, tm_setting, g_setting);

					if(safety == SAFE)
					{
						flowpipe.safety = SAFE;
						flowpipes.tmv_flowpipes.push_back(flowpipe);
					}
					else if(safety == UNSAFE)
					{
						if(!bContracted)
						{
							flowpipe.safety = UNSAFE;
							flowpipes.tmv_flowpipes.push_back(flowpipe);
							return COMPLETED_UNSAFE;
						}
						else
						{
							flowpipe.safety = UNKNOWN;
							flowpipes.tmv_flowpipes.push_back(flowpipe);
							return COMPLETED_UNKNOWN;
						}
					}
					else
					{
						flowpipe.safety = UNKNOWN;
						flowpipes.tmv_flowpipes.push_back(flowpipe);

						if(checking_result == COMPLETED_SAFE)
						{
							checking_result = COMPLETED_UNKNOWN;
						}
					}
				}
				else
				{
					flowpipe.safety = SAFE;
					flowpipes.tmv_flowpipes.push_back(flowpipe);
				}

				current_flowpipe = new_flowpipe;
				break;
			}
			case CONTRACTED: 	// the domain is contracted but the time interval is not
			{
				bContracted = true;

				TaylorModelFlowpipe flowpipe;
				flowpipe.tmv_flowpipe = tmv_flowpipe;
				flowpipe.domain = contracted_domain;
				flowpipe.bConstrained = true;

				new_flowpipe.domain = contracted_domain;
				new_flowpipe.normalize(tm_setting.cutoff_threshold);

				symbolic_remainder.reset(new_flowpipe.tmvPre.tms.size());

				if(bSafetyChecking)
				{
					int safety = safetyChecking(tmv_flowpipe, contracted_domain, safeSet, tm_setting, g_setting);

					if(safety == SAFE)
					{
						flowpipe.safety = SAFE;
						flowpipes.tmv_flowpipes.push_back(flowpipe);
					}
					else if(safety == UNSAFE)
					{
						flowpipe.safety = UNKNOWN;
						flowpipes.tmv_flowpipes.push_back(flowpipe);
						return COMPLETED_UNKNOWN;
					}
					else
					{
						flowpipe.safety = UNKNOWN;
						flowpipes.tmv_flowpipes.push_back(flowpipe);

						if(checking_result == COMPLETED_SAFE)
						{
							checking_result = COMPLETED_UNKNOWN;
						}
					}
				}
				else
				{
					flowpipe.safety = SAFE;
					flowpipes.tmv_flowpipes.push_back(flowpipe);
				}

				current_flowpipe = new_flowpipe;

				break;
			}
			case TIME_RANGE_CONTRACTED: 	// time interval is contracted
			{
				Real zero(0);

				if(contracted_domain[0].greaterThan(zero))
				{
					return checking_result;
				}
				else
				{
					bContracted = true;

					TaylorModelFlowpipe flowpipe;
					flowpipe.tmv_flowpipe = tmv_flowpipe;
					flowpipe.domain = contracted_domain;
					flowpipe.bConstrained = true;

					new_flowpipe.domain = contracted_domain;
					new_flowpipe.normalize(tm_setting.cutoff_threshold);

					symbolic_remainder.reset(new_flowpipe.tmvPre.tms.size());

					if(bSafetyChecking)
					{
						int safety = safetyChecking(tmv_flowpipe, contracted_domain, safeSet, tm_setting, g_setting);

						if(safety == SAFE)
						{
							flowpipe.safety = SAFE;
							flowpipes.tmv_flowpipes.push_back(flowpipe);
						}
						else
						{
							flowpipe.safety = UNKNOWN;
							flowpipes.tmv_flowpipes.push_back(flowpipe);

							if(checking_result == COMPLETED_SAFE)
							{
								checking_result = COMPLETED_UNKNOWN;
							}
						}
					}
					else
					{
						flowpipe.safety = SAFE;
						flowpipes.tmv_flowpipes.push_back(flowpipe);
					}

					t += contracted_domain[0].sup();

					if(bPrint)
					{
						printf("time = %f,\t", t);
						printf("step = %f,\t", contracted_domain[0].sup());
						printf("order = %d\n", tm_setting.order);
					}

					last_flowpipe = new_flowpipe;

					return checking_result;
				}
			}
			}


			t += contracted_domain[0].sup();

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", step);
				printf("order = %d\n", tm_setting.order);
			}

			if(tm_setting.order > tm_setting.order_min)
			{
				--tm_setting.order;
			}

			if(symbolic_remainder.J.size() >= symbolic_remainder.max_size)
			{
				symbolic_remainder.reset(current_flowpipe.tmvPre.tms.size());
			}
		}
		else
		{
			switch(checking_result)
			{
			case COMPLETED_SAFE:
				return UNCOMPLETED_SAFE;
			case COMPLETED_UNSAFE:
				return UNCOMPLETED_UNSAFE;
			case COMPLETED_UNKNOWN:
				return UNCOMPLETED_UNKNOWN;
			}
		}
	}

	last_flowpipe = current_flowpipe;

	return checking_result;
}

template <class DATA_TYPE>
void ODE<DATA_TYPE>::reach_inv(Result_of_Reachability & result, const Flowpipe & initialSet, const double T, const std::vector<Constraint> & invariant, Computational_Setting & setting, const std::vector<Constraint> & safeSet) const
{
	if(T <= 0)
	{
		printf("Time horizon should be positive.\n");
		return;
	}

	bool bSafetyChecking = false;

	if(safeSet.size() > 0)
	{
		bSafetyChecking = true;
	}

	if(setting.tm_setting.step_min > 0)
	{
		// adaptive stepsizes
		result.status = reach_inv_adaptive_stepsize(result.tmv_flowpipes, result.fp_end_of_time, T, initialSet, invariant, setting.tm_setting, setting.g_setting,
				setting.bPrint, safeSet, bSafetyChecking);
	}
	else if(setting.tm_setting.order_max > 0)
	{
		// adaptive orders
		result.status = reach_inv_adaptive_order(result.tmv_flowpipes, result.fp_end_of_time, T, initialSet, invariant, setting.tm_setting, setting.g_setting,
				setting.bPrint, safeSet, bSafetyChecking);
	}
	else
	{
		// fixed stepsizes and orders
		result.status = reach_inv(result.tmv_flowpipes, result.fp_end_of_time, T, initialSet, invariant, setting.tm_setting, setting.g_setting,
				setting.bPrint, safeSet, bSafetyChecking);
	}

	if(result.tmv_flowpipes.size() > 0)
	{
		TaylorModelVec<Real> last_flowpipe = result.tmv_flowpipes.tmv_flowpipes.back().tmv_flowpipe;

		Real t;
		(result.tmv_flowpipes.tmv_flowpipes.back().domain)[0].sup(t);

		std::vector<Real> realVec;
		realVec.push_back(1);
		realVec.push_back(t);

		Real tmp = t;

		for(unsigned int i=2; i<=setting.tm_setting.step_end_exp_table.size(); ++i)
		{
			tmp *= t;
			realVec.push_back(tmp);
		}

		last_flowpipe.evaluate_time(result.tmv_fp_end_of_time.tmv_flowpipe, realVec);

		result.fp_end_of_time.tmvPre.evaluate_time(result.fp_end_of_time.tmvPre, realVec);
		result.fp_end_of_time.domain[0] = 0;

		result.tmv_fp_end_of_time.domain = result.fp_end_of_time.domain;
	}
}

template <class DATA_TYPE>
void ODE<DATA_TYPE>::reach_inv(Result_of_Reachability & result, const Flowpipe & initialSet, const double T, const std::vector<Constraint> & invariant, Computational_Setting & setting, const std::vector<Constraint> & safeSet, Symbolic_Remainder & symbolic_remainder) const
{
	if(T <= 0)
	{
		printf("Time horizon should be positive.\n");
		return;
	}

	bool bSafetyChecking = false;

	if(safeSet.size() > 0)
	{
		bSafetyChecking = true;
	}

	if(setting.tm_setting.step_min > 0)
	{
		// adaptive stepsizes
		result.status = reach_inv_symbolic_remainder_adaptive_stepsize(result.tmv_flowpipes, result.fp_end_of_time, T, initialSet, invariant, symbolic_remainder, setting.tm_setting, setting.g_setting,
				setting.bPrint, safeSet, bSafetyChecking);
	}
	else if(setting.tm_setting.order_max > 0)
	{
		// adaptive orders
		result.status = reach_inv_symbolic_remainder_adaptive_order(result.tmv_flowpipes, result.fp_end_of_time, T, initialSet, invariant, symbolic_remainder, setting.tm_setting, setting.g_setting,
				setting.bPrint, safeSet, bSafetyChecking);
	}
	else
	{
		// fixed stepsizes and orders
		result.status = reach_inv_symbolic_remainder(result.tmv_flowpipes, result.fp_end_of_time, T, initialSet, invariant, symbolic_remainder, setting.tm_setting, setting.g_setting,
				setting.bPrint, safeSet, bSafetyChecking);
	}

	if(result.tmv_flowpipes.size() > 0)
	{
		TaylorModelVec<Real> last_flowpipe = result.tmv_flowpipes.tmv_flowpipes.back().tmv_flowpipe;

		Real t;
		(result.tmv_flowpipes.tmv_flowpipes.back().domain)[0].sup(t);

		std::vector<Real> realVec;
		realVec.push_back(1);
		realVec.push_back(t);

		Real tmp = t;

		for(unsigned int i=2; i<=setting.tm_setting.step_end_exp_table.size(); ++i)
		{
			tmp *= t;
			realVec.push_back(tmp);
		}

		last_flowpipe.evaluate_time(result.tmv_fp_end_of_time.tmv_flowpipe, realVec);

		result.fp_end_of_time.tmvPre.evaluate_time(result.fp_end_of_time.tmvPre, realVec);
		result.fp_end_of_time.domain[0] = 0;

		result.tmv_fp_end_of_time.domain = result.fp_end_of_time.domain;
	}
}














class Plot_Setting
{
public:
	Variables variables;
	std::vector<Expression<Real> > outputDims;
	std::vector<std::string> labels;
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

//	void plot_2D(const std::string & path, const std::string & fileName, const TaylorModelFlowpipes & flowpipes, Computational_Setting & setting) const;

	void plot_2D_MATLAB(const std::string & path, const std::string & fileName, const TaylorModelFlowpipes & flowpipes, Computational_Setting & setting) const;
	void plot_2D_interval_MATLAB(const std::string & path, const std::string & fileName, const TaylorModelFlowpipes & flowpipes, Computational_Setting & setting) const;
	void plot_2D_octagon_MATLAB(const std::string & path, const std::string & fileName, const TaylorModelFlowpipes & flowpipes, Computational_Setting & setting) const;
	void plot_2D_grids_MATLAB(const std::string & path, const std::string & fileName, const unsigned int num, const TaylorModelFlowpipes & flowpipes, Computational_Setting & setting) const;

	void plot_2D_GNUPLOT(const std::string & path, const std::string & fileName, const TaylorModelFlowpipes & flowpipes, Computational_Setting & setting) const;
	void plot_2D_interval_GNUPLOT(const std::string & path, const std::string & fileName, const TaylorModelFlowpipes & flowpipes, Computational_Setting & setting) const;
	void plot_2D_octagon_GNUPLOT(const std::string & path, const std::string & fileName, const TaylorModelFlowpipes & flowpipes, Computational_Setting & setting) const;
	void plot_2D_grids_GNUPLOT(const std::string & path, const std::string & fileName, const unsigned int num, const TaylorModelFlowpipes & flowpipes, Computational_Setting & setting) const;
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
	std::vector<TaylorModelFlowpipes> flowpipes;					// contracted flowpipes
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





void gridBox(std::list<std::vector<Interval> > & grids, const std::vector<Interval> & box, const unsigned int num);

unsigned int findProperOrder(Real & error, const Real & max, const Real & min, const Real & tolerance, const unsigned int start_order);

void check_connectivities(Matrix<bool> & result, Matrix<bool> & adjMatrix);

void compute_one_step_trans(Matrix<UnivariateTaylorModel<Real> > & utm_Phi_t, Matrix<UnivariateTaylorModel<Real> > & utm_Phi_end, Matrix<UnivariateTaylorModel<Real> > & utm_Psi_t,
		Matrix<UnivariateTaylorModel<Real> > & utm_Psi_end, Matrix<UnivariateTaylorModel<Real> > & utm_Omega_t, Matrix<UnivariateTaylorModel<Real> > & utm_Omega_end, Matrix<Interval> & tv_part,
		LTV_ODE & ltv_ode, const UnivariateTaylorModel<Real> & utm_t0, const unsigned int order, const double step, const Global_Setting & g_setting);


// flowpipe/guard intersections
void intersect_a_guard(Intersected_Flowpipes & result, const TaylorModelFlowpipes & flowpipes, const std::vector<Constraint> & guard, const bool boundary_of_invariant, const Computational_Setting & setting);


void compute_weight_matrix(Matrix<double> & weightMat, const std::vector<std::vector<double> > & candidate_vectors);
bool check_validity(Matrix<double> & matTemplate, const std::vector<double> & vec, const int rank);
bool select_a_vector(FactorTab & lst_selected, std::list<FactorTab> & lst_unselected, Matrix<double> & matTemplate, const std::vector<std::vector<double> > & candidate_vectors, int & rank);



void eliminate_t(TaylorModelFlowpipe & flowpipe);

void create_initial_set(Flowpipe & initial_set, const TaylorModelFlowpipe & flowpipe);

void merge_consecutive_flowpipes(Flowpipe & result, const TaylorModelFlowpipes & flowpipes, LTI_ODE & lti_ode,
		const std::vector<Constraint> & invariant, const Computational_Setting & setting);





// only for testing
//void test_domain_contraction(Result_of_Reachability & contraction_result, Result_of_Reachability & reachability_result, const std::vector<Constraint> & constraints, const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting);

}

#endif /* CONTINUOUS_H_ */
