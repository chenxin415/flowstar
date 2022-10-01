/*---
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.
  The code is released as is under the GNU General Public License (GPL).
---*/

#ifndef DISCRETE_H_
#define DISCRETE_H_


#include "Continuous.h"



namespace flowstar
{


/*
 * class of linear time-invariant DDEs: x_{k+1} = A x_k + B + C p + D v where
 * x are the state variables, p are the time-invariant parameters
 * B is a constant matrix
 * D is the constant matrix for the time-varying uncertainties in [-1,1]
 */

class LTI_DDE : public Dynamics
{
protected:
	Matrix<Real> rm_dyn_A;
	Matrix<Real> rm_dyn_B;
	Matrix<Real> rm_dyn_C;
	Matrix<Real> rm_dyn_D;

public:
	LTI_DDE(const Matrix<Real> & A);
	LTI_DDE(const Matrix<Real> & A, const Matrix<Real> & B);
	LTI_DDE(const Matrix<Real> & A, const Matrix<Real> & B, const Matrix<Real> & C);
	LTI_DDE(const Matrix<Real> & A, const Matrix<Real> & B, const Matrix<Real> & C, const Matrix<Real> & D);
	LTI_DDE(const LTI_DDE & lti_dde);
	~LTI_DDE();

	LTI_DDE & operator = (const LTI_DDE & lti_dde);

	void abstract(FlowmapAbstraction & abstraction, const int N, const int zono_order = -1);
};










/*
 * Class for the dynamics defined by Discrete Difference Equation (DDE)
 * which may have interval coefficients representing time-varying uncertainties.
 */

template <class DATA_TYPE>
class DDE : public Dynamics
{
protected:
	Variables stateVars;
	std::vector<Expression<DATA_TYPE> >	expressions;

public:
	DDE(const Variables & vars);
	DDE(const std::vector<std::string> & str_dde, Variables & vars);
	DDE(const DDE<DATA_TYPE> & dde);
	virtual ~DDE() {} ;

	DDE<DATA_TYPE> & operator = (const DDE<DATA_TYPE> & dde);

	bool defNextState(const std::string & varName, const Expression<DATA_TYPE> & nextState);

protected:
	int reach(std::list<Flowpipe> & flowpipes, const unsigned int n, const Flowpipe & initialSet, const Taylor_Model_Setting & tm_setting,
			const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking) const;

	int reach_inv(TaylorModelFlowpipes & flowpipes, Flowpipe & last_flowpipe, const unsigned int n, const Flowpipe & initialSet, const std::vector<Constraint> & invariant, const Taylor_Model_Setting & tm_setting,
			const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking) const;


	int reach_symbolic_remainder(std::list<Flowpipe> & flowpipes, const unsigned int n, const Flowpipe & initialSet, Symbolic_Remainder & symbolic_remainder, const Taylor_Model_Setting & tm_setting,
			const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking) const;

	int reach_inv_symbolic_remainder(TaylorModelFlowpipes & flowpipes, Flowpipe & last_flowpipe, const unsigned int n, const Flowpipe & initialSet, const std::vector<Constraint> & invariant, Symbolic_Remainder & symbolic_remainder, const Taylor_Model_Setting & tm_setting,
			const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking) const;


public:
	void reach(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const unsigned int n, const std::vector<Constraint> & unsafeSet) const;
	void reach_inv(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const unsigned int n, const std::vector<Constraint> & invariant, const std::vector<Constraint> & unsafeSet) const;

	void reach(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const unsigned int n, const std::vector<Constraint> & unsafeSet, Symbolic_Remainder & symbolic_remainder) const;
	void reach_inv(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const unsigned int n, const std::vector<Constraint> & invariant, const std::vector<Constraint> & unsafeSet, Symbolic_Remainder & symbolic_remainder) const;
};



template <class DATA_TYPE>
DDE<DATA_TYPE>::DDE(const Variables & vars)
{
	stateVars = vars;
}

template <class DATA_TYPE>
DDE<DATA_TYPE>::DDE(const std::vector<std::string> & str_dde, Variables & vars)
{
	if(str_dde.size() > vars.size())
	{
		printf("DDE: There are more derivatives than the state variables.\n");

		Expression<DATA_TYPE> zero(DATA_TYPE(0));
		expressions.resize(vars.size(), zero);
	}
	else
	{
		Expression<DATA_TYPE> zero(DATA_TYPE(0));
		expressions.resize(vars.size(), zero);

		for(int i=0; i<str_dde.size(); ++i)
		{
			Expression<DATA_TYPE> next_state(str_dde[i], vars);
			expressions[i] = next_state;
		}
	}
}

template <class DATA_TYPE>
DDE<DATA_TYPE>::DDE(const DDE<DATA_TYPE> & dde)
{
	stateVars	= dde.stateVars;
	expressions	= dde.expressions;
}

template <class DATA_TYPE>
DDE<DATA_TYPE> & DDE<DATA_TYPE>::operator = (const DDE<DATA_TYPE> & dde)
{
	if(this == &dde)
		return *this;

	stateVars	= dde.stateVars;
	expressions	= dde.expressions;

	return *this;
}

template <class DATA_TYPE>
bool DDE<DATA_TYPE>::defNextState(const std::string & varName, const Expression<DATA_TYPE> & nextState)
{
	int varID = stateVars.getIDForVar(varName);

	if(varID < 0)
	{
		printf("State variable %s is not declared.\n", varName.c_str());
		return false;
	}

	expressions[varID] = nextState;

	return true;
}

template <class DATA_TYPE>
int DDE<DATA_TYPE>::reach(std::list<Flowpipe> & flowpipes, const unsigned int n, const Flowpipe & initialSet, const Taylor_Model_Setting & tm_setting,
		const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking) const
{
	int checking_result = COMPLETED_SAFE;
	unsigned int rangeDim = expressions.size();
	Interval intZero, intUnit(-1,1);

	flowpipes.push_back(initialSet);

	Flowpipe current_flowpipe = initialSet;


	for(int i=0; i<n; ++i)
	{
		Flowpipe new_flowpipe;
		new_flowpipe.domain = current_flowpipe.domain;

		TaylorModelVec<Real> tmv_of_x0 = current_flowpipe.tmvPre;

		// the center point of x0's polynomial part
		std::vector<Real> const_of_x0;
		tmv_of_x0.constant(const_of_x0);
/*
		for(unsigned int i=0; i<rangeDim; ++i)
		{
			Real c;
			tmv_of_x0.tms[i].remainder.remove_midpoint(c);
			const_of_x0[i] += c;
		}
*/
		TaylorModelVec<Real> tmv_c0(const_of_x0, rangeDim + 1);

		// introduce a new variable r0 such that x0 = c0 + A*r0, then r0 is origin-centered
		tmv_of_x0.rmConstant();

		std::vector<Interval> range_of_x0;

		std::vector<Interval> tmvPolyRange;
		current_flowpipe.tmv.polyRange(tmvPolyRange, current_flowpipe.domain);
		tmv_of_x0.insert_ctrunc_normal(new_flowpipe.tmv, current_flowpipe.tmv, tmvPolyRange, tm_setting.step_end_exp_table, rangeDim + 1, tm_setting.order, tm_setting.cutoff_threshold);

		new_flowpipe.tmv.intEvalNormal(range_of_x0, tm_setting.step_end_exp_table);



		// Compute the scaling matrix S.
		std::vector<Real> S, invS;

		for(int i=0; i<rangeDim; ++i)
		{
			Real sup;
			range_of_x0[i].mag(sup);

			if(sup == 0)
			{
				S.push_back(0);
				invS.push_back(1);
			}
			else
			{
				S.push_back(sup);
				Real tmp = 1/sup;
				invS.push_back(tmp);
				range_of_x0[i] = intUnit;
			}
		}

		range_of_x0.insert(range_of_x0.begin(), intZero);

		new_flowpipe.tmv.scale_assign(invS);

		TaylorModelVec<Real> new_x0(S);
		new_x0 += tmv_c0;
		TaylorModelVec<Real> x = new_x0;

		for(int i=0; i<expressions.size(); ++i)
		{
			TaylorModel<Real> tmTemp;
			expressions[i].evaluate(tmTemp, x.tms, tm_setting.order, range_of_x0, tm_setting.cutoff_threshold, g_setting);
			new_flowpipe.tmvPre.tms.push_back(tmTemp);
		}


		if(bPrint)
		{
			printf("Step %d\n", i+1);
		}


		if(bSafetyChecking)
		{
			int safety = new_flowpipe.safetyChecking(unsafeSet, tm_setting, g_setting);

			flowpipes.push_back(new_flowpipe);

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
			flowpipes.push_back(new_flowpipe);
		}

		current_flowpipe = new_flowpipe;
	}

	return checking_result;
}

template <class DATA_TYPE>
int DDE<DATA_TYPE>::reach_inv(TaylorModelFlowpipes & flowpipes, Flowpipe & last_flowpipe, const unsigned int n, const Flowpipe & initialSet, const std::vector<Constraint> & invariant, const Taylor_Model_Setting & tm_setting,
		const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking) const
{
	int checking_result = COMPLETED_SAFE;
	unsigned int rangeDim = expressions.size();
	Interval intZero, intUnit(-1,1);

	Flowpipe current_flowpipe = initialSet;

	TaylorModelVec<Real> tmv_flowpipe;
	initialSet.compose(tmv_flowpipe, tm_setting.order, tm_setting.cutoff_threshold);

	TaylorModelFlowpipe init_flowpipe;
	init_flowpipe.tmv_flowpipe = tmv_flowpipe;
	init_flowpipe.domain = initialSet.domain;
	flowpipes.tmv_flowpipes.push_back(init_flowpipe);

	bool bContracted = false;

	for(int i=0; i<n; ++i)
	{
		Flowpipe new_flowpipe;
		new_flowpipe.domain = current_flowpipe.domain;

		TaylorModelVec<Real> tmv_of_x0 = current_flowpipe.tmvPre;

		// the center point of x0's polynomial part
		std::vector<Real> const_of_x0;
		tmv_of_x0.constant(const_of_x0);
/*
		for(unsigned int i=0; i<rangeDim; ++i)
		{
			Real c;
			tmv_of_x0.tms[i].remainder.remove_midpoint(c);
			const_of_x0[i] += c;
		}
*/
		TaylorModelVec<Real> tmv_c0(const_of_x0, rangeDim + 1);

		// introduce a new variable r0 such that x0 = c0 + A*r0, then r0 is origin-centered
		tmv_of_x0.rmConstant();

		std::vector<Interval> range_of_x0;

		std::vector<Interval> tmvPolyRange;
		current_flowpipe.tmv.polyRange(tmvPolyRange, current_flowpipe.domain);
		tmv_of_x0.insert_ctrunc_normal(new_flowpipe.tmv, current_flowpipe.tmv, tmvPolyRange, tm_setting.step_end_exp_table, rangeDim + 1, tm_setting.order, tm_setting.cutoff_threshold);

		new_flowpipe.tmv.intEvalNormal(range_of_x0, tm_setting.step_end_exp_table);



		// Compute the scaling matrix S.
		std::vector<Real> S, invS;

		for(int i=0; i<rangeDim; ++i)
		{
			Real sup;
			range_of_x0[i].mag(sup);

			if(sup == 0)
			{
				S.push_back(0);
				invS.push_back(1);
			}
			else
			{
				S.push_back(sup);
				Real tmp = 1/sup;
				invS.push_back(tmp);
				range_of_x0[i] = intUnit;
			}
		}

		range_of_x0.insert(range_of_x0.begin(), intZero);

		new_flowpipe.tmv.scale_assign(invS);

		TaylorModelVec<Real> new_x0(S);
		new_x0 += tmv_c0;
		TaylorModelVec<Real> x = new_x0;

		for(int i=0; i<expressions.size(); ++i)
		{
			TaylorModel<Real> tmTemp;
			expressions[i].evaluate(tmTemp, x.tms, tm_setting.order, range_of_x0, tm_setting.cutoff_threshold, g_setting);
			new_flowpipe.tmvPre.tms.push_back(tmTemp);
		}



		// constrain the flowpipe by the invariant
		TaylorModelVec<Real> tmv_flowpipe;
		new_flowpipe.compose_normal(tmv_flowpipe, tm_setting.step_exp_table, tm_setting.order, tm_setting.cutoff_threshold);

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
				int safety = safetyChecking(tmv_flowpipe, contracted_domain, unsafeSet, tm_setting, g_setting);

				if(safety == SAFE)
				{
					flowpipe.safety = SAFE;
					flowpipes.tmv_flowpipes.push_back(flowpipe);
				}
				else if(safety == UNSAFE && !bContracted)
				{
					flowpipe.safety = UNSAFE;
					flowpipes.tmv_flowpipes.push_back(flowpipe);
					return COMPLETED_UNSAFE;
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
				int safety = safetyChecking(tmv_flowpipe, contracted_domain, unsafeSet, tm_setting, g_setting);

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
		}

		if(bPrint)
		{
			printf("Step %d\n", i+1);
		}
	}

	last_flowpipe = current_flowpipe;

	return checking_result;
}

template <class DATA_TYPE>
int DDE<DATA_TYPE>::reach_symbolic_remainder(std::list<Flowpipe> & flowpipes, const unsigned int n, const Flowpipe & initialSet, Symbolic_Remainder & symbolic_remainder, const Taylor_Model_Setting & tm_setting,
		const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking) const
{
	int checking_result = COMPLETED_SAFE;
	unsigned int rangeDim = expressions.size();
	Interval intZero, intUnit(-1,1);

	flowpipes.push_back(initialSet);

	Flowpipe current_flowpipe = initialSet;

	for(int k=0; k<n; ++k)
	{
		Flowpipe new_flowpipe;
		new_flowpipe.domain = current_flowpipe.domain;

		TaylorModelVec<Real> tmv_of_x0 = current_flowpipe.tmvPre;

		// the center point of x0's polynomial part
		std::vector<Real> const_of_x0;
		tmv_of_x0.constant(const_of_x0);
/*
		for(unsigned int i=0; i<rangeDim; ++i)
		{
			Real c;
			tmv_of_x0.tms[i].remainder.remove_midpoint(c);
			const_of_x0[i] += c;
		}
*/
		TaylorModelVec<Real> tmv_c0(const_of_x0, rangeDim + 1);

		// introduce a new variable r0 such that x0 = c0 + A*r0, then r0 is origin-centered
		tmv_of_x0.rmConstant();


		// decompose the linear and nonlinear part
		TaylorModelVec<Real> x0_linear, x0_other;
		tmv_of_x0.decompose(x0_linear, x0_other);

		Matrix<Real> Phi_L_i(rangeDim, rangeDim);

		x0_linear.linearCoefficients(Phi_L_i);

		Matrix<Real> local_trans_linear = Phi_L_i;

		Phi_L_i.right_scale_assign(symbolic_remainder.scalars);


		// compute the remainder part under the linear transformation
		Matrix<Interval> J_i(rangeDim, 1);

		for(unsigned int i=1; i<symbolic_remainder.Phi_L.size(); ++i)
		{
			symbolic_remainder.Phi_L[i] = Phi_L_i * symbolic_remainder.Phi_L[i];
		}

		symbolic_remainder.Phi_L.push_back(Phi_L_i);

		for(unsigned int i=1; i<symbolic_remainder.Phi_L.size(); ++i)
		{
			J_i += symbolic_remainder.Phi_L[i] * symbolic_remainder.J[i-1];
		}

		Matrix<Interval> J_ip1(rangeDim, 1);

		std::vector<Interval> range_of_x0;


		// compute the local initial set
		if(symbolic_remainder.J.size() > 0)
		{
			// compute the other part
			std::vector<Interval> tmvPolyRange;
			current_flowpipe.tmv.polyRange(tmvPolyRange, current_flowpipe.domain);
			x0_other.insert_ctrunc_normal(new_flowpipe.tmv, current_flowpipe.tmv, tmvPolyRange, tm_setting.step_end_exp_table, rangeDim + 1, tm_setting.order, tm_setting.cutoff_threshold);

			new_flowpipe.tmv.Remainder(J_ip1);

			std::vector<Polynomial<Real> > poly_tmv;
			current_flowpipe.tmv.Expansion(poly_tmv);
			std::vector<Polynomial<Real> > linear_part = local_trans_linear * poly_tmv;

			for(int i=0; i<rangeDim; ++i)
			{
				new_flowpipe.tmv.tms[i].expansion += linear_part[i];
			}

			for(int i=0; i<rangeDim; ++i)
			{
				new_flowpipe.tmv.tms[i].remainder = J_ip1[i][0] + J_i[i][0];
			}

			new_flowpipe.tmv.intEvalNormal(range_of_x0, tm_setting.step_end_exp_table);
		}
		else
		{
			std::vector<Interval> tmvPolyRange;
			current_flowpipe.tmv.polyRange(tmvPolyRange, current_flowpipe.domain);
			tmv_of_x0.insert_ctrunc_normal(new_flowpipe.tmv, current_flowpipe.tmv, tmvPolyRange, tm_setting.step_end_exp_table, rangeDim + 1, tm_setting.order, tm_setting.cutoff_threshold);

			new_flowpipe.tmv.Remainder(J_ip1);

			new_flowpipe.tmv.intEvalNormal(range_of_x0, tm_setting.step_end_exp_table);
		}

		symbolic_remainder.J.push_back(J_ip1);

		// Compute the scaling matrix S.
		std::vector<Real> S, invS;

		for(int i=0; i<rangeDim; ++i)
		{
			Real sup;
			range_of_x0[i].mag(sup);

			if(sup == 0)
			{
				S.push_back(0);
				invS.push_back(1);
				symbolic_remainder.scalars[i] = 0;
			}
			else
			{
				S.push_back(sup);
				Real tmp = 1/sup;
				invS.push_back(tmp);
				symbolic_remainder.scalars[i] = tmp;
				range_of_x0[i] = intUnit;
			}
		}

		range_of_x0.insert(range_of_x0.begin(), intZero);

		new_flowpipe.tmv.scale_assign(invS);

		Interval init_cft(-INITIAL_SIMP, INITIAL_SIMP);
		new_flowpipe.tmv.cutoff_normal(tm_setting.step_end_exp_table, init_cft);

		TaylorModelVec<Real> new_x0(S);
		new_x0 += tmv_c0;
		TaylorModelVec<Real> x = new_x0;

		for(int i=0; i<expressions.size(); ++i)
		{
			TaylorModel<Real> tmTemp;
			expressions[i].evaluate(tmTemp, x.tms, tm_setting.order, range_of_x0, tm_setting.cutoff_threshold, g_setting);
			new_flowpipe.tmvPre.tms.push_back(tmTemp);
		}


		if(bPrint)
		{
			printf("Step %d\n", k+1);
		}

		if(bSafetyChecking)
		{
			int safety = new_flowpipe.safetyChecking(unsafeSet, tm_setting, g_setting);

			new_flowpipe.safety = safety;
			flowpipes.push_back(new_flowpipe);

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
			new_flowpipe.safety = SAFE;
			flowpipes.push_back(new_flowpipe);
		}

		current_flowpipe = new_flowpipe;


		if(symbolic_remainder.J.size() >= symbolic_remainder.max_size)
		{
			symbolic_remainder.reset(current_flowpipe.tmvPre.tms.size());
		}
	}

	return checking_result;
}

template <class DATA_TYPE>
int DDE<DATA_TYPE>::reach_inv_symbolic_remainder(TaylorModelFlowpipes & flowpipes, Flowpipe & last_flowpipe, const unsigned int n, const Flowpipe & initialSet, const std::vector<Constraint> & invariant, Symbolic_Remainder & symbolic_remainder, const Taylor_Model_Setting & tm_setting,
		const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking) const
{
	int checking_result = COMPLETED_SAFE;
	unsigned int rangeDim = expressions.size();
	Interval intZero, intUnit(-1,1);

	Flowpipe current_flowpipe = initialSet;

	TaylorModelVec<Real> tmv_flowpipe;
	initialSet.compose(tmv_flowpipe, tm_setting.order, tm_setting.cutoff_threshold);

	TaylorModelFlowpipe init_flowpipe;
	init_flowpipe.tmv_flowpipe = tmv_flowpipe;
	init_flowpipe.domain = initialSet.domain;
	init_flowpipe.bConstrained = initialSet.bConstrained;
	flowpipes.tmv_flowpipes.push_back(init_flowpipe);


	bool bContracted = initialSet.bConstrained;


	for(int k=0; k<n; ++k)
	{
		Flowpipe new_flowpipe;
		new_flowpipe.domain = current_flowpipe.domain;

		TaylorModelVec<Real> tmv_of_x0 = current_flowpipe.tmvPre;

		// the center point of x0's polynomial part
		std::vector<Real> const_of_x0;
		tmv_of_x0.constant(const_of_x0);
/*
		for(unsigned int i=0; i<rangeDim; ++i)
		{
			Real c;
			tmv_of_x0.tms[i].remainder.remove_midpoint(c);
			const_of_x0[i] += c;
		}
*/
		TaylorModelVec<Real> tmv_c0(const_of_x0, rangeDim + 1);

		// introduce a new variable r0 such that x0 = c0 + A*r0, then r0 is origin-centered
		tmv_of_x0.rmConstant();


		// decompose the linear and nonlinear part
		TaylorModelVec<Real> x0_linear, x0_other;
		tmv_of_x0.decompose(x0_linear, x0_other);

		Matrix<Real> Phi_L_i(rangeDim, rangeDim);

		x0_linear.linearCoefficients(Phi_L_i);

		Matrix<Real> local_trans_linear = Phi_L_i;

		Phi_L_i.right_scale_assign(symbolic_remainder.scalars);


		// compute the remainder part under the linear transformation
		Matrix<Interval> J_i(rangeDim, 1);

		for(unsigned int i=1; i<symbolic_remainder.Phi_L.size(); ++i)
		{
			symbolic_remainder.Phi_L[i] = Phi_L_i * symbolic_remainder.Phi_L[i];
		}

		symbolic_remainder.Phi_L.push_back(Phi_L_i);

		for(unsigned int i=1; i<symbolic_remainder.Phi_L.size(); ++i)
		{
			J_i += symbolic_remainder.Phi_L[i] * symbolic_remainder.J[i-1];
		}

		Matrix<Interval> J_ip1(rangeDim, 1);

		std::vector<Interval> range_of_x0;


		// compute the local initial set
		if(symbolic_remainder.J.size() > 0)
		{
			// compute the other part
			std::vector<Interval> tmvPolyRange;
			current_flowpipe.tmv.polyRange(tmvPolyRange, current_flowpipe.domain);
			x0_other.insert_ctrunc_normal(new_flowpipe.tmv, current_flowpipe.tmv, tmvPolyRange, tm_setting.step_end_exp_table, rangeDim + 1, tm_setting.order, tm_setting.cutoff_threshold);

			new_flowpipe.tmv.Remainder(J_ip1);

			std::vector<Polynomial<Real> > poly_tmv;
			current_flowpipe.tmv.Expansion(poly_tmv);
			std::vector<Polynomial<Real> > linear_part = local_trans_linear * poly_tmv;

			for(int i=0; i<rangeDim; ++i)
			{
				new_flowpipe.tmv.tms[i].expansion += linear_part[i];
			}

			// contract J_ip1 and J_i
			if(invariant.size() > 0)
			{
				std::vector<Interval> polyRangeOfx0;
				current_flowpipe.tmv.polyRangeNormal(polyRangeOfx0, tm_setting.step_end_exp_table);

				std::vector<Interval> intVecTmp(rangeDim);
				std::vector<Interval> original_remainders(rangeDim);
				std::vector<Interval> contracted_remainders(rangeDim);

				for(int i=0; i<rangeDim; ++i)
				{
					intVecTmp[i] = polyRangeOfx0[i] + const_of_x0[i];
					contracted_remainders[i] = original_remainders[i] = J_ip1[i][0] + J_i[i][0];
				}

				int res = remainder_contraction_int(intVecTmp, contracted_remainders, invariant);

				if(res < 0)
				{
					return -1;
				}

				for(int i=0; i<rangeDim; ++i)
				{
					double lo_diff = contracted_remainders[i].inf() - original_remainders[i].inf();
					J_ip1[i][0].shrink_lo(lo_diff);

					double up_diff = original_remainders[i].sup() - contracted_remainders[i].sup();
					J_ip1[i][0].shrink_up(up_diff);

					new_flowpipe.tmv.tms[i].remainder = contracted_remainders[i];
					range_of_x0.push_back(polyRangeOfx0[i] + new_flowpipe.tmv.tms[i].remainder);
				}
			}
			else
			{
				for(int i=0; i<rangeDim; ++i)
				{
					new_flowpipe.tmv.tms[i].remainder = J_ip1[i][0] + J_i[i][0];
				}

				new_flowpipe.tmv.intEvalNormal(range_of_x0, tm_setting.step_end_exp_table);
			}
		}
		else
		{
			std::vector<Interval> tmvPolyRange;
			current_flowpipe.tmv.polyRange(tmvPolyRange, current_flowpipe.domain);
			tmv_of_x0.insert_ctrunc_normal(new_flowpipe.tmv, current_flowpipe.tmv, tmvPolyRange, tm_setting.step_end_exp_table, rangeDim + 1, tm_setting.order, tm_setting.cutoff_threshold);

			// contract J_ip1
			if(invariant.size() > 0)
			{
				std::vector<Interval> polyRangeOfx0;
				new_flowpipe.tmv.polyRangeNormal(polyRangeOfx0, tm_setting.step_end_exp_table);

				std::vector<Interval> intVecTmp(rangeDim);
				for(int i=0; i<rangeDim; ++i)
				{
					intVecTmp[i] = polyRangeOfx0[i] + const_of_x0[i];
				}

				std::vector<Interval> contracted_remainders(rangeDim);
				for(int i=0; i<rangeDim; ++i)
				{
					contracted_remainders[i] = new_flowpipe.tmv.tms[i].remainder;
				}

				int res = remainder_contraction_int(intVecTmp, contracted_remainders, invariant);

				if(res < 0)
				{
					return -1;
				}

				for(int i=0; i<rangeDim; ++i)
				{
					new_flowpipe.tmv.tms[i].remainder = contracted_remainders[i];
					range_of_x0.push_back(polyRangeOfx0[i] + new_flowpipe.tmv.tms[i].remainder);
				}
			}
			else
			{
				new_flowpipe.tmv.intEvalNormal(range_of_x0, tm_setting.step_end_exp_table);
			}

			new_flowpipe.tmv.Remainder(J_ip1);
		}

		symbolic_remainder.J.push_back(J_ip1);

		// Compute the scaling matrix S.
		std::vector<Real> S, invS;

		for(int i=0; i<rangeDim; ++i)
		{
			Real sup;
			range_of_x0[i].mag(sup);

			if(sup == 0)
			{
				S.push_back(0);
				invS.push_back(1);
				symbolic_remainder.scalars[i] = 0;
			}
			else
			{
				S.push_back(sup);
				Real tmp = 1/sup;
				invS.push_back(tmp);
				symbolic_remainder.scalars[i] = tmp;
				range_of_x0[i] = intUnit;
			}
		}

		range_of_x0.insert(range_of_x0.begin(), intZero);

		new_flowpipe.tmv.scale_assign(invS);

		Interval init_cft(-INITIAL_SIMP, INITIAL_SIMP);
		new_flowpipe.tmv.cutoff_normal(tm_setting.step_end_exp_table, init_cft);

		TaylorModelVec<Real> new_x0(S);
		new_x0 += tmv_c0;
		TaylorModelVec<Real> x = new_x0;

		for(int i=0; i<expressions.size(); ++i)
		{
			TaylorModel<Real> tmTemp;
			expressions[i].evaluate(tmTemp, x.tms, tm_setting.order, range_of_x0, tm_setting.cutoff_threshold, g_setting);
			new_flowpipe.tmvPre.tms.push_back(tmTemp);
		}

		if(bPrint)
		{
			printf("Step %d\n", k+1);
		}


		// constrain the flowpipe by the invariant
		TaylorModelVec<Real> tmv_flowpipe;
		new_flowpipe.compose_normal(tmv_flowpipe, tm_setting.step_exp_table, tm_setting.order, tm_setting.cutoff_threshold);

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
				int safety = safetyChecking(tmv_flowpipe, contracted_domain, unsafeSet, tm_setting, g_setting);

				if(safety == SAFE)
				{
					flowpipe.safety = SAFE;
					flowpipes.tmv_flowpipes.push_back(flowpipe);
				}
				else if(safety == UNSAFE && !bContracted)
				{
					flowpipe.safety = UNSAFE;
					flowpipes.tmv_flowpipes.push_back(flowpipe);
					return COMPLETED_UNSAFE;
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
				int safety = safetyChecking(tmv_flowpipe, contracted_domain, unsafeSet, tm_setting, g_setting);

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
		}

		if(symbolic_remainder.J.size() >= symbolic_remainder.max_size)
		{
			symbolic_remainder.reset(current_flowpipe.tmvPre.tms.size());
		}
	}

	last_flowpipe = current_flowpipe;

	return checking_result;
}

template <class DATA_TYPE>
void DDE<DATA_TYPE>::reach(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const unsigned int n, const std::vector<Constraint> & unsafeSet) const
{
	bool bSafetyChecking = false;

	if(unsafeSet.size() > 0)
	{
		bSafetyChecking = true;
	}

	result.status = reach(result.flowpipes, n, initialSet, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking);

	if(result.flowpipes.size() > 0)
	{
		result.fp_end_of_time = result.flowpipes.back();
	}
}

template <class DATA_TYPE>
void DDE<DATA_TYPE>::reach_inv(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const unsigned int n, const std::vector<Constraint> & invariant, const std::vector<Constraint> & unsafeSet) const
{
	bool bSafetyChecking = false;

	if(unsafeSet.size() > 0)
	{
		bSafetyChecking = true;
	}

	result.status = reach_inv(result.tmv_flowpipes, result.fp_end_of_time, n, initialSet, invariant, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking);

	if(result.tmv_flowpipes.size() > 0)
	{
		result.tmv_fp_end_of_time = result.tmv_flowpipes.tmv_flowpipes.back();
	}
}

template <class DATA_TYPE>
void DDE<DATA_TYPE>::reach(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const unsigned int n, const std::vector<Constraint> & unsafeSet, Symbolic_Remainder & symbolic_remainder) const
{
	bool bSafetyChecking = false;

	if(unsafeSet.size() > 0)
	{
		bSafetyChecking = true;
	}

	result.status = reach_symbolic_remainder(result.flowpipes, n, initialSet, symbolic_remainder, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking);

	if(result.flowpipes.size() > 0)
	{
		result.fp_end_of_time = result.flowpipes.back();
	}
}

template <class DATA_TYPE>
void DDE<DATA_TYPE>::reach_inv(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const unsigned int n, const std::vector<Constraint> & invariant, const std::vector<Constraint> & unsafeSet, Symbolic_Remainder & symbolic_remainder) const
{
	bool bSafetyChecking = false;

	if(unsafeSet.size() > 0)
	{
		bSafetyChecking = true;
	}

	result.status = reach_inv_symbolic_remainder(result.tmv_flowpipes, result.fp_end_of_time, n, initialSet, invariant, symbolic_remainder, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking);

	if(result.tmv_flowpipes.size() > 0)
	{
		result.tmv_fp_end_of_time = result.tmv_flowpipes.tmv_flowpipes.back();
	}
}
















}



#endif /* DISCRETE_H_ */
