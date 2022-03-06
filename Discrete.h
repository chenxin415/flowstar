/*---
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.
  The code is released as is under the GNU General Public License (GPL).
---*/

#ifndef DISCRETE_H_
#define DISCRETE_H_


#include "Continuous.h"



namespace flowstar
{


/* The discrete system should be defined in the form of x_k+1 = A x_k + B,
 * such that A and B are constant matrices which may have interval entries.
 */

class Linear_Discrete_Dynamics : public Dynamics
{
protected:
	Matrix<Real> A;
	Matrix<UnivariateTaylorModel<Real> > B;
	Matrix<Real> C;

public:
	Linear_Discrete_Dynamics();
	Linear_Discrete_Dynamics(const Matrix<Real> & matrix_A, const Matrix<UnivariateTaylorModel<Real> > & matrix_B, const Matrix<Real> & matrix_C);
	Linear_Discrete_Dynamics(const Linear_Discrete_Dynamics & ldd);
	~Linear_Discrete_Dynamics();

	Linear_Discrete_Dynamics & operator = (const Linear_Discrete_Dynamics & ldd);

	int reach(std::list<LinearFlowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const unsigned int n, const Flowpipe & initialSet, const int rem_size, const Taylor_Model_Setting & tm_setting,
			const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput);

	// starting from an initial set defined by a linear flowpipe
	int reach(std::list<LinearFlowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const unsigned int n, const LinearFlowpipe & l_initial_set, const Flowpipe & g_initial_set,
			const int rem_size, const Taylor_Model_Setting & tm_setting,
			const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput);

	int reach_inv(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<std::vector<Interval> > & domains, LinearFlowpipe & last_flowpipe, Flowpipe & contracted_initialSet, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const unsigned int n, const Flowpipe & initialSet, const int rem_size, const std::vector<Constraint> & invariant,
			const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking, const bool bPlot, const bool bTMOutput);

	int reach_inv(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<std::vector<Interval> > & domains, LinearFlowpipe & last_flowpipe, Flowpipe & contracted_initialSet, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const unsigned int n, const LinearFlowpipe & l_initial_set, const Flowpipe & g_initial_set, const int rem_size, const std::vector<Constraint> & invariant,
			const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking, const bool bPlot, const bool bTMOutput);


	void reach(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const unsigned int n, const int rem_size, const std::vector<Constraint> & unsafeSet);
	void reach(Result_of_Reachability & result, Computational_Setting & setting, const LinearFlowpipe & l_initial_set, const Flowpipe & g_initial_set, const unsigned int n, const int rem_size, const std::vector<Constraint> & unsafeSet);

	void reach_inv(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, Flowpipe & contracted_initialSet, const unsigned int n, const int rem_size, const std::vector<Constraint> & invariant, const std::vector<Constraint> & unsafeSet);
	void reach_inv(Result_of_Reachability & result, Computational_Setting & setting, const LinearFlowpipe & l_initial_set, const Flowpipe & g_initial_set, Flowpipe & contracted_initialSet, const unsigned int n, const int rem_size, const std::vector<Constraint> & invariant, const std::vector<Constraint> & unsafeSet);
};










/* Class for the dynamics defined by nonlinear difference equations
 * which may have interval coefficients representing time-varying uncertainties.
 */

class Nonlinear_Discrete_Dynamics : public Dynamics
{
protected:
	std::vector<Expression<Interval> >	expressions;

public:
	Nonlinear_Discrete_Dynamics();
	Nonlinear_Discrete_Dynamics(const std::vector<Expression<Interval> > & dynamics);
	Nonlinear_Discrete_Dynamics(const Nonlinear_Discrete_Dynamics & ndd);
	~Nonlinear_Discrete_Dynamics();

	Nonlinear_Discrete_Dynamics & operator = (const Nonlinear_Discrete_Dynamics & ndd);


	int reach(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const unsigned int n, const Flowpipe & initialSet, const Taylor_Model_Setting & tm_setting,
			const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;

	int reach_inv(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<std::vector<Interval> > & domains, Flowpipe & last_flowpipe, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const unsigned int n, const Flowpipe & initialSet, const std::vector<Constraint> & invariant, const Taylor_Model_Setting & tm_setting,
			const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;


	int reach_symbolic_remainder(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			unsigned long & num_of_flowpipes, const unsigned int n, const Flowpipe & initialSet, Symbolic_Remainder & symbolic_remainder, const Taylor_Model_Setting & tm_setting,
			const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;

	int reach_inv_symbolic_remainder(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<std::vector<Interval> > & domains, Flowpipe & last_flowpipe, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
			std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const unsigned int n, const Flowpipe & initialSet, const std::vector<Constraint> & invariant, Symbolic_Remainder & symbolic_remainder, const Taylor_Model_Setting & tm_setting,
			const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
			const bool bPlot, const bool bTMOutput) const;



	void reach(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const unsigned int n, const std::vector<Constraint> & unsafeSet) const;
	void reach_inv(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const unsigned int n, const std::vector<Constraint> & invariant, const std::vector<Constraint> & unsafeSet) const;

	void reach_sr(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const unsigned int n, Symbolic_Remainder & symbolic_remainder, const std::vector<Constraint> & unsafeSet) const;
	void reach_inv_sr(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const unsigned int n, const std::vector<Constraint> & invariant, Symbolic_Remainder & symbolic_remainder, const std::vector<Constraint> & unsafeSet) const;
};





}



#endif /* DISCRETE_H_ */
