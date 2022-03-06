/*---
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.
  The code is released as is under the GNU General Public License (GPL).
---*/

#include "Discrete.h"


using namespace flowstar;


Linear_Discrete_Dynamics::Linear_Discrete_Dynamics(const Matrix<Real> & matrix_A, const Matrix<UnivariateTaylorModel<Real> > & matrix_B, const Matrix<Real> & matrix_C)
{
	A = matrix_A;
	B = matrix_B;
	C = matrix_C;
}

Linear_Discrete_Dynamics::Linear_Discrete_Dynamics(const Linear_Discrete_Dynamics & ldd)
{
	A = ldd.A;
	B = ldd.B;
	C = ldd.C;
}

Linear_Discrete_Dynamics::~Linear_Discrete_Dynamics()
{
}

Linear_Discrete_Dynamics & Linear_Discrete_Dynamics::operator = (const Linear_Discrete_Dynamics & ldd)
{
	if(this == &ldd)
		return *this;

	A = ldd.A;
	B = ldd.B;
	C = ldd.C;

	return *this;
}

int Linear_Discrete_Dynamics::reach(std::list<LinearFlowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const unsigned int n, const Flowpipe & initialSet, const int rem_size, const Taylor_Model_Setting & tm_setting,
		const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput)
{
	unsigned int rangeDim = A.rows();
	num_of_flowpipes = 0;

	int checking_result = COMPLETED_SAFE;
	std::vector<Interval> polyRangeX0;
	std::vector<Interval> range_of_X0(rangeDim);

	std::vector<Interval> domain = initialSet.domain;

	Matrix<UnivariateTaylorModel<Real> > utm_global_Psi(rangeDim, 1);
	Matrix<UnivariateTaylorModel<Real> > rm_global_Phi(rangeDim);

	Zonotope global_tv_remainder(rangeDim);

	Interval intUnit(-1,1);
	Matrix<Interval> uncertain_range(C.cols(), 1, intUnit);
	Zonotope zono_step;

	if(C.rows() > 0)
	{
		Matrix<Interval> im_temp = A * C * uncertain_range;
		Zonotope zonoTmp(im_temp);
		zono_step = zonoTmp;
	}


	if(bSafetyChecking)
	{
		initialSet.tmvPre.polyRange(polyRangeX0, initialSet.domain);

		for(int k=0; k<rangeDim; ++k)
		{
			range_of_X0[k] = polyRangeX0[k] + initialSet.tmvPre.tms[k].remainder;
		}
	}

	interval_utm_setting.order = 0;

	for(int i=0; i<n; ++i)
	{
		LinearFlowpipe newFlowpipe;

		newFlowpipe.Phi = rm_global_Phi * A;

		if(B.rows() > 0)
		{
			newFlowpipe.Psi = A * utm_global_Psi + B;
			utm_global_Psi += rm_global_Phi * B;
		}

		if(C.rows() > 0)
		{
			if(rem_size >= 0 && global_tv_remainder.numOfGen() > rem_size)
			{
				global_tv_remainder.simplify();
			}

			newFlowpipe.tv_remainder = global_tv_remainder + zono_step;

			global_tv_remainder = A * newFlowpipe.tv_remainder;
		}

		rm_global_Phi = newFlowpipe.Phi;

		++num_of_flowpipes;

		if(bPrint)
		{
			printf("Step %d\n", i+1);
		}

		if(bSafetyChecking)
		{
			if(bTMOutput || bPlot)
			{
				flowpipes.push_back(newFlowpipe);
			}

			int safety;

			safety = newFlowpipe.safetyChecking(unsafeSet, tm_setting, g_setting, initialSet.tmvPre, polyRangeX0, range_of_X0, domain);

			if(bTMOutput || bPlot)
			{
				flowpipes_safety.push_back(safety);
			}

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
			if(bTMOutput || bPlot)
			{
				flowpipes.push_back(newFlowpipe);
				flowpipes_safety.push_back(SAFE);
			}
		}
	}

	return checking_result;
}

int Linear_Discrete_Dynamics::reach(std::list<LinearFlowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const unsigned int n, const LinearFlowpipe & l_initial_set, const Flowpipe & g_initial_set,
		const int rem_size, const Taylor_Model_Setting & tm_setting,
		const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput)
{
	unsigned int rangeDim = A.rows();
	num_of_flowpipes = 0;

	int checking_result = COMPLETED_SAFE;
	std::vector<Interval> polyRangeX0;
	std::vector<Interval> range_of_X0(rangeDim);

	std::vector<Interval> domain = g_initial_set.domain;

	Matrix<UnivariateTaylorModel<Real> > utm_global_Psi = l_initial_set.Psi;
	Matrix<UnivariateTaylorModel<Real> > rm_global_Phi = l_initial_set.Phi;
	Zonotope global_tv_remainder = A * l_initial_set.tv_remainder;

	Interval intUnit(-1,1);
	Matrix<Interval> uncertain_range(C.cols(), 1, intUnit);
	Zonotope zono_step;

	if(C.rows() > 0)
	{
		Matrix<Interval> im_temp = A * C * uncertain_range;
		Zonotope zonoTmp(im_temp);
		zono_step = zonoTmp;
	}



	if(bSafetyChecking)
	{
		g_initial_set.tmvPre.polyRange(polyRangeX0, g_initial_set.domain);

		for(int k=0; k<rangeDim; ++k)
		{
			range_of_X0[k] = polyRangeX0[k] + g_initial_set.tmvPre.tms[k].remainder;
		}
	}

	interval_utm_setting.order = 0;


	for(int i=0; i<n; ++i)
	{
		LinearFlowpipe newFlowpipe;

		newFlowpipe.Phi = A * rm_global_Phi;

		if(B.rows() > 0)
		{
			newFlowpipe.Psi = A * utm_global_Psi + B;
			utm_global_Psi += rm_global_Phi * B;
		}

		if(C.rows() > 0)
		{
			if(rem_size >= 0 && global_tv_remainder.numOfGen() > rem_size)
			{
				global_tv_remainder.simplify();
			}

			newFlowpipe.tv_remainder = global_tv_remainder + zono_step;

			global_tv_remainder = A * newFlowpipe.tv_remainder;
		}

		rm_global_Phi = newFlowpipe.Phi;

		++num_of_flowpipes;

		if(bPrint)
		{
			printf("Step %d\n", i+1);
		}

		if(bSafetyChecking)
		{
			if(bTMOutput || bPlot)
			{
				flowpipes.push_back(newFlowpipe);
			}

			int safety;

			safety = newFlowpipe.safetyChecking(unsafeSet, tm_setting, g_setting, g_initial_set.tmvPre, polyRangeX0, range_of_X0, domain);

			if(bTMOutput || bPlot)
			{
				flowpipes_safety.push_back(safety);
			}

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
			if(bTMOutput || bPlot)
			{
				flowpipes.push_back(newFlowpipe);
				flowpipes_safety.push_back(SAFE);
			}
		}
	}

	return checking_result;
}

int Linear_Discrete_Dynamics::reach_inv(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<std::vector<Interval> > & domains, LinearFlowpipe & last_flowpipe, Flowpipe & contracted_initialSet, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const unsigned int n, const Flowpipe & initialSet, const int rem_size, const std::vector<Constraint> & invariant,
		const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking, const bool bPlot, const bool bTMOutput)
{
	unsigned int rangeDim = A.rows();
	num_of_flowpipes = 0;

	int checking_result = COMPLETED_SAFE;
	std::vector<Interval> polyRangeX0;
	std::vector<Interval> range_of_X0(rangeDim);

	std::vector<Interval> domain = initialSet.domain;

	Matrix<UnivariateTaylorModel<Real> > utm_global_Psi(rangeDim, 1);
	Matrix<UnivariateTaylorModel<Real> > rm_global_Phi(rangeDim);
	Zonotope global_tv_remainder(rangeDim);

	Interval intUnit(-1,1);
	Matrix<Interval> uncertain_range(C.cols(), 1, intUnit);
	Zonotope zono_step;

	if(C.rows() > 0)
	{
		Matrix<Interval> im_temp = A * C * uncertain_range;
		Zonotope zonoTmp(im_temp);
		zono_step = zonoTmp;
	}

	contracted_initialSet = initialSet;

	if(bSafetyChecking)
	{
		initialSet.tmvPre.polyRange(polyRangeX0, initialSet.domain);

		for(int k=0; k<rangeDim; ++k)
		{
			range_of_X0[k] = polyRangeX0[k] + initialSet.tmvPre.tms[k].remainder;
		}
	}

	interval_utm_setting.order = 0;

	bool bContracted = false;
	bool previous_flowpipe_contracted = false;


	Matrix<Interval> im_step_rem(rangeDim, 1);


	// remainder increment in each step
	if(B.rows() > 0)
	{
		for(int i=0; i<rangeDim; ++i)
		{
			im_step_rem[i][0] = B[i][0].remainder;
		}
	}

	if(C.rows() > 0)
	{
		im_step_rem += C * uncertain_range;
	}

	Matrix<Interval> previous_remainder;


	for(int i=0; i<n; ++i)
	{
		LinearFlowpipe new_flowpipe;

		new_flowpipe.Phi = rm_global_Phi * A;

		if(B.rows() > 0)
		{
			new_flowpipe.Psi = A * utm_global_Psi + B;
			utm_global_Psi += rm_global_Phi * B;
		}

		if(C.rows() > 0)
		{
			if(rem_size >= 0 && global_tv_remainder.numOfGen() > rem_size)
			{
				global_tv_remainder.simplify();
			}

			new_flowpipe.tv_remainder = global_tv_remainder + zono_step;

			global_tv_remainder = A * new_flowpipe.tv_remainder;
		}

		rm_global_Phi = new_flowpipe.Phi;




		if(previous_flowpipe_contracted)
		{
			polyRangeX0.clear();
			initialSet.tmvPre.polyRange(polyRangeX0, domain);

			for(int k=0; k<rangeDim; ++k)
			{
				range_of_X0[k] = polyRangeX0[k] + initialSet.tmvPre.tms[k].remainder;
			}

			previous_flowpipe_contracted = false;
		}



		TaylorModelVec<Real> tmv_flowpipe;
		new_flowpipe.evaluate(tmv_flowpipe, initialSet.tmvPre, polyRangeX0, range_of_X0, domain, tm_setting);


		if(previous_remainder.rows() > 0)
		{
			// refine the remainder
			Matrix<Interval> refinement = A * previous_remainder + im_step_rem;

			if(new_flowpipe.remainder_constraints.rows() == 0)
			{
				Matrix<Interval> im_temp(rangeDim, 1);
				new_flowpipe.remainder_constraints = im_temp;
			}

			// contracting the remainder by updating the remainder constraints
			for(int j=0; j<rangeDim; ++j)
			{
				tmv_flowpipe.tms[j].remainder.intersect_assign(refinement[j][0]);
				new_flowpipe.remainder_constraints[j][0] = tmv_flowpipe.tms[j].remainder;
			}
		}


		// contract the remainder firstly
		std::vector<Interval> tmv_flowpipe_polyRange;
		tmv_flowpipe.polyRange(tmv_flowpipe_polyRange, domain);

		std::vector<Interval> contracted_remainders(rangeDim);

		for(int k=0; k<rangeDim; ++k)
		{
			contracted_remainders[k] = tmv_flowpipe.tms[k].remainder;
		}

		// contracting the remainder of the current flowpipe
		// the return value indicates whether the remainder is contracted
		int rem_contraction = remainder_contraction_int(tmv_flowpipe_polyRange, contracted_remainders, invariant);


		// the remainder is contracted
		if(rem_contraction == 1)
		{
			if(new_flowpipe.remainder_constraints.rows() == 0)
			{
				Matrix<Interval> im_temp(rangeDim, 1);
				new_flowpipe.remainder_constraints = im_temp;
			}

			for(int k=0; k<rangeDim; ++k)
			{
				tmv_flowpipe.tms[k].remainder = contracted_remainders[k];

				new_flowpipe.remainder_constraints[k][0] = tmv_flowpipe.tms[k].remainder;
			}
		}



		if(bPrint)
		{
			printf("Step %d\n", i+1);
		}


		std::vector<Interval> contracted_domain = domain;
		int type = domain_contraction_int(tmv_flowpipe, contracted_domain, invariant, tm_setting.order, tm_setting.cutoff_threshold, g_setting);


		switch(type)
		{
		case UNSAT:		// the intersection is empty
			return checking_result;
		case SAT:		// the domain is not contracted
		{
			++num_of_flowpipes;
			flowpipe_orders.push_back(tm_setting.order);
			contraction_of_flowpipes.push_back(bContracted);

			if(bSafetyChecking)
			{
				int safety = safetyChecking(tmv_flowpipe, contracted_domain, unsafeSet, tm_setting, g_setting);

				if(bTMOutput || bPlot)
				{
					tmv_flowpipes.push_back(tmv_flowpipe);
					domains.push_back(domain);
				}

				if(safety == SAFE)
				{
					flowpipes_safety.push_back(SAFE);
				}
				else if(safety == UNSAFE && !bContracted)
				{
					flowpipes_safety.push_back(UNSAFE);
					return COMPLETED_UNSAFE;
				}
				else
				{
					flowpipes_safety.push_back(UNKNOWN);

					if(checking_result == COMPLETED_SAFE)
					{
						checking_result = COMPLETED_UNKNOWN;
					}
				}
			}
			else
			{
				if(bTMOutput || bPlot)
				{
					tmv_flowpipes.push_back(tmv_flowpipe);
					domains.push_back(domain);
					flowpipes_safety.push_back(SAFE);
				}
			}

			break;
		}
		case CONTRACTED: 	// the domain is contracted but the time interval is not
		{
			++num_of_flowpipes;
			flowpipe_orders.push_back(tm_setting.order);

			bContracted = true;
			contraction_of_flowpipes.push_back(true);

			previous_flowpipe_contracted = true;

			domain = contracted_domain;
			contracted_initialSet.domain = contracted_domain;

			if(bSafetyChecking)
			{
				int safety = safetyChecking(tmv_flowpipe, contracted_domain, unsafeSet, tm_setting, g_setting);

				if(bTMOutput || bPlot)
				{
					tmv_flowpipes.push_back(tmv_flowpipe);
					domains.push_back(domain);
				}

				if(safety == SAFE)
				{
					flowpipes_safety.push_back(SAFE);
				}
				else
				{
					flowpipes_safety.push_back(UNKNOWN);

					if(checking_result == COMPLETED_SAFE)
					{
						checking_result = COMPLETED_UNKNOWN;
					}
				}
			}
			else
			{
				if(bTMOutput || bPlot)
				{
					tmv_flowpipes.push_back(tmv_flowpipe);
					domains.push_back(domain);
					flowpipes_safety.push_back(SAFE);
				}
			}

			break;
		}
		}

		if(bContracted)
		{
			if(previous_remainder.rows() == 0)
			{
				Matrix<Interval> im_temp(rangeDim, 1);
				previous_remainder = im_temp;
			}

			for(int k=0; k<rangeDim; ++k)
			{
				previous_remainder[k][0] = tmv_flowpipe.tms[k].remainder;
			}
		}

		last_flowpipe = new_flowpipe;
	}

	return checking_result;
}

int Linear_Discrete_Dynamics::reach_inv(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<std::vector<Interval> > & domains, LinearFlowpipe & last_flowpipe, Flowpipe & contracted_initialSet, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const unsigned int n, const LinearFlowpipe & l_initial_set, const Flowpipe & g_initial_set, const int rem_size, const std::vector<Constraint> & invariant,
		const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking, const bool bPlot, const bool bTMOutput)
{
	unsigned int rangeDim = A.rows();
	num_of_flowpipes = 0;

	int checking_result = COMPLETED_SAFE;
	std::vector<Interval> polyRangeX0;
	std::vector<Interval> range_of_X0(rangeDim);

	std::vector<Interval> domain = g_initial_set.domain;;

	Matrix<UnivariateTaylorModel<Real> > utm_global_Psi = l_initial_set.Psi;
	Matrix<UnivariateTaylorModel<Real> > rm_global_Phi = l_initial_set.Phi;
	Zonotope global_tv_remainder = A * l_initial_set.tv_remainder;

	Interval intUnit(-1,1);
	Matrix<Interval> uncertain_range(C.cols(), 1, intUnit);
	Zonotope zono_step;


	if(C.rows() > 0)
	{
		Matrix<Interval> im_temp = A * C * uncertain_range;
		Zonotope zonoTmp(im_temp);
		zono_step = zonoTmp;
	}

	contracted_initialSet = g_initial_set;

	if(bSafetyChecking)
	{
		g_initial_set.tmvPre.polyRange(polyRangeX0, g_initial_set.domain);

		for(int k=0; k<rangeDim; ++k)
		{
			range_of_X0[k] = polyRangeX0[k] + g_initial_set.tmvPre.tms[k].remainder;
		}
	}

	interval_utm_setting.order = 0;

	bool bContracted = false;
	bool previous_flowpipe_contracted = false;


	Matrix<Interval> im_step_rem(rangeDim, 1);


	// remainder increment in each step
	if(B.rows() > 0)
	{
		for(int i=0; i<rangeDim; ++i)
		{
			im_step_rem[i][0] = B[i][0].remainder;
		}
	}

	if(C.rows() > 0)
	{
		im_step_rem += C * uncertain_range;
	}

	Matrix<Interval> previous_remainder;

	if(l_initial_set.remainder_constraints.rows() > 0)
	{
		previous_remainder = l_initial_set.remainder_constraints;
	}


	for(int i=0; i<n; ++i)
	{
		LinearFlowpipe new_flowpipe;

		new_flowpipe.Phi = rm_global_Phi * A;

		if(B.rows() > 0)
		{
			new_flowpipe.Psi = A * utm_global_Psi + B;
			utm_global_Psi += rm_global_Phi * B;
		}

		if(C.rows() > 0)
		{
			if(rem_size >= 0 && global_tv_remainder.numOfGen() > rem_size)
			{
				global_tv_remainder.simplify();
			}

			new_flowpipe.tv_remainder = global_tv_remainder + zono_step;

			global_tv_remainder = A * new_flowpipe.tv_remainder;
		}

		rm_global_Phi = new_flowpipe.Phi;




		if(previous_flowpipe_contracted)
		{
			polyRangeX0.clear();
			g_initial_set.tmvPre.polyRange(polyRangeX0, domain);

			for(int k=0; k<rangeDim; ++k)
			{
				range_of_X0[k] = polyRangeX0[k] + g_initial_set.tmvPre.tms[k].remainder;
			}

			previous_flowpipe_contracted = false;
		}



		TaylorModelVec<Real> tmv_flowpipe;
		new_flowpipe.evaluate(tmv_flowpipe, g_initial_set.tmvPre, polyRangeX0, range_of_X0, domain, tm_setting);


		if(previous_remainder.rows() > 0)
		{
			// refine the remainder
			Matrix<Interval> refinement = A * previous_remainder + im_step_rem;

			if(new_flowpipe.remainder_constraints.rows() == 0)
			{
				Matrix<Interval> im_temp(rangeDim, 1);
				new_flowpipe.remainder_constraints = im_temp;
			}

			// contracting the remainder by updating the remainder constraints
			for(int j=0; j<rangeDim; ++j)
			{
				tmv_flowpipe.tms[j].remainder.intersect_assign(refinement[j][0]);
				new_flowpipe.remainder_constraints[j][0] = tmv_flowpipe.tms[j].remainder;
			}
		}


		// contract the remainder firstly
		std::vector<Interval> tmv_flowpipe_polyRange;
		tmv_flowpipe.polyRange(tmv_flowpipe_polyRange, domain);

		std::vector<Interval> contracted_remainders(rangeDim);

		for(int k=0; k<rangeDim; ++k)
		{
			contracted_remainders[k] = tmv_flowpipe.tms[k].remainder;
		}

		// contracting the remainder of the current flowpipe
		// the return value indicates whether the remainder is contracted
		int rem_contraction = remainder_contraction_int(tmv_flowpipe_polyRange, contracted_remainders, invariant);


		// the remainder is contracted
		if(rem_contraction == 1)
		{
			if(new_flowpipe.remainder_constraints.rows() == 0)
			{
				Matrix<Interval> im_temp(rangeDim, 1);
				new_flowpipe.remainder_constraints = im_temp;
			}

			for(int k=0; k<rangeDim; ++k)
			{
				tmv_flowpipe.tms[k].remainder = contracted_remainders[k];

				new_flowpipe.remainder_constraints[k][0] = tmv_flowpipe.tms[k].remainder;
			}
		}



		if(bPrint)
		{
			printf("Step %d\n", i+1);
		}


		std::vector<Interval> contracted_domain = domain;
		int type = domain_contraction_int(tmv_flowpipe, contracted_domain, invariant, tm_setting.order, tm_setting.cutoff_threshold, g_setting);


		switch(type)
		{
		case UNSAT:		// the intersection is empty
			return checking_result;
		case SAT:		// the domain is not contracted
		{
			++num_of_flowpipes;
			flowpipe_orders.push_back(tm_setting.order);
			contraction_of_flowpipes.push_back(bContracted);

			if(bSafetyChecking)
			{
				int safety = safetyChecking(tmv_flowpipe, contracted_domain, unsafeSet, tm_setting, g_setting);

				if(bTMOutput || bPlot)
				{
					tmv_flowpipes.push_back(tmv_flowpipe);
					domains.push_back(domain);
				}

				if(safety == SAFE)
				{
					flowpipes_safety.push_back(SAFE);
				}
				else if(safety == UNSAFE && !bContracted)
				{
					flowpipes_safety.push_back(UNSAFE);
					return COMPLETED_UNSAFE;
				}
				else
				{
					flowpipes_safety.push_back(UNKNOWN);

					if(checking_result == COMPLETED_SAFE)
					{
						checking_result = COMPLETED_UNKNOWN;
					}
				}
			}
			else
			{
				if(bTMOutput || bPlot)
				{
					tmv_flowpipes.push_back(tmv_flowpipe);
					domains.push_back(domain);
					flowpipes_safety.push_back(SAFE);
				}
			}

			break;
		}
		case CONTRACTED: 	// the domain is contracted but the time interval is not
		{
			++num_of_flowpipes;
			flowpipe_orders.push_back(tm_setting.order);

			bContracted = true;
			contraction_of_flowpipes.push_back(true);

			previous_flowpipe_contracted = true;

			domain = contracted_domain;
			contracted_initialSet.domain = contracted_domain;

			if(bSafetyChecking)
			{
				int safety = safetyChecking(tmv_flowpipe, contracted_domain, unsafeSet, tm_setting, g_setting);

				if(bTMOutput || bPlot)
				{
					tmv_flowpipes.push_back(tmv_flowpipe);
					domains.push_back(domain);
				}

				if(safety == SAFE)
				{
					flowpipes_safety.push_back(SAFE);
				}
				else
				{
					flowpipes_safety.push_back(UNKNOWN);

					if(checking_result == COMPLETED_SAFE)
					{
						checking_result = COMPLETED_UNKNOWN;
					}
				}
			}
			else
			{
				if(bTMOutput || bPlot)
				{
					tmv_flowpipes.push_back(tmv_flowpipe);
					domains.push_back(domain);
					flowpipes_safety.push_back(SAFE);
				}
			}

			break;
		}
		}

		if(bContracted)
		{
			if(previous_remainder.rows() == 0)
			{
				Matrix<Interval> im_temp(rangeDim, 1);
				previous_remainder = im_temp;
			}

			for(int k=0; k<rangeDim; ++k)
			{
				previous_remainder[k][0] = tmv_flowpipe.tms[k].remainder;
			}
		}

		last_flowpipe = new_flowpipe;
	}

	return checking_result;
}


void Linear_Discrete_Dynamics::reach(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const unsigned int n, const int rem_size, const std::vector<Constraint> & unsafeSet)
{
	bool bSafetyChecking = false;

	if(unsafeSet.size() > 0)
	{
		bSafetyChecking = true;
	}

	result.status = reach(result.linear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes,
			result.num_of_flowpipes, n, initialSet, rem_size, setting.tm_setting,
			setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);

	if(result.linear_flowpipes.size() > 0)
	{
		std::vector<Interval> polyRangeX0;
		initialSet.tmvPre.polyRange(polyRangeX0, initialSet.domain);

		unsigned int rangeDim = initialSet.tmvPre.tms.size();
		std::vector<Interval> range_of_X0(rangeDim);

		for(int k=0; k<rangeDim; ++k)
		{
			range_of_X0[k] = polyRangeX0[k] + initialSet.tmvPre.tms[k].remainder;
		}

		TaylorModelVec<Real> tmvTmp;

		result.linear_flowpipes.back().evaluate(result.tmv_end_of_time, initialSet.tmvPre, polyRangeX0, range_of_X0, initialSet.domain, setting.tm_setting);


		result.lfp_end_of_time = result.linear_flowpipes.back();


		Flowpipe fp(result.tmv_end_of_time, initialSet.domain, setting.tm_setting.cutoff_threshold);
		result.fp_end_of_time = fp;
	}
}

void Linear_Discrete_Dynamics::reach(Result_of_Reachability & result, Computational_Setting & setting, const LinearFlowpipe & l_initial_set, const Flowpipe & g_initial_set, const unsigned int n, const int rem_size, const std::vector<Constraint> & unsafeSet)
{
	bool bSafetyChecking = false;

	if(unsafeSet.size() > 0)
	{
		bSafetyChecking = true;
	}

	result.status = reach(result.linear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes,
			result.num_of_flowpipes, n, l_initial_set, g_initial_set, rem_size, setting.tm_setting,
			setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);

	if(result.linear_flowpipes.size() > 0)
	{
		std::vector<Interval> polyRangeX0;
		g_initial_set.tmvPre.polyRange(polyRangeX0, g_initial_set.domain);

		unsigned int rangeDim = g_initial_set.tmvPre.tms.size();
		std::vector<Interval> range_of_X0(rangeDim);

		for(int k=0; k<rangeDim; ++k)
		{
			range_of_X0[k] = polyRangeX0[k] + g_initial_set.tmvPre.tms[k].remainder;
		}

		TaylorModelVec<Real> tmvTmp;

		result.linear_flowpipes.back().evaluate(result.tmv_end_of_time, g_initial_set.tmvPre, polyRangeX0, range_of_X0, g_initial_set.domain, setting.tm_setting);


		result.lfp_end_of_time = result.linear_flowpipes.back();


		Flowpipe fp(result.tmv_end_of_time, g_initial_set.domain, setting.tm_setting.cutoff_threshold);
		result.fp_end_of_time = fp;
	}
}

void Linear_Discrete_Dynamics::reach_inv(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, Flowpipe & contracted_initialSet, const unsigned int n, const int rem_size, const std::vector<Constraint> & invariant, const std::vector<Constraint> & unsafeSet)
{
	bool bSafetyChecking = false;

	if(unsafeSet.size() > 0)
	{
		bSafetyChecking = true;
	}

	result.status = reach_inv(result.tmv_flowpipes, result.tmv_flowpipes_domains, result.lfp_end_of_time, contracted_initialSet, result.orders_of_flowpipes, result.safety_of_flowpipes,
			result.contraction_of_flowpipes, result.num_of_flowpipes, n, initialSet, rem_size, invariant, setting.tm_setting, setting.g_setting,
			setting.bPrint, unsafeSet, bSafetyChecking, true, true);


	// evaluate the Taylor model overapproximation at the end of the time
	if(result.tmv_flowpipes.size() > 0)
	{
		result.tmv_end_of_time = result.tmv_flowpipes.back();

		Flowpipe fp(result.tmv_end_of_time, result.tmv_flowpipes_domains.back(), setting.tm_setting.cutoff_threshold);
		result.fp_end_of_time = fp;
	}
}

void Linear_Discrete_Dynamics::reach_inv(Result_of_Reachability & result, Computational_Setting & setting, const LinearFlowpipe & l_initial_set, const Flowpipe & g_initial_set, Flowpipe & contracted_initialSet, const unsigned int n, const int rem_size, const std::vector<Constraint> & invariant, const std::vector<Constraint> & unsafeSet)
{
	bool bSafetyChecking = false;

	if(unsafeSet.size() > 0)
	{
		bSafetyChecking = true;
	}

	result.status = reach_inv(result.tmv_flowpipes, result.tmv_flowpipes_domains, result.lfp_end_of_time, contracted_initialSet, result.orders_of_flowpipes, result.safety_of_flowpipes,
			result.contraction_of_flowpipes, result.num_of_flowpipes, n, l_initial_set, g_initial_set, rem_size, invariant, setting.tm_setting, setting.g_setting,
			setting.bPrint, unsafeSet, bSafetyChecking, true, true);


	// evaluate the Taylor model overapproximation at the end of the time
	if(result.tmv_flowpipes.size() > 0)
	{
		result.tmv_end_of_time = result.tmv_flowpipes.back();

		Flowpipe fp(result.tmv_end_of_time, result.tmv_flowpipes_domains.back(), setting.tm_setting.cutoff_threshold);
		result.fp_end_of_time = fp;
	}
}
















Nonlinear_Discrete_Dynamics::Nonlinear_Discrete_Dynamics()
{
}

Nonlinear_Discrete_Dynamics::Nonlinear_Discrete_Dynamics(const std::vector<Expression<Interval> > & dynamics)
{
	expressions = dynamics;
}

Nonlinear_Discrete_Dynamics::Nonlinear_Discrete_Dynamics(const Nonlinear_Discrete_Dynamics & ndd)
{
	expressions = ndd.expressions;
}

Nonlinear_Discrete_Dynamics::~Nonlinear_Discrete_Dynamics()
{
}

Nonlinear_Discrete_Dynamics & Nonlinear_Discrete_Dynamics::operator = (const Nonlinear_Discrete_Dynamics & ndd)
{
	if(this == &ndd)
		return *this;

	expressions = ndd.expressions;

	return *this;
}

int Nonlinear_Discrete_Dynamics::reach(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const unsigned int n, const Flowpipe & initialSet, const Taylor_Model_Setting & tm_setting,
		const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
	int checking_result = COMPLETED_SAFE;
	unsigned int rangeDim = expressions.size();
	Interval intZero, intUnit(-1,1);

	flowpipes.push_back(initialSet);
	flowpipes_safety.push_back(COMPLETED_SAFE);
	flowpipe_orders.push_back(tm_setting.order);

	Flowpipe current_flowpipe = initialSet;


	for(int i=0; i<n; ++i)
	{
		Flowpipe new_flowpipe;
		new_flowpipe.domain = current_flowpipe.domain;

		TaylorModelVec<Real> tmv_of_x0 = current_flowpipe.tmvPre;

		// the center point of x0's polynomial part
		std::vector<Real> const_of_x0;
		tmv_of_x0.constant(const_of_x0);

		for(unsigned int i=0; i<rangeDim; ++i)
		{
			Real c;
			tmv_of_x0.tms[i].remainder.remove_midpoint(c);
			const_of_x0[i] += c;
		}

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

		++num_of_flowpipes;

		if(bPrint)
		{
			printf("Step %d\n", i+1);
		}


		if(bSafetyChecking)
		{
			int safety = new_flowpipe.safetyChecking(unsafeSet, tm_setting, g_setting);

			if(bTMOutput || bPlot)
			{
				flowpipes.push_back(new_flowpipe);
				flowpipes_safety.push_back(safety);
				flowpipe_orders.push_back(tm_setting.order);
			}

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
			if(bTMOutput || bPlot)
			{
				flowpipes.push_back(new_flowpipe);
				flowpipes_safety.push_back(SAFE);
				flowpipe_orders.push_back(tm_setting.order);
			}
		}

		current_flowpipe = new_flowpipe;
	}

	return checking_result;
}

int Nonlinear_Discrete_Dynamics::reach_inv(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<std::vector<Interval> > & domains, Flowpipe & last_flowpipe, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const unsigned int n, const Flowpipe & initialSet, const std::vector<Constraint> & invariant, const Taylor_Model_Setting & tm_setting,
		const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
	int checking_result = COMPLETED_SAFE;
	unsigned int rangeDim = expressions.size();
	Interval intZero, intUnit(-1,1);

	flowpipes_safety.push_back(COMPLETED_SAFE);
	flowpipe_orders.push_back(tm_setting.order);

	Flowpipe current_flowpipe = initialSet;

	TaylorModelVec<Real> tmv_flowpipe;
	initialSet.compose(tmv_flowpipe, tm_setting.order, tm_setting.cutoff_threshold);
	tmv_flowpipes.push_back(tmv_flowpipe);
	domains.push_back(initialSet.domain);

	bool bContracted = false;

	for(int i=0; i<n; ++i)
	{
		Flowpipe new_flowpipe;
		new_flowpipe.domain = current_flowpipe.domain;

		TaylorModelVec<Real> tmv_of_x0 = current_flowpipe.tmvPre;

		// the center point of x0's polynomial part
		std::vector<Real> const_of_x0;
		tmv_of_x0.constant(const_of_x0);

		for(unsigned int i=0; i<rangeDim; ++i)
		{
			Real c;
			tmv_of_x0.tms[i].remainder.remove_midpoint(c);
			const_of_x0[i] += c;
		}

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
			++num_of_flowpipes;
			flowpipe_orders.push_back(tm_setting.order);
			contraction_of_flowpipes.push_back(bContracted);

			if(bSafetyChecking)
			{
				int safety = safetyChecking(tmv_flowpipe, contracted_domain, unsafeSet, tm_setting, g_setting);

				if(bTMOutput || bPlot)
				{
					tmv_flowpipes.push_back(tmv_flowpipe);
					domains.push_back(contracted_domain);
				}

				if(safety == SAFE)
				{
					flowpipes_safety.push_back(SAFE);
				}
				else if(safety == UNSAFE && !bContracted)
				{
					flowpipes_safety.push_back(UNSAFE);
					return COMPLETED_UNSAFE;
				}
				else
				{
					flowpipes_safety.push_back(UNKNOWN);

					if(checking_result == COMPLETED_SAFE)
					{
						checking_result = COMPLETED_UNKNOWN;
					}
				}
			}
			else
			{
				if(bTMOutput || bPlot)
				{
					tmv_flowpipes.push_back(tmv_flowpipe);
					domains.push_back(contracted_domain);
					flowpipes_safety.push_back(SAFE);
				}
			}

			current_flowpipe = new_flowpipe;
			break;
		}
		case CONTRACTED: 	// the domain is contracted but the time interval is not
		{
			++num_of_flowpipes;
			flowpipe_orders.push_back(tm_setting.order);

			bContracted = true;
			contraction_of_flowpipes.push_back(true);

			new_flowpipe.domain = contracted_domain;
			new_flowpipe.normalize(tm_setting.cutoff_threshold);

			if(bSafetyChecking)
			{
				int safety = safetyChecking(tmv_flowpipe, contracted_domain, unsafeSet, tm_setting, g_setting);

				if(bTMOutput || bPlot)
				{
					tmv_flowpipes.push_back(tmv_flowpipe);
					domains.push_back(contracted_domain);
				}

				if(safety == SAFE)
				{
					flowpipes_safety.push_back(SAFE);
				}
				else
				{
					flowpipes_safety.push_back(UNKNOWN);

					if(checking_result == COMPLETED_SAFE)
					{
						checking_result = COMPLETED_UNKNOWN;
					}
				}
			}
			else
			{
				if(bTMOutput || bPlot)
				{
					tmv_flowpipes.push_back(tmv_flowpipe);
					domains.push_back(contracted_domain);
					flowpipes_safety.push_back(SAFE);
				}
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

int Nonlinear_Discrete_Dynamics::reach_symbolic_remainder(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const unsigned int n, const Flowpipe & initialSet, Symbolic_Remainder & symbolic_remainder, const Taylor_Model_Setting & tm_setting,
		const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
	int checking_result = COMPLETED_SAFE;
	unsigned int rangeDim = expressions.size();
	Interval intZero, intUnit(-1,1);


	flowpipes.push_back(initialSet);
	flowpipes_safety.push_back(COMPLETED_SAFE);
	flowpipe_orders.push_back(tm_setting.order);

	Flowpipe current_flowpipe = initialSet;



	for(int k=0; k<n; ++k)
	{
		Flowpipe new_flowpipe;
		new_flowpipe.domain = current_flowpipe.domain;

		TaylorModelVec<Real> tmv_of_x0 = current_flowpipe.tmvPre;

		// the center point of x0's polynomial part
		std::vector<Real> const_of_x0;
		tmv_of_x0.constant(const_of_x0);

		for(unsigned int i=0; i<rangeDim; ++i)
		{
			Real c;
			tmv_of_x0.tms[i].remainder.remove_midpoint(c);
			const_of_x0[i] += c;
		}

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

		for(unsigned int i=0; i<symbolic_remainder.Phi_L.size(); ++i)
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
			// compute the polynomial part for the linear transformation
			std::vector<Polynomial<Real> > initial_linear = symbolic_remainder.Phi_L[0] * symbolic_remainder.polynomial_of_initial_set;

			// compute the other part
			std::vector<Interval> tmvPolyRange;
			current_flowpipe.tmv.polyRange(tmvPolyRange, current_flowpipe.domain);
			x0_other.insert_ctrunc_normal(new_flowpipe.tmv, current_flowpipe.tmv, tmvPolyRange, tm_setting.step_end_exp_table, rangeDim + 1, tm_setting.order, tm_setting.cutoff_threshold);


			new_flowpipe.tmv.Remainder(J_ip1);

			Matrix<Interval> x0_rem(rangeDim, 1);
			tmv_of_x0.Remainder(x0_rem);
			J_ip1 += x0_rem;

			for(int i=0; i<rangeDim; ++i)
			{
				new_flowpipe.tmv.tms[i].expansion += initial_linear[i];
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

		TaylorModelVec<Real> new_x0(S);
		new_x0 += tmv_c0;
		TaylorModelVec<Real> x = new_x0;

		for(int i=0; i<expressions.size(); ++i)
		{
			TaylorModel<Real> tmTemp;
			expressions[i].evaluate(tmTemp, x.tms, tm_setting.order, range_of_x0, tm_setting.cutoff_threshold, g_setting);
			new_flowpipe.tmvPre.tms.push_back(tmTemp);
		}


		++num_of_flowpipes;

		if(bPrint)
		{
			printf("Step %d\n", k+1);
		}

		if(bSafetyChecking)
		{
			int safety = new_flowpipe.safetyChecking(unsafeSet, tm_setting, g_setting);

			if(bTMOutput || bPlot)
			{
				flowpipes.push_back(new_flowpipe);
				flowpipes_safety.push_back(safety);
				flowpipe_orders.push_back(tm_setting.order);
			}

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
			if(bTMOutput || bPlot)
			{
				flowpipes.push_back(new_flowpipe);
				flowpipes_safety.push_back(SAFE);
				flowpipe_orders.push_back(tm_setting.order);
			}
		}

		current_flowpipe = new_flowpipe;


		if(symbolic_remainder.J.size() >= symbolic_remainder.max_size)
		{
			symbolic_remainder.reset(current_flowpipe);
		}
	}

	return checking_result;
}

int Nonlinear_Discrete_Dynamics::reach_inv_symbolic_remainder(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<std::vector<Interval> > & domains, Flowpipe & last_flowpipe, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const unsigned int n, const Flowpipe & initialSet, const std::vector<Constraint> & invariant, Symbolic_Remainder & symbolic_remainder, const Taylor_Model_Setting & tm_setting,
		const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
	int checking_result = COMPLETED_SAFE;
	unsigned int rangeDim = expressions.size();
	Interval intZero, intUnit(-1,1);


	flowpipes_safety.push_back(COMPLETED_SAFE);
	flowpipe_orders.push_back(tm_setting.order);

	Flowpipe current_flowpipe = initialSet;

	TaylorModelVec<Real> tmv_flowpipe;
	initialSet.compose(tmv_flowpipe, tm_setting.order, tm_setting.cutoff_threshold);
	tmv_flowpipes.push_back(tmv_flowpipe);
	domains.push_back(initialSet.domain);

	bool bContracted = false;


	for(int k=0; k<n; ++k)
	{
		Flowpipe new_flowpipe;
		new_flowpipe.domain = current_flowpipe.domain;

		TaylorModelVec<Real> tmv_of_x0 = current_flowpipe.tmvPre;

		// the center point of x0's polynomial part
		std::vector<Real> const_of_x0;
		tmv_of_x0.constant(const_of_x0);

		for(unsigned int i=0; i<rangeDim; ++i)
		{
			Real c;
			tmv_of_x0.tms[i].remainder.remove_midpoint(c);
			const_of_x0[i] += c;
		}

		TaylorModelVec<Real> tmv_c0(const_of_x0, rangeDim + 1);

		// introduce a new variable r0 such that x0 = c0 + A*r0, then r0 is origin-centered
		tmv_of_x0.rmConstant();


		// decompose the linear and nonlinear part
		TaylorModelVec<Real> x0_linear, x0_other;
		tmv_of_x0.decompose(x0_linear, x0_other);

		Matrix<Real> Phi_L_i(rangeDim, rangeDim);

		x0_linear.linearCoefficients(Phi_L_i);


		Phi_L_i.right_scale_assign(symbolic_remainder.scalars);


		// compute the remainder part under the linear transformation
		Matrix<Interval> J_i(rangeDim, 1);

		for(unsigned int i=0; i<symbolic_remainder.Phi_L.size(); ++i)
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
			// compute the polynomial part for the linear transformation
			std::vector<Polynomial<Real> > initial_linear = symbolic_remainder.Phi_L[0] * symbolic_remainder.polynomial_of_initial_set;

			// compute the other part
			std::vector<Interval> tmvPolyRange;
			current_flowpipe.tmv.polyRange(tmvPolyRange, current_flowpipe.domain);
			x0_other.insert_ctrunc_normal(new_flowpipe.tmv, current_flowpipe.tmv, tmvPolyRange, tm_setting.step_end_exp_table, rangeDim + 1, tm_setting.order, tm_setting.cutoff_threshold);


			new_flowpipe.tmv.Remainder(J_ip1);

			Matrix<Interval> x0_rem(rangeDim, 1);
			tmv_of_x0.Remainder(x0_rem);
			J_ip1 += x0_rem;

			for(int i=0; i<rangeDim; ++i)
			{
				new_flowpipe.tmv.tms[i].expansion += initial_linear[i];
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
			++num_of_flowpipes;
			flowpipe_orders.push_back(tm_setting.order);
			contraction_of_flowpipes.push_back(bContracted);

			if(bSafetyChecking)
			{
				int safety = safetyChecking(tmv_flowpipe, contracted_domain, unsafeSet, tm_setting, g_setting);

				if(bTMOutput || bPlot)
				{
					tmv_flowpipes.push_back(tmv_flowpipe);
					domains.push_back(contracted_domain);
				}

				if(safety == SAFE)
				{
					flowpipes_safety.push_back(SAFE);
				}
				else if(safety == UNSAFE && !bContracted)
				{
					flowpipes_safety.push_back(UNSAFE);
					return COMPLETED_UNSAFE;
				}
				else
				{
					flowpipes_safety.push_back(UNKNOWN);

					if(checking_result == COMPLETED_SAFE)
					{
						checking_result = COMPLETED_UNKNOWN;
					}
				}
			}
			else
			{
				if(bTMOutput || bPlot)
				{
					tmv_flowpipes.push_back(tmv_flowpipe);
					domains.push_back(contracted_domain);
					flowpipes_safety.push_back(SAFE);
				}
			}

			current_flowpipe = new_flowpipe;
			break;
		}
		case CONTRACTED: 	// the domain is contracted but the time interval is not
		{
			++num_of_flowpipes;
			flowpipe_orders.push_back(tm_setting.order);

			bContracted = true;
			contraction_of_flowpipes.push_back(true);

			new_flowpipe.domain = contracted_domain;
			new_flowpipe.normalize(tm_setting.cutoff_threshold);

			symbolic_remainder.reset(new_flowpipe);

			if(bSafetyChecking)
			{
				int safety = safetyChecking(tmv_flowpipe, contracted_domain, unsafeSet, tm_setting, g_setting);

				if(bTMOutput || bPlot)
				{
					tmv_flowpipes.push_back(tmv_flowpipe);
					domains.push_back(contracted_domain);
				}

				if(safety == SAFE)
				{
					flowpipes_safety.push_back(SAFE);
				}
				else
				{
					flowpipes_safety.push_back(UNKNOWN);

					if(checking_result == COMPLETED_SAFE)
					{
						checking_result = COMPLETED_UNKNOWN;
					}
				}
			}
			else
			{
				if(bTMOutput || bPlot)
				{
					tmv_flowpipes.push_back(tmv_flowpipe);
					domains.push_back(contracted_domain);
					flowpipes_safety.push_back(SAFE);
				}
			}

			current_flowpipe = new_flowpipe;

			break;
		}
		}

		if(symbolic_remainder.J.size() >= symbolic_remainder.max_size)
		{
			symbolic_remainder.reset(current_flowpipe);
		}
	}

	last_flowpipe = current_flowpipe;

	return checking_result;
}

void Nonlinear_Discrete_Dynamics::reach(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const unsigned int n, const std::vector<Constraint> & unsafeSet) const
{
	bool bSafetyChecking = false;

	if(unsafeSet.size() > 0)
	{
		bSafetyChecking = true;
	}

	result.status = reach(result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes, result.num_of_flowpipes,
			n, initialSet, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);

	if(result.tmv_flowpipes.size() > 0)
	{
		result.fp_end_of_time = result.nonlinear_flowpipes.back();
	}
}

void Nonlinear_Discrete_Dynamics::reach_inv(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const unsigned int n, const std::vector<Constraint> & invariant, const std::vector<Constraint> & unsafeSet) const
{
	bool bSafetyChecking = false;

	if(unsafeSet.size() > 0)
	{
		bSafetyChecking = true;
	}

	result.status = reach_inv(result.tmv_flowpipes, result.tmv_flowpipes_domains, result.fp_end_of_time, result.orders_of_flowpipes, result.safety_of_flowpipes, result.contraction_of_flowpipes, result.num_of_flowpipes,
			n, initialSet, invariant, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);

	if(result.tmv_flowpipes.size() > 0)
	{
		result.tmv_end_of_time = result.tmv_flowpipes.back();
	}
}

void Nonlinear_Discrete_Dynamics::reach_sr(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const unsigned int n, Symbolic_Remainder & symbolic_remainder, const std::vector<Constraint> & unsafeSet) const
{
	bool bSafetyChecking = false;

	if(unsafeSet.size() > 0)
	{
		bSafetyChecking = true;
	}

	result.status = reach_symbolic_remainder(result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes, result.num_of_flowpipes,
			n, initialSet, symbolic_remainder, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);

	if(result.tmv_flowpipes.size() > 0)
	{
		result.fp_end_of_time = result.nonlinear_flowpipes.back();
	}
}

void Nonlinear_Discrete_Dynamics::reach_inv_sr(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const unsigned int n, const std::vector<Constraint> & invariant, Symbolic_Remainder & symbolic_remainder, const std::vector<Constraint> & unsafeSet) const
{
	bool bSafetyChecking = false;

	if(unsafeSet.size() > 0)
	{
		bSafetyChecking = true;
	}

	result.status = reach_inv_symbolic_remainder(result.tmv_flowpipes, result.tmv_flowpipes_domains, result.fp_end_of_time, result.orders_of_flowpipes, result.safety_of_flowpipes, result.contraction_of_flowpipes, result.num_of_flowpipes,
			n, initialSet, invariant, symbolic_remainder, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);

	if(result.tmv_flowpipes.size() > 0)
	{
		result.tmv_end_of_time = result.tmv_flowpipes.back();
	}
}








