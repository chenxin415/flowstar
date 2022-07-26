/*---
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.
  The code is released as is under the GNU General Public License (GPL).
---*/

#include "Discrete.h"


using namespace flowstar;


LTI_DDE::LTI_DDE(const Matrix<Real> & A)
{
	rm_dyn_A = A;
}

LTI_DDE::LTI_DDE(const Matrix<Real> & A, const Matrix<Real> & B)
{
	rm_dyn_A = A;
	rm_dyn_B = B;
}

LTI_DDE::LTI_DDE(const Matrix<Real> & A, const Matrix<Real> & B, const Matrix<Real> & C)
{
	rm_dyn_A = A;
	rm_dyn_B = B;
	rm_dyn_C = C;
}

LTI_DDE::LTI_DDE(const Matrix<Real> & A, const Matrix<Real> & B, const Matrix<Real> & C, const Matrix<Real> & D)
{
	rm_dyn_A = A;
	rm_dyn_B = B;
	rm_dyn_C = C;
	rm_dyn_D = D;
}

LTI_DDE::LTI_DDE(const LTI_DDE & lti_dde)
{
	rm_dyn_A = lti_dde.rm_dyn_A;
	rm_dyn_B = lti_dde.rm_dyn_B;
	rm_dyn_C = lti_dde.rm_dyn_C;
	rm_dyn_D = lti_dde.rm_dyn_D;
}

LTI_DDE::~LTI_DDE()
{
}

LTI_DDE & LTI_DDE::operator = (const LTI_DDE & lti_dde)
{
	if(this == &lti_dde)
		return *this;

	rm_dyn_A = lti_dde.rm_dyn_A;
	rm_dyn_B = lti_dde.rm_dyn_B;
	rm_dyn_C = lti_dde.rm_dyn_C;
	rm_dyn_D = lti_dde.rm_dyn_D;

	return *this;
}

void LTI_DDE::abstract(FlowmapAbstraction & abstraction, const int N, const int zono_order)
{
	unsigned int rangeDim = rm_dyn_A.rows();

	Matrix<Real> rm_global_Phi(rangeDim);
	Matrix<Real> rm_global_Psi(rangeDim, 1);
	Matrix<Real> rm_global_Omega(rangeDim, rm_dyn_C.cols());

	Zonotope global_tv_remainder(rangeDim);

	Interval intUnit(-1,1);
	Matrix<Interval> uncertain_range(rm_dyn_D.cols(), 1, intUnit);
	Zonotope zono_step;

	if(rm_dyn_D.cols() > 0)
	{
		Matrix<Interval> im_temp = rm_dyn_D * uncertain_range;
		Zonotope zonoTmp(im_temp);
		zono_step = zonoTmp;
	}


	for(int i=0; i<N; ++i)
	{
		LinearFlowmap flowmap;

		Matrix<Real> rm_temp = rm_dyn_A * rm_global_Phi;
		flowmap.Phi = rm_temp;

		if(rm_dyn_B.cols() > 0)
		{
			rm_global_Psi = rm_global_Psi + rm_global_Phi * rm_dyn_B;
			flowmap.Psi = rm_global_Psi;
		}

		if(rm_dyn_C.cols() > 0)
		{
			rm_global_Omega = rm_global_Omega + rm_global_Omega * rm_dyn_C;
			flowmap.Omega = rm_global_Omega;
		}

		if(rm_dyn_D.cols() > 0)
		{
			if(zono_order >= 0 && global_tv_remainder.numOfGen() > zono_order)
			{
				global_tv_remainder.simplify();
			}


			global_tv_remainder = global_tv_remainder + rm_global_Phi * zono_step;
			flowmap.tv_remainder = global_tv_remainder;
		}

		rm_global_Phi = rm_temp;

		abstraction.flowmaps.push_back(flowmap);
	}
}



/*
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
*/



















