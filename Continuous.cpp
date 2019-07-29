/*---
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.
  The code is released as is under the GNU General Public License (GPL).
---*/

#include "Continuous.h"

using namespace flowstar;

Variables stateVars;
Variables tmVars;
Variables tvPars;

Continuous_Reachability_Problem_Description problem_description;
Continuous_Reachability reachability_for_outputFile;



Symbolic_Remainder::Symbolic_Remainder()
{
}

Symbolic_Remainder::Symbolic_Remainder(const Flowpipe & initial_set)
{
	scalars.resize(initial_set.tmvPre.tms.size(), 1);
	initial_set.tmv.Expansion(polynomial_of_initial_set);
}

Symbolic_Remainder::Symbolic_Remainder(const Symbolic_Remainder & symbolic_remainder)
{
	J							= symbolic_remainder.J;
	Phi_L						= symbolic_remainder.Phi_L;
	scalars						= symbolic_remainder.scalars;
	polynomial_of_initial_set	= symbolic_remainder.polynomial_of_initial_set;
}

Symbolic_Remainder::~Symbolic_Remainder()
{
}

void Symbolic_Remainder::reset(const Flowpipe & initial_set)
{
	scalars.clear();
	scalars.resize(initial_set.tmvPre.tms.size(), 1);

	initial_set.tmv.Expansion(polynomial_of_initial_set);

	J.clear();
	Phi_L.clear();
}

Symbolic_Remainder & Symbolic_Remainder::operator = (const Symbolic_Remainder & symbolic_remainder)
{
	if(this == &symbolic_remainder)
		return *this;

	J							= symbolic_remainder.J;
	Phi_L						= symbolic_remainder.Phi_L;
	scalars						= symbolic_remainder.scalars;
	polynomial_of_initial_set	= symbolic_remainder.polynomial_of_initial_set;

	return *this;
}



Computational_Setting::Computational_Setting()
{
	time = 0;
	bPrint = false;
}

Computational_Setting::Computational_Setting(const Computational_Setting & setting)
{
	tm_setting			= setting.tm_setting;
	g_setting			= setting.g_setting;
	symbolic_remainder	= setting.symbolic_remainder;
	time				= setting.time;
	bPrint				= setting.bPrint;
}

Computational_Setting::~Computational_Setting()
{
}

void Computational_Setting::clear()
{
	tm_setting.clear();
}

bool Computational_Setting::setTime(const double t)
{
	if(t < 0)
	{
		return false;
	}

	time = t;

	return true;
}

bool Computational_Setting::setFixedStepsize(const double step, const unsigned int order)
{
	tm_setting.step_min = -1;
	return tm_setting.setStepsize(step, order);
}

bool Computational_Setting::setFixedStepsize(const double step, const unsigned int order_min, const unsigned int order_max)
{
	if(step <= 0 || order_min <= 1 || order_max <= 1 || order_min > order_max)
	{
		return false;
	}

	tm_setting.setStepsize(step, order_max);

	tm_setting.order = order_min;
	tm_setting.order_min = order_min;
	tm_setting.order_max = order_max;

	return true;
}

bool Computational_Setting::setAdaptiveStepsize(const double step_min, const double step_max, const unsigned int order)
{
	if(step_min <= 0 || step_max < 0 || step_min > step_max || order <= 1)
	{
		return false;
	}

	tm_setting.step_min = step_min;
	tm_setting.step_max = step_max;
	tm_setting.order = order;

	tm_setting.setStepsize(step_max, order);

	return true;
}

bool Computational_Setting::setCutoffThreshold(const double threshold)
{
	if(threshold < 0)
	{
		return false;
	}

	Interval I(-threshold, threshold);
	tm_setting.cutoff_threshold = I;

	return true;
}

void Computational_Setting::setRemainderEstimation(const std::vector<Interval> & estimation)
{
	tm_setting.remainder_estimation = estimation;
}

void Computational_Setting::setQueueSize(const unsigned int m)
{
	tm_setting.queue_size = m;
}

bool Computational_Setting::resetOrder(const unsigned int order)
{
	bool bValid = tm_setting.resetOrder(order);

	if(!bValid)
	{
		return false;
	}

	bValid = g_setting.resetOrder(order);
	return bValid;
}

bool Computational_Setting::resetOrder(const unsigned int order_min, const unsigned int order_max)
{
	if(order_min < 2 || order_max < 2 || order_min > order_max)
	{
		return false;
	}

	tm_setting.resetOrder(order_max);
	g_setting.resetOrder(order_max);

	return true;
}

void Computational_Setting::printOn()
{
	bPrint = true;
}

void Computational_Setting::printOff()
{
	bPrint = false;
}

void Computational_Setting::prepare()
{
	unsigned int maxOrder = tm_setting.order_min > tm_setting.order_max ? tm_setting.order_min : tm_setting.order_max;

	maxOrder = tm_setting.order > maxOrder ? tm_setting.order : maxOrder;

	g_setting.prepareForReachability(maxOrder);
}

Computational_Setting & Computational_Setting::operator = (const Computational_Setting & setting)
{
	if(this == &setting)
		return *this;

	tm_setting			= setting.tm_setting;
	g_setting			= setting.g_setting;
	symbolic_remainder	= setting.symbolic_remainder;
	time				= setting.time;
	bPrint				= setting.bPrint;

	return *this;
}






Flowpipe::Flowpipe()
{
}

Flowpipe::Flowpipe(const std::vector<Interval> & box)
{
	TaylorModelVec<Real> tmv1(box, domain);
	TaylorModelVec<Real> tmv2(box.size());

	tmvPre = tmv1;
	tmv = tmv2;
}

Flowpipe::Flowpipe(const TaylorModelVec<Real> & tmv_flowpipe, const std::vector<Interval> & flowpipe_domain, const Interval & cutoff_threshold)
{
	tmvPre = tmv_flowpipe;
	domain = flowpipe_domain;
	tmvPre.normalize(domain, cutoff_threshold);

	TaylorModelVec<Real> tmvTmp(domain.size() - 1);
	tmv = tmvTmp;
}

Flowpipe::Flowpipe(const Flowpipe & flowpipe)
{
	tmvPre = flowpipe.tmvPre;
	tmv = flowpipe.tmv;
	domain = flowpipe.domain;
}

Flowpipe::~Flowpipe()
{
}

void Flowpipe::clear()
{
	tmvPre.clear();
	tmv.clear();
	domain.clear();
}

void Flowpipe::compose(TaylorModelVec<Real> & result, const unsigned int order, const Interval & cutoff_threshold) const
{
	std::vector<Interval> tmvPolyRange;
	tmv.polyRange(tmvPolyRange, domain);
	tmvPre.insert_ctrunc(result, tmv, tmvPolyRange, domain, order, cutoff_threshold);
}

void Flowpipe::compose(TaylorModelVec<Real> & result, const std::vector<unsigned int> & orders, const Interval & cutoff_threshold) const
{
	std::vector<Interval> tmvPolyRange;
	tmv.polyRange(tmvPolyRange, domain);
	tmvPre.insert_ctrunc(result, tmv, tmvPolyRange, domain, orders, cutoff_threshold);
}

void Flowpipe::compose(TaylorModelVec<Real> & result, const std::vector<unsigned int> & outputAxes, const unsigned int order, const Interval & cutoff_threshold) const
{
	std::vector<Interval> tmvPolyRange;
	tmv.polyRange(tmvPolyRange, domain);

	result.clear();

	for(unsigned int i=0; i<outputAxes.size(); ++i)
	{
		TaylorModel<Real> tmTmp;
		tmvPre.tms[outputAxes[i]].insert_ctrunc(tmTmp, tmv, tmvPolyRange, domain, order, cutoff_threshold);
		result.tms.push_back(tmTmp);
	}
}

void Flowpipe::compose_normal(TaylorModelVec<Real> & result, const std::vector<Interval> & step_exp_table, const unsigned int order, const Interval & cutoff_threshold) const
{
	std::vector<Interval> tmvPolyRange;
	tmv.polyRangeNormal(tmvPolyRange, step_exp_table);
	tmvPre.insert_ctrunc_normal(result, tmv, tmvPolyRange, step_exp_table, domain.size(), order, cutoff_threshold);
}

void Flowpipe::compose_normal(TaylorModelVec<Real> & result, const std::vector<Interval> & step_exp_table, const std::vector<unsigned int> & orders, const Interval & cutoff_threshold) const
{
	std::vector<Interval> tmvPolyRange;
	tmv.polyRangeNormal(tmvPolyRange, step_exp_table);
	tmvPre.insert_ctrunc_normal(result, tmv, tmvPolyRange, step_exp_table, domain.size(), orders, cutoff_threshold);
}
void Flowpipe::compose_normal(TaylorModelVec<Real> & result, const std::vector<unsigned int> & outputAxes, const unsigned int order, const std::vector<Interval> & step_exp_table, const Interval & cutoff_threshold) const
{
	std::vector<Interval> tmvPolyRange;
	tmv.polyRange(tmvPolyRange, domain);

	result.clear();

	for(unsigned int i=0; i<outputAxes.size(); ++i)
	{
		TaylorModel<Real> tmTmp;
		tmvPre.tms[outputAxes[i]].insert_ctrunc_normal(tmTmp, tmv, tmvPolyRange, step_exp_table, domain.size(), order, cutoff_threshold);
		result.tms.push_back(tmTmp);
	}
}

void Flowpipe::intEval(std::vector<Interval> & result, const unsigned int order, const Interval & cutoff_threshold) const
{
	TaylorModelVec<Real> tmvTmp;
	compose(tmvTmp, order, cutoff_threshold);
	tmvTmp.intEval(result, domain);
}

void Flowpipe::intEvalNormal(std::vector<Interval> & result, const std::vector<Interval> & step_exp_table, const unsigned int order, const Interval & cutoff_threshold) const
{
	TaylorModelVec<Real> tmvTmp;
	compose_normal(tmvTmp, step_exp_table, order, cutoff_threshold);
	tmvTmp.intEvalNormal(result, step_exp_table);
}

void Flowpipe::normalize(const Interval & cutoff_threshold)
{
	// tmv is always normalized, we only need to normalize tmvPre
	tmvPre.normalize(domain, cutoff_threshold);
}

int Flowpipe::safetyChecking(const std::vector<Constraint> & unsafeSet, const Taylor_Model_Computation_Setting & tm_setting, const Global_Computation_Setting & g_setting) const
{
	if(unsafeSet.size() == 0)
	{
		return SAFE;
	}

	unsigned int rangeDim = tmvPre.tms.size();
	int result = UNKNOWN;
	bool bContained = true;

	std::vector<Interval> tmvRange;
	tmvPre.intEvalNormal(tmvRange, tm_setting.step_exp_table);

	for(unsigned int i=0; i<unsafeSet.size(); ++i)
	{
		Interval I;

		// interval evaluation on the constraint
		unsafeSet[i].expression.evaluate(I, tmvRange);

		if(unsafeSet[i].bound < I.inf())
		{
			// no intersection with the unsafe set
			result = SAFE;
			break;
		}
		else
		{
			if(!(unsafeSet[i].bound >= I.sup()) && bContained)
			{
				bContained = false;
			}
		}
	}

	if(result == UNKNOWN)
	{
		if(bContained)
		{
			return UNSAFE;
		}
		else
		{
			// do a simple branch & bound for safety checking
			TaylorModelVec<Real> tmvFlowpipe;
			compose(tmvFlowpipe, tm_setting.order, tm_setting.cutoff_threshold);

			std::vector<HornerForm<Real> > obj_hfs;
			std::vector<Interval> obj_rems;

			result = SAFE;

			for(unsigned int i=0; i<unsafeSet.size(); ++i)
			{
				TaylorModel<Real> tmTmp;

				// interval evaluation on the constraint
				unsafeSet[i].expression.evaluate(tmTmp, tmvFlowpipe.tms, tm_setting.order, domain, tm_setting.cutoff_threshold, g_setting);

				HornerForm<Real> obj_hf;
				tmTmp.expansion.toHornerForm(obj_hf);
				obj_hfs.push_back(obj_hf);
				obj_rems.push_back(tmTmp.remainder);
			}

			std::vector<Interval> refined_domain = domain;

			std::list<Interval> subdivisions;

			if(domain[0].width() > REFINEMENT_PREC)
			{
				subdivisions.push_back(domain[0]);
			}

			for(; subdivisions.size() > 0; )
			{
				Interval subdivision = subdivisions.front();
				subdivisions.pop_front();

				int result_iter = UNKNOWN;
				bool bContained_iter = true;

				refined_domain[0] = subdivision;

				for(int i=0; i<unsafeSet.size(); ++i)
				{
					Interval I;
					obj_hfs[i].evaluate(I, refined_domain);

					I += obj_rems[i];

					if(unsafeSet[i].bound < I.inf())
					{
						// no intersection with the unsafe set
						result_iter = SAFE;
						break;
					}
					else
					{
						if(!(unsafeSet[i].bound >= I.sup()) && bContained_iter)
						{
							bContained_iter = false;
						}
					}
				}

				if(result_iter == UNKNOWN)
				{
					if(bContained_iter)
					{
						return UNSAFE;
					}
					else
					{
						if(subdivision.width() <= REFINEMENT_PREC)
						{
							return UNKNOWN;
						}

						// split the domain
						Interval I1, I2;
						subdivision.split(I1, I2);

						if(I1.width() <= REFINEMENT_PREC)
						{
							if(result == SAFE)
								result = UNKNOWN;
						}
						else
						{
							subdivisions.push_back(I1);
						}

						if(I2.width() <= REFINEMENT_PREC)
						{
							if(result == SAFE)
								result = UNKNOWN;
						}
						else
						{
							subdivisions.push_back(I2);
						}
					}
				}
			}

			return result;
		}
	}
	else
	{
		return SAFE;
	}
}

Flowpipe & Flowpipe::operator = (const Flowpipe & flowpipe)
{
	if(this == &flowpipe)
		return *this;

	tmvPre = flowpipe.tmvPre;
	tmv = flowpipe.tmv;
	domain = flowpipe.domain;

	return *this;
}

int Flowpipe::advance_deterministic(Flowpipe & result, const std::vector<Expression_AST<Real> > & ode, const Taylor_Model_Computation_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Computation_Setting & g_setting) const
{
	unsigned int rangeDim = ode.size();
	unsigned int rangeDimExt = rangeDim + 1;
	Interval intUnit(-1,1);
	result.clear();

	// evaluate the the initial set x0
	TaylorModelVec<Real> tmv_of_x0;
	tmvPre.evaluate_time(tmv_of_x0, tm_setting.step_end_exp_table);

	// the center point of x0's polynomial part
	std::vector<Real> const_of_x0;
	tmv_of_x0.constant(const_of_x0);

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		Real c;
		tmv_of_x0.tms[i].remainder.remove_midpoint(c);
		const_of_x0[i] += c;
	}

	TaylorModelVec<Real> tmv_c0(const_of_x0, rangeDimExt);

	// introduce a new variable r0 such that x0 = c0 + A*r0, then r0 is origin-centered
	tmv_of_x0.rmConstant();

	std::vector<Interval> tmvPolyRange;
	tmv.polyRangeNormal(tmvPolyRange, tm_setting.step_end_exp_table);
	tmv_of_x0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, tm_setting.step_end_exp_table, domain.size(), tm_setting.order, tm_setting.cutoff_threshold);

	std::vector<Interval> range_of_x0;

	// contract the remainder part of the initial set
	if(invariant.size() > 0)
	{
		std::vector<Interval> polyRangeOfx0;
		result.tmv.polyRangeNormal(polyRangeOfx0, tm_setting.step_end_exp_table);

		std::vector<Interval> intVecTmp(rangeDim);
		for(int i=0; i<rangeDim; ++i)
		{
			intVecTmp[i] = polyRangeOfx0[i] + const_of_x0[i];
		}

		std::vector<Interval> contracted_remainders(rangeDim);
		for(int i=0; i<rangeDim; ++i)
		{
			contracted_remainders[i] = result.tmv.tms[i].remainder;
		}

		int res = contract_remainder(intVecTmp, contracted_remainders, invariant);

		if(res < 0)
		{
			return -1;
		}

		for(int i=0; i<rangeDim; ++i)
		{
			result.tmv.tms[i].remainder = contracted_remainders[i];
			range_of_x0.push_back(polyRangeOfx0[i] + result.tmv.tms[i].remainder);
		}
	}
	else
	{
		result.tmv.intEvalNormal(range_of_x0, tm_setting.step_end_exp_table);
	}


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
			invS.push_back(1/sup);
			range_of_x0[i] = intUnit;
		}
	}

	result.tmv.scale_assign(invS);
//	result.tmv.cutoff_normal(tm_setting.step_end_exp_table, tm_setting.cutoff_threshold);

	TaylorModelVec<Real> new_x0(S);
	new_x0 += tmv_c0;
	TaylorModelVec<Real> x = new_x0;

	for(unsigned int i=1; i<=tm_setting.order; ++i)
	{
		x.Picard_no_remainder_assign(new_x0, ode, rangeDimExt, i, tm_setting.cutoff_threshold);
	}

//	x.cutoff(tm_setting.cutoff_threshold);

	bool bfound = true;

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tm_setting.remainder_estimation[i];
	}

	TaylorModelVec<Real> tmvTmp;
	std::list<Interval> intermediate_ranges;

	x.Picard_ctrunc_normal(tmvTmp, new_x0, ode, tm_setting.step_exp_table, rangeDimExt, tm_setting.order, tm_setting.cutoff_threshold, intermediate_ranges, g_setting);

	// compute the interval evaluation of the polynomial difference due to the roundoff error
	std::vector<Interval> intDifferences;
	for(unsigned int i=0; i<rangeDim; ++i)
	{
		Polynomial<Real> polyTmp;
		polyTmp = tmvTmp.tms[i].expansion - x.tms[i].expansion;

		Interval I;
		polyTmp.intEvalNormal(I, tm_setting.step_exp_table);

		intDifferences.push_back(I);
	}

	// add the uncertainties and the cutoff intervals onto the result
	for(unsigned int i=0; i<rangeDim; ++i)
	{
		tmvTmp.tms[i].remainder += intDifferences[i];

		if( ! tmvTmp.tms[i].remainder.subseteq(x.tms[i].remainder) )
		{
			bfound = false;
			break;
		}
	}

	if(!bfound)
	{
		return 0;
	}
	else
	{
		for(unsigned int i=0; i<rangeDim; ++i)
		{
			x.tms[i].remainder = tmvTmp.tms[i].remainder;
		}
	}

	bool bfinished = false;
	for(int rSteps = 0; !bfinished && (rSteps <= MAX_REFINEMENT_STEPS); ++rSteps)
	{
		bfinished = true;

		std::vector<Interval> newRemainders;
		x.Picard_ctrunc_normal_remainder(newRemainders, ode, tm_setting.step_exp_table[1], tm_setting.order, intermediate_ranges, g_setting);

		// add the uncertainties and the cutoff intervals onto the result
		for(unsigned int i=0; i<rangeDim; ++i)
		{
			newRemainders[i] += intDifferences[i];

			if(newRemainders[i].subseteq(x.tms[i].remainder))
			{
				if(x.tms[i].remainder.widthRatio(newRemainders[i]) <= STOP_RATIO)
				{
					bfinished = false;
				}

				x.tms[i].remainder = newRemainders[i];
			}
			else
			{
				bfinished = true;
				break;
			}
		}
	}

	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = tm_setting.step_exp_table[1];

	return 1;
}

int Flowpipe::advance_nondeterministic(Flowpipe & result, const std::vector<Expression_AST<Interval> > & ode, const Taylor_Model_Computation_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Computation_Setting & g_setting) const
{
	unsigned int rangeDim = ode.size();
	unsigned int rangeDimExt = rangeDim + 1;
	Interval intUnit(-1,1);
	result.clear();

	// evaluate the the initial set x0
	TaylorModelVec<Real> tmv_of_x0;
	tmvPre.evaluate_time(tmv_of_x0, tm_setting.step_end_exp_table);

	// the center point of x0's polynomial part
	std::vector<Real> const_of_x0;
	tmv_of_x0.constant(const_of_x0);

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		Real c;
		tmv_of_x0.tms[i].remainder.remove_midpoint(c);
		const_of_x0[i] += c;
	}

	TaylorModelVec<Real> tmv_c0(const_of_x0, rangeDimExt);

	// introduce a new variable r0 such that x0 = c0 + A*r0, then r0 is origin-centered
	tmv_of_x0.rmConstant();

	std::vector<Interval> tmvPolyRange;
	tmv.polyRangeNormal(tmvPolyRange, tm_setting.step_end_exp_table);
	tmv_of_x0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, tm_setting.step_end_exp_table, domain.size(), tm_setting.order, tm_setting.cutoff_threshold);

	std::vector<Interval> range_of_x0;

	// contract the remainder part of the initial set
	if(invariant.size() > 0)
	{
		std::vector<Interval> polyRangeOfx0;
		result.tmv.polyRangeNormal(polyRangeOfx0, tm_setting.step_end_exp_table);

		std::vector<Interval> intVecTmp(rangeDim);
		for(int i=0; i<rangeDim; ++i)
		{
			intVecTmp[i] = polyRangeOfx0[i] + const_of_x0[i];
		}

		std::vector<Interval> contracted_remainders(rangeDim);
		for(int i=0; i<rangeDim; ++i)
		{
			contracted_remainders[i] = result.tmv.tms[i].remainder;
		}

		int res = contract_remainder(intVecTmp, contracted_remainders, invariant);

		if(res < 0)
		{
			return -1;
		}

		for(int i=0; i<rangeDim; ++i)
		{
			result.tmv.tms[i].remainder = contracted_remainders[i];
			range_of_x0.push_back(polyRangeOfx0[i] + result.tmv.tms[i].remainder);
		}
	}
	else
	{
		result.tmv.intEvalNormal(range_of_x0, tm_setting.step_end_exp_table);
	}


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
			invS.push_back(1/sup);
			range_of_x0[i] = intUnit;
		}
	}

	result.tmv.scale_assign(invS);
//	result.tmv.cutoff_normal(tm_setting.step_end_exp_table, tm_setting.cutoff_threshold);

	TaylorModelVec<Real> new_x0(S);
	new_x0 += tmv_c0;
	TaylorModelVec<Real> x = new_x0;

	for(unsigned int i=1; i<=tm_setting.order; ++i)
	{
		x.Picard_no_remainder_assign(new_x0, ode, rangeDimExt, i, tm_setting.cutoff_threshold);
	}

//	x.cutoff(tm_setting.cutoff_threshold);

	bool bfound = true;

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tm_setting.remainder_estimation[i];
	}

	TaylorModelVec<Real> tmvTmp;
	std::list<Interval> intermediate_ranges;

	x.Picard_ctrunc_normal(tmvTmp, new_x0, ode, tm_setting.step_exp_table, rangeDimExt, tm_setting.order, tm_setting.cutoff_threshold, intermediate_ranges, g_setting);

	// compute the interval evaluation of the polynomial difference due to the roundoff error
	std::vector<Interval> intDifferences;
	for(unsigned int i=0; i<rangeDim; ++i)
	{
		Polynomial<Real> polyTmp;
		polyTmp = tmvTmp.tms[i].expansion - x.tms[i].expansion;

		Interval I;
		polyTmp.intEvalNormal(I, tm_setting.step_exp_table);

		intDifferences.push_back(I);
	}

	// add the uncertainties and the cutoff intervals onto the result
	for(unsigned int i=0; i<rangeDim; ++i)
	{
		tmvTmp.tms[i].remainder += intDifferences[i];

		if( ! tmvTmp.tms[i].remainder.subseteq(x.tms[i].remainder) )
		{
			bfound = false;
			break;
		}
	}

	if(!bfound)
	{
		return 0;
	}
	else
	{
		for(unsigned int i=0; i<rangeDim; ++i)
		{
			x.tms[i].remainder = tmvTmp.tms[i].remainder;
		}
	}

	bool bfinished = false;
	for(int rSteps = 0; !bfinished && (rSteps <= MAX_REFINEMENT_STEPS); ++rSteps)
	{
		bfinished = true;

		std::vector<Interval> newRemainders;
		x.Picard_ctrunc_normal_remainder(newRemainders, ode, tm_setting.step_exp_table[1], tm_setting.order, intermediate_ranges, g_setting);

		// add the uncertainties and the cutoff intervals onto the result
		for(unsigned int i=0; i<rangeDim; ++i)
		{
			newRemainders[i] += intDifferences[i];

			if(newRemainders[i].subseteq(x.tms[i].remainder))
			{
				if(x.tms[i].remainder.widthRatio(newRemainders[i]) <= STOP_RATIO)
				{
					bfinished = false;
				}

				x.tms[i].remainder = newRemainders[i];
			}
			else
			{
				bfinished = true;
				break;
			}
		}
	}

	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = tm_setting.step_exp_table[1];

	return 1;
}

int Flowpipe::advance_deterministic_adaptive_stepsize(Flowpipe & result, const std::vector<Expression_AST<Real> > & ode, Taylor_Model_Computation_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Computation_Setting & g_setting) const
{
	unsigned int rangeDim = ode.size();
	unsigned int rangeDimExt = rangeDim + 1;
	Interval intUnit(-1,1);
	result.clear();

	// evaluate the the initial set x0
	TaylorModelVec<Real> tmv_of_x0;
	tmvPre.evaluate_time(tmv_of_x0, tm_setting.step_end_exp_table);

	// the center point of x0's polynomial part
	std::vector<Real> const_of_x0;
	tmv_of_x0.constant(const_of_x0);

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		Real c;
		tmv_of_x0.tms[i].remainder.remove_midpoint(c);
		const_of_x0[i] += c;
	}

	TaylorModelVec<Real> tmv_c0(const_of_x0, rangeDimExt);

	// introduce a new variable r0 such that x0 = c0 + A*r0, then r0 is origin-centered
	tmv_of_x0.rmConstant();

	std::vector<Interval> tmvPolyRange;
	tmv.polyRangeNormal(tmvPolyRange, tm_setting.step_end_exp_table);
	tmv_of_x0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, tm_setting.step_end_exp_table, domain.size(), tm_setting.order, tm_setting.cutoff_threshold);

	std::vector<Interval> range_of_x0;

	// contract the remainder part of the initial set
	if(invariant.size() > 0)
	{
		std::vector<Interval> polyRangeOfx0;
		result.tmv.polyRangeNormal(polyRangeOfx0, tm_setting.step_end_exp_table);

		std::vector<Interval> intVecTmp(rangeDim);
		for(int i=0; i<rangeDim; ++i)
		{
			intVecTmp[i] = polyRangeOfx0[i] + const_of_x0[i];
		}

		std::vector<Interval> contracted_remainders(rangeDim);
		for(int i=0; i<rangeDim; ++i)
		{
			contracted_remainders[i] = result.tmv.tms[i].remainder;
		}

		int res = contract_remainder(intVecTmp, contracted_remainders, invariant);

		if(res < 0)
		{
			return -1;
		}

		for(int i=0; i<rangeDim; ++i)
		{
			result.tmv.tms[i].remainder = contracted_remainders[i];
			range_of_x0.push_back(polyRangeOfx0[i] + result.tmv.tms[i].remainder);
		}
	}
	else
	{
		result.tmv.intEvalNormal(range_of_x0, tm_setting.step_end_exp_table);
	}


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
			invS.push_back(1/sup);
			range_of_x0[i] = intUnit;
		}
	}

	result.tmv.scale_assign(invS);
//	result.tmv.cutoff_normal(tm_setting.step_end_exp_table, tm_setting.cutoff_threshold);

	TaylorModelVec<Real> new_x0(S);
	new_x0 += tmv_c0;
	TaylorModelVec<Real> x = new_x0;

	for(unsigned int i=1; i<=tm_setting.order; ++i)
	{
		x.Picard_no_remainder_assign(new_x0, ode, rangeDimExt, i, tm_setting.cutoff_threshold);
	}

//	x.cutoff(tm_setting.cutoff_threshold);

	bool bfound = true;

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tm_setting.remainder_estimation[i];
	}

	TaylorModelVec<Real> tmvTmp;
	std::list<Interval> intermediate_ranges;

	x.Picard_ctrunc_normal(tmvTmp, new_x0, ode, tm_setting.step_exp_table, rangeDimExt, tm_setting.order, tm_setting.cutoff_threshold, intermediate_ranges, g_setting);

	// compute the interval evaluation of the polynomial difference due to the roundoff error
	std::vector<Interval> intDifferences;
	std::vector<Polynomial<Real> > polyDifferences;
	for(unsigned int i=0; i<rangeDim; ++i)
	{
		Polynomial<Real> polyTmp;
		polyTmp = tmvTmp.tms[i].expansion - x.tms[i].expansion;
		polyDifferences.push_back(polyTmp);

		Interval I;
		polyTmp.intEvalNormal(I, tm_setting.step_exp_table);

		intDifferences.push_back(I);
	}

	// add the uncertainties and the cutoff intervals onto the result
	for(unsigned int i=0; i<rangeDim; ++i)
	{
		tmvTmp.tms[i].remainder += intDifferences[i];

		if( ! tmvTmp.tms[i].remainder.subseteq(x.tms[i].remainder) )
		{
			bfound = false;
			break;
		}
	}

	for(; !bfound;)
	{
		bfound = true;
		double newStep = tm_setting.step_exp_table[1].sup() * LAMBDA_DOWN;	// reduce the time step size

		if(newStep < tm_setting.step_min)
		{
			return 0;
		}

		tm_setting.setStepsize(newStep, tm_setting.order);

		intermediate_ranges.clear();
		x.Picard_ctrunc_normal(tmvTmp, new_x0, ode, tm_setting.step_exp_table, rangeDimExt, tm_setting.order, tm_setting.cutoff_threshold, intermediate_ranges, g_setting);

		// recompute the interval evaluation of the polynomial differences
		for(int i=0; i<rangeDim; ++i)
		{
			polyDifferences[i].intEvalNormal(intDifferences[i], tm_setting.step_exp_table);

			tmvTmp.tms[i].remainder += intDifferences[i];

			if( ! tmvTmp.tms[i].remainder.subseteq(x.tms[i].remainder) )
			{
				bfound = false;
				break;
			}
		}
	}

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tmvTmp.tms[i].remainder;
	}

	bool bfinished = false;
	for(int rSteps = 0; !bfinished && (rSteps <= MAX_REFINEMENT_STEPS); ++rSteps)
	{
		bfinished = true;

		std::vector<Interval> newRemainders;
		x.Picard_ctrunc_normal_remainder(newRemainders, ode, tm_setting.step_exp_table[1], tm_setting.order, intermediate_ranges, g_setting);

		// add the uncertainties and the cutoff intervals onto the result
		for(unsigned int i=0; i<rangeDim; ++i)
		{
			newRemainders[i] += intDifferences[i];

			if(newRemainders[i].subseteq(x.tms[i].remainder))
			{
				if(x.tms[i].remainder.widthRatio(newRemainders[i]) <= STOP_RATIO)
				{
					bfinished = false;
				}

				x.tms[i].remainder = newRemainders[i];
			}
			else
			{
				bfinished = true;
				break;
			}
		}
	}

	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = tm_setting.step_exp_table[1];

	return 1;
}

int Flowpipe::advance_nondeterministic_adaptive_stepsize(Flowpipe & result, const std::vector<Expression_AST<Interval> > & ode, Taylor_Model_Computation_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Computation_Setting & g_setting) const
{
	unsigned int rangeDim = ode.size();
	unsigned int rangeDimExt = rangeDim + 1;
	Interval intUnit(-1,1);
	result.clear();

	// evaluate the the initial set x0
	TaylorModelVec<Real> tmv_of_x0;
	tmvPre.evaluate_time(tmv_of_x0, tm_setting.step_end_exp_table);

	// the center point of x0's polynomial part
	std::vector<Real> const_of_x0;
	tmv_of_x0.constant(const_of_x0);

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		Real c;
		tmv_of_x0.tms[i].remainder.remove_midpoint(c);
		const_of_x0[i] += c;
	}

	TaylorModelVec<Real> tmv_c0(const_of_x0, rangeDimExt);

	// introduce a new variable r0 such that x0 = c0 + A*r0, then r0 is origin-centered
	tmv_of_x0.rmConstant();

	std::vector<Interval> tmvPolyRange;
	tmv.polyRangeNormal(tmvPolyRange, tm_setting.step_end_exp_table);
	tmv_of_x0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, tm_setting.step_end_exp_table, domain.size(), tm_setting.order, tm_setting.cutoff_threshold);

	std::vector<Interval> range_of_x0;

	// contract the remainder part of the initial set
	if(invariant.size() > 0)
	{
		std::vector<Interval> polyRangeOfx0;
		result.tmv.polyRangeNormal(polyRangeOfx0, tm_setting.step_end_exp_table);

		std::vector<Interval> intVecTmp(rangeDim);
		for(int i=0; i<rangeDim; ++i)
		{
			intVecTmp[i] = polyRangeOfx0[i] + const_of_x0[i];
		}

		std::vector<Interval> contracted_remainders(rangeDim);
		for(int i=0; i<rangeDim; ++i)
		{
			contracted_remainders[i] = result.tmv.tms[i].remainder;
		}

		int res = contract_remainder(intVecTmp, contracted_remainders, invariant);

		if(res < 0)
		{
			return -1;
		}

		for(int i=0; i<rangeDim; ++i)
		{
			result.tmv.tms[i].remainder = contracted_remainders[i];
			range_of_x0.push_back(polyRangeOfx0[i] + result.tmv.tms[i].remainder);
		}
	}
	else
	{
		result.tmv.intEvalNormal(range_of_x0, tm_setting.step_end_exp_table);
	}


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
			invS.push_back(1/sup);
			range_of_x0[i] = intUnit;
		}
	}

	result.tmv.scale_assign(invS);

	TaylorModelVec<Real> new_x0(S);
	new_x0 += tmv_c0;
	TaylorModelVec<Real> x = new_x0;

	for(unsigned int i=1; i<=tm_setting.order; ++i)
	{
		x.Picard_no_remainder_assign(new_x0, ode, rangeDimExt, i, tm_setting.cutoff_threshold);
	}

	bool bfound = true;

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tm_setting.remainder_estimation[i];
	}

	TaylorModelVec<Real> tmvTmp;
	std::list<Interval> intermediate_ranges;

	x.Picard_ctrunc_normal(tmvTmp, new_x0, ode, tm_setting.step_exp_table, rangeDimExt, tm_setting.order, tm_setting.cutoff_threshold, intermediate_ranges, g_setting);

	// compute the interval evaluation of the polynomial difference due to the roundoff error
	std::vector<Interval> intDifferences;
	std::vector<Polynomial<Real> > polyDifferences;
	for(unsigned int i=0; i<rangeDim; ++i)
	{
		Polynomial<Real> polyTmp;
		polyTmp = tmvTmp.tms[i].expansion - x.tms[i].expansion;
		polyDifferences.push_back(polyTmp);

		Interval I;
		polyTmp.intEvalNormal(I, tm_setting.step_exp_table);

		intDifferences.push_back(I);
	}

	// add the uncertainties and the cutoff intervals onto the result
	for(unsigned int i=0; i<rangeDim; ++i)
	{
		tmvTmp.tms[i].remainder += intDifferences[i];

		if( ! tmvTmp.tms[i].remainder.subseteq(x.tms[i].remainder) )
		{
			bfound = false;
			break;
		}
	}

	for(; !bfound;)
	{
		bfound = true;
		double newStep = tm_setting.step_exp_table[1].sup() * LAMBDA_DOWN;	// reduce the time step size

		if(newStep < tm_setting.step_min)
		{
			return 0;
		}

		tm_setting.setStepsize(newStep, tm_setting.order);

		intermediate_ranges.clear();
		x.Picard_ctrunc_normal(tmvTmp, new_x0, ode, tm_setting.step_exp_table, rangeDimExt, tm_setting.order, tm_setting.cutoff_threshold, intermediate_ranges, g_setting);

		// recompute the interval evaluation of the polynomial differences
		for(int i=0; i<rangeDim; ++i)
		{
			polyDifferences[i].intEvalNormal(intDifferences[i], tm_setting.step_exp_table);

			tmvTmp.tms[i].remainder += intDifferences[i];

			if( ! tmvTmp.tms[i].remainder.subseteq(x.tms[i].remainder) )
			{
				bfound = false;
				break;
			}
		}
	}

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tmvTmp.tms[i].remainder;
	}

	bool bfinished = false;
	for(int rSteps = 0; !bfinished && (rSteps <= MAX_REFINEMENT_STEPS); ++rSteps)
	{
		bfinished = true;

		std::vector<Interval> newRemainders;
		x.Picard_ctrunc_normal_remainder(newRemainders, ode, tm_setting.step_exp_table[1], tm_setting.order, intermediate_ranges, g_setting);

		// add the uncertainties and the cutoff intervals onto the result
		for(unsigned int i=0; i<rangeDim; ++i)
		{
			newRemainders[i] += intDifferences[i];

			if(newRemainders[i].subseteq(x.tms[i].remainder))
			{
				if(x.tms[i].remainder.widthRatio(newRemainders[i]) <= STOP_RATIO)
				{
					bfinished = false;
				}

				x.tms[i].remainder = newRemainders[i];
			}
			else
			{
				bfinished = true;
				break;
			}
		}
	}

	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = tm_setting.step_exp_table[1];

	return 1;
}

int Flowpipe::advance_deterministic_adaptive_order(Flowpipe & result, const std::vector<Expression_AST<Real> > & ode, Taylor_Model_Computation_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Computation_Setting & g_setting) const
{
	unsigned int rangeDim = ode.size();
	unsigned int rangeDimExt = rangeDim + 1;
	Interval intUnit(-1,1);
	result.clear();

	// evaluate the the initial set x0
	TaylorModelVec<Real> tmv_of_x0;
	tmvPre.evaluate_time(tmv_of_x0, tm_setting.step_end_exp_table);

	// the center point of x0's polynomial part
	std::vector<Real> const_of_x0;
	tmv_of_x0.constant(const_of_x0);

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		Real c;
		tmv_of_x0.tms[i].remainder.remove_midpoint(c);
		const_of_x0[i] += c;
	}

	TaylorModelVec<Real> tmv_c0(const_of_x0, rangeDimExt);

	// introduce a new variable r0 such that x0 = c0 + A*r0, then r0 is origin-centered
	tmv_of_x0.rmConstant();

	std::vector<Interval> tmvPolyRange;
	tmv.polyRangeNormal(tmvPolyRange, tm_setting.step_end_exp_table);
	tmv_of_x0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, tm_setting.step_end_exp_table, domain.size(), tm_setting.order, tm_setting.cutoff_threshold);

	std::vector<Interval> range_of_x0;

	// contract the remainder part of the initial set
	if(invariant.size() > 0)
	{
		std::vector<Interval> polyRangeOfx0;
		result.tmv.polyRangeNormal(polyRangeOfx0, tm_setting.step_end_exp_table);

		std::vector<Interval> intVecTmp(rangeDim);
		for(int i=0; i<rangeDim; ++i)
		{
			intVecTmp[i] = polyRangeOfx0[i] + const_of_x0[i];
		}

		std::vector<Interval> contracted_remainders(rangeDim);
		for(int i=0; i<rangeDim; ++i)
		{
			contracted_remainders[i] = result.tmv.tms[i].remainder;
		}

		int res = contract_remainder(intVecTmp, contracted_remainders, invariant);

		if(res < 0)
		{
			return -1;
		}

		for(int i=0; i<rangeDim; ++i)
		{
			result.tmv.tms[i].remainder = contracted_remainders[i];
			range_of_x0.push_back(polyRangeOfx0[i] + result.tmv.tms[i].remainder);
		}
	}
	else
	{
		result.tmv.intEvalNormal(range_of_x0, tm_setting.step_end_exp_table);
	}


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
			invS.push_back(1/sup);
			range_of_x0[i] = intUnit;
		}
	}

	result.tmv.scale_assign(invS);

	TaylorModelVec<Real> new_x0(S);
	new_x0 += tmv_c0;
	TaylorModelVec<Real> x = new_x0;

	for(unsigned int i=1; i<=tm_setting.order; ++i)
	{
		x.Picard_no_remainder_assign(new_x0, ode, rangeDimExt, i, tm_setting.cutoff_threshold);
	}

	bool bfound = true;

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tm_setting.remainder_estimation[i];
	}

	TaylorModelVec<Real> tmvTmp;
	std::list<Interval> intermediate_ranges;

	x.Picard_ctrunc_normal(tmvTmp, new_x0, ode, tm_setting.step_exp_table, rangeDimExt, tm_setting.order, tm_setting.cutoff_threshold, intermediate_ranges, g_setting);

	// compute the interval evaluation of the polynomial difference due to the roundoff error
	std::vector<Interval> intDifferences;

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		Polynomial<Real> polyTmp;
		polyTmp = tmvTmp.tms[i].expansion - x.tms[i].expansion;

		Interval I;
		polyTmp.intEvalNormal(I, tm_setting.step_exp_table);

		intDifferences.push_back(I);
	}

	// add the uncertainties and the cutoff intervals onto the result
	for(unsigned int i=0; i<rangeDim; ++i)
	{
		tmvTmp.tms[i].remainder += intDifferences[i];

		if( ! tmvTmp.tms[i].remainder.subseteq(x.tms[i].remainder) )
		{
			bfound = false;
			break;
		}
	}

	for(; !bfound;)
	{
		bfound = true;

		++tm_setting.order;

		if(tm_setting.order > tm_setting.order_max)
		{
			return 0;
		}

		x.Picard_no_remainder_assign(new_x0, ode, rangeDimExt, tm_setting.order, tm_setting.cutoff_threshold);

		for(int i=0; i<rangeDim; ++i)
		{
			x.tms[i].remainder = tm_setting.remainder_estimation[i];
		}

		intermediate_ranges.clear();
		x.Picard_ctrunc_normal(tmvTmp, new_x0, ode, tm_setting.step_exp_table, rangeDimExt, tm_setting.order, tm_setting.cutoff_threshold, intermediate_ranges, g_setting);

		// recompute the interval evaluation of the polynomial differences
		for(int i=0; i<rangeDim; ++i)
		{
			Polynomial<Real> polyTmp;
			polyTmp = tmvTmp.tms[i].expansion - x.tms[i].expansion;

			Interval I;
			polyTmp.intEvalNormal(I, tm_setting.step_exp_table);

			intDifferences[i] = I;

			tmvTmp.tms[i].remainder += intDifferences[i];

			if( ! tmvTmp.tms[i].remainder.subseteq(x.tms[i].remainder) )
			{
				bfound = false;
				break;
			}
		}
	}

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tmvTmp.tms[i].remainder;
	}

	bool bfinished = false;
	for(int rSteps = 0; !bfinished && (rSteps <= MAX_REFINEMENT_STEPS); ++rSteps)
	{
		bfinished = true;

		std::vector<Interval> newRemainders;
		x.Picard_ctrunc_normal_remainder(newRemainders, ode, tm_setting.step_exp_table[1], tm_setting.order, intermediate_ranges, g_setting);

		// add the uncertainties and the cutoff intervals onto the result
		for(unsigned int i=0; i<rangeDim; ++i)
		{
			newRemainders[i] += intDifferences[i];

			if(newRemainders[i].subseteq(x.tms[i].remainder))
			{
				if(x.tms[i].remainder.widthRatio(newRemainders[i]) <= STOP_RATIO)
				{
					bfinished = false;
				}

				x.tms[i].remainder = newRemainders[i];
			}
			else
			{
				bfinished = true;
				break;
			}
		}
	}

	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = tm_setting.step_exp_table[1];

	return 1;
}

int Flowpipe::advance_nondeterministic_adaptive_order(Flowpipe & result, const std::vector<Expression_AST<Interval> > & ode, Taylor_Model_Computation_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Computation_Setting & g_setting) const
{
	unsigned int rangeDim = ode.size();
	unsigned int rangeDimExt = rangeDim + 1;
	Interval intUnit(-1,1);
	result.clear();

	// evaluate the the initial set x0
	TaylorModelVec<Real> tmv_of_x0;
	tmvPre.evaluate_time(tmv_of_x0, tm_setting.step_end_exp_table);

	// the center point of x0's polynomial part
	std::vector<Real> const_of_x0;
	tmv_of_x0.constant(const_of_x0);

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		Real c;
		tmv_of_x0.tms[i].remainder.remove_midpoint(c);
		const_of_x0[i] += c;
	}

	TaylorModelVec<Real> tmv_c0(const_of_x0, rangeDimExt);

	// introduce a new variable r0 such that x0 = c0 + A*r0, then r0 is origin-centered
	tmv_of_x0.rmConstant();

	std::vector<Interval> tmvPolyRange;
	tmv.polyRangeNormal(tmvPolyRange, tm_setting.step_end_exp_table);
	tmv_of_x0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, tm_setting.step_end_exp_table, domain.size(), tm_setting.order, tm_setting.cutoff_threshold);

	std::vector<Interval> range_of_x0;

	// contract the remainder part of the initial set
	if(invariant.size() > 0)
	{
		std::vector<Interval> polyRangeOfx0;
		result.tmv.polyRangeNormal(polyRangeOfx0, tm_setting.step_end_exp_table);

		std::vector<Interval> intVecTmp(rangeDim);
		for(int i=0; i<rangeDim; ++i)
		{
			intVecTmp[i] = polyRangeOfx0[i] + const_of_x0[i];
		}

		std::vector<Interval> contracted_remainders(rangeDim);
		for(int i=0; i<rangeDim; ++i)
		{
			contracted_remainders[i] = result.tmv.tms[i].remainder;
		}

		int res = contract_remainder(intVecTmp, contracted_remainders, invariant);

		if(res < 0)
		{
			return -1;
		}

		for(int i=0; i<rangeDim; ++i)
		{
			result.tmv.tms[i].remainder = contracted_remainders[i];
			range_of_x0.push_back(polyRangeOfx0[i] + result.tmv.tms[i].remainder);
		}
	}
	else
	{
		result.tmv.intEvalNormal(range_of_x0, tm_setting.step_end_exp_table);
	}


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
			invS.push_back(1/sup);
			range_of_x0[i] = intUnit;
		}
	}

	result.tmv.scale_assign(invS);

	TaylorModelVec<Real> new_x0(S);
	new_x0 += tmv_c0;
	TaylorModelVec<Real> x = new_x0;

	for(unsigned int i=1; i<=tm_setting.order; ++i)
	{
		x.Picard_no_remainder_assign(new_x0, ode, rangeDimExt, i, tm_setting.cutoff_threshold);
	}

	bool bfound = true;

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tm_setting.remainder_estimation[i];
	}

	TaylorModelVec<Real> tmvTmp;
	std::list<Interval> intermediate_ranges;

	x.Picard_ctrunc_normal(tmvTmp, new_x0, ode, tm_setting.step_exp_table, rangeDimExt, tm_setting.order, tm_setting.cutoff_threshold, intermediate_ranges, g_setting);

	// compute the interval evaluation of the polynomial difference due to the roundoff error
	std::vector<Interval> intDifferences;

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		Polynomial<Real> polyTmp;
		polyTmp = tmvTmp.tms[i].expansion - x.tms[i].expansion;

		Interval I;
		polyTmp.intEvalNormal(I, tm_setting.step_exp_table);

		intDifferences.push_back(I);
	}

	// add the uncertainties and the cutoff intervals onto the result
	for(unsigned int i=0; i<rangeDim; ++i)
	{
		tmvTmp.tms[i].remainder += intDifferences[i];

		if( ! tmvTmp.tms[i].remainder.subseteq(x.tms[i].remainder) )
		{
			bfound = false;
			break;
		}
	}

	for(; !bfound;)
	{
		bfound = true;

		++tm_setting.order;

		if(tm_setting.order > tm_setting.order_max)
		{
			return 0;
		}

		x.Picard_no_remainder_assign(new_x0, ode, rangeDimExt, tm_setting.order, tm_setting.cutoff_threshold);

		for(int i=0; i<rangeDim; ++i)
		{
			x.tms[i].remainder = tm_setting.remainder_estimation[i];
		}

		intermediate_ranges.clear();
		x.Picard_ctrunc_normal(tmvTmp, new_x0, ode, tm_setting.step_exp_table, rangeDimExt, tm_setting.order, tm_setting.cutoff_threshold, intermediate_ranges, g_setting);

		// recompute the interval evaluation of the polynomial differences
		for(int i=0; i<rangeDim; ++i)
		{
			Polynomial<Real> polyTmp;
			polyTmp = tmvTmp.tms[i].expansion - x.tms[i].expansion;

			Interval I;
			polyTmp.intEvalNormal(I, tm_setting.step_exp_table);

			intDifferences[i] = I;

			tmvTmp.tms[i].remainder += intDifferences[i];

			if( ! tmvTmp.tms[i].remainder.subseteq(x.tms[i].remainder) )
			{
				bfound = false;
				break;
			}
		}
	}

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tmvTmp.tms[i].remainder;
	}

	bool bfinished = false;
	for(int rSteps = 0; !bfinished && (rSteps <= MAX_REFINEMENT_STEPS); ++rSteps)
	{
		bfinished = true;

		std::vector<Interval> newRemainders;
		x.Picard_ctrunc_normal_remainder(newRemainders, ode, tm_setting.step_exp_table[1], tm_setting.order, intermediate_ranges, g_setting);

		// add the uncertainties and the cutoff intervals onto the result
		for(unsigned int i=0; i<rangeDim; ++i)
		{
			newRemainders[i] += intDifferences[i];

			if(newRemainders[i].subseteq(x.tms[i].remainder))
			{
				if(x.tms[i].remainder.widthRatio(newRemainders[i]) <= STOP_RATIO)
				{
					bfinished = false;
				}

				x.tms[i].remainder = newRemainders[i];
			}
			else
			{
				bfinished = true;
				break;
			}
		}
	}

	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = tm_setting.step_exp_table[1];

	return 1;
}






int Flowpipe::advance_deterministic(Flowpipe & result, const std::vector<Expression_AST<Real> > & ode, const Taylor_Model_Computation_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Computation_Setting & g_setting, Symbolic_Remainder & symbolic_remainder) const
{
	unsigned int rangeDim = ode.size();
	unsigned int rangeDimExt = rangeDim + 1;
	Interval intUnit(-1,1);
	result.clear();

	// evaluate the the initial set x0
	TaylorModelVec<Real> tmv_of_x0;
	tmvPre.evaluate_time(tmv_of_x0, tm_setting.step_end_exp_table);

	// the center point of x0's polynomial part
	std::vector<Real> const_of_x0;
	tmv_of_x0.constant(const_of_x0);

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		Real c;
		tmv_of_x0.tms[i].remainder.remove_midpoint(c);
		const_of_x0[i] += c;
	}

	TaylorModelVec<Real> tmv_c0(const_of_x0, rangeDimExt);

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
		// compute the polynomial part under the linear transformation
		std::vector<Polynomial<Real> > initial_linear = symbolic_remainder.Phi_L[0] * symbolic_remainder.polynomial_of_initial_set;

		// compute the other part
		std::vector<Interval> tmvPolyRange;
		tmv.polyRangeNormal(tmvPolyRange, tm_setting.step_end_exp_table);
		x0_other.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, tm_setting.step_end_exp_table, domain.size(), tm_setting.order, tm_setting.cutoff_threshold);

		result.tmv.Remainder(J_ip1);

		Matrix<Interval> x0_rem(rangeDim, 1);
		tmv_of_x0.Remainder(x0_rem);
		J_ip1 += x0_rem;

		for(int i=0; i<rangeDim; ++i)
		{
			result.tmv.tms[i].expansion += initial_linear[i];
		}

		// contract J_ip1 and J_i
		if(invariant.size() > 0)
		{
			std::vector<Interval> polyRangeOfx0;
			result.tmv.polyRangeNormal(polyRangeOfx0, tm_setting.step_end_exp_table);

			std::vector<Interval> intVecTmp(rangeDim);
			std::vector<Interval> original_remainders(rangeDim);
			std::vector<Interval> contracted_remainders(rangeDim);

			for(int i=0; i<rangeDim; ++i)
			{
				intVecTmp[i] = polyRangeOfx0[i] + const_of_x0[i];
				contracted_remainders[i] = original_remainders[i] = J_ip1[i][0] + J_i[i][0];
			}

			int res = contract_remainder(intVecTmp, contracted_remainders, invariant);

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

				result.tmv.tms[i].remainder = contracted_remainders[i];
				range_of_x0.push_back(polyRangeOfx0[i] + result.tmv.tms[i].remainder);
			}
		}
		else
		{
			for(int i=0; i<rangeDim; ++i)
			{
				result.tmv.tms[i].remainder = J_ip1[i][0] + J_i[i][0];
			}

			result.tmv.intEvalNormal(range_of_x0, tm_setting.step_end_exp_table);
		}
	}
	else
	{
		std::vector<Interval> tmvPolyRange;
		tmv.polyRangeNormal(tmvPolyRange, tm_setting.step_end_exp_table);
		tmv_of_x0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, tm_setting.step_end_exp_table, domain.size(), tm_setting.order, tm_setting.cutoff_threshold);

		// contract J_ip1
		if(invariant.size() > 0)
		{
			std::vector<Interval> polyRangeOfx0;
			result.tmv.polyRangeNormal(polyRangeOfx0, tm_setting.step_end_exp_table);

			std::vector<Interval> intVecTmp(rangeDim);
			for(int i=0; i<rangeDim; ++i)
			{
				intVecTmp[i] = polyRangeOfx0[i] + const_of_x0[i];
			}

			std::vector<Interval> contracted_remainders(rangeDim);
			for(int i=0; i<rangeDim; ++i)
			{
				contracted_remainders[i] = result.tmv.tms[i].remainder;
			}

			int res = contract_remainder(intVecTmp, contracted_remainders, invariant);

			if(res < 0)
			{
				return -1;
			}

			for(int i=0; i<rangeDim; ++i)
			{
				result.tmv.tms[i].remainder = contracted_remainders[i];
				range_of_x0.push_back(polyRangeOfx0[i] + result.tmv.tms[i].remainder);
			}
		}
		else
		{
			result.tmv.intEvalNormal(range_of_x0, tm_setting.step_end_exp_table);
		}

		result.tmv.Remainder(J_ip1);
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

	result.tmv.scale_assign(invS);

	TaylorModelVec<Real> new_x0(S);
	new_x0 += tmv_c0;
	TaylorModelVec<Real> x = new_x0;

	for(unsigned int i=1; i<=tm_setting.order; ++i)
	{
		x.Picard_no_remainder_assign(new_x0, ode, rangeDimExt, i, tm_setting.cutoff_threshold);
	}

	bool bfound = true;

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tm_setting.remainder_estimation[i];
	}

	TaylorModelVec<Real> tmvTmp;
	std::list<Interval> intermediate_ranges;

	x.Picard_ctrunc_normal(tmvTmp, new_x0, ode, tm_setting.step_exp_table, rangeDimExt, tm_setting.order, tm_setting.cutoff_threshold, intermediate_ranges, g_setting);

	// compute the interval evaluation of the polynomial difference due to the roundoff error
	std::vector<Interval> intDifferences;
	for(unsigned int i=0; i<rangeDim; ++i)
	{
		Polynomial<Real> polyTmp;
		polyTmp = tmvTmp.tms[i].expansion - x.tms[i].expansion;

		Interval I;
		polyTmp.intEvalNormal(I, tm_setting.step_exp_table);

		intDifferences.push_back(I);
	}

	// add the uncertainties and the cutoff intervals onto the result
	for(unsigned int i=0; i<rangeDim; ++i)
	{
		tmvTmp.tms[i].remainder += intDifferences[i];

		if( ! tmvTmp.tms[i].remainder.subseteq(x.tms[i].remainder) )
		{
			bfound = false;
			break;
		}
	}

	if(!bfound)
	{
		return 0;
	}
	else
	{
		for(unsigned int i=0; i<rangeDim; ++i)
		{
			x.tms[i].remainder = tmvTmp.tms[i].remainder;
		}
	}

	bool bfinished = false;
	for(int rSteps = 0; !bfinished && (rSteps <= MAX_REFINEMENT_STEPS); ++rSteps)
	{
		bfinished = true;

		std::vector<Interval> newRemainders;
		x.Picard_ctrunc_normal_remainder(newRemainders, ode, tm_setting.step_exp_table[1], tm_setting.order, intermediate_ranges, g_setting);

		// add the uncertainties and the cutoff intervals onto the result
		for(unsigned int i=0; i<rangeDim; ++i)
		{
			newRemainders[i] += intDifferences[i];

			if(newRemainders[i].subseteq(x.tms[i].remainder))
			{
				if(x.tms[i].remainder.widthRatio(newRemainders[i]) <= STOP_RATIO)
				{
					bfinished = false;
				}

				x.tms[i].remainder = newRemainders[i];
			}
			else
			{
				bfinished = true;
				break;
			}
		}
	}

	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = tm_setting.step_exp_table[1];

	return 1;
}

int Flowpipe::advance_nondeterministic(Flowpipe & result, const std::vector<Expression_AST<Interval> > & ode, const Taylor_Model_Computation_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Computation_Setting & g_setting, Symbolic_Remainder & symbolic_remainder) const
{
	unsigned int rangeDim = ode.size();
	unsigned int rangeDimExt = rangeDim + 1;
	Interval intUnit(-1,1);
	result.clear();

	// evaluate the the initial set x0
	TaylorModelVec<Real> tmv_of_x0;
	tmvPre.evaluate_time(tmv_of_x0, tm_setting.step_end_exp_table);

	// the center point of x0's polynomial part
	std::vector<Real> const_of_x0;
	tmv_of_x0.constant(const_of_x0);

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		Real c;
		tmv_of_x0.tms[i].remainder.remove_midpoint(c);
		const_of_x0[i] += c;
	}

	TaylorModelVec<Real> tmv_c0(const_of_x0, rangeDimExt);

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
		// compute the polynomial part under the linear transformation
		std::vector<Polynomial<Real> > initial_linear = symbolic_remainder.Phi_L[0] * symbolic_remainder.polynomial_of_initial_set;

		// compute the other part
		std::vector<Interval> tmvPolyRange;
		tmv.polyRangeNormal(tmvPolyRange, tm_setting.step_end_exp_table);
		x0_other.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, tm_setting.step_end_exp_table, domain.size(), tm_setting.order, tm_setting.cutoff_threshold);

		result.tmv.Remainder(J_ip1);

		Matrix<Interval> x0_rem(rangeDim, 1);
		tmv_of_x0.Remainder(x0_rem);
		J_ip1 += x0_rem;

		for(int i=0; i<rangeDim; ++i)
		{
			result.tmv.tms[i].expansion += initial_linear[i];
		}

		// contract J_ip1 and J_i
		if(invariant.size() > 0)
		{
			std::vector<Interval> polyRangeOfx0;
			result.tmv.polyRangeNormal(polyRangeOfx0, tm_setting.step_end_exp_table);

			std::vector<Interval> intVecTmp(rangeDim);
			std::vector<Interval> original_remainders(rangeDim);
			std::vector<Interval> contracted_remainders(rangeDim);

			for(int i=0; i<rangeDim; ++i)
			{
				intVecTmp[i] = polyRangeOfx0[i] + const_of_x0[i];
				contracted_remainders[i] = original_remainders[i] = J_ip1[i][0] + J_i[i][0];
			}

			int res = contract_remainder(intVecTmp, contracted_remainders, invariant);

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

				result.tmv.tms[i].remainder = contracted_remainders[i];
				range_of_x0.push_back(polyRangeOfx0[i] + result.tmv.tms[i].remainder);
			}
		}
		else
		{
			for(int i=0; i<rangeDim; ++i)
			{
				result.tmv.tms[i].remainder = J_ip1[i][0] + J_i[i][0];
			}

			result.tmv.intEvalNormal(range_of_x0, tm_setting.step_end_exp_table);
		}
	}
	else
	{
		std::vector<Interval> tmvPolyRange;
		tmv.polyRangeNormal(tmvPolyRange, tm_setting.step_end_exp_table);
		tmv_of_x0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, tm_setting.step_end_exp_table, domain.size(), tm_setting.order, tm_setting.cutoff_threshold);

		// contract J_ip1
		if(invariant.size() > 0)
		{
			std::vector<Interval> polyRangeOfx0;
			result.tmv.polyRangeNormal(polyRangeOfx0, tm_setting.step_end_exp_table);

			std::vector<Interval> intVecTmp(rangeDim);
			for(int i=0; i<rangeDim; ++i)
			{
				intVecTmp[i] = polyRangeOfx0[i] + const_of_x0[i];
			}

			std::vector<Interval> contracted_remainders(rangeDim);
			for(int i=0; i<rangeDim; ++i)
			{
				contracted_remainders[i] = result.tmv.tms[i].remainder;
			}

			int res = contract_remainder(intVecTmp, contracted_remainders, invariant);

			if(res < 0)
			{
				return -1;
			}

			for(int i=0; i<rangeDim; ++i)
			{
				result.tmv.tms[i].remainder = contracted_remainders[i];
				range_of_x0.push_back(polyRangeOfx0[i] + result.tmv.tms[i].remainder);
			}
		}
		else
		{
			result.tmv.intEvalNormal(range_of_x0, tm_setting.step_end_exp_table);
		}

		result.tmv.Remainder(J_ip1);
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

	result.tmv.scale_assign(invS);

	TaylorModelVec<Real> new_x0(S);
	new_x0 += tmv_c0;
	TaylorModelVec<Real> x = new_x0;

	for(unsigned int i=1; i<=tm_setting.order; ++i)
	{
		x.Picard_no_remainder_assign(new_x0, ode, rangeDimExt, i, tm_setting.cutoff_threshold);
	}

	bool bfound = true;

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tm_setting.remainder_estimation[i];
	}

	TaylorModelVec<Real> tmvTmp;
	std::list<Interval> intermediate_ranges;

	x.Picard_ctrunc_normal(tmvTmp, new_x0, ode, tm_setting.step_exp_table, rangeDimExt, tm_setting.order, tm_setting.cutoff_threshold, intermediate_ranges, g_setting);

	// compute the interval evaluation of the polynomial difference due to the roundoff error
	std::vector<Interval> intDifferences;
	for(unsigned int i=0; i<rangeDim; ++i)
	{
		Polynomial<Real> polyTmp;
		polyTmp = tmvTmp.tms[i].expansion - x.tms[i].expansion;

		Interval I;
		polyTmp.intEvalNormal(I, tm_setting.step_exp_table);

		intDifferences.push_back(I);
	}

	// add the uncertainties and the cutoff intervals onto the result
	for(unsigned int i=0; i<rangeDim; ++i)
	{
		tmvTmp.tms[i].remainder += intDifferences[i];

		if( ! tmvTmp.tms[i].remainder.subseteq(x.tms[i].remainder) )
		{
			bfound = false;
			break;
		}
	}

	if(!bfound)
	{
		return 0;
	}
	else
	{
		for(unsigned int i=0; i<rangeDim; ++i)
		{
			x.tms[i].remainder = tmvTmp.tms[i].remainder;
		}
	}

	bool bfinished = false;
	for(int rSteps = 0; !bfinished && (rSteps <= MAX_REFINEMENT_STEPS); ++rSteps)
	{
		bfinished = true;

		std::vector<Interval> newRemainders;
		x.Picard_ctrunc_normal_remainder(newRemainders, ode, tm_setting.step_exp_table[1], tm_setting.order, intermediate_ranges, g_setting);

		// add the uncertainties and the cutoff intervals onto the result
		for(unsigned int i=0; i<rangeDim; ++i)
		{
			newRemainders[i] += intDifferences[i];

			if(newRemainders[i].subseteq(x.tms[i].remainder))
			{
				if(x.tms[i].remainder.widthRatio(newRemainders[i]) <= STOP_RATIO)
				{
					bfinished = false;
				}

				x.tms[i].remainder = newRemainders[i];
			}
			else
			{
				bfinished = true;
				break;
			}
		}
	}

	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = tm_setting.step_exp_table[1];

	return 1;
}

int Flowpipe::advance_deterministic_adaptive_stepsize(Flowpipe & result, const std::vector<Expression_AST<Real> > & ode, Taylor_Model_Computation_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Computation_Setting & g_setting, Symbolic_Remainder & symbolic_remainder) const
{
	unsigned int rangeDim = ode.size();
	unsigned int rangeDimExt = rangeDim + 1;
	Interval intUnit(-1,1);
	result.clear();

	// evaluate the the initial set x0
	TaylorModelVec<Real> tmv_of_x0;
	tmvPre.evaluate_time(tmv_of_x0, tm_setting.step_end_exp_table);

	// the center point of x0's polynomial part
	std::vector<Real> const_of_x0;
	tmv_of_x0.constant(const_of_x0);

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		Real c;
		tmv_of_x0.tms[i].remainder.remove_midpoint(c);
		const_of_x0[i] += c;
	}

	TaylorModelVec<Real> tmv_c0(const_of_x0, rangeDimExt);

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
		// compute the polynomial part under the linear transformation
		std::vector<Polynomial<Real> > initial_linear = symbolic_remainder.Phi_L[0] * symbolic_remainder.polynomial_of_initial_set;

		// compute the other part
		std::vector<Interval> tmvPolyRange;
		tmv.polyRangeNormal(tmvPolyRange, tm_setting.step_end_exp_table);
		x0_other.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, tm_setting.step_end_exp_table, domain.size(), tm_setting.order, tm_setting.cutoff_threshold);

		result.tmv.Remainder(J_ip1);

		Matrix<Interval> x0_rem(rangeDim, 1);
		tmv_of_x0.Remainder(x0_rem);
		J_ip1 += x0_rem;

		for(int i=0; i<rangeDim; ++i)
		{
			result.tmv.tms[i].expansion += initial_linear[i];
		}

		// contract J_ip1 and J_i
		if(invariant.size() > 0)
		{
			std::vector<Interval> polyRangeOfx0;
			result.tmv.polyRangeNormal(polyRangeOfx0, tm_setting.step_end_exp_table);

			std::vector<Interval> intVecTmp(rangeDim);
			std::vector<Interval> original_remainders(rangeDim);
			std::vector<Interval> contracted_remainders(rangeDim);

			for(int i=0; i<rangeDim; ++i)
			{
				intVecTmp[i] = polyRangeOfx0[i] + const_of_x0[i];
				contracted_remainders[i] = original_remainders[i] = J_ip1[i][0] + J_i[i][0];
			}

			int res = contract_remainder(intVecTmp, contracted_remainders, invariant);

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

				result.tmv.tms[i].remainder = contracted_remainders[i];
				range_of_x0.push_back(polyRangeOfx0[i] + result.tmv.tms[i].remainder);
			}
		}
		else
		{
			for(int i=0; i<rangeDim; ++i)
			{
				result.tmv.tms[i].remainder = J_ip1[i][0] + J_i[i][0];
			}

			result.tmv.intEvalNormal(range_of_x0, tm_setting.step_end_exp_table);
		}
	}
	else
	{
		std::vector<Interval> tmvPolyRange;
		tmv.polyRangeNormal(tmvPolyRange, tm_setting.step_end_exp_table);
		tmv_of_x0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, tm_setting.step_end_exp_table, domain.size(), tm_setting.order, tm_setting.cutoff_threshold);

		// contract J_ip1
		if(invariant.size() > 0)
		{
			std::vector<Interval> polyRangeOfx0;
			result.tmv.polyRangeNormal(polyRangeOfx0, tm_setting.step_end_exp_table);

			std::vector<Interval> intVecTmp(rangeDim);
			for(int i=0; i<rangeDim; ++i)
			{
				intVecTmp[i] = polyRangeOfx0[i] + const_of_x0[i];
			}

			std::vector<Interval> contracted_remainders(rangeDim);
			for(int i=0; i<rangeDim; ++i)
			{
				contracted_remainders[i] = result.tmv.tms[i].remainder;
			}

			int res = contract_remainder(intVecTmp, contracted_remainders, invariant);

			if(res < 0)
			{
				return -1;
			}

			for(int i=0; i<rangeDim; ++i)
			{
				result.tmv.tms[i].remainder = contracted_remainders[i];
				range_of_x0.push_back(polyRangeOfx0[i] + result.tmv.tms[i].remainder);
			}
		}
		else
		{
			result.tmv.intEvalNormal(range_of_x0, tm_setting.step_end_exp_table);
		}

		result.tmv.Remainder(J_ip1);
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

	result.tmv.scale_assign(invS);

	TaylorModelVec<Real> new_x0(S);
	new_x0 += tmv_c0;
	TaylorModelVec<Real> x = new_x0;

	for(unsigned int i=1; i<=tm_setting.order; ++i)
	{
		x.Picard_no_remainder_assign(new_x0, ode, rangeDimExt, i, tm_setting.cutoff_threshold);
	}

	bool bfound = true;

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tm_setting.remainder_estimation[i];
	}

	TaylorModelVec<Real> tmvTmp;
	std::list<Interval> intermediate_ranges;

	x.Picard_ctrunc_normal(tmvTmp, new_x0, ode, tm_setting.step_exp_table, rangeDimExt, tm_setting.order, tm_setting.cutoff_threshold, intermediate_ranges, g_setting);

	// compute the interval evaluation of the polynomial difference due to the roundoff error
	std::vector<Interval> intDifferences;
	std::vector<Polynomial<Real> > polyDifferences;
	for(unsigned int i=0; i<rangeDim; ++i)
	{
		Polynomial<Real> polyTmp;
		polyTmp = tmvTmp.tms[i].expansion - x.tms[i].expansion;
		polyDifferences.push_back(polyTmp);

		Interval I;
		polyTmp.intEvalNormal(I, tm_setting.step_exp_table);

		intDifferences.push_back(I);
	}

	// add the uncertainties and the cutoff intervals onto the result
	for(unsigned int i=0; i<rangeDim; ++i)
	{
		tmvTmp.tms[i].remainder += intDifferences[i];

		if( ! tmvTmp.tms[i].remainder.subseteq(x.tms[i].remainder) )
		{
			bfound = false;
			break;
		}
	}

	for(; !bfound;)
	{
		bfound = true;
		double newStep = tm_setting.step_exp_table[1].sup() * LAMBDA_DOWN;	// reduce the time step size

		if(newStep < tm_setting.step_min)
		{
			return 0;
		}

		tm_setting.setStepsize(newStep, tm_setting.order);

		intermediate_ranges.clear();
		x.Picard_ctrunc_normal(tmvTmp, new_x0, ode, tm_setting.step_exp_table, rangeDimExt, tm_setting.order, tm_setting.cutoff_threshold, intermediate_ranges, g_setting);

		// recompute the interval evaluation of the polynomial differences
		for(int i=0; i<rangeDim; ++i)
		{
			polyDifferences[i].intEvalNormal(intDifferences[i], tm_setting.step_exp_table);

			tmvTmp.tms[i].remainder += intDifferences[i];

			if( ! tmvTmp.tms[i].remainder.subseteq(x.tms[i].remainder) )
			{
				bfound = false;
				break;
			}
		}
	}

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tmvTmp.tms[i].remainder;
	}

	bool bfinished = false;
	for(int rSteps = 0; !bfinished && (rSteps <= MAX_REFINEMENT_STEPS); ++rSteps)
	{
		bfinished = true;

		std::vector<Interval> newRemainders;
		x.Picard_ctrunc_normal_remainder(newRemainders, ode, tm_setting.step_exp_table[1], tm_setting.order, intermediate_ranges, g_setting);

		// add the uncertainties and the cutoff intervals onto the result
		for(unsigned int i=0; i<rangeDim; ++i)
		{
			newRemainders[i] += intDifferences[i];

			if(newRemainders[i].subseteq(x.tms[i].remainder))
			{
				if(x.tms[i].remainder.widthRatio(newRemainders[i]) <= STOP_RATIO)
				{
					bfinished = false;
				}

				x.tms[i].remainder = newRemainders[i];
			}
			else
			{
				bfinished = true;
				break;
			}
		}
	}

	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = tm_setting.step_exp_table[1];

	return 1;
}

int Flowpipe::advance_nondeterministic_adaptive_stepsize(Flowpipe & result, const std::vector<Expression_AST<Interval> > & ode, Taylor_Model_Computation_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Computation_Setting & g_setting, Symbolic_Remainder & symbolic_remainder) const
{
	unsigned int rangeDim = ode.size();
	unsigned int rangeDimExt = rangeDim + 1;
	Interval intUnit(-1,1);
	result.clear();

	// evaluate the the initial set x0
	TaylorModelVec<Real> tmv_of_x0;
	tmvPre.evaluate_time(tmv_of_x0, tm_setting.step_end_exp_table);

	// the center point of x0's polynomial part
	std::vector<Real> const_of_x0;
	tmv_of_x0.constant(const_of_x0);

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		Real c;
		tmv_of_x0.tms[i].remainder.remove_midpoint(c);
		const_of_x0[i] += c;
	}

	TaylorModelVec<Real> tmv_c0(const_of_x0, rangeDimExt);

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
		// compute the polynomial part under the linear transformation
		std::vector<Polynomial<Real> > initial_linear = symbolic_remainder.Phi_L[0] * symbolic_remainder.polynomial_of_initial_set;

		// compute the other part
		std::vector<Interval> tmvPolyRange;
		tmv.polyRangeNormal(tmvPolyRange, tm_setting.step_end_exp_table);
		x0_other.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, tm_setting.step_end_exp_table, domain.size(), tm_setting.order, tm_setting.cutoff_threshold);

		result.tmv.Remainder(J_ip1);

		Matrix<Interval> x0_rem(rangeDim, 1);
		tmv_of_x0.Remainder(x0_rem);
		J_ip1 += x0_rem;

		for(int i=0; i<rangeDim; ++i)
		{
			result.tmv.tms[i].expansion += initial_linear[i];
		}

		// contract J_ip1 and J_i
		if(invariant.size() > 0)
		{
			std::vector<Interval> polyRangeOfx0;
			result.tmv.polyRangeNormal(polyRangeOfx0, tm_setting.step_end_exp_table);

			std::vector<Interval> intVecTmp(rangeDim);
			std::vector<Interval> original_remainders(rangeDim);
			std::vector<Interval> contracted_remainders(rangeDim);

			for(int i=0; i<rangeDim; ++i)
			{
				intVecTmp[i] = polyRangeOfx0[i] + const_of_x0[i];
				contracted_remainders[i] = original_remainders[i] = J_ip1[i][0] + J_i[i][0];
			}

			int res = contract_remainder(intVecTmp, contracted_remainders, invariant);

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

				result.tmv.tms[i].remainder = contracted_remainders[i];
				range_of_x0.push_back(polyRangeOfx0[i] + result.tmv.tms[i].remainder);
			}
		}
		else
		{
			for(int i=0; i<rangeDim; ++i)
			{
				result.tmv.tms[i].remainder = J_ip1[i][0] + J_i[i][0];
			}

			result.tmv.intEvalNormal(range_of_x0, tm_setting.step_end_exp_table);
		}
	}
	else
	{
		std::vector<Interval> tmvPolyRange;
		tmv.polyRangeNormal(tmvPolyRange, tm_setting.step_end_exp_table);
		tmv_of_x0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, tm_setting.step_end_exp_table, domain.size(), tm_setting.order, tm_setting.cutoff_threshold);

		// contract J_ip1
		if(invariant.size() > 0)
		{
			std::vector<Interval> polyRangeOfx0;
			result.tmv.polyRangeNormal(polyRangeOfx0, tm_setting.step_end_exp_table);

			std::vector<Interval> intVecTmp(rangeDim);
			for(int i=0; i<rangeDim; ++i)
			{
				intVecTmp[i] = polyRangeOfx0[i] + const_of_x0[i];
			}

			std::vector<Interval> contracted_remainders(rangeDim);
			for(int i=0; i<rangeDim; ++i)
			{
				contracted_remainders[i] = result.tmv.tms[i].remainder;
			}

			int res = contract_remainder(intVecTmp, contracted_remainders, invariant);

			if(res < 0)
			{
				return -1;
			}

			for(int i=0; i<rangeDim; ++i)
			{
				result.tmv.tms[i].remainder = contracted_remainders[i];
				range_of_x0.push_back(polyRangeOfx0[i] + result.tmv.tms[i].remainder);
			}
		}
		else
		{
			result.tmv.intEvalNormal(range_of_x0, tm_setting.step_end_exp_table);
		}

		result.tmv.Remainder(J_ip1);
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

	result.tmv.scale_assign(invS);

	TaylorModelVec<Real> new_x0(S);
	new_x0 += tmv_c0;
	TaylorModelVec<Real> x = new_x0;

	for(unsigned int i=1; i<=tm_setting.order; ++i)
	{
		x.Picard_no_remainder_assign(new_x0, ode, rangeDimExt, i, tm_setting.cutoff_threshold);
	}

	bool bfound = true;

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tm_setting.remainder_estimation[i];
	}

	TaylorModelVec<Real> tmvTmp;
	std::list<Interval> intermediate_ranges;

	x.Picard_ctrunc_normal(tmvTmp, new_x0, ode, tm_setting.step_exp_table, rangeDimExt, tm_setting.order, tm_setting.cutoff_threshold, intermediate_ranges, g_setting);

	// compute the interval evaluation of the polynomial difference due to the roundoff error
	std::vector<Interval> intDifferences;
	std::vector<Polynomial<Real> > polyDifferences;
	for(unsigned int i=0; i<rangeDim; ++i)
	{
		Polynomial<Real> polyTmp;
		polyTmp = tmvTmp.tms[i].expansion - x.tms[i].expansion;
		polyDifferences.push_back(polyTmp);

		Interval I;
		polyTmp.intEvalNormal(I, tm_setting.step_exp_table);

		intDifferences.push_back(I);
	}

	// add the uncertainties and the cutoff intervals onto the result
	for(unsigned int i=0; i<rangeDim; ++i)
	{
		tmvTmp.tms[i].remainder += intDifferences[i];

		if( ! tmvTmp.tms[i].remainder.subseteq(x.tms[i].remainder) )
		{
			bfound = false;
			break;
		}
	}

	for(; !bfound;)
	{
		bfound = true;
		double newStep = tm_setting.step_exp_table[1].sup() * LAMBDA_DOWN;	// reduce the time step size

		if(newStep < tm_setting.step_min)
		{
			return 0;
		}

		tm_setting.setStepsize(newStep, tm_setting.order);

		intermediate_ranges.clear();
		x.Picard_ctrunc_normal(tmvTmp, new_x0, ode, tm_setting.step_exp_table, rangeDimExt, tm_setting.order, tm_setting.cutoff_threshold, intermediate_ranges, g_setting);

		// recompute the interval evaluation of the polynomial differences
		for(int i=0; i<rangeDim; ++i)
		{
			polyDifferences[i].intEvalNormal(intDifferences[i], tm_setting.step_exp_table);

			tmvTmp.tms[i].remainder += intDifferences[i];

			if( ! tmvTmp.tms[i].remainder.subseteq(x.tms[i].remainder) )
			{
				bfound = false;
				break;
			}
		}
	}

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tmvTmp.tms[i].remainder;
	}

	bool bfinished = false;
	for(int rSteps = 0; !bfinished && (rSteps <= MAX_REFINEMENT_STEPS); ++rSteps)
	{
		bfinished = true;

		std::vector<Interval> newRemainders;
		x.Picard_ctrunc_normal_remainder(newRemainders, ode, tm_setting.step_exp_table[1], tm_setting.order, intermediate_ranges, g_setting);

		// add the uncertainties and the cutoff intervals onto the result
		for(unsigned int i=0; i<rangeDim; ++i)
		{
			newRemainders[i] += intDifferences[i];

			if(newRemainders[i].subseteq(x.tms[i].remainder))
			{
				if(x.tms[i].remainder.widthRatio(newRemainders[i]) <= STOP_RATIO)
				{
					bfinished = false;
				}

				x.tms[i].remainder = newRemainders[i];
			}
			else
			{
				bfinished = true;
				break;
			}
		}
	}

	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = tm_setting.step_exp_table[1];

	return 1;
}

int Flowpipe::advance_deterministic_adaptive_order(Flowpipe & result, const std::vector<Expression_AST<Real> > & ode, Taylor_Model_Computation_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Computation_Setting & g_setting, Symbolic_Remainder & symbolic_remainder) const
{
	unsigned int rangeDim = ode.size();
	unsigned int rangeDimExt = rangeDim + 1;
	Interval intUnit(-1,1);
	result.clear();

	// evaluate the the initial set x0
	TaylorModelVec<Real> tmv_of_x0;
	tmvPre.evaluate_time(tmv_of_x0, tm_setting.step_end_exp_table);

	// the center point of x0's polynomial part
	std::vector<Real> const_of_x0;
	tmv_of_x0.constant(const_of_x0);

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		Real c;
		tmv_of_x0.tms[i].remainder.remove_midpoint(c);
		const_of_x0[i] += c;
	}

	TaylorModelVec<Real> tmv_c0(const_of_x0, rangeDimExt);

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
		// compute the polynomial part under the linear transformation
		std::vector<Polynomial<Real> > initial_linear = symbolic_remainder.Phi_L[0] * symbolic_remainder.polynomial_of_initial_set;

		// compute the other part
		std::vector<Interval> tmvPolyRange;
		tmv.polyRangeNormal(tmvPolyRange, tm_setting.step_end_exp_table);
		x0_other.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, tm_setting.step_end_exp_table, domain.size(), tm_setting.order, tm_setting.cutoff_threshold);

		result.tmv.Remainder(J_ip1);

		Matrix<Interval> x0_rem(rangeDim, 1);
		tmv_of_x0.Remainder(x0_rem);
		J_ip1 += x0_rem;

		for(int i=0; i<rangeDim; ++i)
		{
			result.tmv.tms[i].expansion += initial_linear[i];
		}

		// contract J_ip1 and J_i
		if(invariant.size() > 0)
		{
			std::vector<Interval> polyRangeOfx0;
			result.tmv.polyRangeNormal(polyRangeOfx0, tm_setting.step_end_exp_table);

			std::vector<Interval> intVecTmp(rangeDim);
			std::vector<Interval> original_remainders(rangeDim);
			std::vector<Interval> contracted_remainders(rangeDim);

			for(int i=0; i<rangeDim; ++i)
			{
				intVecTmp[i] = polyRangeOfx0[i] + const_of_x0[i];
				contracted_remainders[i] = original_remainders[i] = J_ip1[i][0] + J_i[i][0];
			}

			int res = contract_remainder(intVecTmp, contracted_remainders, invariant);

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

				result.tmv.tms[i].remainder = contracted_remainders[i];
				range_of_x0.push_back(polyRangeOfx0[i] + result.tmv.tms[i].remainder);
			}
		}
		else
		{
			for(int i=0; i<rangeDim; ++i)
			{
				result.tmv.tms[i].remainder = J_ip1[i][0] + J_i[i][0];
			}

			result.tmv.intEvalNormal(range_of_x0, tm_setting.step_end_exp_table);
		}
	}
	else
	{
		std::vector<Interval> tmvPolyRange;
		tmv.polyRangeNormal(tmvPolyRange, tm_setting.step_end_exp_table);
		tmv_of_x0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, tm_setting.step_end_exp_table, domain.size(), tm_setting.order, tm_setting.cutoff_threshold);

		// contract J_ip1
		if(invariant.size() > 0)
		{
			std::vector<Interval> polyRangeOfx0;
			result.tmv.polyRangeNormal(polyRangeOfx0, tm_setting.step_end_exp_table);

			std::vector<Interval> intVecTmp(rangeDim);
			for(int i=0; i<rangeDim; ++i)
			{
				intVecTmp[i] = polyRangeOfx0[i] + const_of_x0[i];
			}

			std::vector<Interval> contracted_remainders(rangeDim);
			for(int i=0; i<rangeDim; ++i)
			{
				contracted_remainders[i] = result.tmv.tms[i].remainder;
			}

			int res = contract_remainder(intVecTmp, contracted_remainders, invariant);

			if(res < 0)
			{
				return -1;
			}

			for(int i=0; i<rangeDim; ++i)
			{
				result.tmv.tms[i].remainder = contracted_remainders[i];
				range_of_x0.push_back(polyRangeOfx0[i] + result.tmv.tms[i].remainder);
			}
		}
		else
		{
			result.tmv.intEvalNormal(range_of_x0, tm_setting.step_end_exp_table);
		}

		result.tmv.Remainder(J_ip1);
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

	result.tmv.scale_assign(invS);

	TaylorModelVec<Real> new_x0(S);
	new_x0 += tmv_c0;
	TaylorModelVec<Real> x = new_x0;

	for(unsigned int i=1; i<=tm_setting.order; ++i)
	{
		x.Picard_no_remainder_assign(new_x0, ode, rangeDimExt, i, tm_setting.cutoff_threshold);
	}

	bool bfound = true;

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tm_setting.remainder_estimation[i];
	}

	TaylorModelVec<Real> tmvTmp;
	std::list<Interval> intermediate_ranges;

	x.Picard_ctrunc_normal(tmvTmp, new_x0, ode, tm_setting.step_exp_table, rangeDimExt, tm_setting.order, tm_setting.cutoff_threshold, intermediate_ranges, g_setting);

	// compute the interval evaluation of the polynomial difference due to the roundoff error
	std::vector<Interval> intDifferences;
	for(unsigned int i=0; i<rangeDim; ++i)
	{
		Polynomial<Real> polyTmp;
		polyTmp = tmvTmp.tms[i].expansion - x.tms[i].expansion;

		Interval I;
		polyTmp.intEvalNormal(I, tm_setting.step_exp_table);

		intDifferences.push_back(I);
	}

	// add the uncertainties and the cutoff intervals onto the result
	for(unsigned int i=0; i<rangeDim; ++i)
	{
		tmvTmp.tms[i].remainder += intDifferences[i];

		if( ! tmvTmp.tms[i].remainder.subseteq(x.tms[i].remainder) )
		{
			bfound = false;
			break;
		}
	}

	for(; !bfound;)
	{
		bfound = true;

		++tm_setting.order;

		if(tm_setting.order > tm_setting.order_max)
		{
			return 0;
		}

		x.Picard_no_remainder_assign(new_x0, ode, rangeDimExt, tm_setting.order, tm_setting.cutoff_threshold);

		for(int i=0; i<rangeDim; ++i)
		{
			x.tms[i].remainder = tm_setting.remainder_estimation[i];
		}

		intermediate_ranges.clear();
		x.Picard_ctrunc_normal(tmvTmp, new_x0, ode, tm_setting.step_exp_table, rangeDimExt, tm_setting.order, tm_setting.cutoff_threshold, intermediate_ranges, g_setting);

		// recompute the interval evaluation of the polynomial differences
		for(int i=0; i<rangeDim; ++i)
		{
			Polynomial<Real> polyTmp;
			polyTmp = tmvTmp.tms[i].expansion - x.tms[i].expansion;

			Interval I;
			polyTmp.intEvalNormal(I, tm_setting.step_exp_table);

			intDifferences[i] = I;

			tmvTmp.tms[i].remainder += intDifferences[i];

			if( ! tmvTmp.tms[i].remainder.subseteq(x.tms[i].remainder) )
			{
				bfound = false;
				break;
			}
		}
	}

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tmvTmp.tms[i].remainder;
	}

	bool bfinished = false;
	for(int rSteps = 0; !bfinished && (rSteps <= MAX_REFINEMENT_STEPS); ++rSteps)
	{
		bfinished = true;

		std::vector<Interval> newRemainders;
		x.Picard_ctrunc_normal_remainder(newRemainders, ode, tm_setting.step_exp_table[1], tm_setting.order, intermediate_ranges, g_setting);

		// add the uncertainties and the cutoff intervals onto the result
		for(unsigned int i=0; i<rangeDim; ++i)
		{
			newRemainders[i] += intDifferences[i];

			if(newRemainders[i].subseteq(x.tms[i].remainder))
			{
				if(x.tms[i].remainder.widthRatio(newRemainders[i]) <= STOP_RATIO)
				{
					bfinished = false;
				}

				x.tms[i].remainder = newRemainders[i];
			}
			else
			{
				bfinished = true;
				break;
			}
		}
	}

	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = tm_setting.step_exp_table[1];

	return 1;
}

int Flowpipe::advance_nondeterministic_adaptive_order(Flowpipe & result, const std::vector<Expression_AST<Interval> > & ode, Taylor_Model_Computation_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Computation_Setting & g_setting, Symbolic_Remainder & symbolic_remainder) const
{
	unsigned int rangeDim = ode.size();
	unsigned int rangeDimExt = rangeDim + 1;
	Interval intUnit(-1,1);
	result.clear();

	// evaluate the the initial set x0
	TaylorModelVec<Real> tmv_of_x0;
	tmvPre.evaluate_time(tmv_of_x0, tm_setting.step_end_exp_table);

	// the center point of x0's polynomial part
	std::vector<Real> const_of_x0;
	tmv_of_x0.constant(const_of_x0);

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		Real c;
		tmv_of_x0.tms[i].remainder.remove_midpoint(c);
		const_of_x0[i] += c;
	}

	TaylorModelVec<Real> tmv_c0(const_of_x0, rangeDimExt);

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
		// compute the polynomial part under the linear transformation
		std::vector<Polynomial<Real> > initial_linear = symbolic_remainder.Phi_L[0] * symbolic_remainder.polynomial_of_initial_set;

		// compute the other part
		std::vector<Interval> tmvPolyRange;
		tmv.polyRangeNormal(tmvPolyRange, tm_setting.step_end_exp_table);
		x0_other.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, tm_setting.step_end_exp_table, domain.size(), tm_setting.order, tm_setting.cutoff_threshold);

		result.tmv.Remainder(J_ip1);

		Matrix<Interval> x0_rem(rangeDim, 1);
		tmv_of_x0.Remainder(x0_rem);
		J_ip1 += x0_rem;

		for(int i=0; i<rangeDim; ++i)
		{
			result.tmv.tms[i].expansion += initial_linear[i];
		}

		// contract J_ip1 and J_i
		if(invariant.size() > 0)
		{
			std::vector<Interval> polyRangeOfx0;
			result.tmv.polyRangeNormal(polyRangeOfx0, tm_setting.step_end_exp_table);

			std::vector<Interval> intVecTmp(rangeDim);
			std::vector<Interval> original_remainders(rangeDim);
			std::vector<Interval> contracted_remainders(rangeDim);

			for(int i=0; i<rangeDim; ++i)
			{
				intVecTmp[i] = polyRangeOfx0[i] + const_of_x0[i];
				contracted_remainders[i] = original_remainders[i] = J_ip1[i][0] + J_i[i][0];
			}

			int res = contract_remainder(intVecTmp, contracted_remainders, invariant);

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

				result.tmv.tms[i].remainder = contracted_remainders[i];
				range_of_x0.push_back(polyRangeOfx0[i] + result.tmv.tms[i].remainder);
			}
		}
		else
		{
			for(int i=0; i<rangeDim; ++i)
			{
				result.tmv.tms[i].remainder = J_ip1[i][0] + J_i[i][0];
			}

			result.tmv.intEvalNormal(range_of_x0, tm_setting.step_end_exp_table);
		}
	}
	else
	{
		std::vector<Interval> tmvPolyRange;
		tmv.polyRangeNormal(tmvPolyRange, tm_setting.step_end_exp_table);
		tmv_of_x0.insert_ctrunc_normal(result.tmv, tmv, tmvPolyRange, tm_setting.step_end_exp_table, domain.size(), tm_setting.order, tm_setting.cutoff_threshold);

		// contract J_ip1
		if(invariant.size() > 0)
		{
			std::vector<Interval> polyRangeOfx0;
			result.tmv.polyRangeNormal(polyRangeOfx0, tm_setting.step_end_exp_table);

			std::vector<Interval> intVecTmp(rangeDim);
			for(int i=0; i<rangeDim; ++i)
			{
				intVecTmp[i] = polyRangeOfx0[i] + const_of_x0[i];
			}

			std::vector<Interval> contracted_remainders(rangeDim);
			for(int i=0; i<rangeDim; ++i)
			{
				contracted_remainders[i] = result.tmv.tms[i].remainder;
			}

			int res = contract_remainder(intVecTmp, contracted_remainders, invariant);

			if(res < 0)
			{
				return -1;
			}

			for(int i=0; i<rangeDim; ++i)
			{
				result.tmv.tms[i].remainder = contracted_remainders[i];
				range_of_x0.push_back(polyRangeOfx0[i] + result.tmv.tms[i].remainder);
			}
		}
		else
		{
			result.tmv.intEvalNormal(range_of_x0, tm_setting.step_end_exp_table);
		}

		result.tmv.Remainder(J_ip1);
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

	result.tmv.scale_assign(invS);

	TaylorModelVec<Real> new_x0(S);
	new_x0 += tmv_c0;
	TaylorModelVec<Real> x = new_x0;

	for(unsigned int i=1; i<=tm_setting.order; ++i)
	{
		x.Picard_no_remainder_assign(new_x0, ode, rangeDimExt, i, tm_setting.cutoff_threshold);
	}

	bool bfound = true;

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tm_setting.remainder_estimation[i];
	}

	TaylorModelVec<Real> tmvTmp;
	std::list<Interval> intermediate_ranges;

	x.Picard_ctrunc_normal(tmvTmp, new_x0, ode, tm_setting.step_exp_table, rangeDimExt, tm_setting.order, tm_setting.cutoff_threshold, intermediate_ranges, g_setting);

	// compute the interval evaluation of the polynomial difference due to the roundoff error
	std::vector<Interval> intDifferences;
	for(unsigned int i=0; i<rangeDim; ++i)
	{
		Polynomial<Real> polyTmp;
		polyTmp = tmvTmp.tms[i].expansion - x.tms[i].expansion;

		Interval I;
		polyTmp.intEvalNormal(I, tm_setting.step_exp_table);

		intDifferences.push_back(I);
	}

	// add the uncertainties and the cutoff intervals onto the result
	for(unsigned int i=0; i<rangeDim; ++i)
	{
		tmvTmp.tms[i].remainder += intDifferences[i];

		if( ! tmvTmp.tms[i].remainder.subseteq(x.tms[i].remainder) )
		{
			bfound = false;
			break;
		}
	}

	for(; !bfound;)
	{
		bfound = true;

		++tm_setting.order;

		if(tm_setting.order > tm_setting.order_max)
		{
			return 0;
		}

		x.Picard_no_remainder_assign(new_x0, ode, rangeDimExt, tm_setting.order, tm_setting.cutoff_threshold);

		for(int i=0; i<rangeDim; ++i)
		{
			x.tms[i].remainder = tm_setting.remainder_estimation[i];
		}

		intermediate_ranges.clear();
		x.Picard_ctrunc_normal(tmvTmp, new_x0, ode, tm_setting.step_exp_table, rangeDimExt, tm_setting.order, tm_setting.cutoff_threshold, intermediate_ranges, g_setting);

		// recompute the interval evaluation of the polynomial differences
		for(int i=0; i<rangeDim; ++i)
		{
			Polynomial<Real> polyTmp;
			polyTmp = tmvTmp.tms[i].expansion - x.tms[i].expansion;

			Interval I;
			polyTmp.intEvalNormal(I, tm_setting.step_exp_table);

			intDifferences[i] = I;

			tmvTmp.tms[i].remainder += intDifferences[i];

			if( ! tmvTmp.tms[i].remainder.subseteq(x.tms[i].remainder) )
			{
				bfound = false;
				break;
			}
		}
	}

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		x.tms[i].remainder = tmvTmp.tms[i].remainder;
	}

	bool bfinished = false;
	for(int rSteps = 0; !bfinished && (rSteps <= MAX_REFINEMENT_STEPS); ++rSteps)
	{
		bfinished = true;

		std::vector<Interval> newRemainders;
		x.Picard_ctrunc_normal_remainder(newRemainders, ode, tm_setting.step_exp_table[1], tm_setting.order, intermediate_ranges, g_setting);

		// add the uncertainties and the cutoff intervals onto the result
		for(unsigned int i=0; i<rangeDim; ++i)
		{
			newRemainders[i] += intDifferences[i];

			if(newRemainders[i].subseteq(x.tms[i].remainder))
			{
				if(x.tms[i].remainder.widthRatio(newRemainders[i]) <= STOP_RATIO)
				{
					bfinished = false;
				}

				x.tms[i].remainder = newRemainders[i];
			}
			else
			{
				bfinished = true;
				break;
			}
		}
	}

	result.tmvPre = x;
	result.domain = domain;
	result.domain[0] = tm_setting.step_exp_table[1];

	return 1;
}













LinearFlowpipe::LinearFlowpipe()
{
}

LinearFlowpipe::LinearFlowpipe(const LinearFlowpipe & flowpipe)
{
	Phi		= flowpipe.Phi;
	Psi		= flowpipe.Psi;

	tv_remainder	= flowpipe.tv_remainder;
}

LinearFlowpipe::~LinearFlowpipe()
{
}

int LinearFlowpipe::safetyChecking(const std::vector<Constraint> & unsafeSet, const Taylor_Model_Computation_Setting & tm_setting, const Global_Computation_Setting & g_setting,
		const TaylorModelVec<Real> & tmv_of_X0, const std::vector<Interval> & polyRangeX0, const std::vector<Interval> & range_of_X0, const std::vector<Interval> & domain)
{
	if(unsafeSet.size() == 0)
	{
		return SAFE;
	}

	unsigned int rangeDim = Phi.rows();
	int result = UNKNOWN;
	bool bContained = true;

	Matrix<Interval> range_of_Phi(rangeDim, rangeDim);
	Phi.evaluate(range_of_Phi, interval_utm_setting.val_exp_table);

	Matrix<Interval> range_of_Psi(rangeDim, 1);
	Psi.evaluate(range_of_Psi, interval_utm_setting.val_exp_table);

	std::vector<Interval> range_of_x = range_of_Phi * range_of_X0;

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		range_of_x[i] += range_of_Psi[i][0];
	}

	for(unsigned int i=0; i<unsafeSet.size(); ++i)
	{
		Interval I;

		// interval evaluation on the constraint
		unsafeSet[i].expression.evaluate(I, range_of_x);

		if(unsafeSet[i].bound < I.inf())
		{
			// no intersection with the unsafe set
			result = SAFE;
			break;
		}
		else
		{
			if(!(unsafeSet[i].bound >= I.sup()) && bContained)
			{
				bContained = false;
			}
		}
	}

	if(result == UNKNOWN)
	{
		if(bContained)
		{
			return UNSAFE;
		}
		else
		{
			// do a simple branch & bound for safety checking
			TaylorModelVec<Real> tmvFlowpipe;
			evaluate(tmvFlowpipe, tmv_of_X0, polyRangeX0, range_of_X0, domain, tm_setting);

			std::vector<HornerForm<Real> > obj_hfs;
			std::vector<Interval> obj_rems;

			result = SAFE;

			for(unsigned int i=0; i<unsafeSet.size(); ++i)
			{
				TaylorModel<Real> tmTmp;

				// interval evaluation on the constraint
				unsafeSet[i].expression.evaluate(tmTmp, tmvFlowpipe.tms, tm_setting.order, domain, tm_setting.cutoff_threshold, g_setting);

				HornerForm<Real> obj_hf;
				tmTmp.expansion.toHornerForm(obj_hf);
				obj_hfs.push_back(obj_hf);
				obj_rems.push_back(tmTmp.remainder);
			}

			std::vector<Interval> refined_domain = domain;

			std::list<Interval> subdivisions;

			if(domain[0].width() > REFINEMENT_PREC)
			{
				subdivisions.push_back(domain[0]);
			}

			for(; subdivisions.size() > 0; )
			{
				Interval subdivision = subdivisions.front();
				subdivisions.pop_front();

				int result_iter = UNKNOWN;
				bool bContained_iter = true;

				refined_domain[0] = subdivision;

				for(int i=0; i<unsafeSet.size(); ++i)
				{
					Interval I;
					obj_hfs[i].evaluate(I, refined_domain);

					I += obj_rems[i];

					if(unsafeSet[i].bound < I.inf())
					{
						// no intersection with the unsafe set
						result_iter = SAFE;
						break;
					}
					else
					{
						if(!(unsafeSet[i].bound >= I.sup()) && bContained_iter)
						{
							bContained_iter = false;
						}
					}
				}

				if(result_iter == UNKNOWN)
				{
					if(bContained_iter)
					{
						return UNSAFE;
					}
					else
					{
						if(subdivision.width() <= REFINEMENT_PREC)
						{
							return UNKNOWN;
						}

						// split the domain
						Interval I1, I2;
						subdivision.split(I1, I2);

						if(I1.width() <= REFINEMENT_PREC)
						{
							if(result == SAFE)
								result = UNKNOWN;
						}
						else
						{
							subdivisions.push_back(I1);
						}

						if(I2.width() <= REFINEMENT_PREC)
						{
							if(result == SAFE)
								result = UNKNOWN;
						}
						else
						{
							subdivisions.push_back(I2);
						}
					}
				}
			}

			return result;
		}
	}
	else
	{
		return SAFE;
	}
}

void LinearFlowpipe::evaluate(TaylorModelVec<Real> & result, const TaylorModelVec<Real> & tmv_of_X0, const std::vector<Interval> & polyRangeX0, const std::vector<Interval> & range_of_X0, const std::vector<Interval> & domain, const Taylor_Model_Computation_Setting & tm_setting)
{
	unsigned int rangeDim = Phi.rows();
	unsigned int domainDim = domain.size();

	result.clear();

	TaylorModelVec<Real> tmvTmp;

	for(int i=0; i<rangeDim; ++i)
	{
		TaylorModel<Real> tmTmp1;

		for(int j=0; j<rangeDim; ++j)
		{
			TaylorModel<Real> tmTmp2(Phi[i][j], domainDim, true);
			tmTmp2.mul_assign(j+1, range_of_X0[j]);

			tmTmp1 += tmTmp2;
		}

		tmvTmp.tms.push_back(tmTmp1);
	}

	tmvTmp.insert_ctrunc(result, tmv_of_X0, polyRangeX0, domain, tm_setting.order, tm_setting.cutoff_threshold);


	if(Psi.rows() > 0)
	{
		for(int i=0; i<rangeDim; ++i)
		{
			TaylorModel<Real> tmTmp(Psi[i][0], domainDim, true);
			result.tms[i] += tmTmp;
		}
	}

	if(!tv_remainder.isEmpty())
	{
		Matrix<Interval> im_tv_remainder(rangeDim, 1);
		tv_remainder.intEval(im_tv_remainder);

		for(int i=0; i<rangeDim; ++i)
		{
			result.tms[i].remainder += im_tv_remainder[i][0];
		}
	}
}

void LinearFlowpipe::evaluate(TaylorModelVec<Real> & result, const std::vector<unsigned int> & outputAxes, const TaylorModelVec<Real> & tmv_of_X0, const std::vector<Interval> & polyRangeX0, const std::vector<Interval> & range_of_X0, const std::vector<Interval> & domain, const Taylor_Model_Computation_Setting & tm_setting)
{
	unsigned int rangeDim = Phi.rows();
	unsigned int domainDim = domain.size();

	result.clear();

	TaylorModelVec<Real> tmvTmp;

	for(int i=0; i<outputAxes.size(); ++i)
	{
		TaylorModel<Real> tmTmp1;

		for(int j=0; j<rangeDim; ++j)
		{
			TaylorModel<Real> tmTmp2(Phi[outputAxes[i]][j], domainDim, true);
			tmTmp2.mul_assign(j+1, range_of_X0[j]);

			tmTmp1 += tmTmp2;
		}

		tmvTmp.tms.push_back(tmTmp1);
	}

	tmvTmp.insert_ctrunc(result, tmv_of_X0, polyRangeX0, domain, tm_setting.order, tm_setting.cutoff_threshold);


	if(Psi.rows() > 0)
	{
		for(int i=0; i<outputAxes.size(); ++i)
		{
			TaylorModel<Real> tmTmp(Psi[outputAxes[i]][0], domainDim, true);
			result.tms[i] += tmTmp;
		}
	}

	if(!tv_remainder.isEmpty())
	{
		Matrix<Interval> im_tv_remainder(rangeDim, 1);
		tv_remainder.intEval(im_tv_remainder);

		for(int i=0; i<outputAxes.size(); ++i)
		{
			result.tms[i].remainder += im_tv_remainder[outputAxes[i]][0];
		}
	}
}


LinearFlowpipe & LinearFlowpipe::operator = (const LinearFlowpipe & flowpipe)
{
	if(this == &flowpipe)
		return *this;

	Phi				= flowpipe.Phi;
	Psi				= flowpipe.Psi;
	tv_remainder	= flowpipe.tv_remainder;

	return *this;
}










Result_of_Reachability::Result_of_Reachability()
{
	status = -1;
	num_of_flowpipes = 0;
}

Result_of_Reachability::Result_of_Reachability(const Result_of_Reachability & result)
{
	status				= result.status;
	num_of_flowpipes	= result.num_of_flowpipes;
	fp_end_of_time		= result.fp_end_of_time;
	nonlinear_flowpipes	= result.nonlinear_flowpipes;
	orders_of_flowpipes	= result.orders_of_flowpipes;
	safety_of_flowpipes	= result.safety_of_flowpipes;
}

Result_of_Reachability::~Result_of_Reachability()
{
}

void Result_of_Reachability::clear()
{
	status = -1;
	num_of_flowpipes = 0;

	nonlinear_flowpipes.clear();
	tmv_flowpipes.clear();
	orders_of_flowpipes.clear();
	safety_of_flowpipes.clear();
}

void Result_of_Reachability::transformToTaylorModels(const Taylor_Model_Computation_Setting & tm_setting, const bool bPrint)
{
	unsigned int prog = 0, total_size = nonlinear_flowpipes.size();

	if(bPrint)
	{
		printf("Translating the flowpipes...\n");
	}

	std::list<Flowpipe>::const_iterator fpIter = nonlinear_flowpipes.begin();
	std::list<unsigned int>::const_iterator orderIter = orders_of_flowpipes.begin();

	for(; fpIter != nonlinear_flowpipes.end(); ++fpIter, ++orderIter)
	{
		TaylorModelVec<Real> tmvTmp;

		fpIter->compose(tmvTmp, *orderIter, tm_setting.cutoff_threshold);

		tmv_flowpipes.push_back(tmvTmp);

		if(bPrint)
		{
			++prog;
			printf("\b\b\b");
			printf(BOLD_FONT "%%" RESET_COLOR);
			printf(BOLD_FONT "%2d" RESET_COLOR, (int)(prog*100/total_size));
			fflush(stdout);
		}
	}

	if(bPrint)
	{
		printf("\nDone.\n");
	}
}

void Result_of_Reachability::transformToTaylorModels(const Computational_Setting & c_setting)
{
	transformToTaylorModels(c_setting.tm_setting, c_setting.bPrint);
}

Result_of_Reachability & Result_of_Reachability::operator = (const Result_of_Reachability & result)
{
	if(this == &result)
		return *this;

	status				= result.status;
	num_of_flowpipes	= result.num_of_flowpipes;
	fp_end_of_time		= result.fp_end_of_time;
	nonlinear_flowpipes	= result.nonlinear_flowpipes;
	orders_of_flowpipes	= result.orders_of_flowpipes;
	safety_of_flowpipes	= result.safety_of_flowpipes;

	return *this;
}











Linear_Time_Invariant_Dynamics::Linear_Time_Invariant_Dynamics(const Matrix<Real> & A, const Matrix<UnivariateTaylorModel<Real> > & B)
{
	rm_dyn_A = A;
	utm_dyn_B = B;

	if(B.isZero())
	{
		bAuto = true;
	}
	else
	{
		bAuto = false;
	}

	unsigned int n = A.rows();
	Matrix<bool> conMatrix(n, n), adjMatrix(n, n);

	for(unsigned int i=0; i<n; ++i)
	{
		for(unsigned int j=0; j<n; ++j)
		{
			if(rm_dyn_A[i][j] != 0)
			{
				adjMatrix[i][j] = true;
			}
		}
	}

	check_connectivities(conMatrix, adjMatrix);
	connectivity = conMatrix;
}

Linear_Time_Invariant_Dynamics::Linear_Time_Invariant_Dynamics(const Linear_Time_Invariant_Dynamics & dynamics)
{
	rm_dyn_A		= dynamics.rm_dyn_A;
	utm_dyn_B		= dynamics.utm_dyn_B;
	bAuto			= dynamics.bAuto;
	connectivity	= dynamics.connectivity;
}

Linear_Time_Invariant_Dynamics::~Linear_Time_Invariant_Dynamics()
{
}

Linear_Time_Invariant_Dynamics & Linear_Time_Invariant_Dynamics::operator = (const Linear_Time_Invariant_Dynamics & dynamics)
{
	if(this == &dynamics)
		return *this;

	rm_dyn_A		= dynamics.rm_dyn_A;
	utm_dyn_B		= dynamics.utm_dyn_B;
	bAuto			= dynamics.bAuto;
	connectivity	= dynamics.connectivity;

	return *this;
}

int Linear_Time_Invariant_Dynamics::reach_LTI(std::list<LinearFlowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, const Taylor_Model_Computation_Setting & tm_setting,
		const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput)
{
	// find a proper parameter r = 2^n such that |A*delta/r| < 0.1
	Real A_max = rm_dyn_A.max_norm();
	Real threshold = 0.1;
	double step = tm_setting.step_max;
	Real rStep = tm_setting.step_max;

	unsigned int r = 1, n = 0;
	while(A_max >= threshold)
	{
		r *= 2;
		++n;
		A_max /= r;
	}

	// find the proper order
	Real A_min = rm_dyn_A.min_entry();
	A_min /= r;
	Real tolerance = APPROX_TOLERANCE;
	Real error;

	unsigned int approx_order = findProperOrder(error, A_max, A_min, tolerance, tm_setting.order);


	unsigned int nec_order = 2*approx_order + 1;
	std::vector<Real> step_end_exp_table(nec_order + 1, 1);
	step_end_exp_table[1] = step;

	for(unsigned int i=2; i<=nec_order; ++i)
	{
		step_end_exp_table[i] = step_end_exp_table[i-1] * rStep;
	}


	Interval intStep(0, step);
	interval_utm_setting.setValue(intStep, approx_order);

	int rangeDim = rm_dyn_A.rows();

	// identity matrix
	Matrix<Real> identity(rangeDim);

	Matrix<Real> A_scaled = rm_dyn_A / r;

	// compute the Taylor series to the order of approx_order
	std::vector<Matrix<Real> > A_exp_table;
	compute_mat_pow(A_exp_table, A_scaled, approx_order);

	Matrix<UnivariatePolynomial<Real> > expansion_exp_A_t_k = (Matrix<UnivariatePolynomial<Real> >)identity;

	Real tmp = 1;

	for(unsigned int i=1; i<=approx_order; ++i)
	{
		Matrix<UnivariatePolynomial<Real> > A_t_i = (Matrix<UnivariatePolynomial<Real> >)A_exp_table[i];
		A_t_i.times_x(i);

		tmp *= i;
		A_t_i /= tmp;

		expansion_exp_A_t_k += A_t_i;
	}

	Matrix<UnivariateTaylorModel<Real> > utm_Phi_0 = expansion_exp_A_t_k;


	// compute a proper remainder
	Interval intErr;
	error.to_sym_int(intErr);

	Matrix<Interval> im_error(rangeDim, rangeDim);
	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			if(connectivity[i][j])
			{
				im_error[i][j] = intErr;
			}
		}
	}

	utm_Phi_0.setRemainder(im_error);

	for(unsigned int i=0; i<n; ++i)
	{
		utm_Phi_0 *= utm_Phi_0;
	}

	// evaluate the one-step mapping matrix
	Matrix<Real> rm_Phi(rangeDim, rangeDim);
	utm_Phi_0.evaluate(rm_Phi, step_end_exp_table);


	// compute the linear mapping for the constant part
	Matrix<UnivariateTaylorModel<Real> > utm_Psi_0;
	Matrix<UnivariateTaylorModel<Real> > utm_Psi(rangeDim, 1);

	LinearFlowpipe flowpipe;

	if(!bAuto)
	{
		utm_Psi_0 = utm_Phi_0 * utm_dyn_B;
		utm_Psi_0.integral(intStep);
		utm_Psi_0.evaluate(utm_Psi, step_end_exp_table);
		utm_Psi_0.ctrunc(tm_setting.order);
		flowpipe.Psi = utm_Psi_0;
	}

	utm_Phi_0.ctrunc(tm_setting.order);
	flowpipe.Phi = utm_Phi_0;

	interval_utm_setting.setOrder(tm_setting.order);


	flowpipes.clear();
	flowpipes_safety.clear();

	num_of_flowpipes = 1;


	// perform the safety checking on the first flowpipe
	int checking_result = COMPLETED_SAFE;
	std::vector<std::vector<Interval> > polyRangeX0;
	std::vector<std::vector<Interval> > range_of_X0;

	std::vector<Interval> domain = initialSets[0].domain;
	domain[0] = intStep;


	if(bSafetyChecking)
	{
		for(int m=0; m<initialSets.size(); ++m)
		{
			std::vector<Interval> intVecTemp;
			initialSets[m].tmvPre.polyRangeNormal(intVecTemp, interval_utm_setting.val_exp_table);
			polyRangeX0.push_back(intVecTemp);

			std::vector<Interval> rangeX0(rangeDim);
			for(int k=0; k<rangeDim; ++k)
			{
				rangeX0[k] = intVecTemp[k] + initialSets[m].tmvPre.tms[k].remainder;
			}

			range_of_X0.push_back(rangeX0);
		}

		if(bTMOutput || bPlot)
		{
			flowpipes.push_back(flowpipe);
		}

		for(int m=0; m<initialSets.size(); ++m)
		{
			int safety;

			safety = flowpipe.safetyChecking(unsafeSet, tm_setting, g_setting, initialSets[m].tmvPre, polyRangeX0[m], range_of_X0[m], domain);

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
	}
	else
	{
		if(bTMOutput || bPlot)
		{
			flowpipes.push_back(flowpipe);

			for(int m=0; m<initialSets.size(); ++m)
			{
				flowpipes_safety.push_back(SAFE);
			}
		}
	}


	Matrix<UnivariateTaylorModel<Real> > utm_global_Psi = utm_Psi;
	Matrix<Real> rm_global_Phi(rangeDim);

	int N = (int)ceil(time/step);

	if(bPrint)
	{
		printf("time = %f,\t", step);
		printf("step = %f,\t", step);
		printf("order = %d\n", tm_setting.order);
	}


	for(int i=1; i<N; ++i)
	{
		LinearFlowpipe newFlowpipe;

//		rm_global_Phi = rm_Phi * rm_global_Phi;
		newFlowpipe.Phi = rm_global_Phi * utm_Phi_0;
		newFlowpipe.Phi.evaluate(rm_global_Phi, step_end_exp_table);

		if(!bAuto)
		{
			newFlowpipe.Psi = utm_Phi_0 * utm_global_Psi + utm_Psi_0;
			utm_global_Psi += rm_global_Phi * utm_Psi;
		}

		++num_of_flowpipes;

		if(bSafetyChecking)
		{
			if(bTMOutput || bPlot)
			{
				flowpipes.push_back(newFlowpipe);
			}

			for(int m=0; m<initialSets.size(); ++m)
			{
				int safety;

				safety = newFlowpipe.safetyChecking(unsafeSet, tm_setting, g_setting, initialSets[m].tmvPre, polyRangeX0[m], range_of_X0[m], domain);

				if(bTMOutput || bPlot)
				{
					flowpipes_safety.push_back(safety);
				}

				if(safety == UNSAFE)
				{
					if(true)
					{
						printf("time = %f,\t", (i+1)*step);
						printf("step = %f,\t", step);
						printf("order = %d\n", tm_setting.order);
					}

					return COMPLETED_UNSAFE;
				}
				else if(safety == UNKNOWN && checking_result == COMPLETED_SAFE)
				{
					checking_result = COMPLETED_UNKNOWN;
				}
			}
		}
		else
		{
			if(bTMOutput || bPlot)
			{
				flowpipes.push_back(newFlowpipe);

				for(int m=0; m<initialSets.size(); ++m)
				{
					flowpipes_safety.push_back(SAFE);
				}
			}
		}

		if(bPrint)
		{
			printf("time = %f,\t", (i+1)*step);
			printf("step = %f,\t", step);
			printf("order = %d\n", tm_setting.order);
		}
	}

	return checking_result;
}

int Linear_Time_Invariant_Dynamics::reach_LTV(std::list<LinearFlowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, const Taylor_Model_Computation_Setting & tm_setting,
		const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput)
{
	return 0;
}

int Linear_Time_Invariant_Dynamics::reach(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, const Taylor_Model_Computation_Setting & tm_setting,
		const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
	return 0;
}

int Linear_Time_Invariant_Dynamics::reach_adaptive_stepsize(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, Taylor_Model_Computation_Setting & tm_setting,
		const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
	return 0;
}

int Linear_Time_Invariant_Dynamics::reach_adaptive_order(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, Taylor_Model_Computation_Setting & tm_setting,
		const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
	return 0;
}

int Linear_Time_Invariant_Dynamics::reach_symbolic_remainder(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, const Taylor_Model_Computation_Setting & tm_setting,
		const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
	return 0;
}

int Linear_Time_Invariant_Dynamics::reach_symbolic_remainder_adaptive_stepsize(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, Taylor_Model_Computation_Setting & tm_setting,
		const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
	return 0;
}

int Linear_Time_Invariant_Dynamics::reach_symbolic_remainder_adaptive_order(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, Taylor_Model_Computation_Setting & tm_setting,
		const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
	return 0;
}










Linear_Time_Varying_Dynamics::Linear_Time_Varying_Dynamics(const Matrix<UnivariatePolynomial<Real> > & A, const Matrix<UnivariatePolynomial<Real> > & B, const Matrix<UnivariatePolynomial<Real> > & C)
{
	upm_dyn_A			= A;
	upm_dyn_B			= B;
	upm_dyn_tv			= C;

	if(B.isZero())
	{
		bAuto = true;
	}
	else
	{
		bAuto = false;
	}

	Interval intUnit(-1,1);
	Matrix<Interval> tmp(tvPars.size(), 1, intUnit);
	uncertain_range	= tmp;

	unsigned int n = A.rows();
	Matrix<bool> conMatrix(n, n), adjMatrix(n, n);

	for(unsigned int i=0; i<n; ++i)
	{
		for(unsigned int j=0; j<n; ++j)
		{
			if(!(upm_dyn_A[i][j].isZero()))
			{
				adjMatrix[i][j] = true;
			}
		}
	}

	check_connectivities(conMatrix, adjMatrix);
	connectivity = conMatrix;
}

Linear_Time_Varying_Dynamics::Linear_Time_Varying_Dynamics(const Linear_Time_Varying_Dynamics & dynamics)
{
	upm_dyn_A			= dynamics.upm_dyn_A;
	upm_dyn_B			= dynamics.upm_dyn_B;
	upm_dyn_tv			= dynamics.upm_dyn_tv;
	uncertain_range		= dynamics.uncertain_range;
	bAuto				= dynamics.bAuto;
}

Linear_Time_Varying_Dynamics::~Linear_Time_Varying_Dynamics()
{
}

Linear_Time_Varying_Dynamics & Linear_Time_Varying_Dynamics::operator = (const Linear_Time_Varying_Dynamics & dynamics)
{
	if(this == &dynamics)
		return *this;

	upm_dyn_A			= dynamics.upm_dyn_A;
	upm_dyn_B			= dynamics.upm_dyn_B;
	upm_dyn_tv			= dynamics.upm_dyn_tv;
	uncertain_range		= dynamics.uncertain_range;
	bAuto				= dynamics.bAuto;

	return *this;
}

int Linear_Time_Varying_Dynamics::reach_LTI(std::list<LinearFlowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, const Taylor_Model_Computation_Setting & tm_setting,
		const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput)
{
	return 0;
}

int Linear_Time_Varying_Dynamics::reach_LTV(std::list<LinearFlowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, const Taylor_Model_Computation_Setting & tm_setting,
		const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput)
{
	double step = tm_setting.step_max;
	Real rStep = tm_setting.step_max;
	const int rangeDim = upm_dyn_A.rows();
	Interval intStep(0, step);

	unsigned int maxOrder = upm_dyn_A.degree();
	unsigned int maxOrder_B = upm_dyn_B.degree();

	if(maxOrder < maxOrder_B)
	{
		maxOrder = maxOrder_B;
	}

	unsigned int maxOrder_tv = upm_dyn_tv.degree();

	if(maxOrder < maxOrder_tv)
	{
		maxOrder = maxOrder_tv;
	}

	unsigned int nec_order = 2 * (tm_setting.order + 1 + maxOrder) + 1;

	interval_utm_setting.setValue(intStep, nec_order);

	std::vector<Real> step_end_exp_table(nec_order + 1, 1);
	step_end_exp_table[1] = step;

	for(unsigned int i=2; i<=nec_order; ++i)
	{
		step_end_exp_table[i] = step_end_exp_table[i-1] * rStep;
	}


	unsigned int numTVPar = upm_dyn_tv.cols();

	Matrix<Real> Phi_t_0(rangeDim), Psi_t_0(rangeDim, 1);
	Zonotope global_tv_remainder(rangeDim);

	flowpipes.clear();
	flowpipes_safety.clear();

	int checking_result = COMPLETED_SAFE;
	std::vector<std::vector<Interval> > polyRangeX0;
	std::vector<std::vector<Interval> > range_of_X0;


	if(bSafetyChecking)
	{
		for(int m=0; m<initialSets.size(); ++m)
		{
			std::vector<Interval> intVecTemp;
			initialSets[m].tmvPre.polyRangeNormal(intVecTemp, interval_utm_setting.val_exp_table);
			polyRangeX0.push_back(intVecTemp);

			std::vector<Interval> range(rangeDim);
			for(int i1=0; i1<rangeDim; ++i1)
			{
				range[i1] = intVecTemp[i1] + initialSets[m].tmvPre.tms[i1].remainder;
			}

			range_of_X0.push_back(range);
		}
	}

	std::vector<Interval> domain = initialSets[0].domain;
	domain[0] = intStep;

	unsigned int num = 0;

	num_of_flowpipes = 0;

	Matrix<Interval> tv_part(rangeDim, numTVPar);

	for(double t0 = 0; t0 < time - THRESHOLD_HIGH; )
	{
		std::vector<Real> t0_coefficients;
		t0_coefficients.push_back(t0);
		t0_coefficients.push_back(1);
		UnivariatePolynomial<Real> up_t0(t0_coefficients);

		LinearFlowpipe flowpipe;
		Matrix<Real> Phi_step_end(rangeDim, rangeDim);
		Matrix<Real> Psi_step_end(rangeDim, 1);

		compute_one_step_trans(flowpipe.Phi, Phi_step_end, flowpipe.Psi, Psi_step_end,
				tv_part, upm_dyn_A, upm_dyn_B, upm_dyn_tv, connectivity, bAuto, up_t0,
				tm_setting.order, step_end_exp_table);


		if(numTVPar > 0)
		{
			if(tm_setting.queue_size >= 0 && num > tm_setting.queue_size)
			{
				num = 0;
				global_tv_remainder.simplify();
			}

			Matrix<Interval> im_temp = tv_part * uncertain_range;
			im_temp *= intStep;
			Zonotope zonoTmp(im_temp);
			flowpipe.tv_remainder = global_tv_remainder + zonoTmp;

			++num;
		}

		if(!bAuto)
		{
			Matrix<UnivariateTaylorModel<Real> > utm_tmp = flowpipe.Phi * Psi_t_0;
			flowpipe.Psi = utm_tmp + flowpipe.Psi;
			utm_tmp.evaluate(Psi_t_0, step_end_exp_table);
			Psi_t_0 += Psi_step_end;
		}

		if(numTVPar > 0)
		{
			global_tv_remainder = Phi_step_end * flowpipe.tv_remainder;
		}

		flowpipe.Phi *= Phi_t_0;
		Phi_t_0 = Phi_step_end * Phi_t_0;

		++num_of_flowpipes;

		if(bSafetyChecking)
		{
			if(bTMOutput || bPlot)
			{
				flowpipes.push_back(flowpipe);
			}

			for(int m=0; m<initialSets.size(); ++m)
			{
				int safety;

				std::vector<Matrix<Real> > constraints;
				std::vector<Matrix<Interval> > precond_Phi, precond_Psi;

				safety = flowpipe.safetyChecking(unsafeSet, tm_setting, g_setting, initialSets[m].tmvPre, polyRangeX0[m], range_of_X0[m], domain);

				if(bTMOutput || bPlot)
				{
					flowpipes_safety.push_back(safety);
				}

				if(safety == UNSAFE)
				{
					if(bPrint)
					{
						printf("time = %f,\t", t0 + step);
						printf("step = %f,\t", step);
						printf("order = %d\n", tm_setting.order);
					}

					return COMPLETED_UNSAFE;
				}
				else if(safety == UNKNOWN && checking_result == COMPLETED_SAFE)
				{
					checking_result = COMPLETED_UNKNOWN;
				}
			}
		}
		else
		{
			if(bTMOutput || bPlot)
			{
				flowpipes.push_back(flowpipe);
				flowpipes_safety.push_back(SAFE);
			}
		}

		t0 += step;

		if(bPrint)
		{
			printf("time = %f,\t", t0);
			printf("step = %f,\t", step);
			printf("order = %d\n", tm_setting.order);
		}
	}

	return checking_result;
}

int Linear_Time_Varying_Dynamics::reach(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, const Taylor_Model_Computation_Setting & tm_setting,
		const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
	return 0;
}

int Linear_Time_Varying_Dynamics::reach_adaptive_stepsize(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, Taylor_Model_Computation_Setting & tm_setting,
		const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
	return 0;
}

int Linear_Time_Varying_Dynamics::reach_adaptive_order(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, Taylor_Model_Computation_Setting & tm_setting,
		const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
	return 0;
}

int Linear_Time_Varying_Dynamics::reach_symbolic_remainder(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, const Taylor_Model_Computation_Setting & tm_setting,
		const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
	return 0;
}

int Linear_Time_Varying_Dynamics::reach_symbolic_remainder_adaptive_stepsize(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, Taylor_Model_Computation_Setting & tm_setting,
		const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
	return 0;
}

int Linear_Time_Varying_Dynamics::reach_symbolic_remainder_adaptive_order(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, Taylor_Model_Computation_Setting & tm_setting,
		const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
	return 0;
}






Deterministic_Continuous_Dynamics::Deterministic_Continuous_Dynamics(const std::vector<Expression_AST<Real> > & dynamics)
{
	for(unsigned int i=0; i<dynamics.size(); ++i)
	{
		expressions.push_back(dynamics[i]);
	}
}

Deterministic_Continuous_Dynamics::Deterministic_Continuous_Dynamics(const Deterministic_Continuous_Dynamics & dynamics)
{
	expressions		= dynamics.expressions;
}

Deterministic_Continuous_Dynamics::~Deterministic_Continuous_Dynamics()
{
}

Deterministic_Continuous_Dynamics & Deterministic_Continuous_Dynamics::operator = (const Deterministic_Continuous_Dynamics & dynamics)
{
	if(this == &dynamics)
		return *this;

	expressions		= dynamics.expressions;

	return *this;
}

int Deterministic_Continuous_Dynamics::reach_LTI(std::list<LinearFlowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, const Taylor_Model_Computation_Setting & tm_setting,
		const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput)
{
	return 0;
}

int Deterministic_Continuous_Dynamics::reach_LTV(std::list<LinearFlowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, const Taylor_Model_Computation_Setting & tm_setting,
		const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput)
{
	return 0;
}

int Deterministic_Continuous_Dynamics::reach(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, const Taylor_Model_Computation_Setting & tm_setting,
		const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
//	flowpipes.clear();
//	flowpipe_orders.clear();
//	flowpipes_safety.clear();

//	num_of_flowpipes = 0;
	std::vector<Constraint> dummy_invariant;

	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	for(int m=0; m<initialSets.size(); ++m)
	{
		Flowpipe newFlowpipe, currentFlowpipe = initialSets[m];

		for(double t=THRESHOLD_HIGH; t < time;)
		{
			int res = currentFlowpipe.advance_deterministic(newFlowpipe, expressions, tm_setting, dummy_invariant, g_setting);

			if(res == 1)
			{
				++num_of_flowpipes;
				flowpipe_orders.push_back(tm_setting.order);

				if(bSafetyChecking)
				{
					int safety = newFlowpipe.safetyChecking(unsafeSet, tm_setting, g_setting);

					if(bTMOutput || bPlot)
					{
						flowpipes.push_back(newFlowpipe);
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
	}

	return checking_result;
}

int Deterministic_Continuous_Dynamics::reach_adaptive_stepsize(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, Taylor_Model_Computation_Setting & tm_setting,
		const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
//	flowpipes.clear();
//	flowpipe_orders.clear();
//	flowpipes_safety.clear();

//	num_of_flowpipes = 0;
	std::vector<Constraint> dummy_invariant;

	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	for(int m=0; m<initialSets.size(); ++m)
	{
		Flowpipe newFlowpipe, currentFlowpipe = initialSets[m];

		for(double t=THRESHOLD_HIGH; t < time;)
		{
			int res = currentFlowpipe.advance_deterministic_adaptive_stepsize(newFlowpipe, expressions, tm_setting, dummy_invariant, g_setting);

			if(res == 1)
			{
				++num_of_flowpipes;
				flowpipe_orders.push_back(tm_setting.order);

				if(bSafetyChecking)
				{
					int safety = newFlowpipe.safetyChecking(unsafeSet, tm_setting, g_setting);

					if(bTMOutput || bPlot)
					{
						flowpipes.push_back(newFlowpipe);
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

				currentFlowpipe = newFlowpipe;

				double current_stepsize = tm_setting.step_exp_table[1].sup();
				t += current_stepsize;

				if(bPrint)
				{
					printf("time = %f,\t", t);
					printf("step = %f,\t", current_stepsize);
					printf("order = %d\n", tm_setting.order);
				}

				double new_stepsize = current_stepsize * LAMBDA_UP;
				double last_step = time - t;

				if(new_stepsize > last_step)
				{
					new_stepsize = last_step + THRESHOLD_LOW;
				}

				if(new_stepsize <= tm_setting.step_max - THRESHOLD_LOW)
				{
					tm_setting.setStepsize(new_stepsize, tm_setting.order);
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
	}

	return checking_result;
}

int Deterministic_Continuous_Dynamics::reach_adaptive_order(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, Taylor_Model_Computation_Setting & tm_setting,
		const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
//	flowpipes.clear();
//	flowpipe_orders.clear();
//	flowpipes_safety.clear();

//	num_of_flowpipes = 0;
	std::vector<Constraint> dummy_invariant;

	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	for(int m=0; m<initialSets.size(); ++m)
	{
		Flowpipe newFlowpipe, currentFlowpipe = initialSets[m];

		for(double t=THRESHOLD_HIGH; t < time;)
		{
			int res = currentFlowpipe.advance_deterministic_adaptive_order(newFlowpipe, expressions, tm_setting, dummy_invariant, g_setting);

			if(res == 1)
			{
				++num_of_flowpipes;
				flowpipe_orders.push_back(tm_setting.order);

				if(bSafetyChecking)
				{
					int safety = newFlowpipe.safetyChecking(unsafeSet, tm_setting, g_setting);

					if(bTMOutput || bPlot)
					{
						flowpipes.push_back(newFlowpipe);
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
	}

	return checking_result;
}

int Deterministic_Continuous_Dynamics::reach_symbolic_remainder(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, const Taylor_Model_Computation_Setting & tm_setting,
		const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
//	flowpipes.clear();
//	flowpipe_orders.clear();
//	flowpipes_safety.clear();

//	num_of_flowpipes = 0;
	std::vector<Constraint> dummy_invariant;

	std::vector<Real> initial_scalars(expressions.size(), 1);

	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	for(int m=0; m<initialSets.size(); ++m)
	{
		Flowpipe newFlowpipe, currentFlowpipe = initialSets[m];

		Symbolic_Remainder symbolic_remainder(currentFlowpipe);

		for(double t=THRESHOLD_HIGH; t < time;)
		{
			int res = currentFlowpipe.advance_deterministic(newFlowpipe, expressions, tm_setting, dummy_invariant, g_setting, symbolic_remainder);

			if(res == 1)
			{
				++num_of_flowpipes;
				flowpipe_orders.push_back(tm_setting.order);

				if(bSafetyChecking)
				{
					int safety = newFlowpipe.safetyChecking(unsafeSet, tm_setting, g_setting);

					if(bTMOutput || bPlot)
					{
						flowpipes.push_back(newFlowpipe);
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

				currentFlowpipe = newFlowpipe;

				t += step;

				if(bPrint)
				{
					printf("time = %f,\t", t);
					printf("step = %f,\t", step);
					printf("order = %d\n", tm_setting.order);
				}

				if(symbolic_remainder.J.size() >= tm_setting.queue_size)
				{
					symbolic_remainder.reset(currentFlowpipe);
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
	}

	return checking_result;
}

int Deterministic_Continuous_Dynamics::reach_symbolic_remainder_adaptive_stepsize(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, Taylor_Model_Computation_Setting & tm_setting,
		const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
//	flowpipes.clear();
//	flowpipe_orders.clear();
//	flowpipes_safety.clear();

//	num_of_flowpipes = 0;
	std::vector<Constraint> dummy_invariant;

	std::vector<Real> initial_scalars(expressions.size(), 1);

	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	for(int m=0; m<initialSets.size(); ++m)
	{
		Flowpipe newFlowpipe, currentFlowpipe = initialSets[m];

		tm_setting.setStepsize(tm_setting.step_max, tm_setting.order);

		Symbolic_Remainder symbolic_remainder(currentFlowpipe);

		for(double t=THRESHOLD_HIGH; t < time;)
		{
			int res = currentFlowpipe.advance_deterministic_adaptive_stepsize(newFlowpipe, expressions, tm_setting, dummy_invariant, g_setting, symbolic_remainder);

			if(res == 1)
			{
				++num_of_flowpipes;
				flowpipe_orders.push_back(tm_setting.order);

				if(bSafetyChecking)
				{
					int safety = newFlowpipe.safetyChecking(unsafeSet, tm_setting, g_setting);

					if(bTMOutput || bPlot)
					{
						flowpipes.push_back(newFlowpipe);
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

				currentFlowpipe = newFlowpipe;

				if(symbolic_remainder.J.size() >= tm_setting.queue_size)
				{
					symbolic_remainder.reset(currentFlowpipe);
				}

				double current_stepsize = tm_setting.step_exp_table[1].sup();
				t += current_stepsize;

				if(bPrint)
				{
					printf("time = %f,\t", t);
					printf("step = %f,\t", current_stepsize);
					printf("order = %d\n", tm_setting.order);
				}

				double new_stepsize = current_stepsize * LAMBDA_UP;
				double last_step = time - t;

				if(new_stepsize > last_step)
				{
					new_stepsize = last_step + THRESHOLD_LOW;
				}

				if(new_stepsize <= tm_setting.step_max - THRESHOLD_LOW)
				{
					tm_setting.setStepsize(new_stepsize, tm_setting.order);
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
	}

	return checking_result;
}

int Deterministic_Continuous_Dynamics::reach_symbolic_remainder_adaptive_order(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, Taylor_Model_Computation_Setting & tm_setting,
		const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
//	flowpipes.clear();
//	flowpipe_orders.clear();
//	flowpipes_safety.clear();

//	num_of_flowpipes = 0;
	std::vector<Constraint> dummy_invariant;

	std::vector<Real> initial_scalars(expressions.size(), 1);

	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	for(int m=0; m<initialSets.size(); ++m)
	{
		Flowpipe newFlowpipe, currentFlowpipe = initialSets[m];

		Symbolic_Remainder symbolic_remainder(currentFlowpipe);

		for(double t=THRESHOLD_HIGH; t < time;)
		{
			int res = currentFlowpipe.advance_deterministic_adaptive_order(newFlowpipe, expressions, tm_setting, dummy_invariant, g_setting, symbolic_remainder);

			if(res == 1)
			{
				++num_of_flowpipes;
				flowpipe_orders.push_back(tm_setting.order);

				if(bSafetyChecking)
				{
					int safety = newFlowpipe.safetyChecking(unsafeSet, tm_setting, g_setting);

					if(bTMOutput || bPlot)
					{
						flowpipes.push_back(newFlowpipe);
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

				currentFlowpipe = newFlowpipe;

				t += step;

				if(bPrint)
				{
					printf("time = %f,\t", t);
					printf("step = %f,\t", step);
					printf("order = %d\n", tm_setting.order);
				}

				if(symbolic_remainder.J.size() >= tm_setting.queue_size)
				{
					symbolic_remainder.reset(currentFlowpipe);
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
	}

	return checking_result;
}

void Deterministic_Continuous_Dynamics::reach(Result_of_Reachability & result, Computational_Setting & setting, const std::vector<Flowpipe> & initialSets, const std::vector<Constraint> & unsafeSet) const
{
	bool bSafetyChecking = false;

	if(unsafeSet.size() > 0)
	{
		bSafetyChecking = true;
	}

	if(setting.tm_setting.queue_size > 0)
	{
		// symbolic remainder

		if(setting.tm_setting.step_min > 0)
		{
			// adaptive stepsizes
			result.status = reach_symbolic_remainder_adaptive_stepsize(result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes, result.num_of_flowpipes,
					setting.time, initialSets, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);
		}
		else if(setting.tm_setting.order_max > 0)
		{
			// adaptive orders
			result.status = reach_symbolic_remainder_adaptive_order(result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes, result.num_of_flowpipes,
					setting.time, initialSets, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);
		}
		else
		{
			// fixed stepsizes and orders
			result.status = reach_symbolic_remainder(result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes, result.num_of_flowpipes,
					setting.time, initialSets, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);
		}
	}
	else
	{
		if(setting.tm_setting.step_min > 0)
		{
			// adaptive stepsizes
			result.status = reach_adaptive_stepsize(result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes, result.num_of_flowpipes,
					setting.time, initialSets, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);
		}
		else if(setting.tm_setting.order_max > 0)
		{
			// adaptive orders
			result.status = reach_adaptive_order(result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes, result.num_of_flowpipes,
					setting.time, initialSets, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);
		}
		else
		{
			// fixed stepsizes and orders
			result.status = reach(result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes, result.num_of_flowpipes,
					setting.time, initialSets, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);
		}
	}
}

void Deterministic_Continuous_Dynamics::reach(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const std::vector<Constraint> & unsafeSet) const
{
	bool bSafetyChecking = false;

	if(unsafeSet.size() > 0)
	{
		bSafetyChecking = true;
	}

	std::vector<Flowpipe> initialSets;
	initialSets.push_back(initialSet);

	if(setting.tm_setting.queue_size > 0)
	{
		// symbolic remainder

		if(setting.tm_setting.step_min > 0)
		{
			// adaptive stepsizes
			result.status = reach_symbolic_remainder_adaptive_stepsize(result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes, result.num_of_flowpipes,
					setting.time, initialSets, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);
		}
		else if(setting.tm_setting.order_max > 0)
		{
			// adaptive orders
			result.status = reach_symbolic_remainder_adaptive_order(result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes, result.num_of_flowpipes,
					setting.time, initialSets, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);
		}
		else
		{
			// fixed stepsizes and orders
			result.status = reach_symbolic_remainder(result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes, result.num_of_flowpipes,
					setting.time, initialSets, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);
		}
	}
	else
	{
		if(setting.tm_setting.step_min > 0)
		{
			// adaptive stepsizes
			result.status = reach_adaptive_stepsize(result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes, result.num_of_flowpipes,
					setting.time, initialSets, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);
		}
		else if(setting.tm_setting.order_max > 0)
		{
			// adaptive orders
			result.status = reach_adaptive_order(result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes, result.num_of_flowpipes,
					setting.time, initialSets, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);
		}
		else
		{
			// fixed stepsizes and orders
			result.status = reach(result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes, result.num_of_flowpipes,
					setting.time, initialSets, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);
		}
	}

	Flowpipe fpTmp = result.nonlinear_flowpipes.back();
	result.fp_end_of_time = fpTmp;

	fpTmp.tmvPre.evaluate_time(result.fp_end_of_time.tmvPre, setting.tm_setting.step_end_exp_table);
}






Nondeterministic_Continuous_Dynamics::Nondeterministic_Continuous_Dynamics(const std::vector<Expression_AST<Interval> > & dynamics)
{
	for(unsigned int i=0; i<dynamics.size(); ++i)
	{
		expressions.push_back(dynamics[i]);
	}
}

Nondeterministic_Continuous_Dynamics::Nondeterministic_Continuous_Dynamics(const Nondeterministic_Continuous_Dynamics & dynamics)
{
	expressions		= dynamics.expressions;
}

Nondeterministic_Continuous_Dynamics::~Nondeterministic_Continuous_Dynamics()
{
}

Nondeterministic_Continuous_Dynamics & Nondeterministic_Continuous_Dynamics::operator = (const Nondeterministic_Continuous_Dynamics & dynamics)
{
	if(this == &dynamics)
		return *this;

	expressions		= dynamics.expressions;

	return *this;
}

int Nondeterministic_Continuous_Dynamics::reach_LTI(std::list<LinearFlowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, const Taylor_Model_Computation_Setting & tm_setting,
		const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput)
{
	return 0;
}

int Nondeterministic_Continuous_Dynamics::reach_LTV(std::list<LinearFlowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, const Taylor_Model_Computation_Setting & tm_setting,
		const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput)
{
	return 0;
}

int Nondeterministic_Continuous_Dynamics::reach(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, const Taylor_Model_Computation_Setting & tm_setting,
		const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
//	flowpipes.clear();
//	flowpipe_orders.clear();
//	flowpipes_safety.clear();

//	num_of_flowpipes = 0;
	std::vector<Constraint> dummy_invariant;

	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	for(int m=0; m<initialSets.size(); ++m)
	{
		Flowpipe newFlowpipe, currentFlowpipe = initialSets[m];

		for(double t=THRESHOLD_HIGH; t < time;)
		{
			int res = currentFlowpipe.advance_nondeterministic(newFlowpipe, expressions, tm_setting, dummy_invariant, g_setting);

			if(res == 1)
			{
				++num_of_flowpipes;
				flowpipe_orders.push_back(tm_setting.order);

				if(bSafetyChecking)
				{
					int safety = newFlowpipe.safetyChecking(unsafeSet, tm_setting, g_setting);

					if(bTMOutput || bPlot)
					{
						flowpipes.push_back(newFlowpipe);
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
	}

	return checking_result;
}

int Nondeterministic_Continuous_Dynamics::reach_adaptive_stepsize(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, Taylor_Model_Computation_Setting & tm_setting,
		const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
//	flowpipes.clear();
//	flowpipe_orders.clear();
//	flowpipes_safety.clear();

//	num_of_flowpipes = 0;
	std::vector<Constraint> dummy_invariant;

	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	for(int m=0; m<initialSets.size(); ++m)
	{
		Flowpipe newFlowpipe, currentFlowpipe = initialSets[m];

		for(double t=THRESHOLD_HIGH; t < time;)
		{
			int res = currentFlowpipe.advance_nondeterministic_adaptive_stepsize(newFlowpipe, expressions, tm_setting, dummy_invariant, g_setting);

			if(res == 1)
			{
				++num_of_flowpipes;
				flowpipe_orders.push_back(tm_setting.order);

				if(bSafetyChecking)
				{
					int safety = newFlowpipe.safetyChecking(unsafeSet, tm_setting, g_setting);

					if(bTMOutput || bPlot)
					{
						flowpipes.push_back(newFlowpipe);
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

				currentFlowpipe = newFlowpipe;

				double current_stepsize = tm_setting.step_exp_table[1].sup();
				t += current_stepsize;

				if(bPrint)
				{
					printf("time = %f,\t", t);
					printf("step = %f,\t", current_stepsize);
					printf("order = %d\n", tm_setting.order);
				}

				double new_stepsize = current_stepsize * LAMBDA_UP;
				double last_step = time - t;

				if(new_stepsize > last_step)
				{
					new_stepsize = last_step + THRESHOLD_LOW;
				}

				if(new_stepsize <= tm_setting.step_max - THRESHOLD_LOW)
				{
					tm_setting.setStepsize(new_stepsize, tm_setting.order);
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
	}

	return checking_result;
}

int Nondeterministic_Continuous_Dynamics::reach_adaptive_order(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, Taylor_Model_Computation_Setting & tm_setting,
		const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
//	flowpipes.clear();
//	flowpipe_orders.clear();
//	flowpipes_safety.clear();

//	num_of_flowpipes = 0;
	std::vector<Constraint> dummy_invariant;

	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	for(int m=0; m<initialSets.size(); ++m)
	{
		Flowpipe newFlowpipe, currentFlowpipe = initialSets[m];

		for(double t=THRESHOLD_HIGH; t < time;)
		{
			int res = currentFlowpipe.advance_nondeterministic_adaptive_order(newFlowpipe, expressions, tm_setting, dummy_invariant, g_setting);

			if(res == 1)
			{
				++num_of_flowpipes;
				flowpipe_orders.push_back(tm_setting.order);

				if(bSafetyChecking)
				{
					int safety = newFlowpipe.safetyChecking(unsafeSet, tm_setting, g_setting);

					if(bTMOutput || bPlot)
					{
						flowpipes.push_back(newFlowpipe);
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
	}

	return checking_result;
}

int Nondeterministic_Continuous_Dynamics::reach_symbolic_remainder(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, const Taylor_Model_Computation_Setting & tm_setting,
		const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
//	flowpipes.clear();
//	flowpipe_orders.clear();
//	flowpipes_safety.clear();

//	num_of_flowpipes = 0;
	std::vector<Constraint> dummy_invariant;

	std::vector<Real> initial_scalars(expressions.size(), 1);

	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	for(int m=0; m<initialSets.size(); ++m)
	{
		Flowpipe newFlowpipe, currentFlowpipe = initialSets[m];

		Symbolic_Remainder symbolic_remainder(currentFlowpipe);

		for(double t=THRESHOLD_HIGH; t < time;)
		{
			int res = currentFlowpipe.advance_nondeterministic(newFlowpipe, expressions, tm_setting, dummy_invariant, g_setting, symbolic_remainder);

			if(res == 1)
			{
				++num_of_flowpipes;
				flowpipe_orders.push_back(tm_setting.order);

				if(bSafetyChecking)
				{
					int safety = newFlowpipe.safetyChecking(unsafeSet, tm_setting, g_setting);

					if(bTMOutput || bPlot)
					{
						flowpipes.push_back(newFlowpipe);
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

				currentFlowpipe = newFlowpipe;

				t += step;

				if(bPrint)
				{
					printf("time = %f,\t", t);
					printf("step = %f,\t", step);
					printf("order = %d\n", tm_setting.order);
				}

				if(symbolic_remainder.J.size() >= tm_setting.queue_size)
				{
					symbolic_remainder.reset(currentFlowpipe);
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
	}

	return checking_result;
}

int Nondeterministic_Continuous_Dynamics::reach_symbolic_remainder_adaptive_stepsize(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, Taylor_Model_Computation_Setting & tm_setting,
		const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
//	flowpipes.clear();
//	flowpipe_orders.clear();
//	flowpipes_safety.clear();

//	num_of_flowpipes = 0;
	std::vector<Constraint> dummy_invariant;

	std::vector<Real> initial_scalars(expressions.size(), 1);

	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	for(int m=0; m<initialSets.size(); ++m)
	{
		Flowpipe newFlowpipe, currentFlowpipe = initialSets[m];

		Symbolic_Remainder symbolic_remainder(currentFlowpipe);

		for(double t=THRESHOLD_HIGH; t < time;)
		{
			int res = currentFlowpipe.advance_nondeterministic_adaptive_stepsize(newFlowpipe, expressions, tm_setting, dummy_invariant, g_setting, symbolic_remainder);

			if(res == 1)
			{
				++num_of_flowpipes;
				flowpipe_orders.push_back(tm_setting.order);

				if(bSafetyChecking)
				{
					int safety = newFlowpipe.safetyChecking(unsafeSet, tm_setting, g_setting);

					if(bTMOutput || bPlot)
					{
						flowpipes.push_back(newFlowpipe);
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

				currentFlowpipe = newFlowpipe;

				if(symbolic_remainder.J.size() >= tm_setting.queue_size)
				{
					symbolic_remainder.reset(currentFlowpipe);
				}

				double current_stepsize = tm_setting.step_exp_table[1].sup();
				t += current_stepsize;

				if(bPrint)
				{
					printf("time = %f,\t", t);
					printf("step = %f,\t", current_stepsize);
					printf("order = %d\n", tm_setting.order);
				}

				double new_stepsize = current_stepsize * LAMBDA_UP;
				double last_step = time - t;

				if(new_stepsize > last_step)
				{
					new_stepsize = last_step + THRESHOLD_LOW;
				}

				if(new_stepsize <= tm_setting.step_max - THRESHOLD_LOW)
				{
					tm_setting.setStepsize(new_stepsize, tm_setting.order);
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
	}

	return checking_result;
}

int Nondeterministic_Continuous_Dynamics::reach_symbolic_remainder_adaptive_order(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const std::vector<Flowpipe> & initialSets, Taylor_Model_Computation_Setting & tm_setting,
		const Global_Computation_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
//	flowpipes.clear();
//	flowpipe_orders.clear();
//	flowpipes_safety.clear();

//	num_of_flowpipes = 0;
	std::vector<Constraint> dummy_invariant;

	std::vector<Real> initial_scalars(expressions.size(), 1);

	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	for(int m=0; m<initialSets.size(); ++m)
	{
		Flowpipe newFlowpipe, currentFlowpipe = initialSets[m];

		Symbolic_Remainder symbolic_remainder(currentFlowpipe);

		for(double t=THRESHOLD_HIGH; t < time;)
		{
			int res = currentFlowpipe.advance_nondeterministic_adaptive_order(newFlowpipe, expressions, tm_setting, dummy_invariant, g_setting, symbolic_remainder);

			if(res == 1)
			{
				++num_of_flowpipes;
				flowpipe_orders.push_back(tm_setting.order);

				if(bSafetyChecking)
				{
					int safety = newFlowpipe.safetyChecking(unsafeSet, tm_setting, g_setting);

					if(bTMOutput || bPlot)
					{
						flowpipes.push_back(newFlowpipe);
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

				currentFlowpipe = newFlowpipe;

				t += step;

				if(bPrint)
				{
					printf("time = %f,\t", t);
					printf("step = %f,\t", step);
					printf("order = %d\n", tm_setting.order);
				}

				if(symbolic_remainder.J.size() >= tm_setting.queue_size)
				{
					symbolic_remainder.reset(currentFlowpipe);
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
	}

	return checking_result;
}

void Nondeterministic_Continuous_Dynamics::reach(Result_of_Reachability & result, Computational_Setting & setting, const std::vector<Flowpipe> & initialSets, const std::vector<Constraint> & unsafeSet) const
{
	bool bSafetyChecking = false;

	if(unsafeSet.size() > 0)
	{
		bSafetyChecking = true;
	}

	if(setting.tm_setting.queue_size > 0)
	{
		// symbolic remainder

		if(setting.tm_setting.step_min > 0)
		{
			// adaptive stepsizes
			result.status = reach_symbolic_remainder_adaptive_stepsize(result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes, result.num_of_flowpipes,
					setting.time, initialSets, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);
		}
		else if(setting.tm_setting.order_max > 0)
		{
			// adaptive orders
			result.status = reach_symbolic_remainder_adaptive_order(result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes, result.num_of_flowpipes,
					setting.time, initialSets, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);
		}
		else
		{
			// fixed stepsizes and orders
			result.status = reach_symbolic_remainder(result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes, result.num_of_flowpipes,
					setting.time, initialSets, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);
		}
	}
	else
	{
		if(setting.tm_setting.step_min > 0)
		{
			// adaptive stepsizes
			result.status = reach_adaptive_stepsize(result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes, result.num_of_flowpipes,
					setting.time, initialSets, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);
		}
		else if(setting.tm_setting.order_max > 0)
		{
			// adaptive orders
			result.status = reach_adaptive_order(result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes, result.num_of_flowpipes,
					setting.time, initialSets, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);
		}
		else
		{
			// fixed stepsizes and orders
			result.status = reach(result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes, result.num_of_flowpipes,
					setting.time, initialSets, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);
		}
	}
}

void Nondeterministic_Continuous_Dynamics::reach(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const std::vector<Constraint> & unsafeSet) const
{
	bool bSafetyChecking = false;

	if(unsafeSet.size() > 0)
	{
		bSafetyChecking = true;
	}

	std::vector<Flowpipe> initialSets;
	initialSets.push_back(initialSet);

	if(setting.tm_setting.queue_size > 0)
	{
		// symbolic remainder

		if(setting.tm_setting.step_min > 0)
		{
			// adaptive stepsizes
			result.status = reach_symbolic_remainder_adaptive_stepsize(result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes, result.num_of_flowpipes,
					setting.time, initialSets, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);
		}
		else if(setting.tm_setting.order_max > 0)
		{
			// adaptive orders
			result.status = reach_symbolic_remainder_adaptive_order(result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes, result.num_of_flowpipes,
					setting.time, initialSets, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);
		}
		else
		{
			// fixed stepsizes and orders
			result.status = reach_symbolic_remainder(result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes, result.num_of_flowpipes,
					setting.time, initialSets, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);
		}
	}
	else
	{
		if(setting.tm_setting.step_min > 0)
		{
			// adaptive stepsizes
			result.status = reach_adaptive_stepsize(result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes, result.num_of_flowpipes,
					setting.time, initialSets, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);
		}
		else if(setting.tm_setting.order_max > 0)
		{
			// adaptive orders
			result.status = reach_adaptive_order(result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes, result.num_of_flowpipes,
					setting.time, initialSets, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);
		}
		else
		{
			// fixed stepsizes and orders
			result.status = reach(result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes, result.num_of_flowpipes,
					setting.time, initialSets, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);
		}
	}

	Flowpipe fpTmp = result.nonlinear_flowpipes.back();
	result.fp_end_of_time = fpTmp;

	fpTmp.tmvPre.evaluate_time(result.fp_end_of_time.tmvPre, setting.tm_setting.step_end_exp_table);
}







Plot_Setting::Plot_Setting()
{
	type_of_file	= 0;
	type_of_object	= 0;
	num_of_pieces	= 0;
	bProjected		= false;
	bPrint			= true;
}

Plot_Setting::Plot_Setting(const Plot_Setting & setting)
{
	outputDims		= setting.outputDims;
	type_of_file	= setting.type_of_file;
	type_of_object	= setting.type_of_object;
	num_of_pieces	= setting.num_of_pieces;
	bProjected		= setting.bProjected;
	bPrint			= setting.bPrint;
}

Plot_Setting::~Plot_Setting()
{
}

Plot_Setting & Plot_Setting::operator = (const Plot_Setting & setting)
{
	if(this == &setting)
		return *this;

	outputDims		= setting.outputDims;
	type_of_file	= setting.type_of_file;
	type_of_object	= setting.type_of_object;
	num_of_pieces	= setting.num_of_pieces;
	bProjected		= setting.bProjected;
	bPrint			= setting.bPrint;

	return *this;
}

void Plot_Setting::setOutputDims(const unsigned int x, const unsigned int y)
{
	outputDims.clear();
	outputDims.push_back(x);
	outputDims.push_back(y);
}

void Plot_Setting::setOutputDims(const std::vector<unsigned int> & dims)
{
	outputDims = dims;
}

void Plot_Setting::setFileType(const unsigned int type)
{
	type_of_file = type;
}

void Plot_Setting::setObjectType(const unsigned int type)
{
	type_of_object = type;
}

void Plot_Setting::setNumOfPieces(const unsigned int n)
{
	num_of_pieces = n;
}

void Plot_Setting::printOn()
{
	bPrint = true;
}

void Plot_Setting::printOff()
{
	bPrint = false;
}

void Plot_Setting::plot_2D(const std::string & fileName, const Result_of_Reachability & result) const
{
	switch(type_of_file)
	{
	case PLOT_GNUPLOT:
		plot_2D_GNUPLOT(fileName, result);
		break;
	case PLOT_MATLAB:
		plot_2D_MATLAB(fileName, result);
		break;
	}
}

void Plot_Setting::plot_2D_MATLAB(const std::string & fileName, const Result_of_Reachability & result) const
{
	switch(type_of_object)
	{
	case PLOT_INTERVAL:
		plot_2D_interval_MATLAB(fileName, result);
		break;
	case PLOT_OCTAGON:
		plot_2D_octagon_MATLAB(fileName, result);
		break;
	case PLOT_GRID:
		plot_2D_grids_MATLAB(fileName, num_of_pieces, result);
		break;
	}
}

void Plot_Setting::plot_2D_interval_MATLAB(const std::string & fileName, const Result_of_Reachability & result) const
{
	if(bPrint)
	{
		printf("Generating the plot file...\n");
	}

	std::string name = outputDir + fileName + ".m";
	FILE *plotFile = fopen(name.c_str(), "w");

	if(plotFile == NULL)
	{
		printf("Can not create the output file.\n");
		exit(1);
	}

	std::vector<unsigned int> varIDs;
	if(bProjected)
	{
		varIDs.push_back(0);
		varIDs.push_back(1);
	}
	else
	{
		varIDs.push_back(outputDims[0]);
		varIDs.push_back(outputDims[1]);
	}

	unsigned int prog = 0;

	std::list<TaylorModelVec<Real> >::const_iterator tmvIter = result.tmv_flowpipes.begin();
	std::list<Flowpipe>::const_iterator fpIter = result.nonlinear_flowpipes.begin();
	std::list<int>::const_iterator safetyIter = result.safety_of_flowpipes.begin();

	unsigned int total_size = result.safety_of_flowpipes.size();

	if(total_size > 0)
	{
		for(; safetyIter != result.safety_of_flowpipes.end() ; ++tmvIter, ++fpIter, ++safetyIter)
		{
			std::vector<Interval> box;
			tmvIter->intEval(box, fpIter->domain, varIDs);

			Interval X = box[0], Y = box[1];

			switch(*safetyIter)
			{
			case SAFE:
				fprintf(plotFile,"plot( [%e,%e,%e,%e,%e] , [%e,%e,%e,%e,%e] , 'color' , '[0 0.4 0]');\nhold on;\nclear;\n",
						X.inf(), X.sup(), X.sup(), X.inf(), X.inf(), Y.inf(), Y.inf(), Y.sup(), Y.sup(), Y.inf());
				break;
			case UNSAFE:
				fprintf(plotFile,"plot( [%e,%e,%e,%e,%e] , [%e,%e,%e,%e,%e] , 'color' , '[1 0 0]');\nhold on;\nclear;\n",
						X.inf(), X.sup(), X.sup(), X.inf(), X.inf(), Y.inf(), Y.inf(), Y.sup(), Y.sup(), Y.inf());
				break;
			case UNKNOWN:
				fprintf(plotFile,"plot( [%e,%e,%e,%e,%e] , [%e,%e,%e,%e,%e] , 'color' , '[0 0 1]');\nhold on;\nclear;\n",
						X.inf(), X.sup(), X.sup(), X.inf(), X.inf(), Y.inf(), Y.inf(), Y.sup(), Y.sup(), Y.inf());
				break;
			}

			++prog;
			printf("\b\b\b\b");
			printf(BOLD_FONT "%%" RESET_COLOR);
			printf(BOLD_FONT "%3d" RESET_COLOR, (int)(prog*100/total_size));
			fflush(stdout);

			if(*safetyIter == UNSAFE)
			{
				break;
			}
		}

		printf("\b\b\b\b");
		printf(BOLD_FONT "%%100\n" RESET_COLOR);
		fflush(stdout);
	}
	else
	{
		total_size = result.tmv_flowpipes.size();

		for(; tmvIter != result.tmv_flowpipes.end() ; ++tmvIter, ++fpIter)
		{
			std::vector<Interval> box;
			tmvIter->intEval(box, fpIter->domain, varIDs);

			Interval X = box[0], Y = box[1];

			fprintf(plotFile,"plot( [%e,%e,%e,%e,%e] , [%e,%e,%e,%e,%e] , 'color' , '[0 0.4 0]');\nhold on;\nclear;\n",
					X.inf(), X.sup(), X.sup(), X.inf(), X.inf(), Y.inf(), Y.inf(), Y.sup(), Y.sup(), Y.inf());

			if(bPrint)
			{
				++prog;
				printf("\b\b\b\b");
				printf(BOLD_FONT "%%" RESET_COLOR);
				printf(BOLD_FONT "%3d" RESET_COLOR, (int)(prog*100/total_size));
				fflush(stdout);
			}
		}

		if(bPrint)
		{
			printf("\b\b\b\b");
			printf(BOLD_FONT "%%100\n" RESET_COLOR);
			fflush(stdout);
		}
	}

	fclose(plotFile);

	if(bPrint)
	{
		printf("Done.\n");
	}
}

void Plot_Setting::plot_2D_octagon_MATLAB(const std::string & fileName, const Result_of_Reachability & result) const
{
	if(bPrint)
	{
		printf("Generating the plot file...\n");
	}

	std::string name = outputDir + fileName + ".m";
	FILE *plotFile = fopen(name.c_str(), "w");

	if(plotFile == NULL)
	{
		printf("Can not create the output file.\n");
		exit(1);
	}

	int x, y;

	int rangeDim;

	if(bProjected)
	{
		x = 0;
		y = 1;
		rangeDim = 2;
	}
	else
	{
		x = outputDims[0];
		y = outputDims[1];
		rangeDim = stateVars.size();
	}

	std::vector<std::vector<Real> > output_poly_temp(8, std::vector<Real>(rangeDim, 0));

	output_poly_temp[0][x] = 1;
	output_poly_temp[1][y] = 1;
	output_poly_temp[2][x] = -1;
	output_poly_temp[3][y] = -1;
	output_poly_temp[4][x] = 1/sqrt(2);
	output_poly_temp[4][y] = 1/sqrt(2);
	output_poly_temp[5][x] = 1/sqrt(2);
	output_poly_temp[5][y] = -1/sqrt(2);
	output_poly_temp[6][x] = -1/sqrt(2);
	output_poly_temp[6][y] = 1/sqrt(2);
	output_poly_temp[7][x] = -1/sqrt(2);
	output_poly_temp[7][y] = -1/sqrt(2);

	// Construct the 2D template matrix.
	int rows = 8;
	int cols = rangeDim;

	std::vector<std::vector<Real> > sortedTemplate(rows, std::vector<Real>(cols, 0));
	std::vector<double> rowVec(cols);

	std::vector<std::vector<Real> > sortedRows;
	std::vector<std::vector<Real> >::iterator iterp, iterq;

	sortedRows.push_back(output_poly_temp[0]);

	bool bInserted;

	// Sort the row vectors in the template by anti-clockwise order (only in the x-y space).
	for(int i=1; i<rows; ++i)
	{
		iterp = sortedRows.begin();
		iterq = iterp;
		++iterq;
		bInserted = false;

		for(; iterq != sortedRows.end();)
		{
			Real tmp1 = output_poly_temp[i][x] * (*iterp)[y] - output_poly_temp[i][y] * (*iterp)[x];
			Real tmp2 = output_poly_temp[i][x] * (*iterq)[y] - output_poly_temp[i][y] * (*iterq)[x];

			if(tmp1 < 0 && tmp2 > 0)
			{
				sortedRows.insert(iterq, output_poly_temp[i]);
				bInserted = true;
				break;
			}
			else
			{
				++iterp;
				++iterq;
			}
		}

		if(!bInserted)
		{
			sortedRows.push_back(output_poly_temp[i]);
		}
	}

	iterp = sortedRows.begin();
	for(int i=0; i<rows; ++i, ++iterp)
	{
		for(int j=0; j<cols; ++j)
		{
			sortedTemplate[i][j] = (*iterp)[j];
		}
	}

	// Compute the intersections of two facets.
	// The vertices are ordered clockwisely.

	gsl_matrix *C = gsl_matrix_alloc(2,2);
	gsl_vector *d = gsl_vector_alloc(2);
	gsl_vector *vertex = gsl_vector_alloc(2);

	unsigned int prog = 0;

	std::list<TaylorModelVec<Real> >::const_iterator tmvIter = result.tmv_flowpipes.begin();
	std::list<Flowpipe>::const_iterator fpIter = result.nonlinear_flowpipes.begin();
	std::list<int>::const_iterator safetyIter = result.safety_of_flowpipes.begin();

	unsigned int total_size = result.safety_of_flowpipes.size();

	if(total_size > 0)
	{
		for(; safetyIter != result.safety_of_flowpipes.end(); ++tmvIter, ++fpIter, ++safetyIter)
		{
			Polyhedron polyTemplate(sortedTemplate, *tmvIter, fpIter->domain);

			double f1, f2;

			std::vector<LinearConstraint>::iterator iterp, iterq;
			iterp = iterq = polyTemplate.constraints.begin();
			++iterq;

			std::vector<double> vertices_x, vertices_y;

			for(; iterq != polyTemplate.constraints.end(); ++iterp, ++iterq)
			{
				gsl_matrix_set(C, 0, 0, iterp->A[x].toDouble());
				gsl_matrix_set(C, 0, 1, iterp->A[y].toDouble());
				gsl_matrix_set(C, 1, 0, iterq->A[x].toDouble());
				gsl_matrix_set(C, 1, 1, iterq->A[y].toDouble());

				gsl_vector_set(d, 0, iterp->B.toDouble());
				gsl_vector_set(d, 1, iterq->B.toDouble());

				gsl_linalg_HH_solve(C, d, vertex);

				double v1 = gsl_vector_get(vertex, 0);
				double v2 = gsl_vector_get(vertex, 1);

				if(iterp == polyTemplate.constraints.begin())
				{
					f1 = v1;
					f2 = v2;
				}

				vertices_x.push_back(v1);
				vertices_y.push_back(v2);
			}

			iterp = polyTemplate.constraints.begin();
			--iterq;

			gsl_matrix_set(C, 0, 0, iterp->A[x].toDouble());
			gsl_matrix_set(C, 0, 1, iterp->A[y].toDouble());
			gsl_matrix_set(C, 1, 0, iterq->A[x].toDouble());
			gsl_matrix_set(C, 1, 1, iterq->A[y].toDouble());

			gsl_vector_set(d, 0, iterp->B.toDouble());
			gsl_vector_set(d, 1, iterq->B.toDouble());

			gsl_linalg_HH_solve(C, d, vertex);

			double v1 = gsl_vector_get(vertex, 0);
			double v2 = gsl_vector_get(vertex, 1);

			vertices_x.push_back(v1);
			vertices_y.push_back(v2);
			vertices_x.push_back(f1);
			vertices_y.push_back(f2);

			fprintf(plotFile, "plot( ");

			fprintf(plotFile, "[ ");
			for(int i=0; i<vertices_x.size()-1; ++i)
			{
				fprintf(plotFile, "%e , ", vertices_x[i]);
			}
			fprintf(plotFile, "%e ] , ", vertices_x.back());

			fprintf(plotFile, "[ ");
			for(int i=0; i<vertices_y.size()-1; ++i)
			{
				fprintf(plotFile, "%e , ", vertices_y[i]);
			}
			fprintf(plotFile, "%e ] , ", vertices_y.back());

			switch(*safetyIter)
			{
			case SAFE:
				fprintf(plotFile, "'color' , '[0 0.4 0]');\nhold on;\nclear;\n");
				break;
			case UNSAFE:
				fprintf(plotFile, "'color' , '[1 0 0]');\nhold on;\nclear;\n");
				break;
			case UNKNOWN:
				fprintf(plotFile, "'color' , '[0 0 1]');\nhold on;\nclear;\n");
				break;
			}

			if(bPrint)
			{
				++prog;
				printf("\b\b\b\b");
				printf(BOLD_FONT "%%" RESET_COLOR);
				printf(BOLD_FONT "%3d" RESET_COLOR, (int)(prog*100/total_size));
				fflush(stdout);
			}

			if(*safetyIter == UNSAFE)
			{
				break;
			}
		}

		if(bPrint)
		{
			printf("\b\b\b\b");
			printf(BOLD_FONT "%%100\n" RESET_COLOR);
			fflush(stdout);
		}
	}
	else
	{
		total_size = result.tmv_flowpipes.size();

		for(; tmvIter != result.tmv_flowpipes.end(); ++tmvIter, ++fpIter)
		{
			Polyhedron polyTemplate(sortedTemplate, *tmvIter, fpIter->domain);

			double f1, f2;

			std::vector<LinearConstraint>::iterator iterp, iterq;
			iterp = iterq = polyTemplate.constraints.begin();
			++iterq;

			std::vector<double> vertices_x, vertices_y;

			for(; iterq != polyTemplate.constraints.end(); ++iterp, ++iterq)
			{
				gsl_matrix_set(C, 0, 0, iterp->A[x].toDouble());
				gsl_matrix_set(C, 0, 1, iterp->A[y].toDouble());
				gsl_matrix_set(C, 1, 0, iterq->A[x].toDouble());
				gsl_matrix_set(C, 1, 1, iterq->A[y].toDouble());

				gsl_vector_set(d, 0, iterp->B.toDouble());
				gsl_vector_set(d, 1, iterq->B.toDouble());

				gsl_linalg_HH_solve(C, d, vertex);

				double v1 = gsl_vector_get(vertex, 0);
				double v2 = gsl_vector_get(vertex, 1);

				if(iterp == polyTemplate.constraints.begin())
				{
					f1 = v1;
					f2 = v2;
				}

				vertices_x.push_back(v1);
				vertices_y.push_back(v2);
			}

			iterp = polyTemplate.constraints.begin();
			--iterq;

			gsl_matrix_set(C, 0, 0, iterp->A[x].toDouble());
			gsl_matrix_set(C, 0, 1, iterp->A[y].toDouble());
			gsl_matrix_set(C, 1, 0, iterq->A[x].toDouble());
			gsl_matrix_set(C, 1, 1, iterq->A[y].toDouble());

			gsl_vector_set(d, 0, iterp->B.toDouble());
			gsl_vector_set(d, 1, iterq->B.toDouble());

			gsl_linalg_HH_solve(C, d, vertex);

			double v1 = gsl_vector_get(vertex, 0);
			double v2 = gsl_vector_get(vertex, 1);

			vertices_x.push_back(v1);
			vertices_y.push_back(v2);
			vertices_x.push_back(f1);
			vertices_y.push_back(f2);

			fprintf(plotFile, "plot( ");

			fprintf(plotFile, "[ ");
			for(int i=0; i<vertices_x.size()-1; ++i)
			{
				fprintf(plotFile, "%le , ", vertices_x[i]);
			}
			fprintf(plotFile, "%le ] , ", vertices_x.back());

			fprintf(plotFile, "[ ");
			for(int i=0; i<vertices_y.size()-1; ++i)
			{
				fprintf(plotFile, "%e , ", vertices_y[i]);
			}
			fprintf(plotFile, "%e ] , ", vertices_y.back());


			fprintf(plotFile, "'color' , '[0 0.4 0]');\nhold on;\nclear;\n");

			if(bPrint)
			{
				++prog;
				printf("\b\b\b\b");
				printf(BOLD_FONT "%%" RESET_COLOR);
				printf(BOLD_FONT "%3d" RESET_COLOR, (int)(prog*100/total_size));
				fflush(stdout);
			}
		}

		if(bPrint)
		{
			printf("\b\b\b\b");
			printf(BOLD_FONT "%%100\n" RESET_COLOR);
			fflush(stdout);
		}
	}

	gsl_matrix_free(C);
	gsl_vector_free(d);
	gsl_vector_free(vertex);

	fclose(plotFile);

	if(bPrint)
	{
		printf("Done.\n");
	}
}

void Plot_Setting::plot_2D_grids_MATLAB(const std::string & fileName, const unsigned int num, const Result_of_Reachability & result) const
{
	if(bPrint)
	{
		printf("Generating the plot file...\n");
	}

	std::string name = outputDir + fileName + ".m";
	FILE *plotFile = fopen(name.c_str(), "w");

	if(plotFile == NULL)
	{
		printf("Can not create the output file.\n");
		exit(1);
	}

	int x, y;
	if(bProjected)
	{
		x = 0;
		y = 1;
	}
	else
	{
		x = outputDims[0];
		y = outputDims[1];
	}

	unsigned int prog = 0;

	std::list<TaylorModelVec<Real> >::const_iterator tmvIter = result.tmv_flowpipes.begin();
	std::list<Flowpipe>::const_iterator fpIter = result.nonlinear_flowpipes.begin();
	std::list<int>::const_iterator safetyIter = result.safety_of_flowpipes.begin();

	unsigned int total_size = result.safety_of_flowpipes.size();

	if(total_size > 0)
	{
		for(; safetyIter != result.safety_of_flowpipes.end(); ++tmvIter, ++fpIter, ++safetyIter)
		{
			// decompose the domain
			std::list<std::vector<Interval> > grids;

			gridBox(grids, fpIter->domain, num);

			// we only consider the output dimensions
			HornerForm<Real> hfOutputX;
			Interval remainderX;

			tmvIter->tms[x].toHornerForm(hfOutputX, remainderX);


			HornerForm<Real> hfOutputY;
			Interval remainderY;

			tmvIter->tms[y].toHornerForm(hfOutputY, remainderY);


			// evaluate the images from all of the grids
			std::list<std::vector<Interval> >::const_iterator gIter = grids.begin();
			for(; gIter!=grids.end(); ++gIter)
			{
				Interval X;
				hfOutputX.evaluate(X, *gIter);
				X += remainderX;

				Interval Y;
				hfOutputY.evaluate(Y, *gIter);
				Y += remainderY;

				switch(*safetyIter)
				{
				case SAFE:
					fprintf(plotFile,"plot( [%e,%e,%e,%e,%e] , [%e,%e,%e,%e,%e] , 'color' , '[0 0.4 0]');\nhold on;\nclear;\n",
							X.inf(), X.sup(), X.sup(), X.inf(), X.inf(), Y.inf(), Y.inf(), Y.sup(), Y.sup(), Y.inf());
					break;
				case UNSAFE:
					fprintf(plotFile,"plot( [%e,%e,%e,%e,%e] , [%e,%e,%e,%e,%e] , 'color' , '[1 0 0]');\nhold on;\nclear;\n",
							X.inf(), X.sup(), X.sup(), X.inf(), X.inf(), Y.inf(), Y.inf(), Y.sup(), Y.sup(), Y.inf());
					break;
				case UNKNOWN:
					fprintf(plotFile,"plot( [%e,%e,%e,%e,%e] , [%e,%e,%e,%e,%e] , 'color' , '[0 0 1]');\nhold on;\nclear;\n",
							X.inf(), X.sup(), X.sup(), X.inf(), X.inf(), Y.inf(), Y.inf(), Y.sup(), Y.sup(), Y.inf());
					break;
				}
			}

			if(bPrint)
			{
				++prog;
				printf("\b\b\b\b");
				printf(BOLD_FONT "%%" RESET_COLOR);
				printf(BOLD_FONT "%3d" RESET_COLOR, (int)(prog*100/total_size));
				fflush(stdout);
			}

			if(*safetyIter == UNSAFE)
			{
				break;
			}
		}

		if(bPrint)
		{
			printf("\b\b\b\b");
			printf(BOLD_FONT "%%100\n" RESET_COLOR);
			fflush(stdout);
		}
	}
	else
	{
		total_size = result.tmv_flowpipes.size();

		for(; tmvIter != result.tmv_flowpipes.end(); ++tmvIter, ++fpIter)
		{
			// decompose the domain
			std::list<std::vector<Interval> > grids;

			gridBox(grids, fpIter->domain, num);

			// we only consider the output dimensions
			HornerForm<Real> hfOutputX;
			Interval remainderX;

			tmvIter->tms[x].toHornerForm(hfOutputX, remainderX);


			HornerForm<Real> hfOutputY;
			Interval remainderY;

			tmvIter->tms[y].toHornerForm(hfOutputY, remainderY);


			// evaluate the images from all of the grids
			std::list<std::vector<Interval> >::const_iterator gIter = grids.begin();
			for(; gIter!=grids.end(); ++gIter)
			{
				Interval X;
				hfOutputX.evaluate(X, *gIter);
				X += remainderX;

				Interval Y;
				hfOutputY.evaluate(Y, *gIter);
				Y += remainderY;

				fprintf(plotFile, "plot( [%e,%e,%e,%e,%e] , [%e,%e,%e,%e,%e] , 'color' , '[0 0.4 0]');\nhold on;\nclear;\n",
						X.inf(), X.sup(), X.sup(), X.inf(), X.inf(), Y.inf(), Y.inf(), Y.sup(), Y.sup(), Y.inf());
			}

			if(bPrint)
			{
				++prog;
				printf("\b\b\b\b");
				printf(BOLD_FONT "%%" RESET_COLOR);
				printf(BOLD_FONT "%3d" RESET_COLOR, (int)(prog*100/total_size));
				fflush(stdout);
			}
		}

		if(bPrint)
		{
			printf("\b\b\b\b");
			printf(BOLD_FONT "%%100\n" RESET_COLOR);
			fflush(stdout);
		}
	}

	fclose(plotFile);

	if(bPrint)
	{
		printf("Done.\n");
	}
}

void Plot_Setting::plot_2D_GNUPLOT(const std::string & fileName, const Result_of_Reachability & result) const
{
	switch(type_of_object)
	{
	case PLOT_INTERVAL:
		plot_2D_interval_GNUPLOT(fileName, result);
		break;
	case PLOT_OCTAGON:
		plot_2D_octagon_GNUPLOT(fileName, result);
		break;
	case PLOT_GRID:
		plot_2D_grids_GNUPLOT(fileName, num_of_pieces, result);
		break;
	}
}

void Plot_Setting::plot_2D_interval_GNUPLOT(const std::string & fileName, const Result_of_Reachability & result) const
{
	if(bPrint)
	{
		printf("Generating the plot file...\n");
	}

	std::string plot_file_name = outputDir + fileName + ".plt";
	FILE *plotFile = fopen(plot_file_name.c_str(), "w");

	if(plotFile == NULL)
	{
		printf("Can not create the output file.\n");
		exit(1);
	}

	std::string image_file_name = imageDir + fileName + ".eps";

	fprintf(plotFile, "set terminal postscript enhanced color\n");

	fprintf(plotFile, "set output '%s'\n", image_file_name.c_str());

	fprintf(plotFile, "set style line 1 linecolor rgb \"blue\"\n");
	fprintf(plotFile, "set autoscale\n");
	fprintf(plotFile, "unset label\n");
	fprintf(plotFile, "set xtic auto\n");
	fprintf(plotFile, "set ytic auto\n");
	fprintf(plotFile, "set xlabel \"%s\"\n", stateVars.varNames[outputDims[0]].c_str());
	fprintf(plotFile, "set ylabel \"%s\"\n", stateVars.varNames[outputDims[1]].c_str());
	fprintf(plotFile, "plot '-' notitle with lines ls 1\n");

	std::vector<unsigned int> varIDs;
	if(bProjected)
	{
		varIDs.push_back(0);
		varIDs.push_back(1);
	}
	else
	{
		varIDs.push_back(outputDims[0]);
		varIDs.push_back(outputDims[1]);
	}

	unsigned int prog = 0;

	std::list<TaylorModelVec<Real> >::const_iterator tmvIter = result.tmv_flowpipes.begin();
	std::list<Flowpipe>::const_iterator fpIter = result.nonlinear_flowpipes.begin();
	std::list<int>::const_iterator safetyIter = result.safety_of_flowpipes.begin();

	unsigned int total_size = result.safety_of_flowpipes.size();

	if(total_size > 0)
	{
		for(; safetyIter != result.safety_of_flowpipes.end(); ++tmvIter, ++fpIter, ++safetyIter)
		{
			std::vector<Interval> box;
			tmvIter->intEval(box, fpIter->domain, varIDs);

			// output the vertices
			fprintf(plotFile, "%e %e\n", box[0].inf(), box[1].inf());
			fprintf(plotFile, "%e %e\n", box[0].sup(), box[1].inf());
			fprintf(plotFile, "%e %e\n", box[0].sup(), box[1].sup());
			fprintf(plotFile, "%e %e\n", box[0].inf(), box[1].sup());
			fprintf(plotFile, "%e %e\n", box[0].inf(), box[1].inf());
			fprintf(plotFile, "\n\n");

			if(bPrint)
			{
				++prog;
				printf("\b\b\b\b");
				printf(BOLD_FONT "%%" RESET_COLOR);
				printf(BOLD_FONT "%3d" RESET_COLOR, (int)(prog*100/total_size));
				fflush(stdout);
			}

			if(*safetyIter == UNSAFE)
			{
				break;
			}
		}

		if(bPrint)
		{
			fprintf(plotFile, "e\n");
			printf("\b\b\b\b");
			printf(BOLD_FONT "%%100\n" RESET_COLOR);
			fflush(stdout);
		}
	}
	else
	{
		total_size = result.tmv_flowpipes.size();

		for(; tmvIter != result.tmv_flowpipes.end(); ++tmvIter, ++fpIter)
		{
			std::vector<Interval> box;
			tmvIter->intEval(box, fpIter->domain, varIDs);

			// output the vertices
			fprintf(plotFile, "%e %e\n", box[0].inf(), box[1].inf());
			fprintf(plotFile, "%e %e\n", box[0].sup(), box[1].inf());
			fprintf(plotFile, "%e %e\n", box[0].sup(), box[1].sup());
			fprintf(plotFile, "%e %e\n", box[0].inf(), box[1].sup());
			fprintf(plotFile, "%e %e\n", box[0].inf(), box[1].inf());
			fprintf(plotFile, "\n\n");

			if(bPrint)
			{
				++prog;
				printf("\b\b\b\b");
				printf(BOLD_FONT "%%" RESET_COLOR);
				printf(BOLD_FONT "%3d" RESET_COLOR, (int)(prog*100/total_size));
				fflush(stdout);
			}
		}

		if(bPrint)
		{
			fprintf(plotFile, "e\n");
			printf("\b\b\b\b");
			printf(BOLD_FONT "%%100\n" RESET_COLOR);
			fflush(stdout);
		}
	}

	fclose(plotFile);

	if(bPrint)
	{
		printf("Done.\n");
	}
}

void Plot_Setting::plot_2D_octagon_GNUPLOT(const std::string & fileName, const Result_of_Reachability & result) const
{
	if(bPrint)
	{
		printf("Generating the plot file...\n");
	}

	std::string name = outputDir + fileName + ".plt";
	FILE *plotFile = fopen(name.c_str(), "w");

	if(plotFile == NULL)
	{
		printf("Can not create the output file.\n");
		exit(1);
	}

	int x, y;

	int rangeDim;

	if(bProjected)
	{
		x = 0;
		y = 1;
		rangeDim = 2;
	}
	else
	{
		x = outputDims[0];
		y = outputDims[1];
		rangeDim = stateVars.size();
	}

	std::vector<std::vector<Real> > output_poly_temp(8, std::vector<Real>(rangeDim, 0));

	output_poly_temp[0][x] = 1;
	output_poly_temp[1][y] = 1;
	output_poly_temp[2][x] = -1;
	output_poly_temp[3][y] = -1;
	output_poly_temp[4][x] = 1/sqrt(2);
	output_poly_temp[4][y] = 1/sqrt(2);
	output_poly_temp[5][x] = 1/sqrt(2);
	output_poly_temp[5][y] = -1/sqrt(2);
	output_poly_temp[6][x] = -1/sqrt(2);
	output_poly_temp[6][y] = 1/sqrt(2);
	output_poly_temp[7][x] = -1/sqrt(2);
	output_poly_temp[7][y] = -1/sqrt(2);

	// Construct the 2D template matrix.
	int rows = 8;
	int cols = rangeDim;

	std::vector<std::vector<Real> > sortedTemplate(rows, std::vector<Real>(cols, 0));
	std::vector<Real> rowVec(cols);

	std::vector<std::vector<Real> > sortedRows;
	std::vector<std::vector<Real> >::iterator iterp, iterq;

	sortedRows.push_back(output_poly_temp[0]);

	bool bInserted;

	// Sort the row vectors in the template by anti-clockwise order (only in the x-y space).
	for(int i=1; i<rows; ++i)
	{
		iterp = sortedRows.begin();
		iterq = iterp;
		++iterq;
		bInserted = false;

		for(; iterq != sortedRows.end();)
		{
			Real tmp1 = output_poly_temp[i][x] * (*iterp)[y] - output_poly_temp[i][y] * (*iterp)[x];
			Real tmp2 = output_poly_temp[i][x] * (*iterq)[y] - output_poly_temp[i][y] * (*iterq)[x];

			if(tmp1 < 0 && tmp2 > 0)
			{
				sortedRows.insert(iterq, output_poly_temp[i]);
				bInserted = true;
				break;
			}
			else
			{
				++iterp;
				++iterq;
			}
		}

		if(!bInserted)
		{
			sortedRows.push_back(output_poly_temp[i]);
		}
	}

	iterp = sortedRows.begin();
	for(int i=0; i<rows; ++i, ++iterp)
	{
		for(int j=0; j<cols; ++j)
		{
			sortedTemplate[i][j] = (*iterp)[j];
		}
	}

	std::string image_file_name = imageDir + fileName + ".eps";

	fprintf(plotFile, "set terminal postscript enhanced color\n");

	fprintf(plotFile, "set output '%s'\n", image_file_name.c_str());

	fprintf(plotFile, "set style line 1 linecolor rgb \"blue\"\n");
	fprintf(plotFile, "set autoscale\n");
	fprintf(plotFile, "unset label\n");
	fprintf(plotFile, "set xtic auto\n");
	fprintf(plotFile, "set ytic auto\n");
	fprintf(plotFile, "set xlabel \"%s\"\n", stateVars.varNames[outputDims[0]].c_str());
	fprintf(plotFile, "set ylabel \"%s\"\n", stateVars.varNames[outputDims[1]].c_str());
	fprintf(plotFile, "plot '-' notitle with lines ls 1\n");

	// Compute the intersections of two facets.
	// The vertices are ordered clockwisely.

	gsl_matrix *C = gsl_matrix_alloc(2,2);
	gsl_vector *d = gsl_vector_alloc(2);
	gsl_vector *vertex = gsl_vector_alloc(2);

	unsigned int prog = 0;

	std::list<TaylorModelVec<Real> >::const_iterator tmvIter = result.tmv_flowpipes.begin();
	std::list<Flowpipe>::const_iterator fpIter = result.nonlinear_flowpipes.begin();
	std::list<int>::const_iterator safetyIter = result.safety_of_flowpipes.begin();

	unsigned int total_size = result.safety_of_flowpipes.size();

	if(total_size > 0)
	{
		for(; safetyIter != result.safety_of_flowpipes.end(); ++tmvIter, ++fpIter, ++safetyIter)
		{
			Polyhedron polyTemplate(sortedTemplate, *tmvIter, fpIter->domain);

			double f1, f2;

			std::vector<LinearConstraint>::iterator iterp, iterq;
			iterp = iterq = polyTemplate.constraints.begin();
			++iterq;

			for(; iterq != polyTemplate.constraints.end(); ++iterp, ++iterq)
			{
				gsl_matrix_set(C, 0, 0, iterp->A[x].toDouble());
				gsl_matrix_set(C, 0, 1, iterp->A[y].toDouble());
				gsl_matrix_set(C, 1, 0, iterq->A[x].toDouble());
				gsl_matrix_set(C, 1, 1, iterq->A[y].toDouble());

				gsl_vector_set(d, 0, iterp->B.toDouble());
				gsl_vector_set(d, 1, iterq->B.toDouble());

				gsl_linalg_HH_solve(C, d, vertex);

				double v1 = gsl_vector_get(vertex, 0);
				double v2 = gsl_vector_get(vertex, 1);

				if(iterp == polyTemplate.constraints.begin())
				{
					f1 = v1;
					f2 = v2;
				}

				fprintf(plotFile, "%e %e\n", v1, v2);
			}

			iterp = polyTemplate.constraints.begin();
			--iterq;

			gsl_matrix_set(C, 0, 0, iterp->A[x].toDouble());
			gsl_matrix_set(C, 0, 1, iterp->A[y].toDouble());
			gsl_matrix_set(C, 1, 0, iterq->A[x].toDouble());
			gsl_matrix_set(C, 1, 1, iterq->A[y].toDouble());

			gsl_vector_set(d, 0, iterp->B.toDouble());
			gsl_vector_set(d, 1, iterq->B.toDouble());

			gsl_linalg_HH_solve(C, d, vertex);

			double v1 = gsl_vector_get(vertex, 0);
			double v2 = gsl_vector_get(vertex, 1);

			fprintf(plotFile, "%e %e\n", v1, v2);

			fprintf(plotFile, "%e %e\n", f1, f2);
			fprintf(plotFile, "\n\n");

			if(bPrint)
			{
				++prog;
				printf("\b\b\b\b");
				printf(BOLD_FONT "%%" RESET_COLOR);
				printf(BOLD_FONT "%3d" RESET_COLOR, (int)(prog*100/total_size));
				fflush(stdout);
			}

			if(*safetyIter == UNSAFE)
			{
				break;
			}
		}

		if(bPrint)
		{
			fprintf(plotFile, "e\n");
			printf("\b\b\b\b");
			printf(BOLD_FONT "%%100\n" RESET_COLOR);
			fflush(stdout);
		}
	}
	else
	{
		total_size = result.tmv_flowpipes.size();

		for(; tmvIter != result.tmv_flowpipes.end(); ++tmvIter, ++fpIter)
		{
			Polyhedron polyTemplate(sortedTemplate, *tmvIter, fpIter->domain);

			double f1, f2;

			std::vector<LinearConstraint>::iterator iterp, iterq;
			iterp = iterq = polyTemplate.constraints.begin();
			++iterq;

			for(; iterq != polyTemplate.constraints.end(); ++iterp, ++iterq)
			{
				gsl_matrix_set(C, 0, 0, iterp->A[x].toDouble());
				gsl_matrix_set(C, 0, 1, iterp->A[y].toDouble());
				gsl_matrix_set(C, 1, 0, iterq->A[x].toDouble());
				gsl_matrix_set(C, 1, 1, iterq->A[y].toDouble());

				gsl_vector_set(d, 0, iterp->B.toDouble());
				gsl_vector_set(d, 1, iterq->B.toDouble());

				gsl_linalg_HH_solve(C, d, vertex);

				double v1 = gsl_vector_get(vertex, 0);
				double v2 = gsl_vector_get(vertex, 1);

				if(iterp == polyTemplate.constraints.begin())
				{
					f1 = v1;
					f2 = v2;
				}

				fprintf(plotFile, "%e %e\n", v1, v2);
			}

			iterp = polyTemplate.constraints.begin();
			--iterq;

			gsl_matrix_set(C, 0, 0, iterp->A[x].toDouble());
			gsl_matrix_set(C, 0, 1, iterp->A[y].toDouble());
			gsl_matrix_set(C, 1, 0, iterq->A[x].toDouble());
			gsl_matrix_set(C, 1, 1, iterq->A[y].toDouble());

			gsl_vector_set(d, 0, iterp->B.toDouble());
			gsl_vector_set(d, 1, iterq->B.toDouble());

			gsl_linalg_HH_solve(C, d, vertex);

			double v1 = gsl_vector_get(vertex, 0);
			double v2 = gsl_vector_get(vertex, 1);

			fprintf(plotFile, "%e %e\n", v1, v2);

			fprintf(plotFile, "%e %e\n", f1, f2);
			fprintf(plotFile, "\n\n");

			if(bPrint)
			{
				++prog;
				printf("\b\b\b\b");
				printf(BOLD_FONT "%%" RESET_COLOR);
				printf(BOLD_FONT "%3d" RESET_COLOR, (int)(prog*100/total_size));
				fflush(stdout);
			}

			if(*safetyIter == UNSAFE)
			{
				break;
			}
		}

		if(bPrint)
		{
			fprintf(plotFile, "e\n");
			printf("\b\b\b\b");
			printf(BOLD_FONT "%%100\n" RESET_COLOR);
			fflush(stdout);
		}
	}

	gsl_matrix_free(C);
	gsl_vector_free(d);
	gsl_vector_free(vertex);

	fclose(plotFile);

	if(bPrint)
	{
		printf("Done.\n");
	}
}

void Plot_Setting::plot_2D_grids_GNUPLOT(const std::string & fileName, const unsigned int num, const Result_of_Reachability & result) const
{
	if(bPrint)
	{
		printf("Generating the plot file...\n");
	}

	std::string name = outputDir + fileName + ".plt";
	FILE *plotFile = fopen(name.c_str(), "w");

	if(plotFile == NULL)
	{
		printf("Can not create the output file.\n");
		exit(1);
	}

	int x, y;
	if(bProjected)
	{
		x = 0;
		y = 1;
	}
	else
	{
		x = outputDims[0];
		y = outputDims[1];
	}

	std::string image_file_name = imageDir + fileName + ".eps";

	fprintf(plotFile, "set terminal postscript enhanced color\n");

	fprintf(plotFile, "set output '%s'\n", image_file_name.c_str());

	fprintf(plotFile, "set style line 1 linecolor rgb \"blue\"\n");
	fprintf(plotFile, "set autoscale\n");
	fprintf(plotFile, "unset label\n");
	fprintf(plotFile, "set xtic auto\n");
	fprintf(plotFile, "set ytic auto\n");
	fprintf(plotFile, "set xlabel \"%s\"\n", stateVars.varNames[outputDims[0]].c_str());
	fprintf(plotFile, "set ylabel \"%s\"\n", stateVars.varNames[outputDims[1]].c_str());
	fprintf(plotFile, "plot '-' notitle with lines ls 1\n");

	unsigned int prog = 0;

	std::list<TaylorModelVec<Real> >::const_iterator tmvIter = result.tmv_flowpipes.begin();
	std::list<Flowpipe>::const_iterator fpIter = result.nonlinear_flowpipes.begin();
	std::list<int>::const_iterator safetyIter = result.safety_of_flowpipes.begin();

	unsigned int total_size = result.safety_of_flowpipes.size();

	if(total_size > 0)
	{
		for(; safetyIter != result.safety_of_flowpipes.end(); ++tmvIter, ++fpIter, ++safetyIter)
		{
			// decompose the domain
			std::list<std::vector<Interval> > grids;

			gridBox(grids, fpIter->domain, num);

			// we only consider the output dimensions
			HornerForm<Real> hfOutputX;
			Interval remainderX;

			tmvIter->tms[x].toHornerForm(hfOutputX, remainderX);

			HornerForm<Real> hfOutputY;
			Interval remainderY;

			tmvIter->tms[y].toHornerForm(hfOutputY, remainderY);

			// evaluate the images from all of the grids
			std::list<std::vector<Interval> >::const_iterator gIter = grids.begin();
			for(; gIter!=grids.end(); ++gIter)
			{
				Interval X;
				hfOutputX.evaluate(X, *gIter);
				X += remainderX;

				Interval Y;
				hfOutputY.evaluate(Y, *gIter);
				Y += remainderY;

				// output the vertices
				fprintf(plotFile, "%e %e\n", X.inf(), Y.inf());
				fprintf(plotFile, "%e %e\n", X.sup(), Y.inf());
				fprintf(plotFile, "%e %e\n", X.sup(), Y.sup());
				fprintf(plotFile, "%e %e\n", X.inf(), Y.sup());
				fprintf(plotFile, "%e %e\n", X.inf(), Y.inf());
				fprintf(plotFile, "\n\n");
			}

			if(bPrint)
			{
				++prog;
				printf("\b\b\b\b");
				printf(BOLD_FONT "%%" RESET_COLOR);
				printf(BOLD_FONT "%3d" RESET_COLOR, (int)(prog*100/total_size));
				fflush(stdout);
			}

			if(*safetyIter == UNSAFE)
			{
				break;
			}
		}

		if(bPrint)
		{
			fprintf(plotFile, "e\n");
			printf("\b\b\b\b");
			printf(BOLD_FONT "%%100\n" RESET_COLOR);
			fflush(stdout);
		}
	}
	else
	{
		total_size = result.tmv_flowpipes.size();
		for(; tmvIter != result.tmv_flowpipes.end(); ++tmvIter, ++fpIter)
		{
			// decompose the domain
			std::list<std::vector<Interval> > grids;

			gridBox(grids, fpIter->domain, num);

			// we only consider the output dimensions
			HornerForm<Real> hfOutputX;
			Interval remainderX;

			tmvIter->tms[x].toHornerForm(hfOutputX, remainderX);

			HornerForm<Real> hfOutputY;
			Interval remainderY;

			tmvIter->tms[y].toHornerForm(hfOutputY, remainderY);

			// evaluate the images from all of the grids
			std::list<std::vector<Interval> >::const_iterator gIter = grids.begin();
			for(; gIter!=grids.end(); ++gIter)
			{
				Interval X;
				hfOutputX.evaluate(X, *gIter);
				X += remainderX;

				Interval Y;
				hfOutputY.evaluate(Y, *gIter);
				Y += remainderY;

				// output the vertices
				fprintf(plotFile, "%e %e\n", X.inf(), Y.inf());
				fprintf(plotFile, "%e %e\n", X.sup(), Y.inf());
				fprintf(plotFile, "%e %e\n", X.sup(), Y.sup());
				fprintf(plotFile, "%e %e\n", X.inf(), Y.sup());
				fprintf(plotFile, "%e %e\n", X.inf(), Y.inf());
				fprintf(plotFile, "\n\n");
			}

			if(bPrint)
			{
				++prog;
				printf("\b\b\b\b");
				printf(BOLD_FONT "%%" RESET_COLOR);
				printf(BOLD_FONT "%3d" RESET_COLOR, (int)(prog*100/total_size));
				fflush(stdout);
			}
		}

		if(bPrint)
		{
			fprintf(plotFile, "e\n");
			printf("\b\b\b\b");
			printf(BOLD_FONT "%%100\n" RESET_COLOR);
			fflush(stdout);
		}
	}

	fclose(plotFile);

	if(bPrint)
	{
		printf("Done.\n");
	}
}














Continuous_Reachability_Problem_Description::Continuous_Reachability_Problem_Description()
{
	step_max = 0;
	step_min = 0;
	order_min = 0;
	order_max = 0;
	time = 0;

	bSafetyChecking = false;
	bPlot = false;
	bTMOutput = false;
	bPrint = false;

	type_of_dynamics = -1;
	bDeterministic = true;
	bSymbolicRemainder = false;
	queue_size = 0;
}

Continuous_Reachability_Problem_Description::Continuous_Reachability_Problem_Description(const Continuous_Reachability_Problem_Description & description)
{
	step_max							= description.step_max;
	step_min							= description.step_min;
	order_min							= description.order_min;
	order_max							= description.order_max;
	time								= description.time;

	bSafetyChecking						= description.bSafetyChecking;
	bPlot								= description.bPlot;
	bTMOutput							= description.bTMOutput;
	bPrint								= description.bPrint;

	remainder_estimation				= description.remainder_estimation;
	cutoff_threshold					= description.cutoff_threshold;

	unsafeSet							= description.unsafeSet;

	type_of_dynamics					= description.type_of_dynamics;
	bDeterministic						= description.bDeterministic;
	fileName							= description.fileName;

	deterministic_dynamics				= description.deterministic_dynamics;
	nondeterministic_dynamics			= description.nondeterministic_dynamics;
	bSymbolicRemainder					= description.bSymbolicRemainder;
	queue_size							= description.queue_size;
}

Continuous_Reachability_Problem_Description::~Continuous_Reachability_Problem_Description()
{
}

Continuous_Reachability_Problem_Description & Continuous_Reachability_Problem_Description::operator = (const Continuous_Reachability_Problem_Description & description)
{
	if(this == &description)
		return *this;

	step_max							= description.step_max;
	step_min							= description.step_min;
	order_min							= description.order_min;
	order_max							= description.order_max;
	time								= description.time;

	bSafetyChecking						= description.bSafetyChecking;
	bPlot								= description.bPlot;
	bTMOutput							= description.bTMOutput;
	bPrint								= description.bPrint;

	remainder_estimation				= description.remainder_estimation;
	cutoff_threshold					= description.cutoff_threshold;

	unsafeSet							= description.unsafeSet;

	type_of_dynamics					= description.type_of_dynamics;
	bDeterministic						= description.bDeterministic;
	fileName							= description.fileName;

	deterministic_dynamics				= description.deterministic_dynamics;
	nondeterministic_dynamics			= description.nondeterministic_dynamics;
	bSymbolicRemainder					= description.bSymbolicRemainder;
	queue_size							= description.queue_size;

	return *this;
}

void Continuous_Reachability_Problem_Description::setStateVars(const Variables & vars)
{
	stateVars = vars;
}

void Continuous_Reachability_Problem_Description::setTMVars(const Variables & vars)
{
	tmVars = vars;
}

bool Continuous_Reachability_Problem_Description::setTimeHorizon(const double t)
{
	if(t <= 0)
	{
		std::cout << "The time horizon should be positive." << std::endl;
		return false;
	}
	else
	{
		time = t;
		return true;
	}
}

bool Continuous_Reachability_Problem_Description::setFixedStepsize(const double delta)
{
	if(delta <= 0)
	{
		std::cout << "The stepsize should be positive." << std::endl;
		return false;
	}
	else
	{
		step_max = delta;
		return true;
	}
}

bool Continuous_Reachability_Problem_Description::setAdaptiveStepsize(const double delta_min, const double delta_max)
{
	if(delta_min <= 0 || delta_max <= 0)
	{
		std::cout << "The stepsize should be positive." << std::endl;
		return false;
	}
	else if(delta_min >= delta_max)
	{
		std::cout << "The minimum stepsize should be less than the maximum one." << std::endl;
		return false;
	}
	else
	{
		step_min = delta_min;
		step_max = delta_max;
		return true;
	}
}

bool Continuous_Reachability_Problem_Description::setFixedOrder(const unsigned int k)
{
	if(k < 2)
	{
		std::cout << "The order should be larger than 1." << std::endl;
		return false;
	}
	else
	{
		order_min = k;
		return true;
	}
}

bool Continuous_Reachability_Problem_Description::setAdaptiveOrder(const unsigned int k_min, const unsigned int k_max)
{
	if(k_min < 2 || k_max < 2)
	{
		std::cout << "The order should be larger than 1." << std::endl;
		return false;
	}
	else if(k_min >= k_max)
	{
		std::cout << "The lower order should be less than the higher one." << std::endl;
		return false;
	}
	else
	{
		order_min = k_min;
		order_max = k_max;
		return true;
	}
}

bool Continuous_Reachability_Problem_Description::setRemainderEstimation(const std::vector<Interval> & intVec)
{
	remainder_estimation = intVec;
	return true;
}

bool Continuous_Reachability_Problem_Description::setCutoff(const Interval & cutoff)
{
	cutoff_threshold = cutoff;
	return true;
}

bool Continuous_Reachability_Problem_Description::setPrecision(const unsigned int prec)
{
	if(prec < 53)
	{
		std::cout << "The precision should be at least 53." << std::endl;
		return false;
	}
	else
	{
		intervalNumPrecision = prec;
		return true;
	}
}

void Continuous_Reachability_Problem_Description::printOn()
{
	bPrint = true;
}

void Continuous_Reachability_Problem_Description::printOff()
{
	bPrint = false;
}

void Continuous_Reachability_Problem_Description::safetyCheckingOn()
{
	bSafetyChecking = true;
}

void Continuous_Reachability_Problem_Description::safetyCheckingOff()
{
	bSafetyChecking = false;
}

void Continuous_Reachability_Problem_Description::setOutputDims(const unsigned int x, const unsigned int y)
{
	std::vector<unsigned int> outputDims;
	outputDims.push_back(x);
	outputDims.push_back(y);

	plot_setting.setOutputDims(outputDims);
}

void Continuous_Reachability_Problem_Description::setFileType(const unsigned int type)
{
	plot_setting.setFileType(type);
}

void Continuous_Reachability_Problem_Description::setObjectType(const unsigned int type)
{
	plot_setting.setObjectType(type);
}

void Continuous_Reachability_Problem_Description::setNumOfPieces(const unsigned int n)
{
	plot_setting.setNumOfPieces(n);
}

void Continuous_Reachability_Problem_Description::plotOn()
{
	bPlot = true;
}

void Continuous_Reachability_Problem_Description::plotOff()
{
	bPlot = false;
}

void Continuous_Reachability_Problem_Description::tmOutputOn()
{
	bTMOutput = true;
}

void Continuous_Reachability_Problem_Description::tmOutputOff()
{
	bTMOutput = false;
}

void Continuous_Reachability_Problem_Description::setUnsafe(const std::vector<Constraint> & unsafe_constraints)
{
	unsafeSet = unsafe_constraints;
}

void Continuous_Reachability_Problem_Description::setInitialSets(const std::vector<Flowpipe> & flowpipes)
{
	initialSets = flowpipes;
}

void Continuous_Reachability_Problem_Description::setFileName(const std::string & str)
{
	fileName = str;
}












Continuous_Reachability::Continuous_Reachability()
{
	pDynamics			= NULL;
	p_tm_setting		= NULL;
	p_g_setting			= NULL;
	p_p_setting			= NULL;

	type_of_dynamics	= 0;
	time				= 0;
	bSafetyChecking		= false;
	bPlot				= false;
	bTMOutput			= false;
	bPrint				= false;
	bSymbolicRemainder	= false;
}

Continuous_Reachability::Continuous_Reachability(const Continuous_Reachability_Problem_Description & problem_description)
{
	type_of_dynamics = problem_description.type_of_dynamics;

	switch(problem_description.type_of_dynamics)
	{
	case LINEAR_TIME_INVARIANT:
		pDynamics = new Linear_Time_Invariant_Dynamics(problem_description.rm_dyn_A, problem_description.utm_dyn_B);
		break;

	case LINEAR_TIME_VARYING:
		pDynamics = new Linear_Time_Varying_Dynamics(problem_description.upm_dyn_A, problem_description.upm_dyn_B, problem_description.upm_dyn_tv);
		break;

	case DETERMINISTIC_DYN:
		pDynamics = new Deterministic_Continuous_Dynamics(problem_description.deterministic_dynamics);
		break;

	case NONDETERMINISTIC_DYN:
		pDynamics = new Nondeterministic_Continuous_Dynamics(problem_description.nondeterministic_dynamics);
		break;
	}

	bSymbolicRemainder = problem_description.bSymbolicRemainder;

	p_tm_setting = new Taylor_Model_Computation_Setting(problem_description.stateVars);

	if(bSymbolicRemainder)
	{
		p_tm_setting->queue_size = problem_description.queue_size;
	}

	// initialize the TM computation setting
	p_tm_setting->initializeAdaptiveSettings(problem_description.step_min, problem_description.step_max, problem_description.order_min, problem_description.order_max);

	if(problem_description.order_max > 0)
	{
		p_tm_setting->setStepsize(problem_description.step_max, problem_description.order_max);
		p_tm_setting->order = problem_description.order_min;
	}
	else
	{
		p_tm_setting->setStepsize(problem_description.step_max, problem_description.order_min);
	}

	p_tm_setting->setCutoff(problem_description.cutoff_threshold);
	p_tm_setting->setRemainderEstimation(problem_description.remainder_estimation);

	p_g_setting = new Global_Computation_Setting();

	// initialize the global computation setting
	p_g_setting->prepareForReachability(problem_description.order_min > problem_description.order_max ? problem_description.order_min : problem_description.order_max);

	p_p_setting = new Plot_Setting(problem_description.plot_setting);

	time			= problem_description.time;
	bSafetyChecking	= problem_description.bSafetyChecking;
	bPlot			= problem_description.bPlot;
	bTMOutput		= problem_description.bTMOutput;
	bPrint			= problem_description.bPrint;

	unsafeSet		= problem_description.unsafeSet;
	initialSets		= problem_description.initialSets;
	fileName		= problem_description.fileName;
}

Continuous_Reachability::~Continuous_Reachability()
{
	delete pDynamics;
	delete p_tm_setting;
	delete p_g_setting;
	delete p_p_setting;
}

void Continuous_Reachability::setup(const Continuous_Reachability_Problem_Description & problem_description)
{
	p_tm_setting = new Taylor_Model_Computation_Setting(problem_description.stateVars);

	// initialize the TM computation setting
	p_tm_setting->setCutoff(problem_description.cutoff_threshold);

	unsigned int order = problem_description.order_min > problem_description.order_max ? problem_description.order_min : problem_description.order_max;
	p_tm_setting->order = order;

	p_g_setting = new Global_Computation_Setting();

	// initialize the global computation setting
	p_g_setting->prepareForReachability(order);

	p_p_setting = new Plot_Setting(problem_description.plot_setting);

	bSafetyChecking	= problem_description.bSafetyChecking;
	bPlot			= problem_description.bPlot;
	bTMOutput		= problem_description.bTMOutput;

	unsafeSet		= problem_description.unsafeSet;
	fileName		= problem_description.fileName;
}

unsigned long Continuous_Reachability::run()
{
	int result;

	if(type_of_dynamics == LINEAR_TIME_INVARIANT)
	{
		result = pDynamics->reach_LTI(linear_flowpipes, result_of_reachability.orders_of_flowpipes, result_of_reachability.safety_of_flowpipes,
				result_of_reachability.num_of_flowpipes, time, initialSets, *p_tm_setting, *p_g_setting, bPrint, unsafeSet, bSafetyChecking, bPlot, bTMOutput);
	}
	else if(type_of_dynamics == LINEAR_TIME_VARYING)
	{
		result = pDynamics->reach_LTV(linear_flowpipes, result_of_reachability.orders_of_flowpipes, result_of_reachability.safety_of_flowpipes,
				result_of_reachability.num_of_flowpipes, time, initialSets, *p_tm_setting, *p_g_setting, bPrint, unsafeSet, bSafetyChecking, bPlot, bTMOutput);
	}
	else
	{
		if(bSymbolicRemainder)
		{
			if(p_tm_setting->bAdaptiveStepSize())
			{
				result = pDynamics->reach_symbolic_remainder_adaptive_stepsize(result_of_reachability.nonlinear_flowpipes, result_of_reachability.orders_of_flowpipes, result_of_reachability.safety_of_flowpipes,
						result_of_reachability.num_of_flowpipes, time, initialSets, *p_tm_setting, *p_g_setting, bPrint, unsafeSet, bSafetyChecking, bPlot, bTMOutput);
			}
			else if(p_tm_setting->bAdaptiveOrder())
			{
				result = pDynamics->reach_symbolic_remainder_adaptive_order(result_of_reachability.nonlinear_flowpipes, result_of_reachability.orders_of_flowpipes, result_of_reachability.safety_of_flowpipes,
						result_of_reachability.num_of_flowpipes, time, initialSets, *p_tm_setting, *p_g_setting, bPrint, unsafeSet, bSafetyChecking, bPlot, bTMOutput);
			}
			else
			{
				result = pDynamics->reach_symbolic_remainder(result_of_reachability.nonlinear_flowpipes, result_of_reachability.orders_of_flowpipes, result_of_reachability.safety_of_flowpipes,
						result_of_reachability.num_of_flowpipes, time, initialSets, *p_tm_setting, *p_g_setting, bPrint, unsafeSet, bSafetyChecking, bPlot, bTMOutput);
			}
		}
		else
		{
			if(p_tm_setting->bAdaptiveStepSize())
			{
				result = pDynamics->reach_adaptive_stepsize(result_of_reachability.nonlinear_flowpipes, result_of_reachability.orders_of_flowpipes, result_of_reachability.safety_of_flowpipes,
						result_of_reachability.num_of_flowpipes, time, initialSets, *p_tm_setting, *p_g_setting, bPrint, unsafeSet, bSafetyChecking, bPlot, bTMOutput);
			}
			else if(p_tm_setting->bAdaptiveOrder())
			{
				result = pDynamics->reach_adaptive_order(result_of_reachability.nonlinear_flowpipes, result_of_reachability.orders_of_flowpipes, result_of_reachability.safety_of_flowpipes,
						result_of_reachability.num_of_flowpipes, time, initialSets, *p_tm_setting, *p_g_setting, bPrint, unsafeSet, bSafetyChecking, bPlot, bTMOutput);
			}
			else
			{
				result = pDynamics->reach(result_of_reachability.nonlinear_flowpipes, result_of_reachability.orders_of_flowpipes, result_of_reachability.safety_of_flowpipes,
						result_of_reachability.num_of_flowpipes, time, initialSets, *p_tm_setting, *p_g_setting, bPrint, unsafeSet, bSafetyChecking, bPlot, bTMOutput);
			}
		}
	}
/*
	Flowpipe fpTmp = result_of_reachability.nonlinear_flowpipes.back();
	result_of_reachability.fp_end_of_time = fpTmp;

	fpTmp.tmvPre.evaluate_time(result_of_reachability.fp_end_of_time.tmvPre, p_tm_setting->step_end_exp_table);

	std::vector<Interval> range;
	result_of_reachability.fp_end_of_time.intEvalNormal(range, p_tm_setting->step_exp_table, p_tm_setting->order_min, p_tm_setting->cutoff_threshold);
	std::cout << range[3].width() << std::endl;
*/
	return result;
}

int Continuous_Reachability::safetyChecking()
{
	if(unsafeSet.size() == 0)
	{
		return UNSAFE;	// since the whole state space is unsafe, the system is not safe
	}

	int checking_result = SAFE;

	std::list<TaylorModelVec<Real> > unsafe_tmv_flowpipes;
	std::list<std::vector<Interval> > unsafe_flowpipe_domains;

	std::list<TaylorModelVec<Real> > unknown_tmv_flowpipes;
	std::list<std::vector<Interval> > unknown_flowpipe_domains;

	result_of_reachability.safety_of_flowpipes.clear();

	std::list<TaylorModelVec<Real> >::const_iterator tmvIter = result_of_reachability.tmv_flowpipes.begin();
	std::list<Flowpipe>::const_iterator fpIter = result_of_reachability.nonlinear_flowpipes.begin();

	unsigned int prog = 0, total_size = result_of_reachability.tmv_flowpipes.size();

	for(; tmvIter != result_of_reachability.tmv_flowpipes.end(); ++tmvIter, ++fpIter)
	{
		int safety = flowstar::safetyChecking(*tmvIter, fpIter->domain, unsafeSet, *p_tm_setting, *p_g_setting);

		if(safety == UNSAFE)
		{
			result_of_reachability.safety_of_flowpipes.push_back(UNSAFE);
			checking_result = UNSAFE;

			if(bTMOutput)
			{
				unsafe_tmv_flowpipes.push_back(*tmvIter);
				unsafe_flowpipe_domains.push_back(fpIter->domain);
			}

			break;
		}
		else if(safety == UNKNOWN)
		{
			result_of_reachability.safety_of_flowpipes.push_back(UNKNOWN);

			if(checking_result == SAFE)
			{
				checking_result = UNKNOWN;
			}

			if(bTMOutput)
			{
				unknown_tmv_flowpipes.push_back(*tmvIter);
				unknown_flowpipe_domains.push_back(fpIter->domain);
			}
		}
		else
		{
			result_of_reachability.safety_of_flowpipes.push_back(SAFE);
		}

		++prog;
		printf("\b\b\b\b");
		printf(BOLD_FONT "%%" RESET_COLOR);
		printf(BOLD_FONT "%3d" RESET_COLOR, (int)(prog*100/total_size));
		fflush(stdout);
	}

	printf("\b\b\b\b");
	printf(BOLD_FONT "%%100\n" RESET_COLOR);
	fflush(stdout);

	if(bTMOutput)
	{
		int mkres = mkdir(counterexampleDir.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
		if(mkres < 0 && errno != EEXIST)
		{
			printf("Can not create the directory for counterexamples.\n");
			return checking_result;
		}

		std::ofstream os_counterexamples(counterexampleDir + fileName + str_counterexample_dumping_name_suffix, std::ofstream::out);

		os_counterexamples << "Unsafe flowpipes:\n\n";

		tmvIter = unsafe_tmv_flowpipes.begin();
		std::list<std::vector<Interval> >::const_iterator doIter = unsafe_flowpipe_domains.begin();

		for(; tmvIter!=unsafe_tmv_flowpipes.end(); ++tmvIter, ++doIter)
		{
			os_counterexamples << "{\n";

			tmvIter->output(os_counterexamples, stateVars, tmVars);

			for(int i=0; i<doIter->size(); ++i)
			{
				os_counterexamples << tmVars.varNames[i] << " in " << (*doIter)[i] << "\n";
			}

			os_counterexamples << "}\n\n\n";
		}

		os_counterexamples << "Unknown flowpipes:\n\n";

		tmvIter = unknown_tmv_flowpipes.begin();
		doIter = unknown_flowpipe_domains.begin();

		for(; tmvIter!=unknown_tmv_flowpipes.end(); ++tmvIter, ++doIter)
		{
			os_counterexamples << "{\n";

			tmvIter->output(os_counterexamples, stateVars, tmVars);

			for(int i=0; i<doIter->size(); ++i)
			{
				os_counterexamples << tmVars.varNames[i] << " in " << (*doIter)[i] << "\n";
			}

			os_counterexamples << "}\n\n\n";
		}

		os_counterexamples.close();
	}

	return checking_result;
}

void Continuous_Reachability::prepareForPlotting()
{
	result_of_reachability.tmv_flowpipes.clear();

	if(bPrint)
	{
		printf("Preparing for plotting...\n");
	}

	if(type_of_dynamics == LINEAR_TIME_INVARIANT || type_of_dynamics == LINEAR_TIME_VARYING)
	{
		Interval intStep(0, p_tm_setting->step_max), intUnit(-1,1);

		unsigned int prog = 0, total_size = initialSets.size() * linear_flowpipes.size();
		int rangeDim = initialSets[0].tmvPre.tms.size();

		for(unsigned int m=0; m<initialSets.size(); ++m)
		{
			std::vector<Interval> newDomain = initialSets[0].domain;
			newDomain[0] = intStep;

			std::vector<Interval> polyRangeX0;
			initialSets[m].tmvPre.polyRange(polyRangeX0, initialSets[m].domain);

			std::vector<Interval> range_of_X0(rangeDim);
			for(int k=0; k<rangeDim; ++k)
			{
				range_of_X0[k] = polyRangeX0[k] + initialSets[m].tmvPre.tms[k].remainder;
			}

			std::list<LinearFlowpipe>::iterator iter;

			for(iter = linear_flowpipes.begin(); iter != linear_flowpipes.end(); ++iter)
			{
				TaylorModelVec<Real> tmvFlowpipe;

				iter->evaluate(tmvFlowpipe, p_p_setting->outputDims, initialSets[m].tmvPre, polyRangeX0, range_of_X0, newDomain, *p_tm_setting);

				result_of_reachability.tmv_flowpipes.push_back(tmvFlowpipe);

				Flowpipe flowpipe;
				flowpipe.domain = newDomain;
				result_of_reachability.nonlinear_flowpipes.push_back(flowpipe);

				if(bPrint)
				{
					++prog;
					printf("\b\b\b");
					printf(BOLD_FONT "%%" RESET_COLOR);
					printf(BOLD_FONT "%2d" RESET_COLOR, (int)(prog*100/total_size));
					fflush(stdout);
				}
			}
		}

		linear_flowpipes.clear();

		if(bPrint)
		{
			printf("\n");
		}
	}
	else
	{
		unsigned int prog = 0, total_size = result_of_reachability.nonlinear_flowpipes.size();

		std::list<Flowpipe>::const_iterator fpIter = result_of_reachability.nonlinear_flowpipes.begin();
		std::list<unsigned int>::const_iterator orderIter = result_of_reachability.orders_of_flowpipes.begin();

		for(; fpIter != result_of_reachability.nonlinear_flowpipes.end(); ++fpIter, ++orderIter)
		{
			TaylorModelVec<Real> tmvTmp;

			fpIter->compose(tmvTmp, p_p_setting->outputDims, *orderIter, p_tm_setting->cutoff_threshold);

			result_of_reachability.tmv_flowpipes.push_back(tmvTmp);

			if(bPrint)
			{
				++prog;
				printf("\b\b\b");
				printf(BOLD_FONT "%%" RESET_COLOR);
				printf(BOLD_FONT "%2d" RESET_COLOR, (int)(prog*100/total_size));
				fflush(stdout);
			}
		}

		if(bPrint)
		{
			printf("\n");
		}
	}

	if(bPrint)
	{
		printf("Done.\n");
	}
}

void Continuous_Reachability::prepareForTMOutput()
{
	result_of_reachability.tmv_flowpipes.clear();

	if(type_of_dynamics == LINEAR_TIME_INVARIANT || type_of_dynamics == LINEAR_TIME_VARYING)
	{
		Interval intStep(0, p_tm_setting->step_max), intUnit(-1,1);

		int i = 0, total_size = initialSets.size() * linear_flowpipes.size();
		int rangeDim = initialSets[0].tmvPre.tms.size();

		for(unsigned int m=0; m<initialSets.size(); ++m)
		{
			std::vector<Interval> newDomain = initialSets[0].domain;
			newDomain[0] = intStep;

			std::vector<Interval> polyRangeX0;
			initialSets[m].tmvPre.polyRange(polyRangeX0, initialSets[m].domain);

			std::vector<Interval> range_of_X0(rangeDim);
			for(int k=0; k<rangeDim; ++k)
			{
				range_of_X0[k] = polyRangeX0[k] + initialSets[m].tmvPre.tms[k].remainder;
			}

			std::list<LinearFlowpipe>::iterator iter;

			for(iter = linear_flowpipes.begin(); iter != linear_flowpipes.end(); ++iter)
			{
				TaylorModelVec<Real> tmvFlowpipe;

				iter->evaluate(tmvFlowpipe, initialSets[m].tmvPre, polyRangeX0, range_of_X0, newDomain, *p_tm_setting);

				result_of_reachability.tmv_flowpipes.push_back(tmvFlowpipe);

				Flowpipe flowpipe;
				flowpipe.domain = newDomain;
				result_of_reachability.nonlinear_flowpipes.push_back(flowpipe);

				if(bPrint)
				{
					++i;
					printf("\b\b\b");
					printf(BOLD_FONT "%%" RESET_COLOR);
					printf(BOLD_FONT "%2d" RESET_COLOR, (int)(i*100/total_size));
					fflush(stdout);
				}
			}
		}

		linear_flowpipes.clear();

		if(bPrint)
		{
			printf("\n");
		}
	}
	else
	{
		result_of_reachability.transformToTaylorModels(*p_tm_setting, bPrint);
	}
}

void Continuous_Reachability::plot_2D() const
{
	p_p_setting->plot_2D(fileName, result_of_reachability);
}

void Continuous_Reachability::tmOutput(std::ostream & os) const
{
	if(bPrint)
	{
		printf("Writing the flowpipe(s)...\n");
	}

	os << "state var ";

	for(unsigned int i=0; i<stateVars.varNames.size() - 1; ++i)
	{
		os << stateVars.varNames[i] << ", ";
	}

	os << stateVars.varNames.back() << "\n\n";


	switch(p_p_setting->type_of_file)
	{
	case PLOT_GNUPLOT:
		switch(p_p_setting->type_of_object)
		{
		case PLOT_INTERVAL:
			os << "gnuplot interval " << stateVars.varNames[p_p_setting->outputDims[0]] << " , " << stateVars.varNames[p_p_setting->outputDims[1]] << "\n\n";
			break;
		case PLOT_OCTAGON:
			os << "gnuplot octagon " << stateVars.varNames[p_p_setting->outputDims[0]] << " , " << stateVars.varNames[p_p_setting->outputDims[1]] << "\n\n";
			break;
		case PLOT_GRID:
			os << "gnuplot grid " << p_p_setting->num_of_pieces << " " << stateVars.varNames[p_p_setting->outputDims[0]] << " , " << stateVars.varNames[p_p_setting->outputDims[1]] << "\n\n";
			break;
		}
		break;
	case PLOT_MATLAB:
		switch(p_p_setting->type_of_object)
		{
		case PLOT_INTERVAL:
			os << "matlab interval " << stateVars.varNames[p_p_setting->outputDims[0]] << " , " << stateVars.varNames[p_p_setting->outputDims[1]] << "\n\n";
			break;
		case PLOT_OCTAGON:
			os << "matlab octagon " << stateVars.varNames[p_p_setting->outputDims[0]] << " , " << stateVars.varNames[p_p_setting->outputDims[1]] << "\n\n";
			break;
		case PLOT_GRID:
			os << "matlab grid " << p_p_setting->num_of_pieces << " " << stateVars.varNames[p_p_setting->outputDims[0]] << " , " << stateVars.varNames[p_p_setting->outputDims[1]] << "\n\n";
			break;
		}
		break;
	}

	if(p_tm_setting->order_max == 0)
	{
		os << "order " << p_tm_setting->order_min << "\n\n";
	}
	else
	{
		os << "order " << p_tm_setting->order_max << "\n\n";
	}

	os << "cutoff " << p_tm_setting->cutoff_threshold.sup() << "\n\n";
	os << "output " << fileName << "\n\n";

	if(bSafetyChecking)
	{
		// output the unsafe set
		os << "unsafe\n{\n";

		for(int i=0; i<unsafeSet.size(); ++i)
		{
			unsafeSet[i].output(os, stateVars);
		}

		os << "}\n\n";
	}

	os << "continuous flowpipes\n{\n";

	os << "tm var ";

	for(int i=0; i<tmVars.varNames.size()-1; ++i)
	{
		os << tmVars.varNames[i] << ", ";
	}

	os << tmVars.varNames.back() << "\n\n";

	std::list<TaylorModelVec<Real> > unsafe_tmv_flowpipes;
	std::list<std::vector<Interval> > unsafe_flowpipe_domains;

	std::list<TaylorModelVec<Real> > unknown_tmv_flowpipes;
	std::list<std::vector<Interval> > unknown_flowpipe_domains;

	std::list<TaylorModelVec<Real> >::const_iterator tmvIter = result_of_reachability.tmv_flowpipes.begin();
	std::list<Flowpipe>::const_iterator fpIter = result_of_reachability.nonlinear_flowpipes.begin();
	std::list<int>::const_iterator safetyIter = result_of_reachability.safety_of_flowpipes.begin();

	for(; safetyIter != result_of_reachability.safety_of_flowpipes.end(); ++tmvIter, ++fpIter, ++safetyIter)
	{
		if(*safetyIter == UNSAFE)
		{
			unsafe_tmv_flowpipes.push_back(*tmvIter);
			unsafe_flowpipe_domains.push_back(fpIter->domain);
		}
		else if(*safetyIter == UNKNOWN)
		{
			unknown_tmv_flowpipes.push_back(*tmvIter);
			unknown_flowpipe_domains.push_back(fpIter->domain);
		}

		os << "{\n";
		tmvIter->output(os, stateVars, tmVars);

		for(int i=0; i<fpIter->domain.size(); ++i)
		{
			os << tmVars.varNames[i] << " in " << fpIter->domain[i] << "\n";
		}

		os << "}\n\n";
	}

	os << "}\n";

	bool bDumpCounterexamples = true;

	int mkres = mkdir(counterexampleDir.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
	if(mkres < 0 && errno != EEXIST)
	{
		printf("Can not create the directory for counterexamples.\n");
		bDumpCounterexamples = false;
	}

	if(bDumpCounterexamples)
	{
		std::ofstream os_counterexamples(counterexampleDir + fileName + str_counterexample_dumping_name_suffix, std::ofstream::out);

		os_counterexamples << "Unsafe flowpipes:\n\n";

		tmvIter = unsafe_tmv_flowpipes.begin();
		std::list<std::vector<Interval> >::const_iterator doIter = unsafe_flowpipe_domains.begin();

		for(; tmvIter!=unsafe_tmv_flowpipes.end(); ++tmvIter, ++doIter)
		{
			os_counterexamples << "{\n";

			tmvIter->output(os_counterexamples, stateVars, tmVars);

			for(int i=0; i<doIter->size(); ++i)
			{
				os_counterexamples << tmVars.varNames[i] << " in " << (*doIter)[i] << "\n";
			}

			os_counterexamples << "}\n\n\n";
		}

		os_counterexamples << "Unknown flowpipes:\n\n";

		tmvIter = unknown_tmv_flowpipes.begin();
		doIter = unknown_flowpipe_domains.begin();

		for(; tmvIter!=unknown_tmv_flowpipes.end(); ++tmvIter, ++doIter)
		{
			os_counterexamples << "{\n";

			tmvIter->output(os_counterexamples, stateVars, tmVars);

			for(int i=0; i<doIter->size(); ++i)
			{
				os_counterexamples << tmVars.varNames[i] << " in " << (*doIter)[i] << "\n";
			}

			os_counterexamples << "}\n\n\n";
		}

		os_counterexamples.close();
	}

	if(bPrint)
	{
		printf("Done.\n");
	}
}














namespace flowstar
{

int safetyChecking(const TaylorModelVec<Real> & tmv, const std::vector<Interval> & domain, const std::vector<Constraint> & unsafeSet, const Taylor_Model_Computation_Setting & tm_setting, const Global_Computation_Setting & g_setting)
{
	if(unsafeSet.size() == 0)
	{
		return SAFE;
	}

	unsigned int rangeDim = tmv.tms.size();
	int result = UNKNOWN;
	bool bContained = true;

	std::vector<Interval> tmvRange;
	tmv.intEval(tmvRange, domain);

	for(unsigned int i=0; i<unsafeSet.size(); ++i)
	{
		Interval I;

		// interval evaluation on the constraint
		unsafeSet[i].expression.evaluate(I, tmvRange);

		if(unsafeSet[i].bound < I.inf())
		{
			// no intersection with the unsafe set
			result = SAFE;
			break;
		}
		else
		{
			if(!(unsafeSet[i].bound >= I.sup()) && bContained)
			{
				bContained = false;
			}
		}
	}

	if(result == UNKNOWN)
	{
		if(bContained)
		{
			return UNSAFE;
		}
		else
		{
			// do a simple branch & bound for safety checking
			std::vector<HornerForm<Real> > obj_hfs;
			std::vector<Interval> obj_rems;

			result = SAFE;

			for(unsigned int i=0; i<unsafeSet.size(); ++i)
			{
				TaylorModel<Real> tmTmp;

				// interval evaluation on the constraint
				unsafeSet[i].expression.evaluate(tmTmp, tmv.tms, tm_setting.order, domain, tm_setting.cutoff_threshold, g_setting);

				HornerForm<Real> obj_hf;
				tmTmp.expansion.toHornerForm(obj_hf);
				obj_hfs.push_back(obj_hf);
				obj_rems.push_back(tmTmp.remainder);
			}

			std::vector<Interval> refined_domain = domain;

			std::list<Interval> subdivisions;

			if(domain[0].width() > REFINEMENT_PREC)
			{
				subdivisions.push_back(domain[0]);
			}

			for(; subdivisions.size() > 0; )
			{
				Interval subdivision = subdivisions.front();
				subdivisions.pop_front();

				int result_iter = UNKNOWN;
				bool bContained_iter = true;

				refined_domain[0] = subdivision;

				for(int i=0; i<unsafeSet.size(); ++i)
				{
					Interval I;
					obj_hfs[i].evaluate(I, refined_domain);

					I += obj_rems[i];

					if(unsafeSet[i].bound < I.inf())
					{
						// no intersection with the unsafe set
						result_iter = SAFE;
						break;
					}
					else
					{
						if(!(unsafeSet[i].bound >= I.sup()) && bContained_iter)
						{
							bContained_iter = false;
						}
					}
				}

				if(result_iter == UNKNOWN)
				{
					if(bContained_iter)
					{
						return UNSAFE;
					}
					else
					{
						if(subdivision.width() <= REFINEMENT_PREC)
						{
							return UNKNOWN;
						}

						// split the domain
						Interval I1, I2;
						subdivision.split(I1, I2);

						if(I1.width() <= REFINEMENT_PREC)
						{
							if(result == SAFE)
								result = UNKNOWN;
						}
						else
						{
							subdivisions.push_back(I1);
						}

						if(I2.width() <= REFINEMENT_PREC)
						{
							if(result == SAFE)
								result = UNKNOWN;
						}
						else
						{
							subdivisions.push_back(I2);
						}
					}
				}
			}

			return result;
		}
	}
	else
	{
		return SAFE;
	}
}

void gridBox(std::list<std::vector<Interval> > & grids, const std::vector<Interval> & box, const unsigned int num)
{
	grids.clear();
	grids.push_back(box);

	for(unsigned int i=0; i<box.size(); ++i)
	{
		std::list<std::vector<Interval> >::iterator gridIter;
		std::list<std::vector<Interval> > newGrids;

		for(; grids.size() > 0;)
		{
			gridIter = grids.begin();

			std::list<Interval> queue;
			(*gridIter)[i].split(queue, num);

			std::list<Interval>::iterator iterComponent = queue.begin();
			for(; iterComponent != queue.end(); ++iterComponent)
			{
				std::vector<Interval> tmpBox = *gridIter;
				tmpBox[i] = *iterComponent;
				newGrids.push_back(tmpBox);
			}

			grids.pop_front();
		}

		grids = newGrids;
	}
}




int contract_remainder(const std::vector<Interval> & polyRange, std::vector<Interval> & remainders, const std::vector<Constraint> & constraints)
{
	if(constraints.size() == 0)
	{
		return 0;
	}

	bool bvalid = true;
	bool bcontinue = true;

	unsigned int counter = 0;
	unsigned int num = constraints.size();
	std::vector<bool> bNec(num, true);

	unsigned int rangeDim = polyRange.size();
	unsigned int domainDim = rangeDim + 1;

	std::vector<Interval> intVecTemp = polyRange;


	// 1: we check the intersection with every constraint
	for(unsigned int i=0; i<rangeDim; ++i)
	{
		intVecTemp[i] = polyRange[i] + remainders[i];
	}

	for(unsigned int i=0; i<constraints.size(); ++i)
	{
		Interval intTemp;

		constraints[i].expression.evaluate(intTemp, intVecTemp);

		if(constraints[i].bound < intTemp.inf())
		{
			// no intersection on the left half
			bvalid = false;
			break;
		}
		else if(constraints[i].bound > intTemp.sup())
		{
			// do not need to apply domain contraction w.r.t. the current constraint
			bNec[i] = false;
			++counter;
		}
		else
		{
			bNec[i] = true;
			continue;
		}
	}

	if(!bvalid)
	{
		return -1;	// no intersection is detected
	}
	else if(counter == num)
	{
		return 0;	// no need to do contraction
	}

	// 2: contract the remainder
	for(; bcontinue; )
	{
		std::vector<Interval> oldRemainders = remainders;

		for(int i=0; i<rangeDim; ++i)
		{
			Interval newInt = remainders[i];
			std::vector<bool> localNec = bNec;
			int localCounter = counter;

			for(int k=0; k<rangeDim; ++k)
			{
				if(k != i)
				{
					intVecTemp[k] = polyRange[k] + remainders[k];
				}
				else
				{
					intVecTemp[k] = polyRange[k];
				}
			}

			double w = newInt.width();

			// search an approximation for the lower bound
			for(; w > DC_THRESHOLD_SEARCH;)
			{
				Interval intLeft;
				Interval intRight;
				newInt.split(intLeft, intRight);

				for(int j=0; j<constraints.size(); ++j)
				{
					if(localNec[j])
					{
						Interval intTemp;
						std::vector<Interval> newIntVecTemp = intVecTemp;
						newIntVecTemp[i] += intLeft;

						constraints[j].expression.evaluate(intTemp, newIntVecTemp);

						if(constraints[j].bound < intTemp.inf())
						{
							// no intersection on the left half
							newInt = intRight;
							w = newInt.width();
							break;
						}
						else if(constraints[j].bound >= intTemp.sup())
						{
							// do not need to apply domain contraction w.r.t. the current constraint
							newInt = intLeft;
							w = newInt.width();
							localNec[j] = false;
							++localCounter;
						}
						else
						{
							// refine the interval
							newInt = intLeft;
							w = newInt.width();

							continue;
						}
					}
				}

				if(localCounter == constraints.size())
				{
					break;
				}
			}

			// set the lower bound
			remainders[i].setInf(newInt.inf());

			newInt = remainders[i];
			w = newInt.width();

			localNec = bNec;
			localCounter = counter;

			// search an approximation for the upper bound
			for(; w > DC_THRESHOLD_SEARCH;)
			{
				Interval intLeft;
				Interval intRight;
				newInt.split(intLeft, intRight);

				for(int j=0; j<constraints.size(); ++j)
				{
					if(localNec[j])
					{
						Interval intTemp;
						std::vector<Interval> newIntVecTemp = intVecTemp;
						newIntVecTemp[i] += intRight;

						constraints[j].expression.evaluate(intTemp, newIntVecTemp);

						if(constraints[j].bound < intTemp.inf())
						{
							// no intersection on the right half
							newInt = intLeft;
							w = newInt.width();
							break;
						}
						else if(constraints[j].bound > intTemp.sup())
						{
							// do not need to apply domain contraction w.r.t. the current constraint
							newInt = intRight;
							w = newInt.width();
							localNec[j] = false;
							++localCounter;
						}
						else
						{
							// refine the interval
							newInt = intRight;
							w = newInt.width();
							continue;
						}
					}
				}

				if(localCounter == constraints.size())
				{
					break;
				}
			}

			remainders[i].setSup(newInt.sup());	// set the upper bound

			if(!remainders[i].valid())
			{
				bvalid = false;
				break;
			}
		}

		if(!bvalid)
		{
			break;
		}

		bcontinue = false;
		for(int i=0; i<rangeDim; ++i)
		{
			if(oldRemainders[i].widthRatio(remainders[i]) <= DC_THRESHOLD_IMPROV)
			{
				bcontinue = true;
				break;
			}
		}
	}

	if(!bvalid)
	{
		return -1;	// no intersection is detected
	}
	else
	{
		return 1;
	}
}

unsigned int findProperOrder(Real & error, const Real & max, const Real & min, const Real & tolerance, const unsigned int start_order)
{
	unsigned int order = start_order;
	unsigned int inc = ceil((double)(order)/2);

	while(true)
	{
		error = 1 / (1 - max/(order+2));

		if(error >= 0)
		{
			unsigned int k = order + 1;

			Real tmp1 = max;
			tmp1.pow_assign(k);

			Real tmp2;
			tmp2.factorial(k);
			tmp1 /= tmp2;

			error *= tmp1;
			Real tmp3 = min * tolerance;

			if(error < min * tolerance)
			{
				break;
			}
			else
			{
				order += inc;
			}
		}
		else
		{
			order += inc;
		}
	}

	return order;
}

void check_connectivities(Matrix<bool> & result, Matrix<bool> & adjMatrix)
{
	int n = adjMatrix.rows();
	result = adjMatrix;

	std::vector<bool> bselected;
	for(int i=0; i<n; ++i)
	{
		bselected.push_back(false);
	}

	for(int i=0; i<n; ++i)
	{
		// BFS is used to check the connectivity of two nodes
		for(int j=0; j<n; ++j)
		{
			bselected[j] = false;
		}

		std::list<int> unvisited;
		unvisited.push_back(i);

		while(unvisited.size() > 0)
		{
			int j = unvisited.front();
			unvisited.pop_front();

			for(int k=0; k<n; ++k)
			{
				if(bselected[k] == false && adjMatrix[j][k] == true)
				{
					bselected[k] = true;
					unvisited.push_back(k);
					result[i][k] = true;
				}
			}
		}
	}
}

void compute_one_step_trans(Matrix<UnivariateTaylorModel<Real> > & utm_Phi_t, Matrix<Real> & rm_Phi_t, Matrix<UnivariateTaylorModel<Real> > & utm_Psi_t, Matrix<Real> & rm_Psi_t,
		Matrix<Interval> & tv_part, const Matrix<UnivariatePolynomial<Real> > & A_t, const Matrix<UnivariatePolynomial<Real> > & B_t, const Matrix<UnivariatePolynomial<Real> > & tv_t,
		Matrix<bool> & connectivity, const bool bAuto, const UnivariatePolynomial<Real> & up_t, const unsigned int order, std::vector<Real> & step_end_exp_table)
{
	unsigned int rangeDim = A_t.rows(), numTVPar = tv_t.cols();
	double step = step_end_exp_table[1].toDouble();
	Real rStep = step_end_exp_table[1];

	Matrix<Real> identity(rangeDim);
	Matrix<Interval> im_zero_Phi(rangeDim, rangeDim), im_zero_Psi(rangeDim, 1);


	// evaluate a guaranteed remainder interval
	Matrix<UnivariatePolynomial<Real> > local_A_t(rangeDim, rangeDim);
	A_t.substitute(local_A_t, up_t);
	utm_Phi_t = identity;

	Matrix<Interval> im_A_t(rangeDim, rangeDim);
	local_A_t.evaluate(im_A_t, interval_utm_setting.val_exp_table);

	Real A_max = im_A_t.max_norm();

	// find the proper order
	Real A_min = im_A_t.min_entry();

	Real tolerance = APPROX_TOLERANCE;
	Real error;

	unsigned int approx_order = findProperOrder(error, A_max, A_min, tolerance, order);

	// reconstruct step_exp_table
	unsigned int nec_order = 2*approx_order + 1;

	unsigned int currentOrder = step_end_exp_table.size();

	if(currentOrder < nec_order)
	{
		Real tmp = step_end_exp_table.back();

		for(unsigned int i=currentOrder; i<=nec_order; ++i)
		{
			tmp *= rStep;
			step_end_exp_table.push_back(tmp);
		}
	}


	Interval intStep(0, step);
	interval_utm_setting.resetOrder(intStep, approx_order);

	Matrix<UnivariateTaylorModel<Real> > utm_tmp_Psi;
	Matrix<UnivariatePolynomial<Real> > local_B_t(rangeDim, 1);
	if(!bAuto)
	{
		B_t.substitute(local_B_t, up_t);
		utm_Psi_t = local_B_t;
		utm_Psi_t.integral(intStep);
		utm_tmp_Psi = utm_Psi_t;
	}

	Matrix<UnivariateTaylorModel<Real> > utm_tmp_Phi = utm_Phi_t;


	// compute the polynomial approximations
	for(int i=1; i<=order; ++i)
	{
		utm_tmp_Phi = local_A_t * utm_tmp_Phi;
		utm_tmp_Phi.integral(intStep);

		if(i < order)
		{
			utm_tmp_Phi.ctrunc(order);
		}

		utm_Phi_t += utm_tmp_Phi;

		if(!bAuto)
		{
			utm_tmp_Psi = local_A_t * utm_tmp_Psi;
			utm_tmp_Psi.integral(intStep);
			if(i < order)
			{
				utm_tmp_Psi.ctrunc(order);
			}

			utm_Psi_t += utm_tmp_Psi;
		}
	}

	Interval intErr;
	error.to_sym_int(intErr);

	Matrix<Interval> im_error(rangeDim, rangeDim);
	for(int i=0; i<rangeDim; ++i)
	{
		for(int j=0; j<rangeDim; ++j)
		{
			if(connectivity[i][j])
			{
				im_error[i][j] = intErr;
			}
		}
	}

	utm_Phi_t.addRemainder(im_error);
	utm_Phi_t.evaluate(rm_Phi_t, step_end_exp_table);

	if(!bAuto)
	{
		Matrix<Interval> im_B_t(rangeDim, 1);
		local_B_t.evaluate(im_B_t);
		im_error *= im_B_t;
		im_error *= intStep;

		utm_Psi_t.addRemainder(im_error);
		utm_Psi_t.evaluate(rm_Psi_t, step_end_exp_table);
	}

	utm_Phi_t.ctrunc(order);
	utm_Psi_t.ctrunc(order);

	if(numTVPar > 0)
	{
		Matrix<UnivariatePolynomial<Real> > local_tv_t(rangeDim, numTVPar);
		tv_t.substitute(local_tv_t, up_t);
		local_tv_t.evaluate(tv_part, interval_utm_setting.val_exp_table);
	}
}

}
















