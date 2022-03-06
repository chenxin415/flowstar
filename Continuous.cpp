/*---
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.
  The code is released as is under the GNU General Public License (GPL).
---*/

#include "Continuous.h"

using namespace flowstar;



Symbolic_Remainder::Symbolic_Remainder()
{
	max_size = 0;
}

Symbolic_Remainder::Symbolic_Remainder(const Flowpipe & initial_set, const unsigned int s)
{
	scalars.resize(initial_set.tmv.tms.size(), 1);
	initial_set.tmv.Expansion(polynomial_of_initial_set);
	max_size = s;
}

Symbolic_Remainder::Symbolic_Remainder(const Symbolic_Remainder & symbolic_remainder)
{
	J							= symbolic_remainder.J;
	Phi_L						= symbolic_remainder.Phi_L;
	scalars						= symbolic_remainder.scalars;
	polynomial_of_initial_set	= symbolic_remainder.polynomial_of_initial_set;
	max_size					= symbolic_remainder.max_size;
}

Symbolic_Remainder::~Symbolic_Remainder()
{
}

void Symbolic_Remainder::reset(const Flowpipe & initial_set)
{
	scalars.clear();
	scalars.resize(initial_set.tmv.tms.size(), 1);

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
	max_size					= symbolic_remainder.max_size;

	return *this;
}



Computational_Setting::Computational_Setting()
{
	time = 0;
	bPrint = false;

	// default size of the cutoff threshold
	Interval cutoff_threshold(-1e-10,1e-10);
	tm_setting.setCutoff(cutoff_threshold);
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

Flowpipe::Flowpipe(const TaylorModelVec<Real> & tmv_flowpipe, const std::vector<Interval> & flowpipe_domain)
{
	tmvPre = tmv_flowpipe;
	domain = flowpipe_domain;

	Interval cutoff_threshold(-1e-10,1e-10);

	tmvPre.normalize(domain, cutoff_threshold);

	TaylorModelVec<Real> tmvTmp(domain.size() - 1);
	tmv = tmvTmp;
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
	// we first normalize the Taylor model tmv
	tmv.normalize(domain, cutoff_threshold);

	int rangeDim = tmv.tms.size();

	// compute the center point of tmv
	std::vector<Real> const_part;
	tmv.constant(const_part);
	tmv.rmConstant();

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		Real m;
		tmv.tms[i].remainder.remove_midpoint(m);
		const_part[i] += m;
	}

	std::vector<Interval> tmvRange;
	tmv.intEval(tmvRange, domain);

	Matrix<Real> coefficients(rangeDim, rangeDim+1);
	Real zero(0);

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		Real sup;
		tmvRange[i].mag(sup);

		if(sup == zero)
		{
			coefficients[i][i+1] = 0;
		}
		else
		{
			coefficients[i][i+1] = sup;
			tmv.tms[i] /= sup;
		}
	}

	TaylorModelVec<Real> newVars(coefficients);
	for(unsigned int i=0; i<rangeDim; ++i)
	{
		TaylorModel<Real> tmTemp(const_part[i], rangeDim+1);
		newVars.tms[i] += tmTemp;
	}

	std::vector<Interval> polyRange;
	newVars.polyRange(polyRange, domain);

	for(unsigned int i=0; i<tmvPre.tms.size(); ++i)
	{
		TaylorModel<Real> tmTemp;
		tmvPre.tms[i].insert_ctrunc(tmTemp, newVars, polyRange, domain, tmvPre.tms[i].expansion.degree(), cutoff_threshold);
		tmvPre.tms[i] = tmTemp;
	}
}

int Flowpipe::safetyChecking(const std::vector<Constraint> & unsafeSet, const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting) const
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
			if(domain[0].width() <= REFINEMENT_PREC)
				return UNKNOWN;


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

bool Flowpipe::isInTarget(const std::vector<Constraint> & targetSet, const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting) const
{
	unsigned int rangeDim = tmvPre.tms.size();

	std::vector<Interval> tmvRange;
	tmvPre.intEvalNormal(tmvRange, tm_setting.step_exp_table);

	for(unsigned int i=0; i<targetSet.size(); ++i)
	{
		Interval I;

		// interval evaluation on the constraint
		targetSet[i].expression.evaluate(I, tmvRange);

		if(!(targetSet[i].bound >= I.sup()))
		{
			// not entirely contained in the target set
			return false;
		}
	}

	return true;
}

bool Flowpipe::isInTarget(const std::vector<Constraint> & targetSet, const Computational_Setting & setting) const
{
	return isInTarget(targetSet, setting.tm_setting, setting.g_setting);
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

int Flowpipe::advance_deterministic(Flowpipe & result, const std::vector<Expression<Real> > & ode, const Taylor_Model_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Setting & g_setting) const
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

		int res = remainder_contraction_int(intVecTmp, contracted_remainders, invariant);

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

int Flowpipe::advance_nondeterministic(Flowpipe & result, const std::vector<Expression<Interval> > & ode, const Taylor_Model_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Setting & g_setting) const
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

		int res = remainder_contraction_int(intVecTmp, contracted_remainders, invariant);

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

int Flowpipe::advance_deterministic_adaptive_stepsize(Flowpipe & result, const std::vector<Expression<Real> > & ode, const double new_stepsize, Taylor_Model_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Setting & g_setting) const
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

		int res = remainder_contraction_int(intVecTmp, contracted_remainders, invariant);

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


	if(new_stepsize > 0)
	{
		tm_setting.setStepsize(new_stepsize, tm_setting.order);
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

int Flowpipe::advance_nondeterministic_adaptive_stepsize(Flowpipe & result, const std::vector<Expression<Interval> > & ode, const double new_stepsize, Taylor_Model_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Setting & g_setting) const
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

		int res = remainder_contraction_int(intVecTmp, contracted_remainders, invariant);

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

	if(new_stepsize > 0)
	{
		tm_setting.setStepsize(new_stepsize, tm_setting.order);
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

int Flowpipe::advance_deterministic_adaptive_order(Flowpipe & result, const std::vector<Expression<Real> > & ode, Taylor_Model_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Setting & g_setting) const
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

		int res = remainder_contraction_int(intVecTmp, contracted_remainders, invariant);

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

int Flowpipe::advance_nondeterministic_adaptive_order(Flowpipe & result, const std::vector<Expression<Interval> > & ode, Taylor_Model_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Setting & g_setting) const
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

		int res = remainder_contraction_int(intVecTmp, contracted_remainders, invariant);

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






int Flowpipe::advance_deterministic(Flowpipe & result, const std::vector<Expression<Real> > & ode, const Taylor_Model_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Setting & g_setting, Symbolic_Remainder & symbolic_remainder) const
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
/*
		Matrix<Interval> x0_rem(rangeDim, 1);
		tmv_of_x0.Remainder(x0_rem);
		J_ip1 += x0_rem;
*/
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

			int res = remainder_contraction_int(intVecTmp, contracted_remainders, invariant);

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
//			range_of_x0[i] = intUnit;
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

int Flowpipe::advance_nondeterministic(Flowpipe & result, const std::vector<Expression<Interval> > & ode, const Taylor_Model_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Setting & g_setting, Symbolic_Remainder & symbolic_remainder) const
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
/*
		Matrix<Interval> x0_rem(rangeDim, 1);
		tmv_of_x0.Remainder(x0_rem);
		J_ip1 += x0_rem;
*/
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

			int res = remainder_contraction_int(intVecTmp, contracted_remainders, invariant);

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
//			range_of_x0[i] = intUnit;
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

int Flowpipe::advance_deterministic_adaptive_stepsize(Flowpipe & result, const std::vector<Expression<Real> > & ode, const double new_stepsize, Taylor_Model_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Setting & g_setting, Symbolic_Remainder & symbolic_remainder) const
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
/*
		Matrix<Interval> x0_rem(rangeDim, 1);
		tmv_of_x0.Remainder(x0_rem);
		J_ip1 += x0_rem;
*/
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

			int res = remainder_contraction_int(intVecTmp, contracted_remainders, invariant);

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
//			range_of_x0[i] = intUnit;
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

	if(new_stepsize > 0)
	{
		tm_setting.setStepsize(new_stepsize, tm_setting.order);
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

int Flowpipe::advance_nondeterministic_adaptive_stepsize(Flowpipe & result, const std::vector<Expression<Interval> > & ode, const double new_stepsize, Taylor_Model_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Setting & g_setting, Symbolic_Remainder & symbolic_remainder) const
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
/*
		Matrix<Interval> x0_rem(rangeDim, 1);
		tmv_of_x0.Remainder(x0_rem);
		J_ip1 += x0_rem;
*/
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

			int res = remainder_contraction_int(intVecTmp, contracted_remainders, invariant);

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
//			range_of_x0[i] = intUnit;
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

	if(new_stepsize > 0)
	{
		tm_setting.setStepsize(new_stepsize, tm_setting.order);
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

int Flowpipe::advance_deterministic_adaptive_order(Flowpipe & result, const std::vector<Expression<Real> > & ode, Taylor_Model_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Setting & g_setting, Symbolic_Remainder & symbolic_remainder) const
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
/*
		Matrix<Interval> x0_rem(rangeDim, 1);
		tmv_of_x0.Remainder(x0_rem);
		J_ip1 += x0_rem;
*/
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

			int res = remainder_contraction_int(intVecTmp, contracted_remainders, invariant);

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
//			range_of_x0[i] = intUnit;
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

int Flowpipe::advance_nondeterministic_adaptive_order(Flowpipe & result, const std::vector<Expression<Interval> > & ode, Taylor_Model_Setting & tm_setting, const std::vector<Constraint> & invariant, const Global_Setting & g_setting, Symbolic_Remainder & symbolic_remainder) const
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
/*
		Matrix<Interval> x0_rem(rangeDim, 1);
		tmv_of_x0.Remainder(x0_rem);
		J_ip1 += x0_rem;
*/
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

			int res = remainder_contraction_int(intVecTmp, contracted_remainders, invariant);

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
//			range_of_x0[i] = intUnit;
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

	tv_remainder = flowpipe.tv_remainder;
	remainder_constraints = flowpipe.remainder_constraints;
}

LinearFlowpipe::~LinearFlowpipe()
{
}

int LinearFlowpipe::safetyChecking(const std::vector<Constraint> & unsafeSet, const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting,
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


	Matrix<Interval> im_remainder(rangeDim, 1);

	if(!tv_remainder.isEmpty())
	{
		tv_remainder.intEval(im_remainder);

		for(int i=0; i<rangeDim; ++i)
		{
			range_of_x[i] += im_remainder[i][0];
		}
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
			if(domain[0].width() <= REFINEMENT_PREC)
				return UNKNOWN;

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

void LinearFlowpipe::evaluate(TaylorModelVec<Real> & result, const TaylorModelVec<Real> & tmv_of_X0, const std::vector<Interval> & polyRangeX0, const std::vector<Interval> & range_of_X0, const std::vector<Interval> & domain, const Taylor_Model_Setting & tm_setting)
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
			TaylorModel<Real> tmTmp2(Phi[i][j], rangeDim+1, true);
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

void LinearFlowpipe::evaluate(TaylorModelVec<Real> & result, const TaylorModelVec<Real> & initialSet, const std::vector<Interval> & domain, const unsigned int order, const Interval & cutoff_threshold)
{
	unsigned int rangeDim = Phi.rows();
	unsigned int domainDim = domain.size();

	result.clear();

	std::vector<Interval> polyRangeX0;
	std::vector<Interval> range_of_X0(rangeDim);

	initialSet.polyRangeNormal(polyRangeX0, interval_utm_setting.val_exp_table);

	for(int k=0; k<rangeDim; ++k)
	{
		range_of_X0[k] = polyRangeX0[k] + initialSet.tms[k].remainder;
	}


	TaylorModelVec<Real> tmvTmp;

	for(int i=0; i<rangeDim; ++i)
	{
		TaylorModel<Real> tmTmp1;

		for(int j=0; j<rangeDim; ++j)
		{
			TaylorModel<Real> tmTmp2(Phi[i][j], rangeDim+1, true);
			tmTmp2.mul_assign(j+1, range_of_X0[j]);

			tmTmp1 += tmTmp2;
		}

		tmvTmp.tms.push_back(tmTmp1);
	}


	tmvTmp.insert_ctrunc(result, initialSet, polyRangeX0, domain, order, cutoff_threshold);


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

void LinearFlowpipe::evaluate(TaylorModelVec<Real> & result, const std::vector<unsigned int> & outputAxes, const TaylorModelVec<Real> & tmv_of_X0, const std::vector<Interval> & polyRangeX0, const std::vector<Interval> & range_of_X0, const std::vector<Interval> & domain, const Taylor_Model_Setting & tm_setting)
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
			TaylorModel<Real> tmTmp2(Phi[outputAxes[i]][j], rangeDim+1, true);
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

	Phi = flowpipe.Phi;
	Psi = flowpipe.Psi;
	tv_remainder = flowpipe.tv_remainder;
	remainder_constraints = flowpipe.remainder_constraints;

	return *this;
}










Result_of_Reachability::Result_of_Reachability()
{
	status = -1;
	num_of_flowpipes = 0;
}

Result_of_Reachability::Result_of_Reachability(const Result_of_Reachability & result)
{
	status					= result.status;
	num_of_flowpipes		= result.num_of_flowpipes;
	fp_end_of_time			= result.fp_end_of_time;
	tmv_end_of_time			= result.tmv_end_of_time;
	linear_flowpipes		= result.linear_flowpipes;
	nonlinear_flowpipes		= result.nonlinear_flowpipes;
	orders_of_flowpipes		= result.orders_of_flowpipes;
	safety_of_flowpipes		= result.safety_of_flowpipes;
	contraction_of_flowpipes= result.contraction_of_flowpipes;
	tmv_flowpipes			= result.tmv_flowpipes;
	tmv_flowpipes_domains	= result.tmv_flowpipes_domains;
}

Result_of_Reachability::~Result_of_Reachability()
{
}

void Result_of_Reachability::clear()
{
	status = -1;
	num_of_flowpipes = 0;

	linear_flowpipes.clear();
	nonlinear_flowpipes.clear();
	tmv_flowpipes.clear();
	tmv_flowpipes_domains.clear();
	orders_of_flowpipes.clear();
	safety_of_flowpipes.clear();
	contraction_of_flowpipes.clear();
}

void Result_of_Reachability::merge(const Result_of_Reachability & result)
{
	tmv_flowpipes.insert(tmv_flowpipes.end(), result.tmv_flowpipes.begin(), result.tmv_flowpipes.end());
	tmv_flowpipes_domains.insert(tmv_flowpipes_domains.end(), result.tmv_flowpipes_domains.begin(), result.tmv_flowpipes_domains.end());
	orders_of_flowpipes.insert(orders_of_flowpipes.end(), result.orders_of_flowpipes.begin(), result.orders_of_flowpipes.end());
	safety_of_flowpipes.insert(safety_of_flowpipes.end(), result.safety_of_flowpipes.begin(), result.safety_of_flowpipes.end());
}

void Result_of_Reachability::transformToTaylorModels(const Taylor_Model_Setting & tm_setting, const bool bPrint)
{
	if(nonlinear_flowpipes.size() == 0)
		return;

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
		tmv_flowpipes_domains.push_back(fpIter->domain);

		if(bPrint)
		{
			++prog;
			printf("\b\b\b");
			printf(BOLD_FONT "%%" RESET_COLOR);
			printf(BOLD_FONT "%2d" RESET_COLOR, (int)(prog*100/total_size));
			fflush(stdout);
		}
	}

	nonlinear_flowpipes.clear();

	if(bPrint)
	{
		printf("\nDone.\n");
	}
}

void Result_of_Reachability::transformToTaylorModels(const Computational_Setting & c_setting)
{
	transformToTaylorModels(c_setting.tm_setting, c_setting.bPrint);
}

void Result_of_Reachability::computeBoxOverapproximations(std::list<std::vector<Interval> > & boxes, const Taylor_Model_Setting & tm_setting, const bool bPrint)
{
	unsigned int prog = 0, total_size = nonlinear_flowpipes.size();

	if(bPrint)
	{
		printf("Computing box overapproximations...\n");
	}

	std::list<TaylorModelVec<Real> >::const_iterator tmvIter = tmv_flowpipes.begin();
	std::list<std::vector<Interval> >::const_iterator domainIter = tmv_flowpipes_domains.begin();

	for(; tmvIter != tmv_flowpipes.end(); ++tmvIter, ++domainIter)
	{
		std::vector<Interval> box;
		tmvIter->intEval(box, *domainIter);

		boxes.push_back(box);

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

void Result_of_Reachability::computeBoxOverapproximations(std::list<std::vector<Interval> > & boxes, const Computational_Setting & c_setting)
{
	computeBoxOverapproximations(boxes, c_setting.tm_setting, c_setting.bPrint);
}

void Result_of_Reachability::computeDiscreteBoxOverapproximations(std::list<std::vector<Interval> > & boxes, const Taylor_Model_Setting & tm_setting, const bool bPrint)
{
	unsigned int prog = 0, total_size = nonlinear_flowpipes.size();

	if(bPrint)
	{
		printf("Computing discrete box overapproximations...\n");
	}

	std::list<TaylorModelVec<Real> >::const_iterator tmvIter = tmv_flowpipes.begin();
	std::list<std::vector<Interval> >::const_iterator domainIter = tmv_flowpipes_domains.begin();

	for(; tmvIter != tmv_flowpipes.end(); ++tmvIter, ++domainIter)
	{
		std::vector<Interval> box;

		std::vector<Interval> newDomain = *domainIter;
		newDomain[0] = (*domainIter)[0].sup();
		tmvIter->intEval(box, newDomain);

		boxes.push_back(box);

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

void Result_of_Reachability::computeDiscreteBoxOverapproximations(std::list<std::vector<Interval> > & boxes, const Computational_Setting & c_setting)
{
	computeDiscreteBoxOverapproximations(boxes, c_setting.tm_setting, c_setting.bPrint);
}


void Result_of_Reachability::transformToTaylorModels(const Taylor_Model_Setting & tm_setting, const bool bPrint, const Flowpipe & initialSet)
{
	if(linear_flowpipes.size() == 0)
		return;

	unsigned int prog = 0, total_size = linear_flowpipes.size();

	if(bPrint)
	{
		printf("Translating the flowpipes...\n");
	}

	std::vector<Interval> newDomain = initialSet.domain;

	if(tm_setting.step_exp_table.size() > 0)
	{
		newDomain[0] = tm_setting.step_exp_table[1];
	}

	std::vector<Interval> polyRangeX0;
	initialSet.tmvPre.polyRange(polyRangeX0, initialSet.domain);

	unsigned int rangeDim = initialSet.tmvPre.tms.size();
	std::vector<Interval> range_of_X0(rangeDim);

	for(int k=0; k<rangeDim; ++k)
	{
		range_of_X0[k] = polyRangeX0[k] + initialSet.tmvPre.tms[k].remainder;
	}

	std::list<LinearFlowpipe>::iterator iter;

	for(iter = linear_flowpipes.begin(); iter != linear_flowpipes.end(); ++iter)
	{
		TaylorModelVec<Real> tmvTmp;

		iter->evaluate(tmvTmp, initialSet.tmvPre, polyRangeX0, range_of_X0, newDomain, tm_setting);

		tmv_flowpipes.push_back(tmvTmp);
		tmv_flowpipes_domains.push_back(newDomain);

		if(bPrint)
		{
			++prog;
			printf("\b\b\b");
			printf(BOLD_FONT "%%" RESET_COLOR);
			printf(BOLD_FONT "%2d" RESET_COLOR, (int)(prog*100/total_size));
			fflush(stdout);
		}
	}

	linear_flowpipes.clear();

	if(bPrint)
	{
		printf("\nDone.\n");
	}
}

void Result_of_Reachability::transformToTaylorModels(const Computational_Setting & c_setting, const Flowpipe & initialSet)
{
	transformToTaylorModels(c_setting.tm_setting, c_setting.bPrint, initialSet);
}

Result_of_Reachability & Result_of_Reachability::operator = (const Result_of_Reachability & result)
{
	if(this == &result)
		return *this;

	status					= result.status;
	num_of_flowpipes		= result.num_of_flowpipes;
	fp_end_of_time			= result.fp_end_of_time;
	tmv_end_of_time			= result.tmv_end_of_time;
	linear_flowpipes		= result.linear_flowpipes;
	nonlinear_flowpipes		= result.nonlinear_flowpipes;
	orders_of_flowpipes		= result.orders_of_flowpipes;
	safety_of_flowpipes		= result.safety_of_flowpipes;
	contraction_of_flowpipes= result.contraction_of_flowpipes;
	tmv_flowpipes			= result.tmv_flowpipes;
	tmv_flowpipes_domains	= result.tmv_flowpipes_domains;

	return *this;
}









Linear_Time_Invariant_Dynamics::Linear_Time_Invariant_Dynamics(const Matrix<Real> & A, const Matrix<UnivariateTaylorModel<Real> > & B)
{
	rm_dyn_A = A;
	utm_dyn_B = B;

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

Linear_Time_Invariant_Dynamics::Linear_Time_Invariant_Dynamics(const Matrix<Real> & A, const Matrix<UnivariateTaylorModel<Real> > & B, const Matrix<Real> & C)
{
	rm_dyn_A = A;
	utm_dyn_B = B;
	rm_dyn_C = C;

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
	rm_dyn_C		= dynamics.rm_dyn_C;
	connectivity	= dynamics.connectivity;

	return *this;
}

void Linear_Time_Invariant_Dynamics::evaluate(Interval & result, const unsigned int varID, const TaylorModelVec<Real> & tmv_range, const std::vector<Interval> & domain)
{
	TaylorModel<Real> tmTemp;

	for(unsigned int j=0; j<rm_dyn_A.cols(); ++j)
	{
		tmTemp += tmv_range.tms[j] * rm_dyn_A[varID][j];
	}

	tmTemp.intEval(result, domain);

	if(utm_dyn_B.rows() > 0)
	{
		result += utm_dyn_B[varID][0].expansion.coefficients[0];
		result += utm_dyn_B[varID][0].remainder;
	}

	if(rm_dyn_C.rows() > 0)
	{
		double bound = rm_dyn_C[varID][0].toDouble();

		for(unsigned int j=1; j<rm_dyn_C.cols(); ++j)
		{
			bound += rm_dyn_C[varID][j].abs();
		}

		Interval I(-bound, bound);
		result += I;
	}
}

void Linear_Time_Invariant_Dynamics::evaluate(Interval & result, const std::vector<Expression<Real> > & coeff_of_Lie_deriv, const TaylorModelVec<Real> & tmv_range, const std::vector<Interval> & domain, const unsigned int order, const Computational_Setting & setting)
{
	TaylorModelVec<Real> tmDerivRng;

	for(unsigned int i=0; i<rm_dyn_A.rows(); ++i)
	{
		TaylorModel<Real> tmTemp;

		for(unsigned int j=0; j<rm_dyn_A.cols(); ++j)
		{
			tmTemp += tmv_range.tms[j] * rm_dyn_A[i][j];
		}

		if(utm_dyn_B.rows() > 0)
		{
			Polynomial<Real> constant(utm_dyn_B[i][0].expansion.coefficients[0], tmv_range.tms[0].numOfVars());
			tmTemp.expansion += constant;
			tmTemp.remainder += utm_dyn_B[i][0].remainder;
		}

		if(rm_dyn_C.rows() > 0)
		{
			double bound = rm_dyn_C[i][0].toDouble();

			for(unsigned int j=1; j<rm_dyn_C.cols(); ++j)
			{
				bound += rm_dyn_C[i][j].abs();
			}

			Interval I(-bound, bound);
			tmTemp.remainder += I;
		}

		tmDerivRng.tms.push_back(tmTemp);
	}

	result = 0;

	for(unsigned int i=0; i<coeff_of_Lie_deriv.size(); ++i)
	{
		TaylorModel<Real> tmTemp;

		coeff_of_Lie_deriv[i].evaluate(tmTemp, tmv_range.tms, order, domain, setting.tm_setting.cutoff_threshold, setting.g_setting);

		tmTemp.mul_ctrunc_assign(tmDerivRng.tms[i], domain, order, setting.tm_setting.cutoff_threshold);

		Interval I;
		tmTemp.intEval(I, domain);

		result += I;
	}
}


int Linear_Time_Invariant_Dynamics::reach(std::list<LinearFlowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const int zono_order, const Flowpipe & initialSet, const Taylor_Model_Setting & tm_setting,
		const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput)
{
	// find a proper parameter r = 2^n such that |A*delta/r| < 0.1
	Real A_max = rm_dyn_A.max_norm();
	Real threshold = 0.1;
	double step = tm_setting.step_exp_table[1].sup();
	Real rStep = step;

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

	if(utm_dyn_B.rows() > 0)
	{
		utm_Psi_0 = utm_Phi_0 * utm_dyn_B;
		utm_Psi_0.integral(intStep);
		utm_Psi_0.evaluate(utm_Psi, step_end_exp_table);
		utm_Psi_0.ctrunc(tm_setting.order);
		flowpipe.Psi = utm_Psi_0;
	}

	Zonotope global_tv_remainder(rangeDim);
	Matrix<UnivariateTaylorModel<Real> > utm_Omega;

	Interval intUnit(-1,1);
	Matrix<Interval> uncertain_range(rm_dyn_C.cols(), 1, intUnit);
	Zonotope zono_step;

	if(rm_dyn_C.rows() > 0)
	{
		utm_Omega = utm_Phi_0 * rm_dyn_C;
		utm_Omega.integral(intStep);

		Matrix<Interval> im_temp(rangeDim, utm_Omega.cols());
		utm_Omega.evaluate(im_temp, step_end_exp_table);
		im_temp *= uncertain_range;

		Zonotope zonoTmp(im_temp);
		zono_step = zonoTmp;

		flowpipe.tv_remainder = global_tv_remainder + zono_step;

		global_tv_remainder = rm_Phi * flowpipe.tv_remainder;
	}



	utm_Phi_0.ctrunc(tm_setting.order);
	flowpipe.Phi = utm_Phi_0;

	interval_utm_setting.setOrder(tm_setting.order);


	// perform the safety checking on the first flowpipe
	int checking_result = COMPLETED_SAFE;
	std::vector<Interval> polyRangeX0;
	std::vector<Interval> range_of_X0(rangeDim);

	std::vector<Interval> domain = initialSet.domain;
	domain[0] = intStep;


	if(bSafetyChecking)
	{
		initialSet.tmvPre.polyRangeNormal(polyRangeX0, interval_utm_setting.val_exp_table);

		for(int k=0; k<rangeDim; ++k)
		{
			range_of_X0[k] = polyRangeX0[k] + initialSet.tmvPre.tms[k].remainder;
		}

		if(bTMOutput || bPlot)
		{
			flowpipes.push_back(flowpipe);
		}

		int safety;

		safety = flowpipe.safetyChecking(unsafeSet, tm_setting, g_setting, initialSet.tmvPre, polyRangeX0, range_of_X0, domain);

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
			flowpipes.push_back(flowpipe);
			flowpipes_safety.push_back(SAFE);
		}
	}


	Matrix<UnivariateTaylorModel<Real> > utm_global_Psi = utm_Psi;
	Matrix<Real> rm_global_Phi = rm_Phi;


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

		newFlowpipe.Phi = rm_global_Phi * utm_Phi_0;

		if(utm_dyn_B.rows() > 0)
		{
			newFlowpipe.Psi = utm_Phi_0 * utm_global_Psi + utm_Psi_0;
			utm_global_Psi += rm_global_Phi * utm_Psi;
		}

		if(rm_dyn_C.rows() > 0)
		{
			if(zono_order >= 0 && global_tv_remainder.numOfGen()/rangeDim > zono_order)
			{
				global_tv_remainder.simplify();
			}

			newFlowpipe.tv_remainder = global_tv_remainder + zono_step;

			global_tv_remainder = rm_Phi * newFlowpipe.tv_remainder;
		}

		newFlowpipe.Phi.evaluate(rm_global_Phi, step_end_exp_table);

		++num_of_flowpipes;

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
		else
		{
			if(bTMOutput || bPlot)
			{
				flowpipes.push_back(newFlowpipe);
				flowpipes_safety.push_back(SAFE);
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

int Linear_Time_Invariant_Dynamics::reach_one_step(TaylorModelVec<Real> & result, const double stepsize, const TaylorModelVec<Real> & initialSet,
		std::vector<Interval> & domain, const unsigned int order, const Interval & cutoff_threshold)
{
	// find a proper parameter r = 2^n such that |A*delta/r| < 0.1
	Real A_max = rm_dyn_A.max_norm();
	Real threshold = 0.1;
	double step = stepsize;
	Real rStep = step;

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

	unsigned int approx_order = findProperOrder(error, A_max, A_min, tolerance, order);


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

	if(utm_dyn_B.rows() > 0)
	{
		utm_Psi_0 = utm_Phi_0 * utm_dyn_B;
		utm_Psi_0.integral(intStep);
		utm_Psi_0.evaluate(utm_Psi, step_end_exp_table);
		utm_Psi_0.ctrunc(order);
		flowpipe.Psi = utm_Psi_0;
	}

	Zonotope global_tv_remainder(rangeDim);
	Matrix<UnivariateTaylorModel<Real> > utm_Omega;

	Interval intUnit(-1,1);
	Matrix<Interval> uncertain_range(rm_dyn_C.cols(), 1, intUnit);
	Zonotope zono_step;

	if(rm_dyn_C.rows() > 0)
	{
		utm_Omega = utm_Phi_0 * rm_dyn_C;
		utm_Omega.integral(intStep);

		Matrix<Interval> im_temp(rangeDim, utm_Omega.cols());
		utm_Omega.evaluate(im_temp, step_end_exp_table);
		im_temp *= uncertain_range;

		Zonotope zonoTmp(im_temp);
		zono_step = zonoTmp;

		flowpipe.tv_remainder = global_tv_remainder + zono_step;

		global_tv_remainder = rm_Phi * flowpipe.tv_remainder;
	}



	utm_Phi_0.ctrunc(order);
	flowpipe.Phi = utm_Phi_0;

	domain[0] = intStep;
	flowpipe.evaluate(result, initialSet, domain, order, cutoff_threshold);
/*
	std::vector<Interval> range_x;
	initialSet.intEval(range_x, domain);

	for(int i=0; i<range_x.size(); ++i)
	{
		std::cout << range_x[i] << std::endl;
	}
*/

	return 0;
}

int Linear_Time_Invariant_Dynamics::reach(std::list<LinearFlowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const int zono_order, LinearFlowpipe & l_initial_set, const Flowpipe & g_initial_set,
		const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput)
{
	// find a proper parameter r = 2^n such that |A*delta/r| < 0.1
	Real A_max = rm_dyn_A.max_norm();
	Real threshold = 0.1;
	double step = tm_setting.step_exp_table[1].sup();
	Real rStep = step;

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
	utm_Phi_0.evaluate(rm_Phi, tm_setting.step_end_exp_table);


	// compute the linear mapping for the constant part
	Matrix<UnivariateTaylorModel<Real> > utm_Psi_0;
	Matrix<UnivariateTaylorModel<Real> > utm_Psi(rangeDim, 1);

	LinearFlowpipe flowpipe;

	if(utm_dyn_B.rows() > 0)
	{
		utm_Psi_0 = utm_Phi_0 * utm_dyn_B;
		utm_Psi_0.integral(intStep);
		utm_Psi_0.evaluate(utm_Psi, tm_setting.step_end_exp_table);
		utm_Psi_0.ctrunc(tm_setting.order);

		if(l_initial_set.Psi.rows() > 0)
			flowpipe.Psi = utm_Psi_0 + utm_Phi_0 * l_initial_set.Psi;
		else
			flowpipe.Psi = utm_Psi_0;
	}
	else if(l_initial_set.Psi.rows() > 0)
	{
		flowpipe.Psi = utm_Phi_0 * l_initial_set.Psi;
	}


	Zonotope global_tv_remainder = rm_Phi * l_initial_set.tv_remainder;
	Matrix<UnivariateTaylorModel<Real> > utm_Omega;

	Interval intUnit(-1,1);
	Matrix<Interval> uncertain_range(rm_dyn_C.cols(), 1, intUnit);
	Zonotope zono_step;

	if(rm_dyn_C.rows() > 0)
	{
		utm_Omega = utm_Phi_0 * rm_dyn_C;
		utm_Omega.integral(intStep);

		Matrix<Interval> im_temp(rangeDim, utm_Omega.cols());
		utm_Omega.evaluate(im_temp, tm_setting.step_end_exp_table);
		im_temp *= uncertain_range;

		Zonotope zonoTmp(im_temp);
		zono_step = zonoTmp;

		flowpipe.tv_remainder = global_tv_remainder + zono_step;

		global_tv_remainder = rm_Phi * flowpipe.tv_remainder;
	}
	else if(!global_tv_remainder.isEmpty())
	{
		flowpipe.tv_remainder = global_tv_remainder;
		global_tv_remainder = rm_Phi * flowpipe.tv_remainder;
	}



	utm_Phi_0.ctrunc(tm_setting.order);
	flowpipe.Phi = utm_Phi_0 * l_initial_set.Phi;

	interval_utm_setting.setOrder(tm_setting.order);


	// perform the safety checking on the first flowpipe
	int checking_result = COMPLETED_SAFE;
	std::vector<Interval> polyRangeX0;
	std::vector<Interval> range_of_X0(rangeDim);

	std::vector<Interval> domain = g_initial_set.domain;
	domain[0] = intStep;


	if(bSafetyChecking)
	{
		g_initial_set.tmvPre.polyRangeNormal(polyRangeX0, interval_utm_setting.val_exp_table);

		for(int k=0; k<rangeDim; ++k)
		{
			range_of_X0[k] = polyRangeX0[k] + g_initial_set.tmvPre.tms[k].remainder;
		}

		if(bTMOutput || bPlot)
		{
			flowpipes.push_back(flowpipe);
		}

		int safety;

		safety = flowpipe.safetyChecking(unsafeSet, tm_setting, g_setting, g_initial_set.tmvPre, polyRangeX0, range_of_X0, domain);

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
			flowpipes.push_back(flowpipe);
			flowpipes_safety.push_back(SAFE);
		}
	}


	Matrix<UnivariateTaylorModel<Real> > utm_global_Psi = utm_Psi;

	Matrix<Real> rm_global_Phi(rangeDim, rangeDim);
	flowpipe.Phi.evaluate(rm_global_Phi, tm_setting.step_end_exp_table);

	Matrix<Real> rm_local_Phi = rm_Phi;


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

		newFlowpipe.Phi = utm_Phi_0 * rm_global_Phi;

		if(utm_dyn_B.rows() > 0)
		{
			if(l_initial_set.Psi.rows() > 0)
				newFlowpipe.Psi = utm_Phi_0 * (utm_global_Psi + rm_local_Phi * l_initial_set.Psi) + utm_Psi_0;
			else
				newFlowpipe.Psi = utm_Phi_0 * utm_global_Psi + utm_Psi_0;

			utm_global_Psi += rm_local_Phi * utm_Psi;
			rm_local_Phi = rm_Phi * rm_local_Phi;
		}
		else if(l_initial_set.Psi.rows() > 0)
		{
			newFlowpipe.Psi = utm_Phi_0 * (rm_local_Phi * l_initial_set.Psi);
			rm_local_Phi = rm_Phi * rm_local_Phi;
		}

		if(rm_dyn_C.rows() > 0)
		{
			if(zono_order >= 0 && global_tv_remainder.numOfGen()/rangeDim > zono_order)
			{
				global_tv_remainder.simplify();
			}

			newFlowpipe.tv_remainder = global_tv_remainder + zono_step;

			global_tv_remainder = rm_Phi * newFlowpipe.tv_remainder;
		}
		else if(!global_tv_remainder.isEmpty())
		{
			newFlowpipe.tv_remainder = global_tv_remainder;
			global_tv_remainder = rm_Phi * newFlowpipe.tv_remainder;
		}


		newFlowpipe.Phi.evaluate(rm_global_Phi, tm_setting.step_end_exp_table);


		++num_of_flowpipes;

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
		else
		{
			if(bTMOutput || bPlot)
			{
				flowpipes.push_back(newFlowpipe);
				flowpipes_safety.push_back(SAFE);
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

int Linear_Time_Invariant_Dynamics::reach_inv(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<std::vector<Interval> > & domains, LinearFlowpipe & last_flowpipe, Flowpipe & contracted_initialSet, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const double time, const int zono_order, const Flowpipe & initialSet, const std::vector<Constraint> & invariant,
		const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking, const bool bPlot, const bool bTMOutput)
{
	// find a proper parameter r = 2^n such that |A*delta/r| < 0.1
	Real A_max = rm_dyn_A.max_norm();
	Real threshold = 0.1;
	double step = tm_setting.step_exp_table[1].sup();
	Real rStep = step;

	contracted_initialSet = initialSet;

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

	if(utm_dyn_B.rows() > 0)
	{
		utm_Psi_0 = utm_Phi_0 * utm_dyn_B;
		utm_Psi_0.integral(intStep);
		utm_Psi_0.evaluate(utm_Psi, step_end_exp_table);
		utm_Psi_0.ctrunc(tm_setting.order);
		flowpipe.Psi = utm_Psi_0;
	}

	Zonotope global_tv_remainder(rangeDim);
	Matrix<UnivariateTaylorModel<Real> > utm_Omega;

	Interval intUnit(-1,1);
	Matrix<Interval> uncertain_range(rm_dyn_C.cols(), 1, intUnit);
	Zonotope zono_step;

	if(rm_dyn_C.rows() > 0)
	{
		utm_Omega = utm_Phi_0 * rm_dyn_C;
		utm_Omega.integral(intStep);

		Matrix<Interval> im_temp(rangeDim, utm_Omega.cols());
		utm_Omega.evaluate(im_temp, step_end_exp_table);
		im_temp *= uncertain_range;

		Zonotope zonoTmp(im_temp);
		zono_step = zonoTmp;

		flowpipe.tv_remainder = global_tv_remainder + zono_step;

		global_tv_remainder = rm_Phi * flowpipe.tv_remainder;
	}


	utm_Phi_0.ctrunc(tm_setting.order);
	flowpipe.Phi = utm_Phi_0;

	interval_utm_setting.setOrder(tm_setting.order);


	// perform the safety checking on the first flowpipe
	int checking_result = COMPLETED_SAFE;


	// transform the first linear flowpipe to a Taylor model flowpipe
	std::vector<Interval> domain = initialSet.domain;
	domain[0] = intStep;

	std::vector<Interval> polyRangeX0;
	initialSet.tmvPre.polyRange(polyRangeX0, initialSet.domain);

	std::vector<Interval> range_of_X0(rangeDim);

	for(int k=0; k<rangeDim; ++k)
	{
		range_of_X0[k] = polyRangeX0[k] + initialSet.tmvPre.tms[k].remainder;
	}

	TaylorModelVec<Real> tmv_flowpipe;
	flowpipe.evaluate(tmv_flowpipe, initialSet.tmvPre, polyRangeX0, range_of_X0, domain, tm_setting);

	// contract the remainder firstly
	std::vector<Interval> tmv_flowpipe_polyRange;
	tmv_flowpipe.polyRange(tmv_flowpipe_polyRange, domain);

	std::vector<Interval> contracted_remainders(rangeDim);
	for(int k=0; k<rangeDim; ++k)
	{
		contracted_remainders[k] = tmv_flowpipe.tms[k].remainder;
	}

	remainder_contraction_int(tmv_flowpipe_polyRange, contracted_remainders, invariant);

	for(int k=0; k<rangeDim; ++k)
	{
		tmv_flowpipe.tms[k].remainder = contracted_remainders[k];
	}

	std::vector<Interval> contracted_domain = domain;
	int type = domain_contraction_int(tmv_flowpipe, contracted_domain, invariant, tm_setting.order, tm_setting.cutoff_threshold, g_setting);

	bool bContracted = false;
	bool previous_flowpipe_contracted = false;

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
	case 2: 	// time interval is contracted
	{
		Real zero(0);

		if(contracted_domain[0].greaterThan(zero))
		{
			return checking_result;
		}
		else
		{
			++num_of_flowpipes;

			bContracted = true;
			contraction_of_flowpipes.push_back(true);

			domain = contracted_domain;

			for(int k=1; k<domain.size(); ++k)
			{
				contracted_initialSet.domain[k] = contracted_domain[k];
			}


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

			double t = contracted_domain[0].sup();

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", t);
				printf("order = %d\n", tm_setting.order);
			}

			last_flowpipe = flowpipe;

			return checking_result;
		}
	}
	}


	last_flowpipe = flowpipe;


	Matrix<UnivariateTaylorModel<Real> > utm_global_Psi = utm_Psi;
	Matrix<Real> rm_global_Phi = rm_Phi;


	Matrix<Interval> im_step_rem(rangeDim, 1);
	for(int i=0; i<rangeDim; ++i)
	{
		im_step_rem[i][0] = utm_Psi[i][0].remainder;
	}

	Matrix<Interval> previous_remainder(rangeDim, 1);
	for(int k=0; k<rangeDim; ++k)
	{
		previous_remainder[k][0] = tmv_flowpipe.tms[k].remainder;
	}


	int N = (int)ceil(time/step);

	if(bPrint)
	{
		printf("time = %f,\t", step);
		printf("step = %f,\t", step);
		printf("order = %d\n", tm_setting.order);
	}


	for(int i=1; i<N; ++i)
	{
		LinearFlowpipe new_flowpipe;

		new_flowpipe.Phi = rm_global_Phi * utm_Phi_0;

		if(utm_dyn_B.rows() > 0)
		{
			new_flowpipe.Psi = utm_Phi_0 * utm_global_Psi + utm_Psi_0;
			utm_global_Psi += rm_global_Phi * utm_Psi;
		}

		if(rm_dyn_C.rows() > 0)
		{
			if(zono_order >= 0 && global_tv_remainder.numOfGen()/rangeDim > zono_order)
			{
				global_tv_remainder.simplify();
			}

			new_flowpipe.tv_remainder = global_tv_remainder + zono_step;

			global_tv_remainder = rm_Phi * new_flowpipe.tv_remainder;
		}

		new_flowpipe.Phi.evaluate(rm_global_Phi, step_end_exp_table);



		// intersect the flowpipe with the invariant and check the safety

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


		if(bContracted)
		{
			// refine the remainder
			Matrix<Interval> refinement = rm_Phi * previous_remainder + im_step_rem;

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
		tmv_flowpipe_polyRange.clear();
		tmv_flowpipe.polyRange(tmv_flowpipe_polyRange, domain);

		for(int k=0; k<rangeDim; ++k)
		{
			contracted_remainders[k] = tmv_flowpipe.tms[k].remainder;
		}

		remainder_contraction_int(tmv_flowpipe_polyRange, contracted_remainders, invariant);

		for(int k=0; k<rangeDim; ++k)
		{
			tmv_flowpipe.tms[k].remainder = contracted_remainders[k];
		}




		contracted_domain = domain;
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
		case 2: 	// time interval is contracted
		{
			Real zero(0);

			if(contracted_domain[0].greaterThan(zero))
			{
				return checking_result;
			}
			else
			{
				++num_of_flowpipes;

				bContracted = true;
				contraction_of_flowpipes.push_back(true);

				domain = contracted_domain;

				for(int k=1; k<domain.size(); ++k)
				{
					contracted_initialSet.domain[k] = contracted_domain[k];
				}

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

				double t = contracted_domain[0].sup();

				if(bPrint)
				{
					printf("time = %f,\t", i*step + t);
					printf("step = %f,\t", t);
					printf("order = %d\n", tm_setting.order);
				}

				last_flowpipe = new_flowpipe;

				return checking_result;
			}
		}
		}

		if(bContracted)
		{
			for(int k=0; k<rangeDim; ++k)
			{
				previous_remainder[k][0] = tmv_flowpipe.tms[k].remainder;
			}
		}

		if(bPrint)
		{
			printf("time = %f,\t", (i+1)*step);
			printf("step = %f,\t", step);
			printf("order = %d\n", tm_setting.order);
		}

		last_flowpipe = new_flowpipe;
	}

	return checking_result;
}

int Linear_Time_Invariant_Dynamics::reach_inv(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<std::vector<Interval> > & domains, LinearFlowpipe & last_flowpipe, Flowpipe & contracted_initialSet, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const double time, const int zono_order, LinearFlowpipe & l_initial_set, const Flowpipe & g_initial_set,
		const std::vector<Constraint> & invariant, const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput)
{
	// find a proper parameter r = 2^n such that |A*delta/r| < 0.1
	Real A_max = rm_dyn_A.max_norm();
	Real threshold = 0.1;
	double step = tm_setting.step_exp_table[1].sup();
	Real rStep = step;

	contracted_initialSet = g_initial_set;

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

	if(utm_dyn_B.rows() > 0)
	{
		utm_Psi_0 = utm_Phi_0 * utm_dyn_B;
		utm_Psi_0.integral(intStep);
		utm_Psi_0.evaluate(utm_Psi, step_end_exp_table);
		utm_Psi_0.ctrunc(tm_setting.order);

		if(l_initial_set.Psi.rows() > 0)
			flowpipe.Psi = utm_Psi_0 + utm_Phi_0 * l_initial_set.Psi;
		else
			flowpipe.Psi = utm_Psi_0;
	}
	else if(l_initial_set.Psi.rows() > 0)
	{
		flowpipe.Psi = utm_Phi_0 * l_initial_set.Psi;
	}


	Zonotope global_tv_remainder = rm_Phi * l_initial_set.tv_remainder;
	Matrix<UnivariateTaylorModel<Real> > utm_Omega;

	Interval intUnit(-1,1);
	Matrix<Interval> uncertain_range(rm_dyn_C.cols(), 1, intUnit);
	Zonotope zono_step;

	if(rm_dyn_C.rows() > 0)
	{
		utm_Omega = utm_Phi_0 * rm_dyn_C;
		utm_Omega.integral(intStep);

		Matrix<Interval> im_temp(rangeDim, utm_Omega.cols());
		utm_Omega.evaluate(im_temp, step_end_exp_table);
		im_temp *= uncertain_range;

		Zonotope zonoTmp(im_temp);
		zono_step = zonoTmp;

		flowpipe.tv_remainder = global_tv_remainder + zono_step;

		global_tv_remainder = rm_Phi * flowpipe.tv_remainder;
	}
	else if(!global_tv_remainder.isEmpty())
	{
		flowpipe.tv_remainder = global_tv_remainder;
		global_tv_remainder = rm_Phi * flowpipe.tv_remainder;
	}


	utm_Phi_0.ctrunc(tm_setting.order);
	flowpipe.Phi = utm_Phi_0 * l_initial_set.Phi;

	interval_utm_setting.setOrder(tm_setting.order);


	// perform the safety checking on the first flowpipe
	int checking_result = COMPLETED_SAFE;


	// transform the first linear flowpipe to a Taylor model flowpipe
	std::vector<Interval> domain = g_initial_set.domain;
	domain[0] = intStep;

	std::vector<Interval> polyRangeX0;
	g_initial_set.tmvPre.polyRange(polyRangeX0, g_initial_set.domain);

	std::vector<Interval> range_of_X0(rangeDim);

	for(int k=0; k<rangeDim; ++k)
	{
		range_of_X0[k] = polyRangeX0[k] + g_initial_set.tmvPre.tms[k].remainder;
	}

	TaylorModelVec<Real> tmv_flowpipe;
	flowpipe.evaluate(tmv_flowpipe, g_initial_set.tmvPre, polyRangeX0, range_of_X0, domain, tm_setting);

	// contract the remainder firstly
	std::vector<Interval> tmv_flowpipe_polyRange;
	tmv_flowpipe.polyRange(tmv_flowpipe_polyRange, domain);

	std::vector<Interval> contracted_remainders(rangeDim);
	for(int k=0; k<rangeDim; ++k)
	{
		contracted_remainders[k] = tmv_flowpipe.tms[k].remainder;
	}

	remainder_contraction_int(tmv_flowpipe_polyRange, contracted_remainders, invariant);

	for(int k=0; k<rangeDim; ++k)
	{
		tmv_flowpipe.tms[k].remainder = contracted_remainders[k];
	}

	std::vector<Interval> contracted_domain = domain;
	int type = domain_contraction_int(tmv_flowpipe, contracted_domain, invariant, tm_setting.order, tm_setting.cutoff_threshold, g_setting);

	bool bContracted = false;
	bool previous_flowpipe_contracted = false;


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
	case 2: 	// time interval is contracted
	{
		Real zero(0);

		if(contracted_domain[0].greaterThan(zero))
		{
			return checking_result;
		}
		else
		{
			++num_of_flowpipes;

			bContracted = true;
			contraction_of_flowpipes.push_back(true);

			domain = contracted_domain;

			for(int k=1; k<domain.size(); ++k)
			{
				contracted_initialSet.domain[k] = contracted_domain[k];
			}


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

			double t = contracted_domain[0].sup();

			if(bPrint)
			{
				printf("time = %f,\t", t);
				printf("step = %f,\t", t);
				printf("order = %d\n", tm_setting.order);
			}

			last_flowpipe = flowpipe;

			return checking_result;
		}
	}
	}


	last_flowpipe = flowpipe;


	Matrix<UnivariateTaylorModel<Real> > utm_global_Psi = utm_Psi;

	Matrix<Real> rm_global_Phi(rangeDim, rangeDim);
	flowpipe.Phi.evaluate(rm_global_Phi, tm_setting.step_end_exp_table);

	Matrix<Real> rm_local_Phi = rm_Phi;


	Matrix<Interval> im_step_rem(rangeDim, 1);
	for(int i=0; i<rangeDim; ++i)
	{
		im_step_rem[i][0] = utm_Psi[i][0].remainder;
	}

	Matrix<Interval> previous_remainder(rangeDim, 1);
	for(int k=0; k<rangeDim; ++k)
	{
		previous_remainder[k][0] = tmv_flowpipe.tms[k].remainder;
	}


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

		newFlowpipe.Phi = utm_Phi_0 * rm_global_Phi;

		if(utm_dyn_B.rows() > 0)
		{
			if(l_initial_set.Psi.rows() > 0)
				newFlowpipe.Psi = utm_Phi_0 * (utm_global_Psi + rm_local_Phi * l_initial_set.Psi) + utm_Psi_0;
			else
				newFlowpipe.Psi = utm_Phi_0 * utm_global_Psi + utm_Psi_0;

			utm_global_Psi += rm_local_Phi * utm_Psi;
			rm_local_Phi = rm_Phi * rm_local_Phi;
		}
		else if(l_initial_set.Psi.rows() > 0)
		{
			newFlowpipe.Psi = utm_Phi_0 * (rm_local_Phi * l_initial_set.Psi);
			rm_local_Phi = rm_Phi * rm_local_Phi;
		}

		if(rm_dyn_C.rows() > 0)
		{
			if(zono_order >= 0 && global_tv_remainder.numOfGen()/rangeDim > zono_order)
			{
				global_tv_remainder.simplify();
			}

			newFlowpipe.tv_remainder = global_tv_remainder + zono_step;

			global_tv_remainder = rm_Phi * newFlowpipe.tv_remainder;
		}
		else if(!global_tv_remainder.isEmpty())
		{
			newFlowpipe.tv_remainder = global_tv_remainder;
			global_tv_remainder = rm_Phi * newFlowpipe.tv_remainder;
		}


		// intersect the flowpipe with the invariant and check the safety

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
		newFlowpipe.evaluate(tmv_flowpipe, g_initial_set.tmvPre, polyRangeX0, range_of_X0, domain, tm_setting);


		if(bContracted)
		{
			// refine the remainder
			Matrix<Interval> refinement = rm_Phi * previous_remainder + im_step_rem;

			if(newFlowpipe.remainder_constraints.rows() == 0)
			{
				Matrix<Interval> im_temp(rangeDim, 1);
				newFlowpipe.remainder_constraints = im_temp;
			}

			// contracting the remainder by updating the remainder constraints
			for(int j=0; j<rangeDim; ++j)
			{
				tmv_flowpipe.tms[j].remainder.intersect_assign(refinement[j][0]);
				newFlowpipe.remainder_constraints[j][0] = tmv_flowpipe.tms[j].remainder;
			}
		}


		// contract the remainder firstly
		tmv_flowpipe_polyRange.clear();
		tmv_flowpipe.polyRange(tmv_flowpipe_polyRange, domain);

		for(int k=0; k<rangeDim; ++k)
		{
			contracted_remainders[k] = tmv_flowpipe.tms[k].remainder;
		}

		remainder_contraction_int(tmv_flowpipe_polyRange, contracted_remainders, invariant);

		for(int k=0; k<rangeDim; ++k)
		{
			tmv_flowpipe.tms[k].remainder = contracted_remainders[k];
		}




		contracted_domain = domain;
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
		case 2: 	// time interval is contracted
		{
			Real zero(0);

			if(contracted_domain[0].greaterThan(zero))
			{
				return checking_result;
			}
			else
			{
				++num_of_flowpipes;

				bContracted = true;
				contraction_of_flowpipes.push_back(true);

				domain = contracted_domain;

				for(int k=1; k<domain.size(); ++k)
				{
					contracted_initialSet.domain[k] = contracted_domain[k];
				}

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

				double t = contracted_domain[0].sup();

				if(bPrint)
				{
					printf("time = %f,\t", i*step + t);
					printf("step = %f,\t", t);
					printf("order = %d\n", tm_setting.order);
				}

				last_flowpipe = newFlowpipe;

				return checking_result;
			}
		}
		}

		if(bContracted)
		{
			for(int k=0; k<rangeDim; ++k)
			{
				previous_remainder[k][0] = tmv_flowpipe.tms[k].remainder;
			}
		}

		if(bPrint)
		{
			printf("time = %f,\t", (i+1)*step);
			printf("step = %f,\t", step);
			printf("order = %d\n", tm_setting.order);
		}

		last_flowpipe = newFlowpipe;
	}

	return checking_result;
}

void Linear_Time_Invariant_Dynamics::reach(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const std::vector<Constraint> & unsafeSet, const int zono_order)
{
	bool bSafetyChecking = false;

	if(unsafeSet.size() > 0)
	{
		bSafetyChecking = true;
	}

	result.status = reach(result.linear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes,
			result.num_of_flowpipes, setting.time, zono_order, initialSet, setting.tm_setting,
			setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);


	// evaluate the Taylor model overapproximation at the end of the time
	if(result.linear_flowpipes.size() > 0)
	{
		std::vector<Interval> domain = initialSet.domain;
		domain[0] = setting.tm_setting.step_exp_table[1];

		std::vector<Interval> polyRangeX0;
		initialSet.tmvPre.polyRange(polyRangeX0, initialSet.domain);

		unsigned int rangeDim = initialSet.tmvPre.tms.size();
		std::vector<Interval> range_of_X0(rangeDim);

		for(int k=0; k<rangeDim; ++k)
		{
			range_of_X0[k] = polyRangeX0[k] + initialSet.tmvPre.tms[k].remainder;
		}

		TaylorModelVec<Real> tmvTmp;

		result.linear_flowpipes.back().evaluate(tmvTmp, initialSet.tmvPre, polyRangeX0, range_of_X0, domain, setting.tm_setting);

		tmvTmp.evaluate_time(result.tmv_end_of_time, setting.tm_setting.step_end_exp_table);


		// Flowpipe format for the flowpipe at the end of the time
		domain[0] = 0;
		Flowpipe fp_end_of_time(result.tmv_end_of_time, domain, setting.tm_setting.cutoff_threshold);
		result.fp_end_of_time = fp_end_of_time;




		result.lfp_end_of_time = result.linear_flowpipes.back();
		for(int i=0; i<result.lfp_end_of_time.Phi.rows(); ++i)
		{
			for(int j=0; j<result.lfp_end_of_time.Phi.cols(); ++j)
			{
				Real c;
				result.lfp_end_of_time.Phi[i][j].expansion.evaluate(c, setting.tm_setting.step_end_exp_table);
				result.lfp_end_of_time.Phi[i][j].expansion = c;
			}
		}

		for(int i=0; i<result.lfp_end_of_time.Psi.rows(); ++i)
		{
			for(int j=0; j<result.lfp_end_of_time.Psi.cols(); ++j)
			{
				Real c;
				result.lfp_end_of_time.Psi[i][j].expansion.evaluate(c, setting.tm_setting.step_end_exp_table);
				result.lfp_end_of_time.Psi[i][j].expansion = c;
			}
		}
	}
}

void Linear_Time_Invariant_Dynamics::reach(Result_of_Reachability & result, Computational_Setting & setting, LinearFlowpipe & l_initial_set, const Flowpipe & g_initial_set, const std::vector<Constraint> & unsafeSet, const int zono_order)
{
	bool bSafetyChecking = false;

	if(unsafeSet.size() > 0)
	{
		bSafetyChecking = true;
	}

	result.status = reach(result.linear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes,
			result.num_of_flowpipes, setting.time, zono_order, l_initial_set, g_initial_set, setting.tm_setting,
			setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);


	// evaluate the Taylor model overapproximation at the end of the time
	if(result.linear_flowpipes.size() > 0)
	{
		std::vector<Interval> domain = g_initial_set.domain;
		domain[0] = setting.tm_setting.step_exp_table[1];

		std::vector<Interval> polyRangeX0;
		g_initial_set.tmvPre.polyRange(polyRangeX0, g_initial_set.domain);

		unsigned int rangeDim = g_initial_set.tmvPre.tms.size();
		std::vector<Interval> range_of_X0(rangeDim);

		for(int k=0; k<rangeDim; ++k)
		{
			range_of_X0[k] = polyRangeX0[k] + g_initial_set.tmvPre.tms[k].remainder;
		}

		TaylorModelVec<Real> tmvTmp;

		result.linear_flowpipes.back().evaluate(tmvTmp, g_initial_set.tmvPre, polyRangeX0, range_of_X0, domain, setting.tm_setting);
		tmvTmp.evaluate_time(result.tmv_end_of_time, setting.tm_setting.step_end_exp_table);


		// Flowpipe format for the flowpipe at the end of the time
		domain[0] = 0;
		Flowpipe fp_end_of_time(result.tmv_end_of_time, domain, setting.tm_setting.cutoff_threshold);
		result.fp_end_of_time = fp_end_of_time;



		// test
		result.lfp_end_of_time = result.linear_flowpipes.back();
		for(int i=0; i<result.lfp_end_of_time.Phi.rows(); ++i)
		{
			for(int j=0; j<result.lfp_end_of_time.Phi.cols(); ++j)
			{
				Real c;
				result.lfp_end_of_time.Phi[i][j].expansion.evaluate(c, setting.tm_setting.step_end_exp_table);
				result.lfp_end_of_time.Phi[i][j].expansion = c;
			}
		}

		for(int i=0; i<result.lfp_end_of_time.Psi.rows(); ++i)
		{
			for(int j=0; j<result.lfp_end_of_time.Psi.cols(); ++j)
			{
				Real c;
				result.lfp_end_of_time.Psi[i][j].expansion.evaluate(c, setting.tm_setting.step_end_exp_table);
				result.lfp_end_of_time.Psi[i][j].expansion = c;
			}
		}
	}
}

void Linear_Time_Invariant_Dynamics::reach_inv(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, Flowpipe & contracted_initialSet, const std::vector<Constraint> & invariant, const std::vector<Constraint> & unsafeSet, const int zono_order)
{
	bool bSafetyChecking = false;

	if(unsafeSet.size() > 0)
	{
		bSafetyChecking = true;
	}

	result.status = reach_inv(result.tmv_flowpipes, result.tmv_flowpipes_domains, result.lfp_end_of_time, contracted_initialSet, result.orders_of_flowpipes, result.safety_of_flowpipes,
			result.contraction_of_flowpipes, result.num_of_flowpipes, setting.time, zono_order, initialSet, invariant, setting.tm_setting, setting.g_setting,
			setting.bPrint, unsafeSet, bSafetyChecking, true, true);


	// evaluate the Taylor model overapproximation at the end of the time
	if(result.tmv_flowpipes.size() > 0)
	{
		TaylorModelVec<Real> last_flowpipe = result.tmv_flowpipes.back();

		std::vector<Interval> domain = result.tmv_flowpipes_domains.back();

		Real t;
		domain[0].sup(t);

		std::vector<Real> realVec;
		realVec.push_back(1);
		realVec.push_back(t);

		Real tmp = t;

		for(unsigned int i=2; i<=setting.tm_setting.step_end_exp_table.size(); ++i)
		{
			tmp *= t;
			realVec.push_back(tmp);
		}

		last_flowpipe.evaluate_time(result.tmv_end_of_time, realVec);



		// Flowpipe format for the flowpipe at the end of the time
		domain[0] = 0;
		Flowpipe fp_end_of_time(result.tmv_end_of_time, domain, setting.tm_setting.cutoff_threshold);
		result.fp_end_of_time = fp_end_of_time;




		// test
		for(int i=0; i<result.lfp_end_of_time.Phi.rows(); ++i)
		{
			for(int j=0; j<result.lfp_end_of_time.Phi.cols(); ++j)
			{
				Real c;
				result.lfp_end_of_time.Phi[i][j].expansion.evaluate(c, realVec);
				result.lfp_end_of_time.Phi[i][j].expansion = c;
			}
		}

		for(int i=0; i<result.lfp_end_of_time.Psi.rows(); ++i)
		{
			for(int j=0; j<result.lfp_end_of_time.Psi.cols(); ++j)
			{
				Real c;
				result.lfp_end_of_time.Psi[i][j].expansion.evaluate(c, realVec);
				result.lfp_end_of_time.Psi[i][j].expansion = c;
			}
		}
	}
}

void Linear_Time_Invariant_Dynamics::reach_inv(Result_of_Reachability & result, Computational_Setting & setting, LinearFlowpipe & l_initial_set, const Flowpipe & g_initial_set, Flowpipe & contracted_initialSet, const std::vector<Constraint> & invariant, const std::vector<Constraint> & unsafeSet, const int zono_order)
{
	bool bSafetyChecking = false;

	if(unsafeSet.size() > 0)
	{
		bSafetyChecking = true;
	}

	result.status = reach_inv(result.tmv_flowpipes, result.tmv_flowpipes_domains, result.lfp_end_of_time, contracted_initialSet, result.orders_of_flowpipes, result.safety_of_flowpipes,
			result.contraction_of_flowpipes, result.num_of_flowpipes, setting.time, zono_order, l_initial_set, g_initial_set, invariant, setting.tm_setting, setting.g_setting,
			setting.bPrint, unsafeSet, bSafetyChecking, true, true);


	// evaluate the Taylor model overapproximation at the end of the time
	if(result.tmv_flowpipes.size() > 0)
	{
		TaylorModelVec<Real> last_flowpipe = result.tmv_flowpipes.back();

		std::vector<Interval> domain = result.tmv_flowpipes_domains.back();

		Real t;
		domain[0].sup(t);

		std::vector<Real> realVec;
		realVec.push_back(1);
		realVec.push_back(t);

		Real tmp = t;

		for(unsigned int i=2; i<=setting.tm_setting.step_end_exp_table.size(); ++i)
		{
			tmp *= t;
			realVec.push_back(tmp);
		}

		last_flowpipe.evaluate_time(result.tmv_end_of_time, realVec);


		// Flowpipe format for the flowpipe at the end of the time
		domain[0] = 0;
		Flowpipe fp_end_of_time(result.tmv_end_of_time, domain, setting.tm_setting.cutoff_threshold);
		result.fp_end_of_time = fp_end_of_time;




		// test
		for(int i=0; i<result.lfp_end_of_time.Phi.rows(); ++i)
		{
			for(int j=0; j<result.lfp_end_of_time.Phi.cols(); ++j)
			{
				Real c;
				result.lfp_end_of_time.Phi[i][j].expansion.evaluate(c, realVec);
				result.lfp_end_of_time.Phi[i][j].expansion = c;
			}
		}

		for(int i=0; i<result.lfp_end_of_time.Psi.rows(); ++i)
		{
			for(int j=0; j<result.lfp_end_of_time.Psi.cols(); ++j)
			{
				Real c;
				result.lfp_end_of_time.Psi[i][j].expansion.evaluate(c, realVec);
				result.lfp_end_of_time.Psi[i][j].expansion = c;
			}
		}
	}
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
	Matrix<Interval> tmp(C.cols(), 1, intUnit);
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

void Linear_Time_Varying_Dynamics::evaluate(Interval & result, const unsigned int varID, const TaylorModelVec<Real> & tmv_range, const std::vector<Interval> & domain, const Real & t_lb, const unsigned int order, const Interval & cutoff_threshold)
{
	TaylorModel<Real> tmTemp;
	unsigned int numVars = tmv_range.tms[0].numOfVars();

	for(unsigned int j=0; j<upm_dyn_A.cols(); ++j)
	{
		Polynomial<Real> tmp;
		upm_dyn_A[varID][j].toPolynomial(tmp, numVars, t_lb);

		TaylorModel<Real> tmp2(tmp);
		tmp2.mul_ctrunc_assign(tmv_range.tms[j], domain, order, cutoff_threshold);

		tmTemp += tmp2;
	}

	tmTemp.intEval(result, domain);

	if(upm_dyn_B.rows() > 0)
	{
		Interval t_range(t_lb, t_lb + domain[0].inf());
		Interval I;

		upm_dyn_B[varID][0].evaluate(I, t_range);
		result += I;
	}

	if(upm_dyn_tv.rows() > 0)
	{
		Interval t_range(t_lb, t_lb + domain[0].inf());
		Interval I;

		double bound = 0;

		for(unsigned int j=0; j<upm_dyn_tv.cols(); ++j)
		{
			upm_dyn_tv[varID][j].evaluate(I, t_range);
			bound += I.mag();
		}

		Interval J(-bound, bound);
		result += J;
	}
}
/*
void Linear_Time_Varying_Dynamics::evaluate(Interval & result, const std::vector<Expression<Real> > & coeff_of_Lie_deriv, const TaylorModelVec<Real> & tmv_range, const std::vector<Interval> & domain, const Interval & t, const unsigned int order, const Computational_Setting & setting)
{

}
*/

int Linear_Time_Varying_Dynamics::reach(std::list<LinearFlowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, const int zono_order, const Taylor_Model_Setting & tm_setting,
		const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput)
{
	double step = tm_setting.step_exp_table[1].sup();
	Real rStep = step;
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
	std::vector<Interval> polyRangeX0;
	std::vector<Interval> range_of_X0;


	if(bSafetyChecking)
	{
		initialSet.tmvPre.polyRangeNormal(polyRangeX0, interval_utm_setting.val_exp_table);

		std::vector<Interval> range(rangeDim);
		for(int i1=0; i1<rangeDim; ++i1)
		{
			range_of_X0[i1] = polyRangeX0[i1] + initialSet.tmvPre.tms[i1].remainder;
		}
	}

	std::vector<Interval> domain = initialSet.domain;
	domain[0] = intStep;

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
			if(zono_order >= 0 && global_tv_remainder.numOfGen()/rangeDim > zono_order)
			{
				global_tv_remainder.simplify();
			}

			Matrix<Interval> im_Phi(rangeDim, rangeDim);
			flowpipe.Phi.evaluate(im_Phi, tm_setting.step_exp_table);

			Matrix<Interval> im_temp = im_Phi * tv_part * uncertain_range;
			im_temp *= intStep;
			Zonotope zonoTmp(im_temp);
			flowpipe.tv_remainder = global_tv_remainder + zonoTmp;
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

			int safety;

			std::vector<Matrix<Real> > constraints;
			std::vector<Matrix<Interval> > precond_Phi, precond_Psi;

			safety = flowpipe.safetyChecking(unsafeSet, tm_setting, g_setting, initialSet.tmvPre, polyRangeX0, range_of_X0, domain);

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

int Linear_Time_Varying_Dynamics::reach_inv(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, const int zono_order, const std::vector<Constraint> & invariant,
		const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking, const bool bPlot, const bool bTMOutput)
{
	double step = tm_setting.step_exp_table[1].sup();
	Real rStep = step;
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

	num_of_flowpipes = 0;

	Matrix<Interval> tv_part(rangeDim, numTVPar);

	Matrix<Interval> previous_remainder(rangeDim, 1);

	std::vector<Interval> polyRangeX0;
	initialSet.tmvPre.polyRange(polyRangeX0, initialSet.domain);

	std::vector<Interval> range_of_X0(rangeDim);

	for(int k=0; k<rangeDim; ++k)
	{
		range_of_X0[k] = polyRangeX0[k] + initialSet.tmvPre.tms[k].remainder;
	}

	bool bContracted = false;
	bool previous_flowpipe_contracted = false;
	Flowpipe fp_dummy;
	fp_dummy.domain = initialSet.domain;
	fp_dummy.domain[0] = intStep;


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


		Matrix<Interval> im_temp;

		if(numTVPar > 0)
		{
			if(zono_order >= 0 && global_tv_remainder.numOfGen()/rangeDim > zono_order)
			{
				global_tv_remainder.simplify();
			}

			Matrix<Interval> im_Phi(rangeDim, rangeDim);
			flowpipe.Phi.evaluate(im_Phi, tm_setting.step_exp_table);

			im_temp = im_Phi * tv_part * uncertain_range;
			im_temp *= intStep;
			Zonotope zonoTmp(im_temp);
			flowpipe.tv_remainder = global_tv_remainder + zonoTmp;
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



		// perform the contraction on the remainder and the domain
		if(previous_flowpipe_contracted)
		{
			polyRangeX0.clear();
			initialSet.tmvPre.polyRange(polyRangeX0, fp_dummy.domain);

			for(int k=0; k<rangeDim; ++k)
			{
				range_of_X0[k] = polyRangeX0[k] + initialSet.tmvPre.tms[k].remainder;
			}

			previous_flowpipe_contracted = false;
		}



		TaylorModelVec<Real> tmv_flowpipe;
		flowpipe.evaluate(tmv_flowpipe, initialSet.tmvPre, polyRangeX0, range_of_X0, fp_dummy.domain, tm_setting);


		if(bContracted && num_of_flowpipes > 0)
		{
			// refine the remainder
			Matrix<Interval> refinement = Phi_step_end * previous_remainder;

			if(numTVPar > 0)
			{
				refinement += im_temp;
			}

			for(int j=0; j<rangeDim; ++j)
			{
				tmv_flowpipe.tms[j].remainder.intersect_assign(refinement[j][0]);
			}
		}

		// contract the remainder firstly
		std::vector<Interval> tmv_flowpipe_polyRange;
		tmv_flowpipe.polyRange(tmv_flowpipe_polyRange, fp_dummy.domain);

		std::vector<Interval> contracted_remainders(rangeDim);

		for(int k=0; k<rangeDim; ++k)
		{
			contracted_remainders[k] = tmv_flowpipe.tms[k].remainder;
		}

		remainder_contraction_int(tmv_flowpipe_polyRange, contracted_remainders, invariant);

		for(int k=0; k<rangeDim; ++k)
		{
			tmv_flowpipe.tms[k].remainder = contracted_remainders[k];
		}



		std::vector<Interval> contracted_domain = fp_dummy.domain;
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
			fp_dummy.domain = contracted_domain;

			if(bSafetyChecking)
			{
				int safety = safetyChecking(tmv_flowpipe, contracted_domain, unsafeSet, tm_setting, g_setting);

				if(bTMOutput || bPlot)
				{
					tmv_flowpipes.push_back(tmv_flowpipe);
					flowpipes.push_back(fp_dummy);
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
					flowpipes.push_back(fp_dummy);
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

			fp_dummy.domain = contracted_domain;

			if(bSafetyChecking)
			{
				int safety = safetyChecking(tmv_flowpipe, contracted_domain, unsafeSet, tm_setting, g_setting);

				if(bTMOutput || bPlot)
				{
					tmv_flowpipes.push_back(tmv_flowpipe);
					flowpipes.push_back(fp_dummy);
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
					flowpipes.push_back(fp_dummy);
					flowpipes_safety.push_back(SAFE);
				}
			}

			break;
		}
		case 2: 	// time interval is contracted
		{
			Real zero(0);

			if(contracted_domain[0].greaterThan(zero))
			{
				return checking_result;
			}
			else
			{
				++num_of_flowpipes;

				bContracted = true;
				contraction_of_flowpipes.push_back(true);

				fp_dummy.domain = contracted_domain;

				if(bSafetyChecking)
				{
					int safety = safetyChecking(tmv_flowpipe, contracted_domain, unsafeSet, tm_setting, g_setting);

					if(bTMOutput || bPlot)
					{
						tmv_flowpipes.push_back(tmv_flowpipe);
						flowpipes.push_back(fp_dummy);
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
						flowpipes.push_back(fp_dummy);
						flowpipes_safety.push_back(SAFE);
					}
				}

				double t = contracted_domain[0].sup();

				if(bPrint)
				{
					printf("time = %f,\t", t0 + t);
					printf("step = %f,\t", t);
					printf("order = %d\n", tm_setting.order);
				}

				return checking_result;
			}
		}
		}

		if(bContracted)
		{
			for(int k=0; k<rangeDim; ++k)
			{
				previous_remainder[k][0] = tmv_flowpipe.tms[k].remainder;
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

void Linear_Time_Varying_Dynamics::reach(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const std::vector<Constraint> & unsafeSet, const int zono_order)
{
	bool bSafetyChecking = false;

	if(unsafeSet.size() > 0)
	{
		bSafetyChecking = true;
	}

	result.status = reach(result.linear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes,
			result.num_of_flowpipes, setting.time, initialSet, zono_order, setting.tm_setting,
			setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);

	// evaluate the Taylor model overapproximation at the end of the time
	if(result.linear_flowpipes.size() > 0)
	{
		Flowpipe fp_dummy;
		fp_dummy.domain = initialSet.domain;
		fp_dummy.domain[0] = setting.tm_setting.step_exp_table[1];

		std::vector<Interval> polyRangeX0;
		initialSet.tmvPre.polyRange(polyRangeX0, initialSet.domain);

		unsigned int rangeDim = initialSet.tmvPre.tms.size();
		std::vector<Interval> range_of_X0(rangeDim);

		for(int k=0; k<rangeDim; ++k)
		{
			range_of_X0[k] = polyRangeX0[k] + initialSet.tmvPre.tms[k].remainder;
		}

		TaylorModelVec<Real> tmvTmp;

		result.linear_flowpipes.back().evaluate(tmvTmp, initialSet.tmvPre, polyRangeX0, range_of_X0, fp_dummy.domain, setting.tm_setting);
		tmvTmp.evaluate_time(result.tmv_end_of_time, setting.tm_setting.step_end_exp_table);

		fp_dummy.domain[0] = setting.tm_setting.step_end_exp_table[1];
		result.fp_end_of_time = fp_dummy;
	}
}

void Linear_Time_Varying_Dynamics::reach_inv(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const std::vector<Constraint> & invariant, const std::vector<Constraint> & unsafeSet, const int zono_order)
{
	bool bSafetyChecking = false;

	if(unsafeSet.size() > 0)
	{
		bSafetyChecking = true;
	}

	result.status = reach_inv(result.tmv_flowpipes, result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes,
			result.contraction_of_flowpipes, result.num_of_flowpipes, setting.time, initialSet, zono_order, invariant, setting.tm_setting, setting.g_setting,
			setting.bPrint, unsafeSet, bSafetyChecking, true, true);


	// evaluate the Taylor model overapproximation at the end of the time
	if(result.tmv_flowpipes.size() > 0)
	{
		TaylorModelVec<Real> last_flowpipe = result.tmv_flowpipes.back();

		Flowpipe fp_dummy;
		fp_dummy.domain = result.nonlinear_flowpipes.back().domain;

		Real t;
		fp_dummy.domain[0].sup(t);

		fp_dummy.domain[0] = t;
		result.fp_end_of_time = fp_dummy;

		std::vector<Real> realVec;
		realVec.push_back(1);
		realVec.push_back(t);

		Real tmp = t;

		for(unsigned int i=2; i<=setting.tm_setting.step_end_exp_table.size(); ++i)
		{
			tmp *= t;
			realVec.push_back(tmp);
		}

		last_flowpipe.evaluate_time(result.tmv_end_of_time, realVec);
	}
}





Deterministic_Continuous_Dynamics::Deterministic_Continuous_Dynamics(const std::vector<Expression<Real> > & dynamics)
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

void Deterministic_Continuous_Dynamics::evaluate(Interval & result, const unsigned int varID, const TaylorModelVec<Real> & tmv_range, const std::vector<Interval> & domain, const unsigned int order, const Computational_Setting & setting)
{
	TaylorModel<Real> tm;
	expressions[varID].evaluate(tm, tmv_range.tms, order, domain, setting.tm_setting.cutoff_threshold, setting.g_setting);

	tm.intEval(result, domain);
}

int Deterministic_Continuous_Dynamics::reach(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, const Taylor_Model_Setting & tm_setting,
		const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
	std::vector<Constraint> dummy_invariant;

	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

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

	return checking_result;
}

int Deterministic_Continuous_Dynamics::reach_adaptive_stepsize(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, Taylor_Model_Setting & tm_setting,
		const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
	std::vector<Constraint> dummy_invariant;

//	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	Flowpipe newFlowpipe, currentFlowpipe = initialSet;
	double new_stepsize = -1;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = currentFlowpipe.advance_deterministic_adaptive_stepsize(newFlowpipe, expressions, new_stepsize, tm_setting, dummy_invariant, g_setting);

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

int Deterministic_Continuous_Dynamics::reach_adaptive_order(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, Taylor_Model_Setting & tm_setting,
		const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
	std::vector<Constraint> dummy_invariant;

	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

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


	return checking_result;
}

int Deterministic_Continuous_Dynamics::reach_symbolic_remainder(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, Symbolic_Remainder & symbolic_remainder, const Taylor_Model_Setting & tm_setting,
		const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
	std::vector<Constraint> dummy_invariant;

//	std::vector<Real> initial_scalars(expressions.size(), 1);

	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

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

			if(symbolic_remainder.J.size() >= symbolic_remainder.max_size)
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

	return checking_result;
}

int Deterministic_Continuous_Dynamics::reach_symbolic_remainder_adaptive_stepsize(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, Symbolic_Remainder & symbolic_remainder, Taylor_Model_Setting & tm_setting,
		const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
	std::vector<Constraint> dummy_invariant;

//	std::vector<Real> initial_scalars(expressions.size(), 1);

//	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

	tm_setting.setStepsize(tm_setting.step_max, tm_setting.order);
	double new_stepsize = -1;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = currentFlowpipe.advance_deterministic_adaptive_stepsize(newFlowpipe, expressions, new_stepsize, tm_setting, dummy_invariant, g_setting, symbolic_remainder);

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

			if(symbolic_remainder.J.size() >= symbolic_remainder.max_size)
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
/*
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
			*/

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

int Deterministic_Continuous_Dynamics::reach_symbolic_remainder_adaptive_order(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, Symbolic_Remainder & symbolic_remainder, Taylor_Model_Setting & tm_setting,
		const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
	std::vector<Constraint> dummy_invariant;

//	std::vector<Real> initial_scalars(expressions.size(), 1);

	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

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

			if(symbolic_remainder.J.size() >= symbolic_remainder.max_size)
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

	return checking_result;
}

void Deterministic_Continuous_Dynamics::reach(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const std::vector<Constraint> & unsafeSet) const
{
	bool bSafetyChecking = false;

	if(unsafeSet.size() > 0)
	{
		bSafetyChecking = true;
	}

	if(setting.tm_setting.step_min > 0)
	{
		// adaptive stepsizes
		result.status = reach_adaptive_stepsize(result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes, result.num_of_flowpipes,
				setting.time, initialSet, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);
	}
	else if(setting.tm_setting.order_max > 0)
	{
		// adaptive orders
		result.status = reach_adaptive_order(result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes, result.num_of_flowpipes,
				setting.time, initialSet, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);
	}
	else
	{
		// fixed stepsizes and orders
		result.status = reach(result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes, result.num_of_flowpipes,
				setting.time, initialSet, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);
	}


	if(result.nonlinear_flowpipes.size() > 0)
	{
		Flowpipe fpTmp = result.nonlinear_flowpipes.back();
		result.fp_end_of_time = fpTmp;

		fpTmp.tmvPre.evaluate_time(result.fp_end_of_time.tmvPre, setting.tm_setting.step_end_exp_table);
	}
}

void Deterministic_Continuous_Dynamics::reach_sr(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const std::vector<Constraint> & unsafeSet, Symbolic_Remainder & symbolic_remainder) const
{
	bool bSafetyChecking = false;

	if(unsafeSet.size() > 0)
	{
		bSafetyChecking = true;
	}

	if(setting.tm_setting.step_min > 0)
	{
		// adaptive stepsizes
		result.status = reach_symbolic_remainder_adaptive_stepsize(result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes, result.num_of_flowpipes,
				setting.time, initialSet, symbolic_remainder, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);
	}
	else if(setting.tm_setting.order_max > 0)
	{
		// adaptive orders
		result.status = reach_symbolic_remainder_adaptive_order(result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes, result.num_of_flowpipes,
				setting.time, initialSet, symbolic_remainder, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);
	}
	else
	{
		// fixed stepsizes and orders
		result.status = reach_symbolic_remainder(result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes, result.num_of_flowpipes,
				setting.time, initialSet, symbolic_remainder, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);
	}

	if(result.nonlinear_flowpipes.size() > 0)
	{
		Flowpipe fpTmp = result.nonlinear_flowpipes.back();
		result.fp_end_of_time = fpTmp;

		fpTmp.tmvPre.evaluate_time(result.fp_end_of_time.tmvPre, setting.tm_setting.step_end_exp_table);
	}
}


int Deterministic_Continuous_Dynamics::reach_inv(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<std::vector<Interval> > & domains, Flowpipe & last_flowpipe, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, const std::vector<Constraint> & invariant,
		const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking, const bool bPlot, const bool bTMOutput) const
{
	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	Flowpipe new_flowpipe, current_flowpipe = initialSet;
	bool bContracted = false;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = current_flowpipe.advance_deterministic(new_flowpipe, expressions, tm_setting, invariant, g_setting);

		if(res == 1)
		{
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
			case 2: 	// time interval is contracted
			{
				Real zero(0);

				if(contracted_domain[0].greaterThan(zero))
				{
					return checking_result;
				}
				else
				{
					++num_of_flowpipes;

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

int Deterministic_Continuous_Dynamics::reach_inv_adaptive_stepsize(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<std::vector<Interval> > & domains, Flowpipe & last_flowpipe, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, const std::vector<Constraint> & invariant,
		Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking, const bool bPlot, const bool bTMOutput) const
{
//	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	Flowpipe new_flowpipe, current_flowpipe = initialSet;
	bool bContracted = false;
	double new_stepsize = -1;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = current_flowpipe.advance_deterministic_adaptive_stepsize(new_flowpipe, expressions, new_stepsize, tm_setting, invariant, g_setting);

		if(res == 1)
		{
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
			case 2: 	// time interval is contracted
			{
				Real zero(0);

				if(contracted_domain[0].greaterThan(zero))
				{
					return checking_result;
				}
				else
				{
					++num_of_flowpipes;

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

			double current_stepsize = tm_setting.step_exp_table[1].sup();
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

int Deterministic_Continuous_Dynamics::reach_inv_adaptive_order(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<std::vector<Interval> > & domains, Flowpipe & last_flowpipe, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, const std::vector<Constraint> & invariant,
		Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking, const bool bPlot, const bool bTMOutput) const
{
	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	Flowpipe new_flowpipe, current_flowpipe = initialSet;
	bool bContracted = false;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = current_flowpipe.advance_deterministic_adaptive_order(new_flowpipe, expressions, tm_setting, invariant, g_setting);

		if(res == 1)
		{
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
			case 2: 	// time interval is contracted
			{
				Real zero(0);

				if(contracted_domain[0].greaterThan(zero))
				{
					return checking_result;
				}
				else
				{
					++num_of_flowpipes;

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

int Deterministic_Continuous_Dynamics::reach_inv_symbolic_remainder(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<std::vector<Interval> > & domains, Flowpipe & last_flowpipe, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, const std::vector<Constraint> & invariant, Symbolic_Remainder & symbolic_remainder,
		const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking, const bool bPlot, const bool bTMOutput) const
{
	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	Flowpipe new_flowpipe, current_flowpipe = initialSet;
	bool bContracted = false;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = current_flowpipe.advance_deterministic(new_flowpipe, expressions, tm_setting, invariant, g_setting, symbolic_remainder);

		if(res == 1)
		{
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
			case 2: 	// time interval is contracted
			{
				Real zero(0);

				if(contracted_domain[0].greaterThan(zero))
				{
					return checking_result;
				}
				else
				{
					++num_of_flowpipes;

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
				symbolic_remainder.reset(current_flowpipe);
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

int Deterministic_Continuous_Dynamics::reach_inv_symbolic_remainder_adaptive_stepsize(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<std::vector<Interval> > & domains, Flowpipe & last_flowpipe, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, const std::vector<Constraint> & invariant, Symbolic_Remainder & symbolic_remainder,
		Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking, const bool bPlot, const bool bTMOutput) const
{
	int checking_result = COMPLETED_SAFE;

	Flowpipe new_flowpipe, current_flowpipe = initialSet;
	bool bContracted = false;
	double new_stepsize = -1;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = current_flowpipe.advance_deterministic_adaptive_stepsize(new_flowpipe, expressions, new_stepsize, tm_setting, invariant, g_setting, symbolic_remainder);

		if(res == 1)
		{
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
			case 2: 	// time interval is contracted
			{
				Real zero(0);

				if(contracted_domain[0].greaterThan(zero))
				{
					return checking_result;
				}
				else
				{
					++num_of_flowpipes;

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

			double current_stepsize = tm_setting.step_exp_table[1].sup();
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
				symbolic_remainder.reset(current_flowpipe);
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

int Deterministic_Continuous_Dynamics::reach_inv_symbolic_remainder_adaptive_order(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<std::vector<Interval> > & domains, Flowpipe & last_flowpipe, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, const std::vector<Constraint> & invariant, Symbolic_Remainder & symbolic_remainder,
		Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking, const bool bPlot, const bool bTMOutput) const
{
	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	Flowpipe new_flowpipe, current_flowpipe = initialSet;
	bool bContracted = false;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = current_flowpipe.advance_deterministic_adaptive_order(new_flowpipe, expressions, tm_setting, invariant, g_setting, symbolic_remainder);

		if(res == 1)
		{
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
			case 2: 	// time interval is contracted
			{
				Real zero(0);

				if(contracted_domain[0].greaterThan(zero))
				{
					return checking_result;
				}
				else
				{
					++num_of_flowpipes;

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
				symbolic_remainder.reset(current_flowpipe);
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


void Deterministic_Continuous_Dynamics::reach_inv(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const std::vector<Constraint> & invariant, const std::vector<Constraint> & unsafeSet) const
{
	bool bSafetyChecking = false;

	if(unsafeSet.size() > 0)
	{
		bSafetyChecking = true;
	}

	if(setting.tm_setting.step_min > 0)
	{
		// adaptive stepsizes
		result.status = reach_inv_adaptive_stepsize(result.tmv_flowpipes, result.tmv_flowpipes_domains, result.fp_end_of_time, result.orders_of_flowpipes, result.safety_of_flowpipes,
				result.contraction_of_flowpipes, result.num_of_flowpipes, setting.time, initialSet, invariant, setting.tm_setting, setting.g_setting,
				setting.bPrint, unsafeSet, bSafetyChecking, true, true);
	}
	else if(setting.tm_setting.order_max > 0)
	{
		// adaptive orders
		result.status = reach_inv_adaptive_order(result.tmv_flowpipes, result.tmv_flowpipes_domains, result.fp_end_of_time, result.orders_of_flowpipes, result.safety_of_flowpipes,
				result.contraction_of_flowpipes, result.num_of_flowpipes, setting.time, initialSet, invariant, setting.tm_setting, setting.g_setting,
				setting.bPrint, unsafeSet, bSafetyChecking, true, true);
	}
	else
	{
		// fixed stepsizes and orders
		result.status = reach_inv(result.tmv_flowpipes, result.tmv_flowpipes_domains, result.fp_end_of_time, result.orders_of_flowpipes, result.safety_of_flowpipes,
				result.contraction_of_flowpipes, result.num_of_flowpipes, setting.time, initialSet, invariant, setting.tm_setting, setting.g_setting,
				setting.bPrint, unsafeSet, bSafetyChecking, true, true);
	}

	if(result.tmv_flowpipes.size() > 0)
	{
		TaylorModelVec<Real> last_flowpipe = result.tmv_flowpipes.back();

		Real t;
		(result.tmv_flowpipes_domains.back())[0].sup(t);

		std::vector<Real> realVec;
		realVec.push_back(1);
		realVec.push_back(t);

		Real tmp = t;

		for(unsigned int i=2; i<=setting.tm_setting.step_end_exp_table.size(); ++i)
		{
			tmp *= t;
			realVec.push_back(tmp);
		}

		last_flowpipe.evaluate_time(result.tmv_end_of_time, realVec);

		result.fp_end_of_time.tmvPre.evaluate_time(result.fp_end_of_time.tmvPre, realVec);
	}
}

void Deterministic_Continuous_Dynamics::reach_inv_sr(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const std::vector<Constraint> & invariant, const std::vector<Constraint> & unsafeSet, Symbolic_Remainder & symbolic_remainder) const
{
	bool bSafetyChecking = false;

	if(unsafeSet.size() > 0)
	{
		bSafetyChecking = true;
	}

	if(setting.tm_setting.step_min > 0)
	{
		// adaptive stepsizes
		result.status = reach_inv_symbolic_remainder_adaptive_stepsize(result.tmv_flowpipes, result.tmv_flowpipes_domains, result.fp_end_of_time, result.orders_of_flowpipes, result.safety_of_flowpipes,
				result.contraction_of_flowpipes, result.num_of_flowpipes, setting.time, initialSet, invariant, symbolic_remainder, setting.tm_setting, setting.g_setting,
				setting.bPrint, unsafeSet, bSafetyChecking, true, true);
	}
	else if(setting.tm_setting.order_max > 0)
	{
		// adaptive orders
		result.status = reach_inv_symbolic_remainder_adaptive_order(result.tmv_flowpipes, result.tmv_flowpipes_domains, result.fp_end_of_time, result.orders_of_flowpipes, result.safety_of_flowpipes,
				result.contraction_of_flowpipes, result.num_of_flowpipes, setting.time, initialSet, invariant, symbolic_remainder, setting.tm_setting, setting.g_setting,
				setting.bPrint, unsafeSet, bSafetyChecking, true, true);
	}
	else
	{
		// fixed stepsizes and orders
		result.status = reach_inv_symbolic_remainder(result.tmv_flowpipes, result.tmv_flowpipes_domains, result.fp_end_of_time, result.orders_of_flowpipes, result.safety_of_flowpipes,
				result.contraction_of_flowpipes, result.num_of_flowpipes, setting.time, initialSet, invariant, symbolic_remainder, setting.tm_setting, setting.g_setting,
				setting.bPrint, unsafeSet, bSafetyChecking, true, true);
	}

	if(result.tmv_flowpipes.size() > 0)
	{
		TaylorModelVec<Real> last_flowpipe = result.tmv_flowpipes.back();

		Real t;
		(result.tmv_flowpipes_domains.back())[0].sup(t);

		std::vector<Real> realVec;
		realVec.push_back(1);
		realVec.push_back(t);

		Real tmp = t;

		for(unsigned int i=2; i<=setting.tm_setting.step_end_exp_table.size(); ++i)
		{
			tmp *= t;
			realVec.push_back(tmp);
		}

		last_flowpipe.evaluate_time(result.tmv_end_of_time, realVec);



		result.fp_end_of_time.tmvPre.evaluate_time(result.fp_end_of_time.tmvPre, realVec);
	}
}





Nondeterministic_Continuous_Dynamics::Nondeterministic_Continuous_Dynamics(const std::vector<Expression<Interval> > & dynamics)
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

void Nondeterministic_Continuous_Dynamics::evaluate(Interval & result, const unsigned int varID, const TaylorModelVec<Real> & tmv_range, const std::vector<Interval> & domain, const unsigned int order, const Computational_Setting & setting)
{
	TaylorModel<Real> tm;
	expressions[varID].evaluate(tm, tmv_range.tms, order, domain, setting.tm_setting.cutoff_threshold, setting.g_setting);

	tm.intEval(result, domain);
}

int Nondeterministic_Continuous_Dynamics::reach(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, const Taylor_Model_Setting & tm_setting,
		const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
	std::vector<Constraint> dummy_invariant;

	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

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

	return checking_result;
}

int Nondeterministic_Continuous_Dynamics::reach_adaptive_stepsize(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, Taylor_Model_Setting & tm_setting,
		const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
	std::vector<Constraint> dummy_invariant;

//	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	Flowpipe newFlowpipe, currentFlowpipe = initialSet;
	double new_stepsize = -1;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = currentFlowpipe.advance_nondeterministic_adaptive_stepsize(newFlowpipe, expressions, new_stepsize, tm_setting, dummy_invariant, g_setting);

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
/*
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
			*/

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

int Nondeterministic_Continuous_Dynamics::reach_adaptive_order(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, Taylor_Model_Setting & tm_setting,
		const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
	std::vector<Constraint> dummy_invariant;

	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

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

	return checking_result;
}

int Nondeterministic_Continuous_Dynamics::reach_symbolic_remainder(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, const Taylor_Model_Setting & tm_setting, Symbolic_Remainder & symbolic_remainder,
		const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
	std::vector<Constraint> dummy_invariant;

	std::vector<Real> initial_scalars(expressions.size(), 1);

	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

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

			if(symbolic_remainder.J.size() >= symbolic_remainder.max_size)
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

	return checking_result;
}

int Nondeterministic_Continuous_Dynamics::reach_symbolic_remainder_adaptive_stepsize(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, Taylor_Model_Setting & tm_setting, Symbolic_Remainder & symbolic_remainder,
		const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
	std::vector<Constraint> dummy_invariant;

	std::vector<Real> initial_scalars(expressions.size(), 1);

//	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	Flowpipe newFlowpipe, currentFlowpipe = initialSet;
	double new_stepsize = -1;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = currentFlowpipe.advance_nondeterministic_adaptive_stepsize(newFlowpipe, expressions, new_stepsize, tm_setting, dummy_invariant, g_setting, symbolic_remainder);

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

			if(symbolic_remainder.J.size() >= symbolic_remainder.max_size)
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
/*
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
			*/

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

int Nondeterministic_Continuous_Dynamics::reach_symbolic_remainder_adaptive_order(std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, Taylor_Model_Setting & tm_setting,  Symbolic_Remainder & symbolic_remainder,
		const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking,
		const bool bPlot, const bool bTMOutput) const
{
	std::vector<Constraint> dummy_invariant;

	std::vector<Real> initial_scalars(expressions.size(), 1);

	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	Flowpipe newFlowpipe, currentFlowpipe = initialSet;

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

			if(symbolic_remainder.J.size() >= symbolic_remainder.max_size)
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

	return checking_result;
}

void Nondeterministic_Continuous_Dynamics::reach(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const std::vector<Constraint> & unsafeSet) const
{
	bool bSafetyChecking = false;

	if(unsafeSet.size() > 0)
	{
		bSafetyChecking = true;
	}

	if(setting.tm_setting.step_min > 0)
	{
		// adaptive stepsizes
		result.status = reach_adaptive_stepsize(result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes, result.num_of_flowpipes,
				setting.time, initialSet, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);
	}
	else if(setting.tm_setting.order_max > 0)
	{
		// adaptive orders
		result.status = reach_adaptive_order(result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes, result.num_of_flowpipes,
				setting.time, initialSet, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);
	}
	else
	{
		// fixed stepsizes and orders
		result.status = reach(result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes, result.num_of_flowpipes,
				setting.time, initialSet, setting.tm_setting, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);
	}


	if(result.nonlinear_flowpipes.size() > 0)
	{
		Flowpipe fpTmp = result.nonlinear_flowpipes.back();
		result.fp_end_of_time = fpTmp;

		fpTmp.tmvPre.evaluate_time(result.fp_end_of_time.tmvPre, setting.tm_setting.step_end_exp_table);
	}
}

void Nondeterministic_Continuous_Dynamics::reach_sr(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const std::vector<Constraint> & unsafeSet, Symbolic_Remainder & symbolic_remainder) const
{
	bool bSafetyChecking = false;

	if(unsafeSet.size() > 0)
	{
		bSafetyChecking = true;
	}

	if(setting.tm_setting.step_min > 0)
	{
		// adaptive stepsizes
		result.status = reach_symbolic_remainder_adaptive_stepsize(result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes, result.num_of_flowpipes,
				setting.time, initialSet, setting.tm_setting, symbolic_remainder, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);
	}
	else if(setting.tm_setting.order_max > 0)
	{
		// adaptive orders
		result.status = reach_symbolic_remainder_adaptive_order(result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes, result.num_of_flowpipes,
				setting.time, initialSet, setting.tm_setting, symbolic_remainder, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);
	}
	else
	{
		// fixed stepsizes and orders
		result.status = reach_symbolic_remainder(result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes, result.num_of_flowpipes,
				setting.time, initialSet, setting.tm_setting, symbolic_remainder, setting.g_setting, setting.bPrint, unsafeSet, bSafetyChecking, true, true);
	}

	if(result.nonlinear_flowpipes.size() > 0)
	{
		Flowpipe fpTmp = result.nonlinear_flowpipes.back();
		result.fp_end_of_time = fpTmp;

		fpTmp.tmvPre.evaluate_time(result.fp_end_of_time.tmvPre, setting.tm_setting.step_end_exp_table);
	}
}

int Nondeterministic_Continuous_Dynamics::reach_inv(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, const std::vector<Constraint> & invariant,
		const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking, const bool bPlot, const bool bTMOutput) const
{
	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	Flowpipe new_flowpipe, current_flowpipe = initialSet;
	bool bContracted = false;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = current_flowpipe.advance_nondeterministic(new_flowpipe, expressions, tm_setting, invariant, g_setting);

		if(res == 1)
		{
			// constrain the flowpipe by the invariant
			TaylorModelVec<Real> tmv_flowpipe;
			new_flowpipe.compose_normal(tmv_flowpipe, tm_setting.step_exp_table, tm_setting.order, tm_setting.cutoff_threshold);

			std::vector<Interval> contracted_domain = new_flowpipe.domain;
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
						flowpipes.push_back(new_flowpipe);
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
						flowpipes.push_back(new_flowpipe);
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
						flowpipes.push_back(new_flowpipe);
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
						flowpipes.push_back(new_flowpipe);
						flowpipes_safety.push_back(SAFE);
					}
				}

				current_flowpipe = new_flowpipe;

				break;
			}
			case 2: 	// time interval is contracted
			{
				Real zero(0);

				if(contracted_domain[0].greaterThan(zero))
				{
					return checking_result;
				}
				else
				{
					++num_of_flowpipes;

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
							flowpipes.push_back(new_flowpipe);
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
							flowpipes.push_back(new_flowpipe);
							flowpipes_safety.push_back(SAFE);
						}
					}

					t += contracted_domain[0].sup();

					if(bPrint)
					{
						printf("time = %f,\t", t);
						printf("step = %f,\t", contracted_domain[0].sup());
						printf("order = %d\n", tm_setting.order);
					}

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

	return checking_result;
}

int Nondeterministic_Continuous_Dynamics::reach_inv_adaptive_stepsize(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, const std::vector<Constraint> & invariant,
		Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking, const bool bPlot, const bool bTMOutput) const
{
	int checking_result = COMPLETED_SAFE;

	Flowpipe new_flowpipe, current_flowpipe = initialSet;
	bool bContracted = false;
	double new_stepsize = -1;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = current_flowpipe.advance_nondeterministic_adaptive_stepsize(new_flowpipe, expressions, new_stepsize, tm_setting, invariant, g_setting);

		if(res == 1)
		{
			// constrain the flowpipe by the invariant
			TaylorModelVec<Real> tmv_flowpipe;
			new_flowpipe.compose_normal(tmv_flowpipe, tm_setting.step_exp_table, tm_setting.order, tm_setting.cutoff_threshold);

			std::vector<Interval> contracted_domain = new_flowpipe.domain;
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
						flowpipes.push_back(new_flowpipe);
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
						flowpipes.push_back(new_flowpipe);
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
						flowpipes.push_back(new_flowpipe);
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
						flowpipes.push_back(new_flowpipe);
						flowpipes_safety.push_back(SAFE);
					}
				}

				current_flowpipe = new_flowpipe;

				break;
			}
			case 2: 	// time interval is contracted
			{
				Real zero(0);

				if(contracted_domain[0].greaterThan(zero))
				{
					return checking_result;
				}
				else
				{
					++num_of_flowpipes;

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
							flowpipes.push_back(new_flowpipe);
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
							flowpipes.push_back(new_flowpipe);
							flowpipes_safety.push_back(SAFE);
						}
					}

					t += contracted_domain[0].sup();

					if(bPrint)
					{
						printf("time = %f,\t", t);
						printf("step = %f,\t", contracted_domain[0].sup());
						printf("order = %d\n", tm_setting.order);
					}

					return checking_result;
				}
			}
			}

			double current_stepsize = tm_setting.step_exp_table[1].sup();
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

int Nondeterministic_Continuous_Dynamics::reach_inv_adaptive_order(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, const std::vector<Constraint> & invariant,
		Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking, const bool bPlot, const bool bTMOutput) const
{
	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	Flowpipe new_flowpipe, current_flowpipe = initialSet;
	bool bContracted = false;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = current_flowpipe.advance_nondeterministic_adaptive_order(new_flowpipe, expressions, tm_setting, invariant, g_setting);

		if(res == 1)
		{
			// constrain the flowpipe by the invariant
			TaylorModelVec<Real> tmv_flowpipe;
			new_flowpipe.compose_normal(tmv_flowpipe, tm_setting.step_exp_table, tm_setting.order, tm_setting.cutoff_threshold);

			std::vector<Interval> contracted_domain = new_flowpipe.domain;
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
						flowpipes.push_back(new_flowpipe);
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
						flowpipes.push_back(new_flowpipe);
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
						flowpipes.push_back(new_flowpipe);
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
						flowpipes.push_back(new_flowpipe);
						flowpipes_safety.push_back(SAFE);
					}
				}

				current_flowpipe = new_flowpipe;

				break;
			}
			case 2: 	// time interval is contracted
			{
				Real zero(0);

				if(contracted_domain[0].greaterThan(zero))
				{
					return checking_result;
				}
				else
				{
					++num_of_flowpipes;

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
							flowpipes.push_back(new_flowpipe);
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
							flowpipes.push_back(new_flowpipe);
							flowpipes_safety.push_back(SAFE);
						}
					}

					t += contracted_domain[0].sup();

					if(bPrint)
					{
						printf("time = %f,\t", t);
						printf("step = %f,\t", contracted_domain[0].sup());
						printf("order = %d\n", tm_setting.order);
					}

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

	return checking_result;
}

int Nondeterministic_Continuous_Dynamics::reach_inv_symbolic_remainder(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, const std::vector<Constraint> & invariant, Symbolic_Remainder & symbolic_remainder,
		const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking, const bool bPlot, const bool bTMOutput) const
{
	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	Flowpipe new_flowpipe, current_flowpipe = initialSet;
	bool bContracted = false;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = current_flowpipe.advance_nondeterministic(new_flowpipe, expressions, tm_setting, invariant, g_setting, symbolic_remainder);

		if(res == 1)
		{
			// constrain the flowpipe by the invariant
			TaylorModelVec<Real> tmv_flowpipe;
			new_flowpipe.compose_normal(tmv_flowpipe, tm_setting.step_exp_table, tm_setting.order, tm_setting.cutoff_threshold);

			std::vector<Interval> contracted_domain = new_flowpipe.domain;
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
						flowpipes.push_back(new_flowpipe);
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
						flowpipes.push_back(new_flowpipe);
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
						flowpipes.push_back(new_flowpipe);
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
						flowpipes.push_back(new_flowpipe);
						flowpipes_safety.push_back(SAFE);
					}
				}

				current_flowpipe = new_flowpipe;

				break;
			}
			case 2: 	// time interval is contracted
			{
				Real zero(0);

				if(contracted_domain[0].greaterThan(zero))
				{
					return checking_result;
				}
				else
				{
					++num_of_flowpipes;

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
							flowpipes.push_back(new_flowpipe);
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
							flowpipes.push_back(new_flowpipe);
							flowpipes_safety.push_back(SAFE);
						}
					}

					t += contracted_domain[0].sup();

					if(bPrint)
					{
						printf("time = %f,\t", t);
						printf("step = %f,\t", contracted_domain[0].sup());
						printf("order = %d\n", tm_setting.order);
					}

					return checking_result;
				}
			}
			}


			t += contracted_domain[0].sup();

			if(symbolic_remainder.J.size() >= symbolic_remainder.max_size)
			{
				symbolic_remainder.reset(current_flowpipe);
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

	return checking_result;
}

int Nondeterministic_Continuous_Dynamics::reach_inv_symbolic_remainder_adaptive_stepsize(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, const std::vector<Constraint> & invariant, Symbolic_Remainder & symbolic_remainder,
		Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking, const bool bPlot, const bool bTMOutput) const
{
	int checking_result = COMPLETED_SAFE;

	Flowpipe new_flowpipe, current_flowpipe = initialSet;
	bool bContracted = false;
	double new_stepsize = -1;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = current_flowpipe.advance_nondeterministic_adaptive_stepsize(new_flowpipe, expressions, new_stepsize, tm_setting, invariant, g_setting, symbolic_remainder);

		if(res == 1)
		{
			// constrain the flowpipe by the invariant
			TaylorModelVec<Real> tmv_flowpipe;
			new_flowpipe.compose_normal(tmv_flowpipe, tm_setting.step_exp_table, tm_setting.order, tm_setting.cutoff_threshold);

			std::vector<Interval> contracted_domain = new_flowpipe.domain;
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
						flowpipes.push_back(new_flowpipe);
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
						flowpipes.push_back(new_flowpipe);
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
						flowpipes.push_back(new_flowpipe);
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
						flowpipes.push_back(new_flowpipe);
						flowpipes_safety.push_back(SAFE);
					}
				}

				current_flowpipe = new_flowpipe;

				break;
			}
			case 2: 	// time interval is contracted
			{
				Real zero(0);

				if(contracted_domain[0].greaterThan(zero))
				{
					return checking_result;
				}
				else
				{
					++num_of_flowpipes;

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
							flowpipes.push_back(new_flowpipe);
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
							flowpipes.push_back(new_flowpipe);
							flowpipes_safety.push_back(SAFE);
						}
					}

					t += contracted_domain[0].sup();

					if(bPrint)
					{
						printf("time = %f,\t", t);
						printf("step = %f,\t", contracted_domain[0].sup());
						printf("order = %d\n", tm_setting.order);
					}

					return checking_result;
				}
			}
			}

			double current_stepsize = tm_setting.step_exp_table[1].sup();
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

int Nondeterministic_Continuous_Dynamics::reach_inv_symbolic_remainder_adaptive_order(std::list<TaylorModelVec<Real> > & tmv_flowpipes, std::list<Flowpipe> & flowpipes, std::list<unsigned int> & flowpipe_orders, std::list<int> & flowpipes_safety,
		std::list<bool> & contraction_of_flowpipes, unsigned long & num_of_flowpipes, const double time, const Flowpipe & initialSet, const std::vector<Constraint> & invariant, Symbolic_Remainder & symbolic_remainder,
		Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting, const bool bPrint, const std::vector<Constraint> & unsafeSet, const bool bSafetyChecking, const bool bPlot, const bool bTMOutput) const
{
	double step = tm_setting.step_exp_table[1].sup();

	int checking_result = COMPLETED_SAFE;

	Flowpipe new_flowpipe, current_flowpipe = initialSet;
	bool bContracted = false;

	for(double t=THRESHOLD_HIGH; t < time;)
	{
		int res = current_flowpipe.advance_nondeterministic_adaptive_order(new_flowpipe, expressions, tm_setting, invariant, g_setting, symbolic_remainder);

		if(res == 1)
		{
			// constrain the flowpipe by the invariant
			TaylorModelVec<Real> tmv_flowpipe;
			new_flowpipe.compose_normal(tmv_flowpipe, tm_setting.step_exp_table, tm_setting.order, tm_setting.cutoff_threshold);

			std::vector<Interval> contracted_domain = new_flowpipe.domain;
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
						flowpipes.push_back(new_flowpipe);
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
						flowpipes.push_back(new_flowpipe);
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
						flowpipes.push_back(new_flowpipe);
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
						flowpipes.push_back(new_flowpipe);
						flowpipes_safety.push_back(SAFE);
					}
				}

				current_flowpipe = new_flowpipe;

				break;
			}
			case 2: 	// time interval is contracted
			{
				Real zero(0);

				if(contracted_domain[0].greaterThan(zero))
				{
					return checking_result;
				}
				else
				{
					++num_of_flowpipes;

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
							flowpipes.push_back(new_flowpipe);
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
							flowpipes.push_back(new_flowpipe);
							flowpipes_safety.push_back(SAFE);
						}
					}

					t += contracted_domain[0].sup();

					if(bPrint)
					{
						printf("time = %f,\t", t);
						printf("step = %f,\t", contracted_domain[0].sup());
						printf("order = %d\n", tm_setting.order);
					}

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

	return checking_result;
}

void Nondeterministic_Continuous_Dynamics::reach_inv(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const std::vector<Constraint> & invariant, const std::vector<Constraint> & unsafeSet) const
{
	bool bSafetyChecking = false;

	if(unsafeSet.size() > 0)
	{
		bSafetyChecking = true;
	}

	if(setting.tm_setting.step_min > 0)
	{
		// adaptive stepsizes
		result.status = reach_inv_adaptive_stepsize(result.tmv_flowpipes, result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes,
				result.contraction_of_flowpipes, result.num_of_flowpipes, setting.time, initialSet, invariant, setting.tm_setting, setting.g_setting,
				setting.bPrint, unsafeSet, bSafetyChecking, true, true);
	}
	else if(setting.tm_setting.order_max > 0)
	{
		// adaptive orders
		result.status = reach_inv_adaptive_order(result.tmv_flowpipes, result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes,
				result.contraction_of_flowpipes, result.num_of_flowpipes, setting.time, initialSet, invariant, setting.tm_setting, setting.g_setting,
				setting.bPrint, unsafeSet, bSafetyChecking, true, true);
	}
	else
	{
		// fixed stepsizes and orders
		result.status = reach_inv(result.tmv_flowpipes, result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes,
				result.contraction_of_flowpipes, result.num_of_flowpipes, setting.time, initialSet, invariant, setting.tm_setting, setting.g_setting,
				setting.bPrint, unsafeSet, bSafetyChecking, true, true);
	}

	if(result.tmv_flowpipes.size() > 0)
	{
/*
		TaylorModelVec<Real> last_flowpipe = result.tmv_flowpipes.back();

		Real t;
		result.nonlinear_flowpipes.back().domain[0].sup(t);

		std::vector<Real> realVec;
		realVec.push_back(1);
		realVec.push_back(t);

		Real tmp = t;

		for(unsigned int i=2; i<=setting.tm_setting.step_end_exp_table.size(); ++i)
		{
			tmp *= t;
			realVec.push_back(tmp);
		}

		last_flowpipe.evaluate_time(result.tmv_end_of_time, realVec);
*/






		TaylorModelVec<Real> last_flowpipe = result.tmv_flowpipes.back();

		Real t;
		(result.tmv_flowpipes_domains.back())[0].sup(t);

		std::vector<Real> realVec;
		realVec.push_back(1);
		realVec.push_back(t);

		Real tmp = t;

		for(unsigned int i=2; i<=setting.tm_setting.step_end_exp_table.size(); ++i)
		{
			tmp *= t;
			realVec.push_back(tmp);
		}

		last_flowpipe.evaluate_time(result.tmv_end_of_time, realVec);


		result.fp_end_of_time.tmvPre.evaluate_time(result.fp_end_of_time.tmvPre, realVec);
	}
}

void Nondeterministic_Continuous_Dynamics::reach_inv_sr(Result_of_Reachability & result, Computational_Setting & setting, const Flowpipe & initialSet, const std::vector<Constraint> & invariant, const std::vector<Constraint> & unsafeSet, Symbolic_Remainder & symbolic_remainder) const
{
	bool bSafetyChecking = false;

	if(unsafeSet.size() > 0)
	{
		bSafetyChecking = true;
	}

	if(setting.tm_setting.step_min > 0)
	{
		// adaptive stepsizes
		result.status = reach_inv_symbolic_remainder_adaptive_stepsize(result.tmv_flowpipes, result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes,
				result.contraction_of_flowpipes, result.num_of_flowpipes, setting.time, initialSet, invariant, symbolic_remainder, setting.tm_setting, setting.g_setting,
				setting.bPrint, unsafeSet, bSafetyChecking, true, true);
	}
	else if(setting.tm_setting.order_max > 0)
	{
		// adaptive orders
		result.status = reach_inv_symbolic_remainder_adaptive_order(result.tmv_flowpipes, result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes,
				result.contraction_of_flowpipes, result.num_of_flowpipes, setting.time, initialSet, invariant, symbolic_remainder, setting.tm_setting, setting.g_setting,
				setting.bPrint, unsafeSet, bSafetyChecking, true, true);
	}
	else
	{
		// fixed stepsizes and orders
		result.status = reach_inv_symbolic_remainder(result.tmv_flowpipes, result.nonlinear_flowpipes, result.orders_of_flowpipes, result.safety_of_flowpipes,
				result.contraction_of_flowpipes, result.num_of_flowpipes, setting.time, initialSet, invariant, symbolic_remainder, setting.tm_setting, setting.g_setting,
				setting.bPrint, unsafeSet, bSafetyChecking, true, true);
	}

	if(result.tmv_flowpipes.size() > 0)
	{
/*
		TaylorModelVec<Real> last_flowpipe = result.tmv_flowpipes.back();

		Real t;
		result.nonlinear_flowpipes.back().domain[0].sup(t);

		std::vector<Real> realVec;
		realVec.push_back(1);
		realVec.push_back(t);

		Real tmp = t;

		for(unsigned int i=2; i<=setting.tm_setting.step_end_exp_table.size(); ++i)
		{
			tmp *= t;
			realVec.push_back(tmp);
		}

		last_flowpipe.evaluate_time(result.tmv_end_of_time, realVec);
*/




		TaylorModelVec<Real> last_flowpipe = result.tmv_flowpipes.back();

		Real t;
		(result.tmv_flowpipes_domains.back())[0].sup(t);

		std::vector<Real> realVec;
		realVec.push_back(1);
		realVec.push_back(t);

		Real tmp = t;

		for(unsigned int i=2; i<=setting.tm_setting.step_end_exp_table.size(); ++i)
		{
			tmp *= t;
			realVec.push_back(tmp);
		}

		last_flowpipe.evaluate_time(result.tmv_end_of_time, realVec);


		result.fp_end_of_time.tmvPre.evaluate_time(result.fp_end_of_time.tmvPre, realVec);
	}
}











Plot_Setting::Plot_Setting()
{
	type_of_file	= 0;
	type_of_object	= 0;
	num_of_pieces	= 0;
	bProjected		= false;
	bPrint			= true;
	bDiscrete		= false;
}

Plot_Setting::Plot_Setting(const Variables & vars)
{
	variables		= vars;
	type_of_file	= 0;
	type_of_object	= 0;
	num_of_pieces	= 0;
	bProjected		= false;
	bPrint			= true;
	bDiscrete		= false;
}

Plot_Setting::Plot_Setting(const Plot_Setting & setting)
{
	variables		= setting.variables;
	outputDims		= setting.outputDims;
	type_of_file	= setting.type_of_file;
	type_of_object	= setting.type_of_object;
	num_of_pieces	= setting.num_of_pieces;
	bProjected		= setting.bProjected;
	bPrint			= setting.bPrint;
	bDiscrete		= setting.bDiscrete;
}

Plot_Setting::~Plot_Setting()
{
}

Plot_Setting & Plot_Setting::operator = (const Plot_Setting & setting)
{
	if(this == &setting)
		return *this;

	variables		= setting.variables;
	outputDims		= setting.outputDims;
	type_of_file	= setting.type_of_file;
	type_of_object	= setting.type_of_object;
	num_of_pieces	= setting.num_of_pieces;
	bProjected		= setting.bProjected;
	bPrint			= setting.bPrint;

	return *this;
}

void Plot_Setting::setOutputDims(const std::string & x, const std::string & y)
{
	outputDims.clear();

	int x_id = variables.getIDForVar(x);

	if(x_id < 0)
	{
		std::cout << "The variable " << x << " is not declared!" << std::endl;
		exit(0);
	}

	int y_id = variables.getIDForVar(y);

	if(y_id < 0)
	{
		std::cout << "The variable " << y << " is not declared!" << std::endl;
		exit(0);
	}

	outputDims.push_back(x_id);
	outputDims.push_back(y_id);
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

void Plot_Setting::discreteOutput()
{
	bDiscrete = true;
}

void Plot_Setting::continuousOutput()
{
	bDiscrete = false;
}

void Plot_Setting::plot_2D(const std::string & path, const std::string & fileName, const Result_of_Reachability & result) const
{
	switch(type_of_file)
	{
	case PLOT_GNUPLOT:
		plot_2D_GNUPLOT(path, fileName, result);
		break;
	case PLOT_MATLAB:
		plot_2D_MATLAB(path, fileName, result);
		break;
	}
}

void Plot_Setting::plot_2D_MATLAB(const std::string & path, const std::string & fileName, const Result_of_Reachability & result) const
{
	switch(type_of_object)
	{
	case PLOT_INTERVAL:
		plot_2D_interval_MATLAB(path, fileName, result);
		break;
	case PLOT_OCTAGON:
		plot_2D_octagon_MATLAB(path, fileName, result);
		break;
	case PLOT_GRID:
		plot_2D_grids_MATLAB(path, fileName, num_of_pieces, result);
		break;
	}
}

void Plot_Setting::plot_2D_interval_MATLAB(const std::string & path, const std::string & fileName, const Result_of_Reachability & result) const
{
	if(bPrint)
	{
		printf("Generating the plot file...\n");
	}

	std::string name = path + fileName + ".m";
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
	std::list<std::vector<Interval> >::const_iterator domainIter = result.tmv_flowpipes_domains.begin();
	std::list<int>::const_iterator safetyIter = result.safety_of_flowpipes.begin();

	unsigned int total_size = result.safety_of_flowpipes.size();

	if(total_size > 0)
	{
		for(; safetyIter != result.safety_of_flowpipes.end() ; ++tmvIter, ++domainIter, ++safetyIter)
		{
			std::vector<Interval> box;

			if(bDiscrete)
			{
				std::vector<Interval> newDomain = *domainIter;
				newDomain[0] = (*domainIter)[0].sup();
				tmvIter->intEval(box, newDomain, varIDs);
			}
			else
			{
				tmvIter->intEval(box, *domainIter, varIDs);
			}

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

		for(; tmvIter != result.tmv_flowpipes.end() ; ++tmvIter, ++domainIter)
		{
			std::vector<Interval> box;

			if(bDiscrete)
			{
				std::vector<Interval> newDomain = *domainIter;
				newDomain[0] = (*domainIter)[0].sup();
				tmvIter->intEval(box, newDomain, varIDs);
			}
			else
			{
				tmvIter->intEval(box, *domainIter, varIDs);
			}

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

void Plot_Setting::plot_2D_octagon_MATLAB(const std::string & path, const std::string & fileName, const Result_of_Reachability & result) const
{
	if(bPrint)
	{
		printf("Generating the plot file...\n");
	}

	std::string name = path + fileName + ".m";
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
		rangeDim = variables.size();
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
	std::list<std::vector<Interval> >::const_iterator domainIter = result.tmv_flowpipes_domains.begin();
	std::list<int>::const_iterator safetyIter = result.safety_of_flowpipes.begin();

	unsigned int total_size = result.safety_of_flowpipes.size();

	if(total_size > 0)
	{
		for(; safetyIter != result.safety_of_flowpipes.end(); ++tmvIter, ++domainIter, ++safetyIter)
		{
			std::vector<Interval> newDomain = *domainIter;

			if(bDiscrete)
			{
				newDomain[0] = (*domainIter)[0].sup();
			}

			Polyhedron polyTemplate(sortedTemplate, *tmvIter, newDomain);

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

		for(; tmvIter != result.tmv_flowpipes.end(); ++tmvIter, ++domainIter)
		{
			std::vector<Interval> newDomain = *domainIter;

			if(bDiscrete)
			{
				newDomain[0] = (*domainIter)[0].sup();
			}

			Polyhedron polyTemplate(sortedTemplate, *tmvIter, newDomain);

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

void Plot_Setting::plot_2D_grids_MATLAB(const std::string & path, const std::string & fileName, const unsigned int num, const Result_of_Reachability & result) const
{
	if(bPrint)
	{
		printf("Generating the plot file...\n");
	}

	std::string name = path + fileName + ".m";
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
	std::list<std::vector<Interval> >::const_iterator domainIter = result.tmv_flowpipes_domains.begin();
	std::list<int>::const_iterator safetyIter = result.safety_of_flowpipes.begin();

	unsigned int total_size = result.safety_of_flowpipes.size();

	if(total_size > 0)
	{
		for(; safetyIter != result.safety_of_flowpipes.end(); ++tmvIter, ++domainIter, ++safetyIter)
		{
			// decompose the domain
			std::list<std::vector<Interval> > grids;

			if(bDiscrete)
			{
				std::vector<Interval> newDomain = *domainIter;
				newDomain[0] = (*domainIter)[0].sup();
				gridBox(grids, newDomain, num);
			}
			else
			{
				gridBox(grids, *domainIter, num);
			}

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

		for(; tmvIter != result.tmv_flowpipes.end(); ++tmvIter, ++domainIter)
		{
			// decompose the domain
			std::list<std::vector<Interval> > grids;

			if(bDiscrete)
			{
				std::vector<Interval> newDomain = *domainIter;
				newDomain[0] = (*domainIter)[0].sup();
				gridBox(grids, newDomain, num);
			}
			else
			{
				gridBox(grids, *domainIter, num);
			}


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

void Plot_Setting::plot_2D_GNUPLOT(const std::string & path, const std::string & fileName, const Result_of_Reachability & result) const
{
	switch(type_of_object)
	{
	case PLOT_INTERVAL:
		plot_2D_interval_GNUPLOT(path, fileName, result);
		break;
	case PLOT_OCTAGON:
		plot_2D_octagon_GNUPLOT(path, fileName, result);
		break;
	case PLOT_GRID:
		plot_2D_grids_GNUPLOT(path, fileName, num_of_pieces, result);
		break;
	}
}

void Plot_Setting::plot_2D_interval_GNUPLOT(const std::string & path, const std::string & fileName, const Result_of_Reachability & result) const
{
	if(bPrint)
	{
		printf("Generating the plot file...\n");
	}

	std::string plot_file_name = path + fileName + ".plt";
	FILE *plotFile = fopen(plot_file_name.c_str(), "w");

	if(plotFile == NULL)
	{
		printf("Can not create the output file.\n");
		exit(1);
	}

	std::string image_file_name = fileName + ".eps";

	fprintf(plotFile, "set terminal postscript enhanced color\n");

	fprintf(plotFile, "set output '%s'\n", image_file_name.c_str());

	fprintf(plotFile, "set style line 1 linecolor rgb \"blue\"\n");
	fprintf(plotFile, "set autoscale\n");
	fprintf(plotFile, "unset label\n");
	fprintf(plotFile, "set xtic auto\n");
	fprintf(plotFile, "set ytic auto\n");
	fprintf(plotFile, "set xlabel \"%s\"\n", variables.varNames[outputDims[0]].c_str());
	fprintf(plotFile, "set ylabel \"%s\"\n", variables.varNames[outputDims[1]].c_str());
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
	std::list<std::vector<Interval> >::const_iterator domainIter = result.tmv_flowpipes_domains.begin();
	std::list<int>::const_iterator safetyIter = result.safety_of_flowpipes.begin();

	unsigned int total_size = result.safety_of_flowpipes.size();

	if(total_size > 0)
	{
		for(; safetyIter != result.safety_of_flowpipes.end(); ++tmvIter, ++domainIter, ++safetyIter)
		{
			std::vector<Interval> box;

			if(bDiscrete)
			{
				std::vector<Interval> newDomain = *domainIter;
				newDomain[0] = (*domainIter)[0].sup();
				tmvIter->intEval(box, newDomain, varIDs);
			}
			else
			{
				tmvIter->intEval(box, *domainIter, varIDs);
			}


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

		for(; tmvIter != result.tmv_flowpipes.end(); ++tmvIter, ++domainIter)
		{
			std::vector<Interval> box;

			if(bDiscrete)
			{
				std::vector<Interval> newDomain = *domainIter;
				newDomain[0] = (*domainIter)[0].sup();
				tmvIter->intEval(box, newDomain, varIDs);
			}
			else
			{
				tmvIter->intEval(box, *domainIter, varIDs);
			}


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

void Plot_Setting::plot_2D_octagon_GNUPLOT(const std::string & path, const std::string & fileName, const Result_of_Reachability & result) const
{
	if(bPrint)
	{
		printf("Generating the plot file...\n");
	}

	std::string name = path + fileName + ".plt";
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
		rangeDim = variables.size();
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

	std::string image_file_name = fileName + ".eps";

	fprintf(plotFile, "set terminal postscript enhanced color\n");

	fprintf(plotFile, "set output '%s'\n", image_file_name.c_str());

	fprintf(plotFile, "set style line 1 linecolor rgb \"blue\"\n");
	fprintf(plotFile, "set autoscale\n");
	fprintf(plotFile, "unset label\n");
	fprintf(plotFile, "set xtic auto\n");
	fprintf(plotFile, "set ytic auto\n");
	fprintf(plotFile, "set xlabel \"%s\"\n", variables.varNames[outputDims[0]].c_str());
	fprintf(plotFile, "set ylabel \"%s\"\n", variables.varNames[outputDims[1]].c_str());
	fprintf(plotFile, "plot '-' notitle with lines ls 1\n");

	// Compute the intersections of two facets.
	// The vertices are ordered clockwisely.

	gsl_matrix *C = gsl_matrix_alloc(2,2);
	gsl_vector *d = gsl_vector_alloc(2);
	gsl_vector *vertex = gsl_vector_alloc(2);

	unsigned int prog = 0;

	std::list<TaylorModelVec<Real> >::const_iterator tmvIter = result.tmv_flowpipes.begin();
	std::list<std::vector<Interval> >::const_iterator domainIter = result.tmv_flowpipes_domains.begin();
	std::list<int>::const_iterator safetyIter = result.safety_of_flowpipes.begin();

	unsigned int total_size = result.safety_of_flowpipes.size();

	if(total_size > 0)
	{
		for(; safetyIter != result.safety_of_flowpipes.end(); ++tmvIter, ++domainIter, ++safetyIter)
		{
			std::vector<Interval> newDomain = *domainIter;

			if(bDiscrete)
			{
				newDomain[0] = (*domainIter)[0].sup();
			}

			Polyhedron polyTemplate(sortedTemplate, *tmvIter, newDomain);

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

		for(; tmvIter != result.tmv_flowpipes.end(); ++tmvIter, ++domainIter)
		{
			std::vector<Interval> newDomain = *domainIter;

			if(bDiscrete)
			{
				newDomain[0] = (*domainIter)[0].sup();
			}

			Polyhedron polyTemplate(sortedTemplate, *tmvIter, newDomain);

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

void Plot_Setting::plot_2D_grids_GNUPLOT(const std::string & path, const std::string & fileName, const unsigned int num, const Result_of_Reachability & result) const
{
	if(bPrint)
	{
		printf("Generating the plot file...\n");
	}

	std::string name = path + fileName + ".plt";
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

	std::string image_file_name = fileName + ".eps";

	fprintf(plotFile, "set terminal postscript enhanced color\n");

	fprintf(plotFile, "set output '%s'\n", image_file_name.c_str());

	fprintf(plotFile, "set style line 1 linecolor rgb \"blue\"\n");
	fprintf(plotFile, "set autoscale\n");
	fprintf(plotFile, "unset label\n");
	fprintf(plotFile, "set xtic auto\n");
	fprintf(plotFile, "set ytic auto\n");
	fprintf(plotFile, "set xlabel \"%s\"\n", variables.varNames[outputDims[0]].c_str());
	fprintf(plotFile, "set ylabel \"%s\"\n", variables.varNames[outputDims[1]].c_str());
	fprintf(plotFile, "plot '-' notitle with lines ls 1\n");

	unsigned int prog = 0;

	std::list<TaylorModelVec<Real> >::const_iterator tmvIter = result.tmv_flowpipes.begin();
	std::list<std::vector<Interval> >::const_iterator domainIter = result.tmv_flowpipes_domains.begin();
	std::list<int>::const_iterator safetyIter = result.safety_of_flowpipes.begin();

	unsigned int total_size = result.safety_of_flowpipes.size();

	if(total_size > 0)
	{
		for(; safetyIter != result.safety_of_flowpipes.end(); ++tmvIter, ++domainIter, ++safetyIter)
		{
			// decompose the domain
			std::list<std::vector<Interval> > grids;

			if(bDiscrete)
			{
				std::vector<Interval> newDomain = *domainIter;
				newDomain[0] = (*domainIter)[0].sup();
				gridBox(grids, newDomain, num);
			}
			else
			{
				gridBox(grids, *domainIter, num);
			}


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
		for(; tmvIter != result.tmv_flowpipes.end(); ++tmvIter, ++domainIter)
		{
			// decompose the domain
			std::list<std::vector<Interval> > grids;

			if(bDiscrete)
			{
				std::vector<Interval> newDomain = *domainIter;
				newDomain[0] = (*domainIter)[0].sup();
				gridBox(grids, newDomain, num);
			}
			else
			{
				gridBox(grids, *domainIter, num);
			}


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













// class CandidateVectorPool

CandidateVectorPool::CandidateVectorPool()
{
}

CandidateVectorPool::CandidateVectorPool(const std::vector<std::vector<double> > & candidates)
{
	vectors = candidates;
	unsigned int n = vectors.size();

	Matrix<double> wm(n, n);
	compute_weight_matrix(wm, vectors);
	weightMat = wm;
}

CandidateVectorPool::CandidateVectorPool(const CandidateVectorPool & pool)
{
	vectors = pool.vectors;
	weightMat = pool.weightMat;
}

CandidateVectorPool::~CandidateVectorPool()
{
}

unsigned int CandidateVectorPool::size() const
{
	return vectors.size();
}

CandidateVectorPool & CandidateVectorPool::operator = (const CandidateVectorPool & pool)
{
	if(this == &pool)
		return *this;

	vectors = pool.vectors;
	weightMat = pool.weightMat;
	return *this;
}





// class FactorTab

FactorTab::FactorTab()
{
	index		= 0;
	factor		= 0;
	intercept	= 0;
}

FactorTab::FactorTab(const int i, const double & f, const double & interc)
{
	index		= i;
	factor		= f;
	intercept	= interc;
}

FactorTab::~FactorTab()
{
}











Intersected_Flowpipes::Intersected_Flowpipes()
{
}

Intersected_Flowpipes::Intersected_Flowpipes(const Intersected_Flowpipes & intersections)
{
	flowpipes = intersections.flowpipes;
	domains = intersections.domains;
	start_t = intersections.start_t;
	durations = intersections.durations;
}

Intersected_Flowpipes::~Intersected_Flowpipes()
{
}

Intersected_Flowpipes & Intersected_Flowpipes::operator = (const Intersected_Flowpipes & intersections)
{
	if(this == &intersections)
		return *this;

	flowpipes = intersections.flowpipes;
	domains = intersections.domains;
	start_t = intersections.start_t;
	durations = intersections.durations;

	return *this;
}

void Intersected_Flowpipes::clear()
{
	flowpipes.clear();
	domains.clear();
	start_t.clear();
	durations.clear();
}

unsigned int Intersected_Flowpipes::numOfGroups() const
{
	return flowpipes.size();
}

void Intersected_Flowpipes::interval_aggregation(Flowpipe & result) const
{
	if(flowpipes.size() == 0)
		return;

	unsigned int rangeDim = flowpipes[0][0].tms.size();
	std::vector<Interval> int_aggregation(rangeDim);

	std::vector<Real> up(rangeDim, -UNBOUNDED);
	std::vector<Real> lo(rangeDim, UNBOUNDED);

	for(unsigned int k=0; k<flowpipes.size(); ++k)
	{
		std::vector<TaylorModelVec<Real> >::const_iterator tmvIter = flowpipes[k].begin();
		std::vector<std::vector<Interval> >::const_iterator doIter = domains[k].begin();

		for(; tmvIter != flowpipes[k].end(); ++tmvIter, ++doIter)
		{
			std::vector<Interval> box;
			tmvIter->intEval(box, *doIter);

			for(unsigned int i=0; i<rangeDim; ++i)
			{
				Real tmp1;
				box[i].inf(tmp1);

				if(tmp1 < lo[i])
					lo[i] = tmp1;

				Real tmp2;
				box[i].sup(tmp2);

				if(tmp2 > up[i])
					up[i] = tmp2;
			}
		}
	}

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		Interval I(lo[i], up[i], 0);
		int_aggregation[i] = I;
	}

	Flowpipe fp_aggregation(int_aggregation);
	result = fp_aggregation;
}

void Intersected_Flowpipes::interval_aggregation(std::vector<Flowpipe> & result) const
{
	if(flowpipes.size() == 0)
		return;

	unsigned int rangeDim = flowpipes[0][0].tms.size();

	for(unsigned int k=0; k<flowpipes.size(); ++k)
	{
		std::vector<Interval> int_aggregation(rangeDim);

		std::vector<Real> up(rangeDim, -UNBOUNDED);
		std::vector<Real> lo(rangeDim, UNBOUNDED);

		std::vector<TaylorModelVec<Real> >::const_iterator tmvIter = flowpipes[k].begin();
		std::vector<std::vector<Interval> >::const_iterator doIter = domains[k].begin();

		for(; tmvIter != flowpipes[k].end(); ++tmvIter, ++doIter)
		{
			std::vector<Interval> box;
			tmvIter->intEval(box, *doIter);

			for(unsigned int i=0; i<rangeDim; ++i)
			{
				Real tmp1;
				box[i].inf(tmp1);

				if(tmp1 < lo[i])
					lo[i] = tmp1;

				Real tmp2;
				box[i].sup(tmp2);

				if(tmp2 > up[i])
					up[i] = tmp2;
			}
		}

		for(unsigned int i=0; i<rangeDim; ++i)
		{
			Interval I(lo[i], up[i], 0);
			int_aggregation[i] = I;
		}

		Flowpipe fp_aggregation(int_aggregation);
		result.push_back(fp_aggregation);
	}
}

void Intersected_Flowpipes::parallelotope_aggregation(Flowpipe & result, CandidateVectorPool & pool) const
{
	// if there is no flowpipe to aggregate, then we return immediately
	if(flowpipes.size() == 0)
		return;


	unsigned int num_of_cand = pool.size();
	unsigned int rangeDim = flowpipes[0][0].tms.size();

	std::vector<Interval> ranges(rangeDim);
	Matrix<double> paraTemplate(rangeDim, rangeDim);
	std::vector<double> b(2*rangeDim);


	// if the number of candidate vectors are not enough, we also return immediately
	if(num_of_cand < rangeDim)
	{
		std::cout << "The number of candidate vectors are not enough. At least D vectors are required where D is the dimension of the state space." << std::endl;
		return;
	}
	else if(num_of_cand == rangeDim)	// no need to perform selection, since the number of candidate vectors is exactly enough
	{
		for(unsigned int k=0; k<flowpipes.size(); ++k)
		{
			for(unsigned int i=0; i<flowpipes[k].size(); ++i)
			{
				if(ranges.size() == 0)
				{
					for(unsigned int j=0; j<rangeDim; ++j)
					{
						Interval I;
						flowpipes[k][i].rho(I, pool.vectors[j], domains[k][i]);
						ranges.push_back(I);
					}
				}
				else
				{
					for(unsigned int j=0; j<rangeDim; ++j)
					{
						Interval I;
						flowpipes[k][i].rho(I, pool.vectors[j], domains[k][i]);
						ranges[j].hull_assign(I);
					}
				}
			}
		}

		for(unsigned int i=0; i<rangeDim; ++i)
			for(unsigned int j=0; j<rangeDim; ++j)
				paraTemplate[i][j] = pool.vectors[i][j];

		for(unsigned int i=0; i<rangeDim; ++i)
		{
			b[i] = ranges[i].sup();
			b[i+rangeDim] = -ranges[i].inf();
		}
	}
	else		// we select D vectors from the candidate vector pool
	{
		std::vector<double> rhoPos;
		std::vector<double> rhoNeg;

		std::list<FactorTab> lst_unselected;
		std::list<FactorTab> lst_selected;

		int num_selected = 0;

		for(unsigned int i=0; i<num_of_cand; ++i)
		{
			FactorTab facT(i, 0, INVALID);
			lst_unselected.push_back(facT);
			rhoPos.push_back(INVALID);
			rhoNeg.push_back(INVALID);
		}


		// compute the intercepts

		for(unsigned int k=0; k<flowpipes.size(); ++k)
		{
			for(unsigned int i=0; i<flowpipes[k].size(); ++i)
			{
				for(unsigned int j=0; j<num_of_cand; ++j)
				{
					Interval I;
					flowpipes[k][i].rho(I, pool.vectors[j], domains[k][i]);

					double tmp1 = I.sup();
					double tmp2 = -I.inf();

					if(tmp1 >= rhoPos[j])
					{
						rhoPos[j] = tmp1;
					}

					if(tmp2 >= rhoNeg[j])
					{
						rhoNeg[j] = tmp2;
					}
				}
			}
		}

		// select the vector with the smallest intercept
		std::list<FactorTab>::iterator facIter = lst_unselected.begin();
		std::list<FactorTab>::iterator min_facIter = facIter;

		facIter->intercept = rhoPos[0] + rhoNeg[0];
		++facIter;

		for(unsigned int i=1; facIter!=lst_unselected.end(); ++facIter, ++i)
		{
			facIter->intercept = rhoPos[i] + rhoNeg[i];

			if(facIter->intercept < min_facIter->intercept)
			{
				min_facIter = facIter;
			}
		}

		int index_selected = min_facIter->index;


		lst_selected.push_back(*min_facIter);

		b[num_selected] = rhoPos[index_selected];
		b[num_selected+rangeDim] = rhoNeg[index_selected];

		lst_unselected.erase(min_facIter);

		for(unsigned int i=0; i<rangeDim; ++i)
		{
			paraTemplate[num_selected][i] = pool.vectors[index_selected][i];
		}

		++num_selected;

		if(num_selected < rangeDim)
		{
			// update the factor table according to the selected vectors

			for(facIter=lst_unselected.begin(); facIter!=lst_unselected.end(); ++facIter)
			{
				facIter->factor = pool.weightMat[index_selected][facIter->index];
			}

			// consider the remaining vectors

			for(; lst_unselected.size() != 0 && num_selected < rangeDim; )
			{
				FactorTab vector_selected;
				bool bselected = select_a_vector(vector_selected, lst_unselected, paraTemplate, pool.vectors, num_selected);

				// update the factors
				if(bselected)
				{
					b[num_selected-1] = rhoPos[vector_selected.index];
					b[num_selected-1+rangeDim] = rhoNeg[vector_selected.index];

					lst_selected.push_back(vector_selected);

					facIter = lst_unselected.begin();
					for(; facIter!=lst_unselected.end(); ++facIter)
					{
						facIter->factor *= pool.weightMat[facIter->index][vector_selected.index];
					}
				}
				else
				{
					break;
				}
			}
		}
	}


	// we use the template parallelotope to over-approximate the flowpipe union

	Matrix<double> colVecCenter(rangeDim, 1);

	int d = paraTemplate.cols();

	gsl_vector *r = gsl_vector_alloc(d);
	for(int i=0; i<d; ++i)
		gsl_vector_set( r, i, (b[i] - b[i+rangeDim])/2 );

	// We use GSL to solve the linear equations B x = r.

	gsl_matrix *B = gsl_matrix_alloc(d, d);

	for(int i=0; i<d; ++i)
	{
		for(int j=0; j<d; ++j)
		{
			gsl_matrix_set(B, i, j, paraTemplate[i][j]);
		}
	}

	gsl_vector *x = gsl_vector_alloc(d);

	gsl_linalg_HH_solve(B, r, x);

	for(int i=0; i<d; ++i)
	{
		colVecCenter[i][0] = gsl_vector_get(x,i);
	}

	gsl_vector_free(r);
	gsl_matrix_free(B);
	gsl_vector_free(x);


	std::vector<Real> coefficients;
	for(int i=0; i<rangeDim; ++i)
	{
		coefficients.push_back(colVecCenter[i][0]);
	}

	TaylorModelVec<Real> tmvCenter(coefficients, rangeDim+1);

	// 2: we center the parallelotope at 0
	Matrix<double> colVecDiff = paraTemplate * colVecCenter;

	// since a parallelotope is symmetric, we only need to consider half of the intercepts
	Matrix<double> new_b(rangeDim, 1);

	for(unsigned int i=0; i<rangeDim; ++i)
	{
		new_b[i][0] = b[i] - colVecDiff[i][0];
	}

	// 3: compute the generators.
	Matrix<double> generators(rangeDim, rangeDim);
	std::vector<int> zeroRows;	// the row indices for zero intercepts

	for(int i=0; i<rangeDim; ++i)
	{
		if(new_b[i][0] <= THRESHOLD_LOW && new_b[i][0] >= -THRESHOLD_LOW)	// zero
		{
			zeroRows.push_back(i);

			for(int j=0; j<rangeDim; ++j)
			{
				generators[i][j] = paraTemplate[i][j];
			}
		}
		else
		{
			for(int j=0; j<rangeDim; ++j)
			{
				generators[i][j] = paraTemplate[i][j] / new_b[i][0];
			}
		}
	}

	generators.inverse_assign();

	Matrix<Real> tmv_coefficients(rangeDim, rangeDim+1);

	for(int j=0, k=0; j<rangeDim; ++j)
	{
		if(k < zeroRows.size() && j == zeroRows[k])	// neglect the zero length generators
		{
			++k;
		}
		else
		{
			for(int i=0; i<rangeDim; ++i)
			{
				tmv_coefficients[i][j+1] = generators[i][j];
			}
		}
	}

	TaylorModelVec<Real> tmvParallelotope(tmv_coefficients);
	tmvParallelotope += tmvCenter;

	Interval intUnit(-1,1);
	std::vector<Interval> intTmp(rangeDim + 1, intUnit);
	intTmp[0] = 0;

	Flowpipe fp_aggregation(tmvParallelotope, intTmp);
	result = fp_aggregation;
}

void Intersected_Flowpipes::parallelotope_aggregation(std::vector<Flowpipe> & result, CandidateVectorPool & pool) const
{
	// if there is no flowpipe to aggregate, then we return immediately
	if(flowpipes.size() == 0)
		return;


	unsigned int num_of_cand = pool.size();
	unsigned int rangeDim = flowpipes[0][0].tms.size();

	unsigned int num_of_groups = flowpipes.size();

	std::vector<std::vector<Interval> > ranges(num_of_groups, std::vector<Interval>(rangeDim));
	std::vector<Matrix<double> > paraTemplates(num_of_groups, Matrix<double>(rangeDim, rangeDim));
	std::vector<std::vector<double> > b(num_of_groups, std::vector<double>(2*rangeDim));


	// if the number of candidate vectors are not enough, we also return immediately
	if(num_of_cand < rangeDim)
	{
		std::cout << "The number of candidate vectors are not enough. At least D vectors are required where D is the dimension of the state space." << std::endl;
		return;
	}
	else if(num_of_cand == rangeDim)	// no need to perform selection, since the number of candidate vectors is exactly enough
	{
		for(unsigned int k=0; k<num_of_groups; ++k)
		{
			for(unsigned int i=0; i<flowpipes[k].size(); ++i)
			{
				if(ranges[k].size() == 0)
				{
					for(unsigned int j=0; j<rangeDim; ++j)
					{
						Interval I;
						flowpipes[k][i].rho(I, pool.vectors[j], domains[k][i]);
						ranges[k].push_back(I);
					}
				}
				else
				{
					for(unsigned int j=0; j<rangeDim; ++j)
					{
						Interval I;
						flowpipes[k][i].rho(I, pool.vectors[j], domains[k][i]);
						ranges[k][j].hull_assign(I);
					}
				}
			}
		}

		for(unsigned int i=0; i<rangeDim; ++i)
			for(unsigned int j=0; j<rangeDim; ++j)
				paraTemplates[0][i][j] = pool.vectors[i][j];

		for(unsigned int i=1; i<num_of_groups; ++i)
			paraTemplates[i] = paraTemplates[0];

		for(unsigned int k=0; k<num_of_groups; ++k)
		{
			for(unsigned int i=0; i<rangeDim; ++i)
			{
				b[k][i] = ranges[k][i].sup();
				b[k][i+rangeDim] = -ranges[k][i].inf();
			}
		}
	}
	else		// we select D vectors from the candidate vector pool
	{
		for(unsigned int k=0; k<num_of_groups; ++k)
		{
			std::vector<double> rhoPos;
			std::vector<double> rhoNeg;

			std::list<FactorTab> lst_unselected;
			std::list<FactorTab> lst_selected;

			int num_selected = 0;

			for(unsigned int i=0; i<num_of_cand; ++i)
			{
				FactorTab facT(i, 0, INVALID);
				lst_unselected.push_back(facT);
				rhoPos.push_back(INVALID);
				rhoNeg.push_back(INVALID);
			}


			// compute the intercepts

			for(unsigned int i=0; i<flowpipes[k].size(); ++i)
			{
				for(unsigned int j=0; j<num_of_cand; ++j)
				{
					Interval I;
					flowpipes[k][i].rho(I, pool.vectors[j], domains[k][i]);

					double tmp1 = I.sup();
					double tmp2 = -I.inf();

					if(tmp1 >= rhoPos[j])
					{
						rhoPos[j] = tmp1;
					}

					if(tmp2 >= rhoNeg[j])
					{
						rhoNeg[j] = tmp2;
					}
				}
			}

			// select the vector with the smallest intercept
			std::list<FactorTab>::iterator facIter = lst_unselected.begin();
			std::list<FactorTab>::iterator min_facIter = facIter;

			facIter->intercept = rhoPos[0] + rhoNeg[0];
			++facIter;

			for(unsigned int i=1; facIter!=lst_unselected.end(); ++facIter, ++i)
			{
				facIter->intercept = rhoPos[i] + rhoNeg[i];

				if(facIter->intercept < min_facIter->intercept)
				{
					min_facIter = facIter;
				}
			}

			int index_selected = min_facIter->index;


			lst_selected.push_back(*min_facIter);

			b[k][num_selected] = rhoPos[index_selected];
			b[k][num_selected+rangeDim] = rhoNeg[index_selected];

			lst_unselected.erase(min_facIter);

			for(unsigned int i=0; i<rangeDim; ++i)
			{
				paraTemplates[k][num_selected][i] = pool.vectors[index_selected][i];
			}

			++num_selected;

			if(num_selected < rangeDim)
			{
				// update the factor table according to the selected vectors

				for(facIter=lst_unselected.begin(); facIter!=lst_unselected.end(); ++facIter)
				{
					facIter->factor = pool.weightMat[index_selected][facIter->index];
				}

				// consider the remaining vectors

				for(; lst_unselected.size() != 0 && num_selected < rangeDim; )
				{
					FactorTab vector_selected;
					bool bselected = select_a_vector(vector_selected, lst_unselected, paraTemplates[k], pool.vectors, num_selected);

					// update the factors
					if(bselected)
					{
						b[k][num_selected-1] = rhoPos[vector_selected.index];
						b[k][num_selected-1+rangeDim] = rhoNeg[vector_selected.index];

						lst_selected.push_back(vector_selected);

						facIter = lst_unselected.begin();
						for(; facIter!=lst_unselected.end(); ++facIter)
						{
							facIter->factor *= pool.weightMat[facIter->index][vector_selected.index];
						}
					}
					else
					{
						break;
					}
				}
			}
		}
	}


	// we use the template parallelotope to over-approximate the flowpipe union

	for(unsigned int k=0; k<num_of_groups; ++k)
	{
		Matrix<double> colVecCenter(rangeDim, 1);

		int d = paraTemplates[k].cols();

		gsl_vector *r = gsl_vector_alloc(d);
		for(int i=0; i<d; ++i)
			gsl_vector_set( r, i, (b[k][i] - b[k][i+rangeDim])/2 );

		// We use GSL to solve the linear equations B x = r.

		gsl_matrix *B = gsl_matrix_alloc(d, d);

		for(int i=0; i<d; ++i)
		{
			for(int j=0; j<d; ++j)
			{
				gsl_matrix_set(B, i, j, paraTemplates[k][i][j]);
			}
		}

		gsl_vector *x = gsl_vector_alloc(d);

		gsl_linalg_HH_solve(B, r, x);

		for(int i=0; i<d; ++i)
		{
			colVecCenter[i][0] = gsl_vector_get(x,i);
		}

		gsl_vector_free(r);
		gsl_matrix_free(B);
		gsl_vector_free(x);


		std::vector<Real> coefficients;
		for(int i=0; i<rangeDim; ++i)
		{
			coefficients.push_back(colVecCenter[i][0]);
		}

		TaylorModelVec<Real> tmvCenter(coefficients, rangeDim+1);

		// 2: we center the parallelotope at 0
		Matrix<double> colVecDiff = paraTemplates[k] * colVecCenter;

		// since a parallelotope is symmetric, we only need to consider half of the intercepts
		Matrix<double> new_b(rangeDim, 1);

		for(unsigned int i=0; i<rangeDim; ++i)
		{
			new_b[i][0] = b[k][i] - colVecDiff[i][0];
		}

		// 3: compute the generators.
		Matrix<double> generators(rangeDim, rangeDim);
		std::vector<int> zeroRows;	// the row indices for zero intercepts

		for(int i=0; i<rangeDim; ++i)
		{
			if(new_b[i][0] <= THRESHOLD_LOW && new_b[i][0] >= -THRESHOLD_LOW)	// zero
			{
				zeroRows.push_back(i);

				for(int j=0; j<rangeDim; ++j)
				{
					generators[i][j] = paraTemplates[k][i][j];
				}
			}
			else
			{
				for(int j=0; j<rangeDim; ++j)
				{
					generators[i][j] = paraTemplates[k][i][j] / new_b[i][0];
				}
			}
		}

		generators.inverse_assign();

		Matrix<Real> tmv_coefficients(rangeDim, rangeDim+1);

		for(int j=0, k=0; j<rangeDim; ++j)
		{
			if(k < zeroRows.size() && j == zeroRows[k])	// neglect the zero length generators
			{
				++k;
			}
			else
			{
				for(int i=0; i<rangeDim; ++i)
				{
					tmv_coefficients[i][j+1] = generators[i][j];
				}
			}
		}

		TaylorModelVec<Real> tmvParallelotope(tmv_coefficients);
		tmvParallelotope += tmvCenter;

		Interval intUnit(-1,1);
		std::vector<Interval> intTmp(rangeDim + 1, intUnit);
		intTmp[0] = 0;

		Flowpipe fp_aggregation(tmvParallelotope, intTmp);
		result.push_back(fp_aggregation);
	}
}








namespace flowstar
{

bool compareFactor(const FactorTab & a, const FactorTab & b)
{
	if(a.factor > b.factor)
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool compareIntercept(const FactorTab & a, const FactorTab & b)
{
	if(a.intercept < b.intercept)
	{
		return true;
	}
	else
	{
		return false;
	}
}


int safetyChecking(const TaylorModelVec<Real> & tmv, const std::vector<Interval> & domain, const std::vector<Constraint> & unsafeSet, const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting)
{
	if(unsafeSet.size() == 0)
	{
		return SAFE;
	}

//	unsigned int rangeDim = tmv.tms.size();
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
			if(domain[0].width() <= REFINEMENT_PREC)
				return UNKNOWN;


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




int remainder_contraction_int(const std::vector<Interval> & polyRange, std::vector<Interval> & remainders, const std::vector<Constraint> & constraints)
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
//	unsigned int domainDim = rangeDim + 1;

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

int domain_contraction_int(const TaylorModelVec<Real> & tmv_flowpipe, std::vector<Interval> & domain, const std::vector<Constraint> & constraints, const unsigned int order, const Interval & cutoff_threshold, const Global_Setting & g_setting)
{
	if(constraints.size() == 0)
	{
		return 0;
	}

//	int rangeDim = tmv_flowpipe.tms.size();
	int domainDim = domain.size();

	// the Horner forms of p(T(x))
	std::vector<HornerForm<Real> > objHF;

	std::vector<Interval> hf_remainders;

	/*
	 * Given a Taylor model flowpipe x = p(x0) + I and a constraint q(x) <= b,
	 * we firstly translate the function q(p(x0) + I) to a Horner form.
	 */

	std::vector<bool> bChecking(constraints.size(), true);	// used to track the necessity to check a constraint
	unsigned int counter = 0; // counting the number of constraints which are no need to check

	bool bValid = true;

	for(int i=0; i<constraints.size(); ++i)
	{
		TaylorModel<Real> tmTemp;
		constraints[i].expression.evaluate(tmTemp, tmv_flowpipe.tms, order, domain, cutoff_threshold, g_setting);

		HornerForm<Real> hf;
		Interval remainder;
		tmTemp.toHornerForm(hf, remainder);
		objHF.push_back(hf);
		hf_remainders.push_back(remainder);

		// perform a preliminary checking
		Interval I;
		hf.evaluate(I, domain);
		I += remainder;

		if(I.greaterThan(constraints[i].bound))
		{
			// the constraint is not satisfied
			bValid = false;
			break;
		}
		else if(I.lessThanEq(constraints[i].bound))
		{
			// the constraint is satisfied and will not be used for contracting the domain
			bChecking[i] = false;
			++counter;
		}
		else
		{
			bChecking[i] = true;
			continue;
		}
	}

	if(!bValid)
	{
		return UNSAT;	// at least one of the constraints is not satisfied
	}
	else if(counter == constraints.size())
	{
		return SAT;		// all constraints are satisfied, no need to contract the domain
	}


	Interval intTime = domain[0];

	bool bContinue = true;

	/*
	 * Domain contraction:
	 * In every iteration, we look for a safe lower bound and a safe upper bound
	 * for each domain dimension. The loop terminates until there is no great improvement
	 * on the domain.
	 */

	for(; bContinue; )
	{
		std::vector<Interval> original_domain = domain;

		// contract the domain
		for(unsigned int i=0; i<domainDim; ++i)
		{
			Interval intRange = domain[i];
			std::vector<bool> bChecking_local = bChecking;
			unsigned int localCounter = counter;
			double w = intRange.width();

			// looking for a safe lower bound
			for(; w > DC_THRESHOLD_SEARCH;)
			{
				Interval intLeft;
				Interval intRight;
				intRange.split(intLeft, intRight);

				for(unsigned int j=0; j<constraints.size(); ++j)
				{
					if(bChecking_local[j])
					{
						std::vector<Interval> tmpDomain = domain;
						tmpDomain[i] = intLeft;

						Interval I;
						objHF[j].evaluate(I, tmpDomain);
						I += hf_remainders[j];

						if(I.greaterThan(constraints[j].bound))
						{
							// the constraint is not satisfied by the left half
							intRange = intRight;
							w = intRange.width();
							break;
						}
						else if(I.lessThanEq(constraints[j].bound))
						{
							// the constraint is satisfied by the left half
							intRange = intLeft;
							w = intRange.width();
							bChecking_local[j] = false;
							++localCounter;
						}
						else
						{
							// refine the interval
							intRange = intLeft;
							w = intRange.width();
						}
					}
				}

				if(localCounter == constraints.size())
				{
					break;
				}
			}

			// set the lower bound
			Real lo;
			intRange.inf(lo);
			domain[i].setInf(lo);


			intRange = domain[i];
			bChecking_local = bChecking;
			localCounter = counter;
			w = intRange.width();


			// looking for a safe upper bound
			for(; w > DC_THRESHOLD_SEARCH;)
			{
				Interval intLeft;
				Interval intRight;
				intRange.split(intLeft, intRight);

				for(unsigned int j=0; j<constraints.size(); ++j)
				{
					if(bChecking_local[j])
					{
						std::vector<Interval> tmpDomain = domain;
						tmpDomain[i] = intRight;

						Interval I;
						objHF[j].evaluate(I, tmpDomain);
						I += hf_remainders[j];

						if(I.greaterThan(constraints[j].bound))
						{
							// the constraint is not satisfied by the right half
							intRange = intLeft;
							w = intRange.width();
							break;
						}
						else if(I.lessThanEq(constraints[j].bound))
						{
							// the constraint is satisfied by the right half
							intRange = intRight;
							w = intRange.width();
							bChecking_local[j] = false;
							++localCounter;
						}
						else
						{
							// refine the interval
							intRange = intRight;
							w = intRange.width();
						}
					}
				}

				if(localCounter == constraints.size())
				{
					break;
				}
			}

			// set the upper bound
			Real up;
			intRange.sup(up);
			domain[i].setSup(up);

			if(!domain[i].valid())
			{
				bValid = false;
				break;
			}
		}

		if(!bValid)
		{
			break;
		}

		bContinue = false;

		for(int i=0; i<domainDim; ++i)
		{
			if(original_domain[i].widthRatio(domain[i]) <= DC_THRESHOLD_IMPROV)
			{
				bContinue = true;
				break;
			}
		}

		if(bContinue)
		{
			for(unsigned int i=0; i<constraints.size(); ++i)
			{
				TaylorModel<Real> tmTemp;

				constraints[i].expression.evaluate(tmTemp, tmv_flowpipe.tms, order, domain, cutoff_threshold, g_setting);

				tmTemp.toHornerForm(objHF[i], hf_remainders[i]);
			}
		}
	}

	if(!bValid)
	{
		return UNSAT;
	}

	// checking the satisfiability on the contracted domain
	for(unsigned int i=0; i<constraints.size(); ++i)
	{
		Interval I;
		objHF[i].evaluate(I, domain);
		I += hf_remainders[i];

		if(I.greaterThan(constraints[i].bound))
			return UNSAT;
	}

	if(intTime != domain[0])
	{
		return TIME_RANGE_CONTRACTED;
	}
	else
	{
		return CONTRACTED;
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

void intersect_a_guard(Intersected_Flowpipes & result, const Result_of_Reachability & flowpipes, const std::vector<Constraint> & guard, const bool boundary_of_invariant, const Computational_Setting & setting)
{
	result.clear();

	Real start_time, duration, g_time;
	bool newSection = true;
	bool firstIntersection = true;

	std::list<TaylorModelVec<Real> >::const_iterator tmvIter = flowpipes.tmv_flowpipes.begin();
	std::list<std::vector<Interval> >::const_iterator doIter = flowpipes.tmv_flowpipes_domains.begin();
	std::list<unsigned int>::const_iterator orderIter = flowpipes.orders_of_flowpipes.begin();
	std::list<bool>::const_iterator contractionIter = flowpipes.contraction_of_flowpipes.begin();

	std::vector<TaylorModelVec<Real> > tmvVec;
	std::vector<std::vector<Interval> > domainVec;


	for(; tmvIter != flowpipes.tmv_flowpipes.end(); ++tmvIter, ++doIter, ++orderIter, ++contractionIter)
	{
		if(boundary_of_invariant && !(*contractionIter))
		{
			Real r;
			(*doIter)[0].sup(r);
			g_time += r;
			continue;
		}

		std::vector<Interval> contracted_domain = *doIter;
		int res = domain_contraction_int(*tmvIter, contracted_domain, guard, *orderIter, setting.tm_setting.cutoff_threshold, setting.g_setting);

		if(res != UNSAT)
		{
			tmvVec.push_back(*tmvIter);
			domainVec.push_back(contracted_domain);

			Real r;

			if(firstIntersection)
			{
				firstIntersection = false;
				newSection = false;

				contracted_domain[0].inf(r);
				start_time = g_time + r;
			}

			Real t;
			contracted_domain[0].width(t);
			duration += t;

			(*doIter)[0].sup(r);
			g_time += r;
		}
		else
		{
			if(!newSection)
			{
				result.flowpipes.push_back(tmvVec);
				result.domains.push_back(domainVec);
				result.start_t.push_back(start_time);
				result.durations.push_back(duration);

				duration = 0;
				tmvVec.clear();
				domainVec.clear();

				firstIntersection = true;
				newSection = true;
			}

			Real r;
			(*doIter)[0].sup(r);
			g_time += r;
		}
	}

	if(tmvVec.size() > 0)
	{
		result.flowpipes.push_back(tmvVec);
		result.domains.push_back(domainVec);
		result.start_t.push_back(start_time);
		result.durations.push_back(duration);
	}
}


void compute_weight_matrix(Matrix<double> & weightMat, const std::vector<std::vector<double> > & candidate_vectors)
{
	unsigned int matrix_size = candidate_vectors.size();
	unsigned int vec_size = candidate_vectors[0].size();

	for(unsigned int i=0; i<matrix_size; ++i)
	{
		for(unsigned int j=i+1; j<matrix_size; ++j)
		{
			double d = 1.0;

			for(unsigned int k=0; k<vec_size; ++k)
			{
				d -= candidate_vectors[i][k] * candidate_vectors[j][k];
			}

			weightMat[i][j] = d;
			weightMat[j][i] = d;
		}
	}
}

bool check_validity(Matrix<double> & matTemplate, const std::vector<double> & vec, const int rank)
{
	unsigned int num = vec.size();

	for(unsigned int i=0; i<num; ++i)
	{
		matTemplate[rank][i] = vec[i];
	}

	int r = matTemplate.rank();

	if(r == rank+1)
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool select_a_vector(FactorTab & lst_selected, std::list<FactorTab> & lst_unselected, Matrix<double> & matTemplate, const std::vector<std::vector<double> > & candidate_vectors, int & rank)
{
	lst_unselected.sort(compareFactor);

	std::list<FactorTab> candidates;
	bool bvalid = false;


/*
	// =============== test begin ==================
	printf("Candidates:\n");
	list<FactorTab>::iterator testIter = lst_unselected.begin();
	for(; testIter!=lst_unselected.end(); ++testIter)
	{
		printf("vector: ");
		rowVecs[testIter->index].dump(stdout);
		printf("\tintercept: %lf, factor: %lf\n\n", testIter->intercept.midpoint(), testIter->factor.midpoint());
	}
	// =============== test end ==================
*/


	Interval intZero;

	for(; !bvalid && lst_unselected.size()!=0;)
	{
		std::list<FactorTab>::iterator facIter = lst_unselected.begin();

		double factor = facIter->factor;
		candidates.push_back(*facIter);
		facIter = lst_unselected.erase(facIter);

		for(; facIter!=lst_unselected.end(); )
		{
			if(facIter->factor < 0)
			{
				facIter = lst_unselected.erase(facIter);
			}
			else if(facIter->factor >= factor - THRESHOLD_LOW && facIter->factor <= factor + THRESHOLD_LOW)
			{
				candidates.push_back(*facIter);
				facIter = lst_unselected.erase(facIter);
			}
			else
			{
				break;
			}
		}

		for(; !bvalid && candidates.size()!=0; )
		{
			facIter = candidates.begin();
			double min_intercept = facIter->intercept;
			std::list<FactorTab>::iterator iter_selected = facIter;

			++facIter;

			for(; facIter!=candidates.end(); ++facIter)
			{
				if(min_intercept > facIter->intercept)
				{
					min_intercept = facIter->intercept;
					iter_selected = facIter;
				}
			}

			lst_selected = *iter_selected;
			candidates.erase(iter_selected);

			bvalid = check_validity(matTemplate, candidate_vectors[lst_selected.index], rank);
		}
	}

	if(bvalid)
	{
		++rank;

		// insert the unselected elements back
		std::list<FactorTab>::iterator facIter = candidates.begin();

		for(; facIter!=candidates.end(); ++facIter)
		{
			lst_unselected.push_back(*facIter);
		}


/*
		// =============== test begin ==================
		printf("selected: ");
		rowVecs[lst_selected.index].dump(stdout);
		printf("\tintercept: %lf, factor: %lf\n\n\n", lst_selected.intercept.midpoint(), lst_selected.factor.midpoint());
		// =============== test end ==================
*/




		return true;	// one element is selected
	}
	else
	{
		return false;	// nothing in the unselected list is selectable
	}
}


void eliminate_t(TaylorModelVec<Real> & tmv_flowpipe, std::vector<Interval> & domain)
{
	if(domain[0].isZero())
	{
		for(unsigned int i=0; i<tmv_flowpipe.tms.size(); ++i)
		{
			std::list<Term<Real> >::iterator term_iter = tmv_flowpipe.tms[i].expansion.terms.begin();

			Polynomial<Real> result;

			for(; term_iter != tmv_flowpipe.tms[i].expansion.terms.end(); ++term_iter)
			{
				if(term_iter->degrees[0] == 0)
				{
					result += *term_iter;
				}
			}

			tmv_flowpipe.tms[i].expansion = result;
		}
	}
	else
	{
		std::vector<Interval> newDomain = domain;

		// translate the time variable to the new added variable
		newDomain.push_back(domain[0]);
		newDomain[0] = 0;

		for(unsigned int i=0; i<tmv_flowpipe.tms.size(); ++i)
		{
			std::list<Term<Real> >::iterator term_iter = tmv_flowpipe.tms[i].expansion.terms.begin();

			Polynomial<Real> translation;

			for(; term_iter != tmv_flowpipe.tms[i].expansion.terms.end(); ++term_iter)
			{
				Term<Real> term = *term_iter;
				term.degrees.push_back(term.degrees[0]);
				term.degrees[0] = 0;

				translation += term;
			}

			tmv_flowpipe.tms[i].expansion = translation;
		}

		domain = newDomain;
	}
}

void create_initial_set(Flowpipe & initial_set, TaylorModelVec<Real> & tmv_flowpipe, std::vector<Interval> & domain, const unsigned int order, const Interval & cutoff_threshold)
{
	if(!tmv_flowpipe.isFreeOfT())
		eliminate_t(tmv_flowpipe, domain);

	Flowpipe fp(tmv_flowpipe, domain, cutoff_threshold);

	initial_set = fp;
}

void merge_consecutive_flowpipes(Flowpipe & result, std::vector<TaylorModelVec<Real> > & flowpipes, std::vector<std::vector<Interval> > & flowpipe_domains, Linear_Time_Invariant_Dynamics & LTI_dynamics,
		const std::vector<Constraint> & invariant, const unsigned int order, const Computational_Setting & setting)
{
	Real startTime;
	flowpipe_domains[0][0].inf(startTime);

	std::vector<Real> time_exp_table;
	time_exp_table.push_back(1);
	time_exp_table.push_back(startTime);

	Real tmp = startTime;

	for(int i=0; i<=order; ++i)
	{
		tmp *= startTime;
		time_exp_table.push_back(tmp);
	}

	// computing the initial set
	TaylorModelVec<Real> initialSet;
	flowpipes[0].evaluate_time(initialSet, time_exp_table);


	// evaluating the time duration
	double t = flowpipe_domains[0][0].width();

	for(int i=1; i<flowpipe_domains.size(); ++i)
	{
		t += flowpipe_domains[i][0].width();
	}

	std::vector<Interval> domain = flowpipe_domains[0];
	domain[0] = 0;

	for(int i=1; i<flowpipe_domains.size(); ++i)
	{
		for(int j=1; j<flowpipe_domains[i].size(); ++j)
		{
			domain[j].hull_assign(flowpipe_domains[i][j]);
		}
	}

	TaylorModelVec<Real> merged_flowpipe;
	LTI_dynamics.reach_one_step(merged_flowpipe, t, initialSet, domain, order, setting.tm_setting.cutoff_threshold);

	domain_contraction_int(merged_flowpipe, domain, invariant, order, setting.tm_setting.cutoff_threshold, setting.g_setting);

	create_initial_set(result, merged_flowpipe, domain, order, setting.tm_setting.cutoff_threshold);
}










// only for testing
void test_domain_contraction(Result_of_Reachability & contraction_result, Result_of_Reachability & reachability_result, const std::vector<Constraint> & constraints, const Taylor_Model_Setting & tm_setting, const Global_Setting & g_setting)
{
	std::list<TaylorModelVec<Real> >::iterator tmvIter = reachability_result.tmv_flowpipes.begin();
	std::list<Flowpipe>::const_iterator fpIter = reachability_result.nonlinear_flowpipes.begin();
	std::list<int>::const_iterator safetyIter = reachability_result.safety_of_flowpipes.begin();
	std::list<unsigned int>::const_iterator orderIter = reachability_result.orders_of_flowpipes.begin();


	for(; safetyIter != reachability_result.safety_of_flowpipes.end(); ++tmvIter, ++fpIter, ++safetyIter, ++orderIter)
	{
		std::vector<Interval> domain = fpIter->domain;
		int result = domain_contraction_int(*tmvIter, domain, constraints, *orderIter, tm_setting.cutoff_threshold, g_setting);

		if(result != UNSAT)
		{
			Flowpipe fpTemp = *fpIter;
			fpTemp.domain = domain;
			contraction_result.nonlinear_flowpipes.push_back(fpTemp);

			contraction_result.tmv_flowpipes.push_back(*tmvIter);
			contraction_result.orders_of_flowpipes.push_back(*orderIter);
			contraction_result.safety_of_flowpipes.push_back(*safetyIter);
			contraction_result.num_of_flowpipes++;
		}
	}
}

}
















