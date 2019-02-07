/*---
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.
  The code is released as is under the GNU General Public License (GPL).
---*/

#ifndef EXPRESSION_H_
#define EXPRESSION_H_

#include "TaylorModel.h"

#define OPT_PLUS		0
#define OPT_MINU		1
#define OPT_MULT		2
#define OPT_DIV			3
#define OPT_POW			4
#define OPT_NEG			5

#define OPT_SIN			6
#define OPT_COS			7
#define OPT_EXP			8
#define OPT_LOG			9
#define OPT_SQRT		10

#define NODE_BIN_OPT	0
#define NODE_UNA_OPT	1
#define NODE_VAR		2
#define NODE_CONST		3

#define VAR_ID			0
#define PAR_ID			1

void parseExpression();

namespace flowstar
{

inline void exp_taylor_only_remainder(Interval & result, const Interval & remainder, std::list<Interval>::iterator & iterRange, const unsigned int order, const Global_Computation_Setting & setting)
{
	result = 0;

	if(!iterRange->valid())
	{
		++iterRange;
		return;
	}

	Interval const_part = *iterRange;
	++iterRange;

	for(int i=order; i>0; --i)
	{
		result /= i;

		Interval intTemp;
		intTemp = (*iterRange) * remainder;		// P1 x I2
		++iterRange;
		intTemp += (*iterRange) * result;		// P2 x I1
		intTemp += remainder * result;			// I2 x I1
		++iterRange;
		intTemp += (*iterRange);				// truncation
		++iterRange;

		result = intTemp;
	}

	result *= const_part;

	result += (*iterRange);						// cutoff error
	++iterRange;

	Interval tmRange = (*iterRange) + remainder;
	++iterRange;

	Interval rem;
	exp_taylor_remainder(rem, tmRange, order+1, setting);
	result += const_part * rem;
}

inline void rec_taylor_only_remainder(Interval & result, const Interval & remainder, std::list<Interval>::iterator & iterRange, const unsigned int order, const Global_Computation_Setting & setting)
{
	result = 0;

	if(!iterRange->valid())
	{
		++iterRange;
		return;
	}

	Interval const_part = *iterRange;
	++iterRange;

	Interval tmF_c_remainder = remainder * const_part;

	for(int i=order; i>0; --i)
	{
		result *= -1;

		Interval intTemp;
		intTemp = (*iterRange) * tmF_c_remainder;	// P1 x I2
		++iterRange;
		intTemp += (*iterRange) * result;			// P2 x I1
		intTemp += tmF_c_remainder * result;		// I2 x I1
		++iterRange;
		intTemp += (*iterRange);					// truncation
		++iterRange;

		result = intTemp;
	}

	result *= const_part;

	result += (*iterRange);		// cutoff error
	++iterRange;

	Interval rem, tmF_cRange;
	tmF_cRange = (*iterRange) + tmF_c_remainder;
	++iterRange;

	rec_taylor_remainder(rem, tmF_cRange, order+1, setting);

	result += rem * const_part;
}

inline void sin_taylor_only_remainder(Interval & result, const Interval & remainder, std::list<Interval>::iterator & iterRange, const unsigned int order, const Global_Computation_Setting & setting)
{
	result = 0;

	if(!iterRange->valid())
	{
		++iterRange;
		return;
	}

	Interval const_part = *iterRange;
	++iterRange;

	Interval tmPowerTmF_remainder;

	for(int i=1; i<=order; ++i)
	{
		Interval intTemp;
		intTemp = (*iterRange) * remainder;					// P1 x I2
		++iterRange;
		intTemp += (*iterRange) * tmPowerTmF_remainder;		// P2 x I1
		intTemp += remainder * tmPowerTmF_remainder;		// I2 x I1
		++iterRange;
		intTemp += (*iterRange);							// truncation
		++iterRange;

		tmPowerTmF_remainder = intTemp;

		Interval intTemp2 = tmPowerTmF_remainder;

		intTemp2 *= (*iterRange);
		++iterRange;

		result += intTemp2;
	}

	result += (*iterRange);		// cutoff error
	++iterRange;

	Interval tmRange, rem;
	tmRange = (*iterRange) + remainder;
	++iterRange;

	sin_taylor_remainder(rem, const_part, tmRange, order+1, setting);

	result += rem;
}

inline void cos_taylor_only_remainder(Interval & result, const Interval & remainder, std::list<Interval>::iterator & iterRange, const unsigned int order, const Global_Computation_Setting & setting)
{
	result = 0;

	if(!iterRange->valid())
	{
		++iterRange;
		return;
	}

	Interval const_part = *iterRange;
	++iterRange;

	Interval tmPowerTmF_remainder;

	for(int i=1; i<=order; ++i)
	{
		Interval intTemp;
		intTemp = (*iterRange) * remainder;					// P1 x I2
		++iterRange;
		intTemp += (*iterRange) * tmPowerTmF_remainder;		// P2 x I1
		intTemp += remainder * tmPowerTmF_remainder;		// I2 x I1
		++iterRange;
		intTemp += (*iterRange);							// truncation
		++iterRange;

		tmPowerTmF_remainder = intTemp;

		Interval intTemp2 = tmPowerTmF_remainder;

		intTemp2 *= (*iterRange);
		++iterRange;

		result += intTemp2;
	}

	result += (*iterRange);		// cutoff error
	++iterRange;

	Interval tmRange, rem;
	tmRange = (*iterRange) + remainder;
	++iterRange;

	cos_taylor_remainder(rem, const_part, tmRange, order+1, setting);

	result += rem;
}

inline void log_taylor_only_remainder(Interval & result, const Interval & remainder, std::list<Interval>::iterator & iterRange, const unsigned int order)
{
	result = 0;

	if(!iterRange->valid())
	{
		++iterRange;
		return;
	}

	Interval C = *iterRange;
	++iterRange;

	Interval const_part = C;

	const_part.log_assign();

	Interval tmF_c_remainder = remainder / C;

	result = tmF_c_remainder;
	result.div_assign((double)order);

	for(int i=order; i>=2; --i)
	{
		result.inv_assign();

		Interval intTemp;
		intTemp = (*iterRange) * tmF_c_remainder;	// P1 x I2
		++iterRange;
		intTemp += (*iterRange) * result;			// P2 x I1
		intTemp += tmF_c_remainder * result;		// I2 x I1
		++iterRange;
		intTemp += (*iterRange);					// truncation
		++iterRange;

		result = intTemp;
	}

	result += (*iterRange);		// cutoff error
	++iterRange;

	Interval rem, tmF_cRange;
	tmF_cRange = (*iterRange) + tmF_c_remainder;
	++iterRange;

	log_taylor_remainder(rem, tmF_cRange, order+1);

	result += rem;
}

inline void sqrt_taylor_only_remainder(Interval & result, const Interval & remainder, std::list<Interval>::iterator & iterRange, const unsigned int order, const Global_Computation_Setting & setting)
{
	result = 0;

	if(!iterRange->valid())
	{
		++iterRange;
		return;
	}

	Interval C = *iterRange;
	++iterRange;

	Interval const_part = C;
	const_part.sqrt_assign();

	Interval tmF_2c_remainder = (remainder / C) / 2;

	result = tmF_2c_remainder;

	Interval K(1), J(1);

	for(int i=order, j=2*order-3; i>=2; --i, j-=2)
	{
		result /= -i;
		result *= j;

		Interval intTemp;
		intTemp = (*iterRange) * tmF_2c_remainder;	// P1 x I2
		++iterRange;
		intTemp += (*iterRange) * result;			// P2 x I1
		intTemp += tmF_2c_remainder * result;		// I2 x I1
		++iterRange;
		intTemp += (*iterRange);					// truncation
		++iterRange;

		result = intTemp;
	}

	result *= const_part;

	result += (*iterRange);		// cutoff error
	++iterRange;

	Interval rem, tmF_cRange;
	tmF_cRange = (*iterRange);
	++iterRange;

	tmF_cRange += tmF_2c_remainder * 2;

	flowstar::sqrt_taylor_remainder(rem, tmF_cRange, order+1, setting);

	result += rem * const_part;
}



template <class DATA_TYPE>
class TaylorModel;

template <class DATA_TYPE>
class AST_Node;

template <class DATA_TYPE>
class Expression_AST;

template <class DATA_TYPE>
class Node_Operator
{
protected:
	int type;
	std::shared_ptr<AST_Node<DATA_TYPE> > left_operand;		// only the left operand is used if the operator is unary
	std::shared_ptr<AST_Node<DATA_TYPE> > right_operand;

public:
	Node_Operator()
	{
		type = -1;
		left_operand = nullptr;
		right_operand = nullptr;
	}

	Node_Operator(const int opt_type, const std::shared_ptr<AST_Node<DATA_TYPE> > & left)
	{
		type = opt_type;
		left_operand = left;
		right_operand = nullptr;
	}

	Node_Operator(const int opt_type, const std::shared_ptr<AST_Node<DATA_TYPE> > & left,const std::shared_ptr<AST_Node<DATA_TYPE> > & right)
	{
		type = opt_type;
		left_operand = left;
		right_operand = right;
	}

	~Node_Operator()
	{
		left_operand.reset();
		right_operand.reset();
	}

	template <class DATA_TYPE2>
	friend class AST_Node;

	template <class DATA_TYPE2>
	friend class Expression_AST;
};

class Node_Variable
{
protected:
	int type;		// can only be a state variable in the current version
	int id;

public:
	Node_Variable()
	{
		type = -1;
		id = -1;
	}

	Node_Variable(const int var_type, const int var_id)
	{
		type = var_type;
		id = var_id;
	}

	~Node_Variable()
	{
	}

	template <class DATA_TYPE>
	friend class AST_Node;

	template <class DATA_TYPE>
	friend class Expression_AST;
};

// node of the abstract syntax tree
template <class DATA_TYPE>
class AST_Node
{
protected:
	int node_type;

	struct Node_Value
	{
		Node_Operator<DATA_TYPE> opt;
		Node_Variable var;
		DATA_TYPE constant;

		Node_Value()
		{
		}

		~Node_Value()
		{
		}
	} node_value;

public:
	AST_Node();
	AST_Node(const int opt_type, const std::shared_ptr<AST_Node<DATA_TYPE> > & left);
	AST_Node(const int opt_type, const std::shared_ptr<AST_Node<DATA_TYPE> > & left,const std::shared_ptr<AST_Node<DATA_TYPE> > & right);
	AST_Node(const int var_type, const unsigned int var_id);
	AST_Node(const DATA_TYPE & c);
	~AST_Node();

	void evaluate(Interval & result, const std::vector<Interval> & domain) const;

	template <class DATA_TYPE2>
	void evaluate(TaylorModel<DATA_TYPE2> & result, const std::vector<TaylorModel<DATA_TYPE2> > & tms_of_vars, const unsigned int order, const std::vector<Interval> & step_exp_table, const Interval & cutoff_threshold, const unsigned int numVars, const Global_Computation_Setting & setting) const;

	template <class DATA_TYPE2>
	void evaluate(TaylorModel<DATA_TYPE2> & result, const std::vector<TaylorModel<DATA_TYPE2> > & tms_of_vars, const unsigned int order, const std::vector<Interval> & domain, const Interval & cutoff_threshold, const Global_Computation_Setting & setting) const;

//	void evaluate(Real & result, const std::vector<Real> & values_of_vars) const;
//	void evaluate(Interval & result, const std::vector<Interval> & values_of_vars) const;

//	void evaluate(TaylorModel<DATA_TYPE> & result, const std::vector<TaylorModel> & tms_of_vars, const int order, const Taylor_Model_Computation_Setting & setting) const;
//	template <class DATA_TYPE2>
//	void evaluate(TaylorModel<DATA_TYPE2> & result, const std::vector<TaylorModel<DATA_TYPE2> > & tms_of_vars, const unsigned int order, const std::vector<Interval> & domain, const Interval & cutoff_threshold, const Global_Computation_Setting & setting) const;

//	template <class DATA_TYPE2>
//	void evaluate(TaylorModel<DATA_TYPE2> & result, const std::vector<TaylorModel<DATA_TYPE2> > & tms_of_vars, const unsigned int order, const std::vector<Interval> & step_exp_table, const Interval & cutoff_threshold, const unsigned int numVars, const Global_Computation_Setting & setting) const;

	// for internal use
	template <class DATA_TYPE2>
	void evaluate_no_remainder(TaylorModel<DATA_TYPE2> & result, const std::vector<TaylorModel<DATA_TYPE2> > & tms_of_vars, const unsigned int order, const Interval & cutoff_threshold, const unsigned int numVars) const;

	template <class DATA_TYPE2>
	void evaluate(TaylorModel<DATA_TYPE2> & result, const std::vector<TaylorModel<DATA_TYPE2> > & tms_of_vars, const unsigned int order, const std::vector<Interval> & step_exp_table, const Interval & cutoff_threshold, const unsigned int numVars, std::list<Interval> & intermediate_ranges, const Global_Computation_Setting & setting) const;

	template <class DATA_TYPE2>
	void evaluate_remainder(Interval & result, const std::vector<TaylorModel<DATA_TYPE2> > & tms_of_vars, const unsigned int order, std::list<Interval>::iterator & iter, const Global_Computation_Setting & setting) const;


//	void output(std::string & expression, const Taylor_Model_Computation_Setting & setting) const;
	void output(std::string & expression, const Variables & variables) const;

	void toReal(std::shared_ptr<AST_Node<Real> > & pNode) const;

	template <class DATA_TYPE2>
	friend class Expression_AST;
};


template <class DATA_TYPE>
AST_Node<DATA_TYPE>::AST_Node()
{
	node_type = -1;
}

template <class DATA_TYPE>
AST_Node<DATA_TYPE>::AST_Node(const int opt_type, const std::shared_ptr<AST_Node<DATA_TYPE> > & left)
{
	node_type = NODE_UNA_OPT;

	node_value.opt.type = opt_type;
	node_value.opt.left_operand = left;
	node_value.opt.right_operand = nullptr;
}

template <class DATA_TYPE>
AST_Node<DATA_TYPE>::AST_Node(const int opt_type, const std::shared_ptr<AST_Node<DATA_TYPE> > & left, const std::shared_ptr<AST_Node<DATA_TYPE> > & right)
{
	node_type = NODE_BIN_OPT;

	node_value.opt.type = opt_type;
	node_value.opt.left_operand = left;
	node_value.opt.right_operand = right;
}

template <class DATA_TYPE>
AST_Node<DATA_TYPE>::AST_Node(const int var_type, const unsigned int var_id)
{
	node_type = NODE_VAR;

	node_value.var.type = var_type;
	node_value.var.id = var_id;
}

template <class DATA_TYPE>
AST_Node<DATA_TYPE>::AST_Node(const DATA_TYPE & c)
{
	node_type = NODE_CONST;
	node_value.constant = c;
}

template <class DATA_TYPE>
AST_Node<DATA_TYPE>::~AST_Node()
{
}

template <class DATA_TYPE>
void AST_Node<DATA_TYPE>::evaluate(Interval & result, const std::vector<Interval> & domain) const
{
	switch(node_type)
	{
	case NODE_UNA_OPT:
	{
		node_value.opt.left_operand->evaluate(result, domain);

		switch(node_value.opt.type)
		{
		case OPT_NEG:
			result *= -1;
			break;

		case OPT_SIN:
			result.sin_assign();
			break;

		case OPT_COS:
			result.cos_assign();
			break;

		case OPT_EXP:
			result.exp_assign();
			break;

		case OPT_LOG:
			result.log_assign();
			break;

		case OPT_SQRT:
			result.sqrt_assign();
			break;
		}

		break;
	}
	case NODE_BIN_OPT:
	{
		Interval I1, I2;

		switch(node_value.opt.type)
		{
		case OPT_PLUS:
			node_value.opt.left_operand->evaluate(I1, domain);
			node_value.opt.right_operand->evaluate(I2, domain);
			result = I1 + I2;
			break;
		case OPT_MINU:
			node_value.opt.left_operand->evaluate(I1, domain);
			node_value.opt.right_operand->evaluate(I2, domain);
			result = I1 - I2;
			break;
		case OPT_MULT:
			node_value.opt.left_operand->evaluate(I1, domain);
			node_value.opt.right_operand->evaluate(I2, domain);
			result = I1 * I2;
			break;
		case OPT_DIV:
			node_value.opt.left_operand->evaluate(I1, domain);
			node_value.opt.right_operand->evaluate(I2, domain);
			result = I1 / I2;
			break;
		case OPT_POW:
			node_value.opt.left_operand->evaluate(result, domain);
			result.pow_assign((int)node_value.opt.right_operand->node_value.constant.toDouble());
			break;
		}

		break;
	}

	case NODE_VAR:
		if(node_value.var.type == VAR_ID)
		{
			result = domain[node_value.var.id];
		}
		else
		{
			// for the other variable types
		}

		break;

	case NODE_CONST:
		result = node_value.constant;
		break;
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void AST_Node<DATA_TYPE>::evaluate(TaylorModel<DATA_TYPE2> & result, const std::vector<TaylorModel<DATA_TYPE2> > & tms_of_vars, const unsigned int order, const std::vector<Interval> & step_exp_table, const Interval & cutoff_threshold, const unsigned int numVars, const Global_Computation_Setting & setting) const
{
	switch(node_type)
	{
	case NODE_UNA_OPT:
	{
		TaylorModel<DATA_TYPE2> tmTemp;
		node_value.opt.left_operand->evaluate(tmTemp, tms_of_vars, order, step_exp_table, cutoff_threshold, numVars, setting);

		switch(node_value.opt.type)
		{
		case OPT_NEG:
			tmTemp *= -1;
			result = tmTemp;
			break;

		case OPT_SIN:
			tmTemp.sin_taylor(result, step_exp_table, numVars, order, cutoff_threshold, setting);
			break;

		case OPT_COS:
			tmTemp.cos_taylor(result, step_exp_table, numVars, order, cutoff_threshold, setting);
			break;

		case OPT_EXP:
			tmTemp.exp_taylor(result, step_exp_table, numVars, order, cutoff_threshold, setting);
			break;

		case OPT_LOG:
			tmTemp.log_taylor(result, step_exp_table, numVars, order, cutoff_threshold);
			break;

		case OPT_SQRT:
			tmTemp.sqrt_taylor(result, step_exp_table, numVars, order, cutoff_threshold, setting);
			break;
		}
		break;
	}
	case NODE_BIN_OPT:
	{
		TaylorModel<DATA_TYPE2> tm1, tm2;

		switch(node_value.opt.type)
		{
		case OPT_PLUS:
			node_value.opt.left_operand->evaluate(tm1, tms_of_vars, order, step_exp_table, cutoff_threshold, numVars, setting);
			node_value.opt.right_operand->evaluate(tm2, tms_of_vars, order, step_exp_table, cutoff_threshold, numVars, setting);
			result = tm1 + tm2;
			break;

		case OPT_MINU:
			node_value.opt.left_operand->evaluate(tm1, tms_of_vars, order, step_exp_table, cutoff_threshold, numVars, setting);
			node_value.opt.right_operand->evaluate(tm2, tms_of_vars, order, step_exp_table, cutoff_threshold, numVars, setting);
			result = tm1 - tm2;
			break;

		case OPT_MULT:
		{
			node_value.opt.left_operand->evaluate(tm1, tms_of_vars, order, step_exp_table, cutoff_threshold, numVars, setting);
			node_value.opt.right_operand->evaluate(tm2, tms_of_vars, order, step_exp_table, cutoff_threshold, numVars, setting);

			tm1.mul_ctrunc_normal(result, tm2, step_exp_table, order, cutoff_threshold);

			break;
		}

		case OPT_DIV:
		{
			node_value.opt.left_operand->evaluate(result, tms_of_vars, order, step_exp_table, cutoff_threshold, numVars, setting);
			node_value.opt.right_operand->evaluate(tm2, tms_of_vars, order, step_exp_table, cutoff_threshold, numVars, setting);

			TaylorModel<DATA_TYPE2> tmTemp;
			tm2.rec_taylor(tmTemp, step_exp_table, numVars, order, cutoff_threshold, setting);

			result.mul_ctrunc_normal_assign(tmTemp, step_exp_table, order, cutoff_threshold);

			break;
		}

		case OPT_POW:
		{
			node_value.opt.left_operand->evaluate(result, tms_of_vars, order, step_exp_table, cutoff_threshold, numVars, setting);

			int degree = (int)node_value.opt.right_operand->node_value.constant.toDouble();

			if(degree == 0)
			{
				TaylorModel<DATA_TYPE2> tm(1, numVars);
				result = tm;
			}
			else if(degree > 1)
			{
				TaylorModel<DATA_TYPE2> temp = result;

				for(int i = degree - 1; i > 0;)
				{
					Interval intPoly2;
					temp.polyRangeNormal(intPoly2, step_exp_table);

					if(i & 1)
					{
						// result *= temp
						result.mul_insert_ctrunc_normal_assign(temp, intPoly2, step_exp_table, order, cutoff_threshold);
					}

					i >>= 1;

					if(i > 0)
					{
						// temp *= temp
						temp.mul_insert_ctrunc_normal_assign(temp, intPoly2, step_exp_table, order, cutoff_threshold);
					}
				}
			}

			break;
		}
		}

		break;
	}

	case NODE_VAR:
		if(node_value.var.type == VAR_ID)
		{
			result = tms_of_vars[node_value.var.id];
			result.ctrunc_normal(step_exp_table, order);
		}
		else
		{
			// for the other variable types
		}

		break;

	case NODE_CONST:
	{
		TaylorModel<DATA_TYPE2> temp(node_value.constant, numVars);
		result = temp;
		break;
	}
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void AST_Node<DATA_TYPE>::evaluate(TaylorModel<DATA_TYPE2> & result, const std::vector<TaylorModel<DATA_TYPE2> > & tms_of_vars, const unsigned int order, const std::vector<Interval> & domain, const Interval & cutoff_threshold, const Global_Computation_Setting & setting) const
{
	switch(node_type)
	{
	case NODE_UNA_OPT:
	{
		TaylorModel<DATA_TYPE2> tmTemp;
		node_value.opt.left_operand->evaluate(tmTemp, tms_of_vars, order, domain, cutoff_threshold, setting);

		switch(node_value.opt.type)
		{
		case OPT_NEG:
			tmTemp *= -1;
			result = tmTemp;
			break;

		case OPT_SIN:
			tmTemp.sin_taylor(result, domain, order, cutoff_threshold, setting);
			break;

		case OPT_COS:
			tmTemp.cos_taylor(result, domain, order, cutoff_threshold, setting);
			break;

		case OPT_EXP:
			tmTemp.exp_taylor(result, domain, order, cutoff_threshold, setting);
			break;

		case OPT_LOG:
			tmTemp.log_taylor(result, domain, order, cutoff_threshold);
			break;

		case OPT_SQRT:
			tmTemp.sqrt_taylor(result, domain, order, cutoff_threshold, setting);
			break;
		}
		break;
	}
	case NODE_BIN_OPT:
	{
		TaylorModel<DATA_TYPE2> tm1, tm2;

		switch(node_value.opt.type)
		{
		case OPT_PLUS:
			node_value.opt.left_operand->evaluate(tm1, tms_of_vars, order, domain, cutoff_threshold, setting);
			node_value.opt.right_operand->evaluate(tm2, tms_of_vars, order, domain, cutoff_threshold, setting);
			result = tm1 + tm2;
			break;

		case OPT_MINU:
			node_value.opt.left_operand->evaluate(tm1, tms_of_vars, order, domain, cutoff_threshold, setting);
			node_value.opt.right_operand->evaluate(tm2, tms_of_vars, order, domain, cutoff_threshold, setting);
			result = tm1 - tm2;
			break;

		case OPT_MULT:
		{
			node_value.opt.left_operand->evaluate(tm1, tms_of_vars, order, domain, cutoff_threshold, setting);
			node_value.opt.right_operand->evaluate(tm2, tms_of_vars, order, domain, cutoff_threshold, setting);

			tm1.mul_ctrunc(result, tm2, domain, order, cutoff_threshold);

			break;
		}

		case OPT_DIV:
		{
			node_value.opt.left_operand->evaluate(result, tms_of_vars, order, domain, cutoff_threshold, setting);
			node_value.opt.right_operand->evaluate(tm2, tms_of_vars, order, domain, cutoff_threshold, setting);

			TaylorModel<DATA_TYPE2> tmTemp;
			tm2.rec_taylor(tmTemp, domain, order, cutoff_threshold, setting);

			result.mul_ctrunc_assign(tmTemp, domain, order, cutoff_threshold);

			break;
		}

		case OPT_POW:
		{
			node_value.opt.left_operand->evaluate(result, tms_of_vars, order, domain, cutoff_threshold, setting);

			int degree = (int)node_value.opt.right_operand->node_value.constant.toDouble();
			unsigned int numVars = domain.size();

			if(degree == 0)
			{
				TaylorModel<DATA_TYPE2> tm(1, numVars);
				result = tm;
			}
			else if(degree > 1)
			{
				TaylorModel<DATA_TYPE2> temp = result;

				for(int i = degree - 1; i > 0;)
				{
					Interval intPoly2;
					temp.polyRange(intPoly2, domain);

					if(i & 1)
					{
						// result *= temp
						result.mul_insert_ctrunc_assign(temp, intPoly2, domain, order, cutoff_threshold);
					}

					i >>= 1;

					if(i > 0)
					{
						// temp *= temp
						temp.mul_insert_ctrunc_assign(temp, intPoly2, domain, order, cutoff_threshold);
					}
				}
			}

			break;
		}
		}

		break;
	}

	case NODE_VAR:
		if(node_value.var.type == VAR_ID)
		{
			result = tms_of_vars[node_value.var.id];
			result.ctrunc(domain, order);
		}
		else
		{
			// for the other variable types
		}

		break;

	case NODE_CONST:
	{
		TaylorModel<DATA_TYPE2> temp(node_value.constant, domain.size());
		result = temp;
		break;
	}
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void AST_Node<DATA_TYPE>::evaluate_no_remainder(TaylorModel<DATA_TYPE2> & result, const std::vector<TaylorModel<DATA_TYPE2> > & tms_of_vars, const unsigned int order, const Interval & cutoff_threshold, const unsigned int numVars) const
{
	result.remainder = 0;

	switch(node_type)
	{
	case NODE_UNA_OPT:
	{
		TaylorModel<DATA_TYPE2> tmTemp;
		node_value.opt.left_operand->evaluate_no_remainder(tmTemp, tms_of_vars, order, cutoff_threshold, numVars);

		switch(node_value.opt.type)
		{
		case OPT_NEG:
			tmTemp.expansion *= -1;
			result.expansion = tmTemp.expansion;
			break;

		case OPT_SIN:
			tmTemp.expansion.sin_taylor(result.expansion, numVars, order, cutoff_threshold);
			break;

		case OPT_COS:
			tmTemp.expansion.cos_taylor(result.expansion, numVars, order, cutoff_threshold);
			break;

		case OPT_EXP:
			tmTemp.expansion.exp_taylor(result.expansion, numVars, order, cutoff_threshold);
			break;

		case OPT_LOG:
			tmTemp.expansion.log_taylor(result.expansion, numVars, order, cutoff_threshold);
			break;

		case OPT_SQRT:
			tmTemp.expansion.sqrt_taylor(result.expansion, numVars, order, cutoff_threshold);
			break;
		}
		break;
	}
	case NODE_BIN_OPT:
	{
		TaylorModel<DATA_TYPE2> tm1, tm2;

		switch(node_value.opt.type)
		{
		case OPT_PLUS:
			node_value.opt.left_operand->evaluate_no_remainder(tm1, tms_of_vars, order, cutoff_threshold, numVars);
			node_value.opt.right_operand->evaluate_no_remainder(tm2, tms_of_vars, order, cutoff_threshold, numVars);
			result = tm1 + tm2;
			break;
		case OPT_MINU:
			node_value.opt.left_operand->evaluate_no_remainder(tm1, tms_of_vars, order, cutoff_threshold, numVars);
			node_value.opt.right_operand->evaluate_no_remainder(tm2, tms_of_vars, order, cutoff_threshold, numVars);
			result = tm1 - tm2;
			break;
		case OPT_MULT:
			node_value.opt.left_operand->evaluate_no_remainder(tm1, tms_of_vars, order, cutoff_threshold, numVars);
			node_value.opt.right_operand->evaluate_no_remainder(tm2, tms_of_vars, order, cutoff_threshold, numVars);
			result.expansion = tm1.expansion * tm2.expansion;
			result.expansion.nctrunc(order);
			result.expansion.cutoff(cutoff_threshold);
			break;

		case OPT_DIV:
		{
			node_value.opt.left_operand->evaluate_no_remainder(tm1, tms_of_vars, order, cutoff_threshold, numVars);
			node_value.opt.right_operand->evaluate_no_remainder(tm2, tms_of_vars, order, cutoff_threshold, numVars);

			Polynomial<DATA_TYPE2> polyTemp;
			tm2.expansion.rec_taylor(polyTemp, numVars, order, cutoff_threshold);

			result.expansion = tm1.expansion * polyTemp;
			result.expansion.nctrunc(order);
			result.expansion.cutoff(cutoff_threshold);

			break;
		}

		case OPT_POW:
		{
			node_value.opt.left_operand->evaluate_no_remainder(result, tms_of_vars, order, cutoff_threshold, numVars);

			int degree = (int)node_value.opt.right_operand->node_value.constant.toDouble();

			if(degree == 0)
			{
				TaylorModel<DATA_TYPE2> tm(1, numVars);
				result = tm;
			}
			else if(degree > 1)
			{
				TaylorModel<DATA_TYPE2> temp = result;

				for(int i = degree - 1; i > 0;)
				{
					if(i & 1)
					{
						// result *= temp
						result.expansion *= temp.expansion;
						result.expansion.nctrunc(order);
						result.expansion.cutoff(cutoff_threshold);
					}

					i >>= 1;

					if(i > 0)
					{
						// temp *= temp
						temp.expansion *= temp.expansion;
						temp.expansion.nctrunc(order);
						temp.expansion.cutoff(cutoff_threshold);
					}
				}
			}

			break;
		}
		}

		break;
	}

	case NODE_VAR:
		if(node_value.var.type == VAR_ID)
		{
			result.expansion = tms_of_vars[node_value.var.id].expansion;
			result.expansion.nctrunc(order);
		}
		else
		{
			// for the other variable types
		}

		break;

	case NODE_CONST:
	{
		TaylorModel<DATA_TYPE2> temp(node_value.constant, numVars);
		result = temp;
		break;
	}
	}
}

template <>
template <>
inline void AST_Node<Interval>::evaluate_no_remainder<Real>(TaylorModel<Real> & result, const std::vector<TaylorModel<Real> > & tms_of_vars, const unsigned int order, const Interval & cutoff_threshold, const unsigned int numVars) const
{
	result.remainder = 0;

	switch(node_type)
	{
	case NODE_UNA_OPT:
	{
		TaylorModel<Real> tmTemp;
		node_value.opt.left_operand->evaluate_no_remainder(tmTemp, tms_of_vars, order, cutoff_threshold, numVars);

		switch(node_value.opt.type)
		{
		case OPT_NEG:
			tmTemp.expansion *= -1;
			result.expansion = tmTemp.expansion;
			break;

		case OPT_SIN:
			tmTemp.expansion.sin_taylor(result.expansion, numVars, order, cutoff_threshold);
			break;

		case OPT_COS:
			tmTemp.expansion.cos_taylor(result.expansion, numVars, order, cutoff_threshold);
			break;

		case OPT_EXP:
			tmTemp.expansion.exp_taylor(result.expansion, numVars, order, cutoff_threshold);
			break;

		case OPT_LOG:
			tmTemp.expansion.log_taylor(result.expansion, numVars, order, cutoff_threshold);
			break;

		case OPT_SQRT:
			tmTemp.expansion.sqrt_taylor(result.expansion, numVars, order, cutoff_threshold);
			break;
		}
		break;
	}
	case NODE_BIN_OPT:
	{
		TaylorModel<Real> tm1, tm2;

		switch(node_value.opt.type)
		{
		case OPT_PLUS:
			node_value.opt.left_operand->evaluate_no_remainder(tm1, tms_of_vars, order, cutoff_threshold, numVars);
			node_value.opt.right_operand->evaluate_no_remainder(tm2, tms_of_vars, order, cutoff_threshold, numVars);
			result = tm1 + tm2;
			break;
		case OPT_MINU:
			node_value.opt.left_operand->evaluate_no_remainder(tm1, tms_of_vars, order, cutoff_threshold, numVars);
			node_value.opt.right_operand->evaluate_no_remainder(tm2, tms_of_vars, order, cutoff_threshold, numVars);
			result = tm1 - tm2;
			break;
		case OPT_MULT:
			node_value.opt.left_operand->evaluate_no_remainder(tm1, tms_of_vars, order, cutoff_threshold, numVars);
			node_value.opt.right_operand->evaluate_no_remainder(tm2, tms_of_vars, order, cutoff_threshold, numVars);
			result.expansion = tm1.expansion * tm2.expansion;
			result.expansion.nctrunc(order);
			result.expansion.cutoff(cutoff_threshold);
			break;

		case OPT_DIV:
		{
			node_value.opt.left_operand->evaluate_no_remainder(tm1, tms_of_vars, order, cutoff_threshold, numVars);
			node_value.opt.right_operand->evaluate_no_remainder(tm2, tms_of_vars, order, cutoff_threshold, numVars);

			Polynomial<Real> polyTemp;
			tm2.expansion.rec_taylor(polyTemp, numVars, order, cutoff_threshold);

			result.expansion = tm1.expansion * polyTemp;
			result.expansion.nctrunc(order);
			result.expansion.cutoff(cutoff_threshold);

			break;
		}

		case OPT_POW:
		{
			node_value.opt.left_operand->evaluate_no_remainder(result, tms_of_vars, order, cutoff_threshold, numVars);

			int degree = (int)node_value.opt.right_operand->node_value.constant.toDouble();

			if(degree == 0)
			{
				TaylorModel<Real> tm(1, numVars);
				result = tm;
			}
			else if(degree > 1)
			{
				TaylorModel<Real> temp = result;

				for(int i = degree - 1; i > 0;)
				{
					if(i & 1)
					{
						// result *= temp
						result.expansion *= temp.expansion;
						result.expansion.nctrunc(order);
						result.expansion.cutoff(cutoff_threshold);
					}

					i >>= 1;

					if(i > 0)
					{
						// temp *= temp
						temp.expansion *= temp.expansion;
						temp.expansion.nctrunc(order);
						temp.expansion.cutoff(cutoff_threshold);
					}
				}
			}

			break;
		}
		}

		break;
	}

	case NODE_VAR:
		if(node_value.var.type == VAR_ID)
		{
			result.expansion = tms_of_vars[node_value.var.id].expansion;
			result.expansion.nctrunc(order);
		}
		else
		{
			// for the other variable types
		}

		break;

	case NODE_CONST:
	{
		TaylorModel<Real> temp(node_value.constant.toReal(), numVars);
		result = temp;
		break;
	}
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void AST_Node<DATA_TYPE>::evaluate(TaylorModel<DATA_TYPE2> & result, const std::vector<TaylorModel<DATA_TYPE2> > & tms_of_vars, const unsigned int order, const std::vector<Interval> & step_exp_table, const Interval & cutoff_threshold, const unsigned int numVars, std::list<Interval> & intermediate_ranges, const Global_Computation_Setting & setting) const
{
	switch(node_type)
	{
	case NODE_UNA_OPT:
	{
		TaylorModel<DATA_TYPE2> tmTemp;
		node_value.opt.left_operand->evaluate(tmTemp, tms_of_vars, order, step_exp_table, cutoff_threshold, numVars, intermediate_ranges, setting);

		switch(node_value.opt.type)
		{
		case OPT_NEG:
			tmTemp *= -1;
			result = tmTemp;
			break;

		case OPT_SIN:
			tmTemp.sin_taylor(result, intermediate_ranges, step_exp_table, numVars, order, cutoff_threshold, setting);
			break;

		case OPT_COS:
			tmTemp.cos_taylor(result, intermediate_ranges, step_exp_table, numVars, order, cutoff_threshold, setting);
			break;

		case OPT_EXP:
			tmTemp.exp_taylor(result, intermediate_ranges, step_exp_table, numVars, order, cutoff_threshold, setting);
			break;

		case OPT_LOG:
			tmTemp.log_taylor(result, intermediate_ranges, step_exp_table, numVars, order, cutoff_threshold);
			break;

		case OPT_SQRT:
			tmTemp.sqrt_taylor(result, intermediate_ranges, step_exp_table, numVars, order, cutoff_threshold, setting);
			break;
		}
		break;
	}
	case NODE_BIN_OPT:
	{
		TaylorModel<DATA_TYPE2> tm1, tm2;

		switch(node_value.opt.type)
		{
		case OPT_PLUS:
			node_value.opt.left_operand->evaluate(tm1, tms_of_vars, order, step_exp_table, cutoff_threshold, numVars, intermediate_ranges, setting);
			node_value.opt.right_operand->evaluate(tm2, tms_of_vars, order, step_exp_table, cutoff_threshold, numVars, intermediate_ranges, setting);
			result = tm1 + tm2;
			break;

		case OPT_MINU:
			node_value.opt.left_operand->evaluate(tm1, tms_of_vars, order, step_exp_table, cutoff_threshold, numVars, intermediate_ranges, setting);
			node_value.opt.right_operand->evaluate(tm2, tms_of_vars, order, step_exp_table, cutoff_threshold, numVars, intermediate_ranges, setting);
			result = tm1 - tm2;
			break;

		case OPT_MULT:
		{
			node_value.opt.left_operand->evaluate(tm1, tms_of_vars, order, step_exp_table, cutoff_threshold, numVars, intermediate_ranges, setting);
			node_value.opt.right_operand->evaluate(tm2, tms_of_vars, order, step_exp_table, cutoff_threshold, numVars, intermediate_ranges, setting);

			Interval intPoly1, intPoly2, intTrunc;

			tm2.polyRangeNormal(intPoly2, step_exp_table);
			tm1.mul_insert_ctrunc_normal(result, intPoly1, intTrunc, tm2, intPoly2, step_exp_table, order, cutoff_threshold);

			intermediate_ranges.push_back(intPoly1);
			intermediate_ranges.push_back(intPoly2);
			intermediate_ranges.push_back(intTrunc);

			break;
		}

		case OPT_DIV:
		{
			node_value.opt.left_operand->evaluate(result, tms_of_vars, order, step_exp_table, cutoff_threshold, numVars, intermediate_ranges, setting);
			node_value.opt.right_operand->evaluate(tm2, tms_of_vars, order, step_exp_table, cutoff_threshold, numVars, intermediate_ranges, setting);

			TaylorModel<DATA_TYPE2> tmTemp;
			tm2.rec_taylor(tmTemp, intermediate_ranges, step_exp_table, numVars, order, cutoff_threshold, setting);

			Interval intPoly1, intPoly2, intTrunc;

			tmTemp.polyRangeNormal(intPoly2, step_exp_table);
			result.mul_insert_ctrunc_normal_assign(intPoly1, intTrunc, tmTemp, intPoly2, step_exp_table, order, cutoff_threshold);

			intermediate_ranges.push_back(intPoly1);
			intermediate_ranges.push_back(intPoly2);
			intermediate_ranges.push_back(intTrunc);

			break;
		}

		case OPT_POW:
		{
			node_value.opt.left_operand->evaluate(result, tms_of_vars, order, step_exp_table, cutoff_threshold, numVars, intermediate_ranges, setting);

			int degree = (int)node_value.opt.right_operand->node_value.constant.toDouble();

			if(degree == 0)
			{
				TaylorModel<DATA_TYPE2> tm(1, numVars);
				result = tm;
			}
			else if(degree > 1)
			{
				TaylorModel<DATA_TYPE2> temp = result;
				Interval intPoly1, intPoly2, intTrunc;

				for(int i = degree - 1; i > 0;)
				{
					temp.polyRangeNormal(intPoly2, step_exp_table);

					if(i & 1)
					{
						// result *= temp
						result.mul_insert_ctrunc_normal_assign(intPoly1, intTrunc, temp, intPoly2, step_exp_table, order, cutoff_threshold);

						intermediate_ranges.push_back(intPoly1);
						intermediate_ranges.push_back(intPoly2);
						intermediate_ranges.push_back(intTrunc);
					}

					i >>= 1;

					if(i > 0)
					{
						// temp *= temp
						temp.mul_insert_ctrunc_normal_assign(intPoly1, intTrunc, temp, intPoly2, step_exp_table, order, cutoff_threshold);

						intermediate_ranges.push_back(intPoly1);
						intermediate_ranges.push_back(intPoly2);
						intermediate_ranges.push_back(intTrunc);
					}
				}
			}

			break;
		}
		}

		break;
	}

	case NODE_VAR:
		if(node_value.var.type == VAR_ID)
		{
			result = tms_of_vars[node_value.var.id];
			result.ctrunc_normal(step_exp_table, order);
		}
		else
		{
			// for the other variable types
		}

		break;

	case NODE_CONST:
	{
		TaylorModel<DATA_TYPE2> temp(node_value.constant, numVars);
		result = temp;
		break;
	}
	}
}

template <>
template <>
inline void AST_Node<Interval>::evaluate<Real>(TaylorModel<Real> & result, const std::vector<TaylorModel<Real> > & tms_of_vars, const unsigned int order, const std::vector<Interval> & step_exp_table, const Interval & cutoff_threshold, const unsigned int numVars, std::list<Interval> & intermediate_ranges, const Global_Computation_Setting & setting) const
{
	switch(node_type)
	{
	case NODE_UNA_OPT:
	{
		TaylorModel<Real> tmTemp;
		node_value.opt.left_operand->evaluate(tmTemp, tms_of_vars, order, step_exp_table, cutoff_threshold, numVars, intermediate_ranges, setting);

		switch(node_value.opt.type)
		{
		case OPT_NEG:
			tmTemp *= -1;
			result = tmTemp;
			break;

		case OPT_SIN:
			tmTemp.sin_taylor(result, intermediate_ranges, step_exp_table, numVars, order, cutoff_threshold, setting);
			break;

		case OPT_COS:
			tmTemp.cos_taylor(result, intermediate_ranges, step_exp_table, numVars, order, cutoff_threshold, setting);
			break;

		case OPT_EXP:
			tmTemp.exp_taylor(result, intermediate_ranges, step_exp_table, numVars, order, cutoff_threshold, setting);
			break;

		case OPT_LOG:
			tmTemp.log_taylor(result, intermediate_ranges, step_exp_table, numVars, order, cutoff_threshold);
			break;

		case OPT_SQRT:
			tmTemp.sqrt_taylor(result, intermediate_ranges, step_exp_table, numVars, order, cutoff_threshold, setting);
			break;
		}
		break;
	}
	case NODE_BIN_OPT:
	{
		TaylorModel<Real> tm1, tm2;

		switch(node_value.opt.type)
		{
		case OPT_PLUS:
			node_value.opt.left_operand->evaluate(tm1, tms_of_vars, order, step_exp_table, cutoff_threshold, numVars, intermediate_ranges, setting);
			node_value.opt.right_operand->evaluate(tm2, tms_of_vars, order, step_exp_table, cutoff_threshold, numVars, intermediate_ranges, setting);
			result = tm1 + tm2;
			break;

		case OPT_MINU:
			node_value.opt.left_operand->evaluate(tm1, tms_of_vars, order, step_exp_table, cutoff_threshold, numVars, intermediate_ranges, setting);
			node_value.opt.right_operand->evaluate(tm2, tms_of_vars, order, step_exp_table, cutoff_threshold, numVars, intermediate_ranges, setting);
			result = tm1 - tm2;
			break;

		case OPT_MULT:
		{
			node_value.opt.left_operand->evaluate(tm1, tms_of_vars, order, step_exp_table, cutoff_threshold, numVars, intermediate_ranges, setting);
			node_value.opt.right_operand->evaluate(tm2, tms_of_vars, order, step_exp_table, cutoff_threshold, numVars, intermediate_ranges, setting);

			Interval intPoly1, intPoly2, intTrunc;

			tm2.polyRangeNormal(intPoly2, step_exp_table);
			tm1.mul_insert_ctrunc_normal(result, intPoly1, intTrunc, tm2, intPoly2, step_exp_table, order, cutoff_threshold);

			intermediate_ranges.push_back(intPoly1);
			intermediate_ranges.push_back(intPoly2);
			intermediate_ranges.push_back(intTrunc);

			break;
		}

		case OPT_DIV:
		{
			node_value.opt.left_operand->evaluate(result, tms_of_vars, order, step_exp_table, cutoff_threshold, numVars, intermediate_ranges, setting);
			node_value.opt.right_operand->evaluate(tm2, tms_of_vars, order, step_exp_table, cutoff_threshold, numVars, intermediate_ranges, setting);

			TaylorModel<Real> tmTemp;
			tm2.rec_taylor(tmTemp, intermediate_ranges, step_exp_table, numVars, order, cutoff_threshold, setting);

			Interval intPoly1, intPoly2, intTrunc;

			tmTemp.polyRangeNormal(intPoly2, step_exp_table);
			result.mul_insert_ctrunc_normal_assign(intPoly1, intTrunc, tmTemp, intPoly2, step_exp_table, order, cutoff_threshold);

			intermediate_ranges.push_back(intPoly1);
			intermediate_ranges.push_back(intPoly2);
			intermediate_ranges.push_back(intTrunc);

			break;
		}

		case OPT_POW:
		{
			node_value.opt.left_operand->evaluate(result, tms_of_vars, order, step_exp_table, cutoff_threshold, numVars, intermediate_ranges, setting);

			int degree = (int)node_value.opt.right_operand->node_value.constant.toDouble();

			if(degree == 0)
			{
				TaylorModel<Real> tm(1, numVars);
				result = tm;
			}
			else if(degree > 1)
			{
				TaylorModel<Real> temp = result;
				Interval intPoly1, intPoly2, intTrunc;

				for(int i = degree - 1; i > 0;)
				{
					temp.polyRangeNormal(intPoly2, step_exp_table);

					if(i & 1)
					{
						// result *= temp
						result.mul_insert_ctrunc_normal_assign(intPoly1, intTrunc, temp, intPoly2, step_exp_table, order, cutoff_threshold);

						intermediate_ranges.push_back(intPoly1);
						intermediate_ranges.push_back(intPoly2);
						intermediate_ranges.push_back(intTrunc);
					}

					i >>= 1;

					if(i > 0)
					{
						// temp *= temp
						temp.mul_insert_ctrunc_normal_assign(intPoly1, intTrunc, temp, intPoly2, step_exp_table, order, cutoff_threshold);

						intermediate_ranges.push_back(intPoly1);
						intermediate_ranges.push_back(intPoly2);
						intermediate_ranges.push_back(intTrunc);
					}
				}
			}

			break;
		}
		}

		break;
	}

	case NODE_VAR:
		if(node_value.var.type == VAR_ID)
		{
			result = tms_of_vars[node_value.var.id];
			result.ctrunc_normal(step_exp_table, order);
		}
		else
		{
			// for the other variable types
		}

		break;

	case NODE_CONST:
	{
		Real c;
		Interval I = node_value.constant;

		I.remove_midpoint(c);
		TaylorModel<Real> temp(c, numVars);
		temp.remainder = I;

		intermediate_ranges.push_back(I);

		result = temp;
		break;
	}
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void AST_Node<DATA_TYPE>::evaluate_remainder(Interval & result, const std::vector<TaylorModel<DATA_TYPE2> > & tms_of_vars, const unsigned int order, std::list<Interval>::iterator & iter, const Global_Computation_Setting & setting) const
{
	switch(node_type)
	{
	case NODE_UNA_OPT:
	{
		Interval intTemp;
		node_value.opt.left_operand->evaluate_remainder(intTemp, tms_of_vars, order, iter, setting);

		switch(node_value.opt.type)
		{
		case OPT_NEG:
			intTemp.inv(result);
			break;

		case OPT_SIN:
			sin_taylor_only_remainder(result, intTemp, iter, order, setting);
			break;

		case OPT_COS:
			cos_taylor_only_remainder(result, intTemp, iter, order, setting);
			break;

		case OPT_EXP:
			exp_taylor_only_remainder(result, intTemp, iter, order, setting);
			break;

		case OPT_LOG:
			log_taylor_only_remainder(result, intTemp, iter, order);
			break;

		case OPT_SQRT:
			sqrt_taylor_only_remainder(result, intTemp, iter, order, setting);
			break;
		}
		break;
	}
	case NODE_BIN_OPT:
	{
		Interval remainder1, remainder2;

		switch(node_value.opt.type)
		{
		case OPT_PLUS:
			node_value.opt.left_operand->evaluate_remainder(remainder1, tms_of_vars, order, iter, setting);
			node_value.opt.right_operand->evaluate_remainder(remainder2, tms_of_vars, order, iter, setting);
			result = remainder1 + remainder2;
			break;
		case OPT_MINU:
			node_value.opt.left_operand->evaluate_remainder(remainder1, tms_of_vars, order, iter, setting);
			node_value.opt.right_operand->evaluate_remainder(remainder2, tms_of_vars, order, iter, setting);
			result = remainder1 - remainder2;
			break;
		case OPT_MULT:
			node_value.opt.left_operand->evaluate_remainder(remainder1, tms_of_vars, order, iter, setting);
			node_value.opt.right_operand->evaluate_remainder(remainder2, tms_of_vars, order, iter, setting);

			result = (*iter) * remainder2;
			++iter;
			result += (*iter) * remainder1;
			result += remainder1 * remainder2;
			++iter;
			result += (*iter);
			++iter;

			break;

		case OPT_DIV:
		{
			node_value.opt.left_operand->evaluate_remainder(remainder1, tms_of_vars, order, iter, setting);
			node_value.opt.right_operand->evaluate_remainder(remainder2, tms_of_vars, order, iter, setting);

			Interval intTemp;
			rec_taylor_only_remainder(intTemp, remainder2, iter, order, setting);

			result = (*iter) * intTemp;
			++iter;
			result += (*iter) * remainder1;
			result += remainder1 * intTemp;
			++iter;
			result += (*iter);
			++iter;
			break;
		}

		case OPT_POW:
		{
			node_value.opt.left_operand->evaluate_remainder(result, tms_of_vars, order, iter, setting);

			int degree = (int)node_value.opt.right_operand->node_value.constant.toDouble();

			if(degree == 0)
			{
				result = 0;
			}
			else if(degree > 1)
			{
				Interval temp = result;

				for(int i = degree - 1; i > 0;)
				{
					if(i & 1)
					{
						Interval temp2;
						temp2 = (*iter) * temp;
						++iter;
						temp2 += (*iter) * result;
						temp2 += temp * result;
						++iter;
						temp2 += (*iter);
						++iter;

						result = temp2;
					}

					i >>= 1;

					if(i > 0)
					{
						Interval temp2;
						temp2 = (*iter) * temp;
						++iter;
						temp2 += (*iter) * temp;
						temp2 += temp * temp;
						++iter;
						temp2 += (*iter);
						++iter;

						temp = temp2;
					}
				}
			}

			break;
		}
		}

		break;
	}

	case NODE_VAR:
		if(node_value.var.type == VAR_ID)
		{
			result = tms_of_vars[node_value.var.id].remainder;
		}
		else
		{
			// for the other variable types
		}

		break;

	case NODE_CONST:
	{
		result = 0;
		break;
	}
	}
}

template <>
template <>
inline void AST_Node<Interval>::evaluate_remainder<Real>(Interval & result, const std::vector<TaylorModel<Real> > & tms_of_vars, const unsigned int order, std::list<Interval>::iterator & iter, const Global_Computation_Setting & setting) const
{
	switch(node_type)
	{
	case NODE_UNA_OPT:
	{
		Interval intTemp;
		node_value.opt.left_operand->evaluate_remainder(intTemp, tms_of_vars, order, iter, setting);

		switch(node_value.opt.type)
		{
		case OPT_NEG:
			intTemp.inv(result);
			break;

		case OPT_SIN:
			sin_taylor_only_remainder(result, intTemp, iter, order, setting);
			break;

		case OPT_COS:
			cos_taylor_only_remainder(result, intTemp, iter, order, setting);
			break;

		case OPT_EXP:
			exp_taylor_only_remainder(result, intTemp, iter, order, setting);
			break;

		case OPT_LOG:
			log_taylor_only_remainder(result, intTemp, iter, order);
			break;

		case OPT_SQRT:
			sqrt_taylor_only_remainder(result, intTemp, iter, order, setting);
			break;
		}
		break;
	}
	case NODE_BIN_OPT:
	{
		Interval remainder1, remainder2;

		switch(node_value.opt.type)
		{
		case OPT_PLUS:
			node_value.opt.left_operand->evaluate_remainder(remainder1, tms_of_vars, order, iter, setting);
			node_value.opt.right_operand->evaluate_remainder(remainder2, tms_of_vars, order, iter, setting);
			result = remainder1 + remainder2;
			break;
		case OPT_MINU:
			node_value.opt.left_operand->evaluate_remainder(remainder1, tms_of_vars, order, iter, setting);
			node_value.opt.right_operand->evaluate_remainder(remainder2, tms_of_vars, order, iter, setting);
			result = remainder1 - remainder2;
			break;
		case OPT_MULT:
			node_value.opt.left_operand->evaluate_remainder(remainder1, tms_of_vars, order, iter, setting);
			node_value.opt.right_operand->evaluate_remainder(remainder2, tms_of_vars, order, iter, setting);

			result = (*iter) * remainder2;
			++iter;
			result += (*iter) * remainder1;
			result += remainder1 * remainder2;
			++iter;
			result += (*iter);
			++iter;

			break;

		case OPT_DIV:
		{
			node_value.opt.left_operand->evaluate_remainder(remainder1, tms_of_vars, order, iter, setting);
			node_value.opt.right_operand->evaluate_remainder(remainder2, tms_of_vars, order, iter, setting);

			Interval intTemp;
			rec_taylor_only_remainder(intTemp, remainder2, iter, order, setting);

			result = (*iter) * intTemp;
			++iter;
			result += (*iter) * remainder1;
			result += remainder1 * intTemp;
			++iter;
			result += (*iter);
			++iter;
			break;
		}

		case OPT_POW:
		{
			node_value.opt.left_operand->evaluate_remainder(result, tms_of_vars, order, iter, setting);

			int degree = (int)node_value.opt.right_operand->node_value.constant.toDouble();

			if(degree == 0)
			{
				result = 0;
			}
			else if(degree > 1)
			{
				Interval temp = result;

				for(int i = degree - 1; i > 0;)
				{
					if(i & 1)
					{
						Interval temp2;
						temp2 = (*iter) * temp;
						++iter;
						temp2 += (*iter) * result;
						temp2 += temp * result;
						++iter;
						temp2 += (*iter);
						++iter;

						result = temp2;
					}

					i >>= 1;

					if(i > 0)
					{
						Interval temp2;
						temp2 = (*iter) * temp;
						++iter;
						temp2 += (*iter) * temp;
						temp2 += temp * temp;
						++iter;
						temp2 += (*iter);
						++iter;

						temp = temp2;
					}
				}
			}

			break;
		}
		}

		break;
	}

	case NODE_VAR:
		if(node_value.var.type == VAR_ID)
		{
			result = tms_of_vars[node_value.var.id].remainder;
		}
		else
		{
			// for the other variable types
		}

		break;

	case NODE_CONST:
	{
		result = *iter;
		++iter;
		break;
	}
	}
}

template <class DATA_TYPE>
void AST_Node<DATA_TYPE>::output(std::string & expression, const Variables & variables) const
{
	switch(node_type)
	{
	case NODE_UNA_OPT:
	{
		std::string temp;
		node_value.opt.left_operand->output(temp, variables);

		switch(node_value.opt.type)
		{
		case OPT_NEG:
			expression = "(-" + temp + ")";
			break;

		case OPT_SIN:
			expression = "(SIN(" + temp + "))";
			break;

		case OPT_COS:
			expression = "(COS(" + temp + "))";
			break;

		case OPT_EXP:
			expression = "(EXP(" + temp + "))";
			break;

		case OPT_LOG:
			expression = "(LOG(" + temp + "))";
			break;

		case OPT_SQRT:
			expression = "(SQRT(" + temp + "))";
			break;
		}
		break;
	}
	case NODE_BIN_OPT:
	{
		std::string temp1, temp2;

		switch(node_value.opt.type)
		{
		case OPT_PLUS:
			node_value.opt.left_operand->output(temp1, variables);
			node_value.opt.right_operand->output(temp2, variables);
			expression = "(" + temp1 + "+" + temp2 + ")";
			break;
		case OPT_MINU:
			node_value.opt.left_operand->output(temp1, variables);
			node_value.opt.right_operand->output(temp2, variables);
			expression = "(" + temp1 + "-" + temp2 + ")";
			break;
		case OPT_MULT:
			node_value.opt.left_operand->output(temp1, variables);
			node_value.opt.right_operand->output(temp2, variables);
			expression = "(" + temp1 + "*" + temp2 + ")";
			break;
		case OPT_DIV:
			node_value.opt.left_operand->output(temp1, variables);
			node_value.opt.right_operand->output(temp2, variables);
			expression = "(" + temp1 + "/" + temp2 + ")";
			break;
		case OPT_POW:
			node_value.opt.left_operand->output(temp1, variables);
			expression = temp1 + "^" + std::to_string((int)node_value.opt.right_operand->node_value.constant.toDouble());
			break;
		}

		break;
	}

	case NODE_VAR:
		if(node_value.var.type == VAR_ID)
		{
			variables.getVarName(expression, node_value.var.id);
		}
		else
		{
			variables.getVarName(expression, node_value.var.id);
		}

		break;

	case NODE_CONST:
		expression = node_value.constant.toString();
		break;
	}
}

template <>
inline void AST_Node<Interval>::toReal(std::shared_ptr<AST_Node<Real> > & pNode) const
{
	switch(node_type)
	{
	case NODE_CONST:
		pNode = std::shared_ptr<AST_Node<Real> > (new AST_Node<Real>(node_value.constant.toReal()));
		break;

	case NODE_VAR:
		pNode = std::shared_ptr<AST_Node<Real> > (new AST_Node<Real>(node_value.var.type, node_value.var.id));
		break;

	case NODE_UNA_OPT:
	{
		std::shared_ptr<AST_Node<Real> > p_left_node;
		node_value.opt.left_operand->toReal(p_left_node);
		pNode = std::shared_ptr<AST_Node<Real> > (new AST_Node<Real>(node_value.opt.type, p_left_node));
		break;
	}

	case NODE_BIN_OPT:
	{
		std::shared_ptr<AST_Node<Real> > p_left_node, p_right_node;
		node_value.opt.left_operand->toReal(p_left_node);
		node_value.opt.right_operand->toReal(p_right_node);
		pNode = std::shared_ptr<AST_Node<Real> > (new AST_Node<Real>(node_value.opt.type, p_left_node, p_right_node));
		break;
	}
	}
}





// abstract syntax tree
template <class DATA_TYPE>
class Expression_AST
{
protected:
	std::shared_ptr<AST_Node<DATA_TYPE> > root;

public:
	Expression_AST();
//	Expression_AST(const std::string & varName, const Taylor_Model_Computation_Setting & setting);
	Expression_AST(const std::string & varName, const Variables & variables);
	Expression_AST(const DATA_TYPE & c);
	~Expression_AST();

	void toReal(Expression_AST<Real> & expression) const;

	// using lex
	Expression_AST(const std::string & strExpression);
//	Expression_AST(const std::string & strExpression, const Variables & variables, const Parameters & parameters);

	void clear();

//	void output(FILE *fp, const Taylor_Model_Computation_Setting & setting) const;
	void output(std::ostream & os, const Variables & variables) const;

	void evaluate(Interval & result, const std::vector<Interval> & domain) const;

	template <class DATA_TYPE2>
	void evaluate(TaylorModel<DATA_TYPE2> & result, const std::vector<TaylorModel<DATA_TYPE2> > & tms_of_vars, const unsigned int order, const std::vector<Interval> & step_exp_table, const Interval & cutoff_threshold, const unsigned int numVars, const Global_Computation_Setting & setting) const;

	template <class DATA_TYPE2>
	void evaluate(TaylorModel<DATA_TYPE2> & result, const std::vector<TaylorModel<DATA_TYPE2> > & tms_of_vars, const unsigned int order, const std::vector<Interval> & domain, const Interval & cutoff_threshold, const Global_Computation_Setting & setting) const;

	template <class DATA_TYPE2>
	void evaluate_no_remainder(TaylorModel<DATA_TYPE2> & result, const std::vector<TaylorModel<DATA_TYPE2> > & tms_of_vars, const unsigned int order, const Interval & cutoff_threshold, const unsigned int numVars) const;

	template <class DATA_TYPE2>
	void evaluate(TaylorModel<DATA_TYPE2> & result, const std::vector<TaylorModel<DATA_TYPE2> > & tms_of_vars, const unsigned int order, const std::vector<Interval> & step_exp_table, const Interval & cutoff_threshold, const unsigned int numVars, std::list<Interval> & intermediate_ranges, const Global_Computation_Setting & setting) const;

	template <class DATA_TYPE2>
	void evaluate_remainder(Interval & result, const std::vector<TaylorModel<DATA_TYPE2> > & tms_of_vars, const unsigned int order, std::list<Interval>::iterator & iter, const Global_Computation_Setting & setting) const;

	bool isConstant(DATA_TYPE & c) const;

	Expression_AST & operator = (const Expression_AST & expression);
	Expression_AST & operator += (const Expression_AST & expression);
	Expression_AST & operator -= (const Expression_AST & expression);
	Expression_AST & operator *= (const Expression_AST & expression);
	Expression_AST & operator /= (const Expression_AST & expression);

	void pow_assign(const unsigned int n);
	void inv_assign();
	void sin_assign();
	void cos_assign();
	void exp_assign();
	void log_assign();
	void sqrt_assign();

	template <class DATA_TYPE2>
	friend class Expression_AST;
};


template <class DATA_TYPE>
Expression_AST<DATA_TYPE>::Expression_AST()
{
	root = nullptr;
}

template <class DATA_TYPE>
Expression_AST<DATA_TYPE>::Expression_AST(const std::string & varName, const Variables & variables)
{
	int id = variables.getIDForVar(varName);

	if(id < 0)
	{
		printf("%s is not declared.\n", varName.c_str());
		exit(0);
	}
	else
	{
		root = std::shared_ptr<AST_Node<DATA_TYPE> > (new AST_Node<DATA_TYPE>(VAR_ID, id));
	}
}

template <class DATA_TYPE>
Expression_AST<DATA_TYPE>::Expression_AST(const DATA_TYPE & c)
{
	root = std::shared_ptr<AST_Node<DATA_TYPE> > (new AST_Node<DATA_TYPE>(c));
}

template <class DATA_TYPE>
Expression_AST<DATA_TYPE>::~Expression_AST()
{
	root.reset();
}

template <>
inline void Expression_AST<Interval>::toReal(Expression_AST<Real> & expression) const
{
	root->toReal(expression.root);
}

template <>
inline Expression_AST<Interval>::Expression_AST(const std::string & strExpression)
{
	expression_ast_setting.clear();

	std::string prefix(str_prefix_expression_ast);
	std::string suffix(str_suffix);

	expression_ast_setting.strExpression = prefix + strExpression + suffix;

	parseExpression();

	*this = expression_ast_setting.result;
}

template <>
inline Expression_AST<Real>::Expression_AST(const std::string & strExpression)
{
	expression_ast_setting.clear();

	std::string prefix(str_prefix_expression_ast);
	std::string suffix(str_suffix);

	expression_ast_setting.strExpression = prefix + strExpression + suffix;

	parseExpression();

	expression_ast_setting.result.toReal(*this);
}

template <class DATA_TYPE>
void Expression_AST<DATA_TYPE>::clear()
{
	root.reset();
}

template <class DATA_TYPE>
void Expression_AST<DATA_TYPE>::output(std::ostream & os, const Variables & variables) const
{
	std::string expression;
	root->output(expression, variables);

	os << expression;
}

template <class DATA_TYPE>
void Expression_AST<DATA_TYPE>::evaluate(Interval & result, const std::vector<Interval> & domain) const
{
	root->evaluate(result, domain);
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void Expression_AST<DATA_TYPE>::evaluate(TaylorModel<DATA_TYPE2> & result, const std::vector<TaylorModel<DATA_TYPE2> > & tms_of_vars, const unsigned int order, const std::vector<Interval> & step_exp_table, const Interval & cutoff_threshold, const unsigned int numVars, const Global_Computation_Setting & setting) const
{
	root->evaluate(result, tms_of_vars, order, step_exp_table, cutoff_threshold, numVars, setting);
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void Expression_AST<DATA_TYPE>::evaluate(TaylorModel<DATA_TYPE2> & result, const std::vector<TaylorModel<DATA_TYPE2> > & tms_of_vars, const unsigned int order, const std::vector<Interval> & domain, const Interval & cutoff_threshold, const Global_Computation_Setting & setting) const
{
	root->evaluate(result, tms_of_vars, order, domain, cutoff_threshold, setting);
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void Expression_AST<DATA_TYPE>::evaluate_no_remainder(TaylorModel<DATA_TYPE2> & result, const std::vector<TaylorModel<DATA_TYPE2> > & tms_of_vars, const unsigned int order, const Interval & cutoff_threshold, const unsigned int numVars) const
{
	root->evaluate_no_remainder(result, tms_of_vars, order, cutoff_threshold, numVars);
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void Expression_AST<DATA_TYPE>::evaluate(TaylorModel<DATA_TYPE2> & result, const std::vector<TaylorModel<DATA_TYPE2> > & tms_of_vars, const unsigned int order, const std::vector<Interval> & step_exp_table, const Interval & cutoff_threshold, const unsigned int numVars, std::list<Interval> & intermediate_ranges, const Global_Computation_Setting & setting) const
{
	root->evaluate(result, tms_of_vars, order, step_exp_table, cutoff_threshold, numVars, intermediate_ranges, setting);
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void Expression_AST<DATA_TYPE>::evaluate_remainder(Interval & result, const std::vector<TaylorModel<DATA_TYPE2> > & tms_of_vars, const unsigned int order, std::list<Interval>::iterator & iter, const Global_Computation_Setting & setting) const
{
	root->evaluate_remainder(result, tms_of_vars, order, iter, setting);
}

template <class DATA_TYPE>
bool Expression_AST<DATA_TYPE>::isConstant(DATA_TYPE & c) const
{
	if(root->node_type == NODE_CONST)
	{
		c = root->node_value.constant;
		return true;
	}
	else
	{
		return false;
	}

}

template <class DATA_TYPE>
Expression_AST<DATA_TYPE> & Expression_AST<DATA_TYPE>::operator = (const Expression_AST<DATA_TYPE> & expression)
{
	if(this == &expression)
		return *this;

	root = expression.root;

	return *this;
}

template <class DATA_TYPE>
Expression_AST<DATA_TYPE> & Expression_AST<DATA_TYPE>::operator += (const Expression_AST<DATA_TYPE> & expression)
{
	std::shared_ptr<AST_Node<DATA_TYPE> > tmp(new AST_Node<DATA_TYPE>(OPT_PLUS, root, expression.root));
	root = tmp;

	return *this;
}

template <class DATA_TYPE>
Expression_AST<DATA_TYPE> & Expression_AST<DATA_TYPE>::operator -= (const Expression_AST<DATA_TYPE> & expression)
{
	std::shared_ptr<AST_Node<DATA_TYPE> > tmp(new AST_Node<DATA_TYPE>(OPT_MINU, root, expression.root));
	root = tmp;

	return *this;
}

template <class DATA_TYPE>
Expression_AST<DATA_TYPE> & Expression_AST<DATA_TYPE>::operator *= (const Expression_AST<DATA_TYPE> & expression)
{
	std::shared_ptr<AST_Node<DATA_TYPE> > tmp(new AST_Node<DATA_TYPE>(OPT_MULT, root, expression.root));
	root = tmp;

	return *this;
}

template <class DATA_TYPE>
Expression_AST<DATA_TYPE> & Expression_AST<DATA_TYPE>::operator /= (const Expression_AST<DATA_TYPE> & expression)
{
	std::shared_ptr<AST_Node<DATA_TYPE> > tmp(new AST_Node<DATA_TYPE>(OPT_DIV, root, expression.root));
	root = tmp;

	return *this;
}

template <class DATA_TYPE>
void Expression_AST<DATA_TYPE>::pow_assign(const unsigned int n)
{
	std::shared_ptr<AST_Node<DATA_TYPE> > exponent(new AST_Node<DATA_TYPE>(n));

	std::shared_ptr<AST_Node<DATA_TYPE> > tmp(new AST_Node<DATA_TYPE>(OPT_POW, root, exponent));
	root = tmp;
}

template <class DATA_TYPE>
void Expression_AST<DATA_TYPE>::inv_assign()
{
	std::shared_ptr<AST_Node<DATA_TYPE> > tmp(new AST_Node<DATA_TYPE>(OPT_NEG, root));
	root = tmp;
}

template <class DATA_TYPE>
void Expression_AST<DATA_TYPE>::sin_assign()
{
	std::shared_ptr<AST_Node<DATA_TYPE> > tmp(new AST_Node<DATA_TYPE>(OPT_SIN, root));
	root = tmp;
}

template <class DATA_TYPE>
void Expression_AST<DATA_TYPE>::cos_assign()
{
	std::shared_ptr<AST_Node<DATA_TYPE> > tmp(new AST_Node<DATA_TYPE>(OPT_COS, root));
	root = tmp;
}

template <class DATA_TYPE>
void Expression_AST<DATA_TYPE>::exp_assign()
{
	std::shared_ptr<AST_Node<DATA_TYPE> > tmp(new AST_Node<DATA_TYPE>(OPT_EXP, root));
	root = tmp;
}

template <class DATA_TYPE>
void Expression_AST<DATA_TYPE>::log_assign()
{
	std::shared_ptr<AST_Node<DATA_TYPE> > tmp(new AST_Node<DATA_TYPE>(OPT_LOG, root));
	root = tmp;
}

template <class DATA_TYPE>
void Expression_AST<DATA_TYPE>::sqrt_assign()
{
	std::shared_ptr<AST_Node<DATA_TYPE> > tmp(new AST_Node<DATA_TYPE>(OPT_SQRT, root));
	root = tmp;
}


}



#endif /* EXPRESSION_H_ */
