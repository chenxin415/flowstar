/*---
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.
  The code is released as is under the GNU General Public License (GPL).
---*/

#ifndef HYBRID_H_
#define HYBRID_H_

#include "Continuous.h"


class Reset
{
public:
	std::vector<Expression_AST<Real> > resetMapping;
	std::vector<Interval> uncertainties;
	std::vector<bool> isIdentity;

public:
	Reset();
	Reset(const std::vector<Expression_AST<Real> > & expressions);
	Reset(const std::vector<Expression_AST<Real> > & expressions, const std::vector<bool> & ids);
	Reset(const std::vector<Expression_AST<Real> > & expressions, const std::vector<Interval> & uncs);
	Reset(const std::vector<Expression_AST<Real> > & expressions, const std::vector<Interval> & uncs, const std::vector<bool> & ids);
	Reset(const Reset & reset);
	~Reset();

	void reset(TaylorModelVec<Real> & result, const TaylorModelVec<Real> & tmv_input, const std::vector<Interval> & domain, const Taylor_Model_Computation_Setting & t_setting) const;

	Reset & operator = (const Reset & reset);
};



class DiscreteTransition
{
public:
	int jumpID;
	int startID;
	int targetID;

	std::vector<Constraint> guard;
	Reset reset;
public:
	DiscreteTransition();
	DiscreteTransition(const int id, const int start, const int target, const std::vector<Constraint> & constraints, const Reset & mapping);
	DiscreteTransition(const DiscreteTransition & jump);
	~DiscreteTransition();

	DiscreteTransition & operator = (const DiscreteTransition & jump);
};






class Hybrid_Dynamics
{
protected:
	std::vector<Dynamics> dynamics;
	std::vector

public:

};




#endif /* HYBRID_H_ */
