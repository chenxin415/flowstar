/*---
  Flow*: A Verification Tool for Cyber-Physical Systems.
  Authors: Xin Chen, Sriram Sankaranarayanan, and Erika Abraham.
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.
  
  The code is released as is under the GNU General Public License (GPL).
---*/

#ifndef MODELPARSER_H_
#define MODELPARSER_H_

#include "Continuous.h"

using namespace flowstar;

extern int lineNum;
extern Continuous_Reachability_Problem_Description problem_description;
extern Continuous_Reachability reachability_for_outputFile;
extern Variables stateVars;
extern Variables tmVars;
extern Variables tvPars;


extern int yyparse();

void parseError(const char *str, int lnum);




class LTI_Term
{
public:
	flowstar::Interval coefficient;
	int varID;
	int parID;

public:
	LTI_Term()
	{
		varID = -1;
		parID = -1;
	}

	LTI_Term(const Interval & I, const int var_id, const int par_id)
	{
		coefficient = I;
		varID = var_id;
		parID = par_id;
	}

	~LTI_Term()
	{
	}

	LTI_Term & operator = (const LTI_Term & term)
	{
		if(this == &term)
			return *this;

		coefficient = term.coefficient;
		varID = term.varID;
		parID = term.parID;

		return *this;
	}
};



class LTI_ODE_description
{
public:
	Matrix<Real> 							dyn_A;
	Matrix<UnivariateTaylorModel<Real> >	dyn_B;

public:
	LTI_ODE_description(const unsigned int d)
	{
		Matrix<Real> A(d, d);
		dyn_A = A;

		Matrix<UnivariateTaylorModel<Real> > B(d, 1);
		dyn_B = B;
	}

	LTI_ODE_description(const LTI_ODE_description & description)
	{
		dyn_A = description.dyn_A;
		dyn_B = description.dyn_B;
	}

	~LTI_ODE_description()
	{
	}

	LTI_ODE_description & operator = (const LTI_ODE_description & description)
	{
		if(this == &description)
			return *this;

		dyn_A = description.dyn_A;
		dyn_B = description.dyn_B;

		return *this;
	}
};




class LTV_Term
{
public:
	UnivariatePolynomial<Real> coefficient;
	int varID;
	int tiParID;
	int tvParID;

public:
	LTV_Term()
	{
		varID	= -1;
		tiParID	= -1;
		tvParID	= -1;
	}

	LTV_Term(const UnivariatePolynomial<Real> & p, const int var_id, const int ti_par_id, const int tv_par_id)
	{
		coefficient	= p;
		varID		= var_id;
		tiParID		= ti_par_id;
		tvParID		= tv_par_id;
	}

	~LTV_Term()
	{
	}

	LTV_Term & operator = (const LTV_Term & term)
	{
		if(this == &term)
			return *this;

		coefficient	= term.coefficient;
		varID		= term.varID;
		tiParID		= term.tiParID;
		tvParID		= term.tvParID;

		return *this;
	}
};





class LTV_ODE_description
{
public:
	Matrix<UnivariatePolynomial<Real> >	dyn_A;
	Matrix<UnivariatePolynomial<Real> >	dyn_B;
	Matrix<UnivariatePolynomial<Real> >	dyn_tv;

public:
	LTV_ODE_description(const unsigned int numStateVars, const unsigned int numTVPars)
	{
		Matrix<UnivariatePolynomial<Real> > A(numStateVars, numStateVars), B(numStateVars, 1), C(numStateVars, numTVPars);

		dyn_A = A;
		dyn_B = B;
		dyn_tv = C;
	}

	LTV_ODE_description(const LTV_ODE_description & description)
	{
		dyn_A = description.dyn_A;
		dyn_B = description.dyn_B;
		dyn_tv = description.dyn_tv;
	}

	~LTV_ODE_description()
	{
	}

	LTV_ODE_description & operator = (const LTV_ODE_description & description)
	{
		if(this == &description)
			return *this;

		dyn_A = description.dyn_A;
		dyn_B = description.dyn_B;
		dyn_tv = description.dyn_tv;

		return *this;
	}
};





#endif /* MODELPARSER_H_ */
