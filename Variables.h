/*---
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.
  The code is released as is under the GNU General Public License (GPL).
---*/

#ifndef VARIABLES_H_
#define VARIABLES_H_

#include "Interval.h"

namespace flowstar
{

class Variables
{
public:
	std::map<std::string,int> varTab;
	std::vector<std::string> varNames;

public:
	Variables();
	~Variables();
	Variables(const Variables & variables);

	Variables & operator = (const Variables & variables);

	int declareVar(const std::string & vName);
	int getIDForVar(const std::string & vName) const;
	bool getVarName(std::string & vName, const int id) const;
	unsigned int size() const;

	void output(std::ostream & os) const;

	void clear();
};



class Parameters
{
public:
	std::map<std::string,int> parTab;
	std::vector<std::string> parNames;
	std::vector<Real> parValues;

public:
	Parameters();
	~Parameters();
	Parameters(const Parameters & parameters);

	bool declarePar(const std::string & pName, const Real & value);
	int getIDForPar(const std::string & pName) const;

	bool getParName(std::string & pName, const int id) const;

	bool getParValue(Real & pValue, const std::string & pName) const;
	bool getParValue(Real & pValue, const int id) const;

	Parameters & operator = (const Parameters & parameters);

	unsigned int size() const;
	void clear();
};

}

#endif /* VARIABLES_H_ */
