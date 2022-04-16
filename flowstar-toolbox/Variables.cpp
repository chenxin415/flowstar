/*---
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.
  The code is released as is under the GNU General Public License (GPL).
---*/

#include "Variables.h"

using namespace flowstar;

namespace flowstar
{

Variables::Variables()
{
}

Variables::~Variables()
{
	varTab.clear();
	varNames.clear();
}

Variables::Variables(const Variables & variables)
{
	varTab		= variables.varTab;
	varNames	= variables.varNames;
}

Variables & Variables::operator = (const Variables & variables)
{
	if(this == &variables)
		return *this;

	varTab		= variables.varTab;
	varNames	= variables.varNames;

	return *this;
}

int Variables::declareVar(const std::string & vName)
{
	std::map<std::string,int>::const_iterator iter;

	if((iter = varTab.find(vName)) == varTab.end())
	{
		int varID = varNames.size();
		varTab[vName] = varID;
		varNames.push_back(vName);
		return varID;
	}
	else
	{
		return -1;
	}
}

int Variables::getIDForVar(const std::string & vName) const
{
	std::map<std::string,int>::const_iterator iter;
	if((iter = varTab.find(vName)) == varTab.end())
	{
		return -1;
	}

	return iter->second;
}

bool Variables::getVarName(std::string & vName, const int id) const
{
	if(id >= 0 && id < varNames.size())
	{
		vName = varNames[id];
		return true;
	}
	else
	{
		return false;
	}
}

unsigned int Variables::size() const
{
	return varNames.size();
}

void Variables::output(std::ostream & os) const
{
	for(unsigned int i=0; i<varNames.size(); ++i)
	{
		os << i << "\t" << varNames[i] << std::endl;
	}
}

void Variables::clear()
{
	varTab.clear();
	varNames.clear();
}
















Parameters::Parameters()
{
}

Parameters::~Parameters()
{
	parTab.clear();
	parNames.clear();
	parValues.clear();
}

Parameters::Parameters(const Parameters & parameters)
{
	parTab = parameters.parTab;
	parNames = parameters.parNames;
	parValues = parameters.parValues;
}

bool Parameters::declarePar(const std::string & pName, const Real & value)
{
	std::map<std::string,int>::const_iterator iter;

	if((iter = parTab.find(pName)) == parTab.end())
	{
		parTab[pName] = parNames.size();
		parNames.push_back(pName);
		parValues.push_back(value);
		return true;
	}
	else
	{
		return false;
	}
}

int Parameters::getIDForPar(const std::string & pName) const
{
	std::map<std::string,int>::const_iterator iter;
	if((iter = parTab.find(pName)) == parTab.end())
	{
		return -1;
	}

	return iter->second;
}

bool Parameters::getParName(std::string & pName, const int id) const
{
	if(id >= 0 && id < parNames.size())
	{
		pName = parNames[id];
		return true;
	}
	else
	{
		return false;
	}
}

bool Parameters::getParValue(Real & pValue, const std::string & pName) const
{
	int id = getIDForPar(pName);

	if(id >= 0)
	{
		pValue = parValues[id];
		return true;
	}
	else
	{
		return false;
	}
}

bool Parameters::getParValue(Real & pValue, const int id) const
{
	if(id >= 0 && id < parNames.size())
	{
		pValue = parValues[id];
		return true;
	}
	else
	{
		return false;
	}
}

Parameters & Parameters::operator = (const Parameters & parameters)
{
	if(this == &parameters)
		return *this;

	parTab = parameters.parTab;
	parNames = parameters.parNames;
	parValues = parameters.parValues;

	return *this;
}

unsigned int Parameters::size() const
{
	return parNames.size();
}

void Parameters::clear()
{
	parTab.clear();
	parNames.clear();
	parValues.clear();
}


}
