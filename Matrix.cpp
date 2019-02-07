/*---
  Flow*: A Verification Tool for Cyber-Physical Systems.
  Authors: Xin Chen, Sriram Sankaranarayanan, and Erika Abraham.
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.
  
  The code is released as is under the GNU General Public License (GPL).
---*/

#include "Matrix.h"

using namespace flowstar;

namespace flowstar
{

MatrixParseSetting matrixParseSetting;


}









/*
// matrix for multivariate polynomials

mpMatrix::mpMatrix()
{
	size1 = 0;
	size2 = 0;
	data = NULL;
}

mpMatrix::mpMatrix(const int m, const int n)
{
	size1 = m;
	size2 = n;
	data = new Polynomial[size1 * size2];
}

mpMatrix::mpMatrix(const int n)
{
	size1 = n;
	size2 = n;
	data = new Polynomial[size1 * size2];
}

mpMatrix::mpMatrix(const mpMatrix & mpm)
{
	size1 = mpm.size1;
	size2 = mpm.size2;

	int size_total = size1 * size2;
	data = new Polynomial[size_total];

	std::copy(mpm.data, mpm.data + size_total, data);
}

mpMatrix::~mpMatrix()
{
	delete [] data;
}

int mpMatrix::rows() const
{
	return size1;
}

int mpMatrix::cols() const
{
	return size2;
}

void mpMatrix::intEval(mpMatrix & result, const std::vector<Interval> & val_exp_table) const
{
	delete [] result.data;
	result.size1 = size1;
	result.size2 = size2;

	int size_total = size1 * size2;
	result.data = new Polynomial[size_total];

	for(int i=0; i<size_total; ++i)
	{
		data[i].evaluate_t(result.data[i], val_exp_table);
	}
}

void mpMatrix::intEval(iMatrix & result, const std::vector<Interval> & domain) const
{
	delete [] result.data;
	result.size1 = size1;
	result.size2 = size2;

	int size_total = size1 * size2;
	result.data = new Interval[size_total];

	for(int i=0; i<size_total; ++i)
	{
		Interval tmp;
		data[i].intEval(tmp, domain);
		result.data[i] = tmp;
	}
}

void mpMatrix::output(FILE *fp, const std::vector<std::string> & varNames) const
{
	for(int i=0; i<size1; ++i)
	{
		for(int j=0; j<size2; ++j)
		{
			data[i*size2 + j].dump_interval(fp, varNames);
			fprintf(fp, "\t\t");
		}

		fprintf(fp, "\n");
	}
}

mpMatrix & mpMatrix::operator += (const mpMatrix & mpm)
{
	if(size1 != mpm.size1 || size2 != mpm.size2)
	{
		printf("Univariate polynomial matrix addition: Dimensions do not match.\n");
		exit(1);
	}

	int size_total = size1 * size2;

	for(int i=0; i<size_total; ++i)
	{
		data[i] += mpm.data[i];
	}

	return *this;
}

mpMatrix mpMatrix::operator + (const mpMatrix & mpm) const
{
	if(size1 != mpm.size1 || size2 != mpm.size2)
	{
		printf("Univariate polynomial matrix addition: Dimensions do not match.\n");
		exit(1);
	}

	mpMatrix result(size1, size2);
	int size_total = size1 * size2;

	for(int i=0; i<size_total; ++i)
	{
		result.data[i] = data[i] + mpm.data[i];
	}

	return result;
}

Polynomial * mpMatrix::operator [] (const int i)
{
	return &data[i * size2];
}

mpMatrix & mpMatrix::operator = (const mpMatrix & mpm)
{
	if(this == &mpm)
		return *this;

	size1 = mpm.size1;
	size2 = mpm.size2;

	int size_total = size1 * size2;
	delete [] data;

	if(size_total > 0)
	{
		data = new Polynomial[size_total];
		std::copy(mpm.data, mpm.data + size_total, data);
	}
	else
	{
		data = NULL;
	}

	return *this;
}
*/




MatrixParseSetting::MatrixParseSetting()
{
}

MatrixParseSetting::MatrixParseSetting(const MatrixParseSetting & setting)
{
	strExpression = setting.strExpression;
	result = setting.result;
}

MatrixParseSetting::~MatrixParseSetting()
{
}

MatrixParseSetting & MatrixParseSetting::operator = (const MatrixParseSetting & setting)
{
	if(this == &setting)
		return *this;

	strExpression = setting.strExpression;
	result = setting.result;

	return *this;
}





