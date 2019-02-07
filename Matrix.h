/*---
  Flow*: A Verification Tool for Cyber-Physical Systems.
  Authors: Xin Chen, Sriram Sankaranarayanan, and Erika Abraham.
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.
  
  The code is released as is under the GNU General Public License (GPL).
---*/

#ifndef MATRIX_H_
#define MATRIX_H_

#include "include.h"
#include "Interval.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include "UnivariateTaylorModel.h"

namespace flowstar
{

class Real;
class Interval;

template <class DATA_TYPE>
class TaylorModelVec;

template <class DATA_TYPE>
class Matrix;

template <class DATA_TYPE>
std::ostream & operator << (std::ostream & os, const Matrix<DATA_TYPE> & A);

template <class DATA_TYPE>
class Matrix
{
protected:
	DATA_TYPE *data;
	unsigned int size1;
	unsigned int size2;

public:
	Matrix();
	Matrix(const unsigned int m, const unsigned int n);		// Create an m x n matrix, all of the entries are 0.
	Matrix(const unsigned int n);							// Create an n x n identity matrix.
	Matrix(const unsigned int m, const unsigned int n, const DATA_TYPE & value);
	Matrix(const Matrix<DATA_TYPE> & A);
	~Matrix();

	unsigned int rows() const;
	unsigned int cols() const;

	bool getRowVec(Matrix<DATA_TYPE> & vec, const unsigned int index) const;
	bool getColVec(Matrix<DATA_TYPE> & vec, const unsigned int index) const;

	DATA_TYPE * getRowVecRef(const unsigned int index);

	bool isZero() const;

	void sortColumns();										// Sort the columns by size in descending order.
	unsigned int rank() const;

	void inverse(Matrix<DATA_TYPE> & result) const;
	void inverse_assign();

	void transpose(Matrix<DATA_TYPE> & result) const;
	void svd(Matrix<DATA_TYPE> & U) const;

	template <class DATA_TYPE2>
	void right_scale_assign(const std::vector<DATA_TYPE2> & scalars);

	// the matrix should also be a row vector
	double EuclideanNorm() const;
	void normalize();
	void max_norm(Real & norm) const;
	double max_norm() const;
	double min_entry() const;
	double innerProd(const Matrix<DATA_TYPE> & vec) const;


	// should be a matrix of intervals
	double width() const;
	bool isSingle() const;
	void toReal(Matrix<Real> & result) const;


	// should be a matrix of polynomials
	void integral();
	void times_x(const unsigned int order);
	unsigned int degree() const;

	template <class DATA_TYPE2>
	void ctrunc(Matrix<DATA_TYPE2> & rem1, Matrix<DATA_TYPE2> & rem2, const unsigned int order, const std::vector<DATA_TYPE2> & val1_exp_table, const std::vector<DATA_TYPE2> & val2_exp_table);

	template <class DATA_TYPE2, class DATA_TYPE3>
	void evaluate(Matrix<DATA_TYPE2> & result, const std::vector<DATA_TYPE3> & val_exp_table) const;

	void substitute(Matrix<DATA_TYPE> & result, const DATA_TYPE & x) const;


	// should be a matrix of Taylor models
	template <class DATA_TYPE2>
	void ctrunc(const unsigned int order, const std::vector<DATA_TYPE2> & val_exp_table);

	void ctrunc(const unsigned int order);
	void nctrunc(const unsigned int order);
	bool setRemainder(const Matrix<Interval> & remainder);
	bool getRemainder(Matrix<Interval> & remainder) const;
	bool addRemainder(const Matrix<Interval> & remainder);

	template <class DATA_TYPE2>
	void integral(const DATA_TYPE2 & val);

	template <class DATA_TYPE2>
	void evaluate(Matrix<DATA_TYPE2> & result) const;

	double remainderSize() const;

	DATA_TYPE * operator [] (const int i);

	template <class DATA_TYPE1, class DATA_TYPE2>
	friend Matrix<DATA_TYPE1> & operator += (Matrix<DATA_TYPE1> & A, const Matrix<DATA_TYPE2> & B);

	template <class DATA_TYPE1, class DATA_TYPE2>
	friend Matrix<DATA_TYPE1> & operator -= (Matrix<DATA_TYPE1> & A, const Matrix<DATA_TYPE2> & B);

	template <class DATA_TYPE1, class DATA_TYPE2>
	friend Matrix<DATA_TYPE1> & operator *= (Matrix<DATA_TYPE1> & A, const Matrix<DATA_TYPE2> & B);



	template <class DATA_TYPE1, class DATA_TYPE2>
	friend Matrix<DATA_TYPE1> & operator *= (Matrix<DATA_TYPE1> & A, const DATA_TYPE2 & c);

	template <class DATA_TYPE1, class DATA_TYPE2>
	friend Matrix<DATA_TYPE1> & operator /= (Matrix<DATA_TYPE1> & A, const DATA_TYPE2 & c);



	template <class DATA_TYPE1, class DATA_TYPE2>
	friend Matrix<DATA_TYPE1> operator + (const Matrix<DATA_TYPE1> & A, const Matrix<DATA_TYPE2> & B);

//	template <class DATA_TYPE1, class DATA_TYPE2>
//	friend Matrix<DATA_TYPE2> operator + (const Matrix<DATA_TYPE1> & A, const Matrix<DATA_TYPE2> & B);

	template <class DATA_TYPE1, class DATA_TYPE2>
	friend Matrix<DATA_TYPE1> operator - (const Matrix<DATA_TYPE1> & A, const Matrix<DATA_TYPE2> & B);

//	template <class DATA_TYPE1, class DATA_TYPE2>
//	friend Matrix<DATA_TYPE2> operator - (const Matrix<DATA_TYPE1> & A, const Matrix<DATA_TYPE2> & B);


	template <class DATA_TYPE1, class DATA_TYPE2>
	friend Matrix<DATA_TYPE1> operator * (const Matrix<DATA_TYPE1> & A, const Matrix<DATA_TYPE2> & B);

	template <class DATA_TYPE1, class DATA_TYPE2>
	friend std::vector<DATA_TYPE1> operator * (const Matrix<DATA_TYPE2> & A, const std::vector<DATA_TYPE1> & vec);

	friend Matrix<Interval> operator * (const Matrix<Real> & A, const Matrix<Interval> & B);
	friend Matrix<Interval> operator * (const Matrix<Real> & A, const Matrix<Interval> & B);

//	template <class DATA_TYPE1, class DATA_TYPE2>
//	friend Matrix<DATA_TYPE2> operator * (const Matrix<DATA_TYPE1> & A, const Matrix<DATA_TYPE2> & B);



	template <class DATA_TYPE1, class DATA_TYPE2>
	friend Matrix<DATA_TYPE1> operator * (const Matrix<DATA_TYPE1> & A, const DATA_TYPE2 & c);

	template <class DATA_TYPE1, class DATA_TYPE2>
	friend Matrix<DATA_TYPE1> operator * (const DATA_TYPE2 & c, const Matrix<DATA_TYPE1> & A);

	friend Matrix<UnivariateTaylorModel<Real> > operator * (const Matrix<Real> & A, const Matrix<UnivariateTaylorModel<Real> > & B);
	friend Matrix<UnivariateTaylorModel<Real> > operator * (const Matrix<UnivariatePolynomial<Real> > & A, const Matrix<UnivariateTaylorModel<Real> > & B);

	template <class DATA_TYPE1, class DATA_TYPE2>
	friend Matrix<DATA_TYPE1> operator / (const Matrix<DATA_TYPE1> & A, const DATA_TYPE2 & c);



	Matrix<DATA_TYPE> & operator = (const Matrix<DATA_TYPE> & A);

	template <class DATA_TYPE2>
	operator Matrix<DATA_TYPE2> () const;

	friend std::ostream & operator << <DATA_TYPE> (std::ostream & os, const Matrix<DATA_TYPE> & A);

	friend TaylorModelVec<DATA_TYPE> operator * (const Matrix<DATA_TYPE> & A, const TaylorModelVec<DATA_TYPE> & tmv);

	template <class DATA_TYPE2>
	friend class Matrix;
};


template <class DATA_TYPE>
Matrix<DATA_TYPE>::Matrix()
{
	data = NULL;
	size1 = 0;
	size2 = 0;
}

template <class DATA_TYPE>
Matrix<DATA_TYPE>::Matrix(const unsigned int m, const unsigned int n)
{
	size1 = m;
	size2 = n;

	data = new DATA_TYPE[m * n]();
}

template <class DATA_TYPE>
Matrix<DATA_TYPE>::Matrix(const unsigned int n)
{
	size1 = n;
	size2 = n;

	data = new DATA_TYPE[n * n]();

	for(int i=0, pos=0; i<n; ++i, pos+=n)
	{
		data[pos + i] = (DATA_TYPE)1;
	}
}

template <class DATA_TYPE>
Matrix<DATA_TYPE>::Matrix(const unsigned int m, const unsigned int n, const DATA_TYPE & value)
{
	size1 = m;
	size2 = n;

	unsigned int wholeSize = size1 * size2;

	data = new DATA_TYPE[wholeSize];

	for(unsigned int i=0; i<wholeSize; ++i)
	{
		data[i] = value;
	}
}

template <class DATA_TYPE>
Matrix<DATA_TYPE>::Matrix(const Matrix<DATA_TYPE> & A)
{
	size1 = A.size1;
	size2 = A.size2;

	unsigned int wholeSize = size1 * size2;
	data = new DATA_TYPE[wholeSize];
	std::copy(A.data, A.data + wholeSize, data);
}

template <class DATA_TYPE>
Matrix<DATA_TYPE>::~Matrix()
{
	delete [] data;
}

template <class DATA_TYPE>
unsigned int Matrix<DATA_TYPE>::rows() const
{
	return size1;
}

template <class DATA_TYPE>
unsigned int Matrix<DATA_TYPE>::cols() const
{
	return size2;
}

template <class DATA_TYPE>
bool Matrix<DATA_TYPE>::getRowVec(Matrix<DATA_TYPE> & vec, const unsigned int index) const
{
	if(index < 0 || index >= size1)
	{
		return false;
	}
	else
	{
		unsigned int pos = index * size2;
		for(unsigned int j=0; j<size2; ++j)
		{
			vec[0][j] = data[pos + j];
		}

		return true;
	}
}

template <class DATA_TYPE>
bool Matrix<DATA_TYPE>::getColVec(Matrix<DATA_TYPE> & vec, const unsigned int index) const
{
	if(index < 0 || index >= size2)
	{
		return false;
	}
	else
	{
		for(unsigned int i=0, pos=index; i<size1; ++i, pos+=size2)
		{
			vec[i][0] = data[pos];
		}

		return true;
	}
}

template <class DATA_TYPE>
DATA_TYPE * Matrix<DATA_TYPE>::getRowVecRef(const unsigned int index)
{
	if(index < 0 || index >= size1)
	{
		return NULL;
	}
	else
	{
		return data + index * size2;
	}
}

template <class DATA_TYPE>
bool Matrix<DATA_TYPE>::isZero() const
{
	bool result = true;

	unsigned int wholeSize = size1 * size2;

	for(unsigned int i=0; i<wholeSize; ++i)
	{
		if(data[i] != 0)
		{
			result = false;
			break;
		}
	}

	return result;
}

template <>
bool inline Matrix<UnivariatePolynomial<Real> >::isZero() const
{
	bool result = true;

	unsigned int wholeSize = size1 * size2;

	for(unsigned int i=0; i<wholeSize; ++i)
	{
		if(!(data[i].isZero()))
		{
			result = false;
			break;
		}
	}

	return result;
}

template <class DATA_TYPE>
void Matrix<DATA_TYPE>::sortColumns()
{
	unsigned int wholeSize = data->size1 * data->size2;

	double *sizes = new double[size2];

	// Compute the sizes of the columns
	for(unsigned int j=0; j<size2; ++j)
	{
		double size = 0;

		for(unsigned int pos=0; pos<wholeSize; pos+=size2)
		{
			size += (double)(data[pos + j] * data[pos + j]);
		}

		sizes[j] = size;
	}

	// Selection sort
	unsigned int iMax;

	for(unsigned int i=0; i<size2-1; ++i)
	{
		iMax = i;
		for(unsigned int j=i+1; j<size2; ++j)
		{
			if(sizes[j] > sizes[iMax])
			{
				iMax = j;
			}
		}

		//Exchange the columns
		if(iMax != i)
		{
			unsigned int pos_i = i * size2;
			unsigned int pos_iMax = iMax * size2;

			// swap the columns
			for(unsigned int j=0; j<size2; ++j)
			{
				DATA_TYPE tmp = data[pos_i + j];
				data[pos_i + j] = data[pos_iMax + j];
				data[pos_iMax + j] = tmp;
			}

			// swap the sizes
			double tmp = sizes[i];
			sizes[i] = sizes[iMax];
			sizes[iMax] = tmp;
		}
	}

	delete[] sizes;
}

template <class DATA_TYPE>
unsigned int Matrix<DATA_TYPE>::rank() const
{
	// the matrix is assumed to be a square matrix
	// we use the GSL library

	gsl_matrix *temp = gsl_matrix_alloc(size1, size2);

	for(unsigned int i=0, pos=0; i<size1; ++i, pos+=size2)
	{
		for(unsigned int j=0; j<size2; ++j)
		{
			gsl_matrix_set(temp, i, j, (double)data[pos + j]);
		}
	}

	gsl_vector *work = gsl_vector_alloc(size2);
	gsl_vector *S = gsl_vector_alloc(size2);
	gsl_matrix *V = gsl_matrix_alloc(size2, size2);

	gsl_linalg_SV_decomp(temp, V, S, work);

	unsigned int r = 0;
	double tmp;
	for(int i=0; i<size2; ++i)
	{
		tmp = gsl_vector_get(S, i);
		if(tmp < THRESHOLD_HIGH)
			break;
		else
			++r;
	}

	gsl_matrix_free(temp);
	gsl_vector_free(S);
	gsl_vector_free(work);
	gsl_matrix_free(V);

	return r;
}

template <class DATA_TYPE>
void Matrix<DATA_TYPE>::inverse(Matrix<DATA_TYPE> & result) const
{
	if(size1 != size2)
	{
		printf("Not a square matrix.\n");
		return;
	}

	// use the GSL library.
	gsl_matrix *A = gsl_matrix_alloc(size1, size1);
	gsl_permutation *p = gsl_permutation_alloc(size1);
	gsl_matrix *invA = gsl_matrix_alloc(size1, size1);

	// make a copy the matrix
	for(unsigned int i=0, pos=0; i<size1; ++i, pos+=size2)
	{
		for(unsigned int j=0; j<size2; ++j)
		{
			gsl_matrix_set(A, i, j, (double)data[pos + j]);
		}
	}

	int *signum = new int[size1];

	gsl_linalg_LU_decomp(A, p, signum);
	gsl_linalg_LU_invert(A, p, invA);

	Matrix<DATA_TYPE> matTmp(size1, size2);

	for(unsigned int i=0, pos=0; i<size1; ++i, pos+=size2)
	{
		for(unsigned int j=0; j<size2; ++j)
		{
			matTmp.data[pos + j] = gsl_matrix_get(invA, i, j);
		}
	}

	result = matTmp;

	gsl_matrix_free(A);
	gsl_permutation_free(p);
	gsl_matrix_free(invA);
	delete[] signum;
}

template <class DATA_TYPE>
void Matrix<DATA_TYPE>::inverse_assign()
{
	if(size1 != size2)
	{
		printf("Not a square matrix.\n");
		return;
	}

	// use the GSL library.
	gsl_matrix *A = gsl_matrix_alloc(size1, size1);
	gsl_permutation *p = gsl_permutation_alloc(size1);
	gsl_matrix *invA = gsl_matrix_alloc(size1, size1);

	// make a copy the matrix
	for(unsigned int i=0, pos=0; i<size1; ++i, pos+=size2)
	{
		for(unsigned int j=0; j<size2; ++j)
		{
			gsl_matrix_set(A, i, j, (double)data[pos + j]);
		}
	}

	int *signum = new int[size1];

	gsl_linalg_LU_decomp(A, p, signum);
	gsl_linalg_LU_invert(A, p, invA);

	for(unsigned int i=0, pos=0; i<size1; ++i, pos+=size2)
	{
		for(unsigned int j=0; j<size2; ++j)
		{
			data[pos + j] = gsl_matrix_get(invA, i, j);
		}
	}

	gsl_matrix_free(A);
	gsl_permutation_free(p);
	gsl_matrix_free(invA);
	delete[] signum;
}

template <class DATA_TYPE>
void Matrix<DATA_TYPE>::transpose(Matrix<DATA_TYPE> & result) const
{
	Matrix<DATA_TYPE> matTmp(size2, size1);

	for(unsigned int i=0, pos_i=0; i<size1; ++i, pos_i+=size2)
	{
		for(unsigned int j=0, pos_j=0; j<size2; ++j, pos_j+=size1)
		{
			matTmp.data[pos_j + i] = data[pos_i + j];
		}
	}

	result = matTmp;
}

template <class DATA_TYPE>
void Matrix<DATA_TYPE>::svd(Matrix<DATA_TYPE> & U) const
{
	gsl_matrix *temp = gsl_matrix_alloc(size1, size2);

	for(unsigned int i=0, pos=0; i<size1; ++i, pos+=size2)
	{
		for(unsigned int j=0; j<size2; ++j)
		{
			gsl_matrix_set(temp, i, j, (double)data[pos + j]);
		}
	}

	gsl_vector *work = gsl_vector_alloc(size2);
	gsl_vector *S = gsl_vector_alloc(size2);
	gsl_matrix *V = gsl_matrix_alloc(size2, size2);

	gsl_linalg_SV_decomp(temp, V, S, work);

	Matrix<DATA_TYPE> matTmp(size1, size2);

	for(unsigned int i=0, pos=0; i<size1; ++i, pos+=size2)
	{
		for(unsigned int j=0; j<size2; ++j)
		{
			matTmp.data[pos + j] = gsl_matrix_get(temp, i, j);
		}
	}

	U = matTmp;

	gsl_matrix_free(temp);
	gsl_vector_free(S);
	gsl_vector_free(work);
	gsl_matrix_free(V);
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void Matrix<DATA_TYPE>::right_scale_assign(const std::vector<DATA_TYPE2> & scalars)
{
	if(size2 != scalars.size())
	{
		printf("Matrix Multiplication: Dimensions do not match.\n");
		exit(1);
	}

	for(unsigned int i=0, pos=0; i<size1; ++i, pos+=size2)
	{
		for(unsigned int j=0; j<size2; ++j)
		{
			data[pos + j] *= scalars[j];
		}
	}
}

template <class DATA_TYPE>
double Matrix<DATA_TYPE>::EuclideanNorm() const
{
	double result = 0;

	for(unsigned int i=0; i<size2; ++i)
	{
		result += data[i] * data[i];
	}

	return sqrt(result);
}

template <class DATA_TYPE>
void Matrix<DATA_TYPE>::normalize()
{
	double norm = EuclideanNorm();

	(*this) /= norm;
}

template <class DATA_TYPE>
void Matrix<DATA_TYPE>::max_norm(Real & norm) const
{
	double max = 0;

	unsigned int wholeSize = size1 * size2;

	for(unsigned int i=0; i<wholeSize; ++i)
	{
		double tmp = data[i].mag();

		if(tmp > max)
		{
			max = tmp;
		}
	}

	norm.set(max);
}

template <class DATA_TYPE>
double Matrix<DATA_TYPE>::max_norm() const
{
	double max = 0;

	unsigned int wholeSize = size1 * size2;

	for(unsigned int i=0; i<wholeSize; ++i)
	{
		double tmp = data[i].mag();

		if(tmp > max)
		{
			max = tmp;
		}
	}

	return max;
}

template <class DATA_TYPE>
double Matrix<DATA_TYPE>::min_entry() const
{
	double min = 1e20;

	unsigned int wholeSize = size1 * size2;

	for(unsigned int i=0; i<wholeSize; ++i)
	{
		if(data[i] != 0)
		{
			double tmp = data[i].mag();

			if(tmp < min)
			{
				min = tmp;
			}
		}
	}

	return min;
}

template <class DATA_TYPE>
double Matrix<DATA_TYPE>::innerProd(const Matrix<DATA_TYPE> & vec) const
{
	unsigned int n = vec.cols();

	if(n != vec.cols())
	{
		printf("Vector dimensions do not match.\n");
		return INVALID;
	}

	double result = 0;
	for(unsigned int i=0; i<n; ++i)
	{
		result += data[i] * vec.data[i];;
	}

	return result;
}

template <class DATA_TYPE>
double Matrix<DATA_TYPE>::width() const
{
	unsigned int wholeSize = size1 * size2;
	double max_width = 0;

	for(unsigned int i=0; i<wholeSize; ++i)
	{
		double tmp = data[i].width();

		if(tmp > max_width)
		{
			max_width = tmp;
		}
	}

	return max_width;
}

template <class DATA_TYPE>
bool Matrix<DATA_TYPE>::isSingle() const
{
	unsigned int wholeSize = size1 * size2;

	for(unsigned int i=0; i<wholeSize; ++i)
	{
		if(!data[i].isSingle())
		{
			return false;
		}
	}

	return true;
}

template <class DATA_TYPE>
void Matrix<DATA_TYPE>::toReal(Matrix<Real> & result) const
{
	unsigned int wholeSize = size1 * size2;

	for(unsigned int i=0; i<wholeSize; ++i)
	{
		result.data[i] = data[i].toReal();
	}
}

template <class DATA_TYPE>
void Matrix<DATA_TYPE>::integral()
{
	unsigned int wholeSize = size1 * size2;

	for(unsigned int i=0; i<wholeSize; ++i)
	{
		data[i].integral();
	}
}

template <class DATA_TYPE>
void Matrix<DATA_TYPE>::times_x(const unsigned int order)
{
	unsigned int wholeSize = size1 * size2;

	for(unsigned int i=0; i<wholeSize; ++i)
	{
		data[i].times_x(order);
	}
}

template <class DATA_TYPE>
unsigned int Matrix<DATA_TYPE>::degree() const
{
	unsigned int wholeSize = size1 * size2;
	unsigned int max = 0;

	for(unsigned int i=0; i<wholeSize; ++i)
	{
		unsigned int d = data[i].degree();

		if(max < d)
			max = d;
	}

	return max;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void Matrix<DATA_TYPE>::ctrunc(Matrix<DATA_TYPE2> & rem1, Matrix<DATA_TYPE2> & rem2, const unsigned int order, const std::vector<DATA_TYPE2> & val1_exp_table, const std::vector<DATA_TYPE2> & val2_exp_table)
{
	unsigned int wholeSize = size1 * size2;

	delete [] rem1.data;
	rem1.size1 = size1;
	rem1.size2 = size2;
	rem1.data = new DATA_TYPE2[wholeSize]();

	delete [] rem2.data;
	rem2.size1 = size1;
	rem2.size2 = size2;
	rem2.data = new DATA_TYPE2[wholeSize]();

	for(unsigned int i=0; i<wholeSize; ++i)
	{
		data[i].ctrunc(rem1.data[i], rem2.data[i], order, val1_exp_table, val2_exp_table);
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2, class DATA_TYPE3>
void Matrix<DATA_TYPE>::evaluate(Matrix<DATA_TYPE2> & result, const std::vector<DATA_TYPE3> & val_exp_table) const
{
	unsigned int wholeSize = size1 * size2;

	for(unsigned int i=0; i<wholeSize; ++i)
	{
		DATA_TYPE2 tmp;
		data[i].evaluate(result.data[i], val_exp_table);
	}
}

template <class DATA_TYPE>
void Matrix<DATA_TYPE>::substitute(Matrix<DATA_TYPE> & result, const DATA_TYPE & x) const
{
	delete [] result.data;
	result.size1 = size1;
	result.size2 = size2;

	unsigned int wholeSize = size1 * size2;
	result.data = new UnivariatePolynomial<Real>[wholeSize];

	for(unsigned int i=0; i<wholeSize; ++i)
	{
		data[i].substitute(result.data[i], x);
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void Matrix<DATA_TYPE>::ctrunc(const unsigned int order, const std::vector<DATA_TYPE2> & val_exp_table)
{
	unsigned int wholeSize = size1 * size2;

	for(unsigned int i=0; i<wholeSize; ++i)
	{
		data[i].ctrunc(order, val_exp_table);
	}
}

template <class DATA_TYPE>
void Matrix<DATA_TYPE>::ctrunc(const unsigned int order)
{
	unsigned int wholeSize = size1 * size2;

	for(unsigned int i=0; i<wholeSize; ++i)
	{
		data[i].ctrunc(order);
	}
}

template <class DATA_TYPE>
void Matrix<DATA_TYPE>::nctrunc(const unsigned int order)
{
	unsigned int wholeSize = size1 * size2;

	for(unsigned int i=0; i<wholeSize; ++i)
	{
		data[i].nctrunc(order);
	}
}

template <class DATA_TYPE>
bool Matrix<DATA_TYPE>::setRemainder(const Matrix<Interval> & remainder)
{
	if(size1 != remainder.size1 || size2 != remainder.size2)
	{
		printf("Set Remainder: Dimensions do not match.\n");
		return false;
	}
	else
	{
		unsigned int wholeSize = size1 * size2;

		for(unsigned int i=0; i<wholeSize; ++i)
		{
			data[i].setRemainder(remainder.data[i]);
		}

		return true;
	}
}

template <class DATA_TYPE>
bool Matrix<DATA_TYPE>::getRemainder(Matrix<Interval> & remainder) const
{
	if(size1 != remainder.size1 || size2 != remainder.size2)
	{
		printf("Get Remainder: Dimensions do not match.\n");
		return false;
	}
	else
	{
		unsigned int wholeSize = size1 * size2;

		for(unsigned int i=0; i<wholeSize; ++i)
		{
			data[i].getRemainder(remainder.data[i]);
		}

		return true;
	}
}

template <class DATA_TYPE>
bool Matrix<DATA_TYPE>::addRemainder(const Matrix<Interval> & remainder)
{
	if(size1 != remainder.size1 || size2 != remainder.size2)
	{
		printf("Set Remainder: Dimensions do not match.\n");
		return false;
	}
	else
	{
		unsigned int wholeSize = size1 * size2;

		for(unsigned int i=0; i<wholeSize; ++i)
		{
			data[i].addRemainder(remainder.data[i]);
		}

		return true;
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void Matrix<DATA_TYPE>::integral(const DATA_TYPE2 & val)
{
	unsigned int wholeSize = size1 * size2;

	for(unsigned int i=0; i<wholeSize; ++i)
	{
		data[i].integral(val);
	}
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
void Matrix<DATA_TYPE>::evaluate(Matrix<DATA_TYPE2> & result) const
{
	unsigned int wholeSize = size1 * size2;

	for(unsigned int i=0; i<wholeSize; ++i)
	{
		DATA_TYPE2 tmp;
		data[i].evaluate(result.data[i], interval_utm_setting.val_exp_table);
	}
}

template <class DATA_TYPE>
double Matrix<DATA_TYPE>::remainderSize() const
{
	double max = 0;

	unsigned int wholeSize = size1 * size2;

	for(unsigned int i=0; i<wholeSize; ++i)
	{
		double tmp = data[i].remainderSize();

		if(tmp > max)
		{
			max = tmp;
		}
	}

	return max;
}

template <class DATA_TYPE>
DATA_TYPE * Matrix<DATA_TYPE>::operator [] (const int i)
{
	return &data[i * size2];
}


template <class DATA_TYPE1, class DATA_TYPE2>
Matrix<DATA_TYPE1> & operator += (Matrix<DATA_TYPE1> & A, const Matrix<DATA_TYPE2> & B)
{
	if(A.size1 != B.size1 || A.size2 != B.size2)
	{
		printf("Matrix Addition: Dimensions do not match.\n");
		exit(1);
	}

	unsigned int wholeSize = A.size1 * A.size2;

	for(unsigned int i=0; i<wholeSize; ++i)
	{
		A.data[i] += B.data[i];
	}

	return A;
}

template <class DATA_TYPE1, class DATA_TYPE2>
Matrix<DATA_TYPE1> & operator -= (Matrix<DATA_TYPE1> & A, const Matrix<DATA_TYPE2> & B)
{
	if(A.size1 != B.size1 || A.size2 != B.size2)
	{
		printf("Matrix Subtraction: Dimensions do not match.\n");
		exit(1);
	}

	unsigned int wholeSize = A.size1 * A.size2;

	for(unsigned int i=0; i<wholeSize; ++i)
	{
		A.data[i] -= B.data[i];
	}

	return A;
}

template <class DATA_TYPE1, class DATA_TYPE2>
Matrix<DATA_TYPE1> & operator *= (Matrix<DATA_TYPE1> & A, const Matrix<DATA_TYPE2> & B)
{
	if(A.size2 != B.size1)
	{
		printf("Matrix Multiplication: Dimensions do not match.\n");
		exit(1);
	}

	Matrix<DATA_TYPE1> result(A.size1, B.size2);

	for(unsigned int i=0, pos1=0, pos2=0; i<A.size1; ++i, pos1+=A.size2, pos2+=B.size2)
	{
		for(unsigned int j=0; j<B.size2; ++j)
		{
			DATA_TYPE1 tmp = 0;

			for(unsigned int k=0, pos3=0; k<A.size2; ++k, pos3+=B.size2)
			{
				tmp += A.data[pos1 + k] * B.data[pos3 + j];
			}

			result.data[pos2 + j] = tmp;
		}
	}

	A = result;
	return A;
}




template <class DATA_TYPE1, class DATA_TYPE2>
Matrix<DATA_TYPE1> & operator *= (Matrix<DATA_TYPE1> & A, const DATA_TYPE2 & c)
{
	unsigned int wholeSize = A.size1 * A.size2;

	for(unsigned int i=0; i<wholeSize; ++i)
	{
		A.data[i] *= c;
	}

	return A;
}

template <class DATA_TYPE1, class DATA_TYPE2>
Matrix<DATA_TYPE1> & operator /= (Matrix<DATA_TYPE1> & A, const DATA_TYPE2 & c)
{
	unsigned int wholeSize = A.size1 * A.size2;

	for(unsigned int i=0; i<wholeSize; ++i)
	{
		A.data[i] /= c;
	}

	return A;
}




template <class DATA_TYPE1, class DATA_TYPE2>
Matrix<DATA_TYPE1> operator + (const Matrix<DATA_TYPE1> & A, const Matrix<DATA_TYPE2> & B)
{
	Matrix<DATA_TYPE1> result = A;
	result += B;
	return result;
}
/*
template <class DATA_TYPE1, class DATA_TYPE2>
Matrix<DATA_TYPE2> operator + (const Matrix<DATA_TYPE1> & A, const Matrix<DATA_TYPE2> & B)
{
	Matrix<DATA_TYPE2> result = B;
	result += A;
	return result;
}
*/
template <class DATA_TYPE1, class DATA_TYPE2>
Matrix<DATA_TYPE1> operator - (const Matrix<DATA_TYPE1> & A, const Matrix<DATA_TYPE2> & B)
{
	Matrix<DATA_TYPE1> result = A;
	result -= B;
	return result;
}
/*
template <class DATA_TYPE1, class DATA_TYPE2>
Matrix<DATA_TYPE2> operator - (const Matrix<DATA_TYPE1> & A, const Matrix<DATA_TYPE2> & B)
{
	if(A.size1 != B.size1 || A.size2 != B.size2)
	{
		printf("Matrix Subtraction: Dimensions do not match.\n");
		exit(1);
	}

	Matrix<DATA_TYPE2> result;

	unsigned int wholeSize = A.size1 * A.size2;

	for(unsigned int i=0; i<wholeSize; ++i)
	{
		result.data[i] = A.data[i] - B.data[i];
	}

	return result;
}
*/

template <class DATA_TYPE1, class DATA_TYPE2>
Matrix<DATA_TYPE1> operator * (const Matrix<DATA_TYPE1> & A, const Matrix<DATA_TYPE2> & B)
{
	if(A.size2 != B.size1)
	{
		printf("Matrix Multiplication: Dimensions do not match.\n");
		exit(1);
	}

	Matrix<DATA_TYPE1> result(A.size1, B.size2);

	for(unsigned int i=0, pos1=0, pos2=0; i<A.size1; ++i, pos1+=A.size2, pos2+=B.size2)
	{
		for(unsigned int j=0; j<B.size2; ++j)
		{
			DATA_TYPE1 tmp = 0;

			for(unsigned int k=0, pos3=0; k<A.size2; ++k, pos3+=B.size2)
			{
				tmp += A.data[pos1 + k] * B.data[pos3 + j];
			}

			result.data[pos2 + j] = tmp;
		}
	}

	return result;
}

template <class DATA_TYPE1, class DATA_TYPE2>
std::vector<DATA_TYPE1> operator * (const Matrix<DATA_TYPE2> & A, const std::vector<DATA_TYPE1> & vec)
{
	if(A.size2 != vec.size())
	{
		printf("Matrix Multiplication: Dimensions do not match.\n");
		exit(1);
	}

	std::vector<DATA_TYPE1> result;

	for(unsigned int i=0, pos=0; i<A.size1; ++i, pos+=A.size2)
	{
		DATA_TYPE1 tmp;

		for(unsigned int j=0; j<A.size2; ++j)
		{
			tmp += vec[j] * A.data[pos + j];
		}

		result.push_back(tmp);
	}

	return result;
}

inline Matrix<Interval> operator * (const Matrix<Real> & A, const Matrix<Interval> & B)
{
	if(A.size2 != B.size1)
	{
		printf("Matrix Multiplication: Dimensions do not match.\n");
		exit(1);
	}

	Matrix<Interval> result(A.size1, B.size2);

	for(unsigned int i=0, pos1=0, pos2=0; i<A.size1; ++i, pos1+=A.size2, pos2+=B.size2)
	{
		for(unsigned int j=0; j<B.size2; ++j)
		{
			Interval tmp;

			for(unsigned int k=0, pos3=0; k<A.size2; ++k, pos3+=B.size2)
			{
				tmp += A.data[pos1 + k] * B.data[pos3 + j];
			}

			result.data[pos2 + j] = tmp;
		}
	}

	return result;
}

/*
template <class DATA_TYPE1, class DATA_TYPE2>
Matrix<DATA_TYPE2> operator * (const Matrix<DATA_TYPE1> & A, const Matrix<DATA_TYPE2> & B)
{
	if(A.size2 != B.size1)
	{
		printf("Matrix Multiplication: Dimensions do not match.\n");
		exit(1);
	}

	Matrix<DATA_TYPE2> result(A.size1, B.size2);

	for(unsigned int i=0, pos1=0, pos2=0; i<A.size1; ++i, pos1+=A.size2, pos2+=B.size2)
	{
		for(unsigned int j=0; j<B.size2; ++j)
		{
			DATA_TYPE2 tmp = 0;

			for(unsigned int k=0, pos3=0; k<A.size2; ++k, pos3+=B.size2)
			{
				tmp += A.data[pos1 + k] * B.data[pos3 + j];
			}

			result.data[pos2 + j] = tmp;
		}
	}

	return result;
}
*/


template <class DATA_TYPE1, class DATA_TYPE2>
Matrix<DATA_TYPE1> operator * (const Matrix<DATA_TYPE1> & A, const DATA_TYPE2 & c)
{
	Matrix<DATA_TYPE1> result = A;
	result *= c;
	return result;
}

template <class DATA_TYPE1, class DATA_TYPE2>
Matrix<DATA_TYPE1> operator * (const DATA_TYPE2 & c, const Matrix<DATA_TYPE1> & A)
{
	Matrix<DATA_TYPE1> result = A;
	result *= c;
	return result;
}

inline Matrix<UnivariateTaylorModel<Real> > operator * (const Matrix<Real> & A, const Matrix<UnivariateTaylorModel<Real> > & B)
{
	if(A.size2 != B.size1)
	{
		printf("Matrix Multiplication: Dimensions do not match.\n");
		exit(1);
	}

	Matrix<UnivariateTaylorModel<Real> > result(A.size1, B.size2);

	for(unsigned int i=0, pos1=0, pos2=0; i<A.size1; ++i, pos1+=A.size2, pos2+=B.size2)
	{
		for(unsigned int j=0; j<B.size2; ++j)
		{
			UnivariateTaylorModel<Real> tmp;

			for(unsigned int k=0, pos3=0; k<A.size2; ++k, pos3+=B.size2)
			{
				tmp += B.data[pos3 + j] * A.data[pos1 + k];
			}

			result.data[pos2 + j] = tmp;
		}
	}

	return result;
}

inline Matrix<UnivariateTaylorModel<Real> > operator * (const Matrix<UnivariatePolynomial<Real> > & A, const Matrix<UnivariateTaylorModel<Real> > & B)
{
	if(A.size2 != B.size1)
	{
		printf("Matrix Multiplication: Dimensions do not match.\n");
		exit(1);
	}

	Matrix<UnivariateTaylorModel<Real> > result(A.size1, B.size2);

	for(unsigned int i=0, pos1=0, pos2=0; i<A.size1; ++i, pos1+=A.size2, pos2+=B.size2)
	{
		for(unsigned int j=0; j<B.size2; ++j)
		{
			UnivariateTaylorModel<Real> tmp;

			for(unsigned int k=0, pos3=0; k<A.size2; ++k, pos3+=B.size2)
			{
				tmp += B.data[pos3 + j] * A.data[pos1 + k];
			}

			result.data[pos2 + j] = tmp;
		}
	}

	return result;
}

template <class DATA_TYPE1, class DATA_TYPE2>
Matrix<DATA_TYPE1> operator / (const Matrix<DATA_TYPE1> & A, const DATA_TYPE2 & c)
{
	Matrix<DATA_TYPE1> result = A;
	result /= c;
	return result;
}

template <class DATA_TYPE>
Matrix<DATA_TYPE> & Matrix<DATA_TYPE>::operator = (const Matrix<DATA_TYPE> & A)
{
	if(this == &A)
		return *this;

	size1 = A.size1;
	size2 = A.size2;

	unsigned int wholeSize = size1 * size2;
	delete [] data;

	if(wholeSize > 0)
	{
		data = new DATA_TYPE[wholeSize];
		std::copy(A.data, A.data + wholeSize, data);
	}
	else
	{
		data = NULL;
	}

	return *this;
}

template <class DATA_TYPE>
template <class DATA_TYPE2>
Matrix<DATA_TYPE>::operator Matrix<DATA_TYPE2> () const
{
	Matrix<DATA_TYPE2> result(size1, size2);

	unsigned int wholeSize = size1 * size2;

	if(wholeSize == 0)
	{
		return result;
	}
	else
	{
		for(unsigned int i=0; i<wholeSize; ++i)
		{
			DATA_TYPE2 tmp(data[i]);
			result.data[i] = tmp;
		}

		return result;
	}
}

template <class DATA_TYPE>
std::ostream & operator << (std::ostream & os, const Matrix<DATA_TYPE> & A)
{
	for(int i=0, pos=0; i<A.size1; ++i, pos+=A.size2)
	{
		for(int j=0; j<A.size2; ++j)
		{
			os << A.data[pos + j] << '\t';
		}

		os << std::endl;
	}

	return os;
}






template<class DATA_TYPE>
void compute_mat_pow(std::vector<Matrix<DATA_TYPE> > & result, const Matrix<DATA_TYPE> & A, const int order)
{
	unsigned int d = A.rows();

	if(result.size() == 0)
	{
		Matrix<DATA_TYPE> identity(d);

		std::vector<bool> pow_A_computed;

		Matrix<DATA_TYPE> empty;

		for(int i=0; i<=order; ++i)
		{
			pow_A_computed.push_back(false);
			result.push_back(empty);
		}

		pow_A_computed[0] = true;
		result[0] = identity;

		if(order < 1)
		{
			return;
		}

		pow_A_computed[1] = true;
		result[1] = A;

		for(unsigned int i=order; i>1; --i)
		{
			if(!pow_A_computed[i])
			{
				Matrix<DATA_TYPE> temp = A;
				Matrix<DATA_TYPE> temp2 = A;

				unsigned int pos_temp = 1;
				unsigned int pos_result = 1;

				for(unsigned int d=i-1; d > 0;)
				{
					if(d & 1)
					{
						pos_result += pos_temp;

						if(pow_A_computed[pos_result])
						{
							temp2 = result[pos_result];
						}
						else
						{
							temp2 *= temp;
							pow_A_computed[pos_result] = true;
							result[pos_result] = temp2;
						}
					}

					d >>= 1;

					if(d > 0)
					{
						pos_temp <<= 1;

						if(pow_A_computed[pos_temp])
						{
							temp = result[pos_temp];
						}
						else
						{
							temp *= temp;
							pow_A_computed[pos_temp] = true;
							result[pos_temp] = temp;
						}
					}
				}
			}
		}
	}
	else
	{
		unsigned int prev_order = result.size() - 1;

		std::vector<bool> pow_A_computed(order+1);
		for(unsigned int i=0; i<=prev_order; ++i)
		{
			pow_A_computed[i] = true;
		}

		Matrix<DATA_TYPE> empty;
		for(unsigned int i=order; i>prev_order; --i)
		{
			result.push_back(empty);
		}

		for(unsigned int i=order; i>prev_order; --i)
		{
			if(!pow_A_computed[i])
			{
				Matrix<DATA_TYPE> temp = A;
				Matrix<DATA_TYPE> temp2 = A;

				unsigned int pos_temp = 1;
				unsigned int pos_result = 1;

				for(unsigned int d=i-1; d > 0;)
				{
					if(d & 1)
					{
						pos_result += pos_temp;

						if(pow_A_computed[pos_result])
						{
							temp2 = result[pos_result];
						}
						else
						{
							temp2 *= temp;
							pow_A_computed[pos_result] = true;
							result[pos_result] = temp2;
						}
					}

					d >>= 1;

					if(d > 0)
					{
						pos_temp <<= 1;

						if(pow_A_computed[pos_temp])
						{
							temp = result[pos_temp];
						}
						else
						{
							temp *= temp;
							pow_A_computed[pos_temp] = true;
							result[pos_temp] = temp;
						}
					}
				}
			}
		}
	}
}


/*
// matrix for multivariate polynomials
class mpMatrix
{
protected:
	Polynomial *data;
	int size1;
	int size2;

public:
	mpMatrix();
	mpMatrix(const int m, const int n);
	mpMatrix(const int n);
	mpMatrix(const mpMatrix & mpm);
	~mpMatrix();

	int rows() const;
	int cols() const;

	void intEval(mpMatrix & result, const std::vector<Interval> & val_exp_table) const;
	void intEval(iMatrix & result, const std::vector<Interval> & domain) const;

	void output(FILE *fp, const std::vector<std::string> & varNames) const;

	mpMatrix & operator += (const mpMatrix & mpm);

	mpMatrix operator + (const mpMatrix & mpm) const;

	Polynomial * operator [] (const int i);
	mpMatrix & operator = (const mpMatrix & mpm);

	friend class iMatrix;
	friend class upMatrix;
};
*/

class MatrixParseSetting
{
public:
	std::string strExpression;
	Matrix<Real> result;

public:
	MatrixParseSetting();
	MatrixParseSetting(const MatrixParseSetting & setting);
	~MatrixParseSetting();

	MatrixParseSetting & operator = (const MatrixParseSetting & setting);
};


extern MatrixParseSetting matrixParseSetting;

}

//void parse_Matrix();

#endif /* MATRIX_H_ */
