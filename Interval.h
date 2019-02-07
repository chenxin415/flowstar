/*---
  Flow*: A Verification Tool for Cyber-Physical Systems.
  Authors: Xin Chen, Sriram Sankaranarayanan, and Erika Abraham.
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.
  
  The code is released as is under the GNU General Public License (GPL).
---*/

#ifndef INTERVAL_H_
#define INTERVAL_H_

#include "include.h"

extern mpfr_prec_t intervalNumPrecision;

namespace flowstar
{

class Interval;

class Real
{
protected:
	mpfr_t value;
public:
	Real();
	Real(const double d);
	Real(const Real & real);
	~Real();

	void set(const double c);
	bool isZero() const;
	bool belongsTo(const Interval & I) const;

	double getValue_RNDD() const;
	double getValue_RNDU() const;
	void mag(Real & r) const;
	void abs(Real & r) const;
	double abs() const;
	double mag() const;
	void abs_assign();

	void to_sym_int(Interval & I) const;		// to a symmetric interval

	void exp_RNDU(Real & result) const;
	void exp_assign_RNDU();

	void pow_assign_RNDU(const int n);
	void pow_assign(const int n);

	void factorial(const unsigned int n);

	void rec(Real & result) const;
	void rec_assign();

	void add_RNDD(Real & result, const Real & real) const;
	void add_assign_RNDD(const Real & real);
	void add_RNDU(Real & result, const Real & real) const;
	void add_assign_RNDU(const Real & real);
	void add_RNDN(Real & result, const Real & real) const;
	void add_assign_RNDN(const Real & real);

	void sub_RNDD(Real & result, const Real & real) const;
	void sub_assign_RNDD(const Real & real);
	void sub_RNDU(Real & result, const Real & real) const;
	void sub_assign_RNDU(const Real & real);

	void mul_RNDD(Real & result, const Real & real) const;
	void mul_assign_RNDD(const Real & real);
	void mul_RNDU(Real & result, const Real & real) const;
	void mul_assign_RNDU(const Real & real);

	void mul_RNDD(Real & result, const int n) const;
	void mul_assign_RNDD(const int n);
	void mul_RNDU(Real & result, const int n) const;
	void mul_assign_RNDU(const int n);

	void div_RNDD(Real & result, const Real & real) const;
	void div_assign_RNDD(const Real & real);
	void div_RNDU(Real & result, const Real & real) const;
	void div_assign_RNDU(const Real & real);

	void div_RNDD(Real & result, const int n) const;
	void div_assign_RNDD(const int n);
	void div_RNDU(Real & result, const int n) const;
	void div_assign_RNDU(const int n);

	void output(FILE *fp) const;
	void dump(FILE *fp) const;

	void sin_assign();
	void cos_assign();
	void exp_assign();
	void log_assign();
	void sqrt_assign();

	void sin(Real & c) const;
	void cos(Real & c) const;
	void exp(Real & c) const;
	void log(Real & c) const;
	void sqrt(Real & c) const;

	double toDouble() const;
	std::string toString() const;

	Interval operator * (const Interval & I) const;

	Real & operator += (const Real & r);
	Real & operator -= (const Real & r);
	Real & operator *= (const Real & r);
	Real & operator /= (const Real & r);

	Real & operator += (const double d);
	Real & operator -= (const double d);
	Real & operator *= (const double d);
	Real & operator /= (const double d);

	Real operator - () const;
	Real operator + (const Real & r) const;
	Real operator - (const Real & r) const;
	Real operator * (const Real & r) const;
	Real operator / (const Real & r) const;


	bool operator == (const Real & r) const;
	bool operator != (const Real & r) const;
	bool operator >= (const Real & r) const;
	bool operator <= (const Real & r) const;
	bool operator > (const Real & r) const;
	bool operator < (const Real & r) const;

	bool operator == (const double c) const;
	bool operator != (const double c) const;
	bool operator >= (const double c) const;
	bool operator <= (const double c) const;
	bool operator > (const double c) const;
	bool operator < (const double c) const;


	Real & operator = (const Real & r);		// should always be the same precision
	Real & operator = (const double c);

	friend std::ostream & operator << (std::ostream & output, const Real & r);


	friend Real operator + (const double d, const Real & r);
	friend Real operator - (const double d, const Real & r);
	friend Real operator * (const double d, const Real & r);
	friend Real operator / (const double d, const Real & r);

	friend class Interval;
};





class Interval
{
protected:
	mpfr_t lo;		// the lower bound
	mpfr_t up;		// the upper bound

public:
	Interval();
	Interval(const double c);
	Interval(const Real & r);
	Interval(const double l, const double u);
	Interval(const Real & c, const Real & r);
	Interval(const Real & l, const Real & u, const int n);	// n is reserved
	Interval(const char *strLo, const char *strUp);
	Interval(const Interval & I);
	~Interval();

	bool isZero() const;
	bool isSingle() const;

	void set(const double l, const double u);
	void set(const double c);
	void set(const Real & r);

	void setInf(const double l);
	void setInf(const Interval & I);

	void setSup(const double u);
	void setSup(const Interval & S);

	void split(Interval & left, Interval & right) const;			// split the interval at the midpoint
	void split(std::list<Interval> & result, const int n) const;	// split the interval uniformly by n parts

	void set_inf();

	double sup() const;
	double inf() const;

	void sup(Interval & S) const;
	void inf(Interval & I) const;

	void sup(Real & u) const;
	void inf(Real & l) const;

	double midpoint() const;
	void midpoint(Interval & M) const;
	void midpoint(Real & mid) const;

	void toCenterForm(Real & center, Real & radius) const;

	void remove_midpoint(Interval & M);
	void remove_midpoint(Real & c);
	double remove_midpoint();

	Interval intersect(const Interval & I) const;
	void intersect_assign(const Interval & I);

	void bloat(const double e);		// e >= 0
	void bloat(const Real & e);		// e >= 0
	bool within(const Interval & I, const double e) const;

	double width() const;
	void width(Interval & W) const;

	double mag() const;		// max{|lo|,|up|}
	void mag(Real & m) const;
	void mag(Interval & M) const;

	void abs(Interval & result) const;
	void abs_assign();		// absolute value

	bool subseteq(const Interval & I) const;	// returns true if the interval is a subset of I
	bool supseteq(const Interval & I) const;	// returns true if the interval is a superset of I
	bool valid() const;

	bool lessThan(const Interval & I) const;
	bool greaterThan(const Interval & I) const;
	bool lessThan(const Real & r) const;
	bool greaterThan(const Real & r) const;
	bool lessThan(const double r) const;
	bool greaterThan(const double r) const;

	bool lessThanEq(const Interval & I) const;
	bool lessThanEq(const Real & r) const;

	bool operator == (const Interval & I) const;
	bool operator != (const Interval & I) const;

	double toDouble() const;
	Real toReal() const;
	std::string toString() const;

	Interval & operator = (const Interval & I);
	Interval & operator = (const Real & r);
	Interval & operator = (const double d);

	Interval & operator += (const Interval & I);
	Interval & operator += (const Real & r);
	Interval & operator += (const double c);

	Interval & operator -= (const Interval & I);
	Interval & operator -= (const Real & r);
	Interval & operator -= (const double c);

	Interval & operator *= (const Interval & I);
	Interval & operator *= (const Real & r);
	Interval & operator *= (const double c);

	Interval & operator /= (const Interval & I);
	Interval & operator /= (const Real & r);
	Interval & operator /= (const double c);

	Interval & operator ++ ();
	Interval & operator -- ();

	Interval operator - () const;
	const Interval operator + (const Interval & I) const;
	const Interval operator + (const Real & r) const;
	const Interval operator + (const double c) const;

	const Interval operator - (const Interval & I) const;
	const Interval operator - (const Real & r) const;
	const Interval operator - (const double c) const;

	const Interval operator * (const Interval & I) const;
	const Interval operator * (const Real & r) const;
	const Interval operator * (const double c) const;

	const Interval operator / (const Interval & I) const;
	const Interval operator / (const double c) const;

	void sqrt(Interval & result) const;		// square root
	void inv(Interval & result) const;		// additive inverse
	void rec(Interval & result) const;		// reciprocal
	void sqrt_assign();
	void inv_assign();
	void rec_assign();

	void add_assign(const double c);
	void sub_assign(const double c);
	void mul_assign(const double c);
	void div_assign(const double c);

	void mul_add(Interval *result, const Interval *intVec, const int size);

	Interval pow(const int n) const;
	Interval exp() const;
	Interval sin() const;
	Interval cos() const;
	Interval log() const;

	void pow_assign(const int n);
	void exp_assign();
	void sin_assign();
	void cos_assign();
	void log_assign();

	double widthRatio(const Interval & I) const;

	void hull_assign(const Interval & I);

	void dump(FILE *fp) const;
	void output(FILE *fp, const char * msg, const char * msg2) const;
	void output_midpoint(FILE *fp, const int n) const;

	void round(Interval & remainder);

	void shrink_up(const double d);
	void shrink_lo(const double d);

	friend std::ostream & operator << (std::ostream & output, const Interval & I);

	friend class Real;
};


std::ostream & operator << (std::ostream & output, const Real & r);
std::ostream & operator << (std::ostream & output, const Interval & I);



Real operator + (const double d, const Real & r);
Real operator - (const double d, const Real & r);
Real operator * (const double d, const Real & r);
Real operator / (const double d, const Real & r);

}

#endif /* INTERVAL_H_ */
