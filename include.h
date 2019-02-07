/*---
  Flow*: A Verification Tool for Cyber-Physical Systems.
  Authors: Xin Chen, Sriram Sankaranarayanan, and Erika Abraham.
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.
  
  The code is released as is under the GNU General Public License (GPL).
---*/

#ifndef INCLUDE_H_
#define INCLUDE_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <cmath>
#include <mpfr.h>
#include <vector>
#include <string>
#include <list>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cassert>
#include <memory>
#include <map>
#include <time.h>
#include <algorithm>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <glpk.h>


const int normal_precision	=	53;
const int high_precision	=	256;

#define MAX_REFINEMENT_STEPS	49
#define MAX_DYN_REF				49

#define MAX_WIDTH				1e-12

#define THRESHOLD_HIGH			1e-12
#define THRESHOLD_LOW			1e-20

#define STOP_RATIO				0.99
#define ABS_STOP_RATIO			0.99

#define PN 						15 			// the number of digits printed
#define INVALID 				-1e10
#define UNBOUNDED				1e10

#define DC_THRESHOLD_SEARCH 	1e-6
#define DC_THRESHOLD_IMPROV		0.9

#define UNIFORM			0
#define MULTI			1

#define LAMBDA_DOWN		0.5
#define LAMBDA_UP		1.1

#define NAME_SIZE		100
#define PRINT_SIZE		100

#define UNSAFE			-1
#define SAFE			0
#define UNKNOWN			1

#define LINEAR_TIME_INVARIANT		0
#define LINEAR_TIME_VARYING			1
#define DETERMINISTIC_DYN			2
#define NONDETERMINISTIC_DYN		3


#define COMPLETED_UNSAFE		1
#define COMPLETED_SAFE			2
#define COMPLETED_UNKNOWN		3
#define UNCOMPLETED_SAFE		4
#define UNCOMPLETED_UNSAFE		5
#define UNCOMPLETED_UNKNOWN		6


#define PLOT_GNUPLOT	0
#define PLOT_MATLAB		1

#define PLOT_INTERVAL	0
#define PLOT_OCTAGON	1
#define PLOT_GRID		2

#define INTERVAL_AGGREG	0
#define PARA_AGGREG		1
#define PCA_AGGREG		2

#define MSG_SIZE		100
#define NUM_LENGTH		50

#define REFINEMENT_PREC	1e-5

#define APPROX_TOLERANCE 1e-12

#define RESET_COLOR		"\033[0m"
#define BLACK_COLOR		"\033[30m"
#define RED_COLOR		"\033[31m"
#define GREEN_COLOR		"\033[32m"
#define BLUE_COLOR		"\033[34m"
#define BOLD_FONT		"\e[1m"



const char str_pi_up[]	=	"3.14159265358979323846264338327950288419716939937511";
const char str_pi_lo[]	=	"3.14159265358979323846264338327950288419716939937510";

const std::string outputDir = "./outputs/";
const std::string imageDir = "./images/";
const std::string counterexampleDir = "./counterexamples/";
const char local_var_name[] = "local_var_";

const char str_prefix_taylor_picard[] = "taylor picard { ";
const char str_prefix_taylor_remainder[] = "taylor remainder { ";
const char str_prefix_taylor_polynomial[] = "taylor polynomial { ";
const char str_prefix_center[] = "nonpolynomial center { ";

const char str_prefix_combination_picard[] = "combination picard { ";
const char str_prefix_combination_remainder[] = "combination remainder { ";
const char str_prefix_combination_polynomial[] = "combination polynomial { ";

const char str_prefix_univariate_polynomial[] = "univariate polynomial { ";
const char str_prefix_multivariate_polynomial[] = "multivariate polynomial { ";
const char str_prefix_expression_ast[] = "expression ast { ";

const char str_prefix_expression[] = "expression { ";

const char str_prefix_matrix[] = "matrix { ";

const char str_suffix[] = " }";

const std::string str_counterexample_dumping_name_suffix = ".counterexample";


extern int lineNum;

extern void parseODE();

#endif /* INCLUDE_H_ */
