%{
	/*---
	Flow*: A Verification Tool for Cyber-Physical Systems.
	Authors: Xin Chen, Sriram Sankaranarayanan, and Erika Abraham.
	Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.

	The code is released as is under the GNU General Public License (GPL).
	---*/


	#include "modelParser.h"

	extern int yyerror(const char *);
	extern int yyerror(std::string);
	extern int yylex();
	extern int yyparse();
	bool err;

	int lineNum = 1;

	std::vector<flowstar::Flowpipe> initialSets;

	void parseError(const char *str, int lnum)
	{
		std::cerr << "Error @line " << lineNum << ":" << std::string(str) << std::endl;
		exit(1);
	}
%}

%union
{
	double															dblVal;
	int																intVal;
	std::string														*identifier;
	std::vector<flowstar::Interval>									*intVec;
	std::vector<std::vector<flowstar::Interval> >					*intsVec;
	flowstar::Polynomial<flowstar::Interval>						*intPoly;
	flowstar::Polynomial<flowstar::Real>							*realPoly;
	flowstar::UnivariatePolynomial<flowstar::Real>					*uniPoly;
	std::vector<flowstar::Polynomial<flowstar::Interval> >			*intPolyVec;
	std::vector<flowstar::UnivariatePolynomial<flowstar::Real> >	*uniPolyVec;
	std::vector<flowstar::Flowpipe>									*flowpipeVec;
	std::vector<flowstar::Constraint>								*constraintVec;
	flowstar::TaylorModelVec<flowstar::Real>						*taylorModelVec;
	flowstar::Expression_AST<flowstar::Interval>					*pIntExpression;
	std::vector<flowstar::Expression_AST<flowstar::Interval> >		*intExpressionVec;
	LTI_Term														*pLTITerm;
	LTI_ODE_description												*pLTIDescription;
	LTV_Term														*pLTVTerm;
	LTV_ODE_description												*pLTVDescription;
}


%token <dblVal> NUM
%token <identifier> IDENT
%token STATEVAR TMVAR TM EQ GEQ LEQ ASSIGN END
%token MODE INIT BELONGSTO
%token PARAAGGREG INTAGGREG TMAGGREG
%token OUTPUT NOOUTPUT
%token CONTINUOUS HYBRID
%token SETTING
%token FIXEDST FIXEDORD ADAPTIVEST ADAPTIVEORD ORDER
%token MIN MAX
%token REMEST
%token INTERVAL OCTAGON GRID PLOT
%token QRPRECOND IDPRECOND
%token TIME
%token MODES JUMPS INV GUARD RESET START MAXJMPS
%token PRINTON PRINTOFF UNSAFESET
%token CONTINUOUSFLOW HYBRIDFLOW
%token TAYLOR_PICARD TAYLOR_REMAINDER TAYLOR_POLYNOMIAL
%token EXP SIN COS LOG SQRT
%token ODE CUTOFF PRECISION
%token GNUPLOT MATLAB COMPUTATIONPATHS
%token LTIODE LTVODE PAR UNC
%token UNIVARIATE_POLYNOMIAL MULTIVARIATE_POLYNOMIAL
%token TIME_INV TIME_VAR STEP
%token TRUE FALSE
%token LINEARCONTINUOUSFLOW
%token EXPRESSION
%token MATRIX



%type <intPoly>						multivariate_polynomial
%type <uniPoly>						univariate_polynomial
%type <intVec>						remainders
%type <intVec>						interval_vector
%type <intsVec>						set_of_interval_vectors
%type <flowpipeVec>					init
%type <constraintVec>				constraints
%type <taylorModelVec>				vector_of_Taylor_models
%type <intVec>						domain_setting
%type <realPoly>					Taylor_model_polynomial
%type <pIntExpression>				expression
%type <intExpressionVec>			ode
%type <intVec>						lti_polynomial
%type <pLTITerm>					lti_term
%type <pLTIDescription>				lti_ode
%type <uniPolyVec>					ltv_polynomial
%type <pLTVTerm>					ltv_term
%type <pLTVDescription>				ltv_ode


%left GEQ LEQ EQ 
%left '+' '-'
%left '*' '/'
%nonassoc uminus
%right '^'

%start model

%%

model: CONTINUOUS '{' continuous '}' unsafe_setting
{
	int mkres = mkdir(outputDir.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
	if(mkres < 0 && errno != EEXIST)
	{
		printf("Can not create the directory for output files.\n");
		exit(1);
	}

	mkres = mkdir(imageDir.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
	if(mkres < 0 && errno != EEXIST)
	{
		printf("Can not create the directory for images.\n");
		exit(1);
	}

	tmVars.declareVar("local_t");
	for(unsigned int i=0; i<stateVars.varNames.size(); ++i)
	{
		tmVars.declareVar(stateVars.varNames[i] + "_0");
	}

	int result;
	Continuous_Reachability reachability(problem_description);

	clock_t begin, end;
	begin = clock();
	result = reachability.run();
	end = clock();

	printf("\n");

	if(reachability.bSafetyChecking)
	{
		switch(result)
		{
		case COMPLETED_UNSAFE:
			printf("Computation completed: %ld flowpipe(s) computed.\n\n", reachability.result_of_reachability.num_of_flowpipes);
			printf(BOLD_FONT "Result of the safety verification on the computed flowpipes: " RESET_COLOR);
			printf(BOLD_FONT RED_COLOR "UNSAFE\n\n" RESET_COLOR);
			break;
		case COMPLETED_SAFE:
			printf("Computation completed: %ld flowpipe(s) computed.\n\n", reachability.result_of_reachability.num_of_flowpipes);
			printf(BOLD_FONT "Result of the safety verification on the computed flowpipes: " RESET_COLOR);
			printf(BOLD_FONT GREEN_COLOR "SAFE\n\n" RESET_COLOR);
			break;
		case COMPLETED_UNKNOWN:
			printf("Computation completed: %ld flowpipe(s) computed.\n\n", reachability.result_of_reachability.num_of_flowpipes);
			printf(BOLD_FONT "Result of the safety verification on the computed flowpipes: " RESET_COLOR);
			printf(BOLD_FONT BLUE_COLOR "UNKNOWN\n\n" RESET_COLOR);
			break;
		case UNCOMPLETED_UNSAFE:
			printf(BOLD_FONT RED_COLOR "Computation not completed: %ld flowpipe(s) computed.\n" RESET_COLOR, reachability.result_of_reachability.num_of_flowpipes);
			printf(BOLD_FONT "Please try smaller step sizes or larger Taylor model orders.\n\n" RESET_COLOR);
			printf(BOLD_FONT "Result of the safety verification on the computed flowpipes: " RESET_COLOR);
			printf(BOLD_FONT RED_COLOR "UNSAFE\n\n" RESET_COLOR);
			break;
		case UNCOMPLETED_SAFE:
			printf(BOLD_FONT RED_COLOR "Computation not completed: %ld flowpipe(s) computed.\n" RESET_COLOR, reachability.result_of_reachability.num_of_flowpipes);
			printf(BOLD_FONT "Please try smaller step sizes or larger Taylor model orders.\n\n" RESET_COLOR);
			printf(BOLD_FONT "Result of the safety verification on the computed flowpipes: " RESET_COLOR);
			printf(BOLD_FONT GREEN_COLOR "SAFE\n\n" RESET_COLOR);
			break;
		case UNCOMPLETED_UNKNOWN:
			printf(BOLD_FONT RED_COLOR "Computation not completed: %ld flowpipe(s) computed.\n" RESET_COLOR, reachability.result_of_reachability.num_of_flowpipes);
			printf(BOLD_FONT "Please try smaller step sizes or larger Taylor model orders.\n\n" RESET_COLOR);
			printf(BOLD_FONT "Result of the safety verification on the computed flowpipes: " RESET_COLOR);
			printf(BOLD_FONT BLUE_COLOR "UNKNOWN\n\n" RESET_COLOR);
			break;
		}

		printf("Total time cost:" BOLD_FONT " %lf" RESET_COLOR " seconds.\n\n", (double)(end - begin) / CLOCKS_PER_SEC);
	}
	else
	{
		if(result == COMPLETED_SAFE)
		{
			printf("Computation completed: %ld flowpipe(s) computed.\n", reachability.result_of_reachability.num_of_flowpipes);
		}
		else
		{
			printf(BOLD_FONT RED_COLOR "Computation not completed: %ld flowpipe(s) computed.\n" RESET_COLOR, reachability.result_of_reachability.num_of_flowpipes);
			printf(BOLD_FONT "Please try smaller step sizes or larger Taylor model orders.\n" RESET_COLOR);
		}

		printf("Total time cost:" BOLD_FONT " %lf" RESET_COLOR " seconds.\n\n", (double)(end - begin) / CLOCKS_PER_SEC);
	}

	if(problem_description.bPlot)
	{
		if(problem_description.bTMOutput)
		{
			reachability.prepareForTMOutput();

			reachability.p_p_setting->bProjected = false;
			reachability.plot_2D();

			std::ofstream os(outputDir + problem_description.fileName + ".flow", std::ofstream::out);

			reachability.tmOutput(os);

			os.close();
		}
		else
		{
			reachability.prepareForPlotting();
		
			reachability.p_p_setting->bProjected = true;
			reachability.plot_2D();
		}
	}
	else
	{
		if(problem_description.bTMOutput)
		{
			reachability.prepareForTMOutput();
			
			std::ofstream os(outputDir + problem_description.fileName + ".flow", std::ofstream::out);
			reachability.tmOutput(os);

			os.close();
		}
	}
}
|
stateVarDecls plotting ORDER NUM CUTOFF NUM output_env unsafe_setting CONTINUOUSFLOW '{' tmVarDecls continuous_flowpipes '}'
{
	if($4 <= 0)
	{
		parseError("The cutoff threshold should be a positive number.", lineNum);
		exit(1);
	}

	problem_description.order_min = (int)$4;

	flowstar::Interval I(-$6,$6);
	problem_description.setCutoff(I);

	reachability_for_outputFile.setup(problem_description);

	clock_t begin, end;
	int checkingResult;
	printf("Safety checking ...\n");
	begin = clock();

	if(reachability_for_outputFile.bSafetyChecking)
	{
		checkingResult = reachability_for_outputFile.safetyChecking();
	}

	end = clock();
	printf("Done.\n");
	printf("Time cost of safety verification:" BOLD_FONT " %lf" RESET_COLOR " seconds.\n\n", (double)(end - begin) / CLOCKS_PER_SEC);
	printf(BOLD_FONT "Result of safety verification: " RESET_COLOR);

	switch(checkingResult)
	{
	case UNSAFE:
		printf(BOLD_FONT RED_COLOR "UNSAFE\n\n" RESET_COLOR);
		break;
	case SAFE:
		printf(BOLD_FONT GREEN_COLOR "SAFE\n\n" RESET_COLOR);
		break;
	case UNKNOWN:
		printf(BOLD_FONT BLUE_COLOR "UNKNOWN\n\n" RESET_COLOR);
		break;
	}

	if(problem_description.bPlot)
	{
		reachability_for_outputFile.p_p_setting->bProjected = false;
		reachability_for_outputFile.plot_2D();
	}
}
|
MULTIVARIATE_POLYNOMIAL '{' multivariate_polynomial '}'
{
	flowstar::multivariate_polynomial_setting.result = *$3;
	delete $3;
}
|
EXPRESSION '{' expression '}'
{
	flowstar::expression_ast_setting.result = *$3;
	delete $3;
}
;


multivariate_polynomial: multivariate_polynomial '+' multivariate_polynomial
{
	$$ = $1;
	*$$ += *$3;

	delete $3;
}
|
multivariate_polynomial '-' multivariate_polynomial
{
	$$ = $1;
	*$$ -= *$3;

	delete $3;
}
|
'(' multivariate_polynomial ')'
{
	$$ = $2; 
}
|
multivariate_polynomial '*' multivariate_polynomial
{
	$$ = $1;
	*$$ *= *$3;

	delete $3;
}
|
multivariate_polynomial '^' NUM
{
	int exp = (int) $3;

	if(exp == 0)
	{
		$$ = new flowstar::Polynomial<flowstar::Interval>(1, stateVars.size() + 1);
	}
	else
	{
		$$ = new flowstar::Polynomial<flowstar::Interval>;
		(*$1).pow(*$$, exp);
	}

	delete $1;
}
|
'-' multivariate_polynomial %prec uminus
{
	$$ = $2;
	*$$ *= -1;
}
|
IDENT
{
	int id = stateVars.getIDForVar(*$1);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Variable %s is not declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	unsigned int numVars = stateVars.size() + 1;
	$$ = new flowstar::Polynomial<flowstar::Interval>(1, numVars);
	$$->mul_assign(id + 1, 1);

	delete $1;
}
|
'[' NUM ',' NUM ']'
{
	unsigned int numVars = stateVars.size() + 1;
	flowstar::Interval I($2, $4);
	$$ = new flowstar::Polynomial<flowstar::Interval>(I, numVars);
	flowstar::multivariate_polynomial_setting.bDeterministic = false;
	problem_description.bDeterministic = false;
}
|
NUM
{
	unsigned int numVars = stateVars.size() + 1;
	$$ = new flowstar::Polynomial<flowstar::Interval>($1, numVars);
}
;


univariate_polynomial: univariate_polynomial '+' univariate_polynomial
{
	$$ = $1;
	(*$$) += (*$3);

	delete $3;
}
|
univariate_polynomial '-' univariate_polynomial
{
	$$ = $1;
	(*$$) -= (*$3);

	delete $3;
}
|
'(' univariate_polynomial ')'
{
	$$ = $2; 
}
|
univariate_polynomial '*' univariate_polynomial
{
	$$ = $1;
	(*$$) *= (*$3);

	delete $3;
}
|
univariate_polynomial '^' NUM
{
	int exp = (int) $3;

	if(exp == 0)
	{
		$$ = new flowstar::UnivariatePolynomial<flowstar::Real>(1);
	}
	else
	{
		$$ = new flowstar::UnivariatePolynomial<flowstar::Real>;
		(*$1).pow(*$$, exp);
	}

	delete $1;
}
|
'-' univariate_polynomial %prec uminus
{
	$$ = $2;
	(*$$) *= -1;
}
|
IDENT
{
	std::string tVar("t");
	if($1->compare(tVar) == 0)
	{
		$$ = new flowstar::UnivariatePolynomial<flowstar::Real>(1, 1);
	}
	else
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "The time variable should be denoted by t.");
		parseError(errMsg, lineNum);
		exit(1);
	}
}
|
NUM
{
	$$ = new flowstar::UnivariatePolynomial<flowstar::Real>($1);
}
;


expression: expression '+' expression
{
	$$ = $1;
	(*$$) += (*$3);

	delete $3;
}
|
expression '-' expression
{
	$$ = $1;
	(*$$) -= (*$3);

	delete $3;
}
|
'(' expression ')'
{
	$$ = $2; 
}
|
expression '*' expression
{
	$$ = $1;
	(*$$) *= (*$3);

	delete $3;
}
|
expression '^' NUM
{
	$$ = $1;
	$$->pow_assign((int)$3);
}
|
'-' expression %prec uminus
{
	$$ = $2;
	$$->inv_assign();
}
|
IDENT
{
	$$ = new flowstar::Expression_AST<flowstar::Interval>(*$1, stateVars);

	delete $1;
}
|
'[' NUM ',' NUM ']'
{
	flowstar::Interval I($2, $4);
	$$ = new flowstar::Expression_AST<flowstar::Interval>(I);
	problem_description.bDeterministic = false;
}
|
NUM
{
	$$ = new flowstar::Expression_AST<flowstar::Interval>($1);
}
|
expression '/' expression
{
	$$ = $1;
	(*$$) /= (*$3);

	delete $3;
}
|
EXP '(' expression ')'
{
	$$ = $3;
	$$->exp_assign();
}
|
SIN '(' expression ')'
{
	$$ = $3;
	$$->sin_assign();
}
|
COS '(' expression ')'
{
	$$ = $3;
	$$->cos_assign();
}
|
LOG '(' expression ')'
{
	$$ = $3;
	$$->log_assign();
}
|
SQRT '(' expression ')'
{
	$$ = $3;
	$$->sqrt_assign();
}
;




continuous: stateVarDecls SETTING '{' setting '}' ODE '{' ode '}' INIT '{' init '}'
{
	if(problem_description.bDeterministic)
	{
		problem_description.type_of_dynamics = DETERMINISTIC_DYN;
		problem_description.deterministic_dynamics.clear();

		for(unsigned int i=0; i<$8->size(); ++i)
		{
			flowstar::Expression_AST<Real> tmp;
			(*$8)[i].toReal(tmp);
			problem_description.deterministic_dynamics.push_back(tmp);
		}
	}
	else
	{
		problem_description.type_of_dynamics = NONDETERMINISTIC_DYN;
		problem_description.nondeterministic_dynamics = *$8;
	}

	problem_description.initialSets = *$12;

	delete $8;
	delete $12;
}
|
stateVarDecls SETTING '{' setting '}' ODE '{' NUM '}' '{' ode '}' INIT '{' init '}'
{
	if(problem_description.bDeterministic)
	{
		problem_description.type_of_dynamics = DETERMINISTIC_DYN;
		problem_description.deterministic_dynamics.clear();

		for(unsigned int i=0; i<$11->size(); ++i)
		{
			flowstar::Expression_AST<Real> tmp;
			(*$11)[i].toReal(tmp);
			problem_description.deterministic_dynamics.push_back(tmp);
		}
	}
	else
	{
		problem_description.type_of_dynamics = NONDETERMINISTIC_DYN;
		problem_description.nondeterministic_dynamics = *$11;
	}

	problem_description.bSymbolicRemainder = true;
	problem_description.queue_size = (int)$8;
	problem_description.initialSets = *$15;

	delete $11;
	delete $15;
}
|
stateVarDecls SETTING '{' setting '}' LTIODE '{' lti_ode '}' INIT '{' init '}'
{
	problem_description.type_of_dynamics = LINEAR_TIME_INVARIANT;
	problem_description.rm_dyn_A = $8->dyn_A;
	problem_description.utm_dyn_B = $8->dyn_B;

	problem_description.initialSets = *$12;

	delete $8;
	delete $12;
}
|
stateVarDecls SETTING '{' setting '}' LTVODE '{' ltv_ode '}' INIT '{' init '}'
{
	problem_description.type_of_dynamics = LINEAR_TIME_VARYING;
	problem_description.upm_dyn_A = $8->dyn_A;
	problem_description.upm_dyn_B = $8->dyn_B;

	problem_description.initialSets = *$12;

	delete $8;
	delete $12;
}
;
|
stateVarDecls SETTING '{' setting '}' LTVODE '{' tv_par_list '}' '{' NUM '}' '{' ltv_ode '}' INIT '{' init '}'
{
	problem_description.type_of_dynamics = LINEAR_TIME_VARYING;
	problem_description.upm_dyn_A = $14->dyn_A;
	problem_description.upm_dyn_B = $14->dyn_B;
	problem_description.upm_dyn_tv = $14->dyn_tv;

	problem_description.queue_size = (unsigned int)$11;

	problem_description.initialSets = *$18;

	delete $14;
	delete $18;
}
;




stateVarDecls: STATEVAR stateIdDeclList
{
}
;

stateIdDeclList: stateIdDeclList ',' IDENT
{
	if(stateVars.declareVar(*$3) < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s has already been declared.", (*$3).c_str());
		parseError(errMsg, lineNum);
		exit(0);
	}

	delete $3;
}
|
IDENT
{
	if(stateVars.declareVar(*$1) < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s has already been declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(0);
	}

	delete $1;
}
;


setting: FIXEDST NUM TIME NUM remainder_estimation plotting FIXEDORD NUM CUTOFF NUM PRECISION NUM output_env printing
{
	bool bValid = problem_description.setFixedStepsize($2);

	if(!bValid)
		exit(0);

	bValid = problem_description.setFixedOrder((unsigned int)$8);

	if(!bValid)
		exit(0);

	bValid = problem_description.setTimeHorizon($4);

	if(!bValid)
		exit(0);

	if($10 <= 0)
	{
		parseError("The cutoff threshold should be positive", lineNum);
		exit(1);
	}

	Interval cutoff_threshold(-$10,$10);
	problem_description.setCutoff(cutoff_threshold);

	bValid = problem_description.setPrecision($12);

	if(!bValid)
		exit(0);
}
|
FIXEDST NUM TIME NUM remainder_estimation plotting ADAPTIVEORD '{' MIN NUM ',' MAX NUM '}' CUTOFF NUM PRECISION NUM output_env printing
{
	bool bValid = problem_description.setFixedStepsize($2);

	if(!bValid)
		exit(0);

	bValid = problem_description.setAdaptiveOrder((unsigned int)$10, (unsigned int)$13);

	if(!bValid)
		exit(0);

	bValid = problem_description.setTimeHorizon($4);

	if(!bValid)
		exit(0);

	if($16 <= 0)
	{
		parseError("The cutoff threshold should be positive", lineNum);
		exit(1);
	}

	Interval cutoff_threshold(-$16,$16);
	problem_description.setCutoff(cutoff_threshold);

	bValid = problem_description.setPrecision($18);

	if(!bValid)
		exit(0);
}
|
ADAPTIVEST '{' MIN NUM ',' MAX NUM '}' TIME NUM remainder_estimation plotting FIXEDORD NUM CUTOFF NUM PRECISION NUM output_env printing
{
	bool bValid = problem_description.setAdaptiveStepsize($4, $7);

	if(!bValid)
		exit(0);

	bValid = problem_description.setFixedOrder((unsigned int)$14);

	if(!bValid)
		exit(0);

	bValid = problem_description.setTimeHorizon($10);

	if(!bValid)
		exit(0);

	if($16 <= 0)
	{
		parseError("The cutoff threshold should be positive", lineNum);
		exit(1);
	}

	Interval cutoff_threshold(-$16,$16);
	problem_description.setCutoff(cutoff_threshold);

	bValid = problem_description.setPrecision($18);

	if(!bValid)
		exit(0);
}
;


remainder_estimation: REMEST NUM
{
	if($2 <= 0)
	{
		parseError("Remainder estimation should be a positive number.", lineNum);
		exit(1);
	}

	Interval I(-$2, $2);
	problem_description.remainder_estimation.clear();

	for(unsigned int i=0; i<stateVars.size(); ++i)
	{
		problem_description.remainder_estimation.push_back(I);
	}
}
|
REMEST '{' remainders '}'
{
	for(int i=0; i<$3->size(); ++i)
	{
		if((*$3)[i].inf() >= (*$3)[i].sup() - THRESHOLD_LOW)
		{
			parseError("Invalid remainder estimation.", lineNum);
			exit(0);
		}
	}

	problem_description.setRemainderEstimation(*$3);
	delete $3;
}
;

remainders: remainders ',' IDENT ':' '[' NUM ',' NUM ']'
{
	$$ = $1;
	int id = stateVars.getIDForVar(*$3);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$3).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if($6 >= $8)
	{
		parseError("Invalid remainder estimation.", lineNum);
		exit(1);
	}

	flowstar::Interval I($6,$8);
	(*$$)[id] = I;
	delete $3;
}
|
IDENT ':' '[' NUM ',' NUM ']'
{
	int numVars = stateVars.size();
	$$ = new std::vector<flowstar::Interval>(numVars);

	int id = stateVars.getIDForVar(*$1);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if($4 >= $6)
	{
		parseError("Invalid remainder estimation.", lineNum);
		exit(0);
	}

	flowstar::Interval I($4,$6);
	(*$$)[id] = I;
	delete $1;
}
;


plotting: GNUPLOT INTERVAL IDENT ',' IDENT
{
	int x = stateVars.getIDForVar(*$3);
	int y = stateVars.getIDForVar(*$5);

	if(x < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$3).c_str());
		parseError(errMsg, lineNum);
		exit(0);
	}

	if(y < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$5).c_str());
		parseError(errMsg, lineNum);
		exit(0);
	}

	problem_description.setOutputDims(x, y);
	problem_description.setFileType(PLOT_GNUPLOT);
	problem_description.setObjectType(PLOT_INTERVAL);

	delete $3;
	delete $5;
}
|
GNUPLOT OCTAGON IDENT ',' IDENT
{
	int x = stateVars.getIDForVar(*$3);
	int y = stateVars.getIDForVar(*$5);

	if(x < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$3).c_str());
		parseError(errMsg, lineNum);
		exit(0);
	}

	if(y < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$5).c_str());
		parseError(errMsg, lineNum);
		exit(0);
	}

	problem_description.setOutputDims(x, y);
	problem_description.setFileType(PLOT_GNUPLOT);
	problem_description.setObjectType(PLOT_OCTAGON);

	delete $3;
	delete $5;
}
|
GNUPLOT GRID NUM IDENT ',' IDENT
{
	int x = stateVars.getIDForVar(*$4);
	int y = stateVars.getIDForVar(*$6);

	if(x < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$4).c_str());
		parseError(errMsg, lineNum);
		exit(0);
	}

	if(y < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$6).c_str());
		parseError(errMsg, lineNum);
		exit(0);
	}

	problem_description.setOutputDims(x, y);
	problem_description.setFileType(PLOT_GNUPLOT);
	problem_description.setObjectType(PLOT_GRID);
	problem_description.setNumOfPieces((unsigned int)$3);

	delete $4;
	delete $6;
}
|
MATLAB INTERVAL IDENT ',' IDENT
{
	int x = stateVars.getIDForVar(*$3);
	int y = stateVars.getIDForVar(*$5);

	if(x < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$3).c_str());
		parseError(errMsg, lineNum);
		exit(0);
	}

	if(y < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$5).c_str());
		parseError(errMsg, lineNum);
		exit(0);
	}

	problem_description.setOutputDims(x, y);
	problem_description.setFileType(PLOT_MATLAB);
	problem_description.setObjectType(PLOT_INTERVAL);

	delete $3;
	delete $5;
}
|
MATLAB OCTAGON IDENT ',' IDENT
{
	int x = stateVars.getIDForVar(*$3);
	int y = stateVars.getIDForVar(*$5);

	if(x < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$3).c_str());
		parseError(errMsg, lineNum);
		exit(0);
	}

	if(y < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$5).c_str());
		parseError(errMsg, lineNum);
		exit(0);
	}

	problem_description.setOutputDims(x, y);
	problem_description.setFileType(PLOT_MATLAB);
	problem_description.setObjectType(PLOT_OCTAGON);

	delete $3;
	delete $5;
}
|
MATLAB GRID NUM IDENT ',' IDENT
{
	int x = stateVars.getIDForVar(*$4);
	int y = stateVars.getIDForVar(*$6);

	if(x < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$4).c_str());
		parseError(errMsg, lineNum);
		exit(0);
	}

	if(y < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$6).c_str());
		parseError(errMsg, lineNum);
		exit(0);
	}

	problem_description.setOutputDims(x, y);
	problem_description.setFileType(PLOT_MATLAB);
	problem_description.setObjectType(PLOT_GRID);
	problem_description.setNumOfPieces((unsigned int)$3);

	delete $4;
	delete $6;
}
;


output_env: OUTPUT IDENT
{
	problem_description.setFileName(*$2);
	problem_description.plotOn();
	problem_description.tmOutputOn();
}
|
PLOT OUTPUT IDENT
{
	problem_description.setFileName(*$3);
	problem_description.plotOn();
	problem_description.tmOutputOff();
}
|
TM OUTPUT IDENT
{
	problem_description.setFileName(*$3);
	problem_description.plotOff();
	problem_description.tmOutputOn();
}
|
NOOUTPUT
{
	problem_description.plotOff();
	problem_description.tmOutputOff();
}
;


printing: PRINTON
{
	problem_description.printOn();
}
|
PRINTOFF
{
	problem_description.printOff();
}
;


ode: ode IDENT '\'' EQ expression
{
	$$ = $1;

	int id = stateVars.getIDForVar(*$2);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$2).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	(*$$)[id] = *$5;

	delete $2;
	delete $5;
}
|
{
	unsigned int numVars = stateVars.size();
	$$ = new std::vector<flowstar::Expression_AST<flowstar::Interval> >(numVars);
}
;


init: set_of_interval_vectors
{
	$$ = new std::vector<flowstar::Flowpipe>();
	for(unsigned int i=0; i<$1->size(); ++i)
	{
		flowstar::Flowpipe initialSet((*$1)[i]);
		$$->push_back(initialSet);
	}
}
|
interval_vector
{
	$$ = new std::vector<flowstar::Flowpipe>(1);
	flowstar::Flowpipe initialSet(*$1);
	(*$$)[0] = initialSet;
	delete $1;
}
;


interval_vector: interval_vector IDENT BELONGSTO '[' NUM ',' NUM ']'
{
	int id = stateVars.getIDForVar(*$2);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$2).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if($5 > $7)
	{
		parseError("Invalid interval.", lineNum);
		exit(1);
	}

	flowstar::Interval I($5,$7);
	$$ = $1;
	(*$$)[id] = I;

	delete $2;
}
|
{
	int numVars = stateVars.size();
	$$ = new std::vector<flowstar::Interval>(numVars);
}
;


set_of_interval_vectors: set_of_interval_vectors '{' interval_vector '}'
{
	$$->push_back(*$3);

	delete $3;
}
|
'{' interval_vector '}'
{
	$$ = new std::vector<std::vector<flowstar::Interval> >(0);
	$$->push_back(*$2);

	delete $2;
}
;

unsafe_setting: UNSAFESET '{' constraints '}'
{
	problem_description.safetyCheckingOn();
	problem_description.setUnsafe(*$3);
	delete $3;
}
|
{
	problem_description.safetyCheckingOff();
}
;

constraints: constraints expression LEQ NUM
{
	$$ = $1;
	flowstar::Expression_AST<flowstar::Real> realExp;
	$2->toReal(realExp);
	flowstar::Constraint constraint(realExp, $4);
	$$->push_back(constraint);

	delete $2;
}
|
{
	$$ = new std::vector<Constraint>();
}
;

tmVarDecls: TMVAR tmIdDeclList
{
}
;

tmIdDeclList: tmIdDeclList ',' IDENT
{
	if(tmVars.declareVar(*$3) < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "TM variable %s has already been declared.", (*$3).c_str());
		parseError(errMsg, lineNum);
		exit(0);
	}

	delete $3;
}
|
IDENT
{
	if(tmVars.declareVar(*$1) < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "TM variable %s has already been declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(0);
	}

	delete $1;
}
;


continuous_flowpipes: continuous_flowpipes '{' vector_of_Taylor_models domain_setting '}'
{
	reachability_for_outputFile.result_of_reachability.tmv_flowpipes.push_back(*$3);

	Flowpipe flowpipe;
	flowpipe.domain = *$4;
	reachability_for_outputFile.result_of_reachability.nonlinear_flowpipes.push_back(flowpipe);

	delete $3;
	delete $4;
}
|
'{' vector_of_Taylor_models domain_setting '}'
{
	reachability_for_outputFile.result_of_reachability.tmv_flowpipes.push_back(*$2);

	Flowpipe flowpipe;
	flowpipe.domain = *$3;
	reachability_for_outputFile.result_of_reachability.nonlinear_flowpipes.push_back(flowpipe);

	delete $2;
	delete $3;
}
;


vector_of_Taylor_models: vector_of_Taylor_models IDENT EQ Taylor_model_polynomial '+' '[' NUM ',' NUM ']'
{
	int id = stateVars.getIDForVar(*$2);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$2).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if($7 > $9)
	{
		parseError("Invalid interval remainder.", lineNum);
		exit(1);
	}

	flowstar::Interval I($7,$9);
	flowstar::TaylorModel<Real> tm(*$4, I);
	$$ = $1;
	$$->tms[id] = tm;

	delete $2;
	delete $4;
}
|
{
	std::vector<TaylorModel<Real> > tmv(stateVars.size());
	$$ = new TaylorModelVec<Real>();
	$$->tms = tmv;
}
;

domain_setting: domain_setting IDENT BELONGSTO '[' NUM ',' NUM ']'
{
	int id = tmVars.getIDForVar(*$2);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Taylor model variable %s is not declared.", (*$2).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if($5 > $7)
	{
		parseError("Invalid interval.", lineNum);
		exit(1);
	}

	flowstar::Interval I($5,$7);
	$$ = $1;
	(*$$)[id] = I;

	delete $2;
}
|
IDENT BELONGSTO '[' NUM ',' NUM ']'
{
	int id = tmVars.getIDForVar(*$1);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Taylor model variable %s is not declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	if($4 > $6)
	{
		parseError("Invalid interval.", lineNum);
		exit(1);
	}

	flowstar::Interval I($4,$6);
	$$ = new std::vector<Interval>(tmVars.size());
	(*$$)[id] = I;

	delete $1;
}
;


Taylor_model_polynomial: Taylor_model_polynomial '+' Taylor_model_polynomial
{
	$$ = $1;
	*$$ += *$3;

	delete $3;
}
|
Taylor_model_polynomial '-' Taylor_model_polynomial
{
	$$ = $1;
	*$$ -= *$3;

	delete $3;
}
|
'(' Taylor_model_polynomial ')'
{
	$$ = $2; 
}
|
Taylor_model_polynomial '*' Taylor_model_polynomial
{
	$$ = $1;
	*$$ *= *$3;

	delete $3;
}
|
Taylor_model_polynomial '^' NUM
{
	int exp = (int) $3;

	if(exp == 0)
	{
		$$ = new flowstar::Polynomial<flowstar::Real>(1, tmVars.size());
	}
	else
	{
		$$ = new flowstar::Polynomial<flowstar::Real>();
		(*$1).pow(*$$, exp);
	}

	delete $1;
}
|
'-' Taylor_model_polynomial %prec uminus
{
	$$ = $2;
	*$$ *= -1;
}
|
IDENT
{
	int id = tmVars.getIDForVar(*$1);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Taylor model variable %s is not declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	$$ = new flowstar::Polynomial<flowstar::Real>(1, tmVars.size());
	$$->mul_assign(id, 1);

	delete $1;
}
|
NUM
{
	$$ = new flowstar::Polynomial<flowstar::Real>($1, tmVars.size());
}
;

lti_ode: lti_ode IDENT '\'' EQ lti_polynomial
{
	int id = stateVars.getIDForVar(*$2);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$2).c_str());
		parseError(errMsg, lineNum);
		exit(0);
	}

	unsigned int n = $5->size();
	for(unsigned int i=0; i<n-1; ++i)
	{
		$$->dyn_A[id][i] = (*$5)[i].toReal();
	}

	UnivariateTaylorModel<Real> tmp((*$5).back());
	$$->dyn_B[id][0] = tmp;

	delete $2;
	delete $5;
}
|
{
	unsigned int numVars = stateVars.size();
	$$ = new LTI_ODE_description(numVars);
}
;

lti_polynomial: lti_polynomial '+' lti_term
{
	if($3->varID >= 0)
	{
		(*$$)[$3->varID] += $3->coefficient;
	}
	else
	{
		(*$$).back() += $3->coefficient;
	}

	delete $3;
}
|
lti_polynomial '-' lti_term
{
	if($3->varID >= 0)
	{
		(*$$)[$3->varID] -= $3->coefficient;
	}
	else
	{
		(*$$).back() -= $3->coefficient;
	}

	delete $3;
}
|
'-' lti_term
{
	$$ = new std::vector<Interval>(stateVars.size() + 1);

	if($2->varID >= 0)
	{
		(*$$)[$2->varID] -= $2->coefficient;
	}
	else
	{
		(*$$).back() -= $2->coefficient;
	}

	delete $2;
}
|
lti_term
{
	$$ = new std::vector<Interval>(stateVars.size() + 1);

	if($1->varID >= 0)
	{
		(*$$)[$1->varID] += $1->coefficient;
	}
	else
	{
		(*$$).back() += $1->coefficient;
	}

	delete $1;
}
;

lti_term: NUM '*' IDENT
{
	int id = stateVars.getIDForVar(*$3);

	if(id >= 0)
	{
		$$ = new LTI_Term($1, id, -1);
	}
	else
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Symbol %s is not defined.", (*$3).c_str());
		parseError(errMsg, lineNum);
		exit(0);
	}

	delete $3;
}
|
IDENT
{
	int id = stateVars.getIDForVar(*$1);

	if(id >= 0)
	{
		$$ = new LTI_Term(1, id, -1);
	}
	else
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Symbol %s is not defined.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(0);
	}

	delete $1;
}
|
'[' NUM ',' NUM ']'
{
	Interval I($2, $4);
	$$ = new LTI_Term(I, -1, -1);
}
|
NUM
{
	$$ = new LTI_Term($1, -1, -1);
}
;




tv_par_list: tv_par_list ',' IDENT
{
	if(tvPars.declareVar(*$3) < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Symbol %s has already been declared.", (*$3).c_str());
		parseError(errMsg, lineNum);
		exit(0);
	}

	delete $3;
}
|
IDENT
{
	if(tvPars.declareVar(*$1) < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Symbol %s has already been declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(0);
	}

	delete $1;
}
;



ltv_ode: ltv_ode IDENT '\'' EQ ltv_polynomial
{
	int id = stateVars.getIDForVar(*$2);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$2).c_str());
		parseError(errMsg, lineNum);
		exit(0);
	}

	unsigned int numVars = stateVars.size();
	unsigned int numTVPars = tvPars.size();

	for(unsigned int i=0; i<numVars; ++i)
	{
		$$->dyn_A[id][i] = (*$5)[i];
	}

	$$->dyn_B[id][0] = (*$5)[numVars];

	for(unsigned int i=0; i<numTVPars; ++i)
	{
		$$->dyn_tv[id][i] = (*$5)[i + numTVPars + 1];
	}

	delete $2;
	delete $5;
}
|
IDENT '\'' EQ ltv_polynomial
{
	unsigned int numVars = stateVars.size();
	unsigned int numTVPars = tvPars.size();

	$$ = new LTV_ODE_description(numVars, numTVPars);

	int id = stateVars.getIDForVar(*$1);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "State variable %s is not declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(0);
	}

	for(unsigned int i=0; i<numVars; ++i)
	{
		$$->dyn_A[id][i] = (*$4)[i];
	}

	$$->dyn_B[id][0] = (*$4)[numVars];

	for(unsigned int i=0; i<numTVPars; ++i)
	{
		$$->dyn_tv[id][i] = (*$4)[i + numVars + 1];
	}

	delete $1;
	delete $4;
}
;

ltv_polynomial: ltv_polynomial '+' ltv_term
{
	if($3->varID >= 0)
	{
		(*$$)[$3->varID] += $3->coefficient;
	}
	else if($3->tvParID >= 0)
	{
		(*$$)[$3->tvParID + stateVars.size() + 1] += $3->coefficient;
	}
	else
	{
		(*$$)[stateVars.size()] += $3->coefficient;
	}

	delete $3;
}
|
ltv_polynomial '-' ltv_term
{
	if($3->varID >= 0)
	{
		(*$$)[$3->varID] -= $3->coefficient;
	}
	else if($3->tvParID >= 0)
	{
		(*$$)[$3->tvParID + stateVars.size() + 1] -= $3->coefficient;
	}
	else
	{
		(*$$)[stateVars.size()] -= $3->coefficient;
	}

	delete $3;
}
|
'-' ltv_term
{
	$$ = new std::vector<UnivariatePolynomial<Real> >(stateVars.size() + 1 + tvPars.size());

	if($2->varID >= 0)
	{
		(*$$)[$2->varID] -= $2->coefficient;
	}
	else if($2->tvParID >= 0)
	{
		(*$$)[$2->tvParID + stateVars.size() + 1] -= $2->coefficient;
	}
	else
	{
		(*$$)[stateVars.size()] -= $2->coefficient;
	}

	delete $2;
}
|
ltv_term
{
	$$ = new std::vector<UnivariatePolynomial<Real> >(stateVars.size() + 1 + tvPars.size());

	if($1->varID >= 0)
	{
		(*$$)[$1->varID] += $1->coefficient;
	}
	else if($1->tvParID >= 0)
	{
		(*$$)[$1->tvParID + stateVars.size() + 1] += $1->coefficient;
	}
	else
	{
		(*$$)[stateVars.size()] += $1->coefficient;
	}

	delete $1;
}
;

ltv_term: '(' univariate_polynomial ')' '*' IDENT
{
	int id = stateVars.getIDForVar(*$5);

	if(id >= 0)
	{
		$$ = new LTV_Term(*$2, id, -1, -1);
	}
	else
	{
		id = tvPars.getIDForVar(*$5);

		if(id >= 0)
		{
			$$ = new LTV_Term(*$2, -1, -1, id);
		}
		else
		{
			char errMsg[MSG_SIZE];
			sprintf(errMsg, "Symbol %s is not defined.", (*$5).c_str());
			parseError(errMsg, lineNum);
			exit(0);
		}
	}

	delete $2;
	delete $5;
}
|
IDENT
{
	int id = stateVars.getIDForVar(*$1);

	if(id >= 0)
	{
		$$ = new LTV_Term(1, id, -1, -1);
	}
	else
	{
		id = tvPars.getIDForVar(*$1);

		if(id >= 0)
		{
			$$ = new LTV_Term(1, -1, -1, id);
		}
		else
		{
			char errMsg[MSG_SIZE];
			sprintf(errMsg, "Symbol %s is not defined.", (*$1).c_str());
			parseError(errMsg, lineNum);
			exit(0);
		}
	}

	delete $1;
}
|
'(' univariate_polynomial ')'
{
	$$ = new LTV_Term(*$2, -1, -1, -1);
	delete $2;
}
;










%%

int yyerror(const char * what)
{
	fprintf(stderr, "Error line %d: %s\n", lineNum, what);
	err = true;
	return 1;
}

int yyerror(std::string what)
{
	std::cerr << "Error line "<< lineNum << " " << what << std::endl;
	err = true;
	return 1;
}
