%{
/*---
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
	flowstar::UnivariatePolynomial<flowstar::Real>					*uniPoly;
	flowstar::Polynomial<flowstar::Interval>						*intPoly;
	flowstar::Expression<flowstar::Interval>						*pIntExpression;
}


%token <dblVal> NUM
%token <identifier> IDENT
%token EXP SIN COS LOG SQRT
%token UNIVARIATE_POLYNOMIAL MULTIVARIATE_POLYNOMIAL
%token EXPRESSION



%type <intPoly>						multivariate_polynomial
%type <uniPoly>						univariate_polynomial
%type <pIntExpression>				expression



%left GEQ LEQ EQ 
%left '+' '-'
%left '*' '/'
%nonassoc uminus
%right '^'

%start input_content

%%

input_content: MULTIVARIATE_POLYNOMIAL '{' multivariate_polynomial '}'
{
	flowstar::multivariate_polynomial_setting.result = *$3;
	delete $3;
}
|
UNIVARIATE_POLYNOMIAL '{' univariate_polynomial '}'
{
	flowstar::up_parseresult = *$3;
	delete $3;
}
|
EXPRESSION '{' expression '}'
{
	flowstar::expression_setting.result = *$3;
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
		$$ = new flowstar::Polynomial<flowstar::Interval>(1, flowstar::multivariate_polynomial_setting.pVars->size());
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
	int id = flowstar::multivariate_polynomial_setting.pVars->getIDForVar(*$1);

	if(id < 0)
	{
		char errMsg[MSG_SIZE];
		sprintf(errMsg, "Variable %s is not declared.", (*$1).c_str());
		parseError(errMsg, lineNum);
		exit(1);
	}

	unsigned int numVars = flowstar::multivariate_polynomial_setting.pVars->size();
	$$ = new flowstar::Polynomial<flowstar::Interval>(1, numVars);
	$$->mul_assign(id, 1);

	delete $1;
}
|
'[' NUM ',' NUM ']'
{
	unsigned int numVars = flowstar::multivariate_polynomial_setting.pVars->size();
	flowstar::Interval I($2, $4);
	$$ = new flowstar::Polynomial<flowstar::Interval>(I, numVars);
}
|
NUM
{
	unsigned int numVars = flowstar::multivariate_polynomial_setting.pVars->size();
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
	if(flowstar::expression_setting.pVars == NULL)
	{
		$$ = new flowstar::Expression<flowstar::Interval>();
		$$->root = std::shared_ptr<AST_Node<flowstar::Interval> > (new AST_Node<flowstar::Interval>(VAR_ID, -1));
	}
	else
	{
		$$ = new flowstar::Expression<flowstar::Interval>(*$1, flowstar::expression_setting.pVars);
	}

	delete $1;
}
|
'[' NUM ',' NUM ']'
{
	flowstar::Interval I($2, $4);
	$$ = new flowstar::Expression<flowstar::Interval>(I);
}
|
NUM
{
	$$ = new flowstar::Expression<flowstar::Interval>((Interval)$1);
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
