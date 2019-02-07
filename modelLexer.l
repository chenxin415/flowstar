%{
/*---
  Flow*: A Verification Tool for Cyber-Physical Systems.
  Authors: Xin Chen, Sriram Sankaranarayanan, and Erika Abraham.
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.

  The code is released as is under the GNU General Public License (GPL).
---*/

#include "modelParser.h"
#include "modelParser.tab.h"
%}

delim [ \t\r]
line [\n]
whitespace {delim}+
digit [0-9]
letter [a-zA-Z_]
ident {letter}({letter}|{digit}|".")*
number ("-"?)(({digit}+)|({digit}*"."{digit}*)|({digit}+"e"(({digit}+)|("-"{digit}+)))|({digit}*"."{digit}*"e"(({digit}+)|("-"{digit}+)))|({digit}*"."{digit}*"e"("+"({digit}+)|("-"{digit}+))))


%%
"\'" { return '\'';}
"+" {return '+';}
"*" {return '*';}
"-" {return '-';}
"=" {return EQ;}
">=" {return GEQ;}
"<=" {return LEQ;}
"," {return ',';}
";" {return ';';}
":" {return ':';}
"(" {return '(';}
")" {return ')';}
"{" {return '{';}
"}" {return '}';}
"[" {return '[';}
"]" {return ']';}
":=" {return ASSIGN;}
"^" {return '^';}
"/" {return '/';}
">" {return '>';}
"mode" {return MODE;}
"init" {return INIT;}
"in" {return BELONGSTO;}
"end" {return END;}
"setting" {return SETTING;}
"continuous reachability" {return CONTINUOUS;}
"hybrid reachability" {return HYBRID;}
"fixed steps" {return FIXEDST;}
"adaptive steps" {return ADAPTIVEST;}
"fixed orders" {return FIXEDORD;}
"adaptive orders" {return ADAPTIVEORD;}
"order" {return ORDER;}
"min" {return MIN;}
"max" {return MAX;}
"remainder estimation" {return REMEST;}
"interval" {return INTERVAL;}
"octagon" {return OCTAGON;}
"grid" {return GRID;}
"plot" {return PLOT;}
"tm" {return TM;}
"time" {return TIME;}
"modes" {return MODES;}
"jumps" {return JUMPS;}
"inv" {return INV;}
"guard" {return GUARD;}
"reset" {return RESET;}
"start" {return START;}
"max jumps" {return MAXJMPS;}
"output" {return OUTPUT;}
"no output" {return NOOUTPUT;}
"print on" {return PRINTON;}
"print off" {return PRINTOFF;}
"unsafe" {return UNSAFESET;}
"state var" {return STATEVAR;}
"tm var" {return TMVAR;}
"parallelotope aggregation" {return PARAAGGREG;}
"interval aggregation" {return INTAGGREG;}
"taylor model aggregation" {return TMAGGREG;}
"continuous flowpipes" {return CONTINUOUSFLOW;}
"linear continuous flowpipes" {return LINEARCONTINUOUSFLOW;}
"hybrid flowpipes" {return HYBRIDFLOW;}
"taylor picard" {return TAYLOR_PICARD;}
"taylor remainder" {return TAYLOR_REMAINDER;}
"taylor polynomial" {return TAYLOR_POLYNOMIAL;}
"exp" {return EXP;}
"sin" {return SIN;}
"cos" {return COS;}
"log" {return LOG;}
"sqrt" {return SQRT;}
"ode" {return ODE;}
"cutoff" {return CUTOFF;}
"precision" {return PRECISION;}
"gnuplot" {return GNUPLOT;}
"matlab" {return MATLAB;}
"computation paths" {return COMPUTATIONPATHS;}
"par" {return PAR;}
"uncertain" {return UNC;}
"lti ode" {return LTIODE;}
"ltv ode" {return LTVODE;}
"univariate polynomial" {return UNIVARIATE_POLYNOMIAL;}
"multivariate polynomial" {return MULTIVARIATE_POLYNOMIAL;}
"time-inv" {return TIME_INV;}
"time-var" {return TIME_VAR;}
"step" {return STEP;}
"true" {return TRUE;}
"false" {return FALSE;}
"expression ast" {return EXPRESSION;}
"matrix" {return MATRIX;}




{number} { yylval.dblVal = atof( (char *)yytext ); return NUM; }

{ident}	{ yylval.identifier = new std::string(yytext); return IDENT; }

{whitespace} 

{line} { lineNum++; }

"#" {	/* Comment line */
	int c;
	c = yyinput();
	while(c!='\n' && c!=EOF)
	{
		c = yyinput();
	}
	
	if(c == '\n')
	{
		++lineNum;
	}
}

.  { printf("Unknown Character in line %d : %s -- Ignored\n", lineNum, yytext); }

%%


int yywrap()
{
	return 1;
}
/*
void parseUnivariatePolynomial(const std::string & strPolynomial)
{
	int length = strPolynomial.size();
	char *strInput = new char[length+5];

	sprintf(strInput, "%s", strPolynomial.c_str());

	YY_BUFFER_STATE current = YY_CURRENT_BUFFER;

	YY_BUFFER_STATE yybs = yy_scan_string(strInput);

	yy_switch_to_buffer(yybs);
	yyparse();
    yy_delete_buffer(yybs);

	yy_switch_to_buffer(current);

	delete [] strInput;
}
*/

void parseMultivariatePolynomial()
{
	int length = flowstar::multivariate_polynomial_setting.strPolynomial.size();
	char *strPolynomial = new char[length+5];

	sprintf(strPolynomial, "%s", flowstar::multivariate_polynomial_setting.strPolynomial.c_str());

	YY_BUFFER_STATE current = YY_CURRENT_BUFFER;

	YY_BUFFER_STATE yybs = yy_scan_string(strPolynomial);

	yy_switch_to_buffer(yybs);
	yyparse();
    yy_delete_buffer(yybs);

	yy_switch_to_buffer(current);

	delete [] strPolynomial;
}

void parseExpression()
{
	int length = flowstar::expression_ast_setting.strExpression.size();
	char *strExpression = new char[length+5];

	sprintf(strExpression, "%s", flowstar::expression_ast_setting.strExpression.c_str());

	YY_BUFFER_STATE current = YY_CURRENT_BUFFER;

	YY_BUFFER_STATE yybs = yy_scan_string(strExpression);

	yy_switch_to_buffer(yybs);
	yyparse();
    yy_delete_buffer(yybs);

	yy_switch_to_buffer(current);

	delete [] strExpression;
}

