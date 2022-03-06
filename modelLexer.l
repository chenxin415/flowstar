%{
/*---
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
"^" {return '^';}
"/" {return '/';}
">" {return '>';}
"exp" {return EXP;}
"sin" {return SIN;}
"cos" {return COS;}
"log" {return LOG;}
"sqrt" {return SQRT;}
"univariate polynomial" {return UNIVARIATE_POLYNOMIAL;}
"multivariate polynomial" {return MULTIVARIATE_POLYNOMIAL;}
"expression ast" {return EXPRESSION;}




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
	int length = flowstar::expression_setting.strExpression.size();
	char *strExpression = new char[length+5];

	sprintf(strExpression, "%s", flowstar::expression_setting.strExpression.c_str());

	YY_BUFFER_STATE current = YY_CURRENT_BUFFER;

	YY_BUFFER_STATE yybs = yy_scan_string(strExpression);

	yy_switch_to_buffer(yybs);
	yyparse();
    yy_delete_buffer(yybs);

	yy_switch_to_buffer(current);

	delete [] strExpression;
}

