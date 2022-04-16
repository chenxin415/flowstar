/*---
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.
  The code is released as is under the GNU General Public License (GPL).
---*/

#include "modelParser.h"

extern int yyparse();

int main(int argc, const char *argv[])
{
	yyparse();

	return 0;
}




