/*---
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.
  The code is released as is under the GNU General Public License (GPL).
---*/

#ifndef MODELPARSER_H_
#define MODELPARSER_H_

#include "Continuous.h"

using namespace flowstar;

extern int lineNum;


extern int yyparse();

void parseError(const char *str, int lnum);





#endif /* MODELPARSER_H_ */
