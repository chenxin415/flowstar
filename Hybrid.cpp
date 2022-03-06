/*---
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.
  The code is released as is under the GNU General Public License (GPL).
---*/

#include "Hybrid.h"

using namespace flowstar;


Initial_Set_Configuration::Initial_Set_Configuration()
{
	remaining_time = 0;
	previous_mode = 0;
}

Initial_Set_Configuration::Initial_Set_Configuration(const Flowpipe & fp, const double t, const int mode)
{
	initialSet = fp;
	remaining_time = t;
	previous_mode = mode;
}

Initial_Set_Configuration::~Initial_Set_Configuration()
{
}

Initial_Set_Configuration & Initial_Set_Configuration::operator = (const Initial_Set_Configuration & isc)
{
	if(this == &isc)
		return *this;

	initialSet = isc.initialSet;
	remaining_time = isc.remaining_time;
	previous_mode = isc.previous_mode;

	return *this;
}










Queue_of_Initial_Set::Queue_of_Initial_Set()
{
}

Queue_of_Initial_Set::~Queue_of_Initial_Set()
{
	iscs.clear();
}

bool Queue_of_Initial_Set::isEmpty() const
{
	return iscs.size() == 0 ? true : false;
}

void Queue_of_Initial_Set::enqueue(const Initial_Set_Configuration & isc)
{
	iscs.push_back(isc);
}

void Queue_of_Initial_Set::dequeue(Initial_Set_Configuration & isc)
{
	isc = iscs.front();
	iscs.pop_front();
}

