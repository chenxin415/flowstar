/*---
  Email: Xin Chen <chenxin415@gmail.com> if you have questions or comments.
  The code is released as is under the GNU General Public License (GPL).
---*/


#ifndef HYBRID_H_
#define HYBRID_H_

#include "Continuous.h"


class Initial_Set_Configuration
{
public:
	Flowpipe initialSet;
	double remaining_time;
	int previous_mode;

public:
	Initial_Set_Configuration();
	Initial_Set_Configuration(const Flowpipe & fp, const double t, const int mode);
	~Initial_Set_Configuration();

	Initial_Set_Configuration & operator = (const Initial_Set_Configuration & isc);
};




class Queue_of_Initial_Set
{
protected:
	std::list<Initial_Set_Configuration> iscs;

public:
	Queue_of_Initial_Set();
	~Queue_of_Initial_Set();

	bool isEmpty() const;
	void enqueue(const Initial_Set_Configuration & isc);
	void dequeue(Initial_Set_Configuration & isc);
};



#endif /* HYBRID_H_ */
