#include "../../../flowstar-toolbox/Continuous.h"

using namespace flowstar;
using namespace std;


int main()
{
	Variables vars;

	int x1_id = vars.declareVar("x1");
	int y1_id = vars.declareVar("y1");
	int x2_id = vars.declareVar("x2");
	int y2_id = vars.declareVar("y2");
	int t_id = vars.declareVar("t");


	ODE<Real> ode({"y1", "2*(1 - x1^2) * y1 - 2*x1 + x2", "y2", "2*(1 - x2^2) * y2 - 2*x2 + x1", "1"}, vars);


	Computational_Setting setting(vars);

	setting.setAdaptiveStepsize(0.001, 0.03, 5);

	setting.setCutoffThreshold(2e-8);

	Interval remainder(-5e-6,5e-6);
	vector<Interval> estimation(vars.size(), remainder);
	setting.setRemainderEstimation(estimation);

	setting.printOff();	// do not display information during the reachability analysis

	// define the initial set which is a box
	Interval init_x(1.55, 1.85), init_y(2.35, 2.45);


	// subdividing the initial set to smaller boxes
	// we only subdivide the set to 8 pieces in the dimensions of x1 and x2
	list<Interval> subdiv_x1;
	init_x.split(subdiv_x1, 8);

	list<Interval> subdiv_x2 = subdiv_x1;

	list<Interval>::iterator iter1 = subdiv_x1.begin();

	vector<Flowpipe> initial_sets;

	for(; iter1 != subdiv_x1.end(); ++iter1)
	{
		list<Interval>::iterator iter2 = subdiv_x2.begin();

		for(; iter2 != subdiv_x2.end(); ++iter2)
		{
			vector<Interval> box(vars.size());
			box[x1_id] = *iter1;
			box[y1_id] = init_y;
			box[x2_id] = *iter2;
			box[y2_id] = init_y;

			Flowpipe initialSet(box);

			initial_sets.push_back(initialSet);
		}
	}


	// safety: y1 <= 4 and y2 <= 4
	vector<Constraint> safeSet = {Constraint("y1 - 4", vars), Constraint("y2 - 4", vars)};


	Result_of_Reachability result;	// data structure to keep the reachability analysis result


	clock_t begin, end;
	begin = clock();

	double T = 8;	// time horizon

	int safety = SAFE;

	for(int i=0; i<initial_sets.size(); ++i)
	{
		Symbolic_Remainder sr(initial_sets[i], 1000);
		ode.reach(result, initial_sets[i], T, setting, safeSet, sr);

		if(!result.isCompleted()) // if the flowpipes are successfully computed
		{
			safety = -1;
			break;
		}

		if(result.isUnsafe()) // if there is an unsafe flowpipe detected, the computation terminates immediately
		{
			safety = UNSAFE;
			break;
		}
		else if(!result.isSafe())// there is no unsafe flowpipe found, but some of the flowpipes intersect the unsafe set
		{
			safety = UNKNOWN;
		}
	}


	end = clock();
	printf("time cost: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);


	if(safety == -1)
	{
		printf("Flowpipe computation is terminated due to the large overestimation.\n");
	}
	else if(safety == SAFE) // if all of the flowpipes are safe
	{
		printf("All flowpipes are safe.\n");
	}
	else if(safety == UNSAFE) // if there is an unsafe flowpipe detected, the computation terminates immediately
	{
		printf("The last flowpipe is unsafe.\n");
	}
	else // there is no unsafe flowpipe found, but some of the flowpipes intersect the unsafe set
	{
		printf("The safety is unknown.\n");
	}


	return 0;
}

