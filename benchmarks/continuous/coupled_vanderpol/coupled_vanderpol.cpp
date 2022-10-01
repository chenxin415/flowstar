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

	setting.setCutoffThreshold(1e-8);

	setting.setAdaptiveStepsize(0.0001, 0.02, 4);

	Interval remainder(-3e-5,3e-5);
	vector<Interval> estimation(vars.size(), remainder);
	setting.setRemainderEstimation(estimation);

	setting.printOff();	// do not display information during the reachability analysis

	// safety: y1 <= 4 and y2 <= 4
	vector<Constraint> safeSet = {Constraint("y1 - 4", vars), Constraint("y2 - 4", vars)};

	// define the initial set which is a box
	Interval init_x(1.55, 1.85), init_y(2.35, 2.45);

	vector<Interval> global_box(vars.size());
	global_box[x1_id] = init_x;
	global_box[y1_id] = init_y;
	global_box[x2_id] = init_x;
	global_box[y2_id] = init_y;



	Result_of_Reachability result;	// data structure to keep the reachability analysis result


	// subdividing the initial set to smaller boxes 16 x 1 x 16 x 1
	list<Interval> subdiv_x0;
	global_box[0].split(subdiv_x0, 16);

	list<Interval> subdiv_x1;
	global_box[1].split(subdiv_x1, 1);

	list<Interval> subdiv_x2;
	global_box[2].split(subdiv_x2, 16);

	list<Interval> subdiv_x3;
	global_box[3].split(subdiv_x3, 1);

	list<Interval>::iterator iter0 = subdiv_x0.begin();

	vector<Flowpipe> initial_sets;
	vector<vector<Interval> > initial_boxes;

	for(; iter0 != subdiv_x0.end(); ++iter0)
	{
		list<Interval>::iterator iter1 = subdiv_x1.begin();

		for(; iter1 != subdiv_x1.end(); ++iter1)
		{
			list<Interval>::iterator iter2 = subdiv_x2.begin();

			for(; iter2 != subdiv_x2.end(); ++iter2)
			{
				list<Interval>::iterator iter3 = subdiv_x3.begin();

				for(; iter3 != subdiv_x3.end(); ++iter3)
				{
					vector<Interval> box(vars.size());
					box[0] = *iter0;
					box[1] = *iter1;
					box[2] = *iter2;
					box[3] = *iter3;
					box[4] = 0;

					Flowpipe initialSet(box);

					initial_sets.push_back(initialSet);
					initial_boxes.push_back(box);
				}
			}
		}
	}

	clock_t begin, end;
	begin = clock();

	int safety = SAFE;
	double T = 8;

	for(int i=0; i<initial_sets.size(); ++i)
	{
		Flowpipe local_init = initial_sets[i];

		Symbolic_Remainder sr(local_init, 400);

		ode.reach(result, local_init, T, setting, safeSet, sr);

		if(!result.isCompleted()) // if the flowpipes are successfully computed
		{
			safety = -1;

			cout << "Failed: " << initial_boxes[i][0] << "\t" << initial_boxes[i][1] << "\t" << initial_boxes[i][2] << "\t" << initial_boxes[i][3] << endl;

			break;
		}

		if(result.isUnsafe()) // if there is an unsafe flowpipe detected, the computation terminates immediately
		{
			safety = UNSAFE;
			break;
		}
		else if(!result.isSafe())// there is no unsafe flowpipe found, but some of the flowpipes intersect the unsafe set
		{
			cout << "Unknown: " << initial_boxes[i][0] << "\t" << initial_boxes[i][1] << "\t" << initial_boxes[i][2] << "\t" << initial_boxes[i][3] << endl;

			safety = UNKNOWN;
		}

		cout << "Partition " << i+1 << " completed." << endl;
	}


	end = clock();
	printf("time cost: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);

	return 0;
}

