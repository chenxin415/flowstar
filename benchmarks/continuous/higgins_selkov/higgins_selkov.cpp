#include "../../../flowstar-toolbox/Continuous.h"

using namespace flowstar;
using namespace std;


int main()
{
	Variables vars;

	int S_id = vars.declareVar("S");
	int P_id = vars.declareVar("P");

	ODE<Interval> ode({"1 - [0.9995,1.0005]*S*P^2", "[0.9995,1.0005]*S*P^2 - [0.99951, 1.00051]*P"}, vars);


	Computational_Setting setting(vars);

	setting.setFixedStepsize(0.05, 5);

	// set up the remainder estimation
	Interval I(-1e-2, 1e-2);
	vector<Interval> remainder_estimation(vars.size(), I);
	setting.setRemainderEstimation(remainder_estimation);


	Interval init_S(1.99,2.01), init_P(0.99,1.01);

	vector<Interval> initial_box(vars.size());
	initial_box[S_id] = init_S;
	initial_box[P_id] = init_P;


	Flowpipe initialSet(initial_box);


	// no safety specification
	vector<Constraint> safeSet;


	Result_of_Reachability result;

	clock_t begin, end;
	begin = clock();

	Symbolic_Remainder sr(initialSet, 200);

	double T = 10;

	ode.reach(result, initialSet, T, setting, safeSet, sr);

	end = clock();
	printf("time cost: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);


	if(!result.isCompleted())
	{
		printf("Flowpipe computation is terminated due to the large overestimation.\n");
	}


	result.transformToTaylorModels(setting);


	Plot_Setting plot_setting(vars);
	plot_setting.printOn();
	plot_setting.setOutputDims("S", "P");

	plot_setting.plot_2D_octagon_GNUPLOT("./", "higgins_selkov", result.tmv_flowpipes, setting);

	return 0;
}

