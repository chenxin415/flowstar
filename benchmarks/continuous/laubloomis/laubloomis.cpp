#include "../../../flowstar-toolbox/Continuous.h"

using namespace flowstar;
using namespace std;


int main()
{
	Variables vars;

	int x1_id = vars.declareVar("x1");
	int x2_id = vars.declareVar("x2");
	int x3_id = vars.declareVar("x3");
	int x4_id = vars.declareVar("x4");
	int x5_id = vars.declareVar("x5");
	int x6_id = vars.declareVar("x6");
	int x7_id = vars.declareVar("x7");
	int t_id = vars.declareVar("t");


	ODE<Real> ode({	"1.4 * x3 - 0.9 * x1",
					"2.5 * x5 - 1.5 * x2",
					"0.6 * x7 - 0.8 * x3 * x2",
					"2 - 1.3 * x4 * x3",
					"0.7 * x1 - x4 * x5",
					"0.3 * x1 - 3.1 * x6",
					"1.8 * x6 - 1.5 * x7 * x2",
					"1"}, vars);


	Computational_Setting setting(vars);

	// adaptive stepsize from 0.001 to 0.2
	// the TM order is fixed at 5
	setting.setAdaptiveStepsize(0.001, 0.2, 5);

	// set the cutoff threshold
	setting.setCutoffThreshold(1e-7);

	// set up the remainder estimation
	Interval I(-1e-4, 1e-4);
	vector<Interval> remainder_estimation(vars.size(), I);
	setting.setRemainderEstimation(remainder_estimation);


	double w = 0.2;	// radius of the initial set

	// define the initial set which is a box
	Interval init_x1(1.2-w, 1.2+w), init_x2(1.05-w, 1.05+w), init_x3(1.5-w, 1.5+w), init_x4(2.4-w, 2.4+w),
			init_x5(1-w, 1+w), init_x6(0.1-w, 0.1+w), init_x7(0.45-w, 0.45+w), init_t;

	vector<Interval> initial_box(vars.size());
	initial_box[x1_id] = init_x1;
	initial_box[x2_id] = init_x2;
	initial_box[x3_id] = init_x3;
	initial_box[x4_id] = init_x4;
	initial_box[x5_id] = init_x5;
	initial_box[x6_id] = init_x6;
	initial_box[x7_id] = init_x7;
	initial_box[t_id] = init_t;


	Flowpipe initialSet(initial_box);


	vector<Constraint> safeSet;

	Result_of_Reachability result;

	clock_t begin, end;
	begin = clock();

	Symbolic_Remainder symbolic_remainder(initialSet, 500);

	double T = 20;

	ode.reach(result, initialSet, T, setting, safeSet, symbolic_remainder);


	end = clock();
	printf("time cost: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);



	if(!result.isCompleted())
	{
		printf("Flowpipe computation is terminated due to the large overestimation.\n");
	}

	if(result.isSafe())
	{
		printf("All flowpipes are safe.\n");
	}
	else if(result.isUnsafe())
	{
		printf("The last flowpipe is unsafe.\n");
	}
	else
	{
		printf("The safety is unknown.\n");
	}


	result.transformToTaylorModels(setting);

	Plot_Setting plot_setting(vars);
	plot_setting.printOn();
	plot_setting.setOutputDims("t", "x4");

	plot_setting.plot_2D_interval_GNUPLOT("./", "laubloomis", result.tmv_flowpipes);

	return 0;
}
