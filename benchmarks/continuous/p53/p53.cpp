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


	ODE<Real> ode({	"1800 - 0.0359*x1*x5 - 0.0693*x1",
					"5.4 + 54*(x1^2/(547600 + x1^2)) - 2.88*x2",
					"2.88*x2 - 0.5198*x3",
					"59.76*x3 - 3.24*x4",
					"3.24*x4 - 5.976e-4*x4^2 - 0.0359*x5*x6",
					"1800 - 0.1155*x6 - 0.0359*x5*x6"}, vars);


	Computational_Setting setting(vars);

	setting.setFixedStepsize(0.01, 5);

	// set the cutoff threshold
	setting.setCutoffThreshold(1e-8);

	// set up the remainder estimation
	Interval I(-1e-1, 1e-1);
	vector<Interval> remainder_estimation(vars.size(), I);
	setting.setRemainderEstimation(remainder_estimation);


	double w = 0.2;	// radius of the initial set

	// define the initial set which is a box
	Interval init_x(20-w, 20+w);

	vector<Interval> initial_box(vars.size());
	initial_box[x1_id] = init_x;
	initial_box[x2_id] = init_x;
	initial_box[x3_id] = init_x;
	initial_box[x4_id] = init_x;
	initial_box[x5_id] = init_x;
	initial_box[x6_id] = init_x;


	Flowpipe initialSet(initial_box);


	vector<Constraint> safeSet;

	Result_of_Reachability result;

	clock_t begin, end;
	begin = clock();

	Symbolic_Remainder symbolic_remainder(initialSet, 500);

	double T = 10;

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

	plot_setting.setOutputDims("x1", "x2");
	plot_setting.plot_2D_octagon_GNUPLOT("./", "p53_x1_x2", result.tmv_flowpipes, setting);

	plot_setting.setOutputDims("x3", "x6");
	plot_setting.plot_2D_octagon_GNUPLOT("./", "p53_x3_x6", result.tmv_flowpipes, setting);


	return 0;
}

