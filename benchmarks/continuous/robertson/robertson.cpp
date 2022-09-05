#include "../../../flowstar-toolbox/Continuous.h"

using namespace flowstar;
using namespace std;


int main()
{
	Variables vars;

	int x_id = vars.declareVar("x");
	int y_id = vars.declareVar("y");
	int z_id = vars.declareVar("z");
	int t_id = vars.declareVar("t");


	ODE<Real> ode({	"-0.4*x + 1000*y*z",
			"0.4*x - 1000*y*z - 1e7*y^2",
			"1e7*y^2",
			"1"}, vars);

	// set the reachability parameters
	Computational_Setting setting(vars);


	setting.setAdaptiveStepsize(0.00001, 0.1, 3);
	setting.setCutoffThreshold(1e-15);
	Interval I(-1e-8, 1e-8);
	vector<Interval> remainder_estimation(vars.size(), I);
	setting.setRemainderEstimation(remainder_estimation);



	double w = 0;

	// a larger initial set
//	double w = 1e-5;

	Interval init_x(1-w, 1+w), init_y(-w, w), init_z(-w, w);

	// define the initial set
	vector<Interval> box(vars.size());
	box[x_id] = init_x;
	box[y_id] = init_y;
	box[z_id] = init_z;


	Flowpipe initialSet(box);


	vector<Constraint> safeSet;

	/*
	 * The structure of the class Result_of_Reachability is defined as below:
	 * nonlinear_flowpipes: the list of computed flowpipes
	 * tmv_flowpipes: translation of the flowpipes, they will be used for further analysis
	 * fp_end_of_time: the flowpipe at the time T
	 */
	Result_of_Reachability result;

	// run the reachability computation
	clock_t begin, end;
	begin = clock();

	ode.reach(result, initialSet, 40, setting, safeSet);

	end = clock();
	printf("time cost: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);

	if(!result.isCompleted()) // if the flowpipes are not successfully computed
	{
		printf("Flowpipe computation is terminated due to the large overestimation.\n");
	}

	if(result.isSafe()) // if all of the flowpipes are safe
	{
		printf("All flowpipes are safe.\n");
	}
	else if(result.isUnsafe()) // if there is an unsafe flowpipe detected, the computation terminates immediately
	{
		printf("The last flowpipe is unsafe.\n");
	}
	else // there is no unsafe flowpipe found, but some of the flowpipes intersect the unsafe set
	{
		printf("The safety is unknown.\n");
	}


	// flowpipes should be translated to single Taylor model vectors before plotting
	result.transformToTaylorModels(setting);


	Plot_Setting plot_setting(vars);
	plot_setting.printOn();
	plot_setting.setOutputDims("t", "x");
	plot_setting.plot_2D_interval_GNUPLOT("./", "robertson_t_x", result.tmv_flowpipes, setting);

	plot_setting.setOutputDims("t", "y");
	plot_setting.plot_2D_interval_GNUPLOT("./", "robertson_t_y", result.tmv_flowpipes, setting);

	plot_setting.setOutputDims("t", "z");
	plot_setting.plot_2D_interval_GNUPLOT("./", "robertson_t_z", result.tmv_flowpipes, setting);




	return 0;
}

