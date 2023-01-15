#include "../../../flowstar-toolbox/Continuous.h"

using namespace flowstar;
using namespace std;


int main()
{
	Variables vars;

	int r_id = vars.declareVar("r");
	int theta_id = vars.declareVar("theta");
	int dr_id = vars.declareVar("dr");
	int dtheta_id = vars.declareVar("dtheta");

	ODE<Real> ode({"dr", "dtheta", "r*dtheta^2 + 9.8*cos(theta) - 2*(r - 1)", "(-2*dr*dtheta - 9.8*sin(theta))/r"}, vars);


	Computational_Setting setting(vars);

	setting.setFixedStepsize(0.05, 5);	// fixed stepsize is 0.05, fixed order is 5

	// remainder estimation
	Interval I(-1e-1, 1e-1);
	vector<Interval> remainder_estimation(vars.size(), I);
	setting.setRemainderEstimation(remainder_estimation);

	double w = 0.01; //radius of the initial set

	// define the initial set which is a box
	Interval init_r(1.2-w, 1.2+w), init_theta(0.5-w, 0.5+w);

	// ignored dimensions will be associated with zero range
	vector<Interval> box(vars.size());
	box[r_id] = init_r;
	box[theta_id] = init_theta;



	Flowpipe initialSet(box);


	// no safety specification
	vector<Constraint> safeSet;


	Result_of_Reachability result;

	// run the reachability computation
	clock_t begin, end;
	begin = clock();

	Symbolic_Remainder sr(initialSet, 200);

	double T = 30;	// time horizon

	ode.reach(result, initialSet, T, setting, safeSet, sr);


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


	// the following content is only used for plotting
	// verification does not require to translate flowpipes to Taylor models

	// flowpipes should be translated to single Taylor model vectors before plotting
	result.transformToTaylorModels(setting);


	Plot_Setting plot_setting(vars);
	plot_setting.printOn();
	plot_setting.setOutputDims("theta", "r");
	plot_setting.plot_2D_octagon_GNUPLOT("./", "spring_pendulum_1", result.tmv_flowpipes, setting);

	plot_setting.setOutputDims("dtheta", "dr");
	plot_setting.plot_2D_octagon_GNUPLOT("./", "spring_pendulum_2", result.tmv_flowpipes, setting);

	return 0;
}

