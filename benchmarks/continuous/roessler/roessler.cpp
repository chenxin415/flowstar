#include "../../../flowstar-toolbox/Continuous.h"

using namespace flowstar;
using namespace std;


int main()
{
	// pool of the state variables
	Variables vars;

	int x_id = vars.declareVar("x");
	int y_id = vars.declareVar("y");
	int z_id = vars.declareVar("z");

	ODE<Real> ode({"-y - z", "x + 0.2*y", "0.2 + z*(x - 5.7)"}, vars);


	Computational_Setting setting(vars);

	setting.setFixedStepsize(0.025, 8);

	setting.setCutoffThreshold(1e-8);

	// set up the remainder estimation
	Interval I(-1e-1, 1e-1);
	vector<Interval> remainder_estimation(vars.size(), I);
	setting.setRemainderEstimation(remainder_estimation);

	// define the initial set which is a box
	Interval init_x(-0.2, 0.2), init_y(-8.6, -8.2), init_z(-0.2, 0.2);

	// a box set should be defined as a vector of intervals
	// ignored dimensions will be associated with the range 0
	vector<Interval> box(vars.size());
	box[x_id] = init_x;
	box[y_id] = init_y;
	box[z_id] = init_z;


	// translating the box initial set to a flowpipe which is a preconditioned Taylor model
	Flowpipe initialSet(box);


	// no safety specification
	vector<Constraint> safeSet;


	Result_of_Reachability result;	// data structure to keep the reachability analysis result


	clock_t begin, end;
	begin = clock();


	Symbolic_Remainder sr(initialSet, 500);

	double T = 10;	// time horizon

	ode.reach(result, initialSet, T, setting, safeSet, sr);


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


	// the following content is only used for plotting
	// verification does not require the translation from flowpipes to Taylor models

	// flowpipes should be translated to single Taylor model vectors before plotting
	result.transformToTaylorModels(setting);


	Plot_Setting plot_setting(vars);
	plot_setting.printOn();
	plot_setting.setOutputDims("x", "y");
	plot_setting.plot_2D_octagon_MATLAB("./", "roessler_x_y", result.tmv_flowpipes, setting);

	plot_setting.setOutputDims("y", "z");
	plot_setting.plot_2D_octagon_MATLAB("./", "roessler_y_z", result.tmv_flowpipes, setting);

	return 0;
}

