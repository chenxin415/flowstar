#include "../../../flowstar-toolbox/Continuous.h"

using namespace flowstar;
using namespace std;


int main()
{
	Variables vars;

	int x_id = vars.declareVar("x");
	int y_id = vars.declareVar("y");
	int z_id = vars.declareVar("z");

	ODE<Real> ode({"10*(y - x)", "x*(28 - z) - y", "x*y - 2.67*z"}, vars);

	Computational_Setting setting(vars);

	// using the following statement to specify the stepsize as 0.01 and the order as 7
	setting.setFixedStepsize(0.007, 8);

	setting.setCutoffThreshold(1e-12);


	double w = 0.001; // radius of the initial set

	// define the initial set which is a box
	Interval init_x(15-w,15+w), init_y(15-w,15+w), init_z(36-w,36+w);

	// ignored dimensions will be associated with the range 0
	vector<Interval> box(vars.size());
	box[x_id] = init_x;
	box[y_id] = init_y;
	box[z_id] = init_z;


	// creating the initial set using the box
	Flowpipe initialSet(box);


	// no unsafe constraints
	vector<Constraint> unsafeSet;


	Result_of_Reachability result;


	clock_t begin, end;
	begin = clock();


	Symbolic_Remainder sr(initialSet, 1000);

	double T = 7;

	ode.reach(result, initialSet, T, setting, unsafeSet, sr);


	end = clock();
	printf("time cost: %lf\n", (double)(end - begin) / CLOCKS_PER_SEC);

	if(!result.isCompleted()) // if the flowpipes are successfully computed
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

	// flowpipes should be translated to single Taylor model vectors before plotting
	result.transformToTaylorModels(setting);


	Plot_Setting plot_setting(vars);
	plot_setting.printOn();
	plot_setting.discreteOutput();	// only plotting the overapproximations at the end of time steps
	plot_setting.setOutputDims("x", "y");
	plot_setting.plot_2D_octagon_GNUPLOT("./", "lorenz_x_y", result.tmv_flowpipes, setting);

	plot_setting.setOutputDims("y", "z");
	plot_setting.plot_2D_octagon_GNUPLOT("./", "lorenz_y_z", result.tmv_flowpipes, setting);

	return 0;
}

