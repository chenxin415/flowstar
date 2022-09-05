#include "../../../flowstar-toolbox/Continuous.h"

using namespace flowstar;
using namespace std;


int main()
{
	// precision of the float-point numbers
	// the default setting is double
//	intervalNumPrecision = 53;


	// pool of the state variables
	Variables vars;


	// declaring the state variables
	// returned values are their IDs
	int x_id = vars.declareVar("x");
	int t_id = vars.declareVar("t");	// unnecessary if the safety property is not time related


	// creating the ODE object
	// The derivatives should be provided according to the variable declaration order:
	ODE<Real> ode({"1 - sin(x) * sqrt(log(x)) / exp(cos(x))", "1"}, vars);


	// creating the setting (necessary)
	Computational_Setting setting(vars);

	// default setting for reachability computation
	// adaptive stepsize: 0.002 ~ 0.1
	// fixed TM order: 4

	// using the following statement to specify the stepsize as 0.02 and the order as 5
//	setting.setFixedStepsize(0.02, 5);


	// define the initial set which is a box
	Interval init_x(4.8,5.2);

	// a box set should be defined as a vector of intervals
	// ignored dimensions will be associated with the range 0
	vector<Interval> box(vars.size());
	box[x_id] = init_x;


	// translating the box initial set to a flowpipe which is a preconditioned Taylor model
	Flowpipe initialSet(box);


	// there is no safe specification in this example
	vector<Constraint> safeSet;


	Result_of_Reachability result;

	// run the reachability computation
	clock_t begin, end;
	begin = clock();

	/**** Data Structure of Symbolic Remainders ***
	 * The first parameter should be the flowpipe initial set,
	 * the second one should be the size of the queues.
	 * Details are described in our RTSS 2016 paper.
	 */
	Symbolic_Remainder sr(initialSet, 100);

	double T = 5;	// time horizon

	/**** Calling the reachability analysis function ***
	 * 1st arg: Data structure for keeping the result;
	 * 2nd arg: Initial set which should be a flowpipe;
	 * 3rd arg: Time horizon;
	 * 4th arg: Computational setting;
	 * 5th arg: Safe set;
	 * 6th arg: Data structure of the symbolic remainder;
	 *
	 * The symbolic remainder can be ignored in the function call,
	 * if so, then no symbolic remainder will be used in reachability
	 * computation.
	 */
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
	// verification does not require to translate flowpipes to Taylor models

	// flowpipes should be translated to single Taylor model vectors before plotting
	result.transformToTaylorModels(setting);


	Plot_Setting plot_setting(vars);
	plot_setting.printOn();
	plot_setting.setOutputDims("t", "x");
	plot_setting.plot_2D_interval_GNUPLOT("./", "simple", result.tmv_flowpipes, setting);	// producing a GNUPLOT file
//	plot_setting.plot_2D_interval_MATLAB("./", "simple", result.tmv_flowpipes, setting);	// producing a MATLAB file

	return 0;
}

