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
	int x8_id = vars.declareVar("x8");
	int x9_id = vars.declareVar("x9");
	int x10_id = vars.declareVar("x10");
	int x11_id = vars.declareVar("x11");
	int x12_id = vars.declareVar("x12");
	int t_id = vars.declareVar("t");


	ODE<Real> ode({	"cos(x8)*cos(x9)*x4 + (sin(x7)*sin(x8)*cos(x9) - cos(x7)*sin(x9))*x5 + (cos(x7)*sin(x8)*cos(x9) + sin(x7)*sin(x9))*x6",
					"cos(x8)*sin(x9)*x4 + (sin(x7)*sin(x8)*sin(x9) + cos(x7)*cos(x9))*x5 + (cos(x7)*sin(x8)*sin(x9) - sin(x7)*cos(x9))*x6",
					"sin(x8)*x4 - sin(x7)*cos(x8)*x5 - cos(x7)*cos(x8)*x6",
					"x12*x5 - x11*x6 - 9.81*sin(x8)",
					"x10*x6 - x12*x4 + 9.81*cos(x8)*sin(x7)",
					"x11*x4 - x10*x5 + 9.81*cos(x8)*cos(x7) - 9.81 + 7.14285714285714*(x3 - 1) - 2.14285714285714*x6",
					"x10 + (sin(x7)*(sin(x8)/cos(x8)))*x11 + (cos(x7)*(sin(x8)/cos(x8)))*x12",
					"cos(x7)*x11 - sin(x7)*x12",
					"(sin(x7)/cos(x8))*x11 + (cos(x7)/cos(x8))*x12",
					"-0.92592592592593*x11*x12 - 18.51851851851852*(x7 + x10)",
					"0.92592592592593*x10*x12 - 18.51851851851852*(x8 + x11)",
					"0",
					"1"}, vars);


	Computational_Setting setting(vars);

	setting.setFixedStepsize(0.025, 4);

	setting.setCutoffThreshold(1e-6);

	Interval remainder(-1e-2,1e-2);
	vector<Interval> estimation(vars.size(), remainder);
	setting.setRemainderEstimation(estimation);

	Interval I(-0.5,0.5);

	vector<Interval> initial_box(vars.size());
	initial_box[x1_id] = I;
	initial_box[x2_id] = I;
	initial_box[x3_id] = I;
	initial_box[x4_id] = I;
	initial_box[x5_id] = I;
	initial_box[x6_id] = I;


	Flowpipe initialSet(initial_box);

	vector<Constraint> safeSet = {Constraint("x3 - 1.4", vars)}; // x3 <= 1.4

	Result_of_Reachability result;

	// run the reachability computation
	clock_t begin, end;
	begin = clock();

	Symbolic_Remainder symbolic_remainder(initialSet, 20);

	double T = 5;

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
	plot_setting.setOutputDims("t", "x3");

	plot_setting.plot_2D_interval_GNUPLOT("./", "quadrotor", result.tmv_flowpipes, setting);

	return 0;
}

