## <img src="images/flowstar.png" alt="flowstar" width='50'> Flow* Toolbox - <br> *A Platform for Modeling and Analysis of Cyber-Physical Systems*

# Introduction


This is the homepage of the toolbox version of Flow\*. The first version of Flow\* was released in the year of 2013, and improved in 2015 (version 1.2.0) and 2017 (version 2.1.0). The purpose of releasing a toolbox version is to provide a more flexible way to model and analyze cyber-physical systems (CPS), and expose the key functions to the tools for verifying more complex systems, such as the CPS with machine learning components. The main data structures in the toolbox version are completely re-designed and implemented such that the performance is **at least 10x faster** than the version 2.1.0.


### How does the toolbox work?

The Flow\* toolbox does not have a specific interface, it is compiled as a static library. A verification or reachability task should be described as a C++ file which can be compiled with the Flow\* libarary.


<img src='images/usingflowstar.png' width='500'>


### Main Functionalities


**Time-bounded reachability computation.** Flow* computes Taylor Model (TM) flowpipes for a discrete, continuous or hybrid dynamics over a bounded time interval. TMs are *function* rather than pure range overapproximations. As it is illustrated by the following figure, a TM overapproximates the flowmap of a deterministic dynamics, and its range forms an overapproximation of the reachable set.


<img src='images/tmflowpipe.png' width='500'>





**Safety verification.** Conservatively checking the intersection of the TM flowpipe ranges with a given safe or unsafe set.




**Configuration independent relational abstraction for linear dynamics.** The Flow\* toolbox does not directly compute reachable set overapproximations for Linear Time-Invariant (LTI) or Linear Time-Varying (LTV) ODEs. Instead, it generates TM overapproximations for the flowmap functions that maps any initial state to its reachable state at a time. Such abstractions may be reused in verification tasks for different system settings.




The tool can also be used to find potential counterexamples when the overapproximated flowmap is deterministic. E.g., when an unsafe intersection is detected for a TM $p(x_0,t)$, you may use the domain contraction function to reduce $x_0$'s range, which is the initial set, and find an overapproximation of the unsafe initial states.

--
**More content will be added.**
