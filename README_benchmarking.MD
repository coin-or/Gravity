Instructions to reproduce results in "Benchmarking Large-Scale ACOPF Solutions and Optimality Bounds", submitted to PESGM 2022.
-------

To run the SDP-BT algorithm, build target `obbt` by running `make obbt` after `cmake`. To run the MPI version of the code set `-DOpenMPI=ON` while calling `cmake`. Code compile and build instructions are avaiable [here](https://github.com/coin-or/Gravity/blob/PES2022_Largescale/INSTALL.md).

To run the DSP and DSP-RLT relaxations, go to this branch: [PES2022_Largescale](https://github.com/coin-or/Gravity/tree/PES2022_Largescale). 

To run the upper bound you may use the [OPOMO](https://github.com/hhijazi/OPOMO) solver.

Scripts for benchmarking are available in [here](https://github.com/coin-or/Gravity/tree/PESGM2022/scripts_benchmarking).




