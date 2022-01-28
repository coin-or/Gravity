Instructions to reproduce results in "Benchmarking Large-Scale ACOPF Solutions and Optimality Bounds", submitted to PESGM 2022.
-------

To select the DSP relaxation set `current=false` in your local copy of this file: https://github.com/coin-or/Gravity/blob/aaa253403f0b5fa004ddea9e7c4d26d5b94924de/examples/Optimization/NonLinear/Power/OBBT/obbt.cpp#L218
To select the DSP-RLT relaxation, leave `current=true`.

To run the DSP and DSP-RLT relaxations, build the target `obbt` by running `make obbt` after `cmake`. More detailed instructions can be found [here](https://github.com/coin-or/Gravity/blob/PES2022_Largescale/INSTALL.md).

To run the SDP-BT algorithm, go to this [branch](https://github.com/coin-or/Gravity/tree/PES2022_Largescale).

To run the upper bound, you may use the executables provided under the bin folder.

Scripts for benchmarking are available in [here](https://github.com/coin-or/Gravity/tree/PESGM2022/scripts_benchmarking).



