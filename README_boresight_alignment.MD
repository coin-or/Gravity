Instructions to reproduce results in "GLOBALLY OPTIMAL BORESIGHT CALIBRATION OF UAV-LIDAR SYSTEMS", submitted to ICIP 2022.
-------

To run boresight alignment algorithms, build target `lidar` by running `make obbt` after `cmake`. To run the nsBB version of the code set `-DGurobi=ON` while calling `cmake`. To run the MPI version of the code set `-DOpenMPI=ON` while calling `cmake`. Code compile and build instructions are avaiable [here](https://github.com/coin-or/Gravity/blob/Align/README.md)

Data sets used in the paper may be found [here](https://github.com/coin-or/Gravity/tree/Align/data_sets/LiDAR).

To select the aGS algorithm set `algo="ub"` in your local copy of this file: https://github.com/coin-or/Gravity/blob/Align/examples/MachineLearning/Supervised/LiDAR/LiDAR_main.cpp#L82

To select the nsBB algorithm leave `algo="bb"`.

