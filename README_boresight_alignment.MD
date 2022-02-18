Instructions to reproduce results in "GLOBALLY OPTIMAL BORESIGHT CALIBRATION OF UAV-LIDAR SYSTEMS", submitted to ICIP 2022.
-------
Boresight alignment of Car data-set

<img src="https://raw.githubusercontent.com/coin-or/Gravity/Align/data_sets/LiDAR/Cars_before.png" style="max-width:100%;"/>
<img src="https://raw.githubusercontent.com/coin-or/Gravity/Align/data_sets/LiDAR/Cars_after.png" style="max-width:100%;"/>

Boresight alignment of Tent data-set

<img src="https://raw.githubusercontent.com/coin-or/Gravity/Align/data_sets/LiDAR/Tent_before.png" style="max-width:100%;"/>
<img src="https://raw.githubusercontent.com/coin-or/Gravity/Align/data_sets/LiDAR/Tent_after.png" style="max-width:100%;"/>

Boresight alignment of Truck data-set

<img src="https://raw.githubusercontent.com/coin-or/Gravity/Align/data_sets/LiDAR/Truck_before.png" style="max-width:100%;"/>
<img src="https://raw.githubusercontent.com/coin-or/Gravity/Align/data_sets/LiDAR/Truck_after.png" style="max-width:100%;"/>

Code Tested on Mac and Linux

Algorithm options: nsBB, aGS, gurobi

External Dependencies:

Please install [LASLIB](https://github.com/LAStools/LAStools)

`export LASlib_ROOT_DIR="/Users/yourname/Dev/LAStools/LASlib"`

For nsBB ONLY

Please install [openGJK](https://github.com/MattiaMontanari/openGJK) 

`export GJK_ROOT_DIR="/Users/yourname/Dev/openGJK"`

For nsBB ONLY

Please download and extract (`tar -xf archivename`) [Eigen3](https://gitlab.com/libeigen/eigen/-/releases/3.3.9) Eigen3 is header-only and does not need installation 

`export EIGEN3_ROOT_DIR ="/Users/yourname/Dev/eigen-3.3.9"`

For nsBB and gurobi ONLY

Please install [Gurobi 9.5](https://www.gurobi.com/) (not required for aGS) For troubleshooting see [here](https://support.gurobi.com/hc/en-us/articles/360039093112-How-do-I-resolve-undefined-reference-errors-while-linking-Gurobi-in-C-)

`export GRB_LICENSE_FILE="/Users/yourname/Dev/gurobi.lic"`

For Mac,

`export GUROBI_HOME="/Library/gurobi950/macos_universal2"`

For Linux,

`export GUROBI_HOME="$HOME/Dev/gurobi950/linux64"`

`export PATH=$PATH:"$HOME/Dev/gurobi950/linux64/bin"`   

`export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/Dev/gurobi950/linux64/lib"`


To run boresight alignment algorithms, build target `lidar` by running `make lidar` after `cmake`. Code compile and build instructions are avaiable [here](https://github.com/coin-or/Gravity/blob/Align/README.md)

CMAKE OPTIONS

To run the nsBB set `-DGurobi=ON -DGJK=ON -DEigen3=ON` while calling `cmake`. To run the MPI version of nsBB also set `-DOpenMPI=ON` while calling `cmake`. 

To run gurobi set `-DGurobi=ON` while calling `cmake`. To run the MPI version of nsBB also set `-DOpenMPI=ON` while calling `cmake`. 

Data sets used in the paper may be found [here](https://github.com/coin-or/Gravity/tree/Align/data_sets/LiDAR).

Arguments

The first argument is the .laz file name

The second argument is the algorithm: `aGS`(default) for aGS algorithm, `nsBB`, for nsBB, `gurobi` to solve the boresight MIQCQP as it is with gurobi

Alternatively, set algorithm set `algo="nsBB"` in your local copy of this [file](https://github.com/coin-or/Gravity/blob/Align/examples/MachineLearning/Supervised/LiDAR/LiDAR_main.cpp#L82/)


