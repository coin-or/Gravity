

# ARMO

## Supported platforms and Dependencies

We provide a copy of the open-source code of ARMO for review purposes.
ARMO was only tested on Mac OS and Linux. ARMO's only dependency is [Cmake](https://cmake.org/) 3.2 or higher and a C++ compiler.
ARMO includes [Go-ICP](https://github.com/yangjiaolong/Go-ICP) for comparison purposes.

## Compiling


ARMO will automatically download and link to [Ipopt](https://github.com/coin-or/Ipopt), [LASlib](https://github.com/LAStools/LAStools) and [Voro++](http://math.lbl.gov/voro++/).

For a faster Ipopt, preferably build it with HSL libraries (see instructions [here](https://coin-or.github.io/Ipopt/INSTALL.html)) and make sure to specify the path to `IPOPT_ROOT_DIR` in your `bash_profile` file, e.g.:

`export IPOPT_ROOT_DIR="/Users/yourname/Dev/CoinIpopt/build"`

To run the MIQCPs with [Gurobi](http://www.gurobi.com), please install it on your system first (Gurobi offers free academic licenses [here](https://www.gurobi.com/downloads/end-user-license-agreement-academic/)).

Then enter:
`mkdir build`
`cd build`
`cmake -DGurobi=ON ..` if Gurobi is installed in your system, otherwise just enter:
`cmake ..`

Finally enter:
`make -j4`

## Running

For Registration:

`Gravity/bin/Release/ARMO Reg /path_to_toy_model.txt
/path_to_toy_data.txt MIQCP scale 0`
(To run the MIQCP model with non-proportional scaling and 0% outliers)
 or

`Gravity/bin/Release/ARMO Reg /path_to_toy_model.txt
 /path_to_toy_data.txt MIQCP noscale 10`
 (To run the MIQCP model without non-proportional scaling and 10% outliers)

or 

`Gravity/bin/Release/ARMO Reg /path_to_toy_model.txt /path_to_toy_data.txt GoICP` (To run the Go-ICP algorithm for comparison purposes)

For Boresight Alignment:

`Gravity/bin/Release/ARMO Align /path_to_Cars_model.txt /path_to_Cars_data.txt /path_to_Cars_model_sub.txt /path_to_Cars_data_sub.txt`

The first two command line arguments should point to the full point clouds, the second two should point to the subsampled point clouds.
Once finsihed, ARMO will produce LAZ files (in the same directory of the input files) that can be visualized online at [plas.io](https://plas.io/).

Registration datasets can be found under the `datasets/Registration` directory.
LiDAR datasets can be downloaded [here](https://c6cff554-9579-44a7-959e-fab75fd5d22a.usrfiles.com/archives/c6cff5_402c21969b5d4bc49a340f97607027b1.zip) (full point clouds) and [here](https://c6cff554-9579-44a7-959e-fab75fd5d22a.usrfiles.com/archives/c6cff5_e271c09cc9824d0686aed597678615ec.zip) (subsampled point clouds)

