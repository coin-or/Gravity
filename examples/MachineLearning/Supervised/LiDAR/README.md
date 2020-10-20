

## ARMO
Follow the instructions to compile Gravity [here](https://github.com/coin-or/Gravity) then run using the following command:

For Registration:

`Gravity/bin/Release/lidar Reg toy_model.txt toy_data.txt ARMO global` (This will run the nonconvex MIQCP model)

or 

`Gravity/bin/Release/lidar Reg toy_model.txt toy_data.txt ARMO global convex` (This will run the convex MIQCP relaxation model)

or

`Gravity/bin/Release/lidar Reg toy_model.txt toy_data.txt GoICP` (This will run the Go-ICP algorithm)

For Boresight Alignment:

`Gravity/bin/Release/lidar Align Cars_model.txt Cars_data.txt Cars_model_sub.txt Cars_data_sub.txt`

The first two command line arguments should point to the full point clouds, the second two should point to the subsampled point clouds.

Datasets can be downloaded [here](https://c6cff554-9579-44a7-959e-fab75fd5d22a.usrfiles.com/archives/c6cff5_402c21969b5d4bc49a340f97607027b1.zip) (full point clouds) and [here](https://c6cff554-9579-44a7-959e-fab75fd5d22a.usrfiles.com/archives/c6cff5_e271c09cc9824d0686aed597678615ec.zip) (subsampled point clouds)

