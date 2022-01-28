Intructions to reproduce results in 
"Benchmarking Large-Scale ACOPF Solutions and Optimality Bounds", submitted to PESGM 2022

To run the SDP-BT algorithm highlighted in this work:
Please build target obbt. To run the MPI version of the code set -DOpenMPI=ON while compiling the code. Code compile and build instructions are avaiable in the README.md file.

To run the DSP and DSP-RLT relaxations (lower bounds) highlighted in this work:
Please build target obbt of branch PES2022_Largescale (https://github.com/coin-or/Gravity/tree/PES2022_Largescale). 
To select the DSP-RLT relaxation set current=true in obbt.cpp (line 218)
To select the DSP relaxation set current=false in obbt.cpp (line 218) 

To run the upper bound you may use the executable .

Scripts for benchmarking are available in the folder https://github.com/coin-or/Gravity/tree/PESGM2022/scripts_benchmarking




