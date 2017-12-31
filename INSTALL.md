DEPENDENCIES
-------
cmake: http://www.cmake.org (Version 3.2 or better)

OPTIONAL
-------

Ipopt: https://projects.coin-or.org/Ipopt

Bonmin: https://projects.coin-or.org/Bonmin

Gurobi: http://www.gurobi.com

Cplex: https://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/

Mosek: https://www.mosek.com

-------

Follow these simple instructions:
1) cd Gravity

2) If your solver is not installed in a default location, you need to add the following to your bash_profile:
export SOLVERNAME_ROOT_DIR="your_location", where your location contains include/headers and lib/libraries.

For instance, I have the following: 
export IPOPT_ROOT_DIR="/Users/hh/Dev/CoinIpopt/build"

export MOSEK_ROOT_DIR="/Users/hh/mosek/8/tools/platform/osx64x86"

Note: If Cplex and Guraobi were installed using the standard installer, Cmake should automatically locate the latest version on your system, no need to add anything to your bash_profile.

3) You're now ready to compile everything, just enter:

mkdir build

cd build

cmake -DENABLE_IPOPT=ON -DENABLE_GUROBI=ON -DENABLE_CPLEX=ON -DENABLE_BONMIN=ON ..
(All dependencies are switched off by default, to enable a solver that is installed on your system, add the flag -DENABLE_$SOLVERNAME$=ON as shown above)

Note: To build an Xcode project append -G Xcode to the command above

make -j4

This will install all the models found under Gravity/examples.

The corresponding binaries will then appear under Gravity/bin.

All executables can be called with the flag --help or -h to get help on running instructions.
