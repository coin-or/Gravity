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
2) open the file CMakeLists.txt and turn the solvers you have to on. For instance, if you would like to use ipopt you should go to line 11 and change it to: option(ENABLE_IPOPT "Link to IPOPT libraries" ON)
3) If your solver is not installed in a default location, you need to add the following to your bash_profile:
export SOLVERNAME_ROOT_DIR="your_location". For instance, I have the following: export IPOPT_ROOT_DIR="/Users/hlh/Dev/CoinIpopt/build"
4) You're now ready to compile everything, just enter:

mkdir build

cd build

cmake ..

make -j4

This will install all the models found under Gravity/examples.
The corresponding binaries will then appear under Gravity/bin.

