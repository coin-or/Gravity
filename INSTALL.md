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
1) * `cd Gravity`

2) If your solver is not installed in a default location, you need to add the following to your `~/.bash_profile` (edit or create the file named `bash_profile` under your home and add the following lines):
* `export SOLVERNAME_ROOT_DIR="your_location"`, where `your_location` contains the corresponding include/headers and lib/libraries.

For instance, I have the following: 

* `export IPOPT_ROOT_DIR="/Users/hh/Dev/CoinIpopt/build"`

* `export MOSEK_ROOT_DIR="/Users/hh/mosek/8/tools/platform/osx64x86"`

Note: If Cplex and Gurobi were installed using the standard installer, Cmake should automatically locate the latest version on your system, no need to add anything to your bash_profile.

Note: For [Mosek](https://www.mosek.com/downloads/details/10/), refer to [mosek-setup](https://www.google.com.au/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&cad=rja&uact=8&ved=0ahUKEwjy0sja6oLWAhXEwLwKHQR_A5YQFggoMAA&url=http%3A%2F%2Fdocs.mosek.com%2F8.1%2Finstall%2Finstallation.html&usg=AFQjCNGEiUPE05E_5_UedXe1mmpCYOimrQ) for general setup and license information. 
 

You also need to install and setup Mosek Fusion API for C++ by following instructions in [fusion-url](https://www.google.com.au/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&cad=rja&uact=8&ved=0ahUKEwjx1deH34LWAhWLw7wKHWi2An8QFggoMAA&url=http%3A%2F%2Fdocs.mosek.com%2F8.0%2Fcxxfusion%2Finstall.html&usg=AFQjCNFwhQErdOsuD8iSIcDbMo3IERbhdA
). 

* Additionally, you need to set the enviroment variable `MOSEK_ROOT_DIR` and `DYLD_LIBRARY_PATH`for your operating system. 

For instance, if you install MOSEK 8 under `DIR/mosek8/` on a linux machine, you need to add the following to your `.bash_profile` file:

* `export MOSEK_ROOT_DIR = DIR/mosek/8/tools/platform/linux64x86/`  
* `export DYLD_LIBRARY_PATH= DIR/mosek/8/tools/platform/linux64x86/bin:$DYLD_LIBRARY_PATH`

3) You're now ready to compile everything, just enter:

* `mkdir build`

* `cd build`

* `cmake -DENABLE_IPOPT=ON -DENABLE_GUROBI=ON -DENABLE_CPLEX=ON -DENABLE_BONMIN=ON ..`
(All dependencies are switched off by default, to enable a solver that is installed on your system, add the flag -DENABLE_$SOLVERNAME$=ON as shown above)

Note: To build an Xcode project append -G Xcode to the command above

* `make -j4`

This will install all the models found under `Gravity/examples`.

The corresponding binaries will then appear under `Gravity/bin`.

All executables can be called with the flag --help or -h to get help on running instructions.
