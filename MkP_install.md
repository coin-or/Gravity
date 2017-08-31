# Installation Guide

Authors: Guanglei Wang and Hassan Hijazi 
Affiliations: The Australian National University, ACTON 2601, Canberra, Australia; Los Alamos National Laboratory, New Mexico, USA.

## Platform
We have successfully compiled and tested our minimum-k partition (MkP) code on both Mac OS and Ubuntu. So users are advised to use our code on Mac OS or Ubuntu.
   
## Solvers
For numerical tests of the MkP formulations in the paper `Exploiting sparsity for the min k-partition problem`, you need to install and set the following solvers: 

* CPLEX: Academic users can download it from [cplex-url][].

[cplex-url]: https://ibm.onthehub.com/WebStore/OfferingDetails.aspx?o=6fcc1096-7169-e611-9420-b8ca3a5db7a1&cm_mc_uid=02767726627915009646862&cm_mc_sid_50200000=1504137572 

* MOSEK: Download Mosek 8 from [mosek-url].

[mosek-url]: https://www.mosek.com/downloads/details/10/

* Install and setup Mosek Fusion API for C++ by following instructions in [fusion-url][]. Additionally, you need to set the enviroment variable `MOSEK_ROOT_DIR` for your operating system. For instance, if you install MOSEK 8 under `DIR/mosek8/` in a linux platform, you need to specify `export MOSEK_ROOT_DIR = DIR/mosek8/8/tools/platform/linux64x86/` in your `.bashrc` file. 

[fusion-url]: https://docs.mosek.com/8.0/cxxfusion/install.html 

To avoid compilation problem, one may need to copy Mosek libraries `libfusion64.*, libmosek64.*, libcilkrts.*`under the folder `MOSEK_ROOT_DIR/bin/` to `usr/local/lib/`. 

## Install min\_k\_part
Follow the steps below to compile the executable file. 

1. Unzip the file `Gravity.zip` and go to the folder `Gravity`. 
2. Create a build directory by typing: `mkdir build`. 
3. Enter this directory in order to perform the compilation: `cd build`
4. Type: `cmake ..`. 
5. Compile the MkP executable program by typing: `make min_k_part`.
6. Find the executable file `min_k_part` in the folder `../bin/`. 

## Tests
You are now ready to perform numerical tests stated in the paper.

### Test instances

* All test instances mentioned in the paper can be found in the folder `Gravity/data_sets/Minkcut/` and they are in `.txt` format.  

* For each problem instance, they are generally named by `graph_type + the number of vertices in the graph`. For instance, `spinglass2g_1010.txt` refers to the `spinglass2g` problem instances with `10 x 10` vertices. 

	* Since the size of the generated `band` graphs are also related to the partition parameter `k`. We additionally add `_3` or `_4` to specify the instance is used for min-3 partition or min-4 partition test instances.  

### Syntax
The executable file `min_k_part` takes four ordered input arguments, namely:

1. the input graph data `*.txt`
2. the integer parameter `k` (k > 2, we tested k = 3, 4 in our numerical tests) 
3. the formulation name (`MIP` or `Node_edge` or `MIP_tree` or `SDP` or `SDP_tree`) 
4. the boolean type for continuous relaxation (`true` for relax, `false` for integer program). 

Note that `MIP` corresponds to Model 1 in the paper, `Node_edge` to Model 2, `SDP` to Model 3, `MIP_tree` to Model 4, `SDP_tree` to Model 5. 

All the numerical tests results will also be outputed to the file `MkP_result.txt`. 

### Illustration
For example, if you are in the folder `Gravity/build` and want to valid the result of `spinglass2g` with `k = 3` and `|V|= 6 x 6`in Table 1, you need to type the following to get results for Model 1

`../bin/min_k_part ../data_sets/Minkcut/spinglass2g_66.txt 3 MIP true`
 `../bin/min_k_part ../data_sets/Minkcut/spinglass2g_66.txt 3 MIP false`

 
 and  
  
 `../bin/min_k_part ../data_sets/Minkcut/spinglass2g_66.txt 3 MIP_tree true`

 `../bin/min_k_part ../data_sets/Minkcut/spinglass2g_66.txt 3 MIP_tree false`

to get results for Model 4. 

You should get similar results as below in `Gravity/build/MkP_result.txt`:
> 3, 36, 12.1644, -2.14607e+06  
> 3, 36, 12.9085, -2.14607e+06  
> 3, 36, 0.159482, -2.14607e+06  
> 3, 36, 0.368595, -2.14607e+06

The first column is parameter k; the second is the number of nodes which is `6x6=36`; the third column is the cpu time in seconds; the fourth is the optimal value, by which you can calculate the optimality gap (0 for this test instance) for their continuous relaxations. 

Since the data in `MkP_result.txt`is separated by comma,  you can open this `csv (comma-separated-value)` file using MS Excel and calculate the optimality gap easily. 

### Run multiple tests using the bash script
If you want to run multiple the problem instances in the paper, we provide you a bash file called `exp.sh`. 

1. Copy the `exp.sh` file to your build folder. 

2. Modify the content of the `exp.sh` file to run corresponding set of numerical experiments. 

3. In your build folder, type `./exp.sh`. 

4. check the result in `MkP_result.txt`. 


## Contact 

Guanglei Wang: <guanglei.wang@anu.edu.au>

Hassan Hijazi: <hassan.hijazi@anu.edu.au>


