//
//  Stable_Set.cpp
//  Gravity
//
//  Created by Hijazi, Hassan (Data61, Canberra City) on 6/12/17.
//
//
#include <stdio.h>
#include <iostream>
#include <string>
#include <stdio.h>
#include <cstring>
#include <fstream>
#include <gravity/Net.h>
#include <gravity/model.h>
#include <gravity/solver.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;
using namespace gravity;
int main (int argc, const char * argv[])
{
    //  Start Timers
    std::cout << "WELCOME, THIS IS AN IMPLEMENTATION OF THE POOLING PROBLEM AND SOME OF ITS RELAXATIONS IN GRAVITY\n";
    std::cout << "Understanding the numerical limits of your machine:" << endl;
    std::cout << "type\tlowest\thighest\n";
    
    std::cout << "bool\t"
    << std::numeric_limits<bool>::lowest() << '\t'
    << std::numeric_limits<bool>::max() << '\n';
    std::cout << "short\t"
    << std::numeric_limits<short>::lowest() << '\t'
    << std::numeric_limits<short>::max() << '\n';
    std::cout << "unsigned\t"
    << std::numeric_limits<unsigned>::lowest() << '\t'
    << std::numeric_limits<unsigned>::max() << '\n';
    std::cout << "int\t"
    << std::numeric_limits<int>::lowest() << '\t'
    << std::numeric_limits<int>::max() << '\n';
    std::cout << "long int\t"
    << std::numeric_limits<long int>::lowest() << '\t'
    << std::numeric_limits<long int>::max() << '\n';
    std::cout << "double\t"
    << std::numeric_limits<double>::lowest() << '\t'
    << std::numeric_limits<double>::max() << '\n';
    std::cout << "long double\t"
    << std::numeric_limits<long double>::lowest() << '\t' << std::numeric_limits<long double>::max() << '\n';
    unsigned i, j, i1, i2, i3;
    
    SolverType stype = cplex;
    double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();
    cout << "Running the NONCONVEX model\n";
    //s.run();
    double wall1 = get_wall_time();
    double cpu1  = get_cpu_time();
    cout << "Done running the NONCONVEX model\n";
    cout << "Wall clock computing time =  " << wall1 - wall0 << "\n";
    cout << "CPU computing time =  " << cpu1 - cpu0 << "\n";
    return 0;
}
