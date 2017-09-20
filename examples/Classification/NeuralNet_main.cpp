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
#include "DataSet.h"

using namespace std;
using namespace gravity;
#define EPS 0.00001
#define DebugOn(x) cout << x
#define DebugOff(x)
//  Windows
#ifdef _WIN32
#include <Windows.h>
double get_wall_time(){
    LARGE_INTEGER time,freq;
    if (!QueryPerformanceFrequency(&freq)){
        //  Handle error
        return 0;
    }
    if (!QueryPerformanceCounter(&time)){
        //  Handle error
        return 0;
    }
    return (double)time.QuadPart / freq.QuadPart;
}
double get_cpu_time(){
    FILETIME a,b,c,d;
    if (GetProcessTimes(GetCurrentProcess(),&a,&b,&c,&d) != 0){
        //  Returns total user time.
        //  Can be tweaked to include kernel times as well.
        return
        (double)(d.dwLowDateTime |
                 ((unsigned long long)d.dwHighDateTime << 32)) * 0.0000001;
    }else{
        //  Handle error
        return 0;
    }
}

//  Posix/Linux
#else
#include <time.h>
#include <sys/time.h>
double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}
#endif

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
    
    
    string file = "../data_sets/classification/pendigits";
    if (argc>1) {
        file = argv[argc-1];
    }
    
    DataSet<> training_set;
    training_set.parse(file);
    training_set.print_stats();
    
    Model SoftPlus("SoftPlus");
    param<> feature;
    unsigned bound = 10;
    vector<unsigned> nb_nn = {1,2,1}; /*<Number of neurons per layer */
    unsigned nb_layers = nb_nn.size();
    var<> x("x", 0, 1);
    var<> alpha("alpha", -bound, bound);
    var<> beta("beta", -bound, bound);
    var<> teta("beta", -bound, bound);
//    x.indexed_in({1..5});
//    SoftPlus.add_var(x^(x_dim));
    
    
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
