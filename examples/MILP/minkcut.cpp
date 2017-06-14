//
//  MinKpartition.cpp
//  Gravity
//
//  Created by Guanglei Wang on 13/6/17.
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
    auto k = 2; // input
    auto n = 3; // number of nodes.
    
    /** MLP model by Chopra and Rao (1995)**/
    Model MIP;
    var<bool> zij("z");
    MIP.add_var(zij^(n*(n-1)/2));
    constant<int> ones(1);
    func_ obj_MIP = ones.tr()*zij;
    
    /** constraints **/
    for (auto i=1; i<n-2; i++)
        for (auto h=i+1; h<n-1; h++)
            for (auto j=h+1; j<n;j++)
                


    /**  relaxation model for Minmum k-cut probelm **/
    
    Model relax;
    var<float> Xii("Xii", 1, 1);
    var<int> Xij("Xij", -1/(k-1),INFINITY); // i<j
    relax.add_var(Xii^n);
    relax.add_var(Xij^(n*(n-1)/2));
    
    
    constant<float> weight(1*(k-1)/k);
    func_ obj = weight.tr()*Xij+n*(n-1)/(2*k);

    relax.set_objective(min(obj));
    solver s(relax,cplex);
    
    
    double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();
    cout << "Running the LP relaxation\n";
    s.run();
    double wall1 = get_wall_time();
    double cpu1  = get_cpu_time();
    relax.print_solution(); //
    
    
    cout << "Done running the LP relaxation\n";
    cout << "\nWall clock computing time =  " << wall1 - wall0 << "\n";
    cout << "CPU computing time =  " << cpu1 - cpu0 << "\n";
}

