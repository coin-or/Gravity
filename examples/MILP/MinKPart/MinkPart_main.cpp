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
#include "Minkmodel.hpp"
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
    double k = 2;
    ModelType mt = MIP_tree;

    //string fname = "../../data_sets/Minkcut/spinglass2g_44.txt";
    //string fname = "../../data_sets/Minkcut/toy.txt";
    //string fname = "../../data_sets/Minkcut/toy_kojima.txt";
    //graph->get_tree_decomp_bags(true);
    
    const char* fname;
    const char* type;
    if (argc > 2)
    {
        fname = argv[1];
        k = atoi(argv[2]);
        type = argv[3];
        if (type == "MIP")
            mt = MIP;
        else
            mt = MIP_tree;
    }
    else{
        fname = "../data_sets/Minkcut/grid2d_33.txt";
        k = 2;
        mt = MIP_tree;
    }
    
    Net* graph = new Net();
    graph->readrudy(fname);
    graph->get_clique_tree();

    //graph->get_tree_decomp_bags();
    //print bags
    //    for (int i=0; i<graph->_bags.size();i++){
    //        auto bag = graph->_bags[i];
    //        int size=bag.size();
    //        cout << "bag" << i << "= {";
    //        for (int j=0; j< size; j++)
    //            cout << bag[j]->ID << ",";
    //        cout << "}" << endl;
    //    }
    //ModelType mt = MIP;
    SolverType solver= cplex;
    
    Minkmodel mymodel(mt,graph,k,solver);
    mymodel.build();
    
    double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();
    bool relax = false;
    int output = 0;
   
    mymodel.solve(output,relax);
    
    //mymodel.zij.param<bool>::print(true);
    auto sol = (var<bool> *) mymodel._model.get_var("zij");
    sol->param<bool>::print(true);
    
    double wall1 = get_wall_time();
    double cpu1  = get_cpu_time();
    cout << "\nWall clock computing time =  " << wall1 - wall0 << "\n";
    cout << "CPU computing time =  " << cpu1 - cpu0 << "\n";
    cout << "The constructed feasible solution is: " << endl;
    mymodel.construct_fsol();
}
