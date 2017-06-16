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
    
    // instance generation.
    Net graph;
    string fname = "../../data_sets/Minkcut/grid2d_22.txt";
    graph.readrudy(fname);
    cout<< "Num_nodes: " << graph.nodes.size() << endl;
    int n = graph.nodes.size();
    
    /** MLP model by Chopra and Rao (1995)**/
    Model MIP;
    var<bool> zij("z");
    MIP.add_var(zij^(n*(n-1)/2));
    constant<int> ones(1);
    func_ obj_MIP;
    int i=0, j=0;
    for (auto a: graph.arcs){
        i = (a->src)->ID;
        j = (a->dest)->ID;
        if (i <= j)
            obj_MIP += (a->weight)*zij(i,j);
        else
            obj_MIP += (a->weight)*zij(j,i);
    }
    obj_MIP.print();
    
    /** constraints **/
    for (auto i=0; i<n-1; i++)
        for (auto h=i+1; h<n; h++)
            for (auto j=h+1; j<n    ;j++){
                Constraint Triangle1("Triangle1("+to_string(i)+","+to_string(h)+ ","+to_string(j)+")");
                Triangle1 = zij(i,h)+zij(h,j)-zij(i,j);
                Constraint Triangle2("Triangle2("+to_string(i)+","+to_string(h)+ ","+to_string(j)+")");
                Triangle2 = zij(i,h)+zij(i,j)-zij(h,j);
                Constraint Triangle3("Triangle3("+to_string(i)+","+to_string(h)+ ","+to_string(j)+")");
                Triangle3 = zij(i,j)+zij(h,j)- zij(i,h);
                MIP.add_constraint(Triangle1 <=1);
                MIP.add_constraint(Triangle2 <=1);
                MIP.add_constraint(Triangle3 <=1);
            }
    
    // K+1 subsets.
 //   MIP.print_constraints();
    
    for (auto i=0; i<n-1; i++)
        for (auto h=i+1; h<n; h++)
            for (auto j=h+1; j<n;j++){
                Constraint Clique("Clique("+to_string(i)+","+to_string(h)+ ","+to_string(j)+")");
                Clique = zij(i,h) + zij(h,j) + zij(i,j);
                MIP.add_constraint(Clique >=1);
            }
    MIP.print_constraints();
    MIP.set_objective(min(obj_MIP));
    
    solver s_mip(MIP,cplex);
    s_mip.run();
    
    /**  relaxation model for Minmum k-cut probelm **/
    Model relax;
    var<double> Xij("Xij", -1/(k-1), 1); // i<j
    relax.add_var(Xij^(n*(n-1)/2));
    
    graph.get_tree_decomp_bags();
    cout << "total bags: " << graph._bags.size() << endl;
    func_ obj;
    
    for (auto a: graph.arcs){
        i = (a->src)->ID;
        j = (a->dest)->ID;
        if (i <= j)
            obj += a->weight*((k-1)*Xij(i,j) + 1)/k;
        else
            obj += a->weight*((k-1)*Xij(j,i) + 1)/k;
    }
    
    obj.print();
    relax.set_objective(min(obj));
    
    unsigned i1, i2, i3;
    set<tuple<int,int,int>> ids;
    for (i = 0; i < graph._bags.size(); i++){
        auto bag = graph._bags.at(i);
        if (bag.size()<3) {
            continue;
        }
        for (int j = 0; j < bag.size()-2; j++){
            i1 = bag[j]->ID;
            i2 = bag[j+1]->ID;
            i3 = bag[j+2]->ID;
            assert(i2>i1 && i3>i2);
            if(ids.count(make_tuple(i1, i2, i3))==0){
                ids.insert(make_tuple(i1, i2, i3));
            }
            else {
                continue;
            }
            Constraint SDP3("SDP3("+to_string(i1)+","+to_string(i2)+","+to_string(i3)+")");
            SDP3 = -2*Xij(i1,i2)*Xij(i2,i3)*Xij(i1,i3);
            SDP3 -= 1;
            SDP3 += power(Xij(i1,i2),2);
            SDP3 += power(Xij(i1,i3),2);
            SDP3 += power(Xij(i2,i3),2);
     //       SDP3.print();
     //       relax.add_constraint(SDP3);
        }
    }
    
    /** constraints **/
    for (auto i=0; i<n; i++)
        for (auto h=i+1; h<n; h++)
            for (auto j=h+1; j<n    ;j++){
                Constraint Triangle1("Triangle1("+to_string(i)+","+to_string(h)+ ","+to_string(j)+")");
                Triangle1 = Xij(i,h)+Xij(h,j)-Xij(i,j);
                Constraint Triangle2("Triangle2("+to_string(i)+","+to_string(h)+ ","+to_string(j)+")");
                Triangle2 = Xij(i,h)+Xij(i,j)-Xij(h,j);
                Constraint Triangle3("Triangle3("+to_string(i)+","+to_string(h)+ ","+to_string(j)+")");
                Triangle3 = Xij(i,j)+Xij(h,j)- Xij(i,h);
                Triangle3.print();
                relax.add_constraint(Triangle1 <= 1);
                relax.add_constraint(Triangle2 <=1);
                relax.add_constraint(Triangle3 <=1);
            }
    
  
    for (auto i=0; i<n-1; i++)
        for (auto h=i+1; h<n; h++)
            for (auto j=h+1; j<n;j++){
                Constraint Clique("Clique("+to_string(i)+","+to_string(h)+ ","+to_string(j)+")");
                Clique = Xij(i,h) + Xij(h,j) + Xij(i,j);
                Clique.print();
                relax.add_constraint(Clique >=-k/2);
            }

    
    /* Constraints declaration */
    for (int i = 0; i < n; i++){
        for (int j = i+1; j < n; j++){
            Constraint SOCP("SOCP("+to_string(i)+","+to_string(j)+")");
//            SOCP =  Xij(i,j)*Xij(i,j) - Xii(i)*Xii(j);
            SOCP =  Xij(i,j)*Xij(i,j) -1 ;
  //          SOCP.print();
 //           relax.add_constraint(SOCP<=0);
        }
    }

    solver s_relax(relax,cplex);
    double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();
    cout << "Running the SOCP+SDP cut relaxation\n";
    s_relax.run();
    double wall1 = get_wall_time();
    double cpu1  = get_cpu_time();
    relax.print_solution();
    
    cout << "Done running the SOCP relaxation\n";
    cout << "\nWall clock computing time =  " << wall1 - wall0 << "\n";
    cout << "CPU computing time =  " << cpu1 - cpu0 << "\n";
}

