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

void mosek_reduce(Net* _graph, double _K) {
    mosek::fusion::Model:: t M  = new mosek::fusion::Model("mink_reduce");
    auto _M = monty::finally([&](){M->dispose();});
    M->setLogHandler([=](const std::string & msg){std::cout << msg << std::flush;});

    //vector<mosek::fusion::Variable::t> _mosek_vars;
    mosek::fusion::Variable::t Y = M->variable("Y", mosek::fusion::Domain::inPSDCone(_graph->nodes.size()));
//
//    for (auto b: _graph->_bags){
//        auto t = M->variable(" ", mosek::fusion::Domain::inPSDCone(b.size()));
//        _mosek_vars.push_back(t);
//    }
    
int i = 0, j =0;
 M->constraint(Y->diag(), mosek::fusion::Domain::equalsTo(1.0));
    for (auto a: _graph->_chordalextension->arcs){
        i = (a->src)->ID;
        j = (a->dest)->ID;
        M->constraint("", Y->index(i, j), mosek::fusion::Domain::greaterThan(-1/(_K-1)));
            //M->constraint("", Y->index(i, j), mosek::fusion::Domain::lessThan((1.0)));
    }
    
    monty::rc_ptr< ::mosek::fusion::Expression >  expr= mosek::fusion::Expr::constTerm(_graph->arcs.size()/_K);
// expr is a pointer to the Expression.
//    int i, j;
    for (auto a: _graph->arcs) {
        i = (a->src)->ID;
        j = (a->dest)->ID;
        if (i <= j)
            expr = mosek::fusion::Expr::add(expr,mosek::fusion::Expr::mul(a->weight*(_K-1)/_K,Y->index(i,j)));
        else
            expr = mosek::fusion::Expr::add(expr,mosek::fusion::Expr::mul(a->weight*(_K-1)/_K,Y->index(j,i)));
    }
    
    M->objective("obj", mosek::fusion::ObjectiveSense::Minimize, expr);
    M->solve();
    //std::cout << Y->toString() << endl;
    std::cout << expr->toString() <<endl;
    
    std::cout << "Cost = " << M->primalObjValue() << std::endl;
}



void mosekcode(Net* _graph, double _K) {
    mosek::fusion::Model:: t M  = new mosek::fusion::Model("mink");
    auto _M = monty::finally([&](){M->dispose();});
    M->setLogHandler([=](const std::string & msg){std::cout << msg << std::flush;});
    
    
    mosek::fusion::Variable::t Y = M->variable("Y", mosek::fusion::Domain::inPSDCone(_graph->nodes.size()));
    
    int i = 0, j =0;
    M->constraint(Y->diag(), mosek::fusion::Domain::equalsTo(1.0));
    

    for (i= 0; i < _graph->nodes.size(); i ++)
        for (j = i+1; j< _graph->nodes.size(); j++){
            M->constraint("", Y->index(i, j), mosek::fusion::Domain::greaterThan(-1/(_K-1)));
        }
    
    monty::rc_ptr< ::mosek::fusion::Expression >  expr= mosek::fusion::Expr::constTerm(_graph->arcs.size()/_K);
    // expr is a pointer to the Expression.
    for (auto a: _graph->arcs) {
        i = (a->src)->ID;
        j = (a->dest)->ID;
        if (i <= j)
            expr = mosek::fusion::Expr::add(expr,mosek::fusion::Expr::mul(a->weight*(_K-1)/_K,Y->index(i,j)));
        else
            expr = mosek::fusion::Expr::add(expr,mosek::fusion::Expr::mul(a->weight*(_K-1)/_K,Y->index(j,i)));
    }
    
    M->objective("obj", mosek::fusion::ObjectiveSense::Minimize, expr);
    M->solve();
    std::cout << Y->toString() << endl;
    std::cout << expr->toString() <<endl;
    
    std::cout << "Cost = " << M->primalObjValue() << std::endl;
}

int main (int argc, const char * argv[])
{
    double k = 2;
    ModelType mt = MIP;
    bool relax = false;
    
    //string fname = "../../data_sets/Minkcut/toy.txt";
    //string fname = "../../data_sets/Minkcut/toy_kojima.txt";
    //graph->get_tree_decomp_bags(true);
    
    const char* fname;
    const char* type;
    const char* relaxation;
    if (argc > 3)
    {
        fname = argv[1];
        k = atoi(argv[2]);
        type = argv[3];
        relaxation = argv[4];
        cout << "type" << type << endl;
        if (strcmp(type,"MIP")==0)
            mt = MIP;
        else if (strcmp(type,"SDP")==0)
            mt = SDP;
        else if (strcmp(type,"MIP_tree")==0)
            mt = MIP_tree;
        else{
            cout << "invalid input" << endl;
            exit(1);
        }

        if (strcmp(relaxation,"true")==0)
            relax = true;
        else{
            relax = false;
        }
    }
    else{
        fname = "../../data_sets/Minkcut/random_55.txt";
        //fname = "../../data_sets/Minkcut/spinglass2g_55.txt";

        //fname = "../../data_sets/Minkcut/toy_kojima.txt";
        //fname = "../../data_sets/Minkcut/spinglass3g_234.txt";
        k = 3;
        relax = true;
        mt = MIP_tree;
    }
    
    Net* graph = new Net();
    graph->readrudy(fname);
    graph->get_clique_tree();

    SolverType solver= cplex;
    //SolverType solver= mosek_;

    
    Minkmodel mymodel(mt,graph,k,solver);
    double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();
    
    //mosekcode(graph, k);
    // SDP 
    mosek_reduce(graph, k);

    // CPLEX
    mymodel.build();
    int output = 0;
    mymodel.solve(output,relax);
    
    //mymodel.zij.param<bool>::print(true);
    //auto sol = (var<bool> *) mymodel._model.get_var("zij");
    //sol->param<bool>::print(true);
    
    double wall1 = get_wall_time();
    double cpu1  = get_cpu_time();
    cout << "\nWall clock computing time =  " << wall1 - wall0 << "\n";
    cout << "CPU computing time =  " << cpu1 - cpu0 << "\n";
    //mymodel.construct_fsol();
}
