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
using namespace gravity;

#ifdef USE_MOSEK
double mosek_reduce(Net* _graph, double _K) {
    
    // generate chordal extension.
    Net* chordal_extension = _graph->get_chordal_extension();
    
    // build the model
    mosek::fusion::Model:: t M  = new mosek::fusion::Model("mink_reduce");
    auto _M = monty::finally([&](){M->dispose();});
    M->setLogHandler([=](const std::string & msg){std::cout << msg << std::flush;});

    mosek::fusion::Variable::t Y = M->variable("Y", mosek::fusion::Domain::inPSDCone(_graph->nodes.size()));

    int i = 0, j =0;
    M->constraint(Y->diag(), mosek::fusion::Domain::equalsTo(1.0));
    for (auto a: chordal_extension->arcs){
        i = (a->_src)->_id;
        j = (a->_dest)->_id;
        M->constraint("", Y->index(i, j), mosek::fusion::Domain::greaterThan(-1/(_K-1)));
    }
    
    monty::rc_ptr< ::mosek::fusion::Expression >  expr= mosek::fusion::Expr::constTerm(0.0);
    for (auto a: _graph->arcs) {
        i = (a->_src)->_id;
        j = (a->_dest)->_id;
        if (i <= j){
            expr = mosek::fusion::Expr::add(expr,mosek::fusion::Expr::mul(a->_weight*(_K-1)/_K,Y->index(i,j)));
        }
        else{
            expr = mosek::fusion::Expr::add(expr,mosek::fusion::Expr::mul(a->_weight*(_K-1)/_K,Y->index(j,i)));
        }
        expr = mosek::fusion::Expr::add(expr,a->_weight/_K);
    }
    
    M->objective("obj", mosek::fusion::ObjectiveSense::Minimize, expr);
    M->solve();
    //std::cout << Y->toString() << endl;
    std::cout << expr->toString() <<endl;
    
    std::cout << "Cost = " << M->primalObjValue() << std::endl;
    return M->primalObjValue();
}


double mosekcode(Net* _graph, double _K) {
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
    
    monty::rc_ptr< ::mosek::fusion::Expression >  expr= mosek::fusion::Expr::constTerm(0.0);
    // expr is a pointer to the Expression.
    for (auto a: _graph->arcs) {
        i = (a->_src)->_id;
        j = (a->_dest)->_id;
        if (i <= j) {
            expr = mosek::fusion::Expr::add(expr,mosek::fusion::Expr::mul(a->_weight*(_K-1)/_K,Y->index(i,j)));
        }
        else {
            expr = mosek::fusion::Expr::add(expr,mosek::fusion::Expr::mul(a->_weight*(_K-1)/_K,Y->index(j,i)));
        }
        expr = mosek::fusion::Expr::add(expr,a->_weight/_K);

    }
    
    
    M->objective("obj", mosek::fusion::ObjectiveSense::Minimize, expr);
    M->solve();
    std::cout << Y->toString() << endl;
    std::cout << expr->toString() <<endl;
    
    std::cout << "Cost = " << M->primalObjValue() << std::endl;
    return M->primalObjValue();
}
#endif



int main (int argc, const char * argv[])
{
    double k = 3;
    ModelType mt = MIP;
    bool relax = false;
    int output = 0;
    SolverType solver= cplex;

    
    //string fname = "../../data_sets/Minkcut/toy.txt";
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
        if (strcmp(type,"MIP")==0){
            mt = MIP;
            solver = cplex;
        }
        else if (strcmp(type,"SDP")==0){
            mt = SDP;
            solver = Mosek;
        }
        else if (strcmp(type,"MIP_tree")==0){
            mt = MIP_tree;
            solver = cplex;
        }
        else if (strcmp(type,"SDP_tree")==0){
            mt = SDP_tree;
            solver = Mosek;
        }
        else if (strcmp(type,"Node_edge")==0){
            mt = Node_edge;
            solver = cplex;
        }
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
        //fname = "../../data_sets/Minkcut/random10_100.txt";
        fname = "../../data_sets/Minkcut/band100_3.txt";

        k = 3;
        relax = false;
        mt = MIP_tree;
        //mt = Node_edge;
        solver= cplex;
    }
    
    // Graph generation and clique-tree generation.
    Net* graph = new Net();
    graph->readrudy(fname);
    Minkmodel mymodel(mt,graph,k,solver);
    
    if (mt == MIP_tree || mt == SDP_tree)
    {
        mymodel._chordal_extension = graph->get_chordal_extension();
        graph->get_clique_tree_prim();
    }
    
    double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();

    mymodel.build();
    mymodel.solve(output,relax);
    
    
    double wall1 = get_wall_time();
    double cpu1  = get_cpu_time();
    
    cout << "\nWall clock computing time =  " << wall1 - wall0 << "\n";
    cout << "CPU computing time =  " << cpu1 -cpu0<< "\n";
    //mymodel.construct_fsol();
    
    ofstream outfile("MkP_result.txt", ios_base::app);
    if (!outfile)
        cerr << "Oops! Uable to save session data! \n";
    else{
//      outfile << "Instance,  CPU, Value" << endl;
        outfile << k << ","
                << graph->nodes.size() <<","
                << (cpu1 - cpu0)<< ", "
                << mymodel._model._obj_val
                << endl;
    }
}
