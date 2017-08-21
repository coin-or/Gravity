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



typedef IloArray<IloNumVarArray>   NumVarMatrix;
typedef IloArray<IloBoolVarArray> BoolVarMatrix;


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

double mosek_reduce(Net* _graph, double _K) {
    mosek::fusion::Model:: t M  = new mosek::fusion::Model("mink_reduce");
    auto _M = monty::finally([&](){M->dispose();});
    M->setLogHandler([=](const std::string & msg){std::cout << msg << std::flush;});

    mosek::fusion::Variable::t Y = M->variable("Y", mosek::fusion::Domain::inPSDCone(_graph->nodes.size()));

    int i = 0, j =0;
    M->constraint(Y->diag(), mosek::fusion::Domain::equalsTo(1.0));
    for (auto a: _graph->_chordalextension->arcs){
        i = (a->src)->ID;
        j = (a->dest)->ID;
        M->constraint("", Y->index(i, j), mosek::fusion::Domain::greaterThan(-1/(_K-1)));
    }
    
    monty::rc_ptr< ::mosek::fusion::Expression >  expr= mosek::fusion::Expr::constTerm(0.0);
    for (auto a: _graph->arcs) {
        i = (a->src)->ID;
        j = (a->dest)->ID;
        if (i <= j){
            expr = mosek::fusion::Expr::add(expr,mosek::fusion::Expr::mul(a->weight*(_K-1)/_K,Y->index(i,j)));
        }
        else{
            expr = mosek::fusion::Expr::add(expr,mosek::fusion::Expr::mul(a->weight*(_K-1)/_K,Y->index(j,i)));
        }
        expr = mosek::fusion::Expr::add(expr,a->weight/_K);
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
        i = (a->src)->ID;
        j = (a->dest)->ID;
        if (i <= j) {
            expr = mosek::fusion::Expr::add(expr,mosek::fusion::Expr::mul(a->weight*(_K-1)/_K,Y->index(i,j)));
        }
        else {
            expr = mosek::fusion::Expr::add(expr,mosek::fusion::Expr::mul(a->weight*(_K-1)/_K,Y->index(j,i)));
        }
        expr = mosek::fusion::Expr::add(expr,a->weight/_K);

    }
    
    
    M->objective("obj", mosek::fusion::ObjectiveSense::Minimize, expr);
    M->solve();
    std::cout << Y->toString() << endl;
    std::cout << expr->toString() <<endl;
    
    std::cout << "Cost = " << M->primalObjValue() << std::endl;
    return M->primalObjValue();
}

ILOLAZYCONSTRAINTCALLBACK3(CtCallback, IloExprArray, lhs, IloNumArray, rhs, IloNum, eps) {
    IloInt n = lhs.getSize();
    for (IloInt i = 0; i < n; i++) {
        //IloNum xrhs = rhs[i];
        //cout << "LHS Val: " << getValue(lhs[i]) << endl;
        if (getValue(lhs[i]) > rhs[i] + eps ) {
            IloRange cut;
            try {
                cut = (lhs[i] <= rhs[i]);
                add(cut).end();
                //rhs[i] = IloInfinity; // this avoids checking this cut again.
            }
            catch (...) {
                cut.end();
                throw;
            }
        }
    }
}

int mink_tree_lazycut(Minkmodel& mymodel){
    ILOSTLBEGIN
    IloEnv env;
    try{
        IloModel model(env);
        // create variables
        BoolVarMatrix x(env, mymodel._graph->nodes.size());
        for (int r = 0; r <  mymodel._graph->nodes.size(); r++) {
            x[r] = IloBoolVarArray(env, mymodel._graph->nodes.size());
        }
        
        
        // define triangle inequalities
        int i1,i2,i3;
        for (auto it: mymodel._ids) {
            i1 = get<0>(it);
            i2 = get<1>(it);
            i3 = get<2>(it);
            //cout << "(i1, i2, i3): " << i1 << ", " << i2 << ", " << i3 << endl;
            model.add(x[i1][i2]+x[i1][i3]-x[i2][i3] <=1);
            model.add(x[i1][i3]+x[i2][i3]-x[i1][i2] <=1);
            model.add(x[i1][i2]+x[i2][i3]-x[i1][i3] <=1);
        }

        
        // objective function
        IloExpr obj(env);
        int i=0, j=0;
        for (auto a: mymodel._graph->arcs) {
            i = (a->src)->ID;
            j = (a->dest)->ID;
            if (i <= j)
                obj += (a->weight)*x[i][j];
            else
                obj += (a->weight)*x[j][i];
        }
        model.add(IloMinimize(env,obj));
        obj.end();
        
        // Lazy cut callback
        IloExprArray lhs(env);
        IloNumArray  rhs(env);
        for (auto it: (*mymodel._cliqueid)) {
            auto key = it.first;
            auto value = it.second;
            IloExpr temp(env);
            for (int i = 0; i < value.size()-1; i++){
                auto id1 = value[i];
                for (int j = i+1 ; j< value.size(); j++){
                    auto id2 = value[j];
                    if (id1 <= id2)
                        temp += (-x[id1][id2]);
                    else
                        temp += (-x[id2][id1]);
                }
            }
            lhs.add(temp);
            rhs.add(-1);
            //temp.end();
        }

        cout << "LHS size: " << lhs.getSize() << endl;
        cout << "RHS size: " << rhs.getSize()<< endl;
        cout << "LHS 0: " << lhs[0] << endl;

        IloCplex cplex(env);
        cplex.extract(model);
        // CtCallback constructs an instance of our callback class CtCallbackI and
        // returns a handle.
        cplex.use(CtCallback(env, lhs, rhs, cplex.getParam(IloCplex::Param::Simplex::Tolerances::Feasibility)));
    
        if (cplex.solve()) {
            cout << "The number of binary variables: " << cplex.getNbinVars() << endl;
            cout << "Total number of variables: " << cplex.getNcols() << endl;
            cout << "Total number of constraints: " << cplex.getNrows() << endl;
            cout << "Solution status: " << cplex.getStatus() << endl;
            cout << " Optimal Value = " << cplex.getObjValue() << endl;
            cout << " CPLEX solution time = " << cplex.getCplexTime() << endl;
            IloNumArray vals(env);
        }
    }
    catch (IloException& ex) {
        cerr << "Error: " << ex << endl;
    }
    catch (...) {
        cerr << "Error" << endl;
    }
    env.end();
    return 0;
}

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
            solver = mosek_;
        }
        else if (strcmp(type,"MIP_tree")==0){
            mt = MIP_tree;
            solver = cplex;
        }
        else if (strcmp(type,"SDP_tree")==0){
            mt = SDP_tree;
            solver = mosek_;
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
        fname = "../../data_sets/Minkcut/band100.txt";

        k = 3;
        relax = false;
        mt = MIP_tree;
        //mt = Node_edge;
        solver= cplex;
    }
    
    // Graph generation and clique-tree generation.
    Net* graph = new Net();
    graph->readrudy(fname);
    graph->get_clique_tree();
   // mosek_reduce(graph,k); 
    Minkmodel mymodel(mt,graph,k,solver);
  //  mymodel.cliquetree_decompose();
 //   mink_tree_lazycut(mymodel);
    
    
    double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();

    mymodel.build();
    mymodel.solve(output,relax);
    
    
    double wall1 = get_wall_time();
    double cpu1  = get_cpu_time();
    
    cout << "\nWall clock computing time =  " << wall1 - wall0 << "\n";
    cout << "CPU computing time =  " << cpu1 -cpu0<< "\n";
    //mymodel.construct_fsol();
    
    ofstream outfile("SDP.txt", ios_base::app);
    if (!outfile)
        cerr << "Oops! Uable to save session data! \n";
    else{
//      outfile << "Instance,  CPU, Value" << endl;
        outfile << graph->nodes.size() <<","<< (cpu1 - cpu0)
                << ", "<< mymodel._model._obj_val
                << endl;
    }
}
