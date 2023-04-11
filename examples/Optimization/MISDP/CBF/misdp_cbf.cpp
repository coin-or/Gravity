//
//  misdp.cpp
//  misdp
//
//  Created by Smitha on 7/7/22.
//

#include <stdio.h>
#include "read_misdp.h"
#include <gravity/solver.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <fstream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;




/*returns matrix of dimension n-1, without row rowno and col colno*/
vector<vector<double>> get_minor(vector<vector<double>>& X, int rowno, int colno, int n){
    if(n==1)
        throw invalid_argument ("n=1");
    vector<vector<double>> Y(n-1);
    int count=0;
    for(auto i=0; i<n;i++){
        for(auto j=0; j<n && i!=rowno;j++){
            if(j!=colno)
                Y[count].push_back(X[i][j]);
        }
        if(i!=rowno){
            count++;
        }
    }
    return Y;
}

double determinant(vector<vector<double>> X, int n){
    double det;
    if(n==1){
        det=X[0][0];
    }
    else{
        det=0;
        for(auto i=0;i<n;i++){
            auto Y=get_minor(X, 0, i, n);
            det+=std::pow(-1, i)*X[0][i]*determinant(Y,n-1);
        }
    }
    return det;
}
using namespace gravity;
using namespace std;

int main(int argc, char * argv[]){

    string fname=string(prj_dir)+"/data_sets/MISDP/2x7_3bars.cbf";
    bool root_refine = false, add_soc=false, add_threed=false, hierarc=false;
    if(argc>=2){
        fname=argv[1];
    }
    if(argc>=3){
        root_refine = argv[2];
    }
    if(argc>=4){
        add_soc = argv[3];
    }
    if(argc>=5){
        add_threed = argv[4];
    }
    if(argc>=6){
        hierarc = argv[5];
    }
    auto m=make_shared<Model<double>>("misdp_test");
    m->add_soc=add_soc;
    m->add_threed=add_threed;
    m->add_hierarc=hierarc;
    
    auto g=CBF_read(fname.c_str(), m);
    m->print();
    double ub_solver_tol=1e-6, lb_solver_tol=1e-6, range_tol=1e-4, opt_rel_tol=1e-2, opt_abs_tol=1e6;
    unsigned max_iter=1e3;
    SolverType ub_solver_type = ipopt, lb_solver_type = gurobi;
    
    double max_time = 300;
    DebugOn("Instance "<<fname<<endl);

    solver<> sc(m,gurobi);
    bool relax = false;
    auto ts=get_wall_time();
    if(root_refine){
        //auto m2=m->copy();
        sc.run(relax = true);
       // m->print_solution();
//      upper_bound_heur=true;
//        if(upper_bound_heur){
//            auto y2d=m->get_var<int>("y");
//            auto y2=m2->get_var<int>("y");
//
//            for(auto k:*y2d.get_keys()){
//                if(y2d.eval(k)<0.5){
//                    y2.set_lb(k,0);
//                    y2.set_ub(k,0);
//
//                }
//                else{
//                    y2.set_lb(k,1);
//                    y2.set_ub(k,1);
//                }
//            }
//            solver<> s2(m2,gurobi);
//            s2.run();
//            //check_PSD(m1);
//            if(m2->_status==0){
//                vector<double> sol(m->_nb_vars);
//                m2->get_solution(sol);
//                m->set_solution(sol);
//            }
//
//        }
    }
    sc.run();
    auto tf=get_wall_time();
    m->print_solution();
    m->print_constraints_stats(1e-9);
    
    auto eig_value=m->check_PSD();
    
    string out_file_name=fname;
    auto pos=out_file_name.find_last_of("/");
    out_file_name=out_file_name.substr(pos+1);
    pos=out_file_name.find_first_of(".");
    auto out_file_name1=out_file_name.substr(0,pos);
    out_file_name=string(prj_dir)+"/results_misdp/"+out_file_name1+".txt";
    ofstream fout(out_file_name.c_str());
    fout<<out_file_name1<<"&"<<m->get_obj_val()<<"&"<<tf-ts<<"&"<<eig_value<<"&"<<m->_rel_obj_val<<"&"<<m->num_cuts[0]<<"&"<<m->num_cuts[1]<<"&"<<m->num_cuts[2]<<"&"<<m->num_cuts[3]<<"&"<<m->num_cuts[4]<<"&"<<m->num_cuts[5]<<m->num_cuts[6]<<"&"<<m->num_cuts[7]<<"&"<<m->num_cuts[8]<<"&"<<m->num_cuts[9]<<"&"<<m->num_cuts[10]<<"&"<<m->num_cuts[11]<<"\n";
    fout.close();
    cout<<out_file_name1<<" obj "<<m->get_obj_val()<<" time "<<tf-ts<<" smallest eig "<<eig_value<<" lower bound  "<<m->_rel_obj_val<<endl;
    cout<<"soc incumbent "<<m->num_cuts[0]<<" threed incumbent "<<m->num_cuts[1]<<" eig bags incumbent "<<m->num_cuts[2]<<" eig full incumbent "<<m->num_cuts[3]<<" soc mipnode "<<m->num_cuts[4]<<" threed mipnode "<<m->num_cuts[5]<<" eig bags mipnode "<<m->num_cuts[6]<<" eig full mipnode "<<m->num_cuts[7]<<"\n";
}
