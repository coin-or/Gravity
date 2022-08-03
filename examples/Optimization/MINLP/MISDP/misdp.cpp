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

int main(){
//string fname=string(prj_dir)+"/data_sets/MISDP/2x3_3bars.cbf";
   // string fname=string(prj_dir)+"/data_sets/MISDP/2x7_3bars.cbf";
    string fname=string(prj_dir)+"/data_sets/MISDP/2x4_2scen_3bars.cbf";
auto m=make_shared<Model<double>>("misdp_test");
auto g=CBF_read(fname.c_str(), m);
    g.print();
   m->print();
//    neg=true;
//    while(res.size()==0){
//        solver<> s(m,gurobi);
//        s.run(5, 1e-9);
//        auto res= m->cuts_eigen(1e-9);
//        DebugOn("cuts found "<<res.size()<<endl);
//    }
    
    
   solver<> s(m,gurobi);
   s.run(5, 1e-6);

   // auto lin_model=m->buildOA();
   // auto interior_model=lin_model->add_outer_app_solution(*m);

//    int soc_viol=0,soc_found=0,soc_added=0,det_viol=0,det_found=0,det_added=0;
//    auto res=m->cutting_planes_solution(interior_model, 1e-9,soc_viol, soc_found,soc_added,det_viol, det_found, det_added);
    
    
//    solver<> s(m, gurobi);
//    s.run(5,1e-9);
    m->print_solution();
    m->cuts_eigen(1e-10);
    m->print_constraints_stats(1e-9);
    

    
    var<double> X=m->get_var<double>("X");
    var<double> Xij=m->get_var<double>("Xij");
    int dim_full=X._indices->_keys->size();
    
    Eigen::MatrixXd mat_full(dim_full,dim_full);
    vector<vector<double>> mat_Xf(dim_full, std::vector<double>(dim_full, 0));
    int count=0;
    vector<string> all_names;
    for(auto k:*X._indices->_keys){
        mat_full(count, count)=X.eval(k);
        all_names.push_back(k);
        count++;
    }
    
    for(auto i=0;i<all_names.size()-1;i++){
        for(auto j=i+1;j<all_names.size();j++){
            auto k=all_names[i]+","+all_names[j];
            if (Xij._indices->has_key(k)){
                mat_full(i, j)=Xij.eval(k);
                mat_full(j, i)=Xij.eval(k);
            }
            else{
                mat_full(i, j)=0;
                mat_full(j, i)=0;
            }
        }
    }
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es1;
    es1.compute(mat_full);
    
    for(auto m=0;m<dim_full;m++){
        DebugOn(std::setprecision(12)<<es1.eigenvalues()[m]<<" ");
    }
    DebugOn(endl<<"full"<<endl);
    
      for(auto b:g._bags){
          auto dim=b.second.size();
          vector<string> node_names;
          vector<vector<double>> mat_X(dim, std::vector<double>(dim, 0));
          Eigen::MatrixXd mat(dim,dim);
          int count=0;
          for(auto n:b.second){
              node_names.push_back(n->_name);
              mat_X[count][count]=X.eval(n->_name);
              mat(count,count)=mat_X[count][count];
              count++;
              DebugOn(n->_name<<" ");
          }
          DebugOn(endl);
          for(auto i=0;i<node_names.size()-1;i++){
              for(auto j=i+1;j<node_names.size();j++){
                  mat_X[i][j]=Xij.eval(node_names[i]+","+node_names[j]);
                  mat_X[j][i]=Xij.eval(node_names[i]+","+node_names[j]);
                  mat(i,j)=mat_X[i][j];
                  mat(j,i)=mat_X[i][j];
              }
          }
         // double det=determinant(mat_X, dim);
          Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
          es.compute(mat);
          DebugOn("clique "<<dim<<endl);
          for(auto m=0;m<dim;m++){
              //                           cout<<es.eigenvalues()[m].real();
              DebugOn(std::setprecision(12)<<es.eigenvalues()[m]<<" ");
          }
          DebugOn(endl);
          //DebugOn("Determinant "<<std::setprecision(12)<<det<<" clique size "<<dim<<endl);
      }

//    double lower_bound=numeric_limits<double>::min(),upper_bound=numeric_limits<double>::max(), lower_bound_nonlin_init=numeric_limits<double>::min(),total_time=numeric_limits<double>::min();
//    double ub_solver_tol=1e-6, lb_solver_tol=1e-6, range_tol=1e-3, opt_rel_tol=1e-2, opt_abs_tol=1e-2, zero_tol=1e-12;
//    unsigned max_iter=1e3;
//        int oacuts=0, oacuts_init=0, fail=0;
//        SolverType ub_solver_type = ipopt, lb_solver_type = ipopt;
//        auto nonlin_obj=false;
//        std::vector<double> vrbasis;
//        std::map<string,double> crbasis;
//        lower_bound=0;
//        double active_root_tol=1e-6;
//        double lb_scale_value=1;
//        int nb_root_refine=10;
//        bool initialize_primal=false;
//    lin_model->print();
//
//    auto close=m->root_refine(interior_model, lin_model, lb_solver_type, nb_root_refine, upper_bound, lower_bound, lb_scale_value, lb_solver_tol, active_root_tol, oacuts,  opt_abs_tol, opt_rel_tol, zero_tol, "ma27", 10000, 2000, vrbasis, crbasis, initialize_primal);
//    DebugOn("oa cuts "<<oacuts<<endl);
//    lin_model->print();
}
