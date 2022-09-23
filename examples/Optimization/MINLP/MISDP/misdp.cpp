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

int main(int argc, char * argv[]){
//string fname=string(prj_dir)+"/data_sets/MISDP/2x3_3bars.cbf";
//    string fname=string(prj_dir)+"/data_sets/MISDP/2x7_3bars.cbf";
    string fname=string(prj_dir)+"/data_sets/MISDP/2x4_2scen_3bars.cbf";
    //string fname=string(prj_dir)+"/data_sets/MISDP/coloncancer_1_100_5.cbf";
    //string fname=string(prj_dir)+"/data_sets/MISDP/2g_4_164_k3_5_6.cbf";
    bool root_refine = false;
    if(argc>=2){
        fname=argv[1];
    }
    if(argc>=3){
        root_refine = true;
    }
auto m=make_shared<Model<double>>("misdp_test");
auto g=CBF_read(fname.c_str(), m);
//    g.print();
  // m->print();
    auto rel=make_shared<Model<double>>("misdp_rel");
    auto g2=CBF_read(fname.c_str(), rel);
    double ub_solver_tol=1e-6, lb_solver_tol=1e-6, range_tol=1e-4, opt_rel_tol=1e-2, opt_abs_tol=1e6;
    unsigned max_iter=1e3;
    SolverType ub_solver_type = ipopt, lb_solver_type = gurobi;

    double max_time = 300;
    rel->replace_integers();
    m->replace_integers();
    m->set_name("Before");
    m->write();
    auto res=m->run_obbt(rel, max_time, max_iter, opt_rel_tol, opt_abs_tol, 3, ub_solver_type, lb_solver_type, ub_solver_tol, lb_solver_tol, range_tol);
    m->set_name("After");
    m->write();
//    rel->print();
    int soc_viol=0, soc_added=0;
    //m->print_solution();
    //m->cutting_planes_soc(1e-9, soc_viol,soc_added);
    
    //m->cuts_eigen_bags(1e-9);
//    neg=true;
//    while(res.size()==0){
//        solver<> s(m,gurobi);
//        s.run(5, 1e-9);
//        auto res= m->cuts_eigen(1e-9);
//        DebugOn("cuts found "<<res.size()<<endl);
//    }
    
//    auto m=make_shared<Model<double>>("misdp_test");
//    auto g=CBF_read(fname.c_str(), m);
//    m->copy_bounds(m_init);
//    m->copy_solution(m_init);
//    m->print();
//    auto y = m->get_var<int>("y");
//    y.set_ub(0);
//    y.set_lb("3", 1);
//    y.set_ub("3", 1);
//    y.set_lb("21", 1);
//    y.set_ub("21", 1);
//    y.set_lb("30", 1);
//    y.set_ub("30", 1);
//    y.set_lb("48", 1);
//    y.set_ub("48", 1);
//    y.set_lb("54", 1);
//    y.set_ub("54", 1);
//    y.set_lb("75", 1);
//    y.set_ub("75", 1);
//    y.set_lb("78", 1);
//    y.set_ub("78", 1);
//    y.set_lb("87", 1);
//    y.set_ub("87", 1);
//    y.set_lb("93", 1);
//    y.set_ub("93", 1);
//    y.set_lb("111", 1);
//    y.set_ub("111", 1);
//    y.set_lb("132", 1);
//    y.set_ub("132", 1);
//    y.set_lb("141", 1);
//    y.set_ub("141", 1);
//    y.set_lb("165", 1);
//    y.set_ub("165", 1);
//    auto x = m->get_var<double>("x");
//    x.set_lb(0.1);
    m->reset();
   solver<> sc(m,gurobi);
    bool relax = false;
    if(root_refine){
        sc.run(relax = true);
    }
    sc.run(relax=false);
    m->print_solution();
  /*  m->round_solution();
    m->_obj->uneval();
    DebugOn("Rounded solution objective = " << m->_obj->eval() << endl);*/
   // m->print_constraints_stats(1e-5);
    
   // auto lin_model=m->buildOA();
   // auto interior_model=lin_model->add_outer_app_solution(*m);

//    int soc_viol=0,soc_found=0,soc_added=0,det_viol=0,det_found=0,det_added=0;
//    auto res=m->cutting_planes_solution(interior_model, 1e-9,soc_viol, soc_found,soc_added,det_viol, det_found, det_added);
    
    
//    solver<> s(m, gurobi);
//    s.run(5,1e-9);
   // m->print_solution();
    //m->cuts_eigen(1e-10);
    m->print_constraints_stats(1e-9);
    

    
    var<double> X=m->get_var<double>("X");
    var<double> Xij=m->get_var<double>("Xij");
    int dim_full=X._indices->_keys->size();
    
    Eigen::MatrixXd mat_full(dim_full,dim_full);
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

}
