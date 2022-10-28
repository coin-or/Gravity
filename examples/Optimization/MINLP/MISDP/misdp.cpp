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
    string fname=string(prj_dir)+"/data_sets/MISDP/2x7_3bars.cbf";
    //string fname=string(prj_dir)+"/data_sets/MISDP/2x4_2scen_3bars.cbf";
    //string fname=string(prj_dir)+"/data_sets/MISDP/coloncancer_1_100_5.cbf";
    //string fname=string(prj_dir)+"/data_sets/MISDP/2g_4_164_k3_5_6.cbf";

    bool root_refine = true;
    if(argc>=2){
        fname=argv[1];
    }
    if(argc>=3){
        root_refine = true;
    }
    auto m=make_shared<Model<double>>("misdp_test");
    
    auto g=CBF_read(fname.c_str(), m);
    m->print();
    //    g.print();
    // m->print();
//    auto rel=make_shared<Model<double>>("misdp_rel");
//    auto g2=CBF_read(fname.c_str(), rel);
    double ub_solver_tol=1e-6, lb_solver_tol=1e-6, range_tol=1e-4, opt_rel_tol=1e-2, opt_abs_tol=1e6;
    unsigned max_iter=1e3;
    SolverType ub_solver_type = ipopt, lb_solver_type = gurobi;
    
    double max_time = 300;
    DebugOn("Instance "<<fname<<endl);
    //    rel->replace_integers();
    //    m->replace_integers();
    //    m->set_name("Before");
    //    m->write();
    //    auto res=m->run_obbt(rel, max_time, max_iter, opt_rel_tol, opt_abs_tol, 3, ub_solver_type, lb_solver_type, ub_solver_tol, lb_solver_tol, range_tol);
    //    m->set_name("After");
    //    m->write();
    ////    rel->print();
    //    int soc_viol=0, soc_added=0;
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
    bool upper_bound_heur=false;
    
    
    if(upper_bound_heur){
        vector<double> sol(m->_nb_vars);
        auto m2=m->copy();
        auto m1=m->copy();
        var<double> Xa=m2->get_var<double>("X");
        var<double> Xija=m2->get_var<double>("Xij");
        // m2->min(sum(Xija)-sum(Xa));
        auto o=25*(*m2->_obj)+sum(Xija)-sum(Xa);
        m2->min(o);
        Constraint<> diag_cut("diag_cut");
        diag_cut=sum(Xija)-sum(Xa);
        m2->add(diag_cut<=0);
        //m2->print();
        solver<> s(m2,ipopt);
        s.run();
        m2->print_solution();
        auto y2d=m2->get_var<double>("y");
        auto y1=m1->get_var<int>("y");
        
        for(auto k:*y2d.get_keys()){
            if(y2d.eval(k)<0.5){
                y1.set_lb(k,0);
                y1.set_ub(k,0);
                
            }
            else{
                y1.set_lb(k,1);
                y1.set_ub(k,1);
            }
        }
        solver<> s1(m1,gurobi);
        s1.run();
        //check_PSD(m1);
        if(m1->_status==0){
            m1->get_solution(sol);
        }
    }
    
    //m->set_solution(sol);
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
    sc.run(relax=false);
    auto tf=get_wall_time();
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
    
    auto eig_value=m->check_PSD();
    
    string out_file_name=fname;
    auto pos=out_file_name.find_last_of("/");
    out_file_name=out_file_name.substr(pos+1);
    pos=out_file_name.find_first_of(".");
    auto out_file_name1=out_file_name.substr(0,pos);
    out_file_name=string(prj_dir)+"/results_misdp/"+out_file_name1+".txt";
    ofstream fout(out_file_name.c_str());
    fout<<out_file_name1<<"&"<<m->get_obj_val()<<"&"<<tf-ts<<"&"<<eig_value<<"&"<<m->_rel_obj_val<<"&"<<m->num_cuts[0]<<"&"<<m->num_cuts[1]<<"&"<<m->num_cuts[2]<<"&"<<m->num_cuts[3]<<"&"<<m->num_cuts[4]<<"&"<<m->num_cuts[5]<<"\n";
    fout.close();
    cout<<out_file_name1<<" obj "<<m->get_obj_val()<<" time "<<tf-ts<<" smallest eig "<<eig_value<<" lower bound  "<<m->_rel_obj_val<<endl;
    cout<<"soc incumbent "<<m->num_cuts[0]<<" eig bags incumbent "<<m->num_cuts[1]<<" eig full incumbent "<<m->num_cuts[2]<<" soc mipnode "<<m->num_cuts[3]<<" eig bags mipnode "<<m->num_cuts[4]<<" eig full mipnode "<<m->num_cuts[5]<<"\n";
}
