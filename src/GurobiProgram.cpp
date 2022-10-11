#include <gravity/GurobiProgram.h>
#include <thread>         // std::this_thread::sleep_for
#include <chrono>         // std::chrono::seconds
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
class cuts: public GRBCallback
{
public:
    vector<GRBVar> vars;
    int n;
    Model<>* m;
    Model<> interior;
    vector<GRBLinExpr> vec_expi;
    int soc_viol=0, soc_found=0,soc_added=0,det_viol=0, det_found=0, det_added=0;
    int soc_viol_user=0, soc_found_user=0,soc_added_user=0,det_viol_user=0, det_found_user=0, det_added_user=0;
    //cuts(vector<GRBVar> _grb_vars, int xn, Model<>* mn, Model<> interiorn, vector<GRBLinExpr>& vec_exp) {
    cuts(vector<GRBVar> _grb_vars, int xn, Model<>* mn, Model<> interiorn) {
        vars = _grb_vars;
        n    = xn;
        m=mn;
        interior=interiorn;
        //vec_expi=vec_exp;
        m->num_cuts.resize(6,0);
    }
    ~cuts(){
        DebugOn("soc_viol_user "<<soc_viol_user<<" "<<soc_viol<<endl);
        DebugOn("soc_found_user "<<soc_found_user<<" "<<soc_found<<endl);
        DebugOn("soc_added_user "<<soc_added_user<<" "<<soc_added<<endl);
        DebugOn("det_viol_user "<<det_viol_user<<" "<<det_viol<<endl);
        DebugOn("det_found_user "<<det_found_user<<" "<<det_found<<endl);
        DebugOn("det_added_user "<<det_added_user<<" "<<det_added<<endl);
    }
protected:
    void callback() {
        try {
            bool incumbent=true;
            bool mipnode=true;
            bool hierarc = true;
            if(incumbent){
                if (where == GRB_CB_MIPSOL) {
                    /* Found an integer feasible solution */
                    double *x = new double[n];
                    vector<double> vec_x;
                    // double obj=getDoubleInfo(GRB_CB_MIPSOL_OBJ);
                    int i,j;
                    x=getSolution(vars.data(),n);
                    for(i=0;i<n;i++){
                        vec_x.push_back(x[i]);
                    }
                    m->set_solution(vec_x);
                    //                    auto res=m->cutting_planes_solution(interior, 1e-9, soc_viol, soc_found, soc_added,  det_viol, det_found, det_added);
                    //                    if(res.size()>=1){
                    //                        for(auto i=0;i<res.size();i++){
                    //                            GRBLinExpr expr = 0;
                    //                            int j=0;
                    //                            for(j=0;j<res[i].size()-1;j+=2){
                    //                                int c=res[i][j];
                    //                                DebugOff(res[i][j]<<" ");
                    //                                expr += res[i][j+1]*vars[c];
                    //                                DebugOff(std::setprecision(12)<<res[i][j+1]<<" ");
                    //                            }
                    //                            DebugOff(endl);
                    //                            if(abs(res[i][j])>=1e-6){
                    //                                DebugOff("pos resij");
                    //                            }
                    //                            expr+=res[i][j];
                    //
                    //                            vec_expi.push_back(expr);
                    //                           addLazy(expr, GRB_LESS_EQUAL, 0);
                    //                            DebugOff("added soc"<<endl);
                    //                        }
                    //                    }
                    //                    m->set_solution(vec_x);
                    if(  true ){
                        
                        //                        auto res0= m->cutting_planes_solution(interior, 1e-9, soc_viol,soc_found,soc_added,  det_viol, det_found,det_added);
                        //                        //auto res0= m->cutting_planes_soc(1e-9, soc_found, soc_added);
                        //                          if(res0.size()>=1){
                        //                              for(auto i=0;i<res0.size();i++){
                        //                                  GRBLinExpr expr = 0;
                        //                                  int j;
                        //                                  for(j=0;j<res0[i].size()-1;j+=2){
                        //                                      int c=res0[i][j];
                        //                                      expr += res0[i][j+1]*vars[c];
                        //                                  }
                        //                                  addLazy(expr, GRB_LESS_EQUAL, 0);
                        //                                  //vec_expi.push_back(expr);
                        //                              }
                        //                          }
                        //                        else{
                        m->set_solution(vec_x);
                        auto res= m->cutting_planes_soc(1e-9, soc_found, soc_added);
                        if(res.size()>=1){
                            for(auto i=0;i<res.size();i++){
                                GRBLinExpr expr = 0;
                                int j;
                                DebugOff("soc cut at ");
                                for(j=0;j<res[i].size()-1;j+=2){
                                    int c=res[i][j];
                                    if(std::abs(c)>=0)
                                        expr += res[i][j+1]*vars[c];
                                    DebugOff(to_string_with_precision(vec_x[c],10)<<" "<<to_string_with_precision(c,10)<<" "<<to_string_with_precision(res[i][j+1],10)<<" ");
                                }
                                expr += res[i][j];
                                DebugOff(endl);
                                addLazy(expr, GRB_LESS_EQUAL, 0);
                                //                                  addCut(expr, GRB_LESS_EQUAL, 0);
                                
                                //vec_expi.push_back(expr);
                                m->num_cuts[0]++;
                            }
                        }
                        
                        
                        if(res.size()==0 || !hierarc){
                            m->set_solution(vec_x);
                            //                            if(true || res.size()==0){
                            auto res1=m->cuts_eigen_bags(1e-9);
                            if(res1.size()>=1){
                                for(auto i=0;i<res1.size();i++){
                                    GRBLinExpr expr = 0;
                                    int j;
                                    DebugOff("eig cut at");
                                    for(j=0;j<res1[i].size()-1;j+=2){
                                        int c=res1[i][j];
                                        expr += res1[i][j+1]*vars[c];
                                        DebugOff(to_string_with_precision(vec_x[c],10)<<" "<<to_string_with_precision(c,10)<<" "<<to_string_with_precision(res1[i][j+1],10)<<" ");
                                    }
                                    DebugOff(to_string_with_precision(res1[i][j],10));
                                    DebugOff(endl);
                                    if(std::abs(res1[i][j])>=1e-12)
                                        expr += res1[i][j];
                                    addLazy(expr, GRB_LESS_EQUAL, 0);
                                    m->num_cuts[1]++;
                                    //                                addCut(expr, GRB_LESS_EQUAL, 0);
                                    //vec_expi.push_back(expr);
                                }
                            }
                            
                            
                            m->set_solution(vec_x);
                            
                            if(res1.size()==0 || !hierarc){
                                auto res2=m->cuts_eigen_full(1e-9);
                                if(res2.size()>=1){
                                    for(auto i=0;i<res2.size();i++){
                                        GRBLinExpr expr = 0;
                                        int j;
                                        for(j=0;j<res2[i].size()-1;j+=2){
                                            int c=res2[i][j];
                                            expr += res2[i][j+1]*vars[c];
                                        }
                                        if(std::abs(res2[i][j])>=1e-12)
                                            expr += res2[i][j];
                                        addLazy(expr, GRB_LESS_EQUAL, 0);
                                        m->num_cuts[2]++;
                                        //                            addCut(expr, GRB_LESS_EQUAL, 0);
                                        //vec_expi.push_back(expr);
                                    }
                                }
                            }
                        }
                        
                        //auto res=m->cutting_planes_solution(interior, 1e-6);
                        /*auto res=m->cutting_planes_solution(interior, 1e-9, soc_viol, soc_found, soc_added, det_viol, det_found, det_added);
                         if(res.size()>=1){
                         for(i=0;i<res.size();i++){
                         GRBLinExpr expr = 0;
                         for(j=0;j<res[i].size()-1;j+=2){
                         int c=res[i][j];
                         expr += res[i][j+1]*vars[c];
                         }
                         if(abs(res[i][j])>=1e-6){
                         DebugOff("pos resij");
                         }
                         expr+=res[i][j];
                         // vec_expi.push_back(expr);
                         //addLazy(expr, GRB_LESS_EQUAL, 0);
                         }
                         }
                         auto res1=m->cutting_planes_eigen(1e-9);
                         if(res1.size()>=1){
                         for(i=0;i<res1.size();i++){
                         GRBLinExpr expr = 0;
                         for(j=0;j<res1[i].size()-1;j+=2){
                         int c=res1[i][j];
                         expr += res1[i][j+1]*vars[c];
                         }
                         if(abs(res1[i][j])>=1e-6){
                         DebugOff("pos resij");
                         }
                         expr+=res1[i][j];
                         //addLazy(expr, GRB_LESS_EQUAL, 0);
                         }
                         }
                         m->set_solution(vec_x);
                         auto res2=m->cuts_eigen(1e-10);
                         if(res2.size()>=1){
                         for(i=0;i<res2.size();i++){
                         GRBLinExpr expr = 0;
                         for(j=0;j<res2[i].size()-1;j+=2){
                         int c=res2[i][j];
                         expr += res2[i][j+1]*vars[c];
                         }
                         if(abs(res2[i][j])>=1e-6){
                         DebugOff("pos resij");
                         }
                         expr+=res2[i][j];
                         addLazy(expr, GRB_LESS_EQUAL, 0);
                         vec_expi.push_back(expr);
                         }*/
                        // }
                        
                        
                    }
                }
            }
            if(true){
                if (where == GRB_CB_MIPNODE) {
                    if(getIntInfo(GRB_CB_MIPNODE_STATUS)==2 ){
                        int nc= getDoubleInfo(GRB_CB_MIPNODE_NODCNT);
                        //                        if(nc%100==0){
                        // Found an integer feasible solution - does it visit every node?
                        double *x = new double[n];
                        vector<double> vec_x;
                        // double obj=getDoubleInfo(GRB_CB_MIPSOL_OBJ);
                        int i,j;
                        x=getNodeRel(vars.data(),n);
                        //x=getSolution(vars.data(),n);
                        for(i=0;i<n;i++){
                            vec_x.push_back(x[i]);
                        }
                        m->set_solution(vec_x);
                        
                        auto res= m->cutting_planes_soc(1e-9, soc_found, soc_added);
                        if(res.size()>=1){
                            for(auto i=0;i<res.size();i++){
                                GRBLinExpr expr = 0;
                                int j;
                                DebugOff("soc cut at ");
                                for(j=0;j<res[i].size()-1;j+=2){
                                    int c=res[i][j];
                                    if(std::abs(c)>=0)
                                        expr += res[i][j+1]*vars[c];
                                    DebugOff(to_string_with_precision(vec_x[c],10)<<" "<<to_string_with_precision(c,10)<<" "<<to_string_with_precision(res[i][j+1],10)<<" ");
                                }
                                expr += res[i][j];
                                DebugOff(endl);
                                addLazy(expr, GRB_LESS_EQUAL, 0);
//                                addCut(expr, GRB_LESS_EQUAL, 0);
                                m->num_cuts[3]++;
                                //vec_expi.push_back(expr);
                            }
                        }
                        if(res.size()==0 || !hierarc){
                            m->set_solution(vec_x);
                            auto res1=m->cuts_eigen_bags(1e-9);
                            if(res1.size()>=1){
                                for(auto i=0;i<res1.size();i++){
                                    GRBLinExpr expr = 0;
                                    int j;
                                    DebugOff("eig cut at");
                                    for(j=0;j<res1[i].size()-1;j+=2){
                                        int c=res1[i][j];
                                        expr += res1[i][j+1]*vars[c];
                                        DebugOff(to_string_with_precision(vec_x[c],10)<<" "<<to_string_with_precision(c,10)<<" "<<to_string_with_precision(res1[i][j+1],10)<<" ");
                                    }
                                    DebugOff(to_string_with_precision(res1[i][j],10));
                                    DebugOff(endl);
                                    if(std::abs(res1[i][j])>=1e-12)
                                        expr += res1[i][j];
                                    addLazy(expr, GRB_LESS_EQUAL, 0);
                                    m->num_cuts[4]++;
//                                    addCut(expr, GRB_LESS_EQUAL, 0);
                                    //vec_expi.push_back(expr);
                                }
                            }
                            
                            
                            if(res1.size()==0 || !hierarc){
                                m->set_solution(vec_x);
                                auto res2=m->cuts_eigen_full(1e-9);
                                if(res2.size()>=1){
                                    for(auto i=0;i<res2.size();i++){
                                        GRBLinExpr expr = 0;
                                        int j;
                                        for(j=0;j<res2[i].size()-1;j+=2){
                                            int c=res2[i][j];
                                            expr += res2[i][j+1]*vars[c];
                                        }
                                        if(std::abs(res2[i][j])>=1e-12)
                                            expr += res2[i][j];
                                        addLazy(expr, GRB_LESS_EQUAL, 0);
                                        m->num_cuts[5]++;
                                        //                                    addCut(expr, GRB_LESS_EQUAL, 0);
                                        //vec_expi.push_back(expr);
                                    }
                                }
                            }
                        }
                    }
                }
                //auto res=m->cutting_planes_solution(interior, 1e-6);
                /* auto res=m->cutting_planes_solution(interior, 1e-9,soc_viol_user, soc_found_user,soc_added_user,det_viol_user, det_found_user, det_added_user);
                 if(res.size()>=1){
                 for(i=0;i<res.size();i++){
                 GRBLinExpr expr = 0;
                 for(j=0;j<res[i].size()-1;j+=2){
                 int c=res[i][j];
                 expr += res[i][j+1]*vars[c];
                 }
                 if(abs(res[i][j])>=1e-6){
                 DebugOff("pos resij");
                 }
                 expr+=res[i][j];
                 addLazy(expr, GRB_LESS_EQUAL, 0);
                 }
                 }*/
                //                    }}
            }
            
        } catch (GRBException e) {
            cout << "Error number: " << e.getErrorCode() << endl;
            cout << e.getMessage() << endl;
        } catch (...) {
            cout << "Error during callback" << endl;
        }
    }
};

GurobiProgram::GurobiProgram(){
    //    model = m;
    grb_env = new GRBEnv();
    //    grb_env->set(GRB_IntParam_Presolve,0);
    //grb_env->set(GRB_DoubleParam_NodeLimit,1);
//    grb_env->set(GRB_DoubleParam_TimeLimit,9000);
    //    grb_env->set(GRB_DoubleParam_MIPGap,0.01);
    //   grb_env->set(GRB_IntParam_Threads,8);
    //    grb_env->set(GRB_IntParam_Presolve,0);
    //   grb_env->set(GRB_IntParam_NumericFocus,3);
//    grb_env->set(GRB_IntParam_NonConvex,0);
//    grb_env->set(GRB_DoubleParam_FeasibilityTol, 1e-9);
//    grb_env->set(GRB_DoubleParam_OptimalityTol, 1e-9);
    
    
    grb_env->set(GRB_IntParam_OutputFlag,1);
    //    grb_mod = new GRBModel(*grb_env);
    grb_mod = NULL;
}


GurobiProgram::GurobiProgram(Model<>* m) {
    bool found_token = false;
    while (!found_token) {
        try{
            grb_env = new GRBEnv();
            //    grb_env->set(GRB_IntParam_Presolve,0);
            //grb_env->set(GRB_DoubleParam_NodeLimit,1);
//            grb_env->set(GRB_DoubleParam_TimeLimit,9000);
            //   grb_env->set(GRB_IntParam_Threads,8);
            //    grb_env->set(GRB_DoubleParam_MIPGap,0.01);
            //    grb_env->set(GRB_IntParam_Threads,1);
            //    grb_env->set(GRB_IntParam_Presolve,0);
            //     grb_env->set(GRB_IntParam_NumericFocus,3);
           // grb_env->set(GRB_IntParam_NonConvex,0);
           // grb_env->set(GRB_DoubleParam_FeasibilityTol, 1e-9);
           // grb_env->set(GRB_DoubleParam_OptimalityTol, 1e-9);
            
            //grb_env->set(GRB_IntParam_OutputFlag,1);
            grb_mod = new GRBModel(*grb_env);
            //    grb_env->set(GRB_IntParam_OutputFlag,2);
            found_token = true;
        }
        catch(GRBException e) {
            //            cerr << "\nWas not able to create Gurobi environment or model, Error code = " << e.getErrorCode() << endl;
            //            cerr << e.getMessage() << endl;
            //            exit(-1);
            found_token = false;
            this_thread::sleep_for (chrono::seconds(1));
        }
    }
    _model = m;
    m->fill_in_maps();
    m->compute_funcs();
}

GurobiProgram::GurobiProgram(const shared_ptr<Model<>>& m) {
    bool found_token = false;
    while (!found_token) {
        try{
            grb_env = new GRBEnv();
            //    grb_env->set(GRB_IntParam_Presolve,0);
            //grb_env->set(GRB_DoubleParam_NodeLimit,1);
            //   grb_env->set(GRB_IntParam_Threads,8);
//            grb_env->set(GRB_DoubleParam_TimeLimit,9000);
            //    grb_env->set(GRB_DoubleParam_MIPGap,0.01);
            //    grb_env->set(GRB_IntParam_Threads,1);
            //    grb_env->set(GRB_IntParam_Presolve,0);
            //   grb_env->set(GRB_IntParam_NumericFocus,3);
            //grb_env->set(GRB_IntParam_NonConvex,0);
            //grb_env->set(GRB_DoubleParam_FeasibilityTol, 1e-9);
            //grb_env->set(GRB_DoubleParam_OptimalityTol, 1e-9);
            
            // grb_env->set(GRB_IntParam_OutputFlag,1);
            grb_mod = new GRBModel(*grb_env);
            grb_mod->set(GRB_IntParam_LazyConstraints, 1);
            
            
            found_token = true;
        }
        catch(GRBException e) {
            //            cerr << "\nWas not able to create Gurobi environment or model, Error code = " << e.getErrorCode() << endl;
            //            cerr << e.getMessage() << endl;
            //            exit(-1);
            found_token = false;
            this_thread::sleep_for (chrono::seconds(1));
        }
    }
    //    grb_env->set(GRB_IntParam_OutputFlag,2);
    //    _model = m;
    m->fill_in_maps();
    m->compute_funcs();
}


GurobiProgram::~GurobiProgram() {
    //    for (auto p : _grb_vars) delete p.second;
    if (grb_mod) delete grb_mod;
    delete grb_env;
}

void GurobiProgram::reset_model(){
    if (grb_mod != NULL) delete grb_mod;
    _grb_vars.clear();
    //    grb_env->set(GRB_IntParam_OutputFlag,_output);
    grb_mod = new GRBModel(*grb_env);
}

bool GurobiProgram::solve(bool relax, double mipgap, double time_limit){
    //cout << "\n Presolve = " << grb_env->get(GRB_IntParam_Presolve) << endl;
    //    print_constraints();
        if (relax) relax_model();
        else
            unrelax_model();
    //    relax_model();
    grb_mod->set(GRB_DoubleParam_MIPGap, 1e-6);
    grb_mod->set(GRB_DoubleParam_FeasibilityTol, 1e-9);
    grb_mod->set(GRB_DoubleParam_OptimalityTol, 1e-6);
    // grb_mod->set(GRB_IntParam_StartNodeLimit, -3);
    //    grb_mod->getEnv().set(GRB_IntParam_DualReductions, 0);
    //    grb_mod->getEnv().set(GRB_IntParam_PreCrush, 1);
    //        grb_mod->getEnv().set(GRB_IntParam_Method, 1);
    //    grb_mod->getEnv().set(GRB_IntParam_NodeMethod, 1);
    grb_mod->getEnv().set(GRB_IntParam_LazyConstraints, 1);
    grb_mod->set(GRB_IntParam_Threads, 8);
        grb_mod->set(GRB_DoubleParam_IntFeasTol, 1e-9);
    //   grb_mod->set(GRB_IntParam_NumericFocus,3);
   //  grb_mod->set(GRB_IntParam_Presolve,2);
    grb_mod->set(GRB_IntParam_MIPFocus,3);
    //    grb_mod->set(GRB_IntParam_IntegralityFocus,1);
    //    grb_mod->set(GRB_IntParam_MIPFocus,2);
    //    grb_mod->set(GRB_IntParam_PumpPasses,50);
    //    grb_mod->set(GRB_IntParam_RINS,1000);
    //    grb_mod->set(GRB_IntParam_Cuts,0);
    
    grb_mod->set(GRB_DoubleParam_TimeLimit,14300);
    grb_mod->set(GRB_IntParam_OutputFlag,1);
    //grb_mod->set(GRB_DoubleParam_Cutoff,5.33);
    //  grb_mod->set(GRB_IntParam_MinRelNodes,0);
    //    grb_mod->set(GRB_DoubleParam_Heuristics, 1);
    //    grb_mod->set(GRB_DoubleParam_NoRelHeurTime, 5);
    //    grb_mod->set(GRB_DoubleParam_NoRelHeurWork, 5);
    //    grb_mod->set(GRB_IntParam_CutPasses,10000);
    grb_mod->update();
    //grb_env2 = new GRBEnv();
    //auto mod2=GRBModel(grb_mod);
    int n=grb_mod->get(GRB_IntAttr_NumVars);
    Model<> interior;
    //    auto lin=_model->buildOA();
    //    interior=lin->get_interior(*_model);
    
    //interior=lin->add_outer_app_solution(*_model);
    //interior.print_solution();
    //cuts cb = cuts(_grb_vars, n, _model, interior);
    //vector<GRBLinExpr> vec_expr;
        cuts cb(_grb_vars, n, _model, interior);
        grb_mod->setCallback(&cb);
        grb_mod->update();
    
    //grb_mod->set(GRB_IntParam_RINS,1);
    // grb_mod->set(GRB_DoubleParam_Heuristics, 0.5);
    //grb_mod->update();
    grb_mod->optimize();
//    grb_mod->write("out.sol");
//    _model->print_solution();
    //    grb_mod->set(GRB_IntParam_MIPFocus,3);
    //   // grb_mod->set(GRB_DoubleParam_ImproveStartTime,9000);
    //   // grb_mod->set(GRB_DoubleParam_Heuristics, 0.05);
    //    grb_mod->set(GRB_DoubleParam_TimeLimit,9000);
    //    grb_mod->set(GRB_IntParam_RINS,10000);
    // grb_mod->update();
    // grb_mod->optimize();
    if(grb_mod->get(GRB_IntAttr_SolCount)>0)
        update_solution();
    bool not_sdp=false;
    if(_model->check_PSD()<=-1e-9  && grb_mod->get(GRB_IntAttr_Status)==2){
        not_sdp=true;
        grb_mod->setCallback(NULL);
    }
       // grb_mod->set(GRB_IntParam_OutputFlag,0);
    
    int count=0;
    auto ts=get_wall_time();
//    double sol_new=grb_mod->get(GRB_DoubleAttr_ObjVal);
//    double sol_old=-1e6;
    while(not_sdp && count<=1000){
//        var<double> X=_model->get_var<double>("X");
//        var<double> Xij=_model->get_var<double>("Xij");
//        int dim_full=X._indices->_keys->size();
//
//        Eigen::MatrixXd mat_full(dim_full,dim_full);
//        int counti=0;
//        vector<string> all_names;
//        for(auto k:*X._indices->_keys){
//            mat_full(count, count)=X.eval(k);
//            all_names.push_back(k);
//            counti++;
//        }
//
//        for(auto i=0;i<all_names.size()-1;i++){
//            for(auto j=i+1;j<all_names.size();j++){
//                auto k=all_names[i]+","+all_names[j];
//                if (Xij._indices->has_key(k)){
//                    mat_full(i, j)=Xij.eval(k);
//                    mat_full(j, i)=Xij.eval(k);
//                }
//                else{
//                    mat_full(i, j)=0;
//                    mat_full(j, i)=0;
//                }
//            }
//        }
//        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es1;
//        es1.compute(mat_full);
        
        //        for(auto m=0;m<dim_full;m++){
        //            DebugOn(std::setprecision(12)<<es1.eigenvalues()[m]<<" ");
        //        }
        //        DebugOn(endl<<"full"<<endl);
        
        //grb_mod->setCallback(NULL);
        int soc_viol=0, soc_added=0;
        auto res=_model->cutting_planes_soc(1e-6,soc_viol, soc_added);
        if(res.size()>=1){
            for(auto i=0;i<res.size();i++){
                GRBLinExpr expr = 0;
                int j=0;
                for(j=0;j<res[i].size()-1;j+=2){
                    int c=res[i][j];
                    expr += res[i][j+1]*_grb_vars[c];
                    DebugOff(res[i][j+1]<<" ");
                }
                 expr+=res[i][j];
                grb_mod->addConstr(expr, GRB_LESS_EQUAL, 0);
            }
        }
        if(res.size()==0){
            auto res2=_model->cuts_eigen_bags(1e-6);
            if(res2.size()>=1){
                for(auto i=0;i<res2.size();i++){
                    GRBLinExpr expr = 0;
                    int j;
                    for(j=0;j<res2[i].size()-1;j+=2){
                        int c=res2[i][j];
                        if(std::abs(c)!=0)
                            expr += res2[i][j+1]*_grb_vars[c];
                    }
                    expr += res2[i][j];
                    grb_mod->addConstr(expr, GRB_LESS_EQUAL, 0);
                }
            }
            if(res2.size()==0){
            auto res3=_model->cuts_eigen_full(1e-6);
            if(res3.size()>=1){
                for(auto i=0;i<res3.size();i++){
                    GRBLinExpr expr = 0;
                    int j;
                    for(j=0;j<res3[i].size()-1;j+=2){
                        int c=res3[i][j];
                        expr += res3[i][j+1]*_grb_vars[c];
                    }
                        expr += res3[i][j];
                    grb_mod->addConstr(expr, GRB_LESS_EQUAL, 0);
                }
            }
            if(res3.size()==0){
                not_sdp=false;
            }
            }
        }
            if(not_sdp){
                grb_mod->update();
                grb_mod->optimize();
                if(grb_mod->get(GRB_IntAttr_Status)!=2)
                    break;
                if(grb_mod->get(GRB_IntAttr_SolCount)>0)
                    update_solution();
//                sol_old=sol_new;
//                sol_new=grb_mod->get(GRB_DoubleAttr_ObjVal);
            }

        count++;
    }
    auto tf=get_wall_time();
    DebugOn("While loop Count "<<count<<" in time "<<tf-ts<<endl);
    // grb_mod->update();
    // grb_mod->optimize();
    //        for(auto i=0;i<cb.vec_expi.size();i++){
    //            grb_mod->addConstr(cb.vec_expi[i],GRB_LESS_EQUAL, 0);
    //        }
    //    grb_mod->setCallback(NULL);
    //    grb_mod->update();
    //    grb_mod->optimize();
    //    grb_mod->write("mod.lp");
   
    //    GRBVar* gvars = grb_mod->getVars();
    //    for(int i = 0; i < grb_mod->get(GRB_IntAttr_NumVars); ++i) {
    ////        cout << gvars[i].get(GRB_StringAttr_VarName) << "  " << gvars[i].get(GRB_DoubleAttr_X) << endl;
    //        if (gvars[i].get(GRB_CharAttr_VType)==GRB_BINARY) {
    //            cout << gvars[i].get(GRB_StringAttr_VarName) << "  ";
    //            cout << gvars[i].get(GRB_DoubleAttr_X);
    //            cout << "\n";
    //        }
    //    }
    _model->_obj->set_val(grb_mod->get(GRB_DoubleAttr_ObjVal));
    _model->_rel_obj_val=grb_mod->get(GRB_DoubleAttr_ObjBound);
    cout << "\n***** Optimal Objective = " << _model->get_obj_val() << " *****\n";
    if (grb_mod->get(GRB_IntAttr_IsMIP)) {
        cout.setf(ios::fixed);
        cout.precision(3);
        cout << "Results: " << grb_mod->get(GRB_DoubleAttr_ObjVal) << " & ";
        cout.precision(4);
        cout << (grb_mod->get(GRB_DoubleAttr_MIPGap))*100 << "% & ";
        cout.precision(0);
        cout << grb_mod->get(GRB_DoubleAttr_NodeCount) << " & ";
        cout.precision(2);
        cout << grb_mod->get(GRB_DoubleAttr_Runtime) << " & " << endl;
    }
    //    delete[] gvars;
    if (grb_mod->get(GRB_IntAttr_Status) != 2) {
        cerr << "\nModel has not been solved to optimality, error code = " << grb_mod->get(GRB_IntAttr_Status) << endl;
               return false;
    }
    return true;
}

void GurobiProgram::prepare_model(){
    _model->fill_in_maps();
    _model->compute_funcs();
    fill_in_grb_vmap();
    create_grb_constraints();
    set_grb_objective();
    grb_mod->write("gurobiprint.lp");
    //    print_constraints();
}
void GurobiProgram::update_model(){
    _model->fill_in_maps();
    _model->compute_funcs();
    fill_in_grb_vmap();
    create_grb_constraints();
    set_grb_objective();
    
}

void GurobiProgram::update_solution(){
    size_t vid, vid_inst;
    GRBVar gvar;
    param_* v;
    for(auto& v_p: _model->_vars)
    {
        v = v_p.second.get();
        auto idx = v->get_id();
        auto dim = v->_dim[0];
        for (auto i = 0; i < dim; i++) {
            auto vid = idx + v->get_id_inst(i);
            gvar = _grb_vars.at(vid);
            if(v->is_integer()||v->is_binary()||v->_is_relaxed)
                v->set_double_val(i,std::round(gvar.get(GRB_DoubleAttr_X)));
            else
                v->set_double_val(i,gvar.get(GRB_DoubleAttr_X));
        }
    }
}

void GurobiProgram::relax_model(){
    GRBVar* gvars = grb_mod->getVars();
    for(int i = 0; i < grb_mod->get(GRB_IntAttr_NumVars); ++i) {
        if (gvars[i].get(GRB_CharAttr_VType) == 'B' || gvars[i].get(GRB_CharAttr_VType) == 'I')
            gvars[i].set(GRB_CharAttr_VType,'C');
    }
}

void GurobiProgram::unrelax_model(){
    GRBVar gvar;
    param_* v;
    for(auto& v_p: _model->_vars)
    {
        v = v_p.second.get();
        if(v->is_continuous())
            continue;
        auto idx = v->get_id();
        auto dim = v->_dim[0];
        for (auto i = 0; i < dim; i++) {
            auto vid = idx + v->get_id_inst(i);
            gvar = _grb_vars.at(vid);
            gvar.set(GRB_CharAttr_VType, 'I');
        }
    }
}

void GurobiProgram::fill_in_grb_vmap(){
    param_* v;
    _grb_vars.resize(_model->get_nb_vars());
    for(auto& v_p: _model->_vars)
    {
        v = v_p.second.get();
        if (!v->_new) {
            continue;
        }
        v->_new = false;
        auto idx = v->get_id();
        switch (v->get_intype()) {
            case float_: {
                auto real_var = (var<float>*)v;
                for (int i = 0; i < real_var->_dim[0]; i++) {
                    auto vid = idx + i;
                    _grb_vars.at(vid) = (GRBVar(grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_CONTINUOUS, v->get_name(true,true)+"("+v->_indices->_keys->at(i)+")")));
                }
                break;
            }
            case long_:{
                auto real_var = (var<long double>*)v;
                for (int i = 0; i < real_var->_dim[0]; i++) {
                    auto vid = idx + i;
                    _grb_vars.at(vid) = (GRBVar(grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_CONTINUOUS, v->get_name(true,true)+"("+v->_indices->_keys->at(i)+")")));
                }
                break;
            }
            case double_:{
                auto real_var = (var<double>*)v;
                for (size_t i = 0; i < real_var->_dim[0]; i++) {
                    auto vid = idx + i;
                    if(real_var->_is_relaxed){
                        _grb_vars.at(vid) = (GRBVar(grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_INTEGER, v->get_name(true,true)+"("+v->_indices->_keys->at(i)+")")));
                    }
                    else {
                        _grb_vars.at(vid) = (GRBVar(grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_CONTINUOUS, v->get_name(true,true)+"("+v->_indices->_keys->at(i)+")")));
                        _grb_vars.at(vid).set(GRB_DoubleAttr_Start, real_var->eval(i));
                    }
                }
                break;
            }
            case integer_:{
                auto real_var = (var<int>*)v;
                for (int i = 0; i < real_var->_dim[0]; i++) {
                    auto vid = idx + i;
                    _grb_vars.at(vid) = (GRBVar(grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_INTEGER, v->get_name(true,true)+"("+v->_indices->_keys->at(i)+")")));
                    _grb_vars.at(vid).set(GRB_DoubleAttr_Start, real_var->eval(i));
                }
                break;
            }
            case short_:{
                auto real_var = (var<short>*)v;
                for (int i = 0; i < real_var->_dim[0]; i++) {
                    auto vid = idx + i;
                    _grb_vars.at(vid) = (GRBVar(grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_INTEGER, v->get_name(true,true)+"("+v->_indices->_keys->at(i)+")")));
                }
                break;
            }
            case binary_:{
                auto real_var = (var<bool>*)v;
                for (int i = 0; i < real_var->_dim[0]; i++) {
                    auto vid = idx + i;
                    _grb_vars.at(vid) = (GRBVar(grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_BINARY, v->get_name(true,true)+"("+v->_indices->_keys->at(i)+")")));
                }
                break;
            }
            default:
                break;
        }
        
    }
    //    for(auto& v_p: _model->_vars)
    //    {
    //        v = v_p.second;
    //        auto real_var = (var<double>*)v;
    //        for (int i = 0; i < real_var->_dim[0]; i++) {
    //            auto vid = v->get_id() + v->get_id_inst(i);
    //            DebugOn("VID = "<< vid <<" : " << _grb_vars.at(vid).get(GRB_StringAttr_VarName) << " in [" << _grb_vars.at(vid).get(GRB_DoubleAttr_LB) << "," << _grb_vars.at(vid).get(GRB_DoubleAttr_UB) << "]\n" );
    //        }
    //    }
}

void GurobiProgram::create_grb_constraints(){
    char sense;
    size_t idx = 0, idx_inst = 0, idx1 = 0, idx2 = 0, idx_inst1 = 0, idx_inst2 = 0, nb_inst = 0, inst = 0;
    GRBLinExpr lterm, linlhs;
    GRBQuadExpr quadlhs;
    GRBVar gvar1, gvar2;
    double coeff;
    for(auto& p: _model->_cons){
        auto c = p.second;
        if (!c->_new && c->_all_satisfied) {
            continue;
        }
        if(c->_callback)
        {
            DebugOn(c->_name<<"  added as callback"<<endl);
            continue;
        }
        c->_new = false;
        
        if (c->is_nonlinear() && (!(c->_expr->is_uexpr() && c->get_nb_vars()==2) && !c->_expr->is_mexpr())) {
            throw invalid_argument("Gurobi cannot handle nonlinear constraints with more than two variables, try decomposing your constraints by introducing auxiliary variables.\n");
        }
        switch(c->get_ctype()) {
            case geq:
                sense = GRB_GREATER_EQUAL;
                break;
            case leq:
                sense = GRB_LESS_EQUAL;
                break;
            case eq:
                sense = GRB_EQUAL;
                break;
            default:
                break;
        }
        nb_inst = c->get_nb_inst();
        inst = 0;
        if (c->is_linear()) {
            for (size_t i = 0; i< nb_inst; i++){
                //    if (!c->_lazy[i]) {
                linlhs = 0;
                for (auto& it1: c->get_lterms()) {
                    lterm = 0;
                    if (it1.second._p->_is_vector || it1.second._p->is_matrix_indexed() || it1.second._coef->is_matrix()) {
                        auto dim =it1.second._p->get_dim(i);
                        for (int j = 0; j<dim; j++) {
                            coeff = c->eval(it1.second._coef,i,j);
                            gvar1 = _grb_vars[it1.second._p->get_id() + it1.second._p->get_id_inst(i,j)];
                            lterm += coeff*gvar1;
                        }
                    }
                    else {
                        coeff = c->eval(it1.second._coef,i);
                        gvar1 = _grb_vars[it1.second._p->get_id() + it1.second._p->get_id_inst(i)];
                        lterm += coeff*gvar1;
                    }
                    if (!it1.second._sign) {
                        lterm *= -1;
                    }
                    linlhs += lterm;
                }
                linlhs += c->eval(c->get_cst(), i);
                if(c->_indices)
                    grb_mod->addConstr(linlhs,sense,0,c->get_name()+"("+c->_indices->_keys->at(i)+")");
                else
                    grb_mod->addConstr(linlhs,sense,0,c->get_name());
            }
            
            // }
        }
        else if(c->is_quadratic()){
            for (size_t i = 0; i< nb_inst; i++){
                //    if (!c->_lazy[i]) {
                quadlhs = 0;
                for (auto& it1: c->get_lterms()) {
                    lterm = 0;
                    if (it1.second._coef->_is_transposed || it1.second._coef->is_matrix() || it1.second._p->is_matrix_indexed()) {
                        auto dim = it1.second._p->get_dim(i);
                        for (size_t j = 0; j<dim; j++) {
                            coeff = c->eval(it1.second._coef,i,j);
                            gvar1 = _grb_vars[it1.second._p->get_id() + it1.second._p->get_id_inst(i,j)];
                            lterm += coeff*gvar1;
                        }
                    }
                    else {
                        coeff = c->eval(it1.second._coef,i);
                        gvar1 = _grb_vars[it1.second._p->get_id() + it1.second._p->get_id_inst(i)];
                        lterm += coeff*gvar1;
                    }
                    if (!it1.second._sign) {
                        lterm *= -1;
                    }
                    quadlhs += lterm;
                }
                for (auto& it1: c->get_qterms()) {
                    if (it1.second._coef_p1_tr) { // qterm = (coef*p1)^T*p2
                        for (auto i = 0; i<it1.second._p->first->get_dim(); i++) {
                            for (auto j = 0; j<it1.second._p->first->get_dim(); j++) {
                                coeff = _model->_obj->eval(it1.second._coef,i,j);
                                gvar1 = _grb_vars[it1.second._p->first->get_id() + it1.second._p->first->get_id_inst(j)];
                                gvar2 = _grb_vars[it1.second._p->second->get_id() + it1.second._p->second->get_id_inst(i)];
                                if (!it1.second._sign) {
                                    quadlhs -= coeff*gvar1*gvar2;
                                }
                                else {
                                    quadlhs += coeff*gvar1*gvar2;
                                }
                            }
                        }
                    }
                    else if (it1.second._coef->_is_transposed || it1.second._coef->is_matrix_indexed() || it1.second._p->first->is_matrix_indexed()) {
                        auto dim =it1.second._p->first->get_dim(i);
                        for (int j = 0; j<dim; j++) {
                            coeff = c->eval(it1.second._coef,i,j);
                            gvar1 = _grb_vars[it1.second._p->first->get_id() + it1.second._p->first->get_id_inst(i,j)];
                            gvar2 = _grb_vars[it1.second._p->second->get_id() + it1.second._p->second->get_id_inst(i,j)];
                            if (!it1.second._sign) {
                                quadlhs += -1*coeff*gvar1*gvar2;
                            }
                            else {
                                quadlhs += coeff*gvar1*gvar2;
                            }
                        }
                    }
                    else {
                        gvar1 = _grb_vars[it1.second._p->first->get_id() + it1.second._p->first->get_id_inst(i)];
                        gvar2 = _grb_vars[it1.second._p->second->get_id() + it1.second._p->second->get_id_inst(i)];
                        coeff = c->eval(it1.second._coef,i);
                        if (!it1.second._sign) {
                            quadlhs += -1*coeff*gvar1*gvar2;
                        }
                        else {
                            quadlhs += coeff*gvar1*gvar2;
                        }
                    }
                }
                quadlhs += c->eval(c->get_cst(), i);
                
                if(c->_indices)
                    grb_mod->addQConstr(quadlhs,sense,0,c->get_name()+"("+c->_indices->_keys->at(i)+")");
                else
                    grb_mod->addQConstr(quadlhs,sense,0,c->get_name());
                //                grb_mod->re
                //    }
                
            }
        }
        else{/* This is a General constraint with a unary expression and only two vars, e.g., y = cos(x) or y = min/max(x1,x2..). Refer to https://www.gurobi.com/documentation/9.0/refman/constraints.html#subsubsection:GenConstrFunction
              We currently only support trigonometric functions, log function, and min/max */
            for (size_t i = 0; i< nb_inst; i++){
                linlhs = 0;
                if (c->_lterms->size()!=1) {
                    throw invalid_argument("Gurobi does not support this type of nonlinear constraints");
                }
                auto lt = c->_lterms->begin()->second;
                if (lt._p->_is_vector || lt._p->is_matrix_indexed() || lt._coef->is_matrix()) {
                    throw invalid_argument("Gurobi does not support this type of nonlinear constraints");
                }
                else {
                    coeff = c->eval(lt._coef,i);
                    if (!((coeff==1 && lt._sign) || (coeff==-1 && !lt._sign))) {
                        throw invalid_argument("Gurobi does not support this type of nonlinear constraints");
                    }
                    gvar1 = _grb_vars[lt._p->get_id() + lt._p->get_id_inst(i)];
                }
                if (c->eval(c->get_cst(), i)!=0) {
                    throw invalid_argument("Gurobi does not support this type of nonlinear constraints");
                }
                if(c->_expr->is_uexpr()){
                    if (c->_expr->_coef!=-1) {
                        throw invalid_argument("Gurobi does not support this type of nonlinear constraints");
                    }
                    auto uexp = static_pointer_cast<uexpr<>>(c->_expr);
                    if (uexp->_son->is_function()) {
                        auto f = static_pointer_cast<func<>>(uexp->_son);
                        auto p = f->_vars->begin()->second.first;
                        gvar2 = _grb_vars[p->get_id() + p->get_id_inst(i)];
                    }
                    else if(uexp->_son->is_var()) {
                        auto p = static_pointer_cast<param_>(uexp->_son);
                        gvar2 = _grb_vars[p->get_id() + p->get_id_inst(i)];
                    }
                    else{
                        throw invalid_argument("Error in expression construction");
                    }
                    switch (uexp->_otype) {
                        case gravity::sin_:{
                            if(c->_indices)
                                grb_mod->addGenConstrSin(gvar2, gvar1,c->get_name()+"("+c->_indices->_keys->at(i)+")");
                            else
                                grb_mod->addGenConstrSin(gvar2, gvar1);
                            break;
                        }
                        case gravity::cos_:{
                            if(c->_indices)
                                grb_mod->addGenConstrCos(gvar2, gvar1,c->get_name()+"("+c->_indices->_keys->at(i)+")");
                            else
                                grb_mod->addGenConstrCos(gvar2, gvar1);
                            break;
                        }
                        case gravity::tan_:{
                            if(c->_indices)
                                grb_mod->addGenConstrTan(gvar2, gvar1,c->get_name()+"("+c->_indices->_keys->at(i)+")");
                            else
                                grb_mod->addGenConstrTan(gvar2, gvar1);
                            break;
                        }
                        case gravity::log_:{
                            if(c->_indices)
                                grb_mod->addGenConstrLog(gvar2, gvar1,c->get_name()+"("+c->_indices->_keys->at(i)+")");
                            else
                                grb_mod->addGenConstrLog(gvar2, gvar1);
                            break;
                        }
                        default:
                            break;
                    }
                }
                else{/* Multi expression*/
                    
                    if (c->_expr->_coef!=-1) {
                        throw invalid_argument("Gurobi does not support this type of nonlinear constraints");
                    }
                    auto mexp = static_pointer_cast<mexpr<>>(c->_expr);
                    size_t nb_vars = mexp->_children->size();
                    GRBVar gvars[nb_vars];
                    for (int k = 0; k<nb_vars; k++) {
                        gvars[k] = _grb_vars[mexp->_children->at(k).get_id() + mexp->_children->at(k).get_id_inst(i)];
                    }
                    if(mexp->_otype==min_){
                        if(c->_indices)
                            grb_mod->addGenConstrMin(gvar1, gvars, nb_vars,GRB_INFINITY, c->get_name()+"("+c->_indices->_keys->at(i)+")");
                        else
                            grb_mod->addGenConstrMin(gvar1, gvars, nb_vars);
                    }
                    else {
                        if(c->_indices)
                            grb_mod->addGenConstrMax(gvar1, gvars, nb_vars,GRB_INFINITY, c->get_name()+"("+c->_indices->_keys->at(i)+")");
                        else
                            grb_mod->addGenConstrMax(gvar1, gvars,nb_vars);
                    }
                }
            }
        }
    }
}

void GurobiProgram::set_grb_objective(){
    //    size_t idx = 0;
    GRBLinExpr lterm;
    GRBQuadExpr qobj;
    GRBVar gvar1, gvar2;
    int objt;
    double coeff;
    if (!_model->_obj->_new) {
        return;
    }
    _model->_obj->_new = false;
    if (_model->_objt == minimize) objt = GRB_MINIMIZE;
    else objt = GRB_MAXIMIZE;
    qobj = 0;
    for (auto& it1: _model->_obj->get_lterms()) {
        lterm = 0;
        if (it1.second._coef->_is_transposed || it1.second._coef->is_matrix() || it1.second._p->is_matrix_indexed()) {
            auto dim = it1.second._p->get_dim(0);
            for (size_t j = 0; j<dim; j++) {
                coeff = _model->_obj->eval(it1.second._coef,0,j);
                gvar1 = _grb_vars[it1.second._p->get_id() + it1.second._p->get_id_inst(0,j)];
                lterm += coeff*gvar1;
            }
        }
        else {
            coeff = _model->_obj->eval(it1.second._coef);
            gvar1 = _grb_vars[it1.second._p->get_id() + it1.second._p->get_id_inst()];
            lterm += coeff*gvar1;
        }
        if (!it1.second._sign) {
            lterm *= -1;
        }
        qobj += lterm;
    }
    for (auto& it1: _model->_obj->get_qterms()) {
        if (it1.second._coef_p1_tr) { // qterm = (coef*p1)^T*p2
            for (auto i = 0; i<it1.second._p->first->get_dim(); i++) {
                for (auto j = 0; j<it1.second._p->first->get_dim(); j++) {
                    coeff = _model->_obj->eval(it1.second._coef,i,j);
                    gvar1 = _grb_vars[it1.second._p->first->get_id() + it1.second._p->first->get_id_inst(j)];
                    gvar2 = _grb_vars[it1.second._p->second->get_id() + it1.second._p->second->get_id_inst(i)];
                    if (!it1.second._sign) {
                        qobj -= coeff*gvar1*gvar2;
                    }
                    else {
                        qobj += coeff*gvar1*gvar2;
                    }
                }
            }
        }
        else if (it1.second._coef->_is_transposed) {
            auto dim =it1.second._p->first->get_dim();
            for (int j = 0; j<dim; j++) {
                coeff = _model->_obj->eval(it1.second._coef,j);
                gvar1 = _grb_vars[it1.second._p->first->get_id() + it1.second._p->first->get_id_inst(j)];
                gvar2 = _grb_vars[it1.second._p->second->get_id() + it1.second._p->second->get_id_inst(j)];
                if (!it1.second._sign) {
                    qobj += -1*coeff*gvar1*gvar2;
                }
                else {
                    qobj += coeff*gvar1*gvar2;
                }
            }
        }
        else {
            coeff = _model->_obj->eval(it1.second._coef);
            gvar1 = _grb_vars[it1.second._p->first->get_id() + it1.second._p->first->get_id_inst()];
            gvar2 = _grb_vars[it1.second._p->second->get_id() + it1.second._p->second->get_id_inst()];
            if (!it1.second._sign) {
                qobj += -1*coeff*gvar1*gvar2;
            }
            else {
                qobj += coeff*gvar1*gvar2;
            }
        }
    }
    qobj += _model->_obj->eval(_model->_obj->get_cst());
    grb_mod->setObjective(qobj,objt);
    //    grb_mod->update();
}

void GurobiProgram::print_constraints(){
    GRBConstr* gconstrs = grb_mod->getConstrs();
    for(int i = 0; i < grb_mod->get(GRB_IntAttr_NumConstrs); ++i) {
        if(gconstrs[i].get(GRB_CharAttr_Sense)!='=') {
            cout << gconstrs[i].get(GRB_StringAttr_ConstrName) << "  ";
            //            cout << gconstrs[i].get(GRB_DoubleAttr_Slack);
            cout << "\n";
        }
    }
    GRBQConstr* gqconstrs = grb_mod->getQConstrs();
    for(int i = 0; i < grb_mod->get(GRB_IntAttr_NumQConstrs); ++i) {
        if(gqconstrs[i].get(GRB_CharAttr_QCSense)!='=') {
            cout << gqconstrs[i].get(GRB_StringAttr_QCName) << "  ";
            //            cout << gqconstrs[i].get(GRB_DoubleAttr_Slack);
            //            cout << gqconstrs[i].get(GRB_DoubleAttr_QCSlack);
            cout << "\n";
        }
    }
}
