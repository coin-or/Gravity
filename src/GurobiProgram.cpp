#include <gravity/GurobiProgram.h>
#include <thread>         // std::this_thread::sleep_for
#include <chrono>         // std::chrono::seconds



class cuts: public GRBCallback
{
public:
    vector<GRBVar> vars;
    int nb_vars, nb_spatial, nb_data, nb_model;
    Model<>* m;
    double *x;
    vector<double> cont_x, int_x;
    vector<double> rot_trans;
    vector<int> matching;
    shared_ptr<vector<vector<size_t>>> bin_ids;
    shared_ptr<param<>> x1, x2, y1, y2, z1, z2;
    shared_ptr<param<>> new_x1, new_y1, new_z1;
    shared_ptr<param<>> angle_lb, angle_ub, cos_lb, cos_ub, sin_lb, sin_ub, t_lb, t_ub;
    shared_ptr<param<>> model_voronoi_out_radius;
    shared_ptr<var<>> bin, sbin_roll, sbin_pitch, sbin_yaw, sbin_tx, sbin_ty, sbin_tz;
    shared_ptr<var<>> theta11, theta12, theta13, theta21, theta22, theta23, theta31, theta32, theta33, x_shift, y_shift, z_shift;
    shared_ptr<var<>> x_diff, y_diff, z_diff;
    double objc=999;

    
    int soc_viol, soc_found,soc_added,det_viol, det_found, det_added;
    int soc_viol_user=0, soc_found_user=0,soc_added_user=0,det_viol_user=0, det_found_user=0, det_added_user=0;
    Model<> interior;
    cuts(const vector<GRBVar>& _grb_vars, int xn, Model<>* mod, Model<>& mod_int, int& soc_violn, int& soc_foundn, int& soc_addedn, int& det_violn, int& det_foundn, int& det_addedn) {
        vars = _grb_vars;
        nb_vars = xn;
        x = new double[nb_vars];
        m = mod;
        cont_x.resize(nb_vars);
        int_x.resize(nb_vars);
        soc_viol=soc_violn;
        soc_found=soc_foundn;
        soc_added=soc_addedn;
        det_viol=det_violn;
        det_found=det_foundn;
        det_added=det_addedn;
        interior=mod_int;
        rot_trans.resize(12);
    }
    ~cuts(){
        DebugOff("soc_viol "<<soc_viol_user<<endl);
        DebugOff("soc_found "<<soc_found_user<<endl);
        DebugOff("soc_added "<<soc_added_user<<endl);
        DebugOff("det_viol "<<det_viol_user<<endl);
        DebugOff("det_found "<<det_found_user<<endl);
        DebugOff("det_added "<<det_added_user<<endl);

     
        delete [] x;
    }
protected:
    void callback() {
        try {
            bool incumbent=true;
            bool mipnode=true;
            bool get_lb=true;
            if(incumbent){
                if (where == GRB_CB_MIPSOL) {
                        // Found an integer feasible solution - does it visit every node?
                    int i,j;
                    x=getSolution(vars.data(),nb_vars);
                    for(i=0;i<nb_vars;i++){
                        int_x[i] = x[i];
                    }
                    m->set_solution(int_x);
                    auto res=m->cutting_planes_solution(interior, 1e-9, soc_viol, soc_found, soc_added, det_viol, det_found, det_added);
                    if(res.size()>=1){
                        for(i=0;i<res.size();i++){
                            GRBLinExpr expr = 0;
                            for(j=0;j<res[i].size()-1;j+=2){
                                int c=res[i][j];
                                expr += res[i][j+1]*vars[c];
                            }
                            expr+=res[i][j];
                            addLazy(expr, GRB_LESS_EQUAL, 0);
                        }
                    }
//                    m->set_solution(int_x);
                }
            }
            if(mipnode){
                if (where == GRB_CB_MIPNODE){
                    int stat=getIntInfo(GRB_CB_MIPNODE_STATUS);
                    if(stat==2){
                        DebugOff(getIntInfo(GRB_CB_MIPNODE_STATUS)<<endl);
                        int nct=getDoubleInfo(GRB_CB_MIPNODE_NODCNT);
                        if(nct%100==0){
                            double obj=getDoubleInfo(GRB_CB_MIPNODE_OBJBST);
                            double obj1=getDoubleInfo(GRB_CB_MIPNODE_OBJBND);
                            DebugOff(obj<<"\t"<<obj1<<"\t"<<endl);
                            int i,j;
                            m->get_solution(int_x);
                            x=getNodeRel(vars.data(),nb_vars);
                            for(i=0;i<nb_vars;i++){
                                cont_x[i] = x[i];
                            }
                            m->set_solution(cont_x);
                            auto res=m->cutting_planes_solution(interior, 1e-9,soc_viol_user, soc_found_user,soc_added_user,det_viol_user, det_found_user, det_added_user);
                            if(res.size()>=1){
                                for(i=0;i<res.size();i++){
                                    GRBLinExpr expr = 0;
                                    for(j=0;j<res[i].size()-1;j+=2){
                                        int c=res[i][j];
                                        expr += res[i][j+1]*vars[c];
                                    }
                                    expr+=res[i][j];
                                    addCut(expr, GRB_LESS_EQUAL, 0);
                                }
                            }
                        }
                        m->set_solution(int_x);
                    }
                    
                }
            }
           // if(get_lb)
//            {
//                if (where == GRB_CB_MIPNODE){
//                    int stat=getIntInfo(GRB_CB_MIPNODE_STATUS);
//                    if(stat==2){
//                        m->get_solution(int_x);
//                        double* xa = getNodeRel(vars.data(), nb_vars);
//                        for(auto i=0;i<nb_vars;i++){
//                            cont_x[i]=xa[i];
//                        }
//                        setSolution(vars.data(), xa, nb_vars);
//                        delete[] xa;
//                        m->set_solution(cont_x);
//                       // m->reset();
//                        m->_obj->eval();
//                        objc=m->get_obj_val();
//                        m->set_solution(int_x);
//                    }
//                }
//            }
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
    grb_env->set(GRB_DoubleParam_TimeLimit,9000);
        //    grb_env->set(GRB_DoubleParam_MIPGap,0.01);
        //   grb_env->set(GRB_IntParam_Threads,8);
        //    grb_env->set(GRB_IntParam_Presolve,0);
        //   grb_env->set(GRB_IntParam_NumericFocus,3);
    //grb_env->set(GRB_IntParam_NonConvex,2);
        // grb_env->set(GRB_DoubleParam_FeasibilityTol, 1e-8);
        //    grb_env->set(GRB_DoubleParam_OptimalityTol, 1e-8);
    
    
    //grb_env->set(GRB_IntParam_OutputFlag,0);
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
            grb_env->set(GRB_DoubleParam_TimeLimit,9000);
                //   grb_env->set(GRB_IntParam_Threads,8);
                //    grb_env->set(GRB_DoubleParam_MIPGap,0.01);
                //    grb_env->set(GRB_IntParam_Threads,1);
                //    grb_env->set(GRB_IntParam_Presolve,0);
                //     grb_env->set(GRB_IntParam_NumericFocus,3);
            //grb_env->set(GRB_IntParam_NonConvex,2);
                //        grb_env->set(GRB_DoubleParam_FeasibilityTol, 1e-8);
                //            grb_env->set(GRB_DoubleParam_OptimalityTol, 1e-8);
            
            //grb_env->set(GRB_IntParam_OutputFlag,0);
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
            grb_env->set(GRB_DoubleParam_TimeLimit,9000);
                //    grb_env->set(GRB_DoubleParam_MIPGap,0.01);
                //    grb_env->set(GRB_IntParam_Threads,1);
                //    grb_env->set(GRB_IntParam_Presolve,0);
                //   grb_env->set(GRB_IntParam_NumericFocus,3);
            //grb_env->set(GRB_IntParam_NonConvex,2);
                //          grb_env->set(GRB_DoubleParam_FeasibilityTol, 1e-8);
                //            grb_env->set(GRB_DoubleParam_OptimalityTol, 1e-8);
            
                // grb_env->set(GRB_IntParam_OutputFlag,1);
           // grb_env->set(GRB_IntParam_OutputFlag,0);
            grb_mod = new GRBModel(*grb_env);
            //grb_mod->set(GRB_IntParam_LazyConstraints, 1);
            
            
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

bool GurobiProgram::solve(bool relax, double mipgap, bool use_callback, double max_time, double cut_off, int threads){
        //cout << "\n Presolve = " << grb_env->get(GRB_IntParam_Presolve) << endl;
        //    print_constraints();
    if (relax) relax_model();
        //    relax_model();
//    grb_mod->set(GRB_DoubleParam_MIPGap, 1e-6);
    grb_mod->set(GRB_DoubleParam_FeasibilityTol, 1e-6);
    grb_mod->set(GRB_DoubleParam_OptimalityTol, 1e-6);
//    grb_mod->set(GRB_DoubleParam_BarConvTol, 1e-6);
//    grb_mod->set(GRB_DoubleParam_BarQCPConvTol, 1e-6);
    grb_mod->set(GRB_IntParam_StartNodeLimit,-3);
    grb_mod->set(GRB_IntParam_Threads, threads);
    grb_mod->set(GRB_IntParam_OutputFlag,1);
        //    if(use_callback){
//    grb_mod->set(GRB_DoubleParam_NodefileStart,0.1);
    grb_mod->set(GRB_IntParam_NonConvex,2);
    //grb_mod->set(GRB_IntParam_MIPFocus,2);
    grb_mod->set(GRB_IntParam_BranchDir, 1);
    grb_mod->set(GRB_IntParam_CutPasses, 5);
   // grb_mod->set(GRB_IntParam_PrePasses, 2);
    grb_mod->set(GRB_DoubleParam_TimeLimit,max_time);
    
    grb_mod->set(GRB_DoubleParam_Cutoff,cut_off);
   // grb_mod->set(GRB_DoubleParam_Heuristics,0);
    //grb_mod->set(GRB_IntParam_Cuts,3);
    //grb_mod->set(GRB_DoubleParam_BestBdStop,cut_off);
//    if(use_callback){
        grb_mod->getEnv().set(GRB_IntParam_DualReductions, 0);
        grb_mod->getEnv().set(GRB_IntParam_PreCrush, 1);
        grb_mod->getEnv().set(GRB_IntParam_LazyConstraints, 1);
//    }
    grb_mod->update();
    int n=grb_mod->get(GRB_IntAttr_NumVars);
    if(n==0)
        cout << "Gurobi model has zero variables!\n";
    if(n!=_model->get_nb_vars())
        throw invalid_argument("Number of variables in Gurobi model does not match Gravity!");
    Model<> interior;
    auto lin=_model->buildOA();
    int soc_viol=0,soc_found=0,soc_added=0,det_viol=0,det_found=0,det_added=0;
    vector<int> stats;
    stats.resize(6,0);
    //if(use_callback){
        interior=lin->add_outer_app_solution(*_model);
//    //}
    cuts cb(_grb_vars, n, _model, interior, soc_viol,soc_found,soc_added,det_viol,det_found,det_added);
    grb_mod->setCallback(&cb);
    
    
    grb_mod->optimize();
        //            grb_mod->write("~/mod.mps");
    if (grb_mod->get(GRB_IntAttr_Status) == 3 || grb_mod->get(GRB_IntAttr_Status) == 6) {
        cerr << "Infeasible model!\n";
//        _model->print_var_bounds();
        return false;
//        cerr << "\nModel has not been solved to optimality, error code = " << grb_mod->get(GRB_IntAttr_Status) << endl;
            //        return false;
    }
 
    
        //    GRBVar* gvars = grb_mod->getVars();
        //    for(int i = 0; i < grb_mod->get(GRB_IntAttr_NumVars); ++i) {
        ////        cout << gvars[i].get(GRB_StringAttr_VarName) << "  " << gvars[i].get(GRB_DoubleAttr_X) << endl;
        //        if (gvars[i].get(GRB_CharAttr_VType)==GRB_BINARY) {
        //            cout << gvars[i].get(GRB_StringAttr_VarName) << "  ";
        //            cout << gvars[i].get(GRB_DoubleAttr_X);
        //            cout << "\n";
        //        }
        //    }
    double res=-999.0;
    if(grb_mod->get(GRB_IntAttr_SolCount) >= 1){
        _model->_obj->set_val(grb_mod->get(GRB_DoubleAttr_ObjVal));
        update_solution();
    }
    else{
        _model->_obj->set_val(res);
    }
    _model->_rel_obj_val = grb_mod->get(GRB_DoubleAttr_ObjBound);
//    cout << "\n***** Optimal Objective = " << _model->get_obj_val() << " *****\n";
//    if (grb_mod->get(GRB_IntAttr_IsMIP)) {
//        cout.setf(ios::fixed);
//        cout.precision(3);
//        cout << "Results: " << grb_mod->get(GRB_DoubleAttr_ObjVal) << " & ";
//        cout.precision(4);
//        cout << (grb_mod->get(GRB_DoubleAttr_MIPGap))*100 << "% & ";
//        cout.precision(0);
//        cout << grb_mod->get(GRB_DoubleAttr_NodeCount) << " & ";
//        cout.precision(2);
//        cout << grb_mod->get(GRB_DoubleAttr_Runtime) << " & " << endl;
//    }
        //    delete[] gvars;
    return true;
}

void GurobiProgram::prepare_model(){
    _model->fill_in_maps();
    _model->compute_funcs();
    fill_in_grb_vmap();
    create_grb_constraints();
    set_grb_objective();
//            grb_mod->write("gurobiprint.lp");
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
        //    for (auto i = 0; i < _grb_vars.size(); i++) {
        //        gvar = _grb_vars.at(i);
        //        auto dim = _model->_vars[i]->get_dim();
        //        for (auto j = 0; j < _model->_vars[i]->get_dim(); j++) {
        //            poly_set_val(j, gvar.get(GRB_DoubleAttr_X), _model->_vars[i]);
        //        }
        //    }
    for(auto& v_p: _model->_vars)
    {
        v = v_p.second.get();
        auto idx = v->get_id();
        auto dim = v->_dim[0];
        for (auto i = 0; i < dim; i++) {
            auto vid = idx + v->get_id_inst(i);
            gvar = _grb_vars.at(vid);
            v->set_double_val(i,gvar.get(GRB_DoubleAttr_X));
        }
    }
}

void GurobiProgram::relax_model(){
    GRBVar* gvars = grb_mod->getVars();
    for(int i = 0; i < grb_mod->get(GRB_IntAttr_NumVars); ++i) {
        if (gvars[i].get(GRB_CharAttr_VType) == 'B') gvars[i].set(GRB_CharAttr_VType,'C');
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
//        v->_new = false;
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
                        //_grb_vars.at(vid).set(GRB_DoubleAttr_Start, real_var->eval(i));
                        //_grb_vars.at(vid).set(GRB_IntAttr_BranchPriority, real_var->_priority);/* branching priority for integers */
                    }
                    else {
                        _grb_vars.at(vid) = (GRBVar(grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_CONTINUOUS, v->get_name(true,true)+"("+v->_indices->_keys->at(i)+")")));
                        //_grb_vars.at(vid).set(GRB_DoubleAttr_Start, real_var->eval(i));
                    }
                }
                break;
            }
            case integer_:{
                auto real_var = (var<int>*)v;
                for (int i = 0; i < real_var->_dim[0]; i++) {
                    auto vid = idx + i;
                    //if(i<=500){
                    _grb_vars.at(vid) = (GRBVar(grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_BINARY, v->get_name(true,true)+"("+v->_indices->_keys->at(i)+")")));
//                    }else{
//                        _grb_vars.at(vid) = (GRBVar(grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_CONTINUOUS, v->get_name(true,true)+"("+v->_indices->_keys->at(i)+")")));
//                    }
                    //_grb_vars.at(vid).set(GRB_DoubleAttr_Start, real_var->eval(i));
                    //_grb_vars.at(vid).set(GRB_IntAttr_BranchPriority, real_var->_priority);/* branching priority for integers */
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
        if (!c->_new && c->_all_satisfied){
            continue;
        }
        if(c->_callback)
        {
            DebugOff(c->_name<<"  lazy"<<endl);
            continue;
        }
//        c->_new = false;
        
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
//    _model->_obj->_new = false;
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

