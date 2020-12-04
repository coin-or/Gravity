#include <gravity/GurobiProgram.h>
#include <thread>         // std::this_thread::sleep_for
#include <chrono>         // std::chrono::seconds

//GurobiProgram::GurobiProgram(){
////    model = m;
//    grb_env = new GRBEnv();
////    grb_env->set(GRB_IntParam_Presolve,0);
//    //grb_env->set(GRB_DoubleParam_NodeLimit,1);
//    grb_env->set(GRB_DoubleParam_TimeLimit,7200);
////    grb_env->set(GRB_DoubleParam_MIPGap,0.01);
//    grb_env->set(GRB_IntParam_Threads,1);
////       grb_env->set(GRB_IntParam_Presolve,0);
////      grb_env->set(GRB_IntParam_NumericFocus,3);
////     grb_env->set(GRB_IntParam_NonConvex,2);
//    //grb_env->set(GRB_DoubleParam_FeasibilityTol, 1E-6);
//     grb_env->set(GRB_DoubleParam_OptimalityTol, 1E-6);
//
//    grb_env->set(GRB_IntParam_OutputFlag,1);
////    grb_mod = new GRBModel(*grb_env);
//    grb_mod = NULL;
//}


GurobiProgram::GurobiProgram(Model<>* m) {
    bool found_token = false;
    while (!found_token) {
        try{
            grb_env = new GRBEnv();
        //    grb_env->set(GRB_IntParam_Presolve,0);
            //grb_env->set(GRB_DoubleParam_NodeLimit,1);
            //grb_env->set(GRB_DoubleParam_TimeLimit,7200);
            //    grb_env->set(GRB_DoubleParam_MIPGap,0.01);
            //grb_env->set(GRB_IntParam_Threads,1);
        //    grb_env->set(GRB_IntParam_Presolve,0);
        //    grb_env->set(GRB_IntParam_NumericFocus,3);
        //    grb_env->set(GRB_IntParam_NonConvex,2);
         //   grb_env->set(GRB_DoubleParam_FeasibilityTol, 1E-6);
         //   grb_env->set(GRB_DoubleParam_OptimalityTol, 1E-6);
            
            //grb_env->set(GRB_IntParam_OutputFlag,0);
            grb_mod = new GRBModel(*grb_env);
	    //grb_mod->set(GRB_IntParam_Method, 1);
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
    _model = m;
    m->fill_in_maps();
    m->compute_funcs();
    //DebugOn("found token "<<endl);
}

//GurobiProgram::GurobiProgram(const shared_ptr<Model<>>& m) {
//    grb_env = new GRBEnv();
//    //grb_env->set(GRB_IntParam_Presolve,0);
//    //grb_env->set(GRB_DoubleParam_NodeLimit,1);
//    grb_env->set(GRB_DoubleParam_TimeLimit,7200);
//    //    grb_env->set(GRB_DoubleParam_MIPGap,0.01);
//        grb_env->set(GRB_IntParam_Threads,1);
////    grb_env->set(GRB_IntParam_Presolve,0);
////    grb_env->set(GRB_IntParam_NumericFocus,3);
////    grb_env->set(GRB_IntParam_NonConvex,2);
//    //grb_env->set(GRB_DoubleParam_FeasibilityTol, 1E-6);
//    grb_env->set(GRB_DoubleParam_OptimalityTol, 1E-6);
//
//    grb_env->set(GRB_IntParam_OutputFlag,1);
//    grb_mod = new GRBModel(*grb_env);
////    grb_env->set(GRB_IntParam_ƒFlag,2);
////    _model = m;
//    m->fill_in_maps();
//    m->compute_funcs();
//}


GurobiProgram::~GurobiProgram() {
//    for (auto p : _grb_vars) delete p.second;
    if (grb_mod) {
        delete grb_mod;
        grb_mod = nullptr;
    }
    if (grb_env){
        delete grb_env;
        grb_env = nullptr;
    }
   
}

void GurobiProgram::reset_model(){
    if (grb_mod) delete grb_mod;
    _grb_vars.clear();
//    grb_env->set(GRB_IntParam_OutputFlag,_output);
    grb_mod = new GRBModel(*grb_env);
}

bool GurobiProgram::solve(int output, bool relax, double tol, double mipgap, bool gurobi_crossover){
    //cout << "\n Presolve = " << grb_env->get(GRB_IntParam_Presolve) << endl;
//    print_constraints();
 //   if (relax) relax_model();
//    relax_model();
    //grb_mod->set(GRB_DoubleParam_BarConvTol, 1e-8);
      //grb_mod->set(GRB_IntParam_ScaleFlag, 1);
      grb_mod->set(GRB_DoubleParam_FeasibilityTol, tol);
      grb_mod->set(GRB_DoubleParam_OptimalityTol, tol);
    //grb_mod->set(GRB_DoubleParam_MIPGap, mipgap);
      grb_mod->set(GRB_IntParam_Threads, 1);
      grb_mod->set(GRB_DoubleParam_TimeLimit, 2000.0);
    ///grb_mod->set(GRB_IntParam_NumericFocus, 1);
      grb_mod->set(GRB_IntParam_Method, 1);
//    if(!gurobi_crossover){
//        grb_mod->set(GRB_IntParam_Crossover, 0);
//    }
    grb_mod->set(GRB_IntParam_OutputFlag, 0);
//    warm_start(); // No need to reset variables if Gurobi model has not changed.
    //grb_mod->write("gurobiprint.lp");
    try{
        grb_mod->optimize();
    }
    catch(GRBException e) {
        cerr << "\nWas not able to optimize Gurobi model, Error code = " << e.getErrorCode() << endl;
        cerr << e.getMessage() << endl;
        exit(-1);
    }
    //cout<<"Status "<<grb_mod->get(GRB_IntAttr_Status)<<endl;

    // cout<<"BoundVio "<<grb_mod->get(GRB_DoubleAttr_BoundVio)<<endl;
    // cout<<"BoundSVio "<<grb_mod->get(GRB_DoubleAttr_BoundSVio)<<endl;
    // cout<<"ConstrVio "<<grb_mod->get(GRB_DoubleAttr_ConstrVio)<<endl;
    // cout<<"ConstrSVio "<<grb_mod->get(GRB_DoubleAttr_ConstrSVio)<<endl;
    // cout<<"ConstrResidual "<<grb_mod->get(GRB_DoubleAttr_ConstrResidual)<<endl;
    // cout<<"ConstrSResidual "<<grb_mod->get(GRB_DoubleAttr_ConstrSResidual)<<endl;
 
    //grb_mod->write("~/mod.mps");
    grb_first_run=false;
    if (grb_mod->get(GRB_IntAttr_Status) != 2) {
        cerr << "\nModel has not been solved to optimality, error code = " << grb_mod->get(GRB_IntAttr_Status) << endl;
        return false;
    }
    update_solution();
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
//    cout << "\n***** Optimal Objective = " << _model->get_obj_val() << " *****\n";
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
    return true;
}

void GurobiProgram::prepare_model(){
    _model->fill_in_maps();
    _model->compute_funcs();
//    DebugOn("going to fill in vmap"<<endl);
    fill_in_grb_vmap();
//    DebugOn("going to fill in cmap"<<endl);
    create_grb_constraints();
//    DebugOn("going to fill in omap"<<endl);
    set_grb_objective();

//    grb_mod->write("gurobiprint.lp");

//    print_constraints();
}
void GurobiProgram::initialize_basis(std::vector<int> vbasis, std::vector<int> cbasis){
	grb_mod->update();    
	int nv=grb_mod->get(GRB_IntAttr_NumVars);
    int nc=grb_mod->get(GRB_IntAttr_NumConstrs);
    int count=0;
    for(auto &gv:_grb_vars){
        gv.set(GRB_IntAttr_VBasis, vbasis[count++]);
    }

    GRBConstr* gcons= grb_mod->getConstrs();
  
//DebugOn("nc "<<nc<<endl);
//DebugOn("cbasis "<<cbasis.size()<<endl);
    for(auto i=0;i<cbasis.size();i++){
        gcons[i].set(GRB_IntAttr_CBasis, cbasis[i]);
    }
  for(auto i=cbasis.size();i<nc;i++){                        
        gcons[i].set(GRB_IntAttr_CBasis, 0);          
    }  
delete[] gcons;
}

void GurobiProgram::get_basis(std::vector<int>& vbasis, std::vector<int>& cbasis){
    int nv=grb_mod->get(GRB_IntAttr_NumVars);
    int nc=grb_mod->get(GRB_IntAttr_NumConstrs);
    vbasis.resize(nv);
    cbasis.resize(nc);
    GRBVar* gvars= grb_mod->getVars();
    for(auto i=0;i<nv;i++){
        vbasis[i]=gvars[i].get(GRB_IntAttr_VBasis);
    }

    GRBConstr* gcons= grb_mod->getConstrs();

    for(auto i=0;i<nc;i++){
        cbasis[i]=gcons[i].get(GRB_IntAttr_CBasis);
    }
delete[] gcons;
delete[] gvars;
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

void GurobiProgram::warm_start(){
    GRBVar gvar;
    param_* v;
    double value;
    for(auto& v_p: _model->_vars)
    {
        v = v_p.second.get();
        auto idx = v->get_id();
        auto dim = v->_dim[0];
        for (auto i = 0; i < dim; i++) {
            auto vid = idx + i;
            v->get_double_val(i,value);
            _grb_vars.at(vid).set(GRB_DoubleAttr_PStart, value);
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
    auto size_init=_grb_vars.size();
    _grb_vars.resize(_model->get_nb_vars());
    bool add;
    for(auto& v_p: _model->_vars)
    {
        v = v_p.second.get();
        if (!v->_new) {
            continue;
        }
        if(v->_new && size_init>v->get_id())
        {
            add=false;
        }
        else{
            add=true;
        }
        //DebugOn("to add v"<<endl);
        v->_new = false;
        auto idx = v->get_id();
        switch (v->get_intype()) {
            case float_: {
                auto real_var = (var<float>*)v;
                for (int i = 0; i < real_var->_dim[0]; i++) {
                    auto vid = idx + i;
                    if(add){
                    _grb_vars.at(vid) = (GRBVar(grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_CONTINUOUS, v->get_name(true,true)+"("+v->_indices->_keys->at(i)+")")));
                }
                else{
                    _grb_vars.at(vid).set(GRB_DoubleAttr_LB, real_var->get_lb(i));
                    _grb_vars.at(vid).set(GRB_DoubleAttr_UB, real_var->get_ub(i));
                }
                }
                break;
            }
            case long_:{
                auto real_var = (var<long double>*)v;
                for (int i = 0; i < real_var->_dim[0]; i++) {
                    auto vid = idx + i;
                    if(add){
                    _grb_vars.at(vid) = (GRBVar(grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_CONTINUOUS, v->get_name(true,true)+"("+v->_indices->_keys->at(i)+")")));
                }
                else{
                    _grb_vars.at(vid).set(GRB_DoubleAttr_LB, real_var->get_lb(i));
                    _grb_vars.at(vid).set(GRB_DoubleAttr_UB, real_var->get_ub(i));
                }
                }
                break;
            }
            case double_:{
                auto real_var = (var<double>*)v;
                for (size_t i = 0; i < real_var->_dim[0]; i++) {
                    auto vid = idx + i;
                    if(add){
                    _grb_vars.at(vid) = (GRBVar(grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_CONTINUOUS, v->get_name(true,true)+"("+v->_indices->_keys->at(i)+")")));
                    //DebugOn("added var"<<endl);
                }
                    else{
                    _grb_vars.at(vid).set(GRB_DoubleAttr_LB, real_var->get_lb(i));
                    _grb_vars.at(vid).set(GRB_DoubleAttr_UB, real_var->get_ub(i));
                        DebugOff("updated bounds"<<endl);
                }
                }
                break;
            }
            case integer_:{
                auto real_var = (var<int>*)v;
                for (int i = 0; i < real_var->_dim[0]; i++) {
                    auto vid = idx + i;
                    if(add){
                    _grb_vars.at(vid) = (GRBVar(grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_INTEGER, v->get_name(true,true)+"("+v->_indices->_keys->at(i)+")")));
                }
                else{
                    _grb_vars.at(vid).set(GRB_DoubleAttr_LB, real_var->get_lb(i));
                    _grb_vars.at(vid).set(GRB_DoubleAttr_UB, real_var->get_ub(i));
                }

            }
                break;
            }
            case short_:{
                auto real_var = (var<short>*)v;
                for (int i = 0; i < real_var->_dim[0]; i++) {
                    auto vid = idx + i;
                    if(add){
                    _grb_vars.at(vid) = (GRBVar(grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_INTEGER, v->get_name(true,true)+"("+v->_indices->_keys->at(i)+")")));
                }
                else{
                    _grb_vars.at(vid).set(GRB_DoubleAttr_LB, real_var->get_lb(i));
                    _grb_vars.at(vid).set(GRB_DoubleAttr_UB, real_var->get_ub(i));
                }
                }
                break;
            }
            case binary_:{
                auto real_var = (var<bool>*)v;
                for (int i = 0; i < real_var->_dim[0]; i++) {
                    auto vid = idx + i;
                    if(add){
                    _grb_vars.at(vid) = (GRBVar(grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_BINARY, v->get_name(true,true)+"("+v->_indices->_keys->at(i)+")")));
                   }
                   else{
                    _grb_vars.at(vid).set(GRB_DoubleAttr_LB, real_var->get_lb(i));
                    _grb_vars.at(vid).set(GRB_DoubleAttr_UB, real_var->get_ub(i));
                }
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
        c->_new = false;
        //DebugOn("to add constraint"<<endl);
        if (c->is_nonlinear()) {
            DebugOn("nonlinear"<<endl);
            throw invalid_argument("Gurobi cannot handle nonlinear constraints that are not convex quadratic.\n");
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
                if (c->_violated[i]) {
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
                    //DebugOn("added constraint"<<endl);
                }
            }
        }
        else {
            for (size_t i = 0; i< nb_inst; i++){
                if (c->_violated[i]) {
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
//            assert(it_qterm.second._p->first->_dim[1]==1 && it_qterm.second._coef->_dim[0]==it_qterm.second._p->second->_dim[0]);
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
