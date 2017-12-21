#include <Gravity/GurobiProgram.h>

GurobiProgram::GurobiProgram(){
//    model = m;
    grb_env = new GRBEnv();
//    grb_env->set(GRB_IntParam_Presolve,0);
    //grb_env->set(GRB_DoubleParam_NodeLimit,1);
    grb_env->set(GRB_DoubleParam_TimeLimit,7200);
//    grb_env->set(GRB_DoubleParam_MIPGap,0.01);
//    grb_env->set(GRB_IntParam_Threads,1);
//    grb_env->set(GRB_IntParam_OutputFlag,0);
//    grb_mod = new GRBModel(*grb_env);
    grb_mod = NULL;
}

GurobiProgram::GurobiProgram(Model* m):GurobiProgram(){
    grb_mod = new GRBModel(*grb_env);
    _model = m;
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

bool GurobiProgram::solve(bool relax){
    //cout << "\n Presolve = " << grb_env->get(GRB_IntParam_Presolve) << endl;
//    print_constraints();
    if (relax) relax_model();
//    relax_model();
    grb_mod->optimize();
//    grb_mod->write("~/mod.mps");
    if (grb_mod->get(GRB_IntAttr_Status) != 2) {
        cerr << "\nModel has not been solved to optimality, error code = " << grb_mod->get(GRB_IntAttr_Status) << endl;
        return false;
    }
    update_model();
    GRBVar* gvars = grb_mod->getVars();
    for(int i = 0; i < grb_mod->get(GRB_IntAttr_NumVars); ++i) {
//        cout << gvars[i].get(GRB_StringAttr_VarName) << "  " << gvars[i].get(GRB_DoubleAttr_X) << endl;
        if (gvars[i].get(GRB_CharAttr_VType)==GRB_BINARY) {
            cout << gvars[i].get(GRB_StringAttr_VarName) << "  ";
            cout << gvars[i].get(GRB_DoubleAttr_X);
            cout << "\n";
        }
    }
    cout << "\n***** Optimal Objective = " << grb_mod->get(GRB_DoubleAttr_ObjVal) << " *****\n";
    _model->_obj_val = grb_mod->get(GRB_DoubleAttr_ObjVal);
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
    delete[] gvars;
    return true;
}

void GurobiProgram::prepare_model(){
    fill_in_grb_vmap();
    create_grb_constraints();
    set_grb_objective();
//    print_constraints();
}

void GurobiProgram::update_model(){
    size_t vid, vid_inst;
    GRBVar gvar;
    for (auto& it: _model->_vars){
        switch (it.second->get_type()) {
            case binary_c:{
                auto vb =  (var<bool>*)it.second;
                vid = it.second->get_id();
                auto dim = vb->get_dim();
                for (int i = 0; i < dim; i++) {
                    vid_inst = vid + it.second->get_id_inst(i);
                    gvar = _grb_vars.at(vid_inst);
                    vb->set_val(vid_inst, static_cast<int>(gvar.get(GRB_DoubleAttr_X) + 0.5) == 1);
                }
                break;
            }
            case integer_c:{
                auto vi =  (var<int>*)it.second;
                vid = it.second->get_id();
                auto dim = vi->get_dim();
                for (int i = 0; i < dim; i++) {
                    vid_inst = vid + it.second->get_id_inst(i);
                    gvar = _grb_vars.at(vid_inst);
                    vi->set_val(vid_inst, static_cast<int>(gvar.get(GRB_DoubleAttr_X) + 0.5) == 1);
                }
                break;
            }
            case short_c:{
                auto vs =  (var<short>*)it.second;
                vid = it.second->get_id();
                auto dim = vs->get_dim();
                for (int i = 0; i < dim; i++) {
                    vid_inst = vid + it.second->get_id_inst(i);
                    gvar = _grb_vars.at(vid_inst);
                    vs->set_val(vid_inst, static_cast<short>(gvar.get(GRB_DoubleAttr_X) + 0.5) == 1);
                }
                break;
            }
            case float_c:{
                auto vf =  (var<float>*)it.second;
                vid = it.second->get_id();
                auto dim = vf->get_dim();
                for (int i = 0; i < dim; i++) {
                    vid_inst = vid + it.second->get_id_inst(i);
                    gvar = _grb_vars.at(vid_inst);
                    vf->set_val(vid_inst, gvar.get(GRB_DoubleAttr_X));
                }
                break;
            }
            case double_c:{
                auto vd =  (var<double>*)it.second;
                vid = it.second->get_id();
                auto dim = vd->get_dim();
                for (int i = 0; i < dim; i++) {
                    vid_inst = vid + it.second->get_id_inst(i);
                    gvar = _grb_vars.at(vid_inst);
                    vd->set_val(vid_inst,gvar.get(GRB_DoubleAttr_X));
                    DebugOn("\ngrbvar name = " << gvar.get(GRB_StringAttr_VarName) << " grbvar value = " << gvar.get(GRB_DoubleAttr_X));
                }
                break;
            }
            case long_c:{
                auto vl =  (var<long double>*)it.second;
                vid = it.second->get_id();
                auto dim = vl->get_dim();
                for (int i = 0; i < dim; i++) {
                    gvar = _grb_vars.at(vid+i);
                    vl->set_val(vid+i,gvar.get(GRB_DoubleAttr_X));
                    DebugOff("\ngrbvar name = " << gvarit->second->get(GRB_StringAttr_VarName) << " grbvar value = " << gvarit->second->get(GRB_DoubleAttr_X));
                }
                break;
            }
            default:
                break;
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
        v = v_p.second;
        auto idx = v->get_id();
        switch (v->get_intype()) {
            case float_: {
                auto real_var = (var<float>*)v;
                for (int i = 0; i < real_var->get_dim(); i++) {
                    auto vid = idx + v->get_id_inst(i);
                    _grb_vars.at(vid) = (GRBVar(grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_CONTINUOUS, v->get_name()+"_"+to_string(i))));
                }
                break;
            }
            case long_:{
                auto real_var = (var<long double>*)v;
                for (int i = 0; i < real_var->get_dim(); i++) {
                    auto vid = idx + v->get_id_inst(i);
                    _grb_vars.at(vid) = (GRBVar(grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_CONTINUOUS, v->get_name()+"_"+to_string(i))));
                }
                break;
            }
            case double_:{
                auto real_var = (var<double>*)v;
                for (int i = 0; i < real_var->get_dim(); i++) {
                    auto vid = idx + v->get_id_inst(i);
                    _grb_vars.at(vid) = (grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_CONTINUOUS, v->get_name()+"_"+to_string(i)));
                }
                break;
            }
            case integer_:{
                auto real_var = (var<int>*)v;
                for (int i = 0; i < real_var->get_dim(); i++) {
                    auto vid = idx + v->get_id_inst(i);
                    _grb_vars.at(vid) = (GRBVar(grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_INTEGER, v->get_name()+"_"+to_string(i))));
                }
                break;
            }
            case short_:{
                auto real_var = (var<short>*)v;
                for (int i = 0; i < real_var->get_dim(); i++) {
                    auto vid = idx + v->get_id_inst(i);
                    _grb_vars.at(vid) = (GRBVar(grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_INTEGER, v->get_name()+"_"+to_string(i))));
                }
                break;
            }
            case binary_:{
                auto real_var = (var<bool>*)v;
                for (int i = 0; i < real_var->get_dim(); i++) {
                    auto vid = idx + v->get_id_inst(i);
                    _grb_vars.at(vid) = (GRBVar(grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_BINARY, v->get_name()+"_"+to_string(i))));
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
//        for (int i = 0; i < real_var->get_dim(); i++) {
//            auto vid = v->get_id() + v->get_id_inst(i);
//            DebugOn("VID = "<< vid <<" : " << _grb_vars.at(vid).get(GRB_StringAttr_VarName) << " in [" << _grb_vars.at(vid).get(GRB_DoubleAttr_LB) << "," << _grb_vars.at(vid).get(GRB_DoubleAttr_UB) << "]\n" );
//        }
//    }
}

void GurobiProgram::create_grb_constraints(){
    char sense;
    size_t nb_inst = 0;
    GRBLinExpr lterm, linlhs;
    GRBQuadExpr quadlhs;
    GRBVar gvar1, gvar2;
    double coeff;    
    for(auto& p: _model->_cons){
        auto c = p.second;
        if (c->is_nonlinear()) {
            throw invalid_argument("Gurobi cannot handle nonlinear constraints that are not convex quadratic.\n");
        }
        switch(c->get_type()) {
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
        nb_inst = c->get_nb_instances();
        if (c->is_linear()) {
            for (int i = 0; i< nb_inst; i++){
                linlhs = 0;
                for (auto& it1: c->get_lterms()) {
                    lterm = 0;
                    if (it1.second._coef->_is_transposed) {
                        auto dim =it1.second._p->get_dim();
                        for (int j = 0; j<dim; j++) {
                            coeff = poly_eval(it1.second._coef,j);
                            gvar1 = _grb_vars[it1.second._p->get_id() + it1.second._p->get_id_inst(j)];
                            lterm += coeff*gvar1;
                        }
                    }
                    else {
                        coeff = poly_eval(it1.second._coef,i);
                        gvar1 = _grb_vars[it1.second._p->get_id() + it1.second._p->get_id_inst(i)];
                        lterm += coeff*gvar1;
                    }
                    if (!it1.second._sign) {
                        lterm *= -1;
                    }
                    linlhs += lterm;
                }
                linlhs += poly_eval(c->get_cst(), i);                
                grb_mod->addConstr(linlhs,sense,c->get_rhs(),c->get_name()+"_"+to_string(i));
            }
        }
        else {
            for (int i = 0; i< nb_inst; i++){
                quadlhs = 0;
                for (auto& it1: c->get_lterms()) {
                    lterm = 0;
                    if (it1.second._coef->_is_transposed) {
                        auto dim =it1.second._p->get_dim();
                        for (int j = 0; j<dim; j++) {
                            coeff = poly_eval(it1.second._coef,j);
                            gvar1 = _grb_vars[it1.second._p->get_id() + it1.second._p->get_id_inst(j)];
                            lterm += coeff*gvar1;
                        }
                    }
                    else {
                        coeff = poly_eval(it1.second._coef,i);
                        gvar1 = _grb_vars[it1.second._p->get_id() + it1.second._p->get_id_inst(i)];
                        lterm += coeff*gvar1;
                    }
                    if (!it1.second._sign) {
                        lterm *= -1;
                    }
                    quadlhs += lterm;
                }
                for (auto& it1: c->get_qterms()) {
                    gvar1 = _grb_vars[it1.second._p->first->get_id() + it1.second._p->first->get_id_inst(i)];
                    gvar2 = _grb_vars[it1.second._p->second->get_id() + it1.second._p->second->get_id_inst(i)];
                    if (it1.second._coef->_is_transposed) {
                        throw invalid_argument("unsupported operation");
                    }
                    else {
                        coeff = poly_eval(it1.second._coef,i);
                    }
                    if (!it1.second._sign) {
                        quadlhs += -1*coeff*gvar1*gvar2;
                    }
                    else {
                        quadlhs += coeff*gvar1*gvar2;
                    }
                }
                quadlhs += poly_eval(c->get_cst(), i);
                grb_mod->addQConstr(quadlhs,sense,c->get_rhs(),c->get_name()+"_"+to_string(i));
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
    if (_model->_objt == minimize) objt = GRB_MINIMIZE;
    else objt = GRB_MAXIMIZE;
        qobj = 0;
        for (auto& it1: _model->_obj.get_lterms()) {
            lterm = 0;
//            idx = it1.second._p->get_id();
            if (it1.second._coef->_is_transposed) {
                auto dim = it1.second._p->get_dim();
                auto idx = it1.second._p->get_id();
                for (int j = 0; j<dim; j++) {
                    coeff = poly_eval(it1.second._coef,j);
                    gvar1 = _grb_vars[idx + it1.second._p->get_id_inst(j)];
                    lterm += coeff*gvar1;
                }
            }
            else {
                coeff = poly_eval(it1.second._coef);
                gvar1 = _grb_vars[it1.second._p->get_id() + it1.second._p->get_id_inst()];
                lterm += coeff*gvar1;
            }
            if (!it1.second._sign) {
                lterm *= -1;
            }
            qobj += lterm;
        }
    for (auto& it1: _model->_obj.get_qterms()) {
//        idx = it1.second._p->first->get_id();
        gvar1 = _grb_vars[it1.second._p->first->get_id() + it1.second._p->first->get_id_inst()];
//        idx = it1.second._p->second->get_id();
        gvar2 = _grb_vars[it1.second._p->second->get_id() + it1.second._p->second->get_id_inst()];
        if (it1.second._coef->_is_transposed) {
            throw invalid_argument("unsupported operation");
        }
        else {
            coeff = poly_eval(it1.second._coef);
        }
        if (!it1.second._sign) {
            qobj += -1*coeff*gvar1*gvar2;
        }
        else {
            qobj += coeff*gvar1*gvar2;
        }
    }
    qobj += poly_eval(_model->_obj.get_cst());
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
