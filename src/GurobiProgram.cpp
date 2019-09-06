#include <gravity/GurobiProgram.h>

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


GurobiProgram::GurobiProgram(const shared_ptr<Model<>>& m) {
    grb_env = new GRBEnv();
    grb_mod = new GRBModel(*grb_env);
//    grb_env->set(GRB_IntParam_OutputFlag,2);
    _model = m;
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

bool GurobiProgram::solve(bool relax, double mipgap){
    //cout << "\n Presolve = " << grb_env->get(GRB_IntParam_Presolve) << endl;
//    print_constraints();
    if (relax) relax_model();
//    relax_model();
    grb_mod->set(GRB_DoubleParam_MIPGap, mipgap);
    grb_mod->set(GRB_IntParam_Threads, 1);
    grb_mod->optimize();
//    grb_mod->write("~/mod.mps");
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
    cout << "\n***** Optimal Objective = " << grb_mod->get(GRB_DoubleAttr_ObjVal) << " *****\n";
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
    fill_in_grb_vmap();
    create_grb_constraints();
    set_grb_objective();
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
        auto idx = v->get_vec_id();
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
        v->_new = false;
        auto idx = v->get_vec_id();
        switch (v->get_intype()) {
            case float_: {
                auto real_var = (var<float>*)v;
                for (int i = 0; i < real_var->_dim[0]; i++) {
                    auto vid = idx + i;
                    _grb_vars.at(vid) = (GRBVar(grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_CONTINUOUS, v->get_name(true,true)+"_"+to_string(i))));
                }
                break;
            }
            case long_:{
                auto real_var = (var<long double>*)v;
                for (int i = 0; i < real_var->_dim[0]; i++) {
                    auto vid = idx + i;
                    _grb_vars.at(vid) = (GRBVar(grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_CONTINUOUS, v->get_name(true,true)+"_"+to_string(i))));
                }
                break;
            }
            case double_:{
                auto real_var = (var<double>*)v;
                for (size_t i = 0; i < real_var->_dim[0]; i++) {
                    auto vid = idx + i;
                    _grb_vars.at(vid) = (grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_CONTINUOUS, v->get_name(true,true)+"_"+to_string(i)));
                }
                break;
            }
            case integer_:{
                auto real_var = (var<int>*)v;
                for (int i = 0; i < real_var->_dim[0]; i++) {
                    auto vid = idx + i;
                    _grb_vars.at(vid) = (GRBVar(grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_INTEGER, v->get_name(true,true)+"_"+to_string(i))));
                }
                break;
            }
            case short_:{
                auto real_var = (var<short>*)v;
                for (int i = 0; i < real_var->_dim[0]; i++) {
                    auto vid = idx + i;
                    _grb_vars.at(vid) = (GRBVar(grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_INTEGER, v->get_name(true,true)+"_"+to_string(i))));
                }
                break;
            }
            case binary_:{
                auto real_var = (var<bool>*)v;
                for (int i = 0; i < real_var->_dim[0]; i++) {
                    auto vid = idx + i;
                    _grb_vars.at(vid) = (GRBVar(grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_BINARY, v->get_name(true,true)+"_"+to_string(i))));
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
        if (c->is_nonlinear()) {
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
        nb_inst = c->_dim[0];
        inst = 0;
        if (c->is_linear()) {
            for (size_t i = 0; i< nb_inst; i++){
//                if (c->_violated[i]) {
                    linlhs = 0;
                    for (auto& it1: c->get_lterms()) {
                        lterm = 0;
                        if (it1.second._coef->_is_transposed) {
                            auto dim =it1.second._p->get_dim(i);
                            for (int j = 0; j<dim; j++) {
                                coeff = c->eval(it1.second._coef,i,j);
                                gvar1 = _grb_vars[it1.second._p->get_vec_id() + it1.second._p->get_id_inst(i,j)];
                                lterm += coeff*gvar1;
                            }
                        }
                        else {
                            coeff = c->eval(it1.second._coef,i);
                            gvar1 = _grb_vars[it1.second._p->get_vec_id() + it1.second._p->get_id_inst(i)];
                            lterm += coeff*gvar1;
                        }
                        if (!it1.second._sign) {
                            lterm *= -1;
                        }
                        linlhs += lterm;
                    }
                    linlhs += c->eval(c->get_cst(), i);
                    grb_mod->addConstr(linlhs,sense,0,c->get_name()+"_"+to_string(i));
//                }
            }
        }
        else {
            for (size_t i = 0; i< nb_inst; i++){
//                if (c->_violated[i]) {
                    quadlhs = 0;
                    for (auto& it1: c->get_lterms()) {
                        lterm = 0;
                        if (it1.second._p->_is_vector || it1.second._coef->is_matrix()) {
                            auto dim =it1.second._p->get_dim(i);
                            for (size_t j = 0; j<dim; j++) {
                                coeff = c->eval(it1.second._coef,i,j);
                                gvar1 = _grb_vars[it1.second._p->get_vec_id() + it1.second._p->get_id_inst(i,j)];
                                lterm += coeff*gvar1;
                            }
                        }
                        else {
                            coeff = c->eval(it1.second._coef,i);
                            gvar1 = _grb_vars[it1.second._p->get_vec_id() + it1.second._p->get_id_inst(i)];
                            lterm += coeff*gvar1;
                        }
                        if (!it1.second._sign) {
                            lterm *= -1;
                        }
                        quadlhs += lterm;
                    }
                    for (auto& it1: c->get_qterms()) {
                        gvar1 = _grb_vars[it1.second._p->first->get_vec_id() + it1.second._p->first->get_id_inst(i)];
                        gvar2 = _grb_vars[it1.second._p->second->get_vec_id() + it1.second._p->second->get_id_inst(i)];
                        if (it1.second._p->first->_is_vector) {
                            auto dim =it1.second._p->first->get_dim(i);
                            for (int j = 0; j<dim; j++) {
                                coeff = c->eval(it1.second._coef,i,j);
                                gvar1 = _grb_vars[it1.second._p->first->get_vec_id() + it1.second._p->first->get_id_inst(i,j)];
                                gvar2 = _grb_vars[it1.second._p->second->get_vec_id() + it1.second._p->second->get_id_inst(i,j)];
                                if (!it1.second._sign) {
                                    quadlhs += -1*coeff*gvar1*gvar2;
                                }
                                else {
                                    quadlhs += coeff*gvar1*gvar2;
                                }
                            }
                        }
                        else {
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
                    grb_mod->addQConstr(quadlhs,sense,0,c->get_name()+"_"+to_string(i));
                grb_mod->re
//                }
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
//            idx = it1.second._p->get_id();
            if (it1.second._coef->_is_transposed) {
                auto dim = it1.second._p->_dim[0];
                auto idx = it1.second._p->get_vec_id();
                for (int j = 0; j<dim; j++) {
                    coeff = _model->_obj->eval(it1.second._coef,j);
                    gvar1 = _grb_vars[idx + it1.second._p->get_id_inst(j)];
                    lterm += coeff*gvar1;
                }
            }
            else {
                coeff = _model->_obj->eval(it1.second._coef);
                gvar1 = _grb_vars[it1.second._p->get_vec_id() + it1.second._p->get_id_inst()];
                lterm += coeff*gvar1;
            }
            if (!it1.second._sign) {
                lterm *= -1;
            }
            qobj += lterm;
        }
    for (auto& it1: _model->_obj->get_qterms()) {
//        idx = it1.second._p->first->get_id();
        gvar1 = _grb_vars[it1.second._p->first->get_vec_id() + it1.second._p->first->get_id_inst()];
//        idx = it1.second._p->second->get_id();
        gvar2 = _grb_vars[it1.second._p->second->get_vec_id() + it1.second._p->second->get_id_inst()];
        if (it1.second._coef->_is_transposed) {
            auto dim =it1.second._p->first->get_dim();
            for (int j = 0; j<dim; j++) {
                coeff = _model->_obj->eval(it1.second._coef,j);
                gvar1 = _grb_vars[it1.second._p->first->get_vec_id() + it1.second._p->first->get_id_inst(j)];
                gvar2 = _grb_vars[it1.second._p->second->get_vec_id() + it1.second._p->second->get_id_inst(j)];
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
