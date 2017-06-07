#include <Gravity/GurobiProgram.h>

GurobiProgram::GurobiProgram(){
//    model = m;
    grb_env = new GRBEnv();
//    grb_env->set(GRB_IntParam_Presolve,0);
    //grb_env->set(GRB_DoubleParam_NodeLimit,1);
    grb_env->set(GRB_DoubleParam_TimeLimit,7200);
    grb_env->set(GRB_DoubleParam_MIPGap,0.01);
//    grb_env->set(GRB_IntParam_Threads,1);
    grb_env->set(GRB_IntParam_OutputFlag,0);
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
//    grb_mod->write("~/Gqc_ots.lp");
    if (grb_mod->get(GRB_IntAttr_Status) != 2) {
        cerr << "\nModel has not been solved to optimality, error code = " << grb_mod->get(GRB_IntAttr_Status) << endl;
        return false;
    }
    update_model();
    //grb_mod->write("/home/kbestuzheva/PowerTools--/qc_ots.lp");
    GRBVar* gvars = grb_mod->getVars();
    for(int i = 0; i < grb_mod->get(GRB_IntAttr_NumVars); ++i) {
//        cout << gvars[i].get(GRB_StringAttr_VarName) << "  " << gvars[i].get(GRB_DoubleAttr_X) << endl;
        if (gvars[i].get(GRB_CharAttr_VType)==GRB_BINARY) {
            cout << gvars[i].get(GRB_StringAttr_VarName) << "  ";
            cout << gvars[i].get(GRB_DoubleAttr_X);
            cout << "\n";
        }
    }
//    cout << "\n***** Optimal Objective = " << grb_mod->get(GRB_DoubleAttr_ObjVal) << " *****\n";
//    model->_opt = grb_mod->get(GRB_DoubleAttr_ObjVal);
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
    var<bool>* vb;
    var<double>* vr;
    var<int>* vi;
    std::map<int, GRBVar*>::iterator gvarit;
    for (auto& it: _model->_vars){
//        switch (it->get_type()) {
//            case binary:
//                vb =  (var<bool>*)it;
//                gvarit = _grb_vars.find(it->get_idx());
//                vb->set_val(static_cast<int>(gvarit->second->get(GRB_DoubleAttr_X) + 0.5) == 1);
//                break;
//            case integ:
//                vi =  (var<int>*)it;
//                gvarit = _grb_vars.find(it->get_idx());
//                vi->set_val(static_cast<int>(gvarit->second->get(GRB_DoubleAttr_X) + 0.5));
//                break;
//            case real:
//            case longreal:
//                vr =  (var<double>*)it;
//                gvarit = _grb_vars.find(it->get_idx());
//                vr->set_val(gvarit->second->get(GRB_DoubleAttr_X));
////                cout << "\ngrbvar name = " << gvarit->second->get(GRB_StringAttr_VarName) << " grbvar value = " << gvarit->second->get(GRB_DoubleAttr_X);
//                break;
//            default:
//                break;
//        }
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
    for(auto& v_p: _model->_vars)
    {
        v = v_p.second;
        switch (v->get_intype()) {
            case float_: {
                auto real_var = (var<float>*)v;
                for (int i = 0; i < real_var->get_dim(); i++) {
                    _grb_vars.push_back(GRBVar(grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_CONTINUOUS, v->get_name()+"_"+to_string(i))));
                }
                break;
            }
            case long_:{
                auto real_var = (var<long double>*)v;
                for (int i = 0; i < real_var->get_dim(); i++) {
                    _grb_vars.push_back(GRBVar(grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_CONTINUOUS, v->get_name()+"_"+to_string(i))));
                }
                break;
            }
            case double_:{
                auto real_var = (var<double>*)v;
                for (int i = 0; i < real_var->get_dim(); i++) {
                    _grb_vars.push_back(GRBVar(grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_CONTINUOUS, v->get_name()+"_"+to_string(i))));
                }
                break;
            }
            case integer_:{
                auto real_var = (var<int>*)v;
                for (int i = 0; i < real_var->get_dim(); i++) {
                    _grb_vars.push_back(GRBVar(grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_INTEGER, v->get_name()+"_"+to_string(i))));
                }
                break;
            }
            case short_:{
                auto real_var = (var<short>*)v;
                for (int i = 0; i < real_var->get_dim(); i++) {
                    _grb_vars.push_back(GRBVar(grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_INTEGER, v->get_name()+"_"+to_string(i))));
                }
                break;
            }
            case binary_:{
                auto real_var = (var<bool>*)v;
                for (int i = 0; i < real_var->get_dim(); i++) {
                    _grb_vars.push_back(GRBVar(grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_BINARY, v->get_name()+"_"+to_string(i))));
                }
                break;
            }
            default:
                break;
        }
    }
    grb_mod->update();
}

void GurobiProgram::create_grb_constraints(){
    char sense;
    size_t idx = 0, inst = 0, nb_inst = 0;
    GRBLinExpr linlhs;
    GRBLinExpr lterm;
    GRBQuadExpr quadlhs;
    GRBVar gvar1, gvar2;
    double coeff;
    Constraint* c;
    for(auto& p: _model->_cons){
        c = p.second;
        
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
        switch(c->get_ftype()) {
            case lin_:
                inst = 0;
                for (int i = 0; i< nb_inst; i++){
                    linlhs = 0;
                    for (auto& it1: c->get_lterms()) {
                        lterm = 0;
                        idx = it1.second._p->get_id();
                        gvar1 = _grb_vars[idx];
                        if (it1.second._coef->_is_transposed) {
                            for (int j = 0; j<it1.second._p->get_dim(); j++) {
                                coeff = poly_eval(it1.second._coef,j);
                                gvar1 = _grb_vars[idx+j];
                                lterm += coeff*gvar1;
                            }
                        }
                        else {
                            coeff = poly_eval(it1.second._coef,i);
                            gvar1 = _grb_vars[idx+i];
                            lterm += coeff*gvar1;
                        }
                        if (!it1.second._sign) {
                            lterm *= -1;
                        }
                        linlhs += lterm;
                    }
                    linlhs += poly_eval(c->get_cst(), inst);
                    cout << c->get_name() << linlhs << endl;
                    grb_mod->addConstr(linlhs,sense,c->get_rhs(),c->get_name());
                }
                break;
            case quad_:
//                q = c->get_quad();
                quadlhs = 0;
//                for (auto& it1: q->_qmatrix) {
//                    auto gvit = _grb_vars.find(it1.first);
//                    gvar1 = *(gvit->second);
//                    for (auto& it2: *(it1.second)) {
//                        gvit = _grb_vars.find(it2.first);
//                        gvar2 = *(gvit->second);
//                        coeff = it2.second;
//                        quadlhs += gvar1*gvar2*coeff;
//                    }
//                }
//                for (auto& it1: q->_coefs) {
//                    auto gvit = _grb_vars.find(it1.first);
//                    gvar1 = *(gvit->second);
//                    coeff = it1.second;
//                    quadlhs += gvar1*coeff;
//                }
//                quadlhs += q->get_const();
//                cout << c->get_name() << quadlhs << endl;
                grb_mod->addQConstr(quadlhs,sense,c->get_rhs(),c->get_name());
                break;
            case nlin_:
                cerr << "Gurobi cannot handle nonlinear constraints that are not convex quadratic.\n";
                exit(1);
            default:
                break;
        }
    }
}

void GurobiProgram::set_grb_objective(){
    GRBLinExpr lobj;
    GRBQuadExpr qobj;
    GRBVar gvar1, gvar2;
    int objt;
    double coeff;
    const func_* q;
    if (_model->_objt == minimize) objt = GRB_MINIMIZE;
    else objt = GRB_MAXIMIZE;
    switch(_model->_obj.get_type()){
        case lin_:
//            q = model->_obj->get_quad();
//            for (auto& it1: q->_coefs) {
//                auto gvit = _grb_vars.find(it1.first);
//                gvar1 = *(gvit->second);
//                coeff = it1.second;
//                lobj += gvar1*coeff;
//            }
//            lobj += q->get_const();
            grb_mod->setObjective(lobj,objt);
            break;
        case quad_:
//            q = model->_obj->get_quad();
//            for (auto& it1: q->_qmatrix) {
//                auto gvit = _grb_vars.find(it1.first);
//                gvar1 = *(gvit->second);
//                for (auto& it2: *(it1.second)) {
//                   gvit = _grb_vars.find(it2.first);
//                   gvar2 = *(gvit->second);
//                   coeff = it2.second;
//                   qobj += gvar1*gvar2*coeff;
////                   cout << coeff << "*" << gvar1.get(GRB_StringAttr_VarName) << "*" << gvar2.get(GRB_StringAttr_VarName) << " + ";
//                 }
//             }
//             for (auto& it1: q->_coefs) {
//                 auto gvit = _grb_vars.find(it1.first);
//                 gvar1 = *(gvit->second);
//                 coeff = it1.second;
//                 qobj += gvar1*coeff;
////                 cout << coeff << "*" << gvar1.get(GRB_StringAttr_VarName) << " + ";
//             }
////            cout << "\n";
//            qobj += model->_obj->get_const();
            grb_mod->setObjective(qobj,objt);
            break;
        default:
            break;
    }
}

void GurobiProgram::print_constraints(){
    GRBConstr* gconstrs = grb_mod->getConstrs();
    for(int i = 0; i < grb_mod->get(GRB_IntAttr_NumConstrs); ++i) {
        if(gconstrs[i].get(GRB_CharAttr_Sense)!='=') {
            cout << gconstrs[i].get(GRB_StringAttr_ConstrName) << "  ";
            cout << gconstrs[i].get(GRB_DoubleAttr_Slack);
            cout << "\n";
        }
    }
    GRBQConstr* gqconstrs = grb_mod->getQConstrs();
    for(int i = 0; i < grb_mod->get(GRB_IntAttr_NumQConstrs); ++i) {
        if(gqconstrs[i].get(GRB_CharAttr_QCSense)!='=') {
            cout << gqconstrs[i].get(GRB_StringAttr_QCName) << "  ";
            cout << gqconstrs[i].get(GRB_DoubleAttr_QCSlack);
            cout << "\n";
        }
    }
}
