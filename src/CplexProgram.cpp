#include <Gravity/CplexProgram.h>

CplexProgram::CplexProgram(Model* m){
    _cplex_env = new IloEnv();
    _cplex_model = new IloModel(*_cplex_env);
    _model = m;
}


CplexProgram::~CplexProgram() {    
    delete _cplex_model;
    delete _cplex_env;
}


bool CplexProgram::solve(bool relax){
    //cout << "\n Presolve = " << grb_env->get(GRB_IntParam_Presolve) << endl;
    //    print_constraints();
//    if (relax) relax_model();
    //    relax_model();
    try {
        IloCplex cplex(*_cplex_model);
        cplex.setOut(_cplex_env->getNullStream());
        cplex.setWarning(_cplex_env->getNullStream());
        cplex.solve();
        if (cplex.getStatus() == IloAlgorithm::Infeasible)
            _cplex_env->out() << "No Solution" << endl;
        
        _cplex_env->out() << "Solution status: " << cplex.getStatus() << endl;
        
        // Print results
        _cplex_env->out() << "Cost:" << cplex.getObjValue() << endl;
    }
    catch (IloException& ex) {
        cerr << "Error: " << ex << endl;
    }
    catch (...) {
        cerr << "Error" << endl;
    }
    _cplex_env->end();
    return 0;
}

void CplexProgram::fill_in_cplex_vars(){
    param_* v;
    for(auto& v_p: _model->_vars)
    {
        v = v_p.second;
        switch (v->get_intype()) {
            case float_: {
                auto real_var = (var<float>*)v;
                for (int i = 0; i < real_var->get_dim(); i++) {
                    _cplex_vars.push_back(IloNumVar(*_cplex_env,real_var->get_lb(i), real_var->get_ub(i)));
                }
                break;
            }
            case long_:{
                auto real_var = (var<long double>*)v;
                for (int i = 0; i < real_var->get_dim(); i++) {
                    _cplex_vars.push_back(IloNumVar(*_cplex_env,real_var->get_lb(i), real_var->get_ub(i)));
                }
                break;
            }
            case double_:{
                auto real_var = (var<double>*)v;
                for (int i = 0; i < real_var->get_dim(); i++) {
                    _cplex_vars.push_back(IloNumVar(*_cplex_env,real_var->get_lb(i), real_var->get_ub(i)));
                }
                break;
            }
            case integer_:{
                auto real_var = (var<int>*)v;
                for (int i = 0; i < real_var->get_dim(); i++) {
                    _cplex_vars.push_back(IloNumVar(*_cplex_env,real_var->get_lb(i), real_var->get_ub(i), ILOINT));
                }
                break;
            }
            case short_:{
                auto real_var = (var<short>*)v;
                for (int i = 0; i < real_var->get_dim(); i++) {
                    _cplex_vars.push_back(IloNumVar(*_cplex_env,real_var->get_lb(i), real_var->get_ub(i), ILOINT));
                }
                break;
            }
            case binary_:{
                auto real_var = (var<bool>*)v;
                for (int i = 0; i < real_var->get_dim(); i++) {
                    _cplex_vars.push_back(IloNumVar(*_cplex_env,real_var->get_lb(i), real_var->get_ub(i), ILOBOOL));
                }
                break;
            }
            default:
                break;
        }
    }    
}



void CplexProgram::set_cplex_objective(){
    size_t idx = 0;
    IloExpr lterm(*_cplex_env);
    IloExpr lobj(*_cplex_env);
//    GRBQuadExpr qobj;
    double coeff;
    //    const func_* q;
//    lobj = 0;
    for (auto& it1: _model->_obj.get_lterms()) {
        lterm = 0;
        idx = it1.second._p->get_id();
        if (it1.second._coef->_is_transposed) {
            for (int j = 0; j<it1.second._p->get_dim(); j++) {
                coeff = poly_eval(it1.second._coef,j);
                lterm = lterm + coeff*_cplex_vars[idx+j];
            }
        }
        else {
            coeff = poly_eval(it1.second._coef);
            lterm += coeff*_cplex_vars[idx];
        }
        if (!it1.second._sign) {
            lterm *= -1;
        }
        lobj += lterm;
    }
    lobj += poly_eval(_model->_obj.get_cst());
    cout << lobj << endl;
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
//    grb_mod->setObjective(lobj,objt);
    if (_model->_objt == minimize){
        _cplex_model->add(IloMinimize(*_cplex_env,lobj));
    }
    else {
        _cplex_model->add(IloMaximize(*_cplex_env,lobj));
    }
    lterm.end();
    //    else objt = GRB_MAXIMIZE;

}



void CplexProgram::create_cplex_constraints(){
    size_t idx = 0, inst = 0, nb_inst = 0;
    IloExpr linlhs(*_cplex_env);
    IloExpr lterm(*_cplex_env);
    IloExpr quadlhs(*_cplex_env);
    double coeff;
    Constraint* c;
    for(auto& p: _model->_cons){
        c = p.second;
        if (c->is_nonlinear()) {
            throw invalid_argument("Cplex cannot handle nonlinear constraints that are not convex quadratic.\n");
        }
        nb_inst = c->get_nb_instances();
        inst = 0;
        for (int i = 0; i< nb_inst; i++){
            linlhs = 0;
            for (auto& it1: c->get_lterms()) {
                lterm = 0;
                idx = it1.second._p->get_id();
                if (it1.second._coef->_is_transposed) {
                    for (int j = 0; j<it1.second._p->get_dim(); j++) {
                        coeff = poly_eval(it1.second._coef,j);
                        lterm += coeff*_cplex_vars[idx+j];
                    }
                }
                else {
                    coeff = poly_eval(it1.second._coef,i);
                    lterm += coeff*_cplex_vars[idx+i];
                }
                if (!it1.second._sign) {
                    lterm *= -1;
                }
                linlhs += lterm;
            }
            linlhs += poly_eval(c->get_cst(), inst);
            cout << c->get_name() << linlhs << endl;
            if(c->get_type()==geq){
                _cplex_model->add(linlhs >= c->get_rhs());
            }
            else if(c->get_type()==leq){
                _cplex_model->add(linlhs <= c->get_rhs());
            }
            else if(c->get_type()==eq){
                _cplex_model->add(linlhs == c->get_rhs());
            }
        }
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
//        grb_mod->addQConstr(quadlhs,sense,c->get_rhs(),c->get_name());
    }
}


