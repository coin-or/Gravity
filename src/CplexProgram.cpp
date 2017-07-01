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
//        _cplex_model->M
        IloCplex cplex(*_cplex_env);
        
//        cplex.setOut(_cplex_env->getNullStream());
//        cplex.setWarning(_cplex_env->getNullStream());
//        _cplex_model->setIntProperty("Threads", 1);
//        cplex.setParam(IloCplex::Threads, 1);
//        cplex.setParam(IloCplex::EpRHS, 1e-8);
//        cplex.setParam(IloCplex::EpOpt, 1e-8);
//        cplex.setParam(IloCplex::EpInt, 1e-8);
        if(relax){
            IloModel relax(*_cplex_env);
            relax.add(*_cplex_model);
            for (auto &vv: _cplex_vars) {
                relax.add(IloConversion(*_cplex_env, vv, ILOFLOAT));
            }
            cplex.extract(relax);
        }
        else {
            cplex.extract(*_cplex_model);
        }
        cplex.solve();
        if (cplex.getStatus() == IloAlgorithm::Infeasible)
            _cplex_env->out() << "No Solution" << endl;
        
        _cplex_env->out() << "Solution status: " << cplex.getStatus() << endl;
        
        // Print results
        _cplex_env->out() << "Cost:" << cplex.getObjValue() << endl;
        for (auto i = 0; i < _cplex_vars.size(); i++) {
            IloNumArray vals(*_cplex_env);
            cplex.getValues(vals,_cplex_vars[i]);
//            cout << "  x" << i << " = " <<  vals << endl;
            for (auto j = 0; j < _model->_vars[i]->get_dim(); j++) {
                poly_set_val(j, vals[j], _model->_vars[i]);
            }
        }
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
        if(v->get_vec_id()==-1){
            throw invalid_argument("Variable needs to be added to model first: use add_var(v) function:" + v->get_name());
        }
        switch (v->get_intype()) {
            case float_: {
                auto real_var = (var<float>*)v;
                auto lb = IloNumArray(*_cplex_env, real_var->get_dim());
                auto ub = IloNumArray(*_cplex_env, real_var->get_dim());
                for (int i = 0; i < real_var->get_dim(); i++) {
                    lb[i] = real_var->get_lb(i);
                    ub[i] = real_var->get_ub(i);
                }
                _cplex_vars.push_back(IloNumVarArray(*_cplex_env,lb,ub));
                break;
            }
            case long_:{
                auto real_var = (var<long double>*)v;
                auto lb = IloNumArray(*_cplex_env, real_var->get_dim());
                auto ub = IloNumArray(*_cplex_env, real_var->get_dim());
                for (int i = 0; i < real_var->get_dim(); i++) {
                    lb[i] = real_var->get_lb(i);
                    ub[i] = real_var->get_ub(i);
                }
                _cplex_vars.push_back(IloNumVarArray(*_cplex_env,lb,ub));
                break;
            }
            case double_:{
                auto real_var = (var<double>*)v;
                auto lb = IloNumArray(*_cplex_env, real_var->get_dim());
                auto ub = IloNumArray(*_cplex_env, real_var->get_dim());
                for (int i = 0; i < real_var->get_dim(); i++) {
                    lb[i] = real_var->get_lb(i);
                    ub[i] = real_var->get_ub(i);
                }
                _cplex_vars.push_back(IloNumVarArray(*_cplex_env,lb,ub));
                break;
            }
            case integer_:{
                auto real_var = (var<int>*)v;
                auto lb = IloNumArray(*_cplex_env, real_var->get_dim());
                auto ub = IloNumArray(*_cplex_env, real_var->get_dim());
                for (int i = 0; i < real_var->get_dim(); i++) {
                    lb[i] = real_var->get_lb(i);
                    ub[i] = real_var->get_ub(i);
                }
                _cplex_vars.push_back(IloNumVarArray(*_cplex_env,lb,ub, ILOINT));
                break;
            }
            case short_:{
                auto real_var = (var<short>*)v;
                auto lb = IloNumArray(*_cplex_env, real_var->get_dim());
                auto ub = IloNumArray(*_cplex_env, real_var->get_dim());
                for (int i = 0; i < real_var->get_dim(); i++) {
                    lb[i] = real_var->get_lb(i);
                    ub[i] = real_var->get_ub(i);
                }
                _cplex_vars.push_back(IloNumVarArray(*_cplex_env,lb,ub, ILOINT));
                break;
            }
            case binary_:{
                auto real_var = (var<bool>*)v;
                auto lb = IloNumArray(*_cplex_env, real_var->get_dim());
                auto ub = IloNumArray(*_cplex_env, real_var->get_dim());
                for (int i = 0; i < real_var->get_dim(); i++) {
                    lb[i] = real_var->get_lb(i);
                    ub[i] = real_var->get_ub(i);
                }
                _cplex_vars.push_back(IloNumVarArray(*_cplex_env,lb,ub, ILOBOOL));
                break;
            }
            default:
                break;
        }
    }    
}



void CplexProgram::set_cplex_objective(){
    size_t idx = 0, idx_inst = 0, idx1 = 0, idx2 = 0, idx_inst1 = 0, idx_inst2 = 0;
    size_t c_idx_inst = 0;
    IloNumExpr obj(*_cplex_env);
    for (auto& it1: _model->_obj.get_qterms()) {
        idx1 = it1.second._p->first->get_vec_id();
        idx2 = it1.second._p->second->get_vec_id();
        IloNumExpr lterm(*_cplex_env);
        if (it1.second._coef->_is_transposed) {
            IloNumArray coefs(*_cplex_env,it1.second._p->first->get_dim());
            for (int j = 0; j<it1.second._p->first->get_dim(); j++) {
                coefs[j] = poly_eval(it1.second._coef,j);
            }
            obj += IloQuadProd(_cplex_vars[idx1], _cplex_vars[idx2], coefs);
        }
        else {
            IloNumExpr qterm(*_cplex_env);
            idx_inst1 = it1.second._p->first->get_id_inst();
            idx_inst2 = it1.second._p->second->get_id_inst();
            c_idx_inst = get_poly_id_inst(it1.second._coef);
            qterm += poly_eval(it1.second._coef, c_idx_inst)*_cplex_vars[idx1][idx_inst1]*_cplex_vars[idx2][idx_inst2];
            if (!it1.second._sign) {
                qterm *= -1;
            }
            obj += qterm;
            qterm.end();
        }
    }
    
    for (auto& it1: _model->_obj.get_lterms()) {
        idx = it1.second._p->get_vec_id();
        if (it1.second._coef->_is_transposed) {
            IloNumArray coefs(*_cplex_env,it1.second._p->get_dim());
            for (int j = 0; j<it1.second._p->get_dim(); j++) {
                coefs[j] = poly_eval(it1.second._coef,j);
            }
            obj += IloScalProd(coefs, _cplex_vars[idx]);
        }
        else {
            IloNumExpr lterm(*_cplex_env);
            idx_inst = it1.second._p->get_id_inst();
            c_idx_inst = get_poly_id_inst(it1.second._coef);
            lterm += poly_eval(it1.second._coef, c_idx_inst)*_cplex_vars[idx][idx_inst];
            if (!it1.second._sign) {
                lterm *= -1;
            }
            obj += lterm;
            lterm.end();
        }
    }
    obj += poly_eval(_model->_obj.get_cst());
    
    if (_model->_objt == maximize){
        _cplex_obj = IloMaximize(*_cplex_env,obj);
    }
    else {
        _cplex_obj = IloMinimize(*_cplex_env,obj);
    }
    _cplex_model->add(_cplex_obj);
    obj.end();
}



void CplexProgram::create_cplex_constraints(){
    size_t idx = 0, idx_inst = 0, idx1 = 0, idx2 = 0, idx_inst1 = 0, idx_inst2 = 0, nb_inst = 0, inst = 0;
    size_t c_idx_inst = 0;
    Constraint* c;
    for(auto& p: _model->_cons){
        c = p.second;
//        c->print();
        if (c->is_nonlinear()) {
            throw invalid_argument("Cplex cannot handle nonlinear constraints that are not convex quadratic.\n");
        }
        nb_inst = c->get_nb_instances();
        inst = 0;
        for (int i = 0; i< nb_inst; i++){
            IloNumExpr cc(*_cplex_env);
            for (auto& it1: c->get_qterms()) {
                idx1 = it1.second._p->first->get_vec_id();
                idx2 = it1.second._p->second->get_vec_id();
                IloNumExpr lterm(*_cplex_env);
                if (it1.second._coef->_is_transposed) {
                    IloNumArray coefs(*_cplex_env,it1.second._p->first->get_dim());
                    for (int j = 0; j<it1.second._p->first->get_dim(); j++) {
                        coefs[j] = poly_eval(it1.second._coef,j);
                    }
                    cc += IloQuadProd(_cplex_vars[idx1], _cplex_vars[idx2], coefs);
                }
                else {                    
                    IloNumExpr qterm(*_cplex_env);
                    if (is_indexed(it1.second._p->first)) {
                        idx_inst1 = it1.second._p->first->get_id_inst();
                    }
                    else {
                        idx_inst1 = inst;
                    }
                    if (is_indexed(it1.second._p->second)) {
                        idx_inst2 = it1.second._p->second->get_id_inst();
                    }
                    else {
                        idx_inst2 = inst;
                    }
                    if (is_indexed(it1.second._coef)) {
                        c_idx_inst = get_poly_id_inst(it1.second._coef);
                    }
                    else {
                        c_idx_inst = inst;
                    }
                    qterm += poly_eval(it1.second._coef, c_idx_inst)*_cplex_vars[idx1][idx_inst1]*_cplex_vars[idx2][idx_inst2];
                    if (!it1.second._sign) {
                        qterm *= -1;
                    }
                    cc += qterm;
                    qterm.end();
                }
            }
            
            for (auto& it1: c->get_lterms()) {
                idx = it1.second._p->get_vec_id();
                if (it1.second._coef->_is_transposed) {
                    IloNumArray coefs(*_cplex_env,it1.second._p->get_dim());
                    for (int j = 0; j<it1.second._p->get_dim(); j++) {
                        coefs[j] = poly_eval(it1.second._coef,j);
                    }
                    cc += IloScalProd(coefs, _cplex_vars[idx]);
                }
                else {
                    IloNumExpr lterm(*_cplex_env);
                    if (is_indexed(it1.second._p)) {
                        idx_inst = it1.second._p->get_id_inst();
                    }
                    else {
                        idx_inst = inst;
                    }
                    if (is_indexed(it1.second._coef)) {
                        c_idx_inst = get_poly_id_inst(it1.second._coef);
                    }
                    else {
                        c_idx_inst = inst;
                    }
                    lterm += poly_eval(it1.second._coef, c_idx_inst)*_cplex_vars[idx][idx_inst];
                    if (!it1.second._sign) {
                        lterm *= -1;
                    }
                    cc += lterm;
                    lterm.end();
                }
            }
            cc += poly_eval(c->get_cst());
            
            if(c->get_type()==geq){
                _cplex_model->add(cc >= c->get_rhs());
            }
            else if(c->get_type()==leq){
                _cplex_model->add(cc <= c->get_rhs());
            }
            else if(c->get_type()==eq){
                _cplex_model->add(cc == c->get_rhs());
            }
            inst++;
        }
    }
}

void CplexProgram::prepare_model(){
    fill_in_cplex_vars();
    create_cplex_constraints();
    set_cplex_objective();    
    //    print_constraints();
}
