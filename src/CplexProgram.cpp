#include <gravity/CplexProgram.h>
CplexProgram::CplexProgram(Model* m) {
    _cplex_env = new IloEnv();
    _cplex_model = new IloModel(*_cplex_env);
    _model = m;
    m->fill_in_maps();
    m->compute_funcs();
}

void CplexProgram::update_model(){
    _model->fill_in_maps();
    _model->reset_funcs();
    _model->compute_funcs();
    fill_in_cplex_vars();
    create_cplex_constraints();
    set_cplex_objective();
}

CplexProgram::~CplexProgram() {
    delete _cplex_model;
    delete _cplex_env;
}

bool CplexProgram::solve(bool relax) {
    //cout << "\n Presolve = " << grb_env->get(GRB_IntParam_Presolve) << endl;
    //    print_constraints();
    //if (relax) relax_model();
    //    relax_model();
    try {
        IloCplex cplex(*_cplex_env);

        if(relax) {
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

        // set the optimal value.
        _model->_obj_val = cplex.getObjValue();

        for (auto i = 0; i < _cplex_vars.size(); i++) {
            for (auto j = 0; j < _model->_vars[i]->get_dim(); j++) {
                if(cplex.isExtracted(_cplex_vars[i][j])){
                    poly_set_val(j, cplex.getValue(_cplex_vars[i][j]), _model->_vars[i]);
                }
                else {
                    poly_set_val(j, 0, _model->_vars[i]);
                }
            }
        }
        // populate the dual vars.
        unsigned idx = 0; 
        for (auto &cp: _model->_cons) {
            cp.second->_dual.resize(cp.second->_nb_instances);
            assert(_cplex_constraints[idx].getSize() == cp.second->_nb_instances);
            for (unsigned inst = 0; inst < cp.second->_nb_instances; inst++) {
                if (!*cp.second->_all_lazy || !cp.second->_lazy[inst]) {
                    if (cplex.isExtracted(_cplex_constraints[idx][inst])){
                        cp.second->_dual[inst] = cplex.getDual(_cplex_constraints[idx][inst]);
                    }
                    else{
                        cp.second->_dual[inst] = 0;
                    }
                }
            }
            idx++;
        }
    }
    catch (IloException& ex) {
        cerr << "Error: " << ex << endl;
    }
    catch (...) {
        cerr << "Error" << endl;
    }
//    _cplex_env->end();
    return 0;
}

void CplexProgram::fill_in_cplex_vars() {
    _cplex_vars.resize(_model->_vars.size());
    param_* v;
    unsigned vid = -1;
    for(auto& v_p: _model->_vars)
    {
        v = v_p.second;
        if (!v->_new) {
            continue;//Variable already added to the program
        }
        v->_new = false;
        vid = v->get_vec_id();
        if( vid == -1) {
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
            _cplex_vars.at(vid) = IloNumVarArray(*_cplex_env,lb,ub);
            break;
        }
        case long_: {
            auto real_var = (var<long double>*)v;
            auto lb = IloNumArray(*_cplex_env, real_var->get_dim());
            auto ub = IloNumArray(*_cplex_env, real_var->get_dim());
            for (int i = 0; i < real_var->get_dim(); i++) {
                lb[i] = real_var->get_lb(i);
                ub[i] = real_var->get_ub(i);
            }
            _cplex_vars.at(vid) = IloNumVarArray(*_cplex_env,lb,ub);
            break;
        }
        case double_: {
            auto real_var = (var<double>*)v;
            auto lb = IloNumArray(*_cplex_env, real_var->get_dim());
            auto ub = IloNumArray(*_cplex_env, real_var->get_dim());
            for (int i = 0; i < real_var->get_dim(); i++) {
                lb[i] = real_var->get_lb(i);
                ub[i] = real_var->get_ub(i);
            }
            _cplex_vars.at(vid) = IloNumVarArray(*_cplex_env,lb,ub);
            break;
        }
        case integer_: {
            auto real_var = (var<int>*)v;
            auto lb = IloNumArray(*_cplex_env, real_var->get_dim());
            auto ub = IloNumArray(*_cplex_env, real_var->get_dim());
            for (int i = 0; i < real_var->get_dim(); i++) {
                lb[i] = real_var->get_lb(i);
                ub[i] = real_var->get_ub(i);
            }
            _cplex_vars.at(vid) = IloNumVarArray(*_cplex_env,lb,ub, ILOINT);
            break;
        }
        case short_: {
            auto real_var = (var<short>*)v;
            auto lb = IloNumArray(*_cplex_env, real_var->get_dim());
            auto ub = IloNumArray(*_cplex_env, real_var->get_dim());
            for (int i = 0; i < real_var->get_dim(); i++) {
                lb[i] = real_var->get_lb(i);
                ub[i] = real_var->get_ub(i);
            }
            _cplex_vars.at(vid) = IloNumVarArray(*_cplex_env,lb,ub, ILOINT);
            break;
        }
        case binary_: {
            auto real_var = (var<bool>*)v;
            auto lb = IloNumArray(*_cplex_env, real_var->get_dim());
            auto ub = IloNumArray(*_cplex_env, real_var->get_dim());
            for (int i = 0; i < real_var->get_dim(); i++) {
                lb[i] = real_var->get_lb(i);
                ub[i] = real_var->get_ub(i);
            }
            _cplex_vars.at(vid) = IloNumVarArray(*_cplex_env,lb,ub, ILOINT);
            break;
        }
        default:
            break;
        }
    }
}

void CplexProgram::set_cplex_objective() {
    if (!_model->_obj._new) {
        return;//Objective already added to the program
    }
    _model->_obj._new = false;
    size_t idx = 0, idx_inst = 0, idx1 = 0, idx2 = 0, idx_inst1 = 0, idx_inst2 = 0;
    IloNumExpr obj(*_cplex_env);
    for (auto& it_qterm: _model->_obj.get_qterms()) {
        IloNumExpr qterm(*_cplex_env);
        idx1 = it_qterm.second._p->first->get_vec_id();
        idx2 = it_qterm.second._p->second->get_vec_id();
        if (it_qterm.second._coef->_is_transposed) {
            auto dim = it_qterm.second._p->first->get_dim();
            for (int j = 0; j<dim; j++) {
                qterm += poly_eval(it_qterm.second._coef,j)*_cplex_vars[idx1][it_qterm.second._p->first->get_id_inst(j)]*_cplex_vars[idx2][it_qterm.second._p->second->get_id_inst(j)];
            }            
        }
        else {
            idx_inst1 = it_qterm.second._p->first->get_id_inst();
            idx_inst2 = it_qterm.second._p->second->get_id_inst();
            qterm += poly_eval(it_qterm.second._coef)*_cplex_vars[idx1][idx_inst1]*_cplex_vars[idx2][idx_inst2];
        }
        if (!it_qterm.second._sign) {
            qterm *= -1;
        }
        obj += qterm;
        qterm.end();
    }

    for (auto& it_lterm: _model->_obj.get_lterms()) {
        IloNumExpr lterm(*_cplex_env);
        idx = it_lterm.second._p->get_vec_id();
        if (it_lterm.second._coef->_is_transposed) {
            auto dim = it_lterm.second._p->get_dim();
            for (int j = 0; j<dim; j++) {
                lterm += poly_eval(it_lterm.second._coef,j)*_cplex_vars[idx][it_lterm.second._p->get_id_inst(j)];
            }
        }
        else {
            idx_inst = it_lterm.second._p->get_id_inst();
            lterm += poly_eval(it_lterm.second._coef)*_cplex_vars[idx][idx_inst];
        }
        if (!it_lterm.second._sign) {
            lterm *= -1;
        }
        obj += lterm;
        lterm.end();
    }

    obj += poly_eval(_model->_obj.get_cst());

    if (_model->_objt == maximize) {
        _cplex_obj = IloMaximize(*_cplex_env,obj);
    }
    else {
        _cplex_obj = IloMinimize(*_cplex_env,obj);
    }
    _cplex_model->add(_cplex_obj);
    obj.end();
}

void CplexProgram::create_cplex_constraints() {
    size_t idx = 0, idx_inst = 0, idx1 = 0, idx2 = 0, idx_inst1 = 0, idx_inst2 = 0, nb_inst = 0, inst = 0;
//    size_t c_idx_inst = 0;
    Constraint* c;
    for(auto& p: _model->_cons) {
        c = p.second.get();
        if (!c->_new) {
            continue;//Constraint already added to the program
        }
        c->_new = false;
        if (c->is_nonlinear()) {
            throw invalid_argument("Cplex cannot handle nonlinear constraints that are not convex quadratic.\n");
        }
        nb_inst = c->_nb_instances;
        IloRangeArray _cplex_temp(*_cplex_env, nb_inst);
        inst = 0;
        for (int i = 0; i< nb_inst; i++) {
            IloNumExpr cc(*_cplex_env);
            for (auto& it_qterm: c->get_qterms()) {
                IloNumExpr qterm(*_cplex_env);
                idx1 = it_qterm.second._p->first->get_vec_id();
                idx2 = it_qterm.second._p->second->get_vec_id();
                if (it_qterm.second._coef->_is_transposed) {
                    auto dim = it_qterm.second._p->first->get_dim(i);
                    for (int j = 0; j<dim; j++) {
                        qterm += poly_eval(it_qterm.second._coef,i,j)*_cplex_vars[idx1][it_qterm.second._p->first->get_id_inst(i,j)]*_cplex_vars[idx2][it_qterm.second._p->second->get_id_inst(i,j)];
                    }
                }
                else {                    
                    idx_inst1 = it_qterm.second._p->first->get_id_inst(inst);
                    idx_inst2 = it_qterm.second._p->second->get_id_inst(inst);
                    qterm += poly_eval(it_qterm.second._coef, inst)*_cplex_vars[idx1][idx_inst1]*_cplex_vars[idx2][idx_inst2];
                }
                if (!it_qterm.second._sign) {
                    qterm *= -1;
                }
                cc += qterm;
                qterm.end();
            }
            for (auto& it_lterm: c->get_lterms()) {
                IloNumExpr lterm(*_cplex_env);
                idx = it_lterm.second._p->get_vec_id();
                if (it_lterm.second._coef->_is_transposed) {
                    auto dim = it_lterm.second._p->get_dim(i);
                    for (int j = 0; j < dim; j++) {
                        Debug("polyeval: " << i << ","  << j << " "<< poly_eval(it_lterm.second._coef,i,j) <<endl);
                        Debug("_ids size " << it_lterm.second._p->_ids->size() << endl);
                        Debug("cplex_var: " <<it_lterm.second._p->get_id_inst(i,j) << endl);
                        lterm += poly_eval(it_lterm.second._coef,i,j)*_cplex_vars[idx][it_lterm.second._p->get_id_inst(i,j)];
                    }                    
                }
                else {
                    idx_inst = it_lterm.second._p->get_id_inst(inst);
                    lterm += poly_eval(it_lterm.second._coef, inst)*_cplex_vars[idx][idx_inst];
                }
                if (!it_lterm.second._sign) {
                    lterm *= -1;
                }
                cc += lterm;
                lterm.end();
            }
            cc += poly_eval(c->get_cst(), inst);

            if(c->get_type()==geq) {
                //_cplex_model->add(cc >= c->get_rhs());
                _cplex_temp[i] = (cc >= c->get_rhs());
            }
            else if(c->get_type()==leq) {
                //_cplex_model->add(cc <= c->get_rhs());
                _cplex_temp[i] = cc <= c->get_rhs();
            }
            else if(c->get_type()==eq) {
                //_cplex_model->add(cc == c->get_rhs());
                _cplex_temp[i]= (cc == c->get_rhs());
            }
            inst++;
        }
        _cplex_model->add(_cplex_temp);
        _cplex_constraints.push_back(_cplex_temp);
    }    
    _cplex_constraints.resize(_model->_cons.size());
}

void CplexProgram::prepare_model() {
    fill_in_cplex_vars();
    create_cplex_constraints();
    set_cplex_objective();
    //    print_constraints();
}
