#include <gravity/CplexProgram.h>

CplexProgram::CplexProgram(const shared_ptr<Model<>>& m) {
    _cplex_env = make_shared<IloEnv>();
    _cplex_model = make_shared<IloModel>(*_cplex_env);
    _model = m;
    m->fill_in_maps();
    m->compute_funcs();
}

void CplexProgram::update_model(){
    _model->fill_in_maps();
    _model->compute_funcs();
    fill_in_cplex_vars();
    create_cplex_constraints();
    set_cplex_objective();
}


void CplexProgram::warm_start(){
    IloCplex cplex(*_cplex_env);
    IloNumArray vals(*_cplex_env);
    IloNumVarArray vars(*_cplex_env);
    double val = 0;
    for (auto i = 0; i < _cplex_vars.size(); i++) {
        for (auto j = 0; j < _model->_vars[i]->get_dim(); j++) {
            vars.add(_cplex_vars[i][j]);
            _model->_vars[i]->set_double_val(j,val);
            vals.add(val);
        }
    }
    cplex.setStart(vals, 0, vars, 0, 0, 0);
}

bool CplexProgram::solve(bool relax, double mipgap) {
    //cout << "\n Presolve = " << grb_env->get(GRB_IntParam_Presolve) << endl;
    //    print_constraints();
    //if (relax) relax_model();
    //    relax_model();
    int return_status = -1;
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
        cplex.setParam(IloCplex::EpGap, mipgap);
        cplex.setParam(IloCplex::Param::OptimalityTarget, 2);
        cplex.setParam(IloCplex::Param::Threads, 1);
//        cplex.setParam(IloCplex::BarDisplay, 2);
//        cplex.setParam(IloCplex::RootAlg, 1);
//        cplex.setParam(IloCplex::AdvInd, 1);
        cplex.setParam(IloCplex::MIPDisplay, 2);
//        cplex.setParam(IloCplex::SimDisplay, 2);
//        cplex.setParam(IloCplex::PreInd, 0);
        IloNumArray vals(*_cplex_env);
        IloNumVarArray vars(*_cplex_env);
        double val = 0;
        for (auto i = 0; i < _cplex_vars.size(); i++) {
            for (auto j = 0; j < _model->_vars[i]->get_dim(); j++) {
                vars.add(_cplex_vars[i][j]);
                _model->_vars[i]->set_double_val(j,val);
                vals.add(val);
            }
        }
        cplex.setStart(vals, 0, vars, 0, 0, 0);
        cplex.solve();
        if (cplex.getStatus() == IloAlgorithm::Infeasible) {
            _cplex_env->out() << "No Solution" << endl;
        }
        else if(cplex.getStatus() == IloAlgorithm::Optimal){
            return_status = 100;
        }
        _cplex_env->out() << "Solution status: " << cplex.getStatus() << endl;

        // Print results
        _cplex_env->out() << "Cost:" << cplex.getObjValue() << endl;

        // set the optimal value.
        _model->_obj->set_val(cplex.getObjValue());

        for (auto i = 0; i < _cplex_vars.size(); i++) {
            for (auto j = 0; j < _model->_vars[i]->get_dim(); j++) {
                if(cplex.isExtracted(_cplex_vars[i][j])){
                    _model->_vars[i]->get_double_val(j,cplex.getValue(_cplex_vars[i][j]));
                }
                else {
                    _model->_vars[i]->get_double_val(0, j);
                }
            }
        }
    }
    catch (IloException& ex) {
        cerr << "Error: " << ex << endl;
    }
    catch (...) {
        cerr << "Error" << endl;
    }
//    _cplex_env->end();
    return return_status==100;
}

void CplexProgram::fill_in_cplex_vars() {
    _cplex_vars.resize(_model->_vars.size());
    unsigned vid = 0;
    for(auto& v_p: _model->_vars)
    {
        auto v = v_p.second;
//        if (!v->_new) {
//            continue;//Variable already added to the program
//        }
        v->_new = false;
        vid = v->get_vec_id();
        switch (v->get_intype()) {
        case float_: {
            auto real_var = (var<float>*)v.get();
            auto lb = IloNumArray(*_cplex_env, real_var->get_dim());
            auto ub = IloNumArray(*_cplex_env, real_var->get_dim());
            for (auto i = 0; i < real_var->get_dim(); i++) {
                lb[i] = real_var->get_lb(i);
                ub[i] = real_var->get_ub(i);
            }
            _cplex_vars.at(vid) = IloNumVarArray(*_cplex_env,lb,ub);
            break;
        }
        case long_: {
            auto real_var = (var<long double>*)v.get();
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
            auto real_var = (var<double>*)v.get();
            auto lb = IloNumArray(*_cplex_env, real_var->get_dim());
            auto ub = IloNumArray(*_cplex_env, real_var->get_dim());
//            real_var->print();
//            cout << ": ";
            for (int i = 0; i < real_var->get_dim(); i++) {
                lb[i] = real_var->get_lb(i);
                ub[i] = real_var->get_ub(i);
            }
            
            if (real_var->_is_relaxed) {
                _cplex_vars.at(vid) = IloNumVarArray(*_cplex_env,lb,ub, ILOINT);
            }
            else {
                _cplex_vars.at(vid) = IloNumVarArray(*_cplex_env,lb,ub);
            }
//            for (int i = 0; i < real_var->get_dim(); i++) {
//                cout << real_var->_indices->_keys->at(i) << " : ";
//                cout << to_string(_cplex_vars.at(vid)[i].getId()) << " in [";
//                cout << lb[i] << "," << ub[i]<< "]\n";
//            }
//            cout << endl;
            break;
        }
        case integer_: {
            auto real_var = (var<int>*)v.get();
            auto lb = IloNumArray(*_cplex_env, real_var->get_dim());
            auto ub = IloNumArray(*_cplex_env, real_var->get_dim());
            for (int i = 0; i < real_var->get_dim(); i++) {
                lb[i] = real_var->get_lb(i);
                ub[i] = real_var->get_ub(i);
            }
            _cplex_vars.at(vid) = IloNumVarArray(*_cplex_env,lb,ub, ILOINT);
//            for (int i = 0; i < real_var->get_dim(); i++) {
//                cout << real_var->_indices->_keys->at(i) << " : ";
//                cout << to_string(_cplex_vars.at(vid)[i].getId()) << " in [";
//                cout << lb[i] << "," << ub[i]<< "]\n";
//            }
//            cout << endl;
            break;
        }
        case short_: {
            auto real_var = (var<short>*)v.get();
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
            _cplex_vars.at(vid) = IloNumVarArray(*_cplex_env,ILOBOOL);
            _cplex_vars.at(vid).setSize(v->get_dim());
            break;
        }
        default:
            break;
        }
    }
}

void CplexProgram::set_cplex_objective() {
//    if (!_model->_obj->_new) {
//        return;//Objective already added to the program
//    }
    _model->_obj->_new = false;
    size_t idx = 0, idx_inst = 0, idx1 = 0, idx2 = 0, idx_inst1 = 0, idx_inst2 = 0;
    IloNumExpr obj(*_cplex_env);
    for (auto& it_qterm: _model->_obj->get_qterms()) {
        IloNumExpr qterm(*_cplex_env);
        idx1 = it_qterm.second._p->first->get_vec_id();
        idx2 = it_qterm.second._p->second->get_vec_id();
        if (it_qterm.second._p->first->_is_vector) {//Vectorial/Matrix product
            if (it_qterm.second._coef_p1_tr) { // qterm = (coef*p1)^T*p2
                assert(it_qterm.second._p->first->_dim[1]==1 && it_qterm.second._coef->_dim[0]==it_qterm.second._p->second->_dim[0]);
                for (auto i = 0; i<it_qterm.second._p->first->_dim[0]; i++) {
                    for (auto j = 0; j<it_qterm.second._p->first->_dim[0]; j++) {
                        qterm += _model->_obj->eval(it_qterm.second._coef,i,j)*_cplex_vars[idx1][it_qterm.second._p->first->get_id_inst(i)]*_cplex_vars[idx2][it_qterm.second._p->second->get_id_inst(j)];
                    }
                }
            }
            else {//TODO fix this
                auto dim = it_qterm.second._p->first->get_dim();
                for (int j = 0; j<dim; j++) {
                    qterm += _model->_obj->eval(it_qterm.second._coef,j)*_cplex_vars[idx1][it_qterm.second._p->first->get_id_inst(j)]*_cplex_vars[idx2][it_qterm.second._p->second->get_id_inst(j)];
                }
            }
        }
        else {
            idx_inst1 = it_qterm.second._p->first->get_id_inst();
            idx_inst2 = it_qterm.second._p->second->get_id_inst();
            qterm += _model->_obj->eval(it_qterm.second._coef)*_cplex_vars[idx1][idx_inst1]*_cplex_vars[idx2][idx_inst2];
        }
        if (!it_qterm.second._sign) {
            qterm *= -1;
        }
        obj += qterm;
        qterm.end();
    }

    for (auto& it_lterm: _model->_obj->get_lterms()) {
        IloNumExpr lterm(*_cplex_env);
        idx = it_lterm.second._p->get_vec_id();
        if (it_lterm.second._p->_is_vector) {//Vectorial/Matrix product
            assert(it_lterm.second._p->_is_vector && it_lterm.second._coef->_dim[0]==1 && it_lterm.second._p->_dim[1]==1); // We're in the objective dimensions should reduce to one.
            auto dim = it_lterm.second._p->get_dim();
            for (int j = 0; j<dim; j++) {
                lterm += _model->_obj->eval(it_lterm.second._coef,j)*_cplex_vars[idx][it_lterm.second._p->get_id_inst(j)];
            }
        }
        else {
            idx_inst = it_lterm.second._p->get_id_inst();
            lterm += _model->_obj->eval(it_lterm.second._coef)*_cplex_vars[idx][idx_inst];
        }
        if (!it_lterm.second._sign) {
            lterm *= -1;
        }
        obj += lterm;
        lterm.end();
    }

    obj += _model->_obj->eval(_model->_obj->get_cst());

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
    Constraint<>* c;
    for(auto& p: _model->_cons) {
        c = p.second.get();
//        if (!c->_new) {
//            continue;//Constraint already added to the program
//        }
        c->_new = false;
        if (c->is_nonlinear()) {
            throw invalid_argument("Cplex cannot handle nonlinear constraints that are not convex quadratic.\n");
        }
        nb_inst = c->_dim[0];
        inst = 0;
        for (size_t i = 0; i< nb_inst; i++) {
            IloNumExpr cc(*_cplex_env);
            for (auto& it_qterm: c->get_qterms()) {
                IloNumExpr qterm(*_cplex_env);
                idx1 = it_qterm.second._p->first->get_vec_id();
                idx2 = it_qterm.second._p->second->get_vec_id();
                if (it_qterm.second._p->first->_is_vector) {
                    auto dim = it_qterm.second._p->first->get_dim(i);
                    for (size_t j = 0; j<dim; j++) {
                        qterm += c->eval(it_qterm.second._coef,i,j)*_cplex_vars[idx1][it_qterm.second._p->first->get_id_inst(i,j)]*_cplex_vars[idx2][it_qterm.second._p->second->get_id_inst(i,j)];
                    }
                }
                else {                    
                    idx_inst1 = it_qterm.second._p->first->get_id_inst(inst);
                    idx_inst2 = it_qterm.second._p->second->get_id_inst(inst);
                    qterm += c->eval(it_qterm.second._coef, inst)*_cplex_vars[idx1][idx_inst1]*_cplex_vars[idx2][idx_inst2];
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
                if (it_lterm.second._p->_is_vector || it_lterm.second._coef->is_matrix()) {
                    auto dim = it_lterm.second._p->get_dim(i);
                    for (int j = 0; j<dim; j++) {
                        lterm += c->eval(it_lterm.second._coef,i,j)*_cplex_vars[idx][it_lterm.second._p->get_id_inst(i,j)];
                    }                    
                }
                else {
                    idx_inst = it_lterm.second._p->get_id_inst(inst);
                    lterm += c->eval(it_lterm.second._coef, inst)*_cplex_vars[idx][idx_inst];
                }
                if (!it_lterm.second._sign) {
                    lterm *= -1;
                }
                cc += lterm;
                lterm.end();
            }
            cc += c->eval(c->get_cst(), inst);

            if(c->get_ctype()==geq) {
                _cplex_model->add(cc >= 0);
            }
            else if(c->get_ctype()==leq) {
                _cplex_model->add(cc <= 0);
            }
            else if(c->get_ctype()==eq) {
                _cplex_model->add(cc == 0);
            }
            inst++;
        }
    }    
}

void CplexProgram::prepare_model() {
    fill_in_cplex_vars();
    create_cplex_constraints();
    set_cplex_objective();
//    IloCplex cplex(*_cplex_model);
//    cplex.exportModel("lpex2.lp");

    //    print_constraints();
}
