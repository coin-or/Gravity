////
////  model.cpp
////  Gravity
////
////  Created by Hijazi, Hassan (Data61, Canberra City) on 6/5/17.
////
////
//
//#include <gravity/model.h>
//
//
//using namespace std;
//using namespace gravity;
//
///** Constructor */
////@{
//Model::Model(){
//    _nnz_g = 0;
//    _nnz_h = 0;
//};
////@}
//
///* Destructor */
//Model::~Model(){
//    for (auto &vp:_vars) {
//        delete vp.second;
//    }
////    for (auto &cp:_cons) {
////        delete cp.second;
////    }
//};
//
//
///* Accessors */
//size_t Model::get_nb_vars() const{
//    return _nb_vars;
//};
//
//
//
//
//size_t Model::get_nb_cons() const{
//    size_t n = 0;
//    for (auto &cp:_cons) {
//        n += cp.second->get_nb_instances();
//    }
//    return n;
//};
//
//size_t Model::get_nb_ineq() const{
//    size_t n = 0;
//    for (auto &cp:_cons) {
//        if (cp.second->is_ineq()) {
//            n += cp.second->get_nb_instances();
//        }
//    }
//    return n;
//};
//
//
//size_t Model::get_nb_nnz_g(){
//    _nnz_g = 0;
//    for (auto &cp:_cons) {
//        auto c = cp.second;
//        auto nb_inst = c->_dim[0];
//        for (size_t inst = 0; inst<nb_inst; inst++) {
//            if (!*c->_all_lazy || !c->_lazy[inst]) {
//                _nnz_g += c->get_nb_vars(inst);
//            }
//        }
//    }
//    return _nnz_g;
//};
//
//
///* Return the number of nonzeros in the lower left part of the hessian */
//size_t Model::get_nb_nnz_h(){
//    size_t idx = 0;
//    bool idx_inc = false;
//    Constraint* c = nullptr;
//    auto max_f_idx = 0;
//    for (auto &pairs: _hess_link) {
////        auto f_pair = *pairs.second.begin();
////        max_f_idx = 0;
////        auto f_idx = 0;
//        for (auto &f_pair:pairs.second) {
//            idx_inc = false;
//            auto f = f_pair.first;
//            if (f->_is_constraint) {
//                c = (Constraint*)f;
//            }
//            auto d2f = f_pair.second;
//            size_t nb_inst = f->_dim[0];
////            f_idx = 0;
//            for (size_t inst = 0; inst<nb_inst; inst++) {
//                if (!(f->_is_constraint && *c->_all_lazy && c->_lazy[inst])) {
//                    if (d2f->is_matrix()) {
//                        for (size_t i = 0; i < d2f->_dim[0]; i++) {
//                            for (size_t j = i; j < d2f->_dim[1]; j++) {
//                                idx++;
//                            }
//                        }
//                    }
//                    else if(d2f->_is_vector){
//                        for (size_t j = 0; j < d2f->_dim[0]; j++) {
//                            idx++;
//                        }
//                    }
//                    else {
//                        idx++;
//                    }
//                }
//            }
////            if(max_f_idx<f_idx) {
////                max_f_idx = f_idx;
////            }
//        }
////        idx += max_f_idx;
//    }
////    }
//    _nnz_h = idx;
//    return _nnz_h;
//};
//
//
////void Model::add_indices(const string& constr_name, const node_pairs& indices){
////    auto c = _cons_name.at(constr_name);
////    map<size_t, shared_ptr<Constraint>>  new_cons;
////    c->add_indices_in(indices);
////    for (auto &cc : _cons) {
////        if (cc.second->_id > c->_id) {//Ids of all constraints coming after c will be shifted by the size of indices
////            cc.second->_id += indices._keys.size();
////        }
////        new_cons[cc.second->_id] = cc.second;
////    }
////    _cons = new_cons;
////    _nb_cons = get_nb_cons();
////}
//
//shared_ptr<Constraint> Model::get_constraint(const string& cname) const{
//    return _cons_name.at(cname);
//}
//
//param_* Model::get_var_ptr(const string& vname) const{
//    auto it = _vars_name.find(vname);
//    if (it==_vars_name.end()) {
//        return nullptr;
//    }
//    return it->second;
//}
//
//param_* Model::get_var_ptr(size_t idx) const{
//    return _vars.at(idx);
//}
//
//
//template <typename type>
//var<type> Model::get_var(const string& vname) const{
//    auto it = _vars_name.find(vname);
//    if (it==_vars_name.end()) {
//        throw invalid_argument("In function: Model::get_var(const string& vname) const, unable to find variable with given name");
//    }
//    return *(var<type>*)it->second;
//}
//
//template var<double> Model::get_var(const string& vname) const;
//template var<int> Model::get_var(const string& vname) const;
//template var<short> Model::get_var(const string& vname) const;
//template var<bool> Model::get_var(const string& vname) const;
//template var<Cpx> Model::get_var(const string& vname) const;
//template var<float> Model::get_var(const string& vname) const;
//template var<long double> Model::get_var(const string& vname) const;
//
//
///* Modifiers */
//
//void Model::reindex(){
//    size_t cid = 0, new_cid = 0, nb_inst = 0;
//    shared_ptr<Constraint> c = nullptr;
//    map<size_t, shared_ptr<Constraint>>  new_cons;
//    for(auto& c_p: _cons_name)
//    {
//        c = c_p.second;
//        nb_inst = c->get_nb_instances();
//        if (nb_inst==0) {
//            continue;
//        }
//        if (_type==lin_m && c->is_quadratic()) {
//            _type = quad_m;
//        }
//        if ((_type==lin_m || _type==quad_m) && (c->is_nonlinear() || c->is_polynomial())) {
//            _type = nlin_m;
//        }
//        if (_is_convex && !c->is_convex() && !c->is_soc() && !c->is_rotated_soc()) {
//            _is_convex = false;
//        }
//        cid = c->_id;
//        if (cid!=new_cid) {
//            c->_id = new_cid;
//        }
//        new_cons[c->_id] = c;
//        new_cid = c->_id+nb_inst;
//    }
//    _cons = new_cons;
//    _nb_cons = get_nb_cons();
//}
//
//void Model::init_indices(){// Initialize the indices of all variables involved in the model
//    param_* v= nullptr;
//    size_t idx = 0;
//    for(auto& v_p: _vars)
//    {
//        v = v_p.second;
//        v->set_id(idx);
//        for (size_t i = 0; i < v->get_dim(); i++) {
//            idx++;
//        }
//    }
//}
//
//void Model::add_var(param_* v){
//    if (v->is_indexed()) {
//        throw invalid_argument("Should not add an indexed variable to the model");
//        return;
//    }
//
//    if (_vars_name.count(v->_name) == 0) {
//        v->set_id(_nb_vars);
//        v->set_vec_id(_vars.size());
//        _vars_name[v->_name] = v;
//        _vars[v->get_vec_id()] = v;
//        _nb_vars += v->get_dim();
//    }
//};
//
//
//
//
//void Model::del_var(const param_& v){
//    auto it = _vars.find(v.get_id());
//    if (it!=_vars.end()) {
//        _nb_vars -= v.get_dim();
//        delete it->second;
//        _vars.erase(it);
//    }
//};
//
//
//void Model::add_param(param_* v){
//    if (_params.count(v->get_id())==0) {
//        _nb_params += v->get_dim();
//        v->set_id(_params.size());
//        v->set_vec_id(_vars.size());
//        _params_name[v->get_name()] = v;
//        _params[v->get_vec_id()] = v;
//    }
//};
//
//void Model::add_param(param_& v){
//    if (_params.count(v.get_id())==0) {
//        _nb_params += v.get_dim();
//        auto newv = (param_*)copy(v);
//        v.set_id(_params.size());
//        v.set_vec_id(_vars.size());
//        newv->set_id(_params.size());
//        newv->set_vec_id(_vars.size());
//        _params_name[v.get_name()] = newv;
//        _params[v.get_vec_id()] = newv;
//    }
//};
//
//
//void Model::del_param(const param_& v){
//    auto it = _params.find(v.get_id());
//    if (it!=_params.end()) {
//        _nb_params -= v.get_dim();
//        delete it->second;
//        _params.erase(it);
//    }
//};
//
//
//void Model::add_integrality_cuts(){
//    for (auto &v_pair:_vars) {
//        auto v = v_pair.second;
//        switch (v->get_intype()) {
//            case float_: {
//                break;
//            }
//            case long_:{
//                break;
//            }
//            case double_:{
//                if (v->_is_relaxed) {
//                    auto real_var = (var<double>*)v;
//                    Constraint Integer_Cut("Integer_Cut"+v->_name);
//                    Integer_Cut += power((*real_var),2) - (*real_var);
//                    add(Integer_Cut==0);
//                }
//                break;
//            }
//            case integer_:{
//                auto real_var = (var<int>*)v;
//                Constraint Integer_Cut("Integer_Cut"+v->_name);
//                Integer_Cut += power((*real_var),2) - (*real_var);
//                add(Integer_Cut==0);
//                break;
//            }
//            case short_:{
//                auto real_var = (var<short>*)v;
//                Constraint Integer_Cut("Integer_Cut"+v->_name);
//                Integer_Cut += power((*real_var),2) - (*real_var);
//                add(Integer_Cut==0);
//                break;
//            }
//            case binary_:{
//                auto real_var = (var<bool>*)v;
//                Constraint Integer_Cut("Integer_Cut"+v->_name);
//                Integer_Cut += power((*real_var),2) - (*real_var);
//                add(Integer_Cut==0);
//                break;
//            }
//
//        }
//    }
//}
//
//void Model::add(const Constraint& c){
//    if (c.get_dim()==0) {
//        return;
//    }
//    add_constraint(c);
//}
//void Model::add_lazy(const Constraint& c){
//    if (c.get_dim()==0) {
//        return;
//    }
//    *c._all_lazy = true;
//    auto newc = add_constraint(c);
//    newc->make_lazy();
//    _has_lazy = true;
//}
//
//void Model::replace(param_* v, func_& f){
//    for (auto &c_p: _cons_name) {
//        auto c = c_p.second;
//        if (!c->has_var(*v)) {
//            continue;
//        }
//        c->replace(v, f);
//    }
//    _vars_name.erase(v->_name);
//    auto vid = *v->_vec_id;
//    delete _vars.at(vid);
//    _vars.erase(vid);
//
//}
//
//
//void Model::project() {
//    for (auto& c_pair:_cons_name) {
//        if (!c_pair.second->is_ineq()) {
//            auto &lterms = c_pair.second->get_lterms();
//            if (!lterms.empty()) {
//                auto first = lterms.begin();
//                auto v = first->second._p;
//                if (v->_is_vector) {
//                    continue;
//                }
//                auto f = *c_pair.second - c_pair.second->get_rhs();
//                if (first->second._sign) {
//                    f -= *v;
//                    f *= -1;
//                }
//                else {
//                    f += *v;
//                }
//                DebugOff(f.to_str());
//                _cons.erase(c_pair.second->_id);
//                _cons_name.erase(c_pair.first);
//                replace(v,f);
//                project();
//                return;
//            }
//        }
//    }
//}
//shared_ptr<Constraint> Model::add_constraint(const Constraint& c){
//    if (c.get_dim()==0) {
//        return nullptr;
//    }
//    if (_cons_name.count(c.get_name())==0) {
//        auto newc = make_shared<Constraint>(c);
//        if (newc->is_constant()) {
//            switch (newc->_ctype) {
//                case leq:
//                    if (newc->eval()>c.get_rhs()) {
//                        throw invalid_argument("Adding violated constant constraint!\n");
//                    }
//                    break;
//                case geq:
//                    if (newc->eval()<c.get_rhs()) {
//                        throw invalid_argument("Adding violated constant constraint!\n");
//                    }
//                    break;
//                case eq:
//                    if (newc->eval()!=c.get_rhs()) {
//                        throw invalid_argument("Adding violated constant constraint!\n");
//                    }
//                    break;
//                default:
//                    break;
//            }
//            Warning("WARNING: Adding redundant constant constraint, Gravity will be ignoring it.\n");
//            return newc;
//        }
//
//        newc->_val->resize(newc->get_dim());
//        newc->update_to_str();
////        newc->print();
////        if (newc->is_nonlinear()) {
//            for (auto &p_t: newc->get_lterms()) {
//                if (p_t.second._coef->is_function()) {
//                    auto f = (func_*)p_t.second._coef;
//                    auto exp = f->get_expr();
//                    if (exp) {
//                        embed(*exp);
//                    }
//                }
//            }
//            for (auto &p_t: newc->get_qterms()) {
//                if (p_t.second._coef->is_function()) {
//                    auto f = (func_*)p_t.second._coef;
//                    auto exp = f->get_expr();
//                    if (exp) {
//                        embed(*exp);
//                    }
//                }
//            }
//            for (auto &p_t: newc->get_pterms()) {
//                if (p_t.second._coef->is_function()) {
//                    auto f = (func_*)p_t.second._coef;
//                    auto exp = f->get_expr();
//                    if (exp) {
//                        embed(*exp);
//                    }
//                }
//            }
//            if (newc->get_cst()->is_function()) {
//                auto cc = (func_*) newc->get_cst();
//                auto exp = cc->get_expr();
//                if (exp) {
//                    embed(*exp);
//                }
//            }
//            auto exp = newc->get_expr();
//        if (exp) {
//            if (exp->is_uexpr()) {
//                auto ue = (uexpr*)exp.get();
//                auto f = ue->_son;
//                bool found_cpy = false;
//                auto name = f->to_str();
//                if (name.back()=='T') {
//                    name = name.substr(0,name.size()-2);
//                    if (_nl_funcs_map.count(name)>0) {
//                        auto cpy = _nl_funcs_map.at(name);
//                        f->_val = cpy->_val;
//                        f->_evaluated = true;
//                        found_cpy = true;
//                    }
//                }
//                else {
//                    auto name1 = "[["+name+"]]\u1D40";
//                    if (_nl_funcs_map.count(name1)>0) {
//                        auto cpy = _nl_funcs_map.at(name1);
//                        f->_val = cpy->_val;
//                        f->_evaluated = true;
//                        found_cpy = true;
//                    }
//                    auto name2 = name+"\u1D40";
//                    if (_nl_funcs_map.count(name2)>0) {
//                        auto cpy = _nl_funcs_map.at(name2);
//                        f->_val = cpy->_val;
//                        f->_evaluated = true;
//                        found_cpy = true;
//                    }
//                }
//                if (!found_cpy) {
//                    auto f_p = _nl_funcs_map.insert(make_pair<>(f->to_str(), f));
//                    if (f_p.second) {
//                        embed(f);
//                        _nl_funcs.push_back(f);
//                        DebugOff(f->to_str() << endl);
//                        //                f->_val = make_shared<vector<double>>();
//                        //                f->_val->resize(f->_dim[0]);
//                    }
//                    else {
//                        ue->_son = f_p.first->second;
//                    }
//                }
//            }
//            else {
//                auto be = (bexpr*)exp.get();
//                auto f = be->_lson;
//                bool found_cpy = false;
//                auto name = f->to_str();
//                if (name.back()=='T') {
//                    name = name.substr(0,name.size()-2);
//                    if (_nl_funcs_map.count(name)>0) {
//                        auto cpy = _nl_funcs_map.at(name);
//                        f->_val = cpy->_val;
//                        f->_evaluated = true;
//                        found_cpy = true;
//                    }
//                }
//                else {
//                    auto name1 = "[["+name+"]]\u1D40";
//                    if (_nl_funcs_map.count(name1)>0) {
//                        auto cpy = _nl_funcs_map.at(name1);
//                        f->_val = cpy->_val;
//                        f->_evaluated = true;
//                        found_cpy = true;
//                    }
//                    auto name2 = name+"\u1D40";
//                    if (_nl_funcs_map.count(name2)>0) {
//                        auto cpy = _nl_funcs_map.at(name2);
//                        f->_val = cpy->_val;
//                        f->_evaluated = true;
//                        found_cpy = true;
//                    }
//                }
//                if (!found_cpy) {
//                    auto f_p = _nl_funcs_map.insert(make_pair<>(f->to_str(), f));
//                    if (f_p.second) {
//                        embed(f);
//                        DebugOff(f->to_str() << endl);
//                        _nl_funcs.push_back(f);
//                        //                f->_val = make_shared<vector<double>>();
//                        //                f->_val->resize(f->_dim[0]);
//                    }
//                    else {
//                        be->_lson = f_p.first->second;
//                    }
//                }
//                f = be->_rson;
//                found_cpy = false;
//                name = f->to_str();
//                if (name.back()=='T') {
//                    name = name.substr(0,name.size()-2);
//                    if (_nl_funcs_map.count(name)>0) {
//                        auto cpy = _nl_funcs_map.at(name);
//                        f->_val = cpy->_val;
//                        f->_evaluated = true;
//                        found_cpy = true;
//                    }
//                }
//                else {
//                    auto name1 = "[["+name+"]]\u1D40";
//                    if (_nl_funcs_map.count(name1)>0) {
//                        auto cpy = _nl_funcs_map.at(name1);
//                        f->_val = cpy->_val;
//                        f->_evaluated = true;
//                        found_cpy = true;
//                    }
//                    auto name2 = name+"\u1D40";
//                    if (_nl_funcs_map.count(name2)>0) {
//                        auto cpy = _nl_funcs_map.at(name2);
//                        f->_val = cpy->_val;
//                        f->_evaluated = true;
//                        found_cpy = true;
//                    }
//                }
//                if (!found_cpy) {
//                    auto f_p = _nl_funcs_map.insert(make_pair<>(f->to_str(), f));
//                    if (f_p.second) {
//                        embed(f);
//                        DebugOff(f->to_str() << endl);
//                        _nl_funcs.push_back(f);
//                        //                f->_val = make_shared<vector<double>>();
//                        //                f->_val->resize(f->_dim[0]);
//                    }
//                    else {
//                        be->_rson = f_p.first->second;
//                    }
//                }
//            }
//        }
////            embed();
////        }
//        
//        newc->_violated.resize(newc->_dim[0],true);
//        _cons_name[c.get_name()] = newc;
//        if(*newc->_all_lazy){
//            newc->_lazy.resize(newc->_dim[0],true);
//            newc->allocate_mem();
//            return newc;
//        }
//        size_t nb_inst = c.get_nb_instances();
//        if (nb_inst>0) {
////            if (_cons.empty()) {
//                newc->_id = _nb_cons;
////            }
////            else {
////                auto last_added = *_cons.end();
////                newc->_id = last_added.first + last_added.second->_dim[0];
////            }
//    //        embed(newc);
//            _cons[newc->_id] = newc;
//            _built = false;
//    //        newc->print_expanded();
//            _nb_cons += nb_inst;
//            for (size_t inst = 0; inst<nb_inst; inst++) {
//                _nnz_g += c.get_nb_vars(inst);
//            }
//            if (_type==lin_m && c.is_quadratic()) {
//                _type = quad_m;
//            }
//            if ((_type==lin_m || _type==quad_m) && c.is_polynomial()) {
//                _type = pol_m;
//            }
//            if(c.is_nonlinear()){
//                _type = nlin_m;
//            }
//            if (_is_convex && !c.is_convex()) {
//                _is_convex = false;
//            }
//            newc->allocate_mem();
//        }
//        return newc;
//    }
//    else {
//        return _cons_name[c.get_name()];
////        throw invalid_argument("rename constraint as this name has been used by another one: " + c.to_str());
//    }
//};
//
//
//
//void Model::del_constraint(const Constraint& c){
//    _cons_name.erase(c.get_name());
//    _cons.erase(c._id);
//    _built = false;
//    /* TODO: make sure other infos in model are updated */
//};
//
//void Model::set_objective(const func_& f, ObjectiveType t) {
//    _obj = f;
//    _objt = t;
//    if (_is_convex) {
//        if ((t==maximize && f.is_convex()) || (t==minimize && !f.is_convex())) {
//            _is_convex = false;
//        }
//    }
////    embed(_obj);
//    /* TODO make sure the objective receivs the same treatment as the constraints. */
//}
//
//void Model::min(const func_& f){
//    _obj = f;
//    _obj._dim[0] = 1;
//    _obj.allocate_mem();
//    _objt = minimize;
//}
//
//void Model::max(const func_& f){
//    _obj = f;
//    _obj._dim[0] = 1;
//    _obj.allocate_mem();
//    _objt = maximize;
//}
//
//
//void Model::set_objective(pair<func_*, ObjectiveType> p){
//    _obj = *p.first;
//    _obj.allocate_mem();
//    _objt = p.second;
//}
//
//void Model::set_objective_type(ObjectiveType t) {
//    _objt = t;
//}
//
//bool Model::has_violated_constraints(double tol){
////    if (!_has_lazy) {
////        return false;
////    }
////    int cid = 0;
//    size_t nb_inst = 0, nb_viol = 0, nb_viol_all = 0;
//    size_t nb_active = 0, nb_active_all = 0;
//    double diff = 0;
//    shared_ptr<Constraint> c = nullptr;
//    bool violated = false;
//    for(auto& c_p: _cons_name)
//    {
//        c = c_p.second;
////        cid = c->_id;
//        nb_inst = c->_dim[0];
//        nb_viol = 0;
//        nb_active = 0;
//        c->_all_satisfied = true;
//        c->_violated.resize(nb_inst);
//        c->_active.resize(nb_inst);
//        switch (c->get_type()) {
//            case eq:
//                for (size_t inst=0; inst<nb_inst; inst++) {
//                    diff = abs(c->eval(inst) - c->_rhs);
//                    if(diff > tol) {
//                        DebugOff("Violated equation: ");
////                        c->print(inst);
//                        DebugOff(", violation = "<< diff << endl);
//                        nb_viol++;
////                        violated = true;
//                        if (*c->_all_lazy) {
//                            c->_all_satisfied = false;
//                            c->_violated[inst] = true;
//                            violated = true;
//                            c->_lazy[inst] = false;
//                        }
//                        else {
////                            throw runtime_error("Non-lazy constraint is violated, solution declared optimal by solver!\n" + c->to_str(inst));
//                        }
////                        c->_violated[inst] = true;
//                    }
//                    else {
////                        c->_violated[inst] = false;
//                    }
////                    nb_active++;
//                }
//                break;
//            case leq:
//                for (size_t inst=0; inst<nb_inst; inst++) {
//                    c->_violated[inst] = false;
//                    diff = c->eval(inst) - c->_rhs;
//                    if(diff > tol) {
//                        DebugOn("Violated inequality: ");
//                        c->print(inst);
//                        DebugOn(", violation = "<< diff << endl);
//                        nb_viol++;
////                        violated = true;
//                        if (*c->_all_lazy) {
//                            *c->_all_lazy = false;
//                            c->_all_satisfied = false;
//                            c->_violated[inst] = true;
//                            violated = true;
//                            c->_lazy[inst] = false;
//                        }
//                        else {
////                            throw runtime_error("Non-lazy constraint is violated, solution declared optimal by solver!\n" + c->to_str(inst));
//                        }
//                    }
//                    else if (abs(diff)>tol) {
//                        c->_active[inst] = false;
////                        if (*c->_all_lazy) {
////                            c->_lazy[inst] = true;
////                        }
//                    }
//                    else {
//                        nb_active++;
//                    }
//                }
//                break;
//            case geq:
//                for (size_t inst=0; inst<nb_inst; inst++) {
//                    c->_violated[inst] = false;
//                    diff = c->eval(inst) - c->_rhs;
//                    if(diff < -tol) {
//                        DebugOff("Violated inequality: ");
////                        c->print(inst);
//                        DebugOff(", violation = "<< diff << endl);
//                        nb_viol++;
////                        violated = true;
//                        if (*c->_all_lazy) {
//                            *c->_all_lazy = false;
//                            c->_all_satisfied = false;
//                            c->_violated[inst] = true;
//                            violated = true;
//                            c->_lazy[inst] = false;
//                        }
//                        else {
////                            throw runtime_error("Non-lazy constraint is violated, solution declared optimal by solver!\n" + c->to_str(inst));
//                        }
//                    }
//                    else if (abs(diff)> tol) {
//                        c->_active[inst] = false;
////                        if (*c->_all_lazy) {
////                            c->_lazy[inst] = true;
////                        }
//                    }
//                    else {
//                        nb_active++;
//                    }
//                }
//                break;
//                
//            default:
//                break;
//        }
////        *c->_all_lazy = false;
//        nb_viol_all += nb_viol;
//        nb_active_all += nb_active;
//        if (nb_viol>0) {
//            DebugOn("Percentage of violated constraints for " << c->get_name() << " = " << to_string(100.*nb_viol/nb_inst) << "%\n");
//        }
//        if (c->get_type()!=eq) {
//            DebugOff("Percentage of active constraints for " << c->get_name() << " = " << to_string(100.*nb_active/nb_inst) << "%\n");
//        }
//    }
//    DebugOn("Total percentage of violated constraints = " << to_string(100.*nb_viol_all/_nb_cons) << "%\n");
//    auto nb_ineq = get_nb_ineq();
//    DebugOn("Total percentage of active constraints = " << to_string(100.*nb_active_all/nb_ineq) << "%\n");
//    return violated;
//}
//
//bool Model::is_feasible(double tol){
//    size_t vid;
//    param_* v;
//    double viol = 0;
//    bool feasible = true;
//    for(auto& v_p: _vars)
//    {
//        v = v_p.second;
//        vid = v->get_id();
//        switch (v->get_intype()) {
//            case float_: {
//                auto real_var = (var<float>*)v;
//                for (size_t i = 0; i < real_var->get_dim(); i++) {
//                    viol = real_var->get_lb(i) - real_var->_val->at(i);
//                    if(viol  > tol){
//                        clog << "Violated lower bound for var: " << v->_name << ", index: "<< i << ", violation = " << viol << endl;
//                        feasible = false;
//                        
//                    }
//                    viol = real_var->_val->at(i) - real_var->get_ub(i);
//                    if(viol > tol){
//                        clog << "Violated upper bound for var: " << v->_name << ", index: "<< i << ", violation = " << viol << endl;
//                        feasible = false;
//                    }
//                }
//                break;
//            }
//            case long_:{
//                auto real_var = (var<long double>*)v;
//                for (size_t i = 0; i < real_var->get_dim(); i++) {
//                    viol = real_var->get_lb(i) - real_var->_val->at(i);
//                    if(viol  > tol){
//                        clog << "Violated lower bound for var: " << v->_name << ", index: "<< i << ", violation = " << viol << endl;
//                        feasible = false;
//                    }
//                    viol = real_var->_val->at(i) - real_var->get_ub(i);
//                    if(viol > tol){
//                        clog << "Violated upper bound for var: " << v->_name << ", index: "<< i << ", violation = " << viol << endl;
//                        feasible = false;
//                    }
//                }
//                break;
//            }
//            case double_:{
//                auto real_var = (var<double>*)v;
//                for (size_t i = 0; i < real_var->get_dim(); i++) {
//                    viol = real_var->get_lb(i) - real_var->_val->at(i);
//                    if(viol  > tol){
//                        clog << "Violated lower bound for var: " << v->_name << ", index: "<< i << ", violation = " << viol << endl;
//                        feasible = false;
//                    }
//                    viol = real_var->_val->at(i) - real_var->get_ub(i);
//                    if(viol > tol){
//                        clog << "Violated upper bound for var: " << v->_name << ", index: "<< i << ", violation = " << viol << endl;
//                        feasible = false;
//                    }
//                }
//                break;
//            }
//            case integer_:{
//                auto real_var = (var<int>*)v;
//                for (size_t i = 0; i < real_var->get_dim(); i++) {
//                    viol = real_var->get_lb(i) - real_var->_val->at(i);
//                    if(viol  > tol){
//                        clog << "Violated lower bound for var: " << v->_name << ", index: "<< i << ", violation = " << viol << endl;
//                        feasible = false;
//                    }
//                    viol = real_var->_val->at(i) - real_var->get_ub(i);
//                    if(viol > tol){
//                        clog << "Violated upper bound for var: " << v->_name << ", index: "<< i << ", violation = " << viol << endl;
//                        feasible = false;
//                    }
//                }
//                break;
//            }
//            case short_:{
//                auto real_var = (var<short>*)v;
//                for (size_t i = 0; i < real_var->get_dim(); i++) {
//                    viol = real_var->get_lb(i) - real_var->_val->at(i);
//                    if(viol  > tol){
//                        clog << "Violated lower bound for var: " << v->_name << ", index: "<< i << ", violation = " << viol << endl;
//                        feasible = false;
//                    }
//                    viol = real_var->_val->at(i) - real_var->get_ub(i);
//                    if(viol > tol){
//                        clog << "Violated upper bound for var: " << v->_name << ", index: "<< i << ", violation = " << viol << endl;
//                        feasible = false;
//                    }
//                }
//                break;
//            }
//            case binary_:{
//                break;
//            }
//            default:
//                break;
//        } ;
//    }
//    
//    feasible = feasible && !has_violated_constraints(tol);
//    return feasible;
//}
//
//
//#ifdef USE_BONMIN
//void Model::fill_in_var_types(Bonmin::TMINLP::VariableType* var_types){
//    size_t vid;
//    param_* v;
//    for(auto& v_p: _vars)
//    {
//        v = v_p.second;
//        vid = v->get_id();
//        switch (v->get_intype()) {
//            case float_: {
//                auto real_var = (var<float>*)v;
//                for (size_t i = 0; i < real_var->get_dim(); i++) {
//                   var_types[vid+i] = Bonmin::TMINLP::CONTINUOUS;
//                }
//                break;
//            }
//            case long_:{
//                auto real_var = (var<long double>*)v;
//                for (size_t i = 0; i < real_var->get_dim(); i++) {
//                    var_types[vid+i] = Bonmin::TMINLP::CONTINUOUS;
//                }
//                break;
//            }
//            case double_:{
//                auto real_var = (var<double>*)v;
//                DebugOff(real_var->get_name() << " in:" << endl);
//                if (real_var->_is_relaxed) {
//                    for (size_t i = 0; i < real_var->get_dim(); i++) {
//                        var_types[vid+i] = Bonmin::TMINLP::BINARY;
//                    }
//                }
//                else {
//                    for (size_t i = 0; i < real_var->get_dim(); i++) {
//                        var_types[vid+i] = Bonmin::TMINLP::CONTINUOUS;
//                    }
//                }
//                DebugOff(";" << endl);
//                break;
//            }
//            case integer_:{
//                auto real_var = (var<int>*)v;
//                for (size_t i = 0; i < real_var->get_dim(); i++) {
//                    var_types[vid+i] = Bonmin::TMINLP::INTEGER;
//                }
//                break;
//            }
//            case short_:{
//                auto real_var = (var<short>*)v;
//                for (size_t i = 0; i < real_var->get_dim(); i++) {
//                    var_types[vid+i] = Bonmin::TMINLP::INTEGER;
//                }
//                break;
//            }
//            case binary_:{
//                auto real_var = (var<bool>*)v;
//                for (size_t i = 0; i < real_var->get_dim(); i++) {
//                    var_types[vid+i] = Bonmin::TMINLP::BINARY;
//                }
//                break;
//            }
//            default:
//                break;
//        } ;
//    }
//}
//#endif
//
//
//void Model::fill_in_var_bounds(double* x_l ,double* x_u) {
//    size_t vid;
//    param_* v;
//    for(auto& v_p: _vars)
//    {
//        v = v_p.second;
//        vid = v->get_id();
//        switch (v->get_intype()) {
//            case float_: {
//                auto real_var = (var<float>*)v;
//                for (size_t i = 0; i < real_var->get_dim(); i++) {
//                    x_l[vid+i] = (double)real_var->get_lb(i);
//                    x_u[vid+i] = (double)real_var->get_ub(i);
//                }
//                break;
//            }
//            case long_:{
//                auto real_var = (var<long double>*)v;
//                for (size_t i = 0; i < real_var->get_dim(); i++) {
//                    x_l[vid+i] = (double)real_var->get_lb(i);
//                    x_u[vid+i] = (double)real_var->get_ub(i);
//                }
//                break;
//            }
//            case double_:{
//                auto real_var = (var<double>*)v;
//                DebugOff(real_var->get_name() << " in:" << endl);
//                for (size_t i = 0; i < real_var->get_dim(); i++) {
//                    x_l[vid+i] = real_var->get_lb(i);
//                    x_u[vid+i] = real_var->get_ub(i);
//                    DebugOff("(" << i << ")" << " : [" << x_l[vid+i] << "," << x_u[vid+i] << "]\n");
//                }
//                DebugOff(";" << endl);
//                break;
//            }
//            case integer_:{
//                auto real_var = (var<int>*)v;
//                for (size_t i = 0; i < real_var->get_dim(); i++) {
//                    x_l[vid+i] = (double)real_var->get_lb(i);
//                    x_u[vid+i] = (double)real_var->get_ub(i);
//                }
//                break;
//            }
//            case short_:{
//                auto real_var = (var<short>*)v;
//                for (size_t i = 0; i < real_var->get_dim(); i++) {
//                    x_l[vid+i] = (double)real_var->get_lb(i);
//                    x_u[vid+i] = (double)real_var->get_ub(i);
//                }
//                break;
//            }
//            case binary_:{
//                auto real_var = (var<bool>*)v;
//                for (size_t i = 0; i < real_var->get_dim(); i++) {
//                    x_l[vid+i] = (double)real_var->get_lb(i);
//                    x_u[vid+i] = (double)real_var->get_ub(i);
//                }
//                break;
//            }
//            default:
//                break;
//        } ;
//    }
//    //    cout << "idx = " << idx << endl;
//}
//
//void Model::set_x(const double* x){
//    size_t vid;
//    param_* v;
//    for(auto& v_p: _vars)
//    {
//        v = v_p.second;
//        vid = v->get_id();
//        switch (v->get_intype()) {
//            case float_: {
//                auto real_var = (var<float>*)v;
//                for (size_t i = 0; i < real_var->get_dim(); i++) {
//                    real_var->set_val(i, x[vid+i]);
//                }
//                break;
//            }
//            case long_:{
//                auto real_var = (var<long double>*)v;
//                for (size_t i = 0; i < real_var->get_dim(); i++) {
//                    real_var->set_val(i, x[vid+i]);
//                }
//                break;
//            }
//            case double_:{
//                auto real_var = (var<double>*)v;
//                for (size_t i = 0; i < real_var->get_dim(); i++) {
//                    real_var->set_val(i, x[vid+i]);
//                }
//                break;
//            }
//            case integer_:{
//                auto real_var = (var<int>*)v;
//                for (size_t i = 0; i < real_var->get_dim(); i++) {
//                    real_var->set_val(i, x[vid+i]);
//                }
//                break;
//            }
//            case short_:{
//                auto real_var = (var<short>*)v;
//                for (size_t i = 0; i < real_var->get_dim(); i++) {
//                    real_var->set_val(i, x[vid+i]);
//                }
//                break;
//            }
//            case binary_:{
//                auto real_var = (var<bool>*)v;
//                for (size_t i = 0; i < real_var->get_dim(); i++) {
//                    real_var->set_val(i, x[vid+i]);
//                }
//                break;
//            }
//            default:
//                break;
//        } ;
//    }
//}
//
//void Model::compute_funcs() {
////    if (_type!=nlin_m) {
////        return;
////    }
//    auto it = _nl_funcs.begin();
//    while (it!=_nl_funcs.end()) {
//        auto f = (*it++);
//        DebugOff(f->to_str() << endl);
//        if (f->is_constant() && f->_evaluated) {
//            continue;
//        }
//        if (!f->is_constant() && (_type!=nlin_m)) {//no need to precompute these for polyomial programs
//            continue;
//        }
//        if (!f->is_matrix() && !f->_is_vector) {
//            DebugOff(f->to_str()<<endl);
//            for (size_t inst = 0; inst < f->_dim[0]; inst++) {
//                f->eval(inst);
//            }
//        }
//        else if(!f->is_matrix()){ //vector
//            f->eval_vector();
//        }
//        else {
//            DebugOff(f->to_str()<<endl);
//            f->eval_matrix();
//        }
//        if (f->is_constant()) {
//            f->_evaluated = true;
//        }
//    }
//////        for (int inst = 0; inst < c->_dim[0]; inst++) {
//////            c->eval(inst);
//////        }
////        if (c->is_nonlinear() || (c->is_linear() && !_first_call_jac)) {
//////        if ((c->is_linear() && !_first_call_jac)) {
////            continue;
////        }
////        for (auto &dfdx: *c->get_dfdx()) {
////            if (dfdx.second->is_matrix()) {
////                for (int i = 0; i<dfdx.second->_dim[0]; i++) {
////                    for (size_t j = 0; j<dfdx.second->_dim[1]; j++) {
////                        dfdx.second->eval(i,j);
////                    }
////                }
////            }
//////            else if (dfdx.second->_is_vector) {
//////                for (int i = 0; i<dfdx.second->_dim[0]; i++) {
//////                    dfdx.second->_val->at(i) = dfdx.second->eval(i);
//////                }
//////            }
////            else {
////                for (int inst = 0; inst < dfdx.second->_dim[0]; inst++) {
////                    dfdx.second->eval(inst);
////                }
////            }
////            if (c->is_quadratic() && !_first_call_hess) {
////                continue;
////            }
////            for (auto &dfd2x: *dfdx.second->get_dfdx()) {
////                if (dfd2x.second->is_matrix()) {
////                    for (int i = 0; i<dfd2x.second->_dim[0]; i++) {
////                        for (size_t j = 0; j<dfd2x.second->_dim[1]; j++) {
////                            dfd2x.second->eval(i,j);
////                        }
////                    }
////                }
////                else {
////                    for (int inst = 0; inst < dfd2x.second->_dim[0]; inst++) {
////                        dfd2x.second->eval(inst);
////                    }
////                }
////            }
////        }
////    }
//}
//
//void Model::fill_in_obj(const double* x , double& res, bool new_x){
////    if (!new_x) {// IPOPT SEEMS TO BE INCONSISTENT ON NEW_X HERE!!
//    if (new_x) {
//        set_x(x);
//        compute_funcs();
//    }
////        set_x(x);
//        res = _obj.eval();
//        _obj_val = res;
//    _obj._new = false;
////    }
////    else {
////        res = _obj_val;
////    }
//    DebugOff("Objective = " << to_string(res) << endl);
//}
//
//void Model::fill_in_grad_obj(const double* x , double* res, bool new_x){
//    param_* v;
//    shared_ptr<func_> df;
//    size_t vid, vid_inst, index = 0;
//    unique_id v_unique;
//    if (new_x) {
//        set_x(x);
//        compute_funcs();
//    }
//    for (size_t i = 0; i<_nb_vars; i++) {
//        res[i] = 0;
//    }
//    if (_first_call_gard_obj) {
//        _obj_grad_vals.resize(_obj.get_nb_vars());
//        _first_call_gard_obj = false;
//    }
//    else if (_obj.is_linear()) {
////    else if (false) { /* No need to recompute jacobian for linear objectives */
//        for(auto& vi_p: _obj.get_vars())
//        {
//            v = vi_p.second.first.get();
//            vid = v->get_id();
//            if (v->_is_vector) {
//                for (size_t i = 0; i < v->get_dim(); i++) {
//                    vid_inst = vid + v->get_id_inst(i);
//                    res[vid_inst] = _obj_grad_vals[index++];
//                }
//            }
//            else {
//                vid_inst = vid + v->get_id_inst();
//                res[vid_inst] = _obj_grad_vals[index++];
//            }
//        }
//        return;
//    }
//    for(auto& vi_p: _obj.get_vars()) {
//        v = vi_p.second.first.get();
//        vid = v->get_id();
//        df = _obj.get_stored_derivative(v->_name);
////        if (v->is_matrix()) {
////            for (size_t i = 0; i < v->_dim[0]; i++) {
////                for (size_t j = 0; j < v->_dim[1]; j++) {
////                    vid_inst = vid + v->get_id_inst(i);
////                    res[vid_inst] = df->eval(i);
////                    _obj_grad_vals[index++] =res[vid_inst];
////                }
////            }
////        }
//        if (v->_is_vector) {
//            for (size_t i = 0; i < v->get_dim(); i++) {
//                vid_inst = vid + v->get_id_inst(i);
//                res[vid_inst] = df->eval(i);
//                _obj_grad_vals[index++] =res[vid_inst];
//            }
//        }
//        else {
//            vid_inst = vid + v->get_id_inst();
//            res[vid_inst] = df->eval();
//            _obj_grad_vals[index++] =res[vid_inst];
//        }
//    }
//}
//
//
//void Model::add_McCormick(std::string name, param_* vlift, param_* p1, param_* p2) {
//    var<> v = *(var<double>*)vlift;
//    var<> v1 = *(var<double>*)p1;
//    var<> v2 = *(var<double>*)p2;
//    Constraint MC1(name+"_McCormick1");
//    auto lb1 = v1.get_lb(v1.get_id_inst());
//    auto lb2 = v2.get_lb(v2.get_id_inst());
//    auto ub1 = v1.get_ub(v1.get_id_inst());
//    auto ub2 = v2.get_ub(v2.get_id_inst());
//    bool template_cstr = p1->_dim[0]>1;
//    MC1 += v;
//    if(template_cstr){//Template constraint
//        MC1 -= (*v1._lb)*v2 + (*v2._lb)*v1 - (*v1._lb)*(*v2._lb);
//    }
//    else {
//        MC1 -= lb1*v2 + lb2*v1 - lb1*lb2;
//    }
//    MC1 >= 0;
//    add(MC1);
//    //    MC1.print();
//    Constraint MC2(name+"_McCormick2");
//    MC2 += v;
//    if(template_cstr){//Template constraint
//        MC2 -= (*v1._ub)*v2 + (*v2._ub)*v1 - (*v1._ub)*(*v2._ub);
//    }
//    else {
//        MC2 -= ub1*v2 + ub2*v1 - ub1*ub2;
//    }
//    MC2 >= 0;
//    add(MC2);
////    //    MC2.print();
//    Constraint MC3(name+"_McCormick3");
//    MC3 += v;
//    if(template_cstr){//Template constraint
//        MC3 -= (*v1._lb)*v2 + (*v2._ub)*v1 - (*v1._lb)*(*v2._ub);
//    }
//    else {
//        MC3 -= lb1*v2 + ub2*v1 - lb1*ub2;
//    }
//    MC3 <= 0;
//    add(MC3);
////    //    MC3.print();
//    Constraint MC4(name+"_McCormick4");
//    MC4 += v;
//    if(template_cstr){//Template constraint
//        MC4 -= (*v1._ub)*v2 + (*v2._lb)*v1 - (*v1._ub)*(*v2._lb);
//    }
//    else{
//        MC4 -= ub1*v2 + lb2*v1 - ub1*lb2;
//    }
//    MC4 <= 0;
//    add(MC4);
//    //    MC4.print();
//}
//
//
///** Build the sequential McCormick relaxation for polynomial programs **/
//shared_ptr<Model> Model::build_McCormick(){
//    replace_integers();
//    if (_type==nlin_m) {
//        cerr << "Can only build a McCormick relaxation for polynomial programs, returning null" << endl;
//        return nullptr;
//    }
//    if (_type==lin_m) {
//        cerr << "No need to build a McCormick relaxation for a linear program, returning null" << endl;
//        return nullptr;
//    }
//    shared_ptr<Model> Mc = make_shared<Model>("McCormick Relaxation");
//    shared_ptr<Constraint> cstr;
//    param_* v;
//    for (auto &var_p:_vars) {
//        v = var_p.second;
//        switch (v->get_intype()) {// TODO check for other types
//            case double_:
//                Mc->add(*(var<double>*)v);
//                break;
//            case integer_:
//                Mc->add(*(var<int>*)v);
//                break;
//            case binary_:
//                Mc->add(*(var<bool>*)v);
//                break;
//            default:
//                break;
//        }
//    }
//    if ((_obj.is_convex() && _objt==minimize) || (_obj.is_concave() && _objt==maximize)) {
//        Mc->_obj = _obj;
//    }
//    else{
//        func_ new_obj = *_obj.get_cst();
//        Mc->_objt = _objt;
//        for(auto &lt : _obj.get_lterms()){
//            new_obj.insert(lt.second);
//        }
//        for(auto &qt : _obj.get_qterms()){
//            
//            if ((_objt==maximize && _obj.get_convexity(qt.second)==concave_) || (_objt==minimize && _obj.get_convexity(qt.second)==convex_)) {
//                new_obj.insert(qt.second);
//            }
//            else {//Lift products
//                string new_name;
//                if(qt.second._p->first==qt.second._p->second){
//                    new_name  = qt.second._p->first->_name+"²_lifted";
//                }
//                else {
//                    new_name = qt.second._p->first->_name+"_"+qt.second._p->second->_name+"_lifted";
//                }
//                var<> v(new_name);
//                auto new_v = Mc->get_var_ptr(new_name);
//                if(new_v == nullptr){
//                    Mc->add(v);
//                    new_v = (param_*)&v;
//                }
//                Mc->add_McCormick(new_v->_name, new_v, qt.second._p->first, qt.second._p->second);
//                if(qt.second._sign) {
//                    new_obj += *new_v;
//                }
//                else {
//                    new_obj -= *new_v;
//                }
//            }
//        }
//    }
//    for (auto &cst_p:_cons) {
//        cstr = cst_p.second;
//        if (cstr->is_convex()) {
//            Mc->add(*cstr);
//        }
//        else if(cstr->func_::is_convex()){
//            Constraint new_cstr(*cstr);
//            new_cstr._ctype = leq;
//            Mc->add(new_cstr);
//        }
//        else if(cstr->func_::is_concave()){
//            Constraint new_cstr(*cstr);
//            new_cstr._ctype = geq;
//            Mc->add(new_cstr);
//        }
//        else{
//            Constraint new_cstr(cstr->get_name()+"_lifted");
//            new_cstr += *cstr->get_cst();
//            new_cstr._ctype = cstr->_ctype;
//            new_cstr._rhs = cstr->get_rhs();
//            for(auto &lt : cstr->get_lterms()){
//                new_cstr.insert(lt.second);
//            }
//            for(auto &qt : cstr->get_qterms()){
//                
//                if ((cstr->_ctype==geq && cstr->get_convexity(qt.second)==concave_) || (cstr->_ctype==leq && cstr->get_convexity(qt.second)==convex_)) {
//                    new_cstr.insert(qt.second);
//                }
//                else {//Lift products
//                    string new_name;
//                    if(qt.second._p->first==qt.second._p->second){
//                        new_name  = qt.second._p->first->_name+"²_lifted";
//                    }
//                    else {
//                        new_name = qt.second._p->first->_name+"_"+qt.second._p->second->_name+"_lifted";
//                    }
//                    var<> v(new_name);
//                    auto new_v = Mc->get_var_ptr(new_name);
//                    if(new_v == nullptr){
//                        Mc->add(v);
//                        new_v = (param_*)&v;
//                    }
//                    Mc->add_McCormick(new_v->_name, new_v, qt.second._p->first, qt.second._p->second);
//                    if(qt.second._sign) {
//                        new_cstr += *new_v;
//                    }
//                    else {
//                        new_cstr -= *new_v;
//                    }
//                }
//            }
//            //            for(auto &qt : cstr->get_qterms()){
//            //TODO polynomial part
//            Mc->add(new_cstr);
//        }
//    }
//    return Mc;
//    
//}
//
//void compute_constrs(vector<shared_ptr<Constraint>>& v, double* res, int i, int j){
////    DebugOff("Calling compute_constrts with i =  " << i << "and j = "<< j << endl);
//    for (size_t idx = i; idx < j; idx++) {
//        auto c = v[idx];
//        size_t nb_ins = c->_dim[0];
//        size_t id = 0;
//        for (size_t inst = 0; inst< nb_ins; inst++){
//            if (!*c->_all_lazy || !c->_lazy[inst]) {
//                res[c->_id+id++] = c->eval(inst);
//                DebugOff("Accessing res at position " << c->_id+inst << endl);
//            //                _cons_vals[index++] = res[c->_id+inst];
//                DebugOff("g[" << to_string(c->_id+inst) << "] = " << to_string(res[c->_id+inst]) << endl);
//            }
//        }
//    }
//}
//
//void Model::fill_in_cstr(const double* x , double* res, bool new_x){
//    Constraint* c = nullptr;
////    size_t index = 0;
//    if (new_x) {
//        set_x(x);
//        compute_funcs();
//    }
////    if (_type==nlin_m) {
////        for(auto& c_p: _cons)
////        {
////            c = c_p.second.get();
////            auto nb_ins = c->_dim[0];
////            for (int inst = 0; inst< nb_ins; inst++){
////                res[c->_id+inst] = c->get_val(inst);
//////                _cons_vals[index++] = res[c->_id+inst];
////                DebugOff("g[" << to_string(c->_id+inst) << "] = " << to_string(res[c->_id+inst]) << endl);
////            }
////        }
////    }
////    else {
////    vector<Constraint*> cons;
//        _cons_vec.clear();
//        for(auto& c_p: _cons)
//        {
//            c_p.second->_new = false;
//            _cons_vec.push_back(c_p.second);
////            auto nb_ins = c->_dim[0];
////            for (int inst = 0; inst< nb_ins; inst++){
//////                res[c->_id+inst] = c->get_val(inst);
////                res[c->_id+inst] = c->eval(inst);
////                //                _cons_vals[index++] = res[c->_id+inst];
////                DebugOff("g[" << to_string(c->_id+inst) << "] = " << to_string(res[c->_id+inst]) << endl);
////            }
//        }
////    compute_constrs(_cons_vec, res, 0, _cons_vec.size());
//    unsigned nr_threads = std::thread::hardware_concurrency()/2;
//    if (nr_threads==0) {
//        nr_threads = 1;
//    }
//    vector<thread> threads;
//    /* Split cons into nr_threads parts */
//    vector<int> limits = bounds(nr_threads, _cons_vec.size());
//    DebugOff("limits size = " << limits.size() << endl);
//    for (size_t i = 0; i < limits.size(); ++i) {
//        DebugOff("limits[" << i << "] = " << limits[i] << endl);
//    }
//    /* Launch all threads in parallel */
//    for (unsigned i = 0; i < nr_threads; ++i) {
//        DebugOff("i = " << i << endl);
//        DebugOff("limits[" << i << "] = " << limits[i] << endl);
//        DebugOff("limits[" << i+1 << "] = " << limits[i+1] << endl);
//        threads.push_back(thread(compute_constrs, ref(_cons_vec), res, limits[i], limits[i+1]));
////        threads.push_back(thread(dummy));
//    }
//    /* Join the threads with the main thread */
//    for(auto &t : threads){
//        t.join();
//    }
//    
////    }
////    }
////    else {
////
////        for(auto& c_p: _cons) {
////            c = c_p.second;
////            auto nb_ins = c->_dim[0];
////            for (int inst = 0; inst< nb_ins; inst++){
////                res[c->_id+inst] = _cons_vals[index++];
////            }
////        }
////    }
//}
//
//
//
//void Model::fill_in_jac_nnz(int* iRow , int* jCol){
//    size_t idx=0, id = 0;
//    size_t cid = 0;
//    size_t vid = 0;
//    Constraint* c = NULL;
//    param_* v = NULL;
//    /* return the structure of the jacobian */
//    for(auto& c_p :_cons)
//    {
//        c = c_p.second.get();
//        c->_jac_cstr_idx = idx;
//        auto nb_ins = c->_dim[0];
//        for (auto &v_p: c->get_vars()){
//            v = v_p.second.first.get();
//            vid = v->get_id();
//            id = 0;
//            for (size_t inst = 0; inst< nb_ins; inst++){
//                if (!*c->_all_lazy || !c->_lazy[inst]) {
//                    cid = c->_id+id++;
//                    if (v->_is_vector) {
//                        auto dim = v->get_dim(inst);
//                        for (size_t j = 0; j<dim; j++) {
//                            iRow[idx] = cid;
//                            jCol[idx] = vid + v->get_id_inst(inst, j);
//                            idx++;
//                        }
//                    }
//                    else {
//                        iRow[idx] = cid;
//                        jCol[idx] = vid + v->get_id_inst(inst);
//                        idx++;
//                    }
//                }
//            }
//        }
//    }
//    if (idx!=_nnz_g) {
//        throw invalid_argument("idx!=_nnz_g");
//    }
//}
//
//
//void compute_jac(vector<shared_ptr<Constraint>>& vec, double* res, int i, int j, bool first_call, vector<double>& jac_vals){
//    size_t cid = 0, id_inst = 0;
//    string vid;
//    shared_ptr<Constraint> c = NULL;
//    param_* v = NULL;
//    shared_ptr<func_> dfdx;
//    auto idx = vec[i]->_jac_cstr_idx;
//    for (size_t id = i; id < j; id++) {
//        c = vec[id];
//        auto nb_ins = c->_dim[0];
//        id_inst = 0;
//        if (c->is_linear() && !first_call) {
//            //        if (false) {
//            DebugOff("Linear constraint, using stored jacobian!\n");
//            for (size_t i = 0; i<nb_ins; i++) {
//                if (!*c->_all_lazy || !c->_lazy[i]) {
//                    for (size_t j = 0; j<c->get_nb_vars(i); j++) {
//                        res[idx] = jac_vals[idx];
//                        idx++;
//                    }
//                }
//            }
//        }
//        else {
//            for (auto &v_p: c->get_vars()){
//                v = v_p.second.first.get();
//                vid = v->_name;
//                dfdx = c->get_stored_derivative(vid);
//                id_inst = 0;
//                for (size_t inst = 0; inst< nb_ins; inst++){
//                    if (!*c->_all_lazy || !c->_lazy[inst]) {
//                        cid = c->_id+id_inst++;
//                        if (v->_is_vector) {
//                            auto dim = v->get_dim(inst);
//                            for (size_t j = 0; j<dim; j++) {
//                                res[idx] = dfdx->eval(inst,j);//TODO: res[idx] += .. (account for vectors with repeated identical entries)
//                                jac_vals[idx] = res[idx];
//                                DebugOff("jac_val["<< idx <<"] = " << jac_vals[idx] << endl);
//                                idx++;
//                            }
//                        }
//                        else {
//                            res[idx] = dfdx->eval(inst);
//                            jac_vals[idx] = res[idx];
//                            idx++;
//                        }
//                    }
//                }
//            }
//        }
//    }
//}
//
///* Fill the nonzero values in the jacobian */
//void Model::fill_in_jac(const double* x , double* res, bool new_x){
////    if (!_first_call_jac && (!new_x || _type==lin_m)) { /* No need to recompute jacobian for linear models */
//    if (new_x) {
//        set_x(x);
//        compute_funcs();
//    }
//
//    if (!_first_call_jac && (_type==lin_m)) { /* No need to recompute jacobian for linear models */
//        for (size_t i = 0; i< _nnz_g; i++) {
//            res[i] = _jac_vals[i];
//        }
//        return;
//    }
//    size_t idx=0;
//    size_t cid = 0;
//    string vid;
//    Constraint* c = NULL;
//    param_* v = NULL;
//    shared_ptr<func_> dfdx;
////    vector<Constraint*> cons;
//    if (_type!=nlin_m) {//Polynomial, Quadratic or Linear
//        _cons_vec.clear();
//        for(auto& c_p: _cons)
//        {
//            c_p.second->_new = false;
//            _cons_vec.push_back(c_p.second);
//        }
//        
////        compute_jac(_cons_vec, res, 0, _cons_vec.size(), _first_call_jac, _jac_vals);
//        unsigned nr_threads = std::thread::hardware_concurrency()/2;
//        if (nr_threads==0) {
//            nr_threads = 1;
//        }
//        vector<thread> threads;
//        /* Split cons into nr_threads parts */
//        vector<int> limits = bounds(nr_threads, _cons_vec.size());
//        
//        /* Launch all threads in parallel */
//        for (unsigned i = 0; i < nr_threads; ++i) {
//            threads.push_back(thread(compute_jac, ref(_cons_vec), res, limits[i], limits[i+1], _first_call_jac, ref(_jac_vals)));
//        }
//        /* Join the threads with the main thread */
//        for(auto &t : threads){
//            t.join();
//        }
//        _first_call_jac = false;
//        return;
//    }
//    for(auto& c_p :_cons)
//    {
//        c = c_p.second.get();
//        auto nb_ins = c->_dim[0];
//        size_t id = 0;
//        if (c->is_linear() && !_first_call_jac) {
////        if (false) {
//            DebugOff("Linear constraint, using stored jacobian!\n");
//            for (size_t i = 0; i<nb_ins; i++) {
//                if (!*c->_all_lazy || !c->_lazy[i]) {
//                    for (size_t j = 0; j<c->get_nb_vars(i); j++) {
//                        res[idx] = _jac_vals[idx];
//                        idx++;
//                    }
//                }
//            }
//        }
//        else {
////            if (_type==nlin_m) {
//                for (auto &v_p: c->get_vars()){
//                    v = v_p.second.first.get();
//                    vid = v->_name;
//                    dfdx = c->get_stored_derivative(vid);
//                    if (dfdx->is_number()) {
//                        for (size_t inst = 0; inst< nb_ins; inst++){
//                            if (!*c->_all_lazy || !c->_lazy[inst]) {
//                                cid = c->_id+id++;
//                                if (v->_is_vector) {
//                                    auto dim = v->get_dim(inst);
//                                    for (size_t j = 0; j<dim; j++) {
//                                        res[idx] = dfdx->_val->at(0);
//                                        _jac_vals[idx] = res[idx];
//                                        idx++;
//                                    }
//                                }
//                                else {
//                                    res[idx] = dfdx->_val->at(0);
//                                    _jac_vals[idx] = res[idx];
//                                    idx++;
//                                }
//                            }
//                        }
//                    }
//                    else {
//                        for (size_t inst = 0; inst< nb_ins; inst++){
//                            if (!*c->_all_lazy || !c->_lazy[inst]) {
//                                cid = c->_id+id++;
//                                if (v->_is_vector) {
//                                    auto dim = v->get_dim(inst);
//                                    if (dfdx->is_matrix()) {
//                                        for (size_t j = 0; j<dim; j++) {
//                                            res[idx] = dfdx->get_val(j,inst);
//                                            _jac_vals[idx] = res[idx];
//                                            idx++;
//                                        }
//                                    }
//                                    else {
//                                        for (size_t j = 0; j<dim; j++) {
//                                            res[idx] = dfdx->get_val(j);//TODO check double indexed funcs
//                                            _jac_vals[idx] = res[idx];
//                                            DebugOff("jac_val["<< idx <<"] = " << _jac_vals[idx] << endl);
//                                            idx++;
//                                        }
//                                    }
//                                }
//                                else {
//                                    res[idx] = dfdx->get_val(inst);
//                                    _jac_vals[idx] = res[idx];
//                                    idx++;
//                                }
//                            }
//                        }
//                    }
//                }
////            }
////            else {
////                for (auto &v_p: c->get_vars()){
////                    v = v_p.second.first.get();
////                    vid = v->_name;
////                    dfdx = c->get_stored_derivative(vid);
////                    for (size_t inst = 0; inst< nb_ins; inst++){
////                        cid = c->_id+inst;
////                        if (v->_is_vector) {
////                            auto dim = v->get_dim(inst);
////                            for (size_t j = 0; j<dim; j++) {
////                                res[idx] = dfdx->eval(inst,j);
////                                _jac_vals[idx] = res[idx];
////                                idx++;
////                            }
////                        }
////                        else {
////                            res[idx] = dfdx->eval(inst);
////                            _jac_vals[idx] = res[idx];
////                            idx++;
////                        }
////                    }
////                }
////            }
//        }
//    }
//    _first_call_jac = false;
//}
//
//
//
//
//#ifdef USE_IPOPT
//void Model::fill_in_var_linearity(Ipopt::TNLP::LinearityType* param_types){
//    size_t vid = 0;
//    bool linear = true;
//    for(auto& vi: _vars)
//    {
//        vid = vi.second->get_id();
//        for (size_t i = 0; i < vi.second->get_dim(); i++) {
////            linear = true;
////            for(auto &c: _v_in_cons[vid + i])
////            {
////                if (!c->is_linear()) {
////                    linear=false;
////                }
////            }
//            if (linear) param_types[vid + i]=Ipopt::TNLP::LINEAR;
//            else param_types[vid + i] = Ipopt::TNLP::NON_LINEAR;
//        }
//    }
//}
//
//
//void Model::fill_in_cstr_linearity(Ipopt::TNLP::LinearityType* const_types){
//    Constraint* c = nullptr;
//    bool lin = false;
//    size_t cid = 0, id = 0;
//    for(auto& c_p :_cons)
//    {
//        c = c_p.second.get();
//        id = 0;
//        if (c->is_linear() || c->is_constant()) {
//            lin = true;
//        }
//        else {
//            lin = false;
//        }
//        auto nb_ins = c->_dim[0];
//        for (size_t i = 0; i< nb_ins; i++){
//            if (!*c->_all_lazy || !c->_lazy[i]) {
//                cid = c->_id+id++;
//                if (lin) {
//                    const_types[cid]=Ipopt::TNLP::LINEAR;
//                }
//                else {
//                    const_types[cid] = Ipopt::TNLP::NON_LINEAR;
//                }
//            }
//        }
//    }
//}
//#endif
//
//
//void Model::fill_in_hess_nnz(int* iRow , int* jCol){
//    size_t idx = 0, idx_all=0, idx_pair=0, vid, vjd;
//    string vi_name, vj_name;
//    shared_ptr<param_> vi;
//    shared_ptr<param_> vj;
//    Constraint* c;
//    for (auto &pairs: _hess_link) {
//        vi_name = pairs.first.first;
//        vj_name = pairs.first.second;
//        vi = (pairs.second.begin())->first->get_var(vi_name);
//        vj = (pairs.second.begin())->first->get_var(vj_name);
//        if (vi_name.compare(vj_name) > 0) {//ONLY STORE LOWER TRIANGULAR PART OF HESSIAN
//            throw invalid_argument("SHOULD BE SORTED CORRECTLY IN FILL_MAPS");
//        }
//        vid = vi->get_id();
//        vjd = vj->get_id();
//        idx_pair = idx;
////        auto max_f_idx = 0;
//        for (auto &f_pair:pairs.second) {
////            auto f_idx = 0;
////            idx = idx_pair;
////        auto f_pair = *pairs.second.begin();
//            auto f = f_pair.first;
//            if (f->_is_constraint) {
//                c = (Constraint*)f;
//            }
//            auto d2f = f_pair.second;
//            size_t nb_inst = f->_dim[0];
//            for (size_t inst = 0; inst<nb_inst; inst++) {
//                if (!(f->_is_constraint && *c->_all_lazy && c->_lazy[inst])) {
//                    
//                    if (d2f->is_matrix()) {
//                        for (size_t i = 0; i < d2f->_dim[0]; i++) {
//                            for (size_t j = i; j < d2f->_dim[1]; j++) {
//                                idx_all++;
//                                iRow[idx] = vid + vi->get_id_inst(i);
//                                jCol[idx] = vjd + vj->get_id_inst(j);
//                                idx++;
////                                f_idx++;
//                            }
//                        }
//                    }
//                    else if(d2f->_is_vector){
//    //                    if (d2f->_dim[0] != d2f->_dim[0]) {
//    //                        throw invalid_argument("error");
//    //                    }
//                        for (size_t j = 0; j < d2f->_dim[0]; j++) {
//                            idx_all++;
//                            iRow[idx] = vid + vi->get_id_inst(j);
//                            jCol[idx] = vjd + vj->get_id_inst(j);
//                            idx++;
////                            f_idx++;
//                        }
//                    }
//                    else {
//                        idx_all++;
//                        iRow[idx] = vid + vi->get_id_inst(inst);
//                        jCol[idx] = vjd + vj->get_id_inst(inst);
//                        idx++;
////                        f_idx++;
//                    }
//                }
//            }
////            if(max_f_idx < f_idx){
////                max_f_idx = f_idx;
////            }
//        }
////        idx = idx_pair+max_f_idx;
//    }
//    if (idx!=_nnz_h) {
//        throw invalid_argument("idx!=_nnz_h");
//    }
//    _hess_vals.resize(idx_all);
//}
//
//void Model::fill_in_hess(const double* x , double obj_factor, const double* lambda, double* res, bool new_x){
//    size_t idx = 0, idx_in = 0, c_inst = 0, idx_pair=0;
//    Constraint* c;
//    bool idx_inc = false;
//    double hess = 0;
//    for (size_t i = 0; i<_nnz_h; i++) {
//        res[i] = 0;
//    }
//    if (new_x) {
//        set_x(x);
//        compute_funcs();
//    }
//    if (_first_call_hess) {
//        for (auto &pairs: _hess_link) {
//            idx_pair = idx;
////            auto max_f_idx = 0;
//            for (auto &f_pair:pairs.second) {
////                auto f_idx = 0;
////                idx = idx_pair;
//                idx_inc = false;
//                auto f = f_pair.first;
//                if (f->_is_constraint) {
//                    c = (Constraint*)f;
//                }
//                auto d2f = f_pair.second;
//                size_t nb_inst = f->_dim[0];
//                size_t id_inst = 0;
//                for (size_t inst = 0; inst<nb_inst; inst++) {
//                    if (f->_is_constraint) {
//                        if (!*c->_all_lazy || !c->_lazy[inst]) {
//                            if (c->is_nonlinear()) {
//                                idx_inc = false;
//                                c_inst = c->get_id_inst(id_inst++);
//    //                            if(f->is_nonlinear()){
//                                if (d2f->is_matrix()) {
//                                    for (size_t i = 0; i < d2f->_dim[0]; i++) {
//                                        for (size_t j = i; j < d2f->_dim[1]; j++) {
//                                            hess = d2f->get_val(i,j);
//                                            _hess_vals[idx_in++] = hess;
//                                            res[idx++] += lambda[c->_id + c_inst] * hess;
//    //                                        f_idx++;
//                                            idx_inc = true;
//                                        }
//                                    }
//                                }
//                                else if(d2f->_is_vector){
//                                    for (size_t j = 0; j < d2f->_dim[0]; j++) {
//                                        hess = d2f->get_val(j);
//                                        _hess_vals[idx_in++] = hess;
//                                        res[idx] += lambda[c->_id + c_inst] * hess;
//                                        idx++;
//    //                                    f_idx++;
//                                        idx_inc = true;
//                                        //                    }
//                                    }
//                                }
//                                else {
//                                    if (d2f->is_number()) {
//                                        hess = d2f->_val->at(0);
//                                    }
//                                    else {
//                                        hess = d2f->_val->at(inst);
//                                    }
//                                    _hess_vals[idx_in++] = hess;
//                                    res[idx] += lambda[c->_id + c_inst] * hess;
//                                }
//                            }
//                            else {
//                                idx_inc = false;
//                                c_inst = c->get_id_inst(id_inst++);
//                                //                            if(f->is_nonlinear()){
//                                if (d2f->is_matrix()) {
//                                    for (size_t i = 0; i < d2f->_dim[0]; i++) {
//                                        for (size_t j = i; j < d2f->_dim[1]; j++) {
//                                            hess = d2f->eval(i,j);
//                                            _hess_vals[idx_in++] = hess;
//                                            res[idx++] += lambda[c->_id + c_inst] * hess;
//                                            //                                        f_idx++;
//                                            idx_inc = true;
//                                        }
//                                    }
//                                }
//                                else if(d2f->_is_vector){
//                                    for (size_t j = 0; j < d2f->_dim[0]; j++) {
//                                        hess = d2f->eval(j);
//                                        _hess_vals[idx_in++] = hess;
//                                        res[idx] += lambda[c->_id + c_inst] * hess;
//                                        idx++;
//                                        //                                    f_idx++;
//                                        idx_inc = true;
//                                        //                    }
//                                    }
//                                }
//                                else {
//                                    if (d2f->is_number()) {
//                                        hess = d2f->_val->at(0);
//                                    }
//                                    else {
//                                        hess = d2f->eval(inst);
//                                    }
//                                    _hess_vals[idx_in++] = hess;
//                                    res[idx] += lambda[c->_id + c_inst] * hess;
//                                }
//                            }
//                        }
//                        else {
//                            idx_inc = true;
//                        }
//                    }
//                    else {
//                        if (d2f->is_matrix()) {
//                            for (size_t i = 0; i < d2f->_dim[0]; i++) {
//                                for (size_t j = i; j < d2f->_dim[1]; j++) {
//                                    hess = d2f->eval(i,j);
//                                    _hess_vals[idx_in++] = hess;
//                                    res[idx++] += obj_factor * hess;
//                                    //                                        f_idx++;
//                                    idx_inc = true;
//                                }
//                            }
//                        }
//                        else if(d2f->_is_vector){
//                            for (size_t j = 0; j < d2f->_dim[0]; j++) {
//                                hess = d2f->eval(j);
//                                _hess_vals[idx_in++] = hess;
//                                res[idx] += obj_factor * hess;
//                                idx++;
////                                f_idx++;
//                                idx_inc = true;
//                            }
//                        }
//                        else {
//                            hess = d2f->eval();
//                            _hess_vals[idx_in++] = hess;
//                            res[idx] += obj_factor * hess;
//                        }
//                    }
//                    if (!idx_inc) {
//                        idx++;
////                        f_idx++;
//                    }
//                }
////                if(max_f_idx < f_idx){
////                    max_f_idx = f_idx;
////                }
//            }
////            idx = idx_pair+max_f_idx;
//        }
//        _first_call_hess = false;
//        return;
//    }
//    if ((_type==lin_m || _type==quad_m)) { /* No need to recompute Hessian for quadratic models, used stored values */
//        size_t id_inst = 0;
//        for (auto &pairs: _hess_link) {
//            idx_pair = idx;
////            auto max_f_idx = 0;
//            for (auto &f_pair:pairs.second) {
////                auto f_idx = 0;
////                idx = idx_pair;
//                idx_inc = false;
//                auto f = f_pair.first;
//                if (f->_is_constraint) {
//                    c = (Constraint*)f;
//                }
//                auto d2f = f_pair.second;
//                size_t nb_inst = f->_dim[0];
//                id_inst = 0;
//                for (size_t inst = 0; inst<nb_inst; inst++) {
//                    if (f->_is_constraint) {
//                        if (!*c->_all_lazy || !c->_lazy[inst]) {
//                            idx_inc = false;
//                            c_inst = c->get_id_inst(id_inst++);
//                            if (d2f->is_matrix()) {
//                                for (size_t i = 0; i < d2f->_dim[0]; i++) {
//                                    for (size_t j = i; j < d2f->_dim[1]; j++) {
//                                        res[idx] += lambda[c->_id + c_inst] * _hess_vals[idx_in++];
//                                        idx++;
////                                        f_idx++;
//                                        idx_inc = true;
//                                    }
//                                }
//                            }
//                            else if(d2f->_is_vector){
//                                for (size_t j = 0; j < d2f->_dim[0]; j++) {
//                                    res[idx] += lambda[c->_id + c_inst] * _hess_vals[idx_in++];
//                                    idx++;
////                                    f_idx++;
//                                    idx_inc = true;
//                                }
//                            }
//                            else {
//                                res[idx] += lambda[c->_id + c_inst] * _hess_vals[idx_in++];
//                            }
//                        }
//                        else {
//                            idx_inc = true;
//                        }
//                    }
//                    else {
//                        if(d2f->_is_vector){
//                            for (size_t j = 0; j < d2f->_dim[0]; j++) {
//                                res[idx] += obj_factor * _hess_vals[idx_in++];
//                                idx++;
////                                f_idx++;
//                                idx_inc = true;
//                            }
//                        }
//                        else {
//                            res[idx] += obj_factor * _hess_vals[idx_in++];
//                        }
//                    }
//                    if (!idx_inc) {
//                        idx++;
////                        f_idx++;
//                    }
//                }
////                if(max_f_idx < f_idx){
////                    max_f_idx = f_idx;
////                }
//            }
////            idx = idx_pair+max_f_idx;
//        }
//        return;
//    }
//    for (auto &pairs: _hess_link) {
//        idx_pair = idx;
////        auto max_f_idx = 0;
//        for (auto &f_pair:pairs.second) {
////            auto f_idx = 0;
////            idx = idx_pair;
//            idx_inc = false;
//            auto f = f_pair.first;
//            if (f->_is_constraint) {
//                c = (Constraint*)f;
//            }
//            auto d2f = f_pair.second;
//            size_t nb_inst = f->_dim[0];
//            size_t id_inst = 0;
//            for (size_t inst = 0; inst<nb_inst; inst++) {
//                if (f->_is_constraint) {
//                    if (!*c->_all_lazy || !c->_lazy[inst]) {
//                        idx_inc = false;
//                        c_inst = c->get_id_inst(id_inst++);
//                        if (c->is_quadratic()) {
//                            if (d2f->is_matrix()) {
//                                for (size_t i = 0; i < d2f->_dim[0]; i++) {
//                                    for (size_t j = i; j < d2f->_dim[1]; j++) {
//                                        res[idx++] += lambda[c->_id + c_inst] * _hess_vals[idx_in++];
////                                        f_idx++;
//                                        idx_inc = true;
//                                    }
//                                }
//                            }
//                            else if(d2f->_is_vector){
//                                for (size_t j = 0; j < d2f->_dim[0]; j++) {
//                                    res[idx] += lambda[c->_id + c_inst] * _hess_vals[idx_in++];
//                                    idx++;
////                                    f_idx++;
//                                    idx_inc = true;
//                                }
//                            }
//                            else {
//                                res[idx] += lambda[c->_id + c_inst] * _hess_vals[idx_in++];
//                            }
//                        }
//                        else if (d2f->is_matrix()) {
//                            if (c->is_nonlinear()) {
//                                for (size_t i = 0; i < d2f->_dim[0]; i++) {
//                                    for (size_t j = i; j < d2f->_dim[1]; j++) {
//                                        hess = d2f->get_val(i,j);
//                                        _hess_vals[idx_in++] = hess;
//                                        res[idx++] += lambda[c->_id + c_inst] * hess;
//    //                                    f_idx++;
//                                        idx_inc = true;
//                                    }
//                                }
//                            }
//                            else {
//                                for (size_t i = 0; i < d2f->_dim[0]; i++) {
//                                    for (size_t j = i; j < d2f->_dim[1]; j++) {
//                                        hess = d2f->eval(i,j);
//                                        _hess_vals[idx_in++] = hess;
//                                        res[idx++] += lambda[c->_id + c_inst] * hess;
//                                        //                                    f_idx++;
//                                        idx_inc = true;
//                                    }
//                                }
//                            }
//                        }
//                        else if(d2f->_is_vector){
//                            if (c->is_nonlinear()) {
//                                for (size_t j = 0; j < d2f->_dim[0]; j++) {
//                                    hess = d2f->get_val(j);
//                                    _hess_vals[idx_in++] = hess;
//                                    res[idx] += lambda[c->_id + c_inst] * hess;
//                                    idx++;
//    //                                f_idx++;
//                                    idx_inc = true;
//                                }
//                            }
//                            else {
//                                for (size_t j = 0; j < d2f->_dim[0]; j++) {
//                                    hess = d2f->eval(j);
//                                    _hess_vals[idx_in++] = hess;
//                                    res[idx] += lambda[c->_id + c_inst] * hess;
//                                    idx++;
//                                    //                                f_idx++;
//                                    idx_inc = true;
//                                }
//
//                            }
//                        }
//                        else {
//                            if (d2f->is_number()) {
//                                hess = d2f->_val->at(0);
//                            }
//                            else {
//                                if (c->is_nonlinear()) {
//                                    hess = d2f->_val->at(inst);
//                                }
//                                else {
//                                    hess = d2f->eval(inst);
//                                }
//                            }
//                            _hess_vals[idx_in++] = hess;
//                            res[idx] += lambda[c->_id + c_inst] * hess;
//                        }
//                    }
//                    else {
//                        idx_inc = true;
//                    }
//                }
//                else {
//                    if (_obj.is_quadratic()) {
//                        if(d2f->_is_vector){
//                            for (size_t j = 0; j < d2f->_dim[0]; j++) {
//                                res[idx] += obj_factor * _hess_vals[idx_in++];
//                                idx++;
////                                f_idx++;
//                                idx_inc = true;
//                            }
//                        }
//                        else {
//                            res[idx] += obj_factor * _hess_vals[idx_in++];
//                        }
//                    }
//                    else if(d2f->_is_vector){
//                        for (size_t j = 0; j < d2f->_dim[0]; j++) {
//                            //                    for (size_t j = i; j < (pairs.second.begin())->second->_dim[1]; j++) {
//                            hess = d2f->eval(j);
//                            _hess_vals[idx_in++] = hess;
//                            res[idx] += obj_factor * hess;
//                            idx++;
////                            f_idx++;
//                            idx_inc = true;
//                            //                    }
//                        }
//                    }
//                    else {
//                        hess = d2f->eval();
//                        _hess_vals[idx_in++] = hess;
//                        res[idx] += obj_factor * hess;
//                    }
//                }
//                if (!idx_inc) {
//                    idx++;
////                    f_idx++;
//                }
//            }
////            if(max_f_idx < f_idx){
////                max_f_idx = f_idx;
////            }
//        }
////        idx = idx_pair+max_f_idx;
//    }
//}
//
//
//
//
//
//
//void Model::reset_funcs() {
//    for (auto& f:_nl_funcs) {
//        f->reset_val();
//    }
//    for (auto& c:_cons) {
//        c.second->reset_val();
//    }
//    _obj.reset_val();
//}
//
///* Compute derivatives and nonzeros in Jacobian and Hessian */
//void Model::fill_in_maps() {
//    string vi_name, vj_name;
//    param_* vi;
//    param_* vj;
//
//    _built = true;
//    
//
//    if (_obj._new) {
//        _obj.compute_derivatives();
//        if (!_obj.is_linear()) {
//            for (auto &vi_p: _obj.get_vars()) {
//                vi = vi_p.second.first.get();
//                vi_name = vi_p.first;
//    //            vid = vi->get_id();
//                auto df = _obj.get_stored_derivative(vi->_name);
//                for (auto &vj_p: df->get_vars()) {
//                    vj = vj_p.second.first.get();
//                    vj_name = vj_p.first;
//    //                vjd = vj->get_id();
//                    if (vi_name.compare(vj_name) < 0) {//ONLY STORE LOWER TRIANGULAR PART OF HESSIAN
//                        _hess_link[make_pair<>(vi_name,vj_name)].insert(make_pair<>(&_obj,_obj.get_stored_derivative(vi->_name)->get_stored_derivative(vj->_name).get()));
//                    }
//                    else {
//                        _hess_link[make_pair<>(vj_name,vi_name)].insert(make_pair<>(&_obj,_obj.get_stored_derivative(vj->_name)->get_stored_derivative(vi->_name).get()));
//                    }
//    //                for (int inst = 0; inst<vi->get_dim(); inst++) {
//    //                    vid = vi->get_id();
//    //                    vid_inst = vid + vi->get_id_inst(inst);
//    //                    vjd = vj->get_id();
//    //                    vjd_inst = vjd + vj->get_id_inst(inst);
//    //                    _hess.insert(make_pair<>(vid_inst, vjd_inst));
//    //                }
//                }
//            }
//        }
//    }
//    Constraint* c = NULL;
//    for(auto& c_p :_cons)
//    {
//        c = c_p.second.get();
////        c->print();
//        if (c->_new) {
//            c->compute_derivatives();
////            if (_type==nlin_m) {
//                for (auto &df_p:*c->get_dfdx()) {
//                    auto df = df_p.second;
//                        DebugOff(df->to_str() << endl);
//                        if (df->get_expr() || _type==nlin_m) {
//                            df_p.second = embed(df);
//                        }
//                        else {
//                            embed(df);
//                        }
//                        for (auto &df2_p:*df_p.second->get_dfdx()) {
////                            if (df2_p.second->get_expr()) {
//                                df2_p.second = embed(df2_p.second);
////                            }
//                        }
//                    }
////            }
//            if (!c->is_linear()) {
//                for (auto &vi_p: c->get_vars()) {
//                    vi = vi_p.second.first.get();
//                    vi_name = vi_p.first;
//                    auto df = c->get_stored_derivative(vi->_name);
//                    for (auto &vj_p: df->get_vars()) {
//                        vj = vj_p.second.first.get();
//                        vj_name = vj_p.first;
//                        if (vi_name.compare(vj_name) <= 0) {//ONLY STORE LOWER TRIANGULAR PART OF HESSIAN
//                            _hess_link[make_pair<>(vi_name,vj_name)].insert(make_pair<>(c, c->get_stored_derivative(vi->_name)->get_stored_derivative(vj->_name).get()));
//                        }
//                        else {
//                            _hess_link[make_pair<>(vj_name,vi_name)].insert(make_pair<>(c, c->get_stored_derivative(vj->_name)->get_stored_derivative(vi->_name).get()));
//                        }
//                    }
//                }
//            }
//            DebugOff(c->to_str() << endl);
//        }
//    }
////    print_nl_functions();
//}
//
//
//void Model::fill_in_duals(double* lambda, double* z_L, double* z_U){
//    for (auto &cp: _cons) {
//        size_t idx = 0;
////        for (size_t inst = 0; inst < cp.second->_dim[0]; inst++) {
////            if (!*cp.second->_all_lazy || !cp.second->_lazy[inst]) {
////                lambda[cp.second->_id + idx++] = 0;
////            }
////
////        }
//        idx = 0;
//        for (size_t inst = 0; inst < cp.second->_dual.size(); inst++) {
//            if (!*cp.second->_all_lazy || !cp.second->_lazy[inst]) {
//                auto index = cp.second->_id + idx++;
//                lambda[index] = cp.second->_dual[inst];
//            }
////            else
////            {
////                lambda[index] = 100;
////            }
//        }
//    }
//    for (auto &vp: _vars) {
//        auto nb_inst = vp.second->get_dim();
//        auto vid = vp.second->get_id();
//        for (size_t inst = 0; inst < nb_inst; inst++) {
//            auto id_inst = vp.second->get_id_inst(inst);
//            z_L[vid + id_inst] = vp.second->_l_dual[inst];
//            z_U[vid + id_inst] = vp.second->_u_dual[inst];
////            z_L[vp.second->get_id() + vp.second->get_id_inst(inst)] = 0;
////            z_U[vp.second->get_id() + vp.second->get_id_inst(inst)] = 0;
//        }
//    }
//    
//}
//
//void Model::fill_in_var_init(double* x) {
//    size_t vid;
//    param_* v;
//    for(auto& v_p: _vars)
//    {
//        v = v_p.second;
//        vid = v->get_id();
//        switch (v->get_intype()) {
//            case float_: {
//                auto real_var = (var<float>*)v;
//                for (size_t i = 0; i < v->get_dim(); i++) {
////                    vid_inst = vid + v->get_id_inst(i);
//                    x[vid+i] = (double)real_var->eval(i);
//                }
//                break;
//            }
//            case long_:{
//                auto real_var = (var<long double>*)v;
//                for (size_t i = 0; i < v->get_dim(); i++) {
////                    vid_inst = vid + v->get_id_inst(i);
//                    x[vid+i] = (double)real_var->eval(i);
//                }
//                break;
//            }
//            case double_:{
//                auto real_var = (var<double>*)v;
//                for (size_t i = 0; i < v->get_dim(); i++) {
////                    vid_inst = vid + v->get_id_inst(i);
//                    x[vid+i] = real_var->eval(i);
//                }
//                break;
//            }
//            case integer_:{
//                auto real_var = (var<int>*)v;
//                for (size_t i = 0; i < v->get_dim(); i++) {
////                    vid_inst = vid + v->get_id_inst(i);
//                    x[vid+i] = (double)real_var->eval(i);
//                }
//                break;
//            }
//            case short_:{
//                auto real_var = (var<short>*)v;
//                for (size_t i = 0; i < v->get_dim(); i++) {
////                    vid_inst = vid + v->get_id_inst(i);
//                    x[vid+i] = (double)real_var->eval(i);
//                }
//                break;
//            }
//            case binary_:{
//                auto real_var = (var<bool>*)v;
//                for (size_t i = 0; i < v->get_dim(); i++) {
////                    vid_inst = vid + v->get_id_inst(i);
//                    x[vid+i] = (double)real_var->eval(i);
//                }
//                break;
//            }
//            default:
//                break;
//        } ;
//    }
//}
//
//void Model::fill_in_cstr_bounds(double* g_l ,double* g_u) {
//    size_t cid = 0;
//    Constraint* c = NULL;
//    for(auto& c_p :_cons)
//    {
//        c = c_p.second.get();
//        switch (c->get_type()) {
//            case eq:{
//                auto nb_ins = c->_dim[0];
//                size_t idx= 0;
//                for (size_t i = 0; i< nb_ins; i++){
//                    if (!*c->_all_lazy || !c->_lazy[i]) {
//                        cid = c->_id+idx++;
//                        g_l[cid] = c->_rhs;
//                        g_u[cid] = c->_rhs;
//                    }
//                }
//                break;
//            }
//            case leq:{
//                auto nb_ins = c->_dim[0];
//                size_t idx= 0;
//                for (size_t i = 0; i< nb_ins; i++){
//                    if (!*c->_all_lazy || !c->_lazy[i]) {
//                        cid = c->_id+idx++;
//                        g_l[cid] = numeric_limits<double>::lowest();
//                        g_u[cid] = c->_rhs;
//                    }
//                }
//                break;
//            }
//            case geq:{
//                auto nb_ins = c->_dim[0];
//                size_t idx= 0;
//                for (size_t i = 0; i< nb_ins; i++){
//                    if (!*c->_all_lazy || !c->_lazy[i]) {
//                        cid = c->_id+idx++;
//                        g_l[cid] = c->_rhs;
//                        g_u[cid] = numeric_limits<double>::max();
//                    }
//                }
//                break;
//            }
//            default:
//                throw invalid_argument("Undefined constraint type!\n");
//                exit(-1);
//                break;
//        }
//    }
//}
//
//
//void Model::embed(expr& e){
//    switch (e.get_type()) {
//        case uexp_c:{
//            auto ue = (uexpr*)&e;
//            auto f = ue->_son;
//            bool found_cpy = false;
//            auto name = f->to_str();
//            if (name.back()=='T') {
//                name = name.substr(0,name.size()-2);
//                if (_nl_funcs_map.count(name)>0) {
//                    auto cpy = _nl_funcs_map.at(name);
//                    f->_val = cpy->_val;
//                    f->_evaluated = true;
//                    found_cpy = true;
//                }
//            }
//            else {
//                auto name1 = "[["+name+"]]\u1D40";
//                if (_nl_funcs_map.count(name1)>0) {
//                    auto cpy = _nl_funcs_map.at(name1);
//                    f->_val = cpy->_val;
//                    f->_evaluated = true;
//                    found_cpy = true;
//                }
//                auto name2 = name+"\u1D40";
//                if (_nl_funcs_map.count(name2)>0) {
//                    auto cpy = _nl_funcs_map.at(name2);
//                    f->_val = cpy->_val;
//                    f->_evaluated = true;
//                    found_cpy = true;
//                }
//            }
//            if (!found_cpy) {
//                auto f_p = _nl_funcs_map.insert(make_pair<>(f->to_str(), f));
//                if (f_p.second) {
//                    embed(f);
//                    _nl_funcs.push_back(f);
//                    DebugOff(f->to_str() << endl);
//                    //                f->_val = make_shared<vector<double>>();
//                    //                f->_val->resize(f->_dim[0]);
//                }
//                else {
//                    ue->_son = f_p.first->second;
//                }
//            }
//            break;
//        }
//        case bexp_c:{
//            auto be = (bexpr*)&e;
//            auto f = be->_lson;
//            bool found_cpy = false;
//            auto name = f->to_str();
//            if (name.back()=='T') {
//                name = name.substr(0,name.size()-2);
//                if (_nl_funcs_map.count(name)>0) {
//                    auto cpy = _nl_funcs_map.at(name);
//                    f->_val = cpy->_val;
//                    f->_evaluated = true;
//                    found_cpy = true;
//                }
//            }
//            else {
//                auto name1 = "[["+name+"]]\u1D40";
//                if (_nl_funcs_map.count(name1)>0) {
//                    auto cpy = _nl_funcs_map.at(name1);
//                    f->_val = cpy->_val;
//                    f->_evaluated = true;
//                    found_cpy = true;
//                }
//                auto name2 = name+"\u1D40";
//                if (_nl_funcs_map.count(name2)>0) {
//                    auto cpy = _nl_funcs_map.at(name2);
//                    f->_val = cpy->_val;
//                    f->_evaluated = true;
//                    found_cpy = true;
//                }
//            }
//            if (!found_cpy) {
//                auto f_p = _nl_funcs_map.insert(make_pair<>(f->to_str(), f));
//                if (f_p.second) {
//                    embed(f);
//                    DebugOff(f->to_str() << endl);
//                    _nl_funcs.push_back(f);
//                    //                f->_val = make_shared<vector<double>>();
//                    //                f->_val->resize(f->_dim[0]);
//                }
//                else {
//                    be->_lson = f_p.first->second;
//                }
//            }
//            f = be->_rson;
//            found_cpy = false;
//            name = f->to_str();
//            if (name.back()=='T') {
//                name = name.substr(0,name.size()-2);
//                if (_nl_funcs_map.count(name)>0) {
//                    auto cpy = _nl_funcs_map.at(name);
//                    f->_val = cpy->_val;
//                    f->_evaluated = true;
//                    found_cpy = true;
//                }
//            }
//            else {
//                auto name1 = "[["+name+"]]\u1D40";
//                if (_nl_funcs_map.count(name1)>0) {
//                    auto cpy = _nl_funcs_map.at(name1);
//                    f->_val = cpy->_val;
//                    f->_evaluated = true;
//                    found_cpy = true;
//                }
//                auto name2 = name+"\u1D40";
//                if (_nl_funcs_map.count(name2)>0) {
//                    auto cpy = _nl_funcs_map.at(name2);
//                    f->_val = cpy->_val;
//                    f->_evaluated = true;
//                    found_cpy = true;
//                }
//            }
//            if (!found_cpy) {
//                auto f_p = _nl_funcs_map.insert(make_pair<>(f->to_str(), f));
//                if (f_p.second) {
//                    embed(f);
//                    DebugOff(f->to_str() << endl);
//                    _nl_funcs.push_back(f);
//                    //                f->_val = make_shared<vector<double>>();
//                    //                f->_val->resize(f->_dim[0]);
//                }
//                else {
//                    be->_rson = f_p.first->second;
//                }
//            }
//            break;
//        }
//        default:
//            break;
//    }
//}
//
//shared_ptr<func_> Model::embed(shared_ptr<func_> f){
//    DebugOff(f->to_str() << endl);
//    for (auto &p_t: f->get_lterms()) {
//        if (p_t.second._coef->is_function()) {
//            auto cf = (func_*)p_t.second._coef;
//            auto exp = cf->get_expr();
//            if (exp) {
//                embed(*exp);
//            }
//            if(cf->is_matrix()){
//                auto newf = embed(make_shared<func_>(*cf));
//                cf->_val = newf->_val;
//            }
//        }
//    }
//    for (auto &p_t: f->get_qterms()) {
//        if (p_t.second._coef->is_function()) {
//            auto cf = (func_*)p_t.second._coef;
//            auto exp = cf->get_expr();
//            if (exp) {
//                embed(*exp);
//            }
//            if(cf->is_matrix()){
//                auto newf = embed(make_shared<func_>(*cf));
//                cf->_val = newf->_val;
//            }
//        }
//    }
//    for (auto &p_t: f->get_pterms()) {
//        if (p_t.second._coef->is_function()) {
//            auto cf = (func_*)p_t.second._coef;
//            auto exp = cf->get_expr();
//            if (exp) {
//                embed(*exp);
//            }
//        }
//    }
//    if (f->get_cst()->is_function()) {
//        auto c = (func_*) f->get_cst();
//        auto exp = c->get_expr();
//        if (exp) {
//            embed(*exp);
//        }
//        if(c->is_matrix()){
//            auto newf = embed(make_shared<func_>(*c));
//            c->_val = newf->_val;
//        }
//    }
//    if (f->get_expr()) {
//        embed(*f->get_expr());
//    }
//    bool found_cpy = false;
//    auto name = f->to_str();
//    if (name.back()=='T') {
//        name = name.substr(0,name.size()-2);
//        if (_nl_funcs_map.count(name)>0) {
//            auto cpy = _nl_funcs_map.at(name);
//            f->_val = cpy->_val;
//            f->_evaluated = true;
//            found_cpy = true;
//        }
//    }
//    else {
//        auto name1 = "[["+name+"]]\u1D40";
//        if (_nl_funcs_map.count(name1)>0) {
//            auto cpy = _nl_funcs_map.at(name1);
//            f->_val = cpy->_val;
//            f->_evaluated = true;
//            found_cpy = true;
//        }
//        auto name2 = name+"\u1D40";
//        if (_nl_funcs_map.count(name2)>0) {
//            auto cpy = _nl_funcs_map.at(name2);
//            f->_val = cpy->_val;
//            f->_evaluated = true;
//            found_cpy = true;
//        }
//    }
//    if (!found_cpy) {
//        auto f_p = _nl_funcs_map.insert(make_pair<>(f->to_str(), f));
//        if (f_p.second) {
//            _nl_funcs.push_back(f_p.first->second);
//            DebugOff(f->to_str() << endl);
//            f->allocate_mem();
//            return f;
//            //        f_p.first->second->_val = make_shared<vector<double>>();
//            //        f_p.first->second->_val->resize(f_p.first->second->_dim[0]);
//        }
////        if (f->_new) {
////            f_p.first->second = f;
////            return f;
////        }
//        if (f->_dim[0] > f_p.first->second->_dim[0]) {
//            *f_p.first->second = *f;
//        }
//        else if (f->_dfdx->size()>0) {
//            *f_p.first->second = *f;
//        }
////        f_p.first->second->allocate_mem();
//        return f_p.first->second;
//    }
//    return f;
//}
//
//
//
//void Model::print_nl_functions() const{
//    cout << "Number of atomic functions = " << _nl_funcs.size();
//    cout << endl;
//        for (auto& f: _nl_funcs){
//            f->print_symbolic(false,false);
//            f->print();
//            cout << endl;
//        }
//    cout << endl;    
//}
//
//void Model::print_solution(bool only_discrete) const{
//    for (auto &v_pair:_vars) {
//        auto v = v_pair.second;
//        switch (v->get_intype()) {
//            case float_: {
//                auto real_var = (var<float>*)v;
//                if(!only_discrete || real_var->_is_relaxed)
//                    real_var->param::print(true);
//                break;
//            }
//            case long_:{
//                auto real_var = (var<long double>*)v;
//                if(!only_discrete || real_var->_is_relaxed)
//                    real_var->param::print(true);
//                break;
//            }
//            case double_:{
//                auto real_var = (var<double>*)v;
//                DebugOff(real_var->get_name() << " in:" << endl);
//                if(!only_discrete || real_var->_is_relaxed)
//                    real_var->param::print(true);
//                DebugOff(";" << endl);
//                break;
//            }
//            case integer_:{
//                auto real_var = (var<int>*)v;
//                real_var->param::print(true);
//                break;
//            }
//            case short_:{
//                auto real_var = (var<short>*)v;
//                real_var->param::print(true);
//                break;
//            }
//            case binary_:{
//                auto real_var = (var<bool>*)v;
//                real_var->param::print(true);
//                break;
//            }
//            default:
//                break;
//        } ;
//    }
//}
//
//void Model::round_solution(){
//    for (auto &v_pair:_vars) {
//        auto v = v_pair.second;
//        switch (v->get_intype()) {
//            case float_: {
//                auto real_var = (var<float>*)v;
//                if(real_var->_is_relaxed){
//                    for (size_t i = 0; i < real_var->get_dim(); i++) {
//                        (real_var->_val->at(i)) = round(real_var->_val->at(i));
//                    }
//                }
//                break;
//            }
//            case long_:{
//                auto real_var = (var<long double>*)v;
//                if(real_var->_is_relaxed){
//                    for (size_t i = 0; i < real_var->get_dim(); i++) {
//                        (real_var->_val->at(i)) = round(real_var->_val->at(i));
//                    }
//                }
//                break;
//            }
//            case double_:{
//                auto real_var = (var<double>*)v;
//                if(real_var->_is_relaxed){
//                    for (size_t i = 0; i < real_var->get_dim(); i++) {
//                        if(real_var->_val->at(i)>0.1){
//                            real_var->_val->at(i) = 1;
//                        }
//                        else {
//                            real_var->_val->at(i) = 0;
//                        }
//                    }
//                }
//                break;
//            }
//            case integer_:{
//                break;
//            }
//            case short_:{
//                break;
//            }
//            case binary_:{
//                break;
//            }
//            default:
//                break;
//        } ;
//    }
//    for (auto &v_p:_bin_vars) {
//        auto bin_var = v_p.second;
//        auto real_var = (var<double>*)get_var_ptr(v_p.first);
//        for (size_t i = 0; i < real_var->get_dim(); i++) {
//            if(round(real_var->_val->at(i))==1){
//                bin_var._val->at(i) = true;
//            }
//            else{
//                bin_var._val->at(i) = false;
//            }
//        }
//    }
//}
//
//void Model::add_integrality(){
//    func_ new_obj;
//    for (auto &v_pair:_vars) {
//        auto v = v_pair.second;
//        switch (v->get_intype()) {
//            case float_: {
//                auto real_var = (var<float>*)v;
//                if(real_var->_is_relaxed){
//                    param<int> rhs("rhs");
//                    for (size_t i = 0; i < real_var->get_dim(); i++) {
//                        rhs.set_val(i,round(real_var->_val->at(i)));
//                    }
//                    new_obj += power(*real_var - rhs,2);
//                }
//                break;
//            }
//            case long_:{
//                auto real_var = (var<long double>*)v;
//                if(real_var->_is_relaxed){
//                    param<int> rhs("rhs");
//                    for (size_t i = 0; i < real_var->get_dim(); i++) {
//                        rhs.set_val(i,round(real_var->_val->at(i)));
//                    }
//                    new_obj += power(*real_var - rhs,2);
//                }
//                break;
//            }
//            case double_:{
//                auto real_var = (var<double>*)v;
//                Constraint Fix_int("Fix_int_"+real_var->_name);
//                Fix_int += *real_var - power((*real_var),2);
//                add(Fix_int<=0);
//                break;
//            }
//            case integer_:{
//                //                auto real_var = (var<int>*)v;
//                //                real_var->param::print(true);
//                break;
//            }
//            case short_:{
//                //                auto real_var = (var<short>*)v;
//                //                real_var->param::print(true);
//                break;
//            }
//            case binary_:{
//                //                auto real_var = (var<bool>*)v;
//                //                real_var->param::print(true);
//                break;
//            }
//            default:
//                break;
//        } ;
//    }
//}
//
//void Model::add_round_solution_obj(bool balance_obj){
//    func_ new_obj;
//    for (auto &v_pair:_vars) {
//        auto v = v_pair.second;
//        v->_new = true;
//        switch (v->get_intype()) {
//            case float_: {
//                auto real_var = (var<float>*)v;
//                if(real_var->_is_relaxed){
//                    param<int> rhs("rhs");
//                    for (size_t i = 0; i < real_var->get_dim(); i++) {
//                        rhs.set_val(i,round(real_var->_val->at(i)));
//                    }
//                    new_obj += power(*real_var - rhs,2);
//                }
//                break;
//            }
//            case long_:{
//                auto real_var = (var<long double>*)v;
//                if(real_var->_is_relaxed){
//                    param<int> rhs("rhs");
//                    for (size_t i = 0; i < real_var->get_dim(); i++) {
//                        rhs.set_val(i,round(real_var->_val->at(i)));
//                    }
//                    new_obj += power(*real_var - rhs,2);
//                }
//                break;
//            }
//            case double_:{
//                double rhs = 0;
//                auto real_var = (var<double>*)v;
//                if(real_var->_is_relaxed){
//                    param<double> coef("coef"+v->_name);
////                    auto rev_indices = real_var->get_rev_indices();
//                    for (size_t i = 0; i < real_var->get_dim(); i++) {
////                        auto idx = real_var->get_id_inst(i);
//                        if ((real_var->_val->at(i))>=0.1) {
////                            new_obj -= 500000*(*real_var)(rev_indices->at(i));
////                            coef.set_val(rev_indices->at(i), -500000);
//                            coef = -2;
//                            rhs += 2;
//                        }
//                        else{
////                            new_obj += 500000*(*real_var)(rev_indices->at(i));
////                            coef.set_val(rev_indices->at(i), 500000);
//                            coef = 2;
////                            rhs = 20000;
//                        }
//                        
//                    }
//                    new_obj += rhs + product(coef,(*real_var));
////                    new_obj.print_expanded();
////                    new_obj += power(rhs - *real_var,2);
////                    new_obj += sum(rhs) - sum(*real_var);
//                }
////                else{
//////                    param<int> rhs("rhs");
////                    auto rev_indices = real_var->get_rev_indices();
////                    for (size_t i = 0; i < real_var->get_dim(); i++) {
////                        new_obj += 1e-2*power(real_var->_val->at(i) - (*real_var)(rev_indices->at(i)),2);
//////                        rhs.set_val(rev_indices->at(i),real_var->_val->at(i));
////                    }
//////                    new_obj += constant<>(1e-2).tr()*power(rhs - (*real_var),2);
////                    //                    new_obj += power(rhs - *real_var,2);
////                    //                    new_obj += sum(rhs) - sum(*real_var);
////                }
//                break;
//            }
//            case integer_:{
////                auto real_var = (var<int>*)v;
////                real_var->param::print(true);
//                break;
//            }
//            case short_:{
////                auto real_var = (var<short>*)v;
////                real_var->param::print(true);
//                break;
//            }
//            case binary_:{
////                auto real_var = (var<bool>*)v;
////                real_var->param::print(true);
//                break;
//            }
//            default:
//                break;
//        } ;
//    }
//    if (balance_obj) {
//        _obj = _obj*(1./_obj_val) + new_obj;
//    }
//    else {
//        _obj = new_obj;
////        for(auto& p: _cons){
////            p.second->_new = true;
////            p.second->_dfdx->clear();
////        }
////        _hess_link.clear();
//
//    }
////    _obj = new_obj;
////    _obj.print_expanded();
//    _obj._dfdx->clear();
//    _obj._new = true;
//    for(auto& p: _cons){
//        p.second->_new = true;
//        p.second->_dfdx->clear();
//    }
//    _hess_link.clear();
//}
//
//void Model::add_round_solution_cuts(){
//    for (auto &v_pair:_vars) {
//        auto v = v_pair.second;
//        switch (v->get_intype()) {
//            case float_: {
//                auto real_var = (var<float>*)v;
//                if(real_var->_is_relaxed){
//                    param<int> rhs("rhs");
//                    for (size_t i = 0; i < real_var->get_dim(); i++) {
//                        rhs.set_val(i,round(real_var->_val->at(i)));
//                    }
//                    Constraint Fix_int("Fix_int_"+real_var->_name);
//                    Fix_int += *real_var - rhs;
//                    add(Fix_int==0);
//                }
//                break;
//            }
//            case long_:{
//                auto real_var = (var<long double>*)v;
//                if(real_var->_is_relaxed){
//                    param<int> rhs("rhs");
//                    for (size_t i = 0; i < real_var->get_dim(); i++) {
//                        rhs.set_val(i,round(real_var->_val->at(i)));
//                    }
//                    Constraint Fix_int("Fix_int_"+real_var->_name);
//                    Fix_int += *real_var - rhs;
//                    add(Fix_int==0);
//                }
//                break;
//            }
//            case double_:{
//                auto real_var = (var<double>*)v;
//                if(real_var->_is_relaxed){
//                    param<> rhs("rhs");
//                    for (size_t i = 0; i < real_var->get_dim(); i++) {
//                        rhs.set_val(i,round(real_var->_val->at(i)));
//                    }
//                    Constraint Fix_int("Fix_int_"+real_var->_name);
//                    Fix_int += power((rhs - *real_var),2);
//                    add(Fix_int<=1e-12);
//                }
//                break;
//            }
//            case integer_:{
//                auto real_var = (var<int>*)v;
//                real_var->param::print(true);
//                break;
//            }
//            case short_:{
//                auto real_var = (var<short>*)v;
//                real_var->param::print(true);
//                break;
//            }
//            case binary_:{
//                auto real_var = (var<bool>*)v;
//                real_var->param::print(true);
//                break;
//            }
//            default:
//                break;
//        } ;
//    }
//}
//
//void Model::print_symbolic(){
//    compute_funcs();
//    for (auto &pair:_vars) {
//        auto v = pair.second;
//        switch (v->get_intype()) {
//            case float_: {
//                break;
//            }
//            case long_:{
//                break;
//            }
//            case double_:{
//                auto real_var = (var<double>*)v;
//                real_var->_lb->eval_all();
//                real_var->_ub->eval_all();
//                real_var->print(true);
//                break;
//            }
//            case complex_:{
//                auto cpx_var = (var<Cpx>*)v;
//                cpx_var->_lb->eval_all();
//                cpx_var->_ub->eval_all();
//                cpx_var->print(true);
//                break;
//            }
//            case integer_:{
//                break;
//            }
//            case short_:{
//                break;
//            }
//            case binary_:{
//                break;
//            }
//                
//        }
//    }
//    cout << "Objective ";
//    _obj.print_symbolic();
//    for(auto& p: _cons){
//        p.second->print_symbolic();
//    }
//}
//
//void Model::print(){
//    compute_funcs();
//    for (auto &pair:_vars) {
//        auto v = pair.second;
//        switch (v->get_intype()) {
//            case float_: {
//                break;
//            }
//            case long_:{
//                break;
//            }
//            case double_:{
//                auto real_var = (var<double>*)v;
//                real_var->_lb->eval_all();
//                real_var->_ub->eval_all();
//                real_var->print(true);
//                break;
//            }
//            case integer_:{
//                break;
//            }
//            case short_:{
//                break;
//            }
//            case binary_:{
//                break;
//            }
//                
//            case complex_:{
//                auto cpx_var = (var<Cpx>*)v;
//                cpx_var->_lb->eval_all();
//                cpx_var->_ub->eval_all();
//                cpx_var->print(true);
//                break;
//            }
//        }
//    }
//    cout << "Objective";
//    _obj.print();
////    exit(1);
//    for(auto& p: _cons){
//        p.second->print();
//    }
//}
//
//void Model::print_constraints() const{
//    for(auto& p: _cons){
//        p.second->print_symbolic();
//    }
//}
//
//pair<func_*, ObjectiveType> gravity::max(const func_& f){
//    f._val->resize(1);
//    return make_pair<>((func_*)&f,maximize);
//};
//
//pair<func_*, ObjectiveType> gravity::min(const func_& f){
//    f._val->resize(1);
//    return make_pair<>((func_*)&f,minimize);
//};
//
//
//
//void Model::replace_integers(){
//    bool has_int = false;
//    //this->relax();
//    //this->unrelax();
//    for (auto &v_p:this->_vars_name) {
//        if (v_p.second->is_integer() || v_p.second->is_binary()) {
//            has_int = true;
//            auto v = v_p.second;
//            auto new_v = new var<double>(v_p.second->_name, 0,1);
//            new_v->copy(*v);
//            new_v->_is_relaxed = true;
//            new_v->_val->resize(new_v->get_dim());
//            if (v->get_intype()==integer_) {
//                auto double_var = (var<int>*)v;
//                for (size_t i = 0; i < double_var->get_dim(); i++) {
//                    new_v->_val->at(i) = double_var->_val->at(i);
//                }
//            }
//            if (v->get_intype()==short_) {
//                auto double_var = (var<short>*)v;
//                for (size_t i = 0; i < double_var->get_dim(); i++) {
//                    new_v->_val->at(i) = double_var->_val->at(i);
//                }
//            }
//            if (v->get_intype()==binary_) {
//                auto double_var = (var<bool>*)v;
//                this->_bin_vars[v_p.second->get_vec_id()] = *double_var;
//                for (size_t i = 0; i < double_var->get_dim(); i++) {
//                    new_v->_val->at(i) = double_var->_val->at(i);
//                }
//            }
//            v_p.second = new_v;
//        }
//    }
//    for (auto &v_p:this->_vars) {
//        if (v_p.second->is_integer() || v_p.second->is_binary()) {
//            auto name = v_p.second->_name;
//            delete v_p.second;
//            v_p.second = this->get_var_ptr(name);
//        }
//    }
//    if(has_int){
//        this->_obj.relax(this->_vars);
//        for (auto &c_p: this->_cons) {
//            c_p.second->relax(this->_vars);
//        }
//    }
//    
//}
//
//
//void Model::add_on_off(const Constraint& c, var<bool>& on){
//    if (c.get_ftype() != lin_) {
//        cerr << "Nonlinear constraint.\n";
//        exit(-1);
//    }
//    Constraint res(c.get_name() + "_on/off");
////    double b;
//    //    for(auto it: orig_q->_coefs) {
//    //        v = getparam_<double>(it.first);
//    //        if (!v->is_bounded_below() || !v->is_bounded_above()) {
//    //            cerr << "Variable " << v->_name << " in constraint " << c._name << " does not have finite bounds.\n";
//    //            exit(1);
//    //        }
//    //        if (c.get_type() == leq || c.get_type() == eq) {
//    //            if (it.second < 0) res -= it.second*v->get_lb_off()*(1-on);
//    //            else res -= it.second*v->get_ub_off()*(1-on);
//    //        }
//    //        else{ // geq
//    //            if (it.second < 0) res -= it.second*v->get_ub_off()*(1-on);
//    //            else res -= it.second*v->get_lb_off()*(1-on);
//    //        }
//    //    }
//    //    if (c.get_type() == eq) {
//    //        Constraint res2(c.get_name() + "_on/off2");
//    //        for(auto it: orig_q->_coefs) {
//    //            v = getparam_<double>(it.first);
//    //            if (it.second < 0) res2 -= it.second*v->get_ub_off()*(1-on);
//    //            else res2 -= it.second*v->get_lb_off()*(1-on);
//    //        }
//    //        res2 += *orig_q;
//    //        res2 -= b*on;
//    //        res2 >= 0;
//    //        addConstraint(res2);
//    //    }
//    //    res += *orig_q;
//    //    res -= orig_q->get_const();
//    //    res -= b*on;
//    //    if (c.get_type() == eq or c.get_type() == leq) res <= 0;
//    //    else res >= 0;
//    add_constraint(res);
//}
//
//void Model::add_on_off(var<>& v, var<bool>& on){
//    //    if(v.get_ub() != v.get_ub_off()) {
//    //        Constraint UB(v._name + "_UB_on/off");
//    //        UB += v - v.get_ub() * on - (1 - on) * v.get_ub_off();
//    //        UB <= 0;
//    //        addConstraint(UB);
//    //    }
//    //    if(v.get_lb() != v.get_lb_off()) {
//    //        Constraint LB(v._name + "_LB_on/off");
//    //        LB += v - v.get_lb() * on - (1 - on) * v.get_lb_off();
//    //        LB >= 0;
//    //        addConstraint(LB);
//    //    }
//}
//
//
//
//void Model::add_on_off_McCormick(std::string name, var<>& v, var<>& v1, var<>& v2, var<bool>& on) {
//    //    Constraint MC1(name+"_McCormick1");
//    //    MC1 += v;
//    //    MC1 -= v1.get_lb()*v2 + v2.get_lb()*v1 - v1.get_lb()*v2.get_lb();
//    //    MC1 >= 0;
//    //    add_on_off(MC1, on);
//    //    Constraint MC2(name+"_McCormick2");
//    //    MC2 += v;
//    //    MC2 -= v1.get_ub()*v2 + v2.get_ub()*v1 - v1.get_ub()*v2.get_ub();
//    //    MC2 >= 0;
//    //    add_on_off(MC2, on);
//    //    Constraint MC3(name+"_McCormick3");
//    //    MC3 += v;
//    //    MC3 -= v1.get_lb()*v2 + v2.get_ub()*v1 - v1.get_lb()*v2.get_ub();
//    //    MC3 <= 0;
//    //    add_on_off(MC3, on);
//    //    Constraint MC4(name+"_McCormick4");
//    //    MC4 += v;
//    //    MC4 -= v1.get_ub()*v2 + v2.get_lb()*v1 - v1.get_ub()*v2.get_lb();
//    //    MC4 <= 0;
//    //    add_on_off(MC4, on);
//}
//
//
//
