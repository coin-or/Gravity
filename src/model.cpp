//
//  model.cpp
//  Gravity
//
//  Created by Hijazi, Hassan (Data61, Canberra City) on 6/5/17.
//
//

#include <gravity/model.h>


using namespace std;
using namespace gravity;

/** Constructor */
//@{
Model::Model(){
    _nnz_g = 0;
    _nnz_h = 0;
};
//@}

/* Destructor */
Model::~Model(){
    for (auto &vp:_vars) {
        delete vp.second;
    }
//    for (auto &cp:_cons) {
//        delete cp.second;
//    }
};


/* Accessors */
size_t Model::get_nb_vars() const{
    return _nb_vars;
};




size_t Model::get_nb_cons() const{
    size_t n = 0;
    for (auto &cp:_cons) {
        n += cp.second->_nb_instances;
    }
    return n;
};


size_t Model::get_nb_nnz_g(){
    _nnz_g = 0;
    for (auto &cp:_cons) {
        auto c = cp.second;
        auto nb_inst = c->_nb_instances;
        for (unsigned inst = 0; inst<nb_inst; inst++) {
            _nnz_g += c->get_nb_vars(inst);
        }
    }
    return _nnz_g;
};


/* Return the number of nonzeros in the lower left part of the hessian */
size_t Model::get_nb_nnz_h(){
    size_t idx = 0;
    bool idx_inc = false;
    for (auto &pairs: _hess_link) {
        auto f_pair = *pairs.second.begin();
//        for (auto &f_pair:pairs.second) {
            idx_inc = false;
            auto f = f_pair.first;
            auto d2f = f_pair.second;
            size_t nb_inst = f->_nb_instances;
            for (unsigned inst = 0; inst<nb_inst; inst++) {
                if (d2f->_is_matrix) {
                    for (unsigned i = 0; i < d2f->_dim[0]; i++) {
                        for (unsigned j = i; j < d2f->_dim[1]; j++) {
                            idx++;
                        }
                    }
                }
                else if(d2f->_is_vector){
                    for (unsigned j = 0; j < d2f->_dim[0]; j++) {
                        idx++;
                    }
                }
                else {
                    idx++;
                }
            }
        }
//    }
    _nnz_h = idx;
    return _nnz_h;
};


void Model::add_indices(const string& constr_name, const node_pairs& indices){
    auto c = _cons_name.at(constr_name);
    map<unsigned, shared_ptr<Constraint>>  new_cons;
    c->add_indices_in(indices);
    for (auto &cc : _cons) {
        if (cc.second->_id > c->_id) {
            cc.second->_id += indices._keys.size();
        }
        new_cons[cc.second->_id] = cc.second;
    }
    _cons = new_cons;
    _nb_cons = get_nb_cons();
}

shared_ptr<Constraint> Model::get_constraint(const string& cname) const{
    return _cons_name.at(cname);
}

param_* Model::get_var(const string& vname) const{
        return _vars_name.at(vname);
}



/* Modifiers */

void Model::init_indices(){// Initialize the indices of all variables involved in the model
    param_* v= nullptr;
    size_t idx = 0;
    for(auto& v_p: _vars)
    {
        v = v_p.second;
        v->set_id(idx);
        for (int i = 0; i < v->get_dim(); i++) {
            idx++;
        }
    }
}

void Model::add_var(param_* v){
    if (v->_is_indexed) {
        return;
    }

    if (_vars_name.count(v->_name) == 0) {
        v->set_id(_nb_vars);
        v->set_vec_id(_vars.size());
        _vars_name[v->_name] = v;
        _vars[v->get_vec_id()] = v;
        _nb_vars += v->get_dim();
    }
};




void Model::del_var(const param_& v){
    auto it = _vars.find(v.get_id());
    if (it!=_vars.end()) {
        _nb_vars -= v.get_dim();
        delete it->second;
        _vars.erase(it);
    }
};


void Model::add_param(param_* v){
    if (_params.count(v->get_id())==0) {
        _nb_params += v->get_dim();
        v->set_id(_params.size());
        v->set_vec_id(_vars.size());
        _params_name[v->get_name()] = v;
        _params[v->get_vec_id()] = v;
    }
};

void Model::add_param(param_& v){
    if (_params.count(v.get_id())==0) {
        _nb_params += v.get_dim();
        auto newv = (param_*)copy(v);
        v.set_id(_params.size());
        v.set_vec_id(_vars.size());
        newv->set_id(_params.size());
        newv->set_vec_id(_vars.size());
        _params_name[v.get_name()] = newv;
        _params[v.get_vec_id()] = newv;
    }
};


void Model::del_param(const param_& v){
    auto it = _params.find(v.get_id());
    if (it!=_params.end()) {
        _nb_params -= v.get_dim();
        delete it->second;
        _params.erase(it);
    }
};

void Model::add_constraint(const Constraint& c){
    if (_cons_name.count(c.get_name())==0) {
        auto newc = make_shared<Constraint>(c);
        if (newc->is_constant()) {
            switch (newc->_ctype) {
                case leq:
                    if (newc->eval()>c.get_rhs()) {
                        throw invalid_argument("Adding violated constant constraint!\n");
                    }
                    break;
                case geq:
                    if (newc->eval()<c.get_rhs()) {
                        throw invalid_argument("Adding violated constant constraint!\n");
                    }
                    break;
                case eq:
                    if (newc->eval()!=c.get_rhs()) {
                        throw invalid_argument("Adding violated constant constraint!\n");
                    }
                    break;
                default:
                    break;
            }
            DebugOn("WARNING: Adding redundant constant constraint, Gravity will be ignoring it.\n");
            return;
        }

        newc->update_to_str();
//        newc->print();
//        if (newc->is_nonlinear()) {
            for (auto &p_t: newc->get_lterms()) {
                if (p_t.second._coef->is_function()) {
                    auto f = (func_*)p_t.second._coef;
                    auto exp = f->get_expr();
                    if (exp) {
                        embed(*exp);
                    }
                }
            }
            for (auto &p_t: newc->get_qterms()) {
                if (p_t.second._coef->is_function()) {
                    auto f = (func_*)p_t.second._coef;
                    auto exp = f->get_expr();
                    if (exp) {
                        embed(*exp);
                    }
                }
            }
            for (auto &p_t: newc->get_pterms()) {
                if (p_t.second._coef->is_function()) {
                    auto f = (func_*)p_t.second._coef;
                    auto exp = f->get_expr();
                    if (exp) {
                        embed(*exp);
                    }
                }
            }
            if (newc->get_cst()->is_function()) {
                auto cc = (func_*) newc->get_cst();
                auto exp = cc->get_expr();
                if (exp) {
                    embed(*exp);
                }
            }
            auto exp = newc->get_expr();
        if (exp) {
            if (exp->is_uexpr()) {
                auto ue = (uexpr*)exp.get();
                auto f = ue->_son;
                bool found_cpy = false;
                auto name = f->to_str();
                if (name.back()=='T') {
                    name = name.substr(0,name.size()-2);
                    if (_nl_funcs_map.count(name)>0) {
                        auto cpy = _nl_funcs_map.at(name);
                        f->_val = cpy->_val;
                        f->_evaluated = true;
                        found_cpy = true;
                    }
                }
                else {
                    auto name1 = "[["+name+"]]^T";
                    if (_nl_funcs_map.count(name1)>0) {
                        auto cpy = _nl_funcs_map.at(name1);
                        f->_val = cpy->_val;
                        f->_evaluated = true;
                        found_cpy = true;
                    }
                    auto name2 = name+"^T";
                    if (_nl_funcs_map.count(name2)>0) {
                        auto cpy = _nl_funcs_map.at(name2);
                        f->_val = cpy->_val;
                        f->_evaluated = true;
                        found_cpy = true;
                    }
                }
                if (!found_cpy) {
                    auto f_p = _nl_funcs_map.insert(make_pair<>(f->to_str(), f));
                    if (f_p.second) {
                        embed(f);
                        _nl_funcs.push_back(f);
                        DebugOff(f->to_str() << endl);
                        //                f->_val = make_shared<vector<double>>();
                        //                f->_val->resize(f->_nb_instances);
                    }
                    else {
                        ue->_son = f_p.first->second;
                    }
                }
            }
            else {
                auto be = (bexpr*)exp.get();
                auto f = be->_lson;
                bool found_cpy = false;
                auto name = f->to_str();
                if (name.back()=='T') {
                    name = name.substr(0,name.size()-2);
                    if (_nl_funcs_map.count(name)>0) {
                        auto cpy = _nl_funcs_map.at(name);
                        f->_val = cpy->_val;
                        f->_evaluated = true;
                        found_cpy = true;
                    }
                }
                else {
                    auto name1 = "[["+name+"]]^T";
                    if (_nl_funcs_map.count(name1)>0) {
                        auto cpy = _nl_funcs_map.at(name1);
                        f->_val = cpy->_val;
                        f->_evaluated = true;
                        found_cpy = true;
                    }
                    auto name2 = name+"^T";
                    if (_nl_funcs_map.count(name2)>0) {
                        auto cpy = _nl_funcs_map.at(name2);
                        f->_val = cpy->_val;
                        f->_evaluated = true;
                        found_cpy = true;
                    }
                }
                if (!found_cpy) {
                    auto f_p = _nl_funcs_map.insert(make_pair<>(f->to_str(), f));
                    if (f_p.second) {
                        embed(f);
                        DebugOff(f->to_str() << endl);
                        _nl_funcs.push_back(f);
                        //                f->_val = make_shared<vector<double>>();
                        //                f->_val->resize(f->_nb_instances);
                    }
                    else {
                        be->_lson = f_p.first->second;
                    }
                }
                f = be->_rson;
                found_cpy = false;
                name = f->to_str();
                if (name.back()=='T') {
                    name = name.substr(0,name.size()-2);
                    if (_nl_funcs_map.count(name)>0) {
                        auto cpy = _nl_funcs_map.at(name);
                        f->_val = cpy->_val;
                        f->_evaluated = true;
                        found_cpy = true;
                    }
                }
                else {
                    auto name1 = "[["+name+"]]^T";
                    if (_nl_funcs_map.count(name1)>0) {
                        auto cpy = _nl_funcs_map.at(name1);
                        f->_val = cpy->_val;
                        f->_evaluated = true;
                        found_cpy = true;
                    }
                    auto name2 = name+"^T";
                    if (_nl_funcs_map.count(name2)>0) {
                        auto cpy = _nl_funcs_map.at(name2);
                        f->_val = cpy->_val;
                        f->_evaluated = true;
                        found_cpy = true;
                    }
                }
                if (!found_cpy) {
                    auto f_p = _nl_funcs_map.insert(make_pair<>(f->to_str(), f));
                    if (f_p.second) {
                        embed(f);
                        DebugOff(f->to_str() << endl);
                        _nl_funcs.push_back(f);
                        //                f->_val = make_shared<vector<double>>();
                        //                f->_val->resize(f->_nb_instances);
                    }
                    else {
                        be->_rson = f_p.first->second;
                    }
                }
            }
        }
//            embed();
//        }
        newc->_id = get_nb_cons();
//        embed(newc);
        _cons_name[c.get_name()] = newc;
        _cons[newc->_id] = newc;
        _built = false;
//        newc->print_expanded();
    }
    else {
//        throw invalid_argument("rename constraint as this name has been used by another one: " + c.to_str());
    }
    int nb_inst = c._nb_instances;
    _nb_cons += nb_inst;
    for (unsigned inst = 0; inst<nb_inst; inst++) {
        _nnz_g += c.get_nb_vars(inst);
    }
    if (_type==lin_m && c.is_quadratic()) {
        _type = quad_m;
    }
    if ((_type==lin_m || _type==quad_m) && (c.is_nonlinear() || c.is_polynomial())) {
        _type = nlin_m;
    }
    if (_is_convex && !c.is_convex()) {
            _is_convex = false;
    }
};



void Model::del_constraint(const Constraint& c){
    _cons_name.erase(c.get_name());
    _cons.erase(c._id);
    _built = false;
    /* TODO: make sure other infos in model are updated */
};

void Model::set_objective(const func_& f, ObjectiveType t) {
    _obj = f;
    _objt = t;
    if (_is_convex) {
        if ((t==maximize && f.is_convex()) || (t==minimize && !f.is_convex())) {
            _is_convex = false;
        }
    }
//    embed(_obj);
    /* TODO make sure the objective receivs the same treatment as the constraints. */
}

void Model::min(const func_& f){
    _obj = f;
    _obj._nb_instances = 1;
    _obj._dim.resize(1);
    _obj._dim[0] = 1;
    _objt = minimize;
}

void Model::max(const func_& f){
    _obj = f;
    _obj._nb_instances = 1;
    _obj._dim.resize(1);
    _obj._dim[0] = 1;
    _objt = maximize;
}


void Model::set_objective(pair<func_*, ObjectiveType> p){
    _obj = *p.first;
    _objt = p.second;
}

void Model::set_objective_type(ObjectiveType t) {
    _objt = t;
}

bool Model::has_violated_constraints(double tol){
    int cid = 0;
    unsigned nb_inst = 0;
    double diff = 0;
    shared_ptr<Constraint> c = nullptr;
    for(auto& c_p: _cons)
    {
        c = c_p.second;
        cid = c->_id;
        switch (c->get_type()) {
            case eq:
                nb_inst = c->_nb_instances;
                for (unsigned inst=0; inst<nb_inst; inst++) {
                    diff = fabs(c->eval(inst) - c->_rhs);
                    if(diff > tol) {
                        DebugOff("Violated equation: ");
                        DebugOff(c->to_str(inst));
                        DebugOff(", violation = "<< diff << endl);
                        return true;
                    }
                }
                break;
            case leq:
                nb_inst = c->_nb_instances;
                for (unsigned inst=0; inst<nb_inst; inst++) {
                    diff = fabs(c->eval(inst) - c->_rhs);
                    if(diff > tol) {
                        DebugOff("Violated inequality: ");
                        DebugOff(c->to_str(inst));
                        DebugOff(", violation = "<< diff << endl);
                        return true;
                    }
                }
                break;
            case geq:
                for (unsigned inst=0; inst<nb_inst; inst++) {
                    diff = fabs(c->eval(inst) - c->_rhs);
                    if(diff < -tol) {
                        DebugOff("Violated inequality: ");
                        DebugOff(c->to_str(inst));
                        DebugOff(", violation = "<< diff << endl);
                        return true;
                    }
                }
                break;
                
            default:
                break;
        }
    }
    return false;
}

void Model::check_feasible(const double* x){
//    int vid = 0;
    //    param_* v = NULL;
//    var<>* var = NULL;
    /* return the structure of the hessian */
//    for(auto& v: _vars)
//    {
//        vid = v->get_idx();
//        var = getparam_<double>(vid);
//        if ((x[vid] - var->get_ub())>1e-6) {
//            cerr << "violated upper bound constraint: ";
//            var->print();
//        }
//        if ((x[vid] - var->get_lb())<-1e-6) {
//            cerr << "violated lower bound constraint: ";
//            var->print();
//        }
//    }
//    int cid = 0;
//    for(auto& c: _cons)
//    {
//        cid = c->get_idx();
//        switch (c->get_type()) {
//            case eq:
//                if(fabs(c->eval(x)-c->_rhs) > 1e-6) {
//                    cerr << "violated constraint: ";
//                    c->print();
//                    printf ("violation = %.10f;\n",(c->eval(x)-c->_rhs));
//                }
//                break;
//            case leq:
//                if((c->eval(x)-c->_rhs) > 1e-6) {
//                    cerr << "violated constraint: ";
//                    c->print();
//                    printf ("violation = %.10f;\n",(c->eval(x)-c->_rhs));
//                }
//                break;
//            case geq:
//                if((c->eval(x)-c->_rhs) < -1e-6) {
//                    cerr << "violated constraint: ";
//                    c->print();
//                    printf ("violation = %.10f;\n",(c->eval(x)-c->_rhs));
//                }
//                break;
//                
//            default:
//                break;
//        }
//    }
}


#ifdef USE_BONMIN
void Model::fill_in_var_types(Bonmin::TMINLP::VariableType* var_types){
    unsigned vid;
    param_* v;
    for(auto& v_p: _vars)
    {
        v = v_p.second;
        vid = v->get_id();
        switch (v->get_intype()) {
            case float_: {
                auto real_var = (var<float>*)v;
                for (int i = 0; i < real_var->get_dim(); i++) {
                   var_types[vid+i] = Bonmin::TMINLP::CONTINUOUS;
                }
                break;
            }
            case long_:{
                auto real_var = (var<long double>*)v;
                for (int i = 0; i < real_var->get_dim(); i++) {
                    var_types[vid+i] = Bonmin::TMINLP::CONTINUOUS;
                }
                break;
            }
            case double_:{
                auto real_var = (var<double>*)v;
                DebugOff(real_var->get_name() << " in:" << endl);
                for (int i = 0; i < real_var->get_dim(); i++) {
                    var_types[vid+i] = Bonmin::TMINLP::CONTINUOUS;
                }
                DebugOff(";" << endl);
                break;
            }
            case integer_:{
                auto real_var = (var<int>*)v;
                for (int i = 0; i < real_var->get_dim(); i++) {
                    var_types[vid+i] = Bonmin::TMINLP::INTEGER;
                }
                break;
            }
            case short_:{
                auto real_var = (var<short>*)v;
                for (int i = 0; i < real_var->get_dim(); i++) {
                    var_types[vid+i] = Bonmin::TMINLP::INTEGER;
                }
                break;
            }
            case binary_:{
                auto real_var = (var<bool>*)v;
                for (int i = 0; i < real_var->get_dim(); i++) {
                    var_types[vid+i] = Bonmin::TMINLP::BINARY;
                }
                break;
            }
            default:
                break;
        } ;
    }
}
#endif


void Model::fill_in_var_bounds(double* x_l ,double* x_u) {
    unsigned vid;
    param_* v;
    for(auto& v_p: _vars)
    {
        v = v_p.second;
        vid = v->get_id();
        switch (v->get_intype()) {
            case float_: {
                auto real_var = (var<float>*)v;
                for (int i = 0; i < real_var->get_dim(); i++) {
                    x_l[vid+i] = (double)real_var->get_lb(i);
                    x_u[vid+i] = (double)real_var->get_ub(i);
                }
                break;
            }
            case long_:{
                auto real_var = (var<long double>*)v;
                for (int i = 0; i < real_var->get_dim(); i++) {
                    x_l[vid+i] = (double)real_var->get_lb(i);
                    x_u[vid+i] = (double)real_var->get_ub(i);
                }
                break;
            }
            case double_:{
                auto real_var = (var<double>*)v;
                DebugOff(real_var->get_name() << " in:" << endl);
                for (int i = 0; i < real_var->get_dim(); i++) {
                    x_l[vid+i] = real_var->get_lb(i);
                    x_u[vid+i] = real_var->get_ub(i);
                    DebugOff("(" << i << ")" << " : [" << x_l[vid+i] << "," << x_u[vid+i] << "]\n");
                }
                DebugOff(";" << endl);
                break;
            }
            case integer_:{
                auto real_var = (var<int>*)v;
                for (int i = 0; i < real_var->get_dim(); i++) {
                    x_l[vid+i] = (double)real_var->get_lb(i);
                    x_u[vid+i] = (double)real_var->get_ub(i);
                }
                break;
            }
            case short_:{
                auto real_var = (var<short>*)v;
                for (int i = 0; i < real_var->get_dim(); i++) {
                    x_l[vid+i] = (double)real_var->get_lb(i);
                    x_u[vid+i] = (double)real_var->get_ub(i);
                }
                break;
            }
            case binary_:{
                auto real_var = (var<bool>*)v;
                for (int i = 0; i < real_var->get_dim(); i++) {                    
                    x_l[vid+i] = (double)real_var->get_lb(i);
                    x_u[vid+i] = (double)real_var->get_ub(i);
                }
                break;
            }
            default:
                break;
        } ;
    }
    //    cout << "idx = " << idx << endl;
}

void Model::set_x(const double* x){
    unsigned vid;
    param_* v;
    for(auto& v_p: _vars)
    {
        v = v_p.second;
        vid = v->get_id();
        switch (v->get_intype()) {
            case float_: {
                auto real_var = (var<float>*)v;
                for (int i = 0; i < real_var->get_dim(); i++) {                    
                    real_var->set_val(i, x[vid+i]);
                }
                break;
            }
            case long_:{
                auto real_var = (var<long double>*)v;
                for (int i = 0; i < real_var->get_dim(); i++) {
                    real_var->set_val(i, x[vid+i]);
                }
                break;
            }
            case double_:{
                auto real_var = (var<double>*)v;
                for (int i = 0; i < real_var->get_dim(); i++) {
                    real_var->set_val(i, x[vid+i]);
                }
                break;
            }
            case integer_:{
                auto real_var = (var<int>*)v;
                for (int i = 0; i < real_var->get_dim(); i++) {
                    real_var->set_val(i, x[vid+i]);
                }
                break;
            }
            case short_:{
                auto real_var = (var<short>*)v;
                for (int i = 0; i < real_var->get_dim(); i++) {
                    real_var->set_val(i, x[vid+i]);
                }
                break;
            }
            case binary_:{
                auto real_var = (var<bool>*)v;
                for (int i = 0; i < real_var->get_dim(); i++) {
                    real_var->set_val(i, x[vid+i]);
                }
                break;
            }
            default:
                break;
        } ;
    }
}

void Model::compute_funcs() {
//    if (_type!=nlin_m) {
//        return;
//    }
    auto it = _nl_funcs.begin();
    while (it!=_nl_funcs.end()) {
        auto f = (*it++);
        DebugOff(f->to_str() << endl);
        if (f->is_constant() && f->_evaluated) {
            continue;
        }
        if (!f->is_constant() && (_type==quad_m || _type==lin_m)) {//these are part of the jacobian, no need to precompute them for quadratic programs
            continue;
        }
        if (!f->_is_matrix && !f->_is_vector) {
            DebugOff(f->to_str()<<endl);
            for (int inst = 0; inst < f->_nb_instances; inst++) {
                f->eval(inst);
            }
        }
        else if(!f->_is_matrix){ //vector
            f->eval_vector();
        }
        else {
            DebugOff(f->to_str()<<endl);
            f->eval_matrix();
        }
        if (f->is_constant()) {
            f->_evaluated = true;
        }
    }    
////        for (int inst = 0; inst < c->_nb_instances; inst++) {
////            c->eval(inst);
////        }
//        if (c->is_nonlinear() || (c->is_linear() && !_first_call_jac)) {
////        if ((c->is_linear() && !_first_call_jac)) {
//            continue;
//        }
//        for (auto &dfdx: *c->get_dfdx()) {
//            if (dfdx.second->_is_matrix) {
//                for (int i = 0; i<dfdx.second->_dim[0]; i++) {
//                    for (int j = 0; j<dfdx.second->_dim[1]; j++) {
//                        dfdx.second->eval(i,j);
//                    }
//                }
//            }
////            else if (dfdx.second->_is_vector) {
////                for (int i = 0; i<dfdx.second->_dim[0]; i++) {
////                    dfdx.second->_val->at(i) = dfdx.second->eval(i);
////                }
////            }
//            else {
//                for (int inst = 0; inst < dfdx.second->_nb_instances; inst++) {
//                    dfdx.second->eval(inst);
//                }
//            }
//            if (c->is_quadratic() && !_first_call_hess) {
//                continue;
//            }
//            for (auto &dfd2x: *dfdx.second->get_dfdx()) {
//                if (dfd2x.second->_is_matrix) {
//                    for (int i = 0; i<dfd2x.second->_dim[0]; i++) {
//                        for (int j = 0; j<dfd2x.second->_dim[1]; j++) {
//                            dfd2x.second->eval(i,j);
//                        }
//                    }
//                }
//                else {
//                    for (int inst = 0; inst < dfd2x.second->_nb_instances; inst++) {
//                        dfd2x.second->eval(inst);
//                    }
//                }
//            }
//        }
//    }
}

void Model::fill_in_obj(const double* x , double& res, bool new_x){
//    if (!new_x) {// IPOPT SEEMS TO BE INCONSISTENT ON NEW_X HERE!!
    if (new_x) {
        set_x(x);
        compute_funcs();
    }
//        set_x(x);
        res = _obj.eval();
        _obj_val = res;
    _obj._new = false;
//    }
//    else {
//        res = _obj_val;
//    }
    DebugOff("Objective = " << to_string(res) << endl);
}

void Model::fill_in_grad_obj(const double* x , double* res, bool new_x){
    param_* v;
    shared_ptr<func_> df;
    unsigned vid, vid_inst, index = 0;
    unique_id v_unique;
    if (new_x) {
        set_x(x);
        compute_funcs();
    }
    for (int i = 0; i<_nb_vars; i++) {
        res[i] = 0;
    }
    if (_first_call_gard_obj) {
        _obj_grad_vals.resize(_obj.get_nb_vars());
        _first_call_gard_obj = false;
    }
    else if (_obj.is_linear()) {
//    else if (false) { /* No need to recompute jacobian for linear objectives */
        for(auto& vi_p: _obj.get_vars())
        {
            v = vi_p.second.first.get();
            vid = v->get_id();
            if (v->_is_vector) {
                for (int i = 0; i < v->get_dim(); i++) {
                    vid_inst = vid + v->get_id_inst(i);
                    res[vid_inst] = _obj_grad_vals[index++];
                }
            }
            else {
                vid_inst = vid + v->get_id_inst();
                res[vid_inst] = _obj_grad_vals[index++];
            }
        }
        return;
    }
    for(auto& vi_p: _obj.get_vars()) {
        v = vi_p.second.first.get();
        vid = v->get_id();
        v_unique = v->_unique_id;
        df = _obj.get_stored_derivative(v_unique);
        if (v->_is_vector) {
            for (int i = 0; i < v->get_dim(); i++) {
                vid_inst = vid + v->get_id_inst(i);
                res[vid_inst] = df->eval(i);
                _obj_grad_vals[index++] =res[vid_inst];
            }
        }
        else {
            vid_inst = vid + v->get_id_inst();
            res[vid_inst] = df->eval();
            _obj_grad_vals[index++] =res[vid_inst];
        }
    }
}

void dummy (){
    
    
}
void compute_constrs(vector<shared_ptr<Constraint>>& v, double* res, int i, int j){
//    DebugOff("Calling compute_constrts with i =  " << i << "and j = "<< j << endl);
    for (unsigned idx = i; idx < j; idx++) {
        auto c = v[idx];
        auto nb_ins = c->_nb_instances;
        for (int inst = 0; inst< nb_ins; inst++){
            res[c->_id+inst] = c->eval(inst);
            DebugOff("Accessing res at position " << c->_id+inst << endl);
            //                _cons_vals[index++] = res[c->_id+inst];
            DebugOff("g[" << to_string(c->_id+inst) << "] = " << to_string(res[c->_id+inst]) << endl);
        }
    }
}

void Model::fill_in_cstr(const double* x , double* res, bool new_x){
    Constraint* c = nullptr;
//    unsigned index = 0;
    if (new_x) {
        set_x(x);
        compute_funcs();
    }
//    if (_type==nlin_m) {
//        for(auto& c_p: _cons)
//        {
//            c = c_p.second.get();
//            auto nb_ins = c->_nb_instances;
//            for (int inst = 0; inst< nb_ins; inst++){
//                res[c->_id+inst] = c->get_val(inst);
////                _cons_vals[index++] = res[c->_id+inst];
//                DebugOff("g[" << to_string(c->_id+inst) << "] = " << to_string(res[c->_id+inst]) << endl);
//            }
//        }
//    }
//    else {
//    vector<Constraint*> cons;
        _cons_vec.clear();
        for(auto& c_p: _cons)
        {
            c_p.second->_new = false;
            _cons_vec.push_back(c_p.second);
//            auto nb_ins = c->_nb_instances;
//            for (int inst = 0; inst< nb_ins; inst++){
////                res[c->_id+inst] = c->get_val(inst);
//                res[c->_id+inst] = c->eval(inst);
//                //                _cons_vals[index++] = res[c->_id+inst];
//                DebugOff("g[" << to_string(c->_id+inst) << "] = " << to_string(res[c->_id+inst]) << endl);
//            }
        }
//    compute_constrs(_cons_vec, res, 0, _cons_vec.size());
    unsigned nr_threads = 4;
    vector<thread> threads;
    /* Split cons into nr_threads parts */
    vector<int> limits = bounds(nr_threads, _cons_vec.size());
    DebugOff("limits size = " << limits.size() << endl);
    for (int i = 0; i < limits.size(); ++i) {
        DebugOff("limits[" << i << "] = " << limits[i] << endl);
    }    
    /* Launch all threads in parallel */
    for (int i = 0; i < nr_threads; ++i) {
        DebugOff("i = " << i << endl);
        DebugOff("limits[" << i << "] = " << limits[i] << endl);
        DebugOff("limits[" << i+1 << "] = " << limits[i+1] << endl);
        threads.push_back(thread(compute_constrs, ref(_cons_vec), res, limits[i], limits[i+1]));
//        threads.push_back(thread(dummy));
    }
    /* Join the threads with the main thread */
    for(auto &t : threads){
        t.join();
    }
    
//    }
//    }
//    else {
//        
//        for(auto& c_p: _cons) {
//            c = c_p.second;
//            auto nb_ins = c->_nb_instances;
//            for (int inst = 0; inst< nb_ins; inst++){
//                res[c->_id+inst] = _cons_vals[index++];
//            }
//        }
//    }
}



void Model::fill_in_jac_nnz(int* iRow , int* jCol){
    size_t idx=0;
    size_t cid = 0;
    size_t vid = 0;
    Constraint* c = NULL;
    param_* v = NULL;
    /* return the structure of the jacobian */
    for(auto& c_p :_cons)
    {
        c = c_p.second.get();
        c->_jac_cstr_idx = idx;
        auto nb_ins = c->_nb_instances;
        for (auto &v_p: c->get_vars()){
            v = v_p.second.first.get();
            vid = v->get_id();
            for (int inst = 0; inst< nb_ins; inst++){
                cid = c->_id+inst;
                if (v->_is_vector) {
                    auto dim = v->get_dim(inst);
                    for (int j = 0; j<dim; j++) {
                        iRow[idx] = cid;
                        jCol[idx] = vid + v->get_id_inst(inst, j);
                        idx++;
                    }
                }
                else {
                    iRow[idx] = cid;
                    jCol[idx] = vid + v->get_id_inst(inst);
                    idx++;
                }
            }
        }
    }
    if (idx!=_nnz_g) {
        throw invalid_argument("idx!=_nnz_g");
    }
}


void compute_jac(vector<shared_ptr<Constraint>>& vec, double* res, int i, int j, bool first_call, vector<double>& jac_vals){
    size_t cid = 0;
    unique_id vid;
    shared_ptr<Constraint> c = NULL;
    param_* v = NULL;
    shared_ptr<func_> dfdx;
    auto idx = vec[i]->_jac_cstr_idx;
    for (unsigned id = i; id < j; id++) {
        c = vec[id];
        auto nb_ins = c->_nb_instances;
        if (c->is_linear() && !first_call) {
            //        if (false) {
            DebugOff("Linear constraint, using stored jacobian!\n");
            for (unsigned i = 0; i<nb_ins; i++) {
                for (unsigned j = 0; j<c->get_nb_vars(i); j++) {
                    res[idx] = jac_vals[idx];
                    idx++;
                }
            }
        }
        else {
            for (auto &v_p: c->get_vars()){
                v = v_p.second.first.get();
                vid = v->_unique_id;
                dfdx = c->get_stored_derivative(vid);
                for (int inst = 0; inst< nb_ins; inst++){
                    cid = c->_id+inst;
                    if (v->_is_vector) {
                        auto dim = v->get_dim(inst);
                        for (int j = 0; j<dim; j++) {
                            res[idx] = dfdx->eval(inst,j);
                            jac_vals[idx] = res[idx];
                            idx++;
                        }
                    }
                    else {
                        res[idx] = dfdx->eval(inst);
                        jac_vals[idx] = res[idx];
                        idx++;
                    }
                }
            }
        }
    }
}

/* Fill the nonzero values in the jacobian */
void Model::fill_in_jac(const double* x , double* res, bool new_x){
//    if (!_first_call_jac && (!new_x || _type==lin_m)) { /* No need to recompute jacobian for linear models */
    if (new_x) {
        set_x(x);
        compute_funcs();
    }

    if (!_first_call_jac && (_type==lin_m)) { /* No need to recompute jacobian for linear models */
        for (int i = 0; i< _nnz_g; i++) {
            res[i] = _jac_vals[i];
        }
        return;
    }
    size_t idx=0;
    size_t cid = 0;
    unique_id vid;
    Constraint* c = NULL;
    param_* v = NULL;
    shared_ptr<func_> dfdx;
//    vector<Constraint*> cons;
    if (_type!=nlin_m) {//Quadratic or Linear
        _cons_vec.clear();
        for(auto& c_p: _cons)
        {
            c_p.second->_new = false;
            _cons_vec.push_back(c_p.second);
        }
        
        compute_jac(_cons_vec, res, 0, _cons_vec.size(), _first_call_jac, _jac_vals);
//        unsigned nr_threads = 4;
//        vector<thread> threads;
//        /* Split cons into nr_threads parts */
//        vector<int> limits = bounds(nr_threads, _cons_vec.size());
//        
//        /* Launch all threads in parallel */
//        for (int i = 0; i < nr_threads; ++i) {
//            threads.push_back(thread(compute_jac, ref(_cons_vec), res, limits[i], limits[i+1], _first_call_jac, ref(_jac_vals)));
//        }
//        /* Join the threads with the main thread */
//        for(auto &t : threads){
//            t.join();
//        }
        _first_call_jac = false;
        return;
    }
    for(auto& c_p :_cons)
    {
        c = c_p.second.get();
        auto nb_ins = c->_nb_instances;
        if (c->is_linear() && !_first_call_jac) {
//        if (false) {
            DebugOff("Linear constraint, using stored jacobian!\n");
            for (unsigned i = 0; i<nb_ins; i++) {
                for (unsigned j = 0; j<c->get_nb_vars(i); j++) {
                    res[idx] = _jac_vals[idx];
                    idx++;
                }
            }
        }
        else {
//            if (_type==nlin_m) {
                for (auto &v_p: c->get_vars()){
                    v = v_p.second.first.get();
                    vid = v->_unique_id;
                    dfdx = c->get_stored_derivative(vid);
                    if (dfdx->is_number()) {
                        for (int inst = 0; inst< nb_ins; inst++){
                            cid = c->_id+inst;
                            if (v->_is_vector) {
                                auto dim = v->get_dim(inst);
                                for (int j = 0; j<dim; j++) {
                                    res[idx] = dfdx->_val->at(0);
                                    _jac_vals[idx] = res[idx];
                                    idx++;
                                }                               
                            }
                            else {
                                res[idx] = dfdx->_val->at(0);
                                _jac_vals[idx] = res[idx];
                                idx++;
                            }
                        }
                    }
                    else {
                        for (int inst = 0; inst< nb_ins; inst++){
                            cid = c->_id+inst;
                            if (v->_is_vector) {
                                auto dim = v->get_dim(inst);
                                if (dfdx->_is_matrix) {
                                    for (int j = 0; j<dim; j++) {
                                        res[idx] = dfdx->get_val(j,inst);
                                        _jac_vals[idx] = res[idx];
                                        idx++;
                                    }
                                }
                                else {
                                    for (int j = 0; j<dim; j++) {
                                        res[idx] = dfdx->get_val(j);//TODO check double indexed funcs
                                        _jac_vals[idx] = res[idx];
                                        idx++;
                                    }
                                }
                            }
                            else {
                                res[idx] = dfdx->get_val(inst);
                                _jac_vals[idx] = res[idx];
                                idx++;
                            }
                        }
                    }
                }
//            }
//            else {
//                for (auto &v_p: c->get_vars()){
//                    v = v_p.second.first.get();
//                    vid = v->_unique_id;
//                    dfdx = c->get_stored_derivative(vid);
//                    for (int inst = 0; inst< nb_ins; inst++){
//                        cid = c->_id+inst;
//                        if (v->_is_vector) {
//                            auto dim = v->get_dim(inst);
//                            for (int j = 0; j<dim; j++) {
//                                res[idx] = dfdx->eval(inst,j);
//                                _jac_vals[idx] = res[idx];
//                                idx++;
//                            }
//                        }
//                        else {
//                            res[idx] = dfdx->eval(inst);
//                            _jac_vals[idx] = res[idx];
//                            idx++;
//                        }
//                    }
//                }
//            }
        }
    }
    _first_call_jac = false;
}




#ifdef USE_IPOPT
void Model::fill_in_var_linearity(Ipopt::TNLP::LinearityType* param_types){
    size_t vid = 0;
    bool linear = true;
    for(auto& vi: _vars)
    {
        vid = vi.second->get_id();
        for (int i = 0; i < vi.second->get_dim(); i++) {
//            linear = true;
//            for(auto &c: _v_in_cons[vid + i])
//            {
//                if (!c->is_linear()) {
//                    linear=false;
//                }
//            }
            if (linear) param_types[vid + i]=Ipopt::TNLP::LINEAR;
            else param_types[vid + i] = Ipopt::TNLP::NON_LINEAR;
        }
    }
}


void Model::fill_in_cstr_linearity(Ipopt::TNLP::LinearityType* const_types){
    Constraint* c = nullptr;
    bool lin = false;
    size_t cid = 0;
    for(auto& c_p :_cons)
    {
        c = c_p.second.get();
        if (c->is_linear() || c->is_constant()) {
            lin = true;
        }
        else {
            lin = false;
        }
        auto nb_ins = c->_nb_instances;
        for (int i = 0; i< nb_ins; i++){
            cid = c->_id +i;
            if (lin) {
                const_types[cid]=Ipopt::TNLP::LINEAR;
            }
            else {
                const_types[cid] = Ipopt::TNLP::NON_LINEAR;
            }
        }
    }
}
#endif


void Model::fill_in_hess_nnz(int* iRow , int* jCol){
    size_t idx = 0, idx_all=0, idx_pair=0, vid, vjd;
    string vi_name, vj_name;
    shared_ptr<param_> vi;
    shared_ptr<param_> vj;
    Constraint* c;
    for (auto &pairs: _hess_link) {
        vi_name = pairs.first.first;
        vj_name = pairs.first.second;
        vi = (pairs.second.begin())->first->get_var(vi_name);
        vj = (pairs.second.begin())->first->get_var(vj_name);
        if (vi_name.compare(vj_name) > 0) {//ONLY STORE LOWER TRIANGULAR PART OF HESSIAN
            throw invalid_argument("SHOULD BE SORTED CORRECTLY IN FILL_MAPS");
        }
        vid = vi->get_id();
        vjd = vj->get_id();
        idx_pair = idx;
        for (auto &f_pair:pairs.second) {
            idx = idx_pair;
//        auto f_pair = *pairs.second.begin();
            auto f = f_pair.first;
            if (f->_is_constraint) {
                c = (Constraint*)f;
            }
            auto d2f = f_pair.second;
            size_t nb_inst = f->_nb_instances;
            for (unsigned inst = 0; inst<nb_inst; inst++) {
                if (d2f->_is_matrix) {
                    for (unsigned i = 0; i < d2f->_dim[0]; i++) {
                        for (unsigned j = i; j < d2f->_dim[1]; j++) {
                            idx_all++;
                            iRow[idx] = vid + vi->get_id_inst(i);
                            jCol[idx] = vjd + vj->get_id_inst(j);
                            idx++;
                        }
                    }
                }
                else if(d2f->_is_vector){
//                    if (d2f->_dim[0] != d2f->_nb_instances) {
//                        throw invalid_argument("error");
//                    }
                    for (unsigned j = 0; j < d2f->_dim[0]; j++) {
                        idx_all++;
                        iRow[idx] = vid + vi->get_id_inst(j);
                        jCol[idx] = vjd + vj->get_id_inst(j);
                        idx++;
                    }
                }
                else {
                    idx_all++;
                    iRow[idx] = vid + vi->get_id_inst(inst);
                    jCol[idx] = vjd + vj->get_id_inst(inst);
                    idx++;
                }
            }
        }
    }
    if (idx!=_nnz_h) {
        throw invalid_argument("idx!=_nnz_h");
    }
    _hess_vals.resize(idx_all);
}
 
void Model::fill_in_hess(const double* x , double obj_factor, const double* lambda, double* res, bool new_x){
    size_t idx = 0, idx_in = 0, c_inst = 0, idx_pair=0;
    Constraint* c;
    bool idx_inc = false;
    double hess = 0;
    for (unsigned i = 0; i<_nnz_h; i++) {
        res[i] = 0;
    }
    if (new_x) {
        set_x(x);
        compute_funcs();
    }
    if (_first_call_hess) {
        for (auto &pairs: _hess_link) {
            idx_pair = idx;
            for (auto &f_pair:pairs.second) {
                idx = idx_pair;
                idx_inc = false;
                auto f = f_pair.first;
                if (f->_is_constraint) {
                    c = (Constraint*)f;
                }
                auto d2f = f_pair.second;
                size_t nb_inst = f->_nb_instances;
                for (unsigned inst = 0; inst<nb_inst; inst++) {
                    if (f->_is_constraint) {
                        c_inst = c->get_id_inst(inst);
                        if (d2f->_is_matrix) {
                            for (unsigned i = 0; i < d2f->_dim[0]; i++) {
                                for (unsigned j = i; j < d2f->_dim[1]; j++) {
                                    hess = d2f->get_val(i,j);
                                    _hess_vals[idx_in++] = hess;
                                    res[idx++] += lambda[c->_id + c_inst] * hess;
                                    idx_inc = true;
                                }
                            }
                        }
                        else if(d2f->_is_vector){
                            for (unsigned j = 0; j < d2f->_dim[0]; j++) {
                                hess = d2f->get_val(j);
                                _hess_vals[idx_in++] = hess;
                                res[idx] += lambda[c->_id + c_inst] * hess;
                                idx++;
                                idx_inc = true;
                                //                    }
                            }
                        }
                        else {
                            if (d2f->is_number()) {
                                hess = d2f->_val->at(0);
                            }
                            else {
                                hess = d2f->_val->at(inst);
                            }
                            _hess_vals[idx_in++] = hess;
                            res[idx] += lambda[c->_id + c_inst] * hess;
                        }
                    }
                    else {
                        if(d2f->_is_vector){
                            for (unsigned j = 0; j < d2f->_dim[0]; j++) {
                                hess = d2f->eval(j);
                                _hess_vals[idx_in++] = hess;
                                res[idx] += obj_factor * hess;
                                idx++;
                                idx_inc = true;
                            }
                        }
                        else {
                            hess = d2f->eval();
                            _hess_vals[idx_in++] = hess;
                            res[idx] += obj_factor * hess;
                        }
                    }
                    if (!idx_inc) {
                        idx++;
                    }
                }
            }
        }
        _first_call_hess = false;
        return;
    }
    if ((_type==lin_m || _type==quad_m)) { /* No need to recompute Hessian for quadratic models, used stored values */
        for (auto &pairs: _hess_link) {
            idx_pair = idx;
            for (auto &f_pair:pairs.second) {
                idx = idx_pair;
                idx_inc = false;
                auto f = f_pair.first;
                if (f->_is_constraint) {
                    c = (Constraint*)f;
                }
                auto d2f = f_pair.second;
                size_t nb_inst = f->_nb_instances;
                for (unsigned inst = 0; inst<nb_inst; inst++) {
                    if (f->_is_constraint) {
                        c_inst = c->get_id_inst(inst);
                        if (d2f->_is_matrix) {
                            for (unsigned i = 0; i < d2f->_dim[0]; i++) {
                                for (unsigned j = i; j < d2f->_dim[1]; j++) {
                                    res[idx] += lambda[c->_id + c_inst] * _hess_vals[idx_in++];
                                    idx++;
                                    idx_inc = true;
                                }
                            }
                        }
                        else if(d2f->_is_vector){
                            for (unsigned j = 0; j < d2f->_dim[0]; j++) {
                                res[idx] += lambda[c->_id + c_inst] * _hess_vals[idx_in++];
                                idx++;
                                idx_inc = true;
                            }
                        }
                        else {
                            res[idx] += lambda[c->_id + c_inst] * _hess_vals[idx_in++];
                        }
                    }
                    else {
                        if(d2f->_is_vector){
                            for (unsigned j = 0; j < d2f->_dim[0]; j++) {
                                res[idx] += obj_factor * _hess_vals[idx_in++];
                                idx++;
                                idx_inc = true;
                            }
                        }
                        else {
                            res[idx] += obj_factor * _hess_vals[idx_in++];
                        }
                    }
                    if (!idx_inc) {
                        idx++;
                    }
                }
            }
        }
        return;
    }
    for (auto &pairs: _hess_link) {
        idx_pair = idx;
        for (auto &f_pair:pairs.second) {
            idx = idx_pair;
            idx_inc = false;
            auto f = f_pair.first;
            if (f->_is_constraint) {
                c = (Constraint*)f;
            }
            auto d2f = f_pair.second;
            size_t nb_inst = f->_nb_instances;
            for (unsigned inst = 0; inst<nb_inst; inst++) {
                if (f->_is_constraint) {
                    c_inst = c->get_id_inst(inst);
                    if (c->is_quadratic()) {
                        if (d2f->_is_matrix) {
                            for (unsigned i = 0; i < d2f->_dim[0]; i++) {
                                for (unsigned j = i; j < d2f->_dim[1]; j++) {
                                    res[idx++] += lambda[c->_id + c_inst] * _hess_vals[idx_in++];
                                    idx_inc = true;
                                }
                            }
                        }
                        else if(d2f->_is_vector){
                            for (unsigned j = 0; j < d2f->_dim[0]; j++) {
                                res[idx] += lambda[c->_id + c_inst] * _hess_vals[idx_in++];
                                idx++;
                                idx_inc = true;
                            }
                        }
                        else {
                            res[idx] += lambda[c->_id + c_inst] * _hess_vals[idx_in++];
                        }
                    }
                    else if (d2f->_is_matrix) {
                        for (unsigned i = 0; i < d2f->_dim[0]; i++) {
                            for (unsigned j = i; j < d2f->_dim[1]; j++) {
                                hess = d2f->get_val(i,j);
                                _hess_vals[idx_in++] = hess;
                                res[idx++] += lambda[c->_id + c_inst] * hess;
                                idx_inc = true;
                            }
                        }
                    }
                    else if(d2f->_is_vector){
                        for (unsigned j = 0; j < d2f->_dim[0]; j++) {
                            hess = d2f->get_val(j);
                            _hess_vals[idx_in++] = hess;
                            res[idx] += lambda[c->_id + c_inst] * hess;
                            idx++;
                            idx_inc = true;
                        }
                    }
                    else {
                        if (d2f->is_number()) {
                            hess = d2f->_val->at(0);
                        }
                        else {
                            hess = d2f->_val->at(inst);
                        }
                        _hess_vals[idx_in++] = hess;
                        res[idx] += lambda[c->_id + c_inst] * hess;
                    }
                }
                else {
                    if (_obj.is_quadratic()) {
                        if(d2f->_is_vector){
                            for (unsigned j = 0; j < d2f->_dim[0]; j++) {
                                res[idx] += obj_factor * _hess_vals[idx_in++];
                                idx++;
                                idx_inc = true;
                            }
                        }
                        else {
                            res[idx] += obj_factor * _hess_vals[idx_in++];
                        }
                    }
                    else if(d2f->_is_vector){
                        for (unsigned j = 0; j < d2f->_dim[0]; j++) {
                            //                    for (unsigned j = i; j < (pairs.second.begin())->second->_dim[1]; j++) {
                            hess = d2f->eval(j);
                            _hess_vals[idx_in++] = hess;
                            res[idx] += obj_factor * hess;
                            idx++;
                            idx_inc = true;
                            //                    }
                        }
                    }
                    else {
                        hess = d2f->eval();
                        _hess_vals[idx_in++] = hess;
                        res[idx] += obj_factor * hess;
                    }
                }
                if (!idx_inc) {
                    idx++;
                }
            }
        }
    }
}






void Model::reset_funcs() {
    for (auto& f:_nl_funcs) {
        f->reset_val();
    }
    for (auto& c:_cons) {
        c.second->reset_val();
    }
    _obj.reset_val();
}

/* Compute derivatives and nonzeros in Jacobian and Hessian */
void Model::fill_in_maps() {
    string vi_name, vj_name;
    param_* vi;
    param_* vj;

    _built = true;
    

    if (_obj._new) {
        _obj.compute_derivatives();
        if (!_obj.is_linear()) {
            for (auto &vi_p: _obj.get_vars()) {
                vi = vi_p.second.first.get();
                vi_name = vi_p.first;
    //            vid = vi->get_id();
                auto df = _obj.get_stored_derivative(vi->_unique_id);
                for (auto &vj_p: df->get_vars()) {
                    vj = vj_p.second.first.get();
                    vj_name = vj_p.first;
    //                vjd = vj->get_id();
                    if (vi_name.compare(vj_name) < 0) {//ONLY STORE LOWER TRIANGULAR PART OF HESSIAN
                        _hess_link[make_pair<>(vi_name,vj_name)].insert(make_pair<>(&_obj,_obj.get_stored_derivative(vi->_unique_id)->get_stored_derivative(vj->_unique_id).get()));
                    }
                    else {
                        _hess_link[make_pair<>(vj_name,vi_name)].insert(make_pair<>(&_obj,_obj.get_stored_derivative(vj->_unique_id)->get_stored_derivative(vi->_unique_id).get()));
                    }
    //                for (int inst = 0; inst<vi->get_dim(); inst++) {
    //                    vid = vi->get_id();
    //                    vid_inst = vid + vi->get_id_inst(inst);
    //                    vjd = vj->get_id();
    //                    vjd_inst = vjd + vj->get_id_inst(inst);
    //                    _hess.insert(make_pair<>(vid_inst, vjd_inst));
    //                }
                }
            }
        }
    }
    Constraint* c = NULL;
    for(auto& c_p :_cons)
    {
        c = c_p.second.get();
        if (c->_new) {
            c->compute_derivatives();
//            if (_type==nlin_m) {
                for (auto &df_p:*c->get_dfdx()) {
                    auto df = df_p.second;
                        DebugOff(df->to_str() << endl);
//                        if (df->get_expr()) {
                            df_p.second = embed(df);
//                        }
                        for (auto &df2_p:*df_p.second->get_dfdx()) {
//                            if (df2_p.second->get_expr()) {
                                df2_p.second = embed(df2_p.second);
//                            }
                        }
                    }
//            }
            c->_val = make_shared<vector<double>>();
            c->_val->resize(c->_nb_instances);
            if (!c->is_linear()) {
                for (auto &vi_p: c->get_vars()) {
                    vi = vi_p.second.first.get();
                    vi_name = vi_p.first;
                    auto df = c->get_stored_derivative(vi->_unique_id);
                    for (auto &vj_p: df->get_vars()) {
                        vj = vj_p.second.first.get();
                        vj_name = vj_p.first;
                        if (vi_name.compare(vj_name) <= 0) {//ONLY STORE LOWER TRIANGULAR PART OF HESSIAN
                            _hess_link[make_pair<>(vi_name,vj_name)].insert(make_pair<>(c, c->get_stored_derivative(vi->_unique_id)->get_stored_derivative(vj->_unique_id).get()));
                        }
                        else {
                            _hess_link[make_pair<>(vj_name,vi_name)].insert(make_pair<>(c, c->get_stored_derivative(vj->_unique_id)->get_stored_derivative(vi->_unique_id).get()));
                        }
                    }
                }
            }
            DebugOff(c->to_str() << endl);
        }
    }
}


void Model::fill_in_duals(double* lambda, double* z_L, double* z_U){
    for (auto &cp: _cons) {
        if(cp.second->_dual.size()==0){
            cp.second->_dual.resize(cp.second->_nb_instances);
        }
        for (unsigned inst = 0; inst < cp.second->_nb_instances; inst++) {
            lambda[cp.second->_id + inst] = 0;
            
        }
        for (unsigned inst = 0; inst < cp.second->_dual.size(); inst++) {
            lambda[cp.second->_id + inst] = cp.second->_dual[inst];
        }
    }
    for (auto &vp: _vars) {
        if(vp.second->_l_dual.size()==0 || vp.second->_u_dual.size()==0 ){
            vp.second->_l_dual.resize(vp.second->get_nb_instances());
            vp.second->_u_dual.resize(vp.second->get_nb_instances());
        }
        auto nb_inst = vp.second->get_nb_instances();
        for (unsigned inst = 0; inst < nb_inst; inst++) {
            z_L[vp.second->get_id() + vp.second->get_id_inst(inst)] = vp.second->_l_dual[inst];
            z_U[vp.second->get_id() + vp.second->get_id_inst(inst)] = vp.second->_u_dual[inst];
//            z_L[vp.second->get_id() + vp.second->get_id_inst(inst)] = 0;
//            z_U[vp.second->get_id() + vp.second->get_id_inst(inst)] = 0;
        }
    }
    
}

void Model::fill_in_var_init(double* x) {
    size_t vid;
    param_* v;
    for(auto& v_p: _vars)
    {
        v = v_p.second;
        vid = v->get_id();
        switch (v->get_intype()) {
            case float_: {
                auto real_var = (var<float>*)v;
                for (int i = 0; i < v->get_dim(); i++) {
//                    vid_inst = vid + v->get_id_inst(i);
                    x[vid+i] = (double)real_var->eval(i);
                }
                break;
            }
            case long_:{
                auto real_var = (var<long double>*)v;
                for (int i = 0; i < v->get_dim(); i++) {
//                    vid_inst = vid + v->get_id_inst(i);
                    x[vid+i] = (double)real_var->eval(i);
                }
                break;
            }
            case double_:{
                auto real_var = (var<double>*)v;
                for (int i = 0; i < v->get_dim(); i++) {
//                    vid_inst = vid + v->get_id_inst(i);
                    x[vid+i] = real_var->eval(i);
                }
                break;
            }
            case integer_:{
                auto real_var = (var<int>*)v;
                for (int i = 0; i < v->get_dim(); i++) {
//                    vid_inst = vid + v->get_id_inst(i);
                    x[vid+i] = (double)real_var->eval(i);
                }
                break;
            }
            case short_:{
                auto real_var = (var<short>*)v;
                for (int i = 0; i < v->get_dim(); i++) {
//                    vid_inst = vid + v->get_id_inst(i);
                    x[vid+i] = (double)real_var->eval(i);
                }
                break;
            }
            case binary_:{
                auto real_var = (var<bool>*)v;
                for (int i = 0; i < v->get_dim(); i++) {
//                    vid_inst = vid + v->get_id_inst(i);
                    x[vid+i] = (double)real_var->eval(i);
                }
                break;
            }
            default:
                break;
        } ;
    }
}

void Model::fill_in_cstr_bounds(double* g_l ,double* g_u) {
    size_t cid = 0;
    Constraint* c = NULL;
    for(auto& c_p :_cons)
    {
        c = c_p.second.get();
        switch (c->get_type()) {
            case eq:{
                auto nb_ins = c->_nb_instances;
                for (int i = 0; i< nb_ins; i++){
                    cid = c->_id+i;
                    g_l[cid] = c->_rhs;
                    g_u[cid] = c->_rhs;
                }
                break;
            }
            case leq:{
                auto nb_ins = c->_nb_instances;
                for (int i = 0; i< nb_ins; i++){
                    cid = c->_id+i;
                    g_l[cid] = numeric_limits<double>::lowest();
                    g_u[cid] = c->_rhs;
                }
                break;
            }
            case geq:{
                auto nb_ins = c->_nb_instances;
                for (int i = 0; i< nb_ins; i++){
                    cid = c->_id+i;
                    g_l[cid] = c->_rhs;
                    g_u[cid] = numeric_limits<double>::max();
                }
                break;
            }
            default:
                throw invalid_argument("Undefined constraint type!\n");
                exit(-1);
                break;
        }
    }
}


void Model::embed(expr& e){
    switch (e.get_type()) {
        case uexp_c:{
            auto ue = (uexpr*)&e;
            auto f = ue->_son;
            bool found_cpy = false;
            auto name = f->to_str();
            if (name.back()=='T') {
                name = name.substr(0,name.size()-2);
                if (_nl_funcs_map.count(name)>0) {
                    auto cpy = _nl_funcs_map.at(name);
                    f->_val = cpy->_val;
                    f->_evaluated = true;
                    found_cpy = true;
                }
            }
            else {
                auto name1 = "[["+name+"]]^T";
                if (_nl_funcs_map.count(name1)>0) {
                    auto cpy = _nl_funcs_map.at(name1);
                    f->_val = cpy->_val;
                    f->_evaluated = true;
                    found_cpy = true;
                }
                auto name2 = name+"^T";
                if (_nl_funcs_map.count(name2)>0) {
                    auto cpy = _nl_funcs_map.at(name2);
                    f->_val = cpy->_val;
                    f->_evaluated = true;
                    found_cpy = true;
                }
            }
            if (!found_cpy) {
                auto f_p = _nl_funcs_map.insert(make_pair<>(f->to_str(), f));
                if (f_p.second) {
                    embed(f);
                    _nl_funcs.push_back(f);
                    DebugOff(f->to_str() << endl);
                    //                f->_val = make_shared<vector<double>>();
                    //                f->_val->resize(f->_nb_instances);
                }
                else {
                    ue->_son = f_p.first->second;
                }
            }
            break;
        }
        case bexp_c:{
            auto be = (bexpr*)&e;
            auto f = be->_lson;
            bool found_cpy = false;
            auto name = f->to_str();
            if (name.back()=='T') {
                name = name.substr(0,name.size()-2);
                if (_nl_funcs_map.count(name)>0) {
                    auto cpy = _nl_funcs_map.at(name);
                    f->_val = cpy->_val;
                    f->_evaluated = true;
                    found_cpy = true;
                }
            }
            else {
                auto name1 = "[["+name+"]]^T";
                if (_nl_funcs_map.count(name1)>0) {
                    auto cpy = _nl_funcs_map.at(name1);
                    f->_val = cpy->_val;
                    f->_evaluated = true;
                    found_cpy = true;
                }
                auto name2 = name+"^T";
                if (_nl_funcs_map.count(name2)>0) {
                    auto cpy = _nl_funcs_map.at(name2);
                    f->_val = cpy->_val;
                    f->_evaluated = true;
                    found_cpy = true;
                }
            }
            if (!found_cpy) {
                auto f_p = _nl_funcs_map.insert(make_pair<>(f->to_str(), f));
                if (f_p.second) {
                    embed(f);
                    DebugOff(f->to_str() << endl);
                    _nl_funcs.push_back(f);
                    //                f->_val = make_shared<vector<double>>();
                    //                f->_val->resize(f->_nb_instances);
                }
                else {
                    be->_lson = f_p.first->second;
                }
            }
            f = be->_rson;
            found_cpy = false;
            name = f->to_str();
            if (name.back()=='T') {
                name = name.substr(0,name.size()-2);
                if (_nl_funcs_map.count(name)>0) {
                    auto cpy = _nl_funcs_map.at(name);
                    f->_val = cpy->_val;
                    f->_evaluated = true;
                    found_cpy = true;
                }
            }
            else {
                auto name1 = "[["+name+"]]^T";
                if (_nl_funcs_map.count(name1)>0) {
                    auto cpy = _nl_funcs_map.at(name1);
                    f->_val = cpy->_val;
                    f->_evaluated = true;
                    found_cpy = true;
                }
                auto name2 = name+"^T";
                if (_nl_funcs_map.count(name2)>0) {
                    auto cpy = _nl_funcs_map.at(name2);
                    f->_val = cpy->_val;
                    f->_evaluated = true;
                    found_cpy = true;
                }
            }
            if (!found_cpy) {
                auto f_p = _nl_funcs_map.insert(make_pair<>(f->to_str(), f));
                if (f_p.second) {
                    embed(f);
                    DebugOff(f->to_str() << endl);
                    _nl_funcs.push_back(f);
                    //                f->_val = make_shared<vector<double>>();
                    //                f->_val->resize(f->_nb_instances);
                }
                else {
                    be->_rson = f_p.first->second;
                }
            }
            break;
        }
        default:
            break;
    }
}

shared_ptr<func_> Model::embed(shared_ptr<func_> f){
    DebugOff(f->to_str() << endl);
    for (auto &p_t: f->get_lterms()) {
        if (p_t.second._coef->is_function()) {
            auto cf = (func_*)p_t.second._coef;
            auto exp = cf->get_expr();
            if (exp) {
                embed(*exp);
            }
        }
    }
    for (auto &p_t: f->get_qterms()) {
        if (p_t.second._coef->is_function()) {
            auto cf = (func_*)p_t.second._coef;
            auto exp = cf->get_expr();
            if (exp) {
                embed(*exp);
            }
        }
    }
    for (auto &p_t: f->get_pterms()) {
        if (p_t.second._coef->is_function()) {
            auto cf = (func_*)p_t.second._coef;
            auto exp = cf->get_expr();
            if (exp) {
                embed(*exp);
            }
        }
    }
    if (f->get_cst()->is_function()) {
        auto c = (func_*) f->get_cst();
        auto exp = c->get_expr();
        if (exp) {
            embed(*exp);
        }
    }
    if (f->get_expr()) {
        embed(*f->get_expr());
    }
    bool found_cpy = false;
    auto name = f->to_str();
    if (name.back()=='T') {
        name = name.substr(0,name.size()-2);
        if (_nl_funcs_map.count(name)>0) {
            auto cpy = _nl_funcs_map.at(name);
            f->_val = cpy->_val;
            f->_evaluated = true;
            found_cpy = true;
        }
    }
    else {
        auto name1 = "[["+name+"]]^T";
        if (_nl_funcs_map.count(name1)>0) {
            auto cpy = _nl_funcs_map.at(name1);
            f->_val = cpy->_val;
            f->_evaluated = true;
            found_cpy = true;
        }
        auto name2 = name+"^T";
        if (_nl_funcs_map.count(name2)>0) {
            auto cpy = _nl_funcs_map.at(name2);
            f->_val = cpy->_val;
            f->_evaluated = true;
            found_cpy = true;
        }
    }
    if (!found_cpy) {
        auto f_p = _nl_funcs_map.insert(make_pair<>(f->to_str(), f));
        if (f_p.second) {
            _nl_funcs.push_back(f_p.first->second);
            DebugOff(f->to_str() << endl);
            return f;
            //        f_p.first->second->_val = make_shared<vector<double>>();
            //        f_p.first->second->_val->resize(f_p.first->second->_nb_instances);
        }
//        if (f->_new) {
//            f_p.first->second = f;
//            return f;
//        }
        if (f->_nb_instances > f_p.first->second->_nb_instances) {
            *f_p.first->second = *f;
        }
        else if (f->_dfdx->size()>0) {
            *f_p.first->second = *f;
        }
        return f_p.first->second;
    }
    return f;
}



void Model::print_nl_functions() const{
    cout << "Number of atomic functions = " << _nl_funcs.size();
    cout << endl;
        for (auto& f: _nl_funcs){
            f->print(false,false);
            cout << endl;
        }
    cout << endl;
}

void Model::print_solution() const{
    
}


void Model::print_expanded() const{
    for(auto& p: _cons){
        p.second->print_expanded();
    }
}

void Model::print_constraints() const{
    for(auto& p: _cons){
        p.second->print();
    }
}

pair<func_*, ObjectiveType> gravity::max(const func_& f){
    return make_pair<>((func_*)&f,maximize);
};

pair<func_*, ObjectiveType> gravity::min(const func_& f){
    return make_pair<>((func_*)&f,minimize);
};


void Model::add_on_off(const Constraint& c, var<bool>& on){
    if (c.get_ftype() != lin_) {
        cerr << "Nonlinear constraint.\n";
        exit(-1);
    }
    Constraint res(c.get_name() + "_on/off");
//    double b;
    //    for(auto it: orig_q->_coefs) {
    //        v = getparam_<double>(it.first);
    //        if (!v->is_bounded_below() || !v->is_bounded_above()) {
    //            cerr << "Variable " << v->_name << " in constraint " << c._name << " does not have finite bounds.\n";
    //            exit(1);
    //        }
    //        if (c.get_type() == leq || c.get_type() == eq) {
    //            if (it.second < 0) res -= it.second*v->get_lb_off()*(1-on);
    //            else res -= it.second*v->get_ub_off()*(1-on);
    //        }
    //        else{ // geq
    //            if (it.second < 0) res -= it.second*v->get_ub_off()*(1-on);
    //            else res -= it.second*v->get_lb_off()*(1-on);
    //        }
    //    }
    //    if (c.get_type() == eq) {
    //        Constraint res2(c.get_name() + "_on/off2");
    //        for(auto it: orig_q->_coefs) {
    //            v = getparam_<double>(it.first);
    //            if (it.second < 0) res2 -= it.second*v->get_ub_off()*(1-on);
    //            else res2 -= it.second*v->get_lb_off()*(1-on);
    //        }
    //        res2 += *orig_q;
    //        res2 -= b*on;
    //        res2 >= 0;
    //        addConstraint(res2);
    //    }
    //    res += *orig_q;
    //    res -= orig_q->get_const();
    //    res -= b*on;
    //    if (c.get_type() == eq or c.get_type() == leq) res <= 0;
    //    else res >= 0;
    add_constraint(res);
}

void Model::add_on_off(var<>& v, var<bool>& on){
    //    if(v.get_ub() != v.get_ub_off()) {
    //        Constraint UB(v._name + "_UB_on/off");
    //        UB += v - v.get_ub() * on - (1 - on) * v.get_ub_off();
    //        UB <= 0;
    //        addConstraint(UB);
    //    }
    //    if(v.get_lb() != v.get_lb_off()) {
    //        Constraint LB(v._name + "_LB_on/off");
    //        LB += v - v.get_lb() * on - (1 - on) * v.get_lb_off();
    //        LB >= 0;
    //        addConstraint(LB);
    //    }
}

void Model::add_McCormick(std::string name, var<>& v, var<>& v1, var<>& v2) {
    //    Constraint MC1(name+"_McCormick1");
    //    MC1 += v;
    //    MC1 -= v1.get_lb()*v2 + v2.get_lb()*v1 - v1.get_lb()*v2.get_lb();
    //    MC1 >= 0;
    //    add_constraint(MC1);
    //    //    MC1.print();
    //    Constraint MC2(name+"_McCormick2");
    //    MC2 += v;
    //    MC2 -= v1.get_ub()*v2 + v2.get_ub()*v1 - v1.get_ub()*v2.get_ub();
    //    MC2 >= 0;
    //    add_constraint(MC2);
    //    //    MC2.print();
    //    Constraint MC3(name+"_McCormick3");
    //    MC3 += v;
    //    MC3 -= v1.get_lb()*v2 + v2.get_ub()*v1 - v1.get_lb()*v2.get_ub();
    //    MC3 <= 0;
    //    add_constraint(MC3);
    //    //    MC3.print();
    //    Constraint MC4(name+"_McCormick4");
    //    MC4 += v;
    //    MC4 -= v1.get_ub()*v2 + v2.get_lb()*v1 - v1.get_ub()*v2.get_lb();
    //    MC4 <= 0;
    //    add_constraint(MC4);
    //    MC4.print();
}


void Model::add_on_off_McCormick(std::string name, var<>& v, var<>& v1, var<>& v2, var<bool>& on) {
    //    Constraint MC1(name+"_McCormick1");
    //    MC1 += v;
    //    MC1 -= v1.get_lb()*v2 + v2.get_lb()*v1 - v1.get_lb()*v2.get_lb();
    //    MC1 >= 0;
    //    add_on_off(MC1, on);
    //    Constraint MC2(name+"_McCormick2");
    //    MC2 += v;
    //    MC2 -= v1.get_ub()*v2 + v2.get_ub()*v1 - v1.get_ub()*v2.get_ub();
    //    MC2 >= 0;
    //    add_on_off(MC2, on);
    //    Constraint MC3(name+"_McCormick3");
    //    MC3 += v;
    //    MC3 -= v1.get_lb()*v2 + v2.get_ub()*v1 - v1.get_lb()*v2.get_ub();
    //    MC3 <= 0;
    //    add_on_off(MC3, on);
    //    Constraint MC4(name+"_McCormick4");
    //    MC4 += v;
    //    MC4 -= v1.get_ub()*v2 + v2.get_lb()*v1 - v1.get_ub()*v2.get_lb();
    //    MC4 <= 0;
    //    add_on_off(MC4, on);
}
