//
//  model.cpp
//  Gravity
//
//  Created by Hijazi, Hassan (Data61, Canberra City) on 6/5/17.
//
//

#include <gravity/model.h>
//#define USEDEBUG
#ifdef USEDEBUG
#define Debug(x) cout << x
#else
#define Debug(x)
#endif
#define DebugOn(x) cout << x
#define DebugOff(x)


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
    for (auto &cp:_cons) {
        delete cp.second;
    }
};


/* Accessors */


size_t Model::get_nb_vars() const{
    return _nb_vars;
};

size_t Model::get_nb_cons() const{
    size_t n = 0;
    for (auto &cp:_cons) {
        n += cp.second->get_nb_instances();
    }
    return n;
};


size_t Model::get_nb_nnz_g() const{
    return _nnz_g;
};

//Split "mem" into "parts", e.g. if mem = 10 and parts = 4 you will have: 0,2,4,6,10
//if possible the function will split mem into equal chuncks, if not
//the last chunck will be slightly larger

std::vector<int> bounds(int parts, int mem) {
    std::vector<int>bnd;
    int delta = mem / parts;
    int reminder = mem % parts;
    int N1 = 0, N2 = 0;
    bnd.push_back(N1);
    for (int i = 0; i < parts; ++i) {
        N2 = N1 + delta;
        if (i == parts - 1)
            N2 += reminder;
        bnd.push_back(N2);
        N1 = N2;
    }
    return bnd;
}


/* Return the number of nonzeros in the lower left part of the hessian */
size_t Model::get_nb_nnz_h() const{
    return _nnz_h;
};



Constraint* Model::get_constraint(const string& cname) const{
    return (Constraint*)&_cons_name.at(cname);
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

    if (_vars_name.count(v->get_name()) == 0) {
        v->set_id(_nb_vars);
        v->set_vec_id(_vars.size());
        _vars_name[v->get_name()] = v;
        _vars[v->get_vec_id()] = v;
        _nb_vars += v->get_dim();
    }
};

void Model::add_var(param_& v){
    if (v._is_indexed) { 
        return;
    }
    if (_vars_name.count(v.get_name())==0) {
        v.set_id(_nb_vars);
        v.set_vec_id(_vars.size());
        auto newv = (param_*)copy(v);
        _vars_name[v.get_name()] = newv;
        _vars[v.get_vec_id()] = newv;
        _nb_vars += v.get_dim();
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
        auto newc = new Constraint(c);
//        newc->print();
//        embed(*newc);
        newc->_id = _nb_cons;
        _cons_name[c.get_name()] = newc;
        _cons[newc->_id] = newc;
    }
    else {
        throw invalid_argument("rename constraint as this name has been used by another one: " + c.to_str());
    }
    int nb_inst = c.get_nb_instances();
    _nb_cons += nb_inst;
    _nnz_g += c.get_nb_vars()*nb_inst;
};



void Model::del_constraint(const Constraint& c){
    //    _cons.erase(c->get_idx());
    assert(false);
};

void Model::set_objective(const func_& f, ObjectiveType t) {
    _obj = f;
    _objt = t;
//    embed(_obj);
}

void Model::set_objective(pair<func_*, ObjectiveType> p){
    _obj = *p.first;
    _objt = p.second;
}

void Model::set_objective_type(ObjectiveType t) {
    _objt = t;
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
//                real_var->print(true);
                for (int i = 0; i < real_var->get_dim(); i++) {
                    x_l[vid+i] = (double)real_var->get_lb(i);
                    x_u[vid+i] = (double)real_var->get_ub(i);
                }
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

void Model::fill_in_obj(const double* x , double& res, bool new_x){
    if (new_x) {
        Debug("new_x!" << endl);
        set_x(x);
#ifdef USEDEBUG
        Debug("x = [ ");
        for (int i = 0; i<_nb_vars; i++) {
            Debug(to_string(x[i]) << " ");
        }
        Debug("]\n");
#endif
    }
    res = _obj.eval();
    DebugOff("Objective = " << to_string(res) << endl);
}

void Model::fill_in_cstr(const double* x , double* res, bool new_x){
    if (new_x) {
        set_x(x);
    }
    DebugOff("x = [ ");
    for (int i = 0; i<_nb_vars; i++) {
        DebugOff(to_string(x[i]) << " ");
    }
    DebugOff("]");
    Constraint* c = nullptr;
    for(auto& c_p: _cons)
    {
        c = c_p.second;
        auto nb_ins = c->get_nb_instances();
        for (int inst = 0; inst< nb_ins; inst++){
            res[c->_id+inst] = c->eval(inst);
            DebugOff("g[" << to_string(c->_id+inst) << "] = " << to_string(res[c->_id+inst]) << endl);
        }
    }
}

/* Fill the nonzero values in the jacobian */
void Model::fill_in_jac(const double* x , double* res, bool new_x){
    if (new_x) {
        set_x(x);
    }
    size_t idx=0;
    size_t cid = 0;
    unique_id vid;
    Constraint* c = NULL;
    param_* v = NULL;
    func_* dfdx;
    for(auto& c_p :_cons)
    {
        c = c_p.second;
        auto nb_ins = c->get_nb_instances();
        for (int inst = 0; inst< nb_ins; inst++){
            cid = c->_id+inst;
            for (auto &v_p: c->get_vars()){
                v = v_p.second.first;
                vid = v->_unique_id;
                if (v->_is_vector) {
                    dfdx = c->get_stored_derivative(vid);
                    for (int j = 0; j<v->get_dim(); j++) {
                        res[idx] = dfdx->eval(j);
                        idx++;
                    }
                }
                else {
                    res[idx] = c->get_stored_derivative(vid)->eval(inst);
                    idx++;
                }
            }
        }
    }
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
        c = c_p.second;
        auto nb_ins = c->get_nb_instances();
        for (int inst = 0; inst< nb_ins; inst++){
            cid = c->_id+inst;
            for (auto &v_p: c->get_vars()){
                v = v_p.second.first;
                vid = v->get_id();
                if (v->_is_vector) {
                    for (int j = 0; j<v->get_dim(); j++) {
                        iRow[idx] = cid;
                        jCol[idx] = vid + v->get_id_inst(j);
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
}



#ifdef USE_IPOPT
void Model::fill_in_var_linearity(Ipopt::TNLP::LinearityType* param_types){
    size_t vid = 0;
    bool linear = true;
    for(auto& vi: _vars)
    {
        vid = vi.second->get_id();
        for (int i = 0; i < vi.second->get_dim(); i++) {
            linear = true;
            for(auto &c: _v_in_cons[vid + i])
            {
                if (!c->is_linear()) {
                    linear=false;
                }
            }
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
        c = c_p.second;
        if (c->is_linear() || c->is_constant()) {
            lin = true;
        }
        else {
            lin = false;
        }
        auto nb_ins = c->get_nb_instances();
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
    size_t idx=0;
    for (auto &hess_i: _hess_link) {
//        if (!hess_i.empty()) {
            for (auto &hess_j: hess_i.second) {
                iRow[idx] = hess_i.first.second;
                jCol[idx] = hess_j.first.second;
                idx++;
            }
//        }
    }
}

void Model::fill_in_hess(const double* x , double obj_factor, const double* lambda, double* res, bool new_x){
    size_t idx = 0;
    unique_id vid, vjd;
    double hess = 0;
    Constraint* c;
    for (auto &hess_i: _hess_link) {
        assert(!hess_i.second.empty());
        vid = hess_i.first.first;
        for (auto &hess_j: hess_i.second) {
            vjd = hess_j.first.first;
            res[idx] = 0;
            for (auto cid: hess_j.second){
                if (cid.first==-1) {
                    hess = _obj.get_stored_derivative(vid)->get_stored_derivative(vjd)->eval();
                    res[idx] += obj_factor * hess;
                }
                else {
                    c = _cons[cid.first];
//                    for (int inst = 0; inst< c->get_nb_instances(); inst++){
                        hess = c->get_stored_derivative(vid)->get_stored_derivative(vjd)->eval(cid.second);
//                        if (hess!=1) {
//                            cout << "ok" << endl;
//                        }
                        res[idx] += lambda[cid.first+cid.second] * hess;
//                    }
                }
            }
            idx++;
        }
    }
}






void Model::fill_in_grad_obj(const double* x , double* res, bool new_x){
    if (new_x) {
        set_x(x);
    }
    param_* v;
    func_* df;
    unsigned vid, vid_inst;
    unique_id v_unique;
    for (int i = 0; i<_nb_vars; i++) {
        res[i] = 0;
    }
    
    for(auto& vi_p: _obj.get_vars())
    {
        v = vi_p.second.first;
        vid = v->get_id();
        v_unique = v->_unique_id;
        df = _obj.get_stored_derivative(v_unique);
        if (v->_is_vector) {
            for (int i = 0; i < v->get_dim(); i++) {
                vid_inst = vid + v->get_id_inst(i);
                assert(vid_inst < _nb_vars);
                res[vid_inst] = df->eval(i);
            }
        }
        else {
            vid_inst = vid + v->get_id_inst();
            assert(vid_inst < _nb_vars);
            res[vid_inst] = df->eval();
        }
    }
}

void Model::fill_in_maps() {
    unsigned vid, vjd, expo, temp, vid_inst, vjd_inst;
    unique_id vi_unique, vj_unique, temp_unique;
    param_* vi;
    param_* vj;
    
    _obj.compute_derivatives();
    for (auto &qt_p: _obj.get_qterms()) {
        vi = qt_p.second._p->first;
        vid = vi->get_id();
        vj = qt_p.second._p->second;
        vjd = vj->get_id();
        if (vi->_is_vector) {
            for (int i = 0; i<vi->get_dim(); i++) {
                vid_inst = vid + vi->get_id_inst(i);
                vi_unique = vi->_unique_id;
                vjd_inst = vjd + vj->get_id_inst(i);
                vj_unique = vj->_unique_id;
                if (vid_inst > vjd_inst) {
                    temp = vid_inst;
                    temp_unique = vi_unique;
                    vid_inst = vjd_inst;
                    vi_unique = vj_unique;
                    vjd_inst = temp;
                    vj_unique = temp_unique;
                }
                _hess_link[make_pair<>(vi_unique,vid_inst)][make_pair<>(vj_unique,vjd_inst)].insert(make_pair<>(-1,0));//"-1" indicates a hessian link in the objective
                _obj.get_hess_link()[vid_inst].insert(vjd_inst);
            }
        }
        else {
            vid_inst = vid + vi->get_id_inst();
            vi_unique = vi->_unique_id;
            vjd_inst = vjd + vj->get_id_inst();
            vj_unique = vj->_unique_id;
            if (vid_inst > vjd_inst) {
                temp = vid_inst;
                temp_unique = vi_unique;
                vid_inst = vjd_inst;
                vi_unique = vj_unique;
                vjd_inst = temp;
                vj_unique = temp_unique;
            }
            _hess_link[make_pair<>(vi_unique,vid_inst)][make_pair<>(vj_unique,vjd_inst)].insert(make_pair<>(-1,0));//"-1" indicates a hessian link in the objective
            _obj.get_hess_link()[vid_inst].insert(vjd_inst);

        }
    }
    
    for (auto &pt_p: _obj.get_pterms()) {
        for (auto v_it = pt_p.second._l->begin(); v_it != pt_p.second._l->end(); v_it++) {
            vi = v_it->first;
            vid = vi->get_id();
            vid_inst = vid + vi->get_id_inst();
            vi_unique = vi->_unique_id;
            expo = v_it->second;
            if (expo>1) {
                _hess_link[make_pair<>(vi_unique,vid_inst)][make_pair<>(vi_unique,vid_inst)].insert(make_pair<>(-1,0));//"-1" indicates objective
                _obj.get_hess_link()[vid_inst].insert(vid_inst);
            }
            for (auto v_jt = next(v_it); v_jt != pt_p.second._l->end(); v_jt++) {
                vi = v_it->first;
                vid = vi->get_id();
                vid_inst = vid + vi->get_id_inst();
                vi_unique = vi->_unique_id;
                vj = v_jt->first;
                vjd = vj->get_id();
                vjd_inst = vjd + vj->get_id_inst();
                vj_unique = vj->_unique_id;
                if (vid_inst > vjd_inst) {
                    temp = vid_inst;
                    temp_unique = vi_unique;
                    vid_inst = vjd_inst;
                    vi_unique = vj_unique;
                    vjd_inst = temp;
                    vj_unique = temp_unique;
                }
                _hess_link[make_pair<>(vi_unique,vid_inst)][make_pair<>(vj_unique,vjd_inst)].insert(make_pair<>(-1,0));//"-1" indicates objective
                _obj.get_hess_link()[vid_inst].insert(vjd_inst);
                
            }
        }
    }
    
    Constraint* c = NULL;
    unsigned cid = 0;
    for(auto& c_p :_cons)
    {
        c = c_p.second;
        c->compute_derivatives();
        for (int inst = 0; inst<c->get_nb_instances(); inst++) {
            cid = c->_id + inst;
            for (auto &qt_p: c->get_qterms()) {
                vi = qt_p.second._p->first;
                vid = vi->get_id();
                vid_inst = vid + vi->get_id_inst(inst);
                vi_unique = vi->_unique_id;
                vj = qt_p.second._p->second;
                vjd = vj->get_id();
                vjd_inst = vjd + vj->get_id_inst(inst);
                vj_unique = vj->_unique_id;
                if (vi->_is_vector) {
                    for (int i = 0; i<vi->get_dim(); i++) {
                        vid_inst = vid + vi->get_id_inst(i);
                        vi_unique = vi->_unique_id;
                        vjd_inst = vjd + vj->get_id_inst(i);
                        vj_unique = vj->_unique_id;
                        if (vid_inst > vjd_inst) {
                            temp = vid_inst;
                            temp_unique = vi_unique;
                            vid_inst = vjd_inst;
                            vi_unique = vj_unique;
                            vjd_inst = temp;
                            vj_unique = temp_unique;
                        }
                        _hess_link[make_pair<>(vi_unique,vid_inst)][make_pair<>(vj_unique,vjd_inst)].insert(make_pair<>(c->_id,inst));
                        c->get_hess_link()[vid_inst].insert(vjd_inst);
                    }
                }
                else {
                    if (vid_inst > vjd_inst) {
                        temp = vid_inst;
                        temp_unique = vi_unique;
                        vid_inst = vjd_inst;
                        vi_unique = vj_unique;
                        vjd_inst = temp;
                        vj_unique = temp_unique;
                    }
                    _hess_link[make_pair<>(vi_unique,vid_inst)][make_pair<>(vj_unique,vjd_inst)].insert(make_pair<>(c->_id,inst));
                    c->get_hess_link()[vid_inst].insert(vjd_inst);
                }
            }
            
            for (auto &pt_p: c->get_pterms()) {
                for (auto v_it = pt_p.second._l->begin(); v_it != pt_p.second._l->end(); v_it++) {
                    vi = v_it->first;
                    vid = vi->get_id();
                    vid_inst = vid + vi->get_id_inst(inst);
                    vi_unique = vi->_unique_id;
                    expo = v_it->second;
                    if (expo>1) {
                        _hess_link[make_pair<>(vi_unique,vid_inst)][make_pair<>(vi_unique,vid_inst)].insert(make_pair<>(c->_id,inst));
                        c->get_hess_link()[vid_inst].insert(vid_inst);
                    }
                    for (auto v_jt = next(v_it); v_jt != pt_p.second._l->end(); v_jt++) {
                        vi = v_it->first;
                        vid = vi->get_id();
                        vid_inst = vid + vi->get_id_inst(inst);
                        vi_unique = vi->_unique_id;
                        vj = v_jt->first;
                        vjd = vj->get_id();
                        vjd_inst = vjd + vj->get_id_inst(inst);
                        vj_unique = vj->_unique_id;
                        if (vid_inst==vjd_inst) {
                            continue;
                        }
                        if (vid_inst > vjd_inst) {
                            temp = vid_inst;
                            temp_unique = vi_unique;
                            vid_inst = vjd_inst;
                            vi_unique = vj_unique;
                            vjd_inst = temp;
                            vj_unique = temp_unique;
                        }
                        _hess_link[make_pair<>(vi_unique,vid_inst)][make_pair<>(vj_unique,vjd_inst)].insert(make_pair<>(c->_id,inst));
                        c->get_hess_link()[vid_inst].insert(vjd_inst);
                    }
                }
            }
        }
        //MISSING NONLINEAR PART!!!
        for (auto &v_p: c->get_vars()) {
//            if(!v_p.second.first->_is_indexed){
                for (int i = 0; i < v_p.second.first->get_dim(); i++) {
//                    _v_in_cons[v_p.second.first->get_id()+i].insert(c);
//                }
//            }
//            else {
//                auto ids = v_p.second.first->get_ids();
//                for (int i = 0; i < ids.size(); i++) {
                    _v_in_cons[v_p.second.first->get_id() + v_p.second.first->get_id_inst(i)].insert(c);
//                }                
            }
        }
    }
    _nnz_h = 0;
    for (auto &hess_i: _hess_link) {
        _nnz_h += hess_i.second.size();
    }    
}


void Model::fill_in_var_init(double* x) {
    size_t vid, vid_inst;
    param_* v;
    for(auto& v_p: _vars)
    {
        v = v_p.second;
        vid = v->get_id();
        switch (v->get_intype()) {
            case float_: {
                auto real_var = (var<float>*)v;
                for (int i = 0; i < v->get_dim(); i++) {
                    vid_inst = vid + v->get_id_inst(i);
                    x[vid_inst] = (double)real_var->eval(i);
                }
                break;
            }
            case long_:{
                auto real_var = (var<long double>*)v;
                for (int i = 0; i < v->get_dim(); i++) {
                    vid_inst = vid + v->get_id_inst(i);
                    x[vid_inst] = (double)real_var->eval(i);
                }
                break;
            }
            case double_:{
                auto real_var = (var<double>*)v;
                for (int i = 0; i < v->get_dim(); i++) {
                    vid_inst = vid + v->get_id_inst(i);
                    x[vid_inst] = (double)real_var->eval(i);
                }
                break;
            }
            case integer_:{
                auto real_var = (var<int>*)v;
                for (int i = 0; i < v->get_dim(); i++) {
                    vid_inst = vid + v->get_id_inst(i);
                    x[vid_inst] = (double)real_var->eval(i);
                }
                break;
            }
            case short_:{
                auto real_var = (var<short>*)v;
                for (int i = 0; i < v->get_dim(); i++) {
                    vid_inst = vid + v->get_id_inst(i);
                    x[vid_inst] = (double)real_var->eval(i);
                }
                break;
            }
            case binary_:{
                auto real_var = (var<bool>*)v;
                for (int i = 0; i < v->get_dim(); i++) {
                    vid_inst = vid + v->get_id_inst(i);
                    x[vid_inst] = (double)real_var->eval(i);
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
        c = c_p.second;
        switch (c->get_type()) {
            case eq:{
                auto nb_ins = c->get_nb_instances();
                for (int i = 0; i< nb_ins; i++){
                    cid = c->_id+i;
                    g_l[cid] = c->_rhs;
                    g_u[cid] = c->_rhs;
                }
                break;
            }
            case leq:{
                auto nb_ins = c->get_nb_instances();
                for (int i = 0; i< nb_ins; i++){
                    cid = c->_id+i;
                    g_l[cid] = numeric_limits<double>::lowest();
                    g_u[cid] = c->_rhs;
                }
                break;
            }
            case geq:{
                auto nb_ins = c->get_nb_instances();
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

#ifdef USE_BONMIN
void Model::fill_in_param_types(Bonmin::TMINLP::VariableType* param_types){
    size_t vid;
    param_* v;
    for(auto& v_p: _vars)
    {
        v = v_p.second;
        vid = v->get_id();
        if (v->get_intype()== short_ || v->get_intype() == integer_) {
            for (int i = 0; i < v->get_dim(); i++) {
                param_types[vid + i] = Bonmin::TMINLP::INTEGER;
            }
            
        }
        else if (v->get_intype()== binary_) {
            for (int i = 0; i < v->get_dim(); i++) {
                param_types[vid + i] = Bonmin::TMINLP::BINARY;
            }
        }
        else {
            for (int i = 0; i < v->get_dim(); i++) {
                param_types[vid + i] = Bonmin::TMINLP::CONTINUOUS;
            }
        }
    }
}

#endif


void Model::embed(expr& e){
    switch (e.get_type()) {
        case uexp_c:{
            auto ue = (uexpr*)&e;
            if (ue->_son->is_function()) {
                auto f = (func_*)ue->_son;
                embed(*f);
            }
            else if(ue->_son->is_expr()){
                embed(*(expr*)ue->_son);
            }
            else if (ue->_son->is_var()){
                if (_vars.count(((param_*)ue->_son)->get_id())==0) {
                    add_var((param_*)copy(*ue->_son));
                }
            }
            break;
        }
        case bexp_c:{
            auto be = (bexpr*)&e;
            if (be->_lson->is_function()) {
                auto f = (func_*)be->_lson;
                embed(*f);
            }
            else if(be->_lson->is_expr()){
                embed(*(expr*)be->_lson);
            }
            else if (be->_lson->is_var()){
                if (_vars.count(((param_*)be->_lson)->get_id())==0) {
                    add_var((param_*)copy(*be->_lson));
                }
            }
            if (be->_rson->is_function()) {
                auto f = (func_*)be->_rson;
                embed(*f);
            }
            else if(be->_rson->is_expr()){
                embed(*(expr*)be->_rson);
            }
            else if (be->_rson->is_var()){
                if (_vars.count(((param_*)be->_rson)->get_id())==0) {
                    add_var((param_*)copy(*be->_rson));
                }
            }
            break;
        }
        default:
            break;
    }
}

void Model::embed(func_& f){
    f._embedded = true;
    param_* p = nullptr;
    param_* p1 = nullptr;
    param_* p2 = nullptr;
    for (auto &pair:f.get_lterms()) {
        p = pair.second._p;
        if (p->is_var()) {
            auto it = _vars.find(p->get_id());
            if (it==_vars.end()) {
                add_var(p);
            }
            else{
                p = it->second;
                pair.second._p = p;                
            }
        }
        else {
            auto it = _params.find(p->get_id());
            if (it==_params.end()) {
                add_param(p);
            }
            else{
                p = it->second;
                pair.second._p = p;
            }
        }
    }
    for (auto &pair:f.get_qterms()) {
        p1 = pair.second._p->first;
        p2 = pair.second._p->second;
        if (p1->is_var()) {
            auto it1 = _vars.find(p1->get_id());
            if (it1==_vars.end()) {
                add_var(p1);
            }
            else{
                p1 = it1->second;
                pair.second._p->first = p1;
            }
            auto it2 = _vars.find(p2->get_id());
            if (it2==_vars.end()) {
                add_var(p2);
            }
            else{
                p2 = it2->second;
                pair.second._p->second = p2;
            }
        }
        else {
            auto it1 = _params.find(p1->get_id());
            if (it1==_params.end()) {
                add_param(p1);
            }
            else{
                p1 = it1->second;
                pair.second._p->first = p1;
            }
            auto it2 = _params.find(p2->get_id());
            if (it2==_params.end()) {
                add_param(p2);
            }
            else{
                p2 = it2->second;
                pair.second._p->second = p2;
            }
        }
    }
    for (auto &pair:f.get_pterms()) {
        auto list = pair.second._l;
        for (auto &ppi: *list) {
            p = ppi.first;
            if (p->is_var()) {
                auto it = _vars.find(p->get_id());
                if (it==_vars.end()) {
                    add_var(p);
                }
                else{
                    p = it->second;
                    ppi.first = p;
                }
            }
            else {
                auto it = _params.find(p->get_id());
                if (it==_params.end()) {
                    add_param(p);
                }
                else{
                    p = it->second;
                    ppi.first = p;                    
                }
            }
        }
    }
    if (f.is_nonlinear()) {
        embed(f.get_expr());
    }
    auto old_vars = f.get_vars();
    for (auto &vp: old_vars) {
        auto vv = _vars_name[vp.first];
        if (vv != vp.second.first) {
            delete vp.second.first;
            f.delete_var(vp.first);
            f.add_var(vv);
        }
    }
    auto old_params = f.get_params();
    for (auto &pp: old_params) {
        auto p = _vars_name[pp.first];
        if (p != pp.second.first) {
            delete pp.second.first;
            f.delete_param(pp.first);
            f.add_param(p);
        }
    }
    
    if (f.is_nonlinear()) {
//        f.compute_derivatives();
    }
//    f.embed_derivatives();
    
}



void Model::print_functions() const{
    cout << "Number of atomic functions = " << _functions.size();
    cout << endl;
    //    for (auto& f: _functions){
    //        f->print(false);
    //        cout << endl;
    //    }
    cout << endl;
}

void Model::print_solution() const{
    
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
