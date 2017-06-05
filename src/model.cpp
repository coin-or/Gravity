//
//  model.cpp
//  Gravity
//
//  Created by Hijazi, Hassan (Data61, Canberra City) on 6/5/17.
//
//

#include <Gravity/model.h>


using namespace std;

/** Constructor */
//@{
Model::Model(){
    _name = NULL;
    _obj = NULL;
    _nnz_g = 0;
    _nnz_h = 0;
};
//@}

/* Destructor */
Model::~Model(){
    for (auto it:_cons) {
        delete it;
    }
    delete _obj;
};


/* Accessors */


int Model::get_nb_vars() const{
    return (int)_vars.size();
};

int Model::get_nb_cons() const{
    return (int)_cons.size();
};


int Model::get_nb_nnz_g() const{
    //    return _nnz_g;
    int idx = 0;
    for(auto& c: _cons)
    {
//        idx += c->get_nb_vars();
    }
    //    cout << "Jacobian nnz = " << idx << endl;
    return idx;
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



int Model::get_nb_nnz_h() const{
    int idx=0;
    int vid = 0, vjd = 0;
    /* return the structure of the hessian */
    for(auto& v: _vars)
    {
//        vid = v->get_idx();
//        //        v->print();
//        //        cout << "hessian link :\n";
//        for(auto itv = v->_hess.cbegin(); itv != v->_hess.cend(); ++itv)
//        {
//            vjd = *itv;
//            if (vjd <= vid) { // This is a symmetric matrix, fill the lower left triangle only.
//                //                vj = this->_vars.find(vjd)->second;
//                //                vj->print();
//                //                cout << " , ";
//                idx++;
//            }
//        }
        //        cout << "\n";
    }
    //        cout << "Hessian nnz = " << idx << endl;
    //    exit(-1);
    return idx;
};


bool Model::has_var(const var_& v) const{
//    return (v->get_idx()<_vars.size());
    return false;
};
//

Constraint* Model::get_constraint(const string& cname) const{
//    return _cons[cid];
    return nullptr;
}

var_* Model::get_var(const string& vname) const{
    //    return _cons[cid];
    return nullptr;
}



/* Modifiers */

void Model::add_var(const var_& v){
//    _vars.push_back(&v);
//    v.set_idx(_idx_var++);
    //    var<float>* real_var = NULL;
    //    var<>* long_real_var = NULL;
    //    switch (v.get_type()) {
    //        case real:
    //            real_var = (var<float>*)&v;
    //            real_var->print();
    //            break;
    //        case longreal:
    //            long_real_var = (var<>*)&v;
    //            long_real_var->print();
    //            break;
    //        default:
    //            break;
    //    }
};


void Model::del_var(const var_& v){
    //    _vars.erase(v->get_idx());
    assert(false);
};

void Model::add_constraint(const Constraint& c){
    var_* vi = NULL;
    int vid = 0;
    //var<> v;
    //    Function* dfdx = NULL;
    /*   if(c_.get_ftype()==lin_){
     for(auto it = c_._vars.cbegin(); it != c_._vars.end(); it++) {
     vid = it->first;
     if (fabs(c_.get_coeff(vid)) <= 0.0001) {
     //cout << "\nIgnoring constraint " << c_._name;
     return;
     } // Ignoring constraints with very small coefficients
     }
     }*/
//    Constraint* c = new Constraint(c_);
//    _cons.push_back(c);
//    c->set_idx(_idx_con++);
//    
//    for(auto it = c->_vars.cbegin(); it != c->_vars.end(); it++) {
//        vid = it->first;
//        vi = it->second;
//        vi->addConstraint(c);
//        
//        //        if (c->get_ftype()==nlin_) {
//        //            dfdx = c->get_dfdx(vi);
//        //                    c->print();
//        //                    vi->print();
//        //                    cout << " : ";
//        //                    dfdx->print(true);
//        //            c->set_dfdx(vid, dfdx);
//        //                    if (c->_hess.count(vid)==1) {
//        //                        for (int vjd : *c->_hess[vid]){
//        ////                            if(dfdx->has_var(vjd) && dfdx->_exp && !dfdx->_exp->is_leaf() && dfdx->_exp->has_var(vjd))
//        //                            if(dfdx->has_var(vjd))
//        //                                dfdx->set_dfdx(vjd, dfdx->get_dfdx(vjd));
//        ////                                    c->print();
//        ////                                    vi->print();
//        ////                                    cout << " : ";
//        ////                                    dfdx->print(true);
//        ////                                    getVar(vjd)->print();
//        ////                                    cout << " : ";
//        ////                                    dfdx->get_dfdx(vjd)->print(true);
//        //
//        //
//        //                        }
//        //                    }
//        
//        //        }
//    }
//    _nnz_g+=c->get_nb_vars();
    
    //    c->print();
};




void Model::add_on_off(const Constraint& c, var<bool>& on){
    if (c.get_ftype() != lin_) {
        cerr << "Nonlinear constraint.\n";
        exit(-1);
    }
    Constraint res(c.get_name() + "_on/off");
    double b;
//    for(auto it: orig_q->_coefs) {
//        v = getVar_<double>(it.first);
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
//            v = getVar_<double>(it.first);
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
    Constraint MC1(name+"_McCormick1");
    MC1 += v;
    MC1 -= v1.get_lb()*v2 + v2.get_lb()*v1 - v1.get_lb()*v2.get_lb();
    MC1 >= 0;
    add_constraint(MC1);
    //    MC1.print();
    Constraint MC2(name+"_McCormick2");
    MC2 += v;
    MC2 -= v1.get_ub()*v2 + v2.get_ub()*v1 - v1.get_ub()*v2.get_ub();
    MC2 >= 0;
    add_constraint(MC2);
    //    MC2.print();
    Constraint MC3(name+"_McCormick3");
    MC3 += v;
    MC3 -= v1.get_lb()*v2 + v2.get_ub()*v1 - v1.get_lb()*v2.get_ub();
    MC3 <= 0;
    add_constraint(MC3);
    //    MC3.print();
    Constraint MC4(name+"_McCormick4");
    MC4 += v;
    MC4 -= v1.get_ub()*v2 + v2.get_lb()*v1 - v1.get_ub()*v2.get_lb();
    MC4 <= 0;
    add_constraint(MC4);
    //    MC4.print();
}


void Model::add_on_off_McCormick(std::string name, var<>& v, var<>& v1, var<>& v2, var<bool>& on) {
    Constraint MC1(name+"_McCormick1");
    MC1 += v;
    MC1 -= v1.get_lb()*v2 + v2.get_lb()*v1 - v1.get_lb()*v2.get_lb();
    MC1 >= 0;
    add_on_off(MC1, on);
    Constraint MC2(name+"_McCormick2");
    MC2 += v;
    MC2 -= v1.get_ub()*v2 + v2.get_ub()*v1 - v1.get_ub()*v2.get_ub();
    MC2 >= 0;
    add_on_off(MC2, on);
    Constraint MC3(name+"_McCormick3");
    MC3 += v;
    MC3 -= v1.get_lb()*v2 + v2.get_ub()*v1 - v1.get_lb()*v2.get_ub();
    MC3 <= 0;
    add_on_off(MC3, on);
    Constraint MC4(name+"_McCormick4");
    MC4 += v;
    MC4 -= v1.get_ub()*v2 + v2.get_lb()*v1 - v1.get_ub()*v2.get_lb();
    MC4 <= 0;
    add_on_off(MC4, on);
}


void Model::del_constraint(const Constraint& c){
    //    _cons.erase(c->get_idx());
    assert(false);
};

void Model::set_objective(const func_& f) {
    _obj = new func_(f);
}

void Model::set_objective_type(ObjectiveType t) {
    _objt = t;
}


void Model::check_feasible(const double* x){
    int vid = 0;
    //    var_* v = NULL;
    var<>* var = NULL;
    /* return the structure of the hessian */
//    for(auto& v: _vars)
//    {
//        vid = v->get_idx();
//        var = getVar_<double>(vid);
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
    var<int>* int_var = NULL;
    var<bool>* bin_var = NULL;
    var<float>* real_var = NULL;
    var<>* long_real_var = NULL;
    int idx=0;
    for(auto& v: _vars)
    {
//        switch (v->get_type()) {
//            case real:
//                real_var = (var<float>*)v;
//                //                real_var->print();
//                x_l[idx] = (double)min(real_var->get_lb(), real_var->get_lb_off());
//                x_u[idx] = (double)max(real_var->get_ub(), real_var->get_ub_off());
//                break;
//            case longreal:
//                long_real_var = (var<>*)v;
//                //                long_real_var->print();
//                x_l[idx] = min(long_real_var->get_lb(), long_real_var->get_lb_off());
//                x_u[idx] = max(long_real_var->get_ub(), long_real_var->get_ub_off());
//                break;
//            case integ:
//                int_var = (var<int>*)v;
//                x_l[idx] = (double)int_var->get_lb();
//                x_u[idx] = (double)int_var->get_ub();
//                break;
//            case binary:
//                bin_var = (var<bool>*)v;
//                x_l[idx] = bin_var->get_lb();
//                x_u[idx] = bin_var->get_ub();
//                break;
//            default:
//                break;
//        } ;
        idx++;
    }
    //    cout << "idx = " << idx << endl;
}

void Model::fill_in_obj(const double* x , double& res){
//    res = _obj->eval(x);
}

void Model::fill_in_cstr(const double* x , double* res, bool new_x){
    int idx=0;
    
    for(auto& c: _cons)
    {
//        res[idx] = c->eval(x);
        idx++;
    }
}

void Model::fill_in_jac(const double* x , double* res, bool new_x){
    int idx=0;
    int cid = 0;
    int vid = 0;
    //    Constraint* c = NULL;
    var_* v = NULL;
    func_* dfdx = NULL;
    int meta_link = -1;

    /* return the values of the jacobian of the constraints */
    for(auto& c :_cons)
    {
//        for(auto itv = c->_vars->cbegin(); itv != c->_vars->cend(); ++itv)
//        {
//            vid = itv->first;
//            v = itv->second;
//                res[idx] = c->eval_dfdx(vid, x);
            //                res[idx] = c->eval_dfdx(vid, x);
//            idx++;
//        }
    }
}


void Model::fill_in_jac_nnz(int* iRow , int* jCol){
    int idx=0;
    int cid = 0;
    int vid = 0;
    //    Constraint* c = NULL;
    var_* v = NULL;
    /* return the structure of the jacobian */
    for(auto& c :_cons)
    {
//        cid = c->get_idx();
//        assert(cid==c->get_idx());
//        for(auto itv = c->_vars.cbegin(); itv != c->_vars.cend(); ++itv)
//        {
//            vid = itv->first;
//            v = itv->second;
//            assert(vid==v->get_idx());
//            iRow[idx] = cid;
//            jCol[idx] = vid;
//            idx++;
//        }
    }
}

void Model::fill_in_hess_nnz(int* iRow , int* jCol){
    int idx=0;
    int vid = 0, vjd = 0;
    //    var_* v = NULL;
    /* return the structure of the hessian */
    for(auto& v: _vars)
    {
//        vid = v->get_idx();
//        for(auto itv = v->_hess.cbegin(); itv != v->_hess.cend(); ++itv)
//        {
//            vjd = *itv;
//            if (vjd <= vid) { // This is a symmetric matrix, fill the lower left triangle only.
//                iRow[idx] = vid;
//                jCol[idx] = vjd;
//                v->_hess_id.insert(pair<int, int>(vjd, idx));
//                getVar(vjd)->_hess_id.insert(pair<int, int>(vid, idx));
//                //                _hess_index.insert(pair<string, int>(to_string(vid)+','+to_string(vjd), idx));
//                idx++;
//            }
//        }
    }
}

#ifdef USE_IPOPT
void Model::fill_in_var_linearity(Ipopt::TNLP::LinearityType* var_types){
    int vid = 0, cid = 0;
    Constraint* c = NULL;
    bool linear = true;
    for(auto& vi: _vars)
    {
//        vid = vi->get_idx();
//        linear = true;
//        for(auto itc = vi->_cstrs.cbegin(); itc != vi->_cstrs.cend(); ++itc)
//        {
//            cid = itc->first;
//            c = itc->second;
//            if (c->get_ftype()>lin_) {
//                linear=false;
//            }
//        }
        if (linear) var_types[vid]=Ipopt::TNLP::LINEAR;
        else var_types[vid] = Ipopt::TNLP::NON_LINEAR;
    }
}


void Model::fill_in_cstr_linearity(Ipopt::TNLP::LinearityType* const_types){
    int cid = 0;
    for(auto& c :_cons)
    {
//        cid = c->get_idx();
//        if (c->get_ftype()<=lin_) {
//            const_types[cid]=Ipopt::TNLP::LINEAR;
//        }
//        else {
//            const_types[cid] = Ipopt::TNLP::NON_LINEAR;
//        }
    }
}
#endif


void Model::fill_in_hess(const double* x , double obj_factor, const double* lambda, double* res, bool new_x){
    int vid = 0, vjd = 0, cid = 0, idx = 0;
    var_* vi = NULL;
    func_* obj_dfdx = NULL;
    func_* dfdx = NULL;
    double hess = 0;
    for(auto& v: _vars)
    {
//        vid = v->get_idx();
//        
//        for(auto itv = v->_hess.cbegin(); itv != v->_hess.cend(); ++itv)
//        {
//            vjd = *itv;
//            if (vjd <= vid) { // This is a symmetric matrix, fill the lower left triangle only.
//                res[idx] = 0;
//                idx++;
//            }
//        }
    }
    
//    for (auto &it1:_obj->_vars){
//        vid = it1.first;
//        vi = it1.second;
//        if (_obj->get_type()==nlin_) {
//            //            obj_dfdx = _obj->get_dfdx(vi);
//            //            if (obj_dfdx && !_obj->has_dfdx(vid)) {
//            //                _obj->set_dfdx(vid, obj_dfdx);
//            //            }
//            //        }
//            _obj->compute_dfdx(vi);
//        }
//        for (auto &it2:_obj->_vars){
//            vjd = it2.first;
//            if (vjd<=vid) {
//                if(_obj->_hess.count(vid)==1 && _obj->_hess.count(vjd)==1) {
//                    if(_obj->get_type()==nlin_ && obj_dfdx){
//                        if(obj_dfdx->has_var(vjd))
//                            hess = obj_dfdx->eval_dfdx(vjd, x);
//                        else
//                            hess = 0;
//                    }
//                    else{
//                        hess = _obj->get_q_coeff(vid, vjd);
//                        if (vid==vjd)
//                            hess *= 2;
//                    }
//                    if(hess != 0)
//                        res[get_hess_index(vi, vjd)] = obj_factor * hess;
//                }
//                else
//                    res[get_hess_index(vi, vjd)] = 0;
//                //                else
//                //                    res[get_hess_index(vid, vjd)] = 0;
//            }
//        }
//    }

    for (auto& c:_cons) {
        
//        for (auto it1:c->_vars){
//            vid = it1.first;
//            if(c->_hess.count(vid)==0)
//                continue;
//            vi = it1.second;
//            for (auto& it2:*c->_hess[vid]){
//                vjd = it2;
//                hess = c->get_q_coeff(vid, vjd);
//                
//                if (vid==vjd) {
//                        hess *= 2;
//                }
//                if(hess != 0)
//                    res[get_hess_index(vi, vjd)] += lambda[cid] * hess;
//            }
//        }
    }
}






void Model::fill_in_grad_obj(const double* x , double* res){
    int idx=0;
    //    var_* v = NULL;
    int vid = 0;
    for(auto& vi: _vars)
    {
//        vid = vi->get_idx();
//        
//        if (!_obj->has_var(vid)) {
//            res[idx] = 0;
//            idx++;
//            continue;
//        }
//        res[idx] = _obj->eval_dfdx(vid, x);
        idx++;
    }
}


void Model::fill_in_var_init(double* x) {
    var<int>* int_var = NULL;
    var<bool>* bin_var = NULL;
    var<float>* real_var = NULL;
    var<>* long_real_var = NULL;
    int idx=0;
    int vid = 0;
    for(auto& vi: _vars)
    {
//        vid = vi->get_idx();
//        
//        switch (vi->get_type()) {
//            case real:
//                real_var = (var<float>*)vi;
//                x[idx] = (double)real_var->get_value();
//                break;
//            case longreal:
//                long_real_var = (var<>*)vi;
//                x[idx] = long_real_var->get_value();
//                break;
//            case integ:
//                int_var = (var<int>*)vi;
//                x[idx] = (double)int_var->get_value();
//                break;
//            case binary:
//                bin_var = (var<bool>*)vi;
//                x[idx] = (double)bin_var->get_value();
//                break;
//            default:
//                exit(-1);
//                break;
//        } ;
        idx++;
    }
}

void Model::fill_in_cstr_bounds(double* g_l ,double* g_u) {
    int idx=0;
    int cid = 0;
    //    Constraint* c = NULL;
    for(auto& c : _cons)
    {
//        cid = c->get_idx();
//        switch (c->get_type()) {
//            case eq:
//                g_l[idx] = c->_rhs;
//                g_u[idx] = c->_rhs;
//                //                g_l[idx] = 0;
//                //                g_u[idx] = 0;
//                break;
//            case leq:
//                g_l[idx] = -INFINITY;
//                g_u[idx] = c->_rhs;
//                //                g_u[idx] = -c->get_const();
//                //                cout << "gl[idx] = " <<g_l[idx];
//                //                cout << ", gu[idx] = " <<g_u[idx] << endl;
//                break;
//            case geq:
//                g_l[idx] = c->_rhs;
//                //                g_l[idx] = -c->get_const();
//                g_u[idx] = INFINITY;
//                break;                
//                
//            default:
//                exit(-1);
//                break;
//        }
        idx++;
    }
}

#ifdef USE_BONMIN
void Model::fill_in_var_types(Bonmin::TMINLP::VariableType* var_types){
    int vid = 0;
    for(auto& vi: _vars){
//        vid = vi->get_idx();
//        switch (vi->get_type()) {
//            case real:
//                var_types[vid] = Bonmin::TMINLP::CONTINUOUS;
//                break;
//            case longreal:
//                var_types[vid] = Bonmin::TMINLP::CONTINUOUS;
//                break;
//            case integ:
//                var_types[vid] = Bonmin::TMINLP::INTEGER;
//                break;
//            case binary:
//                var_types[vid] = Bonmin::TMINLP::BINARY;
//                break;
//            default:
//                break;
//        } ;
    }
}

#endif

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
    var<int>* int_var = NULL;
    var<bool>* bin_var = NULL;
    var<float>* real_var = NULL;
    var<double>* long_real_var = NULL;
    int idx=0;
    int vid = 0;
    for(auto& v: _vars ) {
//        vid = v->get_idx();
//        //        assert(vid==idx);
//        switch (v->get_type()) {
//            case real:
//                real_var = (var<float>*)v;
//                real_var->print();
//                break;
//            case longreal:
//                long_real_var = (var<double>*)v;
//                long_real_var->print();
//                break;
//            case integ:
//                int_var = (var<int>*)v;
//                int_var->print();
//                break;
//            case binary:
//                bin_var = (var<bool>*)v;
//                bin_var->print();
//                break;
//            default:
//                break;
//        } ;
        idx++;
    }
}

void Model::print_constraints() const{
    for(auto& c: _cons){
        c->print();
    }
}
