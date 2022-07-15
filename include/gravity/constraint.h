    //
    //  constraint.hpp
    //  Gravity
    //
    //  Created by Hijazi, Hassan (Data61, Canberra City) on 6/5/17.
    //
    //

#ifndef constraint_hpp
#define constraint_hpp

#include <stdio.h>
#include <gravity/func.h>

namespace gravity {
class Constraint_{
    
public:
    size_t                      _id = 0;
    size_t                      _jac_cstr_idx; /* First index of the corresponding non-zero values in the Jacobian */
    ConstraintType              _ctype = leq; /**< Constraint type: leq, geq or eq */
    vector<double>              _dual ; /**< Lagrange multipliers at a KKT point */
    bool                        _relaxed = false; /**< True if this constraint is a relaxed version of an non-convex constraint, i.e. McCormick or from == to <= or >= */
    bool                        _all_active = true;
    vector<bool>                _active;
    shared_ptr<bool>            _all_lazy;
    vector<bool>                _lazy;
    bool                        _all_satisfied = true;
    vector<bool>                _violated;
    param<double>               _onCoef; /** Coefficient vector for on in on/off constraints */
    param<double>               _offCoef; /** Coefficient vector for on in on/off constraints */
    
    
    
    
    /* Accessors */
    ConstraintType get_ctype() const;
    bool is_ineq() const{
        return (_ctype==leq || _ctype==geq);
    }
    bool is_eq() const{
        return (_ctype==eq);
    }
    
    bool is_leq() const{
        return (_ctype==leq);
    }
    
    bool is_geq() const{
        return (_ctype==geq);
    }
    
    
    
};

template<typename type = double>
class Constraint: public Constraint_, public func<type>{
    
public:
    
    
    Constraint():Constraint("noname"){};
    Constraint(const string& name):Constraint(name, leq){};
    Constraint(const string& name, ConstraintType ctype){
        this->_name = name;
        this->_ctype = ctype;
        this->_is_constraint = true;
        this->_all_lazy = make_shared<bool>(false);
        this->_dim[0] = 1;
        this->_onCoef.set_name(name+"_onCoef");
        this->_offCoef.set_name(name+"_offCoef");
        
    };
    
        //        Constraint& operator<=(type rhs) {
        //            _ctype = leq;
        //            *this -= rhs;
        //            return *this;
        //        }
        //
        //        Constraint& operator==(type rhs) {
        //            _ctype = eq;
        //            *this -= rhs;
        //            return *this;
        //        }
        //
        //        Constraint& operator>=(type rhs) {
        //            _ctype = geq;
        //            *this -= rhs;
        //            return *this;
        //        }
    
    Constraint& operator<=(const param<type>& rhs) {
        _ctype = leq;
        *this -= rhs;
        return *this;
    }
    bool equal(const Constraint<type>& c) const{
        return _ctype==c._ctype && this->func<type>::operator==(c);
    }
    
    Constraint& operator==(const param<type>& rhs) {
        _ctype = eq;
        *this -= rhs;
        return *this;
    }
    
    Constraint& operator>=(const param<type>& rhs) {
        _ctype = geq;
        *this -= rhs;
        return *this;
    }
    
    Constraint& operator==(const var<type>& rhs) {
        _ctype = eq;
        *this -= rhs;
        return *this;
    }
    
    Constraint& operator>=(const var<type>& rhs) {
        _ctype = geq;
        *this -= rhs;
        return *this;
    }
    
    Constraint& operator <=(const func<type>& f){
        _ctype = leq;
        (*this) -= f;
        return *this;
    };
    Constraint& operator >=(const func<type>& f){
        _ctype = geq;
        (*this) -= f;
        return *this;
    };
    
    Constraint& operator==(const func<type>& f){
        _ctype = eq;
        (*this) -= f;
        return *this;
    }
    
    /* Modifiers */
    
    void make_lazy() {
        *_all_lazy = true;
        _lazy.resize(this->get_nb_instances(),true);
    }
    
    
    Constraint& in(const vector<Node*>& vec) {
        this->func<type>::in(vec);
        return *this;
    }
    
    Constraint& in(const indices& ids){
        if(ids.empty()){
            this->_dim[0] = 0;
            return *this;
        }
        this->func<type>::in(ids);
        this->_dim[0] = ids.size();
        return *this;
    };
    
    Constraint(const Constraint& c){
        *this = c;
    }
    
    Constraint(Constraint&& c){
        *this = move(c);
    }
    
    void deep_copy(const Constraint& c){
        this->func<type>::deep_copy(c);
        _jac_cstr_idx = c._jac_cstr_idx;
        _id = c._id;
        _ctype = c._ctype;
        _dual = c._dual;
        _all_active = c._all_active;
        _active = c._active;
        this->_all_lazy = make_shared<bool>(*c._all_lazy);
        _lazy = c._lazy;
        _all_satisfied = c._all_satisfied;
        _violated = c._violated;
        _relaxed = c._relaxed;
        this->_name = c._name;
        this->_is_constraint = true;
        _onCoef = c._onCoef.deep_copy();
        _offCoef = c._offCoef.deep_copy();
    }
    
    Constraint& operator=(const Constraint& c){
        _jac_cstr_idx = c._jac_cstr_idx;
        _id = c._id;
        _ctype = c._ctype;
        _dual = c._dual;
        _all_active = c._all_active;
        _active = c._active;
        _all_lazy = c._all_lazy;
        _lazy = c._lazy;
        _all_satisfied = c._all_satisfied;
        _violated = c._violated;
        _relaxed = c._relaxed;
        this->func<type>::operator=(c);
        this->_name = c._name;
        this->_is_constraint = true;
        _onCoef = c._onCoef.deep_copy();
        _offCoef = c._offCoef.deep_copy();
        return *this;
    }
    
    /** Assuming an outer point and an inner point, this function uses a binary line search to find an active point for the current constraint and updates the value of the variables.
     @param[in] x_start: interior point
     @param[in] nb_inst: instance number
     @param[in] ctype: ineq type
     @return True if line search successfully solved
     The function assumes that the current value stored in vars is the outer point.
     Interior and outer point classification depends on constraint type (\geq 0 or \leq 0) as input by ctype
     **/
    bool binary_line_search(const vector<double>& x_start, size_t nb_inst);
    
    Constraint& operator=(const func<type>& c){
        this->func<type>::operator=(c);
        this->_is_constraint = true;
        return *this;
    }
    
    Constraint& operator=(Constraint&& c){
        _jac_cstr_idx = c._jac_cstr_idx;
        _id = c._id;
        _ctype = c._ctype;
        _dual = c._dual;
        _all_active = c._all_active;
        _active = c._active;
        _all_lazy = c._all_lazy;
        _lazy = c._lazy;
        _all_satisfied = c._all_satisfied;
        _violated = c._violated;
        _relaxed = c._relaxed;
        this->func<type>::operator=(move(c));
        this->_name = c._name;
        this->_is_constraint = true;
        _onCoef = move(c._onCoef);
        _offCoef = move(c._offCoef);
        return *this;
    }
    
    template<typename T=type>
    Constraint<type> replace(const var<T>& v, const func<T>& f, int& tag_iter) {/**<  Replace v with function f everywhere it appears */
        Constraint<type> cpy = *this;
        int nb_inst = get_nb_instances();
        cpy.func<type>::operator=(this->func<type>::replace(v,f,tag_iter));
        if(cpy._indices && cpy._indices->size()!=nb_inst){
            cpy._name += "_projected_in("+cpy._indices->get_name()+")";
        }
        return cpy;
    }
    
    size_t get_nb_instances() const{
        size_t nb = 0;
        if (!*_all_lazy) {
            return this->get_nb_inst();
        }
        if(_lazy.size()==0){
            return 0;
        }
        for (size_t i = 0; i<this->get_nb_inst(); i++) {
            if (!_lazy[i]) {
                nb++;
            }
        }
        return nb;
    }
    
    string get_name() const{
        return this->_name;
    };
    
    size_t get_id_inst(size_t ind) const{
        if (this->_dim[0]==1) {
            return 0;
        }
        return ind;
    }
    
    bool is_redundant() {
        if((_ctype==eq && this->is_zero()) || (_ctype==leq && this->_range.second <= 0) || (_ctype==geq && this->_range.first >=0)){
            return true;
        }
        return false;
    }
    
    bool is_convex() const{
        return (this->_all_convexity==linear_ || (this->_all_convexity==convex_ &&_ctype==leq) || (this->_all_convexity==concave_ &&_ctype==geq));
    }
    
    bool is_concave() const{
        return (this->_all_convexity==linear_ || (this->_all_convexity==convex_ &&_ctype==geq) || (this->_all_convexity==concave_ &&_ctype==leq));
    }
    
    /* Output */
    void print_symbolic(){
        cout << " " << this->_name;
        string str;
        if (this->is_constant()) {
            str += " (Constant";
        }
        else if (this->is_linear()) {
            str += " (Linear";
        }
        else if (is_convex()) {
            str += " (Convex";
        }
        else if (is_concave()){
            str += " (Concave";
        }
        else {
            str += " (Unknown";
        }
        if (this->is_complex()) {
            str += " Complex) : ";
        }
        else {
            str += ") : ";
        }
        cout << str << this->func<type>::to_str();
        switch (_ctype) {
            case leq:
                cout << " <= ";
                break;
            case geq:
                cout << " >=  ";
                break;
            case eq:
                cout << " = ";
                break;
            default:
                break;
        }
        cout << 0 << ";\n";
    };
    
    void print(int prec = 30){
        string str;
        str += " " + this->_name;
        if (this->is_constant()) {
            str += " (Constant";
        }
        else if (this->is_linear()) {
            str += " (Linear";
        }
        else if (is_convex()) {
            str += " (Convex";
        }
        else if (is_concave()){
            str += " (Concave";
        }
        else {
            str += " (Unknown";
        }
        if (this->is_complex()) {
            str += " Complex) : ";
        }
        else {
            str += ") : \n";
        }
            //            auto space_size = str.size();
        auto space_size = 0;
        auto nb_inst = this->get_nb_instances();
        this->allocate_mem();
        if (this->is_matrix()) {
            auto max_cell_size = this->get_max_cell_size();
            for (size_t i = 0; i<this->_dim[0]; i++) {
                if (i>0) {
                    str.insert(str.end(), space_size, ' ');
                }
                str += "|";
                for (size_t j = 0; j<this->_dim[1]; j++) {
                    auto cell = this->to_str(i,j,prec);
                    auto cell_size = cell.size();
                    cell.insert(0, floor((max_cell_size - cell_size)/2.), ' ');
                    cell.append(ceil((max_cell_size - cell_size)/2.), ' ');
                    str += cell;
                    if(j!=this->_dim[1]-1){
                        str += " ";
                    }
                }
                str += "|\n";
            }
        }
        else {
            for (auto inst = 0; inst<nb_inst; inst++) {
                if (*_all_lazy && _lazy[inst]) {
                    continue;
                }
                this->eval(inst);
                if (inst>0) {
                    str.insert(str.end(), space_size, ' ');
                }
                if (this->_dim[0]>1) {
                    str += this->_name;
                    if (!this->_indices || this->_indices->empty()) {
                        str += "[" + to_string(inst) + "]: ";
                    }
                    else {
                        str += "[" + this->_indices->_keys->at(inst) + "]: ";
                    }
                }
                str += this->to_str(inst,prec);
                switch (_ctype) {
                    case leq:
                        str += " <= ";
                        break;
                    case geq:
                        str += " >=  ";
                        break;
                    case eq:
                        str += " = ";
                        break;
                    default:
                        break;
                }
                str += "0;\n";
            }
        }
        cout << str;
    }
    
    bool is_active(size_t inst = 0, double tol = 1e-6) const{
        return fabs(this->_val->at(inst)) < tol;
    }
    
    int get_lterm_cont_var_name(int i, int inst) const{
        int idx = 0;
        int nb_terms = this->nb_linear_terms(inst);
        if(i>=nb_terms){
            throw invalid_argument("in call to get_lterm_cont_var_name(), out of bounds index");
        }
        map<int, lterm> ordered_lterms;
        for(const auto &lt: *this->_lterms){
            int unique_id = 0;
            while(!ordered_lterms.insert({lt.second._p->get_id()+lt.second._p->get_id_inst()+unique_id,lterm(lt.second._sign, lt.second._coef,lt.second._p)}).second)
                unique_id++;
        }
        auto iter =ordered_lterms.begin();
        while (iter != ordered_lterms.end()) {
            if(!iter->second._p->_is_relaxed){
                if(iter->second._p->is_matrix_indexed()){
                    auto n = iter->second._p->get_dim(inst);
                    for (int k = 0; k<n; k++) {
                        if(idx==i)
                            return iter->second._p->get_vec_id();
                        idx++;
                    }
                }
                else{
                    if(idx==i)
                        return iter->second._p->get_vec_id();
                    idx++;
                }
            }
            iter++;
        }
        throw invalid_argument("in call to get_lterm_cont_var_name(), out of bounds index");
    }
    
    int get_lterm_cont_var_id(int i, int inst) const{
        int idx = 0;
        int nb_terms = this->nb_linear_terms(inst);
        if(i>=nb_terms){
            throw invalid_argument("in call to get_lterm_cont_var_id(), out of bounds index");
        }
        map<int, lterm> ordered_lterms;
        
        for(const auto &lt: *this->_lterms){
            int unique_id = 0;
            while(!ordered_lterms.insert({lt.second._p->get_id()+lt.second._p->get_id_inst()+unique_id,lterm(lt.second._sign, lt.second._coef,lt.second._p)}).second)
                unique_id++;
        }
        auto iter =ordered_lterms.begin();
        while (iter != ordered_lterms.end()) {
            if(!iter->second._p->_is_relaxed){
                if(iter->second._p->is_matrix_indexed()){
                    auto n = iter->second._p->get_dim(inst);
                    for (int k = 0; k<n; k++) {
                        if(idx==i)
                            return iter->second._p->_indices->_ids->at(inst).at(k);
                        idx++;
                    }
                }
                else{
                    if(idx==i)
                        return iter->second._p->get_id_inst(inst);
                    idx++;
                }
            }
            iter++;
        }
        throw invalid_argument("in get_lterm_cont_var_id(i), i out of bounds");
    }
    
    int get_qterm_cont_var_id1(int i, int inst) const{
        int idx = 0;
        int nb_terms = this->nb_quad_terms();
        if(i>=nb_terms){
            throw invalid_argument("in call to get_qterm_cont_var_id1(), out of bounds index");
        }
        auto iter =this->_qterms->begin();
        while (iter != this->_qterms->end()) {
            if(!iter->second._p->first->_is_relaxed && !iter->second._p->second->_is_relaxed){
                if(idx==i)
                    return iter->second._p->first->get_id_inst(inst);
                idx++;
            }
            iter++;
        }
        throw invalid_argument("in get_qterm_cont_var_id1(i), i out of bounds");
    }
    
    int get_qterm_cont_var_name1(int i, int inst) const{
        int idx = 0;
        int nb_terms = this->nb_quad_terms();
        if(i>=nb_terms){
            throw invalid_argument("in call to get_qterm_cont_var_name1(), out of bounds index");
        }
        auto iter =this->_qterms->begin();
        while (iter != this->_qterms->end()) {
            if(!iter->second._p->first->_is_relaxed && !iter->second._p->second->_is_relaxed){
                if(idx==i)
                    return iter->second._p->first->get_vec_id();
                idx++;
            }
            iter++;
        }
        throw invalid_argument("in get_qterm_cont_var_name1(i), i out of bounds");
    }
    
    int get_qterm_hyb_var_id1(int i, int inst) const{/* Return the continuous variable index appearing in the ith hybrid quadratic term */
        int idx = 0;
        int nb_terms = this->nb_quad_terms();
        if(i>=nb_terms){
            throw invalid_argument("in call to get_qterm_hyb_var_id1(), out of bounds index");
        }
        auto iter =this->_qterms->begin();
        while (iter != this->_qterms->end()) {
            if(iter->second._p->first->_is_relaxed != iter->second._p->second->_is_relaxed){
                if(idx==i){
                    if(!iter->second._p->first->_is_relaxed)
                        return iter->second._p->first->get_id_inst(inst);
                    return iter->second._p->second->get_id_inst(inst);
                }
                idx++;
            }
            iter++;
        }
        return -1;
    }
    
    int get_qterm_hyb_var_name1(int i, int inst) const{/* Return the continuous variable index appearing in the ith hybrid quadratic term */
        int idx = 0;
        int nb_terms = this->nb_quad_terms();
        if(i>=nb_terms){
            throw invalid_argument("in call to get_qterm_hyb_var_name1(), out of bounds index");
        }
        auto iter =this->_qterms->begin();
        while (iter != this->_qterms->end()) {
            if(iter->second._p->first->_is_relaxed != iter->second._p->second->_is_relaxed){
                if(idx==i){
                    if(!iter->second._p->first->_is_relaxed)
                        return iter->second._p->first->get_vec_id();
                    return iter->second._p->second->get_vec_id();
                }
                idx++;
            }
            iter++;
        }
        throw invalid_argument("in call to get_qterm_hyb_var_name1(), out of bounds index");
    }
    
    int get_qterm_hyb_var_id2(int i, int inst) const{/* Return the integer variable index appearing in the ith hybrid quadratic term */
        int idx = 0;
        int nb_terms = this->nb_quad_terms();
        if(i>=nb_terms){
            throw invalid_argument("in call to get_qterm_hyb_var_id2(), out of bounds index");
        }
        auto iter =this->_qterms->begin();
        while (iter != this->_qterms->end()) {
            if(iter->second._p->first->_is_relaxed != iter->second._p->second->_is_relaxed){
                if(idx==i){
                    if(iter->second._p->first->_is_relaxed)
                        return iter->second._p->first->get_id_inst(inst);
                    return iter->second._p->second->get_id_inst(inst);
                }
                idx++;
            }
            iter++;
        }
        return -1;
    }
    
    int get_qterm_hyb_var_name2(int i, int inst) const{/* Return the integer variable index appearing in the ith hybrid quadratic term */
        int idx = 0;
        int nb_terms = this->nb_quad_terms();
        if(i>=nb_terms){
            throw invalid_argument("in call to get_qterm_hyb_var_name2(), out of bounds index");
        }
        auto iter =this->_qterms->begin();
        while (iter != this->_qterms->end()) {
            if(iter->second._p->first->_is_relaxed != iter->second._p->second->_is_relaxed){
                if(idx==i){
                    if(iter->second._p->first->_is_relaxed)
                        return iter->second._p->first->get_vec_id();
                    return iter->second._p->second->get_vec_id();
                }
                idx++;
            }
            iter++;
        }
        throw invalid_argument("in call to get_qterm_hyb_var_name2(), out of bounds index");
    }
    
    int get_qterm_int_var_id1(int i, int inst) const{
        int idx = 0;
        int nb_terms = this->nb_quad_terms();
        if(i>=nb_terms){
            throw invalid_argument("in call to get_qterm_int_var_id1(), out of bounds index");
        }
        auto iter =this->_qterms->begin();
        while (iter != this->_qterms->end()) {
            if(iter->second._p->first->_is_relaxed && iter->second._p->second->_is_relaxed){
                if(idx==i)
                    return iter->second._p->first->get_id_inst(inst);
                idx++;
            }
            iter++;
        }
        return -1;
    }
    
    int get_qterm_int_var_name1(int i, int inst) const{
        int idx = 0;
        int nb_terms = this->nb_quad_terms();
        if(i>=nb_terms){
            throw invalid_argument("in call to get_qterm_int_var_name1(), out of bounds index");
        }
        auto iter =this->_qterms->begin();
        while (iter != this->_qterms->end()) {
            if(iter->second._p->first->_is_relaxed && iter->second._p->second->_is_relaxed){
                if(idx==i)
                    return iter->second._p->first->get_vec_id();
                idx++;
            }
            iter++;
        }
        throw invalid_argument("in call to get_qterm_int_var_name1(), out of bounds index");
    }
    
    int get_qterm_int_var_id2(int i, int inst) const{
        int idx = 0;
        int nb_terms = this->nb_quad_terms();
        if(i>=nb_terms){
            throw invalid_argument("in call to get_qterm_int_var_id2(), out of bounds index");
        }
        auto iter =this->_qterms->begin();
        while (iter != this->_qterms->end()) {
            if(iter->second._p->first->_is_relaxed && iter->second._p->second->_is_relaxed){
                if(idx==i)
                    return iter->second._p->second->get_id_inst(inst);
                idx++;
            }
            iter++;
        }
	return -1;
    }
    
    int get_qterm_int_var_name2(int i, int inst) const{
        int idx = 0;
        int nb_terms = this->nb_quad_terms();
        if(i>=nb_terms){
            throw invalid_argument("in call to get_qterm_int_var_name2(), out of bounds index");
        }
        auto iter =this->_qterms->begin();
        while (iter != this->_qterms->end()) {
            if(iter->second._p->first->_is_relaxed && iter->second._p->second->_is_relaxed){
                if(idx==i)
                    return iter->second._p->second->get_vec_id();
                idx++;
            }
            iter++;
        }
        throw invalid_argument("in call to get_qterm_int_var_name2(), out of bounds index");
    }
    
    int get_qterm_cont_var_id2(int i, int inst) const{
        int idx = 0;
        int nb_terms = this->nb_quad_terms();
        if(i>=nb_terms){
            throw invalid_argument("in call to get_qterm_cont_var_id2(), out of bounds index");
        }
        auto iter =this->_qterms->begin();
        while (iter != this->_qterms->end()) {
            if(!iter->second._p->first->_is_relaxed && !iter->second._p->second->_is_relaxed){
                if(idx==i)
                    return iter->second._p->second->get_id_inst(inst);
                idx++;
            }
            iter++;
        }
	return -1;
    }
    
    int get_qterm_cont_var_name2(int i, int inst) const{
        int idx = 0;
        int nb_terms = this->nb_quad_terms();
        if(i>=nb_terms){
            throw invalid_argument("in call to get_qterm_cont_var_name2(), out of bounds index");
        }
        auto iter =this->_qterms->begin();
        while (iter != this->_qterms->end()) {
            if(!iter->second._p->first->_is_relaxed && !iter->second._p->second->_is_relaxed){
                if(idx==i)
                    return iter->second._p->second->get_vec_id();
                idx++;
            }
            iter++;
        }
        throw invalid_argument("in call to get_qterm_cont_var_name2(), out of bounds index");
    }
    
    
    int get_lterm_int_var_id(int i, int inst) const{
        int idx = 0;
        int nb_terms = this->nb_linear_terms(inst);
        if(i>=nb_terms){
            throw invalid_argument("in call to get_lterm_int_var_id(), out of bounds index");
        }
        map<int, lterm> ordered_lterms;
        for(const auto &lt: *this->_lterms){
            int unique_id = 0;
            while(!ordered_lterms.insert({lt.second._p->get_id()+lt.second._p->get_id_inst()+unique_id,lterm(lt.second._sign, lt.second._coef,lt.second._p)}).second)
                unique_id++;
        }
        auto iter =ordered_lterms.begin();
        while (iter != ordered_lterms.end()) {
            if(iter->second._p->_is_relaxed){
                if(iter->second._p->is_matrix_indexed()){
                    auto n = iter->second._p->get_dim(inst);
                    for (int k = 0; k<n; k++) {
                        if(idx==i)
                            return iter->second._p->_indices->_ids->at(inst).at(k);
                        idx++;
                    }
                }
                else{
                    if(idx==i)
                        return iter->second._p->get_id_inst(inst);
                    idx++;
                }
            }
            iter++;
        }
        return -1;
    }
    
    int get_lterm_int_var_name(int i, int inst) const{
        int idx = 0;
        int nb_terms = this->nb_linear_terms(inst);
        if(i>=nb_terms){
            throw invalid_argument("in call to get_lterm_int_var_name(), out of bounds index");
        }
        map<int, lterm> ordered_lterms;
        for(const auto &lt: *this->_lterms){
            int unique_id = 0;
            while(!ordered_lterms.insert({lt.second._p->get_id()+lt.second._p->get_id_inst()+unique_id,lterm(lt.second._sign, lt.second._coef,lt.second._p)}).second)
                unique_id++;
        }
        auto iter =ordered_lterms.begin();
        while (iter != ordered_lterms.end()) {
            if(iter->second._p->_is_relaxed){
                if(iter->second._p->is_matrix_indexed()){
                    auto n = iter->second._p->get_dim(inst);
                    for (int k = 0; k<n; k++) {
                        if(idx==i)
                            return iter->second._p->get_vec_id();
                        idx++;
                    }
                }
                else{
                    if(idx==i)
                        return iter->second._p->get_vec_id();
                    idx++;
                }
            }
            iter++;
        }
        throw invalid_argument("in call to get_lterm_int_var_name(), out of bounds index");
    }
    
    type eval_qterm_cont_coef(int i, int inst) const{
        int nb_terms = this->nb_quad_terms();
        int idx = 0;
        auto iter =this->_qterms->begin();
        while (iter != this->_qterms->end()) {
            if(!iter->second._p->first->_is_relaxed && !iter->second._p->second->_is_relaxed){
                if(idx==i){
                    if(iter->second._coef->is_param()) {
                        auto p_cst = ((param<>*)(iter->second._coef.get()));
                        if (!iter->second._sign) {
                            return -1*p_cst->eval(inst);
                        }
                        else {
                            return p_cst->eval(inst);
                        }
                    }
                    if(iter->second._coef->is_function()) {
                        auto f = static_pointer_cast<func<>>(iter->second._coef);
                        if (!iter->second._sign) {
                            return -1*f->eval(inst);
                        }
                        else {
                            return f->eval(inst);
                        }
                    }
                    if(iter->second._coef->is_number()) {
                        auto p_cst = ((constant<>*)(iter->second._coef.get()));
                        if (!iter->second._sign) {
                            return -1*p_cst->eval();
                        }
                        else {
                            return p_cst->eval();
                        }
                    }
                }
                idx++;
            }
            iter++;
        }
	return 0;
    }
    
    type eval_lterm_cont_coef(int i, int inst) const{
        int nb_terms = this->nb_linear_terms(inst);
        int idx = 0;
        
        map<int, lterm> ordered_lterms;
        for(const auto &lt: *this->_lterms){
            int unique_id = 0;
            while(!ordered_lterms.insert({lt.second._p->get_id()+lt.second._p->get_id_inst()+unique_id,lterm(lt.second._sign, lt.second._coef,lt.second._p)}).second)
                unique_id++;
        }
        auto iter =ordered_lterms.begin();
        while (iter != ordered_lterms.end()) {
            if(!iter->second._p->_is_relaxed){
                if(iter->second._p->is_matrix_indexed()){
                    auto n = iter->second._p->get_dim(inst);
                    for (int k = 0; k<n; k++) {
                        if(idx==i){
                            if(iter->second._coef->is_param()) {
                                auto p_cst = ((param<>*)(iter->second._coef.get()));
                                if (!iter->second._sign) {
                                    return -1*p_cst->eval(k);
                                }
                                else {
                                    return p_cst->eval(k);
                                }
                            }
                            if(iter->second._coef->is_function()) {
                                auto f = static_pointer_cast<func<>>(iter->second._coef);
                                if (!iter->second._sign) {
                                    return -1*f->eval(k);
                                }
                                else {
                                    return f->eval(k);
                                }
                            }
                            if(iter->second._coef->is_number()) {
                                auto p_cst = ((constant<>*)(iter->second._coef.get()));
                                if (!iter->second._sign) {
                                    return -1*p_cst->eval();
                                }
                                else {
                                    return p_cst->eval();
                                }
                            }
                        }
                        idx++;
                    }
                }
                else{
                    if(idx==i){
                        if(iter->second._coef->is_param()) {
                            auto p_cst = ((param<>*)(iter->second._coef.get()));
                            if (!iter->second._sign) {
                                return -1*p_cst->eval(inst);
                            }
                            else {
                                return p_cst->eval(inst);
                            }
                        }
                        if(iter->second._coef->is_function()) {
                            auto f = static_pointer_cast<func<>>(iter->second._coef);
                            f->eval_all();
                            if (!iter->second._sign) {
                                return -1*f->eval(inst);
                            }
                            else {
                                return f->eval(inst);
                            }
                        }
                        if(iter->second._coef->is_number()) {
                            auto p_cst = ((constant<>*)(iter->second._coef.get()));
                            if (!iter->second._sign) {
                                return -1*p_cst->eval();
                            }
                            else {
                                return p_cst->eval();
                            }
                        }
                    }
                    idx++;
                }
            }
            iter++;
        }
        return 0;
    }
    
    type eval_qterm_int_coef(int i, int inst) const{
        int nb_terms = this->nb_quad_terms();
        int idx = 0;
        auto iter =this->_qterms->begin();
        while (iter != this->_qterms->end()) {
            if(iter->second._p->first->_is_relaxed && iter->second._p->second->_is_relaxed){
                if(idx==i){
                    if(iter->second._coef->is_param()) {
                        auto p_cst = ((param<>*)(iter->second._coef.get()));
                        if (!iter->second._sign) {
                            return -1*p_cst->eval(inst);
                        }
                        else {
                            return p_cst->eval(0);
                        }
                    }
                    if(iter->second._coef->is_function()) {
                        auto f = static_pointer_cast<func<>>(iter->second._coef);
                        if (!iter->second._sign) {
                            return -1*f->eval(inst);
                        }
                        else {
                            return f->eval(inst);
                        }
                    }
                    if(iter->second._coef->is_number()) {
                        auto p_cst = ((constant<>*)(iter->second._coef.get()));
                        if (!iter->second._sign) {
                            return -1*p_cst->eval();
                        }
                        else {
                            return p_cst->eval();
                        }
                    }
                }
                idx++;
            }
            iter++;
        }
	return 0;
    }
    
    type eval_qterm_hyb_coef(int i, int inst) const{
        int nb_terms = this->nb_quad_terms();
        int idx = 0;
        auto iter =this->_qterms->begin();
        while (iter != this->_qterms->end()) {
            if(iter->second._p->first->_is_relaxed != iter->second._p->second->_is_relaxed){
                if(idx==i){
                    if(iter->second._coef->is_param()) {
                        auto p_cst = ((param<>*)(iter->second._coef.get()));
                        if (!iter->second._sign) {
                            return -1*p_cst->eval(inst);
                        }
                        else {
                            return p_cst->eval(inst);
                        }
                    }
                    if(iter->second._coef->is_function()) {
                        auto f = static_pointer_cast<func<>>(iter->second._coef);
                        if (!iter->second._sign) {
                            return -1*f->eval(inst);
                        }
                        else {
                            return f->eval(inst);
                        }
                    }
                    if(iter->second._coef->is_number()) {
                        auto p_cst = ((constant<>*)(iter->second._coef.get()));
                        if (!iter->second._sign) {
                            return -1*p_cst->eval();
                        }
                        else {
                            return p_cst->eval();
                        }
                    }
                }
                idx++;
            }
            iter++;
        }
	return 0;
    }
    
    type eval_lterm_int_coef(int i, int inst) const{
        int nb_terms = this->nb_linear_terms(inst);
        int idx = 0;
        map<int, lterm> ordered_lterms;
        for(const auto &lt: *this->_lterms){
            int unique_id = 0;
            while(!ordered_lterms.insert({lt.second._p->get_id()+lt.second._p->get_id_inst()+unique_id,lterm(lt.second._sign, lt.second._coef,lt.second._p)}).second)
                unique_id++;
        }
        auto iter =ordered_lterms.begin();
        while (iter != ordered_lterms.end()) {
            if(iter->second._p->_is_relaxed){
                if(idx==i){
                    if(iter->second._coef->is_param()) {
                        auto p_cst = ((param<>*)(iter->second._coef.get()));
                        if (!iter->second._sign) {
                            return -1*p_cst->eval(inst);
                        }
                        else {
                            return p_cst->eval(inst);
                        }
                    }
                    if(iter->second._coef->is_function()) {
                        auto f = static_pointer_cast<func<>>(iter->second._coef);
                        if (!iter->second._sign) {
                            return -1*f->eval(inst);
                        }
                        else {
                            return f->eval(inst);
                        }
                    }
                    if(iter->second._coef->is_number()) {
                        auto p_cst = ((constant<>*)(iter->second._coef.get()));
                        if (!iter->second._sign) {
                            return -1*p_cst->eval();
                        }
                        else {
                            return p_cst->eval();
                        }
                    }
                }
                idx++;
            }
            iter++;
        }
	return 0;
    }
    
    /*Adds row(or new instance) in a quadratic constraint
     @param[in] con: quadratic constraint to add in current symbolic constraint
     */
    template<typename T=type>
    void add_quad_row(const shared_ptr<Constraint<type>>& con, int inst){
        if(!(this->is_quadratic() || this->is_linear() || this->is_constant())){
            throw invalid_argument("calling add_quad_row on a nonlinear constraint that is not quadratic!");
        }
        int nb_inst = this->get_nb_instances()+1;
        this->_indices->add("inst_"+to_string(nb_inst));
        this->_dim[0] = this->_indices->_keys->size();
        this->_violated.push_back(true);
        DebugOff("nb inst "<<nb_inst);
        int nb_cont_lin_terms = con->nb_cont_lterms(inst);
        int nb_int_lin_terms = con->nb_int_lterms(inst);
        int nb_cont_quad_terms = con->nb_cont_quad_terms(inst);
        int nb_int_quad_terms = con->nb_int_quad_terms(inst);
        int nb_hyb_quad_terms = con->nb_hyb_quad_terms(inst);
        if(nb_cont_lin_terms!=this->nb_cont_lterms(0) || nb_int_lin_terms!=this->nb_int_lterms(0) || nb_cont_quad_terms!= this->nb_cont_quad_terms(0) || nb_int_quad_terms != this->nb_int_quad_terms(0) || nb_hyb_quad_terms != this->nb_hyb_quad_terms(0)){
            throw invalid_argument("adding row with different sparsity structure");
        }
        map<int, lterm> ordered_lterms;
        for(const auto &lt: *this->_lterms){
            int unique_id = 0;
            while(!ordered_lterms.insert({lt.second._p->get_id()+lt.second._p->get_id_inst()+unique_id,lterm(lt.second._sign, lt.second._coef,lt.second._p)}).second)
                unique_id++;
        }
        auto iter =ordered_lterms.begin();
        for(int i = 0;i < nb_cont_lin_terms; i++){/* terms with continuous variables */
            auto l = iter++;
            if(l->second._coef->is_param()) {
                auto p_cst = ((param<>*)(l->second._coef.get()));
                p_cst->add_val("inst_"+to_string(nb_inst), con->eval_lterm_cont_coef(i,inst));
                DebugOff("added p"<<endl);
            }
            else {
                auto f = static_pointer_cast<func<>>(l->second._coef);
                if(!f->func_is_param()){
                    throw invalid_argument("function should be a param");
                }
                auto p = static_pointer_cast<param<>>(f->_params->begin()->second.first);
                p->add_val("inst_"+to_string(nb_inst), con->eval_lterm_cont_coef(i,inst));
                l->second._coef = p;
            }
            l->second._p->_indices->add_ref(con->get_lterm_cont_var_id(i,inst));
        }
        for(int i = 0;i < nb_int_lin_terms; i++){/* terms with integer variables */
            auto l =  iter++;
            if(l->second._coef->is_param()) {
                auto p_cst = ((param<>*)(l->second._coef.get()));
                p_cst->add_val("inst_"+to_string(nb_inst), con->eval_lterm_int_coef(i,inst));
                DebugOff("added p"<<endl);
            }
            else {
                auto f = static_pointer_cast<func<>>(l->second._coef);
                if(!f->func_is_param()){
                    throw invalid_argument("function should be a param");
                }
                auto p = static_pointer_cast<param<>>(f->_params->begin()->second.first);
                p->add_val("inst_"+to_string(nb_inst), con->eval_lterm_int_coef(i,inst));
                l->second._coef = p;
            }
            l->second._p->_indices->add_ref(con->get_lterm_int_var_id(i,inst));
        }
        auto qiter =this->_qterms->begin();
        for(int i = 0;i < nb_cont_quad_terms; i++){/* terms with continuous variables */
            auto l = qiter++;
            if(l->second._coef->is_param()) {
                auto p_cst = ((param<>*)(l->second._coef.get()));
                p_cst->add_val("inst_"+to_string(nb_inst), con->eval_qterm_cont_coef(i,inst));
                DebugOff("added p"<<endl);
            }
            else {
                auto f = static_pointer_cast<func<>>(l->second._coef);
                if(!f->func_is_param()){
                    throw invalid_argument("function should be a param");
                }
                auto p = static_pointer_cast<param<>>(f->_params->begin()->second.first);
                p->add_val("inst_"+to_string(nb_inst), con->eval_qterm_cont_coef(i,inst));
                l->second._coef = p;
            }
            if(l->second._p->first->_indices->size()<this->_dim[0])/* make sure the index set has not been extended already (variable appearing in another quadratic term) */
                l->second._p->first->_indices->add_ref(con->get_qterm_cont_var_id1(i,inst));
            if(l->second._p->second->_indices->size()<this->_dim[0])
                l->second._p->second->_indices->add_ref(con->get_qterm_cont_var_id2(i,inst));
        }
        for(int i = 0;i < nb_int_quad_terms; i++){/* terms with integer variables */
            auto l =  qiter++;
            if(l->second._coef->is_param()) {
                auto p_cst = ((param<>*)(l->second._coef.get()));
                p_cst->add_val("inst_"+to_string(nb_inst), con->eval_qterm_int_coef(i,inst));
                DebugOff("added p"<<endl);
            }
            else {
                auto f = static_pointer_cast<func<>>(l->second._coef);
                if(!f->func_is_param()){
                    throw invalid_argument("function should be a param");
                }
                auto p = static_pointer_cast<param<>>(f->_params->begin()->second.first);
                p->add_val("inst_"+to_string(nb_inst), con->eval_qterm_int_coef(i,inst));
                l->second._coef = p;
            }
            if(l->second._p->first->_indices->size()<this->_dim[0])/* make sure the index set has not been extended already (variable appearing in another quadratic term) */
                l->second._p->first->_indices->add_ref(con->get_qterm_int_var_id1(i,inst));
            if(l->second._p->second->_indices->size()<this->_dim[0])
                l->second._p->second->_indices->add_ref(con->get_qterm_int_var_id2(i,inst));
        }
        for(int i = 0;i < nb_hyb_quad_terms; i++){/* terms with mixed cont*integer variables */
            auto l =  qiter++;
            if(l->second._coef->is_param()) {
                auto p_cst = ((param<>*)(l->second._coef.get()));
                p_cst->add_val("inst_"+to_string(nb_inst), con->eval_qterm_int_coef(i,inst));
                DebugOff("added p"<<endl);
            }
            else {
                auto f = static_pointer_cast<func<>>(l->second._coef);
                if(!f->func_is_param()){
                    throw invalid_argument("function should be a param");
                }
                auto p = static_pointer_cast<param<>>(f->_params->begin()->second.first);
                p->add_val("inst_"+to_string(nb_inst), con->eval_qterm_int_coef(i,inst));
                l->second._coef = p;
            }
            if(l->second._p->first->_indices->size()<this->_dim[0])/* make sure the index set has not been extended already (variable appearing in another quadratic term) */
                l->second._p->first->_indices->add_ref(con->get_qterm_hyb_var_id1(i,inst));
            if(l->second._p->second->_indices->size()<this->_dim[0])
                l->second._p->second->_indices->add_ref(con->get_qterm_hyb_var_id2(i,inst));
        }
            //Set value of the constant
        if(this->_cst->is_param()){
            auto co_cst = ((param<>*)(this->_cst.get()));
            co_cst->add_val("inst_"+to_string(nb_inst), con->eval_cst(inst));
        }
        else if(this->_cst->is_function()){
            auto rhs_f = static_pointer_cast<func<>>(this->_cst);
            if(!rhs_f->func_is_param()){
                throw invalid_argument("function should be a param");
            }
            auto rhs_p = static_pointer_cast<param<>>(rhs_f->_params->begin()->second.first);
            rhs_p->add_val("inst_"+to_string(nb_inst), con->eval_cst(inst));
            this->_cst = rhs_p;
        }
    }
    
    
    /*Adds row(or new instance) of a linear constraint
     @param[in] con: linear constraint to add in current symbolic constraint
     */
    template<typename T=type>
    void add_linear_row(const shared_ptr<Constraint<type>>& con, int inst){
        if(!(this->is_linear() || this->is_constant())){
            throw invalid_argument("calling add_linear_row on a nonlinear constraint!");
        }
        int nb_inst = this->get_nb_instances()+1;
        this->_indices->add("inst_"+to_string(nb_inst));
        this->_dim[0] = this->_indices->_keys->size();
        this->_violated.push_back(true);
        DebugOff("nb inst "<<nb_inst);
        int nb_int_vars = con->nb_int_lterms(inst);
        int nb_cont_vars = con->nb_cont_lterms(inst);
        if(nb_cont_vars!=this->nb_cont_lterms(inst) || nb_int_vars!=this->nb_int_lterms(inst)){
            throw invalid_argument("adding row with different sparsity structure");
        }
        map<int, lterm> ordered_lterms;
        for(const auto &lt: *this->_lterms){
            int unique_id = 0;
            while(!ordered_lterms.insert({lt.second._p->get_id()+lt.second._p->get_id_inst()+unique_id,lterm(lt.second._sign, lt.second._coef,lt.second._p)}).second)
                unique_id++;
        }
        auto iter =ordered_lterms.begin();
        for(int i = 0;i < nb_cont_vars; i++){/* terms with continuous variables */
            auto l = iter++;
            if(l->second._coef->is_param()) {
                auto p_cst = ((param<>*)(l->second._coef.get()));
                p_cst->add_val("inst_"+to_string(nb_inst), con->eval_lterm_cont_coef(i,inst));
                DebugOff("added p"<<endl);
            }
            else {
                auto f = static_pointer_cast<func<>>(l->second._coef);
                if(!f->func_is_param()){
                    throw invalid_argument("function should be a param");
                }
                auto p = static_pointer_cast<param<>>(f->_params->begin()->second.first);
                p->add_val("inst_"+to_string(nb_inst), con->eval_lterm_cont_coef(i,inst));
                l->second._coef = p;
            }
            l->second._p->_indices->add_ref(con->get_lterm_cont_var_id(i,inst));
        }
        for(int i = 0;i < nb_int_vars; i++){/* terms with integer variables */
            auto l =  iter++;
            if(l->second._coef->is_param()) {
                auto p_cst = ((param<>*)(l->second._coef.get()));
                p_cst->add_val("inst_"+to_string(nb_inst), con->eval_lterm_int_coef(i,inst));
                DebugOff("added p"<<endl);
            }
            else {
                auto f = static_pointer_cast<func<>>(l->second._coef);
                if(!f->func_is_param()){
                    throw invalid_argument("function should be a param");
                }
                auto p = static_pointer_cast<param<>>(f->_params->begin()->second.first);
                p->add_val("inst_"+to_string(nb_inst), con->eval_lterm_int_coef(i,inst));
                l->second._coef = p;
            }
            l->second._p->_indices->add_ref(con->get_lterm_int_var_id(i,inst));
        }
            //Set value of the constant
        if(this->_cst->is_param()){
            auto co_cst = ((param<>*)(this->_cst.get()));
            co_cst->add_val("inst_"+to_string(nb_inst), con->eval_cst(inst));
        }
        else if(this->_cst->is_function()){
            auto rhs_f = static_pointer_cast<func<>>(this->_cst);
            if(!rhs_f->func_is_param()){
                throw invalid_argument("function should be a param");
            }
            auto rhs_p = static_pointer_cast<param<>>(rhs_f->_params->begin()->second.first);
            rhs_p->add_val("inst_"+to_string(nb_inst), con->eval_cst(inst));
            this->_cst = rhs_p;
        }
    }
    
    bool check_convex_region(int inst){
        bool convex_region=true;
        vector<double> xres;
        if(!this->is_convex() && !this->is_rotated_soc() && !this->check_soc()) //For the SDP determinant constraint, check if the point is feasible with respect to to the SOC constraints
        {
            xres=this->get_x(inst);
            auto soc1=std::pow(xres[0],2)+std::pow(xres[3],2)-xres[6]*xres[7];
            auto soc2=std::pow(xres[1],2)+std::pow(xres[4],2)-xres[7]*xres[8];
            auto soc3=std::pow(xres[2],2)+std::pow(xres[5],2)-xres[6]*xres[8];
            if(soc1<=0 && soc2<=0 && soc3<=0){
                convex_region=true;
            }
            else{
                convex_region=false;
            }
        }
        return convex_region;
    }
};
}
#endif /* constraint_hpp */
