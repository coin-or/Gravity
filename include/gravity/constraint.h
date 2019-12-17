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
    void replace(const var<T>& v, const func<T>& f){/**<  Replace v with function f everywhere it appears */
        func<T> new_f = *this;
        if(!v.is_indexed()){/* Replace all instances of v, no need to split the constraint */
             for (auto &lt:this->get_lterms()) {
                 auto vv = lt.second._p;
                 if(v.get_vec_id()==vv->get_vec_id()){
                     if(!vv->is_indexed()){/* No need to select a subset of indices in f */
                         func<T> coef;
                         if (lt.second._coef->is_function()) {
                             coef = *static_pointer_cast<func<T>>(lt.second._coef);
                         }
                         else if(lt.second._coef->is_param()) {
                             coef = *static_pointer_cast<param<T>>(lt.second._coef);
                         }
                         else if(lt.second._coef->is_number()) {
                             coef = *static_pointer_cast<constant<T>>(lt.second._coef);
                         }
                         if(lt.second._sign) {
                             f -= coef*v;
                             f += coef*f;
                         }
                         else {
                             f += coef*v;
                             f -= coef*f;
                         }
                     }
                 }
             }
             for (auto &lt:this->get_qterms()) {
                 if(v.get_vec_id()==lt.second._p->first->get_vec_id()){
                 }
//                     lt.second._p->second = new_vars->at(lt.second._p->second->get_name(false,false)).first;
             }
             for (auto &lt:this->get_pterms()) {
                 for (auto &v_p:*lt.second._l) {
                     if(v.get_vec_id()==v_p.first->get_vec_id()){
                         
                     }
                 }
             }
             if (this->_expr) {
                 if (this->_expr->is_uexpr()) {
                     auto ue = static_pointer_cast<uexpr<type>>(this->_expr);
                 }
                 else {
                     auto be = static_pointer_cast<bexpr<type>>(this->_expr);
                 }
             }
        }
        else {
        }
        *this = f;
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
    
    void print(int prec = 5){
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
};
}
#endif /* constraint_hpp */
