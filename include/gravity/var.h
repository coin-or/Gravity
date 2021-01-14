//
//  var.h
//  Gravity
//
//  Created by Hijazi, Hassan on 21/10/16.
//
//

#ifndef var_h
#define var_h

#include <gravity/param.h>
#include <gravity/Net.h>
#include <stdio.h>
#include <string>
#include <set>
#include <list>
#include <limits>
#include <random>



using namespace std;

namespace gravity {
    
    template<typename type>
    class func;
    
    
    
    /** A variable can be a bool, a short, an int, a float, a double, a long double or a complex<double>*/
    template<typename type = double>
    /** define variable as a parameter with bounds */
    class var: public param<type>{
        
    public:
        shared_ptr<func<type>>   _lb; /**< Lower Bound */
        shared_ptr<func<type>>   _ub; /**< Upper Bound */
        bool _in_q_cone = false; /**< Used by Mosek */
        bool _psd = false; /**< Has to be positive semidefinite */
        bool _lift = false;/*flag to show if variable is a lifted variable*/
        bool _lift_lb = false;/*flag to show if the lifted variables need lower-bounding function*/
        bool _lift_ub = false;/*flag to show if the lifted variables need upper-bounding function*/
        bool _in_SOC_partn = false;/*flag to show if variable appers in a SOC partition*/
        vector<shared_ptr<var>> _original_vars;/*< If this is a lifted variable, pointers to the corresponding original variables */
        /*These should eventually be shared_ptr<int>, or an object with an access to get_id_inst, or eval */
        shared_ptr<int> _num_partns;/*number of partitons*/
        int _cur_partn = 1;/*current partition we are focused on*/
        
        
        /* Constructors */
        //@{
        /** Unbounded variable constructor */
        var(){};
        ~var() {};
        template<typename T=type,
        typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr>
        var(const string& name){
            constant_::set_type(var_c);
            this->_name = name;
            _num_partns = make_shared<int>(1);
            _lb = make_shared<func<type>>(constant<type>(numeric_limits<type>::lowest()));
            _ub = make_shared<func<type>>(constant<type>(numeric_limits<type>::max()));
            this->_range->first = _lb->_range->first;
            this->_range->second = _ub->_range->second;
        }
        
        template<class T=type, class = typename enable_if<is_same<T, Cpx>::value>::type>
        var(const string& name){
            constant_::set_type(var_c);
            this->_name = name;
            _num_partns = make_shared<int>(1);
            _lb = make_shared<func<type>>(constant<type>(Cpx(numeric_limits<double>::lowest(), numeric_limits<double>::lowest())));
            _ub = make_shared<func<type>>(constant<type>(Cpx(numeric_limits<double>::max(), numeric_limits<double>::max())));
            this->_range->first = _lb->_range->first;
            this->_range->second = _ub->_range->second;
        }
        
        var(const string& name, Sign s):var(name){
            if (s==non_neg_ || s==pos_) {
                    add_lb_only(zero<type>().eval());
            }
            else if (s==non_pos_ || s==neg_) {
                    add_ub_only(zero<type>().eval());
            }
        };
        
        var(const var<type>& v);
        var(var<type>&& v);
        
        var deep_copy() const;
        //@}
        
        shared_ptr<param_> pcopy() const{return make_shared<var>(*this);};
        
        shared_ptr<constant_> copy() const{return make_shared<var>(*this);};
        
        shared_ptr<param_> ptr_deep_copy() const{return make_shared<var>(this->deep_copy());};
        //@{
                
//        var(const string& name, type lb, type ub){
//            this->_name = name;
//            constant_::set_type(var_c);
//            _lb = make_shared<func<type>>(constant<type>(lb));
//            _ub = make_shared<func<type>>(constant<type>(ub));
//            param<type>::_range->first = lb;
//            param<type>::_range->second = ub;
//        };
        
        var(const string& name, const func<type>& lb, const func<type>& ub){
            this->_name = name;
            constant_::set_type(var_c);
            _num_partns = make_shared<int>(1);
            _lb = make_shared<func<type>>(lb);
            _ub = make_shared<func<type>>(ub);
            if(_lb->get_dim()==0 ||_ub->get_dim()==0)
            {
                this->_range->first = 0;
                this->_range->second = 0;
            }
            else {
                this->_range->first = _lb->_range->first;
                this->_range->second = _ub->_range->second;
            }

        };
        
        var(const string& name, func<type>&& lb, func<type>&& ub){
            this->_name = name;
            constant_::set_type(var_c);
            _num_partns = make_shared<int>(1);
            _lb = make_shared<func<type>>(move(lb));
            _ub = make_shared<func<type>>(move(ub));
            if(_lb->get_dim()==0 ||_ub->get_dim()==0)
            {
                this->_range->first = 0;
                this->_range->second = 0;
            }
            else {
                this->_range->first = _lb->_range->first;
                this->_range->second = _ub->_range->second;
            }
        };
        
        //@}
        
        
        Sign get_all_sign() const{
            return get_all_sign_();
        }
        
        template<class T=type, class = typename enable_if<is_same<T, Cpx>::value>::type> Sign get_all_sign_() const{
            if (_lb->is_zero() && _ub->is_zero()) {
                return zero_;
            }
            if ((_ub->_range->second.real() < 0 && _ub->_range->second.imag() < 0)) {
                return neg_;
            }
            if ((_lb->_range->first.real() > 0 && _lb->_range->first.imag() > 0)) {
                return pos_;
            }
            if (_lb->is_zero()) {
                return non_pos_;
            }
            if (_ub->is_zero()) {
                return non_neg_;
            }
            return unknown_;
        }
        
        template<typename T=type,
        typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr> Sign get_all_sign_() const {
            if (_lb->is_zero() && _ub->is_zero()) {
                return zero_;
            }
            if (_ub->_range->second < 0) {
                return neg_;
            }
            if (_lb->_range->first > 0) {
                return pos_;
            }
            if (_ub->is_zero()) {
                return non_pos_;
            }
            if (_lb->is_zero()) {
                return non_neg_;
            }
            return unknown_;
        }
        
        Sign get_sign(size_t idx = 0) const{ /**< returns the sign of one instance of the current parameter/variable. **/
            return param<type>::get_sign(idx);
        }
        
        void initialize_all(type v) {
            this->set_val(v);
        }
        
        void set_val(type val) {
            if(this->is_indexed()){
                for(auto &idx: this->_indices->_ids->at(0)){
                    this->_val->at(idx) = val;
                }
            }
            else {
                for (auto i = 0; i<this->_val->size() ;i++) {
                    this->_val->at(i) = val;
                }
            }
        }
        
        /* Retrieve specified indexed variable. */
        template<typename... Args>
        var operator()(size_t i, size_t j) {
            bool indexed = param<type>::_indices!=nullptr;
            var<type> res(*this);
            res.param<type>::operator=(param<type>::operator()(i, j));
            res._type = var_c;
            if(!indexed && !res._lb->is_number()){
                (res._lb->in(*res._indices));
            }
            if(!indexed && !res._ub->is_number()){
                (res._ub->in(*res._indices));
            }
            res._range = make_shared<pair<type,type>>(res._lb->_range->first,res._ub->_range->second);
            return res;
        }
        
        
        
        template<typename... Args>
        var operator()(string t1, Args&&... args) {
            bool indexed = param<type>::_indices!=nullptr;
            var<type> res(*this);
            res.param<type>::operator=(param<type>::operator()(t1, args...));
            res.set_type(var_c);
            if(!indexed & !res._lb->is_number()){
                (res._lb->in(*res._indices));
            }
            if(!indexed & !res._ub->is_number()){
                (res._ub->in(*res._indices));
            }
            res._range = make_shared<pair<type,type>>(res._lb->_range->first,res._ub->_range->second);
            return res;
        }
        
        
        
        
        
        template<typename... Args>
        var operator()(size_t i) {
            return (*this)(to_string(i));
        }
        
        
        template<typename... Args>
        var operator[](size_t i) {
            return (*this)(to_string(i));
        }
        
//        var in_pairs(){
//            var<type> res(*this);
//            res._name += ".in_pairs";
//            res._indices->_type = in_pairs_;
//            return res;
//        }
//
        var from(){
            var<type> res(*this);
            res._name += ".from";
            res._indices->_type = from_;
            return res;
        }
        
        var from(const indices& ids){
            return this->from().in(ids);
        }
        
        
        var to(){
            var<type> res(*this);
            res._name += ".to";
            res._indices->_type = to_;
            return res;
        }
        
        var to(const indices& ids){
            return this->to().in(ids);
        }
//        var out_arcs(){
//            var<type> res(*this);
//            res._name += ".out_arcs";
//            res._indices->_type = out_arcs_;
//            return res;
//        }
//        
//        var in_arcs(){
//            var<type> res(*this);
//            res._name += ".in_arcs";
//            res._indices->_type = in_arcs_;
//            return res;
//        }
//        
//        var in_gens(){
//            var<type> res(*this);
//            res._name += ".in_gens";
//            res._indices->_type = in_gens_;
//            return res;
//        }
        
        
        
        var in_pairs(const indices& ids) {
            bool indexed = param<type>::_indices!=nullptr;
            var<type> res(*this);
            res.param<type>::operator=(param<type>::in_pairs(ids));
            if(!indexed && !res._lb->is_number()){
                (res._lb->in(*res._indices));
            }
            if(!indexed && !res._ub->is_number()){
                (res._ub->in(*res._indices));
            }
            res._range = make_shared<pair<type,type>>(res._lb->_range->first,res._ub->_range->second);
            return res;
        }
        
        template<class T=type, class = typename enable_if<is_same<T, Cpx>::value>::type>
        void real_imag(const var<>& pr, const var<>& pi){
            this->_real = make_shared<var<>>(pr);
            this->_imag = make_shared<var<>>(pi);
            this->_range->first.real(pr._range->first);
            this->_range->first.imag(pi._range->first);
            this->_range->second.real(pr._range->second);
            this->_range->second.imag(pi._range->second);
        }
        
        
        void mag_ang(const var<>& pmag, const var<>& pang){
            this->_mag = make_shared<var<>>(pmag);
            this->_ang = make_shared<var<>>(pang);
            this->_polar = true;
        }
        
        void set_real(const var<>& p){
            this->_real = make_shared<var<>>(p);
        }
        
        void set_imag(const var<>& p){
            this->_imag = make_shared<var<>>(p);
        }

//        void reset_range(){
//            param<type>::reset_range();
//        }

//        void index_in(const indices& ids) {
//            *this = this->in(ids);
//        }
        
        template<typename... Args>
        var in(const indices& vec1, Args&&... args) {
//            bool indexed = param<type>::_indices!=nullptr;
            var<type> res(*this);
            bool first_time_indexing = this->_indices==nullptr;
            res.param<type>::operator=(param<type>::in(vec1, forward<Args>(args)...));//TODO assert lb dim = res dim
            res._type = var_c;
            if(first_time_indexing){/* Create a parameter for each bound and assign it the bounds values. This is important for bound tightening algorithms. */
                _lb->allocate_mem();
                _ub->allocate_mem();
                _lb->eval_all();
                _ub->eval_all();
                if(!_lift){
                    param<type> lb(this->get_name(true,true)+"-lb");
                    lb.index_in(*res._indices);
                    lb.copy_vals(*_lb);
                    *this->_lb = lb;
                    param<type> ub(this->get_name(true,true)+"-ub");
                    ub.index_in(*res._indices);
                    ub.copy_vals(*_ub);
                    *this->_ub = ub;
                    _lb->_val = lb._val;
                    _ub->_val = ub._val;
                    _lb->_evaluated = true;
                    _ub->_evaluated = true;
                }
                else {
                    _lb->index_in(*res._indices);
                    _ub->index_in(*res._indices);
                    _lb->allocate_mem();
                    _ub->allocate_mem();
                }
                res._lb = _lb;
                res._ub = _ub;
            }
            else if(!res._lb->func_is_number() && !res._ub->func_is_number()){
                res._lb->allocate_mem();
                res._ub->allocate_mem();
                auto new_lb(*res._lb);
                auto new_ub(*res._ub);
                new_lb.in(*res._indices);
                new_ub.in(*res._indices);
                res._range = make_shared<pair<type,type>>(new_lb._range->first,new_ub._range->second);
            }
            if(res._real){
                auto real_var = static_pointer_cast<var<>>(res._real);
                res._real = make_shared<var<>>(real_var->in(*res._indices));
            }
            if(res._imag){
                auto imag_var = static_pointer_cast<var<>>(res._imag);
                res._imag = make_shared<var<>>(imag_var->in(*res._indices));
            }
            return res;
        }
        
        /** Index variable in ids, look for the keys starting at the ith position
         @param[in] start_position If ids has keys with additional entries, use the substring starting after the start_position comma separator
         @param[in] ids_ index set
         */
        var from_ith(unsigned start_position, const indices& ids) {
            var<type> res(*this);
            res.param<type>::operator=(param<type>::from_ith(start_position, ids));//TODO assert lb dim = res dim
            res._type = var_c;
            if(res._real){
                auto real_var = static_pointer_cast<var<>>(res._real);
                res._real = make_shared<var<>>(real_var->in(*res._indices));
            }
            if(res._imag){
                auto imag_var = static_pointer_cast<var<>>(res._imag);
                res._imag = make_shared<var<>>(imag_var->in(*res._indices));
            }
            return res;
        }
        
        /** Index parameter/variable in ids, remove keys starting at the ith position and spanning nb_entries
         @param[in] start_position
         @param[in] ids_ index set
         */        
        var in_ignore_ith(unsigned start_position, unsigned nb_entries, const indices& ids_) {
            var<type> res(*this);
            res.param<type>::operator=(param<type>::in_ignore_ith(start_position, nb_entries, ids_));//TODO assert lb dim = res dim
            res._type = var_c;
            if(res._real){
                auto real_var = static_pointer_cast<var<>>(res._real);
                res._real = make_shared<var<>>(real_var->in(*res._indices));
            }
            if(res._imag){
                auto imag_var = static_pointer_cast<var<>>(res._imag);
                res._imag = make_shared<var<>>(imag_var->in(*res._indices));
            }
            return res;
        }
        
//        template<typename... Args>
//        var in_matrix(unsigned start_pos) const{
//            var<type> res(*this);
//            res.param<type>::operator=(param<type>::in_matrix(start_pos));//TODO assert lb dim = res dim
//            res._type = var_c;
//            if(res._real){
//                auto real_var = static_pointer_cast<var<>>(res._real);
//                res._real = make_shared<var<>>(real_var->in(*res._indices));
//            }
//            if(res._imag){
//                auto imag_var = static_pointer_cast<var<>>(res._imag);
//                res._imag = make_shared<var<>>(imag_var->in(*res._indices));
//            }
//            res._range = make_shared<pair<type,type>>(res._lb->_range->first,res._ub->_range->second);
//            return res;
//        }
        
        /** Create a matrix version of variable where each row will be indexed based on the entries starting at start_pos and spanning nb_entries.
         Example:
         dv = {
         [1,8] = 0
         [1,9] = 0
         [1,10] = 0
         [1,11] = 0
         [1,12] = 0
         [2,8] = 0
         [2,9] = 0
         [2,10] = 0
         [2,11] = 0
         [2,12] = 0
         [3,8] = 0
         [3,9] = 0
         [3,10] = 0
         [3,11] = 0
         [3,12] = 0
         };
         sum(dv.in_matrix(0,1)) <= 0 gives: dv[1,8] + dv[2,8] + dv[3,8] <= 0;
         sum(dv.in_matrix(1,1)) <= 0 gives: dv[1,8] + dv[1,9] + dv[1,10] + dv[1,11] + dv[1,12] <= 0;
         */
        template<typename... Args>
        var in_matrix(unsigned start_pos, unsigned nb_entries) const{
            var<type> res(*this);
            res.param<type>::operator=(param<type>::in_matrix(start_pos, nb_entries));//TODO assert lb dim = res dim
            res._type = var_c;
            if(res._real){
                auto real_var = static_pointer_cast<var<>>(res._real);
                res._real = make_shared<var<>>(real_var->in(*res._indices));
            }
            if(res._imag){
                auto imag_var = static_pointer_cast<var<>>(res._imag);
                res._imag = make_shared<var<>>(imag_var->in(*res._indices));
            }
            res._range = make_shared<pair<type,type>>(res._lb->_range->first,res._ub->_range->second);
            return res;
        }
        
        void reset_bounds(){
//            _lb->_evaluated = false;
//            _ub->_evaluated = false;
            _lb->uneval();
            _ub->uneval();
            _lb->eval_all();
            _ub->eval_all();
        }
        
        //    var in(const node_pairs& np){
        //        return this->in(np._keys);
        //    }
        
        var in_arcs(const vector<Node*>& vec) {
            var<type> res(*this);
            res.param<type>::operator=(param<type>::in_arcs(vec));
            res._type = var_c;
            return res;
        }
        
        var out_arcs(const vector<Node*>& vec) {
            var<type> res(*this);
            res.param<type>::operator=(param<type>::out_arcs(vec));
            res._type = var_c;
            return res;
        }
        
        
        var in_aux(const vector<Node*>& vec, const string& aux_type) {
            var<type> res(*this);
            res.param<type>::operator=(param<type>::in_aux(vec,aux_type));
            res._type = var_c;
            return res;
        }
        
        
        var excl(size_t index) {
            var<type> res(*this);
            res._indices->_type = excl_;
            return res;
        }
        
        
        /* Create a vector of variables indexed based on nodes from bags of size bag_size
         e.g., given bag { 1, 5, 7 } index the first variable (1), the second (5) and the last (7)
         */
        vector<var> in_bags(const vector<pair<string,vector<Node*>>>& bags, size_t bag_size);
        
        /* Create a vector of variables indexed as pair of nodes from bags of size bag_size
         e.g., given bag { 1, 5, 7 } index the first variable (1,5), the second (5,7) and the last (1,7)
         */
        vector<var<type>> pairs_in_bags(const vector<pair<string,vector<Node*>>>& bags, size_t bag_size);
        
        /* Querries */
        
        type    get_lb(size_t i) const;
        
        type    get_ub(size_t i) const;
        
        type    get_lb(const string& key) const;
        
        type    get_ub(const string& key) const;
        
        param<type>    get_lb() const;
        param<type>    get_ub() const;
        
        /* If this is a lifted variable lifted(xy)= xy, return the lowerbound on x*/
        shared_ptr<param<type>>    get_bilinear_lb1() const;
        /* If this is a lifted variable lifted(xy)= xy, return the lowerbound on y*/
        shared_ptr<param<type>>    get_bilinear_lb2() const;
        /* If this is a lifted variable lifted(xy)= xy, return the upperbound on x*/
        shared_ptr<param<type>>    get_bilinear_ub1() const;
        /* If this is a lifted variable lifted(xy)= xy, return the upperbound on y*/
        shared_ptr<param<type>>    get_bilinear_ub2() const;
        /* If this is a lifted variable lifted(x^2)= x^2, return the lowerbound on x*/
        shared_ptr<param<type>>    get_square_lb() const;
        /* If this is a lifted variable lifted(x^2)= x^2, return the upperbound on x*/
        shared_ptr<param<type>>    get_square_ub() const;
        
        
        
        template<typename T=type,
        typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr>
        bool is_bounded_above(size_t i = 0) const{
            return (_ub->eval(i)!=numeric_limits<type>::max());
        };
        
        template<typename T=type,
        typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr>
        bool is_bounded_below(size_t i = 0) const{
            return (_lb->eval(i)!=numeric_limits<type>::lowest());
        };
        
        template<typename T=type, typename enable_if<is_same<T, Cpx>::value>::type* = nullptr>
        bool is_bounded_below(size_t i = 0) const{
            return (_lb->eval(i).real()!=numeric_limits<type>::lowest() && _lb->eval(i).imag()!=numeric_limits<type>::lowest());
        };
        
        template<typename T=type, typename enable_if<is_same<T, Cpx>::value>::type* = nullptr>
        bool is_bounded_above(size_t i = 0) const{
            return (_ub->eval(i).real()!=numeric_limits<type>::max() && _ub->eval(i).imag()!=numeric_limits<type>::max());
        };
        
        template<typename T=type,
        typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr>
        bool is_constant(size_t i=0) const;
        
        template<typename T=type,
        typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr>
        Sign get_sign(size_t idx = 0) const;
        
        /* Modifiers */
        void    set_size(vector<size_t>);
        void    set_size(size_t s);
        
        void    add_bounds(type lb, type ub);
        void    add_lb_only(type v); /**< Adds a new lower bound and infinity in the corresponding upper bound*/
        void    add_ub_only(type v); /**< Adds a new upper bound and -infinity in the corresponding lower bound*/
        
        
        void    set_lb(type v);
        void    set_ub(type v);
        
        void    set_lb(const string& key, type v);
        void    set_ub(const string& key, type v);
        
        void in_q_cone(); /**< States that the variable is contained in a quadratic cone, for mosek */
        
        
        /* Operators */
        var& operator=(const var& v);
        var& operator=(var&& v);
        
        var& operator=(type v) {
            this->set_val(v);
            return *this;
        }
        
        /** let this share the bounds of p */
        void share_bounds(shared_ptr<param_> p){
            switch (p->get_intype()) {
                case binary_:{
                    auto pp =  static_pointer_cast<var<bool>>(p);
                    share_bounds_(*pp);
                }
                    break;
                case short_:{
                    auto pp =  static_pointer_cast<var<short>>(p);
                    share_bounds_(*pp);
                }
                    break;
                case integer_:{
                    auto pp =  static_pointer_cast<var<int>>(p);
                    share_bounds_(*pp);
                }
                    break;
                case float_:{
                    auto pp =  static_pointer_cast<var<float>>(p);
                    share_bounds_(*pp);
                }
                    break;
                case double_:{
                    auto pp =  (var<double>*)(p.get());
                    share_bounds_(*pp);
                }
                    break;
                case long_:{
                    auto pp =  static_pointer_cast<var<long double>>(p);
                    share_bounds_(*pp);
                }
                    break;
                case complex_:{
                    auto pp =  static_pointer_cast<var<Cpx>>(p);
                    share_bounds_(*pp);
                }
                    break;
                default:
                    break;
            }
        }
        
        /** let this share the values of p */
        void share_vals(shared_ptr<param_> p){
            switch (p->get_intype()) {
                case binary_:{
                    auto pp =  static_pointer_cast<var<bool>>(p);
                    share_vals_(*pp);
                }
                    break;
                case short_:{
                    auto pp =  static_pointer_cast<var<short>>(p);
                    share_vals_(*pp);
                }
                    break;
                case integer_:{
                    auto pp =  static_pointer_cast<var<int>>(p);
                    share_vals_(*pp);
                }
                    break;
                case float_:{
                    auto pp =  static_pointer_cast<var<float>>(p);
                    share_vals_(*pp);
                }
                    break;
                case double_:{
                    auto pp =  (var<double>*)(p.get());
                    share_vals_(*pp);
                }
                    break;
                case long_:{
                    auto pp =  static_pointer_cast<var<long double>>(p);
                    share_vals_(*pp);
                }
                    break;
                case complex_:{
                    auto pp =  static_pointer_cast<var<Cpx>>(p);
                    share_vals_(*pp);
                }
                    break;
                default:
                    break;
            }
        }
        
        
        /** let this share the bounds of p */
        template<class T2, typename std::enable_if<!is_same<T2, type>::value>::type* = nullptr>
        void share_bounds_(var<T2>& p){
            throw invalid_argument("cannot share bounds with different typed params/vars");
        }
        
        /** let this share the bounds of p */
        template<class T2, typename std::enable_if<is_same<T2, type>::value>::type* = nullptr>
        void share_bounds_(var<T2>& pp){
            this->_lb->_val = pp._lb->_val;
            this->_ub->_val = pp._ub->_val;
        }
        
        /** let this share the values of p */
        template<class T2, typename std::enable_if<!is_same<T2, type>::value>::type* = nullptr>
        void share_vals_(var<T2>& p){
            throw invalid_argument("cannot share vals with different typed params/vars");
        }
        
        /** let this share the values of p */
        template<class T2, typename std::enable_if<is_same<T2, type>::value>::type* = nullptr>
        void share_vals_(var<T2>& pp){
            this->_val = pp._val;
        }
        
        /**
         \brief Update dimensions based on the current indexing-set
         \todo sparse matrix/double indexing
         */
        void update_dim(){
            this->_dim[0] = this->_indices->size();
            this->_val->resize(this->get_dim());
            _lb->_dim[0] = std::max(_lb->_dim[0], _lb->_indices->size());
            _ub->_dim[0] = std::max(_ub->_dim[0], _ub->_indices->size());
            _lb->_val->resize(_lb->_dim[0]);
            _ub->_val->resize(_ub->_dim[0]);
        }
        
        /**
         \brief Add bounds to the current variable
         \warning Assumes identical index sets for lb and ub
         \warning Will also increase the variable index set to include new indices found in lb and ub
         \warning Will update the bounds if the same indices already exist
         @return index set of added indices.
         */
        template<typename T=type>
        indices add_bounds(const param<T>& lb, const param<T>& ub){
            /* Make sure lb and ub indices are identical */
            assert(*ub._indices==*lb._indices);
            auto ids = *lb._indices;
            /* If the variable hasn't been indexed yet */
            if(!this->_indices){
                this->index_in(ids);
                this->_lb = make_shared<func<T>>(lb);
                this->_ub = make_shared<func<T>>(ub);
                return ids;
            }
            else{
                /* Add indices and update dimension on variable and bounds*/
                auto added = this->_indices->add(ids);
                if(!added.empty()){
                    _lb->_indices->add(ids);
                    _ub->_indices->add(ids);
                    this->update_dim();
                    _lb->index_in(indices(*_lb->_indices));/* Update subexpression indices */
                    _ub->index_in(indices(*_ub->_indices));/* Update subexpression indices */
                    _lb->uneval();
                    _ub->uneval();
                    _lb->eval_all();
                    _ub->eval_all();
//                    for (auto i = 0; i< ids.size(); i++) {
//                        auto idx = lb.get_id_inst(i);
//                        auto key = ids._keys->at(idx);
//                        auto lb_idx = _lb->_indices->_keys_map->at(key);
//                        auto lb_val = lb._val->at(idx);
//                        _lb->_val->at(lb_idx) = lb_val; /* Update the bound */
//                        _lb->update_range(lb_val);
//                        this->update_range(lb_val);
//                        auto ub_idx = _ub->_indices->_keys_map->at(key);
//                        auto ub_val = ub._val->at(idx);
//                        _ub->_val->at(ub_idx) = ub_val; /* Update the bound */
//                        _ub->update_range(ub_val);
//                        this->update_range(ub_val);
//                    }
                }
                return added;
            }
        }
        
        bool operator==(const var& v) const;
        bool operator!=(const var& v) const;
        
            var in(const space& s){
                set_size(s._dim);
                if(s._dim.size()==1){ /* We can afford to build indices since this is a 1-d set */
                    return this->in(range(0,s._dim[0]-1));
                }
                _lb->set_size(s._dim);
                _ub->set_size(s._dim);
                this->_off.resize(s._dim[0],false);
                _lb->allocate_mem();
                _ub->allocate_mem();
                return *this;
            }
        
        
        
        void initialize_zero(){
            for (int i = 0; i<this->_val->size(); i++) {
                this->_val->at(i) = 0.;
            }
        };
        
        template<typename T=type, typename enable_if<is_same<T, Cpx>::value>::type* = nullptr> void initialize_binary_() {
            std::random_device dev;
            std::mt19937 engine(dev());
            for (int i = 0; i<param<type>::_val->size(); i++) {
                std::uniform_real_distribution<double> real_distribution(get_lb(i).real(),get_ub(i).real());
                std::uniform_real_distribution<double> imag_distribution(get_lb(i).imag(),get_ub(i).imag());
                if(real_distribution(engine) <= (get_ub(i).real()-get_lb(i).real())/2.){
                    param<type>::_val->at(i) = get_lb(i).real();
                }
                else {
                    param<type>::_val->at(i) = get_ub(i).real();
                }
                if(imag_distribution(engine) <= (get_ub(i).imag()-get_lb(i).imag())/2.){
                    param<type>::_val->at(i) = get_lb(i).imag();
                }
                else {
                    param<type>::_val->at(i) = get_ub(i).imag();
                }
            }
        }
        
        template<class T=type, class = typename enable_if<is_arithmetic<T>::value>::type> void initialize_binary_() {
            std::random_device dev;
            std::mt19937 engine(dev());
            for (int i = 0; i<param<type>::_val->size(); i++) {
                std::uniform_real_distribution<double> distribution(get_lb(i),get_ub(i));
                if(distribution(engine) <= (get_ub(i)-get_lb(i))/2.){
                    param<type>::_val->at(i) = get_lb(i);
                }
                else {
                    param<type>::_val->at(i) = get_ub(i);
                }
            }
        }
        
        void initialize_binary(){initialize_binary_();};
        
        void initialize_uniform(){initialize_uniform_();};
        
        template<typename T=type, typename enable_if<is_same<T, Cpx>::value>::type* = nullptr> void initialize_uniform_() {
            std::random_device dev;
            std::mt19937 engine(dev());
            for (int i = 0; i<param<type>::_val->size(); i++) {
                std::uniform_real_distribution<double> real_distribution(get_lb(i).real(),get_ub(i).real());
                std::uniform_real_distribution<double> imag_distribution(get_lb(i).imag(),get_ub(i).imag());
                param<type>::_val->at(i).real(real_distribution(engine));
                param<type>::_val->at(i).imag(imag_distribution(engine));
            }
        }
        
        template<class T=type, class = typename enable_if<is_arithmetic<T>::value>::type> void initialize_uniform_() {
            std::random_device dev;
            std::mt19937 engine(dev());
            for (int i = 0; i<param<type>::_val->size(); i++) {
                std::uniform_real_distribution<double> distribution(get_lb(i),get_ub(i));
                param<type>::_val->at(i) = distribution(engine);
            }
        }
        
        void initialize_ub(){
            for (int i = 0; i<param<type>::_val->size(); i++) {
                param<type>::_val->at(i) = get_ub(i);
            }
        };
        
        /** Initialize variable using midpoint in bounds
        */
        void initialize_midpoint(){
            for (int i = 0; i<param<type>::_val->size(); i++) {
                param<type>::_val->at(i) = (get_lb(i) + get_ub(i))/2.;
            }
        };
        
        /** Initialize variable using random value in bounds drawn from a uniform distribution
        */
        void initialize_uniform(type lb, type ub){initialize_uniform_(lb,ub);};
        
        template<typename T=type, typename enable_if<is_same<T, Cpx>::value>::type* = nullptr> void initialize_uniform_(type lb, type ub) {
            std::random_device dev;
            std::mt19937 engine(dev());
            for (int i = 0; i<param<type>::_val->size(); i++) {
                std::uniform_real_distribution<double> real_distribution(lb.real(),ub.real());
                std::uniform_real_distribution<double> imag_distribution(lb.imag(),ub.imag());
                param<type>::_val->at(i).real(real_distribution(engine));
                param<type>::_val->at(i).imag(imag_distribution(engine));
            }
        }
        
        template<class T=type, class = typename enable_if<is_arithmetic<T>::value>::type> void initialize_uniform_(type lb, type ub) {
            std::random_device dev;
            std::mt19937 engine(dev());
            for (int i = 0; i<param<type>::_val->size(); i++) {
                std::uniform_real_distribution<double> distribution(lb,ub);
                param<type>::_val->at(i) = distribution(engine);
            }
        }
        
        void vectorize() {
            if(!this->_is_vector){
                this->_name = "["+this->_name+"]";
            }
            this->_is_vector = true;
        }
        
        void transpose() {
            if(!this->_is_vector){
                this->_name = "["+this->_name+"]";
            }
            this->_is_transposed = !this->_is_transposed;
            this->_is_vector = true;
            auto temp = this->_dim[0];
            this->_dim[0] = this->_dim[1];
            this->_dim[1] = temp;
        }
        
        var tr() const {
            auto v = var(*this);
            if(!this->_is_vector){
                v._name = "["+v._name+"]";
            }
            v.constant_::transpose();
            return v;
        }
        
        var vec() const {
            auto v = var(*this);
            v._is_vector = true;
            v._name = "["+v._name+"]";
            return v;
        }
        
        /* Output */
        string to_str(size_t index1, size_t index2, int prec) {
            return this->get_name(index1,index2);
        }
        
        string to_str(){
            return this->get_name(false,false);
        }
        
        string to_str(size_t index, int prec) {
            return this->get_name(index);
        }
        string to_str_bounds(bool bounds=true, int prec = 10);        
        void print();
        void print(int prec);
        void print_vals(int prec){param<type>::print_vals(prec);};
        void print_symbolic() const{
            string str = this->_name;
            str += " ∈ [" + _lb->to_str() +"," + _ub->to_str() +"]^" + to_string(this->get_dim());
            cout << str << endl;
        }
        
        
        /* Return the scaling factor needed to make sure all bounds are in [-unit,unit]*/
        double get_scale_factor(double unit){
            auto absmax = std::max(std::abs(this->_range->first),std::abs(this->_range->second));
            if(absmax>unit)
                return unit/absmax;
            return 1;
        }
        
        /* Make sure all bounds are in [-unit,unit] */
        void scale(double unit){
            _lb->eval_all();
            _ub->eval_all();
            auto dim = this->get_dim();
            auto factor = get_scale_factor(unit);
            if(factor==1)
                return;
            for (size_t i = 0; i < dim; i++) {
                _lb->_val->at(i) = factor*_lb->_val->at(i);
                _ub->_val->at(i) = factor*_ub->_val->at(i);
            }
            _lb->_range->first = factor*_lb->_range->first;
            _ub->_range->first = factor*_ub->_range->first;
            this->_range->first = _lb->_range->first;
            this->_range->second = _ub->_range->second;
        }
        
        template<typename T=type, typename=enable_if<is_arithmetic<T>::value>>
        void copy_bounds(const shared_ptr<param_>& p){
            auto dim = p->get_dim();
            if(dim!=this->get_dim()){
                throw invalid_argument("calling function copy_bounds with non-matching dimensions");
            }
            _lb->_val->resize(dim);
            _ub->_val->resize(dim);
            for (size_t i = 0; i < dim; i++) {
                _lb->_val->at(i) = p->get_double_lb(i);
                _lb->update_range(_lb->_val->at(i));
                _ub->_val->at(i) = p->get_double_ub(i);
                _ub->update_range(_ub->_val->at(i));
            }
            this->_range->first = _lb->_range->first;
            this->_range->second = _ub->_range->second;
        }
        
        
        template<typename T=type, typename T2, typename enable_if<is_convertible<T2, T>::value>::type* = nullptr>
        void copy_bounds(const var<T2>& p){
            auto dim = p.get_dim();
            if(dim!=this->get_dim()){
                throw invalid_argument("calling function copy_bounds with non-matching dimensions");
            }
            for (size_t i = 0; i < dim; i++) {
                _lb->_val->at(i) = p._lb->_val->at(i);
                _lb->update_range(_lb->_val->at(i));
                _ub->_val->at(i) = p._ub->_val->at(i);
                _ub->update_range(_ub->_val->at(i));
            }
        }
        
        
        /** Fill x with the variable's lower bound values */
        void get_double_lb(double* x) const{get_double_lb_(x);};
        /** Fill x with the variable's upper bound values */
        void get_double_ub(double* x) const{get_double_ub_(x);};
        
        
        template<typename T=type, typename enable_if<is_arithmetic<T>::value && is_convertible<T, double>::value>::type* = nullptr>
        void get_double_lb_(double* x) const{
            auto vid = this->get_id();
            for (size_t i = 0; i < this->get_dim(); i++) {
                x[vid+i] = (double)this->get_double_lb_(i);
            }
        };
        
        template<typename T=type, typename enable_if<is_arithmetic<T>::value && is_convertible<T, double>::value>::type* = nullptr>
        void get_double_ub_(double* x) const{
            auto vid = this->get_id();
            for (size_t i = 0; i < this->get_dim(); i++) {
                x[vid+i] = (double)this->get_double_ub_(i);
            }
        };
        
        template<typename T=type, typename enable_if<is_same<T, Cpx>::value>::type* = nullptr>
        void get_double_lb_(double* x) const{
            throw invalid_argument("Cannot call get_double_lb_ with a non-arithmetic type.");
        };
        
        template<typename T=type, typename enable_if<is_same<T, Cpx>::value>::type* = nullptr>
        void get_double_ub_(double* x) const{
            throw invalid_argument("Cannot call get_double_ub_ with a non-arithmetic type.");
        };
        
        
        /** Return lower bound violation */
        double get_lb_violation(size_t i){
            return get_lb_violation_(i);
        };
        
        /** Return upper bound violation */
        double get_ub_violation(size_t i){
            return get_ub_violation_(i);
        };
        
        template<typename T=type, typename enable_if<is_arithmetic<T>::value && is_convertible<T, double>::value>::type* = nullptr>
        double get_lb_violation_(size_t i){
            return get_double_lb_(i) - this->_val->at(i);
        };
        
        template<typename T=type, typename enable_if<is_arithmetic<T>::value && is_convertible<T, double>::value>::type* = nullptr>
        double get_ub_violation_(size_t i){
            return this->_val->at(i) - get_double_ub_(i);
        };

        template<typename T=type, typename enable_if<is_same<T, Cpx>::value>::type* = nullptr>
        double get_lb_violation_(size_t i){
            throw invalid_argument("Cannot call get_lb_violation_ with a non-arithmetic type.");
        };
        
        template<typename T=type, typename enable_if<is_same<T, Cpx>::value>::type* = nullptr>
        double get_ub_violation_(size_t i){
            throw invalid_argument("Cannot call get_ub_violation_ with a non-arithmetic type.");
        };
        
        double get_double_lb(size_t i) const{return get_double_lb_(i);};
        
        double get_double_ub(size_t i) const{return get_double_ub_(i);};
        
        template<typename T=type, typename enable_if<is_same<T, Cpx>::value>::type* = nullptr>
        double get_double_lb_(size_t i) const{throw invalid_argument("Cannot call get_double_lb_ with a non-arithmetic type.");};
        
        template<typename T=type, typename enable_if<is_same<T, Cpx>::value>::type* = nullptr>
        double get_double_ub_(size_t i) const{throw invalid_argument("Cannot call get_double_ub_ with a non-arithmetic type.");};
        
        template<typename T=type, typename enable_if<is_arithmetic<T>::value && is_convertible<T, double>::value>::type* = nullptr>
        double get_double_lb_(size_t i) const{return _lb->eval(i);};
        
        template<typename T=type, typename enable_if<is_arithmetic<T>::value && is_convertible<T, double>::value>::type* = nullptr>
        double get_double_ub_(size_t i) const{return _ub->eval(i);};
        
        // should fix this by considering get_id_inst(i), but in that case _num_partns should be at least a param object
        
        int get_num_partns() const{ return *_num_partns;};
        int get_cur_partn() const{ return _cur_partn;};
        
        vector<shared_ptr<param_>> get_original_vars(){vector<shared_ptr<param_>> res; res.push_back(_original_vars[0]); res.push_back(_original_vars[1]); return res;};
        bool is_lifted() const{return _lift;};
        bool get_in_SOC_partn() const{return _in_SOC_partn;};
        void set_in_SOC_partn(bool in_SOC_partn) {this->_in_SOC_partn = in_SOC_partn;};
        
    };
    
    var<Cpx> conj(const var<Cpx>& p);
    var<Cpx> ang(const var<Cpx>& p);
    var<Cpx> sqrmag(const var<Cpx>& p);
    var<Cpx> real(const var<Cpx>& p);
    var<Cpx> imag(const var<Cpx>& p);
    
}
#endif /* var_h */
