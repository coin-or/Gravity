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
            _lb = make_shared<func<type>>(constant<type>(numeric_limits<type>::lowest()));
            _ub = make_shared<func<type>>(constant<type>(numeric_limits<type>::max()));
            this->_range->first = _lb->_range->first;
            this->_range->second = _ub->_range->second;
        }
        
        template<class T=type, class = typename enable_if<is_same<T, Cpx>::value>::type>
        var(const string& name){
            constant_::set_type(var_c);
            this->_name = name;
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
        //@}
        
        shared_ptr<param_> pcopy() const{return make_shared<var>(*this);};
        
        shared_ptr<constant_> copy() const{return make_shared<var>(*this);};
        
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
            return res;
        }
        
        /** let p share the values and indices of current var */
        void share_vals(param<type>& p){
            if(this->_indices){
                p._indices = this->_indices;
            }
            p._dim[0] = this->_dim[0];
            p._dim[1] = this->_dim[1];
            p._val = this->_val;
            p._range = this->_range;
        }
        
        template<typename... Args>
        var operator()(size_t i) {
            bool indexed = param<type>::_indices!=nullptr;
            var<type> res(*this);
            res.param<type>::operator=(param<type>::operator()(i));
            res.param<type>::set_type(var_c);
            if(!indexed & !res._lb->is_number()){
                (res._lb->in(*res._indices));
            }
            if(!indexed & !res._ub->is_number()){
                (res._ub->in(*res._indices));
            }
            return res;
        }
        
        
        template<typename... Args>
        var operator[](size_t i) {
            auto res = (*this)(i-1);
            res._name = this->_name+"["+to_string(i)+"]";
            return res;
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
            res.param<type>::operator=(param<type>::in(ids));
            if(!indexed && !res._lb->is_number()){
                (res._lb->in(*res._indices));
            }
            if(!indexed && !res._ub->is_number()){
                (res._ub->in(*res._indices));
            }
            return res;
        }
        
        void real_imag(const var<>& pr, const var<>& pi){
            this->_real = make_shared<var<>>(pr);
            this->_imag = make_shared<var<>>(pi);
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


        template<typename... Args>
        var in(const indices& vec1, Args&&... args) {
//            bool indexed = param<type>::_indices!=nullptr;
            var<type> res(*this);
            res.param<type>::operator=(param<type>::in(vec1, forward<Args>(args)...));//TODO assert lb dim = res dim
            res._type = var_c;
//            if(!indexed && !res._lb->is_number()){
//                (res._lb->index_in(*res._indices));
//            }
//            if(!indexed && !res._ub->is_number()){
//                (res._ub->index_in(*res._indices));
//            }
            res._lb->allocate_mem();
            res._ub->allocate_mem();
            return res;
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
        vector<var> in_bags(const vector<vector<Node*>>& bags, size_t bag_size);
        
        /* Create a vector of variables indexed as pair of nodes from bags of size bag_size
         e.g., given bag { 1, 5, 7 } index the first variable (1,5), the second (5,7) and the last (1,7)
         */
        vector<var> pairs_in_bags(const vector<vector<Node*>>& bags, size_t bag_size);
        
        /* Querries */
        
        type    get_lb(size_t i = 0) const;
        
        type    get_ub(size_t i = 0) const;
        
        
        template<typename T=type,
        typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr>
        bool is_bounded_above(size_t i = 0) const;
        
        template<typename T=type,
        typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr>
        bool is_bounded_below(size_t i = 0) const;
        
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
        
        void in_q_cone(); /**< States that the variable is contained in a quadratic cone, for mosek */
        
        
        /* Operators */
        var& operator=(const var& v);
        var& operator=(var&& v);
        
        var& operator=(type v) {
            this->set_val(v);
            return *this;
        }
        
        bool operator==(const var& v) const;
        bool operator!=(const var& v) const;
        
            var in(const space& s){
                set_size(s._dim);
                if(s._dim.size()==1){ /* We can afford to build indices since this is a 1-d set */
                    this->_indices = make_shared<indices>(indices(0,s._dim[0]-1));
                }
                _lb->set_size(s._dim);
                _ub->set_size(s._dim);
                _lb->allocate_mem();
                _ub->allocate_mem();
                return *this;
            }
        
        void initialize_zero(){
            for (int i = 0; i<this->_val->size(); i++) {
                this->_val->at(i) = 0.;
            }
        };
        
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
        
        void initialize_av(){
            for (int i = 0; i<param<type>::_val->size(); i++) {
                param<type>::_val->at(i) = (get_lb(i) + get_ub(i))/2.;
            }
        };
        
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
        
        
        template<typename T=type, typename=enable_if<is_arithmetic<T>::value>>
        void copy_bounds(const shared_ptr<param_>& p){
            auto dim = p->get_dim();
            if(dim!=this->get_dim()){
                throw invalid_argument("calling function copy_bounds with non-matching dimensions");
            }
            _lb->reset();
            _ub->reset();
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
        void set_double_lb(double* x){set_double_lb_(x);};
        /** Fill x with the variable's upper bound values */
        void set_double_ub(double* x){set_double_ub_(x);};
        
        
        template<typename T=type, typename enable_if<is_arithmetic<T>::value && is_convertible<T, double>::value>::type* = nullptr>
        void set_double_lb_(double* x){
            auto vid = this->get_id();
            for (size_t i = 0; i < this->get_dim(); i++) {
                x[vid+i] = (double)this->get_double_lb_(i);
            }
        };
        
        template<typename T=type, typename enable_if<is_arithmetic<T>::value && is_convertible<T, double>::value>::type* = nullptr>
        void set_double_ub_(double* x){
            auto vid = this->get_id();
            for (size_t i = 0; i < this->get_dim(); i++) {
                x[vid+i] = (double)this->get_double_ub_(i);
            }
        };
        
        template<typename T=type, typename enable_if<is_same<T, Cpx>::value>::type* = nullptr>
        void set_double_lb_(double* x){
            throw invalid_argument("Cannot call get_lb_violation_ with a non-arithmetic type.");
        };
        
        template<typename T=type, typename enable_if<is_same<T, Cpx>::value>::type* = nullptr>
        void set_double_ub_(double* x){
            throw invalid_argument("Cannot call get_ub_violation_ with a non-arithmetic type.");
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
        
    };
    
    var<Cpx> conj(const var<Cpx>& p);
    var<Cpx> ang(const var<Cpx>& p);
    var<Cpx> sqrmag(const var<Cpx>& p);
    var<Cpx> real(const var<Cpx>& p);
    var<Cpx> imag(const var<Cpx>& p);
    
}
#endif /* var_h */
