////
////  var.h
////  Gravity
////
////  Created by Hijazi, Hassan on 21/10/16.
////
////
//
//#ifndef var_h
//#define var_h
//
//#include <gravity/param.h>
//#include <gravity/Net.h>
//#include <stdio.h>
//#include <string>
//#include <set>
//#include <list>
//#include <limits>
//#include <random>
//
//
//
//using namespace std;
//
//namespace gravity {
//    
//    
//    class func_;
//    
//    
//    /** A variable can be a bool, a short, an int, a float, a double, a long double or a complex<double>*/
//    template<typename type = double>
//    // define variable as a parameter with bounds
//    class var: public param<type>{
//        
//    public:
//        shared_ptr<func_>   _lb; /**< Lower Bound */
//        shared_ptr<func_>   _ub; /**< Upper Bound */
//        bool _in_q_cone = false; /**< Used by Mosek */
//        bool _psd = false; /**< Has to be positive semidefinite */
//        
//        /* Constructors */
//        //@{
//        /** Unbounded variable constructor */
//        var();
//        ~var() {};
//        var(const string& name);
//        var(const string& name, Sign s);
//        var(const var<type>& v);
//        var(var<type>&& v);
//        //@}
//        
//        //@{
//        /** Bounded variable constructor */
//        template<class T=type, class = typename enable_if<is_same<T, Cpx>::value>::type> var(const string& name, const T& lb, const T& ub):param<T>(name){
//            param<T>::set_type(var_c);
//            _lb = make_shared<func_>(constant<T>(lb));
//            _ub = make_shared<func_>(constant<T>(ub));
//            param<T>::_range->first = lb;
//            param<T>::_range->second = ub;
//        }
//        
//        template<class T=type, class = typename enable_if<is_same<T, Cpx>::value>::type> var(const string& name, T&& lb, T&& ub):param<T>(name){
//            param<T>::set_type(var_c);
//            _lb = make_shared<func_>(constant<T>(lb));
//            _ub = make_shared<func_>(constant<T>(ub));
//            param<T>::_range->first = lb;
//            param<T>::_range->second = ub;
//        }
//        
//        var(const string& name, type lb, type ub):var(name){
//            _lb = make_shared<func_>(constant<type>(lb));
//            _ub = make_shared<func_>(constant<type>(ub));
//            param<type>::_range->first = lb;
//            param<type>::_range->second = ub;
//        };
//        
//        var(const string& name, const param<type>& lb, const param<type>& ub):var(name){
//            _lb = make_shared<func_>(lb);
//            _ub = make_shared<func_>(ub);
//            param<type>::_range->first = lb._range->first;
//            param<type>::_range->second = ub._range->second;
//        };
//        
//        template<class T=type, class = typename enable_if<is_arithmetic<T>::value>::type>
//        var(const string& name, const param<T>& sb):var(name) {
//            _lb = make_shared<func_>(-1*sb);
//            _ub = make_shared<func_>(sb);
//            param<type>::_range->first = min(-1*sb._range->first, -1*sb._range->second);
//            param<type>::_range->second = sb._range->second;
//        }
//        
//        var(const string& name, const func_& lb, const func_& ub):var(name){
//            _lb = make_shared<func_>(move(lb));
//            _ub = make_shared<func_>(move(ub));
//        };
//        
//        var(const string& name, func_&& lb, func_&& ub):var(name) {
//            _lb = make_shared<func_>(lb);
//            _ub = make_shared<func_>(ub);
//        };
//        
//        //@}
//        
//        
//        void initialize_all(type v) {
//            this->set_val_all(v);
//        }
//        
//        /* Retrieve specified indexed variable. */
//        template<typename... Args>
//        var operator()(size_t i, size_t j) {
//            bool indexed = param<type>::_indices!=nullptr;
//            var<type> res(*this);
//            res.param<type>::operator=(param<type>::operator()(i, j));
//            if(!indexed && !res._lb->is_number()){
//                (res._lb->in(*res._indices));
//            }
//            if(!indexed && !res._ub->is_number()){
//                (res._ub->in(*res._indices));
//            }
//            return res;
//        }
//        
//        template<typename... Args>
//        var operator()(string t1, Args&&... args) {
//            bool indexed = param<type>::_indices!=nullptr;
//            var<type> res(*this);
//            res.param<type>::operator=(param<type>::operator()(t1, args...));
//            res.param<type>::set_type(var_c);
//            if(!indexed & !res._lb->is_number()){
//                (res._lb->in(*res._indices));
//            }
//            if(!indexed & !res._ub->is_number()){
//                (res._ub->in(*res._indices));
//            }
//            return res;
//        }
//        
//        template<typename... Args>
//        var operator()(size_t i) {
//            bool indexed = param<type>::_indices!=nullptr;
//            var<type> res(*this);
//            res.param<type>::operator=(param<type>::operator()(i));
//            res.param<type>::set_type(var_c);
//            if(!indexed & !res._lb->is_number()){
//                (res._lb->in(*res._indices));
//            }
//            if(!indexed & !res._ub->is_number()){
//                (res._ub->in(*res._indices));
//            }
//            return res;
//        }
//        
//        var in_pairs(){
//            var<type> res(*this);
//            res._name += ".in_pairs";
//            res._indices->_type = in_pairs_;
//            return res;
//        }
//        
//        var from(){
//            var<type> res(*this);
//            res._name += ".from";
//            res._indices->_type = from_;
//            return res;
//        }
//        
//        
//        var to(){
//            var<type> res(*this);
//            res._name += ".to";
//            res._indices->_type = to_;
//            return res;
//        }
//        
//        
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
//        
//        
//        
//        var in_pairs(const indices& ids) {
//            bool indexed = param<type>::_indices!=nullptr;
//            var<type> res(*this);
//            res.param<type>::operator=(param<type>::in(ids));
//            if(!indexed && !res._lb->is_number()){
//                (res._lb->in(*res._indices));
//            }
//            if(!indexed && !res._ub->is_number()){
//                (res._ub->in(*res._indices));
//            }
//            return res;
//        }
//        
//        template<typename... Args>
//        var in(const indices& vec1, Args&&... args) {
//            bool indexed = param<type>::_indices!=nullptr;
//            var<type> res(*this);
//            res.param<type>::operator=(param<type>::in(vec1, forward<Args>(args)...));
//            if(!indexed && !res._lb->is_number()){
//                (res._lb->in(*res._indices));
//            }
//            if(!indexed && !res._ub->is_number()){
//                (res._ub->in(*res._indices));
//            }
//            return res;
//        }
//        
//        
//        //    var in(const node_pairs& np){
//        //        return this->in(np._keys);
//        //    }
//        
//        var in_arcs(const vector<Node*>& vec) {
//            var<type> res(*this);
//            res.param<type>::operator=(param<type>::in_arcs(vec));
//            return res;
//        }
//        
//        var out_arcs(const vector<Node*>& vec) {
//            var<type> res(*this);
//            res.param<type>::operator=(param<type>::out_arcs(vec));
//            return res;
//        }
//        
//        
//        var in_aux(const vector<Node*>& vec, const string& aux_type) {
//            var<type> res(*this);
//            res.param<type>::operator=(param<type>::in_aux(vec,aux_type));
//            return res;
//        }
//        
//        
//        var excl(size_t index) {
//            var<type> res(*this);
//            res._indices->_type = excl_;
//            return res;
//        }
//        
//        
//        /* Create a vector of variables indexed based on nodes from bags of size bag_size
//         e.g., given bag { 1, 5, 7 } index the first variable (1), the second (5) and the last (7)
//         */
//        vector<var> in_bags(const vector<vector<Node*>>& bags, size_t bag_size);
//        
//        /* Create a vector of variables indexed as pair of nodes from bags of size bag_size
//         e.g., given bag { 1, 5, 7 } index the first variable (1,5), the second (5,7) and the last (1,7)
//         */
//        vector<var> pairs_in_bags(const vector<vector<Node*>>& bags, size_t bag_size);
//        
//        /* Querries */
//        
//        type    get_lb(size_t i = 0) const;
//        
//        type    get_ub(size_t i = 0) const;
//        
//        
//        template<typename T=type,
//        typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr>
//        bool is_bounded_above(size_t i = 0) const;
//        
//        template<typename T=type,
//        typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr>
//        bool is_bounded_below(size_t i = 0) const;
//        
//        template<typename T=type,
//        typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr>
//        bool is_constant(size_t i=0) const;
//        
//        template<typename T=type,
//        typename std::enable_if<is_arithmetic<T>::value>::type* = nullptr>
//        Sign get_sign(size_t idx = 0) const;
//        
//        /* Modifiers */
//        void    set_size(vector<size_t>);
//        void    set_size(size_t s);
//        
//        void    add_bounds(type lb, type ub);
//        void    add_lb_only(type v); /**< Adds a new lower bound and infinity in the corresponding upper bound*/
//        void    add_ub_only(type v); /**< Adds a new upper bound and -infinity in the corresponding lower bound*/
//        
//        
//        void    set_lb(type v);
//        void    set_ub(type v);
//        
//        void in_q_cone(); /**< States that the variable is contained in a quadratic cone, for mosek */
//        
//        
//        /* Operators */
//        var& operator=(const var& v);
//        var& operator=(var&& v);
//        
//        var& operator=(type v) {
//            this->set_val_all(v);
//            return *this;
//        }
//        
//        bool operator==(const var& v) const;
//        bool operator!=(const var& v) const;
//        
//            var in(const space& s){
//                set_size(s._dim);
//                if(s._dim.size()==1){ /* We can afford to build indices since this is a 1-d set */
//                    this->_indices = make_shared<indices>(indices(0,s._dim[0]-1));
//                }
//                return *this;
//            }
//        
//        template<class T=type, class = typename enable_if<is_arithmetic<T>::value>::type> void initialize_uniform() {
//            std::default_random_engine generator;
//            for (int i = 0; i<param<type>::_val->size(); i++) {
//                std::uniform_real_distribution<double> distribution(get_lb(i),get_ub(i));
//                param<type>::_val->at(i) = distribution(generator);
//            }
//        }
//        
//        var tr() const {
//            auto v = var(*this);
//            if(!this->_is_vector){
//                v._name = "["+v._name+"]";
//            }
//            v.constant_::transpose();
//            return v;
//        }
//        
//        var vec() const {
//            auto v = var(*this);
//            v._is_vector = true;
//            v._name = "["+v._name+"]";
//            return v;
//        }
//        
//        /* Output */
//        string to_str(bool bounds=true, int prec = 10) const;
//        void print(bool bounds=true, int prec = 10) const;
//        
//    };
//    
//    var<Cpx> conj(const var<Cpx>& p);
//    var<Cpx> ang(const var<Cpx>& p);
//    var<Cpx> sqrmag(const var<Cpx>& p);
//    var<Cpx> real(const var<Cpx>& p);
//    var<Cpx> imag(const var<Cpx>& p);
//    
//}
//#endif /* var_h */
