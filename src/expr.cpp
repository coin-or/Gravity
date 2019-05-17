////
////  expr.cpp
////  Gravity
////
////  Created by Hijazi, Hassan on 24/10/16.
////
////
#include <cmath>
#include <gravity/expr.h>
#include <gravity/func.h>

using namespace std;
string operator_str(gravity::OperatorType ot){
    switch (ot) {
        case gravity::log_:
            return "log";
        case gravity::exp_:
            return "exp";
        case gravity::cos_:
            return "cos";
        case gravity::sin_:
            return "sin";
        case gravity::tan_:
            return "tan";
        case gravity::atan2_:
            return "atan2";
        case gravity::acos_:
            return "acos";
        case gravity::asin_:
            return "asin";
        case gravity::sqrt_:
            return "sqrt";
        case gravity::relu_:
            return "ReLU";
        case gravity::unit_step_:
            return "UnitStep";
        case gravity::min_:
            return "min";
        case gravity::max_:
            return "max";
        default:
            break;
    }
    throw invalid_argument("Unsupported unitary operator");
}


namespace gravity {
    template<typename type>
    bexpr<type>::bexpr(OperatorType otype, shared_ptr<constant_> lson, shared_ptr<constant_> rson){
        _otype = otype;
        _lson = lson;
        _rson = rson;
        this->_type = bexp_c;
        this->_to_str = to_str();
        if(otype==product_){
            this->_dim[0] = _lson->_dim[0];
            this->_dim[1] = _rson->_dim[1];
            
            /* Instructions above work for matrix and vector dot products, below we check if it's a component-wise vector,matrix product */
            if(otype==product_ && !_lson->is_matrix() && _rson->is_matrix()){
                this->_dim[0] = _rson->_dim[0];
            }
            if(otype==product_ && _lson->is_matrix() && !_rson->is_matrix() && _rson->_is_transposed){
                this->_dim[1] = _lson->_dim[1];
            }
            if(this->is_matrix()){
                this->_is_vector = true;
            }
        }
        else {
            this->_dim[0] = std::max(this->_dim[0], _lson->_dim[0]);
            this->_dim[0] = std::max(this->_dim[0], _rson->_dim[0]);
            this->_dim[1] = std::max(this->_dim[1], _lson->_dim[1]);
            this->_dim[1] = std::max(this->_dim[1], _rson->_dim[1]);
        }
        shared_ptr<pair<type,type>> lson_range, rson_range;
        if (_lson->is_function()) {
            lson_range = static_pointer_cast<func<type>>(_lson)->_range;
        }
        else if(_lson->is_expr()){
            lson_range = static_pointer_cast<expr<type>>(_lson)->_range;
        }
        else if (_lson->is_param() || _lson->is_var() ){
            lson_range = static_pointer_cast<param<type>>(_lson)->_range;
        }
        if (_rson->is_function()) {
            rson_range = static_pointer_cast<func<type>>(_rson)->_range;
        }
        else if(_rson->is_expr()){
            rson_range = static_pointer_cast<expr<type>>(_rson)->_range;
        }
        else if (_rson->is_param() || _rson->is_var() ){
            rson_range = static_pointer_cast<param<type>>(_rson)->_range;
        }
        switch(otype){
            case product_:{
                this->_range = get_product_range(lson_range, rson_range);
            }
                break;
            case plus_:{
                this->_range = get_plus_range(lson_range, rson_range);
            }
                break;
            case minus_:{
                this->_range = get_minus_range(lson_range, rson_range);
            }
                break;
            case div_:{
                this->_range = get_div_range(lson_range, rson_range);
            }
                break;
            default:
                break;
        }
    };
    template class bexpr<bool>;
    template class bexpr<short>;
    template class bexpr<int>;
    template class bexpr<float>;
    template class bexpr<double>;
    template class bexpr<long double>;
    template class bexpr<Cpx>;

}
