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
namespace gravity{


    string operator_str(OperatorType ot){
        switch (ot) {
            case log_:
                return "log";
            case exp_:
                return "exp";
            case cos_:
                return "cos";
            case sin_:
                return "sin";
            case tan_:
                return "tan";
            case sqrt_:
                return "sqrt";
            case relu_:
                return "ReLU";
            case unit_step_:
                return "UnitStep";
            default:
                break;
        }
        throw invalid_argument("Unsupported unitary operator");
    }

    uexpr::uexpr(const uexpr& exp){
        *this = exp;
    }

    uexpr::uexpr(uexpr&& exp){
        *this = move(exp);
    }
    
    uexpr::uexpr(OperatorType ot, shared_ptr<constant_> son){
        _otype = ot;
        _son = son;
        _type = uexp_c;
        _dim[0] = son->_dim[0];
        _dim[1] = son->_dim[1];
        _to_str = to_str();
        _is_vector = son->_is_vector;
        _is_transposed = son->_is_transposed;
    }


    uexpr& uexpr::operator=(uexpr&& exp){
        _type = uexp_c;
        _son = move(exp._son);
        _otype = exp._otype;
        _to_str = exp._to_str;
        _coef = exp._coef;
        _is_vector = exp._is_vector;
        _is_transposed = exp._is_transposed;
        _dim[0] = exp._dim[0]; _dim[1] = exp._dim[1];
        return *this;
    }


    bool uexpr::operator==(const uexpr &c)const{
//        return (_otype == c._otype && equals(_son,c._son));
        return (this->_to_str.compare(c._to_str)==0);
    }



    /* BINARY EXPRESSIONS */

    bool bexpr::operator==(const bexpr &c)const{
//        return (_otype == c._otype && equals(_lson,c._lson) && equals(_rson,c._rson));
        return (this->_to_str.compare(c._to_str)==0);
    }


    bexpr::bexpr(){
        _type = bexp_c;
    }


    bexpr::bexpr(const bexpr& exp){ /**< Copy constructor from binary expression tree */
        *this = exp;
    };

    bexpr::bexpr(bexpr&& exp){ /**< Move constructor from binary expression tree */
        *this = move(exp);
    };
    

    bexpr& bexpr::operator=(bexpr&& exp){
        _type = bexp_c;
        _lson = move(exp._lson);
        _rson = move(exp._rson);
        _otype = exp._otype;
        _to_str = exp._to_str;
        _coef = exp._coef;
        _is_vector = exp._is_vector;
        _is_transposed = exp._is_transposed;
        _dim[0] = exp._dim[0]; _dim[1] = exp._dim[1];
        return *this;
    }


    /* UNARY EXPRESSIONS */

    uexpr::uexpr(){
        _type = uexp_c;
    }

    void uexpr::print() {
        cout << _to_str << endl;
    }

    void bexpr::print() {
        cout << _to_str << endl;
    }
    
    void expr::allocate_mem(){
        if (is_uexpr()) {
            auto ue = (uexpr*)(this);
            ue->_son->allocate_mem();
        }
        else {
            auto be = (bexpr*)(this);
            be->_lson->allocate_mem();
            be->_rson->allocate_mem();
        }
    }
    
    void expr::propagate_dim(size_t d){
        if(_is_transposed){
            _dim[1] = d;
        }
        else {
            _dim[0] = d;
        }
        if (is_uexpr()) {
            auto ue = (uexpr*)(this);
            ue->_son->propagate_dim(d);
        }
        else {
            auto be = (bexpr*)(this);
            be->_lson->propagate_dim(d);
            be->_rson->propagate_dim(d);
        }
    }
    
    
    
    uexpr& uexpr::operator=(const uexpr& exp){
        _type = uexp_c;
        _son = exp._son->copy();
        _otype = exp._otype;
        _to_str = exp._to_str;
        _coef = exp._coef;
        _is_vector = exp._is_vector;
        _is_transposed = exp._is_transposed;
        _dim[0] = exp._dim[0]; _dim[1] = exp._dim[1];
        return *this;
    }
    
    
    string uexpr::to_str(){
        string str;
        if (_coef!=1) {
            if (_coef!=-1) {
                str+= to_string_with_precision(_coef,3);
            }
            else {
                str+= "-";
            }
            str+="(";
        }
        str += operator_str(_otype) +"("+_son->to_str()+")";
        if (_coef!=1) {
            str += ")";
        }
        return str;
    }

        
    string uexpr::to_str(int prec){
        string str;
        if (_coef!=1) {
            if (_coef!=-1) {
                str+= to_string_with_precision(_coef,prec);
            }
            else {
                str+= "-";
            }
            str+="(";
        }
        str += operator_str(_otype) +"("+_son->to_str(prec)+")";
        if (_coef!=1) {
            str += ")";
        }
        return str;
    }
    
    string uexpr::to_str(size_t inst, int prec) {
        string str;
        if (_coef!=1) {
            if (_coef!=-1) {
                str+= to_string_with_precision(_coef,prec);
            }
            else {
                str+= "-";
            }
            str+="(";
        }
        str += operator_str(_otype) +"("+_son->to_str(inst,prec)+")";
        if (_coef!=1) {
            str += ")";
        }
        return str;
    }
    
    string uexpr::to_str(size_t inst1, size_t inst2, int prec) {
        string str;
        if (_coef!=1) {
            if (_coef!=-1) {
                str+= to_string_with_precision(_coef,prec);
            }
            else {
                str+= "-";
            }
            str+="(";
        }
        str += operator_str(_otype) +"("+_son->to_str(inst1,inst2,prec)+")";
        if (_coef!=1) {
            str += ")";
        }
        return str;
    }
//    vector<shared_ptr<param_>> uexpr::get_nl_vars() const{
//        vector<shared_ptr<param_>> res;
//        if (_son->is_function()) {
//            auto vars = ((func_*)_son.get())->get_vars();
//            for (auto &pairs: vars){
//                res.push_back(pairs.second.first);
//            }
//        }
//        else if(_son->is_uexpr()) {
//            return ((uexpr*)_son.get())->get_nl_vars();
//        }
//        else if (_son->is_bexpr()){
//            return ((bexpr*)_son.get())->get_nl_vars();
//        }
//        return res;
//    }
    
//    func_ expr::get_derivative(const param_ &v) const{
//        if(is_uexpr()){
//            return ((uexpr*)this)->get_derivative(v);
//        }
//        else {
//            return ((bexpr*)this)->get_derivative(v);
//        }
//    }
//    func_ uexpr::get_derivative(const param_ &v) const{
//        // Unary operators
//        // f(g(x))' = f'(g(x))*g'(x).
//        switch (_otype) {
//            case cos_:
//                return _coef*-1*_son->get_derivative(v)* sin(*_son);
//                break;
//            case sin_:
//                return _coef*_son->get_derivative(v)*cos(*_son);
//                break;
//            case sqrt_:
//                return _coef*_son->get_derivative(v)/(2*sqrt(*_son));
//                break;
//            case exp_:{
//                auto rson = expo(*_son);
//                if (rson.is_matrix() || rson._is_vector) {
//                    rson.transpose();
//                }
//                return _coef*_son->get_derivative(v)*rson;
//                break;
//            }
//            case log_:
//                return _coef*_son->get_derivative(v)/(*_son);
//                break;
//            default:
//                throw invalid_argument("Unsupported unary operation");
//                break;
//        }
//        return func_();
//    }
//    vector<shared_ptr<param_>> bexpr::get_nl_vars() const{
//        vector<shared_ptr<param_>> res;
//        if (_lson->is_function()) {
//            auto vars = ((func_*)_lson.get())->get_vars();
//            for (auto &pairs: vars){
//                res.push_back(pairs.second.first);
//            }
//        }
//        else if(_lson->is_uexpr()) {
//            res = ((uexpr*)_lson.get())->get_nl_vars();
//        }
//        else if (_lson->is_bexpr()){
//            res = ((bexpr*)_lson.get())->get_nl_vars();
//        }
//        if (_rson->is_function()) {
//            auto vars = ((func_*)_rson.get())->get_vars();
//            for (auto &pairs: vars){
//                res.push_back(pairs.second.first);
//            }
//        }
//        else if(_rson->is_uexpr()) {
//            auto vars = ((uexpr*)_rson.get())->get_nl_vars();
//            res.insert(res.end(), vars.begin(), vars.end() );
//        }
//        else if (_rson->is_bexpr()){
//            auto vars = ((bexpr*)_rson.get())->get_nl_vars();
//            res.insert(res.end(), vars.begin(), vars.end() );
//        }
//        return res;
//    }
    
//    func_ bexpr::get_derivative(const param_ &v) const{
//        switch (_otype) {
//            case plus_:
//                return _coef*(_lson->get_derivative(v) + _rson->get_derivative(v));
//                break;
//            case minus_:
//                return _coef*(_lson->get_derivative(v) - _rson->get_derivative(v));
//                break;
//            case product_:{
//                return _coef*(_lson->get_derivative(v)*(*_rson) + (*_lson)*_rson->get_derivative(v));
//                //                else {
//                //                    auto lson = *_lson;
//                ////                    rson_df._is_transposed = !rson_df._is_transposed ;
//                ////                    lson._is_transposed = !lson._is_transposed;
//                ////                    lson.transpose();
//                //                    return _coef*(_lson->get_derivative(v)*(*_rson) + rson_df*lson);
//                //                }
//                // f'g + fg'
//                break;
//            }
//            case div_:
//                return _coef*((_lson->get_derivative(v)*(*_rson) - (*_lson)*_rson->get_derivative(v))/((*_rson)*(*_rson)));
//                // (f'g - fg')/g^2
//                break;
//            case power_:
//                if (!_rson->is_number()) {
//                    throw invalid_argument("Function in exponent not supported yet.\n");
//                }
//                //                auto exponent = poly_eval(_rson);
//                //                return (exponent*get_poly_derivative(_lson,v)*power(*_lson, exponent-1));// nf'f^n-1
//                break;
//            default:
//                throw invalid_argument("unsupported operation");
//                break;
//        }
//        return func_();
//    }
    
    string bexpr::to_str() {
        string str;
        if (_coef!=1) {
            if (_coef!=-1) {
                str+= to_string_with_precision(_coef, 3);
            }
            else {
                str+= "-";
            }
            str+="(";
        }
        if((_otype==product_ || _otype==div_) && (_lson->get_type()==uexp_c || _lson->get_type()==bexp_c)) {
            str += "(";
            str+= _lson->to_str();
            str += ")";
        }
        else
            str+= _lson->to_str();
        
        if (_otype==plus_) {
            str+= " + ";
        }
        if (_otype==minus_) {
            str+= " - ";
        }
        if (_otype==product_) {
            str+= " * ";
        }
        if (_otype==div_) {
            str+= "/";
        }
        
        if (_otype==power_) {
            str+= "^";
        }
        
        if (_otype==plus_ || (_rson->get_type()!=uexp_c && _rson->get_type()!=bexp_c)) {
            str+= _rson->to_str();
        }
        else {
            str+= "(";
            str+= _rson->to_str();
            str+= ")";
        }
        if (_coef!=1) {
            str += ")";
        }
        return str;
    }
    
    string bexpr::to_str(int prec) {
        string str;
        if (_coef!=1) {
            if (_coef!=-1) {
                str+= to_string_with_precision(_coef, prec);
            }
            else {
                str+= "-";
            }
            str+="(";
        }
        if((_otype==product_ || _otype==div_) && (_lson->get_type()==uexp_c || _lson->get_type()==bexp_c)) {
            str += "(";
            str+= _lson->to_str(prec);
            str += ")";
        }
        else
            str+= _lson->to_str(prec);
        
        if (_otype==plus_) {
            str+= " + ";
        }
        if (_otype==minus_) {
            str+= " - ";
        }
        if (_otype==product_) {
            str+= " * ";
        }
        if (_otype==div_) {
            str+= "/";
        }
        
        if (_otype==power_) {
            str+= "^";
        }
        
        if (_otype==plus_ || (_rson->get_type()!=uexp_c && _rson->get_type()!=bexp_c)) {
            str+= _rson->to_str(prec);
        }
        else {
            str+= "(";
            str+= _rson->to_str(prec);
            str+= ")";
        }
        if (_coef!=1) {
            str += ")";
        }
        return str;
    }
    
    string bexpr::to_str(size_t inst,int prec) {
        string str;
        if (_coef!=1) {
            if (_coef!=-1) {
                str+= to_string_with_precision(_coef, prec);
            }
            else {
                str+= "-";
            }
            str+="(";
        }
        if((_otype==product_ || _otype==div_) && (_lson->get_type()==uexp_c || _lson->get_type()==bexp_c)) {
            str += "(";
            str+= _lson->to_str(inst,prec);
            str += ")";
        }
        else
            str+= _lson->to_str(inst,prec);
        
        if (_otype==plus_) {
            str+= " + ";
        }
        if (_otype==minus_) {
            str+= " - ";
        }
        if (_otype==product_) {
            str+= " * ";
        }
        if (_otype==div_) {
            str+= "/";
        }
        
        if (_otype==power_) {
            str+= "^";
        }
        
        if (_otype==plus_ || (_rson->get_type()!=uexp_c && _rson->get_type()!=bexp_c)) {
            str+= _rson->to_str(inst,prec);
        }
        else {
            str+= "(";
            str+= _rson->to_str(inst,prec);
            str+= ")";
        }
        if (_coef!=1) {
            str += ")";
        }
        return str;
    }
    
    string bexpr::to_str(size_t inst1,size_t inst2,int prec) {
        string str;
        if (_coef!=1) {
            if (_coef!=-1) {
                str+= to_string_with_precision(_coef, prec);
            }
            else {
                str+= "-";
            }
            str+="(";
        }
        if((_otype==product_ || _otype==div_) && (_lson->get_type()==uexp_c || _lson->get_type()==bexp_c)) {
            str += "(";
            str+= _lson->to_str(inst1,inst2,prec);
            str += ")";
        }
        else
            str+= _lson->to_str(inst1,inst2,prec);
        
        if (_otype==plus_) {
            str+= " + ";
        }
        if (_otype==minus_) {
            str+= " - ";
        }
        if (_otype==product_) {
            str+= " * ";
        }
        if (_otype==div_) {
            str+= "/";
        }
        
        if (_otype==power_) {
            str+= "^";
        }
        
        if (_otype==plus_ || (_rson->get_type()!=uexp_c && _rson->get_type()!=bexp_c)) {
            str+= _rson->to_str(inst1,inst2,prec);
        }
        else {
            str+= "(";
            str+= _rson->to_str(inst1,inst2,prec);
            str+= ")";
        }
        if (_coef!=1) {
            str += ")";
        }
        return str;
    }
    
    bool bexpr::is_inner_product() const{
        return _otype==product_ && (_lson->get_dim(1)==_rson->get_dim(0) || (_lson->_is_transposed && _lson->get_dim(0)==_rson->get_dim(0)));
    }
    

    bexpr::bexpr(OperatorType otype, shared_ptr<constant_> lson, shared_ptr<constant_> rson){
        _otype = otype;
        _lson = lson;
        _rson = rson;
        _type = bexp_c;
        _to_str = to_str();
        if(otype==product_){
            _dim[0] = _lson->_dim[0];
            _dim[1] = _rson->_dim[1];
        
            /* Instructions above work for matrix and vector dot products, below we check if it's a component-wise vector,matrix product */
            if(otype==product_ && !_lson->is_matrix() && _rson->is_matrix()){
                _dim[0] = _rson->_dim[0];
            }
            if(otype==product_ && _lson->is_matrix() && !_rson->is_matrix() && _rson->_is_transposed){
                _dim[1] = _lson->_dim[1];
            }
            if(is_matrix()){
                _is_vector = true;
            }
        }
        else {
            _dim[0] = std::max(_dim[0], _lson->_dim[0]);
            _dim[0] = std::max(_dim[0], _rson->_dim[0]);
            _dim[1] = std::max(_dim[1], _lson->_dim[1]);
            _dim[1] = std::max(_dim[1], _rson->_dim[1]);
        }
    };
    
    bexpr& bexpr::operator=(const bexpr& exp){
        _type = bexp_c;
        _lson = exp._lson->copy();
        _rson = exp._rson->copy();
        _otype = exp._otype;
        _to_str = exp._to_str;
        _coef = exp._coef;
        _is_vector = exp._is_vector;
        _is_transposed = exp._is_transposed;
        _dim[0] = exp._dim[0]; _dim[1] = exp._dim[1];
        return *this;
    }
}
