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


//    void expr::reset_val(){
//        if (is_uexpr()) {
//            return ((uexpr*)this)->reset_val();
//        }
//        else {
//            return ((bexpr*)this)->reset_val();
//        }
//
//    }
    
    
    

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
//        _is_matrix = exp._is_matrix;
        _is_transposed = exp._is_transposed;
        _dim[0] = exp._dim[0]; _dim[1] = exp._dim[1];
        return *this;
    }


    bool uexpr::operator==(const uexpr &c)const{
//        return (_otype == c._otype && equals(_son,c._son));
        return (this->to_str().compare(c.to_str())==0);
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
//        _is_matrix = exp._is_matrix;
        _is_transposed = exp._is_transposed;
        _dim[0] = exp._dim[0]; _dim[1] = exp._dim[1];
        return *this;
    }


    /* UNARY EXPRESSIONS */

    uexpr::uexpr(){
        _type = uexp_c;
    }


    void bexpr::print() const {
        cout << _to_str;
    //    if((_otype==product_ || _otype==div_) && (_lson->get_type()==uexp_c || _lson->get_type()==bexp_c)) {
    //        cout << "(";
    //        poly_print(_lson);
    //        cout << ")";
    //    }
    //    else
    //        poly_print(_lson);
    //
    //    if (_otype==plus_) {
    //        cout << " + ";
    //    }
    //    if (_otype==minus_) {
    //        cout << " - ";
    //    }
    //    if (_otype==product_) {
    //        cout << " * ";
    //    }
    //    if (_otype==div_) {
    //        cout << "/";
    //    }
    //
    //    if (_otype==power_) {
    //        cout << "^";
    //    }
    //
    //    if (_otype==plus_ || (_rson->get_type()!=uexp_c && _rson->get_type()!=bexp_c)) {
    //        poly_print(_rson);
    //    }
    //    else {
    //        cout << "(";
    //        poly_print(_rson);
    //        cout << ")";
    //    }
//        if(endline)
            cout << endl;
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
        _son = make_shared<constant_>(*exp._son);
        _otype = exp._otype;
        _to_str = exp._to_str;
        _coef = exp._coef;
        _is_vector = exp._is_vector;
//        _is_matrix = exp._is_matrix;
        _is_transposed = exp._is_transposed;
        _dim[0] = exp._dim[0]; _dim[1] = exp._dim[1];
        return *this;
    }
    
//    void uexpr::reset_val(){
//        _son->reset_val();
//    }
    
    string uexpr::to_str(int prec) const{
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
        switch (_otype) {
            case log_:
                str += "log(";
                str += _son->to_str(prec);
                str += ")";
                break;
                
            case exp_:
                str += "exp(";
                str += _son->to_str(prec);
                str += ")";
                break;
                
            case cos_:
                str += "cos(";
                str += _son->to_str(prec);
                str += ")";
                break;
                
            case sin_:
                str += "sin(";
                str += _son->to_str(prec);
                str += ")";
                break;
                
            case sqrt_:
                str += "sqrt(";
                str += _son->to_str(prec);
                str += ")";
                break;
            default:
                break;
        }
        if (_coef!=1) {
            str += ")";
        }
        return str;
    }
    
    string uexpr::to_str(size_t inst, int prec) const{
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
        switch (_otype) {
            case log_:
                str += "log(";
                str += _son->to_str(inst,prec);
                str += ")";
                break;
                
            case exp_:
                str += "exp(";
                str += _son->to_str(inst,prec);
                str += ")";
                break;
                
            case cos_:
                str += "cos(";
                str += _son->to_str(inst,prec);
                str += ")";
                break;
                
            case sin_:
                str += "sin(";
                str += _son->to_str(inst,prec);
                str += ")";
                break;
                
            case sqrt_:
                str += "sqrt(";
                str += _son->to_str(inst,prec);
                str += ")";
                break;
            default:
                break;
        }
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
    
    string bexpr::to_str(int prec) const{
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
    
    string bexpr::to_str(size_t inst,int prec) const{
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
//    void bexpr::reset_val(){
//        _lson->reset_val();
//        _rson->reset_val();
//    }
    
    bool bexpr::is_inner_product() const{
        return _otype==product_ && (_lson->get_dim(1)==_rson->get_dim(0) || (_lson->_is_transposed && _lson->get_dim(0)==_rson->get_dim(0)));
    }
    
//    double uexpr::eval(size_t i) const{
//        if (_son->is_constant() && !_son->_evaluated) {//TODO what if son is matrix?
//            for (unsigned inst = 0; inst < _son->_val->size(); inst++) {
//                _son->_val->at(inst) = _son->eval(inst);
//            }
//            _son->_evaluated = true;
//        }
//        double val = 0;
//        if (_son->is_number()) {
//            val = _son->_val->at(0);
//        }
//        else {
//            val = _son->get_val(i);
//        }
//        switch (_otype) {
//            case cos_:
//                return _coef*std::cos(val);
//                break;
//            case sin_:
//                return _coef*std::sin(val);
//                break;
//            case sqrt_:
//                return _coef*std::sqrt(val);
//                break;
//            case log_:
//                return _coef*std::log(val);
//                break;
//            case exp_:
//                return _coef*std::exp(val);
//                break;
//            default:
//                throw invalid_argument("Unsupported unary operator");
//                break;
//        }
//
//    }
//
//    double uexpr::eval(size_t i, size_t j) const{
//        if (!is_matrix()) {
//            return eval(j);//TODO what if son is transposed
//        }
//        if (_son->is_constant() && !_son->_evaluated) {
//            unsigned index = 0;
//            if (_son->is_matrix()) {
//                for (unsigned row = 0; row<_son->_dim[0]; row++) {
//                    for (unsigned col = 0; col<_son->_dim[1]; col++) {
//                        if (_is_transposed) {
//                            index = _son->_dim[0]*col + row;
//                        }
//                        else {
//                            index = _son->_dim[1]*row + col;
//                        }
//
//                        _son->_val->at(index) = _son->eval(row,col);
//                    }
//                }
//            }
//            else {
//                for (size_t row = 0; row<_son->_dim[0]; row++) {
//                    _son->_val->at(index) = _son->eval(index);
//                }
//            }
//            _son->_evaluated = true;
//        }
//        double val = 0;
//        if (_son->is_number()) {
//            val = _son->_val->at(0);
//        }
//        else {
//            val = _son->get_val(i,j);
//        }
//        switch (_otype) {
//            case cos_:
//                return _coef*std::cos(val);
//                break;
//            case sin_:
//                return _coef*std::sin(val);
//                break;
//            case sqrt_:
//                return _coef*std::sqrt(val);
//                break;
//            case log_:
//                return _coef*std::log(val);
//                break;
//            case exp_:
//                return _coef*std::exp(val);
//                break;
//            default:
//                throw invalid_argument("Unsupported unary operator");
//                break;
//        }
//
//    }
    
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
            _dim[0] = max(_dim[0], _lson->_dim[0]);
            _dim[0] = max(_dim[0], _rson->_dim[0]);
            _dim[1] = max(_dim[1], _lson->_dim[1]);
            _dim[1] = max(_dim[1], _rson->_dim[1]);
        }
    };
    
    bexpr& bexpr::operator=(const bexpr& exp){
        _type = bexp_c;
        _lson = make_shared<constant_>(*exp._lson);
        _rson = make_shared<constant_>(*exp._rson);
        _otype = exp._otype;
        _to_str = exp._to_str;
        _coef = exp._coef;
        _is_vector = exp._is_vector;
//        _is_matrix = exp._is_matrix;
        _is_transposed = exp._is_transposed;
        _dim[0] = exp._dim[0]; _dim[1] = exp._dim[1];
        return *this;
    }
    
    
    
//    double  bexpr::eval(size_t i, size_t j) const{
//
//        switch (_otype) {
//            case plus_:
//                return _coef*(_lson->get_val(i,j) + _rson->get_val(i,j));
//                break;
//            case minus_:
//                return _coef*(_lson->get_val(i,j) - _rson->get_val(i,j));
//                break;
//            case product_:{
//                double res = 0;
//                if (_lson->is_matrix() && _rson->is_matrix()) {
//                    //matrix product
//                    for (unsigned col = 0; col<_lson->_dim[1]; col++) {
//                        res += _lson->get_val(i,col) * _rson->get_val(col,j);
//                    }
//                    return _coef*res;
//                }
//                if (_lson->is_matrix() && !_rson->is_matrix() && _rson->_is_transposed) {//matrix * transposed vect
//                    return _coef*(_lson->get_val(i,j) * _rson->get_val(j));
//                }
//                if (!_lson->is_matrix() && !_lson->_is_transposed && _rson->is_matrix() ) {//vect * matrix
//                    return _coef*(_lson->get_val(i) * _rson->get_val(i,j));//TODO i ot j?
//                }
//                if (_lson->is_matrix() && _rson->_is_vector && _rson->_is_transposed) {//matrix * vector
//                    return _coef*(_lson->get_val(i,j) * _rson->get_val(j));
//                }
//                throw invalid_argument("eval(i,j) on non matrix function");
//                break;
//            }
//            case div_:
//                return _coef*(_lson->get_val(i,j)/_rson->get_val(i,j));
//                break;
//            case power_:
//                return _coef*(powl(_lson->get_val(i,j),_rson->get_val(i,j)));
//                break;
//            default:
//                throw invalid_argument("Unsupported binary operator");
//                break;
//        }
//
//    }
//    double  bexpr::eval(size_t i) const{
//        //        if (_lson->is_constant() && !_lson->_evaluated) {
//        //            _lson->_val->resize(_lson->get_dim());//TODO what if son is a matrix?
//        ////            _lson->_val = make_shared<vector<double>>();
//        ////            if (_lson->_is_transposed) {
//        ////                _lson->_val->resize(_lson->_dim[0]);
//        ////            }
//        //            for (unsigned inst = 0; inst < _lson->_val->size(); inst++) {
//        //                _lson->_val->at(inst) = _lson->eval(inst);
//        //            }
//        //            _lson->_evaluated = true;
//        //        }
//        //        if (_rson->is_constant() && !_rson->_evaluated) {
//        //            _rson->_val->resize(_rson->get_dim());
//        ////            _rson->_val = make_shared<vector<double>>();
//        ////            _rson->_val->resize(_rson->get_dim());
//        ////            if (_rson->_is_transposed) {
//        ////                _rson->_val->resize(_rson->_dim[0]);
//        ////            }
//        //            for (unsigned inst = 0; inst < _rson->_val->size(); inst++) {
//        //                _rson->_val->at(inst) = _rson->eval(inst);
//        //            }
//        //            _rson->_evaluated = true;
//        //        }
//        double lval = 0, rval = 0;
//        if (_lson->is_number()) {
//            if (_lson->_val->size()==0) {
//                lval = _lson->eval(0);
//            }
//            else {
//                lval = _lson->_val->at(0);
//            }
//        }
//        else if(!_lson->_is_vector){
//            if (_lson->_val->size()<=i) {
//                lval = _lson->eval(i);
//            }
//            else {
//                lval = _lson->get_val(i);
//            }
//        }
//        if (_rson->is_number()) {
//            if (_rson->_val->size()==0) {
//                rval = _rson->eval(0);
//            }
//            else {
//                rval = _rson->_val->at(0);
//            }
//        }
//        else if(!_rson->_is_vector){
//            if (_rson->_val->size()<=i) {
//                rval = _rson->eval(i);
//            }
//            else {
//                rval = _rson->get_val(i);
//            }
//        }
//        switch (_otype) {
//            case plus_:
//                return _coef*(lval + rval);
//                break;
//            case minus_:
//                return _coef*(lval - rval);
//                break;
//            case product_:
//                if (_lson->is_matrix() && !_rson->is_matrix() && !_rson->_is_transposed) {//matrix * vect
//                    double res = 0;
//                    for (size_t j = 0; j<_rson->_dim[0]; j++) {
//                        res += _lson->get_val(i,j) * _rson->get_val(j);
//                    }
//                    return _coef*(res);
//                }
//                if (!_lson->is_matrix() && _lson->_is_transposed && _rson->is_matrix() ) {//transposed vect * matrix
//                    double res = 0;
//                    for (size_t j = 0; j<_lson->_dim[0]; j++) {
//                        res += _lson->get_val(j) * _rson->get_val(j,i);
//                    }
//                    return _coef*(res);
//                }
//                if (!_lson->is_matrix() && _lson->_is_transposed && !_rson->is_matrix() && i==0) {//transposed vect * vec, a dot product of two vectors
//                    double res = 0;
//                    for (size_t j = 0; j<_lson->_dim[1]; j++) {
//                        res += _lson->get_val(j) * _rson->get_val(j);
//                    }
//                    return _coef*(res);
//                }
//                return _coef*(lval*rval);
//                break;
//            case div_:
//                return _coef*(lval/rval);
//                break;
//            case power_:
//                return _coef*(powl(lval,rval));
//                break;
//            default:
//                throw invalid_argument("Unsupported binary operator");
//                break;
//        }
//
//    }
}
