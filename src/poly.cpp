//
//  poly.cpp
//  Gravity
//
//  Created by Hijazi, Hassan on 25 Oct 18.
//
//
#include <cmath>
#include <gravity/poly.h>
#include <gravity/func.h>

using namespace std;
namespace gravity{        

    lterm::lterm(bool sign, constant_* coef, param_* p){
        _coef = coef;
        _p = p;
        _sign = sign;
        if (coef->_is_transposed && p->_is_transposed) {
            throw invalid_argument("in lterm(bool sign, constant_* coef, param_* p), both coef and p are transposed!\n");
        }
    };


    lterm& lterm::operator=(const lterm &l){
        if (_coef) {
            delete _coef;
        }
        if (_p) {
            delete _p;
        }
        _coef = copy(*l._coef);
        _p = (param_*)copy(*l._p);
        _sign = l._sign;
        return *this;
    }

    lterm& lterm::operator=(lterm&& l){
        if (_coef) {
            delete _coef;
        }
        if (_p) {
            delete _p;
        }
        _coef = l._coef;
        l._coef = nullptr;
        _p = l._p;
        l._p = nullptr;
        _sign = l._sign;
        return *this;
    }

    qterm& qterm::operator=(const qterm &q){
        if (_coef) {
            delete _coef;
        }
        if (_p) {
            delete _p;
        }
        _coef = copy(*q._coef);
        _p = new pair<param_*, param_*>(make_pair<>((param_*)copy(*q._p->first), (param_*)copy(*q._p->second)));
        _sign = q._sign;
        return *this;
    }


    qterm& qterm::operator=(qterm&& q){
        if (_coef) {
            delete _coef;
        }
        if (_p) {
            delete _p;
        }
        _coef = q._coef;
        q._coef = nullptr;
        _p = q._p;
        q._p = nullptr;
        _sign = q._sign;
        return *this;
    }



    pterm& pterm::operator=(const pterm &p){
        if (_coef) {
            delete _coef;
        }
        if (_l) {
            delete _l;
        }
        _coef = copy(*p._coef);
        _l = new list<pair<param_*, int>>();
        for (auto &pair : *p._l) {
            _l->push_back(make_pair<>((param_*)copy(*pair.first), pair.second));
        }
        _sign = p._sign;
        return *this;
    }


    pterm& pterm::operator=(pterm&& p){
        if (_coef) {
            delete _coef;
        }
        if (_l) {
            delete _l;
        }
        _coef = p._coef;
        p._coef = nullptr;
        _l = p._l;
        p._l = nullptr;
        _sign = p._sign;
        return *this;
    }

    
    double lterm::eval(size_t i) const{
        double res = 0;
//        if (_p->is_indexed() && _p->_indices->size()>1 && _p->_ids->at(i).size()==0) {
//            return 0;
//        }
        if ((_coef->_is_transposed || _coef->is_matrix() || (_p->is_indexed() && _p->_indices->_ids->size()>1)) && !_p->is_matrix()) {
            auto dim = _p->get_dim(i);
            for (size_t j = 0; j<dim; j++) {
                res += t_eval(_coef,i,j) * t_eval(_p,i,j);
            }
        }
        else {
            res = t_eval(_coef,i) * t_eval(_p, i);
        }
        if (!_sign) {
            res *= -1;
        }
        return res;
    }

    double lterm::eval(size_t i, size_t j) const{
        if (_coef->is_matrix() && _p->is_matrix()) {
            //matrix product
                    auto res = 0.;
                    for (size_t col = 0; col<_coef->_dim[1]; col++) {
                        res += t_eval(_coef, i,col) * t_eval(_p,col,j);
                    }
            if(_sign){
                return res;
            }
            else {
                return -1*res;
            }
        }
        if (_coef->is_matrix() && !_p->is_matrix() && _p->_is_transposed) {//matrix * transposed vect
            if(_sign){
                return t_eval(_coef,i,j) * t_eval(_p,j);
            }
            else {
                return -1*t_eval(_coef,i,j) * t_eval(_p,j);
            }
        }
        
        if (!_coef->is_matrix() && !_coef->_is_transposed && _p->is_matrix()) {//vect * matrix
            if(_sign) {
                return t_eval(_coef,i) * t_eval(_p,i,j);
            }
            else {
                return -1 * t_eval(_coef,i) * t_eval(_p,i,j);
            }
        }
        if (_coef->is_matrix() && _p->_is_vector) {//matrix*vect
            if(_sign) {
                return t_eval(_coef,i,j) * t_eval(_p,i);
            }
            else {
                return -1 * t_eval(_coef,i,j) * t_eval(_p,i);
            }
        }
        double res = t_eval(_coef,i,j) * t_eval(_p,i,j);
        if (!_sign) {
            res *= -1;
        }
        return res;
    }


    double qterm::eval(size_t i) const{
        double res = 0;
        if (_p->first->is_matrix() && !_p->second->is_matrix() && !_p->second->_is_transposed) {//matrix * vect            
            for (size_t j = 0; j<_p->second->_dim[0]; j++) {
                res += t_eval(_p->first,i,j) * t_eval(_p->second,j);
            }
            res *= t_eval(_coef,i);
        }
        else if (!_p->first->is_matrix() && _p->first->_is_transposed && _p->second->is_matrix() ) {//transposed vect * matrix
            for (size_t j = 0; j<_p->first->_dim[0]; j++) {
                res += t_eval(_p->first,j) * t_eval(_p->second,j,i);
            }
            res *= t_eval(_coef,i);
        }
        else if (!_p->first->is_matrix() && _p->first->_is_transposed && !_p->second->is_matrix() && i==0) {//transposed vect * vec, a dot product of two vectors
            for (size_t j = 0; j<_p->first->_dim[1]; j++) {
                res += t_eval(_p->first,j) * t_eval(_p->second,j);
            }
            res *= t_eval(_coef,i);
        }
        else if (!_coef->is_matrix() && _coef->_is_transposed && !_p->first->is_matrix()) {//transposed vect * vec, a dot product of two vectors
            for (size_t j = 0; j<_p->first->_dim[0]; j++) {
                res += t_eval(_coef,j) * t_eval(_p->first,j) * t_eval(_p->second,j);
            }
        }
        else {
            res = t_eval(_coef,i)*t_eval(_p->first,i)*t_eval(_p->second,i);
        }
        if (!_sign) {
            res *= -1;
        }
        return res;
    }

    double qterm::eval(size_t i, size_t j) const{
        double res = 0;
        if (_p->first->is_matrix() && _p->second->is_matrix()) {
            //matrix product
            for (size_t col = 0; col<_p->first->_dim[1]; col++) {
                res += t_eval(_p->first,i,col) * t_eval(_p->second,col,j);
            }
            res *= t_eval(_coef,i);
        }
        else if (_p->first->is_matrix() && !_p->second->is_matrix() && _p->second->_is_transposed) {//matrix * transposed vect
            res = t_eval(_coef,i)*(t_eval(_p->first,i,j) * t_eval(_p->second,j));
        }
        else if (!_p->first->is_matrix() && !_p->first->_is_transposed && _p->second->is_matrix() ) {//vect * matrix
            res = t_eval(_coef,i)*(t_eval(_p->first,i) * t_eval(_p->second,i,j));
        }
        else {
            throw invalid_argument("eval(i,j) on non-matrix function");
        }
        if (!_sign) {
            res *= -1;
        }
        return res;
    }


    double pterm::eval(size_t i) const{
        double res = 0;
        if (_coef->_is_transposed) {
            throw invalid_argument("Unspported operation\n");
        }// TREAT TRANSPOSED VECTORS IN POLYNOMIAL TERMS HERE
        else {
            res =1;
            for (auto &pair: *_l) {
                res *= pow(t_eval(pair.first, i), pair.second);
            }
            res *= t_eval(_coef,i);
        }
        if (!_sign) {
            res *= -1;
        }
        return res;
    }

    double pterm::eval(size_t i, size_t j) const{
        double res = 0;
        if (_coef->_is_transposed) {
            throw invalid_argument("Unspported operation\n");
        }// TREAT TRANSPOSED VECTORS IN POLYNOMIAL TERMS HERE
        else {
            res =1;
            for (auto &pair: *_l) {
                res *= pow(t_eval(pair.first,i,j), pair.second);
            }
            
            res *= t_eval(_coef,i,j);
        }
        if (!_sign) {
            res *= -1;
        }
        return res;
    }


        
    string pterm::to_str(size_t ind) const{
        string str;
        constant_* c_new = _coef;
        if (c_new->is_number()){
            string v = poly_to_str(c_new);
            if (_sign) {
                if (v=="-1") {
                    str += " - ";
                }
                else if (ind>0) {
                    str += " + ";
                    if(v!="1") {
                        str += v;
                    }
                }
                else if(v!="1") {
                    str += v;
                }
            }
            if(!_sign) {
                if (v == "-1" && ind>0) {
                    str += " + ";
                }
                else if (v.front()=='-'){
                    if (ind > 0) {
                        str += " + ";
                    }
                    str += v.substr(1);
                }
                else if (v=="1"){
                    str += " - ";
                }
                else if(v!="-1"){
                    str += " - " + v;
                }
            }
        }
        else{
            if (!_sign) {
                str += " - ";
            }
            if(ind > 0 && _sign) {
                str += " + ";
            }
            str += "(";
            str += poly_to_str(c_new);
            str += ")";
        }
        for (auto& p: *_l) {
            str += poly_to_str(p.first);
            if (p.second != 1) {
                switch (p.second) {
                    case 2:
                        str += "\u00B2";
                        break;
                    case 3:
                        str += "\u00B3";
                        break;
                    case 4:
                        str += "\u00B4";
                        break;
                    case 5:
                        str += "\u00B5";
                        break;
                    case 6:
                        str += "\u00B6";
                        break;
                    case 7:
                        str += "\u00B7";
                        break;
                    case 8:
                        str += "\u00B8";
                        break;
                    case 9:
                        str += "\u00B9";
                        break;
                    default:
                        str += "^" + to_string(p.second);
                        break;
                }
            }
        }
        return str;
    }
    
    string pterm::to_str(size_t ind, size_t inst) const{
        string str;
        constant_* c_new = _coef;
        if (c_new->is_number()){
            string v = poly_to_str(c_new);
            if (_sign) {
                if (v=="-1") {
                    str += " - ";
                }
                else if (ind>0) {
                    str += " + ";
                    if(v!="1") {
                        str += v;
                    }
                }
                else if(v!="1") {
                    str += v;
                }
            }
            if(!_sign) {
                if (v == "-1" && ind>0) {
                    str += " + ";
                }
                else if (v.front()=='-'){
                    if (ind > 0) {
                        str += " + ";
                    }
                    str += v.substr(1);
                }
                else if (v=="1"){
                    str += " - ";
                }
                else if(v!="-1"){
                    str += " - " + v;
                }
            }
        }
        else{
            if (!_sign) {
                str += " - ";
            }
            if(ind > 0 && _sign) {
                str += " + ";
            }
            str += "(";
            str += poly_to_str(c_new, inst);
            str += ")";
        }
        for (auto& p: *_l) {
            str += poly_to_str(p.first, inst);
            if (p.second != 1) {
                switch (p.second) {
                    case 2:
                        str += "\u00B2";
                        break;
                    case 3:
                        str += "\u00B3";
                        break;
                    case 4:
                        str += "\u00B4";
                        break;
                    case 5:
                        str += "\u00B5";
                        break;
                    case 6:
                        str += "\u00B6";
                        break;
                    case 7:
                        str += "\u00B7";
                        break;
                    case 8:
                        str += "\u00B8";
                        break;
                    case 9:
                        str += "\u00B9";
                        break;
                    default:
                        str += "^" + to_string(p.second);
                        break;
                }
            }
        }
        return str;
    }

    void pterm::print(size_t ind) const{
        cout << this->to_str(ind);
    }

    string qterm::to_str(size_t ind) const {
        string str;
        constant_* c_new = _coef;
        param_* p_new1 = (param_*)_p->first;
        param_* p_new2 = (param_*)_p->second;
        if (c_new->is_number()){
            string v = poly_to_str(c_new);
            if (_sign) {
                if (v=="-1") {
                    str += " - ";
                }
                else if (ind>0) {
                    str += " + ";
                    if(v!="1") {
                        str += v;
                    }
                }
                else if(v!="1") {
                    str += v;
                }
            }
            if(!_sign) {
                if (v == "-1" && ind>0) {
                    str += " + ";
                }
                else if (v.front()=='-'){
                    if (ind > 0) {
                        str += " + ";
                    }
                    str += v.substr(1);
                }
                else if (v=="1"){
                    str += " - ";
                }
                else if(v!="-1"){
                    str += " - " + v;
                }
            }
        }
        else{
            if (!_sign) {
                str += " - ";
            }
            if(ind > 0 && _sign) {
                str += " + ";
            }
            str += "(";
            str += poly_to_str(c_new);
            str += ")";
        }
        str += poly_to_str(p_new1);
        if (p_new1==p_new2) {
            str += "²";
        }
        else {
            str += "*";
            str += poly_to_str(p_new2);
        }
        return str;
    }
    
    string qterm::to_str(size_t ind, size_t inst) const {
        string str;
        constant_* c_new = _coef;
        param_* p_new1 = (param_*)_p->first;
        param_* p_new2 = (param_*)_p->second;
        unsigned dim = 1;
        if (c_new->_is_transposed) {
            dim = p_new1->get_dim(inst);
        }
        for (unsigned idx = 0; idx <dim; idx++) {
            if (c_new->is_number()){
                string v;
                if (c_new->_is_transposed) {
                    v = poly_to_str(c_new,idx);
                }
                else{
                    v = poly_to_str(c_new,inst);
                }
                
                if (_sign) {
                    if (v=="-1") {
                        str += " - ";
                    }
                    else if (idx>0 || ind>0) {
                        str += " + ";
                        if(v!="1") {
                            str += v;
                        }
                    }
                    else if(v!="1") {
                        str += v;
                    }
                }
                if(!_sign) {
                    if (v == "-1" && (idx>0 || ind>0)) {
                        str += " + ";
                    }
                    else if (v.front()=='-'){
                        if (ind > 0) {
                            str += " + ";
                        }
                        str += v.substr(1);
                    }
                    else if (v=="1"){
                        str += " - ";
                    }
                    else if(v!="-1"){
                        str += " - " + v;
                    }
                }
            }
            else{
                if (!_sign) {
                    str += " - ";
                }
                if((ind > 0 || idx >0) && _sign) {
                    str += " + ";
                }
                str += "(";
                if (c_new->_is_transposed) {
                    str += poly_to_str(c_new, idx);
                }
                else{
                    str += poly_to_str(c_new, inst);
                }
                str += ")";
            }
            if (c_new->_is_transposed) {
                str += poly_to_str(p_new1, idx);
            }
            else{
                str += poly_to_str(p_new1, inst);
            }
        
            if (p_new1==p_new2) {
                str += "²";
            }
            else {
                str += "*";
                if (c_new->_is_transposed) {
                    str += poly_to_str(p_new2, idx);
                }
                else{
                    str += poly_to_str(p_new2, inst);
                }
            }
        }
        return str;
    }


    void qterm::print(size_t ind) const {
        cout << this->to_str(ind);
    }


    string lterm::to_str(size_t ind) const{
        string str;
        constant_* c_new = _coef;
        param_* p_new = (param_*)_p;
        if (c_new->is_number()){
            string v = poly_to_str(c_new);
            if (_sign) {
                if (v=="-1") {
                    str += " - ";
                }
                else if (ind>0) {
                    str += " + ";
                    if(v!="1") {
                        str += v;
                    }
                }
                else if(v!="1") {
                    str += v;
                }
            }
            if(!_sign) {
                if (v == "-1" && ind>0) {
                    str += " + ";
                }
                else if (v.front()=='-'){
                    if (ind > 0) {
                        str += " + ";
                    }
                    str += v.substr(1);
                }
                else if (v=="1"){
                    str += " - ";
                }
                else if(v!="-1"){
                    str += " - " + v;
                }
            }
        }
        else{
            if (!_sign) {
                str += " - ";
            }
            if(ind > 0 && _sign) {
                str += " + ";
            }
            str += "(";
            str += poly_to_str(c_new);
            str += ")";
        }
        str += poly_to_str(p_new);
        return str;
    }


    void lterm::print(size_t ind) const{
        cout << this->to_str(ind);
    }

    string lterm::to_str(size_t ind, size_t inst) const{
        string str;
        constant_* c_new = _coef;
        param_* p_new = (param_*)_p;        
        unsigned dim = 1;
        if (c_new->_is_transposed) {
            dim = p_new->get_dim(inst);
        }
        for (unsigned idx = 0; idx <dim; idx++) {
            if (c_new->constant_::is_number()){
                string v;
                v = poly_to_str(c_new);
                if (_sign) {
                    if (v=="-1") {
                        str += " - ";
                    }
                    else if (idx>0 || ind>0) {
                        str += " + ";
                        if(v!="1") {
                            str += v;
                        }
                    }
                    else if(v!="1") {
                        str += v;
                    }
                }
                if(!_sign) {
                    if (v == "-1" && (idx>0 || ind>0)) {
                        str += " + ";
                    }
                    else if (v.front()=='-'){
                        if (ind > 0 || idx >0) {
                            str += " + ";
                        }
                        str += v.substr(1);
                    }
                    else if (v=="1"){
                        str += " - ";
                    }
                    else if(v!="-1"){
                        str += " - " + v;
                    }
                }
            }
            else{
                
                if (!_sign) {
                    str += " - ";
                }
                if((ind > 0 || idx >0) && _sign) {
                    str += " + ";
                }
                if (!c_new->is_function() || !((func_*)c_new)->is_unit_instance()) {
                    str += "(";
                    if (c_new->is_matrix()) {
                        if (c_new->is_function() && !((func_*)(c_new))->is_constant()){
                            str += ((func_*)(c_new))->to_str(inst);
                        }
                        else {
                            str += "[\n";
                            for (size_t i = 0; i<c_new->_dim[0]; i++) {
                                for (size_t j = 0; j<c_new->_dim[1]; j++) {
                                    str += poly_to_str(c_new, i, j)+ " ";
                                }
                                str += "\n";
                            }
                            str += "\n] * ";
                        }
                    }
                    else if (c_new->_is_transposed) {
                        str += poly_to_str(c_new, inst, idx);
                    }
                    else{
                        str += poly_to_str(c_new, inst);
                    }
                    str += ")";
                }
            }
            if (p_new->is_matrix()) {
                str += "[\n";
                for (size_t i = 0; i<c_new->_dim[0]; i++) {
                    for (size_t j = 0; j<c_new->_dim[1]; j++) {
                        str += poly_to_str(p_new, i, j)+ " ";
                    }
                    str += "\n";
                }
                str += "\n]";
            }
            else if (c_new->_is_transposed || (p_new->is_indexed() && p_new->_indices->_ids->size()>1)) {
                str += poly_to_str(p_new, inst, idx);
            }
            else{
                str += poly_to_str(p_new, inst);
            }
        }
        return str;
    }

}
