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

    lterm::lterm(bool sign, shared_ptr<constant_> coef, shared_ptr<param_> p){
        _coef = coef;
        _p = p;
        _sign = sign;
        if (coef->_is_transposed && p->_is_transposed) {
            throw invalid_argument("in lterm(bool sign, constant_* coef, param_* p), both coef and p are transposed!\n");
        }
    };


    lterm& lterm::operator=(const lterm &l){
        _coef = l._coef->copy();
        _p = l._p->pcopy();
        _sign = l._sign;
        return *this;
    }

    lterm& lterm::operator=(lterm&& l){
        _coef = move(l._coef);
        _p = move(l._p);
        _sign = l._sign;
        return *this;
    }

    qterm& qterm::operator=(const qterm &q){
        _coef = q._coef->copy();
        _p = make_shared<pair<shared_ptr<param_>,shared_ptr<param_>>>(make_pair<>(q._p->first->pcopy(), q._p->second->pcopy()));
        _sign = q._sign;
        _coef_p1_transposed = q._coef_p1_transposed;
        return *this;
    }


    qterm& qterm::operator=(qterm&& q){
        _coef = move(q._coef);
        _p = move(q._p);
        _sign = q._sign;
        _coef_p1_transposed = q._coef_p1_transposed;
        return *this;
    }



    pterm& pterm::operator=(const pterm &p){
        _coef = p._coef->copy();
        _l = make_shared<list<pair<shared_ptr<param_>, int>>>();
        for (auto &pair : *p._l) {
            _l->push_back(make_pair<>(pair.first->pcopy(), pair.second));
        }
        _sign = p._sign;
        return *this;
    }


    pterm& pterm::operator=(pterm&& p){
        _coef = move(p._coef);
        _l = move(p._l);
        _sign = p._sign;
        return *this;
    }

    
//    double lterm::eval(size_t i) const{
//        double res = 0;
////        if (_p->is_indexed() && _p->_indices->size()>1 && _p->_ids->at(i).size()==0) {
////            return 0;
////        }
//        if ((_coef->_is_transposed || _coef->is_matrix() || (_p->is_indexed() && _p->_indices->_ids->size()>1)) && !_p->is_matrix()) {
//            auto dim = _p->get_dim(i);
//            for (size_t j = 0; j<dim; j++) {
//                res += t_eval(_coef,i,j) * t_eval(_p,i,j);
//            }
//        }
//        else {
//            res = t_eval(_coef,i) * t_eval(_p, i);
//        }
//        if (!_sign) {
//            res *= -1;
//        }
//        return res;
//    }
//
//    double lterm::eval(size_t i, size_t j) const{
//        if (_coef->is_matrix() && _p->is_matrix()) {
//            //matrix product
//                    auto res = 0.;
//                    for (size_t col = 0; col<_coef->_dim[1]; col++) {
//                        res += t_eval(_coef, i,col) * t_eval(_p,col,j);
//                    }
//            if(_sign){
//                return res;
//            }
//            else {
//                return -1*res;
//            }
//        }
//        if (_coef->is_matrix() && !_p->is_matrix() && _p->_is_transposed) {//matrix * transposed vect
//            if(_sign){
//                return t_eval(_coef,i,j) * t_eval(_p,j);
//            }
//            else {
//                return -1*t_eval(_coef,i,j) * t_eval(_p,j);
//            }
//        }
//
//        if (!_coef->is_matrix() && !_coef->_is_transposed && _p->is_matrix()) {//vect * matrix
//            if(_sign) {
//                return t_eval(_coef,i) * t_eval(_p,i,j);
//            }
//            else {
//                return -1 * t_eval(_coef,i) * t_eval(_p,i,j);
//            }
//        }
//        if (_coef->is_matrix() && _p->_is_vector) {//matrix*vect
//            if(_sign) {
//                return t_eval(_coef,i,j) * t_eval(_p,i);
//            }
//            else {
//                return -1 * t_eval(_coef,i,j) * t_eval(_p,i);
//            }
//        }
//        double res = t_eval(_coef,i,j) * t_eval(_p,i,j);
//        if (!_sign) {
//            res *= -1;
//        }
//        return res;
//    }
//
//
//    double qterm::eval(size_t i) const{
//        double res = 0;
//        if (_coef_p1_transposed) { // qterm = (coef*p1)^T*p2
//            assert(_p->first->_dim[1]==1 && _coef->_dim[0]==_p->second->_dim[0]);
//            for (auto i = 0; i<_p->first->_dim[0]; i++) {
//                for (auto j = 0; j<_p->first->_dim[0]; j++) {
//                    res += t_eval(_coef,i,j)* t_eval(_p->first,i) * t_eval(_p->second,j);
//                }
//            }
//            if (!_sign) {
//                res *= -1;
//            }
//            return res;
//
//        }
//        if (_p->first->is_matrix() && !_p->second->is_matrix() && !_p->second->_is_transposed) {//matrix * vect
//            for (size_t j = 0; j<_p->second->_dim[0]; j++) {
//                res += t_eval(_p->first,i,j) * t_eval(_p->second,j);
//            }
//            res *= t_eval(_coef,i);
//        }
//        else if (!_p->first->is_matrix() && _p->first->_is_transposed && _p->second->is_matrix() ) {//transposed vect * matrix
//            for (size_t j = 0; j<_p->first->_dim[0]; j++) {
//                res += t_eval(_p->first,j) * t_eval(_p->second,j,i);
//            }
//            res *= t_eval(_coef,i);
//        }
//        else if (!_p->first->is_matrix() && _p->first->_is_transposed && !_p->second->is_matrix() && i==0) {//transposed vect * vec, a dot product of two vectors
//            for (size_t j = 0; j<_p->first->_dim[1]; j++) {
//                res += t_eval(_p->first,j) * t_eval(_p->second,j);
//            }
//            res *= t_eval(_coef,i);
//        }
//        else if (!_coef->is_matrix() && _coef->_is_transposed && !_p->first->is_matrix()) {//transposed vect * vec, a dot product of two vectors
//            for (size_t j = 0; j<_p->first->_dim[0]; j++) {
//                res += t_eval(_coef,j) * t_eval(_p->first,j) * t_eval(_p->second,j);
//            }
//        }
//        else {
//            res = t_eval(_coef,i)*t_eval(_p->first,i)*t_eval(_p->second,i);
//        }
//        if (!_sign) {
//            res *= -1;
//        }
//        return res;
//    }
//
//    double qterm::eval(size_t i, size_t j) const{
//        double res = 0;
//        if (_p->first->is_matrix() && _p->second->is_matrix()) {
//            //matrix product
//            for (size_t col = 0; col<_p->first->_dim[1]; col++) {
//                res += t_eval(_p->first,i,col) * t_eval(_p->second,col,j);
//            }
//            res *= t_eval(_coef,i);
//        }
//        else if (_p->first->is_matrix() && !_p->second->is_matrix() && _p->second->_is_transposed) {//matrix * transposed vect
//            res = t_eval(_coef,i)*(t_eval(_p->first,i,j) * t_eval(_p->second,j));
//        }
//        else if (!_p->first->is_matrix() && !_p->first->_is_transposed && _p->second->is_matrix() ) {//vect * matrix
//            res = t_eval(_coef,i)*(t_eval(_p->first,i) * t_eval(_p->second,i,j));
//        }
//        else {
//            throw invalid_argument("eval(i,j) on non-matrix function");
//        }
//        if (!_sign) {
//            res *= -1;
//        }
//        return res;
//    }
//
//
//    double pterm::eval(size_t i) const{
//        double res = 0;
//        if (_coef->_is_transposed) {
//            throw invalid_argument("Unspported operation\n");
//        }// TREAT TRANSPOSED VECTORS IN POLYNOMIAL TERMS HERE
//        else {
//            res =1;
//            for (auto &pair: *_l) {
//                res *= pow(t_eval(pair.first, i), pair.second);
//            }
//            res *= t_eval(_coef,i);
//        }
//        if (!_sign) {
//            res *= -1;
//        }
//        return res;
//    }
//
//    double pterm::eval(size_t i, size_t j) const{
//        double res = 0;
//        if (_coef->_is_transposed) {
//            throw invalid_argument("Unspported operation\n");
//        }// TREAT TRANSPOSED VECTORS IN POLYNOMIAL TERMS HERE
//        else {
//            res =1;
//            for (auto &pair: *_l) {
//                res *= pow(t_eval(pair.first,i,j), pair.second);
//            }
//
//            res *= t_eval(_coef,i,j);
//        }
//        if (!_sign) {
//            res *= -1;
//        }
//        return res;
//    }
//
//
//
    string pterm::to_str() const{
        string str;
        auto c_new = _coef;
        if (c_new->is_number()){
            string v = c_new->to_str();
            if (_sign) {
                if (v=="-1") {
                    str += " - ";
                }
                else if(v!="1") {
                    str += v;
                }
            }
            if(!_sign) {
                if (v == "-1") {
                    str += " + ";
                }
                else if (v.front()=='-'){
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
            if(_sign) {
                str += " + ";
            }
            str += "(";
            str += c_new->to_str();
            str += ")";
        }
        for (auto& p: *_l) {
            str += p.first->get_name(true,true);
            if (p.second != 1) {
                switch (p.second) {
                    case 2:
                        str += "²";
                        break;
                    case 3:
                        str += "³";
                        break;
                    case 4:
                        str += "⁴";
                        break;
                    case 5:
                        str += "⁵";
                        break;
                    case 6:
                        str += "⁶";
                        break;
                    case 7:
                        str += "⁷";
                        break;
                    case 8:
                        str += "⁸";
                        break;
                    case 9:
                        str += "⁹";
                        break;
                    default:
                        str += "^" + to_string(p.second);
                        break;
                }
            }
        }
        return str;
    }
    string pterm::to_str(size_t ind, int prec) const{
        string str;
        auto c_new = _coef;
        if (c_new->is_number()){
            string v = c_new->to_str(prec);
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
            str += c_new->to_str(ind,prec);
            str += ")";
        }
        for (auto& p: *_l) {
            str += p.first->get_name(ind);
            if (p.second != 1) {
                switch (p.second) {
                    case 2:
                        str += "²";
                        break;
                    case 3:
                        str += "³";
                        break;
                    case 4:
                        str += "⁴";
                        break;
                    case 5:
                        str += "⁵";
                        break;
                    case 6:
                        str += "⁶";
                        break;
                    case 7:
                        str += "⁷";
                        break;
                    case 8:
                        str += "⁸";
                        break;
                    case 9:
                        str += "⁹";
                        break;
                    default:
                        str += "^" + to_string(p.second);
                        break;
                }
            }
        }
        return str;
    }
    
    string pterm::to_str(size_t ind, size_t inst, int prec) const{
        string str;
        auto c_new = _coef;
        if (c_new->is_number()){
            string v = c_new->to_str(inst,prec);
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
            str += c_new->to_str(inst,prec);
            str += ")";
        }
        for (auto& p: *_l) {
            str += p.first->get_name(inst);
            if (p.second != 1) {
                switch (p.second) {
                    case 2:
                        str += "²";
                        break;
                    case 3:
                        str += "³";
                        break;
                    case 4:
                        str += "⁴";
                        break;
                    case 5:
                        str += "⁵";
                        break;
                    case 6:
                        str += "⁶";
                        break;
                    case 7:
                        str += "⁷";
                        break;
                    case 8:
                        str += "⁸";
                        break;
                    case 9:
                        str += "⁹";
                        break;
                    default:
                        str += "^" + to_string(p.second);
                        break;
                }
            }
        }
        return str;
    }

//    void pterm::print(size_t ind) const{
//        cout << this->to_str(ind);
//    }

    string qterm::to_str() const {
        string str;
        auto c_new = _coef;
        auto p_new1 = _p->first;
        auto p_new2 = _p->second;
        if (c_new->is_number()){
            string v = c_new->to_str();
            if (_sign) {
                if (v=="-1") {
                    str += " - ";
                }
                else if(v!="1") {
                    str += v;
                }
            }
            if(!_sign) {
                if (v == "-1") {
                    str += " + ";
                }
                else if (v.front()=='-'){
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
            if(_sign) {
                str += " + ";
            }
            if(_coef_p1_transposed){
                str += "(";
            }
            str += "(";
            str += c_new->to_str();
            str += ")";
        }
        str += p_new1->get_name(true,true);
        if(_coef_p1_transposed){
            str += ")\u1D40";
        }
        else if (p_new1==p_new2) {
            str += "²";
            return str;
        }
        str += p_new2->get_name(true,true);
        return str;
    }
    
    string qterm::to_str(size_t ind, int prec) const {
        string str;
        auto c_new = _coef;
        auto p_new1 = _p->first;
        auto p_new2 = _p->second;
        if (c_new->is_number()){
            string v = c_new->to_str(prec);
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
            if(_coef_p1_transposed){
                str += "(";
            }
            str += "(";
            str += c_new->to_str(prec);
            str += ")";
        }
        str += p_new1->get_name(ind);
        if(_coef_p1_transposed){
            str += ")\u1D40";
        }
        else if (p_new1==p_new2) {
            str += "²";
            return str;
        }
//            str += ".";
        str += p_new2->get_name(ind);
        return str;
    }
    
    string qterm::to_str(size_t ind, size_t inst, int prec) const {
        string str;
        auto c_new = _coef;
        auto p_new1 = _p->first;
        auto p_new2 = _p->second;
        unsigned dim = 1;
        if (c_new->_is_transposed) {
            dim = p_new1->get_dim(inst);
        }
        for (auto idx = 0; idx <dim; idx++) {
            if (c_new->is_number()){
                string v;
                if (c_new->_is_transposed) {
                    v = c_new->to_str(idx,prec);
                }
                else{
                    v = c_new->to_str(inst,prec);
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
                    str += c_new->to_str(idx);
                }
                else{
                    str += c_new->to_str(inst);
                }
                str += ")";
            }
            if (c_new->_is_transposed) {
                str += p_new1->to_str(idx);
            }
            else{
                str += p_new1->to_str(inst);
            }
        
            if (p_new1==p_new2) {
                str += "²";
            }
            else {
                str += "*";
                if (c_new->_is_transposed) {
                    str += p_new2->get_name(idx);
                }
                else{
                    str += p_new2->get_name(inst);
                }
            }
        }
        return str;
    }


//    void qterm::print(size_t ind) const {
//        cout << this->to_str(ind);
//    }

    string lterm::to_str() const{
        string str;
        auto c_new = _coef;
        auto p_new = _p;
        if (c_new->is_number()){
            string v = c_new->to_str();
            if (_sign) {
                if (v=="-1") {
                    str += " - ";
                }
                else if(v!="1") {
                    str += v;
                }
            }
            if(!_sign) {
                if (v == "-1") {
                    str += " + ";
                }
                else if (v.front()=='-'){
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
            if(_sign) {
                str += " + ";
            }
            str += "(";
            str += c_new->to_str();
            str += ")";
        }
        str += p_new->get_name(true,true);
        return str;
    }

    string lterm::to_str(size_t ind, int prec) const{
        string str;
        auto c_new = _coef;
        auto p_new = _p;
        unsigned dim = 1;
        if (c_new->_is_transposed) {
            dim = p_new->get_dim();
            for (auto inst = 0; inst <dim; inst++) {
                string v;
                if (c_new->is_number()){
                    v = c_new->to_str(prec);
                }
                else {
                    v += c_new->to_str(inst,prec);
                }
                if (_sign) {
                    if (v=="-1") {
                        str += " - ";
                    }
                    else {
                        str += " + ";
                        if(v!="1") {
                            str += v;
                        }
                    }
                }
                if(!_sign) {
                    if (v == "-1") {
                        str += " + ";
                    }
                    else if (v.front()=='-'){
                        str += " + ";
                        str += v.substr(1);
                    }
                    else if (v=="1"){
                        str += " - ";
                    }
                    else if(v!="-1"){
                        str += " - " + v;
                    }
                }
                str += p_new->get_name(inst);
            }
        }
        else{
            string v;
            if (c_new->is_number()){
                v = c_new->to_str(prec);
            }
            else {
                v += c_new->to_str(ind,prec);
            }
            if (_sign) {
                if (v=="-1") {
                    str += " - ";
                }
                else {
                    str += " + ";
                    if(v!="1") {
                        str += v;
                    }
                }
            }
            if(!_sign) {
                if (v == "-1") {
                    str += " + ";
                }
                else if (v.front()=='-'){
                    str += " + ";
                    str += v.substr(1);
                }
                else if (v=="1"){
                    str += " - ";
                }
                else if(v!="-1"){
                    str += " - " + v;
                }
            }
            str += p_new->get_name(ind);
        }
        return str;
    }


//    void lterm::print(size_t ind) const{
//        cout << this->to_str(ind);
//    }

    string lterm::to_str(size_t ind, size_t inst, int prec) const{
        string str;
        auto c_new = _coef;
        auto p_new = _p;
        unsigned dim = 1;
        if (c_new->_is_transposed) {
            dim = p_new->get_dim(inst);
        }
        for (unsigned idx = 0; idx <dim; idx++) {
            if (c_new->constant_::is_number()){
                string v;
                v = c_new->to_str(prec);
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
                if (!c_new->is_function() || !(c_new->is_unit())) {
                    str += "(";
                    if (c_new->is_matrix()) {
                        if (c_new->is_function() && !(c_new)->is_constant()){
                            str += c_new->to_str(inst);
                        }
                        else {
                            str += "[\n";
                            for (size_t i = 0; i<c_new->_dim[0]; i++) {
                                for (size_t j = 0; j<c_new->_dim[1]; j++) {
                                    str += c_new->to_str(i, j)+ " ";
                                }
                                str += "\n";
                            }
                            str += "\n] * ";
                        }
                    }
                    else if (c_new->_is_transposed) {
                        str += c_new->to_str(inst, idx);
                    }
                    else{
                        str += c_new->to_str(inst);
                    }
                    str += ")";
                }
            }
            if (p_new->is_matrix()) {
                str += "[\n";
                for (size_t i = 0; i<c_new->_dim[0]; i++) {
                    for (size_t j = 0; j<c_new->_dim[1]; j++) {
                        str += p_new->to_str(i, j)+ " ";
                    }
                    str += "\n";
                }
                str += "\n]";
            }
            else if (c_new->_is_transposed || (p_new->is_indexed() && p_new->_indices->_ids->size()>1)) {
                str += p_new->to_str(inst, idx);
            }
            else{
                str += p_new->get_name(ind);
            }
        }
        return str;
    }

}
