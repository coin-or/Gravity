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

    lterm::lterm(bool sign, shared_ptr<constant_> coef, shared_ptr<param_> p):lterm(){
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
        _coef_p1_tr = q._coef_p1_tr;
        return *this;
    }


    qterm& qterm::operator=(qterm&& q){
        _coef = move(q._coef);
        _p = move(q._p);
        _sign = q._sign;
        _coef_p1_tr = q._coef_p1_tr;
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
//        if (_coef_p1_tr) { // qterm = (coef*p1)^T*p2
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
//    string clean_print(bool pos, const string& v, bool brackets = false){
//        if(pos){
//            if (v=="-1" || v==" - 1" || v=="(-1,0)") {
//                return " - ";
//            }
//            else if (v.front()=='-'){
//                return " - " + v.substr(1);
//            }
//            else if(v=="1" || v==" + 1" || v=="(1,0)") {
//                return " + ";
//            }
//            else if(brackets){
//                return " + ("+v+")";
//            }
//            else{
//                return " + " + v;
//            }
//        }
//        else {
//            if (v == "-1" || v==" - 1" || v=="(-1,0)") {
//                return " + ";
//            }
//            else if (v.front()=='-'){
//                return " + " + v.substr(1);
//            }
//            else if (v=="1" || v==" + 1" || v=="(1,0)"){
//                return " - ";
//            }
//            else if(brackets){
//                return " - ("+v+")";
//            }
//            else{
//                return " - " + v;
//            }
//        }
//    }
    
    string print_expo(int exp){
        string str;
        if (exp != 1) {
            switch (exp) {
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
                    str += "^" + to_string(exp);
                    break;
            }
        }
        return str;
    }
    
    string pterm::print_poly_vars() const{
        string str;
        for (auto& p: *_l) {
            str += p.first->get_name(true,false);
            str += print_expo(p.second);
        }
        return str;
    }
    
    string pterm::print_poly_vars(size_t inst) const{
        string str;
        for (auto& p: *_l) {
            str += p.first->get_name(inst);
            str += print_expo(p.second);
        }
        return str;
    }
    
    string pterm::print_poly_vars(size_t inst1, size_t inst2) const{
        string str;
        for (auto& p: *_l) {
            str += p.first->get_name(inst1, inst2);
            str += print_expo(p.second);
        }
        return str;
    }
    
    string pterm::to_str() const{
        string str;
        if (_coef->is_number()){
            str += clean_print(_sign,_coef->to_str());
        }
        else{
            str += clean_print(_sign,_coef->to_str(),true);
        }
        str += print_poly_vars();
        return str;
    }
    
    string pterm::to_str(size_t ind, int prec) const{
        string str;
        if (_coef->is_number()){
            str += clean_print(_sign,_coef->to_str(prec));
        }
        else {
            str += clean_print(_sign,_coef->to_str(ind,prec));
        }
        str += print_poly_vars(ind);
        return str;
    }
    
    string pterm::to_str(size_t ind1, size_t ind2, int prec) const{
        string str;
        if (_coef->is_number()){
            str += clean_print(_sign,_coef->to_str(prec));
        }
        else {
            str += clean_print(_sign,_coef->to_str(ind1,ind2,prec));
        }
        str += print_poly_vars(ind1,ind2);
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
            str += clean_print(_sign,c_new->to_str());
        }
        else{
            str += clean_print(_sign,c_new->to_str(),true);
        }
        str += p_new1->get_name(true,false);
        if(_coef_p1_tr){
            str += ")\u1D40";
        }
        else if (p_new1==p_new2 || *p_new1==*p_new2) {
            str += "²";
            return str;
        }
        str += p_new2->get_name(true,false);
        return str;
    }
    
    string qterm::to_str(size_t ind, int prec) const {
        string str;
        auto c_new = _coef;
        auto p_new1 = _p->first;
        auto p_new2 = _p->second;
        string coef;
        if (p_new1->is_matrix_indexed() || c_new->_is_transposed) {
            if (c_new->_is_transposed) {
                auto dim = c_new->get_dim(0);
                string term="";
                for(auto i=0; i<dim;i++){
                    term = print_transposed(i,prec);
                    if(str.size()>0 && term.size()>0 && i!=0){
                        str += " + " + term;
                    }
                    else {
                        str += term;
                    }
                }
            }
            else {
                str += print_row(ind,prec);
            }
        }
        else{
            if (c_new->is_number()){
                coef = c_new->to_str(prec);
            }
            else {
                coef = c_new->to_str(ind,prec);
            }
            if(_coef_p1_tr){
                str += "(";
            }
            str += clean_print(_sign,coef);
            str += p_new1->get_name(ind);
            if(_coef_p1_tr){
                str += ")\u1D40";
            }
            else if (p_new1==p_new2 || *p_new1==*p_new2) {
                str += "²";
                return str;
            }
            str += p_new2->get_name(ind);
        }
        return str;
    }

    /* Return a new polynomial term excluding p */
    pterm pterm::exclude(const shared_ptr<param_>& p) const{
        pterm new_pt(*this);
        auto new_l = make_shared<list<pair<shared_ptr<param_>, int>>>();
        for (auto &pair : *_l) {
            if(pair.first->get_name(true, false)!=p->get_name(true, false))
                new_l->push_back(make_pair<>(pair.first, pair.second));
        }
        new_pt._l = new_l;
        return new_pt;
    }
    
    string qterm::to_str(size_t ind1, size_t ind2, int prec) const {
        string str, coef;
        auto p_new2 = _p->second;
        if (_coef->is_number()){
            coef = _coef->to_str(prec);
        }
        else {
            coef = _coef->to_str(ind1,ind2,prec);
        }
        if(_coef_p1_tr){
            str += "(";
        }
        str += clean_print(_sign,coef);
        str += _p->first->get_name(ind1,ind2);
        if(_coef_p1_tr){
            str += ")\u1D40";
        }
        else if (_p->first==_p->second) {
            str += "²";
            return str;
        }
        str += _p->second->get_name(ind1,ind2);
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
            str += clean_print(_sign,c_new->to_str());
        }
        else{
            str += clean_print(_sign,c_new->to_str(),true);
        }
        str += p_new->get_name(true,false);
        return str;
    }
    
    string lterm::print_transposed(int prec) const{
        auto dim = _p->get_dim();
        string str;
        for (auto idx = 0; idx <dim; idx++) {
            string coef;
            if (_coef->is_number()){
                coef = _coef->to_str(prec);
            }
            else {
                coef = _coef->to_str(idx,prec);
            }
            str += clean_print(_sign,coef);
            str += _p->get_name(idx);
        }
        return str;
    }
    
    string qterm::print_transposed(int prec) const{
        auto dim = _p->first->get_dim(0);
        string str;
        for (auto idx = 0; idx <dim; idx++) {
            string coef;
            if (_coef->is_number()){
                coef = _coef->to_str(prec);
            }
            else {
                coef = _coef->to_str(idx,prec);
            }
            str += clean_print(_sign,coef);
            if (_p->first==_p->second) {
                str += _p->first->get_name(idx)+ "²";
            }
            else {
                str += _p->first->get_name(idx)+ _p->second->get_name(idx);
            }
        }
        return str;
    }

    string qterm::print_row(size_t inst, int prec) const{
        if(!_coef_p1_tr && !_p->first->is_matrix_indexed()){
            return print_transposed(prec);
        }
        string str;
        size_t dim = 0;
        if (_coef_p1_tr) { // qterm = (coef*p1)^T*p2
            for (auto i = 0; i<_p->first->get_dim(); i++) {
                for (auto j = 0; j<_coef->get_dim(i); j++) {
                    string coef;
                    if (_coef->is_number()){
                        coef = _coef->to_str(prec);
                    }
                    else {
                        coef = _coef->to_str(i,j,prec);
                    }
                    str += clean_print(_sign,coef);
                    auto name1 = _p->first->get_name(j);
                    auto name2 = _p->second->get_name(i);
                    if (name1==name2) {
                        str += name1 + "²";
                    }
                    else {
                        str += name1 + name2;
                    }
                }
            }
            return str;
        }

        dim = std::min(_p->first->get_dim(inst),_p->second->get_dim(inst));
        if(dim==0){
            return str;
        }
        
        for (auto idx = 0; idx <dim; idx++) {
            string coef;
            if (_coef->is_number()){
                coef = _coef->to_str(prec);
            }
            else {
                coef = _coef->to_str(inst, idx,prec);
            }
            str += clean_print(_sign,coef);
            string name1, name2;
            if(_p->first->is_matrix_indexed())
                name1 = _p->first->get_name(inst,idx);
            else
                name1 = _p->first->get_name(inst);
            if(_p->second->is_matrix_indexed())
                name2 = _p->second->get_name(inst,idx);
            else
                name2 = _p->second->get_name(inst);
            if (name1==name2) {
                str += name1 + "²";
            }
            else {
                str += name1 + name2;
            }
        }
        return str;
    }
    
    string qterm::print_transposed(size_t inst, int prec) const{
        if(!_coef_p1_tr && !_p->first->is_matrix_indexed()){
            return print_transposed(prec);
        }
        string str;
        size_t dim = 0;
        if (_coef_p1_tr) { // qterm = (coef*p1)^T*p2
            for (auto i = 0; i<_p->first->get_dim(); i++) {
                for (auto j = 0; j<_coef->get_dim(i); j++) {
                    string coef;
                    if (_coef->is_number()){
                        coef = _coef->to_str(prec);
                    }
                    else {
                        coef = _coef->to_str(i,j,prec);
                    }
                    str += clean_print(_sign,coef);
                    auto name1 = _p->first->get_name(j);
                    auto name2 = _p->second->get_name(i);
                    if (name1==name2) {
                        str += name1 + "²";
                    }
                    else {
                        str += name1 + name2;
                    }
                }
            }
            return str;
        }

        dim = std::min(_p->first->get_dim(inst),_p->second->get_dim(inst));
        if(dim==0){
            return str;
        }
        
        for (auto idx = 0; idx <dim; idx++) {
            string coef;
            if (_coef->is_number()){
                coef = _coef->to_str(prec);
            }
            else {
                coef = _coef->to_str(inst, idx,prec);
            }
            str += clean_print(_sign,coef);
            auto name1 = _p->first->get_name(inst,idx);
            auto name2 = _p->second->get_name(inst,idx);
            if (name1==name2) {
                str += name1 + "²";
            }
            else {
                str += name1 + name2;
            }
        }
        return str;
    }
    
    string lterm::print_transposed(size_t inst, int prec) const{
        if(!_p->is_matrix_indexed() && !_coef->is_matrix()){
            return print_transposed(prec);
        }
        string str;
        auto dim = _p->get_dim(inst);
        if(dim==0){
            return str;
        }
        
        for (auto idx = 0; idx <dim; idx++) {
            string coef;
            if (_coef->is_number()){
                coef = _coef->to_str(prec);
            }
            else {
                coef = _coef->to_str(inst, idx,prec);
            }
            str += clean_print(_sign,coef);
            if(_p->is_matrix_indexed()){
                str += _p->get_name(inst,idx);
            }
           else {
               str += _p->get_name(idx);
           }
        }
        return str;
    }
    
    
    string lterm::to_str(size_t ind, int prec) const{
        string str;
        auto c_new = _coef;
        auto p_new = _p;
        if (p_new->_is_vector || p_new->is_matrix_indexed()) {
            str += print_transposed(ind,prec);
        }
        else{
            string coef;
            if (c_new->is_number()){
                coef = c_new->to_str(prec);
            }
            else {
                coef = c_new->to_str(ind,prec);
            }
            str += clean_print(_sign,coef);
            str += p_new->get_name(ind);
        }
        return str;
    }


//    void lterm::print(size_t ind) const{
//        cout << this->to_str(ind);
//    }

    string lterm::to_str(size_t ind1, size_t ind2, int prec) const{
        string str;
        auto c_new = _coef;
        auto p_new = _p;
        if (c_new->_is_transposed) {
            str += print_transposed(ind1, prec);
        }
        else {
            str += clean_print(_sign,c_new->to_str(ind1,ind2,prec));
        }
        str += p_new->get_name(ind1,ind2);
        return str;
    }

}
