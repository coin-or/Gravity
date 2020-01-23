//
//  param.cpp
//
//
//  Created by Hassan on 20/05/2016.
//
//
#include <gravity/param.h>
#include <gravity/func.h>

namespace gravity {
    

/* Deletes the ith column for the matrix indeed param/var */
indices param_::delete_column(size_t i){
    indices deleted_col = indices(_indices->get_name()+"\\col_"+to_string(i));
   if(_indices->_type!=matrix_){
        throw invalid_argument("delete_column(size_t i) can only be called on a matrix indexed param/var");
    }
    for (size_t j = 0; j<_indices->_ids->size(); j++) {
        if(_indices->_ids->at(j).size()<=i){
            throw invalid_argument("delete_column(size_t i) out of range: " + to_string(i) + " at row " + to_string(j));
        }
        deleted_col.add(_indices->_keys->at(_indices->_ids->at(j).at(i)));
        vector<size_t> new_row;
        for (size_t k = 0; k<_indices->_ids->at(j).size(); k++) {
            if(k!=i)
                new_row.push_back(_indices->_ids->at(j).at(k));
        }
        _indices->_ids->at(j) = new_row;
    }
    _name += +"\\col_"+to_string(i);
    return deleted_col;
}

/* Returns a pair of boolean vectors <z_v, nnz_v> indicating which row p has a zero/non-zero coefficient in */
pair<vector<bool>,vector<bool>> param_::get_nnz_rows() const{
    if(_indices->_type!=matrix_){
        throw invalid_argument("get_nnz_rows() can only be called on a matrix indexed param/var");
    }
    pair<vector<bool>,vector<bool>> nnz_rows;
    nnz_rows.first.resize(_indices->get_nb_rows(),false);
    nnz_rows.second.resize(_indices->get_nb_rows(),false);
    for (size_t i = 0; i<_indices->_ids->size(); i++) {
        if(_indices->_ids->at(i).size()>0){
            nnz_rows.second[i] = true;
        }
        else{
            nnz_rows.first[i] = true;
        }
    }
    return nnz_rows;
}

    template<typename T>
    void param<T>::copy_vals(const shared_ptr<param_>& p){
        switch (p->get_intype()) {
            case binary_:{
                auto pp =  static_pointer_cast<param<bool>>(p);
                copy_vals_(*pp);
            }
                break;
            case short_:{
                auto pp =  static_pointer_cast<param<short>>(p);
                copy_vals_(*pp);
            }
                break;
            case integer_:{
                auto pp =  static_pointer_cast<param<int>>(p);
                copy_vals_(*pp);
            }
                break;
            case float_:{
                auto pp =  static_pointer_cast<param<float>>(p);
                copy_vals_(*pp);
            }
                break;
            case double_:{
                auto pp =  (param<double>*)(p.get());
                copy_vals_(*pp);
            }
                break;
            case long_:{
                auto pp =  static_pointer_cast<param<long double>>(p);
                copy_vals_(*pp);
            }
                break;
            case complex_:{
                auto pp =  static_pointer_cast<param<Cpx>>(p);
                copy_vals_(*pp);
            }
                break;
            default:
                break;
        }
    }
    
    template<>
    void param<complex<double>>::update_range(const complex<double>& val) {
        if(_range->first.real()>val.real()){
            _range->first.real(val.real());
        }
        if(_range->second.real()<val.real()){
            _range->second.real(val.real());
        }
        if(_range->first.imag()>val.imag()){
            _range->first.imag(val.imag());
        }
        if(_range->second.imag()<val.imag()){
            _range->second.imag(val.imag());
        }
    }
    
        
    
    param<Cpx> conj(const param<Cpx>& p){
        param<Cpx> newp(p);
        if(newp._is_conjugate){
            newp._name = newp._name.substr(newp._name.find("("),newp._name.find(")"));
        }
        else {
            newp._name = "conj("+newp._name+")";
        }
        newp._is_conjugate = !newp._is_conjugate;
        return newp;
    }
    
    param<Cpx> ang(const param<Cpx>& p){
        param<Cpx> newp(p);
        newp._is_angle = true;
        newp._name = "ang("+newp._name+")";
        return newp;
    }
    
    param<Cpx> sqrmag(const param<Cpx>& p){
        param<Cpx> newp(p);
        newp._is_sqrmag = true;
        newp._name = "|"+newp._name+"|²";
        return newp;
    }
    
    param<Cpx> real(const param<Cpx>& p){
        param<Cpx> newp(p);
        newp._name = "real("+newp._name+")";
        newp._is_real = true;
        return newp;
    }
    
    param<Cpx> imag(const param<Cpx>& p){
        param<Cpx> newp(p);
        newp._name = "imag("+newp._name+")";
        newp._is_imag = true;
        return newp;
    }
    
    template class param<bool>;
    template class param<short>;
    template class param<int>;
    template class param<float>;
    template class param<double>;
    template class param<long double>;
    template class param<Cpx>;
}
