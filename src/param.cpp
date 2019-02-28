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
    
    template<typename T>
    void param<T>::copy_vals(const shared_ptr<param_>& p){
        auto dim = p->get_dim();
        func<T> f;
        _val->resize(dim);
        for (size_t i = 0; i < dim; i++) {
            _val->at(i) = f.eval(p,i);
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
    
    
    //    template<>
    //    void param<complex<double>>::set_vals(const Eigen::SparseMatrix<complex<double>,Eigen::RowMajor>& SM){
    //        if (!is_complex()) {
    //            throw invalid_argument("Function void set_complex_vals(const Eigen::SparseMatrix<complex<double>,Eigen::RowMajor>& SM) is only implemented for complex<double> typed params/vars");
    //        }
    //        for (size_t k=0; k<SM.outerSize(); ++k) {
    //            for (Eigen::SparseMatrix<complex<double>,Eigen::RowMajor>::InnerIterator it(SM,k); it; ++it){
    //                set_val(it.row(), it.col(), it.value());
    //            }
    //        }
    //    }
    
    
    
    param<Cpx> conj(const param<Cpx>& p){
        param<Cpx> newp(p);
        if(newp._is_conjugate){
            newp._name = newp._name.substr(newp._name.find("("),newp._name.find(")"));
        }
        else {
            newp._name = "conj("+newp._name+")";
        }
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
