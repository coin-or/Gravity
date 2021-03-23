    //
    //  constraint.cpp
    //  Gravity
    //
    //
    //
#include <math.h>
#include <gravity/constraint.h>


using namespace gravity;
/** Constructor */
    //@{

/* Accessors */




ConstraintType Constraint_::get_ctype() const{
    return _ctype;
};


double l2norm(const vector<double>& x)
{
    double res=0;
    for(auto i=0;i<x.size();i++)
    {
        res+=x[i]*x[i];
    }
    res=std::sqrt(res);
    return(res);
}

/** Assuming an outer point and an inner point, this function uses a binary line search to find an active point for the current constraint and updates the value of the variables.
 @param[in] x_start: interior point
 @param[in] nb_inst: instance number
 @param[in] ctype: ineq type
 @return True if line search successfully solved
 The function assumes that the current value stored in vars is the outer point.
 Interior and outer point classification depends on constraint type (\geq 0 or \leq 0) as input by ctype
 **/
template<typename type>
bool Constraint<type>::binary_line_search(const vector<double>& x_start, size_t nb_inst)
{
    bool success=false;
    const double int_tol=1e-14, zero_tol=1e-10;
    const int max_iter=1000;
    vector<double> x_f, x_t, x_end, interval, mid;
    double  f_a,f_b,f_f, f_t, f_mid, interval_norm;
    int iter=0;
    x_end=this->get_x(nb_inst);/* outer point */
    this->uneval();
    f_b=this->eval(nb_inst);
    
    this->set_x(nb_inst, x_start);/* inner point */
    this->uneval();
    f_a=this->eval(nb_inst);
    
    if(this->_ctype==leq)
    {
        f_f=f_a;
        f_t=f_b;
        x_f=x_start;
        x_t=x_end;
    }
    else
    {
        f_f=f_b;
        f_t=f_a;
        x_f=x_end;
        x_t=x_start;
    }
    if(x_start.size()!=x_t.size() || x_start.size()!=x_f.size() || x_start.size()!=x_end.size()){
        throw invalid_argument("size mismatch");
    }
    for(auto i=0;i<x_start.size();i++)
    {
        interval.push_back(x_t[i]-x_f[i]);
        mid.push_back((x_end[i]+x_start[i])*0.5);
    }
    
    interval_norm=l2norm(interval);
    
    if(f_f<=(zero_tol*(-1)) && f_t>=zero_tol)
    {
        while(interval_norm>=int_tol && iter<=max_iter)
        {
            for(auto i=0;i<x_start.size();i++)
            {
                mid[i]=(x_f[i]+x_t[i])*0.5;
            }
            this->set_x(nb_inst, mid);
            this->uneval();
            f_mid=this->eval(nb_inst);
                // DebugOn("F_mid "<< f_mid<<endl);
                //DebugOn("xf\t xt\t xmid"<<endl);
                //                    for(auto i=0;i<x_start.size();i++)
                //                    {
                //                        //  DebugOn(x_f[i]<<"\t"<<x_t[i]<<"\t"<<mid[i]<<endl);
                //
                //                    }
                //
            if(f_mid>=zero_tol)
            {
                x_t=mid;
            }
            else if(f_mid<=zero_tol*(-1.))
            {
                x_f=mid;
            }
            else
            {
                if(true||(this->_ctype==leq) && (f_mid>=0)){
                    success=true;
                    
                    
                    DebugOff("F_mid "<<f_mid<<endl);
                    DebugOff("Interval_Norm "<<interval_norm<<endl);
                    DebugOff("Iter "<<iter<<endl);
                    break;
                }
                if(true||(this->_ctype==geq) && (f_mid<=0)){
                    success=true;
                    
                    
                    DebugOff("F_mid "<<f_mid<<endl);
                    DebugOff("Interval_Norm "<<interval_norm<<endl);
                    DebugOff("Iter "<<iter<<endl);
                    break;
                }
                
                else{
                    break;
                }
            }
            for(auto i=0;i<x_start.size();i++)
            {
                interval[i]=x_t[i]-x_f[i];
            }
            
            interval_norm=l2norm(interval);
            iter++;
        }
        
    }
    
    
    return success;
}

template bool Constraint<double>::binary_line_search(const vector<double>& x_start, size_t nb_inst);
/* Operators */









