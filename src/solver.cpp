//
//  solver.cpp
//  Gravity++
//
//  Created by Hassan on 30/01/2015.

//

#include <gravity/solver.h>

#ifdef USE_BONMIN
#include <coin/BonBonminSetup.hpp>
#include <coin/BonCbc.hpp>
#endif
using namespace gravity;
using namespace std;

void gurobiNotAvailable()
{
    cerr << "Can't use Gurobi as a solver: this version of Gravity "
    "was compiled without Gurobi support." << endl;
    exit(1);
}
void cplexNotAvailable()
{
    cerr << "Can't use Cplex as a solver: this version of Gravity "
    "was compiled without Cplex support." << endl;
    exit(1);
}

void bonminNotAvailable()
{
    cerr << "Can't use Bonmin as a solver: this version of Gravity "
    "was compiled without Bonmin support." << endl;
    exit(1);
}
void ipoptNotAvailable()
{
    cerr << "Can't use Ipopt as a solver: this version of Gravity "
    "was compiled without Ipopt support." << endl;
    exit(1);
}


void mosekNotAvailable()
{
    cerr << "Can't use Mosek as a solver: this version of Gravity "
    "was compiled without Mosek support." << endl;
    exit(1);
}
void ClpNotAvailable()
{
    cerr << "Can't use Clp as a solver: this version of Gravity "
    "was compiled without Clp support." << endl;
    exit(1);
}

namespace gravity {
    /** Outer approximation of model. Throws exception if model has nonlinear equality constraints
     @param[in] nb_discr: number of OA cuts per nonlinear constraint
     @return Model with OA cuts. OA cuts are added to the model (for all func instances) in an uniform grid (nb_discr)
     **/
    template<typename type>
    template<typename T>
    shared_ptr<Model<type>> Model<type>::buildOA()
    {
        
        auto cpy = this->copy();
        //        cpy->initialize_all(1);
        //    this->print_solution();
        vector<double> xsolution(_nb_vars);
        vector<double> xinterior(_nb_vars);
        vector<double> xcurrent;
        get_solution(xsolution);
        
        
        auto OA=make_shared<Model<>>(_name+"-OA Model");
        for (auto &it: cpy->_vars)
        {
            auto v = it.second;
            if(!OA->has_var(*v)){
                OA->add_var(v);
            }
        }
        auto obj=*cpy->_obj;
        if(_objt==minimize){
            OA->min(obj);
        }
        else {
            OA->max(obj);
        }
        string cname;
        for (auto &con: cpy->_cons_vec)
        {
            if(!con->is_linear()) {
                if(con->_ctype==eq)
                {
                    throw invalid_argument("Equality constraint is not currently supported");
                }
            }
            else
            {
                OA->add(*con);
            }
        }
        set_solution(xsolution);
        return OA;
    }
    //    template<typename type>
    //    template<typename T>
    //    shared_ptr<Model<type>> Model<type>::buildOA(int nb_discr, int nb_perturb)
    //    {
    //
    //
    //        //    this->print_solution();
    //        vector<double> xsolution(_nb_vars);
    //        vector<double> xinterior(_nb_vars);
    //        vector<double> xcurrent;
    //        get_solution(xsolution);
    //
    //
    //        auto OA=make_shared<Model<>>(_name+"-OA Model");
    //        for (auto &it: _vars)
    //        {
    //            auto v = it.second;
    //            if(!OA->has_var(*v)){
    //                OA->add_var(v);
    //            }
    //        }
    //        auto obj=*_obj;
    //        if(_objt==minimize){
    //            OA->min(obj);
    //        }
    //        else {
    //            OA->max(obj);
    //        }
    //        string cname;
    //        for (auto &con: _cons_vec)
    //        {
    //            if(!con->is_linear()) {
    //                if(con->_ctype==eq)
    //                {
    //                    DebugOn("Exception: Equality constraint is not currently supported"<<endl);
    //                    DebugOn("Throw exception" <<endl);
    //
    //                }
    //                else
    //                {
    //
    //                    OA->add_outer_app_uniform(nb_discr, *con);
    //                }
    //            }
    //            else
    //            {
    //                OA->add(*con);
    //            }
    //        }
    //        set_solution(xsolution);
    //        OA->add_outer_app_active(*this, nb_perturb);
    //
    //        set_solution(xsolution);
    //        return OA;
    //    }
    
    
    /** Returns a model that can compute an interior point to the current model
     **/
    template<typename type>
    template<typename T>
    Model<type> Model<type>::build_model_interior() const
    {
        Model<type> Interior(_name+"Interior");
        
        /* Variables */
        for (auto &it: _vars)
        {
            auto v = it.second;
            Interior.add_var(v);
        }
        indices ind_eta("ind_eta");
        vector<indices> ind_eta_vec;
        for (auto &con: _cons_vec)
        {
            if(!con->is_linear()) {
                if(!con->is_convex() || con->is_rotated_soc() || con->check_soc()){
                    indices ind_eta_c("ind_eta_c");
                    for(auto i=0;i<con->get_nb_instances();i++){
                        ind_eta.add(con->_name+to_string(i));
                        ind_eta_c.add(con->_name+to_string(i));
                    }
                    ind_eta_vec.push_back(ind_eta_c);
                }
            }
        }
        var<> eta_int("eta_interior", -1, 0);
        Interior.add(eta_int.in(ind_eta));
        
        /* Objective */
        Interior.min(sum(eta_int));
        
        /* Constraints */
        int count=0;
        for (auto &con: _cons_vec)
        {
            if(!con->is_linear()) {
                /* We are only interested in an iterior point for constraints defining a convex region but having a non-convex description, e.g., SDP-determinant cuts and SOC constraints.*/
                if(!con->is_convex() || con->is_rotated_soc() || con->check_soc()){
                    auto ind=ind_eta_vec[count++];
                    
                    Constraint<> Inter_con(*con);
                    
                    if(con->_ctype==leq)
                    {
                        Interior.add(Inter_con<=eta_int.in(ind));
                    }
                    else  if(con->_ctype==geq)
                    {
                        Interior.add(Inter_con>=-1*eta_int.in(ind));
                    }
                }
            }
        }
        return *Interior.copy();
    }
    
    /* Runds models stored in the vector in parallel, using solver of stype and tolerance tol */
    int run_parallel_new(const std::vector<std::string> objective_models, std::vector<double>& sol_obj, std::vector<int>& sol_status, const std::vector<shared_ptr<gravity::Model<double>>>& models, gravity::SolverType stype, double tol, unsigned nr_threads, const string& lin_solver, int max_iter){
        std::vector<thread> threads;
        std::vector<bool> feasible;
        std::string mname, msname,vname, key, dir;
        var<> var;
        if(models.size()==0){
            DebugOff("in run_parallel(models...), models is empty, returning");
            return -1;
        }
        for (auto s=0;s<objective_models.size();s++){
            msname=objective_models[s];
            mname=msname;
            std::size_t pos = msname.find("|");
            vname.assign(msname, 0, pos);
            msname=msname.substr(pos+1);
            pos=msname.find("|");
            key.assign(msname, 0, pos);
            dir=msname.substr(pos+1);
            var=models[s]->get_var<double>(vname);
            if(dir=="LB")
            {
                models[s]->min(var(key));
            }
            else
            {
                models[s]->max(var(key));
            
            }
            models[s]->reindex();
        }
        /* Split models into nr_threads parts */
        auto nr_threads_ = std::min((size_t)nr_threads,models.size());
        std::vector<size_t> limits = bounds(nr_threads_, models.size());
        DebugOff("Running on " << nr_threads_ << " threads." << endl);
        DebugOff("limits size = " << limits.size() << endl);
        for (size_t i = 0; i < limits.size(); ++i) {
            DebugOff("limits[" << i << "] = " << limits[i] << endl);
        }
        /* Launch all threads in parallel */
        auto vec = vector<shared_ptr<gravity::Model<double>>>(models);
        for (size_t i = 0; i < nr_threads_; ++i) {
            threads.push_back(thread(run_models<double>, ref(vec), limits[i], limits[i+1], stype, tol, lin_solver, max_iter));
        }
        /* Join the threads with the main thread */
        for(auto &t : threads){
            t.join();
        }
        for(auto &m:models){
            sol_status.push_back(m->_status);
            sol_obj.push_back(m->get_obj_val());
        }
        return 0;
    }
    
    
    
    /* Runds models stored in the vector in parallel, using solver of stype and tolerance tol */
    int run_parallel(const vector<shared_ptr<gravity::Model<double>>>& models, gravity::SolverType stype, double tol, unsigned nr_threads, const string& lin_solver, int max_iter){
        std::vector<thread> threads;
        std::vector<bool> feasible;
        if(models.size()==0){
            DebugOff("in run_parallel(models...), models is empty, returning");
            return -1;
        }
        /* Split models into nr_threads parts */
        auto nr_threads_ = std::min((size_t)nr_threads,models.size());
        std::vector<size_t> limits = bounds(nr_threads_, models.size());
        DebugOff("Running on " << nr_threads_ << " threads." << endl);
        DebugOff("limits size = " << limits.size() << endl);
        for (size_t i = 0; i < limits.size(); ++i) {
            DebugOff("limits[" << i << "] = " << limits[i] << endl);
        }
        /* Launch all threads in parallel */
        auto vec = vector<shared_ptr<gravity::Model<double>>>(models);
        for (size_t i = 0; i < nr_threads_; ++i) {
            threads.push_back(thread(run_models<double>, ref(vec), limits[i], limits[i+1], stype, tol, lin_solver, max_iter));
        }
        /* Join the threads with the main thread */
        for(auto &t : threads){
            t.join();
        }
        return 0;
    }
    
    int run_parallel(const vector<shared_ptr<gravity::Model<double>>>& models, gravity::SolverType stype, double tol, unsigned nr_threads, int max_iter){
        return run_parallel(models,stype,tol,nr_threads,"",max_iter);
    };
    
    /** Discretizes Constraint con and adds OA cuts to the model that calls it. Discretization of squared constraint only currently implemented
     @param[in] nb_discr:
     @param[in] con:
     @return void. OA cuts are added to the model that calls the function (for all func instances)
     **/
    //     template<>
    //    void Model<>::add_outer_app_uniform(int nb_discr, Constraint<> con)
    //    {
    //
    //        func<> res;
    //        double lb,ub;
    //        size_t posv;
    //        bool sign_x, sign_y;
    //        double coef_x, coef_y, xval;
    //        if(con.is_quadratic() && con._lterms->size()==1 && con._qterms->size()==1 && con._qterms->begin()->second._p->first==con._qterms->begin()->second._p->second) //This if is specific to constraints of the form ay- bx^2 \geq 0 or bx^2-ay \leq 0
    //        {
    //
    //            auto x=static_pointer_cast<var<double>>(con._qterms->begin()->second._p->first);
    //            auto y=static_pointer_cast<var<double>>(con._lterms->begin()->second._p);
    //            indices D("D");
    //            for(auto d=0;d<nb_discr;d++){
    //                D.add("D"+to_string(d));
    //            }
    //            indices I("I");
    //            for(auto i=0;i<con.get_nb_inst();i++){
    //                I.add("I"+to_string(i));
    //            }
    //
    //            indices UniDI(D, I);
    //            param<double> cx("Paramcx"+con._name);
    //            cx.in(UniDI);
    //            param<double> cy("Paramcy"+con._name);
    //            cy.in(UniDI);
    //            param<double> c0("Paramc0"+con._name);
    //            c0.in(UniDI);
    //            sign_x=con._qterms->begin()->second._sign;
    //            sign_y=con._lterms->begin()->second._sign;
    //            for(auto i=0;i<con.get_nb_inst();i++)
    //            {
    //                posv=x->get_id_inst(i);
    //                lb=x->get_double_lb(posv);
    //                ub=x->get_double_ub(posv);
    //                coef_x=con.eval_coef(con._qterms->begin()->second._coef, i);
    //                coef_y=con.eval_coef(con._lterms->begin()->second._coef, i);
    //                for(auto d=0;d<nb_discr;d++)
    //                {
    //                    xval= lb+d*(ub-lb)/nb_discr;
    //                    if(sign_y){
    //                        cy.set_val("D"+to_string(d)+",I"+to_string(i), coef_y);
    //                    }
    //                    else{
    //                        cy.set_val("D"+to_string(d)+",I"+to_string(i), coef_y*(-1));
    //                    }
    //                    if(sign_x){
    //                        cx.set_val("D"+to_string(d)+",I"+to_string(i), 2*coef_x*xval);
    //                        c0.set_val("D"+to_string(d)+",I"+to_string(i), coef_x*(-1)*xval*xval);
    //                    }
    //                    else{
    //                        cx.set_val("D"+to_string(d)+",I"+to_string(i), (-2)*coef_x*xval);
    //                        c0.set_val("D"+to_string(d)+",I"+to_string(i), coef_x*xval*xval);
    //                    }
    //
    //                }
    //            }
    //
    //
    //
    //            auto x_ids = indices(D,*x->_indices);
    //            auto y_ids = indices(D,*y->_indices);
    //            Constraint<> OA_uniform("OA_cuts_uniform "+con._name);
    //            OA_uniform=cy*(y->from_ith(1,y_ids))+cx*(x->from_ith(1,x_ids))+c0;
    //            if(con._ctype==leq){
    //                add(OA_uniform.in(UniDI)<=0);
    //            }
    //            else{
    //                add(OA_uniform.in(UniDI)>=0);
    //            }
    //             // OA_uniform.print();
    //        }
    //    }
    template<>
    void Model<>::add_outer_app_uniform(int nb_discr, Constraint<> con)
    {
        
        func<> res;
        double lb,ub;
        size_t posv;
        bool sign_x, sign_y;
        double coef_x, coef_y, xval;
        if(con.is_quadratic() && con._lterms->size()==1 && con._qterms->size()==1 && con._qterms->begin()->second._p->first==con._qterms->begin()->second._p->second) //This if is specific to constraints of the form ay- bx^2 \geq 0 or bx^2-ay \leq 0
        {
            
            auto x=static_pointer_cast<var<double>>(con._qterms->begin()->second._p->first);
            auto y=static_pointer_cast<var<double>>(con._lterms->begin()->second._p);
            indices D("D");
            for(auto d=0;d<nb_discr;d++){
                D.add("D"+to_string(d));
            }
            indices I("I");
            for(auto i=0;i<con.get_nb_inst();i++){
                I.add("I"+to_string(i));
            }
            
            indices UniDI(D, I);
            param<double> dc("d_"+con._name);
            dc.in(UniDI);
            
            double epsc=1e-6;
            
            
            for(auto i=0;i<con.get_nb_inst();i++)
            {
                for(auto d=0;d<nb_discr;d++)
                {
                    //  dc.set_val("D"+to_string(d)+",I"+to_string(i),d);
                    dc.set_val("D"+to_string(d)+",I"+to_string(i),d*1./nb_discr);
                    
                }
            }
            
            
            
            
            
            auto x_ids = indices(D,*x->_indices);
            auto y_ids = indices(D,*y->_indices);
            
            Constraint<> OA_uniform("OA_uniform "+con._name);
            func<> xval =x->get_lb().from_ith(1,x_ids)+dc.from_ith(0,UniDI)*(x->get_ub().from_ith(1,x_ids)-x->get_lb().from_ith(1,x_ids));
            
            OA_uniform=(2*x->from_ith(1,x_ids)*xval -xval*xval-y->from_ith(1,y_ids));
            add(OA_uniform.in(UniDI)<=0);
            //            OA_uniform=(2*(x->from_ith(1,x_ids))*(x->get_lb().from_ith(1,x_ids)+dc.from_ith(0,UniDI)*(x->get_ub().from_ith(1,x_ids)-x->get_lb().from_ith(1,x_ids))/nb_discr) -pow((x->get_lb().from_ith(1,x_ids)+dc.from_ith(0,UniDI)*(x->get_ub().from_ith(1,x_ids)-x->get_lb().from_ith(1,x_ids))/nb_discr),2)-y->from_ith(1,y_ids));
            //            add(OA_uniform.in(UniDI)<=0);
            
            
            //          OA_uniform.print();
            //            OA_uniform.print_symbolic();
        }
        /*TODO Else (discretization for general constraint)*/
        //        else
        //        {
        //            int N_coordinates=con._nb_vars;
        //            int N_dim=N_coordinates-1;
        //
        //            vector<double> xn(N_dim);
        //
        //            int c;
        //
        //            for(auto a=0;a<std::pow(nb_discr,N_dim);a++)
        //            {
        //                c=a;
        //                std::fill (xn.begin(),xn.end(),0);
        //                int pos=N_dim-1;
        //                while(c>0)//
        //                {
        //                    xn[pos--]=c%nb_discr;
        //                    c/=nb_discr;
        //                }
        //                for (auto &i:xn)
        //                    DebugOn(i<<"\t");
        //                DebugOn(endl);
        //            }
        //
        //
        //        }
    }
    
    template<>
    Model<> Model<>::add_outer_app_solution(Model<>& nonlin)
    {
        const double active_tol=1e-6,active_tol_sol=1e-12, perturb_dist=1e-1, zero_tol=1e-6;
        vector<double> xactive, xcurrent, xinterior, xres, xtest;
        bool interior=false;
        double fk, fkk;
        double rhs_tol = 0;
        bool outer;
        int count=0;
        size_t posv;
        vector<double> c_val ;
        double c0_val;
        vector<indices> vec_Pert;
        vector<param<double>> oa_vec_c;
        param<double> oa_c0;
        
        int output=0;
        double tol=1e-8;
        bool convex_region=true;
        double scale=1.0;
        auto Ointerior = nonlin.build_model_interior();
        solver<> modelI(Ointerior, ipopt);
        //Ointerior.print();
        modelI.run(output=0, tol);
        
        vector<double> xsolution(_nb_vars);
        nonlin.get_solution(xsolution);
        
        //    Ointerior.print();
        
        if((Ointerior._status==0||Ointerior._status==1) && Ointerior.get_obj_val() <0)
        {
            interior=true;
            //        Ointerior.print_solution();
        }
        for (auto &con: nonlin._cons_vec)
        {
            if(!con->is_linear()) {
                if(!con->is_convex() || con->is_rotated_soc() || con->check_soc())
                {
                    Constraint<> OA_sol("OA_cuts_"+con->_name);
                    indices activeset("active_"+con->_name);
                    auto keys=con->_indices->_keys;
                    for(auto i=0;i<con->get_nb_inst();i++){
                        //                    con->uneval();
                        /** Generate an OA cut if constraint is active or if it has a convex representation */
                        //if(con->is_active(i,active_tol_sol) || (con->is_convex() && !con->is_rotated_soc() && !con->check_soc())){
                        convex_region=true;
                        if(!con->is_convex()) //For the SDP determinant constraint, check if the point is feasible with repsecto to the SOC constraints
                        {
                            
                            xres=con->get_x(i);
                            auto soc1=std::pow(xres[0],2)+std::pow(xres[3],2)-xres[6]*xres[7];
                            auto soc2=std::pow(xres[1],2)+std::pow(xres[4],2)-xres[7]*xres[8];
                            auto soc3=std::pow(xres[2],2)+std::pow(xres[5],2)-xres[6]*xres[8];
                            if(soc1<=0 && soc2<=0 && soc3<=0){
                                convex_region=true;
                            }
                            else{
                                convex_region=false;
                            }
                        }
                        
                        if(convex_region){
                            con->uneval();
                            con->eval_all();
                            if(con->is_active(i,active_tol_sol)){
                                
                                activeset.add((*keys)[i]);
                            }
                        }
                    }
                    OA_sol=con->get_outer_app(activeset, scale);
                    if(con->_ctype==leq) {
                        add(OA_sol.in(activeset)<=0);
                    }
                    else {
                        add(OA_sol.in(activeset)>=0);
                    }
                }
                else if(con->is_quadratic() && con->_lterms->size()==1 && con->_qterms->size()==1 && con->_qterms->begin()->second._p->first==con->_qterms->begin()->second._p->second) //This if is specific to constraints of the form ay- bx^2 \geq 0 or bx^2-ay \leq 0
                {
                    double xval;
                    auto x=static_pointer_cast<var<double>>(con->_qterms->begin()->second._p->first);
                    Constraint<> OA_sol("OA_cuts_"+con->_name);
                    indices allset("active_"+con->_name);
                    auto keys=con->_indices->_keys;
                     for(auto i=0;i<con->get_nb_inst();i++){
                         posv=x->get_id_inst(i);
                         x->get_double_val(posv, xval);
                         if(std::abs(xval)>=zero_tol){
                             string keyi=(*keys)[i];
                             allset.add(keyi);
                         }
                     }
                    OA_sol=con->get_outer_app(allset, scale);
                    if(con->_ctype==leq) {
                        add(OA_sol.in(allset)<=0);
                    }
                    else {
                        add(OA_sol.in(allset)>=0);
                    }
                }
                    else if(con->is_convex() && !con->is_rotated_soc() && !con->check_soc())
                    {
                        //add_outer_app_uniform(5, *con);
                        //Need to reset after this
                        Constraint<> OA_sol("OA_cuts_"+con->_name);
                        indices allset("active_"+con->_name);
                        auto keys=con->_indices->_keys;
                        for(auto i=0;i<con->get_nb_inst();i++){
                            string keyi=(*keys)[i];
                            allset.add(keyi);
                            
                        }
                        OA_sol=con->get_outer_app(allset, scale);
                        if(con->_ctype==leq) {
                            add(OA_sol.in(allset)<=0);
                        }
                        else {
                            add(OA_sol.in(allset)>=0);
                        }
                    }
            }
        }
        nonlin.set_solution(xsolution);
        reindex();
        DebugOn("Number of constraints in linear model "<<_nb_cons<<endl);
        bool add_new, oa_cut;
        int nb_added_cuts=0;
        int nb_perturb=1;
        int count_var=0;
        if(interior){
            for (auto &con: nonlin._cons_vec)
            {
                if(!con->is_linear() && (!con->is_convex() || con->is_rotated_soc() || con->check_soc())) {
                    auto con_lin_name="OA_cuts_"+con->_name;
                    if(_cons_name.find(con_lin_name)!=_cons_name.end()){
                        add_new=false;
                    }
                    else{
                        add_new=true;
                    }
                    auto cnb_inst=con->get_nb_inst();
                    for(auto i=0;i<cnb_inst;i++){
                        auto cname=con->_name;
                        xcurrent=con->get_x(i);
                        count_var=0;
                        for(auto &it: *(con->_vars))
                        {
                            auto v = it.second.first;
                            auto vname=v->_name;
                            if(v->_is_vector)
                            {
                                throw invalid_argument("vector variables are not supported");
                                break;
                            }
                            else
                            {
                                posv=v->get_id_inst(i);
                                auto ub_v=v->get_double_ub(posv);
                                auto lb_v=v->get_double_lb(posv);
                                for(auto j=1;j<=nb_perturb;j++)
                                {
                                    outer=false;
                                    oa_cut=false;
                                    c0_val=0;
                                    c_val.resize(con->_nb_vars,0);
                                    v->set_double_val(posv, std::min(std::max(xcurrent[count_var]*(1 - j*perturb_dist), lb_v),ub_v)); /** Perturbed point with negative epsilon */
                                    con->uneval();
                                    fk=con->eval(i);
                                    if((fk > active_tol && con->_ctype==leq) || (fk < -active_tol && con->_ctype==geq)){
                                        outer=true;
                                    }
                                    if(!outer){
                                        v->set_double_val(posv, std::min(std::max(xcurrent[count_var]*(1 + j*perturb_dist), lb_v),ub_v)); /** Perturbed point with positive epsilon */
                                        con->uneval();
                                        fk=con->eval(i);
                                        if((fk > active_tol && con->_ctype==leq) || (fk < -active_tol && con->_ctype==geq)){
                                            outer=true;
                                        }
                                    }
                                    if(outer){
                                        // DebugOn("outer"<<endl);
                                        auto con_interior=Ointerior.get_constraint(cname);
                                        xinterior=con_interior->get_x_ignore(i, "eta_interior"); /** ignore the Eta (slack) variable */
                                        xres=con->get_x(i);
                                        //                                            for(auto n=0;n<xres.size();n++){
                                        //                                                DebugOn(xinterior[n]<<"\t"<<xres[n]<<"\t");
                                        //                                            }
                                        // DebugOn(endl);
                                        auto res_search=con->binary_line_search(xinterior, i);
                                        if(res_search){
                                            oa_cut=true;
                                        }
                                    }
                                    if(oa_cut){
                                        //DebugOn("oacut"<<endl);
                                        convex_region=true;
                                        oa_cut=false;
                                        if(!con->is_convex() && !con->is_rotated_soc() && !con->check_soc()) //For the SDP determinant constraint, check if the point is feasible with repsecto to the SOC constraints
                                        {
                                            con->uneval();
                                            fkk=con->eval(i);
                                            xres=con->get_x(i);
                                            auto soc1=std::pow(xres[0],2)+std::pow(xres[3],2)-xres[6]*xres[7];
                                            auto soc2=std::pow(xres[1],2)+std::pow(xres[4],2)-xres[7]*xres[8];
                                            auto soc3=std::pow(xres[2],2)+std::pow(xres[5],2)-xres[6]*xres[8];
                                            if(soc1<=0 && soc2<=0 && soc3<=0){
                                                convex_region=true;
                                            }
                                            else{
                                                convex_region=false;
                                            }
                                        }
                                        if(convex_region){
                                            scale=1.0;
                                            if(add_new){
                                                nb_added_cuts++;
                                                indices activeset("active_"+con->_name);
                                                activeset.add((*con->_indices->_keys)[i]);
                                                Constraint<> OA_cut(con_lin_name);
                                                OA_cut=con->get_outer_app(activeset, scale);
                                                if(con->_ctype==leq) {
                                                    add(OA_cut.in(activeset)<=0);
                                                }
                                                else {
                                                    add(OA_cut.in(activeset)>=0);
                                                }
                                                add_new=false;
                                                oa_cut=false;
                                            }
                                            else{
                                                con->get_outer_coef(i, c_val, c0_val);
                                                vector<int> coefs;
                                                for (auto k = 0; k<c_val.size(); k++) {
                                                    if(c_val[k]!=0 && std::abs(c_val[k])<zero_tol){
                                                        scale=1.0e3;
                                                    }
                                                    coefs.push_back(1e5*c_val[k]);
                                                }
                                                coefs.push_back(1e5*c0_val);
                                                if(_OA_cuts[con->_id*100+i].insert(coefs).second)
                                                    oa_cut=true;
                                            }
                                        }
                                    }
                                    if(oa_cut){
                                        //DebugOn("adding inst"<<endl);
                                        nb_added_cuts++;
                                        auto con_lin=get_constraint("OA_cuts_"+con->_name);
                                        auto nb_inst = con_lin->get_nb_instances();
                                        con_lin->_indices->add("inst_"+to_string(nb_inst));
                                        con_lin->_dim[0] = con_lin->_indices->size();
                                        auto count=0;
                                        for(auto &l: *(con_lin->_lterms)){
                                            auto name=l.first;
                                            if(!l.second._sign){
                                                throw invalid_argument("symbolic negative");
                                            }
                                            if(l.second._coef->is_param()) {
                                                auto p_cst = ((param<>*)(l.second._coef.get()));
                                                p_cst->add_val("inst_"+to_string(p_cst->_indices->_keys->size()), c_val[count]*scale);
                                            }
                                            else {
                                                throw invalid_argument("Coefficient must be parameter");
                                            }
                                            auto parkeys=l.second._p->_indices->_keys;
                                            auto varc = con->get_var(l.second._p->_name);
                                            l.second._p->_indices->add_ref((*parkeys)[varc->get_id_inst(i)]);
                                            count++;
                                        }
                                        //Set value of the constant!!!
                                        if(con_lin->_cst->is_param()){
                                            auto co_cst = ((param<>*)(con_lin->_cst.get()));
                                            co_cst->add_val("inst_"+to_string(co_cst->_indices->_keys->size()), c0_val*scale);
                                        }
                                        else if(con_lin->_cst->is_function()){
                                            auto rhs_f = static_pointer_cast<func<>>(con_lin->_cst);
                                            if(!rhs_f->func_is_param()){
                                                throw invalid_argument("function should be a param");
                                            }
                                            auto p = static_pointer_cast<param<>>(rhs_f->_params->begin()->second.first);
                                            p->add_val("inst_"+to_string(p->_indices->_keys->size()), c0_val*scale);
                                            rhs_f->_indices->add("inst_"+to_string(nb_inst));
                                            rhs_f->_dim[0] = rhs_f->_indices->size();
                                        }
                                    }
                                    con->set_x(i, xcurrent);
                                    xres.clear();
                                }
                            }
                            count_var++;
                        }
                        xcurrent.clear();
                        xinterior.clear();
                    }
                }
            }
        }
        reindex();
        DebugOn("Number of constraints in linear model after perturb "<<_nb_cons<<endl);
        set_solution(xsolution);
        return Ointerior;
    }
    //
    template<>
    void Model<>::add_iterative(const Model<>& interior, vector<double>& obbt_solution, shared_ptr<Model<>> lin, string modelname, int& nb_oacuts, double active_tol)
    {
        
        vector<double> xsolution(_nb_vars);
        vector<double> xinterior(_nb_vars);
        vector<double> xcurrent, xres;
        get_solution(xsolution);
        set_solution(obbt_solution);
        const double active_tol_sol=1e-12, zero_tol=1e-6;
        
        bool interior_solv=true;
        vector<double> c_val ;
        double c0_val;
        bool oa_cut=true;
        bool convex_region=true;
        bool add_new=false;
        int nb_added_cuts = 0;
        string keyv;
        double scale=1.0;
        //    Ointerior.print();
        
        string vkname,keyk,dirk;
        var<> vk;
        shared_ptr<param_> vck;
        if(modelname!="allvar"){
            auto mname=modelname;
            std::size_t pos = mname.find("|");
            vkname.assign(mname, 0, pos);
            //  DebugOn("vkname "<<vkname<<endl);
            mname=mname.substr(pos+1);
            pos=mname.find("|");
            keyk.assign(mname, 0, pos);
            dirk=mname.substr(pos+1);
            vk=this->template get_var<double>(vkname);
        }
        bool var_found=false, key_found=false;
        for (auto &con: _cons_vec)
        {
            if(!con->is_linear()) {
                var_found=false;
                if(modelname=="allvar"){
                    var_found=true;
                }
                else if(con->has_sym_var(vk)){
                    var_found=true;
                }
                if(var_found){
                    auto cname=con->_name;
                    auto con_lin_name="OA_cuts_"+con->_name;
                    if(lin->_cons_name.find(con_lin_name)!=lin->_cons_name.end()){
                        add_new=false;
                    }
                    else{
                        add_new=true;
                    }
                    auto cnb_inst=con->get_nb_inst();
                    for(auto i=0;i<cnb_inst;i++){
                        for(auto &v: *con->_vars){
                            key_found=false;
                            if(modelname=="allvar"){
                                key_found=true;
                            }
                            else{
                                if(v.second.first->_vec_id==vk._vec_id){
                                    vck=con->_vars->at(v.first).first;
                                    if(!vck->is_indexed()){
                                        auto posv=vck->get_id_inst(i);
                                        keyv=(*vck->_indices->_keys)[posv];
                                    }
                                    else{
                                        auto posv=(vck->_indices->_ids->at(0))[i];
                                        keyv=(*vck->_indices->_keys)[posv];
                                    }
                                    if(keyv==keyk){
                                        key_found=true;
                                    }
                                }
                            }
                            // DebugOn(vkname<<""<<keyv<<" "<<keyk<<endl);
                            if(key_found){
                                oa_cut=false;
                                c0_val=0;
                                c_val.resize(con->_nb_vars,0);
                                auto cname=con->_name;
                                xcurrent=con->get_x(i);
                                con->uneval();
                                con->eval_all();
                                if(con->is_active(i,active_tol_sol)){
                                    oa_cut=true;
                                }
                                else
                                {   con->uneval();
                                    auto fk=con->eval(i);
                                    if((fk > active_tol && con->_ctype==leq) || (fk < -active_tol && con->_ctype==geq)){
                                        //if((!con->is_convex()||con->is_rotated_soc() || con->check_soc()) && (interior._status==0||interior._status==1))  {
                                        if((!con->is_convex()||con->is_rotated_soc() || con->check_soc()))  {
                                            auto con_interior=interior.get_constraint(cname);
                                            xinterior=con_interior->get_x_ignore(i, "eta_interior"); /** ignore the Eta (slack) variable */
                                            auto res_search=con->binary_line_search(xinterior, i);
                                            if(res_search){
                                                oa_cut=true;
                                            }
                                        }
                                        else{
                                            if (con->is_convex()){
                                                //  DebugOn(con->_name<<" "<<con->is_convex()<<" "<<con->is_rotated_soc()<<" "<<con->check_soc()<<endl);
                                                oa_cut=true;
                                            }
                                        }
                                        
                                    }
                                }
                                if(oa_cut){
                                    convex_region=true;
                                    oa_cut=false;
                                    
                                    if(!con->is_convex() && !con->is_rotated_soc() && !con->check_soc()) //For the SDP determinant constraint, check if the point is feasible with repsecto to the SOC constraints
                                    {
                                        
                                        xres=con->get_x(i);
                                        auto soc1=std::pow(xres[0],2)+std::pow(xres[3],2)-xres[6]*xres[7];
                                        auto soc2=std::pow(xres[1],2)+std::pow(xres[4],2)-xres[7]*xres[8];
                                        auto soc3=std::pow(xres[2],2)+std::pow(xres[5],2)-xres[6]*xres[8];
                                        if(soc1<=0 && soc2<=0 && soc3<=0){
                                            convex_region=true;
                                        }
                                        else{
                                            convex_region=false;
                                        }
                                    }
                                    if(convex_region){
                                        //                                                con->get_outer_coef(i, c_val, c0_val); /* Get the coefficients of the OA cut corresponding to instance i and store them in c_val and c0_val */
                                        //                                                for(auto l=0;l<c_val.size();l++)
                                        //                                                    DebugOn(c_val[l]<<"\t");
                                        //                                                DebugOn(c0_val<<endl);
                                        scale=1.0;
                                        if(add_new){
                                            nb_added_cuts++;
                                            indices activeset("active_"+con->_name);
                                            activeset.add((*con->_indices->_keys)[i]);
                                            Constraint<> OA_cut(con_lin_name);
                                            OA_cut=con->get_outer_app(activeset, scale);
                                            if(con->_ctype==leq) {
                                                lin->add(OA_cut.in(activeset)<=0);
                                            }
                                            else {
                                                lin->add(OA_cut.in(activeset)>=0);
                                            }
                                            add_new=false;
                                            oa_cut=false;
                                        }
                                        else{
                                            con->get_outer_coef(i, c_val, c0_val);
                                            vector<int> coefs;
                                            for (auto j = 0; j<c_val.size(); j++) {
                                                if(c_val[j]!=0 && std::abs(c_val[j])<zero_tol){
                                                    scale=1.0e3;
                                                }
                                                coefs.push_back(1e5*c_val[j]);
                                            }
                                            coefs.push_back(1e5*c0_val);
                                            if(_OA_cuts[con->_id*100+i].insert(coefs).second)
                                                oa_cut=true;
                                            //                                else {
                                            //                                    DebugOn("discarded OA cut");
                                            //                                }
                                        }
                                        //
                                        //                                                    Constraint<> con_oa(con->_name+to_string(i)+vname+to_string(j));
                                        //                                                    con_oa=con->get_outer_app_insti(i, false);
                                        //                                                    if(con->_ctype==geq)
                                        //                                                        add(con_oa>=0);
                                        //                                                    else
                                        //                                                        add(con_oa<=0);
                                    }
                                    
                                }
                                if(oa_cut){
                                    nb_added_cuts++;
                                    auto con_lin=lin->get_constraint("OA_cuts_"+con->_name);
                                    //                    DebugOn("added "<<con->_name<<endl);
                                    auto nb_inst = con_lin->get_nb_instances();
                                    con_lin->_indices->add("inst_"+to_string(nb_inst));
                                    con_lin->_dim[0] = con_lin->_indices->size();
                                    auto count=0;
                                    for(auto &l: *(con_lin->_lterms)){
                                        auto name=l.first;
                                        if(!l.second._sign){
                                            throw invalid_argument("symbolic negative");
                                        }
                                        if(l.second._coef->is_param()) {
                                            auto p_cst = ((param<>*)(l.second._coef.get()));
                                            //                                    DebugOn(p_cst->_indices->_keys->size());
                                            
                                            p_cst->add_val("inst_"+to_string(p_cst->_indices->_keys->size()), c_val[count]*scale);
                                            
                                            //                                    DebugOn(p_cst->_indices->_keys->size());
                                        }
                                        else {
                                            //                            auto f = static_pointer_cast<func<>>(l.second._coef);
                                            //                            if(!f->func_is_param()){
                                            //                                throw invalid_argument("function should be a param");
                                            //                            }
                                            //                            auto p = static_pointer_cast<param<>>(f->_params->begin()->second.first);
                                            //                            f->_indices->add("inst_"+to_string(p->_indices->_keys->size()));
                                            //                            p->add_val("inst_"+to_string(p->_indices->_keys->size()), c_val[count]);
                                            //                            f->_dim[0] = f->_indices->size();
                                            //                            f->uneval();
                                            //                            f->allocate_mem();
                                            throw invalid_argument("Coefficient must be parameter");
                                        }
                                        auto parkeys=l.second._p->_indices->_keys;
                                        //                                auto vname = l.second._p->_name.substr(0,l.second._p->_name.find_last_of("."));
                                        auto v = con->get_var(l.second._p->_name);
                                        l.second._p->_indices->add_ref((*parkeys)[v->get_id_inst(i)]);
                                        count++;
                                    }
                                    //Set value of the constant!!!
                                    if(con_lin->_cst->is_param()){
                                        auto co_cst = ((param<>*)(con_lin->_cst.get()));
                                        co_cst->add_val("inst_"+to_string(co_cst->_indices->_keys->size()), c0_val*scale);
                                    }
                                    else if(con_lin->_cst->is_function()){
                                        auto rhs_f = static_pointer_cast<func<>>(con_lin->_cst);
                                        if(!rhs_f->func_is_param()){
                                            throw invalid_argument("function should be a param");
                                        }
                                        auto p = static_pointer_cast<param<>>(rhs_f->_params->begin()->second.first);
                                        p->add_val("inst_"+to_string(p->_indices->_keys->size()), c0_val*scale);
                                        rhs_f->_indices->add("inst_"+to_string(nb_inst));
                                        rhs_f->_dim[0] = rhs_f->_indices->size();
                                    }
                                    //                            DebugOn("a"<<endl);
                                }
                                con->set_x(i, xcurrent);
                                xcurrent.clear();
                                xinterior.clear();
                                xres.clear();
                                break;
                            }
                        }
                    }
                }
            }
        }
        
        
        set_solution(xsolution);
        nb_oacuts+=nb_added_cuts;
        DebugOn("Number of constraints in linear model = " << nb_oacuts << endl);
        //OA_iter.print();
    }
    
    
    
    
    
    /** Returns an interior point of a model
     @param[in] nonlin: model for which interior point with respect to nonlinear constraints is computed
     Assuming model has no nonlinear equality constraints
     **/
    template<typename type>
    template<typename T>
    shared_ptr<Model<type>> Model<type>::build_model_IIS()
    {
        auto IIS = make_shared<Model<>>(_name+"IIS");
        int i, cs;
        
        for (auto &it: _vars)
        {
            auto v = it.second;
            if(!IIS->has_var(*v)){
                IIS->add_var(v);
            }
        }
        var<> eta_int("eta_int", -100, 100);
        
        IIS->add(eta_int.in(range(0,0)));
        auto obj=eta_int;
        
        IIS->min(obj);
        
        int n=0;
        for (auto &con: _cons_vec)
        {
            n+=con->get_nb_instances();
        }
        
        
        
        // DebugOn("N"<<n<<endl);
        // param<> eta_length("eta_length");
        //  eta_length.set_val(n);
        
        
        var<> eta_i("eta_i", 0, 1);
        IIS->add(eta_i.in(range(0,n-1)));
        
        Constraint<> sum_eta("sum_eta");
        sum_eta=eta_int-sum(eta_i);
        IIS->add(sum_eta==0);
        
        
        
        int counter=0;
        
        for (auto &con: _cons_vec)
        {
            cs=con->get_nb_instances();
            //            for(i=0;i<cs;i++)
            //            {
            
            // DebugOn(cs<<endl);
            Constraint<> Inter_con(con->_name);
            Inter_con=*con;
            
            if(con->_ctype==leq)
            {
                //  Inter_con -= eta_i.in(range(counter, counter));
                //IIS->add(Inter_con.in(range(0,con->get_nb_instances()-1))<=eta_i.in(range(counter, counter)));
                IIS->add(Inter_con<= eta_i.in(range(counter, counter+con->get_nb_instances()-1)));
            }
            else  if(con->_ctype==geq)
            {
                Inter_con += eta_i.in(range(counter, counter+con->get_nb_instances()-1));
                //IIS->add(Inter_con.in(range(0,con->get_nb_instances()-1))>=eta_i.in(range(counter, counter)));
                //   Inter_con=(*con)+eta_i.in(range(counter, counter));
                IIS->add(Inter_con>=0);
                //IIS->add(Inter_con.in(range(i,i))>=0);
            }
            else  if(con->_ctype==eq)
            {
                // IIS->add(Inter_con.in(range(0,con->get_nb_instances()-1))==0);
                Inter_con=(*con);
                IIS->add(Inter_con==0);
            }
            counter+=con->get_nb_instances();
            //eta_i.in(range(counter, counter))}
            
            
        }
        IIS->print();
        return IIS;
    }
#ifdef USE_MPI
    
    /** Runs models stored in the vector in parallel using MPI
     *      @models vector of models to run in parallel
     *           @stype Solver type
     *                @tol numerical tolerance
     *                     @max_iter max number of iterations per model
     *                          @max_batch_time max wall clock time of each batch
     *                               @nb_threads Number of parallel threads per worker
     *                                    @lin_solver linear system solver
     *                                         @share_all propagate model status and solutions to all workers, if false, only worker 0 has updated solutions and status flags for all models
     *                                              @share_all_obj propagate only objective values and model status to all workers
     *
     */
    int run_MPI_new(const std::vector<std::string> objective_models, std::vector<double>& sol_obj, std::vector<int>& sol_status, const vector<shared_ptr<gravity::Model<double>>>& models, gravity::SolverType stype, double tol, unsigned nr_threads, const string& lin_solver, int max_iter, int max_batch_time, bool share_all, bool share_all_obj){
        int worker_id, nb_workers;
        auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
        auto err_size = MPI_Comm_size(MPI_COMM_WORLD, &nb_workers);
        auto nb_workers_ = std::min((size_t)nb_workers, objective_models.size());
        std::string msname, mname, vname, key, dir;
        var<> var;
        if(objective_models.size()!=0){
            /* Split models into equal loads */
            auto nb_total_threads_ = std::min((size_t)nr_threads*nb_workers, objective_models.size());
            auto nb_threads_per_worker = std::min((size_t)nr_threads, objective_models.size());
            DebugOff("I have " << nb_workers_ << " workers" << endl);
            DebugOff("I will be using  " << nb_total_threads_ << " thread(s) in total" << endl);
            std::vector<size_t> limits = bounds(nb_workers_, objective_models.size());
            DebugOff("I will be splitting " << objective_models.size() << " tasks ");
            DebugOff("among " << nb_workers_ << " worker(s)" << endl);
            DebugOff("limits size = " << limits.size() << endl);
            for (size_t i = 0; i < limits.size(); ++i) {
                DebugOff("limits[" << i << "] = " << limits[i] << endl);
            }
            if(worker_id+1<limits.size()){
                /* Launch all threads in parallel */
                if(limits[worker_id] == limits[worker_id+1]){
                    throw invalid_argument("limits[worker_id]==limits[worker_id+1]");
                }
                DebugOff("I'm worker ID: " << worker_id << ", I will be running models " << limits[worker_id] << " to " << limits[worker_id+1]-1 << endl);
                int count=0;
                auto vec = vector<shared_ptr<gravity::Model<double>>>();
                for (auto i = limits[worker_id]; i < limits[worker_id+1]; i++) {
                    msname=objective_models[i];
                    mname=msname;
                    std::size_t pos = msname.find("|");
                    vname.assign(msname, 0, pos);
                    msname=msname.substr(pos+1);
                    pos=msname.find("|");
                    key.assign(msname, 0, pos);
                    dir=msname.substr(pos+1);
                    var=models[count]->get_var<double>(vname);
                    if(dir=="LB")
                    {
                        models[count]->min(var(key));
                    }
                    else
                    {
                        models[count]->max(var(key));
                        
                    }
                    models[count]->set_name(mname);
                    models[count]->reindex();
                    vec.push_back(models[count++]);
                }
                run_parallel(vec,stype,tol,nr_threads,lin_solver,max_iter);
            }
            sol_status.resize(objective_models.size(),0);
            sol_obj.resize(objective_models.size(),0);
            MPI_Barrier(MPI_COMM_WORLD);
            send_status_new(models,limits, sol_status);
            MPI_Barrier(MPI_COMM_WORLD);
            if(share_all_obj){
                /* We will send the objective value of successful models */
                send_obj_all_new(models,limits, sol_obj);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        return max(err_rank, err_size);
        
    }
    
    int run_MPI(const vector<shared_ptr<gravity::Model<double>>>& models, gravity::SolverType stype, double tol, unsigned nr_threads, const string& lin_solver, int max_iter, int max_batch_time, bool share_all, bool share_all_obj){
        int worker_id, nb_workers;
        auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
        auto err_size = MPI_Comm_size(MPI_COMM_WORLD, &nb_workers);
        auto nb_workers_ = std::min((size_t)nb_workers, models.size());
        MPI_Request send_reqs[nb_workers_*models.size()];
        
        if(models.size()!=0){
            /* Split models into equal loads */
            auto nb_total_threads_ = std::min((size_t)nr_threads*nb_workers, models.size());
            auto nb_threads_per_worker = std::min((size_t)nr_threads, models.size());
            DebugOff("I have " << nb_workers_ << " workers" << endl);
            DebugOff("I will be using  " << nb_total_threads_ << " thread(s) in total" << endl);
            std::vector<size_t> limits = bounds(nb_workers_, models.size());
            DebugOff("I will be splitting " << models.size() << " tasks ");
            DebugOff("among " << nb_workers_ << " worker(s)" << endl);
            DebugOff("limits size = " << limits.size() << endl);
            for (size_t i = 0; i < limits.size(); ++i) {
                DebugOff("limits[" << i << "] = " << limits[i] << endl);
            }
            if(worker_id+1<limits.size()){
                /* Launch all threads in parallel */
                if(limits[worker_id] == limits[worker_id+1]){
                    throw invalid_argument("limits[worker_id]==limits[worker_id+1]");
                }
                DebugOff("I'm worker ID: " << worker_id << ", I will be running models " << limits[worker_id] << " to " << limits[worker_id+1]-1 << endl);
                auto vec = vector<shared_ptr<gravity::Model<double>>>();
                for (auto i = limits[worker_id]; i < limits[worker_id+1]; i++) {
                    vec.push_back(models[i]);
                }
                run_parallel(vec,stype,tol,nr_threads,lin_solver,max_iter);
            }
            send_status(models,limits);
            MPI_Barrier(MPI_COMM_WORLD);
            if(!share_all && !share_all_obj){/* Only share solution with worker 0 */
                if (worker_id == 0){
                    DebugOff("I'm the main worker, I'm waiting for the solutions broadcasted by the other workers " << endl);
                    for (auto w_id = 1; w_id<nb_workers_; w_id++) {
                        for (auto i = limits[w_id]; i < limits[w_id+1]; i++) {
                            auto model = models[i];
                            auto nb_vars = model->get_nb_vars();
                            vector<double> solution;
                            solution.resize(nb_vars);
                            DebugOff("I'm the main worker, I'm waiting for the solution of task " << i << " broadcasted by worker " << w_id << endl);
                            MPI_Recv(&solution[0], nb_vars, MPI_DOUBLE, w_id, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                            DebugOff("I'm the main worker, I received the solution of task " << i << " broadcasted by worker " << w_id << endl);
                            model->set_solution(solution);
                        }
                    }
                }
                else {
                    DebugOff("I'm worker ID: " << worker_id << ", I will be sending my solutions to main worker " << endl);
                    for (auto i = limits[worker_id]; i < limits[worker_id+1]; i++) {
                        auto model = models[i];
                        auto nb_vars = model->get_nb_vars();
                        vector<double> solution;
                        solution.resize(nb_vars);
                        model->get_solution(solution);
                        DebugOff("I'm worker ID: " << worker_id << ", I finished loading solution of task " << i << endl);
                        MPI_Send(&solution[0], nb_vars, MPI_DOUBLE, 0, i, MPI_COMM_WORLD);
                        DebugOff("I'm worker ID: " << worker_id << ", I finished sending solution of task " << i << endl);
                    }
                }
            }
            else if(share_all){
                /* We will send the solution of successful models */
                send_solution_all(models,limits);
            }
            if(share_all_obj){
                /* We will send the objective value of successful models */
                send_obj_all(models,limits);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        return max(err_rank, err_size);
    }
    void run_MPI(const initializer_list<shared_ptr<gravity::Model<double>>>& models, gravity::SolverType stype, double tol, unsigned nr_threads, const string& lin_solver, int max_iter, int max_batch_time, bool share_all, bool share_all_obj){
        run_MPI(vector<shared_ptr<gravity::Model<double>>>(models), stype, tol, nr_threads, lin_solver);}
#endif
    template shared_ptr<Model<double>> Model<double>::buildOA();
    template Model<double> Model<double>::build_model_interior() const;
    template shared_ptr<Model<double>> Model<double>::build_model_IIS();
}
