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
/** Returns copy of current Model  which has all variables in current model, same objective as current model and only linear constraints in current Model. Throws exception if model has nonlinear equality constraints
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
            Constraint<type> temp_c;
            temp_c.deep_copy(*con);
            OA->add(temp_c);
            auto c=OA->get_constraint(con->_name);
            c->_new=true;
            for(auto i=0;i<c->get_nb_inst();i++){
                c->_violated[i]=true;
            }
        }
    }
    set_solution(xsolution);
    return OA;
}
/** Outer approximation of model. Throws exception if model has nonlinear equality constraints
 @param[in] nb_discr: number of OA cuts per nonlinear constraint
 @return Model with OA cuts. OA cuts are added to the model (for all func instances) in an uniform grid (nb_discr)
 **/
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


/** Returns a model which has all variables in current model. If current model has nonlinear constraints g(x), a new model is created with g_i(x) \le eta_i and the objective is to
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



/* Runs models stored in the vector in parallel, using solver of stype and tolerance tol */
int run_parallel_new(const std::vector<std::string> objective_models, std::vector<double>& sol_obj, std::vector<int>& sol_status, const std::vector<shared_ptr<gravity::Model<double>>>& models, gravity::SolverType stype, double tol, unsigned nr_threads, const string& lin_solver, int max_iter, int max_batch_time){
    std::vector<thread> threads;
    std::vector<bool> feasible;
    std::vector<double> solution(models[0]->_nb_vars);
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
        models[s]->set_name(mname);
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
        threads.push_back(thread(run_models<double>, ref(vec), limits[i], limits[i+1], stype, tol, lin_solver, max_iter, max_batch_time));
    }
    /* Join the threads with the main thread */
    for(auto &t : threads){
        t.join();
    }
    int count=0;
    for(auto &m:models){
        if(count<objective_models.size()){
            sol_status.at(count)=m->_status;
            sol_obj.at(count)=m->get_obj_val();
            count++;
        }
    }
    return 0;
}

/* Runs models stored in the vector models in parallel, using the solver object solvers. Objective of models to run is given in objective_models. Solution status and objective of each model is added to sol_status, and sol_obj. Liner constraints are added to each model calling add_iterative.  relaxed_model, interior, cut_type,  nb_oacuts, active_tol are args to add_iterative*/
int run_parallel_new(const std::vector<std::string> objective_models, std::vector<double>& sol_obj, std::vector<int>& sol_status, std::vector<shared_ptr<gravity::Model<double>>>& models, const shared_ptr<gravity::Model<double>>& relaxed_model, const gravity::Model<double>& interior, string cut_type, double active_tol, gravity::SolverType stype, double tol, unsigned nr_threads, const string& lin_solver, int max_iter, int max_batch_time, bool linearize, int nb_refine){
    std::vector<thread> threads;
    std::vector<bool> feasible;
    std::vector<double> solution(models[0]->_nb_vars);
    std::string mname, msname,vname, key, dir, modelname;
    var<> var;
    vector<shared_ptr<solver<double>>> batch_solvers;
    int viol_i=0, viol=0, count=0, ncuts=0;
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
        models[s]->set_name(mname);
        if(dir=="LB")
        {
            models[s]->min(var(key));
            models[s]->_objt=minimize;
        }
        else
        {
            models[s]->max(var(key));
            models[s]->_objt=maximize;
        }
        models[s]->reindex();
        models[s]->reset();
        models[s]->reset_lifted_vars_bounds();
        models[s]->reset_constrs();
        DebugOff(objective_models[s]);
        //models[s]->print();
        auto solverk = make_shared<solver<double>>(models[s], stype);
        batch_solvers.push_back(solverk);
    }
    DebugOn(endl);
    /* Split models into nr_threads parts */
    auto nr_threads_ = std::min((size_t)nr_threads,objective_models.size());
    std::vector<size_t> limits = bounds(nr_threads_, objective_models.size());
    DebugOn("Running on " << nr_threads_ << " threads." << endl);
    DebugOn("limits size = " << limits.size() << endl);
    for (size_t i = 0; i < limits.size(); ++i) {
        DebugOff("limits[" << i << "] = " << limits[i] << endl);
    }
    /* Launch all threads in parallel */
    for(auto l=0;l<nb_refine;l++){
        viol_i=0;
        viol=0;
        count=0;
        if(l>=1){
            DebugOn("resolving"<<endl);
        }
        if(l>=1 && stype==ipopt){
            for (auto s=0;s<objective_models.size();s++){
                if(batch_solvers[s]->_model->_objt==maximize){
                    *batch_solvers[s]->_model->_obj *= -1;
                }
            }
            
        }
        auto vec = vector<shared_ptr<gravity::Model<double>>>(models);
        for (size_t i = 0; i < nr_threads_; ++i) {
            threads.push_back(thread(run_models_solver<double>, ref(vec), ref(batch_solvers), limits[i], limits[i+1], stype, tol, lin_solver, max_iter, max_batch_time));
        }
        /* Join the threads with the main thread */
        for(auto &t : threads){
            t.join();
        }
        if(linearize){
            for(auto &m:models){
                if(count<objective_models.size()){
                    m->get_solution(solution);
                    if(m->_status==0 && linearize){
                        if(cut_type=="modelname"){
                            modelname=m->get_name();
                        }
                        else if(cut_type=="allvar"){
                            modelname="allvar";
                        }
                        DebugOff("active tol "<<active_tol<<endl);
                        viol_i=relaxed_model->add_iterative(interior,  solution, models[count],  modelname,  ncuts,  active_tol);
                    }
                    if(viol_i==1){
                        viol=1;
                    }
                    count++;
                    m->reindex();
                    m->reset();
                    m->reset_constrs();
                    // m->reset_lifted_vars_bounds();
                }
            }
        }
        threads.clear();
        if(viol==0){
            break;
        }
    }
    count=0;
    for(auto &m:models){
        if(count<objective_models.size()){
            sol_status.at(count)=m->_status;
            sol_obj.at(count)=m->get_obj_val();
            count++;
        }
    }
    return viol;
}

/* Runds models stored in the vector in parallel, using solver of stype and tolerance tol */
int run_parallel(const vector<shared_ptr<gravity::Model<double>>>& models, gravity::SolverType stype, double tol, unsigned nr_threads, const string& lin_solver, int max_iter, int max_batch_time){
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
        threads.push_back(thread(run_models<double>, ref(vec), limits[i], limits[i+1], stype, tol, lin_solver, max_iter, max_batch_time));
    }
    /* Join the threads with the main thread */
    for(auto &t : threads){
        t.join();
    }
    return 0;
}
/* Runds models stored in the vector in parallel, using solver solvers. solvers is passed as an argument ot run_model_solvers */
//int run_parallel(const vector<shared_ptr<gravity::Model<double>>>& models, const vector<shared_ptr<solver<>>>& solvers, gravity::SolverType stype, double tol, unsigned nr_threads, const string& lin_solver, int max_iter, int max_batch_time){
//    std::vector<thread> threads;
//    std::vector<bool> feasible;
//    if(models.size()==0){
//        DebugOff("in run_parallel(models...), models is empty, returning");
//        return -1;
//    }
//    /* Split models into nr_threads parts */
//    auto nr_threads_ = std::min((size_t)nr_threads,models.size());
//    std::vector<size_t> limits = bounds(nr_threads_, models.size());
//    DebugOff("Running on " << nr_threads_ << " threads." << endl);
//    DebugOff("limits size = " << limits.size() << endl);
//    for (size_t i = 0; i < limits.size(); ++i) {
//        DebugOff("limits[" << i << "] = " << limits[i] << endl);
//    }
//    /* Launch all threads in parallel */
//    auto vec = vector<shared_ptr<gravity::Model<double>>>(models);
//    for (size_t i = 0; i < nr_threads_; ++i) {
//        threads.push_back(thread(run_models_solver<double>, ref(vec), solvers, limits[i], limits[i+1], stype, tol, lin_solver, max_iter, max_batch_time));
//    }
//    /* Join the threads with the main thread */
//    for(auto &t : threads){
//        t.join();
//    }
//    return 0;
//}


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
/*To add OA constraints using a grid, not used in current linearization approach*/
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
/*Returns true if constraint cname of instance int in the current model is violated by x, as per active_tol. Filtering approach, currently not used.*/
template<>
bool Model<>::linearmodel_violates_x(vector<double>& x, string cname, int inst, double tol){
#ifdef USE_MPI
    int worker_id, nb_workers;
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
    auto err_size = MPI_Comm_size(MPI_COMM_WORLD, &nb_workers);
#endif
    bool violated=false;
    double active_tol=1e-6;
    if(tol<active_tol){
        // active_tol=tol;
    }
    vector<double> xcurr;
    string lin_cname="OA_cuts_"+cname;
    int nb_inst;
    if(_cons_name.find(lin_cname)!=_cons_name.end()){
        auto con_lin=_cons_name.at(lin_cname);
        nb_inst=con_lin->get_nb_instances();
        for(auto i=0;i<nb_inst;i++){
            if(con_lin->_lazy[i]){
                if(con_lin->_indices->_keys->at(i).find("inst_"+to_string(inst)+"_")!=string::npos){
                    xcurr=con_lin->get_x(i);
                    con_lin->set_x(i, x);
                    con_lin->uneval();
                    auto fk=con_lin->eval(i);
                    con_lin->set_x(i, xcurr);
                    if((fk > active_tol && con_lin->_ctype==leq) || (fk < -active_tol && con_lin->_ctype==geq)){
                        violated=true;
#ifdef USE_MPI
                        if(worker_id==0)
                            DebugOn("avoided fcut"<<endl);
#endif
                        return violated;
                    }
                }
            }
        }
    }
    return violated;
}
/*Returns a model that can compute interior point of nonlin model. Generates OA cuts to nonlin model at solution of nonlin model. Adds the generated cuts to the current model. If nb_perturb>0, perturbs the solution of nonlin model and generates cuts at the perturbed solution too*/
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
    //modelI.set_option("bound_relax_factor", tol*1e-2);
    //modelI.set_option("check_violation", true);
    modelI.run(output=0, tol);
    
    
    vector<double> xsolution(_nb_vars);
    nonlin.get_solution(xsolution);
    
    //    Ointerior.print();
    
    if((Ointerior._status==0||Ointerior._status==1) && Ointerior.get_obj_val() <0)
    {
        interior=true;
        //        Ointerior.print_solution();
    }
    else{
        DebugOn("Failed to solve interior model"<<endl);
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
                auto oa_sol_inst=OA_sol.get_nb_instances();
                for(auto o=0;o<oa_sol_inst;o++){
                    OA_sol._lazy.push_back(false);
                }
                
                if(con->_ctype==leq) {
                    add(OA_sol.in(activeset)<=0);
                }
                else {
                    add(OA_sol.in(activeset)>=0);
                }
                //OA_sol.print();
                if(oa_sol_inst>0){
                    auto coa=get_constraint(OA_sol._name);
                    for(auto o=0;o<oa_sol_inst;o++){
                        coa->_violated[o]=true;
                    }
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
                auto oa_sol_inst=OA_sol.get_nb_instances();
                for(auto o=0;o<oa_sol_inst;o++){
                    OA_sol._lazy.push_back(false);
                }
                if(con->_ctype==leq) {
                    add(OA_sol.in(allset)<=0);
                }
                else {
                    add(OA_sol.in(allset)>=0);
                }
                if(oa_sol_inst>0){
                    auto coa=get_constraint(OA_sol._name);
                    for(auto o=0;o<oa_sol_inst;o++){
                        coa->_violated[o]=true;
                    }
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
                auto oa_sol_inst=OA_sol.get_nb_instances();
                for(auto o=0;o<oa_sol_inst;o++){
                    OA_sol._lazy.push_back(false);
                }
                if(con->_ctype==leq) {
                    add(OA_sol.in(allset)<=0);
                }
                else {
                    add(OA_sol.in(allset)>=0);
                }
                if(oa_sol_inst>0){
                    auto coa=get_constraint(OA_sol._name);
                    for(auto o=0;o<oa_sol_inst;o++){
                        coa->_violated[o]=true;
                    }
                }
            }
        }
    }
    nonlin.set_solution(xsolution);
    reindex();
    DebugOff("Number of constraints in linear model "<<_nb_cons<<endl);
    bool add_new, oa_cut;
    int nb_added_cuts=0;
    int nb_perturb=0;
    int count_var=0;
    bool near_zero=false;
    if(interior && nb_perturb>0){
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
                                    if(con->xval_within_bounds(i, xcurrent)){
                                        if(!(linearmodel_violates_x(xcurrent, con->_name, i, active_tol))){
                                            auto res_search=con->binary_line_search(xinterior, i);
                                            if(res_search){
                                                oa_cut=true;
                                            }
                                        }
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
                                            OA_cut._lazy.push_back(false);
                                            if(con->_ctype==leq) {
                                                add(OA_cut.in(activeset)<=0);
                                            }
                                            else {
                                                add(OA_cut.in(activeset)>=0);
                                            }
                                            auto coa=get_constraint(OA_cut._name);
                                            for(auto o=0;o<coa->get_nb_inst();o++){
                                                coa->_violated[o]=true;
                                            }
                                            add_new=false;
                                            oa_cut=false;
                                        }
                                        else{
                                            con->get_outer_coef(i, c_val, c0_val);
                                            vector<int> coefs;
                                            scale=1;
                                            for (auto k = 0; k<c_val.size(); k++) {
                                                near_zero=true;
                                                scale=1.0;
                                                if(c_val[k]!=0 && std::abs(c_val[k])<zero_tol){
                                                    if(zero_tol/std::abs(c_val[k])>scale){
                                                        scale=zero_tol/std::abs(c_val[k]);
                                                    }
                                                }
                                                if(near_zero && c_val[k]!=0 && std::abs(c_val[k])<zero_tol){
                                                    near_zero=true;
                                                }
                                                else{
                                                    near_zero=false;
                                                }
                                                //coefs.push_back(1e5*c_val[j]);
                                            }
                                            // coefs.push_back(1e5*c0_val);
                                            // if(_OA_cuts[con->_id*100+i].insert(coefs).second)
                                            oa_cut=true;
                                            if(near_zero)
                                                oa_cut=false;
                                            
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
                                    con_lin->_lazy.push_back(false);
                                    con_lin->_violated.push_back(true);
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
                                        con_lin->_cst=p;
                                    }
                                    con_lin->eval_all();
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
    DebugOff("Number of constraints in linear model after perturb "<<_nb_cons<<endl);
    set_solution(xsolution);
    //print();
    return Ointerior;
}
/*Generates and adds constraints to model lin. The curent model solution is set to obbt_solution and OA cuts are generated for the current model at the obbt_solution point. These cuts are added to model lin. The cut addition strategy is given by argument modelname. Can be for "allvar" or only "modelname", for varibale whose name is given in the models->_name. When a cut is added nb_oacuts is incremented. active_tol is the tolerance used to generate the cuts. */
template<typename type>
template<typename T>
bool Model<type>::add_iterative(const Model<type>& interior, vector<double>& obbt_solution, shared_ptr<Model<type>>& lin, string modelname, int& nb_oacuts, double active_tol)
{
    vector<double> xsolution(_nb_vars);
    vector<double> xinterior(_nb_vars);
    vector<double> xcurrent, xres;
    get_solution(xsolution);
    set_solution(obbt_solution);
    const double active_tol_sol=1e-12, zero_tol=1e-6;
    bool constr_viol=false;
    bool interior_solv=true;
    vector<double> c_val ;
    double c0_val;
    bool oa_cut=true;
    bool convex_region=true;
    bool add_new=false;
    int nb_added_cuts = 0;
    string keyv;
    double scale=1.0;
    bool near_zero=true;
    string cname, con_lin_name;
    
    string vkname,keyk,dirk;
    var<> vk;
    shared_ptr<param_> vck;
    if(modelname!="allvar"){
        auto mname=modelname;
        std::size_t pos = mname.find("|");
        vkname.assign(mname, 0, pos);
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
                cname=con->_name;
                con_lin_name="OA_cuts_"+con->_name;
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
                                    constr_viol=true;
                                    //ToDo fix interior status
                                    DebugOff("status "<< interior._status);
                                    //                                        if((!con->is_convex()||con->is_rotated_soc() || con->check_soc()) && (interior._status==0||interior._status==1))  {
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
                                    scale=1.0;
                                    if(add_new){
                                        nb_added_cuts++;
                                        indices activeset("active_"+con->_name);
                                        activeset.add((*con->_indices->_keys)[i]);
                                        Constraint<> OA_cut(con_lin_name);
                                        OA_cut=con->get_outer_app(activeset, scale);
                                        //OA_cut._lazy.push_back(false);
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
                                            near_zero=true;
                                            scale=1.0;
                                            if(c_val[j]!=0 && std::abs(c_val[j])<zero_tol){
                                                if(zero_tol/std::abs(c_val[j])>scale){
                                                    scale=zero_tol/std::abs(c_val[j]);
                                                }
                                            }
                                            if(near_zero && c_val[j]!=0 && std::abs(c_val[j])<zero_tol){
                                                near_zero=true;
                                            }
                                            else{
                                                near_zero=false;
                                            }
                                            //coefs.push_back(1e5*c_val[j]);
                                        }
                                        // coefs.push_back(1e5*c0_val);
                                        // if(_OA_cuts[con->_id*100+i].insert(coefs).second)
                                        oa_cut=true;
                                        if(near_zero)
                                            oa_cut=false;
                                    }
                                }
                                
                            }
                            if(oa_cut){
                                nb_added_cuts++;
                                auto con_lin=lin->get_constraint("OA_cuts_"+con->_name);
                                auto nb_inst = con_lin->get_nb_instances();
                                con_lin->_indices->add("inst_"+to_string(i)+"_"+to_string(nb_inst));
                                con_lin->_dim[0] = con_lin->_indices->_keys->size();
                                //con_lin->_lazy.push_back(true);
                                con_lin->_violated.push_back(true);
                                auto count=0;
                                DebugOff("nb inst "<<nb_inst);
                                //auto func_true=false;
                                for(auto &l: *(con_lin->_lterms)){
                                    auto name=l.first;
                                    if(!l.second._sign){
                                        throw invalid_argument("symbolic negative");
                                    }
                                    if(l.second._coef->is_param()) {
                                        auto p_cst = ((param<>*)(l.second._coef.get()));
                                        p_cst->add_val("inst_"+to_string(i)+"_"+to_string(nb_inst), c_val[count]*scale);
                                        DebugOff("added p"<<endl);
                                    }
                                    else {
                                        auto f = static_pointer_cast<func<>>(l.second._coef);
                                        if(!f->func_is_param()){
                                            throw invalid_argument("function should be a param");
                                        }
                                        auto p = static_pointer_cast<param<>>(f->_params->begin()->second.first);
                                        p->add_val("inst_"+to_string(i)+"_"+to_string(nb_inst), c_val[count]*scale);
                                        l.second._coef = p;
                                    }
                                    auto parkeys=l.second._p->_indices->_keys;
                                    //                                auto vname = l.second._p->_name.substr(0,l.second._p->_name.find_last_of("."));
                                    auto v = con->get_var(l.second._p->_name);
                                    l.second._p->_indices->add_ref((*parkeys)[v->get_id_inst(i)]);
                                    count++;
                                }
                                //Set value of the constant
                                if(con_lin->_cst->is_param()){
                                    auto co_cst = ((param<>*)(con_lin->_cst.get()));
                                    co_cst->add_val("inst_"+to_string(i)+"_"+to_string(nb_inst), c0_val*scale);
                                }
                                else if(con_lin->_cst->is_function()){
                                    auto rhs_f = static_pointer_cast<func<>>(con_lin->_cst);
                                    if(!rhs_f->func_is_param()){
                                        throw invalid_argument("function should be a param");
                                    }
                                    auto rhs_p = static_pointer_cast<param<>>(rhs_f->_params->begin()->second.first);
                                    rhs_p->add_val("inst_"+to_string(i)+"_"+to_string(nb_inst), c0_val*scale);
                                    con_lin->_cst = rhs_p;
                                }
                                
                                
                                con_lin->uneval();
                                con_lin->eval_all();
                                //con_lin->eval(nb_inst);
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
    lin->set_solution(obbt_solution);
    nb_oacuts+=nb_added_cuts;
    DebugOff("Number of constraints in linear model = " << nb_oacuts << endl);
    return(constr_viol);
}
/*Generates cuts (just as using the strategy in the add_iterative function above) but does not add them to lin model. Stores new row of existing constraints in lin in vector cuts. If a constraint does not exist in lin, point at which a new constraint must be added is stored in vector cuts*/
template<>
int Model<>::generate_cuts_iterative(const Model<>& interior, vector<double>& obbt_solution, shared_ptr<Model<>> lin, string modelname, int& nb_oacuts, double active_tol, vector<double>& cuts)
{
    vector<double> xsolution(_nb_vars);
    vector<double> xinterior(_nb_vars);
    vector<double> xcurrent, xres;
#ifdef USE_MPI
    int worker_id, nb_workers;
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
    auto err_size = MPI_Comm_size(MPI_COMM_WORLD, &nb_workers);
#endif
    get_solution(xsolution);
    set_solution(obbt_solution);
    const double active_tol_sol=1e-12, zero_tol=1e-6;
    int constr_viol=0;
    bool interior_solv=true;
    vector<double> c_val ;
    double c0_val;
    bool oa_cut=true;
    bool convex_region=true;
    bool add_new=false;
    int nb_added_cuts = 0;
    string keyv;
    double scale=1.0;
    bool near_zero=true;
    string cname, con_lin_name;
    
    string vkname,keyk,dirk;
    var<> vk;
    shared_ptr<param_> vck;
    //        inst_old.clear();
    //        for(auto &clin:lin->_cons_vec){
    //            inst_old[clin->_name]=clin->get_nb_instances();
    //        }
    if(modelname!="allvar"){
        auto mname=modelname;
        std::size_t pos = mname.find("|");
        vkname.assign(mname, 0, pos);
        mname=mname.substr(pos+1);
        pos=mname.find("|");
        keyk.assign(mname, 0, pos);
        dirk=mname.substr(pos+1);
        vk=this->template get_var<double>(vkname);
    }
    bool var_found=false, key_found=false;
    int con_pos=0;
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
                cname=con->_name;
                con_lin_name="OA_cuts_"+con->_name;
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
                                    constr_viol=1;
                                    //ToDo fix interior status
                                    DebugOff("status "<< interior._status);
                                    //                                        if((!con->is_convex()||con->is_rotated_soc() || con->check_soc()) && (interior._status==0||interior._status==1))  {
                                    if((!con->is_convex()||con->is_rotated_soc() || con->check_soc()))  {
                                        if(true){
                                            if(true){
                                                auto con_interior=interior.get_constraint(cname);
                                                xinterior=con_interior->get_x_ignore(i, "eta_interior"); /** ignore the Eta (slack) variable */
                                                auto res_search=con->binary_line_search(xinterior, i);
                                                if(res_search){
                                                    oa_cut=true;
                                                }
                                            }
                                        }
                                        else{
#ifdef USE_MPI
                                            if(worker_id==0)
                                                DebugOn("avoided xcut "<<endl);
#endif
                                        }
                                        
                                    }
                                    else{
                                        if (con->is_convex()){
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
                                    scale=1.0;
                                    if(add_new){
                                        oa_cut=true;
                                        
                                    }
                                    else{
                                        con->get_outer_coef(i, c_val, c0_val);
                                        vector<int> coefs;
                                        for (auto j = 0; j<c_val.size(); j++) {
                                            near_zero=true;
                                            scale=1.0;
                                            if(c_val[j]!=0 && std::abs(c_val[j])<zero_tol){
                                                if(zero_tol/std::abs(c_val[j])>scale){
                                                    scale=zero_tol/std::abs(c_val[j]);
                                                }
                                            }
                                            if(near_zero && c_val[j]!=0 && std::abs(c_val[j])<zero_tol){
                                                near_zero=true;
                                            }
                                            else{
                                                near_zero=false;
                                            }
                                            //coefs.push_back(1e5*c_val[j]);
                                        }
                                        // coefs.push_back(1e5*c0_val);
                                        // if(_OA_cuts[con->_id*100+i].insert(coefs).second)
                                        oa_cut=true;
                                        if(near_zero)
                                            oa_cut=false;
                                    }
                                }
                            }
                            if(oa_cut && !add_new){
                                nb_added_cuts++;
                                auto con_lin=lin->get_constraint("OA_cuts_"+con->_name);
                                auto count_con=0;
                                auto not_found =true;
                                for(auto c:lin->_cons_vec){
                                    if(c->_name=="OA_cuts_"+con->_name){
                                        not_found=false;
                                        break;
                                    }
                                    count_con++;
                                }
                                if(!not_found){
                                    cuts.push_back(count_con);
                                    cuts.push_back(i);
                                    for (auto j = 0; j<c_val.size(); j++) {
                                        cuts.push_back(c_val[j]*scale);
                                    }
                                    cuts.push_back(c0_val*scale);
                                }
                            }
                            else if(oa_cut && add_new){
                                nb_added_cuts++;
                                cuts.push_back(-1);
                                cuts.push_back(con_pos);
                                cuts.push_back(i);
                                xres=con->get_x(i);
                                for(auto &x:xres){
                                    cuts.push_back(x);
                                }
                                
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
        con_pos++;
    }
    
    
    set_solution(xsolution);
    nb_oacuts+=nb_added_cuts;
    DebugOff("Number of constraints in linear model = " << nb_oacuts << endl);
    return constr_viol;
}
/*Adds cuts to the current model*/
template<>
void Model<>::add_cuts_to_model(vector<double>& cuts, Model<>& nonlin, int &oacuts){
    int start=0,i, count_con,count_vec, start_orig;
    bool not_found=true, c_nfound, near_zero;
    vector<double> res1, res2, c_val;
    double c0_val, scale;
    const double zero_tol=1e-6;
    string clin="OA_cuts_";
    int l_clin=clin.length();
    while(start!=cuts.size()){
        not_found=true;
        if(cuts[start]==-1){
            start_orig=start;
            start++;
            for(count_con=0;count_con<nonlin._cons_vec.size();count_con++){
                if(count_con==cuts[start]){
                    not_found=false;
                    break;
                }
            }
            if(!not_found){
                start++;
                i=cuts[start++];
                auto con_nonlin=nonlin._cons_vec.at(count_con);
                res1.resize(0);
                res2.resize(0);
                res1=con_nonlin->get_x(i);
                for (auto r=0;r<res1.size();r++){
                    res2.push_back(cuts[start++]);
                }
                con_nonlin->set_x(i, res2);
                con_nonlin->uneval();
                con_nonlin->eval_all();
                con_nonlin->eval(i);
                auto con_lin_name="OA_cuts_"+con_nonlin->_name;
                if(_cons_name.find(con_lin_name)!=_cons_name.end()){
                    count_vec=0;
                    c_nfound=true;
                    for(auto c:_cons_vec){
                        if(c->_name==con_lin_name){
                            c_nfound=false;
                            break;
                        }
                        count_vec++;
                    }
                    if(!c_nfound){
                        c0_val=0;
                        c_val.resize(con_nonlin->_nb_vars,0);
                        near_zero=true;
                        scale=1.0;
                        con_nonlin->get_outer_coef(i, c_val, c0_val);
                        for (auto j = 0; j<c_val.size(); j++) {
                            if(c_val[j]!=0 && std::abs(c_val[j])<zero_tol){
                                if(zero_tol/std::abs(c_val[j])>scale){
                                    scale=zero_tol/std::abs(c_val[j]);
                                }
                            }
                            if(near_zero && c_val[j]!=0 && std::abs(c_val[j])<zero_tol){
                                near_zero=true;
                            }
                            else{
                                near_zero=false;
                            }
                            //coefs.push_back(1e5*c_val[j]);
                        }
                        if(!near_zero){
                            start=start_orig;
                            cuts[start_orig]=count_vec;
                            cuts[start_orig+1]=i;
                            for(auto s=1;s<=c_val.size();s++){
                                cuts[start_orig+s+1]=c_val[s-1]*scale;
                            }
                            cuts[start_orig+c_val.size()+2]=c0_val*scale;
                            
                        }
                        
                    }
                }
                else{
                    Constraint<> OA_cut("OA_cuts_"+con_nonlin->_name);
                    indices activeset("active_"+con_nonlin->_name);
                    activeset.add((*con_nonlin->_indices->_keys)[i]);
                    OA_cut=con_nonlin->get_outer_app(activeset, 1.0);
                    //OA_cut._lazy.push_back(false);
                    if(con_nonlin->_ctype==leq) {
                        add(OA_cut.in(activeset)<=0);
                    }
                    else {
                        add(OA_cut.in(activeset)>=0);
                    }
                    auto coa=get_constraint(OA_cut._name);
                    for(auto o=0;o<coa->get_nb_instances();o++){
                        coa->_violated.push_back(true);
                    }
                    oacuts++;
                    reindex();
                }
                con_nonlin->set_x(i, res1);
                
            }
        }
        else{
            count_con=0;
            for(count_con=0;count_con<_cons_vec.size();count_con++){
                if(count_con==cuts[start]){
                    not_found=false;
                    break;
                }
            }
            if(!not_found){
                auto con_lin=_cons_vec.at(count_con);
                auto con_lin_name=con_lin->_name;
                auto con_nonlin_name=con_lin_name.substr(l_clin);
                auto con_nonlin=nonlin.get_constraint(con_nonlin_name);
                start++;
                i=cuts[start++];
                
                auto nb_inst = con_lin->get_nb_instances();
                con_lin->_indices->add("inst_"+to_string(i)+"_"+to_string(nb_inst));
                con_lin->_dim[0] = con_lin->_indices->size();
                //con_lin->_lazy.push_back(true);
                con_lin->_violated.push_back(true);
                auto count=0;
                for(auto &l: *(con_lin->_lterms)){
                    auto name=l.first;
                    if(!l.second._sign){
                        throw invalid_argument("symbolic negative");
                    }
                    if(l.second._coef->is_param()) {
                        auto p_cst = ((param<>*)(l.second._coef.get()));
                        p_cst->add_val("inst_"+to_string(p_cst->_indices->_keys->size()), cuts[start++]);
                    }
                    else {
                        throw invalid_argument("Coefficient must be parameter");
                    }
                    auto parkeys=l.second._p->_indices->_keys;
                    //                                auto vname = l.second._p->_name.substr(0,l.second._p->_name.find_last_of("."));
                    auto v = con_nonlin->get_var(l.second._p->_name);
                    l.second._p->_indices->add_ref((*parkeys)[v->get_id_inst(i)]);
                    count++;
                }
                //Set value of the constant
                if(con_lin->_cst->is_param()){
                    auto co_cst = ((param<>*)(con_lin->_cst.get()));
                    co_cst->add_val("inst_"+to_string(co_cst->_indices->_keys->size()), cuts[start++]);
                }
                else if(con_lin->_cst->is_function()){
                    auto rhs_f = static_pointer_cast<func<>>(con_lin->_cst);
                    if(!rhs_f->func_is_param()){
                        throw invalid_argument("function should be a param");
                    }
                    auto p = static_pointer_cast<param<>>(rhs_f->_params->begin()->second.first);
                    p->add_val("inst_"+to_string(p->_indices->_keys->size()), cuts[start++]);
                    con_lin->_cst=p;
                    //  rhs_f->eval_all();
                }
                con_lin->eval_all();
                oacuts++;
            }
            if(not_found){
                DebugOn("Constraint not found in solver.cpp 1382"<<endl);
                break;
            }
        }
    }
}
/** function to resolve the lower bound problem and add cuts to it, either until nb_refine times or if no violated point is found as per tolerance active_tol
 @param[in] interior model: Model which can give interior point of current model
 @param[in] obbt_model: Model which is resolved and cuts are added to
 @param[in] lb_solver_type: solver algorithm
 @param[in] nb_refine: max number of model resolves allowed
 @param[in] upper bound: upper_bound of the nonconvex model
 @param[in] lower bound: objective value of obbt_model
 @param[in] ub_scale_value: scaled value of upper bound objective
 @param[in] lb_solver_tol: tolerance with which obbt_model is solved
 @param[in] active_tol: constraint violation tolerance based on which cuts are added
 @param[in] oacuts: incremented by added cuts
 @param[in] abs_tol, rel_tol, zero_tol: parameters to determine if gap closed
 @param[in] abs_tol, rel_tol, zero_tol: parameters to determine if gap closed
 @param[in] abs_tol, rel_tol, zero_tol: parameters to determine if gap closed
 @param[in] max_iter, amx_time: to solve obbt_model
 @return bool close. returns true iff gap is closed
 **/
template<typename type>
template<typename T>
bool Model<type>::root_refine(const Model<type>& interior_model, shared_ptr<Model<type>>& obbt_model, SolverType lb_solver_type, int nb_refine, const double upper_bound, double& lower_bound, const double ub_scale_value, double lb_solver_tol, double active_tol, int& oacuts,  const double abs_tol, const double rel_tol, const double zero_tol, string lin_solver, int max_iter, int max_time){
    int constr_viol=1, lin_count=0, output;
    solver<> LB_solver(obbt_model, lb_solver_type);
    bool close=false;
    while (constr_viol==1 && lin_count<nb_refine){
        LB_solver.run(output = 0, lb_solver_tol, lin_solver, max_iter, max_time);
        if(obbt_model->_status==0){
            lower_bound=obbt_model->get_obj_val()*upper_bound/ub_scale_value;
            DebugOn("Iter linear gap = "<<(upper_bound- lower_bound)/(std::abs(upper_bound))*100<<"%"<<endl);
            DebugOn("lin count "<<lin_count<<endl);
            if (std::abs(upper_bound- lower_bound)<=abs_tol && ((upper_bound- lower_bound))/(std::abs(upper_bound)+zero_tol)<=rel_tol)
            {
                close= true;
                break;
            }
            vector<double> solution(obbt_model->_nb_vars);
            obbt_model->get_solution(solution);
            constr_viol=add_iterative(interior_model, solution, obbt_model, "allvar", oacuts, active_tol);
            DebugOn("oacuts "<<oacuts<<endl);
            obbt_model->reindex();
            obbt_model->reset();
            obbt_model->reset_constrs();
        }
        else{
            lower_bound=numeric_limits<double>::min();
            obbt_model->reindex();
            obbt_model->reset();
            obbt_model->reset_constrs();
            break;
        }
        lin_count++;
        
    }
    return(close);
}
/** function to update variable bounds of current model and vector models for the OBBT algorithm
 @param[in] objective models: Vec with names of each model in models
 @param[in] sol_obj: Vec with objective values of each model in models
 @param[in] sol_obj: Vec with solution status of each model in models
 @param[in] models: Vec of models
 @param[in] fixed point: for each obbt_problem, set to true if fixed point reached
 @param[in] interval_original: interval of each variable given at the start
 @param[in] ub_original: upper bound of each variable given at the start
 @param[in] lb_original: lower bound of each variable given at the start
 @param[in] terminate: iff obbt algorithm reached fixed point then terminate is true
 @param[in] lb_solver_tol: tolerance with which obbt_model is solved
 @param[in] fail: increment by 1 when a obbt subproblem fails
 @param[in] range_tol: vvariable interval width should be greater than equal to range_tol
 @param[in]  fixed_tol_abs, fixed_tol_rel, zero_tol: parameters to determine if fixed point reached for a variable
 @return returns true
 */
template<typename type>
template<typename T>
bool Model<type>::obbt_update_bounds(const std::vector<std::string> objective_models, const std::vector<double>& sol_obj, const std::vector<int>& sol_status, std::vector<shared_ptr<gravity::Model<type>>>& models,    map<string, bool>& fixed_point,  const map<string, double>& interval_original, const map<string, double>& ub_original, const map<string, double>& lb_original, bool& terminate, int& fail, const double range_tol, const double fixed_tol_abs, const double fixed_tol_rel, const double zero_tol){
    std::string msname, mkname,vkname,keyk,dirk, var_key_k;
    double objk, boundk1, temp, tempa, mid, left, right;
    var<> vk;
    for (auto s=0;s<objective_models.size();s++)
    {
        /* Update bounds only if the model status is solved to optimal */
        if(sol_status.at(s)==0)
        {
            msname=objective_models.at(s);
            mkname=msname;
            std::size_t pos = mkname.find("|");
            vkname.assign(mkname, 0, pos);
            mkname=mkname.substr(pos+1);
            pos=mkname.find("|");
            keyk.assign(mkname, 0, pos);
            dirk=mkname.substr(pos+1);
            vk=this->template get_var<T>(vkname);
            var_key_k=vkname+"|"+keyk;
            
            objk=sol_obj.at(s);
            
            auto update_lb=false;
            auto update_ub=false;
            if(dirk=="LB")
            {
                boundk1=vk.get_lb(keyk);
                //Uncertainty in objk=obk+-solver_tolerance, here we choose lowest possible value in uncertainty interval
                objk=std::max(objk-range_tol, boundk1);
            }
            else
            {
                boundk1=vk.get_ub(keyk);
                //Uncertainty in objk=obk+-solver_tolerance, here we choose highest possible value in uncertainty interval
                objk=std::min(objk+range_tol, boundk1);
            }
            if((std::abs(boundk1-objk) <= fixed_tol_abs || std::abs((boundk1-objk)/(boundk1+zero_tol))<=fixed_tol_rel))
            {
                fixed_point[msname]=true;
            }
            else
            {
                if(dirk=="LB"){
                    vk.set_lb(keyk, objk);
                    update_lb=true;
                }
                else{
                    vk.set_ub(keyk, objk);
                    update_ub=true;
                }
                //If crossover in bounds,just exchange them
                if(vk.get_ub(keyk)<vk.get_lb(keyk))
                {
                    fixed_point[var_key_k+"|LB"]=true;
                    fixed_point[var_key_k+"|UB"]=true;
                    temp=vk.get_ub(keyk);
                    tempa=vk.get_lb(keyk);
                    vk.set_ub(keyk, tempa);
                    vk.set_lb(keyk, temp);
                    update_lb=true;
                    update_ub=true;
                }/*If fixed point not reached for any variable, terminate is false*/
                else if(!vk._lift){
                    fixed_point[msname]=false;
                    terminate=false;
                }
            }
            //If interval becomes smaller than range_tol, reset bounds so that interval=range_tol
            if(std::abs(vk.get_ub(keyk)-vk.get_lb(keyk))<range_tol)
            {
                //If original interval is itself smaller than range_tol, do not have to reset interval
                if(interval_original.at(var_key_k)>=range_tol)
                {
                    DebugOff("Entered reset");
                    //Mid is the midpoint of interval
                    mid=(vk.get_ub(keyk)+vk.get_lb(keyk))/2.0;
                    left=mid-range_tol/2.0;
                    right=mid+range_tol/2.0;
                    //If resized interval does not cross original bounds, reset
                    if(right<=ub_original.at(var_key_k) && left>=lb_original.at(var_key_k))
                    {
                        vk.set_ub(keyk, right);
                        vk.set_lb(keyk, left);
                        update_lb=true;
                        update_ub=true;
                    }
                    //If resized interval crosses original upperbound, set the new bound to upperbound, and lower bound is expanded to upperbound-range_tolerance
                    else if(right>ub_original.at(var_key_k))
                    {
                        
                        vk.set_ub(keyk, ub_original.at(var_key_k));
                        vk.set_lb(keyk, ub_original.at(var_key_k)-range_tol);
                        update_lb=true;
                        update_ub=true;
                    }
                    //If resized interval crosses original lowerbound, set the new bound to lowerbound, and upper bound is expanded to lowerbound+range_tolerance
                    else if(left<lb_original.at(var_key_k))
                    {
                        vk.set_lb(keyk, lb_original.at(var_key_k));
                        vk.set_ub(keyk, lb_original.at(var_key_k)+range_tol);
                        update_lb=true;
                        update_ub=true;
                        
                    }
                    //In the resized interval both original lower and upper bounds can not be crosses, because original interval is greater
                    //than range_tol
                    
                }
            }
            if(update_lb||update_ub){
                auto vkptr=get_var_ptr(vkname);
                vkptr->_new=true;
                for(auto &mod:models){
                    auto vkmod=mod->template get_var<T>(vkname);
                    //                                                            auto vmodptr=mod->get_var_ptr(vkname);
                    //                                                            vmodptr->_new=true;
                    if(update_lb){
                        vkmod.set_lb(keyk, vk.get_lb(keyk));
                    }
                    if(update_ub){
                        vkmod.set_ub(keyk, vk.get_ub(keyk));
                    }
                }
            }
        }
        else
        {
            fail++;
        }
    }
    return 0;
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
template<>
int Model<>::cuts_parallel(vector<shared_ptr<Model<>>> batch_models, int batch_model_count, const Model<>& interior_model, shared_ptr<Model<>> lin, int& oacuts, double active_tol, int run_obbt_iter, double range_tol, string vname){
    std::vector<double> obbt_solution, cut_vec;
    string msname;
    int added_cuts=0;
    int viol=0, viol_i=0;
    for (auto s=0;s<batch_model_count;s++)
    {
        auto m=batch_models[s];
        /* Update bounds only if the model status is solved to optimal */
        if(m->_status==0)
        {
            if(vname=="modelname"){
                msname=m->_name;
            }
            else if(vname=="allvar"){
                msname="allvar";
            }
            obbt_solution.resize(m->_nb_vars);
            m->get_solution(obbt_solution);
            cut_vec.resize(0);
            viol_i=generate_cuts_iterative(interior_model, obbt_solution, lin, msname, oacuts, active_tol, cut_vec);
            lin->add_cuts_to_model(cut_vec, *this, added_cuts);
            if(viol==1 || viol_i==1){
                viol=1;
            }
        }
    }
    return viol;
}
template<>
int Model<>::cuts_parallel(vector<shared_ptr<Model<>>> batch_models, int batch_model_count, const Model<>& interior_model, shared_ptr<Model<>> lin, int& oacuts, double active_tol, int run_obbt_iter, double range_tol, string vname, std::vector<std::string>& repeat_list){
    std::vector<double> obbt_solution, cut_vec;
    string msname;
    int added_cuts=0;
    int viol=0, viol_i=0;
    for (auto s=0;s<batch_model_count;s++)
    {
        auto m=batch_models[s];
        /* Update bounds only if the model status is solved to optimal */
        if(m->_status==0)
        {
            if(vname=="modelname"){
                msname=m->_name;
            }
            else if(vname=="allvar"){
                msname="allvar";
            }
            obbt_solution.resize(m->_nb_vars);
            m->get_solution(obbt_solution);
            cut_vec.resize(0);
            viol_i=generate_cuts_iterative(interior_model, obbt_solution, lin, msname, oacuts, active_tol, cut_vec);
            m->set_solution(obbt_solution);
            m->add_cuts_to_model(cut_vec, *this, added_cuts);
            if(viol_i==1){
                repeat_list.push_back(m->_name);
            }
            if(viol==1 || viol_i==1){
                viol=1;
            }
        }
    }
    return viol;
}

#ifdef USE_MPI

template<>
int Model<>::cuts_MPI(vector<shared_ptr<Model<>>>& batch_models, int batch_model_count, const Model<>& interior_model, shared_ptr<Model<>> lin, int& oacuts, double active_tol, int run_obbt_iter, double range_tol, vector<int>& sol_status, string vname){
    int worker_id, nb_workers;
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
    auto err_size = MPI_Comm_size(MPI_COMM_WORLD, &nb_workers);
    auto nb_workers_ = std::min((size_t)nb_workers, (size_t)batch_model_count);
    std::string msname;
    vector<double> obbt_solution, cut_vec;
    vector<int> bcast_array;
    bcast_array.resize(3,0);
    int added_cuts=0;
    int cut_size=0;
    int viol=0, viol_i=0;
    if(batch_model_count!=0){
        std::vector<size_t> limits = bounds(nb_workers_, batch_model_count);
        for (auto w_id = 0; w_id<nb_workers; w_id++) {
            if(w_id+1<limits.size()){
                for (auto i = limits[w_id]; i < limits[w_id+1]; i++) {
                    if(sol_status[i]==0){
                        if(worker_id==w_id){
                            if(vname=="modelname"){
                                msname=batch_models[i-limits[w_id]]->_name;
                            }
                            else if(vname=="allvar"){
                                msname="allvar";
                            }
                            obbt_solution.resize(batch_models[i-limits[w_id]]->_nb_vars);
                            batch_models[i-limits[w_id]]->get_solution(obbt_solution);
                            cut_vec.resize(0);
                            viol_i=generate_cuts_iterative(interior_model, obbt_solution, lin, msname, oacuts, active_tol, cut_vec);
                            if((viol==1) || (viol_i==1)){
                                viol=1;
                            }
                            
                            cut_size=cut_vec.size();
                            bcast_array[0]=viol;
                            bcast_array[1]=cut_size;
                            bcast_array[2]=oacuts;
                        }
                        MPI_Bcast(&bcast_array[0], 3, MPI_INT, w_id, MPI_COMM_WORLD);
                        MPI_Barrier(MPI_COMM_WORLD);
                        viol=bcast_array[0];
                        cut_size=bcast_array[1];
                        oacuts=bcast_array[2];
                        if(cut_size!=0){
                            if(worker_id!=w_id){
                                cut_vec.resize(cut_size, 0);
                            }
                            MPI_Bcast(&cut_vec[0], cut_size, MPI_DOUBLE, w_id, MPI_COMM_WORLD);
                            lin->add_cuts_to_model(cut_vec, *this, added_cuts);
                            MPI_Barrier(MPI_COMM_WORLD);
                        }
                    }
                }
            }
        }
    }
    return viol;
}

template<>
int Model<>::cuts_MPI(vector<shared_ptr<Model<>>>& batch_models, int batch_model_count, const Model<>& interior_model, shared_ptr<Model<>> lin, int& oacuts, double active_tol, int run_obbt_iter, double range_tol, vector<int>& sol_status, string vname, std::vector<std::string>& repeat_list, const std::vector<size_t>& limits){
    int worker_id, nb_workers;
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
    auto err_size = MPI_Comm_size(MPI_COMM_WORLD, &nb_workers);
    auto nb_workers_ = std::min((size_t)nb_workers, (size_t)batch_model_count);
    std::string msname;
    vector<double> obbt_solution, cut_vec, cut_array_w, cut_array_allw;
    int cut_size=0;
    int viol=0, viol_i, c_old;
    int oacuts_old=oacuts;
    double t_start=0, t=0, tb_time=0, t1, t2;
    vector<int> cut_size_vec, d;
    cut_size_vec.resize(nb_workers, 0);
    //t_start=get_wall_time();
    if(batch_model_count!=0){
        if(worker_id+1<limits.size()){
            for (auto i = limits[worker_id]; i < limits[worker_id+1]; i++) {
                if(sol_status[i]==0){
                    if(vname=="modelname"){
                        msname=batch_models[i-limits[worker_id]]->_name;
                    }
                    else if(vname=="allvar"){
                        msname="allvar";
                    }
                    obbt_solution.resize(batch_models[i-limits[worker_id]]->_nb_vars);
                    batch_models[i-limits[worker_id]]->get_solution(obbt_solution);
                    cut_vec.clear();
                    viol_i=generate_cuts_iterative(interior_model, obbt_solution, lin, msname, oacuts, active_tol, cut_vec);
                    batch_models[i-limits[worker_id]]->set_solution(obbt_solution);
                    for(auto &c:cut_vec){
                        cut_array_w.push_back(c);
                    }
                }
            }
            cut_size=cut_array_w.size();
        }
        //t1=get_wall_time();
        MPI_Allgather(&cut_size, 1, MPI_INT, &cut_size_vec[0], 1, MPI_INT, MPI_COMM_WORLD);
        //t2=get_wall_time();
        //tb_time+=t2-t1;
        int cut_all_size=0;
        c_old=0;
        d.push_back(0);
        for(auto &c:cut_size_vec){
            if((viol==0) && c>0){
                viol=1;
            }
            cut_all_size+=c;
            d.push_back(c_old+c);
            c_old+=c;
        }
        d.pop_back();
        cut_array_allw.resize(cut_all_size);
        //t1=get_wall_time();
        MPI_Allgatherv(&cut_array_w[0], cut_size_vec[worker_id], MPI_DOUBLE,
                       &cut_array_allw[0], &cut_size_vec[0], &d[0], MPI_DOUBLE, MPI_COMM_WORLD);
        //t2=get_wall_time();
        //tb_time+=t2-t1;
        oacuts=oacuts_old;
        lin->add_cuts_to_model(cut_array_allw, *this, oacuts);
    }
    //t=get_wall_time();
    //DebugOff(endl<<endl<<"wid "<<worker_id<<" cuts_MPI "<<(t-t_start)<<" Broad "<<tb_time<<endl<<endl);
    MPI_Barrier(MPI_COMM_WORLD);
    return viol;
}




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

int run_MPI_new(const std::vector<std::string> objective_models, std::vector<double>& sol_obj, std::vector<int>& sol_status, std::vector<shared_ptr<gravity::Model<double>>>& models, const shared_ptr<gravity::Model<double>>& relaxed_model, const gravity::Model<double>& interior, string cut_type, double active_tol, gravity::SolverType stype, double tol, unsigned nr_threads, const string& lin_solver, int max_iter, int max_batch_time, bool linearize, int nb_refine, std::vector<size_t> limits){
    
    
    int worker_id, nb_workers;
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
    auto err_size = MPI_Comm_size(MPI_COMM_WORLD, &nb_workers);
    auto nb_workers_ = std::min((size_t)nb_workers, objective_models.size());
    std::vector<std::string> objective_models_worker;
    std::vector<double> sol_obj_worker;
    std::vector<int> sol_status_worker;
    if(nb_workers_!=limits.size()-1){
        DebugOn("Error4 in computing limits");
    }
    std::string msname, mname, vname, key, dir;
    var<> var;
    if(objective_models.size()!=0){
        /* Split models into equal loads */
        //std::vector<size_t> limits = bounds(nb_workers_, objective_models.size());
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
                objective_models_worker.push_back(objective_models[i]);
            }
            sol_status_worker.resize(limits[worker_id+1]-limits[worker_id],0);
            sol_obj_worker.resize(limits[worker_id+1]-limits[worker_id],0);
            run_parallel_new(objective_models_worker, sol_obj_worker, sol_status_worker, models, relaxed_model, interior, cut_type, active_tol, stype, tol, nr_threads, lin_solver, max_iter, max_batch_time, linearize, nb_refine);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        send_status_new(models,limits, sol_status);
        MPI_Barrier(MPI_COMM_WORLD);
        send_obj_all_new(models,limits, sol_obj);
    }
    //   MPI_Barrier(MPI_COMM_WORLD);
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
template bool Model<double>::add_iterative(const Model<double>& interior, vector<double>& obbt_solution, shared_ptr<Model<double>>& lin, string modelname, int& nb_oacuts, double active_tol);
template bool Model<double>::root_refine(const Model<double>& interior_model, shared_ptr<Model<double>>& obbt_model, SolverType lb_solver_type, int nb_refine, const double upper_bound, double& lower_bound, const double ub_scale_value, double lb_solver_tol, double active_tol, int& oacuts, const double abs_tol, const double rel_tol, const double zero_tol, string lin_solver, int max_iter, int max_time);
template bool Model<double>::obbt_update_bounds(const std::vector<std::string> objective_models, const std::vector<double>& sol_obj, const std::vector<int>& sol_status, std::vector<shared_ptr<gravity::Model<double>>>& models, map<string, bool>& fixed_point,  const map<string, double>& interval_original, const map<string, double>& ub_original, const map<string, double>& lb_original, bool& terminate, int& fail, const double range_tol, const double fixed_tol_abs, const double fixed_tol_rel, const double zero_tol);
}
