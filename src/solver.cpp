    //
    //  solver.cpp
    //  Gravity++
    //
    //  Created by Hassan on 30/01/2015.

    //

#include <gravity/solver.h>
#include <mutex>

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
//        if(con->is_eq())
            
        if(!con->is_linear()) {
            if(con->_ctype==eq)
            {
                continue;
//                throw invalid_argument("Equality constraint is not currently supported");
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
/** /** Returns an interior point model
 @param[in] nonlin: model for which interior point with respect to nonconvex constraints that describe a convex region(SDP,SOC,rotated SOC constraints) is computed
 Assuming model has no nonlinear equality constraints
 Returns a model which has all variables in current model. If current model has SDP, SOC, rotated SOC constraints g(x), a new model is created with g_i(x) \le eta_i and the objective is to min (sum eta_i)
 **/
template<typename type>
template<typename T>
Model<type> Model<type>::build_model_interior() const
{
    Model<type> Interior(_name+"Interior");
    
    /* Variables of current model are added*/
    for (auto &it: _vars)
    {
        auto v = it.second;
        Interior.add_var(v);
    }
    /* Index set for eta variables*/
    indices ind_eta("ind_eta");
    /* Vector of indices, each element corresponds to a nonconvex constraint, and has indices of eta corresponding to all its instances */
    vector<indices> ind_eta_vec;
    for (auto &con: _cons_vec)
    {
        if(!con->is_linear()) {
            if(con->_callback){
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
    }
    /*Add eta variables to model*/
    var<> eta_int("eta_interior", -1, 0);
    Interior.add(eta_int.in(ind_eta));
    
    /* Objective */
    Interior.min(sum(eta_int));
    
    /* Constraints */
    int count=0;
    for (auto &con: _cons_vec)
    {
        if(!con->is_linear()) {
            if(con->_callback){
            /* We are only interested in an interior point for constraints defining a convex region but having a non-convex description, e.g., SDP-determinant cuts and SOC constraints.*/
            if(!con->is_convex() || con->is_rotated_soc() || con->check_soc()){
                /*ind has indices of all eta_int elements corressponding to con*/
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
//        else{
//            Constraint<> Inter_con(*con);
//
//            if(con->_ctype==leq)
//            {
//                Interior.add(Inter_con<=0);
//            }
//            else  if(con->_ctype==geq)
//            {
//                Interior.add(Inter_con>=0);
//            }
//
//        }
    }
    return *Interior.copy();
}

/* Runs models stored in the vector models in parallel, using the solver object solvers. Objective of models to run is given in objective_models. Solution status and objective of each model is added to sol_status, and sol_obj. Liner constraints are added to each model calling add_iterative.  relaxed_model, interior, cut_type,  nb_oacuts, active_tol are args to add_iterative*/
int run_parallel_new(const std::vector<std::string> objective_models, std::vector<double>& sol_obj, std::vector<int>& sol_status, std::vector<shared_ptr<gravity::Model<double>>>& models, const shared_ptr<gravity::Model<double>>& relaxed_model, const gravity::Model<double>& interior, string cut_type, double active_tol, gravity::SolverType stype, double tol, unsigned nr_threads, const string& lin_solver, int max_iter, int max_batch_time, bool linearize, int nb_refine, vector<vector<double>>& vbasis, vector<std::map<string,double>>& cbasis, bool initialize_primal){
    std::vector<thread> threads;
    std::vector<double> solution(models[0]->_nb_vars);
    std::string mname, msname,vname, key, dir, modelname;
    var<> var;
    vector<shared_ptr<solver<double>>> batch_solvers;
    int viol=0, count=0, ncuts=0, nref;
    vector<int> viol_i;
    if(objective_models.size()==0){
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
        models[s]->reset_lifted_vars_bounds();
        models[s]->reset_constrs();
        models[s]->reset();
        models[s]->_status=0;
            //models[s]->print();
        DebugOff("to create solver"<<endl);
        auto solverk = make_shared<solver<double>>(models[s], stype);
        if(stype==gurobi && initialize_primal){
            solverk->initialize_basis(vbasis.at(s), cbasis.at(s));
        }
        batch_solvers.push_back(solverk);
        
        DebugOff("created solver"<<endl);
    }
    /* Split models into nr_threads parts */
    auto nr_threads_ = std::min((size_t)nr_threads,objective_models.size());
    std::vector<size_t> limits = bounds(nr_threads_, objective_models.size());
    viol_i.resize(nr_threads_, 1);
    /* Launch all threads in parallel */
    for(auto l=0;l<nb_refine;l++){
        viol=0;
        count=0;
        if(l>=1){
            for (auto s=0;s<objective_models.size();s++){
                if(stype==ipopt && batch_solvers[s]->_model->_objt==maximize){
                    *batch_solvers[s]->_model->_obj *= -1;
                }
                if(stype==gurobi){
                    batch_solvers[s]->_model->_obj->_new= false;
                }
            }
        }
            //auto vec = vector<shared_ptr<gravity::Model<double>>>(models);
        for (size_t i = 0; i < nr_threads_; ++i) {
            if(viol_i.at(i)==1 && models[i]->_status==0){
                if(l>=1)
                    DebugOff("resolving"<<endl);
                threads.push_back(thread(run_models_solver<double>, ref(models), ref(batch_solvers), limits[i], limits[i+1], stype, tol, lin_solver, max_iter, max_batch_time));
            }
        }
        /* Join the threads with the main thread */
        for(auto &t : threads){
            t.join();
        }
        if(linearize){
            for(auto &m:models){
                if(count<objective_models.size()){
                    if(m->_status==0 && viol_i.at(count)==1){
                        m->get_solution(solution);
                        DebugOff("active tol "<<active_tol<<endl);
                        viol_i.at(count)=relaxed_model->add_iterative(interior,  solution, models[count],  "allvar",  ncuts,  active_tol);
                    }
                    if(viol_i.at(count)==1){
                        viol=1;
                    }
                    count++;
                }
            }
        }
        threads.clear();
        nref=l+1;
        if(viol==0){
            break;
        }
    }
    if(nref>1)
        DebugOn("nrefine "<<nref<<" ncuts "<<ncuts<<endl<<endl);
    count=0;
    sol_status.resize(objective_models.size(),-1);
    sol_obj.resize(objective_models.size(),-1.0);
    for(auto &m:models){
        if(count<objective_models.size()){
            sol_status.at(count)=m->_status;
            sol_obj.at(count)=m->get_obj_val();
            count++;
        }
    }
    count=0;
    if(stype==gurobi && initialize_primal){
        for(auto&s: batch_solvers){
            if(sol_status.at(count)==0){
                s->get_pstart(vbasis.at(count), cbasis.at(count));
            }
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
int run_parallel(const vector<shared_ptr<gravity::Model<double>>>& models, gravity::SolverType stype, double tol, unsigned nr_threads, int max_iter){
    return run_parallel(models,stype,tol,nr_threads,"",max_iter);
};

/*Returns a model that can compute interior point of nonlin model. Generates OA cuts to nonlin model at solution of nonlin model and adds it to current model.
 @param[in] nonlin: Nonlinear model whose nonlinear constraints are linearzied and added to current model
 @return interior model: Model which can give interior point of current model
 */
template<>
Model<> Model<>::add_outer_app_solution(Model<>& nonlin)
{
    const double active_tol_sol=1e-12, zero_tol=1e-6;
    vector<double> xactive, xcurrent, xinterior, xres, xtest;
    bool interior=false, convex_region=true;
    size_t posv;
    int output=0;
    double scale=1.0, tol=1e-8;
    auto Ointerior = nonlin.build_model_interior();
    //Ointerior.print();
    solver<> modelI(Ointerior, ipopt);
    modelI.run(output=0, tol);
    vector<double> xsolution(_nb_vars);
    nonlin.get_solution(xsolution);
    if((Ointerior._status==0||Ointerior._status==1) && Ointerior.get_obj_val()<0)
    {
        interior=true;
    }
    else{
        throw invalid_argument("Failed to solve interior model");
    }
    for (auto &con: nonlin._cons_vec)
    {
        if(con->is_eq())
            continue;
        if(!con->is_linear() && !con->has_int()){
            if(!con->is_convex() || con->is_rotated_soc() || con->check_soc())
            {
                Constraint<> OA_sol("OA_cuts_"+con->_name);
                indices activeset("active_"+con->_name);
                auto keys=con->_indices->_keys;
                for(auto i=0;i<con->get_nb_inst();i++){
                    convex_region=con->check_convex_region(i);
                    if(convex_region){
                        con->uneval();
                        con->eval_all();
                        if(con->is_active(i,active_tol_sol)){
                            activeset.add((*keys)[i]);
                        }
                    }
                }
                
                    //  OA_sol=con->get_outer_app(activeset, scale);
                OA_sol=con->get_outer_app(activeset, scale, true, 1e-3);
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
                    /*Avoiding linearizations too close to zero*/
                    posv=x->get_id_inst(i);
                    x->get_double_val(posv, xval);
                    if(std::abs(xval)>=1e-3){
                        string keyi=(*keys)[i];
                        allset.add(keyi);
                    }
                }
                    // OA_sol=con->get_outer_app(allset, scale);
                OA_sol=con->get_outer_app(allset, scale, true, 1e-3);
                if(con->_ctype==leq) {
                    add(OA_sol.in(allset)<=0);
                }
                else {
                    add(OA_sol.in(allset)>=0);
                }
            }
            else if(con->is_convex() && !con->is_rotated_soc() && !con->check_soc())
            {
               // if(con->_name!="limit_neg"){
                Constraint<> OA_sol("OA_cuts_"+con->_name);
                indices allset("active_"+con->_name);
                auto keys=con->_indices->_keys;
                for(auto i=0;i<con->get_nb_inst();i++){
                    string keyi=(*keys)[i];
                    allset.add(keyi);
                }
                    //OA_sol=con->get_outer_app(allset, scale);
                OA_sol=con->get_outer_app(allset, scale, true, 1e-3);
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
    reset();
    reset_constrs();
    set_solution(xsolution);
    return Ointerior;
}
/*Generates and adds constraints to model lin. The curent model solution is set to obbt_solution and OA cuts are generated for the nonlinear constraints in the current model at the obbt_solution point. These cuts are added to model lin.
 @param[in] interior model: Model which can give interior point of current model
 @param[in] obbt_solution: Point at which constraints are linearized. For non-convex constraints that define a convex region, the point is shifted to an active point of that constraint and its instance
 @param[in] lin: Model to which linear cuts are added
 @param[in] nb_oacuts: When a cut is added nb_oacuts is incremented
 @param[in] active_tol: the obbt_solution x is said to violate a nonlinear constraint g in current model if abs(g(x))> active_tol */
template<typename type>
template<typename T>
bool Model<type>::add_iterative(const Model<type>& interior, vector<double>& obbt_solution, shared_ptr<Model<type>>& lin, string modelname, int& nb_oacuts, double active_tol)
{
    vector<double> xsolution(_nb_vars);
    vector<double> xinterior(_nb_vars);
    vector<double> xcurrent, c_val;
    get_solution(xsolution);
    set_solution(obbt_solution);
    const double active_tol_sol=1e-12, zero_tol=1e-6;
    double c0_val, scale=1.0;
    bool constr_viol=false, oa_cut=true, convex_region=true, add_new=false;
    int nb_added_cuts = 0;
    string cname,con_lin_name;
    for (auto &con: _cons_vec)
    {
        if(!con->is_linear()) {
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
                {
                    auto fk=con->eval(i);
                    if((fk > active_tol && con->_ctype==leq) || (fk < -active_tol && con->_ctype==geq)){
                        constr_viol=true;
                            //ToDo fix interior status and check for it
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
                    oa_cut=false;
                    convex_region=con->check_convex_region(i);
                    if(convex_region){
                        scale=1.0;
                        if(add_new){
                            nb_added_cuts++;
                            indices activeset("active_"+con->_name);
                            activeset.add((*con->_indices->_keys)[i]);
                            Constraint<> OA_cut(con_lin_name);
                                //OA_cut=con->get_outer_app(activeset, scale);
                            OA_cut=con->get_outer_app(activeset, scale, true, 1e-3);
                            if(con->_ctype==leq){
                                lin->add(OA_cut.in(activeset)<=0);
                            }
                            else{
                                lin->add(OA_cut.in(activeset)>=0);
                            }
                            add_new=false;
                            oa_cut=false;
                        }
                        else{
                            con->get_outer_coef(i, c_val, c0_val);
                            get_row_scaling(c_val, scale, oa_cut, zero_tol, 1e-3, 1000);
                        }
                    }
                }
                if(oa_cut){
                    nb_added_cuts++;
                    lin->add_linear_row(*con,i, c_val, c0_val, scale);
                }
                con->set_x(i, xcurrent);
                xcurrent.clear();
                xinterior.clear();
            }
        }
    }
    set_solution(xsolution);
    lin->set_solution(obbt_solution);
    lin->reindex();
    lin->reset();
    lin->reset_constrs();
    nb_oacuts+=nb_added_cuts;
    DebugOff("Number of constraints in linear model = " << nb_oacuts << endl);
    return(constr_viol);
}

/*Generates and creates a vector of cuts when any solution is violated by the model this . The curent model solution is set to obbt_solution and OA cuts are generated for the nonlinear constraints in the current model at the obbt_solution point. These cuts are added to model lin.
 @param[in] interior model: Model which can give interior point of current model
 @param[in] obbt_solution: Point at which constraints are linearized. For non-convex constraints that define a convex region, the point is shifted to an active point of that constraint and its instance
 @param[in] lin: Model to which linear cuts are added
 @param[in] nb_oacuts: When a cut is added nb_oacuts is incremented
 @param[in] active_tol: the obbt_solution x is said to violate a nonlinear constraint g in current model if abs(g(x))> active_tol */
template<typename type>
template<typename T>
vector<vector<double>> Model<type>::cutting_planes_solution(const Model<type>& interior, double active_tol)
{
    vector<double> xsolution(_nb_vars);
    vector<double> xinterior(_nb_vars);
    vector<double> xcurrent, c_val;
    vector<vector<double>> res;
    vector<double> cut;
    const double active_tol_sol=1e-12, zero_tol=1e-6;
    double c0_val, scale=1.0, fk;
    bool constr_viol=false, oa_cut=true, convex_region=true, add_new=false;
    int nb_added_cuts = 0;
    for (auto &con: _cons_vec)
    {
        if(!con->is_linear() && con->_callback) {
           // if(con->_name!="limit_neg"){
            auto cnb_inst=con->get_nb_inst();
            for(auto i=0;i<cnb_inst;i++){
                oa_cut=false;
                c0_val=0;
                c_val.resize(con->_nb_vars,0);
                auto cname=con->_name;
                xcurrent=con->get_x(i);
              
                con->uneval();
                con->eval_all();
                DebugOff(con->_name<<"\t"<<con->eval(i)<<endl);
                if(con->is_active(i,active_tol_sol)){
                    oa_cut=false;
                }
                else
                {
                    con->uneval();
                    fk=con->eval(i);
                    if((fk >= active_tol && con->_ctype==leq) || (fk <= -active_tol && con->_ctype==geq)){
                        constr_viol=true;
                            //ToDo fix interior status and check for it
                        if((!con->is_convex()||con->is_rotated_soc() || con->check_soc()))  {
                            auto con_interior=interior.get_constraint(cname);
                            xinterior=con_interior->get_x_ignore(i, "eta_interior"); /** ignore the Eta (slack) variable */
                            auto res_search=con->binary_line_search(xinterior, i);
                            if(res_search){
                                oa_cut=true;
                                    //oa_cut=false;
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
                    oa_cut=false;
                    convex_region=true;
                    if(cname.find("soc")!=std::string::npos){
                        convex_region=con->check_convex_region_soc_rotation(i);
                    }
                    else if(cname.find("det")!=std::string::npos){
                        convex_region=con->check_convex_region_det_rotation(i);
                    }
                    if(convex_region){
                        scale=1.0;
                        con->get_outer_coef(i, c_val, c0_val);
                        get_row_scaling(c_val, scale, oa_cut, zero_tol, 1e-3, 1000);
                    }
                }
                if(oa_cut){
                    auto j=0;
                    for (auto &v_p: con->get_vars()){
                        auto vid=v_p.second.first->get_id() + v_p.second.first->get_id_inst(i);
                        cut.push_back(vid);
                        cut.push_back(c_val[j++]*scale);
                    }
                    cut.push_back(c0_val*scale);
                    res.push_back(cut);
                    if(con->is_convex()){
                        DebugOff(con->_name<<" "<<fk );
                    }
                }
                con->set_x(i, xcurrent);
                xcurrent.clear();
                xinterior.clear();
                cut.clear();
            }
            
        }
    }
    
    return res;
    
}


/*Adds row(or new instance) of a linear constraint to a model by linearizing a nonlinear constraint con
 @param[in] con: Nonlinear constraint
 @param[in] c_inst: Instance of nonlinear constraint which is to be linearized
 @param[in] c_val: Coefficients of the new row of the linear constraint
 @param[in] c0_val: Constant term of the new row of the linear constraint
 @param[in] scaling factor of the new linear constraint row */
template<typename type>
template<typename T>
void Model<type>::add_linear_row(Constraint<type>& con, int c_inst, const vector<double>& c_val, const double c0_val, const double scale){
    auto con_lin=this->get_constraint("OA_cuts_"+con._name);
    auto nb_inst = con_lin->get_nb_instances();
    con_lin->_indices->add("inst_"+to_string(c_inst)+"_"+to_string(nb_inst));
    con_lin->_dim[0] = con_lin->_indices->_keys->size();
    con_lin->_violated.push_back(true);
    auto count=0;
    DebugOff("nb inst "<<nb_inst);
    for(auto &l: *(con_lin->_lterms)){
        auto name=l.first;
        if(!l.second._sign){
            throw invalid_argument("symbolic negative");
        }
        if(l.second._coef->is_param()) {
            auto p_cst = ((param<>*)(l.second._coef.get()));
            p_cst->add_val("inst_"+to_string(c_inst)+"_"+to_string(nb_inst), c_val[count]*scale);
            DebugOff("added p"<<endl);
        }
        else {
            auto f = static_pointer_cast<func<>>(l.second._coef);
            if(!f->func_is_param()){
                throw invalid_argument("function should be a param");
            }
            auto p = static_pointer_cast<param<>>(f->_params->begin()->second.first);
            p->add_val("inst_"+to_string(c_inst)+"_"+to_string(nb_inst), c_val[count]*scale);
            l.second._coef = p;
        }
        auto parkeys=l.second._p->_indices->_keys;
        auto v = con.get_var(l.second._p->_name);
        l.second._p->_indices->add_ref((*parkeys)[v->get_id_inst(c_inst)]);
        count++;
    }
        //Set value of the constant
    if(con_lin->_cst->is_param()){
        auto co_cst = ((param<>*)(con_lin->_cst.get()));
        co_cst->add_val("inst_"+to_string(c_inst)+"_"+to_string(nb_inst), c0_val*scale);
    }
    else if(con_lin->_cst->is_function()){
        auto rhs_f = static_pointer_cast<func<>>(con_lin->_cst);
        if(!rhs_f->func_is_param()){
            throw invalid_argument("function should be a param");
        }
        auto rhs_p = static_pointer_cast<param<>>(rhs_f->_params->begin()->second.first);
        rhs_p->add_val("inst_"+to_string(c_inst)+"_"+to_string(nb_inst), c0_val*scale);
        con_lin->_cst = rhs_p;
    }
    con_lin->uneval();
    con_lin->eval_all();
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
bool Model<type>::root_refine(const Model<type>& interior_model, shared_ptr<Model<type>>& obbt_model, SolverType lb_solver_type, int nb_refine, const double upper_bound, double& lower_bound, const double lb_scale_value, double lb_solver_tol, double active_tol, int& oacuts,  const double abs_tol, const double rel_tol, const double zero_tol, string lin_solver, int max_iter, int max_time, std::vector<double>& vrbasis, std::map<string,double>& crbasis, bool initialize_primal){
    int constr_viol=1, lin_count=0, output;
    solver<> LB_solver(obbt_model, lb_solver_type);
    if(initialize_primal && lb_solver_type==gurobi){
        LB_solver.initialize_basis(vrbasis, crbasis);
    }
    vector<double> solution(obbt_model->_nb_vars);
    bool close=false;
#ifdef USE_MPI
    int worker_id, nb_workers;
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
    auto err_size = MPI_Comm_size(MPI_COMM_WORLD, &nb_workers);
#endif
    while (constr_viol==1 && lin_count<nb_refine){
        LB_solver.run(output = 0, lb_solver_tol, max_iter, max_time);
        if(obbt_model->_status==0){
            lower_bound=obbt_model->get_obj_val()*lb_scale_value;
                // gap=(upper_bound- lower_bound)/(std::abs(upper_bound))*100;
#ifdef USE_MPI
            if(worker_id==0)
#endif
                DebugOn("Iter linear gap = "<<(upper_bound- lower_bound)/(std::abs(upper_bound)+zero_tol)*100<<"% "<<lin_count<<endl);
            if (std::abs(upper_bound- lower_bound)<=abs_tol && ((upper_bound- lower_bound))/(std::abs(upper_bound)+zero_tol)<=rel_tol)
            {
                close= true;
                break;
            }
            /*    if(gap_old-gap<=0.01){
             #ifdef USE_MPI
             if(worker_id==0)
             #endif
             DebugOn("gap less than 0.01"<<endl);
             break;
             }
             gap_old=gap;*/
            obbt_model->get_solution(solution);
            constr_viol=add_iterative(interior_model, solution, obbt_model, "allvar", oacuts, active_tol);
        }
        else{
#ifdef USE_MPI
            if(worker_id==0)
            {
#endif
                DebugOn("lower bounding failed in root refine"<<lin_count<<endl);
                LB_solver.run(output = 5, lb_solver_tol, max_iter, max_time);
                obbt_model->print();
#ifdef USE_MPI
            }
#endif
            lower_bound=numeric_limits<double>::min();
            break;
        }
        lin_count++;
    }
    if(initialize_primal && lb_solver_type==gurobi && obbt_model->_status==0){
        LB_solver.get_pstart(vrbasis,crbasis);
    }
    obbt_model->reindex();
    obbt_model->reset();
    obbt_model->reset_constrs();
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
bool Model<type>::obbt_batch_update_bounds(const std::vector<std::string> objective_models, const std::vector<double>& sol_obj, const std::vector<int>& sol_status, std::vector<shared_ptr<gravity::Model<type>>>& models, shared_ptr<gravity::Model<type>>& obbt_model,   map<string, bool>& fixed_point,  const map<string, double>& interval_original, const map<string, double>& ub_original, const map<string, double>& lb_original, bool& terminate, int& fail, const double range_tol, const double fixed_tol_abs, const double fixed_tol_rel, const double zero_tol, int run_obbt_iter){
#ifdef USE_MPI
    int worker_id, nb_workers;
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
    auto err_size = MPI_Comm_size(MPI_COMM_WORLD, &nb_workers);
#endif
    std::string msname, mkname,vkname,keyk,dirk, var_key_k;
    double objk;
    bool bound_converge=true;
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
            objk=sol_obj.at(s);
            
            this->obbt_update_bounds(bound_converge,objk, msname,vkname, keyk, dirk,  models,  obbt_model,   fixed_point,  interval_original,  ub_original,  lb_original, terminate, fail, range_tol, fixed_tol_abs, fixed_tol_rel,  zero_tol, run_obbt_iter);
            
            
            
#ifdef USE_MPI
            if(worker_id==0)
#endif
                DebugOff("success "<<objective_models.at(s)<<endl);
        }
        else
        {
#ifdef USE_MPI
            if(worker_id==0)
#endif
                DebugOn("failed "<<objective_models.at(s)<<endl);
            fail++;
        }
    }
    
    return 0;
}

template<typename type>
template<typename T>
bool Model<type>::obbt_update_bounds(bool bound_converge,double objk, std::string msname,std::string vkname, std::string keyk, std::string dirk, std::vector<shared_ptr<gravity::Model<type>>>& models, shared_ptr<gravity::Model<type>>& obbt_model,   std::map<string, bool>& fixed_point, const map<string, double>& interval_original, const map<string, double>& ub_original, const map<string, double>& lb_original, bool& terminate, int& fail, const double range_tol, const double fixed_tol_abs, const double fixed_tol_rel, const double zero_tol, int run_obbt_iter)
{
    double boundk1, temp, tempa, mid, left, right, interval, vi_ub_val;
    bool in_orig_model;
    var<> vk, var_ub;
    auto var_key_k=vkname+"|"+keyk;
    auto update_lb=false;
    auto update_ub=false;
    vk=obbt_model->template get_var<T>(vkname);
        // if interval is less than range_tol, fixed point is reached, bounds are updated below Update bounds
    if(dirk=="LB"){
        boundk1=vk.get_lb(keyk);
        interval=vk.get_ub(keyk)-objk;
    }
    else{
        boundk1=vk.get_ub(keyk);
        interval=objk-vk.get_lb(keyk);
    }
    if(std::abs(interval)<=range_tol)
    {
        fixed_point[var_key_k+"|LB"]=true;
        fixed_point[var_key_k+"|UB"]=true;
    }
    
        // if new bound converges to previous bound, fixed point is reached and no need to update bounds
    if(!vk._is_relaxed && bound_converge && (std::abs(boundk1-objk) <= fixed_tol_abs || std::abs((boundk1-objk)/(boundk1+zero_tol))<=fixed_tol_rel))
    {
        if(run_obbt_iter>1){
            fixed_point[msname]=true;
        }
    }
    else{/*Update bounds*/
        if(dirk=="LB"){
                //Uncertainty in objk=obk+-solver_tolerance, here we choose lowest possible value in uncertainty interval
            if(!vk._is_relaxed){
                objk=std::max(objk-range_tol, boundk1);
            }
            vk.set_lb(keyk, objk);
            update_lb=true;
        }
        else{
                //Uncertainty in objk=obk+-solver_tolerance, here we choose highest possible value in uncertainty interval
            if(!vk._is_relaxed){
                objk=std::min(objk+range_tol, boundk1);
            }
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
        }
        /*If fixed point not reached for any variable, terminate is false*/
        if(!fixed_point[msname]){
            if(!vk._lift){
                terminate=false;
            }
        }
    }
    
        //If interval becomes smaller than range_tol, reset bounds so that interval=range_tol
    if(!vk._is_relaxed && std::abs(vk.get_ub(keyk)-vk.get_lb(keyk))<range_tol)
    {
            //If original interval is itself smaller than range_tol, do not have to reset interval
        if(interval_original.at(var_key_k)>=range_tol)
        {
            DebugOff("Entered reset"<<endl);
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
                //In the resized interval both original lower and upper bounds can not be crossed, because original interval is greater
                //than range_tol
        }
    }
    if(update_lb||update_ub){
        for(auto &mod:models){
            auto vkmod=mod->template get_var<T>(vkname);
            if(update_lb){
                vkmod.set_lb(keyk, vk.get_lb(keyk));
                DebugOff(vk.get_lb(keyk)<<endl);
            }
            if(update_ub){
                vkmod.set_ub(keyk, vk.get_ub(keyk));
                DebugOff(vk.get_ub(keyk)<<endl);
            }
        }
    }
    return true;
}
template<typename type>
template<typename T>
bool Model<type>::obbt_update_lagrange_bounds(std::vector<shared_ptr<gravity::Model<type>>>& models, shared_ptr<gravity::Model<type>>& obbt_model,   map<string, bool>& fixed_point,  const map<string, double>& interval_original, const map<string, double>& ub_original, const map<string, double>& lb_original, bool& terminate, int& fail, const double range_tol, const double fixed_tol_abs, const double fixed_tol_rel, const double zero_tol, int run_obbt_iter, std::map<int, double>& map_lb, std::map<int, double>& map_ub){
#ifdef USE_MPI
    int worker_id, nb_workers;
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
    auto err_size = MPI_Comm_size(MPI_COMM_WORLD, &nb_workers);
#endif
    auto id_old=0;
    auto in_it_start=obbt_model->_vars.begin()++;
    auto v_old=obbt_model->_vars.at(0);
    bool found=false;
    string vname, keyname, dirname;
    for(auto &it:map_ub){
        dirname="UB";
        found=false;
        for(auto in_it=in_it_start;in_it!=obbt_model->_vars.end();in_it++)
        {
            auto v = (*in_it).second;
            auto id=v->get_id();
            if(it.first>=id_old && it.first<id){
                vname=v_old->_name;
                auto inst=it.first-id_old;
                auto keys=v_old->get_keys();
                keyname=(*keys)[inst];
                in_it_start=in_it;
                found=true;
                break;
            }
            id_old=id;
            v_old=v;
        }
        if(found){
            this->obbt_update_bounds(false,it.second, vname+"|"+keyname+"|"+dirname,vname, keyname, dirname,  models,  obbt_model,   fixed_point,  interval_original,  ub_original,  lb_original, terminate, fail, range_tol, fixed_tol_abs, fixed_tol_rel,  zero_tol, run_obbt_iter);
        }else{
#ifdef USE_MPI
            if(worker_id==0)
#endif
                DebugOn("Not found ub" <<it.first<<" "<<it.second<<endl);
        }
        
    }
    id_old=0;
    in_it_start=obbt_model->_vars.begin()++;
    v_old=obbt_model->_vars.at(0);
    for(auto &it:map_lb){
        dirname="LB";
        found=false;
        for(auto in_it=in_it_start;in_it!=obbt_model->_vars.end();in_it++)
        {
            auto v = (*in_it).second;
            auto id=v->get_id();
            if(it.first>=id_old && it.first<id){
                vname=v_old->_name;
                auto inst=it.first-id_old;
                auto keys=v_old->get_keys();
                keyname=(*keys)[inst];
                in_it_start=in_it;
                found=true;
                break;
            }
            id_old=id;
            v_old=v;
        }
        if(found){
            this->obbt_update_bounds(false,it.second, vname+"|"+keyname+"|"+dirname,vname, keyname, dirname,  models,  obbt_model,   fixed_point,  interval_original,  ub_original,  lb_original, terminate, fail, range_tol, fixed_tol_abs, fixed_tol_rel,  zero_tol, run_obbt_iter);
        }else{
#ifdef USE_MPI
            if(worker_id==0)
#endif
                DebugOn("Not found lb" <<it.first<<" "<<it.second<<endl);
        }
        
    }
    return true;
}
template<typename type>
template<typename T>
void Model<type>::generate_lagrange_bounds(const std::vector<std::string> objective_models, std::vector<shared_ptr<gravity::Model<type>>>& models, shared_ptr<gravity::Model<type>>& obbt_model,   std::map<string, bool>& fixed_point,  const double range_tol, const double zero_tol, std::map<int, double>& map_lb, std::map<int, double>& map_ub){
    double vi_ub, vi_lb, vi_int;
    std::string viname, vikey, vidir;
    bool in_orig_model=false;
    for(auto &modeli:models){
        if (modeli->_status==0)
        {
            /* code */
            auto miname=modeli->_name;
            auto   mname=miname;
            std::size_t pos = mname.find("|");
            viname.assign(mname, 0, pos);
            mname=mname.substr(pos+1);
            pos=mname.find("|");
            vikey.assign(mname, 0, pos);
            vidir=mname.substr(pos+1);
            /*if!fixed_point in both dirs implies that the interval is not closed hence we check only in this case*/
            var<> vi=obbt_model->template get_var<T>(viname);
            vi_ub=vi.get_ub(vikey);
            vi_lb=vi.get_lb(vikey);
            vi_int=vi_ub-vi_lb;
            for (auto &vjp: modeli->_vars) {
                auto nb_inst = vjp.second->get_dim();
                auto ldvj= vjp.second->_l_dual;
                auto udvj=vjp.second->_u_dual;
                auto vjkeys=vjp.second->get_keys();
                auto vjname=vjp.second->_name;
                auto vjid=vjp.second->get_id();
                for(auto vpiter=0;vpiter<nb_inst;vpiter++)
                {
                    auto vjkey=(*vjkeys)[vpiter];
                    auto mjname_lb=vjname+"|"+(*vjkeys)[vpiter]+"|LB";
                    auto mjname_ub=vjname+"|"+(*vjkeys)[vpiter]+"|UB";
                    if(!(fixed_point[mjname_lb] && fixed_point[mjname_ub])){
                        if(udvj[vpiter] >= 1E-3)
                        {
                            var<> vj=modeli->template get_var<T>(vjname);
                            auto lb_vj=vj.get_lb(vjkey);
                            auto ub_vj=vj.get_ub(vjkey);
                            auto lb_vj_new=ub_vj-(vi_int)/udvj[vpiter];
                            if((lb_vj_new-lb_vj)>=range_tol){
                                if(map_lb.find(vjid+vpiter)!=map_lb.end()){
                                    if((lb_vj_new-map_lb.at(vjid+vpiter))>=range_tol){
                                        map_lb.at(vjid+vpiter)=lb_vj_new;
                                        DebugOff(vjname<<"\t"<<vjkey<<"\t"<<lb_vj_new<<" LB"<<endl);
                                    }
                                }
                                else{
                                    DebugOff(vjname<<"\t"<<vjkey<<"\t"<<lb_vj_new<<" LB"<<endl);
                                    map_lb[vjid+vpiter]=lb_vj_new;
                                }
                                
                            }
                        }
                        if(ldvj[vpiter] >= 1E-3 && vjname!="Wii")
                        {
                            var<> vj=modeli->template get_var<T>(vjname);
                            auto lb_vj=vj.get_lb(vjkey);
                            auto ub_vj=vj.get_ub(vjkey);
                            auto ub_vj_new=lb_vj+(vi_int)/ldvj[vpiter];
                            if((ub_vj-ub_vj_new)>=range_tol){
                                if(map_ub.find(vjid+vpiter)!=map_ub.end()){
                                    if((map_ub.at(vjid+vpiter)-ub_vj_new)>=range_tol){
                                        DebugOff(vjname<<"\t"<<vjkey<<"\t"<<ub_vj_new<<" UB"<<endl);
                                        map_ub.at(vjid+vpiter)=ub_vj_new;
                                    }
                                }
                                else{
                                    DebugOff(vjname<<"\t"<<vjkey<<"\t"<<ub_vj_new<<" UB"<<endl);
                                    map_ub[vjid+vpiter]=ub_vj_new;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

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
    
    var<> eta_i("eta_i", 0, 1);
    IIS->add(eta_i.in(range(0,n-1)));
    
    Constraint<> sum_eta("sum_eta");
    sum_eta=eta_int-sum(eta_i);
    IIS->add(sum_eta==0);
    
    int counter=0;
    
    for (auto &con: _cons_vec)
    {
        cs=con->get_nb_instances();
        Constraint<> Inter_con(con->_name);
        Inter_con=*con;
        
        if(con->_ctype==leq)
        {
            IIS->add(Inter_con<= eta_i.in(range(counter, counter+con->get_nb_instances()-1)));
        }
        else  if(con->_ctype==geq)
        {
            Inter_con += eta_i.in(range(counter, counter+con->get_nb_instances()-1));
            IIS->add(Inter_con>=0);
        }
        else  if(con->_ctype==eq)
        {
            Inter_con=(*con);
            IIS->add(Inter_con==0);
        }
        counter+=con->get_nb_instances();
    }
    IIS->print();
    return IIS;
}


#ifdef USE_MPI
/** function to run models in parallel across machines. Populates solution status and objective of the runs. If linearize cuts are also added, and each problem is refined upto nb_refine times
 @param[in] objective models: Vec with names of each model in models. If linearize order of elements in this vector is changed at the end of the function
 @param[in] sol_obj: Vec with objective values of each model in models
 @param[in] sol_obj: Vec with solution status of each model in models
 @param[in] models: Vec of models (per machine, of size threads_per_machine nr_threads)
 @param[in] relaxed model: original nonlinear nonconvex model
 @param[in] interioir model: model to compute interior point of relaxed model
 @param[in] cut_type: cut add strategy allvar or modelname
 @param[in] active tol: a constraint is violated and linearized as per this tolerance
 @param[in] stype: solver type
 @param[in] tol: tolerance with which each model in objective_models is solved
 @param[in] nr_threads: threads per machine
 @param[in] lin_solver: linear solver used by ipopt (ma27, ma57, etc)
 @param[in] max_iter, max_batch_time: max iter and batch_time for each time a problem is run
 @param[in] linearize: true if linear obbt algorithm is used
 @param[in] nb_refine: number of refinement steps, used for linearizs only
 @param[in] old_map: map of which worker_id each modelname is assigned to, used for linearize only
 @param[in] vbasis, cbasis: initialization of variable and constraint basis(defined by Gurobi documentation) for each model in models.
 @return returns true
 */

int run_MPI_new(std::vector<std::string>& objective_models, std::vector<double>& sol_obj, std::vector<int>& sol_status, std::vector<shared_ptr<gravity::Model<double>>>& models, const shared_ptr<gravity::Model<double>>& relaxed_model, const gravity::Model<double>& interior, string cut_type, double active_tol, gravity::SolverType stype, double tol, unsigned nr_threads, const string& lin_solver, int max_iter, int max_batch_time, bool linearize, int nb_refine, std::map<string,int>& old_map, vector<vector<double>>& vbasis, vector<std::map<string,double>>& cbasis, bool initialize_primal){
    int worker_id, nb_workers;
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
    auto err_size = MPI_Comm_size(MPI_COMM_WORLD, &nb_workers);
    auto nb_workers_ = std::min((size_t)nb_workers, objective_models.size());
    std::vector<std::string> objective_models_worker;
    std::vector<double> sol_obj_worker;
    std::vector<int> sol_status_worker;
    std::vector<size_t> limits;
    if(objective_models.size()!=0){
        if(!linearize){
            limits=bounds(nb_workers_, objective_models.size());
        }
        else{
            limits=bounds(nb_workers_, objective_models.size());
                //limits=bounds_reassign(nb_workers_, objective_models, old_map);
        }
        if(nb_workers_!=limits.size()-1){
            DebugOn("Error4 in computing limits");
        }
        sol_status.resize(objective_models.size(),-1);
        sol_obj.resize(objective_models.size(),-1.0);
        if(worker_id+1<limits.size()){
            /* Launch all threads in parallel */
            if(limits[worker_id] == limits[worker_id+1]){
                throw invalid_argument("limits[worker_id]==limits[worker_id+1]");
            }
            DebugOff("I'm worker ID: " << worker_id << ", I will be running models " << limits[worker_id] << " to " << limits[worker_id+1]-1 << endl);
            for (auto i = limits[worker_id]; i < limits[worker_id+1]; i++) {
                objective_models_worker.push_back(objective_models[i]);
            }
            sol_status_worker.resize(limits[worker_id+1]-limits[worker_id],0);
            sol_obj_worker.resize(limits[worker_id+1]-limits[worker_id],0);
            run_parallel_new(objective_models_worker, sol_obj_worker, sol_status_worker, models, relaxed_model, interior, cut_type, active_tol, stype, tol, nr_threads, lin_solver, max_iter, max_batch_time, linearize, nb_refine, vbasis, cbasis, initialize_primal);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        send_status_new(models,limits, sol_status);
        MPI_Barrier(MPI_COMM_WORLD);
        send_obj_all_new(models,limits, sol_obj);
    }
        //   MPI_Barrier(MPI_COMM_WORLD);
    return nb_workers_;
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
template bool Model<double>::root_refine(const Model<double>& interior_model, shared_ptr<Model<double>>& obbt_model, SolverType lb_solver_type, int nb_refine, const double upper_bound, double& lower_bound, const double ub_scale_value, double lb_solver_tol, double active_tol, int& oacuts, const double abs_tol, const double rel_tol, const double zero_tol, string lin_solver, int max_iter, int max_time, std::vector<double>& vrbasis, std::map<string,double>& crbasis, bool init);
template bool Model<double>::obbt_update_bounds(bool bound_converge,double objk, std::string msname,std::string vkname, std::string keyk, std::string dirk, std::vector<shared_ptr<gravity::Model<double>>>& models, shared_ptr<gravity::Model<double>>& obbt_model,   std::map<string, bool>& fixed_point, const map<string, double>& interval_original, const map<string, double>& ub_original, const map<string, double>& lb_original, bool& terminate, int& fail, const double range_tol, const double fixed_tol_abs, const double fixed_tol_rel, const double zero_tol, int run_obbt_iter);
template bool Model<double>::obbt_batch_update_bounds(const std::vector<std::string> objective_models, const std::vector<double>& sol_obj, const std::vector<int>& sol_status, std::vector<shared_ptr<gravity::Model<double>>>& models, shared_ptr<gravity::Model<double>>& obbt_model,   map<string, bool>& fixed_point,  const map<string, double>& interval_original, const map<string, double>& ub_original, const map<string, double>& lb_original, bool& terminate, int& fail, const double range_tol, const double fixed_tol_abs, const double fixed_tol_rel, const double zero_tol, int run_obbt_iter);
template void Model<double>::add_linear_row(Constraint<double>& con, int c_inst, const vector<double>& c_val, const double c0_val, const double scale);
template void Model<double>::generate_lagrange_bounds(const std::vector<std::string> objective_models, std::vector<shared_ptr<gravity::Model<double>>>& models, shared_ptr<gravity::Model<double>>& obbt_model,   std::map<string, bool>& fixed_point,  const double range_tol, const double zero_tol, std::map<int, double>& map_lb, std::map<int, double>& map_ub);
template bool Model<double>::obbt_update_lagrange_bounds(std::vector<shared_ptr<gravity::Model<double>>>& models, shared_ptr<gravity::Model<double>>& obbt_model,   map<string, bool>& fixed_point,  const map<string, double>& interval_original, const map<string, double>& ub_original, const map<string, double>& lb_original, bool& terminate, int& fail, const double range_tol, const double fixed_tol_abs, const double fixed_tol_rel, const double zero_tol, int run_obbt_iter, std::map<int, double>& map_lb, std::map<int, double>& map_ub);
template vector<vector<double>> Model<double>::cutting_planes_solution(const Model<double>& interior, double active_tol);
}
