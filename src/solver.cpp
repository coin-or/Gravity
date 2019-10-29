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
    shared_ptr<Model<type>> Model<type>::buildOA(int nb_discr, int nb_perturb)
    {
        
        
        //    this->print_solution();
        vector<double> xsolution(_nb_vars);
        vector<double> xinterior(_nb_vars);
        vector<double> xcurrent;
        get_solution(xsolution);
        
        
        auto OA=make_shared<Model<>>(_name+"-OA Model");
        for (auto &it: _vars)
        {
            auto v = it.second;
            if(!OA->has_var(*v)){
                OA->add_var(v);
            }
        }
        auto obj=*_obj;
        if(_objt==minimize){
            OA->min(obj);
        }
        else {
            OA->max(obj);
        }
        string cname;
        for (auto &con: _cons_vec)
        {
            if(!con->is_linear()) {
                if(con->_ctype==eq)
                {
                    DebugOn("Exception: Equality constraint is not currently supported"<<endl);
                    DebugOn("Throw exception" <<endl);
                    
                }
                else
                {
                    
                    OA->add_outer_app_uniform(nb_discr, *con);
                    
                    set_solution(xsolution);
                }
                
            }
            else
            {
                OA->add(*con);
            }
        }
        
        OA->add_outer_app_active(*this, nb_perturb);
        
        set_solution(xsolution);
        return OA;
    }
    
    
    /** Returns an interior point of a model
     @param[in] nonlin: model for which interior point with respect to nonlinear constraints is computed
     Assuming model has no nonlinear equality constraints
     **/
    template<typename type>
    template<typename T>
    shared_ptr<Model<type>> Model<type>::build_model_interior(Model<type> nonlin)
    {
        auto Interior = make_shared<Model<>>(nonlin._name+"Interior");
        
        
        for (auto &it: nonlin._vars)
        {
            auto v = it.second;
            if(!Interior->has_var(*v)){
                Interior->add_var(v);
            }
        }
        var<> eta_int("eta_int", -1, -0.000001);
        
        Interior->add(eta_int.in(range(0,0)));
        auto obj=eta_int;
        
        Interior->min(obj);
        
        
        
        for (auto &con: nonlin._cons_vec)
        {
            if(!con->is_linear()) {
                 if(!con->is_convex() || con->is_rotated_soc() || con->check_soc()){
            Constraint<> Inter_con(*con);
                
            if(con->_ctype==leq)
            {
                Inter_con -= eta_int;
                Interior->add(Inter_con<=0);
            }
            else  if(con->_ctype==geq)
            {
                Inter_con += eta_int;
                Interior->add(Inter_con>=0);
            }
                   }
            }
            //        else
            //        {
            //            if(con->_ctype==leq||con->_ctype==geq)
            //            Interior->add(*con);
            //        }
        }
        //    Interior->print();
        return Interior;
    }
    
    
    
    
    /** Runds models stored in the vector in parallel, using solver of stype and tolerance tol */
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
    
    /** Discretizes Constraint con and adds OA cuts to the model that calls it. Discretization of squared constraint only currently implemented
     @param[in] nb_discr:
     @param[in] con:
     @return void. OA cuts are added to the model that calls the function (for all func instances)
     **/
    template<>
    void Model<>::add_outer_app_uniform(int nb_discr, Constraint<> con)
    {
        
        func<> res;
        double lb,ub;
        size_t posv;
        if(con.is_quadratic() && con._lterms->size()==1 && con._qterms->size()==1 && con._qterms->begin()->second._p->first==con._qterms->begin()->second._p->second) //This if is specific to constraints of the form ay- x^2 or x^2-ay
        {
            
            auto x=con._qterms->begin()->second._p->first;
            
            for(auto d=0;d<nb_discr;d++)
            {
                for(auto i=0;i<con.get_nb_inst();i++)
                {
                    posv=x->get_id_inst(i);
                    lb=x->get_double_lb(posv);
                    ub=x->get_double_ub(posv);
                    x->set_double_val(posv, lb+d*(ub-lb)/nb_discr);
                }
                Constraint<> OA_uniform("OA_cuts_uniform "+con._name+to_string(d));
                OA_uniform=con.get_outer_app_squared();
                if(con._ctype==leq) {
                    add(OA_uniform<=0);
                }
                else {
                    add(OA_uniform>=0);
                }
                
            }
        }/*TODO Else (discretization for general constraint)*/
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
//                while(c>0)
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
    //
    /** Outer approximation of active (nonlinear) constraints of the model
     @param[in] nonlin: original nonlinear model at whose solution (at the active point) OA cuts are added:
     @param[in] nb_perturb:
     @return void. OA cuts are added to the model that calls the function (for all func instances) at the solution and by perturbing the solution
     Assumption: nonlinear constraint to be perturbed does not have any vector variables
     **/
    template<>
    void Model<>::add_outer_app_active(Model<> nonlin, int nb_perturb)
    {
        const double active_tol=1e-6,active_tol_sol=1e-8, perturb_dist=1e-1;
        vector<double> xsolution(_nb_vars);
        vector<double> xactive, xcurrent, xinterior, xres;
        double a,b,c;
        bool interior=false;
        double fk;
        bool outer;
        int counter=0;
        size_t posv;
        
        int output=0;
        double tol=1e-8;
        bool convex_fr=true;
        bool scale=false;
        auto Ointerior=build_model_interior(nonlin);
        solver<> modelI(Ointerior, ipopt);
        modelI.run(output, tol);
        
        if((Ointerior->_status==0||Ointerior->_status==1) && Ointerior->get_obj_val() <0)
        {
            interior=true;
            //        Ointerior->print_solution();
        }
        
        for (auto &con: nonlin._cons_vec)
        {
            //if(!con->is_linear() && !con->is_convex()) {
                 if(!con->is_linear()) {
                
                con->uneval();
                Constraint<> OA_sol("OA_cuts_solution_"+con->_name);
                indices activeset("active_"+con->_name);
                for(auto i=0;i<con->get_nb_inst();i++){
                    if(std::abs(con->eval(i))<=active_tol_sol || (con->is_convex() && !con->is_rotated_soc() && !con->check_soc())){
                        auto keys=con->_indices->_keys;
                        activeset.add((*keys)[i]);
                    }
                }
                OA_sol=con->get_outer_app();
                if(con->_ctype==leq) {
                    add(OA_sol.in(activeset)<=0);
                }
                else {
                    add(OA_sol.in(activeset)>=0);
                }
            }
        }
       
        if(interior)
        {
            get_solution(xsolution);
            for (auto &con: nonlin._cons_vec)
            {
                if(!con->is_linear()) {
                    for(auto i=0;i<con->get_nb_inst();i++){
                        con->uneval();
                        if (!con->is_convex() || con->is_rotated_soc() || con->check_soc()){
                            auto cname=con->_name;
                            auto con_interior=Ointerior->get_constraint(cname);
                            xinterior=con_interior->get_x(i);
                            xinterior.pop_back();
                            xcurrent=con->get_x(i);
                            if(std::abs(con->eval(i))<=active_tol_sol){
                                xactive=xcurrent;
                            }
                            else
                            {
                                auto res=con->get_any_active_point(i, con->_ctype);
                                if(res.first){
                                    xactive=res.second;
                                    
                                }
                                else{
                                    continue;
                                }
                            }
                            con->set_x(i, xactive);
                            
                            for(auto j=1;j<=nb_perturb;j++)
                            {
                                counter=0;
                                for(auto &it: *(con->_vars))
                                {
                                    auto v = it.second.first;
                                    auto vname=v->_name;
                                    if(v->_is_vector)
                                    {
                                        DebugOn("Exception: Vector variables are not currently supported"<<endl);
                                        DebugOn("Throw exception" <<endl);
                                        break;
                                    }
                                    else
                                    {
                                        outer=false;
                                        for(auto k=-1;k<=1;k+=2)
                                        {
                                            posv=v->get_id_inst(i);
                                            v->set_double_val(posv, xactive[counter]*(1+k*j*perturb_dist));
                                            con->uneval();
                                            fk=con->eval(i);
                                            if((fk>=active_tol && con->_ctype==leq) || (fk<=(active_tol*(-1)) && con->_ctype==geq)){
                                                outer=true;
                                                break;
                                            }
                                        }
                                        if(outer)
                                        {
                                            auto res_search=con->linesearchbinary(xinterior, i, con->_ctype);
                                            if(res_search){
                                                convex_fr=true;
                                                if(!con->is_convex() && !con->is_rotated_soc() && !con->check_soc()) //assuming con is the SDP cut as it is the only nonconvex one
                                                {
                                                    xres=con->get_x(i);
                                                    con->uneval();
                                                    fk=con->eval(i);
                                                    a=std::pow(xres[0],2)+std::pow(xres[3],2)-xres[6]*xres[7];
                                                    b=std::pow(xres[1],2)+std::pow(xres[4],2)-xres[7]*xres[8];
                                                    c=std::pow(xres[2],2)+std::pow(xres[5],2)-xres[6]*xres[8];
                                                    if(a<=0 && b<=0 && c<=0){
                                                        convex_fr=true;
                                                    }
                                                    else{
                                                        convex_fr=false;
                                                    }
                                                   }
                                            
                                                if(convex_fr){
                                                
                                                Constraint<> OA_active("OA_active_"+con->_name+"_"+to_string(i)+"_"+to_string(j)+"_"+v->_name);
                                                OA_active=con->get_outer_app_insti(i, false);
                                                if(con->_ctype==leq) {
                                                    add(OA_active<=0);
                                                }
                                                else {
                                                    add(OA_active>=0);
                                                }
                                                }
                                                
                                            }
                                            
                                        }
                                    }
                                    con->set_x(i, xactive);
                                    counter++;
                                }
                            }
                            con->set_x(i, xcurrent);
                            xcurrent.clear();
                            xactive.clear();
                            xinterior.clear();
                            
                        }
                    }
                }
            }
        }
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
        
        
       
        DebugOn("N"<<n<<endl);
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
    template shared_ptr<Model<double>> Model<double>::buildOA(int,int);
    template shared_ptr<Model<double>> Model<double>::build_model_interior(Model<double>);
    template shared_ptr<Model<double>> Model<double>::build_model_IIS();
}




