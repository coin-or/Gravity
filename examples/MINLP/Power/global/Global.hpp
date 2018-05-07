//
//  Global.hpp
//  projects
//
//  Created by Guanglei Wang on 5/3/18.
//
//
#ifndef Global_hpp
#define Global_hpp
#include "Partition.hpp"


class Global {
public:
    PowerNet* grid; 
    Partition* P_;
    int Num_parts; //  spatial decomposition
    int Num_time; // time decomposition

    // Schedule Parameters
    param<Real> rate_ramp;
    param<Real> rate_switch;
    param<int> min_up;
    param<int> min_down;
    param<Real> cost_up;
    param<Real> cost_down;
    param<Real> Pg_initial;
    param<bool> On_off_initial;
    // Variables 
    vector<var<Real>> Pg;
    vector<var<Real>> Pg2; // new var introduced for the perspective formulation.
    vector<var<Real>> Qg;

    //Lifted variables.
    vector<var<Real>> R_Xij;
    vector<var<Real>> Im_Xij;
    vector<var<Real>> Xii;
    // Commitment variables
    vector<var<bool>> On_off; 
    vector<var<bool>> Start_up;
    vector<var<bool>> Shut_down;
    // sol
    vector<param<bool>> On_off_sol_; 
    vector<param<Real>> Pg_sol_;
    vector<param<bool>> Start_up_sol_;
    vector<param<bool>> Shut_down_sol_;
    // power flow vars;
    vector<var<Real>> Pf_from;
    vector<var<Real>> Qf_from;
    vector<var<Real>> Pf_to;
    vector<var<Real>> Qf_to;

   // multipliers time
    param<Real> lambda_up; // inter temporal: start up and shut down constraints
    param<Real> lambda_down;
    param<Real> zeta_up; // ramping constraints
    param<Real> zeta_down;
    param<Real> mu; // dual of min up down constraints
    param<Real> mu_up;
    param<Real> mu_down;
    bool include_min_updown_ =false;
    
    // multipliers spatial
    param<Real> R_lambda_;
    param<Real> Im_lambda_;
    param<Real> lambda_;
    // vals of each subproblem
    vector<double> Sub_;
    // SOCP constraints
    vector<Constraint> SOC_;
    vector<shared_ptr<Constraint>> SOC_outer_;
    
    //Constructors
    Global();
    Global(PowerNet*, int parts, int T);
    ~Global();
    
    // Accessors
    double getdual_relax_time_(bool include);
    double LR_bound_time_(bool included_min_up_down);
    double Upper_bound_sequence_(bool included_min_up_down);
    double Subproblem_time_(int l);
    double Subproblem_upper_time_(int l);
    void add_var_Sub_time(Model&, int t);
    void add_obj_Sub_time(Model&, int t);
    void add_obj_Sub_upper_time(Model&, int t);
    void add_SOCP_Sub_time(Model&, int t);
    void add_SOCP_Outer_Sub_time(Model&, int t);
    void add_KCL_Sub_time(Model&, int t);
    void add_thermal_Sub_time(Model& , int t);
    void add_perspective_OnOff_Sub_time(Model&, int t);
    void add_MC_upper_Sub_time(Model&, int t);

    void add_MC_intertemporal_Sub_upper_time(Model&, int t);
    void add_OnOff_Sub_upper_time(Model&, int t);
    void add_Ramp_Sub_upper_time(Model&, int t);
    void add_minupdown_Sub_upper_time(Model&, int t);
    
    // check rank 1 constraint.
    void check_rank1_constraint(Model& Sub, int t);

    double getdual_relax_spatial();
    double LR_bound_spatial_();
    double Subproblem_spatial_(int l);
};
#endif /* Global_hpp */

