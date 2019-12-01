
//
//  SpaceCraft.cpp
//  Gravity
//
//  Created by Hassan Hijazi on 27 Nov 2019.
//
//
#include <stdio.h>
#include <gravity/solver.h>
#include <gravity/rapidcsv.h>

using namespace std;
using namespace gravity;

int main (int argc, char * argv[])
{
    int output = 0;
    double tol = 1e-6;
    double solver_time_end, total_time_end, solve_time, total_time;
    string log_level="0";
    string fname_phi = string(prj_dir)+"/data_sets/SpaceCraft/phi_obs.csv";
    string fname_ne0 = string(prj_dir)+"/data_sets/SpaceCraft/ne0.csv";
    
    if(argc>2){
        fname_phi = argv[1];
        fname_ne0 = argv[2];
    }
    
    /* Declaring Parameters */
    auto ne0 = param<>("ne0");
    auto Phi_obs = param<>("ϕ_obs");

    /* Read Phi observation file */
    rapidcsv::Document  in_phi(fname_phi, rapidcsv::LabelParams(-1, -1));
    int n = in_phi.GetRowCount();
    for (int i = 0; i< n; i++) { // Input iterator
        Phi_obs.add_val(in_phi.GetCell<double>(0, i));
    }
    /* Read ne0 file */
    rapidcsv::Document  in_ne0(fname_ne0, rapidcsv::LabelParams(-1, -1));
    if(n!=in_ne0.GetRowCount()){
        throw invalid_argument("ne0 and phi_obs dimension mismatch, check input files");
    }
    for (int i = 0; i< n; i++) { // Input iterator
        ne0.add_val(in_ne0.GetCell<double>(0, i));
    }
    
    double total_time_start = get_wall_time();
    
    /* Declaring Model */
    Model<> M("Spacecraft Model");
    
    /* Declaring Indices */
    auto time = indices("time");
    time = range(1,n);
    Phi_obs.in(time);
    ne0.in(time);

    auto neg_time = Phi_obs.get_negative_ids();
    auto pos_time = Phi_obs.get_non_negative_ids();
    double me = 9.1e-31, e = 1.6e-19, rsc = 1., A = 0.3;
    
    /* Declaring Variables */
    auto Te = var<>("Te", 1/std::sqrt(10), 1/std::sqrt(0.1)); auto Tph = var<>("Tph",1/5,1/0.5); auto Jph = var<>("Jph",6.4e-6,1.6e-4);
    auto Phi = var<>("ϕ",10,20);
    double scale = 1e4;
    
    M.add(Te.in(time), Phi.in(time));
    M.add(Tph.in(R(1)), Jph.in(R(1)));
    /* Initializing variables */
    M.initialize_midpoint();
    Phi.copy_vals(Phi_obs);
    
    auto neg_ids = Tph.repeat_id(neg_time.size());
    auto pos_ids = Tph.repeat_id(pos_time.size());
    
    /* Declaring Constraints */
    Constraint<> Eq_neg("Eq_neg");
    Eq_neg = scale*A*4*pi*pow(rsc,2)*Jph.in(neg_ids)*Te.in(neg_time) - scale*e*4*pi*pow(rsc,2)*ne0.in(neg_time)*std::sqrt(e/(2*pi*me))*exp(Phi.in(neg_time)*pow(Te.in(neg_time),2));
    M.add(Eq_neg.in(neg_time)==0);
    
    Constraint<> Eq_pos("Eq_pos");
    Eq_pos = scale*A*4*pi*pow(rsc,2)*Jph.in(pos_ids)*(1+Phi.in(pos_time)*Tph.in(pos_ids))*exp(-1*Phi.in(pos_time)*Tph.in(pos_ids))*Te.in(pos_time) - scale*e*4*pi*pow(rsc,2)*ne0.in(pos_time)*std::sqrt(e/(2*pi*me))*(1+Phi.in(pos_time)*pow(Te.in(pos_time),2));
    M.add(Eq_pos.in(pos_time)==0);

    
/* Declaring Objective */
    func<> objective;
    for (int i = 0; i<n; i++) {
        objective += pow(Phi_obs(i) - Phi(i),2);
    }
    M.min(objective);
//    M.min(norm2(Phi_obs-Phi));
    
    /* Uncomment next line to print expanded model */
    M.print();
    
    /* Declaring Solver */
    solver<> NLSolver(M,ipopt);
    auto solver_time_start = get_wall_time();
    NLSolver.run(output = 5, tol = 1e-12);
    solver_time_end = get_wall_time();
    total_time_end = get_wall_time();
    solve_time = solver_time_end - solver_time_start;
    total_time = total_time_end - total_time_start;

    int precision = 14;
    Phi.print_vals(precision);
    DebugOn("Jph = " << to_string_with_precision(Jph.eval(),precision) << "."<<endl);
    DebugOn("Tph = " << to_string_with_precision(1./Tph.eval(), precision) << "."<<endl);
    auto Tinv = 1./Te;
    Tinv.eval_all();
    param<> Te_("Te");
    Te_ = Tinv;
    Te_.print_vals(precision);
    DebugOn("Optimal objective = " << to_string(M.get_obj_val()) << "."<<endl);
    DebugOn("Total wall clock time = " << to_string(total_time)<<endl);
    return 0;
}
