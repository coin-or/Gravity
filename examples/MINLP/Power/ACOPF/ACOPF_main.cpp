//
//  ACOPF.cpp
//  Gravity
//
//  Created by Guanglei Wang on 6/9/17.
//
//
#include <stdio.h>
#include <iostream>
#include <string>
#include <stdio.h>
#include <cstring>
#include <fstream>
#include "../PowerNet.h"
#include <gravity/model.h>
#include <gravity/solver.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;
using namespace gravity;
#define EPS 0.00001
#define DebugOn(x) cout << x
#define DebugOff(x)
//  Windows
#ifdef _WIN32
#include <Windows.h>
double get_wall_time() {
    LARGE_INTEGER time,freq;
    if (!QueryPerformanceFrequency(&freq)) {
        return 0;
    }
    if (!QueryPerformanceCounter(&time)) {
        return 0;
    }
    return (double)time.QuadPart / freq.QuadPart;
}
double get_cpu_time() {
    FILETIME a,b,c,d;
    if (GetProcessTimes(GetCurrentProcess(),&a,&b,&c,&d) != 0) {
        return
            (double)(d.dwLowDateTime |
                     ((unsigned long long)d.dwHighDateTime << 32)) * 0.0000001;
    } else {
        return 0;
    }
}
#else
#include <time.h>
#include <sys/time.h>
double get_wall_time() {
    struct timeval time;
    if (gettimeofday(&time,NULL)) {
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
double get_cpu_time() {
    return (double)clock() / CLOCKS_PER_SEC;
}
#endif

int main (int argc, const char * argv[])
{
    // ACOPF
    PowerNet* grid = new PowerNet();
    const char* fname;
//   fname = "../../data_sets/Power/nesta_case5_pjm.m";
   fname = "../../data_sets/Power/nesta_case14_ieee.m";
//   fname = "../../data_sets/Power/nesta_case1354_pegase.m";
//   fname = "../../data_sets/Power/nesta_case2383wp_mp.m";
    grid->readgrid(fname);

    // Grid Parameters
    unsigned nb_gen = grid->gens.size();
    unsigned nb_lines = grid->arcs.size();
    unsigned nb_buses = grid->nodes.size();

    /** build model */

    Model ACOPF("AC-OPF Model");
    /** Variables */
    // power generation
    var<double> Pg("Pg", grid->pg_min, grid->pg_max);
    var<double> Qg ("Qg", grid->qg_min, grid->qg_max);
    ACOPF.add_var(Pg^(nb_gen));
    ACOPF.add_var(Qg^(nb_gen));

    // power flow
    var<double> Pf_from("Pf_from");
    var<double> Qf_from("Qf_from");
    var<double> Pf_to("Pf_to");
    var<double> Qf_to("Qf_to");
    ACOPF.add_var(Pf_from^(nb_lines));
    ACOPF.add_var(Qf_from^(nb_lines));
    ACOPF.add_var(Pf_to^(nb_lines));
    ACOPF.add_var(Qf_to^(nb_lines));

    // voltage related variables.
    //    var<double> vr("vr");
    var<double> vr("vr", grid->v_max);
    var<double> vi("vi", grid->v_max);
    ACOPF.add_var(vr^(nb_buses));
    ACOPF.add_var(vi^(nb_buses));
    vr.initialize_all(1);

    /** Construct the objective function*/
    func_ obj = sum(grid->c0) + sum(grid->c1,Pg) + sum(grid->c2, power(Pg,2));
    obj.print();
    ACOPF.set_objective(min(obj));


    /** Define constraints */

    //KCL
    for (auto b: grid->nodes) {
        Bus* bus = (Bus*) b;
        Constraint KCL_P("KCL_P"+bus->_name);
        Constraint KCL_Q("KCL_Q"+bus->_name);

        /* Power Conservation */
        KCL_P  = sum(Pf_from.in(b->get_out())) + sum(Pf_to.in(b->get_in())) + bus->pl()- sum(Pg.in(bus->_gen));
        KCL_Q  = sum(Qf_from.in(b->get_out())) + sum(Qf_to.in(b->get_in())) + bus->ql()- sum(Qg.in(bus->_gen));
        
        /* Shunts */
        KCL_P +=  bus->gs()*(power(vr(bus->ID+1), 2) + power(vi(bus->ID+1), 2));
        KCL_Q -=  bus->bs()*(power(vr(bus->ID+1), 2) + power(vi(bus->ID+1), 2));

//        KCL_P.print();
        ACOPF.add_constraint(KCL_P = 0);
        ACOPF.add_constraint(KCL_Q = 0);
    }
    //AC Power Flow.
//
    Constraint Flow_P_From("Flow_P_From");
    Flow_P_From += Pf_from.in(grid->arcs);
    Flow_P_From -= grid->g_ff.in(grid->arcs)*(power(vr.from(grid->arcs), 2) + power(vi.from(grid->arcs), 2));
    Flow_P_From -= grid->g_ft.in(grid->arcs)*(vr.from(grid->arcs)*vr.to(grid->arcs) + vi.from(grid->arcs)*vi.to(grid->arcs));
    Flow_P_From -= grid->b_ft.in(grid->arcs)*(vi.from(grid->arcs)*vr.to(grid->arcs) - vr.from(grid->arcs)*vi.to(grid->arcs));
    Flow_P_From = 0;
//    Flow_P_From.print();
    ACOPF.add_constraint(Flow_P_From);

    Constraint Flow_P_To("Flow_P_To");
    Flow_P_To += Pf_to.in(grid->arcs);
    Flow_P_To -= grid->g_tt.in(grid->arcs)*(power(vr.to(grid->arcs), 2) + power(vi.to(grid->arcs), 2));
    Flow_P_To -= grid->g_tf.in(grid->arcs)*(vr.from(grid->arcs)*vr.to(grid->arcs) + vi.from(grid->arcs)*vi.to(grid->arcs));
    Flow_P_To -= grid->b_tf.in(grid->arcs)*(vr.from(grid->arcs)*vi.to(grid->arcs) - vr.to(grid->arcs)*vi.from(grid->arcs));
    Flow_P_To = 0;
    ACOPF.add_constraint(Flow_P_To);

    Constraint Flow_Q_From("Flow_Q_From");
    Flow_Q_From += Qf_from.in(grid->arcs);
    Flow_Q_From += grid->b_ff.in(grid->arcs)*(power(vr.from(grid->arcs), 2) + power(vi.from(grid->arcs), 2));
    Flow_Q_From += grid->b_ft.in(grid->arcs)*(vr.from(grid->arcs)*vr.to(grid->arcs) + vi.from(grid->arcs)*vi.to(grid->arcs));
    Flow_Q_From -= grid->g_ft.in(grid->arcs)*(vi.from(grid->arcs)*vr.to(grid->arcs) - vr.from(grid->arcs)*vi.to(grid->arcs));
    Flow_Q_From = 0;
    ACOPF.add_constraint(Flow_Q_From);

    Constraint Flow_Q_To("Flow_Q_To");
    Flow_Q_To += Qf_to.in(grid->arcs);
    Flow_Q_To += grid->b_tt.in(grid->arcs)*(power(vr.to(grid->arcs), 2) + power(vi.to(grid->arcs), 2));
    Flow_Q_To += grid->b_tf.in(grid->arcs)*(vr.from(grid->arcs)*vr.to(grid->arcs) + vi.from(grid->arcs)*vi.to(grid->arcs));
    Flow_Q_To -= grid->g_tf.in(grid->arcs)*(vr.from(grid->arcs)*vi.to(grid->arcs) - vr.to(grid->arcs)*vi.from(grid->arcs));
    Flow_Q_To = 0;
    ACOPF.add_constraint(Flow_Q_To);

    // AC voltage limit constraints.

    Constraint Vol_limit_UB("Vol_limit_UB");
    Vol_limit_UB = power(vr.in(grid->nodes), 2) + power(vi.in(grid->nodes), 2);
    Vol_limit_UB -= power(grid->v_max.in(grid->nodes), 2);
    ACOPF.add_constraint(Vol_limit_UB <= 0);

    Constraint Vol_limit_LB("Vol_limit_LB");
    Vol_limit_LB = power(vr.in(grid->nodes), 2) + power(vi.in(grid->nodes), 2);
    Vol_limit_LB -= power(grid->v_min.in(grid->nodes),2);
    ACOPF.add_constraint(Vol_limit_LB >= 0);

    
    /* REF BUS */
//    Constraint Ref_Bus("Ref_Bus");
//    Ref_Bus = vi(grid->arcs[grid->ref]);
//    ACOPF.add_constraint(Ref_Bus = 0);
    
    /* Phase Angle Bounds constraints */
    Constraint PAD_UB("PAD_UB");
    PAD_UB = vi.from(grid->arcs)*vr.to(grid->arcs) - vr.from(grid->arcs)*vi.to(grid->arcs);
    PAD_UB -= grid->tan_th_max.in(grid->arcs)*(vr.from(grid->arcs)*vr.to(grid->arcs) + vi.from(grid->arcs)*vi.to(grid->arcs));
    ACOPF.add_constraint(PAD_UB <= 0);
    
    Constraint PAD_LB("PAD_LB:");
    PAD_LB = vi.from(grid->arcs)*vr.to(grid->arcs) - vr.from(grid->arcs)*vi.to(grid->arcs);
    PAD_LB -= grid->tan_th_min.in(grid->arcs)*(vr.from(grid->arcs)*vr.to(grid->arcs) + vi.to(grid->arcs)*vi.from(grid->arcs));
    ACOPF.add_constraint(PAD_LB >= 0);


    /* Thermal Limit Constraints */
    Constraint Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from += power(Pf_from.in(grid->arcs), 2) + power(Qf_from.in(grid->arcs), 2);
    Thermal_Limit_from -= power(grid->S_max.in(grid->arcs), 2);
    ACOPF.add_constraint(Thermal_Limit_from <= 0);

    Constraint Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to += power(Pf_to.in(grid->arcs), 2) + power(Qf_to.in(grid->arcs), 2);
    Thermal_Limit_to -= power(grid->S_max.in(grid->arcs),2);
    ACOPF.add_constraint(Thermal_Limit_to <= 0);

    //solver OPF(ACOPF,cplex);
    solver OPF(ACOPF,ipopt);
    OPF.run();
    return 0;
}
