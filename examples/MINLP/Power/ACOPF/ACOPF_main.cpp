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
    fname = "/Users/guangleiwang/phD/kernel/Gravity/data_sets/Power/nesta_case5_pjm.m";
    grid->readgrid(fname);

    // Additional Parameters
    int nb_gen = grid->gens.size();
    int nb_lines = grid->arcs.size();
    int nb_buses = grid->nodes.size();

    /** build model */

    Model ACOPF("AC-OPF Model");

    /** Variables */
    // power generation
    var<double> Pg("Pg");
    var<double> Qg ("Qg");
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
    var<double> vr("vr");
    var<double> vi("vf");
    ACOPF.add_var(vr^(nb_buses));
    ACOPF.add_var(vi^(nb_buses));

    /** Construct the objective function*/
    func_ obj;
    for (auto g: grid->gens) {
        if (!g->_active)
            continue;
        obj += g->_cost->c0 + (g->_cost->c1)*Pg((g->ID + 1)) + g->_cost->c2*power(Pg(g->ID + 1), 2);
    }
    ACOPF.set_objective(min(obj));


    /** Define constraints */
    constant<int> ones(1);

    //KCL
    for (auto b: grid->nodes) {
        Bus* bus = (Bus*) b;
        Constraint KCL_P("KCL_P"+bus->_name);
        Constraint KCL_Q("KCL_Q"+bus->_name);

        KCL_P  = ones.tr()*Pf_from.in(b->get_out());
        KCL_Q  = ones.tr()*Qf_from.in(b->get_out());
        KCL_P += ones.tr()*Pf_to.in(b->get_in());
        KCL_Q += ones.tr()*Qf_to.in(b->get_in());

        KCL_P += bus->pl() + bus->gs()*(power(vr(bus->ID +1), 2) + power(vi(bus->ID + 1), 2));
        KCL_Q += bus->ql() - bus->bs()*(power(vr(bus->ID +1), 2) + power(vi(bus->ID + 1), 2));

        for (auto g: bus->_gen) {
            KCL_P -= Pg(g->ID + 1);
            KCL_Q -= Qg(g->ID + 1);
        }

        ACOPF.add_constraint(KCL_P = 0);
        ACOPF.add_constraint(KCL_Q = 0);
    }

    //AC Power Flow.
    param<double> g_ff("g_ff");
    param<double> g_ft("g_ft");
    param<double> g_tf("g_tf");
    param<double> g_tt("g_tt");

    param<double> b_ff("b_ff");
    param<double> b_ft("b_ft");
    param<double> b_tf("b_tf");
    param<double> b_tt("b_tt");

    for (auto a: grid->arcs) {
        auto la = (Line*) a;
        g_ff.add_val(la->g/pow(la->tr, 2.));
        g_ft.add_val((-la->g*la->cc + la->b*la->dd)/(pow(la->cc, 2) + pow(la->dd, 2)));
        b_ft.add_val((-la->b*la->cc - la->g*la->dd)/(pow(la->cc, 2) + pow(la->dd, 2)));

        g_tt.add_val(la->g);
        g_tf.add_val((-la->g*la->cc - la->b*la->dd)/(pow(la->cc, 2) + pow(la->dd, 2)));
        b_tf.add_val((-la->b*la->cc + la->g*la->dd)/(pow(la->cc, 2) + pow(la->dd, 2)));

        b_ff.add_val((la->ch/2 + la->b)/pow(la->tr, 2.));
        b_tt.add_val((la->ch/2 + la->b));
    }

    Constraint Flow_P_From("Flow_P_From");
    Flow_P_From += Pf_from.in(grid->arcs);
    Flow_P_From -= g_ff.in(grid->arcs)*(power(vr.from(grid->arcs), 2) + power(vi.from(grid->arcs), 2));
    Flow_P_From -= g_ft.in(grid->arcs)*(vr.from(grid->arcs)*vr.to(grid->arcs) + vi.from(grid->arcs)*vi.to(grid->arcs));
    Flow_P_From -= b_ft.in(grid->arcs)*(vi.from(grid->arcs)*vr.to(grid->arcs) - vr.from(grid->arcs)*vi.to(grid->arcs));
    Flow_P_From = 0;
    //ACOPF.add_constraint(Flow_P_From);

    Constraint Flow_P_To("Flow_P_To");
    Flow_P_To += Pf_to.in(grid->arcs);
    Flow_P_To -= g_tt.in(grid->arcs)*(power(vr.to(grid->arcs), 2) + power(vi.to(grid->arcs), 2));
    Flow_P_To -= g_tf.in(grid->arcs)*(vr.from(grid->arcs)*vr.to(grid->arcs) + vi.from(grid->arcs)*vi.to(grid->arcs));
    Flow_P_To -= b_tf.in(grid->arcs)*(vr.from(grid->arcs)*vi.to(grid->arcs) - vr.to(grid->arcs)*vi.from(grid->arcs));
    Flow_P_To = 0;
    // ACOPF.add_constraint(Flow_P_To);

    Constraint Flow_Q_From("Flow_Q_From");
    Flow_Q_From += Qf_from.in(grid->arcs);
    Flow_Q_From += b_ff.in(grid->arcs)*(power(vr.from(grid->arcs), 2) + power(vi.from(grid->arcs), 2));
    Flow_Q_From += b_ft.in(grid->arcs)*(vr.from(grid->arcs)*vr.to(grid->arcs) + vi.from(grid->arcs)*vi.to(grid->arcs));
    Flow_Q_From -= g_ft.in(grid->arcs)*(vi.from(grid->arcs)*vr.to(grid->arcs) - vr.from(grid->arcs)*vi.to(grid->arcs));
    Flow_Q_From = 0;
    //  ACOPF.add_constraint(Flow_Q_From);

    Constraint Flow_Q_To("Flow_Q_To");
    Flow_Q_To += Qf_to.in(grid->arcs);
    Flow_Q_To += b_tt.in(grid->arcs)*(power(vr.to(grid->arcs), 2) + power(vi.to(grid->arcs), 2));
    Flow_Q_To -= b_tf.in(grid->arcs)*(vr.from(grid->arcs)*vr.to(grid->arcs) + vi.from(grid->arcs)*vi.to(grid->arcs));
    Flow_Q_To -= g_tf.in(grid->arcs)*(vr.from(grid->arcs)*vi.to(grid->arcs) - vr.to(grid->arcs)*vi.from(grid->arcs));
    Flow_Q_To = 0;
//   ACOPF.add_constraint(Flow_Q_To);

    // AC voltage limit constraints.
    param<double> vbound_max_square("vbound_max_square");
    param<double> vbound_min_square("vbound_min_square");

    for (auto b: grid->nodes) {
        vbound_max_square.add_val(pow(((Bus*)b)->vbound.max, 2.));
        vbound_min_square.add_val(pow(((Bus*)b)->vbound.min, 2.));
    }

    Constraint Vol_limit_UB("Vol_limit_UB");
    Vol_limit_UB = power(vr.in(grid->nodes), 2) + power(vi.in(grid->nodes), 2);
    Vol_limit_UB -= vbound_max_square.in(grid->nodes);
    ACOPF.add_constraint(Vol_limit_UB <= 0);

    Constraint Vol_limit_LB("Vol_limit_LB");
    Vol_limit_LB = power(vr.in(grid->nodes), 2) + power(vi.in(grid->nodes), 2);
    Vol_limit_LB -= vbound_max_square.in(grid->nodes);
//   ACOPF.add_constraint(Vol_limit_LB >= 0);

    //AC-PAD constraints
    param<double> tbound_max_tan("vbound_max_tan");
    param<double> tbound_min_tan("vbound_min_tan");

    for (auto a: grid->arcs) {
        tbound_max_tan.add_val(tan(((Line*)a)->tbound.max));
        tbound_min_tan.add_val(tan(((Line*)a)->tbound.min));
    }

    Constraint PAD_UB("PAD_UB");
    PAD_UB = vr.from(grid->arcs)*vi.to(grid->arcs) + vr.to(grid->arcs)*vi.from(grid->arcs);
    PAD_UB -= tbound_max_tan.in(grid->arcs)*(vr.from(grid->arcs)*vr.to(grid->arcs) - vi.to(grid->arcs)*vi.from(grid->arcs));
    //ACOPF.add_constraint(PAD_UB <= 0);


    Constraint PAD_LB("PAD_LB:");
    PAD_LB = vr.from(grid->arcs)*vi.to(grid->arcs) + vr.to(grid->arcs)*vi.from(grid->arcs);
    PAD_LB -= tbound_min_tan.in(grid->arcs)*(vr.from(grid->arcs)*vr.to(grid->arcs) - vi.to(grid->arcs)*vi.from(grid->arcs));
    //ACOPF.add_constraint(PAD_LB >= 0);


    // Thermal_Limit {(l,i,j) in arcs}: p[l,i,j]^2 + q[l,i,j]^2 <= s[l]^2;*/
    param<double> Thermal_limit_square("Thermal_limit_square");
    for (auto a: grid->arcs) {
        Line* la = (Line *) a;
        Thermal_limit_square.add_val(pow(la->limit, 2));
    }

    Constraint Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from += power(Pf_from.in(grid->arcs), 2) + power(Qf_from.in(grid->arcs), 2);
    Thermal_Limit_from -= Thermal_limit_square;
    ACOPF.add_constraint(Thermal_Limit_from <= 0);

    Constraint Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_from += power(Pf_to.in(grid->arcs), 2) + power(Qf_to.in(grid->arcs), 2);
    Thermal_Limit_from -= Thermal_limit_square;
    ACOPF.add_constraint(Thermal_Limit_to <= 0);

    // Power generation constraints.
    param<double> PUB("PUB");
    param<double> PLB("PUB");
    param<double> QUB("QUB");
    param<double> QLB("QUB");

    for (auto g: grid->gens) {
        PUB.add_val(g->pbound.max);
        PLB.add_val(g->pbound.min);
        QUB.add_val(g->qbound.max);
        QLB.add_val(g->qbound.min);
    }

    Constraint Pbound_UB("Pbound_UB");
    Constraint Pbound_LB("Pbound_LB");
    Constraint Qbound_UB("Qbound_UB");
    Constraint Qbound_LB("Qbound_LB");

    Pbound_UB = Pg.in(grid->gens) - PUB.in(grid->gens);
    Pbound_LB = Pg.in(grid->gens) - PLB.in(grid->gens);
    Qbound_UB = Qg.in(grid->gens) - QUB.in(grid->gens);
    Qbound_UB = Qg.in(grid->gens) - QLB.in(grid->gens);

    ACOPF.add_constraint(Pbound_UB <= 0);
    ACOPF.add_constraint(Pbound_LB >= 0);
    ACOPF.add_constraint(Qbound_UB <= 0);
    ACOPF.add_constraint(Qbound_LB >= 0);

    solver OPF(ACOPF,cplex);
    OPF.run();
    return 0;
}
