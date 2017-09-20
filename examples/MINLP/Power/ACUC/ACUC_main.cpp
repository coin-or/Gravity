//
//  ACUC.cpp
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
    // ACUC
    PowerNet* grid = new PowerNet();
    const char* fname;
    fname = "/Users/guangleiwang/phD/kernel/Gravity/data_sets/Power/nesta_case5_pjm.m";
    grid->readgrid(fname);

    // Additional Parameters
    int nb_gen = grid->gens.size();
    int nb_lines = grid->arcs.size();
    int nb_buses = grid->nodes.size();

    /** build model */
    int T = 2;
    Model ACUC("ACUC Model");

    /** Variables */
    // power generation
    var<double> Pg("Pg");
    var<double> Qg ("Qg");
    ACUC.add_var(Pg^(T*nb_gen));
    ACUC.add_var(Qg^(T*nb_gen));

    // power flow
    var<double> Pf_from("Pf_from");
    var<double> Qf_from("Qf_from");
    var<double> Pf_to("Pf_to");
    var<double> Qf_to("Qf_to");
    ACUC.add_var(Pf_from^(T*nb_lines));
    ACUC.add_var(Qf_from^(T*nb_lines));
    ACUC.add_var(Pf_to^(T*nb_lines));
    ACUC.add_var(Qf_to^(T*nb_lines));

    // voltage related variables.
    var<double> vr("vr");
    var<double> vi("vf");
    ACUC.add_var(vr^(T*nb_buses));
    ACUC.add_var(vi^(T*nb_buses));

    /** Construct the objective function*/
    func_ obj;
    for (int t = 0; t < T; t++)
        for (auto g: grid->gens) {
            if (!g->_active)
                continue;
            obj += g->_cost->c0 + (g->_cost->c1)*Pg((g->ID + 1), t) + g->_cost->c2*power(Pg(g->ID + 1, t), 2);
        }
    ACUC.set_objective(min(obj));


    /** Define constraints */
    constant<int> ones(1);

    //KCL
    for (int t  = 0; t < T; t++) {
        for (auto b: grid->nodes) {
            Bus* bus = (Bus*) b;
            Constraint KCL_P("KCL_P"+bus->_name + to_string(t));
            Constraint KCL_Q("KCL_Q"+bus->_name + to_string(t));

            KCL_P  = ones.tr()*Pf_from.in(b->get_out(), t);
            KCL_Q  = ones.tr()*Qf_from.in(b->get_out(), t);
            KCL_P += ones.tr()*Pf_to.in(b->get_in(), t);
            KCL_Q += ones.tr()*Qf_to.in(b->get_in(), t);

            KCL_P += bus->pl() + bus->gs()*(power(vr(bus->ID +1, t), 2) + power(vi(bus->ID + 1, t), 2));
            KCL_Q += bus->ql() - bus->bs()*(power(vr(bus->ID +1, t), 2) + power(vi(bus->ID + 1, t), 2));

            for (auto g: bus->_gen) {
                KCL_P -= Pg(g->ID + 1, t);
                KCL_Q -= Qg(g->ID + 1, t);
            }

            ACUC.add_constraint(KCL_P = 0);
            ACUC.add_constraint(KCL_Q = 0);
        }
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

    for (int t = 0; t < T; t++) {
        Constraint Flow_P_From("Flow_P_From" + to_string(t));
        Flow_P_From += Pf_from.in(grid->arcs, t);
        Flow_P_From -= g_ff.in(grid->arcs)*(power(vr.from(grid->arcs, t), 2) + power(vi.from(grid->arcs, t), 2));
        Flow_P_From -= g_ft.in(grid->arcs)*(vr.from(grid->arcs, t)*vr.to(grid->arcs, t) + vi.from(grid->arcs, t)*vi.to(grid->arcs, t));
        Flow_P_From -= b_ft.in(grid->arcs)*(vi.from(grid->arcs, t)*vr.to(grid->arcs, t) - vr.from(grid->arcs, t)*vi.to(grid->arcs, t));
        Flow_P_From = 0;
        //ACUC.add_constraint(Flow_P_From);
    }

    for (int t = 0; t < T; t++) {
        Constraint Flow_P_To("Flow_P_To"+ to_string(t));
        Flow_P_To += Pf_to.in(grid->arcs, t);
        Flow_P_To -= g_tt.in(grid->arcs)*(power(vr.to(grid->arcs, t), 2) + power(vi.to(grid->arcs, t), 2));
        Flow_P_To -= g_tf.in(grid->arcs)*(vr.from(grid->arcs, t)*vr.to(grid->arcs, t) + vi.from(grid->arcs, t)*vi.to(grid->arcs, t));
        Flow_P_To -= b_tf.in(grid->arcs)*(vr.from(grid->arcs, t)*vi.to(grid->arcs, t) - vr.to(grid->arcs, t)*vi.from(grid->arcs, t));
        Flow_P_To = 0;
        //ACUC.add_constraint(Flow_P_To);
    }

    for (int t = 0; t < T; t++) {
        Constraint Flow_Q_From("Flow_Q_From" + to_string(t));
        Flow_Q_From += Qf_from.in(grid->arcs, t);
        Flow_Q_From += b_ff.in(grid->arcs)*(power(vr.from(grid->arcs, t), 2) + power(vi.from(grid->arcs, t), 2));
        Flow_Q_From += b_ft.in(grid->arcs)*(vr.from(grid->arcs, t)*vr.to(grid->arcs, t) + vi.from(grid->arcs, t)*vi.to(grid->arcs, t));
        Flow_Q_From -= g_ft.in(grid->arcs)*(vi.from(grid->arcs, t)*vr.to(grid->arcs, t) - vr.from(grid->arcs, t)*vi.to(grid->arcs, t));
        Flow_Q_From = 0;
        //ACUC.add_constraint(Flow_Q_From);
    }

    for (int t = 0; t < T; t++) {
        Constraint Flow_Q_To("Flow_Q_To");
        Flow_Q_To += Qf_to.in(grid->arcs, t);
        Flow_Q_To += b_tt.in(grid->arcs)*(power(vr.to(grid->arcs, t), 2) + power(vi.to(grid->arcs, t), 2));
        Flow_Q_To -= b_tf.in(grid->arcs)*(vr.from(grid->arcs, t)*vr.to(grid->arcs, t) + vi.from(grid->arcs, t)*vi.to(grid->arcs, t));
        Flow_Q_To -= g_tf.in(grid->arcs)*(vr.from(grid->arcs, t)*vi.to(grid->arcs, t) - vr.to(grid->arcs, t)*vi.from(grid->arcs, t));
        Flow_Q_To = 0;
        //ACUC.add_constraint(Flow_Q_To);
    }


    // AC voltage limit constraints.
    param<double> vbound_max_square("vbound_max_square");
    param<double> vbound_min_square("vbound_min_square");

    for (auto b: grid->nodes) {
        vbound_max_square.add_val(pow(((Bus*)b)->vbound.max, 2.));
        vbound_min_square.add_val(pow(((Bus*)b)->vbound.min, 2.));
    }

    for (int t = 0; t < T; t++) {
        Constraint Vol_limit_UB("Vol_limit_UB" + to_string(t));
        Vol_limit_UB = power(vr.in(grid->nodes, t), 2) + power(vi.in(grid->nodes, t), 2);
        Vol_limit_UB -= vbound_max_square.in(grid->nodes);
        ACUC.add_constraint(Vol_limit_UB <= 0);
    }

    for (int t = 0; t < T; t++) {
        Constraint Vol_limit_LB("Vol_limit_LB" + to_string(t));
        Vol_limit_LB = power(vr.in(grid->nodes, t), 2) + power(vi.in(grid->nodes, t), 2);
        Vol_limit_LB -= vbound_max_square.in(grid->nodes);
       // ACUC.add_constraint(Vol_limit_LB >= 0);
    }

    //AC-PAD constraints
    param<double> tbound_max_tan("vbound_max_tan");
    param<double> tbound_min_tan("vbound_min_tan");

    for (auto a: grid->arcs) {
        tbound_max_tan.add_val(tan(((Line*)a)->tbound.max));
        tbound_min_tan.add_val(tan(((Line*)a)->tbound.min));
    }

    for (int t = 0; t < T; t++) {
        Constraint PAD_UB("PAD_UB" + to_string(t));
        PAD_UB = vr.from(grid->arcs, t)*vi.to(grid->arcs, t) + vr.to(grid->arcs, t)*vi.from(grid->arcs, t);
        PAD_UB -= tbound_max_tan.in(grid->arcs)*(vr.from(grid->arcs, t)*vr.to(grid->arcs, t) - vi.to(grid->arcs, t)*vi.from(grid->arcs, t));
        //ACUC.add_constraint(PAD_UB <= 0);

        Constraint PAD_LB("PAD_LB:");
        PAD_LB = vr.from(grid->arcs, t)*vi.to(grid->arcs, t) + vr.to(grid->arcs, t)*vi.from(grid->arcs, t);
        PAD_LB -= tbound_min_tan.in(grid->arcs, t)*(vr.from(grid->arcs, t)*vr.to(grid->arcs, t) - vi.to(grid->arcs, t)*vi.from(grid->arcs, t));
       // ACUC.add_constraint(PAD_LB >= 0);
    }

    // Thermal_Limit {(l,i,j) in arcs}: p[l,i,j]^2 + q[l,i,j]^2 <= s[l]^2;*/
    param<double> Thermal_limit_square("Thermal_limit_square");
    for (auto a: grid->arcs) {
        Line* la = (Line *) a;
        Thermal_limit_square.add_val(pow(la->limit, 2));
    }

    for (int t = 0; t < T; t++) {
        Constraint Thermal_Limit_from("Thermal_Limit_from" + to_string(t));
        Thermal_Limit_from += power(Pf_from.in(grid->arcs, t), 2) + power(Qf_from.in(grid->arcs, t), 2);
        Thermal_Limit_from -= Thermal_limit_square;
        ACUC.add_constraint(Thermal_Limit_from <= 0);

        Constraint Thermal_Limit_to("Thermal_Limit_to" + to_string(t));
        Thermal_Limit_from += power(Pf_to.in(grid->arcs, t), 2) + power(Qf_to.in(grid->arcs, t), 2);
        Thermal_Limit_from -= Thermal_limit_square;
        ACUC.add_constraint(Thermal_Limit_to <= 0);
    }

//
//    // Power generation constraints.
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

    for (int t = 0; t < T; t++) {
        Constraint Pbound_UB("Pbound_UB" + to_string(t));
        Constraint Pbound_LB("Pbound_LB" + to_string(t));
        Constraint Qbound_UB("Qbound_UB" + to_string(t));
        Constraint Qbound_LB("Qbound_LB" + to_string(t));

        Pbound_UB = Pg.in(grid->gens, t) - PUB.in(grid->gens);
        Pbound_LB = Pg.in(grid->gens, t) - PLB.in(grid->gens);
        Qbound_UB = Qg.in(grid->gens, t) - QUB.in(grid->gens);
        Qbound_UB = Qg.in(grid->gens, t) - QLB.in(grid->gens);

        ACUC.add_constraint(Pbound_UB <= 0);
        ACUC.add_constraint(Pbound_LB >= 0);
        ACUC.add_constraint(Qbound_UB <= 0);
        ACUC.add_constraint(Qbound_LB >= 0);
    }

    solver OPF(ACUC,cplex);
    OPF.run();
    return 0;
}
