//
//  MinKpartition.cpp
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
        //  Handle error
        return 0;
    }
    if (!QueryPerformanceCounter(&time)) {
        //  Handle error
        return 0;
    }
    return (double)time.QuadPart / freq.QuadPart;
}
double get_cpu_time() {
    FILETIME a,b,c,d;
    if (GetProcessTimes(GetCurrentProcess(),&a,&b,&c,&d) != 0) {
        //  Returns total user time.
        //  Can be tweaked to include kernel times as well.
        return
            (double)(d.dwLowDateTime |
                     ((unsigned long long)d.dwHighDateTime << 32)) * 0.0000001;
    } else {
        //  Handle error
        return 0;
    }
}

//  Posix/Linux
#else
#include <time.h>
#include <sys/time.h>
double get_wall_time() {
    struct timeval time;
    if (gettimeofday(&time,NULL)) {
        //  Handle error
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
    var<double> Pf("Pf");
    var<double> Qf("Qf");
    cout << nb_lines << endl;
    ACOPF.add_var(Pf^(2*nb_lines));
    ACOPF.add_var(Qf^(2*nb_lines));

    // voltage related variables.
    // V_i = vr + jvi
    var<double> vr("vr");
    var<double> vi("vf");
    ACOPF.add_var(vr^(nb_buses));
    ACOPF.add_var(vi^(nb_buses));
    
    // Construct the objective function.
    func_ obj;
    for (auto g: grid->gens){
        if (!g->_active)
            continue;
        obj += g->_cost->c0 + (g->_cost->c1)*Pg(g->ID) + g->_cost->c2*power(Pg(g->ID), 2);
    }
    ACOPF.set_objective(min(obj));
    obj.print();
    
    
    /** Define constraints */
    // AC-KCL constraints
    //-- for each bus.
    /** subject to KCL_P {i in buses}: sum{(l,i,j) in arcs} p[l,i,j] + shunt_g[i]*v[i]^2 + load_p[i] = sum{(i,gen) in bus_gen} pg[gen];
     subject to KCL_Q {i in buses}: sum{(l,i,j) in arcs} q[l,i,j] - shunt_g[i]*v[i]^2 + load_q[i] = sum{(i,gen) in bus_gen} qg[gen];
     */
    constant<int> ones(1);
    for (auto b: grid->nodes) {
        Bus* bus = (Bus*)b;
        Constraint KCL_P("KCL_P"+bus->_name);
        Constraint KCL_Q("KCL_Q"+bus->_name);

        var<double> temp = Pf.in(b->get_out());
        cout << temp.get_dim() << endl;
        temp.print(true); /** bug */
//        for (auto a: b->get_out()){
//            KCL_P += Pf(a->id);
//            KCL_Q += Qf(a->id);
//        }
//        
//        for (auto a: b->get_in()){
//            KCL_P += Pf(a->id);
//            KCL_Q += Qf(a->id);
//        }
//        
//        KCL_P += bus->pl() + bus->gs()*(power(vr(bus->ID),2) + power(vi(bus->ID), 2));
//        KCL_Q += bus->ql()- bus->bs()*(power(vr(bus->ID),2)+power(vi(bus->ID), 2));
//
//        for (auto g: bus->_gen)
//            KCL_P -= Pg(g->ID);
//
//        ACOPF.add_constraint(KCL_P=0);
//        KCL_P.print();
//
//        for (auto g: bus->_gen)
//            KCL_Q -= Qg(g->ID);
//
//        ACOPF.add_constraint(KCL_Q=0);
//        KCL_Q.print();

    }
    //AC Power Flow definitions.
//    for (auto a: grid->arcs) {
//        Line* la = (Line *) a;
//        if (la->status == 1) {
//            Bus* src = (Bus*)(la->src);
//            Bus* dest = (Bus*)(la->dest);
//
//            Constraint Flow_P_From("Flow_P_From: "+ la->_name);
//            Flow_P_From += Pf(la->id);
//            Flow_P_From -= la->g*(power(vr(la->src->ID),2) + power(vi(la->src->ID),2))/pow(la->tr,2);
//            Flow_P_From -= (-la->g*la->cc + la->b*la->dd)/(pow(la->cc,2)+pow(la->dd,2))*(vr(src->ID)*vr(dest->ID) + vi(src->ID)*vi(dest->ID));
//            Flow_P_From -= (-la->b*la->cc - la->g*la->dd)/(pow(la->cc,2)+pow(la->dd,2))*(vr(src->ID)*vi(dest->ID) - vr(dest->ID)*vi(src->ID));
//            Flow_P_From = 0;
//            ACOPF.add_constraint(Flow_P_From);
//
//            Constraint Flow_P_To("Flow_P_To"+la->_name);
//            Flow_P_To += Pf(nb_lines + la->id);
//            Flow_P_To -= la->g*(power(vr(la->dest->ID),2) + power(vi(la->dest->ID),2))/pow(la->tr,2);
//            Flow_P_To -= (-la->g*la->cc - la->b*la->dd)/(pow(la->cc,2)+pow(la->dd,2))*(vr(dest->ID)*vr(src->ID) + vi(dest->ID)*vi(src->ID));
//            Flow_P_To -= (-la->b*la->cc + la->g*la->dd)/(pow(la->cc,2)+pow(la->dd,2))*(vi(dest->ID)*vr(src->ID) - vr(dest->ID)*vi(src->ID));
//            Flow_P_To = 0;
//            ACOPF.add_constraint(Flow_P_To);
//
//            Constraint Flow_Q_From("Flow_Q_From: "+ la->_name);
//            Flow_Q_From += Qf(la->id);
//            Flow_Q_From += (la->ch/2 + la->b)*(power(vr(la->src->ID),2) + power(vi(la->src->ID),2))/pow(la->tr,2);
//            Flow_Q_From += (-la->b*la->cc - la->g*la->dd)/(pow(la->cc,2)+pow(la->dd,2))*(vr(dest->ID)*vr(src->ID) + vi(dest->ID)*vi(src->ID));
//            Flow_Q_From -= (-la->g*la->cc + la->b*la->dd)/(pow(la->cc,2)+pow(la->dd,2))*(vi(src->ID)*vr(dest->ID) - vr(src->ID)*vi(dest->ID));
//            Flow_Q_From = 0;
//            ACOPF.add_constraint(Flow_Q_From);
//
//            Constraint Flow_Q_To("Flow_Q_To"+la->_name);
//            Flow_Q_From += Qf(nb_lines + la->id);
//            Flow_Q_To += (la->ch/2+la->b)*(power(vr(la->src->ID),2) + power(vi(la->src->ID),2))/pow(la->tr,2);
//            Flow_Q_To += (-la->b*la->cc + la->g*la->dd)/(pow(la->cc,2)+pow(la->dd,2))*(vr(dest->ID)*vr(src->ID) + vi(dest->ID)*vi(src->ID));
//            Flow_Q_To -= (-la->g*la->cc - la->b*la->dd)/(pow(la->cc,2)+pow(la->dd,2))*(vi(dest->ID)*vr(src->ID) - vr(dest->ID)*vi(src->ID));
//            Flow_Q_To = 0;
//            ACOPF.add_constraint(Flow_Q_To);
//        }
//    }
//    
//    // AC voltage limit constraints.
//    param<double> vbound_max_square("vbound_max_square");
//    param<double> vbound_min_square("vbound_min_square");
//
//    for (auto b: grid->nodes){
//        vbound_max_square.add_val(pow(((Bus*)b)->vbound.max, 2.));
//        vbound_min_square.add_val(pow(((Bus*)b)->vbound.min, 2.));
//    }
//    
//    Constraint Vol_limit_UB("Vol_limit_UB");
//    Vol_limit_UB = power(vr.in(grid->nodes), 2) + power(vi.in(grid->nodes),2);
//    Vol_limit_UB -= vbound_max_square.in(grid->nodes);
//    ACOPF.add_constraint(Vol_limit_UB <= 0);
//    
//    Constraint Vol_limit_LB("Vol_limit_LB");
//    Vol_limit_LB = power(vr.in(grid->nodes), 2) + power(vi.in(grid->nodes),2);
//    Vol_limit_LB -= vbound_max_square.in(grid->nodes);
//    ACOPF.add_constraint(Vol_limit_LB >= 0);
//
//     //AC-PAD constraints
//    param<double> tbound_max_tan("vbound_max_tan");
//    param<double> tbound_min_tan("vbound_min_tan");
//    
//    for (auto a: grid->arcs){
//        tbound_max_tan.add_val(tan(((Line*)a)->tbound.max));
//        tbound_min_tan.add_val(tan(((Line*)a)->tbound.min));
//    }
//    
//    Constraint PAD_UB("PAD_UB");
//    PAD_UB = vr.from(grid->arcs)*vi.to(grid->arcs) + vr.to(grid->arcs)*vi.from(grid->arcs);
//    PAD_UB -= tbound_max_tan.in(grid->arcs)*(vr.from(grid->arcs)*vr.to(grid->arcs) - vi.to(grid->arcs)*vi.from(grid->arcs));
//    ACOPF.add_constraint(PAD_UB <= 0);
//
//    
//    Constraint PAD_LB("PAD_LB:");
//    PAD_LB = vr.from(grid->arcs)*vi.to(grid->arcs) + vr.to(grid->arcs)*vi.from(grid->arcs);
//    PAD_LB -= tbound_min_tan.in(grid->arcs)*(vr.from(grid->arcs)*vr.to(grid->arcs) - vi.to(grid->arcs)*vi.from(grid->arcs));
//    ACOPF.add_constraint(PAD_LB >= 0);
//
//    
//    // Thermal_Limit {(l,i,j) in arcs}: p[l,i,j]^2 + q[l,i,j]^2 <= s[l]^2;*/
//    param<double> Thermal_limit_square("Thermal_limit_square");
//    for (auto a: grid->arcs){
//        Line* la = (Line *) a;
//        Thermal_limit_square.add_val(pow(la->limit, 2));
//    }
//    
//    Constraint Thermal_Limit_from("Thermal_Limit_from");
//    Thermal_Limit_from += power(Pf.in(grid->arcs), 2) + power(Qf.in(grid->arcs),2);
//    Thermal_Limit_from -= Thermal_limit_square;
//    ACOPF.add_constraint(Thermal_Limit_from <= 0);
//    
//    for (auto a: grid->arcs) {
//        Line* la = (Line *) a;
//        if (la->status == 1){
//            Constraint Thermal_Limit_to("Thermal_Limit_to"+ to_string(la->id));
//            Thermal_Limit_to += power(Pf(la->id + nb_lines), 2) + power(Qf(la->id + nb_lines),2);
//            Thermal_Limit_to <= pow(la->limit, 2);
//            ACOPF.add_constraint(Thermal_Limit_to);
//        }
//    }
//    
//    // Power generation constraints.
//    for (auto g: grid->gens){
//        Constraint Pbound_UB("Pbound_UB" + g->_name);
//        Pbound_UB += Pg(g->ID);
//        ACOPF.add_constraint(Pbound_UB <= g->pbound.max);
//        
//        Constraint Pbound_LB("Pbound_LB" + g->_name);
//        Pbound_LB += Pg(g->ID);
//        ACOPF.add_constraint(Pbound_LB >= g->pbound.min);
//        
//        
//        Constraint Qbound_UB("Qbound_UB" + g->_name);
//        Qbound_UB += Qg(g->ID);
//        ACOPF.add_constraint(Qbound_UB <= g->qbound.max);
//        
//        Constraint Qbound_LB("Qbound_LB" + g->_name);
//        Qbound_LB += Qg(g->ID);
//        ACOPF.add_constraint(Qbound_LB >= g->qbound.min);
//    }
    
    solver OPF(ACOPF,ipopt);
    OPF.run();
    return 0;
}
