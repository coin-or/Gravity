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
    var<double> Pf("Pf");
    var<double> Qf("Qf");
    cout << nb_lines << endl;
    ACOPF.add_var(Pf^(2*nb_lines));
    ACOPF.add_var(Qf^(2*nb_lines));

    // voltage related variables.
    var<double> vr("vr");
    var<double> vi("vf");
    ACOPF.add_var(vr^(nb_buses));
    ACOPF.add_var(vi^(nb_buses));
    
    /** Construct the objective function*/
    func_ obj;
    for (auto g: grid->gens){
        if (!g->_active)
            continue;
        obj += g->_cost->c0 + (g->_cost->c1)*Pg(g->ID) + g->_cost->c2*power(Pg(g->ID), 2);
    }
    ACOPF.set_objective(min(obj));
    obj.print();
    
    
    /** Define constraints */
    constant<int> ones(1);
    
    //KCL
    for (auto b: grid->nodes) {
        Bus* bus = (Bus*)b;
        Constraint KCL_P("KCL_P"+bus->_name);
        Constraint KCL_Q("KCL_Q"+bus->_name);
        var<double> temp = Pf.in(b->get_out());
        KCL_P = ones.tr()*(temp);
        KCL_Q = ones.tr()*Qf.in(b->get_out());
        KCL_P += ones.tr()*Pf.in(b->get_in());
        KCL_Q += ones.tr()*Qf.in(b->get_in());
        
        KCL_P += bus->pl() + bus->gs()*(power(vr(bus->ID),2) + power(vi(bus->ID), 2));
        KCL_Q += bus->ql()- bus->bs()*(power(vr(bus->ID),2)+power(vi(bus->ID), 2));
        
        for (auto g: bus->_gen) {
            KCL_P -= Pg(g->ID);
            KCL_Q -= Qg(g->ID);
        }
        
        ACOPF.add_constraint(KCL_P=0);
        ACOPF.add_constraint(KCL_Q=0);
        break;
    }

    //AC Power Flow.
    param<double> gij("gij");
    param<double> bij("bij");

    for (auto a: grid->arcs) {
        gij.add_val(((Line *)a)->g);
        bij.add_val(((Line *)a)->g);
    }

    Constraint Flow_P_From("Flow_P_From");
    Flow_P_From += Pf.in(grid->arcs);
    Flow_P_From -= gij.in(grid->arcs)*(power(vr.from(grid->arcs),2) + power(vi.from(grid->arcs),2));
    Flow_P_From += gij.in(grid->arcs)*(vr.from(grid->arcs)*vr.to(grid->arcs) + vi.from(grid->arcs)*vi.to(grid->arcs));
    Flow_P_From -= bij.in(grid->arcs)*(vr.from(grid->arcs)*vi.to(grid->arcs) - vr.to(grid->arcs)*vi.from(grid->arcs));
    Flow_P_From = 0;
    ACOPF.add_constraint(Flow_P_From);
    
    //Constraint Flow_P_To("Flow_P_To");
    //Flow_P_To += Pf.in(grid->arcs);
    //Flow_P_To -= gij.in(grid->arcs)*(power(vr.to(grid->arcs),2) + power(vi.to(grid->arcs),2));
    //Flow_P_To += gij.in(grid->arcs)*(vr.from(grid->arcs)*vr.to(grid->arcs) + vi.from(grid->arcs)*vi.to(grid->arcs));
    //Flow_P_To -= bij.in(grid->arcs)*(vr.from(grid->arcs)*vi.to(grid->arcs) - vr.to(grid->arcs)*vi.from(grid->arcs));
    //Flow_P_To = 0;
    //ACOPF.add_constraint(Flow_P_To);

    Constraint Flow_Q_From("Flow_Q_From");
    Flow_Q_From += Qf.in(grid->arcs);
    Flow_Q_From -= bij.in(grid->arcs)*(power(vr.from(grid->arcs),2) + power(vi.from(grid->arcs),2));
    Flow_Q_From -= bij.in(grid->arcs)*(vr.from(grid->arcs)*vr.to(grid->arcs) + vi.from(grid->arcs)*vi.to(grid->arcs));
    Flow_Q_From -= gij.in(grid->arcs)*(vr.from(grid->arcs)*vi.to(grid->arcs) - vr.to(grid->arcs)*vi.from(grid->arcs));
    Flow_Q_From = 0;
    ACOPF.add_constraint(Flow_Q_From);
    
    //Constraint Flow_Q_To("Flow_Q_To");
    //Flow_Q_To += Qf.in(grid->arcs);
    //Flow_Q_To -= bij.in(grid->arcs)*(power(vr.to(grid->arcs),2) + power(vi.to(grid->arcs),2));
    //Flow_Q_To -= bij.in(grid->arcs)*(vr.from(grid->arcs)*vr.to(grid->arcs) + vi.from(grid->arcs)*vi.to(grid->arcs));
    //Flow_Q_To -= gij.in(grid->arcs)*(vr.from(grid->arcs)*vi.to(grid->arcs) - vr.to(grid->arcs)*vi.from(grid->arcs));
    //Flow_Q_To = 0;
    //ACOPF.add_constraint(Flow_Q_To);

//    // AC voltage limit constraints.
    param<double> vbound_max_square("vbound_max_square");
    param<double> vbound_min_square("vbound_min_square");

    for (auto b: grid->nodes){
        vbound_max_square.add_val(pow(((Bus*)b)->vbound.max, 2.));
        vbound_min_square.add_val(pow(((Bus*)b)->vbound.min, 2.));
    }
    
    Constraint Vol_limit_UB("Vol_limit_UB");
    Vol_limit_UB = power(vr.in(grid->nodes), 2) + power(vi.in(grid->nodes),2);
    Vol_limit_UB -= vbound_max_square.in(grid->nodes);
    ACOPF.add_constraint(Vol_limit_UB <= 0);
    
    Constraint Vol_limit_LB("Vol_limit_LB");
    Vol_limit_LB = power(vr.in(grid->nodes), 2) + power(vi.in(grid->nodes),2);
    Vol_limit_LB -= vbound_max_square.in(grid->nodes);
    ACOPF.add_constraint(Vol_limit_LB >= 0);

     //AC-PAD constraints
    param<double> tbound_max_tan("vbound_max_tan");
    param<double> tbound_min_tan("vbound_min_tan");
    
    for (auto a: grid->arcs){
        tbound_max_tan.add_val(tan(((Line*)a)->tbound.max));
        tbound_min_tan.add_val(tan(((Line*)a)->tbound.min));
    }
    
    Constraint PAD_UB("PAD_UB");
    PAD_UB = vr.from(grid->arcs)*vi.to(grid->arcs) + vr.to(grid->arcs)*vi.from(grid->arcs);
    PAD_UB -= tbound_max_tan.in(grid->arcs)*(vr.from(grid->arcs)*vr.to(grid->arcs) - vi.to(grid->arcs)*vi.from(grid->arcs));
    ACOPF.add_constraint(PAD_UB <= 0);

    
    Constraint PAD_LB("PAD_LB:");
    PAD_LB = vr.from(grid->arcs)*vi.to(grid->arcs) + vr.to(grid->arcs)*vi.from(grid->arcs);
    PAD_LB -= tbound_min_tan.in(grid->arcs)*(vr.from(grid->arcs)*vr.to(grid->arcs) - vi.to(grid->arcs)*vi.from(grid->arcs));
    ACOPF.add_constraint(PAD_LB >= 0);

    
    // Thermal_Limit {(l,i,j) in arcs}: p[l,i,j]^2 + q[l,i,j]^2 <= s[l]^2;*/
    param<double> Thermal_limit_square("Thermal_limit_square");
    for (auto a: grid->arcs){
        Line* la = (Line *) a;
        Thermal_limit_square.add_val(pow(la->limit, 2));
    }
    
    Constraint Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from += power(Pf.in(grid->arcs), 2) + power(Qf.in(grid->arcs),2);
    Thermal_Limit_from -= Thermal_limit_square;
    ACOPF.add_constraint(Thermal_Limit_from <= 0);

    // Power generation constraints.
    param<double> PUB("PUB");
    param<double> PLB("PUB");
    param<double> QUB("QUB");
    param<double> QLB("QUB");
    
    for (auto g: grid->gens){
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
