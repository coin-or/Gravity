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
#include <gravity/PowerNet.h>
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
    // ACUC
    PowerNet* grid = new PowerNet();
    const char* fname;
    fname = "../../data_sets/ACUC/nesta_case5_pjm.m";
    grid->readgrid(fname);

    // Time periods
    int T = 24;
    int nb_gen = grid->gens.size();
    int nb_lines = grid->arcs.size();
    int nb_bues = grid->nodes.size();

    /** build model */

    Model ACUC("AC-UC Model");

    /** define variables */

    // power generation
    var<double> Pg("Pg");
    var<double> Qg ("Qg");
    ACUC.add_var(Pg^(nb_gen*T));
    ACUC.add_var(Qg^(nb_gen*T));


    // power flow
    var<double> Pf("p_f");
    var<double> Qf("q_f");
    ACUC.add_var(Pf^(nb_lines*T));
    ACUC.add_var(Qf^(nb_lines*T));

    // voltage related variables.
    var<double> WR("WR");
    var<double> WI("WI");
    ACUC.add_var(WR^(T*nb_lines*(nb_lines- 1)/2));
    ACUC.add_var(WI^(T*nb_lines*(nb_lines- 1)/2));

    // Construct the objective function.
    func_ obj;
    for (int t =0; t < T; t++)
        for (int i = 0; i < nb_gen; i++) {
            auto g = grid->gens.at(i);
            if (!g->_active)
                continue;
            obj += g->_cost->c0 + (g->_cost->c1)*Pg(i,t) + g->_cost->c2*Pg(i,t);
        }

    /** Define constraints */

    // AC-KCL constraints
    //-- for generators
    /** subject to KCL_P {i in buses}: sum{(l,i,j) in arcs} p[l,i,j]  + load_p[i] = sum{(i,gen) in bus_gen} pg[gen];
     subject to KCL_Q {i in buses}: sum{(l,i,j) in arcs} q[l,i,j] + load_q[i] = sum{(i,gen) in bus_gen} qg[gen];
     */
    constant<int> ones(1);
    for (auto b: grid->nodes) {
        Bus* bus = (Bus*)b;
        Constraint KCL_P("KCL_P"+bus->_name);
        // for each t
        //for (auto a: b->get_out()){
        KCL_P += Pf.to(b->get_out())*ones;
        KCL_P += Pf.from(b->get_out())*ones;
        KCL_P += bus->pl();
        for (auto g: bus->_gen)
            KCL_P -= Pg(g->ID);
        KCL_P = 0;

        ACUC.add_constraint(KCL_P);

        Constraint KCL_Q("KCL_Q");
        KCL_Q += Qf.to(b->get_out())*ones;
        KCL_Q += Qf.from(b->get_out())*ones;
        KCL_Q += bus->ql();
        for (auto g: bus->_gen)
            KCL_P -= Qg(g->ID);
        KCL_Q = 0;
        ACUC.add_constraint(KCL_Q);
    }

    //AC-PF definitions.
    for (auto a: grid->arcs) {
        Line* la = (Line *) a;
        if (la->status==1) {
            Bus* src = (Bus*) la->src;
            Bus* dest = (Bus*) la->dest;
            Constraint Flow_P_From("Flow_P_From: "+ la->_name);
            /** subject to Flow_P_From {(l,i,j) in arcs_from}:
             p[l,i,j] = g[l]*(v[i]/tr[l])^2
             + -g[l]*v[i]/tr[l]*v[j]*cos(t[i]-t[j]-as[l])
             + -b[l]*v[i]/tr[l]*v[j]*sin(t[i]-t[j]-as[l]);
             */
            Flow_P_From += la->pi;
            Flow_P_From -= la->g*(src->_V_.square_magnitude())/pow(la->tr,2);
            Flow_P_From -= (-a->g*a->cc + a->b*a->dd)/(pow(a->cc,2)+pow(a->dd,2))*(src->vr*dest->vr + src->vi*dest->vi);
            Flow_P_From -= (-a->b*a->cc - a->g*a->dd)/(pow(a->cc,2)+pow(a->dd,2))*(src->vi*dest->vr - src->vr*dest->vi);
            Flow_P_From = 0;
            ACUC->add_constraint(Flow_P_From);
           

            Constraint Flow_P_To("Flow_P_To"+a->pj._name);
            Flow_P_To += a->pj;
            Flow_P_To -= a->g*(dest->_V_.square_magnitude());
            Flow_P_To -= (-a->g*a->cc - a->b*a->dd)/(pow(a->cc,2)+pow(a->dd,2))*(dest->vr*src->vr + dest->vi*src->vi);
            Flow_P_To -= (-a->b*a->cc + a->g*a->dd)/(pow(a->cc,2)+pow(a->dd,2))*(dest->vi*src->vr - dest->vr*src->vi);
            Flow_P_To = 0;
            _model->addConstraint(Flow_P_To);
            /** subject to Flow_Q_From {(l,i,j) in arcs_from}:
             q[l,i,j] = -(charge[l]/2+b[l])*(v[i]/tr[l])^2
             +  b[l]*v[i]/tr[l]*v[j]*cos(t[i]-t[j]-as[l])
             + -g[l]*v[i]/tr[l]*v[j]*sin(t[i]-t[j]-as[l]);
             */
            Constraint Flow_Q_From("Flow_Q_From"+a->qi._name);
            Flow_Q_From += a->qi;
            Flow_Q_From += (a->ch/2+a->b)*(src->_V_.square_magnitude())/pow(a->tr,2.);
            Flow_Q_From += (-a->b*a->cc - a->g*a->dd)/(pow(a->cc,2)+pow(a->dd,2))*(dest->vr*src->vr + dest->vi*src->vi);
            Flow_Q_From -= (-a->g*a->cc + a->b*a->dd)/(pow(a->cc,2)+pow(a->dd,2))*(src->vi*dest->vr - src->vr*dest->vi);
            Flow_Q_From = 0;
            ACUC->add_constraint(Flow_Q_From);
            /** subject to Flow_Q_To {(l,i,j) in arcs_to}:
             q[l,i,j] = -(charge[l]/2+b[l])*v[i]^2
             +  b[l]*v[i]*v[j]/tr[l]*cos(t[i]-t[j]+as[l])
             + -g[l]*v[i]*v[j]/tr[l]*sin(t[i]-t[j]+as[l]);
             */
            Constraint Flow_Q_To("Flow_Q_To"+a->qj._name);
            Flow_Q_To += a->qj;
            Flow_Q_To += (a->ch/2+a->b)*(dest->_V_.square_magnitude());
            Flow_Q_To += (-a->b*a->cc + a->g*a->dd)/(pow(a->cc,2)+pow(a->dd,2))*(dest->vr*src->vr + dest->vi*src->vi);
            Flow_Q_To -= (-a->g*a->cc - a->b*a->dd)/(pow(a->cc,2)+pow(a->dd,2))*(dest->vi*src->vr - dest->vr*src->vi);
            Flow_Q_To = 0;
            ACUC->add_constraint(Flow_Q_To);
        }
    }


    return 0;
}
