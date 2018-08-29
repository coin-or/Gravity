//
//  PowerNet.cpp
//
//
//  Created by Guanglei Wang on 03/06/2017.
//

#include "PowerNet.h"
#include <algorithm>
#include <map>
#define _USE_MATH_DEFINES
#include <cmath>
#include <list>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <queue>
#include <time.h>
//#include <xlnt/xlnt.hpp>
#include <gravity/solver.h>
//#define USEDEBUG
#ifdef USEDEBUG
#define DebugOn(x) cout << x
#else
#define Debug(x)
#endif
#define DebugOff(x)

using namespace std;

PowerNet::PowerNet() {
    bMVA = 0;
    pg_min.set_name("pg_min");
    pg_max.set_name("pg_max");
    qg_min.set_name("qg_min");
    qg_max.set_name("qg_max");
    pb_min.set_name("pb_min");
    pb_max.set_name("pb_max");
    qb_min.set_name("qb_min");
    qb_max.set_name("qb_max");
    pv_min.set_name("pv_min");
    pv_max.set_name("pv_max");
    qv_min.set_name("qv_min");
    qv_max.set_name("qv_max");
    pw_min.set_name("pw_min");
    pw_max.set_name("pw_max");
    qw_min.set_name("qw_min");
    qw_max.set_name("qw_max");
    pv_out.set_name("pv_out");
    pv_capcost.set_name("pv_capcost");
    pv_varcost.set_name("pv_varcost");
    pg_s.set_name("pg_s");
    qg_s.set_name("qg_s");
    cb_f.set_name("cb_f");
    cb_v.set_name("cb_v");
    c0.set_name("c0");
    c1.set_name("c1");
    c2.set_name("c2");
    p_factor.set_name("p_factor");
    ramp_up.set_name("ramp_up");
    ramp_down.set_name("ramp_down");
    gen_eff.set_name("gen_eff");
    min_ut.set_name("min_ut");
    min_dt.set_name("min_dt");
    min_diesel_invest.set_name("min_diesel_invest");
    max_diesel_invest.set_name("max_diesel_invest");
    min_batt_invest.set_name("min_batt_invest");
    max_batt_invest.set_name("max_batt_invest");
    gen_capcost.set_name("cg");
    expansion_capcost.set_name("ce");
    inverter_capcost.set_name("cb");
    th_min.set_name("th_min");
    th_max.set_name("th_max");
    cphi.set_name("cphi");
    sphi.set_name("sphi");
    cos_d.set_name("cos_d");
    tan_th_min.set_name("tan_th_min");
    tan_th_max.set_name("tan_th_max");
    v_diff_max.set_name("v_diff_max");
    v_min.set_name("v_min");
    v_max.set_name("v_max");
    w_min.set_name("w_min");
    w_max.set_name("w_max");
    wr_min.set_name("wr_min");
    wr_max.set_name("wr_max");
    wi_min.set_name("wi_min");
    wi_max.set_name("wi_max");
    v_s.set_name("v_s");
    pl.set_name("pl");
    pl_ratio.set_name("pl_ratio");
    ql.set_name("ql");
    gs.set_name("gs");
    bs.set_name("bs");
    g.set_name("g");
    b.set_name("b");
    r.set_name("r");
    x.set_name("x");
    ch.set_name("ch");
    as.set_name("as");
    tr.set_name("tr");
    S_max.set_name("S_max");
    eff_a.set_name("eff_a");
    eff_b.set_name("eff_b");
    g_ff.set_name("g_ff");
    g_ft.set_name("g_ft");
    g_tf.set_name("g_tf");
    g_tt.set_name("g_tt");

    b_ff.set_name("b_ff");
    b_ft.set_name("b_ft");
    b_tf.set_name("b_tf");
    b_tt.set_name("b_tt");
    Y.set_name("Y");
    Y_t.set_name("Y_t");
    Y_charge.set_name("Y_charge");
    Y_charge_t.set_name("Y_charge_t");
}

PowerNet::~PowerNet() {
    if(!gens.empty()) {
        for (Gen* g:gens) {
            delete g;
        }
        gens.clear();
    }
    for (Node* n:nodes) {
        delete n;
    }
    nodes.clear();
    for (Arc* a:arcs) {
        delete a;
    }
    arcs.clear();
}
//
// Read a grid
// @discussion line with delimiter ";"

string PowerNet::get_ref_bus() {
    return ref_bus;
}

unsigned PowerNet::get_nb_active_gens() const {
    unsigned nb=0;
    for (auto g: gens) {
        if (g->_active) {
            nb++;
        }
    }
    return nb;
}

unsigned PowerNet::get_nb_active_bus_pairs() const {
    unsigned nb=0;
    for (auto bp: _bus_pairs._keys) {
//        if (bp->_active) {
            nb++;
//        }
    }
    return nb;
}


unsigned PowerNet::get_nb_active_arcs() const {
    unsigned nb=0;
    for (auto a: arcs) {
        if (a->_active) {
            nb++;
        }
    }
    return nb;
}

unsigned PowerNet::get_nb_active_nodes() const {
    unsigned nb=0;
    for (auto n: nodes) {
        if (n->_active) {
            nb++;
        }
        else {
            DebugOff("Inactive Node" << n->_name << endl);
        }
    }
    return nb;
}

void PowerNet::time_expand(const indices& T) {
//    c0.time_expand(T);
//    c1.time_expand(T);
//    c2.time_expand(T);
    pl._time_extended = true;
    pl_ratio._time_extended = true;
    ql._time_extended = true;
    pv_out._time_extended = true;
    pw_min._time_extended = true;
    pw_max._time_extended = true;
//    S_max.time_expand(T);
//    th_min.time_expand(T);
//    th_max.time_expand(T);
//    tan_th_min.time_expand(T);
//    tan_th_max.time_expand(T);
//    r.time_expand(T);
//    x.time_expand(T);
//    pg_min.time_expand(T);
//    pg_max.time_expand(T);
//    pv_min.time_expand(T);
//    pv_max.time_expand(T);
//    qg_min.time_expand(T);
//    qg_max.time_expand(T);
//    pb_min.time_expand(T);
//    pb_max.time_expand(T);
//    qb_min.time_expand(T);
//    qb_max.time_expand(T);
//    w_min.time_expand(T);
//    w_max.time_expand(T);
//    v_min.time_expand(T);
//    v_max.time_expand(T);
//    eff_a.time_expand(T);
//    eff_b.time_expand(T);
}

void PowerNet::time_expand(unsigned T) {
    c0.time_expand(T);
    c1.time_expand(T);
    c2.time_expand(T);
    S_max.time_expand(T);
    th_min.time_expand(T);
    th_max.time_expand(T);
    tan_th_min.time_expand(T);
    tan_th_max.time_expand(T);
    g_tt.time_expand(T);
    g_ff.time_expand(T);
    g_ft.time_expand(T);
    g_tf.time_expand(T);
    b_tt.time_expand(T);
    b_ff.time_expand(T);
    b_ft.time_expand(T);
    b_tf.time_expand(T);
    pg_min.time_expand(T);
    pg_max.time_expand(T);
    qg_min.time_expand(T);
    qg_max.time_expand(T);
    w_min.time_expand(T);
    w_max.time_expand(T);
    gs.time_expand(T);
    bs.time_expand(T);
    pl.time_expand(T);
    ql.time_expand(T);
    v_min.time_expand(T);
    v_max.time_expand(T);
}

int PowerNet::readGAMS(const string& fname) {
    double pi = 4.*atan(1.);
    string name;
    double kvb = 0;
    //    int id = 0;
    unsigned index = 0;
    cout << "Loading file " << fname << endl;
    ifstream file(fname, std::ifstream::in);
    if(!file.is_open()) {
        auto new_fname = "../" + string(fname);
        file.open(new_fname);
        if(!file.is_open()) {
            throw invalid_argument("Could not open file\n");
        }
    }
    string word;
    while (word.compare("BaseMVA")) {
        file >> word;
    }
    getline(file, word);
    getline(file, word,'/');
    bMVA = atoi(word.c_str());
    DebugOn("bMVA = " << bMVA << endl);
    //  cout << _name << endl;
    while (word.compare("branches(*)")) {
        file >> word;
    }
    
    getline(file, word);
    file >> word;
    while (word.compare("/;")) {
        auto arc = new Line(word);
        _exist_arcs.push_back(arc);
        arcMap[word] = arc;
        file >> word;
    }
    DebugOn("Number of lines = " << _exist_arcs.size() << endl);
    
    while (word.compare("buses(*)")) {
        file >> word;
    }
    getline(file, word);
    file >> word;
    while (word.compare("/;")) {
        auto bus = new Bus(word);
        add_node(bus);
        file >> word;
    }
    ref_bus = nodes[0]->_name;
    DebugOn("Number of buses = " << nodes.size() << endl);
    
    while (word.compare("contingencies(*)")) {
        file >> word;
    }
    getline(file, word);
    file >> word;
    while (word.compare("/;")) {
        auto pos = word.find('*');
        if (pos!=string::npos) {
            int loperand= atoi(word.substr(0,pos).c_str());
            int roperand= atoi(word.substr(pos+1,word.size()-1).c_str());
            _nb_conting = loperand*roperand;
        }
        else {
            _nb_conting  = atoi(word.c_str());
        }
        file >> word;
    }
    DebugOn("Number of contingencies = " << _nb_conting << endl);
    v_cont.resize(_nb_conting);
    theta_cont.resize(_nb_conting);
    Qg_cont.resize(_nb_conting);
    v_diff_n_cont.resize(_nb_conting);
    v_diff_p_cont.resize(_nb_conting);
    _delta.resize(_nb_conting);
    p_from.resize(_nb_conting);
    q_from.resize(_nb_conting);
    p_to.resize(_nb_conting);
    q_to.resize(_nb_conting);
    
    /* Generator data */
    while (word.compare("generators(*)")) {
        file >> word;
    }
    getline(file, word);
    file >> word;
    while (word.compare("/;")) {
        auto gen = new Gen(word);
        gens.push_back(gen);
        genMap[word] = gen;
        file >> word;
    }
    DebugOn("Number of generators = " << gens.size() << endl);
    
    /* Branch origin/dest */
    while (word.compare("ijOrigin(i,j)")) {
        file >> word;
    }
    getline(file, word);
    file >> word;
    unsigned idx = 0;
    while (word.compare("/;")) {
        word =word.substr(word.find('.')+1, word.size()-1);
        auto arc = _exist_arcs[idx++];
        arc->_src = get_node(word);
        file >> word;
    }
    while (word.compare("ijDestination(i,j)")) {
        file >> word;
    }
    getline(file, word);
    file >> word;
    idx = 0;
    while (word.compare("/;")) {
        word =word.substr(word.find('.')+1, word.size()-1);
        auto arc = _exist_arcs[idx++];
        arc->_dest = get_node(word);
        file >> word;
    }
    while (word.compare("icMap(i,c)")) {
        file >> word;
    }
    getline(file, word);
    file >> word;
    idx = 0;
    while (word.compare("/;")) {
        auto pos = word.find('.');
        word =word.substr(pos+1, word.size()-1);
        auto arc = _exist_arcs[idx++];
        if (word.find("\"")!=std::string::npos) {
            word = word.substr(1,word.size()-2);
        }
        arc->_circuit_id = word;
        file >> word;
    }
    
    
    while (word.compare("ljMap(l,j)")) {
        file >> word;
    }
    getline(file, word);
    file >> word;
    idx = 0;
    while (word.compare("/;")) {
        word =word.substr(word.find('.')+1, word.size()-1);
        auto gen = gens[idx++];
        gen->_bus = (Bus*)get_node(word);
        gen->_bus->_gen.push_back(gen);
        file >> word;
    }
    
    while (word.compare("luMap(l,u)")) {
        file >> word;
    }
    getline(file, word);
    file >> word;
    idx = 0;
    while (word.compare("/;")) {
        word =word.substr(word.find(".")+1, word.size()-1);
        if (word.find("\"")!=std::string::npos) {
            word = word.substr(1,word.size()-1);
        }
        if (word.find("'")!=std::string::npos) {
            word += " '";
        }
        auto gen = gens[idx++];
        gen->_unit_name = word;
        getline(file, word);
        file >> word;
    }
    
    while (word.compare("ikInactive(i,k)")) {
        file >> word;
    }
    getline(file, word);
    file >> word;
    while (word.compare("/;")) {
        word = word.substr(0, word.find('.'));
        auto arc = arcMap.at(word);
        _conting_lines.push_back(arc);        
        file >> word;
    }
    DebugOn("Number of line contingencies = " << _conting_lines.size() << endl);
    
    while (word.compare("lkInactive(l,k)")) {
        file >> word;
    }
    getline(file, word);
    file >> word;
    while (word.compare("/;")) {
        word = word.substr(0, word.find('.'));
        auto gen = genMap.at(word);
        _conting_gens.push_back(gen);
        file >> word;
    }
    DebugOn("Number of generator contingencies = " << _conting_gens.size() << endl);
    
    while (word.compare("jShuntConductance(j)")) {
        file >> word;
    }
    getline(file, word);
    file >> word;
    idx = 0;
    while (word.compare("/;")) {
        file >> word;
        auto bus = (Bus*)nodes[idx++];
        bus->_cond.at(0)->_gs = atof(word.c_str());
        file >> word;
    }
    while (word.compare("jShuntSusceptance(j)")) {
        file >> word;
    }
    getline(file, word);
    file >> word;
    idx = 0;
    while (word.compare("/;")) {
        file >> word;
        auto bus = (Bus*)nodes[idx++];
        bus->_cond.at(0)->_bs = atof(word.c_str());
        file >> word;
    }
    while (word.compare("jVoltageMagnitudeMin(j)")) {
        file >> word;
    }
    getline(file, word);
    file >> word;
    idx = 0;
    while (word.compare("/;")) {
        file >> word;
        auto bus = (Bus*)nodes[idx++];
        bus->vbound.min = atof(word.c_str());
        file >> word;
    }
    while (word.compare("jVoltageMagnitudeMax(j)")) {
        file >> word;
    }
    getline(file, word);
    file >> word;
    idx = 0;
    while (word.compare("/;")) {
        file >> word;
        auto bus = (Bus*)nodes[idx++];
        bus->vbound.max = atof(word.c_str());
        file >> word;
    }
    while (word.compare("jRealPowerDemand(j)")) {
        file >> word;
    }
    getline(file, word);
    file >> word;
    idx = 0;
    while (word.compare("/;")) {
        file >> word;
        auto bus = (Bus*)nodes[idx++];
        bus->_cond.at(0)->_pl = atof(word.c_str());
        file >> word;
    }
    while (word.compare("jReactivePowerDemand(j)")) {
        file >> word;
    }
    getline(file, word);
    file >> word;
    idx = 0;
    while (word.compare("/;")) {
        file >> word;
        auto bus = (Bus*)nodes[idx++];
        bus->_cond.at(0)->_ql = atof(word.c_str());
        file >> word;
    }
    while (word.compare("iSeriesResistance(i)")) {
        file >> word;
    }
    getline(file, word);
    file >> word;
    idx = 0;
    while (word.compare("/;")) {
        file >> word;
        auto arc = (Line*)_exist_arcs[idx++];
        arc->r = atof(word.c_str());
        file >> word;
    }
    while (word.compare("iSeriesReactance(i)")) {
        file >> word;
    }
    getline(file, word);
    file >> word;
    idx = 0;
    while (word.compare("/;")) {
        file >> word;
        auto arc = (Line*)_exist_arcs[idx++];
        arc->x = atof(word.c_str());
        auto res = pow(arc->r,2) + pow(arc->x,2);
        
        if (res==0) {
            cerr << " line with r = x = 0" << endl;
            exit(-1);
        }
        // define g and b for each conductor.
        arc->g = arc->r/res;
        arc->b = -arc->x/res;
        file >> word;
    }
    while (word.compare("iChargingSusceptance(i)")) {
        file >> word;
    }
    getline(file, word);
    file >> word;
    idx = 0;
    while (word.compare("/;")) {
        file >> word;
        auto arc = (Line*)_exist_arcs[idx++];
        arc->ch = atof(word.c_str());
        file >> word;
    }
    while (word.compare("iTapRatio(i)")) {
        file >> word;
    }
    getline(file, word);
    file >> word;
    idx = 0;
    while (word.compare("/;")) {
        file >> word;
        auto arc = (Line*)_exist_arcs[idx++];
        arc->tr = atof(word.c_str());
        file >> word;
    }
    while (word.compare("iPhaseShift(i)")) {
        file >> word;
    }
    getline(file, word);
    file >> word;    
    idx = 0;
    while (word.compare("/;")) {
        file >> word;
        auto arc = (Line*)_exist_arcs[idx++];
        arc->as = atof(word.c_str());
        file >> word;
    }
    while (word.compare("iPowerMagnitudeMax(i)")) {
        file >> word;
    }
    getline(file, word);
    file >> word;
    idx = 0;
    while (word.compare("/;")) {
        file >> word;
        auto arc = (Line*)_exist_arcs[idx++];
        arc->smax = atof(word.c_str());
        arc->limit = arc->smax/bMVA;
        file >> word;
    }
    while (word.compare("lRealPowerMin(l)")) {
        file >> word;
    }
    getline(file, word);
    file >> word;
    idx = 0;
    while (word.compare("/;")) {
        file >> word;
        auto gen = gens[idx++];
        gen->_pbound.min = atof(word.c_str())/bMVA;
        file >> word;
    }
    while (word.compare("lRealPowerMax(l)")) {
        file >> word;
    }
    getline(file, word);
    file >> word;
    idx = 0;
    while (word.compare("/;")) {
        file >> word;
        auto gen = gens[idx++];
        gen->_pbound.max = atof(word.c_str())/bMVA;
        file >> word;
    }
    while (word.compare("lReactivePowerMin(l)")) {
        file >> word;
    }
    getline(file, word);
    file >> word;
    idx = 0;
    while (word.compare("/;")) {
        file >> word;
        auto gen = gens[idx++];
        gen->_qbound.min = atof(word.c_str())/bMVA;
        file >> word;
    }
    while (word.compare("lReactivePowerMax(l)")) {
        file >> word;
    }
    getline(file, word);
    file >> word;
    idx = 0;
    while (word.compare("/;")) {
        file >> word;
        auto gen = gens[idx++];
        gen->_qbound.max = atof(word.c_str())/bMVA;
        file >> word;
    }
    while (word.compare("lmRealPowerCostCoefficient(l,m)")) {
        file >> word;
    }
    getline(file, word);
    file >> word;
    idx = 0;
    while (word.compare("/;")) {
        file >> word;
        auto gen = gens[idx++];
        gen->_cost->c0 = atof(word.c_str());
        file >> word;
        file >> word;
        gen->_cost->c1 = atof(word.c_str());
        file >> word;
        file >> word;
        gen->_cost->c2 = atof(word.c_str());
        file >> word;
    }
    while (word.compare("lParticipationFactor(l)")) {
        file >> word;
    }
    getline(file, word);
    file >> word;
    idx = 0;
    double tot_part_fact = 0;
    while (word.compare("/;")) {
        file >> word;
        auto gen = gens[idx++];
        gen->_part_factor = atof(word.c_str());
        tot_part_fact +=gen->_part_factor;
        file >> word;
    }
    
    DebugOn("Bus Data : "<<endl);
    for (auto n: nodes) {
        auto b = (Bus*)n;
        b->print();
        name = b->_name;
        pl.set_val(name,b->pl()/bMVA);
        ql.set_val(name,b->ql()/bMVA);
        gs.set_val(name,b->gs()/bMVA);
        bs.set_val(name,b->bs()/bMVA);
        v_max.set_val(name,b->vmax());
        v_min.set_val(name,b->vmin());
        w_min.set_val(name,pow(b->vmin(), 2));
        w_max.set_val(name,pow(b->vmax(), 2));
    }
    DebugOn("Gen Data : "<<endl);
    for (auto gen:gens) {
        gen->print();
        name = gen->_name;
        qg_max.set_val(name,gen->_qbound.max);
        qg_min.set_val(name,gen->_qbound.min);
        pg_max.set_val(name,gen->_pbound.max);
        pg_min.set_val(name,gen->_pbound.min);
        c0.set_val(name, gen->_cost->c0);
        c1.set_val(name, gen->_cost->c1*bMVA);
        c2.set_val(name, gen->_cost->c2*pow(bMVA,2));
        p_factor.set_val(name, gen->_part_factor/tot_part_fact);
//        p_factor.set_val(name, gen->_part_factor);
    }
    
    /* Lines data */
    Bus* bus_s = nullptr;
    Bus* bus_d = nullptr;
    string src, dest;
    set<string> bus_pair_names;
    DebugOn("Line Data : "<<endl);
    for(auto a: _exist_arcs) {
        auto arc = (Line*)a;
        bus_s = (Bus*)arc->_src;
        bus_d = (Bus*)arc->_dest;
        src = arc->_src->_name;
        dest = arc->_dest->_name;
        arc->_name += "," + src + "," + dest;
        name = arc->_name;
        g.set_val(name,arc->g);
        b.set_val(name,arc->b);
        tr.set_val(name,arc->tr);
        as.set_val(name,arc->as);
        arc->cc = arc->tr*cos(arc->as); // Rectangular values for transformer phase shifters
        arc->dd = arc->tr*sin(arc->as);
        g_ff.set_val(name,arc->g/(pow(arc->cc, 2) + pow(arc->dd, 2)));
        g_ft.set_val(name,(-arc->g*arc->cc + arc->b*arc->dd)/(pow(arc->cc, 2) + pow(arc->dd, 2)));
        g_tt.set_val(name,arc->g);
        g_tf.set_val(name,(-arc->g*arc->cc - arc->b*arc->dd)/(pow(arc->cc, 2) + pow(arc->dd, 2)));
        b_ff.set_val(name,(arc->ch*0.5 + arc->b)/(pow(arc->cc, 2) + pow(arc->dd, 2)));
        b_ft.set_val(name,(-arc->b*arc->cc - arc->g*arc->dd)/(pow(arc->cc, 2) + pow(arc->dd, 2)));
        b_tt.set_val(name,(arc->ch*0.5 + arc->b));
        b_tf.set_val(name,(-arc->b*arc->cc + arc->g*arc->dd)/(pow(arc->cc, 2) + pow(arc->dd, 2)));
        ch.set_val(name,arc->ch);
        S_max.set_val(name,arc->limit);
        arc->connect();
        add_arc(arc);
        arc->tbound.min = -30*pi/180;
        arc->tbound.max = 30*pi/180;
        m_theta_lb += arc->tbound.min;
        m_theta_ub += arc->tbound.max;
        arc->print();
        /* Switching to bus_pairs keys */
        name = src + "," + dest;
        if (arc->_active && bus_pair_names.count(name)==0) {
            _bus_pairs._keys.push_back(new index_pair(index_(src), index_(dest), arc->_active));
            bus_pair_names.insert(name);
        }
        if (!arc->_parallel) {
            th_min.set_val(name,arc->tbound.min);
            th_max.set_val(name,arc->tbound.max);
            tan_th_min.set_val(name,tan(arc->tbound.min));
            tan_th_max.set_val(name,tan(arc->tbound.max));
            
        }
        else {
            th_min.set_val(name,max(th_min.eval(name), arc->tbound.min));
            th_max.set_val(name,min(th_max.eval(name), arc->tbound.max));
            tan_th_min.set_val(name,tan(th_min.eval(name)));
            tan_th_max.set_val(name,tan(th_max.eval(name)));
        }
        if (arc->tbound.min >= 0) {
            wr_max.set_val(name,bus_s->vbound.max*bus_d->vbound.max*cos(th_min(name).eval()));
            wr_min.set_val(name,bus_s->vbound.min*bus_d->vbound.min*cos(th_max(name).eval()));
            wi_max.set_val(name,bus_s->vbound.max*bus_d->vbound.max*sin(th_max(name).eval()));
            wi_min.set_val(name,bus_s->vbound.min*bus_d->vbound.min*sin(th_min(name).eval()));
        };
        if (arc->tbound.max <= 0) {
            wr_max.set_val(name,bus_s->vbound.max*bus_d->vbound.max*cos(th_max(name).eval()));
            wr_min.set_val(name,bus_s->vbound.min*bus_d->vbound.min*cos(th_min(name).eval()));
            wi_max.set_val(name,bus_s->vbound.min*bus_d->vbound.min*sin(th_max(name).eval()));
            wi_min.set_val(name,bus_s->vbound.max*bus_d->vbound.max*sin(th_min(name).eval()));
        }
        if (arc->tbound.min < 0 && arc->tbound.max > 0) {
            wr_max.set_val(name,bus_s->vbound.max*bus_d->vbound.max);
            wr_min.set_val(name,bus_s->vbound.min*bus_d->vbound.min*min(cos(th_min(name).eval()), cos(th_max(name).eval())));
            wi_max.set_val(name,bus_s->vbound.max*bus_d->vbound.max*sin(th_max(name).eval()));
            wi_min.set_val(name,bus_s->vbound.max*bus_d->vbound.max*sin(th_min(name).eval()));
        }
        cphi.set_val(name, cos(0.5*(arc->tbound.min+arc->tbound.max)));
        sphi.set_val(name, sin(0.5*(arc->tbound.min+arc->tbound.max)));
        cos_d.set_val(name, cos(0.5*(arc->tbound.max-arc->tbound.min)));
        getline(file, word,'\n');
        file >> word;
    }
    DebugOff(ch.to_str(true) << endl);
    DebugOff(as.to_str(true) << endl);
    DebugOff(tr.to_str(true) << endl);    
    file.close();
    return 0;
}


int PowerNet::readgrid(const char* fname) {
    double pi = 4.*atan(1.);
    string name;
    double kvb = 0;
//    int id = 0;
    unsigned index = 0;
    cout << "Loading file " << fname << endl;
    ifstream file(fname, std::ifstream::in);
    if(!file.is_open()) {
        auto new_fname = "../" + string(fname);
        file.open(new_fname);
        if(!file.is_open()) {
            throw invalid_argument("Could not open file\n");
        }
    }
    string word;
    while (word.compare("function")) {
        file >> word;
    }

    file.ignore(6);
    file >> word;
    _name = word;

//  cout << _name << endl;
    while (word.compare("mpc.baseMVA")) {
        file >> word;
    }

    file.ignore(3);
    getline(file, word,';');
    bMVA = atoi(word.c_str());
    /* Nodes data */
    while (word.compare("mpc.bus")) {
        file >> word;
    }

    getline(file, word);
    Bus* bus = NULL;
//    Bus* bus_clone= NULL;
    file >> word;
    int status;
    while(word.compare("];")) {
        name = word.c_str();
        file >> ws >> word;
        status = atoi(word.c_str());
        if (status==3) {
            ref_bus = name;
            DebugOff("Ref bus = " << ref_bus << endl);
        }
        file >> ws >> word;
        pl.set_val(name,atof(word.c_str())/bMVA);
        file >> word;
        ql.set_val(name,atof(word.c_str())/bMVA);
        file >> word;
        gs.set_val(name,atof(word.c_str())/bMVA);
        file >> word;
        bs.set_val(name,atof(word.c_str())/bMVA);
        file >> ws >> word >> ws >> word;
        v_s.set_val(name,atof(word.c_str()));
        file >> ws >> word >> ws >> word;
        kvb = atof(word.c_str());
        file >> ws >> word >> ws >> word;
        v_max.set_val(name,atof(word.c_str()));
        getline(file, word,';');
        v_min.set_val(name,atof(word.c_str()));
        w_min.set_val(name,pow(v_min(name).eval(), 2));
        w_max.set_val(name,pow(v_max(name).eval(), 2));
        // single phase

        bus = new Bus(name, pl(name).eval(), ql(name).eval(), gs(name).eval(), bs(name).eval(), v_min(name).eval(), v_max(name).eval(), kvb, 1);
//        bus_clone = new Bus(name, pl(name).eval(), ql(name).eval(), gs(name).eval(), bs(name).eval(), v_min(name).eval(), v_max(name).eval(), kvb, 1);
        bus->vs = v_s(name).eval();
//        bus_clone->vs = v_s(name).eval();
        if (status>=4) {
            bus->_active = false;
//            bus_clone->_active = false;
        }

        this->Net::add_node(bus);
        if (status>=4) {
            DebugOff("INACTIVE NODE!\n" << name << endl);
        }
        file >> word;
    }
//    ref_bus = nodes.front()->_name;
    file.seekg (0, file.beg);


    /* Generator data */
    while (word.compare("mpc.gen")) {
        file >> word;
    }
//    double qmin = 0, qmax = 0, pmin = 0, pmax = 0, ps = 0, qs = 0;
//    int status = 0;
    getline(file, word);


    file >> word;
//    std::vector<bool> gen_status;
    index = 0;
    string bus_name;
    while(word.compare("];")) {
        bus_name = word.c_str();
        // name -> node.
        bus = (Bus*)(Net::get_node(bus_name));
        name = to_string(index);
        file >> word;
        pg_s.set_val(name,atof(word.c_str())/bMVA);
        file >> word;
        qg_s.set_val(name,atof(word.c_str())/bMVA);
        file >> word;
        qg_max.set_val(name,atof(word.c_str())/bMVA);
        file >> word;
        qg_min.set_val(name,atof(word.c_str())/bMVA);

        file >> ws >> word >> ws >> word >> ws >> word;
        status = atoi(word.c_str());
        file >> word;
        pg_max.set_val(name,atof(word.c_str())/bMVA);

        file >> word;
        pg_min.set_val(name,atof(word.c_str())/bMVA);
        getline(file, word,'\n');
//        gen_status.push_back(status==1);


        bus->_has_gen = true;
        /** generator name, ID */
        Gen* g = new Gen(bus, name, pg_min.eval(index), pg_max.eval(index), qg_min.eval(index), qg_max.eval(index));
        g->_id = index;
        g->_ps = pg_s.eval();
        g->_qs = qg_s.eval();
        gens.push_back(g);
        bus->_gen.push_back(g);
        if(status!=1 || !bus->_active) {
            DebugOff("INACTIVE GENERATOR!\n" << name << endl);
            g->_active = false;
        }
        index++;
//        getline(file, word);
        file >> word;
    }


    file.seekg (0, file.beg);

    /* Generator costs */
    while (word.compare("mpc.gencost")) {
        file >> word;
    }
//    double c0 = 0, c1 = 0,c2 = 0;
    getline(file, word);

    int gen_counter = 0;
    for (int i = 0; i < gens.size(); ++i) {
        file >> ws >> word >> ws >> word >> ws >> word >> ws >> word >> ws >> word;
        c2.set_val(to_string(i),atof(word.c_str())*pow(bMVA,2));
        file >> word;
        c1.set_val(to_string(i),atof(word.c_str())*bMVA);
        file >> word;
        c0.set_val(to_string(i),atof(word.c_str()));
//        c2(i) = atof(word.c_str())*pow(bMVA,2);
//        file >> word;
//        c1(i) = atof(word.c_str())*bMVA;
//        file >> word;
//        c0(i) = atof(word.c_str());
        gens[gen_counter++]->set_costs(c0.eval(), c1.eval(), c2.eval());
        getline(file, word);
    }
    file.seekg (0, file.beg);

    /* Lines data */
    while (word.compare("mpc.branch")) {
        file >> word;
    }
    getline(file, word);
    double res = 0;
    set<string> bus_pair_names;
    Line* arc = NULL;
    string src,dest,key;
    file >> word;
    index = 0;
    bool reversed = false;
    while(word.compare("];")) {
        src = word;
        file >> dest;
        key = dest+","+src;//Taking care of reversed direction arcs
        reversed = false;
        if(get_node(src)->_id > get_node(dest)->_id) {//Reverse arc direction
            //            if(arcID.find(key)!=arcID.end()) {//Reverse arc direction
            DebugOn("Adding arc linking " +src+" and "+dest);
            DebugOn(" with reversed direction, reversing source and destination.\n");
            reversed = true;
            key = src;
            src = dest;
            dest = key;
        }
        
        arc = new Line(to_string(index) + "," + src + "," + dest); // Name of lines
        arc->_id = index++;
        arc->_src = get_node(src);
        arc->_dest= get_node(dest);
        
        file >> word;
        arc->r = atof(word.c_str());
        file >> word;
        arc->x = atof(word.c_str());
        res = pow(arc->r,2) + pow(arc->x,2);
        
        if (res==0) {
            cerr << " line with r = x = 0" << endl;
            exit(-1);
        }
        // define g and b for each conductor.
        arc->g = arc->r/res;
        arc->b = -arc->x/res;
        file >> word;
        arc->ch = atof(word.c_str());
        file >> word;
        arc->limit = atof(word.c_str())/bMVA;
        
        // skip rate A rate B rate C.
        file >> ws >> word >> ws >> word >> ws >> word;
        if(atof(word.c_str()) == 0)
            arc->tr = 1.0;
        else
            arc->tr = atof(word.c_str());
        file >> ws >> word;
        arc->as = atof(word.c_str())*pi/180.;
        file >> ws >> word;
        
        
        
        arc->status = atoi(word.c_str());
        file >> ws >> word;
        
        arc->tbound.min = atof(word.c_str())*pi/180.;
        //        arc->tbound.min = -30*pi/180;
        m_theta_lb += arc->tbound.min;
        file >>  ws >>word;
        
        arc->tbound.max = atof(word.c_str())*pi/180.;
        if (arc->tbound.min==0 && arc->tbound.max==0) {
            DebugOn("Angle bounds are equal to zero. Setting them to -+60");
            arc->tbound.min = -60*pi/180.;
            arc->tbound.max = 60*pi/180.;
            
        }
        if (reversed) {
            arc->g /= pow(arc->tr,2);
            arc->b /= pow(arc->tr,2);
            arc->tr = 1./arc->tr;
            arc->as *= -1.;
            auto temp = arc->tbound.max;
            arc->tbound.max = -1.*arc->tbound.min;
            arc->tbound.min = -1.*temp;
        }
        arc->cc = arc->tr*cos(arc->as); // Rectangular values for transformer phase shifters
        arc->dd = arc->tr*sin(arc->as);
        //        arc->tbound.max = 30*pi/180;
        m_theta_ub += arc->tbound.max;
        
        Bus* bus_s = (Bus*)(arc->_src);
        Bus* bus_d = (Bus*)(arc->_dest);
        
        arc->smax = max(
                        pow(bus_s->vbound.max,2)*(arc->g*arc->g + arc->b*arc->b)*(pow(bus_s->vbound.max,2) + pow(bus_d->vbound.max,2)),
                        pow(bus_d->vbound.max,2)*(arc->g*arc->g+arc->b*arc->b)*(pow(bus_d->vbound.max,2) + pow(bus_s->vbound.max,2))
                        );
        name = arc->_name;
        g.set_val(name,arc->g);
        b.set_val(name,arc->b);
        tr.set_val(name,arc->tr);
        as.set_val(name,arc->as);
        g_ff.set_val(name,arc->g/(pow(arc->cc, 2) + pow(arc->dd, 2)));
        g_ft.set_val(name,(-arc->g*arc->cc + arc->b*arc->dd)/(pow(arc->cc, 2) + pow(arc->dd, 2)));
        
        g_tt.set_val(name,arc->g);
        g_tf.set_val(name,(-arc->g*arc->cc - arc->b*arc->dd)/(pow(arc->cc, 2) + pow(arc->dd, 2)));
        
        
        b_ff.set_val(name,(arc->ch*0.5 + arc->b)/(pow(arc->cc, 2) + pow(arc->dd, 2)));
        b_ft.set_val(name,(-arc->b*arc->cc - arc->g*arc->dd)/(pow(arc->cc, 2) + pow(arc->dd, 2)));
        
        b_tt.set_val(name,(arc->ch*0.5 + arc->b));
        b_tf.set_val(name,(-arc->b*arc->cc + arc->g*arc->dd)/(pow(arc->cc, 2) + pow(arc->dd, 2)));
        
        Y.set_val(name,sqrt(arc->g*arc->g + arc->b*arc->b));
        if (arc->g!=0) {
            Y_t.set_val(name,atan(arc->b/arc->g));
            if(arc->b < 0) {
                Y_t.set_val(name,atan(arc->b/arc->g) - pi);
            } else {
                Y_t.set_val(name,atan(arc->b/arc->g) + pi);
            }
        }
        else {
            if(arc->b < 0){
                Y_t.set_val(name,-pi*0.5);
            } else {
                Y_t.set_val(name,pi*0.5);
            }
        }
        
        Y_charge.set_val(name,sqrt(pow(arc->g,2) + pow((arc->b+arc->ch*0.5),2)));
        if(arc->g != 0) {
            Y_charge_t.set_val(name,atan((arc->b+arc->ch*0.5)/arc->g));
        } else {
            if(arc->b < 0) {
                Y_charge_t.set_val(name,-pi*0.5);
            } else {
                Y_charge_t.set_val(name,pi*0.5);
            }
        }
        ch.set_val(name,arc->ch);
        S_max.set_val(name,arc->limit);
        //        DebugOn("charge = " << arc->ch << endl);
        //        DebugOn("as = " << arc->as << endl);
        //        DebugOn("tr = " << arc->tr << endl);
        
        
        if(arc->status != 1 || !bus_s->_active || !bus_d->_active) {
            arc->_active = false;
            DebugOff("INACTIVE ARC!\n" << arc->_name << endl);
        }
        arc->connect();
        add_arc(arc);
        /* Switching to bus_pairs keys */
        name = bus_s->_name + "," + bus_d->_name;
        if (arc->_active && bus_pair_names.count(name)==0) {
            _bus_pairs._keys.push_back(new index_pair(index_(bus_s->_name), index_(bus_d->_name), arc->_active));
            bus_pair_names.insert(name);
        }
        if (!arc->_parallel) {
            th_min.set_val(name,arc->tbound.min);
            th_max.set_val(name,arc->tbound.max);
            tan_th_min.set_val(name,tan(arc->tbound.min));
            tan_th_max.set_val(name,tan(arc->tbound.max));
            
        }
        else {
            th_min.set_val(name,max(th_min.eval(name), arc->tbound.min));
            th_max.set_val(name,min(th_max.eval(name), arc->tbound.max));
            tan_th_min.set_val(name,tan(th_min.eval(name)));
            tan_th_max.set_val(name,tan(th_max.eval(name)));
        }
        if (arc->tbound.min >= 0) {
            wr_max.set_val(name,bus_s->vbound.max*bus_d->vbound.max*cos(th_min(name).eval()));
            wr_min.set_val(name,bus_s->vbound.min*bus_d->vbound.min*cos(th_max(name).eval()));
            wi_max.set_val(name,bus_s->vbound.max*bus_d->vbound.max*sin(th_max(name).eval()));
            wi_min.set_val(name,bus_s->vbound.min*bus_d->vbound.min*sin(th_min(name).eval()));
        };
        if (arc->tbound.max <= 0) {
            wr_max.set_val(name,bus_s->vbound.max*bus_d->vbound.max*cos(th_max(name).eval()));
            wr_min.set_val(name,bus_s->vbound.min*bus_d->vbound.min*cos(th_min(name).eval()));
            wi_max.set_val(name,bus_s->vbound.min*bus_d->vbound.min*sin(th_max(name).eval()));
            wi_min.set_val(name,bus_s->vbound.max*bus_d->vbound.max*sin(th_min(name).eval()));
        }
        if (arc->tbound.min < 0 && arc->tbound.max > 0) {
            wr_max.set_val(name,bus_s->vbound.max*bus_d->vbound.max);
            wr_min.set_val(name,bus_s->vbound.min*bus_d->vbound.min*min(cos(th_min(name).eval()), cos(th_max(name).eval())));
            wi_max.set_val(name,bus_s->vbound.max*bus_d->vbound.max*sin(th_max(name).eval()));
            wi_min.set_val(name,bus_s->vbound.max*bus_d->vbound.max*sin(th_min(name).eval()));
        }
        cphi.set_val(name, cos(0.5*(arc->tbound.min+arc->tbound.max)));
        sphi.set_val(name, sin(0.5*(arc->tbound.min+arc->tbound.max)));
        cos_d.set_val(name, cos(0.5*(arc->tbound.max-arc->tbound.min)));
        getline(file, word,'\n');
        file >> word;
    }
    DebugOff(ch.to_str(true) << endl);
    DebugOff(as.to_str(true) << endl);
    DebugOff(tr.to_str(true) << endl);
    
    file.close();
    if (nodes.size()>1000) {
        add_3d_nlin = false;
    }
    return 0;
}

/* Create imaginary lines, fill bus_pairs_chord, set lower and upper bounds */
void PowerNet::update_net(){
    string name;
    double cos_max_, cos_min_, sin_max_, sin_min_;
    double wr_max_, wr_min_, wi_max_, wi_min_, w_max_, w_min_;
    Node *src, *dest, *n;
    Arc *new_arc;
    int fixed = 1, id_sorted = 0; //id of the current bag in bags_sorted
    Arc *a12, *a13, *a32;
    std::vector<std::vector<Node*>> bags_sorted;

    // bags are cliques in the chordal completion graph
    for(auto& b: _bags){
        for(int i = 0; i < b.size()-1; i++) {
            for(int j = i+1; j < b.size(); j++) {
                Arc* a = get_arc(b[i]->_name,b[j]->_name);
                if (a==nullptr) {
                    src = get_node(b[i]->_name);
                    dest = get_node(b[j]->_name);
                    new_arc = new Line(to_string((int) arcs.size() + 1));
                    new_arc->_id = arcs.size();
                    new_arc->_src = src;
                    new_arc->_dest = dest;
                    new_arc->_active = false;
                    new_arc->_imaginary = true;
                    new_arc->_free = true;
                    new_arc->connect();
                    add_undirected_arc(new_arc);
                }
            }
        }
    }

    while (fixed != 0) {
        fixed = 0;
        DebugOff("\nNew iteration");
        for(auto b_it = _bags.begin(); b_it != _bags.end();) {
            std::vector<Node*> b = *b_it;
            if(b.size() == 3) {
                DebugOff("\nBag: " << b[0]->_name << ", " << b[1]->_name << ", " << b[2]->_name);
                a12 = get_arc(b[0], b[1]);
                a13 = get_arc(b[0], b[2]);
                a32 = get_arc(b[2], b[1]);
                if ((a12->_free && a13->_free) || (a12->_free && a32->_free) || (a13->_free && a32->_free) ||
                    (!a12->_free && !a13->_free && !a32->_free)) { // at least two missing lines or all lines real
                    ++b_it;
                    continue;
                }
                if (a12->_free) {
                    a12->_free = false;
                    DebugOff("\nFixing arc a12 (" << a12->_src->_name << ", " << a12->_dest->_name << "), adding bag #" << id_sorted);
                    fixed++;
                }
                if (a13->_free) {
                    a13->_free = false;
                    DebugOff("\nFixing arc a13 (" << a13->_src->_name << ", " << a13->_dest->_name << "), adding bag #" << id_sorted);
                    fixed++;
                }
                if (a32->_free) {
                    a32->_free = false;
                    DebugOff("\nFixing arc a32 (" << a32->_src->_name << ", " << a32->_dest->_name << "), adding bag #" << id_sorted);
                    fixed++;
                }
                bags_sorted.push_back(b);
                _bags.erase(b_it);
                id_sorted++;
            }
            else{ // Bags with size > 3; todo: leave only this as the general case?
                DebugOff("\nBag with size > 3");

                for(int i = 0; i < b.size()-1; i++) {
                    for (int j = i + 1; j < b.size(); j++) {
                        Arc* a = get_arc(b[i]->_name, b[j]->_name);
                        if (!a->_free) continue;
                        n = a->_src;
                        //by now, all arcs in bags should be created
                        for (auto n1: b) {
                            if(n==n1) continue;
                            Arc* a2 = get_arc(n->_name, n1->_name);
                            if (a2->_free) continue;
                            Arc *a1 = get_arc(a->_dest, n1);
                            if (!a1->_free) {
                                a->_free = false;

                                vector<Node *> bag;
                                bag.push_back(get_node(n->_name));
                                bag.push_back(get_node(a->_dest->_name));
                                bag.push_back(get_node(n1->_name));
//                                sort(bag.begin(), bag.end(),
//                                     [](const Node *a, const Node *b) -> bool { return a->_id < b->_id; });

                                fixed++;
                                sort(bag.begin(), bag.end(), [](const Node* a, const Node* b) -> bool{return a->_id < b->_id;});
                                bags_sorted.push_back(bag);
                                id_sorted++;
                                DebugOff("\nFixing arc in a larger bag (" << a->_src->_name << ", " << a->_dest->_name << ")");
                                break;
                            }
                        }
                    } // j
                } // i
                ++b_it;

            } // size > 3
        } // bags loop
    } // while

    //add all remaining bags to bags_sorted
    for(auto b_it = _bags.begin(); b_it != _bags.end();) {
        std::vector<Node*> b = *b_it;
            if(b.size() >= 2) bags_sorted.push_back(b);
            _bags.erase(b_it);
//            id_sorted++;
    }
    _bags = bags_sorted;

    for(auto& a: arcs) {
        if(a->_imaginary) a->_free = true;
    }


    for(auto& k: _bus_pairs._keys){
        _bus_pairs_chord._keys.push_back(new index_pair(*k));
    }

    for(auto& a: arcs){
        if(a->_imaginary){
            Bus* bus_s = (Bus*)(a->_src);
            Bus* bus_d = (Bus*)(a->_dest);

            name = bus_s->_name + "," + bus_d->_name;
            _bus_pairs_chord._keys.push_back(new index_pair(index_(bus_s->_name), index_(bus_d->_name)));

            if (m_theta_lb < -3.14 && m_theta_ub > 3.14) {
                cos_max_ = 1;
                cos_min_ = -1;
            } else if (m_theta_lb < 0 && m_theta_ub > 0){
                cos_max_ = 1;
                cos_min_ = min(cos(m_theta_lb), cos(m_theta_ub));
            } else{
                cos_max_ = max(cos(m_theta_lb),cos(m_theta_ub));
                cos_min_ = min(cos(m_theta_lb), cos(m_theta_ub));
            }
            w_max_ = bus_s->vbound.max*bus_d->vbound.max;
            w_min_ = bus_s->vbound.min*bus_d->vbound.min;

            wr_max_ = cos_max_*w_max_;
            if(cos_min_ < 0) wr_min_ = cos_min_*w_max_;
            else wr_min_ = cos_min_*w_min_;

            if(m_theta_lb < -1.57 && m_theta_ub > 1.57){
                sin_max_ = 1;
                sin_min_ = -1;
            } else{
                sin_max_ = sin(m_theta_ub);
                sin_min_ = sin(m_theta_lb);
            }

            if(sin_max_ > 0) wi_max_ = sin_max_*w_max_;
            else wi_max_ = sin_max_*w_min_;
            if(sin_min_ > 0) wi_min_ = sin_min_*w_min_;
            else wi_min_ = sin_min_*w_max_;

//            cout << "\nImaginary line, bounds: (" << wr_min_ << "," << wr_max_ << "); (" << wi_min_ << "," << wi_max_ << ")";

            wr_max.set_val(name,wr_max_);
            wr_min.set_val(name,wr_min_);
            wi_max.set_val(name,wi_max_);
            wi_min.set_val(name,wi_min_);
        }
    }
    DebugOff("\nBags sorted: " << endl);
    for(auto& b: _bags) {
        DebugOff("bag = {");
        for (int i = 0; i < b.size(); i++) {
            DebugOff(b.at(i)->_name << " ");
        }
        DebugOff("}" << endl);
        if(add_3d_nlin && b.size()==3){
            for(int i = 0; i < 2; i++) {
                for(int j = i+1; j < 3; j++) {
                    Arc* aij = get_arc(b[i],b[j]);
                    aij->_free = false;
                }
            }
        }
    }
}

//void PowerNet::update_status(unique_ptr<Model> model){
//    for (auto &v_p:model->_vars) {
//        if (v_p.second->is_integer() || v_p.second->is_binary()) {
//            auto v = v_p.second;
//            auto new_v = new var<Real>(v_p.second->_name, 0,1);
//            new_v->copy(*v);
//            new_v->_is_relaxed = true;
//            new_v->_val->resize(new_v->get_dim());
//            if (v->get_intype()==integer_) {
//                auto real_var = (var<int>*)v;
//                for (int i = 0; i < real_var->get_dim(); i++) {
//                    new_v->_val->at(i) = real_var->_val->at(i);
//                }
//            }
//            if (v->get_intype()==short_) {
//                auto real_var = (var<short>*)v;
//                for (int i = 0; i < real_var->get_dim(); i++) {
//                    new_v->_val->at(i) = real_var->_val->at(i);
//                }
//            }
//            if (v->get_intype()==binary_) {
//                auto real_var = (var<bool>*)v;
//                _model->_bin_vars[v_p.second->get_vec_id()] = *real_var;
//                for (int i = 0; i < real_var->get_dim(); i++) {
//                    new_v->_val->at(i) = real_var->_val->at(i);
//                }
//            }
//            v_p.second = new_v;
//        }
//    }
//}

shared_ptr<Model> PowerNet::build_SCOPF_gen_contingency(int cont, const string& name, PowerModelType pmt, int output, double tol){
    /** MODEL DECLARATION */
    shared_ptr<Model> CSCOPF(new Model(name));
    return CSCOPF;
}

shared_ptr<Model> PowerNet::build_SCOPF_line_contingency(int cont, const string& name, PowerModelType pmt, int output, double tol){
    /** MODEL DECLARATION */
    shared_ptr<Model> CSCOPF(new Model(name));
    bool polar = (pmt==ACPOL);
    if (polar) {
        DebugOn("Using polar model\n");
    }
    else {
        DebugOn("Using rectangular model\n");
    }
    /** Variables */
    /* Power generation variables */
    var<Real> Pg("Pg_cont", pg_min, pg_max);
    var<Real> Qg ("Qg_cont", qg_min, qg_max);
    var<Real> Sp ("Sp", pos_);
    var<Real> Sn ("Sn", pos_);
    CSCOPF->add(Pg.in(gens), Qg.in(gens), Sp.in(gens), Sn.in(gens));
    Qg_cont[cont] = Qg;
    
    /* Power flow variables */
    var<Real> Pf_from("Pf_from", S_max);
    var<Real> Qf_from("Qf_from", S_max);
    var<Real> Pf_to("Pf_to", S_max);
    var<Real> Qf_to("Qf_to", S_max);
    CSCOPF->add(Pf_from.in(arcs), Qf_from.in(arcs),Pf_to.in(arcs),Qf_to.in(arcs));
    
    p_from[cont] = Pf_from;
    q_from[cont] = Qf_from;
    p_to[cont] = Pf_to;
    q_to[cont] = Qf_to;
    
    var<Real> delta("delta", pos_);
    CSCOPF->add(delta.in(R(1)));
    _delta[cont] = delta;
    
    indices Ng; /*<< Nodes with generators */
    for (auto node: nodes) {
        auto bus = (Bus*)node;
        if (bus->_gen.size()>0) {
            Ng.add(bus->_name);
        }
    }
    var<Real> v_diff("v_diff", pos_);
    CSCOPF->add(v_diff.in(Ng));
    
    /** Voltage related variables */
    var<Real> theta("theta");
    var<Real> v("v_cont", v_min, v_max);
    var<Real> vr("vr");
    var<Real> vi("vi");
    
    if (polar) {
        CSCOPF->add_var(v.in(nodes));
        CSCOPF->add_var(theta.in(nodes));
        v.initialize_all(1);
    }
    else {
        CSCOPF->add_var(vr.in(nodes));
        CSCOPF->add_var(vi.in(nodes));
        vr.initialize_all(1.0);
    }
    v_cont[cont] = v;
    theta_cont[cont] = theta;
    
    /** Construct the objective function */
    CSCOPF->min(sum(v_diff) + sum(Sp) + sum(Sn));
//    CSCOPF->min(1);
    
/** Define constraints */
    
    /* OBJ Constraints */
    
    Constraint V_diff_p("V_diff_p");
    V_diff_p = v_diff - (v - v_base);
    CSCOPF->add(V_diff_p.in(Ng) >= 0);
    
    Constraint V_diff_n("V_diff_n");
    V_diff_n = v_diff - (v_base - v);
    CSCOPF->add(V_diff_n.in(Ng) >= 0);
//    CSCOPF->get_constraint("V_diff_n")->print_expanded();
    
    /* Generator Constraints */
    Constraint Gen_resp("g_resp");
    Gen_resp = Pg - (Pg_base + p_factor*delta) - (Sp - Sn);
    CSCOPF->add(Gen_resp.in(indices(gens)) == 0);
//    CSCOPF->get_constraint("g_resp")->print_expanded();
    /* REF BUS */
    Constraint Ref_Bus("Ref_Bus");
    if (polar) {
        Ref_Bus = theta(get_ref_bus());
    }
    else {
        Ref_Bus = vi(get_ref_bus());
    }
    CSCOPF->add(Ref_Bus == 0);
    
    /** KCL Flow conservation */
    Constraint KCL_P("KCL_P");
    Constraint KCL_Q("KCL_Q");
    KCL_P  = sum(Pf_from.out_arcs()) + sum(Pf_to.in_arcs()) + pl - sum(Pg.in_gens());
    KCL_Q  = sum(Qf_from.out_arcs()) + sum(Qf_to.in_arcs()) + ql - sum(Qg.in_gens());
    /* Shunts */
    if (polar) {
        KCL_P +=  gs*power(v,2);
        KCL_Q -=  bs*power(v,2);
    }
    else {
        KCL_P +=  gs*(power(vr,2)+power(vi,2));
        KCL_Q -=  bs*(power(vr,2)+power(vi,2));
    }
    CSCOPF->add(KCL_P.in(nodes) == 0);
    CSCOPF->add(KCL_Q.in(nodes) == 0);
    
    /** AC Power Flows */
    Constraint Flow_P_From("Flow_P_From");
    Flow_P_From += Pf_from;
    if (polar) {
        Flow_P_From -= g/power(tr,2)*power(v.from(),2);
        Flow_P_From += g/tr*(v.from()*v.to()*cos(theta.from() - theta.to() - as));
        Flow_P_From += b/tr*(v.from()*v.to()*sin(theta.from() - theta.to() - as));
    }
    else {
        Flow_P_From -= g_ff*(power(vr.from(), 2) + power(vi.from(), 2));
        Flow_P_From -= g_ft*(vr.from()*vr.to() + vi.from()*vi.to());
        Flow_P_From -= b_ft*(vi.from()*vr.to() - vr.from()*vi.to());
    }
    CSCOPF->add(Flow_P_From.in(arcs)==0);
    
    Constraint Flow_P_To("Flow_P_To");
    Flow_P_To += Pf_to;
    if (polar) {
        Flow_P_To -= g*power(v.to(), 2);
        Flow_P_To += g/tr*(v.from()*v.to()*cos(theta.to() - theta.from() + as));
        Flow_P_To += b/tr*(v.from()*v.to()*sin(theta.to() - theta.from() + as));
    }
    else {
        Flow_P_To -= g_tt*(power(vr.to(), 2) + power(vi.to(), 2));
        Flow_P_To -= g_tf*(vr.from()*vr.to() + vi.from()*vi.to());
        Flow_P_To -= b_tf*(vi.to()*vr.from() - vr.to()*vi.from());
    }
    CSCOPF->add(Flow_P_To.in(arcs)==0);
    
    Constraint Flow_Q_From("Flow_Q_From");
    Flow_Q_From += Qf_from;
    if (polar) {
        Flow_Q_From += (0.5*ch+b)/power(tr,2)*power(v.from(),2);
        Flow_Q_From -= b/tr*(v.from()*v.to()*cos(theta.from() - theta.to() - as));
        Flow_Q_From += g/tr*(v.from()*v.to()*sin(theta.from() - theta.to() - as));
    }
    else {
        Flow_Q_From += b_ff*(power(vr.from(), 2) + power(vi.from(), 2));
        Flow_Q_From += b_ft*(vr.from()*vr.to() + vi.from()*vi.to());
        Flow_Q_From -= g_ft*(vi.from()*vr.to() - vr.from()*vi.to());
    }
    CSCOPF->add(Flow_Q_From.in(arcs)==0);
    Constraint Flow_Q_To("Flow_Q_To");
    Flow_Q_To += Qf_to;
    if (polar) {
        Flow_Q_To += (0.5*ch+b)*power(v.to(),2);
        Flow_Q_To -= b/tr*(v.from()*v.to()*cos(theta.to() - theta.from() + as));
        Flow_Q_To += g/tr*(v.from()*v.to()*sin(theta.to() - theta.from() + as));
    }
    else {
        Flow_Q_To += b_tt*(power(vr.to(), 2) + power(vi.to(), 2));
        Flow_Q_To += b_tf*(vr.from()*vr.to() + vi.from()*vi.to());
        Flow_Q_To -= g_tf*(vi.to()*vr.from() - vr.to()*vi.from());
    }
    CSCOPF->add(Flow_Q_To.in(arcs)==0);
    
    /** AC voltage limit constraints. */
    if (!polar) {
        Constraint Vol_limit_UB("Vol_limit_UB");
        Vol_limit_UB = power(vr, 2) + power(vi, 2);
        Vol_limit_UB -= power(v_max, 2);
        CSCOPF->add(Vol_limit_UB.in(nodes) <= 0);
        
        Constraint Vol_limit_LB("Vol_limit_LB");
        Vol_limit_LB = power(vr, 2) + power(vi, 2);
        Vol_limit_LB -= power(v_min,2);
        CSCOPF->add(Vol_limit_LB.in(nodes) >= 0);
    }
    
    
    
    /*  Thermal Limit Constraints */
    Constraint Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from += power(Pf_from, 2) + power(Qf_from, 2);
    Thermal_Limit_from -= power(S_max, 2);
    CSCOPF->add(Thermal_Limit_from.in(arcs) <= 0);
    
    Constraint Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to += power(Pf_to, 2) + power(Qf_to, 2);
    Thermal_Limit_to -= power(S_max,2);
    CSCOPF->add(Thermal_Limit_to.in(arcs) <= 0);
//    CSCOPF->print_expanded();
    return CSCOPF;
}


shared_ptr<Model> PowerNet::build_ROMDST_contingency(const string& name, PowerModelType pmt, int output, double tol, int max_nb_hours){
    
    
    auto new_Et = Et_opt;
    auto new_Gt = indices(gens,T);
    auto new_Bt = Bt_opt;
    auto new_Bt1 = Bt1_opt; /**< Excluding first hour */
    auto new_PVt = PVt_opt;
    
    /** MODEL DECLARATION */
    shared_ptr<Model> ROMDST(new Model(name));
    /** VARIABLES */
    
    
    /* Diesel power generation variables */
    var<Real> Pg("Pg", pg_min, pg_max);
    var<Real> Qg ("Qg", qg_min, qg_max);
    var<Real> Pg_ ("Pg_", pg_min, pg_max);/**< Active power generation before losses */
    ROMDST->add(Pg.in(new_Gt),Pg_.in(new_Gt),Qg.in(new_Gt));
    
    
    /* Battery power generation variables */
    var<Real> Pb("Pb", pb_min, pb_max);/**< Active power generation outside the battery */
    var<Real> Pb_diff("Pb_diff", pos_);/**< Difference in active power generation outside and inside the battery */
    var<Real> Qb ("Qb", qb_min, qb_max);/**< Reactive power generation outside the battery */
    var<Real> Pb_("Pb_", pb_min, pb_max);/**< Active power generation in the battery */
//    var<Real> Qb_("Qb_", qb_min, qb_max);/**< Reactive power generation in the battery */
    ROMDST->add(Pb.in(new_Bt), Qb.in(new_Bt), Pb_.in(new_Bt), Pb_diff.in(new_Bt));
    
    
    /* PV power generation variables */
    var<Real> Pv("Pv", pos_);
    ROMDST->add(Pv.in(new_PVt));
    
    /* Battery state of charge variables */
    var<Real> Sc("Sc", pos_);
    ROMDST->add(Sc.in(new_Bt));
    
    /* Wind power generation variables */
    var<Real> Pw("Pw", pw_min, pw_max);
    //    pw_max.print(true);
    ROMDST->add(Pw.in(Wt));
    
    /* Power flow variables */
    var<Real> Pij("Pfrom", S_max);
    var<Real> Qij("Qfrom", S_max);
    var<Real> Pji("Pto", S_max);
    var<Real> Qji("Qto", S_max);
    
    var<Real> P_shed("P_shed", pos_);
    var<Real> P_shed_max("P_shed_max", pos_);
//    var<Real> Q_shed("Q_shed", pos_);
    
    ROMDST->add(Pij.in(new_Et),Qij.in(new_Et), P_shed.in(Nt), P_shed_max.in(R(1)));
    if (pmt!=LDISTF) {
        ROMDST->add(Pji.in(new_Et),Qji.in(new_Et));
    }
    
    /** Voltage magnitude (squared) variables */
    var<Real> v("v", vmin, vmax);
    ROMDST->add(v.in(Nt));
    v.initialize_all(0.64);
    
    /** Power loss variables */
    var<Real> loss("loss", pos_);
    if (pmt==DISTF || pmt==CDISTF) {
        ROMDST->add(loss.in(Nt));
    }
    
    /** OBJECTIVE FUNCTION */
    ROMDST->min(100*P_shed_max);
    
    
    /** CONSTRAINTS **/
    
/** Max Load Shed **/
    
    Constraint max_load_shed("max_load_shed");
//    max_load_shed = bMVA*P_shed/(pl*bMVA-1e-5) - P_shed_max;//TODO fix this
    max_load_shed = P_shed*pl_ratio - P_shed_max;
    ROMDST->add(max_load_shed.in(Nt_load)<=0);
//    ROMDST->get_constraint("max_load_shed")->print_expanded();

    Constraint no_load_shed("no_load_shed");
    no_load_shed = P_shed;
    ROMDST->add(no_load_shed.in(Nt_no_load)==0);
//    ROMDST->get_constraint("no_load_shed")->print_expanded();

/** FLOW CONSERVATION **/
    
    /** KCL P and Q */
    Constraint KCL_P("KCL_P");
    Constraint KCL_Q("KCL_Q");
    if (pmt==LDISTF) {
        KCL_P  = sum(Pij.out_arcs()) - sum(Pij.in_arcs()) + pl - sum(Pg.in_gens()) - sum(Pb.in_bats()) - sum(Pw.in_wind()) - sum(Pv.in_pv()) - P_shed;
        //        KCL_P  = sum(Pij.out_arcs()) - sum(Pij.in_arcs()) + pl - sum(Pg.in_gens()) - sum(Pb.in_bats()) - sum(Pw.in_wind());
        KCL_Q  = sum(Qij.out_arcs()) - sum(Qij.in_arcs()) + ql - sum(Qg.in_gens()) - sum(Qb.in_bats());
    }
    else{
        KCL_P  = sum(Pij.out_arcs()) + sum(Pji.in_arcs()) + pl - sum(Pg.in_gens())- sum(Pb.in_bats()) - sum(Pw.in_wind()) - sum(Pv.in_pv()) - P_shed;
        KCL_Q  = sum(Qij.out_arcs()) + sum(Qji.in_arcs()) + ql - sum(Qg.in_gens()) - sum(Qb.in_bats());
    }
    ROMDST->add(KCL_P.in(nodes, T) == 0);
    ROMDST->add(KCL_Q.in(nodes, T) == 0);
    
/**  THERMAL LIMITS **/
    
    
    
    /*  Thermal Limit Constraints */
    Constraint Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from += power(Pij, 2) + power(Qij, 2);
    Thermal_Limit_from -= power(S_max, 2);
    ROMDST->add(Thermal_Limit_from.in(new_Et) <= 0);
    
    
    
/**  POWER FLOW **/
    
    /** Power Flow Constraints */
    Constraint Flow_P_From("Flow_P_From");
    Flow_P_From = v.from() - 2*(r*Pij + x*Qij) - v.to();
    ROMDST->add(Flow_P_From.in(arcs, T)==0);
    
    
    /**  PV **/
    
    /*  On/Off on newly installed PV */
    Constraint OnOffPV("OnOffPV");
    OnOffPV += Pv - Pv_cap_*pv_out;
    ROMDST->add(OnOffPV.in(new_PVt) <= 0);
    
    
/**  BATTERIES **/
    
    /*  Apparent Power Limit on existing and installed Batteries */
    Constraint Apparent_Limit_Batt_Pot("Apparent_Limit_Batt");
    Apparent_Limit_Batt_Pot += power(Pb, 2) + power(Qb, 2);
    Apparent_Limit_Batt_Pot -= power(pb_max, 2);
    ROMDST->add(Apparent_Limit_Batt_Pot.in(new_Bt) <= 0);
    
    /*  State Of Charge */
    Constraint State_Of_Charge("State_Of_Charge");
    State_Of_Charge = Sc - Sc.prev() + Pb_;
    ROMDST->add(State_Of_Charge.in(new_Bt1) == 0);
    
    /*  State Of Charge 0 */
    auto Bat0 = indices(_battery_inverters,T.start());
    Constraint State_Of_Charge0("State_Of_Charge0");
    State_Of_Charge0 = Sc;
    ROMDST->add(State_Of_Charge0.in(Bat0) == 0);
    Constraint Pb0("Pb0");
    Pb0 = Pb;
    ROMDST->add(Pb0.in(Bat0) == 0);
    
    /*  Efficiencies */
    Constraint DieselEff("DieselEff");
    DieselEff += Pg - gen_eff*Pg_;
    ROMDST->add(DieselEff.in(new_Gt) == 0);
    
    Constraint EfficiencyExist("BatteryEfficiency");
    EfficiencyExist += Pb  - eff_a*Pb_ - eff_b;//TODO without time extending eff_a and eff_b
    ROMDST->add(EfficiencyExist.in(indices(_eff_pieces,new_Bt)) <= 0);
    
//    Constraint PbDiff_P("PbDiff_P");
//    PbDiff_P += Pb - Pb_ - Pb_diff;
//    ROMDST->add(PbDiff_P.in(new_Bt) <= 0);
//    
//    Constraint PbDiff_N("PbDiff_N");
//    PbDiff_N += Pb_ - Pb - Pb_diff;
//    ROMDST->add(PbDiff_N.in(new_Bt) <= 0);
    
//    ROMDST->_obj -= *ROMDST->get_constraint("BatteryEfficiency");

/*  GENERATOR RAMP UP/DOWN */
    Constraint RampDown("RampDown");
    RampDown +=  Pg_ - (Pg_base - ramp_down*pg_max);
    ROMDST->add(RampDown.in(new_Gt) >= 0);
    
    Constraint RampUp("RampUp");
    RampUp +=  Pg_ - (Pg_base + ramp_up*pg_max);
    ROMDST->add(RampUp.in(new_Gt) <= 0);
    
    return ROMDST;
}

void PowerNet::save_base_case_sol(const string& fname){
    ofstream myfile;
    myfile.open (fname);
    myfile << "--generation dispatch\n";
    myfile << "bus id,unit id,pg(MW),qg(MVar)\n";
    auto pg_vals = Pg_base.get_vals();
    auto qg_vals = Qg_base.get_vals();
    unsigned idx = 0;
    for (auto gen:gens) {
        myfile << gen->_bus->_name << "," << gen->_unit_name << "," << to_string_with_precision(pg_vals->at(idx)*bMVA,12) << "," << to_string_with_precision(qg_vals->at(idx)*bMVA,12) << endl;
        idx++;
    }
    myfile << "--end of generation dispatch" << endl;
    myfile.close();
}

void PowerNet::save_all_sol(const string& fname){
    double pi = 4.*atan(1.);
    ofstream myfile;
    myfile.open (fname);
    myfile << "--contingency generator\n";
    myfile << "conID,genID,busID,unitID,q(MW)\n";
    unsigned idx = 0;
    auto cont = 0;
    for (cont = 0; cont<_conting_lines.size(); cont++) {
        auto qg_vals = Qg_cont[cont].get_vals();
        idx = 0;
        for (auto gen:gens) {
            myfile << cont+1 << "," << gen->_name<< "," << gen->_bus->_name << "," << gen->_unit_name << "," << to_string_with_precision(qg_vals->at(idx)*bMVA,12) << endl;
            idx++;
        }
    }
    auto g_idx = 0;
    for (; cont<_nb_conting; cont++) {
        auto qg_vals = Qg_cont[cont].get_vals();
        idx = 0;
        for (auto gen:gens) {
            if (_conting_gens[g_idx]==gen) {
                continue;
            }
            myfile << cont+1 << "," << gen->_name<< "," << gen->_bus->_name << "," << gen->_unit_name << "," << to_string_with_precision(qg_vals->at(idx)*bMVA,12) << endl;
            idx++;
        }
        g_idx++;
    }
    myfile << "--end of contingency generator" << endl;
    myfile << "--bus" << endl;
    myfile << "contingency id,bus id,v(pu),theta(deg)" << endl;
    idx = 0;
    auto v_vals = v_base.get_vals();
    auto theta_vals = theta_base.get_vals();
    for (auto node:nodes) {
        myfile << "0," << node->_name << "," << to_string_with_precision(v_vals->at(idx),12) << "," << to_string_with_precision(theta_vals->at(idx)*180./pi,12) << endl;
        idx++;
    }
    for (auto cont = 0; cont<_nb_conting; cont++) {
        auto v_vals = v_cont[cont].get_vals();
        auto theta_vals = theta_cont[cont].get_vals();
        idx = 0;
        for (auto node:nodes) {
            myfile << cont+1 << "," << node->_name << "," << to_string_with_precision(v_vals->at(idx),12) << "," << to_string_with_precision(theta_vals->at(idx)*180./pi,12) << endl;
            idx++;
        }
    }
    myfile << "--end of bus" << endl;
    myfile << "--Delta" << endl;
    myfile << "contingency id,Delta(MW)" << endl;
    for (auto cont = 0; cont<_nb_conting; cont++) {
        myfile << cont+1 << ","<< to_string_with_precision(_delta[cont].eval()*bMVA,12) << endl;
    }
    myfile << "--end of Delta" << endl;
    myfile << "--line flow" << endl;
    myfile << "contingency id,line id,origin bus id,destination bus id,circuit id,p_origin(MW),q_origin(MVar),p_destination(MW),q_destination(MVar)" << endl;
    auto p_from_vals = p_from_base.get_vals();
    auto q_from_vals = q_from_base.get_vals();
    auto p_to_vals = p_to_base.get_vals();
    auto q_to_vals = q_to_base.get_vals();
    idx = 0;
    for (auto arc:arcs) {
        myfile << "0," << arc->_name << "," << arc->_circuit_id << "," << to_string_with_precision(p_from_vals->at(idx)*bMVA,12) << "," << to_string_with_precision(q_from_vals->at(idx)*bMVA,12) << "," << to_string_with_precision(p_to_vals->at(idx)*bMVA,12) << "," << to_string_with_precision(q_to_vals->at(idx)*bMVA,12) << endl;
        idx++;
    }
    cont = 0;
    /** Line Contingencies **/
    for (cont = 0; cont<_conting_lines.size(); cont++) {
        auto p_from_vals = p_from[cont].get_vals();
        auto q_from_vals = q_from[cont].get_vals();
        auto p_to_vals = p_to[cont].get_vals();
        auto q_to_vals = q_to[cont].get_vals();
        idx = 0;
        for (auto arc:arcs) {
            if (arc == _conting_lines[cont]) {
                continue;
            }
            myfile << cont+1 << "," << arc->_name << "," << arc->_circuit_id << "," << to_string_with_precision(p_from_vals->at(idx)*bMVA,12) << "," << to_string_with_precision(q_from_vals->at(idx)*bMVA,12)<< "," << to_string_with_precision(p_to_vals->at(idx)*bMVA,12) << "," << to_string_with_precision(q_to_vals->at(idx)*bMVA,12) << endl;
            idx++;
        }
    }
    /** Gen Contingencies **/
    for (; cont<_nb_conting; cont++) {
        auto p_from_vals = p_from[cont].get_vals();
        auto q_from_vals = q_from[cont].get_vals();
        auto p_to_vals = p_to[cont].get_vals();
        auto q_to_vals = q_to[cont].get_vals();
        idx = 0;
        for (auto arc:arcs) {
            myfile << cont+1 << "," << arc->_name << "," << arc->_circuit_id << "," << to_string_with_precision(p_from_vals->at(idx)*bMVA,12) << "," << to_string_with_precision(q_from_vals->at(idx)*bMVA,12)<< "," << to_string_with_precision(p_to_vals->at(idx)*bMVA,12) << "," << to_string_with_precision(q_to_vals->at(idx)*bMVA,12) << endl;
            idx++;
        }
    }
    myfile << "--end of line flow" << endl;
    myfile.close();
}

void PowerNet::fix_investment(){
    auto vals = w_g.get_vals();
    for (auto i = 0; i<vals->size(); i++) {
        if (!vals->at(i)) {//This generator was not built
            clog << "Excluding Generator: " << _potential_diesel_gens[i]->_name << endl;
            _potential_diesel_gens[i]->_active = false;
        }
    }
    vals = w_b.get_vals();
    for (auto i = 0; i<vals->size(); i++) {
        if (!vals->at(i)) {//This battery was not built
            clog << "Excluding Battery: " << _potential_battery_inverters[i]->_name << endl;
            _potential_battery_inverters[i]->_active = false;
        }
    }
    vals = w_pv.get_vals();
    for (auto i = 0; i<vals->size(); i++) {
        if (!vals->at(i)) {//This PV was not built
            clog << "Excluding PV: " << _potential_PV_gens[i]->_name << endl;
            _potential_PV_gens[i]->_active = false;
        }
    }
    vals = w_e.get_vals();
    for (auto i = 0; i<vals->size(); i++) {
        if (!vals->at(i)) {//This edge was not built
            clog << "Excluding Edge Expansion: " << _potential_expansion[i]->_name << endl;
            _potential_expansion[i]->_active = false;
        }
    }
    Pv_cap_ = Pv_cap;
    for (auto pv:_existing_PV_gens) {
        Pv_cap_(pv->_name) = pv->_max_cap;
    }
    Pg_base = Pg_;
    Et_opt = indices(arcs, T);
    Gt_opt = indices(gens, T);
    Bt_opt = indices(_battery_inverters, T);
    Bt1_opt = indices(_battery_inverters,T.exclude(T.start())); /**< Excluding first hour */
    PVt_opt = indices(_all_PV_gens,T);
}

unique_ptr<Model> PowerNet::build_SCOPF(PowerModelType pmt, int output, double tol){
    /** MODEL DECLARATION */
    unique_ptr<Model> SCOPF(new Model("SCOPF Model"));
    return SCOPF;
}

unique_ptr<Model> PowerNet::build_ROMDST(PowerModelType pmt, int output, double tol, int max_nb_hours){
    /** Indices Sets */
    hours = time(1,max_nb_hours); /**< Hours */
//        indices months = time("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec"); /**< Months */
//    indices months = time("jan","feb","mar","apr","may","jun"); /**< Months */
//    months = time("apr", "aug", "dec"); /**< Months */
        indices months = time("jan");
    //        indices months = time("jan", "feb", "mar","apr","may","jun");
    typical_days = time("week","peak","weekend");
//        indices typical_days = time("week");
    T = indices(months,typical_days,hours);
    double nT = T.size();
    DebugOn("number of time periods = " << nT << endl);
    time_expand(T);
    Nt = indices(nodes,T);
    Nt_load._time_extended = true;
    Nt_no_load._time_extended = true;
    for (auto &nt:*Nt._indices) {
        if (pl.eval(nt)!=0) {
            Nt_load.add(nt);
        }
        else{
            Nt_no_load.add(nt);
        }
    }
    Et = indices(arcs,T);
    Gt = indices(gens,T);
    exist_Gt = indices(_existing_diesel_gens,T);
    exist_Bt = indices(_existing_battery_inverters,T);
    exist_Et = indices(_exist_arcs,T);
    pot_Gt = indices(_potential_diesel_gens,T);
    pot_Bt = indices(_potential_battery_inverters,T);
    pot_Et = indices(_potential_expansion,T);
    Bt = indices(_battery_inverters,T);
    auto T1 = T.exclude(T.start());/**< Excluding first time step */
    Bt1 = indices(_battery_inverters,T1);
    Gt1 = indices(gens,T1);
    Wt = indices(_wind_gens,T);
    PVt = indices(_all_PV_gens,T);
    PV_pot_t = indices(_potential_PV_gens,T);
    pot_gen = indices(_potential_diesel_gens);
    pot_batt = indices(_potential_battery_inverters);
    pot_edges = indices(_potential_expansion);
    
    
    /** MODEL DECLARATION */
    unique_ptr<Model> ROMDST(new Model("ROMDST Model"));
    /** VARIABLES */
    
    
    /* Investment binaries */
    
//    var<Real> Pv_cap("Pv_cap", 0, pv_max); /**< Real variable indicating the extra capacity of PV to be installed on bus b */
    var<Real> Pv_cap("Pv_cap", pos_); /**< Real variable indicating the extra capacity of PV to be installed on bus b */
    auto pot_pv = indices(_potential_PV_gens);
    ROMDST->add(Pv_cap.in(pot_pv));
//        var<Real> w_g("w_g",0,1); /**< Binary variable indicating if generator g is built on bus */
//        var<Real> w_b("w_b",0,1); /**< Binary variable indicating if battery b is built on bus */
//        var<Real> w_e("w_e",0,1); /**< Binary variable indicating if expansion is selected for edge e */
//        var<Real> w_pv("w_pv",0,1); /**< Binary variable indicating if expansion is selected for edge e */
//        w_g._is_relaxed = true;
//        w_b._is_relaxed = true;
//        w_e._is_relaxed = true;
//        w_pv._is_relaxed = true;
    
    var<bool> w_g("w_g"); /**< Binary variable indicating if generator g is built on bus */
    var<bool> w_b("w_b"); /**< Binary variable indicating if battery b is built on bus */
    var<bool> w_e("w_e"); /**< Binary variable indicating if expansion is selected for edge e */
    var<bool> w_pv("w_pv"); /**< Binary variable indicating if PV is installed on bus b */
    
    
    
    ROMDST->add(w_g.in(pot_gen),w_b.in(pot_batt),w_e.in(pot_edges),w_pv.in(pot_pv));
    w_g.initialize_all(1);
    w_b.initialize_all(1);
    w_e.initialize_all(1);
    w_pv.initialize_all(1);
    
    this->w_g = w_g;
    this->w_b = w_b;
    this->w_e = w_e;
    this->w_pv = w_pv;
    this->Pv_cap = Pv_cap;
    
    DebugOff("size w_g = " << w_g.get_dim() << endl);
    DebugOff("size w_b = " << w_b.get_dim() << endl);
    DebugOff("size w_e = " << w_e.get_dim() << endl);
    DebugOff("size w_pv = " << w_pv.get_dim() << endl);
    DebugOff("size Pv_cap = " << Pv_cap.get_dim() << endl);
    
    
    /* Diesel power generation variables */
    var<Real> Pg("Pg", pg_min, pg_max);
    var<Real> Qg ("Qg", qg_min, qg_max);
    var<Real> Pg_ ("Pg_", pg_min, pg_max);/**< Active power generation before losses */
//    var<Real> Pg2("Pg2", 0, power(pg_max,2));/**< Square of Pg */
    var<Real> Pg2("Pg2", pos_);/**< Square of Pg */
    ROMDST->add(Pg.in(Gt),Pg_.in(Gt),Qg.in(Gt),Pg2.in(pot_Gt));
    DebugOff("size Pg = " << Pg.get_dim() << endl);
    DebugOff("size Pg_ = " << Pg_.get_dim() << endl);
    DebugOff("size Qg = " << Qg.get_dim() << endl);
    DebugOff("size Pg2 = " << Pg2.get_dim() << endl);
//    Pg.initialize_all(0.001);
//    Qg.initialize_all(0.0001);
    
    this->Pg_ = Pg_;    
    
    /* Battery power generation variables */
    var<Real> Pb("Pb", pb_min, pb_max);/**< Active power generation outside the battery */
    var<Real> Qb ("Qb", qb_min, qb_max);/**< Reactive power generation outside the battery */
    var<Real> Pb_("Pb_", pb_min, pb_max);/**< Active power generation in the battery */
//    var<Real> Qb_("Qb_", qb_min, qb_max);/**< Reactive power generation in the battery */
    ROMDST->add(Pb.in(Bt), Qb.in(Bt), Pb_.in(Bt));
    DebugOff("size Pb = " << Pb.get_dim() << endl);
    DebugOff("size Qb = " << Qb.get_dim() << endl);
    
    
    /* PV power generation variables */
//    var<Real> Pv("Pv", 0,pv_max);
    var<Real> Pv("Pv", pos_);
    ROMDST->add(Pv.in(PVt));
    DebugOff("size Pv = " << Pv.get_dim() << endl);
    
    /* Battery state of charge variables */
    var<Real> Sc("Sc", pos_);
    ROMDST->add(Sc.in(Bt));
    DebugOff("size Sc = " << Sc.get_dim() << endl);
    
    /* Wind power generation variables */
    var<Real> Pw("Pw", pw_min, pw_max);
    //    pw_max.print(true);
    ROMDST->add(Pw.in(Wt));
    DebugOff("size Pw = " << Pw.get_dim() << endl);
    
    /* Power flow variables */
    var<Real> Pij("Pfrom", S_max);
    var<Real> Qij("Qfrom", S_max);
    var<Real> Pji("Pto", S_max);
    var<Real> Qji("Qto", S_max);
    
    ROMDST->add(Pij.in(Et),Qij.in(Et));
    DebugOff("size Pij = " << Pij.get_dim() << endl);
    if (pmt!=LDISTF) {
        ROMDST->add(Pji.in(Et),Qji.in(Et));
    }
    
    /** Voltage magnitude (squared) variables */
    var<Real> v("v", vmin, vmax);
    ROMDST->add(v.in(Nt));
    DebugOff("size v = " << v.get_dim() << endl);
//    v.initialize_all(0.64);
    
    /** Power loss variables */
    var<Real> loss("loss", pos_);
    if (pmt==DISTF || pmt==CDISTF) {
        ROMDST->add(loss.in(Nt));
    }
    
    /** OBJECTIVE FUNCTION */
    func_ obj = product(c1.in(exist_Gt), Pg_.in(exist_Gt)) + product(c1.in(pot_Gt), Pg_.in(pot_Gt)) + product(c2.in(exist_Gt), power(Pg_.in(exist_Gt),2)) + product(c2.in(pot_Gt), Pg2.in(pot_Gt)) + sum(c0.in(exist_Gt));
//    obj *= 12./months.size();
    obj += nT*product(c0.in(pot_gen),w_g.in(pot_gen));
    obj += product(gen_capcost.in(pot_gen), w_g.in(pot_gen));
    obj += product(inverter_capcost.in(pot_batt), w_b.in(pot_batt));
    obj += product(expansion_capcost.in(pot_edges), w_e.in(pot_edges));
    obj += product(pv_capcost.in(pot_pv), w_pv.in(pot_pv));
    obj += product(pv_varcost.in(pot_pv), Pv_cap.in(pot_pv));
    obj *= 1e-3;
    ROMDST->min(obj);
    
//    obj.print_expanded();
    
    /** CONSTRAINTS **/
    
/** FLOW CONSERVATION **/
    
    /** KCL P and Q */
    Constraint KCL_P("KCL_P");
    Constraint KCL_Q("KCL_Q");
    if (pmt==LDISTF) {
        KCL_P  = sum(Pij.out_arcs()) - sum(Pij.in_arcs()) + pl - sum(Pg.in_gens()) - sum(Pb.in_bats()) - sum(Pw.in_wind()) - sum(Pv.in_pv());
        //        KCL_P  = sum(Pij.out_arcs()) - sum(Pij.in_arcs()) + pl - sum(Pg.in_gens()) - sum(Pb.in_bats()) - sum(Pw.in_wind());
        KCL_Q  = sum(Qij.out_arcs()) - sum(Qij.in_arcs()) + ql - sum(Qg.in_gens()) - sum(Qb.in_bats());
    }
    else{
        KCL_P  = sum(Pij.out_arcs()) + sum(Pji.in_arcs()) + pl - sum(Pg.in_gens())- sum(Pb.in_bats()) - sum(Pw.in_wind()) - sum(Pv.in_pv());
        KCL_Q  = sum(Qij.out_arcs()) + sum(Qji.in_arcs()) + ql - sum(Qg.in_gens()) - sum(Qb.in_bats());
    }
    ROMDST->add(KCL_P.in(nodes, T) == 0);
    ROMDST->add(KCL_Q.in(nodes, T) == 0);
    
/**  THERMAL LIMITS **/
    
    /*  Thermal Limit Constraints */
    Constraint Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from += power(Pij, 2) + power(Qij, 2);
    Thermal_Limit_from -= power(S_max, 2);
    ROMDST->add(Thermal_Limit_from.in(exist_Et) <= 0);
    
    
    /*  Thermal Limit Constraints for expansion edges */
    Constraint Thermal_Limit_from_exp("Thermal_Limit_From_Exp");
    Thermal_Limit_from_exp += power(Pij, 2) + power(Qij, 2);
    Thermal_Limit_from_exp -= power(w_e,2)*power(S_max, 2);
    ROMDST->add(Thermal_Limit_from_exp.in(pot_Et) <= 0);
//    ROMDST->get_constraint("Thermal_Limit_From_Exp")->print_expanded();
    
/**  GENERATOR INVESTMENT **/
    
    /*  On/Off status */
    Constraint OnOff_maxP("OnOff_maxP");
    OnOff_maxP += Pg_ - pg_max*w_g;
    ROMDST->add(OnOff_maxP.in(pot_Gt) <= 0);
    
    Constraint Perspective_OnOff("Perspective_OnOff");
    Perspective_OnOff += power(Pg_,2) - Pg2*w_g;
    ROMDST->add(Perspective_OnOff.in(pot_Gt) <= 0);
    
    Constraint OnOff_maxQ("OnOff_maxQ");
    OnOff_maxQ += Qg - qg_max*w_g;
    ROMDST->add(OnOff_maxQ.in(pot_Gt) <= 0);

    Constraint OnOff_maxQ_N("OnOff_maxQ_N");
    OnOff_maxQ_N += Qg + qg_max*w_g;
    ROMDST->add(OnOff_maxQ_N.in(pot_Gt) >= 0);

/**  POWER FLOW **/
    
    /** Power Flow Constraints */
    Constraint Flow_P_From("Flow_P_From");
    Flow_P_From = v.from() - 2*(r*Pij + x*Qij) - v.to();
    ROMDST->add(Flow_P_From.in(_exist_arcs, T)==0);
    
    Constraint Flow_P_Expansion_L("Flow_P_Expansion_L");
    Flow_P_Expansion_L = v.from() - 2*(r*Pij + x*Qij) - v.to() - (1-w_e)*(v_diff_max);
    ROMDST->add(Flow_P_Expansion_L.in(_potential_expansion, T)<=0);
    
    Constraint Flow_P_Expansion_U("Flow_P_Expansion_U");
    Flow_P_Expansion_U = v.from() - 2*(r*Pij + x*Qij) - v.to() + (1-w_e)*(v_diff_max);
    ROMDST->add(Flow_P_Expansion_U.in(_potential_expansion, T)>=0);
    
    
/**  PV **/
    
    /*  On/Off on Potential PV */
    Constraint OnOffPV("OnOffPV");
    OnOffPV += Pv_cap - w_pv*pv_max;
    ROMDST->add(OnOffPV.in(indices(_potential_PV_gens)) <= 0);
//    ROMDST->get_constraint("OnOffPV")->print_expanded();
    
    /*  Max Cap on Potential PV */
    Constraint MaxCapPV("MaxCapPV");
    MaxCapPV += Pv - Pv_cap*pv_out;
    ROMDST->add(MaxCapPV.in(PV_pot_t) <= 0);
//    ROMDST->get_constraint("MaxCapPV")->print_expanded();
    
    /*  Existing PV */
    Constraint existPV("existPV");
    existPV += Pv - pv_max*pv_out;
    ROMDST->add(existPV.in(indices(_existing_PV_gens, T)) <= 0);
//    ROMDST->get_constraint("existPV")->print_expanded();
    

/**  BATTERIES **/
    
    /*  Apparent Power Limit on Potential Batteries */
    Constraint Apparent_Limit_Batt_Pot("Apparent_Limit_Batt_Potential");
    Apparent_Limit_Batt_Pot += power(Pb, 2) + power(Qb, 2);
    Apparent_Limit_Batt_Pot -= power(w_b,2)*power(pb_max, 2);
    ROMDST->add(Apparent_Limit_Batt_Pot.in(pot_Bt) <= 0);
    
    /*  Apparent Power Limit on Existing Batteries */
    Constraint Apparent_Limit_Batt("Apparent_Limit_Batt_Existing");
    Apparent_Limit_Batt += power(Pb, 2) + power(Qb, 2);
    Apparent_Limit_Batt -= power(pb_max, 2);
    ROMDST->add(Apparent_Limit_Batt.in(exist_Bt) <= 0);
    
    
    /*  State Of Charge */
    Constraint State_Of_Charge("State_Of_Charge");
    State_Of_Charge = Sc - Sc.prev() + Pb_;
    ROMDST->add(State_Of_Charge.in(Bt1) == 0);
    
    /*  State Of Charge 0 */
    auto Bat0 = indices(_battery_inverters,T.start());
    Constraint State_Of_Charge0("State_Of_Charge0");
    State_Of_Charge0 = Sc;
    ROMDST->add(State_Of_Charge0.in(Bat0) == 0);
    Constraint Pb0("Pb0");
    Pb0 = Pb_;
    ROMDST->add(Pb0.in(Bat0) == 0);
    
/*  EFFICIENCIES */
    Constraint DieselEff("DieselEff");
    DieselEff += Pg - gen_eff*Pg_;
    ROMDST->add(DieselEff.in(Gt) == 0);
    
    Constraint EfficiencyExist("BatteryEfficiencyExisting");
    EfficiencyExist += Pb  - eff_a*Pb_ - eff_b;//TODO without time extending eff_a and eff_b
    ROMDST->add(EfficiencyExist.in(indices(_eff_pieces,exist_Bt)) <= 0);
    
    Constraint EfficiencyPot("BatteryEfficiencyPotential");
    EfficiencyPot += Pb  - eff_a*Pb_ - eff_b*w_b;
    ROMDST->add(EfficiencyPot.in(indices(_eff_pieces,pot_Bt)) <= 0);
    
    
    for (auto n:nodes) {
        auto b = (Bus*)n;
//        b->print();
        for (auto i = 0; i < b->_pot_gen.size(); i++) {
            auto gen = b->_pot_gen[i];
            if(min_diesel_invest.eval(gen->_name)==max_diesel_invest.eval(gen->_name)){
                Constraint FixedDieselInvest("FixedDieselInvest"+gen->_name);
                FixedDieselInvest += w_g(gen->_name);
                ROMDST->add(FixedDieselInvest == 1);
                for (auto j = i+1; j < b->_pot_gen.size(); j++) {
                    auto gen2 = b->_pot_gen[j];
                    if (gen2->_gen_type==gen->_gen_type) {
                        Constraint FixedDieselInvest("FixedDieselInvest"+gen2->_name);
                        FixedDieselInvest += w_g(gen2->_name);
                        ROMDST->add(FixedDieselInvest == 1);
                    }
                }
            }
            else {
                Constraint MinDieselInvest("MinDieselInvest_"+b->_name+"_DG"+to_string(gen->_gen_type));
                MinDieselInvest += w_g(gen->_name);
                for (auto j = i+1; j < b->_pot_gen.size(); j++) {
                    auto gen2 = b->_pot_gen[j];
                    if (gen2->_gen_type==gen->_gen_type) {
                        MinDieselInvest += w_g(gen2->_name);
                    }
                }
                auto rhs = min_diesel_invest.eval(gen->_name);
                if (rhs>0) {
                    ROMDST->add(MinDieselInvest >= rhs);
                }
            }
        }
        for (auto i = 0; i < b->_pot_bat.size(); i++) {
            auto bat = b->_pot_bat[i];
            if(min_batt_invest.eval(bat->_name)==max_batt_invest.eval(bat->_name)){
                Constraint FixedBattInvest("FixedBattInvest"+bat->_name);
                FixedBattInvest += w_b(bat->_name);
                ROMDST->add(FixedBattInvest == 1);
                for (auto j = i+1; j < b->_pot_bat.size(); j++) {
                    auto bat2 = b->_pot_bat[j];
                    if (bat2->_bat_type==bat->_bat_type) {
                        Constraint FixedBattInvest("FixedBattInvest"+bat2->_name);
                        FixedBattInvest += w_b(bat2->_name);
                        ROMDST->add(FixedBattInvest == 1);
                    }
                }
            }
            else {
                Constraint MinBattInvest("MinBattInvest_"+b->_name+"_DG"+to_string(bat->_bat_type));
                MinBattInvest += w_b(bat->_name);
                for (auto j = i+1; j < b->_pot_bat.size(); j++) {
                    auto bat2 = b->_pot_bat[j];
                    if (bat2->_bat_type==bat->_bat_type) {
                        MinBattInvest += w_b(bat2->_name);
                    }
                }
                auto rhs = min_batt_invest.eval(bat->_name);
                if (rhs>0) {
                    ROMDST->add(MinBattInvest >= rhs);
                }
            }
        }
    }
//    ROMDST->print_expanded();
    return ROMDST;
}

unique_ptr<Model> PowerNet::build_fixed_ACOPF_N_1(PowerModelType pmt, int output, double tol, double obj_pen, const vector<indices>& ids_p, const vector<indices>& ids_n){
    bool polar = (pmt==ACPOL);
    if (polar) {
        Debug("Using polar model\n");
    }
    else {
        Debug("Using rectangular model\n");
    }
    unique_ptr<Model> ACOPF(new Model("AC-OPF Model"));
    /** Variables */
    /* Power generation variables */
    var<Real> Pg("Pg", pg_min, pg_max);
    var<Real> Qg ("Qg", qg_min, qg_max);
    ACOPF->add(Pg.in(gens), Qg.in(gens));
    Pg_base = Pg;
    Qg_base = Qg;
    
    
    /* Power flow variables */
    var<Real> Pf_from("Pf_from", S_max);
    var<Real> Qf_from("Qf_from", S_max);
    var<Real> Pf_to("Pf_to", S_max);
    var<Real> Qf_to("Qf_to", S_max);
    ACOPF->add(Pf_from.in(arcs), Qf_from.in(arcs),Pf_to.in(arcs),Qf_to.in(arcs));
    p_from_base = Pf_from;
    p_to_base = Pf_to;
    q_from_base = Qf_from;
    q_to_base = Qf_to;
    
    
    /** Voltage related variables */
    
    var<Real> theta("theta_base", -6.3, 6.3);
    var<Real> v("v_base", v_min, v_max);
    
    //    var<Real> vr("vr");
    //    var<Real> vi("vi");
    var<Real> vr("vr", v_max);
    var<Real> vi("vi", v_max);
    
    if (polar) {
        ACOPF->add(v.in(nodes), theta.in(nodes));
        v.initialize_all(1);
    }
    else {
        ACOPF->add_var(vr.in(nodes));
        ACOPF->add_var(vi.in(nodes));
        vr.initialize_all(1.0);
    }
    v_base = v;
    theta_base = theta;
    
    
    /** Construct the objective function */
    func_ obj = product(c1, Pg) + product(c2, power(Pg,2)) + sum(c0);
    ACOPF->min(obj.in(gens));
    
    /** Define constraints */
    
    
    /* REF BUS */
    Constraint Ref_Bus("Ref_Bus");
    if (polar) {
        Ref_Bus = theta(get_ref_bus());
    }
    else {
        Ref_Bus = vi(get_ref_bus());
    }
    ACOPF->add(Ref_Bus == 0);
    
    /** KCL Flow conservation */
    Constraint KCL_P("KCL_P");
    Constraint KCL_Q("KCL_Q");
    KCL_P  = sum(Pf_from.out_arcs()) + sum(Pf_to.in_arcs()) + pl - sum(Pg.in_gens());
    KCL_Q  = sum(Qf_from.out_arcs()) + sum(Qf_to.in_arcs()) + ql - sum(Qg.in_gens());
    /* Shunts */
    if (polar) {
        KCL_P +=  gs*power(v,2);
        KCL_Q -=  bs*power(v,2);
    }
    else {
        KCL_P +=  gs*(power(vr,2)+power(vi,2));
        KCL_Q -=  bs*(power(vr,2)+power(vi,2));
    }
    ACOPF->add(KCL_P.in(nodes) == 0);
    ACOPF->add(KCL_Q.in(nodes) == 0);
    
    
    /** AC Power Flows */
    /** TODO write the constraints in Complex form */
    Constraint Flow_P_From("Flow_P_From");
    Flow_P_From += Pf_from;
    if (polar) {
        Flow_P_From -= g/power(tr,2)*power(v.from(),2);
        Flow_P_From += g/tr*(v.from()*v.to()*cos(theta.from() - theta.to() - as));
        Flow_P_From += b/tr*(v.from()*v.to()*sin(theta.from() - theta.to() - as));
    }
    else {
        Flow_P_From -= g_ff*(power(vr.from(), 2) + power(vi.from(), 2));
        Flow_P_From -= g_ft*(vr.from()*vr.to() + vi.from()*vi.to());
        Flow_P_From -= b_ft*(vi.from()*vr.to() - vr.from()*vi.to());
    }
    ACOPF->add(Flow_P_From.in(arcs)==0);
    
    Constraint Flow_P_To("Flow_P_To");
    Flow_P_To += Pf_to;
    if (polar) {
        Flow_P_To -= g*power(v.to(), 2);
        Flow_P_To += g/tr*(v.from()*v.to()*cos(theta.to() - theta.from() + as));
        Flow_P_To += b/tr*(v.from()*v.to()*sin(theta.to() - theta.from() + as));
    }
    else {
        Flow_P_To -= g_tt*(power(vr.to(), 2) + power(vi.to(), 2));
        Flow_P_To -= g_tf*(vr.from()*vr.to() + vi.from()*vi.to());
        Flow_P_To -= b_tf*(vi.to()*vr.from() - vr.to()*vi.from());
    }
    ACOPF->add(Flow_P_To.in(arcs)==0);
    
    Constraint Flow_Q_From("Flow_Q_From");
    Flow_Q_From += Qf_from;
    if (polar) {
        Flow_Q_From += (0.5*ch+b)/power(tr,2)*power(v.from(),2);
        Flow_Q_From -= b/tr*(v.from()*v.to()*cos(theta.from() - theta.to() - as));
        Flow_Q_From += g/tr*(v.from()*v.to()*sin(theta.from() - theta.to() - as));
    }
    else {
        Flow_Q_From += b_ff*(power(vr.from(), 2) + power(vi.from(), 2));
        Flow_Q_From += b_ft*(vr.from()*vr.to() + vi.from()*vi.to());
        Flow_Q_From -= g_ft*(vi.from()*vr.to() - vr.from()*vi.to());
    }
    ACOPF->add(Flow_Q_From.in(arcs)==0);
    Constraint Flow_Q_To("Flow_Q_To");
    Flow_Q_To += Qf_to;
    if (polar) {
        Flow_Q_To += (0.5*ch+b)*power(v.to(),2);
        Flow_Q_To -= b/tr*(v.from()*v.to()*cos(theta.to() - theta.from() + as));
        Flow_Q_To += g/tr*(v.from()*v.to()*sin(theta.to() - theta.from() + as));
    }
    else {
        Flow_Q_To += b_tt*(power(vr.to(), 2) + power(vi.to(), 2));
        Flow_Q_To += b_tf*(vr.from()*vr.to() + vi.from()*vi.to());
        Flow_Q_To -= g_tf*(vi.to()*vr.from() - vr.to()*vi.from());
    }
    ACOPF->add(Flow_Q_To.in(arcs)==0);
    
    /** AC voltage limit constraints. */
    if (!polar) {
        Constraint Vol_limit_UB("Vol_limit_UB");
        Vol_limit_UB = power(vr, 2) + power(vi, 2);
        Vol_limit_UB -= power(v_max, 2);
        ACOPF->add(Vol_limit_UB.in(nodes) <= 0);
        
        Constraint Vol_limit_LB("Vol_limit_LB");
        Vol_limit_LB = power(vr, 2) + power(vi, 2);
        Vol_limit_LB -= power(v_min,2);
        ACOPF->add(Vol_limit_LB.in(nodes) >= 0);
        DebugOff(v_min.to_str(true) << endl);
        DebugOff(v_max.to_str(true) << endl);
    }
    
    
    /*  Thermal Limit Constraints */
    Constraint Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from += power(Pf_from, 2) + power(Qf_from, 2);
    Thermal_Limit_from -= power(S_max, 2);
    ACOPF->add(Thermal_Limit_from.in(arcs) <= 0);
    
    Constraint Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to += power(Pf_to, 2) + power(Qf_to, 2);
    Thermal_Limit_to -= power(S_max,2);
    ACOPF->add(Thermal_Limit_to.in(arcs) <= 0);
    //    ACOPF->print_expanded();
    //    ACOPF->_obj*=1e-4;
    
    unsigned nb_added = 0, nb_total = 0, nb_gens_cont = _conting_gens.size(), nb_lines_cont = _conting_lines.size();
    clog << "Building contingency models for lines." << endl;
    auto cont = 0;
    double start_built_time = get_wall_time();
    for (auto i = 0; i<nb_lines_cont; i++) {
        clog << "Switching Off Line: " << _conting_lines[i]->_name << endl;
        _conting_lines[i]->_active = false;
        if (i>0) {
            _conting_lines[i-1]->_active = true;
        }
        var<Real> Pg_c("Pg_cont"+to_string(cont), pg_min, pg_max);
        var<Real> Qg_c ("Qg_cont"+to_string(cont), qg_min, qg_max);
        //    var<Real> Sp ("Sp", pos_);
        //    var<Real> Sn ("Sn", pos_);
        ACOPF->add(Pg_c.in(gens), Qg_c.in(gens));
        
        
        Qg_cont[cont] = Qg_c;
        indices Ng; /*<< Nodes with generators not fixed to upper and lower bound */
        for (auto node: nodes) {
            auto bus = (Bus*)node;
            if (bus->_gen.empty()) {
                continue;
            }
            if (!ids_p[cont].empty() && ids_p[cont]._indices_map->count(bus->_name)!=0) {
                Qg_c.print(true);
                for (auto gen:bus->_gen) {
                    auto idx = qg_min.get_indices()->at(gen->_name);
                    Qg_c._ub->set_val(idx, qg_min.eval(idx));
                }
                Qg_c.print(true);
            }
            else if(!ids_n[cont].empty() && ids_n[cont]._indices_map->count(bus->_name)!=0){
                Qg_c.print(true);
                for (auto gen:bus->_gen) {
                    auto idx = qg_max.get_indices()->at(gen->_name);
                    Qg_c._lb->set_val(idx, qg_max.eval(idx));
                }
                Qg_c.print(true);
            }
            else {
                Ng.add(bus->_name);
            }
        }
        
        var<Real> Pf_from_c("Pf_from_c"+to_string(cont), S_max);
        var<Real> Qf_from_c("Qf_from_c"+to_string(cont), S_max);
        var<Real> Pf_to_c("Pf_to_c"+to_string(cont), S_max);
        var<Real> Qf_to_c("Qf_to_c"+to_string(cont), S_max);
        ACOPF->add(Pf_from_c.in(arcs), Qf_from_c.in(arcs),Pf_to_c.in(arcs),Qf_to_c.in(arcs));
        
        p_from[cont] = Pf_from_c;
        q_from[cont] = Qf_from_c;
        p_to[cont] = Pf_to_c;
        q_to[cont] = Qf_to_c;
        
        var<Real> delta("delta"+to_string(cont), pos_);
        ACOPF->add(delta.in(R(1)));
        _delta[cont] = delta;
        
        var<Real> theta_c("theta_c"+to_string(cont), -6.3, 6.3);
        var<Real> v_c("v_c"+to_string(cont), v_min, v_max);
        
        ACOPF->add(v_c.in(nodes), theta_c.in(nodes));
        v_c.initialize_all(1);
        v_cont[cont] = v_c;
        theta_cont[cont] = theta_c;
        
//        Constraint V_eq_base("V_eq_base"+to_string(cont));// v < v_c
//        V_eq_base = v_c - v_cont[cont];
//        ACOPF->add(V_eq_base.in(Ng) == 0);
//        
//        Constraint T_eq_base("T_eq_base"+to_string(cont));// v < v_c
//        T_eq_base = theta_c - theta_cont[cont];
//        ACOPF->add(T_eq_base.in(Ng) == 0);
        
        Constraint V_eq("V_eq"+to_string(cont));// v < v_c
        V_eq = (v_c - v);
        ACOPF->add(V_eq.in(Ng) == 0);
//        ACOPF->get_constraint("V_eq"+to_string(cont))->print_expanded();
        
        
        
        /* Generator Constraints */
        Constraint Gen_resp("g_resp"+to_string(cont));
        Gen_resp = Pg_c - (Pg + p_factor*delta);
        ACOPF->add(Gen_resp.in(indices(gens)) == 0);
        /** KCL Flow conservation */
        Constraint KCL_P_c("KCL_P_c_"+to_string(cont));
        Constraint KCL_Q_c("KCL_Q_c"+to_string(cont));
        KCL_P_c  = sum(Pf_from_c.out_arcs()) + sum(Pf_to_c.in_arcs()) + pl - sum(Pg_c.in_gens());
        KCL_Q_c  = sum(Qf_from_c.out_arcs()) + sum(Qf_to_c.in_arcs()) + ql - sum(Qg_c.in_gens());
        /* Shunts */
        if (polar) {
            KCL_P_c +=  gs*power(v_c,2);
            KCL_Q_c -=  bs*power(v_c,2);
        }
        ACOPF->add(KCL_P_c.in(nodes) == 0);
        ACOPF->add(KCL_Q_c.in(nodes) == 0);
        
        Constraint Flow_P_From_c("Flow_P_From_c"+to_string(cont));
        Flow_P_From_c += Pf_from_c;
        Flow_P_From_c -= g/power(tr,2)*power(v_c.from(),2);
        Flow_P_From_c += g/tr*(v_c.from()*v_c.to()*cos(theta_c.from() - theta_c.to() - as));
        Flow_P_From_c += b/tr*(v_c.from()*v_c.to()*sin(theta_c.from() - theta_c.to() - as));
        ACOPF->add(Flow_P_From_c.in(arcs)==0);
        
        Constraint Flow_P_To_c("Flow_P_To_c"+to_string(cont));
        Flow_P_To_c += Pf_to_c;
        Flow_P_To_c -= g*power(v_c.to(), 2);
        Flow_P_To_c += g/tr*(v_c.from()*v_c.to()*cos(theta_c.to() - theta_c.from() + as));
        Flow_P_To_c += b/tr*(v_c.from()*v_c.to()*sin(theta_c.to() - theta_c.from() + as));
        ACOPF->add(Flow_P_To_c.in(arcs)==0);
        
        Constraint Flow_Q_From_c("Flow_Q_From_c"+to_string(cont));
        Flow_Q_From_c += Qf_from_c;
        Flow_Q_From_c += (0.5*ch+b)/power(tr,2)*power(v_c.from(),2);
        Flow_Q_From_c -= b/tr*(v_c.from()*v_c.to()*cos(theta_c.from() - theta_c.to() - as));
        Flow_Q_From_c += g/tr*(v_c.from()*v_c.to()*sin(theta_c.from() - theta_c.to() - as));
        ACOPF->add(Flow_Q_From_c.in(arcs)==0);
        Constraint Flow_Q_To_c("Flow_Q_To_c"+to_string(cont));
        Flow_Q_To_c += Qf_to_c;
        Flow_Q_To_c += (0.5*ch+b)*power(v_c.to(),2);
        Flow_Q_To_c -= b/tr*(v_c.from()*v_c.to()*cos(theta_c.to() - theta_c.from() + as));
        Flow_Q_To_c += g/tr*(v_c.from()*v_c.to()*sin(theta_c.to() - theta_c.from() + as));
        ACOPF->add(Flow_Q_To_c.in(arcs)==0);
        
        /*  Thermal Limit Constraints */
        Constraint Thermal_Limit_from_c("Thermal_Limit_from_c"+to_string(cont));
        Thermal_Limit_from_c += power(Pf_from_c, 2) + power(Qf_from_c, 2);
        Thermal_Limit_from_c -= power(S_max, 2);
        ACOPF->add(Thermal_Limit_from_c.in(arcs) <= 0);
        
        Constraint Thermal_Limit_to_c("Thermal_Limit_to_c"+to_string(cont));
        Thermal_Limit_to_c += power(Pf_to_c, 2) + power(Qf_to_c, 2);
        Thermal_Limit_to_c -= power(S_max,2);
        ACOPF->add(Thermal_Limit_to_c.in(arcs) <= 0);
        cont++;
    }
    if (nb_gens_cont>0) {
        clog << "Building contingency models for generators." << endl;
    }
    else {
        clog << "No generator contingencies." << endl;
    }
    for (auto i = 0; i<nb_gens_cont; i++) {
        clog << "Deactivating Generator: " << _conting_gens[i]->_name << endl;
        _conting_gens[i]->_active = false;
        if (i>0) {
            _conting_gens[i-1]->_active = true;
        }
        cont++;
    }
    nb_total = nb_gens_cont + nb_lines_cont;
    double end_build_time = get_wall_time();
    auto model_build_time = end_build_time - start_built_time;
    //    ACOPF->print_expanded();
    DebugOn("Total time for building continency constraints = " << model_build_time << endl);
    DebugOn("nb_vars = " << ACOPF->get_nb_vars()<<endl);
    DebugOn("nb_eqs = " << ACOPF->get_nb_cons() - ACOPF->get_nb_ineq()<<endl);
    DebugOn("nb_ineqs = " << ACOPF->get_nb_ineq()<<endl);
    
    return ACOPF;
}

vector<param<>> PowerNet::signs() {
    vector<param<>> res;
    string key;
    size_t idx;
    res.resize(3);
    for (int i = 0; i<3; i++) {
        res[i].set_name("I_sign_"+to_string(i));
    }
    set<vector<unsigned>> ids;
    for (auto &bag: _bags) {
        if (bag.size() != 3) {
            continue;
        }
        vector<unsigned> ids_bag;
        for (int i = 0; i<3; i++) {
            ids_bag.push_back(bag[i]->_id);
        }
        if(ids.count(ids_bag)==0) {
            ids.insert(ids_bag);
        } else {
            continue;
        }
        for (int i = 0; i< 2; i++) {
            if(has_directed_arc(bag[i], bag[i+1])) {
                key = bag[i]->_name + "," + bag[i+1]->_name;
                idx = res[i].set_val(key,1.0);
                res[i]._ids->at(0).push_back(idx);
            }
            else {
                key = bag[i+1]->_name + "," + bag[i]->_name;
                DebugOff("\nreversed arc " << key);
                idx = res[i].set_val(key,-1.0);
                res[i]._ids->at(0).push_back(idx);
            }
        }
        /* Loop back pair */
        if(has_directed_arc(bag[0], bag[2])) {
            key = bag[0]->_name + "," + bag[2]->_name;
            idx = res[2].set_val(key,1.0);
            res[2]._ids->at(0).push_back(idx);
        }
        else{
            key = bag[2]->_name + "," + bag[0]->_name;
            DebugOff("\nreversed arc " << key);
            idx = res[2].set_val(key,-1.0);
            res[2]._ids->at(0).push_back(idx);
        }
    }
    for (int i = 0; i<3; i++) {
        res[i]._dim[0]=res[i]._ids->at(0).size();
        res[i]._is_indexed = true;
    }
    return res;
}

unique_ptr<Model> PowerNet::build_SOCP_OPF_N_1(PowerModelType pmt, int output, double tol, bool sdp_cuts){
    bool socp = (pmt==SOCP);
    if (socp) {
        Debug("Using SOCP model\n");
    }
    else {
        Debug("Using QC model\n");
    }
    unique_ptr<Model> SOCP(new Model("SOCP Model"));
    /** Variables */
    /* Power generation variables */
    var<Real> Pg("Pg", pg_min, pg_max);
    var<Real> Qg ("Qg", qg_min, qg_max);
    SOCP->add(Pg.in(gens), Qg.in(gens));
    Pg_base = Pg;
    Qg_base = Qg;
    
    
    /* Power flow variables */
    var<Real> Pf_from("Pf_from", S_max);
    var<Real> Qf_from("Qf_from", S_max);
    var<Real> Pf_to("Pf_to", S_max);
    var<Real> Qf_to("Qf_to", S_max);
    SOCP->add(Pf_from.in(arcs), Qf_from.in(arcs),Pf_to.in(arcs),Qf_to.in(arcs));
    p_from_base = Pf_from;
    p_to_base = Pf_to;
    q_from_base = Qf_from;
    q_to_base = Qf_to;
    
    indices Ng; /*<< Nodes with generators */
    for (auto node: nodes) {
        auto bus = (Bus*)node;
        if (bus->_gen.size()>0) {
            Ng.add(bus->_name);
        }
    }
    
    /** Voltage related variables */
    
    /* Real part of Wij = ViVj */
    var<Real>  R_Wij("R_Wij", wr_min, wr_max);
    /* Imaginary part of Wij = ViVj */
    var<Real>  Im_Wij("Im_Wij", wi_min, wi_max);
    /* Magnitude of Wii = Vi^2 */
    var<Real>  Wii("Wii", w_min, w_max);
    SOCP->add_var(Wii.in(nodes));
    auto bus_pairs = get_bus_pairs();
    auto bus_pairs_chord = get_bus_pairs_chord();
    if (sdp_cuts) {
        SOCP->add_var(R_Wij.in(bus_pairs_chord));
        SOCP->add_var(Im_Wij.in(bus_pairs_chord));
    }
    else {
        SOCP->add_var(R_Wij.in(bus_pairs));
        SOCP->add_var(Im_Wij.in(bus_pairs));
    }
    
    /* Initialize variables */
    R_Wij.initialize_all(1.0);
    Wii.initialize_all(1.001);
    
    
    if (sdp_cuts) {
//        auto R_Wij_ = R_Wij.pairs_in_directed(*this, _bags, 3);
//        auto Im_Wij_ = Im_Wij.pairs_in_directed(*this, _bags, 3);
//        auto Wii_ = Wii.in(_bags, 3);
//        auto I_sgn = signs();
//        Constraint SDP3("SDP_3D");
//        SDP3 = 2 * R_Wij_[0] * (R_Wij_[1] * R_Wij_[2] + I_sgn[1] * I_sgn[2] * Im_Wij_[1] * Im_Wij_[2]);
//        SDP3 += 2 * I_sgn[0] * Im_Wij_[0] * (R_Wij_[1] * I_sgn[2] * Im_Wij_[2] - I_sgn[1] * Im_Wij_[1] * R_Wij_[2]);
//        SDP3 -= (power(R_Wij_[0], 2) + power(Im_Wij_[0], 2)) * Wii_[2];
//        SDP3 -= (power(R_Wij_[1], 2) + power(Im_Wij_[1], 2)) * Wii_[0];
//        SDP3 -= (power(R_Wij_[2], 2) + power(Im_Wij_[2], 2)) * Wii_[1];
//        SDP3 += Wii_[0] * Wii_[1] * Wii_[2];
//        DebugOn("\nsdp nb inst = " << SDP3.get_nb_instances() << endl);
//        SOCP->add(SDP3 >= 0);
    }
    
    /* Second-order cone constraints */
    Constraint SOC("SOC");
    SOC = power(R_Wij, 2) + power(Im_Wij, 2) - Wii.from()*Wii.to();
    if (sdp_cuts) {
        SOCP->add(SOC.in(bus_pairs_chord) <= 0);
    }
    else {
        SOCP->add(SOC.in(bus_pairs) <= 0);
    }
    
    
    
    /** Construct the objective function */
    func_ obj = product(c1, Pg) + product(c2, power(Pg,2)) + sum(c0);
    SOCP->min(obj.in(gens));
    
    /** Define constraints */
    
    
    
    /** KCL Flow conservation */
    Constraint KCL_P("KCL_P");
    Constraint KCL_Q("KCL_Q");
    KCL_P  = sum(Pf_from.out_arcs()) + sum(Pf_to.in_arcs()) + pl - sum(Pg.in_gens());
    KCL_Q  = sum(Qf_from.out_arcs()) + sum(Qf_to.in_arcs()) + ql - sum(Qg.in_gens());
    /* Shunts */
    KCL_P +=  gs*Wii;
    KCL_Q -=  bs*Wii;
    
    SOCP->add(KCL_P.in(nodes) == 0);
    SOCP->add(KCL_Q.in(nodes) == 0);
    
    
    /** AC Power Flows */
    /** TODO write the constraints in Complex form */
    Constraint Flow_P_From("Flow_P_From");
    Flow_P_From = Pf_from - (g_ff*Wii.from() + g_ft*R_Wij.in_pairs() + b_ft*Im_Wij.in_pairs());
    SOCP->add(Flow_P_From.in(arcs)==0);
    
    Constraint Flow_P_To("Flow_P_To");
    Flow_P_To = Pf_to - (g_tt*Wii.to() + g_tf*R_Wij.in_pairs() - b_tf*Im_Wij.in_pairs());
    SOCP->add(Flow_P_To.in(arcs)==0);
    
    Constraint Flow_Q_From("Flow_Q_From");
    Flow_Q_From = Qf_from - (g_ft*Im_Wij.in_pairs() - b_ff*Wii.from() - b_ft*R_Wij.in_pairs());
    SOCP->add(Flow_Q_From.in(arcs) == 0);
    
    Constraint Flow_Q_To("Flow_Q_To");
    Flow_Q_To = Qf_to + (b_tt*Wii.to() + b_tf*R_Wij.in_pairs() + g_tf*Im_Wij.in_pairs());
    SOCP->add(Flow_Q_To.in(arcs) == 0);
    
    
    
    /*  Thermal Limit Constraints */
    Constraint Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from += power(Pf_from, 2) + power(Qf_from, 2);
    Thermal_Limit_from -= power(S_max, 2);
    SOCP->add(Thermal_Limit_from.in(arcs) <= 0);
    
    Constraint Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to += power(Pf_to, 2) + power(Qf_to, 2);
    Thermal_Limit_to -= power(S_max,2);
    SOCP->add(Thermal_Limit_to.in(arcs) <= 0);
    //    ACOPF->print_expanded();
    //    ACOPF->_obj*=1e-4;
    
    unsigned nb_added = 0, nb_total = 0, nb_gens_cont = _conting_gens.size(), nb_lines_cont = _conting_lines.size();
    clog << "Building contingency models for lines." << endl;
    auto cont = 0;
    double start_built_time = get_wall_time();
    for (auto i = 0; i<nb_lines_cont; i++) {
        clog << "Switching Off Line: " << _conting_lines[i]->_name << endl;
        _conting_lines[i]->_active = false;
        if (i>0) {
            _conting_lines[i-1]->_active = true;
        }
        var<Real> Pg_c("Pg_cont"+to_string(cont), pg_min, pg_max);
        var<Real> Qg_c ("Qg_cont"+to_string(cont), qg_min, qg_max);
        //    var<Real> Sp ("Sp", pos_);
        //    var<Real> Sn ("Sn", pos_);
        SOCP->add(Pg_c.in(gens), Qg_c.in(gens));
        Qg_cont[cont] = Qg_c;
        
        
        var<Real> Pf_from_c("Pf_from_c"+to_string(cont), S_max);
        var<Real> Qf_from_c("Qf_from_c"+to_string(cont), S_max);
        var<Real> Pf_to_c("Pf_to_c"+to_string(cont), S_max);
        var<Real> Qf_to_c("Qf_to_c"+to_string(cont), S_max);
        SOCP->add(Pf_from_c.in(arcs), Qf_from_c.in(arcs),Pf_to_c.in(arcs),Qf_to_c.in(arcs));
        
        p_from[cont] = Pf_from_c;
        q_from[cont] = Qf_from_c;
        p_to[cont] = Pf_to_c;
        q_to[cont] = Qf_to_c;
        
        var<Real> delta("delta"+to_string(cont), pos_);
        SOCP->add(delta.in(R(1)));
        _delta[cont] = delta;
        /* Real part of Wij = ViVj */
        var<Real>  R_Wij_c("R_Wij"+to_string(cont), wr_min, wr_max);
        /* Imaginary part of Wij = ViVj */
        var<Real>  Im_Wij_c("Im_Wij"+to_string(cont), wi_min, wi_max);
        /* Magnitude of Wii = Vi^2 */
        var<Real>  Wii_c("Wii"+to_string(cont), w_min, w_max);
        SOCP->add_var(Wii_c.in(nodes));
        SOCP->add_var(R_Wij_c.in(bus_pairs));
        SOCP->add_var(Im_Wij_c.in(bus_pairs));
        
        /* Initialize variables */
        R_Wij_c.initialize_all(1.0);
        Wii_c.initialize_all(1.001);
        
        /* Second-order cone constraints */
        Constraint SOC_c("SOC"+to_string(cont));
        SOC_c = power(R_Wij_c, 2) + power(Im_Wij_c, 2) - Wii_c.from()*Wii_c.to();
        SOCP->add(SOC_c.in(bus_pairs) <= 0);
        
        //        auto R_Wij_c_ = R_Wij_c.pairs_in_directed(*this, _bags, 3);
        //        auto Im_Wij_c_ = Im_Wij_c.pairs_in_directed(*this, _bags, 3);
        //        auto Wii_c_ = Wii_c.in_bags(_bags, 3);
        //        Constraint SDP3_c("SDP_3D"+to_string(cont));
        //        SDP3_c = 2 * R_Wij_c_[0] * (R_Wij_c_[1] * R_Wij_c_[2] + I_sgn[1] * I_sgn[2] * Im_Wij_c_[1] * Im_Wij_c_[2]);
        //        SDP3_c += 2 * I_sgn[0] * Im_Wij_c_[0] * (R_Wij_c_[1] * I_sgn[2] * Im_Wij_c_[2] - I_sgn[1] * Im_Wij_c_[1] * R_Wij_c_[2]);
        //        SDP3_c -= (power(R_Wij_c_[0], 2) + power(Im_Wij_c_[0], 2)) * Wii_c_[2];
        //        SDP3_c -= (power(R_Wij_c_[1], 2) + power(Im_Wij_c_[1], 2)) * Wii_c_[0];
        //        SDP3_c -= (power(R_Wij_c_[2], 2) + power(Im_Wij_c_[2], 2)) * Wii_c_[1];
        //        SDP3_c += Wii_c_[0] * Wii_c_[1] * Wii_c_[2];
        ////        DebugOn("\nsdp nb inst = " << SDP3.get_nb_instances() << endl);
        //        SOCP->add(SDP3_c >= 0);
        
        Constraint V_eq("V_eq"+to_string(cont));// v < v_c
        V_eq = (Wii_c - Wii);
        SOCP->add(V_eq.in(Ng) == 0);
        
        
        
        /* Generator Constraints */
        Constraint Gen_resp("g_resp"+to_string(cont));
        Gen_resp = Pg_c - (Pg + p_factor*delta);
        SOCP->add(Gen_resp.in(indices(gens)) == 0);
        
        /** KCL Flow conservation */
        Constraint KCL_P_c("KCL_P_c_"+to_string(cont));
        Constraint KCL_Q_c("KCL_Q_c"+to_string(cont));
        KCL_P_c  = sum(Pf_from_c.out_arcs()) + sum(Pf_to_c.in_arcs()) + pl - sum(Pg_c.in_gens());
        KCL_Q_c  = sum(Qf_from_c.out_arcs()) + sum(Qf_to_c.in_arcs()) + ql - sum(Qg_c.in_gens());
        /* Shunts */
        KCL_P_c +=  gs*Wii_c;
        KCL_Q_c -=  bs*Wii_c;
        
        
        Constraint Flow_P_From_c("Flow_P_From_c"+to_string(cont));
        Flow_P_From_c = Pf_from_c - (g_ff*Wii_c.from() + g_ft*R_Wij_c.in_pairs() + b_ft*Im_Wij_c.in_pairs());
        SOCP->add(Flow_P_From_c.in(arcs)==0);
        
        
        
        Constraint Flow_P_To_c("Flow_P_To_c"+to_string(cont));
        Flow_P_To_c = Pf_to_c - (g_tt*Wii_c.to() + g_tf*R_Wij_c.in_pairs() - b_tf*Im_Wij_c.in_pairs());
        SOCP->add(Flow_P_To_c.in(arcs)==0);
        
        Constraint Flow_Q_From_c("Flow_Q_From_c"+to_string(cont));
        Flow_Q_From_c = Qf_from_c - (g_ft*Im_Wij_c.in_pairs() - b_ff*Wii_c.from() - b_ft*R_Wij_c.in_pairs());
        SOCP->add(Flow_Q_From_c.in(arcs)==0);
        
        Constraint Flow_Q_To_c("Flow_Q_To_c"+to_string(cont));
        Flow_Q_To_c = Qf_to_c + (b_tt*Wii_c.to() + b_tf*R_Wij_c.in_pairs() + g_tf*Im_Wij_c.in_pairs());
        SOCP->add(Flow_Q_To_c.in(arcs)==0);
        
        /*  Thermal Limit Constraints */
        Constraint Thermal_Limit_from_c("Thermal_Limit_from_c"+to_string(cont));
        Thermal_Limit_from_c += power(Pf_from_c, 2) + power(Qf_from_c, 2);
        Thermal_Limit_from_c -= power(S_max, 2);
        SOCP->add(Thermal_Limit_from_c.in(arcs) <= 0);
        
        Constraint Thermal_Limit_to_c("Thermal_Limit_to_c"+to_string(cont));
        Thermal_Limit_to_c += power(Pf_to_c, 2) + power(Qf_to_c, 2);
        Thermal_Limit_to_c -= power(S_max,2);
        SOCP->add(Thermal_Limit_to_c.in(arcs) <= 0);
        cont++;
    }
    _conting_lines.back()->_active = true;
    if (nb_gens_cont>0) {
        clog << "Building contingency models for generators." << endl;
    }
    else {
        clog << "No generator contingencies." << endl;
    }
    for (auto i = 0; i<nb_gens_cont; i++) {
        clog << "Deactivating Generator: " << _conting_gens[i]->_name << endl;
        _conting_gens[i]->_active = false;
        if (i>0) {
            _conting_gens[i-1]->_active = true;
        }
        cont++;
    }
    
    
    
    
    nb_total = nb_gens_cont + nb_lines_cont;
    double end_build_time = get_wall_time();
    auto model_build_time = end_build_time - start_built_time;
    //    ACOPF->print_expanded();
    DebugOn("Total time for building continency constraints = " << model_build_time << endl);
    DebugOn("nb_vars = " << SOCP->get_nb_vars()<<endl);
    DebugOn("nb_eqs = " << SOCP->get_nb_cons() - SOCP->get_nb_ineq()<<endl);
    DebugOn("nb_ineqs = " << SOCP->get_nb_ineq()<<endl);
    
    return SOCP;
}

unique_ptr<Model> PowerNet::build_SOCP_OPF_MINLP(PowerModelType pmt, int output, double tol){
    bool socp = (pmt==SOCP);
    if (socp) {
        Debug("Using SOCP model\n");
    }
    else {
        Debug("Using QC model\n");
    }
    unique_ptr<Model> SOCP(new Model("SOCP Model"));
    /** Variables */
    /* Power generation variables */
    var<Real> Pg("Pg", pg_min, pg_max);
    var<Real> Qg ("Qg", qg_min, qg_max);
    SOCP->add(Pg.in(gens), Qg.in(gens));
    Pg_base = Pg;
    Qg_base = Qg;
    
    
    /* Power flow variables */
    var<Real> Pf_from("Pf_from", S_max);
    var<Real> Qf_from("Qf_from", S_max);
    var<Real> Pf_to("Pf_to", S_max);
    var<Real> Qf_to("Qf_to", S_max);
    SOCP->add(Pf_from.in(arcs), Qf_from.in(arcs),Pf_to.in(arcs),Qf_to.in(arcs));
    p_from_base = Pf_from;
    p_to_base = Pf_to;
    q_from_base = Qf_from;
    q_to_base = Qf_to;
    
    indices Ng; /*<< Nodes with generators */
    for (auto node: nodes) {
        auto bus = (Bus*)node;
        if (bus->_gen.size()>0) {
            Ng.add(bus->_name);
        }
    }
    
    /** Voltage related variables */
    
    
    /* Real part of Wij = ViVj */
    var<Real>  R_Wij("R_Wij", wr_min, wr_max);
    /* Imaginary part of Wij = ViVj */
    var<Real>  Im_Wij("Im_Wij", wi_min, wi_max);
    /* Magnitude of Wii = Vi^2 */
    var<Real>  Wii("Wii", w_min, w_max);
    SOCP->add_var(Wii.in(nodes));
    auto bus_pairs = get_bus_pairs();
    auto bus_pairs_chord = get_bus_pairs_chord();
    SOCP->add_var(R_Wij.in(bus_pairs));
    SOCP->add_var(Im_Wij.in(bus_pairs));
    
    /* Initialize variables */
    R_Wij.initialize_all(1.0);
    Wii.initialize_all(1.001);
    
    
    
    /** Define constraints */
    
    /* Second-order cone constraints */
    Constraint SOC("SOC");
    SOC = power(R_Wij, 2) + power(Im_Wij, 2) - Wii.from()*Wii.to();
    SOCP->add(SOC.in(bus_pairs) <= 0);
    
    
    
    /** Construct the objective function */
    func_ obj = product(c1, Pg) + product(c2, power(Pg,2)) + sum(c0);
    SOCP->min(obj.in(gens));
    
    
    
    
    
    /** KCL Flow conservation */
    Constraint KCL_P("KCL_P");
    Constraint KCL_Q("KCL_Q");
    KCL_P  = sum(Pf_from.out_arcs()) + sum(Pf_to.in_arcs()) + pl - sum(Pg.in_gens());
    KCL_Q  = sum(Qf_from.out_arcs()) + sum(Qf_to.in_arcs()) + ql - sum(Qg.in_gens());
    /* Shunts */
    KCL_P +=  gs*Wii;
    KCL_Q -=  bs*Wii;
    
    SOCP->add(KCL_P.in(nodes) == 0);
    SOCP->add(KCL_Q.in(nodes) == 0);
    
    
    /** AC Power Flows */
    /** TODO write the constraints in Complex form */
    Constraint Flow_P_From("Flow_P_From");
    Flow_P_From = Pf_from - (g_ff*Wii.from() + g_ft*R_Wij.in_pairs() + b_ft*Im_Wij.in_pairs());
    SOCP->add(Flow_P_From.in(arcs)==0);
    
    Constraint Flow_P_To("Flow_P_To");
    Flow_P_To = Pf_to - (g_tt*Wii.to() + g_tf*R_Wij.in_pairs() - b_tf*Im_Wij.in_pairs());
    SOCP->add(Flow_P_To.in(arcs)==0);
    
    Constraint Flow_Q_From("Flow_Q_From");
    Flow_Q_From = Qf_from - (g_ft*Im_Wij.in_pairs() - b_ff*Wii.from() - b_ft*R_Wij.in_pairs());
    SOCP->add(Flow_Q_From.in(arcs) == 0);
    
    Constraint Flow_Q_To("Flow_Q_To");
    Flow_Q_To = Qf_to + (b_tt*Wii.to() + b_tf*R_Wij.in_pairs() + g_tf*Im_Wij.in_pairs());
    SOCP->add(Flow_Q_To.in(arcs) == 0);
    
    
    
    /*  Thermal Limit Constraints */
    Constraint Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from += power(Pf_from, 2) + power(Qf_from, 2);
    Thermal_Limit_from -= power(S_max, 2);
    SOCP->add(Thermal_Limit_from.in(arcs) <= 0);
    
    Constraint Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to += power(Pf_to, 2) + power(Qf_to, 2);
    Thermal_Limit_to -= power(S_max,2);
    SOCP->add(Thermal_Limit_to.in(arcs) <= 0);
    //    ACOPF->print_expanded();
    //    ACOPF->_obj*=1e-4;
    
    unsigned nb_added = 0, nb_total = 0, nb_gens_cont = _conting_gens.size(), nb_lines_cont = _conting_lines.size();
    clog << "Building contingency models for lines." << endl;
    auto cont = 0;
    double start_built_time = get_wall_time();
    for (auto i = 0; i<nb_lines_cont; i++) {
        clog << "Switching Off Line: " << _conting_lines[i]->_name << endl;
        _conting_lines[i]->_active = false;
        if (i>0) {
            _conting_lines[i-1]->_active = true;
        }
        var<Real> Pg_c("Pg_cont"+to_string(cont), pg_min, pg_max);
        var<Real> Qg_c ("Qg_cont"+to_string(cont), qg_min, qg_max);
        //    var<Real> Sp ("Sp", pos_);
        //    var<Real> Sn ("Sn", pos_);
        SOCP->add(Pg_c.in(gens), Qg_c.in(gens));
        Qg_cont[cont] = Qg_c;
        
        
        var<Real> Pf_from_c("Pf_from_c"+to_string(cont), S_max);
        var<Real> Qf_from_c("Qf_from_c"+to_string(cont), S_max);
        var<Real> Pf_to_c("Pf_to_c"+to_string(cont), S_max);
        var<Real> Qf_to_c("Qf_to_c"+to_string(cont), S_max);
        SOCP->add(Pf_from_c.in(arcs), Qf_from_c.in(arcs),Pf_to_c.in(arcs),Qf_to_c.in(arcs));
        
        p_from[cont] = Pf_from_c;
        q_from[cont] = Qf_from_c;
        p_to[cont] = Pf_to_c;
        q_to[cont] = Qf_to_c;
        
        var<Real> delta("delta"+to_string(cont), pos_);
        SOCP->add(delta.in(R(1)));
        _delta[cont] = delta;
        /* Real part of Wij = ViVj */
        var<Real>  R_Wij_c("R_Wij"+to_string(cont), wr_min, wr_max);
        /* Imaginary part of Wij = ViVj */
        var<Real>  Im_Wij_c("Im_Wij"+to_string(cont), wi_min, wi_max);
        /* Magnitude of Wii = Vi^2 */
        var<Real>  Wii_c("Wii"+to_string(cont), w_min, w_max);
        SOCP->add_var(Wii_c.in(nodes));
        SOCP->add_var(R_Wij_c.in(bus_pairs));
        SOCP->add_var(Im_Wij_c.in(bus_pairs));
        
        /* Initialize variables */
        R_Wij_c.initialize_all(1.0);
        Wii_c.initialize_all(1.001);
        
        
        //        auto R_Wij_c_ = R_Wij_c.pairs_in_directed(*this, _bags, 3);
        //        auto Im_Wij_c_ = Im_Wij_c.pairs_in_directed(*this, _bags, 3);
        //        auto Wii_c_ = Wii_c.in_bags(_bags, 3);
        //        Constraint SDP3_c("SDP_3D"+to_string(cont));
        //        SDP3_c = 2 * R_Wij_c_[0] * (R_Wij_c_[1] * R_Wij_c_[2] + I_sgn[1] * I_sgn[2] * Im_Wij_c_[1] * Im_Wij_c_[2]);
        //        SDP3_c += 2 * I_sgn[0] * Im_Wij_c_[0] * (R_Wij_c_[1] * I_sgn[2] * Im_Wij_c_[2] - I_sgn[1] * Im_Wij_c_[1] * R_Wij_c_[2]);
        //        SDP3_c -= (power(R_Wij_c_[0], 2) + power(Im_Wij_c_[0], 2)) * Wii_c_[2];
        //        SDP3_c -= (power(R_Wij_c_[1], 2) + power(Im_Wij_c_[1], 2)) * Wii_c_[0];
        //        SDP3_c -= (power(R_Wij_c_[2], 2) + power(Im_Wij_c_[2], 2)) * Wii_c_[1];
        //        SDP3_c += Wii_c_[0] * Wii_c_[1] * Wii_c_[2];
        ////        DebugOn("\nsdp nb inst = " << SDP3.get_nb_instances() << endl);
        //        SOCP->add(SDP3_c >= 0);
        /* Binary vars for complementarity, on = 1 <=> v_c = v */
        var<bool> on("on"+to_string(cont));
        SOCP->add_var(on.in(Ng));
        /* Binary vars for complementarity, low = 1 <=> qg_c = qg_min forall gens attached to n */
        var<bool> low("low"+to_string(cont));
        SOCP->add_var(low.in(Ng));
        /* Binary vars for complementarity, up = 1 <=> qg_c = qg_max forall gens attached to n */
        var<bool> up("up"+to_string(cont));
        SOCP->add_var(up.in(Ng));
        
        SOCP->_obj -= sum(on);
        
        /* logic constraints */
        Constraint Disj("Disjunction"+to_string(cont));
        Disj = low + up + on;
        SOCP->add(Disj.in(Ng) == 1);
        
        Constraint On_c_u("On_c_u"+to_string(cont));
        On_c_u =  Wii - Wii_c - (1-on)*(w_max - w_min);
        SOCP->add(On_c_u.in(Ng) <= 0);
        
        Constraint On_c_l("On_c_l"+to_string(cont));
        On_c_l =  Wii - Wii_c - (1-on)*(w_min - w_max);
        SOCP->add(On_c_l.in(Ng) >= 0);
        
        //        Constraint Fix_on("Fix_on"+to_string(cont));
        //        Fix_on =  on;
        //        SOCP->add(Fix_on.in(Ng) == 1);
        for (auto n:nodes) {
            auto bus = (Bus*)n;
            if (bus->_gen.empty()) {
                continue;
            }
            for (auto gen: bus->_gen) {
                Constraint Low_c("Low_c_"+gen->_name+"_"+to_string(cont));
                auto gname = gen->_name;
                Low_c =  Qg_c(gname) - qg_min(gname) - (1-low(bus->_name))*(qg_max(gname) - qg_min(gname));
                SOCP->add(Low_c <= 0);
                
                Constraint Up_c("Up_c"+gen->_name+"_"+to_string(cont));
                Up_c =  Qg_c(gname) - qg_max(gname) - (1-up(bus->_name))*(qg_min(gname) - qg_max(gname));
                SOCP->add(Up_c >= 0);
            }
        }
        //
        //        Constraint Low_c_l("Low_c_l"+to_string(cont));
        //        Low_c_l =  Wii - Wii_c - (1-on)*(w_min - w_max);
        //        SOCP->add(Low_c_l.in(Ng) >= 0);
        
        //        Constraint V_eq("V_eq"+to_string(cont));// v < v_c
        //        V_eq = (Wii_c - Wii);
        //        SOCP->add(V_eq.in(Ng) == 0);
        
        Constraint SOC_c("SOC"+to_string(cont));
        SOC_c = power(R_Wij_c, 2) + power(Im_Wij_c, 2) - Wii_c.from()*Wii_c.to();
        SOCP->add(SOC_c.in(bus_pairs) <= 0);
        
        
        /* Generator Constraints */
        Constraint Gen_resp("g_resp"+to_string(cont));
        Gen_resp = Pg_c - (Pg + p_factor*delta);
        SOCP->add(Gen_resp.in(indices(gens)) == 0);
        
        /** KCL Flow conservation */
        Constraint KCL_P_c("KCL_P_c_"+to_string(cont));
        Constraint KCL_Q_c("KCL_Q_c"+to_string(cont));
        KCL_P_c  = sum(Pf_from_c.out_arcs()) + sum(Pf_to_c.in_arcs()) + pl - sum(Pg_c.in_gens());
        KCL_Q_c  = sum(Qf_from_c.out_arcs()) + sum(Qf_to_c.in_arcs()) + ql - sum(Qg_c.in_gens());
        /* Shunts */
        KCL_P_c +=  gs*Wii_c;
        KCL_Q_c -=  bs*Wii_c;
        SOCP->add(KCL_P_c.in(nodes) == 0);
        SOCP->add(KCL_Q_c.in(nodes) == 0);
        
        Constraint Flow_P_From_c("Flow_P_From_c"+to_string(cont));
        Flow_P_From_c = Pf_from_c - (g_ff*Wii_c.from() + g_ft*R_Wij_c.in_pairs() + b_ft*Im_Wij_c.in_pairs());
        SOCP->add(Flow_P_From_c.in(arcs)==0);
        
        
        
        Constraint Flow_P_To_c("Flow_P_To_c"+to_string(cont));
        Flow_P_To_c = Pf_to_c - (g_tt*Wii_c.to() + g_tf*R_Wij_c.in_pairs() - b_tf*Im_Wij_c.in_pairs());
        SOCP->add(Flow_P_To_c.in(arcs)==0);
        
        Constraint Flow_Q_From_c("Flow_Q_From_c"+to_string(cont));
        Flow_Q_From_c = Qf_from_c - (g_ft*Im_Wij_c.in_pairs() - b_ff*Wii_c.from() - b_ft*R_Wij_c.in_pairs());
        SOCP->add(Flow_Q_From_c.in(arcs)==0);
        
        Constraint Flow_Q_To_c("Flow_Q_To_c"+to_string(cont));
        Flow_Q_To_c = Qf_to_c + (b_tt*Wii_c.to() + b_tf*R_Wij_c.in_pairs() + g_tf*Im_Wij_c.in_pairs());
        SOCP->add(Flow_Q_To_c.in(arcs)==0);
        
        /*  Thermal Limit Constraints */
        Constraint Thermal_Limit_from_c("Thermal_Limit_from_c"+to_string(cont));
        Thermal_Limit_from_c += power(Pf_from_c, 2) + power(Qf_from_c, 2);
        Thermal_Limit_from_c -= power(S_max, 2);
        SOCP->add(Thermal_Limit_from_c.in(arcs) <= 0);
        
        Constraint Thermal_Limit_to_c("Thermal_Limit_to_c"+to_string(cont));
        Thermal_Limit_to_c += power(Pf_to_c, 2) + power(Qf_to_c, 2);
        Thermal_Limit_to_c -= power(S_max,2);
        SOCP->add(Thermal_Limit_to_c.in(arcs) <= 0);
        cont++;
    }
    _conting_lines.back()->_active = true;
    if (nb_gens_cont>0) {
        clog << "Building contingency models for generators." << endl;
    }
    else {
        clog << "No generator contingencies." << endl;
    }
    for (auto i = 0; i<nb_gens_cont; i++) {
        clog << "Deactivating Generator: " << _conting_gens[i]->_name << endl;
        _conting_gens[i]->_active = false;
        if (i>0) {
            _conting_gens[i-1]->_active = true;
        }
        cont++;
    }
    
    nb_total = nb_gens_cont + nb_lines_cont;
    double end_build_time = get_wall_time();
    auto model_build_time = end_build_time - start_built_time;
    //    ACOPF->print_expanded();
    DebugOn("Total time for building continency constraints = " << model_build_time << endl);
    DebugOn("nb_vars = " << SOCP->get_nb_vars()<<endl);
    DebugOn("nb_eqs = " << SOCP->get_nb_cons() - SOCP->get_nb_ineq()<<endl);
    DebugOn("nb_ineqs = " << SOCP->get_nb_ineq()<<endl);
    
    return SOCP;
}

unique_ptr<Model> PowerNet::build_ACOPF_MINLP(PowerModelType pmt, int output, double tol, const vector<bool>& cont_in){
    DebugOn("Using MINLP model\n");
    unique_ptr<Model> ACOPF(new Model("ACOPF_MINLP Model"));
    /** Variables */
    /* Power generation variables */
    var<Real> Pg("Pg", pg_min, pg_max);
    var<Real> Qg ("Qg", qg_min, qg_max);
    ACOPF->add(Pg.in(gens), Qg.in(gens));
    Pg_base = Pg;
    Qg_base = Qg;
    
    
    /* Power flow variables */
    var<Real> Pf_from("Pf_from", S_max);
    var<Real> Qf_from("Qf_from", S_max);
    var<Real> Pf_to("Pf_to", S_max);
    var<Real> Qf_to("Qf_to", S_max);
    ACOPF->add(Pf_from.in(arcs), Qf_from.in(arcs),Pf_to.in(arcs),Qf_to.in(arcs));
    p_from_base = Pf_from;
    p_to_base = Pf_to;
    q_from_base = Qf_from;
    q_to_base = Qf_to;
    
    indices Ng; /*<< Nodes with generators */
    for (auto node: nodes) {
        auto bus = (Bus*)node;
        if (bus->_gen.size()>0) {
            Ng.add(bus->_name);
        }
    }
    
    /** Voltage related variables */
    
    
    var<Real> theta("theta_base", -6.3, 6.3);
    var<Real> v("v_base", v_min, v_max);
    
    /** Define constraints */
    
    
    
    /** Construct the objective function */
    func_ obj = product(c1, Pg) + product(c2, power(Pg,2)) + sum(c0);
    ACOPF->min(obj.in(gens));
    
    
    
    
    
    /** KCL Flow conservation */
    Constraint KCL_P("KCL_P");
    Constraint KCL_Q("KCL_Q");
    KCL_P  = sum(Pf_from.out_arcs()) + sum(Pf_to.in_arcs()) + pl - sum(Pg.in_gens());
    KCL_Q  = sum(Qf_from.out_arcs()) + sum(Qf_to.in_arcs()) + ql - sum(Qg.in_gens());
    /* Shunts */
    KCL_P +=  gs*power(v,2);
    KCL_Q -=  bs*power(v,2);
    
    ACOPF->add(KCL_P.in(nodes) == 0);
    ACOPF->add(KCL_Q.in(nodes) == 0);
    
    
    /* REF BUS */
    Constraint Ref_Bus("Ref_Bus");
    Ref_Bus = theta(get_ref_bus());

    ACOPF->add(Ref_Bus == 0);
    
    /** AC Power Flows */
    /** TODO write the constraints in Complex form */
    Constraint Flow_P_From("Flow_P_From");
    Flow_P_From += Pf_from;
    Flow_P_From -= g/power(tr,2)*power(v.from(),2);
    Flow_P_From += g/tr*(v.from()*v.to()*cos(theta.from() - theta.to() - as));
    Flow_P_From += b/tr*(v.from()*v.to()*sin(theta.from() - theta.to() - as));
    ACOPF->add(Flow_P_From.in(arcs)==0);
    
    Constraint Flow_P_To("Flow_P_To");
    Flow_P_To += Pf_to;
    Flow_P_To -= g*power(v.to(), 2);
    Flow_P_To += g/tr*(v.from()*v.to()*cos(theta.to() - theta.from() + as));
    Flow_P_To += b/tr*(v.from()*v.to()*sin(theta.to() - theta.from() + as));
    ACOPF->add(Flow_P_To.in(arcs)==0);
    
    Constraint Flow_Q_From("Flow_Q_From");
    Flow_Q_From += Qf_from;
    Flow_Q_From += (0.5*ch+b)/power(tr,2)*power(v.from(),2);
    Flow_Q_From -= b/tr*(v.from()*v.to()*cos(theta.from() - theta.to() - as));
    Flow_Q_From += g/tr*(v.from()*v.to()*sin(theta.from() - theta.to() - as));
    ACOPF->add(Flow_Q_From.in(arcs) == 0);
    
    Constraint Flow_Q_To("Flow_Q_To");
    Flow_Q_To += Qf_to;
    Flow_Q_To += (0.5*ch+b)*power(v.to(),2);
    Flow_Q_To -= b/tr*(v.from()*v.to()*cos(theta.to() - theta.from() + as));
    Flow_Q_To += g/tr*(v.from()*v.to()*sin(theta.to() - theta.from() + as));
    ACOPF->add(Flow_Q_To.in(arcs) == 0);

    
    
    /*  Thermal Limit Constraints */
    Constraint Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from += power(Pf_from, 2) + power(Qf_from, 2);
    Thermal_Limit_from -= power(S_max, 2);
    ACOPF->add(Thermal_Limit_from.in(arcs) <= 0);
    
    Constraint Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to += power(Pf_to, 2) + power(Qf_to, 2);
    Thermal_Limit_to -= power(S_max,2);
    ACOPF->add(Thermal_Limit_to.in(arcs) <= 0);
    
    unsigned nb_added = 0, nb_total = 0, nb_gens_cont = _conting_gens.size(), nb_lines_cont = _conting_lines.size();
    clog << "Building contingency models for lines." << endl;
    auto cont = 0;
    double start_built_time = get_wall_time();
    for (auto i = 0; i<nb_lines_cont; i++) {
        clog << "Switching Off Line: " << _conting_lines[i]->_name << endl;
        _conting_lines[i]->_active = false;
        if (i>0) {
            _conting_lines[i-1]->_active = true;
        }
        var<Real> Pg_c("Pg_cont"+to_string(cont), pg_min, pg_max);
        var<Real> Qg_c ("Qg_cont"+to_string(cont), qg_min, qg_max);
        //    var<Real> Sp ("Sp", pos_);
        //    var<Real> Sn ("Sn", pos_);
        ACOPF->add(Pg_c.in(gens), Qg_c.in(gens));
        Qg_cont[cont] = Qg_c;
        
        
        var<Real> Pf_from_c("Pf_from_c"+to_string(cont), S_max);
        var<Real> Qf_from_c("Qf_from_c"+to_string(cont), S_max);
        var<Real> Pf_to_c("Pf_to_c"+to_string(cont), S_max);
        var<Real> Qf_to_c("Qf_to_c"+to_string(cont), S_max);
        ACOPF->add(Pf_from_c.in(arcs), Qf_from_c.in(arcs),Pf_to_c.in(arcs),Qf_to_c.in(arcs));
        
        p_from[cont] = Pf_from_c;
        q_from[cont] = Qf_from_c;
        p_to[cont] = Pf_to_c;
        q_to[cont] = Qf_to_c;
        
        var<Real> delta("delta"+to_string(cont), pos_);
        ACOPF->add(delta.in(R(1)));
        _delta[cont] = delta;
        
        var<Real> theta_c("theta_c"+to_string(cont), -6.3, 6.3);
        var<Real> v_c("v_c"+to_string(cont), v_min, v_max);
        
        ACOPF->add(v_c.in(nodes), theta_c.in(nodes));
        v_c.initialize_all(1);
        v_cont[cont] = v_c;
        theta_cont[cont] = theta_c;

        if (cont_in[cont]) {
            /* Binary vars for complementarity, low = 1 <=> qg_c = qg_min forall gens attached to n */
            var<bool> low("low"+to_string(cont));
            ACOPF->add_var(low.in(Ng));
            /* Binary vars for complementarity, up = 1 <=> qg_c = qg_max forall gens attached to n */
            var<bool> up("up"+to_string(cont));
            ACOPF->add_var(up.in(Ng));
            
            ACOPF->_obj += sum(low) + sum(up);
            
            /* logic constraints */
    //        Constraint Disj("Disjunction"+to_string(cont));
    //        Disj = low + up + on;
    //        ACOPF->add(Disj.in(Ng) == 1);

            Constraint Disj("Disjunction"+to_string(cont));
            Disj = low + up;
            ACOPF->add(Disj.in(Ng) <= 1);

            
            Constraint On_c_u("On_c_u"+to_string(cont));
            On_c_u =  v - v_c - (low + up)*(v_max - v_min);
            ACOPF->add(On_c_u.in(Ng) <= 0);
            
            Constraint On_c_l("On_c_l"+to_string(cont));
            On_c_l =  v - v_c - (low + up)*(v_min - v_max);
            ACOPF->add(On_c_l.in(Ng) >= 0);
            
    //        Constraint Fix_on("Fix_on"+to_string(cont));
    //        Fix_on =  on;
    //        SOCP->add(Fix_on.in(Ng) == 1);
            for (auto n:nodes) {
                auto bus = (Bus*)n;
                if (bus->_gen.empty()) {
                    continue;
                }
                for (auto gen: bus->_gen) {
                    Constraint Low_c("Low_c_"+gen->_name+"_"+to_string(cont));
                    auto gname = gen->_name;
                    Low_c =  Qg_c(gname) - qg_min(gname) - (1-low(bus->_name))*(qg_max(gname) - qg_min(gname));
                    ACOPF->add(Low_c <= 0);
                    
                    Constraint Up_c("Up_c"+gen->_name+"_"+to_string(cont));
                    Up_c =  Qg_c(gname) - qg_max(gname) - (1-up(bus->_name))*(qg_min(gname) - qg_max(gname));
                    ACOPF->add(Up_c >= 0);
                }
            }
        }
        else {
            Constraint V_eq("V_eq"+to_string(cont));// v < v_c
            V_eq = (v_c - v);
            ACOPF->add(V_eq.in(Ng) == 0);
        }
        
        /* Generator Constraints */
        Constraint Gen_resp("g_resp"+to_string(cont));
        Gen_resp = Pg_c - (Pg + p_factor*delta);
        ACOPF->add(Gen_resp.in(indices(gens)) == 0);
        
        /** KCL Flow conservation */
        Constraint KCL_P_c("KCL_P_c_"+to_string(cont));
        Constraint KCL_Q_c("KCL_Q_c"+to_string(cont));
        KCL_P_c  = sum(Pf_from_c.out_arcs()) + sum(Pf_to_c.in_arcs()) + pl - sum(Pg_c.in_gens());
        KCL_Q_c  = sum(Qf_from_c.out_arcs()) + sum(Qf_to_c.in_arcs()) + ql - sum(Qg_c.in_gens());
        /* Shunts */
        KCL_P_c +=  gs*power(v_c,2);
        KCL_Q_c -=  bs*power(v_c,2);
        ACOPF->add(KCL_P_c.in(nodes) == 0);
        ACOPF->add(KCL_Q_c.in(nodes) == 0);
        
        Constraint Flow_P_From_c("Flow_P_From_c"+to_string(cont));
        Flow_P_From_c += Pf_from_c;
        Flow_P_From_c -= g/power(tr,2)*power(v_c.from(),2);
        Flow_P_From_c += g/tr*(v_c.from()*v_c.to()*cos(theta_c.from() - theta_c.to() - as));
        Flow_P_From_c += b/tr*(v_c.from()*v_c.to()*sin(theta_c.from() - theta_c.to() - as));
        ACOPF->add(Flow_P_From_c.in(arcs)==0);
        Constraint Flow_P_To_c("Flow_P_To_c"+to_string(cont));
        Flow_P_To_c += Pf_to_c;
        Flow_P_To_c -= g*power(v_c.to(), 2);
        Flow_P_To_c += g/tr*(v_c.from()*v_c.to()*cos(theta_c.to() - theta_c.from() + as));
        Flow_P_To_c += b/tr*(v_c.from()*v_c.to()*sin(theta_c.to() - theta_c.from() + as));
        ACOPF->add(Flow_P_To_c.in(arcs)==0);
        
        Constraint Flow_Q_From_c("Flow_Q_From_c"+to_string(cont));
        Flow_Q_From_c += Qf_from_c;
        Flow_Q_From_c += (0.5*ch+b)/power(tr,2)*power(v_c.from(),2);
        Flow_Q_From_c -= b/tr*(v_c.from()*v_c.to()*cos(theta_c.from() - theta_c.to() - as));
        Flow_Q_From_c += g/tr*(v_c.from()*v_c.to()*sin(theta_c.from() - theta_c.to() - as));
        ACOPF->add(Flow_Q_From_c.in(arcs)==0);
        Constraint Flow_Q_To_c("Flow_Q_To_c"+to_string(cont));
        Flow_Q_To_c += Qf_to_c;
        Flow_Q_To_c += (0.5*ch+b)*power(v_c.to(),2);
        Flow_Q_To_c -= b/tr*(v_c.from()*v_c.to()*cos(theta_c.to() - theta_c.from() + as));
        Flow_Q_To_c += g/tr*(v_c.from()*v_c.to()*sin(theta_c.to() - theta_c.from() + as));
        ACOPF->add(Flow_Q_To_c.in(arcs)==0);
        
        /*  Thermal Limit Constraints */
        Constraint Thermal_Limit_from_c("Thermal_Limit_from_c"+to_string(cont));
        Thermal_Limit_from_c += power(Pf_from_c, 2) + power(Qf_from_c, 2);
        Thermal_Limit_from_c -= power(S_max, 2);
        ACOPF->add(Thermal_Limit_from_c.in(arcs) <= 0);
        
        Constraint Thermal_Limit_to_c("Thermal_Limit_to_c"+to_string(cont));
        Thermal_Limit_to_c += power(Pf_to_c, 2) + power(Qf_to_c, 2);
        Thermal_Limit_to_c -= power(S_max,2);
        ACOPF->add(Thermal_Limit_to_c.in(arcs) <= 0);
        
        cont++;
    }
    _conting_lines.back()->_active = true;
    if (nb_gens_cont>0) {
        clog << "Building contingency models for generators." << endl;
    }
    else {
        clog << "No generator contingencies." << endl;
    }
    for (auto i = 0; i<nb_gens_cont; i++) {
        clog << "Deactivating Generator: " << _conting_gens[i]->_name << endl;
        _conting_gens[i]->_active = false;
        if (i>0) {
            _conting_gens[i-1]->_active = true;
        }
        cont++;
    }
    
    nb_total = nb_gens_cont + nb_lines_cont;
    double end_build_time = get_wall_time();
    auto model_build_time = end_build_time - start_built_time;
    //    ACOPF->print_expanded();
    DebugOn("Total time for building continency constraints = " << model_build_time << endl);
    DebugOn("nb_vars = " << ACOPF->get_nb_vars()<<endl);
    DebugOn("nb_eqs = " << ACOPF->get_nb_cons() - ACOPF->get_nb_ineq()<<endl);
    DebugOn("nb_ineqs = " << ACOPF->get_nb_ineq()<<endl);
    
    return ACOPF;
}

unique_ptr<Model> PowerNet::build_ACOPF_N_1(PowerModelType pmt, int output, double tol){
    bool polar = (pmt==ACPOL);
    if (polar) {
        Debug("Using polar model\n");
    }
    else {
        Debug("Using rectangular model\n");
    }
    unique_ptr<Model> ACOPF(new Model("AC-OPF Model"));
    /** Variables */
    /* Power generation variables */
    var<Real> Pg("Pg", pg_min, pg_max);
    var<Real> Qg ("Qg", qg_min, qg_max);
    ACOPF->add(Pg.in(gens), Qg.in(gens));
    Pg_base = Pg;
    Qg_base = Qg;
    
    
    /* Power flow variables */
    var<Real> Pf_from("Pf_from", S_max);
    var<Real> Qf_from("Qf_from", S_max);
    var<Real> Pf_to("Pf_to", S_max);
    var<Real> Qf_to("Qf_to", S_max);
    ACOPF->add(Pf_from.in(arcs), Qf_from.in(arcs),Pf_to.in(arcs),Qf_to.in(arcs));
    p_from_base = Pf_from;
    p_to_base = Pf_to;
    q_from_base = Qf_from;
    q_to_base = Qf_to;
    
    indices Ng; /*<< Nodes with generators */
    for (auto node: nodes) {
        auto bus = (Bus*)node;
        if (bus->_gen.size()>0) {
            Ng.add(bus->_name);
        }
    }
    
    /** Voltage related variables */
    
    var<Real> theta("theta_base", -6.3, 6.3);
    var<Real> v("v_base", v_min, v_max);
    
    //    var<Real> vr("vr");
    //    var<Real> vi("vi");
    var<Real> vr("vr", v_max);
    var<Real> vi("vi", v_max);
    
    if (polar) {
        ACOPF->add(v.in(nodes), theta.in(nodes));
        v.initialize_all(1);
    }
    else {
        ACOPF->add_var(vr.in(nodes));
        ACOPF->add_var(vi.in(nodes));
        vr.initialize_all(1.0);
    }
    v_base = v;
    theta_base = theta;
    
    
    /** Construct the objective function */
    func_ obj = product(c1, Pg) + product(c2, power(Pg,2)) + sum(c0);
    ACOPF->min(obj.in(gens));
    
    /** Define constraints */
    
    
    /* REF BUS */
    Constraint Ref_Bus("Ref_Bus");
    if (polar) {
        Ref_Bus = theta(get_ref_bus());
    }
    else {
        Ref_Bus = vi(get_ref_bus());
    }
    ACOPF->add(Ref_Bus == 0);
    
    /** KCL Flow conservation */
    Constraint KCL_P("KCL_P");
    Constraint KCL_Q("KCL_Q");
    KCL_P  = sum(Pf_from.out_arcs()) + sum(Pf_to.in_arcs()) + pl - sum(Pg.in_gens());
    KCL_Q  = sum(Qf_from.out_arcs()) + sum(Qf_to.in_arcs()) + ql - sum(Qg.in_gens());
    /* Shunts */
    if (polar) {
        KCL_P +=  gs*power(v,2);
        KCL_Q -=  bs*power(v,2);
    }
    else {
        KCL_P +=  gs*(power(vr,2)+power(vi,2));
        KCL_Q -=  bs*(power(vr,2)+power(vi,2));
    }
    ACOPF->add(KCL_P.in(nodes) == 0);
    ACOPF->add(KCL_Q.in(nodes) == 0);
    
    
    /** AC Power Flows */
    /** TODO write the constraints in Complex form */
    Constraint Flow_P_From("Flow_P_From");
    Flow_P_From += Pf_from;
    if (polar) {
        Flow_P_From -= g/power(tr,2)*power(v.from(),2);
        Flow_P_From += g/tr*(v.from()*v.to()*cos(theta.from() - theta.to() - as));
        Flow_P_From += b/tr*(v.from()*v.to()*sin(theta.from() - theta.to() - as));
    }
    else {
        Flow_P_From -= g_ff*(power(vr.from(), 2) + power(vi.from(), 2));
        Flow_P_From -= g_ft*(vr.from()*vr.to() + vi.from()*vi.to());
        Flow_P_From -= b_ft*(vi.from()*vr.to() - vr.from()*vi.to());
    }
    ACOPF->add(Flow_P_From.in(arcs)==0);
    
    Constraint Flow_P_To("Flow_P_To");
    Flow_P_To += Pf_to;
    if (polar) {
        Flow_P_To -= g*power(v.to(), 2);
        Flow_P_To += g/tr*(v.from()*v.to()*cos(theta.to() - theta.from() + as));
        Flow_P_To += b/tr*(v.from()*v.to()*sin(theta.to() - theta.from() + as));
    }
    else {
        Flow_P_To -= g_tt*(power(vr.to(), 2) + power(vi.to(), 2));
        Flow_P_To -= g_tf*(vr.from()*vr.to() + vi.from()*vi.to());
        Flow_P_To -= b_tf*(vi.to()*vr.from() - vr.to()*vi.from());
    }
    ACOPF->add(Flow_P_To.in(arcs)==0);
    
    Constraint Flow_Q_From("Flow_Q_From");
    Flow_Q_From += Qf_from;
    if (polar) {
        Flow_Q_From += (0.5*ch+b)/power(tr,2)*power(v.from(),2);
        Flow_Q_From -= b/tr*(v.from()*v.to()*cos(theta.from() - theta.to() - as));
        Flow_Q_From += g/tr*(v.from()*v.to()*sin(theta.from() - theta.to() - as));
    }
    else {
        Flow_Q_From += b_ff*(power(vr.from(), 2) + power(vi.from(), 2));
        Flow_Q_From += b_ft*(vr.from()*vr.to() + vi.from()*vi.to());
        Flow_Q_From -= g_ft*(vi.from()*vr.to() - vr.from()*vi.to());
    }
    ACOPF->add(Flow_Q_From.in(arcs)==0);
    Constraint Flow_Q_To("Flow_Q_To");
    Flow_Q_To += Qf_to;
    if (polar) {
        Flow_Q_To += (0.5*ch+b)*power(v.to(),2);
        Flow_Q_To -= b/tr*(v.from()*v.to()*cos(theta.to() - theta.from() + as));
        Flow_Q_To += g/tr*(v.from()*v.to()*sin(theta.to() - theta.from() + as));
    }
    else {
        Flow_Q_To += b_tt*(power(vr.to(), 2) + power(vi.to(), 2));
        Flow_Q_To += b_tf*(vr.from()*vr.to() + vi.from()*vi.to());
        Flow_Q_To -= g_tf*(vi.to()*vr.from() - vr.to()*vi.from());
    }
    ACOPF->add(Flow_Q_To.in(arcs)==0);
    
    /** AC voltage limit constraints. */
    if (!polar) {
        Constraint Vol_limit_UB("Vol_limit_UB");
        Vol_limit_UB = power(vr, 2) + power(vi, 2);
        Vol_limit_UB -= power(v_max, 2);
        ACOPF->add(Vol_limit_UB.in(nodes) <= 0);
        
        Constraint Vol_limit_LB("Vol_limit_LB");
        Vol_limit_LB = power(vr, 2) + power(vi, 2);
        Vol_limit_LB -= power(v_min,2);
        ACOPF->add(Vol_limit_LB.in(nodes) >= 0);
        DebugOff(v_min.to_str(true) << endl);
        DebugOff(v_max.to_str(true) << endl);
    }
    
    
    /*  Thermal Limit Constraints */
    Constraint Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from += power(Pf_from, 2) + power(Qf_from, 2);
    Thermal_Limit_from -= power(S_max, 2);
    ACOPF->add(Thermal_Limit_from.in(arcs) <= 0);
    
    Constraint Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to += power(Pf_to, 2) + power(Qf_to, 2);
    Thermal_Limit_to -= power(S_max,2);
    ACOPF->add(Thermal_Limit_to.in(arcs) <= 0);
    //    ACOPF->print_expanded();
    //    ACOPF->_obj*=1e-4;
    
    unsigned nb_added = 0, nb_total = 0, nb_gens_cont = _conting_gens.size(), nb_lines_cont = _conting_lines.size();
    clog << "Building contingency models for lines." << endl;
    auto cont = 0;
    double start_built_time = get_wall_time();
    for (auto i = 0; i<nb_lines_cont; i++) {
        clog << "Switching Off Line: " << _conting_lines[i]->_name << endl;
        _conting_lines[i]->_active = false;
        if (i>0) {
            _conting_lines[i-1]->_active = true;
        }
        var<Real> Pg_c("Pg_cont"+to_string(cont), pg_min, pg_max);
        var<Real> Qg_c ("Qg_cont"+to_string(cont), qg_min, qg_max);
        //    var<Real> Sp ("Sp", pos_);
        //    var<Real> Sn ("Sn", pos_);
        ACOPF->add(Pg_c.in(gens), Qg_c.in(gens));
        Qg_cont[cont] = Qg_c;
        
        
        var<Real> Pf_from_c("Pf_from_c"+to_string(cont), S_max);
        var<Real> Qf_from_c("Qf_from_c"+to_string(cont), S_max);
        var<Real> Pf_to_c("Pf_to_c"+to_string(cont), S_max);
        var<Real> Qf_to_c("Qf_to_c"+to_string(cont), S_max);
        ACOPF->add(Pf_from_c.in(arcs), Qf_from_c.in(arcs),Pf_to_c.in(arcs),Qf_to_c.in(arcs));
        
        p_from[cont] = Pf_from_c;
        q_from[cont] = Qf_from_c;
        p_to[cont] = Pf_to_c;
        q_to[cont] = Qf_to_c;
        
        var<Real> delta("delta"+to_string(cont), pos_);
        ACOPF->add(delta.in(R(1)));
        _delta[cont] = delta;
        
        var<Real> theta_c("theta_c"+to_string(cont), -6.3, 6.3);
        var<Real> v_c("v_c"+to_string(cont), v_min, v_max);
        
        ACOPF->add(v_c.in(nodes), theta_c.in(nodes));
        v_c.initialize_all(1);
        v_cont[cont] = v_c;
        theta_cont[cont] = theta_c;
        
        Constraint V_eq("V_eq"+to_string(cont));// v < v_c
        V_eq = (v_c - v);
        ACOPF->add(V_eq.in(Ng) == 0);
        
        
        
        /* Generator Constraints */
        Constraint Gen_resp("g_resp"+to_string(cont));
        Gen_resp = Pg_c - (Pg + p_factor*delta);
        ACOPF->add(Gen_resp.in(indices(gens)) == 0);
        /** KCL Flow conservation */
        Constraint KCL_P_c("KCL_P_c_"+to_string(cont));
        Constraint KCL_Q_c("KCL_Q_c"+to_string(cont));
        KCL_P_c  = sum(Pf_from_c.out_arcs()) + sum(Pf_to_c.in_arcs()) + pl - sum(Pg_c.in_gens());
        KCL_Q_c  = sum(Qf_from_c.out_arcs()) + sum(Qf_to_c.in_arcs()) + ql - sum(Qg_c.in_gens());
        /* Shunts */
        if (polar) {
            KCL_P_c +=  gs*power(v_c,2);
            KCL_Q_c -=  bs*power(v_c,2);
        }
        ACOPF->add(KCL_P_c.in(nodes) == 0);
        ACOPF->add(KCL_Q_c.in(nodes) == 0);
        
        Constraint Flow_P_From_c("Flow_P_From_c"+to_string(cont));
        Flow_P_From_c += Pf_from_c;
        Flow_P_From_c -= g/power(tr,2)*power(v_c.from(),2);
        Flow_P_From_c += g/tr*(v_c.from()*v_c.to()*cos(theta_c.from() - theta_c.to() - as));
        Flow_P_From_c += b/tr*(v_c.from()*v_c.to()*sin(theta_c.from() - theta_c.to() - as));
        ACOPF->add(Flow_P_From_c.in(arcs)==0);
        
        Constraint Flow_P_To_c("Flow_P_To_c"+to_string(cont));
        Flow_P_To_c += Pf_to_c;
        Flow_P_To_c -= g*power(v_c.to(), 2);
        Flow_P_To_c += g/tr*(v_c.from()*v_c.to()*cos(theta_c.to() - theta_c.from() + as));
        Flow_P_To_c += b/tr*(v_c.from()*v_c.to()*sin(theta_c.to() - theta_c.from() + as));
        ACOPF->add(Flow_P_To_c.in(arcs)==0);
        
        Constraint Flow_Q_From_c("Flow_Q_From_c"+to_string(cont));
        Flow_Q_From_c += Qf_from_c;
        Flow_Q_From_c += (0.5*ch+b)/power(tr,2)*power(v_c.from(),2);
        Flow_Q_From_c -= b/tr*(v_c.from()*v_c.to()*cos(theta_c.from() - theta_c.to() - as));
        Flow_Q_From_c += g/tr*(v_c.from()*v_c.to()*sin(theta_c.from() - theta_c.to() - as));
        ACOPF->add(Flow_Q_From_c.in(arcs)==0);
        Constraint Flow_Q_To_c("Flow_Q_To_c"+to_string(cont));
        Flow_Q_To_c += Qf_to_c;
        Flow_Q_To_c += (0.5*ch+b)*power(v_c.to(),2);
        Flow_Q_To_c -= b/tr*(v_c.from()*v_c.to()*cos(theta_c.to() - theta_c.from() + as));
        Flow_Q_To_c += g/tr*(v_c.from()*v_c.to()*sin(theta_c.to() - theta_c.from() + as));
        ACOPF->add(Flow_Q_To_c.in(arcs)==0);
        
        /*  Thermal Limit Constraints */
        Constraint Thermal_Limit_from_c("Thermal_Limit_from_c"+to_string(cont));
        Thermal_Limit_from_c += power(Pf_from_c, 2) + power(Qf_from_c, 2);
        Thermal_Limit_from_c -= power(S_max, 2);
        ACOPF->add(Thermal_Limit_from_c.in(arcs) <= 0);
        
        Constraint Thermal_Limit_to_c("Thermal_Limit_to_c"+to_string(cont));
        Thermal_Limit_to_c += power(Pf_to_c, 2) + power(Qf_to_c, 2);
        Thermal_Limit_to_c -= power(S_max,2);
        ACOPF->add(Thermal_Limit_to_c.in(arcs) <= 0);
        cont++;
    }
    _conting_lines.back()->_active = true;
    if (nb_gens_cont>0) {
        clog << "Building contingency models for generators." << endl;
    }
    else {
        clog << "No generator contingencies." << endl;
    }
    for (auto i = 0; i<nb_gens_cont; i++) {
        clog << "Deactivating Generator: " << _conting_gens[i]->_name << endl;
        _conting_gens[i]->_active = false;
        if (i>0) {
            _conting_gens[i-1]->_active = true;
        }
        cont++;
    }
    nb_total = nb_gens_cont + nb_lines_cont;
    double end_build_time = get_wall_time();
    auto model_build_time = end_build_time - start_built_time;
    //    ACOPF->print_expanded();
    DebugOn("Total time for building continency constraints = " << model_build_time << endl);
    DebugOn("nb_vars = " << ACOPF->get_nb_vars()<<endl);
    DebugOn("nb_eqs = " << ACOPF->get_nb_cons() - ACOPF->get_nb_ineq()<<endl);
    DebugOn("nb_ineqs = " << ACOPF->get_nb_ineq()<<endl);
    
    return ACOPF;
}

unique_ptr<Model> PowerNet::build_soft_ACOPF_N_1(PowerModelType pmt, int output, double tol, double obj_pen){
    bool polar = (pmt==ACPOL);
    if (polar) {
        Debug("Using polar model\n");
    }
    else {
        Debug("Using rectangular model\n");
    }
    unique_ptr<Model> ACOPF(new Model("AC-OPF Model"));
    /** Variables */
    /* Power generation variables */
    var<Real> Pg("Pg", pg_min, pg_max);
    var<Real> Qg ("Qg", qg_min, qg_max);
    ACOPF->add(Pg.in(gens), Qg.in(gens));
    Pg_base = Pg;
    Qg_base = Qg;
//    Pg.initialize_uniform();
//    Qg.initialize_uniform();
    
    /* Power flow variables */
    var<Real> Pf_from("Pf_from", S_max);
    var<Real> Qf_from("Qf_from", S_max);
    var<Real> Pf_to("Pf_to", S_max);
    var<Real> Qf_to("Qf_to", S_max);
    ACOPF->add(Pf_from.in(arcs), Qf_from.in(arcs),Pf_to.in(arcs),Qf_to.in(arcs));
    p_from_base = Pf_from;
    p_to_base = Pf_to;
    q_from_base = Qf_from;
    q_to_base = Qf_to;
//    Pf_from.initialize_normal(0, 1);
//    Qf_from.initialize_normal(0, 1);
//    Pf_to.initialize_normal(0, 1);
//    Qf_to.initialize_normal(0, 1);

    indices Ng; /*<< Nodes with generators */
    for (auto node: nodes) {
        auto bus = (Bus*)node;
        if (bus->_gen.size()>0) {
            Ng.add(bus->_name);
        }
    }
    
    /** Voltage related variables */

    var<Real> theta("theta_base", -6.3, 6.3);
    var<Real> v("v_base", v_min, v_max);

    //    var<Real> vr("vr");
    //    var<Real> vi("vi");
    var<Real> vr("vr", v_max);
    var<Real> vi("vi", v_max);
    
    if (polar) {
        ACOPF->add(v.in(nodes), theta.in(nodes));
        v.initialize_all(1);
//        v.initialize_uniform();
//        theta.initialize_normal(0, 1);
    }
    else {
        ACOPF->add_var(vr.in(nodes));
        ACOPF->add_var(vi.in(nodes));
        vr.initialize_all(1.0);
    }
    v_base = v;
    theta_base = theta;
    
    
    /** Construct the objective function */
    func_ obj = product(c1, Pg) + product(c2, power(Pg,2)) + sum(c0);
    ACOPF->min(obj.in(gens));
    
    /** Define constraints */
    
    
    /* REF BUS */
    Constraint Ref_Bus("Ref_Bus");
    if (polar) {
        Ref_Bus = theta(get_ref_bus());
    }
    else {
        Ref_Bus = vi(get_ref_bus());
    }
        ACOPF->add(Ref_Bus == 0);
    
    /** KCL Flow conservation */
    Constraint KCL_P("KCL_P");
    Constraint KCL_Q("KCL_Q");
    KCL_P  = sum(Pf_from.out_arcs()) + sum(Pf_to.in_arcs()) + pl - sum(Pg.in_gens());
    KCL_Q  = sum(Qf_from.out_arcs()) + sum(Qf_to.in_arcs()) + ql - sum(Qg.in_gens());
    /* Shunts */
    if (polar) {
        KCL_P +=  gs*power(v,2);
        KCL_Q -=  bs*power(v,2);
    }
    else {
        KCL_P +=  gs*(power(vr,2)+power(vi,2));
        KCL_Q -=  bs*(power(vr,2)+power(vi,2));
    }
    ACOPF->add(KCL_P.in(nodes) == 0);
    ACOPF->add(KCL_Q.in(nodes) == 0);
    
    
    /** AC Power Flows */
    /** TODO write the constraints in Complex form */
    Constraint Flow_P_From("Flow_P_From");
    Flow_P_From += Pf_from;
    if (polar) {
        Flow_P_From -= g/power(tr,2)*power(v.from(),2);
        Flow_P_From += g/tr*(v.from()*v.to()*cos(theta.from() - theta.to() - as));
        Flow_P_From += b/tr*(v.from()*v.to()*sin(theta.from() - theta.to() - as));
    }
    else {
        Flow_P_From -= g_ff*(power(vr.from(), 2) + power(vi.from(), 2));
        Flow_P_From -= g_ft*(vr.from()*vr.to() + vi.from()*vi.to());
        Flow_P_From -= b_ft*(vi.from()*vr.to() - vr.from()*vi.to());
    }
    ACOPF->add(Flow_P_From.in(arcs)==0);
    
    Constraint Flow_P_To("Flow_P_To");
    Flow_P_To += Pf_to;
    if (polar) {
        Flow_P_To -= g*power(v.to(), 2);
        Flow_P_To += g/tr*(v.from()*v.to()*cos(theta.to() - theta.from() + as));
        Flow_P_To += b/tr*(v.from()*v.to()*sin(theta.to() - theta.from() + as));
    }
    else {
        Flow_P_To -= g_tt*(power(vr.to(), 2) + power(vi.to(), 2));
        Flow_P_To -= g_tf*(vr.from()*vr.to() + vi.from()*vi.to());
        Flow_P_To -= b_tf*(vi.to()*vr.from() - vr.to()*vi.from());
    }
    ACOPF->add(Flow_P_To.in(arcs)==0);
    
    Constraint Flow_Q_From("Flow_Q_From");
    Flow_Q_From += Qf_from;
    if (polar) {
        Flow_Q_From += (0.5*ch+b)/power(tr,2)*power(v.from(),2);
        Flow_Q_From -= b/tr*(v.from()*v.to()*cos(theta.from() - theta.to() - as));
        Flow_Q_From += g/tr*(v.from()*v.to()*sin(theta.from() - theta.to() - as));
    }
    else {
        Flow_Q_From += b_ff*(power(vr.from(), 2) + power(vi.from(), 2));
        Flow_Q_From += b_ft*(vr.from()*vr.to() + vi.from()*vi.to());
        Flow_Q_From -= g_ft*(vi.from()*vr.to() - vr.from()*vi.to());
    }
    ACOPF->add(Flow_Q_From.in(arcs)==0);
    Constraint Flow_Q_To("Flow_Q_To");
    Flow_Q_To += Qf_to;
    if (polar) {
        Flow_Q_To += (0.5*ch+b)*power(v.to(),2);
        Flow_Q_To -= b/tr*(v.from()*v.to()*cos(theta.to() - theta.from() + as));
        Flow_Q_To += g/tr*(v.from()*v.to()*sin(theta.to() - theta.from() + as));
    }
    else {
        Flow_Q_To += b_tt*(power(vr.to(), 2) + power(vi.to(), 2));
        Flow_Q_To += b_tf*(vr.from()*vr.to() + vi.from()*vi.to());
        Flow_Q_To -= g_tf*(vi.to()*vr.from() - vr.to()*vi.from());
    }
    ACOPF->add(Flow_Q_To.in(arcs)==0);
    
    /** AC voltage limit constraints. */
    if (!polar) {
        Constraint Vol_limit_UB("Vol_limit_UB");
        Vol_limit_UB = power(vr, 2) + power(vi, 2);
        Vol_limit_UB -= power(v_max, 2);
        ACOPF->add(Vol_limit_UB.in(nodes) <= 0);
        
        Constraint Vol_limit_LB("Vol_limit_LB");
        Vol_limit_LB = power(vr, 2) + power(vi, 2);
        Vol_limit_LB -= power(v_min,2);
        ACOPF->add(Vol_limit_LB.in(nodes) >= 0);
        DebugOff(v_min.to_str(true) << endl);
        DebugOff(v_max.to_str(true) << endl);
    }
    
    
    /*  Thermal Limit Constraints */
    Constraint Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from += power(Pf_from, 2) + power(Qf_from, 2);
    Thermal_Limit_from -= power(S_max, 2);
    ACOPF->add(Thermal_Limit_from.in(arcs) <= 0);
    
    Constraint Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to += power(Pf_to, 2) + power(Qf_to, 2);
    Thermal_Limit_to -= power(S_max,2);
    ACOPF->add(Thermal_Limit_to.in(arcs) <= 0);
    //    ACOPF->print_expanded();
//    ACOPF->_obj*=1e-4;
    
    unsigned nb_added = 0, nb_total = 0, nb_gens_cont = _conting_gens.size(), nb_lines_cont = _conting_lines.size();
    clog << "Building contingency models for lines." << endl;
    auto cont = 0;
    double start_built_time = get_wall_time();
    for (auto i = 0; i<nb_lines_cont; i++) {
        clog << "Switching Off Line: " << _conting_lines[i]->_name << endl;
        _conting_lines[i]->_active = false;
        if (i>0) {
            _conting_lines[i-1]->_active = true;
        }
        var<Real> Pg_c("Pg_cont"+to_string(cont), pg_min, pg_max);
        var<Real> Qg_c ("Qg_cont"+to_string(cont), qg_min, qg_max);
        //    var<Real> Sp ("Sp", pos_);
        //    var<Real> Sn ("Sn", pos_);
        ACOPF->add(Pg_c.in(gens), Qg_c.in(gens));
        Qg_cont[cont] = Qg_c;

//        Pg_c.initialize_uniform();
//        Qg_c.initialize_uniform();
        
        var<Real> v_diff_p("v_diff_p"+to_string(cont), pos_);
        ACOPF->add(v_diff_p.in(Ng));
        var<Real> v_diff_n("v_diff_n"+to_string(cont), pos_);
        ACOPF->add(v_diff_n.in(Ng));
        
        v_diff_p_cont[cont] = v_diff_p;
        v_diff_n_cont[cont] = v_diff_n;
//
        ACOPF->_obj += obj_pen*sum(v_diff_p) + obj_pen*sum(v_diff_n);
        
//        ACOPF->_obj*=1e-4;
        
        var<Real> Pf_from_c("Pf_from_c"+to_string(cont), S_max);
        var<Real> Qf_from_c("Qf_from_c"+to_string(cont), S_max);
        var<Real> Pf_to_c("Pf_to_c"+to_string(cont), S_max);
        var<Real> Qf_to_c("Qf_to_c"+to_string(cont), S_max);
        ACOPF->add(Pf_from_c.in(arcs), Qf_from_c.in(arcs),Pf_to_c.in(arcs),Qf_to_c.in(arcs));
        
        p_from[cont] = Pf_from_c;
        q_from[cont] = Qf_from_c;
        p_to[cont] = Pf_to_c;
        q_to[cont] = Qf_to_c;
        
//        Pf_from_c.initialize_normal(0, 1);
//        Qf_from_c.initialize_normal(0, 1);
//        Pf_to_c.initialize_normal(0, 1);
//        Qf_to_c.initialize_normal(0, 1);
        
        var<Real> delta("delta"+to_string(cont), pos_);
        ACOPF->add(delta.in(R(1)));
        _delta[cont] = delta;
        
        var<Real> theta_c("theta_c"+to_string(cont), -6.3, 6.3);
        var<Real> v_c("v_c"+to_string(cont), v_min, v_max);
        
        ACOPF->add(v_c.in(nodes), theta_c.in(nodes));
        v_c.initialize_all(1);
        v_cont[cont] = v_c;
        theta_cont[cont] = theta_c;
        
//        v_c.initialize_uniform();
//        theta_c.initialize_normal(0, 1);
        
//        Constraint V_eq("V_eq"+to_string(cont));// v < v_c
//        V_eq = (v_c - v);
//        ACOPF->add(V_eq.in(Ng) == 0);

        Constraint V_diff_def("V_diff_def"+to_string(cont));// v < v_c
        V_diff_def = 1e3*((v_diff_p - v_diff_n) - (v_c - v));
//        V_diff_def = ((v_diff_p - v_diff_n) - (v_c - v));
        ACOPF->add(V_diff_def.in(Ng) == 0);
        
//        Constraint V_diff_p("V_diff_p"+to_string(cont));// v < v_c
//        V_diff_p = v_diff_p - (v_c - v);
//        ACOPF->add(V_diff_p.in(Ng) >= 0);
//        
//        Constraint V_diff_n("V_diff_n"+to_string(cont));// v > v_c
//        V_diff_n = v_diff_n - (v - v_c);
//        ACOPF->add(V_diff_n.in(Ng) >= 0);
//
//        Constraint Complement_p("Complement_p"+to_string(cont));
//        Complement_p = v_diff_p*(Qg_c - qg_min);
//        ACOPF->add(Complement_p.in(Ng) <= 0);
////
//        Constraint Complement_n("Complement_n"+to_string(cont));
//        Complement_n = v_diff_n*(qg_max - Qg_c);
//        ACOPF->add(Complement_n.in(Ng) <= 0);
        
        
        
        /* Generator Constraints */
        Constraint Gen_resp("g_resp"+to_string(cont));
        Gen_resp = Pg_c - (Pg + p_factor*delta);
        ACOPF->add(Gen_resp.in(indices(gens)) == 0);
        /** KCL Flow conservation */
        Constraint KCL_P_c("KCL_P_c_"+to_string(cont));
        Constraint KCL_Q_c("KCL_Q_c"+to_string(cont));
        KCL_P_c  = sum(Pf_from_c.out_arcs()) + sum(Pf_to_c.in_arcs()) + pl - sum(Pg_c.in_gens());
        KCL_Q_c  = sum(Qf_from_c.out_arcs()) + sum(Qf_to_c.in_arcs()) + ql - sum(Qg_c.in_gens());
        /* Shunts */
        if (polar) {
            KCL_P_c +=  gs*power(v_c,2);
            KCL_Q_c -=  bs*power(v_c,2);
        }
        ACOPF->add(KCL_P_c.in(nodes) == 0);
        ACOPF->add(KCL_Q_c.in(nodes) == 0);
        
        Constraint Flow_P_From_c("Flow_P_From_c"+to_string(cont));
        Flow_P_From_c += Pf_from_c;
        Flow_P_From_c -= g/power(tr,2)*power(v_c.from(),2);
        Flow_P_From_c += g/tr*(v_c.from()*v_c.to()*cos(theta_c.from() - theta_c.to() - as));
        Flow_P_From_c += b/tr*(v_c.from()*v_c.to()*sin(theta_c.from() - theta_c.to() - as));
        ACOPF->add(Flow_P_From_c.in(arcs)==0);
        
        Constraint Flow_P_To_c("Flow_P_To_c"+to_string(cont));
        Flow_P_To_c += Pf_to_c;
        Flow_P_To_c -= g*power(v_c.to(), 2);
        Flow_P_To_c += g/tr*(v_c.from()*v_c.to()*cos(theta_c.to() - theta_c.from() + as));
        Flow_P_To_c += b/tr*(v_c.from()*v_c.to()*sin(theta_c.to() - theta_c.from() + as));
        ACOPF->add(Flow_P_To_c.in(arcs)==0);
        
        Constraint Flow_Q_From_c("Flow_Q_From_c"+to_string(cont));
        Flow_Q_From_c += Qf_from_c;
        Flow_Q_From_c += (0.5*ch+b)/power(tr,2)*power(v_c.from(),2);
        Flow_Q_From_c -= b/tr*(v_c.from()*v_c.to()*cos(theta_c.from() - theta_c.to() - as));
        Flow_Q_From_c += g/tr*(v_c.from()*v_c.to()*sin(theta_c.from() - theta_c.to() - as));
        ACOPF->add(Flow_Q_From_c.in(arcs)==0);
        Constraint Flow_Q_To_c("Flow_Q_To_c"+to_string(cont));
        Flow_Q_To_c += Qf_to_c;
        Flow_Q_To_c += (0.5*ch+b)*power(v_c.to(),2);
        Flow_Q_To_c -= b/tr*(v_c.from()*v_c.to()*cos(theta_c.to() - theta_c.from() + as));
        Flow_Q_To_c += g/tr*(v_c.from()*v_c.to()*sin(theta_c.to() - theta_c.from() + as));
        ACOPF->add(Flow_Q_To_c.in(arcs)==0);
        
        /*  Thermal Limit Constraints */
        Constraint Thermal_Limit_from_c("Thermal_Limit_from_c"+to_string(cont));
        Thermal_Limit_from_c += power(Pf_from_c, 2) + power(Qf_from_c, 2);
        Thermal_Limit_from_c -= power(S_max, 2);
        ACOPF->add(Thermal_Limit_from_c.in(arcs) <= 0);
        
        Constraint Thermal_Limit_to_c("Thermal_Limit_to_c"+to_string(cont));
        Thermal_Limit_to_c += power(Pf_to_c, 2) + power(Qf_to_c, 2);
        Thermal_Limit_to_c -= power(S_max,2);
        ACOPF->add(Thermal_Limit_to_c.in(arcs) <= 0);
        cont++;
    }
    _conting_lines.back()->_active = true;
    if (nb_gens_cont>0) {
        clog << "Building contingency models for generators." << endl;
    }
    else {
        clog << "No generator contingencies." << endl;
    }
    for (auto i = 0; i<nb_gens_cont; i++) {
        clog << "Deactivating Generator: " << _conting_gens[i]->_name << endl;
        _conting_gens[i]->_active = false;
        if (i>0) {
            _conting_gens[i-1]->_active = true;
        }
        cont++;
    }
    nb_total = nb_gens_cont + nb_lines_cont;
    double end_build_time = get_wall_time();
    auto model_build_time = end_build_time - start_built_time;
//    ACOPF->print_expanded();
    DebugOn("Total time for building continency constraints = " << model_build_time << endl);
    DebugOn("nb_vars = " << ACOPF->get_nb_vars()<<endl);
    DebugOn("nb_eqs = " << ACOPF->get_nb_cons() - ACOPF->get_nb_ineq()<<endl);
    DebugOn("nb_ineqs = " << ACOPF->get_nb_ineq()<<endl);
    
    return ACOPF;
}

unique_ptr<Model> PowerNet::build_ACOPF(PowerModelType pmt, int output, double tol){
    bool polar = (pmt==ACPOL);
    if (polar) {
        Debug("Using polar model\n");
    }
    else {
        Debug("Using rectangular model\n");
    }
    unique_ptr<Model> ACOPF(new Model("AC-OPF Model"));
    /** Variables */
    /* Power generation variables */
    var<Real> Pg("Pg", pg_min, pg_max);
    var<Real> Qg ("Qg", qg_min, qg_max);
    ACOPF->add(Pg.in(gens), Qg.in(gens));
    Pg_base = Pg;
    Qg_base = Qg;
    
    /* Power flow variables */
    var<Real> Pf_from("Pf_from", S_max);
    var<Real> Qf_from("Qf_from", S_max);
    var<Real> Pf_to("Pf_to", S_max);
    var<Real> Qf_to("Qf_to", S_max);
    ACOPF->add(Pf_from.in(arcs), Qf_from.in(arcs),Pf_to.in(arcs),Qf_to.in(arcs));
    p_from_base = Pf_from;
    p_to_base = Pf_to;
    q_from_base = Qf_from;
    q_to_base = Qf_to;
    /** Voltage related variables */
    var<Real> theta("theta", -6.3, 6.3);
    var<Real> v("|V|", v_min, v_max);
    //    var<Real> vr("vr");
    //    var<Real> vi("vi");
    var<Real> vr("vr", v_max);
    var<Real> vi("vi", v_max);
    
    if (polar) {
        ACOPF->add_var(v.in(nodes));
        ACOPF->add_var(theta.in(nodes));
        v.initialize_all(1);
    }
    else {
        ACOPF->add_var(vr.in(nodes));
        ACOPF->add_var(vi.in(nodes));
        vr.initialize_all(1.0);
    }
    v_base = v;
    theta_base = theta;

    /** Construct the objective function */
    func_ obj = product(c1, Pg) + product(c2, power(Pg,2)) + sum(c0);
    ACOPF->min(obj.in(gens));
    
    /** Define constraints */
    
    /* REF BUS */
    Constraint Ref_Bus("Ref_Bus");
    if (polar) {
        Ref_Bus = theta(get_ref_bus());
    }
    else {
        Ref_Bus = vi(get_ref_bus());
    }
    ACOPF->add(Ref_Bus == 0);
    
    /** KCL Flow conservation */
    Constraint KCL_P("KCL_P");
    Constraint KCL_Q("KCL_Q");
    KCL_P  = sum(Pf_from.out_arcs()) + sum(Pf_to.in_arcs()) + pl - sum(Pg.in_gens());
    KCL_Q  = sum(Qf_from.out_arcs()) + sum(Qf_to.in_arcs()) + ql - sum(Qg.in_gens());
    /* Shunts */
    if (polar) {
        KCL_P +=  gs*power(v,2);
        KCL_Q -=  bs*power(v,2);
    }
    else {
        KCL_P +=  gs*(power(vr,2)+power(vi,2));
        KCL_Q -=  bs*(power(vr,2)+power(vi,2));
    }
    ACOPF->add(KCL_P.in(nodes) == 0);
    ACOPF->add(KCL_Q.in(nodes) == 0);
    
    /** AC Power Flows */
    /** TODO write the constraints in Complex form */
    Constraint Flow_P_From("Flow_P_From");
    Flow_P_From += Pf_from;
    if (polar) {
        Flow_P_From -= g/power(tr,2)*power(v.from(),2);
        Flow_P_From += g/tr*(v.from()*v.to()*cos(theta.from() - theta.to() - as));
        Flow_P_From += b/tr*(v.from()*v.to()*sin(theta.from() - theta.to() - as));
    }
    else {
        Flow_P_From -= g_ff*(power(vr.from(), 2) + power(vi.from(), 2));
        Flow_P_From -= g_ft*(vr.from()*vr.to() + vi.from()*vi.to());
        Flow_P_From -= b_ft*(vi.from()*vr.to() - vr.from()*vi.to());
    }
    ACOPF->add(Flow_P_From.in(arcs)==0);
    
    Constraint Flow_P_To("Flow_P_To");
    Flow_P_To += Pf_to;
    if (polar) {
        Flow_P_To -= g*power(v.to(), 2);
        Flow_P_To += g/tr*(v.from()*v.to()*cos(theta.to() - theta.from() + as));
        Flow_P_To += b/tr*(v.from()*v.to()*sin(theta.to() - theta.from() + as));
    }
    else {
        Flow_P_To -= g_tt*(power(vr.to(), 2) + power(vi.to(), 2));
        Flow_P_To -= g_tf*(vr.from()*vr.to() + vi.from()*vi.to());
        Flow_P_To -= b_tf*(vi.to()*vr.from() - vr.to()*vi.from());
    }
    ACOPF->add(Flow_P_To.in(arcs)==0);
    
    Constraint Flow_Q_From("Flow_Q_From");
    Flow_Q_From += Qf_from;
    if (polar) {
        Flow_Q_From += (0.5*ch+b)/power(tr,2)*power(v.from(),2);
        Flow_Q_From -= b/tr*(v.from()*v.to()*cos(theta.from() - theta.to() - as));
        Flow_Q_From += g/tr*(v.from()*v.to()*sin(theta.from() - theta.to() - as));
    }
    else {
        Flow_Q_From += b_ff*(power(vr.from(), 2) + power(vi.from(), 2));
        Flow_Q_From += b_ft*(vr.from()*vr.to() + vi.from()*vi.to());
        Flow_Q_From -= g_ft*(vi.from()*vr.to() - vr.from()*vi.to());
    }
    ACOPF->add(Flow_Q_From.in(arcs)==0);
    Constraint Flow_Q_To("Flow_Q_To");
    Flow_Q_To += Qf_to;
    if (polar) {
        Flow_Q_To += (0.5*ch+b)*power(v.to(),2);
        Flow_Q_To -= b/tr*(v.from()*v.to()*cos(theta.to() - theta.from() + as));
        Flow_Q_To += g/tr*(v.from()*v.to()*sin(theta.to() - theta.from() + as));
    }
    else {
        Flow_Q_To += b_tt*(power(vr.to(), 2) + power(vi.to(), 2));
        Flow_Q_To += b_tf*(vr.from()*vr.to() + vi.from()*vi.to());
        Flow_Q_To -= g_tf*(vi.to()*vr.from() - vr.to()*vi.from());
    }
    ACOPF->add(Flow_Q_To.in(arcs)==0);
    
    /** AC voltage limit constraints. */
    if (!polar) {
        Constraint Vol_limit_UB("Vol_limit_UB");
        Vol_limit_UB = power(vr, 2) + power(vi, 2);
        Vol_limit_UB -= power(v_max, 2);
        ACOPF->add(Vol_limit_UB.in(nodes) <= 0);
        
        Constraint Vol_limit_LB("Vol_limit_LB");
        Vol_limit_LB = power(vr, 2) + power(vi, 2);
        Vol_limit_LB -= power(v_min,2);
        ACOPF->add(Vol_limit_LB.in(nodes) >= 0);
        DebugOff(v_min.to_str(true) << endl);
        DebugOff(v_max.to_str(true) << endl);
    }
    
    
    /* Phase Angle Bounds constraints */
    Constraint PAD_UB("PAD_UB");
    Constraint PAD_LB("PAD_LB");
    auto bus_pairs = get_bus_pairs();
    if (polar) {
        PAD_UB = theta.from() - theta.to();
        PAD_UB -= th_max;
        PAD_LB = theta.from() - theta.to();
        PAD_LB -= th_min;
        DebugOff(th_min.to_str(true) << endl);
        DebugOff(th_max.to_str(true) << endl);
    }
    else {
        DebugOff("Number of bus_pairs = " << bus_pairs.size() << endl);
        PAD_UB = vi.from()*vr.to() - vr.from()*vi.to();
        PAD_UB -= tan_th_max*(vr.from()*vr.to() + vi.from()*vi.to());
        
        PAD_LB = vi.from()*vr.to() - vr.from()*vi.to();
        PAD_LB -= tan_th_min*(vr.from()*vr.to() + vi.from()*vi.to());
        DebugOff(th_min.to_str(true) << endl);
        DebugOff(th_max.to_str(true) << endl);
    }
//    ACOPF->add(PAD_UB.in(bus_pairs) <= 0);
//    ACOPF->add(PAD_LB.in(bus_pairs) >= 0);
    
    
    /*  Thermal Limit Constraints */
    Constraint Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from += power(Pf_from, 2) + power(Qf_from, 2);
    Thermal_Limit_from -= power(S_max, 2);
    ACOPF->add(Thermal_Limit_from.in(arcs) <= 0);
    
    Constraint Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to += power(Pf_to, 2) + power(Qf_to, 2);
    Thermal_Limit_to -= power(S_max,2);
    ACOPF->add(Thermal_Limit_to.in(arcs) <= 0);
//    ACOPF->print_expanded();
    return ACOPF;
}

double PowerNet::solve_acopf(PowerModelType pmt, int output, double tol){
    
    auto ACOPF = build_ACOPF();
    bool relax;
    solver OPF(*ACOPF,ipopt);
    auto mipgap = 1e-6;
    OPF.run(output, relax = false, tol = 1e-6, mipgap, "ma27");
    return ACOPF->_obj_val;
}


void PowerNet::fill_wbnds(){
    double cos_max_, cos_min_, w_max_, w_min_, wr_max_, wr_min_, sin_max_, sin_min_, wi_max_, wi_min_;
    for(int i = 0; i < nodes.size()-1; i++) {
        for(int j = i+1; j < nodes.size(); j++) {
            Bus *bus_s = (Bus *) (nodes[i]);
            Bus *bus_d = (Bus *) (nodes[j]);
            
            if(get_arc(bus_s, bus_d)) continue;
            
            string name = bus_s->_name + "," + bus_d->_name;
            _bus_pairs_chord._keys.push_back(new index_pair(index_(bus_s->_name), index_(bus_d->_name)));
            
            if (m_theta_lb < -3.14 && m_theta_ub > 3.14) {
                cos_max_ = 1;
                cos_min_ = -1;
            } else if (m_theta_lb < 0 && m_theta_ub > 0) {
                cos_max_ = 1;
                cos_min_ = min(cos(m_theta_lb), cos(m_theta_ub));
            } else {
                cos_max_ = max(cos(m_theta_lb), cos(m_theta_ub));
                cos_min_ = min(cos(m_theta_lb), cos(m_theta_ub));
            }
            w_max_ = bus_s->vbound.max * bus_d->vbound.max;
            w_min_ = bus_s->vbound.min * bus_d->vbound.min;
            
            wr_max_ = cos_max_ * w_max_;
            if (cos_min_ < 0) wr_min_ = cos_min_ * w_max_;
            else wr_min_ = cos_min_ * w_min_;
            
            if (m_theta_lb < -1.57 && m_theta_ub > 1.57) {
                sin_max_ = 1;
                sin_min_ = -1;
            } else {
                sin_max_ = sin(m_theta_ub);
                sin_min_ = sin(m_theta_lb);
            }
            
            if (sin_max_ > 0) wi_max_ = sin_max_ * w_max_;
            else wi_max_ = sin_max_ * w_min_;
            if (sin_min_ > 0) wi_min_ = sin_min_ * w_min_;
            else wi_min_ = sin_min_ * w_max_;
            
            //            cout << "\nImaginary line, bounds: (" << wr_min_ << "," << wr_max_ << "); (" << wi_min_ << "," << wi_max_ << ")";
            
            wr_max.set_val(name, wr_max_);
            wr_min.set_val(name, wr_min_);
            wi_max.set_val(name, wi_max_);
            wi_min.set_val(name, wi_min_);
        }
    }
}



