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
    Sg_min.set_name("Sg_min");
    Sg_max.set_name("Sg_max");
    Sd.set_name("Sd");
    Y.set_name("Y");
    Ysh.set_name("Ysh");
    Ych.set_name("Ych");
    T.set_name("T");
    Sd.set_name("Sd");
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

indices PowerNet::gens_per_node() const{
    indices ids("gens_per_node");
    ids = indices(gens);
    ids._ids = make_shared<vector<vector<size_t>>>();
    ids._ids->resize(get_nb_active_nodes());
    string key;
    size_t inst = 0;
    for (auto n: nodes) {
        if (n->_active) {
            for(auto g: ((Bus*)n)->_gen){
                if (!g->_active) {
                    continue;
                }
                key = g->_name;
                auto it1 = ids._keys_map->find(key);
                if (it1 == ids._keys_map->end()){
                    throw invalid_argument("In function gen_ids(), unknown key.");
                }
                ids._ids->at(inst).push_back(it1->second);
            }
            inst++;
        }
    }
    return ids;
}

indices PowerNet::out_arcs_per_node() const{
    auto ids = indices(arcs);
    ids._name = "out_arcs_per_node";
    ids._ids = make_shared<vector<vector<size_t>>>();
    ids._ids->resize(get_nb_active_nodes());
    string key;
    size_t inst = 0;
    for (auto n: nodes) {
        if (n->_active) {
            for (auto a:n->get_out()) {
                if (!a->_active) {
                    continue;
                }
                key = a->_name;
                auto it1 = ids._keys_map->find(key);
                if (it1 == ids._keys_map->end()){
                    throw invalid_argument("In function out_arcs_per_node(), unknown key.");
                }
                ids._ids->at(inst).push_back(it1->second);
            }
            inst++;
        }
    }
    return ids;
}

indices PowerNet::in_arcs_per_node() const{
    auto ids = indices(arcs);
    ids._name = "in_arcs_per_node";
    ids._ids = make_shared<vector<vector<size_t>>>();
    ids._ids->resize(get_nb_active_nodes());
    string key;
    size_t inst = 0;
    for (auto n: nodes) {
        if (n->_active) {
            for (auto a:n->get_in()) {
                if (!a->_active) {
                    continue;
                }
                key = a->_name;
                auto it1 = ids._keys_map->find(key);
                if (it1 == ids._keys_map->end()){
                    throw invalid_argument("In function in_arcs_per_node(), unknown key.");
                }
                ids._ids->at(inst).push_back(it1->second);
            }
            inst++;
        }
    }
    return ids;
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
        if (bp->_active) {
            nb++;
        }
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

int PowerNet::readGAMS(const string& fname) {
    double pi = 4.*atan(1.);
    string name;
//    double kvb = 0;
    //    int id = 0;
//    unsigned index = 0;
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
        pl.add_val(name,b->pl()/bMVA);
        ql.add_val(name,b->ql()/bMVA);
        gs.add_val(name,b->gs()/bMVA);
        bs.add_val(name,b->bs()/bMVA);
        v_max.add_val(name,b->vmax());
        v_min.add_val(name,b->vmin());
        w_min.add_val(name,pow(b->vmin(), 2));
        w_max.add_val(name,pow(b->vmax(), 2));
    }
    DebugOn("Gen Data : "<<endl);
    for (auto gen:gens) {
        gen->print();
        name = gen->_name;
        qg_max.add_val(name,gen->_qbound.max);
        qg_min.add_val(name,gen->_qbound.min);
        pg_max.add_val(name,gen->_pbound.max);
        pg_min.add_val(name,gen->_pbound.min);
        c0.add_val(name, gen->_cost->c0);
        c1.add_val(name, gen->_cost->c1*bMVA);
        c2.add_val(name, gen->_cost->c2*pow(bMVA,2));
        p_factor.add_val(name, gen->_part_factor/tot_part_fact);
//        p_factor.add_val(name, gen->_part_factor);
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
        g.add_val(name,arc->g);
        b.add_val(name,arc->b);
        tr.add_val(name,arc->tr);
        as.add_val(name,arc->as);
        arc->cc = arc->tr*cos(arc->as); // Rectangular values for transformer phase shifters
        arc->dd = arc->tr*sin(arc->as);
        g_ff.add_val(name,arc->g/(pow(arc->cc, 2) + pow(arc->dd, 2)));
        g_ft.add_val(name,(-arc->g*arc->cc + arc->b*arc->dd)/(pow(arc->cc, 2) + pow(arc->dd, 2)));
        g_tt.add_val(name,arc->g);
        g_tf.add_val(name,(-arc->g*arc->cc - arc->b*arc->dd)/(pow(arc->cc, 2) + pow(arc->dd, 2)));
        b_ff.add_val(name,(arc->ch*0.5 + arc->b)/(pow(arc->cc, 2) + pow(arc->dd, 2)));
        b_ft.add_val(name,(-arc->b*arc->cc - arc->g*arc->dd)/(pow(arc->cc, 2) + pow(arc->dd, 2)));
        b_tt.add_val(name,(arc->ch*0.5 + arc->b));
        b_tf.add_val(name,(-arc->b*arc->cc + arc->g*arc->dd)/(pow(arc->cc, 2) + pow(arc->dd, 2)));
        ch.add_val(name,arc->ch);
        S_max.add_val(name,arc->limit);
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
//            _bus_pairs._keys.push_back(new index_pair(index_(src), index_(dest), arc->_active));
            bus_pair_names.insert(name);
        }
        if (!arc->_parallel) {
            th_min.add_val(name,arc->tbound.min);
            th_max.add_val(name,arc->tbound.max);
            tan_th_min.add_val(name,tan(arc->tbound.min));
            tan_th_max.add_val(name,tan(arc->tbound.max));
            
        }
        else {
            th_min.add_val(name,gravity::max(th_min.eval(name), arc->tbound.min));
            th_max.add_val(name,gravity::min(th_max.eval(name), arc->tbound.max));
            tan_th_min.add_val(name,tan(th_min.eval(name)));
            tan_th_max.add_val(name,tan(th_max.eval(name)));
        }
        if (arc->tbound.min >= 0) {
            wr_max.add_val(name,bus_s->vbound.max*bus_d->vbound.max*cos(th_min.eval(name)));
            wr_min.add_val(name,bus_s->vbound.min*bus_d->vbound.min*cos(th_max.eval(name)));
            wi_max.add_val(name,bus_s->vbound.max*bus_d->vbound.max*sin(th_max.eval(name)));
            wi_min.add_val(name,bus_s->vbound.min*bus_d->vbound.min*sin(th_min.eval(name)));
        };
        if (arc->tbound.max <= 0) {
            wr_max.add_val(name,bus_s->vbound.max*bus_d->vbound.max*cos(th_max.eval(name)));
            wr_min.add_val(name,bus_s->vbound.min*bus_d->vbound.min*cos(th_min.eval(name)));
            wi_max.add_val(name,bus_s->vbound.min*bus_d->vbound.min*sin(th_max.eval(name)));
            wi_min.add_val(name,bus_s->vbound.max*bus_d->vbound.max*sin(th_min.eval(name)));
        }
        if (arc->tbound.min < 0 && arc->tbound.max > 0) {
            wr_max.add_val(name,bus_s->vbound.max*bus_d->vbound.max);
            wr_min.add_val(name,bus_s->vbound.min*bus_d->vbound.min*gravity::min(cos(th_min.eval(name)), cos(th_max.eval(name))));
            wi_max.add_val(name,bus_s->vbound.max*bus_d->vbound.max*sin(th_max.eval(name)));
            wi_min.add_val(name,bus_s->vbound.max*bus_d->vbound.max*sin(th_min.eval(name)));
        }
        cphi.add_val(name, cos(0.5*(arc->tbound.min+arc->tbound.max)));
        sphi.add_val(name, sin(0.5*(arc->tbound.min+arc->tbound.max)));
        cos_d.add_val(name, cos(0.5*(arc->tbound.max-arc->tbound.min)));
        getline(file, word,'\n');
        file >> word;
    }
    DebugOff(ch.to_str(true) << endl);
    DebugOff(as.to_str(true) << endl);
    DebugOff(tr.to_str(true) << endl);
    file.close();
    return 0;
}


int PowerNet::readgrid(const string& fname, bool reverse_arcs) {
    double pi = 4.*atan(1.);
    string name;
    double kvb = 0;
//    int id = 0;
    unsigned index = 0;
    cout << "Loading file " << fname << endl;
    ifstream file(fname.c_str(), std::ifstream::in);
    if(!file.is_open()) {
        throw invalid_argument("Could not open file " + fname);
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
    double total_p_load = 0, total_q_load = 0;
    while(word.compare("];")) {
        name = word.c_str();
        file >> ws >> word;
        status = atoi(word.c_str());
        if (status==3) {
            ref_bus = name;
            DebugOn("Ref bus = " << ref_bus << endl);
        }
        file >> ws >> word;
        pl.add_val(name,atof(word.c_str())/bMVA);
        file >> word;
        ql.add_val(name,atof(word.c_str())/bMVA);
        file >> word;
        gs.add_val(name,atof(word.c_str())/bMVA);
        file >> word;
        bs.add_val(name,atof(word.c_str())/bMVA);
        Ysh.add_val(name, Cpx(gs.eval(),bs.eval()));
        file >> ws >> word >> ws >> word;
        v_s.add_val(name,atof(word.c_str()));
        file >> ws >> word >> ws >> word;
        kvb = atof(word.c_str());
        file >> ws >> word >> ws >> word;
        v_max.add_val(name,atof(word.c_str()));
        getline(file, word,';');
        v_min.add_val(name,atof(word.c_str()));
        w_min.add_val(name,pow(v_min.eval(), 2));
        w_max.add_val(name,pow(v_max.eval(), 2));
        // single phase

        bus = new Bus(name, pl.eval(), ql.eval(), gs.eval(), bs.eval(), v_min.eval(), v_max.eval(), kvb, 1);
//        bus_clone = new Bus(name, pl.eval(), ql.eval(), gs.eval(), bs.eval(), v_min.eval(), v_max.eval(), kvb, 1);
        total_p_load += pl.eval();
        total_q_load += ql.eval();
        bus->vs = v_s.eval();
        V_min.add_val(name,Cpx(-1*bus->vbound.max, -1*bus->vbound.max));
        V_max.add_val(name,Cpx(bus->vbound.max, bus->vbound.max));
        Sd.add_val(name,Cpx(bus->_cond[0]->_pl,bus->_cond[0]->_ql));
//        bus_clone->vs = v_s.eval();
        if (status>=4) {
            bus->_active = false;
//            bus_clone->_active = false;
        }

        this->Net::add_node(bus);
        if (status>=4) {
            DebugOn("INACTIVE NODE!\n" << name << endl);
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
        pg_s.add_val(name,atof(word.c_str())/bMVA);
        file >> word;
        qg_s.add_val(name,atof(word.c_str())/bMVA);
        file >> word;
        qg_max.add_val(name,atof(word.c_str())/bMVA);
        file >> word;
        qg_min.add_val(name,atof(word.c_str())/bMVA);

        file >> ws >> word >> ws >> word >> ws >> word;
        status = atoi(word.c_str());
        file >> word;
        pg_max.add_val(name,atof(word.c_str())/bMVA);

        
        file >> word;
        pg_min.add_val(name,atof(word.c_str())/bMVA);
        getline(file, word,'\n');
//        gen_status.push_back(status==1);


        bus->_has_gen = true;
        /** generator name, ID */
        Gen* g = new Gen(bus, name, pg_min.eval(index), pg_max.eval(index), qg_min.eval(index), qg_max.eval(index));
        g->_id = index;
        g->_ps = pg_s.eval();
        g->_qs = qg_s.eval();
        Sg_min.add_val(name,Cpx(g->_pbound.min,g->_qbound.min));
        Sg_max.add_val(name,Cpx(g->_pbound.max,g->_qbound.max));
        
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
    for (size_t i = 0; i < gens.size(); ++i) {
        file >> ws >> word >> ws >> word >> ws >> word >> ws >> word >> ws >> word;
        c2.add_val(to_string(i),atof(word.c_str())*pow(bMVA,2));
        file >> word;
        c1.add_val(to_string(i),atof(word.c_str())*bMVA);
        file >> word;
        c0.add_val(to_string(i),atof(word.c_str()));
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
//        if(get_node(src)->_id > get_node(dest)->_id) {//Reverse arc direction
        if((reverse_arcs && get_node(src)->_id > get_node(dest)->_id) || arcID.find(key)!=arcID.end()) {//Reverse arc direction
            Warning("Adding arc linking " +src+" and "+dest);
            Warning(" with reversed direction, reversing source and destination.\n");
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
        arc->as = (atof(word.c_str())*pi)/180.;
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
            arc->tbound.min = -60.*pi/180.;
            arc->tbound.max = 60.*pi/180.;
            
        }
        if (reversed && reverse_arcs) {
            arc->g /= pow(arc->tr,2);
            arc->b /= pow(arc->tr,2);
            arc->ch /= pow(arc->tr,2);
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
        
        arc->smax = gravity::max(
                        pow(bus_s->vbound.max,2)*(arc->g*arc->g + arc->b*arc->b)*(pow(bus_s->vbound.max,2) + pow(bus_d->vbound.max,2)),
                        pow(bus_d->vbound.max,2)*(arc->g*arc->g+arc->b*arc->b)*(pow(bus_d->vbound.max,2) + pow(bus_s->vbound.max,2))
                        );
        name = arc->_name;
        g.add_val(name,arc->g);
        b.add_val(name,arc->b);
        tr.add_val(name,arc->tr);
        as.add_val(name,arc->as);
        //(g+g_fr)/tm^2
        g_ff.add_val(name,arc->g/(pow(arc->cc, 2) + pow(arc->dd, 2)));
        g_ft.add_val(name,(-arc->g*arc->cc + arc->b*arc->dd)/(pow(arc->cc, 2) + pow(arc->dd, 2)));
        
        g_tt.add_val(name,arc->g);
        g_tf.add_val(name,(-arc->g*arc->cc - arc->b*arc->dd)/(pow(arc->cc, 2) + pow(arc->dd, 2)));
        
        
        b_ff.add_val(name,(arc->ch*0.5 + arc->b)/(pow(arc->cc, 2) + pow(arc->dd, 2)));
        b_ft.add_val(name,(-arc->b*arc->cc - arc->g*arc->dd)/(pow(arc->cc, 2) + pow(arc->dd, 2)));
        
        b_tt.add_val(name,(arc->ch*0.5 + arc->b));
        b_tf.add_val(name,(-arc->b*arc->cc + arc->g*arc->dd)/(pow(arc->cc, 2) + pow(arc->dd, 2)));
        
        ch.add_val(name,arc->ch);
//        S_max.add_val(name,gravity::min(arc->limit,max(2.*total_p_load, 2.*total_q_load)));
        S_max.add_val(name,arc->limit);
        
        /* Complex params */
        Y.add_val(name,Cpx(arc->g,arc->b));
        Ych.add_val(name,Cpx(0,arc->ch));
        T.add_val(name,Cpx(arc->tr, 0));
        
        if(arc->status != 1 || !bus_s->_active || !bus_d->_active) {
            arc->_active = false;
            DebugOn("INACTIVE ARC!\n" << arc->_name << endl);
        }
        arc->connect();
        add_arc(arc);
        /* Switching to bus_pairs keys */
        name = bus_s->_name + "," + bus_d->_name;
        if (arc->_active && bus_pair_names.count(name)==0) {
//            _bus_pairs._keys.push_back(new index_pair(index_(bus_s->_name), index_(bus_d->_name), arc->_active));
            bus_pair_names.insert(name);
        }
        if (!arc->_parallel) {
            th_min.add_val(name,arc->tbound.min);
            th_max.add_val(name,arc->tbound.max);
            tan_th_min.add_val(name,tan(arc->tbound.min));
            tan_th_max.add_val(name,tan(arc->tbound.max));
            
        }
        else {
            th_min.set_val(name,gravity::max(th_min.eval(name), arc->tbound.min));
            th_max.set_val(name,gravity::min(th_max.eval(name), arc->tbound.max));
            tan_th_min.set_val(name,tan(th_min.eval(name)));
            tan_th_max.set_val(name,tan(th_max.eval(name)));
        }
        if (arc->tbound.min >= 0) {
            wr_max.add_val(name,bus_s->vbound.max*bus_d->vbound.max*cos(th_min.eval(name)));
            wr_min.add_val(name,bus_s->vbound.min*bus_d->vbound.min*cos(th_max.eval(name)));
            wi_max.add_val(name,bus_s->vbound.max*bus_d->vbound.max*sin(th_max.eval(name)));
            wi_min.add_val(name,bus_s->vbound.min*bus_d->vbound.min*sin(th_min.eval(name)));
        };
        if (arc->tbound.max <= 0) {
            wr_max.add_val(name,bus_s->vbound.max*bus_d->vbound.max*cos(th_max.eval(name)));
            wr_min.add_val(name,bus_s->vbound.min*bus_d->vbound.min*cos(th_min.eval(name)));
            wi_max.add_val(name,bus_s->vbound.min*bus_d->vbound.min*sin(th_max.eval(name)));
            wi_min.add_val(name,bus_s->vbound.max*bus_d->vbound.max*sin(th_min.eval(name)));
        }
        if (arc->tbound.min < 0 && arc->tbound.max > 0) {
            wr_max.add_val(name,bus_s->vbound.max*bus_d->vbound.max);
            wr_min.add_val(name,bus_s->vbound.min*bus_d->vbound.min*gravity::min(cos(th_min.eval(name)), cos(th_max.eval(name))));
            wi_max.add_val(name,bus_s->vbound.max*bus_d->vbound.max*sin(th_max.eval(name)));
            wi_min.add_val(name,bus_s->vbound.max*bus_d->vbound.max*sin(th_min.eval(name)));
        }
        cphi.add_val(name, cos(0.5*(arc->tbound.min+arc->tbound.max)));
        sphi.add_val(name, sin(0.5*(arc->tbound.min+arc->tbound.max)));
        cos_d.add_val(name, cos(0.5*(arc->tbound.max-arc->tbound.min)));
        getline(file, word,'\n');
        file >> word;
    }
    DebugOff(ch.to_str(true) << endl);
    DebugOff(as.to_str(true) << endl);
    DebugOff(tr.to_str(true) << endl);
    
    file.close();
//    if (nodes.size()>1000) {
//        add_3d_nlin = false;
//    }
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
                cos_min_ = gravity::min(cos(m_theta_lb), cos(m_theta_ub));
            } else{
                cos_max_ = gravity::max(cos(m_theta_lb),cos(m_theta_ub));
                cos_min_ = gravity::min(cos(m_theta_lb), cos(m_theta_ub));
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

            wr_max.add_val(name,wr_max_);
            wr_min.add_val(name,wr_min_);
            wi_max.add_val(name,wi_max_);
            wi_min.add_val(name,wi_min_);
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


shared_ptr<Model<>> PowerNet::build_SCOPF(PowerModelType pmt, int output, double tol){
    auto bus_pairs = get_bus_pairs();
    /** MODEL DECLARATION */
    shared_ptr<Model<>> SOCPF(new Model<>("SCOPF Model"));
    /** Variables */
    /* power generation variables */
    var<double> Pg("Pg", pg_min, pg_max);
    var<double> Qg ("Qg", qg_min, qg_max);
    SOCPF->add(Pg.in(gens));
    SOCPF->add(Qg.in(gens));


    /* power flow variables */
    var<double> Pf_from("Pf_from", -1*S_max,S_max);
    var<double> Qf_from("Qf_from", -1*S_max,S_max);
    var<double> Pf_to("Pf_to", -1*S_max,S_max);
    var<double> Qf_to("Qf_to", -1*S_max,S_max);
    SOCPF->add(Pf_from.in(arcs));
    SOCPF->add(Qf_from.in(arcs));
    SOCPF->add(Pf_to.in(arcs));
    SOCPF->add(Qf_to.in(arcs));

    /* Real part of Wij = ViVj */
    var<double>  R_Wij("R_Wij", wr_min, wr_max);
    /* Imaginary part of Wij = ViVj */
    var<double>  Im_Wij("Im_Wij", wi_min, wi_max);
    /* Magnitude of Wii = Vi^2 */
    var<double>  Wii("Wii", w_min, w_max);
    SOCPF->add(Wii.in(nodes));
    SOCPF->add(R_Wij.in(bus_pairs));
    SOCPF->add(Im_Wij.in(bus_pairs));

    /* Initialize variables */
    R_Wij.initialize_all(1.0);
    Wii.initialize_all(1.001);
    
    /** Sets */
    auto gens = gens_per_node();
    auto out_arcs = out_arcs_per_node();
    auto in_arcs = in_arcs_per_node();

    /**  Objective */
    auto obj = c1.tr()*Pg + c2.tr()*pow(Pg,2) + sum(c0);
    SOCPF->min(obj);

    /** Constraints */
    /* Second-order cone constraints */
    Constraint<> SOC("SOC");
    SOC = pow(R_Wij, 2) + pow(Im_Wij, 2) - Wii.from(bus_pairs)*Wii.to(bus_pairs);
    SOCPF->add(SOC.in(bus_pairs) <= 0);

    /* Flow conservation */
    Constraint<> KCL_P("KCL_P");
    KCL_P  = sum(Pf_from, out_arcs) + sum(Pf_to, in_arcs) + pl - sum(Pg, gens) + gs*Wii;
    SOCPF->add(KCL_P.in(nodes) == 0);

    Constraint<> KCL_Q("KCL_Q");
    KCL_Q  = sum(Qf_from, out_arcs) + sum(Qf_to, in_arcs) + ql - sum(Qg, gens) - bs*Wii;
    SOCPF->add(KCL_Q.in(nodes) == 0);

    /* AC Power Flow */
    Constraint<> Flow_P_From("Flow_P_From");
    Flow_P_From = Pf_from - (g_ff*Wii.from(arcs) + g_ft*R_Wij + b_ft*Im_Wij);
    SOCPF->add(Flow_P_From.in(arcs) == 0);

    Constraint<> Flow_P_To("Flow_P_To");
    Flow_P_To = Pf_to - (g_tt*Wii.to(arcs) + g_tf*R_Wij - b_tf*Im_Wij);
    SOCPF->add(Flow_P_To.in(arcs) == 0);

    Constraint<> Flow_Q_From("Flow_Q_From");
    Flow_Q_From = Qf_from - (g_ft*Im_Wij - b_ff*Wii.from(arcs) - b_ft*R_Wij);
    SOCPF->add(Flow_Q_From.in(arcs) == 0);

    Constraint<> Flow_Q_To("Flow_Q_To");
    Flow_Q_To = Qf_to + (b_tt*Wii.to(arcs) + b_tf*R_Wij + g_tf*Im_Wij);
    SOCPF->add(Flow_Q_To.in(arcs) == 0);

    /* Phase Angle Bounds constraints */
    Constraint<> PAD_UB("PAD_UB");
    PAD_UB = Im_Wij;
    PAD_UB <= tan_th_max*R_Wij;
    SOCPF->add(PAD_UB);

    Constraint<> PAD_LB("PAD_LB");
    PAD_LB =  Im_Wij;
    PAD_LB >= tan_th_min*R_Wij;
    SOCPF->add(PAD_LB);

    /* Thermal Limit Constraints */
    Constraint<> Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from = pow(Pf_from, 2) + pow(Qf_from, 2);
    Thermal_Limit_from <= pow(S_max,2);
    SOCPF->add(Thermal_Limit_from);


    Constraint<> Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to = pow(Pf_to, 2) + pow(Qf_to, 2);
    Thermal_Limit_to <= pow(S_max,2);
    SOCPF->add(Thermal_Limit_to);

    /* Lifted Nonlinear Cuts */
    Constraint<> LNC1("LNC1");
    LNC1 += (v_min.from(bus_pairs)+v_max.from(bus_pairs))*(v_min.to(bus_pairs)+v_max.to(bus_pairs))*(sphi*Im_Wij + cphi*R_Wij);
    LNC1 -= v_max.to(bus_pairs)*cos_d*(v_min.to(bus_pairs)+v_max.to(bus_pairs))*Wii.from(bus_pairs);
    LNC1 -= v_max.from(bus_pairs)*cos_d*(v_min.from(bus_pairs)+v_max.from(bus_pairs))*Wii.to(bus_pairs);
    LNC1 -= v_max.from(bus_pairs)*v_max.to(bus_pairs)*cos_d*(v_min.from(bus_pairs)*v_min.to(bus_pairs) - v_max.from(bus_pairs)*v_max.to(bus_pairs));
    SOCPF->add(LNC1.in(bus_pairs) >= 0);

    Constraint<> LNC2("LNC2");
    LNC2 += (v_min.from(bus_pairs)+v_max.from(bus_pairs))*(v_min.to(bus_pairs)+v_max.to(bus_pairs))*(sphi*Im_Wij + cphi*R_Wij);
    LNC2 -= v_min.to(bus_pairs)*cos_d*(v_min.to(bus_pairs)+v_max.to(bus_pairs))*Wii.from(bus_pairs);
    LNC2 -= v_min.from(bus_pairs)*cos_d*(v_min.from(bus_pairs)+v_max.from(bus_pairs))*Wii.to(bus_pairs);
    LNC2 += v_min.from(bus_pairs)*v_min.to(bus_pairs)*cos_d*(v_min.from(bus_pairs)*v_min.to(bus_pairs) - v_max.from(bus_pairs)*v_max.to(bus_pairs));
    SOCPF->add(LNC2.in(bus_pairs) >= 0);
    return SOCPF;
}



shared_ptr<Model<>> PowerNet::build_ACOPF(PowerModelType pmt, int output, double tol){
    auto gen_nodes = gens_per_node();
    auto out_arcs = out_arcs_per_node();
    auto in_arcs = in_arcs_per_node();
    bool polar = (pmt==ACPOL);
    if (polar) {
        Debug("Using polar model\n");
    }
    else {
        Debug("Using rectangular model\n");
    }
    shared_ptr<Model<>> ACOPF(new Model<>("AC-OPF Model"));
    /** Variables */
    /* Power generation variables */
    var<double> Pg("Pg", pg_min, pg_max);
    var<double> Qg ("Qg", qg_min, qg_max);
    ACOPF->add(Pg.in(gens), Qg.in(gens));
    Pg_base = Pg;
    Qg_base = Qg;
    
    /* Power flow variables */
    var<double> Pf_from("Pf_from", -1*S_max,S_max);
    var<double> Qf_from("Qf_from", -1*S_max,S_max);
    var<double> Pf_to("Pf_to", -1*S_max,S_max);
    var<double> Qf_to("Qf_to", -1*S_max,S_max);
    ACOPF->add(Pf_from.in(arcs), Qf_from.in(arcs),Pf_to.in(arcs),Qf_to.in(arcs));
    p_from_base = Pf_from;
    p_to_base = Pf_to;
    q_from_base = Qf_from;
    q_to_base = Qf_to;
    /** Voltage related variables */
    var<double> theta("theta", -6.3, 6.3);
    var<double> v("|V|", v_min, v_max);
    //    var<double> vr("vr");
    //    var<double> vi("vi");
    var<double> vr("vr", -1*v_max,v_max);
    var<double> vi("vi", -1*v_max,v_max);
    
    var<> v_from, v_to, theta_from, theta_to;
    var<> vr_from, vr_to, vi_from, vi_to;
    if (polar) {
        ACOPF->add(v.in(nodes));
        ACOPF->add(theta.in(nodes));
        v.initialize_all(1.0);
        v_from = v.from(arcs);
        v_to = v.to(arcs);
        theta_from = theta.from(arcs);
        theta_to = theta.to(arcs);        
    }
    else {
        ACOPF->add_var(vr.in(nodes));
        ACOPF->add_var(vi.in(nodes));
        vr.initialize_all(1.0);
        vr_from = vr.from(arcs);
        vr_to = vr.to(arcs);
        vi_from = vi.from(arcs);
        vi_to = vi.to(arcs);
    }
    v_base = v;
    theta_base = theta;
    
    /** Construct the objective function */
    auto obj = c1.tr()*Pg + c2.tr()*pow(Pg.vec(),2) + sum(c0);
    ACOPF->min(obj);
    
    /** Define constraints */
    
    /* REF BUS */
    Constraint<> Ref_Bus("Ref_Bus");
    if (polar) {
        Ref_Bus = theta(ref_bus);
    }
    else {
        Ref_Bus = vi(ref_bus);
    }
    ACOPF->add(Ref_Bus == 0);
    
    /** KCL Flow conservation */
    Constraint<> KCL_P("KCL_P");
    Constraint<> KCL_Q("KCL_Q");
    KCL_P  = sum(Pf_from, out_arcs) + sum(Pf_to, in_arcs) + pl - sum(Pg, gen_nodes);
    KCL_Q  = sum(Qf_from, out_arcs) + sum(Qf_to, in_arcs) + ql - sum(Qg, gen_nodes);
    /* Shunts */
    if (polar) {
        KCL_P +=  gs*pow(v,2);
        KCL_Q -=  bs*pow(v,2);
    }
    else {
        KCL_P +=  gs*(pow(vr,2)+pow(vi,2));
        KCL_Q -=  bs*(pow(vr,2)+pow(vi,2));
    }
    ACOPF->add(KCL_P.in(nodes) == 0);
    ACOPF->add(KCL_Q.in(nodes) == 0);
    
    /** AC Power Flows */
    /** TODO write the constraints in Complex form */
    Constraint<> Flow_P_From("Flow_P_From");
    Flow_P_From += Pf_from;
    if (polar) {
        Flow_P_From -= g/pow(tr,2)*pow(v_from,2);
        Flow_P_From += g/tr*(v_from*v_to*cos(theta_from - theta_to - as));
        Flow_P_From += b/tr*(v_from*v_to*sin(theta_from - theta_to - as));
    }
    else {
        Flow_P_From -= g_ff*(pow(vr_from, 2) + pow(vi_from, 2));
        Flow_P_From -= g_ft*(vr_from*vr_to + vi_from*vi_to);
        Flow_P_From -= b_ft*(vi_from*vr_to - vr_from*vi_to);
    }
    ACOPF->add(Flow_P_From.in(arcs)==0);
    
    Constraint<> Flow_P_To("Flow_P_To");
    Flow_P_To += Pf_to;
    if (polar) {
        Flow_P_To -= g*pow(v_to, 2);
        Flow_P_To += g/tr*(v_from*v_to*cos(theta_to - theta_from + as));
        Flow_P_To += b/tr*(v_from*v_to*sin(theta_to - theta_from + as));
    }
    else {
        Flow_P_To -= g_tt*(pow(vr_to, 2) + pow(vi_to, 2));
        Flow_P_To -= g_tf*(vr_from*vr_to + vi_from*vi_to);
        Flow_P_To -= b_tf*(vi_to*vr_from - vr_to*vi_from);
    }
    ACOPF->add(Flow_P_To.in(arcs)==0);
    
    Constraint<> Flow_Q_From("Flow_Q_From");
    Flow_Q_From += Qf_from;
    if (polar) {
        Flow_Q_From += (0.5*ch+b)/pow(tr,2)*pow(v_from,2);
        Flow_Q_From -= b/tr*(v_from*v_to*cos(theta_from - theta_to - as));
        Flow_Q_From += g/tr*(v_from*v_to*sin(theta_from - theta_to - as));
    }
    else {
        Flow_Q_From += b_ff*(pow(vr_from, 2) + pow(vi_from, 2));
        Flow_Q_From += b_ft*(vr_from*vr_to + vi_from*vi_to);
        Flow_Q_From -= g_ft*(vi_from*vr_to - vr_from*vi_to);
    }
    ACOPF->add(Flow_Q_From.in(arcs)==0);
    Constraint<> Flow_Q_To("Flow_Q_To");
    Flow_Q_To += Qf_to;
    if (polar) {
        Flow_Q_To += (0.5*ch+b)*pow(v_to,2);
        Flow_Q_To -= b/tr*(v_from*v_to*cos(theta_to - theta_from + as));
        Flow_Q_To += g/tr*(v_from*v_to*sin(theta_to - theta_from + as));
    }
    else {
        Flow_Q_To += b_tt*(pow(vr_to, 2) + pow(vi_to, 2));
        Flow_Q_To += b_tf*(vr_from*vr_to + vi_from*vi_to);
        Flow_Q_To -= g_tf*(vi_to*vr_from - vr_to*vi_from);
    }
    ACOPF->add(Flow_Q_To.in(arcs)==0);
    
    /** AC voltage limit constraints. */
    if (!polar) {
        Constraint<> Vol_limit_UB("Vol_limit_UB");
        Vol_limit_UB = pow(vr, 2) + pow(vi, 2);
        Vol_limit_UB -= pow(v_max, 2);
        ACOPF->add(Vol_limit_UB.in(nodes) <= 0);
        
        Constraint<> Vol_limit_LB("Vol_limit_LB");
        Vol_limit_LB = pow(vr, 2) + pow(vi, 2);
        Vol_limit_LB -= pow(v_min,2);
        ACOPF->add(Vol_limit_LB.in(nodes) >= 0);
        DebugOff(v_min.to_str(true) << endl);
        DebugOff(v_max.to_str(true) << endl);
    }
    
    
    /* Phase Angle Bounds constraints */
    Constraint<> PAD_UB("PAD_UB");
    Constraint<> PAD_LB("PAD_LB");
    auto bus_pairs = get_bus_pairs();
    if (polar) {
        PAD_UB = theta.from(bus_pairs) - theta.to(bus_pairs);
        PAD_UB -= th_max;
        PAD_LB = theta.from(bus_pairs) - theta.to(bus_pairs);
        PAD_LB -= th_min;
    }
    else {
        DebugOff("Number of bus_pairs = " << bus_pairs.size() << endl);
        PAD_UB = vi.from(bus_pairs)*vr.to(bus_pairs) - vr.from(bus_pairs)*vi.to(bus_pairs);
        PAD_UB -= tan_th_max*(vr.from(bus_pairs)*vr.to(bus_pairs) + vi.from(bus_pairs)*vi.to(bus_pairs));
        
        PAD_LB = vi.from(bus_pairs)*vr.to(bus_pairs) - vr.from(bus_pairs)*vi.to(bus_pairs);
        PAD_LB -= tan_th_min*(vr.from(bus_pairs)*vr.to(bus_pairs) + vi.from(bus_pairs)*vi.to(bus_pairs));
    }
    ACOPF->add(PAD_UB.in(bus_pairs) <= 0);
    ACOPF->add(PAD_LB.in(bus_pairs) >= 0);
    
    
    /*  Thermal Limit Constraints */
    Constraint<> Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from += pow(Pf_from, 2) + pow(Qf_from, 2);
    Thermal_Limit_from -= pow(S_max, 2);
    ACOPF->add(Thermal_Limit_from.in(arcs) <= 0);
    
    Constraint<> Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to += pow(Pf_to, 2) + pow(Qf_to, 2);
    Thermal_Limit_to -= pow(S_max,2);
    ACOPF->add(Thermal_Limit_to.in(arcs) <= 0);
    return ACOPF;
}



/** Return the vector of arcs of the chordal completion ignoring parallel lines **/
indices PowerNet::get_bus_pairs_chord(){
    set<pair<Node*,Node*>> unique_pairs;
    indices bpairs("bus_pairs");
    for (auto a: arcs) {
        if (!a->_parallel) {
            unique_pairs.insert({a->_src,a->_dest});
            bpairs.add(a->_src->_name+","+a->_dest->_name);
        }
    }
    string key;
    double cos_max_, cos_min_, sin_max_, sin_min_;
    double wr_max_, wr_min_, wi_max_, wi_min_, w_max_, w_min_;
    if (m_theta_lb < -3.14 && m_theta_ub > 3.14) {
        cos_max_ = 1;
        cos_min_ = -1;
    } else if (m_theta_lb < 0 && m_theta_ub > 0){
        cos_max_ = 1;
        cos_min_ = gravity::min(cos(m_theta_lb), cos(m_theta_ub));
    } else{
        cos_max_ = gravity::max(cos(m_theta_lb),cos(m_theta_ub));
        cos_min_ = gravity::min(cos(m_theta_lb), cos(m_theta_ub));
    }
    if(m_theta_lb < -1.57 && m_theta_ub > 1.57){
        sin_max_ = 1;
        sin_min_ = -1;
    } else{
        sin_max_ = sin(m_theta_ub);
        sin_min_ = sin(m_theta_lb);
    }
    for (auto &bag: _bags) {
        for (size_t i = 0; i< bag.size()-1; i++) {
            if (unique_pairs.insert({bag[i],bag[i+1]}).second) {
                auto bus_s = (Bus*)bag[i];
                auto bus_d = (Bus*)bag[i+1];
                w_max_ = bus_s->vbound.max*bus_d->vbound.max;
                w_min_ = bus_s->vbound.min*bus_d->vbound.min;
                wr_max_ = cos_max_*w_max_;
                if(cos_min_ < 0) wr_min_ = cos_min_*w_max_;
                else wr_min_ = cos_min_*w_min_;
                if(sin_max_ > 0) wi_max_ = sin_max_*w_max_;
                else wi_max_ = sin_max_*w_min_;
                if(sin_min_ > 0) wi_min_ = sin_min_*w_min_;
                else wi_min_ = sin_min_*w_max_;
                auto name = bag[i]->_name + "," + bag[i+1]->_name;
                wr_max.add_val(name,wr_max_);
                wr_min.add_val(name,wr_min_);
                wi_max.add_val(name,wi_max_);
                wi_min.add_val(name,wi_min_);
                bpairs.add(name);
            }
        }
        /* Loop back pair */
        if (unique_pairs.insert({bag[0],bag[bag.size()-1]}).second) {
            auto name = bag[0]->_name + "," + bag[bag.size()-1]->_name;
            auto bus_s = (Bus*)bag[0];
            auto bus_d = (Bus*)bag[bag.size()-1];
            w_max_ = bus_s->vbound.max*bus_d->vbound.max;
            w_min_ = bus_s->vbound.min*bus_d->vbound.min;
            wr_max_ = cos_max_*w_max_;
            if(cos_min_ < 0) wr_min_ = cos_min_*w_max_;
            else wr_min_ = cos_min_*w_min_;
            if(sin_max_ > 0) wi_max_ = sin_max_*w_max_;
            else wi_max_ = sin_max_*w_min_;
            if(sin_min_ > 0) wi_min_ = sin_min_*w_min_;
            else wi_min_ = sin_min_*w_max_;
            wr_max.add_val(name,wr_max_);
            wr_min.add_val(name,wr_min_);
            wi_max.add_val(name,wi_max_);
            wi_min.add_val(name,wi_min_);
            bpairs.add(name);
        }
    }
    return bpairs;
}

double PowerNet::solve_acopf(PowerModelType pmt, int output, double tol){
    
    auto ACOPF = build_ACOPF();
    bool relax;
    solver<> OPF(ACOPF,ipopt);
//    auto mipgap = 1e-6;
    OPF.run(output, relax = false, tol);
    return ACOPF->_obj->get_val();
}


void PowerNet::fill_wbnds(){
    double cos_max_, cos_min_, w_max_, w_min_, wr_max_, wr_min_, sin_max_, sin_min_, wi_max_, wi_min_;
    for(int i = 0; i < nodes.size()-1; i++) {
        for(int j = i+1; j < nodes.size(); j++) {
            Bus *bus_s = (Bus *) (nodes[i]);
            Bus *bus_d = (Bus *) (nodes[j]);
            
            if(get_arc(bus_s, bus_d)) continue;
            
            string name = bus_s->_name + "," + bus_d->_name;
//            _bus_pairs_chord._keys.push_back(new index_pair(index_(bus_s->_name), index_(bus_d->_name)));
            
            if (m_theta_lb < -3.14 && m_theta_ub > 3.14) {
                cos_max_ = 1;
                cos_min_ = -1;
            } else if (m_theta_lb < 0 && m_theta_ub > 0) {
                cos_max_ = 1;
                cos_min_ = gravity::min(cos(m_theta_lb), cos(m_theta_ub));
            } else {
                cos_max_ = gravity::max(cos(m_theta_lb), cos(m_theta_ub));
                cos_min_ = gravity::min(cos(m_theta_lb), cos(m_theta_ub));
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
            
            wr_max.add_val(name, wr_max_);
            wr_min.add_val(name, wr_min_);
            wi_max.add_val(name, wi_max_);
            wi_min.add_val(name, wi_min_);
        }
    }
}



