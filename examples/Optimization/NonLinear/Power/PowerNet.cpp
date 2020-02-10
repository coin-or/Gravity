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
#include <gravity/model.h>

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
    pf_from_min.set_name("pf_from_min");
    pf_from_max.set_name("pf_from_max");
    qf_from_min.set_name("qf_from_min");
    qf_from_max.set_name("qf_from_max");
    pf_to_min.set_name("pf_to_min");
    pf_to_max.set_name("pf_to_max");
    qf_to_min.set_name("qf_to_min");
    qf_to_max.set_name("qf_to_max");
    lij_min.set_name("lij_min");
    lij_max.set_name("lij_max");
    Iij_min.set_name("Iij_min");
    Iij_max.set_name("Iij_max");
    lji_min.set_name("lji_min");
    lji_max.set_name("lji_min");
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
    cc.set_name("cc");
    dd.set_name("dd");
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

void PowerNet::update_pij_bounds()
{
    double Pg_sum_min,Pg_sum_max, p_min,p_max,pload,qload,shunt_wmax, shunt_wmin;
    bool node_src;
    string arc_name;
    for (auto n: nodes) {
        if(n->get_degree()==1)
        {
            DebugOn("Degree 1 node "<<n->_name<<endl);
            if (n->_active)
            {
                Pg_sum_min=0;
                Pg_sum_max=0;
                pload=pl(n->_name).eval();
                qload=pl(n->_name).eval();
                if(((Bus*)n)->get_out().size()>=1)
                {
                    arc_name=((Bus*)n)->get_out()[0]->_name;
                    node_src=true;
                }
                else
                {
                    arc_name=((Bus*)n)->get_in()[0]->_name;
                    node_src=false;
                }
                if(((Bus*)n)->_has_gen)
                {
                    
                    for(auto gen: ((Bus*)n)->_gen)
                    {
                        Pg_sum_min+=pg_min(gen->_name).eval();
                        Pg_sum_max+=pg_max(gen->_name).eval();
                    }
                    //Check if required at all if(((Bus*)n)->gs()>=0)
                    
                    shunt_wmax=((Bus*)n)->gs()*w_max(n->_name).eval();
                    shunt_wmin=((Bus*)n)->gs()*w_min(n->_name).eval();
                    
                    
                    if(node_src)
                    {
                        p_min=std::max(Pg_sum_min-pload-shunt_wmax, pf_from_min(arc_name).eval());
                        pf_from_min.add_val(arc_name, p_min);
                        p_max=std::min(Pg_sum_max-pload-shunt_wmin, pf_from_max(arc_name).eval());
                        pf_from_max.add_val(arc_name, p_max);
                    }
                    else
                    {
                        p_min=std::max(Pg_sum_min-pload-shunt_wmax, pf_to_min(arc_name).eval());
                        pf_to_min.add_val(arc_name,p_min);
                        p_max=std::min(Pg_sum_max-pload-shunt_wmin, pf_to_max(arc_name).eval());
                        pf_to_max.add_val(arc_name,p_max);
                    }
                }
                else
                {
                    DebugOn("No gen for node " <<n->_name<<endl);
                    if(((Bus*)n)->gs()==0)
                    {
                        
                        if(node_src==true)
                        {
                            pf_from_min.add_val(arc_name,pload*(-1));
                            pf_from_max.add_val(arc_name,pload*(-1));
                            //                            if(pload>=0)
                            //                            {
                            //                                pf_to_min.add_val(arc_name,std::max(pload, pf_to_min(arc_name).eval()));
                            //                            }
                        }
                        else
                        {
                            
                            pf_to_min.add_val(arc_name,pload*(-1));
                            pf_to_max.add_val(arc_name,pload*(-1));
                            //                            if(pload>=0)
                            //                            {
                            //                                pf_from_min.add_val(arc_name,std::max(pload, pf_from_min(arc_name).eval()));
                            //                            }
                            
                        }
                    }
                    if(((Bus*)n)->bs()==0)
                    {
                        qload=ql(n->_name).eval();
                        //                        auto q_map=pf_from_min.get_keys_map();
                        //                        auto key_pos=q_map->at(arc_name);
                        if(node_src==true)
                        {
                            qf_from_min.add_val(arc_name,qload*(-1));
                            qf_from_max.add_val(arc_name,qload*(-1));
                        }
                        else
                        {
                            qf_to_min.add_val(arc_name,qload*(-1));
                            qf_to_max.add_val(arc_name,qload*(-1));
                        }
                        
                    }
                }
            }
        }
    }
    for (auto a: arcs)
    {
        arc_name=a->_name;
        Bus* bus_s = (Bus*)(a->_src);
        Bus* bus_d = (Bus*)(a->_dest);
        auto bus_pair_name = bus_s->_name + "," + bus_d->_name;
        
        DebugOn("Arc name\t" <<arc_name<<endl);
        if(pf_from_max(arc_name).eval()<=0 && tr(arc_name).eval()==1 && as(arc_name).eval()==0)
        {
            pf_to_min.add_val(arc_name,std::max(pf_from_max(arc_name).eval()*(-1), pf_to_min(arc_name).eval()));
        }
        if(pf_to_max(arc_name).eval()<=0 && tr(arc_name).eval()==1 && as(arc_name).eval()==0)
        {
            pf_from_min.add_val(arc_name,std::max(pf_to_max(arc_name).eval()*(-1), pf_from_min(arc_name).eval()));
        }
        if(qf_from_max(arc_name).eval()<=0 && tr(arc_name).eval()==1 && as(arc_name).eval()==0 && ch(arc_name).eval()==0)
        {
            qf_to_min.add_val(arc_name,std::max(qf_from_max(arc_name).eval()*(-1), qf_to_min(arc_name).eval()));
        }
        if(qf_to_max(arc_name).eval()<=0 && tr(arc_name).eval()==1 && as(arc_name).eval()==0 && ch(arc_name).eval()==0)
        {
            qf_from_min.add_val(arc_name,std::max(qf_to_max(arc_name).eval()*(-1), qf_from_min(arc_name).eval()));
        }
        
        if(g(arc_name).eval()==0 && b(arc_name).eval()<=0)
        {
            //Set flag to closed for pf_to
            auto p_map=pf_to_min.get_keys_map();
            auto key_pos=p_map->at(arc_name);
            pf_to_min._off[key_pos]=true;
            if(pf_from_min(arc_name).eval()>=0 )
            {
                th_min.add_val(bus_pair_name, std::max(0.0, th_min(bus_pair_name).eval()));
                wi_min.add_val(bus_pair_name, std::max(0.0, wi_min(bus_pair_name).eval()));
            }
            if(pf_from_max(arc_name).eval()<=0 )
            {
                th_max.add_val(bus_pair_name, std::min(0.0, th_max(bus_pair_name).eval()));
                wi_max.add_val(bus_pair_name, std::min(0.0, wi_max(bus_pair_name).eval()));
                
            }
            if(pf_to_min(arc_name).eval()>=0 )
            {
                th_max.add_val(bus_pair_name, std::min(0.0, th_max(bus_pair_name).eval()));
                wi_max.add_val(bus_pair_name, std::min(0.0, wi_max(bus_pair_name).eval()));
            }
            if(pf_to_max(arc_name).eval()<=0 )
            {
                th_min.add_val(bus_pair_name, std::max(0.0, th_min(bus_pair_name).eval()));
                wi_min.add_val(bus_pair_name, std::max(0.0, wi_min(bus_pair_name).eval()));
            }
        }
    }
}
indices PowerNet::gens_per_node() const{
    indices ids("gens_per_node");
    ids = indices(gens);
    ids._type = matrix_;
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
indices PowerNet::get_ref_node_pairs_from()
{
    
    gravity::indices ref_node_pairs_from("ref_node_pairs_from");
    for (auto a: arcs) {
        if (!a->_parallel && a->_dest->_name==ref_bus)
        {
            ref_node_pairs_from.add(a->_src->_name+","+a->_dest->_name);
        }
    }
    return ref_node_pairs_from;
}

indices PowerNet::get_ref_node_pairs_to()
{
    
    gravity::indices ref_node_pairs_to("ref_node_pairs_to");
    for (auto a: arcs) {
        if (!a->_parallel && a->_src->_name==ref_bus) {
            ref_node_pairs_to.add(a->_src->_name+","+a->_dest->_name);
        }
    }
    return ref_node_pairs_to;
}




indices PowerNet::out_arcs_per_node() const{
    auto ids = indices(arcs);
    ids._type = matrix_;
    ids.set_name("out_arcs_per_node");
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
    ids._type = matrix_;
    ids.set_name("in_arcs_per_node");
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

unsigned PowerNet::get_nb_active_node_pairs() const {
    unsigned nb=0;
    for (auto bp: _node_pairs._keys) {
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
        if(pg_min.eval()<=0)
            pg_min_sq.add_val(name, 0.0);
        else
            pg_min_sq.add_val(name, std::min(std::pow(pg_min.eval(),2), std::pow(pg_max.eval(),2)) );
        
        pg_max_sq.add_val(name, std::max(std::pow(pg_min.eval(),2), std::pow(pg_max.eval(),2)) );
        
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
        if(arc->r==0){
            DebugOn("Branch with zero resistance: " << arc->_name << endl);
        }
        file >> word;
        arc->x = atof(word.c_str());
        if(arc->x==0){
            DebugOn("Branch with zero reactance: " << arc->_name << endl);
        }
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
        name = arc->_name;
        cc.add_val(name, arc->cc);
        dd.add_val(name, arc->dd);
        bt.add_val(name,arc->b/pow(arc->tr,2));
        cht.add_val(name, arc->ch/pow(arc->tr,2));
        ch_half.add_val(name, arc->ch*0.5);
        cht_half.add_val(name, arc->ch*0.5/pow(arc->tr,2));
        rty.add_val(name, arc->g*arc->cc-arc->b*arc->dd);
        ity.add_val(name, arc->b*arc->cc+arc->g*arc->dd);
        
        //        arc->tbound.max = 30*pi/180;
        m_theta_ub += arc->tbound.max;
        
        Bus* bus_s = (Bus*)(arc->_src);
        Bus* bus_d = (Bus*)(arc->_dest);
        
        arc->smax = gravity::max(
                                 pow(bus_s->vbound.max,2)*(arc->g*arc->g + arc->b*arc->b)*(pow(bus_s->vbound.max,2) + pow(bus_d->vbound.max,2)),
                                 pow(bus_d->vbound.max,2)*(arc->g*arc->g+arc->b*arc->b)*(pow(bus_d->vbound.max,2) + pow(bus_s->vbound.max,2))
                                 );
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
        Smax.add_val(name,Cpx(arc->limit,arc->limit));
        pf_from_max.add_val(name, arc->limit);
        pf_from_min.add_val(name, (arc->limit)*(-1));
        qf_from_max.add_val(name, arc->limit);
        qf_from_min.add_val(name, (arc->limit)*(-1));
        pf_to_max.add_val(name, arc->limit);
        pf_to_min.add_val(name, (arc->limit)*(-1));
        qf_to_max.add_val(name, arc->limit);
        qf_to_min.add_val(name, (arc->limit)*(-1));
        lij_min.add_val(name, 0);
        lij_max.add_val(name, pow(arc->tr,2)*pow(arc->limit, 2)/pow(bus_s->vbound.min,2));
        lji_min.add_val(name, 0);
        lji_max.add_val(name, pow(arc->limit, 2)/pow(bus_d->vbound.min,2));
        Iij_min.add_val(name, (-1)*(arc->tr)*(arc->limit)/(bus_s->vbound.min));
        Iij_max.add_val(name, (arc->tr)*(arc->limit)/(bus_s->vbound.min));
        
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
        /* Switching to node_pairs keys */
        name = bus_s->_name + "," + bus_d->_name;
        if (arc->_active && bus_pair_names.count(name)==0) {
            //            _node_pairs._keys.push_back(new index_pair(index_(bus_s->_name), index_(bus_d->_name), arc->_active));
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

/* Create imaginary lines, fill node_pairs_chord, set lower and upper bounds */
void PowerNet::update_net(){
    string name;
    double cos_max_, cos_min_, sin_max_, sin_min_;
    double wr_max_, wr_min_, wi_max_, wi_min_, w_max_, w_min_;
    Node *src, *dest, *n;
    Arc *new_arc;
    int fixed = 1, id_sorted = 0; //id of the current bag in bags_sorted
    Arc *a12, *a13, *a32;
    std::vector<pair<string,std::vector<Node*>>> bags_sorted;
    
    // bags are cliques in the chordal completion graph
    for(auto& b: _bags){
        for(int i = 0; i < b.second.size()-1; i++) {
            for(int j = i+1; j < b.second.size(); j++) {
                Arc* a = get_arc(b.second[i]->_name,b.second[j]->_name);
                if (a==nullptr) {
                    src = get_node(b.second[i]->_name);
                    dest = get_node(b.second[j]->_name);
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
            auto b = *(b_it);
            if(b.second.size() == 3) {
                DebugOff("\nBag: " << b.second[0]->_name << ", " << b.second[1]->_name << ", " << b.second[2]->_name);
                a12 = get_arc(b.second[0], b.second[1]);
                a13 = get_arc(b.second[0], b.second[2]);
                a32 = get_arc(b.second[2], b.second[1]);
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
                
                for(int i = 0; i < b.second.size()-1; i++) {
                    for (int j = i + 1; j < b.second.size(); j++) {
                        Arc* a = get_arc(b.second[i]->_name, b.second[j]->_name);
                        if (!a->_free) continue;
                        n = a->_src;
                        //by now, all arcs in bags should be created
                        for (auto n1: b.second) {
                            if(n==n1) continue;
                            Arc* a2 = get_arc(n->_name, n1->_name);
                            if (a2->_free) continue;
                            Arc *a1 = get_arc(a->_dest, n1);
                            if (!a1->_free) {
                                a->_free = false;
                                
                                pair<string,vector<Node*>> bag;
                                bag.second.push_back(get_node(n->_name));
                                bag.second.push_back(get_node(a->_dest->_name));
                                bag.second.push_back(get_node(n1->_name));
                                bag.first =n->_name+","+a->_dest->_name+","+n1->_name;
                                //                                sort(bag.begin(), bag.end(),
                                //                                     [](const Node *a, const Node *b) -> bool { return a->_id < b->_id; });
                                
                                fixed++;
                                sort(bag.second.begin(), bag.second.end(), [](const Node* a, const Node* b) -> bool{return a->_id < b->_id;});
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
        auto b = *b_it;
        if(b.second.size() >= 2) bags_sorted.push_back(b);
        _bags.erase(b_it);
        //            id_sorted++;
    }
    _bags = bags_sorted;
    
    for(auto& a: arcs) {
        if(a->_imaginary) a->_free = true;
    }
    
    
    for(auto& k: _node_pairs._keys){
        _node_pairs_chord._keys.push_back(new index_pair(*k));
    }
    
    for(auto& a: arcs){
        if(a->_imaginary){
            Bus* bus_s = (Bus*)(a->_src);
            Bus* bus_d = (Bus*)(a->_dest);
            
            name = bus_s->_name + "," + bus_d->_name;
            _node_pairs_chord._keys.push_back(new index_pair(index_(bus_s->_name), index_(bus_d->_name)));
            
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
        for (int i = 0; i < b.second.size(); i++) {
            DebugOff(b.at(i)->_name << " ");
        }
        DebugOff("}" << endl);
        if(add_3d_nlin && b.second.size()==3){
            for(int i = 0; i < 2; i++) {
                for(int j = i+1; j < 3; j++) {
                    Arc* aij = get_arc(b.second[i],b.second[j]);
                    aij->_free = false;
                }
            }
        }
    }
}

void PowerNet::update_ref_bus()
{
    string ref_bus1=ref_bus;
    int deg=0, ref_deg=0;
    ref_deg=get_node(ref_bus)->get_degree();
    if (ref_deg!=get_nb_active_nodes()-1)
    {
        for(auto n:nodes)
        {
            deg=(*n).get_degree();
            if(deg>ref_deg)
            {
                ref_bus1=n->_name;
                ref_deg=deg;
            }
        }
        ref_bus=ref_bus1;
    }
    // ref_bus = "1";
}





shared_ptr<Model<>> PowerNet::build_SCOPF(PowerModelType pmt, int output, double tol){
    auto node_pairs = get_node_pairs();
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
    SOCPF->add(R_Wij.in(node_pairs));
    SOCPF->add(Im_Wij.in(node_pairs));
    
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
    SOC = pow(R_Wij, 2) + pow(Im_Wij, 2) - Wii.from(node_pairs)*Wii.to(node_pairs);
    SOCPF->add(SOC.in(node_pairs) <= 0);
    
    /* Flow conservation */
    Constraint<> KCL_P("KCL_P");
    KCL_P  = sum(Pf_from, out_arcs) + sum(Pf_to, in_arcs) + pl - sum(Pg, gens) + gs*Wii;
    SOCPF->add(KCL_P.in(nodes) == 0);
    
    Constraint<> KCL_Q("KCL_Q");
    KCL_Q  = sum(Qf_from, out_arcs) + sum(Qf_to, in_arcs) + ql - sum(Qg, gens) - bs*Wii;
    SOCPF->add(KCL_Q.in(nodes) == 0);
    
    /* AC Power Flow */
    Constraint<> Flow_P_From("Flow_P_From");
    Flow_P_From = Pf_from - (g_ff*Wii.from(arcs) + g_ft*R_Wij.in_pairs(arcs) + b_ft*Im_Wij.in_pairs(arcs) );
    SOCPF->add(Flow_P_From.in(arcs) == 0);
    
    Constraint<> Flow_P_To("Flow_P_To");
    Flow_P_To = Pf_to - (g_tt*Wii.to(arcs) + g_tf*R_Wij.in_pairs(arcs)  - b_tf*Im_Wij.in_pairs(arcs) );
    SOCPF->add(Flow_P_To.in(arcs) == 0);
    
    Constraint<> Flow_Q_From("Flow_Q_From");
    Flow_Q_From = Qf_from - (g_ft*Im_Wij.in_pairs(arcs)  - b_ff*Wii.from(arcs) - b_ft*R_Wij.in_pairs(arcs) );
    SOCPF->add(Flow_Q_From.in(arcs) == 0);
    
    Constraint<> Flow_Q_To("Flow_Q_To");
    Flow_Q_To = Qf_to + (b_tt*Wii.to(arcs) + b_tf*R_Wij.in_pairs(arcs)  + g_tf*Im_Wij.in_pairs(arcs) );
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
    LNC1 += (v_min.from(node_pairs)+v_max.from(node_pairs))*(v_min.to(node_pairs)+v_max.to(node_pairs))*(sphi*Im_Wij + cphi*R_Wij);
    LNC1 -= v_max.to(node_pairs)*cos_d*(v_min.to(node_pairs)+v_max.to(node_pairs))*Wii.from(node_pairs);
    LNC1 -= v_max.from(node_pairs)*cos_d*(v_min.from(node_pairs)+v_max.from(node_pairs))*Wii.to(node_pairs);
    LNC1 -= v_max.from(node_pairs)*v_max.to(node_pairs)*cos_d*(v_min.from(node_pairs)*v_min.to(node_pairs) - v_max.from(node_pairs)*v_max.to(node_pairs));
    SOCPF->add(LNC1.in(node_pairs) >= 0);
    
    Constraint<> LNC2("LNC2");
    LNC2 += (v_min.from(node_pairs)+v_max.from(node_pairs))*(v_min.to(node_pairs)+v_max.to(node_pairs))*(sphi*Im_Wij + cphi*R_Wij);
    LNC2 -= v_min.to(node_pairs)*cos_d*(v_min.to(node_pairs)+v_max.to(node_pairs))*Wii.from(node_pairs);
    LNC2 -= v_min.from(node_pairs)*cos_d*(v_min.from(node_pairs)+v_max.from(node_pairs))*Wii.to(node_pairs);
    LNC2 += v_min.from(node_pairs)*v_min.to(node_pairs)*cos_d*(v_min.from(node_pairs)*v_min.to(node_pairs) - v_max.from(node_pairs)*v_max.to(node_pairs));
    SOCPF->add(LNC2.in(node_pairs) >= 0);
    return SOCPF;
}



shared_ptr<Model<>> build_ACOPF(PowerNet& grid, PowerModelType pmt, int output, double tol){
    if(grid.nodes.empty()){
        throw invalid_argument("empty grid! Make sure to call grid.readgrid(fname);");
    }
    /** Sets */
    auto node_pairs = grid.get_node_pairs();
    auto nodes = indices(grid.nodes);
    auto arcs = indices(grid.arcs);
    
    auto gens = indices(grid.gens);
    auto gen_nodes = grid.gens_per_node();
    auto out_arcs = grid.out_arcs_per_node();
    auto in_arcs = grid.in_arcs_per_node();
    
    /* Grid Parameters */
    auto pg_min = grid.pg_min.in(gens);
    auto pg_max = grid.pg_max.in(gens);
    auto qg_min = grid.qg_min.in(gens);
    auto qg_max = grid.qg_max.in(gens);
    auto c1 = grid.c1.in(gens);
    auto c2 = grid.c2.in(gens);
    auto c0 = grid.c0.in(gens);
    auto pl = grid.pl.in(nodes);
    auto ql = grid.ql.in(nodes);
    auto gs = grid.gs.in(nodes);
    auto bs = grid.bs.in(nodes);
    auto b = grid.b.in(arcs);
    auto g = grid.g.in(arcs);
    auto as = grid.as.in(arcs);
    auto ch = grid.ch.in(arcs);
    auto tr = grid.tr.in(arcs);
    auto th_min = grid.th_min.in(node_pairs);
    auto th_max = grid.th_max.in(node_pairs);
    auto g_ft = grid.g_ft.in(arcs);
    auto g_ff = grid.g_ff.in(arcs);
    auto g_tt = grid.g_tt.in(arcs);
    auto g_tf = grid.g_tf.in(arcs);
    auto b_ft = grid.b_ft.in(arcs);
    auto b_ff = grid.b_ff.in(arcs);
    auto b_tf = grid.b_tf.in(arcs);
    auto b_tt = grid.b_tt.in(arcs);
    auto S_max = grid.S_max.in(arcs);
    auto v_max = grid.v_max.in(nodes);
    auto v_min = grid.v_min.in(nodes);
    auto tan_th_min = grid.tan_th_min.in(node_pairs);
    auto tan_th_max = grid.tan_th_max.in(node_pairs);
    
    grid.update_ref_bus();
    
    bool polar = (pmt==ACPOL);
    if (polar) {
        DebugOn("Using polar model\n");
    }
    else {
        DebugOn("Using rectangular model\n");
    }
    auto ACOPF = make_shared<Model<>>("AC-OPF Model");
    /** Variables */
    /* Power generation variables */
    var<> Pg("Pg", pg_min, pg_max);
    var<> Qg ("Qg", qg_min, qg_max);
    ACOPF->add(Pg.in(gens),Qg.in(gens));
    //    Pg.copy_vals(grid.pg_s);
    //    Pg.initialize_av();
    //    Qg.initialize_uniform();
    /* Power flow variables */
    var<> Pf_from("Pf_from", -1.*S_max,S_max);
    var<> Qf_from("Qf_from", -1.*S_max,S_max);
    var<> Pf_to("Pf_to", -1.*S_max,S_max);
    var<> Qf_to("Qf_to", -1.*S_max,S_max);
    ACOPF->add(Pf_from.in(arcs), Qf_from.in(arcs),Pf_to.in(arcs),Qf_to.in(arcs));
    
    /** Voltage related variables */
    var<> theta("theta");
    var<> v("|V|", v_min, v_max);
    var<> vr("vr", -1.*v_max,v_max);
    var<> vi("vi", -1.*v_max,v_max);
    
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
        ACOPF->add(vr.in(nodes));
        ACOPF->add(vi.in(nodes));
        vr.initialize_all(1);
        vr_from = vr.from(arcs);
        vr_to = vr.to(arcs);
        vi_from = vi.from(arcs);
        vi_to = vi.to(arcs);
        //        vr.initialize_uniform(0.99,1.01);
    }
    
    /** Construct the objective function */
    /**  Objective */
    auto obj = product(c1,Pg) + product(c2,pow(Pg,2)) + sum(c0);
    ACOPF->min(obj);
    
    /** Define constraints */
    
    /* REF BUS */
    Constraint<> Ref_Bus("Ref_Bus");
    if (polar) {
        Ref_Bus = theta(grid.ref_bus);
    }
    else {
        Ref_Bus = vi(grid.ref_bus);
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
    }
    
    
    /* Phase Angle Bounds constraints */
    Constraint<> PAD_UB("PAD_UB");
    Constraint<> PAD_LB("PAD_LB");
    if (polar) {
        PAD_UB = theta.from(node_pairs) - theta.to(node_pairs);
        PAD_UB -= th_max;
        PAD_LB = theta.from(node_pairs) - theta.to(node_pairs);
        PAD_LB -= th_min;
    }
    else {
        DebugOff("Number of node_pairs = " << node_pairs.size() << endl);
        PAD_UB = vi.from(node_pairs)*vr.to(node_pairs) - vr.from(node_pairs)*vi.to(node_pairs);
        PAD_UB -= tan_th_max*(vr.from(node_pairs)*vr.to(node_pairs) + vi.from(node_pairs)*vi.to(node_pairs));
        
        PAD_LB = vi.from(node_pairs)*vr.to(node_pairs) - vr.from(node_pairs)*vi.to(node_pairs);
        PAD_LB -= tan_th_min*(vr.from(node_pairs)*vr.to(node_pairs) + vi.from(node_pairs)*vi.to(node_pairs));
    }
    ACOPF->add(PAD_UB.in(node_pairs) <= 0);
    ACOPF->add(PAD_LB.in(node_pairs) >= 0);
    
    
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

shared_ptr<Model<>> build_SDPOPF_QC(PowerNet& grid, bool loss, double upper_bound, double lower_bound)
{
    bool relax, sdp_cuts = true,  llnc=true, lazy_bool = false, current=true, add_original=true, convexify=true;
    loss=true;
    size_t num_bags = 0;
    string num_bags_s = "100";
    num_bags = atoi(num_bags_s.c_str());
    
    cout << "\nnum bags = " << num_bags << endl;
    
    
    //grid.update_ref_bus();
    
    /* Grid Stats */
    auto nb_gen = grid.get_nb_active_gens();
    auto nb_lines = grid.get_nb_active_arcs();
    auto nb_buses = grid.get_nb_active_nodes();
    DebugOn("nb active gens = " << nb_gen << endl);
    DebugOn("nb active lines = " << nb_lines << endl);
    DebugOn("nb active buses = " << nb_buses << endl);
    
    /** Sets */
    auto bags_3d=grid.decompose_bags_3d();
    auto node_pairs = grid.get_node_pairs();
    auto node_pairs_chord = grid.get_node_pairs_chord(bags_3d);
    if (grid._tree || !grid.add_3d_nlin || !sdp_cuts) {
        node_pairs_chord = node_pairs;
    }
    auto nodes = indices(grid.nodes);
    auto arcs = indices(grid.arcs);
    auto gens = indices(grid.gens);
    auto gen_nodes = grid.gens_per_node();
    auto out_arcs = grid.out_arcs_per_node();
    auto in_arcs = grid.in_arcs_per_node();
    //  grid.update_pij_bounds();
    
    /* Grid Parameters */
    auto pg_min = grid.pg_min.in(gens);
    auto pg_max = grid.pg_max.in(gens);
    auto qg_min = grid.qg_min.in(gens);
    auto qg_max = grid.qg_max.in(gens);
    auto pf_from_min=grid.pf_from_min.in(arcs);
    auto pf_from_max=grid.pf_from_max.in(arcs);
    auto qf_from_min=grid.qf_from_min.in(arcs);
    auto qf_from_max=grid.qf_from_max.in(arcs);
    auto pf_to_min=grid.pf_to_min.in(arcs);
    auto pf_to_max=grid.pf_to_max.in(arcs);
    auto qf_to_min=grid.qf_to_min.in(arcs);
    auto qf_to_max=grid.qf_to_max.in(arcs);
    auto lij_min=grid.lij_min.in(arcs);
    auto lij_max=grid.lij_max.in(arcs);
    auto lji_min=grid.lji_min.in(arcs);
    auto lji_max=grid.lji_max.in(arcs);
    auto Iij_min=grid.Iij_min.in(arcs);
    auto Iij_max=grid.Iij_max.in(arcs);
    auto c1 = grid.c1.in(gens);
    auto c2 = grid.c2.in(gens);
    auto c0 = grid.c0.in(gens);
    auto pl = grid.pl.in(nodes);
    auto ql = grid.ql.in(nodes);
    auto gs = grid.gs.in(nodes);
    auto bs = grid.bs.in(nodes);
    auto b = grid.b.in(arcs);
    auto g = grid.g.in(arcs);
    auto as = grid.as.in(arcs);
    auto ch = grid.ch.in(arcs);
    auto tr = grid.tr.in(arcs);
    auto th_min = grid.th_min.in(node_pairs);
    auto th_max = grid.th_max.in(node_pairs);
    auto g_ft = grid.g_ft.in(arcs);
    auto g_ff = grid.g_ff.in(arcs);
    auto g_tt = grid.g_tt.in(arcs);
    auto g_tf = grid.g_tf.in(arcs);
    auto b_ft = grid.b_ft.in(arcs);
    auto b_ff = grid.b_ff.in(arcs);
    auto b_tf = grid.b_tf.in(arcs);
    auto b_tt = grid.b_tt.in(arcs);
    auto S_max = grid.S_max.in(arcs);
    auto v_max = grid.v_max.in(nodes);
    auto v_min = grid.v_min.in(nodes);
    auto w_max = grid.w_max.in(nodes);
    auto w_min = grid.w_min.in(nodes);
    auto tan_th_min = grid.tan_th_min.in(node_pairs);
    auto tan_th_max = grid.tan_th_max.in(node_pairs);
    auto wr_min = grid.wr_min.in(node_pairs_chord);
    auto wr_max = grid.wr_max.in(node_pairs_chord);
    auto wi_min = grid.wi_min.in(node_pairs_chord);
    auto wi_max = grid.wi_max.in(node_pairs_chord);
    auto cc=grid.cc.in(arcs);
    auto dd=grid.dd.in(arcs);
    auto ch_half=grid.ch_half.in(arcs);
    
    //    auto arcs_charged=grid.arcs_line_charge();
    //    auto arcs_inductive=grid.arcs_inductive_only();
    /** Build model */
    
    auto SDPOPF = make_shared<Model<>>("SDP-OPF Model");
    
    var<> Pg("Pg", pg_min, pg_max);
    var<> Qg ("Qg", qg_min, qg_max);
    SDPOPF->add(Pg.in(gens),Qg.in(gens));
    
    /* Power flow variables */
    var<> Pf_from("Pf_from", -1.*S_max,S_max);
    var<> Qf_from("Qf_from", -1.*S_max,S_max);
    var<> Pf_to("Pf_to", -1.*S_max,S_max);
    var<> Qf_to("Qf_to", -1.*S_max,S_max);
    var<> lij("lij", lij_min,lij_max);
    var<> lji("lji", lji_min,lji_max);
    
    
    //    var<> RIij("RIij", Iij_min,Iij_max);
    //    var<> IIij("IIij", Iij_min,Iij_max);
    
    if(current){
        SDPOPF->add(lij.in(arcs));
        SDPOPF->add(lji.in(arcs));
        //        SDPOPF->add(RIij.in(arcs));
        //        SDPOPF->add(IIij.in(arcs));
    }
    
    SDPOPF->add(Pf_from.in(arcs), Qf_from.in(arcs),Pf_to.in(arcs),Qf_to.in(arcs));
    
    
    /* Real part of Wij = ViVj */
    var<>  R_Wij("R_Wij", wr_min, wr_max);
    /* Imaginary part of Wij = ViVj */
    var<>  Im_Wij("Im_Wij", wi_min, wi_max);
    /* Magnitude of Wii = Vi^2 */
    var<>  Wii("Wii", w_min, w_max);
    SDPOPF->add(Wii.in(nodes));
    SDPOPF->add(R_Wij.in(node_pairs_chord));
    SDPOPF->add(Im_Wij.in(node_pairs_chord));
    
    
    var<>  R_Vi("R_Vi", -1*v_max, v_max);
    var<>  Im_Vi("Im_Vi", -1*v_max, v_max);
    
    if(add_original){
        SDPOPF->add(R_Vi.in(nodes),Im_Vi.in(nodes));
        R_Vi.initialize_all(1);
    }
    
    param<> t_min("t_min"), t_max("t_max"), unit("unit"), zero("zero"), oL("oL"), oU("oU");
    t_min.in(nodes);t_max.in(nodes);unit.in(node_pairs);zero.in(node_pairs);
    t_min.set_val((-1)*pi);
    t_max.set_val(pi);
    zero.set_val(0.0);
    oL=0.0;
    oU=3000.0;
    
    unit.set_val(1.0);
    param<> ViVj_L("ViVj_L"), ViVj_U("ViVj_U");
    ViVj_L.in(node_pairs);ViVj_U.in(node_pairs);
    for(auto &a:*node_pairs._keys)
    {
        auto b1 = a.substr(0, a.find_first_of(","));
        auto b2=  a.substr(a.find_first_of(",")+1);
        ViVj_L.set_val(a, v_min.eval(b1)*v_min.eval(b2));
        ViVj_U.set_val(a, v_max.eval(b1)*v_max.eval(b2));
        
    }
    
    //    var<> objt("objt", oL, oU);
    //    SDPOPF->add(objt);
    
    var<> V_mag("V_mag", v_min, v_max);
    var<> theta("theta", t_min, t_max);
    theta._lift=true;
    var<> theta_ij("theta_ij", th_min, th_max);
    SDPOPF->add(V_mag.in(nodes));
    SDPOPF->add(theta.in(nodes));
    SDPOPF->add(theta_ij.in(node_pairs));
    var<> Vi_Vj("Vi_Vj", ViVj_L, ViVj_U);
    SDPOPF->add(Vi_Vj.in(node_pairs));
    
    var<> costhetaij("costhetaij", min(cos(th_min), cos(th_max)), unit);
    var<> sinthetaij("sinthetaij", sin(th_min), sin(th_max));
    
    SDPOPF->add(sinthetaij.in(node_pairs));
    SDPOPF->add(costhetaij.in(node_pairs));
    
    
    
    /* Initialize variables */
    R_Wij.initialize_all(1.0);
    Wii.initialize_all(1.001);
    
    
    //  Objective
    auto obj = (product(c1,Pg) + product(c2,pow(Pg,2)) + sum(c0));
    SDPOPF->min(obj);
    
    Constraint<> obj_UB("obj_UB");
    obj_UB=(product(c1,Pg) + product(c2,pow(Pg,2)) + sum(c0))-upper_bound;
    SDPOPF->add(obj_UB<=0);
    
    /**  Objective */
    //    auto obj = (product(c1,Pg) + product(c2,pow(Pg,2)) + sum(c0));
    //    SDPOPF->min(obj);
    //    Constraint<> obj_UB("obj_UB");
    //    obj_UB=objt-upper_bound;
    //    SDPOPF->add(obj_UB<=0);
    
    //
    //    Constraint<> obj_LB("obj_LB");
    //    obj_LB=objt-lower_bound;
    //    SDPOPF->add(obj_LB>=0);
    //
    //
    //
    //    var<> Pg2("Pg2", pow(pg_min,2), pow(pg_max,2));
    //    SDPOPF->add(Pg2.in(gens));
    //
    //    Constraint<> PGSquare("PGSquare");
    //    PGSquare = pow(Pg,2) - Pg2;
    //    SDPOPF->add(PGSquare.in(gens)==0,true);
    //
    //
    //    Constraint<> obj_def("obj_def");
    //    obj_def=(product(c1,Pg) + product(c2,Pg2) + sum(c0))-objt;
    //    SDPOPF->add(obj_def==0);
    
    //  SDPOPF->min(objt);
    // SDPOPF->print();
    
    /** Constraints */
    if(!grid._tree && grid.add_3d_nlin && sdp_cuts) {
        
        auto bag_size = bags_3d.size();
        DebugOn("\nNum of bags = " << bag_size << endl);
        DebugOn("Adding 3d determinant polynomial cuts\n");
        auto R_Wij_ = R_Wij.pairs_in_bags(bags_3d, 3);
        auto Im_Wij_ = Im_Wij.pairs_in_bags(bags_3d, 3);
        auto Wii_ = Wii.in_bags(bags_3d, 3);
        
        
        Constraint<> SDP3("SDP_3D");
        SDP3 = 2 * R_Wij_[0] * (R_Wij_[1] * R_Wij_[2] + Im_Wij_[1] * Im_Wij_[2]);
        SDP3 -= 2 * Im_Wij_[0] * (R_Wij_[2] * Im_Wij_[1] - Im_Wij_[2] * R_Wij_[1]);
        SDP3 -= (pow(R_Wij_[0], 2) + pow(Im_Wij_[0], 2)) * Wii_[2];
        SDP3 -= (pow(R_Wij_[1], 2) + pow(Im_Wij_[1], 2)) * Wii_[0];
        SDP3 -= (pow(R_Wij_[2], 2) + pow(Im_Wij_[2], 2)) * Wii_[1];
        SDP3 += Wii_[0] * Wii_[1] * Wii_[2];
        if (lazy_bool) {
            SDPOPF->add_lazy(SDP3 >= 0);
        }
        else {
            SDPOPF->add(SDP3 >= 0);
            DebugOn("Number of 3d determinant cuts = " << SDP3.get_nb_instances() << endl);
        }
        
    }
    
    /** Constraints */
    /* Second-order cone constraints */
    
    Constraint<> SOC_EQ("SOC_EQ");
    SOC_EQ = pow(R_Wij, 2) + pow(Im_Wij, 2) - Wii.from(node_pairs_chord)*Wii.to(node_pairs_chord);
    SDPOPF->add(SOC_EQ.in(node_pairs_chord) == 0, convexify);
    
    /* Flow conservation */
    Constraint<> KCL_P("KCL_P");
    KCL_P  = sum(Pf_from, out_arcs) + sum(Pf_to, in_arcs) + pl - sum(Pg, gen_nodes) + gs*Wii;
    SDPOPF->add(KCL_P.in(nodes) == 0);
    
    Constraint<> KCL_Q("KCL_Q");
    KCL_Q  = sum(Qf_from, out_arcs) + sum(Qf_to, in_arcs) + ql - sum(Qg, gen_nodes) - bs*Wii;
    SDPOPF->add(KCL_Q.in(nodes) == 0);
    
    /* AC Power Flow */
    Constraint<> Flow_P_From("Flow_P_From");
    Flow_P_From = Pf_from - (g_ff*Wii.from(arcs) + g_ft*R_Wij.in_pairs(arcs)  + b_ft*Im_Wij.in_pairs(arcs) );
    SDPOPF->add(Flow_P_From.in(arcs) == 0);
    
    Constraint<> Flow_P_To("Flow_P_To");
    Flow_P_To = Pf_to - (g_tt*Wii.to(arcs) + g_tf*R_Wij.in_pairs(arcs)  - b_tf*Im_Wij.in_pairs(arcs) );
    SDPOPF->add(Flow_P_To.in(arcs) == 0);
    
    Constraint<> Flow_Q_From("Flow_Q_From");
    Flow_Q_From = Qf_from - (g_ft*Im_Wij.in_pairs(arcs)  - b_ff*Wii.from(arcs) - b_ft*R_Wij.in_pairs(arcs) );
    SDPOPF->add(Flow_Q_From.in(arcs) == 0);
    
    Constraint<> Flow_Q_To("Flow_Q_To");
    Flow_Q_To = Qf_to + b_tt*Wii.to(arcs) + b_tf*R_Wij.in_pairs(arcs)  + g_tf*Im_Wij.in_pairs(arcs) ;
    SDPOPF->add(Flow_Q_To.in(arcs) == 0);
    
    /* Phase Angle Bounds constraints */
    Constraint<> PAD_UB("PAD_UB");
    PAD_UB = Im_Wij.in(node_pairs);
    PAD_UB <= tan_th_max*R_Wij.in(node_pairs);
    SDPOPF->add(PAD_UB.in(node_pairs));
    
    Constraint<> PAD_LB("PAD_LB");
    PAD_LB =  Im_Wij.in(node_pairs);
    PAD_LB >= tan_th_min*R_Wij.in(node_pairs);
    SDPOPF->add(PAD_LB.in(node_pairs));
    
    /* Thermal Limit Constraints */
    Constraint<> Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from = pow(Pf_from, 2) + pow(Qf_from, 2);
    Thermal_Limit_from <= pow(S_max,2);
    SDPOPF->add(Thermal_Limit_from.in(arcs));
    
    
    Constraint<> Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to = pow(Pf_to, 2) + pow(Qf_to, 2);
    Thermal_Limit_to <= pow(S_max,2);
    SDPOPF->add(Thermal_Limit_to.in(arcs));
    
    if(add_original){
        Im_Vi.set_lb((grid.ref_bus),0);
        Im_Vi.set_ub((grid.ref_bus),0);
        
        //        theta.set_lb((grid.ref_bus),0);
        //        theta.set_ub((grid.ref_bus),0);
        
        R_Vi.set_lb((grid.ref_bus),v_min(grid.ref_bus).eval());
        R_Vi.set_ub((grid.ref_bus),v_max(grid.ref_bus).eval());
        
        Constraint<> R_V_mag("R_V_mag");
        R_V_mag=R_Vi(grid.ref_bus)-V_mag(grid.ref_bus);
        SDPOPF->add(R_V_mag==0);
        
        
        var<Cpx> Vi("Vi"), Vj("Vj"), Wij("Wij");
        Vi.real_imag(R_Vi.from(node_pairs_chord), Im_Vi.from(node_pairs_chord));
        Vj.real_imag(R_Vi.to(node_pairs_chord), Im_Vi.to(node_pairs_chord));
        Wij.real_imag(R_Wij.in(node_pairs_chord), Im_Wij.in(node_pairs_chord));
        
        Constraint<Cpx> Linking_Wij("Linking_Wij");
        Linking_Wij = Wij - Vi*conj(Vj);
        SDPOPF->add(Linking_Wij.in(node_pairs_chord)==0, convexify);
        
        Vi.real_imag(R_Vi.in(nodes), Im_Vi.in(nodes));
        
        Constraint<Cpx> Linking_Wi("Linking_Wi");
        Linking_Wi = Wii - Vi*conj(Vi);
        SDPOPF->add(Linking_Wi.in(nodes)==0, convexify);
        
        Constraint<> Linking_Wi_V_mag("Linking_Wi_V_mag");
        Linking_Wi_V_mag = Wii - V_mag*V_mag;
        SDPOPF->add(Linking_Wi_V_mag.in(nodes)==0, convexify);
        
        Constraint<Cpx> Linking_V_mag_V("Linking_V_mag_V");
        Linking_V_mag_V =  V_mag*V_mag - pow(R_Vi,2) - pow(Im_Vi,2) ;
        SDPOPF->add(Linking_V_mag_V.in(nodes)>=0);
        
        
        Constraint<> Linking_V_mag_i_V_mag_j("Linking_V_mag_i_V_mag_j");
        Linking_V_mag_i_V_mag_j =  Vi_Vj.in(node_pairs)- V_mag.from(node_pairs)*V_mag.to(node_pairs);
        SDPOPF->add(Linking_V_mag_i_V_mag_j.in(node_pairs)==0, convexify);
        
        
        //
        
        //                func<> thetaij_m;
        //                thetaij_m=max(min(theta.get_ub().from(node_pairs)-theta.get_lb().to(node_pairs), th_max.in(node_pairs)),(max(theta.get_lb().from(node_pairs)-theta.get_ub().to(node_pairs), th_min.in(node_pairs)))*(-1));
        //                thetaij_m.eval_all();
        //                DebugOn("thetaij_m");
        //                thetaij_m.print();
        //
        //                Constraint<> sinenvup("sinenvup");
        //                sinenvup=sinthetaij-cos(thetaij_m*0.5)*(theta.from(node_pairs)-theta.to(node_pairs)-thetaij_m*0.5)-sin(thetaij_m*0.5);
        //                SDPOPF->add(sinenvup.in(node_pairs)<=0);
        //
        //                Constraint<> sinenvlow("sinenvlow");
        //                sinenvlow=sinthetaij-cos(thetaij_m*0.5)*(theta.from(node_pairs)-theta.to(node_pairs)+thetaij_m*0.5)+sin(thetaij_m*0.5);
        //                SDPOPF->add(sinenvlow.in(node_pairs)>=0);
        //
        //                Constraint<> cosenvup("cosenvup");
        //                cosenvup=costhetaij*pow(thetaij_m,2)-1.0*pow(thetaij_m,2)+(1-cos(thetaij_m))*pow((theta.from(node_pairs)-theta.to(node_pairs)), 2);
        //                SDPOPF->add(cosenvup.in(node_pairs)<=0);
        //
        //                Constraint<> cosenvlow("cosenvlow");
        //                cosenvlow=costhetaij-cos(thetaij_m);
        //                SDPOPF->add(cosenvlow.in(node_pairs)>=0);
        
        
        func<> thetaij_m;
        thetaij_m=max(min(theta_ij.get_ub().in(node_pairs), th_max.in(node_pairs)),(max(theta_ij.get_lb().in(node_pairs), th_min.in(node_pairs)))*(-1));
        thetaij_m.eval_all();
        DebugOn("thetaij_m");
        thetaij_m.print();
        
        Constraint<> sinenvup("sinenvup");
        sinenvup=sinthetaij.in(node_pairs)-cos(thetaij_m*0.5)*(theta_ij.in(node_pairs)-thetaij_m*0.5)-sin(thetaij_m*0.5);
        SDPOPF->add(sinenvup.in(node_pairs)<=0);
        
        Constraint<> sinenvlow("sinenvlow");
        sinenvlow=sinthetaij.in(node_pairs)-cos(thetaij_m*0.5)*(theta_ij.in(node_pairs)+thetaij_m*0.5)+sin(thetaij_m*0.5);
        SDPOPF->add(sinenvlow.in(node_pairs)>=0);
        
        Constraint<> cosenvup("cosenvup");
        cosenvup=costhetaij.in(node_pairs)*pow(thetaij_m,2)-1.0*pow(thetaij_m,2)+(1-cos(thetaij_m))*pow(theta_ij.in(node_pairs), 2);
        SDPOPF->add(cosenvup.in(node_pairs)<=0);
        
        Constraint<> cosenvlow("cosenvlow");
        cosenvlow=costhetaij.in(node_pairs)-cos(thetaij_m).in(node_pairs);
        SDPOPF->add(cosenvlow.in(node_pairs)>=0);
        
        
        //        Constraint<> trig("trig");
        //        trig=pow(sinthetaij, 2)+pow(costhetaij,2)-1;
        //
        //        SDPOPF->add(trig.in(node_pairs_chord)==0, true);
        
        
        
        Constraint<> tanU("tanU");
        tanU=Im_Wij-R_Wij*tan(theta_ij.get_ub().in(node_pairs));
        SDPOPF->add(tanU.in(node_pairs)<=0);
        
        Constraint<> tanL("tanL");
        tanL=Im_Wij-R_Wij*tan(theta_ij.get_lb().in(node_pairs));
        SDPOPF->add(tanL.in(node_pairs)>=0);
        
        
        Constraint<> theta_diff("theta_diff");
        theta_diff=theta.from(node_pairs)-theta.to(node_pairs)-theta_ij.in(node_pairs);
        SDPOPF->add(theta_diff.in(node_pairs)==0);
        
        
        Constraint<> Linking_RW_Vtheta("Linking_RW_Vtheta");
        Linking_RW_Vtheta = R_Wij.in(node_pairs) - Vi_Vj.in(node_pairs)*costhetaij.in(node_pairs);
        SDPOPF->add(Linking_RW_Vtheta.in(node_pairs)==0, convexify);
        
        
        Constraint<Cpx> Linking_ImW_Vtheta("Linking_ImW_Vtheta");
        Linking_ImW_Vtheta = Im_Wij.in(node_pairs) - Vi_Vj.in(node_pairs)*sinthetaij.in(node_pairs);
        SDPOPF->add(Linking_ImW_Vtheta.in(node_pairs)==0, convexify);
        
        auto ref_node_pairs_from=grid.get_ref_node_pairs_from();
        auto ref_node_pairs_to=grid.get_ref_node_pairs_to();
        
        
        Constraint<> RV_V("RV_V");
        RV_V =  R_Vi.from(ref_node_pairs_from)-V_mag.from(ref_node_pairs_from)*costhetaij.in(ref_node_pairs_from);
        SDPOPF->add(RV_V.in(ref_node_pairs_from)==0, convexify);
        
        Constraint<> IV_V("IV_V");
        IV_V =  Im_Vi.from(ref_node_pairs_from)-V_mag.from(ref_node_pairs_from)*sinthetaij.in(ref_node_pairs_from);
        SDPOPF->add(IV_V.in(ref_node_pairs_from)==0, convexify);
        
        Constraint<> RV_V1("RV_V1");
        RV_V1 =  R_Vi.to(ref_node_pairs_to)-V_mag.to(ref_node_pairs_to)*costhetaij.in(ref_node_pairs_to);
        SDPOPF->add(RV_V1.in(ref_node_pairs_to)==0, convexify);
        
        Constraint<> IV_V1("IV_V1");
        IV_V1 =  Im_Vi.to(ref_node_pairs_to)+V_mag.to(ref_node_pairs_to)*sinthetaij.in(ref_node_pairs_to);
        SDPOPF->add(IV_V1.in(ref_node_pairs_to)==0, convexify);
        
        
        
        
        
        //                for (auto &p: *ref_node_pairs_from._keys)
        //                {
        //                    auto ngb = p.substr(0,p.find_first_of(","));
        //                    theta.set_lb(ngb, th_min.eval(p));
        //                    theta.set_ub(ngb, th_max.eval(p));
        //
        //                }
        //                for (auto &p: *ref_node_pairs_to._keys)
        //                {
        //                    auto ngb = p.substr(p.find_first_of(",")+1);
        //                    theta.set_lb(ngb, th_min.eval(p));
        //                    theta.set_ub(ngb, th_max.eval(p));;
        //                }
        
        
        
        if(!grid._tree)
        {
            //
            auto Wij_ = Wij.in(node_pairs_chord).pairs_in_bags(grid._bags, 3);
            auto Wii_ = Wii.in_bags(grid._bags, 3);
            auto nb_bags3 = Wij_[0]._indices->size();
            
            //  auto thetaij_=theta_ij.in(node_pairs_chord).pairs_in_bags(grid._bags, 3);
            //auto nb_bagst=thetaij_[0]._indices->size();
            //
            Constraint<Cpx> Rank_type2a("RankType2a");
            Rank_type2a=Wij_[0]*Wij_[1]-Wii_[1]*Wij_[2];
            SDPOPF->add(Rank_type2a.in(range(1,nb_bags3))==0, true);
            
            Constraint<Cpx> Rank_type2b("RankType2b");
            Rank_type2b=Wij_[2]*conj(Wij_[1])-Wii_[2]*Wij_[0];
            SDPOPF->add(Rank_type2b.in(range(1,nb_bags3))==0, true);
            
            Constraint<Cpx> Rank_type2c("RankType2c");
            Rank_type2c=Wij_[2]*conj(Wij_[0])-Wii_[0]*Wij_[1];
            SDPOPF->add(Rank_type2c.in(range(1,nb_bags3))==0, true);
            
            //            Constraint<> Cycle_theta("Cycle_theta");
            //            Cycle_theta=thetaij_[0]+thetaij_[1]-thetaij_[2];
            //            SDPOPF->add(Cycle_theta.in(range(1,nb_bagst))==0);
            
            
            
            
            //            auto ref_node_pairs_ijkl=grid.get_pairsof_node_pairs_ijkl();
            //            DebugOn("firstfirst");
            //            ref_node_pairs_ijkl.first.first.print();
            //            ref_node_pairs_ijkl.first.second.print();
            //            ref_node_pairs_ijkl.second.first.print();
            //            ref_node_pairs_ijkl.second.second.print();
            //            DebugOn("size "<<ref_node_pairs_ijkl.first.first.size());
            //
            //            Constraint<Cpx> Rank_type3("RankType3");
            //            Rank_type3= Wij.in(ref_node_pairs_ijkl.first.first)*Wij.in(ref_node_pairs_ijkl.first.second)-conj(Wij).in(ref_node_pairs_ijkl.second.first)*Wij.in(ref_node_pairs_ijkl.second.second);
            //            SDPOPF->add(Rank_type3.in(indices(1,ref_node_pairs_ijkl.first.first.size()))==0, convexify);
        }
    }
    
    if(current){
        param<Cpx> T("T"), Y("Y"), Ych("Ych");
        var<Cpx> L_from("L_from"), W("W"), Vi("Vi"), Vj("Vj"), I("I");
        T.real_imag(cc.in(arcs), dd.in(arcs));
        Y.real_imag(g.in(arcs), b.in(arcs));
        Ych.set_imag(ch_half.in(arcs));
        
        
        L_from.set_real(lij.in(arcs));
        
        
        W.real_imag(R_Wij.in_pairs(arcs) , Im_Wij.in_pairs(arcs) );
        //        Vi.real_imag(R_Vi.from(arcs), Im_Vi.from(arcs));
        //        Vj.real_imag(R_Vi.to(arcs), Im_Vi.to(arcs));
        //        I.real_imag(RIij.in(arcs), IIij.in(arcs));
        
        //        Constraint<Cpx> I_ij("I_ij");
        //        I_ij=(Y+Ych)*Vi*conj(T)-pow(tr,2)*Y*Vj-pow(tr,2)*I;
        //        SDPOPF->add(I_ij.in(arcs)==0);
        //
        //        Constraint<Cpx> II_conj("II_conj");
        //        II_conj=I*conj(I)-L_from;
        //       SDPOPF->add(II_conj.in(arcs)==0, convexify);
        
        Constraint<Cpx> I_from("I_from");
        I_from=(Y+Ych)*(conj(Y)+conj(Ych))*Wii.from(arcs)-T*Y*(conj(Y)+conj(Ych))*conj(W)-conj(T)*conj(Y)*(Y+Ych)*W+pow(tr,2)*Y*conj(Y)*Wii.to(arcs);
        SDPOPF->add_real(I_from.in(arcs)==pow(tr,2)*L_from);
        
        var<Cpx> L_to("L_to");
        L_to.set_real(lji.in(arcs));
        
        Constraint<Cpx> I_to("I_to");
        I_to=pow(tr,2)*(Y+Ych)*(conj(Y)+conj(Ych))*Wii.to(arcs)-conj(T)*Y*(conj(Y)+conj(Ych))*W-T*conj(Y)*(Y+Ych)*conj(W)+Y*conj(Y)*Wii.from(arcs);
        SDPOPF->add_real(I_to.in(arcs)==pow(tr,2)*L_to);
        
        Constraint<> I_from_Pf("I_from_Pf");
        I_from_Pf=lij*Wii.from(arcs)-pow(tr,2)*(pow(Pf_from,2) + pow(Qf_from,2));
        SDPOPF->add(I_from_Pf.in(arcs)==0, true);
        
        Constraint<> I_to_Pf("I_to_Pf");
        I_to_Pf=lji*Wii.to(arcs)-(pow(Pf_to,2) + pow(Qf_to, 2));
        SDPOPF->add(I_to_Pf.in(arcs)==0, true);
        
        
    }
    
    func<> theta_L = atan(min(Im_Wij.get_lb().in(node_pairs)/R_Wij.get_ub().in(node_pairs),Im_Wij.get_lb().in(node_pairs)/R_Wij.get_lb().in(node_pairs)));
    func<> theta_U = atan(max(Im_Wij.get_ub().in(node_pairs)/R_Wij.get_lb().in(node_pairs),Im_Wij.get_ub().in(node_pairs)/R_Wij.get_ub().in(node_pairs)));
    func<> phi=(theta_U.in(node_pairs)+theta_L.in(node_pairs))/2.0;
    func<> del=(theta_U.in(node_pairs)-theta_L.in(node_pairs))/2.0;
    
    
    Constraint<> LNC1("LNC1");
    LNC1 += (sqrt(Wii.get_lb().from(node_pairs))+sqrt(Wii.get_ub().from(node_pairs)))*(sqrt(Wii.get_lb().to(node_pairs))+sqrt(Wii.get_ub().to(node_pairs)))*(Im_Wij.in(node_pairs)*sin(phi.in(node_pairs)) + R_Wij.in(node_pairs)*cos(phi.in(node_pairs)));
    
    LNC1 -=sqrt(Wii.get_ub().to(node_pairs))*cos(del.in(node_pairs))*(sqrt(Wii.get_lb().to(node_pairs))+sqrt(Wii.get_ub().to(node_pairs)))*Wii.from(node_pairs);
    
    LNC1 -=sqrt(Wii.get_ub().from(node_pairs))*cos(del.in(node_pairs))*(sqrt(Wii.get_lb().from(node_pairs))+sqrt(Wii.get_ub().from(node_pairs)))*Wii.to(node_pairs);
    
    LNC1-=sqrt(Wii.get_ub().from(node_pairs))*sqrt(Wii.get_ub().to(node_pairs))*cos(del)*(sqrt(Wii.get_lb().from(node_pairs))*
                                                                                          sqrt(Wii.get_lb().to(node_pairs)) - sqrt(Wii.get_ub().from(node_pairs))*sqrt(Wii.get_ub().to(node_pairs)));
    SDPOPF->add(LNC1.in(node_pairs) >= 0);
    
    Constraint<> LNC2("LNC2");
    LNC2 += (sqrt(Wii.get_lb().from(node_pairs))+sqrt(Wii.get_ub().from(node_pairs)))*(sqrt(Wii.get_lb().to(node_pairs))+sqrt(Wii.get_ub().to(node_pairs)))*(sin(phi.in(node_pairs))*Im_Wij.in(node_pairs) + cos(phi.in(node_pairs))*R_Wij.in(node_pairs));
    LNC2 -=sqrt(Wii.get_lb().to(node_pairs))*cos(del.in(node_pairs))*(sqrt(Wii.get_lb().to(node_pairs))+sqrt(Wii.get_ub().to(node_pairs)))*Wii.from(node_pairs);
    LNC2 -=sqrt(Wii.get_lb().from(node_pairs))*cos(del.in(node_pairs))*(sqrt(Wii.get_lb().from(node_pairs))+sqrt(Wii.get_ub().from(node_pairs)))*Wii.to(node_pairs);
    LNC2 -=sqrt(Wii.get_lb().from(node_pairs))*sqrt(Wii.get_lb().to(node_pairs))*cos(del.in(node_pairs))*(sqrt(Wii.get_ub().from(node_pairs))*
                                                                                                          sqrt(Wii.get_ub().to(node_pairs))-sqrt(Wii.get_lb().from(node_pairs))*sqrt(Wii.get_lb().to(node_pairs)));
    SDPOPF->add(LNC2.in(node_pairs) >= 0);
    
    if(loss)
    {
        Constraint<> PL("PL");
        PL=Pf_from+Pf_to-g*pow((V_mag.from(arcs)-V_mag.to(arcs)),2);
        SDPOPF->add(PL.in(arcs)>=0);
        
        
        //        Constraint<> PL("PL");
        //        PL=Pf_from+Pf_to-g*(pow(min(cos(theta_ij.in_pairs(arcs).get_lb()),cos(theta_ij.in_pairs(arcs).get_ub())),2)*pow((V_mag.from(arcs)-V_mag.to(arcs)),2)+max(zero.in_pairs(arcs), sin(theta_ij.in_pairs(arcs).get_lb()))*(pow(V_mag.from(arcs),2)+pow(V_mag.to(arcs), 2)));
        //        SDPOPF->add(PL.in(arcs)>=0);
        
        //        Constraint<> PU("PU");
        //        PU=Pf_from+Pf_to-g*pow((V_mag.from(arcs)-V_mag.to(arcs)),2);
        //        SDPOPF->add(PL.in(arcs)>=0);
    }
    
    return SDPOPF;
    
}

shared_ptr<Model<>> build_SDPOPF(PowerNet& grid, bool current, bool nonlin_obj, bool sdp_kim, double upper_bound)
{
    bool relax, sdp_cuts = true,  llnc=false, lazy_bool = false, add_original=false, convexify=true;
    size_t num_bags = 0;
    string num_bags_s = "100";
    num_bags = atoi(num_bags_s.c_str());
    
    cout << "\nnum bags = " << num_bags << endl;
    
    if(!grid._tree && grid._bags.empty()){
        grid.get_tree_decomp_bags();
    }
    //grid.update_ref_bus();
    
    /* Grid Stats */
    auto nb_gen = grid.get_nb_active_gens();
    auto nb_lines = grid.get_nb_active_arcs();
    auto nb_buses = grid.get_nb_active_nodes();
    DebugOn("nb active gens = " << nb_gen << endl);
    DebugOn("nb active lines = " << nb_lines << endl);
    DebugOn("nb active buses = " << nb_buses << endl);
    
    /** Sets */
    auto bags_3d=grid.decompose_bags_3d();
    auto node_pairs = grid.get_node_pairs();
    auto node_pairs_chord = grid.get_node_pairs_chord(bags_3d);
    if (grid._tree || !grid.add_3d_nlin || !sdp_cuts) {
        node_pairs_chord = node_pairs;
    }
    auto nodes = indices(grid.nodes);
    auto arcs = indices(grid.arcs);
    auto arcs_I_to=grid.arcs_not_inductive_only();
    auto gens = indices(grid.gens);
    auto gen_nodes = grid.gens_per_node();
    auto out_arcs = grid.out_arcs_per_node();
    auto in_arcs = grid.in_arcs_per_node();
    //grid.update_pij_bounds();
    
    /* Grid Parameters */
    auto pg_min = grid.pg_min.in(gens);
    auto pg_max = grid.pg_max.in(gens);
    auto qg_min = grid.qg_min.in(gens);
    auto qg_max = grid.qg_max.in(gens);
    auto pg_max_sq = grid.pg_max_sq.in(gens);
    auto pg_min_sq = grid.pg_min_sq.in(gens);
    auto pf_from_min=grid.pf_from_min.in(arcs);
    auto pf_from_max=grid.pf_from_max.in(arcs);
    auto qf_from_min=grid.qf_from_min.in(arcs);
    auto qf_from_max=grid.qf_from_max.in(arcs);
    auto pf_to_min=grid.pf_to_min.in(arcs);
    auto pf_to_max=grid.pf_to_max.in(arcs);
    auto qf_to_min=grid.qf_to_min.in(arcs);
    auto qf_to_max=grid.qf_to_max.in(arcs);
    auto lij_min=grid.lij_min.in(arcs);
    auto lij_max=grid.lij_max.in(arcs);
    auto lji_min=grid.lji_min.in(arcs);
    auto lji_max=grid.lji_max.in(arcs);
    auto Iij_min=grid.Iij_min.in(arcs);
    auto Iij_max=grid.Iij_max.in(arcs);
    auto c1 = grid.c1.in(gens);
    auto c2 = grid.c2.in(gens);
    auto c0 = grid.c0.in(gens);
    auto pl = grid.pl.in(nodes);
    auto ql = grid.ql.in(nodes);
    auto gs = grid.gs.in(nodes);
    auto bs = grid.bs.in(nodes);
    auto b = grid.b.in(arcs);
    auto g = grid.g.in(arcs);
    auto as = grid.as.in(arcs);
    auto ch = grid.ch.in(arcs);
    auto tr = grid.tr.in(arcs);
    auto th_min = grid.th_min.in(node_pairs);
    auto th_max = grid.th_max.in(node_pairs);
    auto g_ft = grid.g_ft.in(arcs);
    auto g_ff = grid.g_ff.in(arcs);
    auto g_tt = grid.g_tt.in(arcs);
    auto g_tf = grid.g_tf.in(arcs);
    auto b_ft = grid.b_ft.in(arcs);
    auto b_ff = grid.b_ff.in(arcs);
    auto b_tf = grid.b_tf.in(arcs);
    auto b_tt = grid.b_tt.in(arcs);
    auto S_max = grid.S_max.in(arcs);
    auto v_max = grid.v_max.in(nodes);
    auto v_min = grid.v_min.in(nodes);
    auto w_max = grid.w_max.in(nodes);
    auto w_min = grid.w_min.in(nodes);
    auto tan_th_min = grid.tan_th_min.in(node_pairs);
    auto tan_th_max = grid.tan_th_max.in(node_pairs);
    auto wr_min = grid.wr_min.in(node_pairs_chord);
    auto wr_max = grid.wr_max.in(node_pairs_chord);
    auto wi_min = grid.wi_min.in(node_pairs_chord);
    auto wi_max = grid.wi_max.in(node_pairs_chord);
    auto cc=grid.cc.in(arcs);
    auto dd=grid.dd.in(arcs);
    auto ch_half=grid.ch_half.in(arcs);
    
    
    
    auto SDPOPF = make_shared<Model<>>("SDP-OPF Model");
    
    /** Variables */
    /* Power generation variables */
    var<> Pg("Pg", pg_min, pg_max);
    var<> Qg ("Qg", qg_min, qg_max);
    SDPOPF->add(Pg.in(gens),Qg.in(gens));
    
    
    var<> Pf_from("Pf_from", pf_from_min,pf_from_max);
    var<> Qf_from("Qf_from", qf_from_min,qf_from_max);
    var<> Pf_to("Pf_to", pf_to_min,pf_to_max);
    var<> Qf_to("Qf_to", qf_to_min,qf_to_max);
    
    
    SDPOPF->add(Pf_from.in(arcs), Qf_from.in(arcs),Pf_to.in(arcs),Qf_to.in(arcs));
    //Pf_to._off=pf_to_min._off;
    
    /* Real part of Wij = ViVj */
    var<>  R_Wij("R_Wij", wr_min, wr_max);
    /* Imaginary part of Wij = ViVj */
    var<>  Im_Wij("Im_Wij", wi_min, wi_max);
    /* Magnitude of Wii = Vi^2 */
    var<>  Wii("Wii", w_min, w_max);
    SDPOPF->add(Wii.in(nodes),R_Wij.in(node_pairs_chord),Im_Wij.in(node_pairs_chord));
    
//    add_original=true;
    if(add_original)
    {
  
        
        var<>  R_Vi("R_Vi", -1*v_max, v_max);
        var<>  Im_Vi("Im_Vi", -1*v_max, v_max);
        
        
        SDPOPF->add(R_Vi.in(nodes),Im_Vi.in(nodes));
        R_Vi.initialize_all(1);
        Im_Vi.set_lb((grid.ref_bus),0);
        Im_Vi.set_ub((grid.ref_bus),0);
        
//        R_Vi.set_lb((grid.ref_bus),v_min(grid.ref_bus).eval());
//        R_Vi.set_ub((grid.ref_bus),v_max(grid.ref_bus).eval());
        
        
        var<Cpx> Vi("Vi"), Vj("Vj"), Wij("Wij");
        Vi.real_imag(R_Vi.from(node_pairs_chord), Im_Vi.from(node_pairs_chord));
        Vj.real_imag(R_Vi.to(node_pairs_chord), Im_Vi.to(node_pairs_chord));
        Wij.real_imag(R_Wij.in(node_pairs_chord), Im_Wij.in(node_pairs_chord));
        
        Constraint<Cpx> Linking_Wij("Linking_Wij");
        Linking_Wij = Wij - Vi*conj(Vj);
        SDPOPF->add(Linking_Wij.in(node_pairs_chord)==0, convexify);
        
        Vi.real_imag(R_Vi.in(nodes), Im_Vi.in(nodes));
        
        Constraint<Cpx> Linking_Wi("Linking_Wi");
        Linking_Wi = Wii - Vi*conj(Vi);
        SDPOPF->add(Linking_Wi.in(nodes)==0, convexify);
        
        //        if(!grid._tree)
        //        {
        //
        //            auto Wij_ = Wij.in(node_pairs_chord).pairs_in_bags(grid._bags, 3);
        //            auto Wii_ = Wii.in_bags(grid._bags, 3);
        //            auto nb_bags3 = Wij_[0]._indices->size();
        //
        //            Constraint<Cpx> Rank_type2a("RankType2a");
        //            Rank_type2a=Wij_[0]*Wij_[1]-Wii_[1]*Wij_[2];
        //            SDPOPF->add(Rank_type2a.in(range(1,nb_bags3))==0, true);
        //
        //            Constraint<Cpx> Rank_type2b("RankType2b");
        //            Rank_type2b=Wij_[2]*conj(Wij_[1])-Wii_[2]*Wij_[0];
        //            SDPOPF->add(Rank_type2b.in(range(1,nb_bags3))==0, true);
        //
        //            Constraint<Cpx> Rank_type2c("RankType2c");
        //            Rank_type2c=Wij_[2]*conj(Wij_[0])-Wii_[0]*Wij_[1];
        //            SDPOPF->add(Rank_type2c.in(range(1,nb_bags3))==0, true);
        //        }
    }
    //
    
    /* Initialize variables */
    R_Wij.initialize_all(1.0);
    Wii.initialize_all(1.0);
    
    var<> lij("lij", lij_min,lij_max);
    var<> lji("lji", lji_min,lji_max);
    //
    //    var<> eta("eta", 0, 1);
    //    SDPOPF->add(eta.in(range(0,0)));
    
    var<> etag("etag", pg_min_sq, pg_max_sq);
    if(!nonlin_obj){
        SDPOPF->add(etag.in(gens));
    }
    
    
    if(current){
        SDPOPF->add(lij.in(arcs));
        SDPOPF->add(lji.in(arcs));
    }
    
    
    
    //    func<> obj = (product(c1,Pg) + product(c2,pow(Pg,2)) + sum(c0));
    //    SDPOPF->min(obj);
    
    //func<> obj=(product(c1,Pg) + product(c2,pow(Pg,2)) + sum(c0))/upper_bound;
    //SDPOPF->min(eta);
    
    /**  Objective */
    
    
    if(nonlin_obj)
    {
        //    Constraint<> obj_UB("obj_UB");
        //    obj_UB=(product(c1,Pg) + product(c2,pow(Pg,2)) + sum(c0))-eta*upper_bound;
        //    SDPOPF->add(obj_UB.in(range(0,0))<=0);
        
        
        auto obj=(product(c1,Pg) + product(c2,pow(Pg,2)) + sum(c0));
        SDPOPF->min(obj);
        
        
//        Constraint<> obj_UB("obj_UB");
//        obj_UB=(product(c1,Pg) + product(c2,pow(Pg,2)) + sum(c0))/upper_bound-1.0;
//        SDPOPF->add(obj_UB.in(range(0,0))<=0);
        //        Constraint<> obj_UB("obj_UB");
        //        obj_UB=(product(c1,Pg) + product(c2,pow(Pg,2)) + sum(c0));
        //        SDPOPF->add(obj_UB.in(range(0,0))<=upper_bound);
    }
    else
    {
        Constraint<> obj_cost("obj_cost");
        obj_cost=etag-pow(Pg,2);
        SDPOPF->add(obj_cost.in(gens)>=0);
        
        
        auto obj=(product(c1,Pg) + product(c2,etag) + sum(c0));
        SDPOPF->min(obj);
        
        
  
        
        
        //        Constraint<> obj_UB("obj_UB");
        //        obj_UB=(product(c1,Pg) + product(c2,etag) + sum(c0));
        //        SDPOPF->add(obj_UB.in(range(0,0))<=upper_bound);
        
        
        
        
        //        Constraint<> obj_UB("obj_UB");
        //        obj_UB=(product(c1,Pg) + product(c2,pow(Pg,2)) + sum(c0))/upper_bound-eta;
        //        SDPOPF->add(obj_UB.in(range(0,0))==0, convexify, "on/off", false);
        //
        //    }
    }
    
    /** Constraints */
    if(!grid._tree && grid.add_3d_nlin && sdp_cuts) {
        
        auto bag_size = bags_3d.size();
        DebugOn("\nNum of bags = " << bag_size << endl);
        DebugOn("Adding 3d determinant polynomial cuts\n");
        auto R_Wij_ = R_Wij.pairs_in_bags(bags_3d, 3);
        auto Im_Wij_ = Im_Wij.pairs_in_bags(bags_3d, 3);
        auto Wii_ = Wii.in_bags(bags_3d, 3);
        
        if(sdp_kim){
            
            
            Constraint<> SDP3("SDP_3D");
            SDP3 = 2 * R_Wij_[0] * (R_Wij_[1] * R_Wij_[2] + Im_Wij_[1] * Im_Wij_[2]);
            SDP3 -= 2 * Im_Wij_[0] * (R_Wij_[2] * Im_Wij_[1] - Im_Wij_[2] * R_Wij_[1]);
            SDP3 -= (pow(R_Wij_[0], 2) + pow(Im_Wij_[0], 2)) * Wii_[2];
            SDP3 -= (pow(R_Wij_[1], 2) + pow(Im_Wij_[1], 2)) * Wii_[0];
            SDP3 -= (pow(R_Wij_[2], 2) + pow(Im_Wij_[2], 2)) * Wii_[1];
            SDP3 += Wii_[0] * Wii_[1] * Wii_[2];
            
            if (lazy_bool) {
                SDPOPF->add_lazy(SDP3.in(range(0,bag_size-1)) >= 0);
                
            }
            else {
                SDPOPF->add(SDP3.in(range(0,bag_size-1)) >= 0);
                DebugOn("Number of 3d determinant cuts = " << SDP3.get_nb_instances() << endl);
            }
        }
        else{
            indices theta_kim("theta_kim");
            param<double> theta;
            theta.in(theta_kim);
            vector<double> theta_val;
            theta_val.push_back(0);
            theta_val.push_back(pi/2.0);
            theta_val.push_back(-pi/2.0);
            
            for(auto i=0;i<theta_val.size();i++){
                theta_kim.add(to_string(i));
                theta.add_val(theta_val[i]);
            }
            indices bags("bags");
            bags=range(1, bag_size);
            
            auto bags_theta=indices(bags,theta_kim);
            R_Wij_[0] = R_Wij_[0].in_ignore_ith(2, 1, indices(*R_Wij_[0]._indices,theta_kim));
            R_Wij_[1] = R_Wij_[1].in_ignore_ith(2, 1, indices(*R_Wij_[1]._indices,theta_kim));
            R_Wij_[2] = R_Wij_[2].in_ignore_ith(2, 1, indices(*R_Wij_[2]._indices,theta_kim));
            
            Im_Wij_[0] = Im_Wij_[0].in_ignore_ith(2, 1, indices(*Im_Wij_[0]._indices,theta_kim));
            Im_Wij_[1] = Im_Wij_[1].in_ignore_ith(2, 1, indices(*Im_Wij_[1]._indices,theta_kim));
            Im_Wij_[2] = Im_Wij_[2].in_ignore_ith(2, 1, indices(*Im_Wij_[2]._indices,theta_kim));
            
            Wii_[0] = Wii_[0].in_ignore_ith(1, 1, indices(*Wii_[0]._indices,theta_kim));
            Wii_[1] = Wii_[1].in_ignore_ith(1, 1, indices(*Wii_[1]._indices,theta_kim));
            Wii_[2] = Wii_[2].in_ignore_ith(1, 1, indices(*Wii_[2]._indices,theta_kim));
            
            theta = theta.in_ignore_ith(0, 2, indices(*R_Wij_[0]._indices,theta_kim));
            
            Constraint<> SOC_Kojima1_theta("SOC_Kojima1_theta");
            SOC_Kojima1_theta = pow(R_Wij_[0] + cos(theta)*R_Wij_[2]-sin(theta)*Im_Wij_[2], 2)+pow(Im_Wij_[0] + cos(theta)*Im_Wij_[2]+sin(theta)*R_Wij_[2], 2)-Wii_[0]*(Wii_[1]+Wii_[2]+2*(cos(theta)*R_Wij_[1]-sin(theta)*Im_Wij_[1]));
            SDPOPF->add(SOC_Kojima1_theta.in(bags_theta) <= 0);
            
            
            Constraint<> SOC_Kojima2_theta("SOC_Kojima2_theta");
            SOC_Kojima2_theta = pow(R_Wij_[0] + cos(theta)*R_Wij_[1]-sin(theta)*Im_Wij_[1], 2)+pow(Im_Wij_[0] - cos(theta)*Im_Wij_[1]-sin(theta)*R_Wij_[1], 2)-Wii_[1]*(Wii_[0]+Wii_[2]+2*(cos(theta)*R_Wij_[2]-sin(theta)*Im_Wij_[2]));
            SDPOPF->add(SOC_Kojima2_theta.in(bags_theta) <= 0);
            
            
            Constraint<> SOC_Kojima3_theta("SOC_Kojima3_theta");
            SOC_Kojima3_theta = pow(R_Wij_[2] + cos(theta)*R_Wij_[1]+sin(theta)*Im_Wij_[1], 2)+pow(Im_Wij_[2] + cos(theta)*Im_Wij_[1]-sin(theta)*R_Wij_[1], 2)-Wii_[2]*(Wii_[0]+Wii_[1]+2*(cos(theta)*R_Wij_[0]-sin(theta)*Im_Wij_[0]));
            SDPOPF->add(SOC_Kojima3_theta.in(bags_theta) <= 0);
            
        }
        
        //        /* Second-order cone constraints */
        //        Constraint<> SOC_Kojima1_0("SOC_Kojima1_0");
        //        SOC_Kojima1_0 = pow(R_Wij_[0] + R_Wij_[2], 2) + pow(Im_Wij_[0] + Im_Wij_[2], 2) - Wii_[0]*(Wii_[1]+Wii_[2]+2*R_Wij_[1]);
        //        SDPOPF->add(SOC_Kojima1_0.in(range(0,bag_size-1)) <= 0);
        ////
        ////
        //        /* Second-order cone constraints */
        //        Constraint<> SOC_Kojima2_0("SOC_Kojima2_0");
        //        SOC_Kojima2_0 = pow(R_Wij_[0] + R_Wij_[1], 2) + pow(Im_Wij_[0] - Im_Wij_[1], 2) - Wii_[1]*(Wii_[0]+Wii_[2]+2*R_Wij_[2]);
        //        SDPOPF->add(SOC_Kojima2_0.in(range(0,bag_size-1)) <= 0);
        //
        //
        //        Constraint<> SOC_Kojima3_0("SOC_Kojima3_0");
        //        SOC_Kojima3_0 = pow(R_Wij_[2] + R_Wij_[1], 2) + pow(Im_Wij_[2] + Im_Wij_[1], 2) - Wii_[2]*(Wii_[0]+Wii_[1]+2*R_Wij_[0]);
        //        SDPOPF->add(SOC_Kojima3_0.in(range(0,bag_size-1)) <= 0);
        //
        //        Constraint<> SOC_Kojima1_90("SOC_Kojima1_90");
        //        SOC_Kojima1_90 = pow(R_Wij_[0] - Im_Wij_[2], 2) + pow(Im_Wij_[0] + R_Wij_[2], 2) - Wii_[0]*(Wii_[1]+Wii_[2]-2*Im_Wij_[1]);
        //        SDPOPF->add(SOC_Kojima1_90.in(range(0,bag_size-1)) <= 0);
        //
        //
        //       Constraint<> SOC_Kojima2_90("SOC_Kojima2_90");
        //       SOC_Kojima2_90 = pow(R_Wij_[0] - Im_Wij_[1], 2) + pow(Im_Wij_[0] - R_Wij_[1], 2) - Wii_[1]*(Wii_[0]+Wii_[2]-2*Im_Wij_[2]);
        //       SDPOPF->add(SOC_Kojima2_90.in(range(0,bag_size-1)) <= 0);
        //
        //        Constraint<> SOC_Kojima3_90("SOC_Kojima3_90");
        //        SOC_Kojima3_90 = pow(R_Wij_[2] + Im_Wij_[1], 2) + pow(Im_Wij_[2] - R_Wij_[1], 2) - Wii_[2]*(Wii_[0]+Wii_[1]-2*Im_Wij_[0]);
        //        SDPOPF->add(SOC_Kojima3_90.in(range(0,bag_size-1)) <= 0);
        //
        //        const double root2=sqrt(2.0);
        //
        //        Constraint<> SOC_Kojima1_45("SOC_Kojima1_45");
        //        SOC_Kojima1_45 = pow(root2*R_Wij_[0] + R_Wij_[2] -Im_Wij_[2], 2) + pow(root2*Im_Wij_[0] + Im_Wij_[2]+ R_Wij_[2], 2) - 2.0*Wii_[0]*(Wii_[1]+Wii_[2]+root2*(R_Wij_[1]-Im_Wij_[1]));
        //        SDPOPF->add(SOC_Kojima1_45.in(range(0,bag_size-1)) <= 0);
        //
        //
        //        /* Second-order cone constraints */
        //        Constraint<> SOC_Kojima2_45("SOC_Kojima2_45");
        //        SOC_Kojima2_45 = pow(root2*R_Wij_[0] + R_Wij_[1] -Im_Wij_[1], 2) + pow(root2*Im_Wij_[0] - Im_Wij_[1]-R_Wij_[1], 2) - 2.0*Wii_[1]*(Wii_[0]+Wii_[2]+root2*(R_Wij_[2]-Im_Wij_[2]));
        //        SDPOPF->add(SOC_Kojima2_45.in(range(0,bag_size-1)) <= 0);
        //
        //
        //        Constraint<> SOC_Kojima3_45("SOC_Kojima3_45");
        //        SOC_Kojima3_45 = pow(root2*R_Wij_[2] + R_Wij_[1] +Im_Wij_[1], 2) + pow(root2*Im_Wij_[2] + Im_Wij_[1]-   R_Wij_[1], 2) - 2.0*Wii_[2]*(Wii_[0]+Wii_[1]+root2*(R_Wij_[0]-Im_Wij_[0]));
        //        SDPOPF->add(SOC_Kojima3_45.in(range(0,bag_size-1)) <= 0);
        
        
        
        
        //        Constraint<> SOC_Kojima1_0_NC("SOC_Kojima1_0_NC");
        //         SOC_Kojima1_0_NC = pow(R_Wij_[0] + R_Wij_[2], 2) + pow(Im_Wij_[0] + Im_Wij_[2], 2) - Wii_[0]*(Wii_[1]+Wii_[2]+2*R_Wij_[1]);
        //         SDPOPF->add(SOC_Kojima1_0_NC.in(range(0,bag_size-1)) >= 0,true);
        //
        //         /* Second-order cone constraints */
        //         Constraint<> SOC_Kojima2_0_NC("SOC_Kojima2_0_NC");
        //         SOC_Kojima2_0_NC = pow(R_Wij_[0] + R_Wij_[1], 2) + pow(Im_Wij_[0] - Im_Wij_[1], 2) - Wii_[1]*(Wii_[0]+Wii_[2]+2*R_Wij_[2]);
        //         SDPOPF->add(SOC_Kojima2_0_NC.in(range(0,bag_size-1)) >= 0,true);
        //
        //
        //         Constraint<> SOC_Kojima3_0_NC("SOC_Kojima3_0_NC");
        //         SOC_Kojima3_0_NC = pow(R_Wij_[2] + R_Wij_[1], 2) + pow(Im_Wij_[2] + Im_Wij_[1], 2) - Wii_[2]*(Wii_[0]+Wii_[1]+2*R_Wij_[0]);
        //         SDPOPF->add(SOC_Kojima3_0_NC.in(range(0,bag_size-1)) >= 0,true);
        //
        //         Constraint<> SOC_Kojima1_90_NC("SOC_Kojima1_90_NC");
        //         SOC_Kojima1_90_NC = pow(R_Wij_[0] - Im_Wij_[2], 2) + pow(Im_Wij_[0] + R_Wij_[2], 2) - Wii_[0]*(Wii_[1]+Wii_[2]-2*Im_Wij_[1]);
        //         SDPOPF->add(SOC_Kojima1_90_NC.in(range(0,bag_size-1)) >= 0,true);
        //
        //        Constraint<> SOC_Kojima2_90_NC("SOC_Kojima2_90_NC");
        //        SOC_Kojima2_90_NC = pow(R_Wij_[0] - Im_Wij_[1], 2) + pow(Im_Wij_[0] - R_Wij_[1], 2) - Wii_[1]*(Wii_[0]+Wii_[2]-2*Im_Wij_[2]);
        //        SDPOPF->add(SOC_Kojima2_90_NC.in(range(0,bag_size-1)) >= 0,true);
        //
        //         Constraint<> SOC_Kojima3_90_NC("SOC_Kojima3_90_NC");
        //         SOC_Kojima3_90_NC = pow(R_Wij_[2] + Im_Wij_[1], 2) + pow(Im_Wij_[2] - R_Wij_[1], 2) - Wii_[2]*(Wii_[0]+Wii_[1]-2*Im_Wij_[0]);
        //         SDPOPF->add(SOC_Kojima3_90_NC.in(range(0,bag_size-1)) >= 0,true);
        
        //        SDPOPF->print();
        //
        //
    }
    
    /** Constraints */
    /* Second-order cone constraints */
    Constraint<> SOC("SOC");
    SOC = pow(R_Wij, 2) + pow(Im_Wij, 2) - Wii.from(node_pairs_chord)*Wii.to(node_pairs_chord);
    SDPOPF->add(SOC.in(node_pairs_chord) == 0,true);
    
    
    
    /* Flow conservation */
    Constraint<> KCL_P("KCL_P");
    KCL_P  = sum(Pf_from, out_arcs) + sum(Pf_to, in_arcs) + pl - sum(Pg, gen_nodes) + gs*Wii;
    SDPOPF->add(KCL_P.in(nodes) == 0);
    
    
    Constraint<> KCL_Q("KCL_Q");
    KCL_Q  = sum(Qf_from, out_arcs) + sum(Qf_to, in_arcs) + ql - sum(Qg, gen_nodes) - bs*Wii;
    SDPOPF->add(KCL_Q.in(nodes) == 0);
    
    
    /* AC Power Flow */
    Constraint<> Flow_P_From("Flow_P_From");
    Flow_P_From = Pf_from - (g_ff*Wii.from(arcs) + g_ft*R_Wij.in_pairs(arcs) + b_ft*Im_Wij.in_pairs(arcs));
    SDPOPF->add(Flow_P_From.in(arcs) == 0);
    
    
    Constraint<> Flow_P_To("Flow_P_To");
    Flow_P_To = Pf_to - (g_tt*Wii.to(arcs) + g_tf*R_Wij.in_pairs(arcs) - b_tf*Im_Wij.in_pairs(arcs));
    SDPOPF->add(Flow_P_To.in(arcs) == 0);
    
    
    Constraint<> Flow_Q_From("Flow_Q_From");
    Flow_Q_From = Qf_from - (g_ft*Im_Wij.in_pairs(arcs) - b_ff*Wii.from(arcs) - b_ft*R_Wij.in_pairs(arcs));
    SDPOPF->add(Flow_Q_From.in(arcs) == 0);
    
    
    Constraint<> Flow_Q_To("Flow_Q_To");
    Flow_Q_To = Qf_to + b_tt*Wii.to(arcs) + b_tf*R_Wij.in_pairs(arcs) + g_tf*Im_Wij.in_pairs(arcs);
    SDPOPF->add(Flow_Q_To.in(arcs) == 0);
    
    
    /* Phase Angle Bounds constraints */
    Constraint<> PAD_UB("PAD_UB");
    PAD_UB = Im_Wij.in(node_pairs);
    PAD_UB <= tan_th_max*R_Wij.in(node_pairs);
    SDPOPF->add(PAD_UB.in(node_pairs));
    
    
    Constraint<> PAD_LB("PAD_LB");
    PAD_LB =  Im_Wij.in(node_pairs);
    PAD_LB >= tan_th_min*R_Wij.in(node_pairs);
    SDPOPF->add(PAD_LB.in(node_pairs));
    
    /* Thermal Limit Constraints */
    Constraint<> Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from = pow(Pf_from, 2) + pow(Qf_from, 2);
    Thermal_Limit_from <= pow(S_max,2);
    // SDPOPF->add(Thermal_Limit_from.in(arcs));
    SDPOPF->add(Thermal_Limit_from.in(arcs), true);
    
    
    Constraint<> Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to = pow(Pf_to, 2) + pow(Qf_to, 2);
    Thermal_Limit_to <= pow(S_max,2);
    //SDPOPF->add(Thermal_Limit_to.in(arcs));
    SDPOPF->add(Thermal_Limit_to.in(arcs), true);
    
    if(llnc)
    {
        
        func<> theta_L = atan(min(Im_Wij.get_lb().in(node_pairs)/R_Wij.get_ub().in(node_pairs),Im_Wij.get_lb().in(node_pairs)/R_Wij.get_lb().in(node_pairs)));
        func<> theta_U = atan(max(Im_Wij.get_ub().in(node_pairs)/R_Wij.get_lb().in(node_pairs),Im_Wij.get_ub().in(node_pairs)/R_Wij.get_ub().in(node_pairs)));
        func<> phi=(theta_U.in(node_pairs)+theta_L.in(node_pairs))/2.0;
        func<> del=(theta_U.in(node_pairs)-theta_L.in(node_pairs))/2.0;
        
        
        Constraint<> LNC1("LNC1");
        LNC1 += (sqrt(Wii.get_lb().from(node_pairs))+sqrt(Wii.get_ub().from(node_pairs)))*(sqrt(Wii.get_lb().to(node_pairs))+sqrt(Wii.get_ub().to(node_pairs)))*(Im_Wij.in(node_pairs)*sin(phi.in(node_pairs)) + R_Wij.in(node_pairs)*cos(phi.in(node_pairs)));
        
        LNC1 -=sqrt(Wii.get_ub().to(node_pairs))*cos(del.in(node_pairs))*(sqrt(Wii.get_lb().to(node_pairs))+sqrt(Wii.get_ub().to(node_pairs)))*Wii.from(node_pairs);
        
        LNC1 -=sqrt(Wii.get_ub().from(node_pairs))*cos(del.in(node_pairs))*(sqrt(Wii.get_lb().from(node_pairs))+sqrt(Wii.get_ub().from(node_pairs)))*Wii.to(node_pairs);
        
        LNC1-=sqrt(Wii.get_ub().from(node_pairs))*sqrt(Wii.get_ub().to(node_pairs))*cos(del)*(sqrt(Wii.get_lb().from(node_pairs))*
                                                                                              sqrt(Wii.get_lb().to(node_pairs)) - sqrt(Wii.get_ub().from(node_pairs))*sqrt(Wii.get_ub().to(node_pairs)));
        SDPOPF->add(LNC1.in(node_pairs) >= 0);
        
        Constraint<> LNC2("LNC2");
        LNC2 += (sqrt(Wii.get_lb().from(node_pairs))+sqrt(Wii.get_ub().from(node_pairs)))*(sqrt(Wii.get_lb().to(node_pairs))+sqrt(Wii.get_ub().to(node_pairs)))*(sin(phi.in(node_pairs))*Im_Wij.in(node_pairs) + cos(phi.in(node_pairs))*R_Wij.in(node_pairs));
        LNC2 -=sqrt(Wii.get_lb().to(node_pairs))*cos(del.in(node_pairs))*(sqrt(Wii.get_lb().to(node_pairs))+sqrt(Wii.get_ub().to(node_pairs)))*Wii.from(node_pairs);
        LNC2 -=sqrt(Wii.get_lb().from(node_pairs))*cos(del.in(node_pairs))*(sqrt(Wii.get_lb().from(node_pairs))+sqrt(Wii.get_ub().from(node_pairs)))*Wii.to(node_pairs);
        LNC2 -=sqrt(Wii.get_lb().from(node_pairs))*sqrt(Wii.get_lb().to(node_pairs))*cos(del.in(node_pairs))*(sqrt(Wii.get_ub().from(node_pairs))*
                                                                                                              sqrt(Wii.get_ub().to(node_pairs))-sqrt(Wii.get_lb().from(node_pairs))*sqrt(Wii.get_lb().to(node_pairs)));
        SDPOPF->add(LNC2.in(node_pairs) >= 0);
    }
    
    if(current){
        param<Cpx> T("T"), Y("Y"), Ych("Ych");
        var<Cpx> L_from("L_from"), W("W"), Vi("Vi"), Vj("Vj"), I("I");
        T.real_imag(cc.in(arcs), dd.in(arcs));
        Y.real_imag(g.in(arcs), b.in(arcs));
        Ych.set_imag(ch_half.in(arcs));
        
        
        L_from.set_real(lij.in(arcs));
        W.real_imag(R_Wij.in_pairs(arcs), Im_Wij.in_pairs(arcs));
        
        
        Constraint<Cpx> I_from("I_from");
        I_from=(Y+Ych)*(conj(Y)+conj(Ych))*Wii.from(arcs)-T*Y*(conj(Y)+conj(Ych))*conj(W)-conj(T)*conj(Y)*(Y+Ych)*W+pow(tr,2)*Y*conj(Y)*Wii.to(arcs)-pow(tr,2)*L_from;
        SDPOPF->add_real(I_from.in(arcs)==0);
        
        
        Constraint<> I_from_Pf("I_from_Pf");
        I_from_Pf=lij*Wii.from(arcs)-pow(tr,2)*(pow(Pf_from,2) + pow(Qf_from,2));
        //SDPOPF->add(I_from_Pf.in(arcs)<=0, true);
        SDPOPF->add(I_from_Pf.in(arcs)==0, true, "on/off", false);
        
        var<Cpx> L_to("L_to");
        L_to.set_real(lji.in(arcs));
        
        Constraint<Cpx> I_to("I_to");
        I_to=pow(tr,2)*(Y+Ych)*(conj(Y)+conj(Ych))*Wii.to(arcs)-conj(T)*Y*(conj(Y)+conj(Ych))*W-T*conj(Y)*(Y+Ych)*conj(W)+Y*conj(Y)*Wii.from(arcs)-pow(tr,2)*L_to;
        //SDPOPF->add_real(I_to.in(arcs_I_to)==0);
        SDPOPF->add_real(I_to.in(arcs)==0);
        
        Constraint<> I_to_Pf("I_to_Pf");
        I_to_Pf=lji*Wii.to(arcs)-(pow(Pf_to,2) + pow(Qf_to, 2));
        //SDPOPF->add(I_to_Pf.in(arcs)<=0, true);
        SDPOPF->add(I_to_Pf.in(arcs)==0, true, "on/off", false);
        
        
    }
    
    return SDPOPF;
    
}

//shared_ptr<Model<>> build_SDPOPF_linear(PowerNet& grid, double upper_bound) {
//    int output = 0;
//    bool sdp_cuts = true;
//
//    bool current_from = true, llnc=true, current_to=true, loss=true, loss_bounds=true;
//
//    size_t num_bags = 0;
//    string num_bags_s = "100";
//    string solver_str = "ipopt";
//    string sdp_cuts_s = "yes";
//    string current_from_s = "yes";
//    string orig_s = "yes";
//    string current_to_s="yes";
//    string lazy_s = "no";
//    bool lazy_bool = false;
//    SolverType solv_type = ipopt;
//    double tol = 1e-6;
//    string mehrotra = "no";
//
//
//
//
//    cout << "\nnum bags = " << num_bags << endl;
//
//    double total_time_start = get_wall_time();
//
//    grid.update_ref_bus();
//
//    grid.get_tree_decomp_bags();
//    auto bags_3d=grid.decompose_bags_3d();
//
//
//    /* Grid Stats */
//    auto nb_gen = grid.get_nb_active_gens();
//    auto nb_lines = grid.get_nb_active_arcs();
//    auto nb_buses = grid.get_nb_active_nodes();
//    DebugOn("nb active gens = " << nb_gen << endl);
//    DebugOn("nb active lines = " << nb_lines << endl);
//    DebugOn("nb active buses = " << nb_buses << endl);
//
//    /** Sets */
//    auto node_pairs = grid.get_node_pairs();
//    auto node_pairs_chord = grid.get_node_pairs_chord(bags_3d);
//    if (grid._tree || !grid.add_3d_nlin || !sdp_cuts) {
//        node_pairs_chord = node_pairs;
//    }
//    auto nodes = indices(grid.nodes);
//    auto arcs = indices(grid.arcs);
//    auto gens = indices(grid.gens);
//    auto gen_nodes = grid.gens_per_node();
//    auto out_arcs = grid.out_arcs_per_node();
//    auto in_arcs = grid.in_arcs_per_node();
//
//    /* Grid Parameters */
//    auto pg_min = grid.pg_min.in(gens);
//    auto pg_max = grid.pg_max.in(gens);
//    auto qg_min = grid.qg_min.in(gens);
//    auto qg_max = grid.qg_max.in(gens);
//    auto c1 = grid.c1.in(gens);
//    auto c2 = grid.c2.in(gens);
//    auto c0 = grid.c0.in(gens);
//    auto pl = grid.pl.in(nodes);
//    auto ql = grid.ql.in(nodes);
//    auto gs = grid.gs.in(nodes);
//    auto bs = grid.bs.in(nodes);
//    auto b = grid.b.in(arcs);
//    auto g = grid.g.in(arcs);
//    auto as = grid.as.in(arcs);
//    auto ch = grid.ch.in(arcs);
//    auto tr = grid.tr.in(arcs);
//    auto th_min = grid.th_min.in(node_pairs);
//    auto th_max = grid.th_max.in(node_pairs);
//    auto g_ft = grid.g_ft.in(arcs);
//    auto g_ff = grid.g_ff.in(arcs);
//    auto g_tt = grid.g_tt.in(arcs);
//    auto g_tf = grid.g_tf.in(arcs);
//    auto b_ft = grid.b_ft.in(arcs);
//    auto b_ff = grid.b_ff.in(arcs);
//    auto b_tf = grid.b_tf.in(arcs);
//    auto b_tt = grid.b_tt.in(arcs);
//    auto S_max = grid.S_max.in(arcs);
//    auto v_max = grid.v_max.in(nodes);
//    auto v_min = grid.v_min.in(nodes);
//    auto w_max = grid.w_max.in(nodes);
//    auto w_min = grid.w_min.in(nodes);
//    auto tan_th_min = grid.tan_th_min.in(node_pairs);
//    auto tan_th_max = grid.tan_th_max.in(node_pairs);
//    auto wr_min = grid.wr_min.in(node_pairs_chord);
//    auto wr_max = grid.wr_max.in(node_pairs_chord);
//    auto wi_min = grid.wi_min.in(node_pairs_chord);
//    auto wi_max = grid.wi_max.in(node_pairs_chord);
//    auto lij_min=grid.lij_min.in(arcs);
//    auto lij_max=grid.lij_max.in(arcs);
//    auto cc=grid.cc.in(arcs);
//    auto dd=grid.dd.in(arcs);
//    auto ch_half=grid.ch_half.in(arcs);
//    auto arcs_inductive=grid.arcs_inductive_only();
//    auto lji_min=grid.lji_min.in(arcs);
//    auto lji_max=grid.lji_max.in(arcs);
//
//
//    /** Build model */
//    Model<> SDP("SDP Model");
//    auto SDPOA = make_shared<Model<>>("SDP-OA Model");
//
//    /** Variables */
//    /* Power generation variables */
//    var<> Pg("Pg", pg_min, pg_max);
//    var<> Qg ("Qg", qg_min, qg_max);
//    SDP.add(Pg.in(gens),Qg.in(gens));
//    SDPOA->add(Pg.in(gens),Qg.in(gens));
//
//    /* Power flow variables */
//    var<> Pf_from("Pf_from", -1.*S_max,S_max);
//    var<> Qf_from("Qf_from", -1.*S_max,S_max);
//    var<> Pf_to("Pf_to", -1.*S_max,S_max);
//    var<> Qf_to("Qf_to", -1.*S_max,S_max);
//
//    SDP.add(Pf_from.in(arcs), Qf_from.in(arcs),Pf_to.in(arcs),Qf_to.in(arcs));
//    SDPOA->add(Pf_from.in(arcs), Qf_from.in(arcs),Pf_to.in(arcs),Qf_to.in(arcs));
//
//
//    /* Real part of Wij = ViVj */
//    var<>  R_Wij("R_Wij", wr_min, wr_max);
//    /* Imaginary part of Wij = ViVj */
//    var<>  Im_Wij("Im_Wij", wi_min, wi_max);
//    /* Magnitude of Wii = Vi^2 */
//    var<>  Wii("Wii", w_min, w_max);
//    SDP.add(Wii.in(nodes),R_Wij.in(node_pairs_chord),Im_Wij.in(node_pairs_chord));
//    SDPOA->add(Wii.in(nodes),R_Wij.in(node_pairs_chord),Im_Wij.in(node_pairs_chord));
//
//
//    /* Initialize variables */
//    R_Wij.initialize_all(1.0);
//    Wii.initialize_all(1.001);
//
//    bool current = true;
//    var<> lij("lij", lij_min,lij_max);
//    var<> lji("lji", lji_min,lji_max);
//    //var<> eta("eta", 0, 10);
//    if(current){
//        SDP.add(lij.in(arcs),lji.in(arcs));
//        SDPOA->add(lij.in(arcs),lji.in(arcs));
//    }
//
//    //SDP.add(eta.in(range(0, 0)));
//    /**  Objective */
//    auto obj = (product(c1,Pg) + product(c2,pow(Pg,2)) + sum(c0));
//    // obj=eta("0")*(-1);
//    SDP.min(obj);
//    SDPOA->min(obj);
//
//
//
//
//    /** Constraints */
//    auto bag_size = bags_3d.size();
//    Constraint<> SDP3("SDP_3D");
//    //      Constraint<> SDPD("SDPD");
//    if(!grid._tree && grid.add_3d_nlin && sdp_cuts)
//    {
//        DebugOn("\nNum of bags = " << bag_size << endl);
//        DebugOn("Adding 3d determinant polynomial cuts\n");
//        auto R_Wij_ = R_Wij.pairs_in_bags(bags_3d, 3);
//        auto Im_Wij_ = Im_Wij.pairs_in_bags(bags_3d, 3);
//        auto Wii_ = Wii.in_bags(bags_3d, 3);
//
//
//
//        SDP3 = 2 * R_Wij_[0] * (R_Wij_[1] * R_Wij_[2] + Im_Wij_[1] * Im_Wij_[2]);
//        SDP3 -= 2 * Im_Wij_[0] * (R_Wij_[2] * Im_Wij_[1] - Im_Wij_[2] * R_Wij_[1]);
//        SDP3 -= (pow(R_Wij_[0], 2) + pow(Im_Wij_[0], 2)) * Wii_[2];
//        SDP3 -= (pow(R_Wij_[1], 2) + pow(Im_Wij_[1], 2)) * Wii_[0];
//        SDP3 -= (pow(R_Wij_[2], 2) + pow(Im_Wij_[2], 2)) * Wii_[1];
//        SDP3 += Wii_[0] * Wii_[1] * Wii_[2];
//        if (lazy_bool) {
//            SDP.add_lazy(SDP3.in(range(0, bag_size-1)) >= 0);
//            //SDPOA->add_lazy(SDP3.in(orig) >= 0);
//        }
//        else {
//            SDP.add(SDP3.in(range(0, bag_size-1)) >= 0);
//            // SDPOA->add(SDP3.in(orig) >= 0);
//            DebugOn("Number of 3d determinant cuts = " << SDP3.get_nb_instances() << endl);
//        }
//
//        //
//        //        SDPD= - 0.002491499038*Im_Wij_[0] - 0.002405078286*Im_Wij_[1] + 0.002576479796*Im_Wij_[2] + 0.0006191834985*R_Wij_[0] - 0.002142755127*R_Wij_[1] - 0.004868374713*R_Wij_[2] + 0.002190591032*Wii_[0] + 0.0007469140267*Wii_[1] + 0.003457224761*Wii_[2] + 1.999982567e-08;
//        //        SDP.add(SDPD.in(orig1) >= 0);
//    }
//
//    /** Constraints */
//    /* Second-order cone constraints */
//    Constraint<> SOC("SOC");
//    SOC = pow(R_Wij, 2) + pow(Im_Wij, 2) - Wii.from(node_pairs_chord)*Wii.to(node_pairs_chord);
//    SDP.add(SOC.in(node_pairs_chord) == 0,true);
//    //SDPOA->add(SOC.in(node_pairs_chord) == 0,true);
//
//    /* Flow conservation */
//    Constraint<> KCL_P("KCL_P");
//    KCL_P  = sum(Pf_from, out_arcs) + sum(Pf_to, in_arcs) + pl - sum(Pg, gen_nodes) + gs*Wii;
//    SDP.add(KCL_P.in(nodes) == 0);
//    SDPOA->add(KCL_P.in(nodes) == 0);
//
//    Constraint<> KCL_Q("KCL_Q");
//    KCL_Q  = sum(Qf_from, out_arcs) + sum(Qf_to, in_arcs) + ql - sum(Qg, gen_nodes) - bs*Wii;
//    SDP.add(KCL_Q.in(nodes) == 0);
//    SDPOA->add(KCL_Q.in(nodes) == 0);
//
//    /* AC Power Flow */
//    Constraint<> Flow_P_From("Flow_P_From");
//    Flow_P_From = Pf_from - (g_ff*Wii.from(arcs) + g_ft*R_Wij.in_pairs(arcs) + b_ft*Im_Wij.in_pairs(arcs));
//    SDP.add(Flow_P_From.in(arcs) == 0);
//    SDPOA->add(Flow_P_From.in(arcs) == 0);
//
//    Constraint<> Flow_P_To("Flow_P_To");
//    Flow_P_To = Pf_to - (g_tt*Wii.to(arcs) + g_tf*R_Wij.in_pairs(arcs) - b_tf*Im_Wij.in_pairs(arcs));
//    SDP.add(Flow_P_To.in(arcs) == 0);
//    SDPOA->add(Flow_P_To.in(arcs) == 0);
//
//    Constraint<> Flow_Q_From("Flow_Q_From");
//    Flow_Q_From = Qf_from - (g_ft*Im_Wij.in_pairs(arcs) - b_ff*Wii.from(arcs) - b_ft*R_Wij.in_pairs(arcs));
//    SDP.add(Flow_Q_From.in(arcs) == 0);
//    SDPOA->add(Flow_Q_From.in(arcs) == 0);
//
//    Constraint<> Flow_Q_To("Flow_Q_To");
//    Flow_Q_To = Qf_to + b_tt*Wii.to(arcs) + b_tf*R_Wij.in_pairs(arcs) + g_tf*Im_Wij.in_pairs(arcs);
//    SDP.add(Flow_Q_To.in(arcs) == 0);
//    SDPOA->add(Flow_Q_To.in(arcs) == 0);
//
//    /* Phase Angle Bounds constraints */
//    Constraint<> PAD_UB("PAD_UB");
//    PAD_UB = Im_Wij.in(node_pairs);
//    PAD_UB <= tan_th_max*R_Wij.in(node_pairs);
//    SDP.add(PAD_UB.in(node_pairs));
//    SDPOA->add(PAD_UB.in(node_pairs));
//
//    Constraint<> PAD_LB("PAD_LB");
//    PAD_LB =  Im_Wij.in(node_pairs);
//    PAD_LB >= tan_th_min*R_Wij.in(node_pairs);
//    SDP.add(PAD_LB.in(node_pairs));
//    SDPOA->add(PAD_LB.in(node_pairs));
//
//    /* Thermal Limit Constraints */
//    Constraint<> Thermal_Limit_from("Thermal_Limit_from");
//    Thermal_Limit_from = pow(Pf_from, 2) + pow(Qf_from, 2);
//    Thermal_Limit_from <= pow(S_max,2);
//    SDP.add(Thermal_Limit_from.in(arcs));
//    //SDPOA->add(Thermal_Limit_from.in(arcs));
//
//
//
//    Constraint<> Thermal_Limit_to("Thermal_Limit_to");
//    Thermal_Limit_to = pow(Pf_to, 2) + pow(Qf_to, 2);
//    Thermal_Limit_to <= pow(S_max,2);
//    SDP.add(Thermal_Limit_to.in(arcs));
//    //SDPOA->add(Thermal_Limit_to.in(arcs));
//
//    func<> theta_L = atan(min(Im_Wij.get_lb().in(node_pairs)/R_Wij.get_ub().in(node_pairs),Im_Wij.get_lb().in(node_pairs)/R_Wij.get_lb().in(node_pairs)));
//    func<> theta_U = atan(max(Im_Wij.get_ub().in(node_pairs)/R_Wij.get_lb().in(node_pairs),Im_Wij.get_ub().in(node_pairs)/R_Wij.get_ub().in(node_pairs)));
//    func<> phi=(theta_U.in(node_pairs)+theta_L.in(node_pairs))/2.0;
//    func<> del=(theta_U.in(node_pairs)-theta_L.in(node_pairs))/2.0;
//
//
//    Constraint<> LNC1("LNC1");
//    LNC1 += (sqrt(Wii.get_lb().from(node_pairs))+sqrt(Wii.get_ub().from(node_pairs)))*(sqrt(Wii.get_lb().to(node_pairs))+sqrt(Wii.get_ub().to(node_pairs)))*(Im_Wij.in(node_pairs)*sin(phi.in(node_pairs)) + R_Wij.in(node_pairs)*cos(phi.in(node_pairs)));
//
//    LNC1 -=sqrt(Wii.get_ub().to(node_pairs))*cos(del.in(node_pairs))*(sqrt(Wii.get_lb().to(node_pairs))+sqrt(Wii.get_ub().to(node_pairs)))*Wii.from(node_pairs);
//
//    LNC1 -=sqrt(Wii.get_ub().from(node_pairs))*cos(del.in(node_pairs))*(sqrt(Wii.get_lb().from(node_pairs))+sqrt(Wii.get_ub().from(node_pairs)))*Wii.to(node_pairs);
//
//    LNC1-=sqrt(Wii.get_ub().from(node_pairs))*sqrt(Wii.get_ub().to(node_pairs))*cos(del)*(sqrt(Wii.get_lb().from(node_pairs))*
//                                                                                        sqrt(Wii.get_lb().to(node_pairs)) - sqrt(Wii.get_ub().from(node_pairs))*sqrt(Wii.get_ub().to(node_pairs)));
//    SDP.add(LNC1.in(node_pairs) >= 0);
//    SDPOA->add(LNC1.in(node_pairs) >= 0);
//
//    Constraint<> LNC2("LNC2");
//    LNC2 += (sqrt(Wii.get_lb().from(node_pairs))+sqrt(Wii.get_ub().from(node_pairs)))*(sqrt(Wii.get_lb().to(node_pairs))+sqrt(Wii.get_ub().to(node_pairs)))*(sin(phi.in(node_pairs))*Im_Wij.in(node_pairs) + cos(phi.in(node_pairs))*R_Wij.in(node_pairs));
//    LNC2 -=sqrt(Wii.get_lb().to(node_pairs))*cos(del.in(node_pairs))*(sqrt(Wii.get_lb().to(node_pairs))+sqrt(Wii.get_ub().to(node_pairs)))*Wii.from(node_pairs);
//    LNC2 -=sqrt(Wii.get_lb().from(node_pairs))*cos(del.in(node_pairs))*(sqrt(Wii.get_lb().from(node_pairs))+sqrt(Wii.get_ub().from(node_pairs)))*Wii.to(node_pairs);
//    LNC2 -=sqrt(Wii.get_lb().from(node_pairs))*sqrt(Wii.get_lb().to(node_pairs))*cos(del.in(node_pairs))*(sqrt(Wii.get_ub().from(node_pairs))*
//                                                                                                       sqrt(Wii.get_ub().to(node_pairs))-sqrt(Wii.get_lb().from(node_pairs))*sqrt(Wii.get_lb().to(node_pairs)));
//    SDP.add(LNC2.in(node_pairs) >= 0);
//    SDPOA->add(LNC2.in(node_pairs) >= 0);
//
//    if(current){
//        param<Cpx> T("T"), Y("Y"), Ych("Ych");
//        var<Cpx> L_from("L_from"), Wij("Wij");
//        T.real_imag(cc.in(arcs), dd.in(arcs));
//        Y.real_imag(g.in(arcs), b.in(arcs));
//        Ych.set_imag(ch_half.in(arcs));
//
//
//        L_from.set_real(lij.in(arcs));
//        Wij.real_imag(R_Wij.in_pairs(arcs), Im_Wij.in_pairs(arcs));
//        var<Cpx> Sij("Sij"), Sji("Sji");
//        Sij.real_imag(Pf_from.in(arcs), Qf_from.in(arcs));
//        Sji.real_imag(Pf_to.in(arcs), Qf_to.in(arcs));
//
//        //        Constraint<> PLoss("PLoss");
//        //        PLoss = pow(Pf_from,2);
//        //        PLoss -= pow((g_ff*Wii.from(arcs) + g_ft*R_Wij.in(arcs) + b_ft*Im_Wij.in(arcs)),2);
//        //        SDP.add(PLoss.in(arcs)==0,true);
//        //        Constraint<> PLoss2("PLoss2");
//        //        PLoss2 = pow(Pf_to,2);
//        //        PLoss2 -= pow((g_tt*Wii.to(arcs) + g_tf*R_Wij.in(arcs) - b_tf*Im_Wij.in(arcs)),2);
//        //        SDP.add(PLoss2.in(arcs)==0,true);
//        //        Constraint<> PLoss("PLoss");
//        //        PLoss = pow(Pf_from,2) - pow(Pf_to,2);
//        //        PLoss -= pow((g_ff*Wii.from(arcs) + g_ft*R_Wij.in(arcs) + b_ft*Im_Wij.in(arcs)),2);
//        //        PLoss += pow((g_tt*Wii.to(arcs) + g_tf*R_Wij.in(arcs) - b_tf*Im_Wij.in(arcs)),2);
//        //        SDP.add(PLoss.in(arcs)==0,true);
//
//        //        Constraint<> QLoss("QLoss");
//        //        QLoss = pow(Qf_from,2) - pow(Qf_to,2);
//        //        QLoss -= pow((g_ft*Im_Wij.in(arcs) - b_ff*Wii.from(arcs) - b_ft*R_Wij.in(arcs)),2);
//        //        QLoss += pow(-1*(b_tt*Wii.to(arcs) + b_tf*R_Wij.in(arcs) + g_tf*Im_Wij.in(arcs)),2);
//        //        SDP.add(QLoss.in(arcs)==0,true);
//        //        SDP.print();
//
//        Constraint<Cpx> I_from("I_from");
//        I_from=(Y+Ych)*(conj(Y)+conj(Ych))*Wii.from(arcs)-T*Y*(conj(Y)+conj(Ych))*conj(Wij)-conj(T)*conj(Y)*(Y+Ych)*Wij+pow(tr,2)*Y*conj(Y)*Wii.to(arcs);
//        //       SDP.add_real(I_from.in(arcs)==pow(tr,2)*L_from);
//        //        SDPOA->add_real(I_from.in(arcs)==pow(tr,2)*L_from);
//
//        var<Cpx> L_to("L_to");
//        L_to.set_real(lji.in(arcs));
//
//        Constraint<Cpx> I_to("I_to");
//        I_to=pow(tr,2)*(Y+Ych)*(conj(Y)+conj(Ych))*Wii.to(arcs)-conj(T)*Y*(conj(Y)+conj(Ych))*Wij-T*conj(Y)*(Y+Ych)*conj(Wij)+Y*conj(Y)*Wii.from(arcs);
//        //SDP.add_real(I_to.in(arcs)==pow(tr,2)*L_to);
//
//        Constraint<> I_from_Pf("I_from_Pf");
//        I_from_Pf=lij*Wii.from(arcs)-pow(tr,2)*(pow(Pf_from,2) + pow(Qf_from,2));
//        //        SDP.add(I_from_Pf.in(arcs)==0, true);
//        //        SDPOA->add(I_from_Pf.in(arcs)==0, true);
//
//        Constraint<> I_to_Pf("I_to_Pf");
//        I_to_Pf=lji*Wii.to(arcs)-(pow(Pf_to,2) + pow(Qf_to, 2));
//        //    SDP.add(I_to_Pf.in(arcs)==0, true);
//
//    }
//
//
//
//    total_time_start = get_wall_time();
//    /* Solver selection */
//    solver<> SDPOPF(SDP,solv_type);
//    double solver_time_start = get_wall_time();
//
//    //    SDP.print();
//    SDPOPF.run(output = 5, tol = 1e-6);
//    //    SDP.print_solution();
//    SDP.print_constraints_stats(tol);
//    SDP.print_nonzero_constraints(tol,true);
//    auto lower_bound = SDP.get_obj_val();
//    SDP.print_solution();
//
//
//    //First iteration
//    //Possible improvements: Find interior point via optimization when point is an active point
//    //Best way to do this is to rewrite the model as g(x)<=\eta y and set y (a flag) for each con
//    //wait for Hassan to do this
//    //generate as many OA iterative cuts as given by num_iter_cuts
//
//    vector<string> con_names={"SOC_convex", "Thermal_Limit_from", "Thermal_Limit_to"};
//
//
//    if(!grid._tree && grid.add_3d_nlin && sdp_cuts)
//    {
//        con_names.push_back("SDP_3D");
//    }
//
//    //const string con_names[]={"SOC_convex", "Thermal_Limit_from"};
//
//    bool interior=false;
//    pair<vector<double>,bool> xactive;
//
//    vector<vector<double>> xouter_array, xactive_array;
//        vector<double> xsolution;
//        int counter;
//        double xv;
//    const double active_tol=1e-6;
//
//
//    for (auto &cname: con_names)
//    {
//        auto con=SDP.get_constraint(cname);
//        for(auto i=0;i<con->get_nb_inst();i++)
//            //  for(auto i=0;i<1;i++)
//        {
//            con->uneval();
//            DebugOn("eval of con "<<con->eval(i)<<endl);
//            con->uneval();
//
//            if(std::abs(con->eval(i))<=active_tol)
//            {
//                con->uneval();
//                func<> oacon=con->get_outer_app_insti(i);
//                oacon.eval_all();
//                Constraint<> OA_sol("OA_cuts_solution"+cname+to_string(i));
//                OA_sol=oacon;
//                if(con->_ctype==leq)
//                    SDPOA->add(OA_sol<=0);
//                else if(con->_ctype==geq)
//                    SDPOA->add(OA_sol>=0);
//
//                oacon.uneval();
//
//                OA_sol.print();
//                DebugOn("OA \t" <<oacon.eval(0));
//
//                DebugOn("Active instant "<<i<<endl);
//            }
//            else //If constraint is not active xsolution is an interior point
//            {
//                xsolution.clear();
//                for (auto &it: *con->_vars)
//                {
//                    auto v = it.second.first;
//                    size_t posv=v->get_id_inst(i);
//                    v->get_double_val(posv, xv);
//                    xsolution.push_back(xv);
//                }
//
//                xactive_array= con->get_any_active_point(i,  con->_ctype);
//
//                for(auto j=0;j<xactive_array.size();j++)
//                {
//                    if(xactive_array[j].size()>0)
//                    {
//                        con->uneval();
//
//                        counter=0;
//                        for (auto &it: *con->_vars)
//                        {
//                            auto v = it.second.first;
//                            size_t posv=v->get_id_inst(i);
//                            v->set_double_val(posv,xactive_array[j][counter++]);
//                        }
//                        con->uneval();
//                        func<> oa_iter=con->get_outer_app_insti(i);
//                        oa_iter.eval_all();
//                        Constraint<> OA_itercon("OA_cuts_iterative "+cname+to_string(i)+","+to_string(j));
//                        OA_itercon=oa_iter;
//                        if(con->_ctype==leq)
//                            SDPOA->add(OA_itercon<=0);
//                        else if(con->_ctype==geq)
//                            SDPOA->add(OA_itercon>=0);
//
//                    }
//                }
//
//
//                counter=0;
//                for (auto &it: *con->_vars)
//                {
//                    auto v = it.second.first;
//                    size_t posv=v->get_id_inst(i);
//                    v->set_double_val(posv, xsolution[counter++]);
//                }
//
//
//                
//
//
//            }
//        }
//    }
//    return SDPOA;
//}
//
/** Return the vector of arcs of the chordal completion ignoring parallel lines **/
indices PowerNet::get_node_pairs_chord(const vector<pair<string,vector<Node*>>>& bags){
    if(!this->node_pairs_chord.empty()){
        return this->node_pairs_chord;
    }
    map<string,pair<Node*,Node*>> unique_pairs;
    for (auto a: arcs) {
        if (!a->_parallel) {
            unique_pairs[a->_src->_name+","+a->_dest->_name] = {a->_src,a->_dest};
            node_pairs_chord.add(a->_src->_name+","+a->_dest->_name);
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
    for (auto &bag: bags) {
        for (size_t i = 0; i< bag.second.size()-1; i++) {
            if (unique_pairs.insert({bag.second[i]->_name+","+bag.second[i+1]->_name,{bag.second[i],bag.second[i+1]}}).second) {
                auto bus_s = (Bus*)bag.second[i];
                auto bus_d = (Bus*)bag.second[i+1];
                w_max_ = bus_s->vbound.max*bus_d->vbound.max;
                w_min_ = bus_s->vbound.min*bus_d->vbound.min;
                wr_max_ = cos_max_*w_max_;
                if(cos_min_ < 0) wr_min_ = cos_min_*w_max_;
                else wr_min_ = cos_min_*w_min_;
                if(sin_max_ > 0) wi_max_ = sin_max_*w_max_;
                else wi_max_ = sin_max_*w_min_;
                if(sin_min_ > 0) wi_min_ = sin_min_*w_min_;
                else wi_min_ = sin_min_*w_max_;
                auto name = bag.second[i]->_name + "," + bag.second[i+1]->_name;
                wr_max.add_val(name,wr_max_);
                wr_min.add_val(name,wr_min_);
                wi_max.add_val(name,wi_max_);
                wi_min.add_val(name,wi_min_);
                node_pairs_chord.add(name);
            }
        }
        /* Loop back pair */
        if (unique_pairs.insert({bag.second[0]->_name+","+bag.second[bag.second.size()-1]->_name,{bag.second[0],bag.second[bag.second.size()-1]}}).second) {
            auto name = bag.second[0]->_name + "," + bag.second[bag.second.size()-1]->_name;
            auto bus_s = (Bus*)bag.second[0];
            auto bus_d = (Bus*)bag.second[bag.second.size()-1];
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
            node_pairs_chord.add(name);
        }
    }
    return node_pairs_chord;
}

indices PowerNet::get_node_pairs_chord_bags(std::vector<pair<string,vector<Node*>>> bags){
    if(!this->node_pairs_chord.empty()){
        return this->node_pairs_chord;
    }
    set<pair<Node*,Node*>> unique_pairs;
    for (auto a: arcs) {
        if (!a->_parallel) {
            unique_pairs.insert({a->_src,a->_dest});
            node_pairs_chord.add(a->_src->_name+","+a->_dest->_name);
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
    for (auto &bag: bags) {
        for (size_t i = 0; i< bag.second.size()-1; i++) {
            if (unique_pairs.insert({bag.second[i],bag.second[i+1]}).second) {
                auto bus_s = (Bus*)bag.second[i];
                auto bus_d = (Bus*)bag.second[i+1];
                w_max_ = bus_s->vbound.max*bus_d->vbound.max;
                w_min_ = bus_s->vbound.min*bus_d->vbound.min;
                wr_max_ = cos_max_*w_max_;
                if(cos_min_ < 0) wr_min_ = cos_min_*w_max_;
                else wr_min_ = cos_min_*w_min_;
                if(sin_max_ > 0) wi_max_ = sin_max_*w_max_;
                else wi_max_ = sin_max_*w_min_;
                if(sin_min_ > 0) wi_min_ = sin_min_*w_min_;
                else wi_min_ = sin_min_*w_max_;
                auto name = bag.second[i]->_name + "," + bag.second[i+1]->_name;
                wr_max.add_val(name,wr_max_);
                wr_min.add_val(name,wr_min_);
                wi_max.add_val(name,wi_max_);
                wi_min.add_val(name,wi_min_);
                node_pairs_chord.add(name);
            }
        }
        /* Loop back pair */
        if (unique_pairs.insert({bag.second[0],bag.second[bag.second.size()-1]}).second) {
            auto name = bag.second[0]->_name + "," + bag.second[bag.second.size()-1]->_name;
            auto bus_s = (Bus*)bag.second[0];
            auto bus_d = (Bus*)bag.second[bag.second.size()-1];
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
            node_pairs_chord.add(name);
        }
    }
    return node_pairs_chord;
}


gravity::indices PowerNet:: node_pairs_no_line_charge()
{
    indices node_pairs_nolinecharge = indices("node_pairs_nolinecharge");
    for(auto &bp:arcs)
    {
        if(ch_half.eval(bp->_name)==0.0 && tr.eval(bp->_name)==1)
        {
            node_pairs_nolinecharge.add(bp->_name);
        }
    }
    return(node_pairs_nolinecharge);
}


gravity::indices PowerNet:: arcs_line_charge()
{
    indices arcs_charged = indices("arcs_charged");
    for(auto &bp:arcs)
    {
        if(ch_half.eval(bp->_name)!=0.0)
        {
            arcs_charged.add(bp->_name);
        }
    }
    return(arcs_charged);
}

gravity::indices PowerNet:: arcs_inductive_only()
{
    indices arcs_inductive = indices("arcs_inductive");
    for(auto &bp:arcs)
    {
        if(ch_half.eval(bp->_name)==0.0 && g.eval(bp->_name)==0 && tr.eval(bp->_name)==1 && cc.eval(bp->_name)==1)
        {
            arcs_inductive.add(bp->_name);
        }
    }
    return(arcs_inductive);
}

gravity::indices PowerNet:: arcs_not_inductive_only()
{
    indices arcs_not_inductive_only = indices("arcs_not_inductive");
    for(auto &bp:arcs)
    {
        if(ch_half.eval(bp->_name)>0.0 && g.eval(bp->_name)>0.0)
        {
            arcs_not_inductive_only.add(bp->_name);
        }
    }
    return(arcs_not_inductive_only);
}

pair<pair<indices,indices>,pair<indices,indices>> PowerNet:: get_pairsof_node_pairs_ijkl()
{
    indices ij = indices("ij");
    indices kl = indices("kl");
    indices il = indices("il");
    indices jk = indices("jk");
    indices p1=node_pairs_chord.exclude(node_pairs_chord.last());
    indices p2=node_pairs_chord.exclude(node_pairs_chord.first());
    pair<pair<indices,indices>,pair<indices,indices>> res;
    for(auto &bp:*p1._keys)
    {
        auto n1 = bp.substr(0,bp.find_first_of(","));
        auto n2 = bp.substr(bp.find_first_of(",")+1);
        for(auto &ip:*p2._keys)
        {
            auto n3=ip.substr(0, ip.find_first_of(","));
            auto n4=ip.substr(ip.find_first_of(",")+1);
            
            if(n1!=n3 && n1!=n4 && n2!=n3 && n2!=n4)
            {
                
                auto il_n=n1+","+n4;
                auto jk_n=n2+","+n3;
                
                // DebugOn(bp<<" "<<ip<<" "<<il_n<<" "<<jk_n<<endl);
                if(node_pairs._keys_map->find(il_n)!=node_pairs._keys_map->end() && node_pairs._keys_map->find(jk_n)!=node_pairs._keys_map->end())
                {
                    
                    
                    
                    ij.add(bp);
                    kl.add(ip);
                    jk.add(jk_n);
                    il.add(il_n);
                    DebugOn("innermost");
                    DebugOn(bp<<" "<<ip<<" "<<il_n<<" "<<jk_n<<endl);
                }
                
                
                
            }
        }
        p2.exclude(p2.first());
    }
    
    auto pair1=make_pair(ij,kl);
    auto pair2=make_pair(jk,il);
    res=make_pair(pair1,pair2);
    return res;
    
}
double PowerNet::solve_acopf(PowerModelType pmt, int output, double tol){
    
    auto ACOPF = build_ACOPF(*this,pmt,output,tol);
    bool relax;
    solver<> OPF(ACOPF,ipopt);
    //    auto mipgap = 1e-6;
    OPF.run(output, tol);
    
    
    //    ACOPF->print_solution();
    
    return ACOPF->_obj->get_val();
}

double PowerNet::solve_sdpopf(bool loss_from, int output, double tol){
    bool relax;
    double upper_bound=solve_acopf(ACPOL, output, tol);
    auto SDPOPF = build_SDPOPF(*this, loss_from, upper_bound);
    solver<> OPF(SDPOPF,ipopt);
    //    auto mipgap = 1e-6;
    OPF.run(output, tol);
    return SDPOPF->_obj->get_val();
}


void PowerNet::fill_wbnds(){
    double cos_max_, cos_min_, w_max_, w_min_, wr_max_, wr_min_, sin_max_, sin_min_, wi_max_, wi_min_;
    for(int i = 0; i < nodes.size()-1; i++) {
        for(int j = i+1; j < nodes.size(); j++) {
            Bus *bus_s = (Bus *) (nodes[i]);
            Bus *bus_d = (Bus *) (nodes[j]);
            
            if(get_arc(bus_s, bus_d)) continue;
            
            string name = bus_s->_name + "," + bus_d->_name;
            //            _node_pairs_chord._keys.push_back(new index_pair(index_(bus_s->_name), index_(bus_d->_name)));
            
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


