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
//#define USEDEBUG
#ifdef USEDEBUG
#define Debug(x) cout << x
#else
#define Debug(x)
#endif
#define DebugOn(x) cout << x
#define DebugOff(x)

using namespace std;

PowerNet::PowerNet() {
    bMVA = 0;
    pg_min.set_name("pg_min");
    pg_max.set_name("pg_max");
    qg_min.set_name("qg_min");
    qg_max.set_name("qg_max");
    pg_s.set_name("pg_s");
    qg_s.set_name("qg_s");
    c0.set_name("c0");
    c1.set_name("c1");
    c2.set_name("c2");
    th_min.set_name("th_min");
    th_max.set_name("th_max");
    tan_th_min.set_name("tan_th_min");
    tan_th_max.set_name("tan_th_max");
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
    ql.set_name("ql");
    gs.set_name("gs");
    bs.set_name("bs");
    g.set_name("g");
    b.set_name("b");
    ch.set_name("ch");
    as.set_name("as");
    tr.set_name("tr");
    S_max.set_name("S_max");

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
        delete (Bus*)n;
    }
    nodes.clear();
    for (Arc* a:arcs) {
        delete (Line*)a;
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
            DebugOn("Inactive Node" << n->_name << endl);
        }
    }
    return nb;
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
    tbound_max_tan.time_expand(T);
    tbound_min_tan.time_expand(T);
    wr_min.time_expand(T);
    wi_min.time_expand(T);
    wr_max.time_expand(T);
    wi_max.time_expand(T);
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
        throw invalid_argument("Could not open file\n");
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

    Line* arc = NULL;
    string src,dest,key;
    file >> word;
    index = 0;
    while(word.compare("];")) {
        src = word;
        file >> dest;
        key = dest+","+src;//Taking care of reversed direction arcs
        if(arcID.find(key)!=arcID.end()) {//Reverse arc direction
            DebugOn("Adding arc linking " +src+" and "+dest);
            DebugOn(" with reversed direction, reversing source and destination.\n");
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


        arc->cc = arc->tr*cos(arc->as); // Rectangular values for transformer phase shifters
        arc->dd = arc->tr*sin(arc->as);
        arc->status = atoi(word.c_str());
        file >> ws >> word;

        arc->tbound.min = atof(word.c_str())*pi/180.;
//        arc->tbound.min = -30*pi/180;
        m_theta_lb += arc->tbound.min;
        file >>  ws >>word;

        arc->tbound.max = atof(word.c_str())*pi/180.;
        if (arc->tbound.min==0 && arc->tbound.max) {
            DebugOn("Angle bounds are equal to zero. Setting them to -+60");
             arc->tbound.min = -60*pi/180;
            arc->tbound.max = 60*pi/180;
            
        }
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
            DebugOn("INACTIVE ARC!\n" << arc->_name << endl);
        }
        arc->connect();
        add_arc(arc);

        /* Switching to bus_pairs keys */
        name = bus_s->_name + "," + bus_d->_name;
        if (!arc->_parallel) {
            th_min.set_val(name,arc->tbound.min);
            th_max.set_val(name,arc->tbound.max);
            tan_th_min.set_val(name,tan(arc->tbound.min));
            tan_th_max.set_val(name,tan(arc->tbound.max));
            _bus_pairs._keys.push_back(new index_pair(index_(bus_s->_name), index_(bus_d->_name), arc->_active));
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
        getline(file, word,'\n');
        file >> word;
    }
    DebugOff(ch.to_str(true) << endl);
    DebugOff(as.to_str(true) << endl);
    DebugOff(tr.to_str(true) << endl);

    file.close();
    return 0;
}

