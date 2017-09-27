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

PowerNet::PowerNet(){
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
    v_min.set_name("v_min"); 
    v_max.set_name("v_max"); 
    v_s.set_name("v_s"); 
    pl.set_name("pl"); 
    ql.set_name("ql"); 
    gs.set_name("gs"); 
    bs.set_name("bs"); 
    g.set_name("g"); 
    b.set_name("b"); 
    ch.set_name("ch"); 
    S_max.set_name("S_max"); 

    g_ff.set_name("g_ff"); 
    g_ft.set_name("g_ft"); 
    g_tf.set_name("g_tf"); 
    g_tt.set_name("g_tt"); 

    b_ff.set_name("b_ff"); 
    b_ft.set_name("b_ft"); 
    b_tf.set_name("b_tf"); 
    b_tt.set_name("b_tt"); 
}

PowerNet::~PowerNet(){
    if(!gens.empty()){
        for (Gen* g:gens){
            delete g;
        }
        gens.clear();
    }
}
//
// Read a grid
// @discussion line with delimiter ";"

int PowerNet::readgrid(const char* fname) {
    
    string name;
    double kvb = 0;
    int id = 0;
    unsigned index = 0;
    cout << "Loading file " << fname << endl;
    ifstream file(fname, std::ifstream::in);
    if(!file.is_open()){
        cout << "Could not open file\n";
        return -1;
    }
    
    string word;
    while (word.compare("function")){
        file >> word;
    }
    
    file.ignore(6);
    file >> word;
    _name = word;
    
//  cout << _name << endl;
    while (word.compare("mpc.baseMVA")){
        file >> word;
    }
    
    file.ignore(3);
    getline(file, word,';');
    bMVA = atoi(word.c_str());
//  cout << "BaseMVA = " << bMVA << endl;
    
    /* Nodes data */
    while (word.compare("mpc.bus")){
        file >> word;
    }
    
    getline(file, word);
    Bus* bus = NULL;
    Bus* bus_clone= NULL;
    file >> word;
    
    while(word.compare("];")){
        name = word.c_str();
        id = atoi(name.c_str());
        file >> ws >> word >> ws >> word;
        pl = atof(word.c_str())/bMVA;
        file >> word;
        ql = atof(word.c_str())/bMVA;
        file >> word;
        gs = atof(word.c_str())/bMVA;
        file >> word;
        bs = atof(word.c_str())/bMVA;
        file >> ws >> word >> ws >> word;
        v_s = atof(word.c_str());
        file >> ws >> word >> ws >> word;
        kvb = atof(word.c_str());
        file >> ws >> word >> ws >> word;
        v_max = atof(word.c_str());
        getline(file, word,';');
        v_min = atof(word.c_str());
        // single phase
        bus = new Bus(name, pl.eval(), ql.eval(), gs.eval(), bs.eval(), v_min.eval(), v_max.eval(), kvb, 1);
        bus_clone = new Bus(name, pl.eval(), ql.eval(), gs.eval(), bs.eval(), v_min.eval(), v_max.eval(), kvb, 1);
        bus->vs = v_s.eval();
        bus_clone->vs = v_s.eval();
        
        this->Net::add_node(bus);
        file >> word;
    }
    file.seekg (0, file.beg);
    
    
    /* Generator data */
    while (word.compare("mpc.gen")){
        file >> word;
    }
//    double qmin = 0, qmax = 0, pmin = 0, pmax = 0, ps = 0, qs = 0;
    int status = 0;
    getline(file, word);
    
    
    file >> word;
    std::vector<bool> gen_status;
    
    while(word.compare("];")){
        name = word.c_str();
        // name -> node. 
        bus = (Bus*)(Net::get_node(name));
        file >> word;
        pg_s = atof(word.c_str())/bMVA;
        file >> word;
        qg_s = atof(word.c_str())/bMVA;
        file >> word;
        qg_max = atof(word.c_str())/bMVA;
        file >> word;
        qg_min = atof(word.c_str())/bMVA;
        
        file >> ws >> word >> ws >> word >> ws >> word;
        status = atof(word.c_str());
        file >> word;
        pg_max = atof(word.c_str())/bMVA;
        
        file >> word;
        pg_min = atof(word.c_str())/bMVA;
        
        getline(file, word,'\n');
        gen_status.push_back(status==1);
        
        if(status==1){
            bus->_has_gen = true;
            /** generator name, ID */
            Gen* g = new Gen(bus, to_string(gens.size()+1), pg_min.eval(), pg_max.eval(), qg_min.eval(), qg_max.eval());
            g->ID = gens.size();
            g->ps = pg_s.eval();
            g->qs = qg_s.eval();
            gens.push_back(g);
            bus->_gen.push_back(g);
        }
//        getline(file, word);
        file >> word;
    }
    
    
    file.seekg (0, file.beg);
    
    /* Generator costs */
    while (word.compare("mpc.gencost")){
        file >> word;
    }
//    double c0 = 0, c1 = 0,c2 = 0;
    getline(file, word);
    
    int gen_counter = 0;
    for (int i = 0; i < gen_status.size(); ++i) {
        file >> ws >> word >> ws >> word >> ws >> word >> ws >> word >> ws >> word;
        c2 = atof(word.c_str())*powf(bMVA,2);
        file >> word;
        c1 = atof(word.c_str())*bMVA;
        file >> word;
        c0 = atof(word.c_str());
        if (gen_status[i]) {
          gens[gen_counter]->set_costs(c0.eval(), c1.eval(), c2.eval());
          gen_counter++;
        }
        getline(file, word);
    }
    file.seekg (0, file.beg);
    
    /* Lines data */
    while (word.compare("mpc.branch")){
        file >> word;
    }
    getline(file, word);
    double res = 0;
    
    Line* arc = NULL;
    string src,dest;
    file >> word;
    while(word.compare("];")){
        src = word;
        file >> dest;
        id = (int)arcs.size();
        
        arc = new Line("(" + src + "," + dest + ")");
        arc->id = id;
        arc->src = get_node(src);
        arc->dest= get_node(dest);
        
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
        file >> word;
        arc->as = atof(word.c_str())*M_PI/180; 
        file >> word;
        
        arc->cc = arc->tr*cos(arc->as); // Rectangular values for transformer phase shifters
        arc->dd = arc->tr*sin(arc->as);
        arc->status = atof(word.c_str());
        file >> word;
        
        arc->tbound.min = atof(word.c_str())*M_PI/180;
        m_theta_lb += arc->tbound.min;
        file >> word;
        
        arc->tbound.max = atof(word.c_str())*M_PI/180;
        m_theta_ub += arc->tbound.max;
        
        Bus* bus_s = (Bus*)(arc->src);
        Bus* bus_d = (Bus*)(arc->dest);

        arc->smax = max(
                        pow(bus_s->vbound.max,2)*(arc->g*arc->g + arc->b*arc->b)*(pow(bus_s->vbound.max,2) + pow(bus_d->vbound.max,2)),
                        pow(bus_d->vbound.max,2)*(arc->g*arc->g+arc->b*arc->b)*(pow(bus_d->vbound.max,2) + pow(bus_s->vbound.max,2))
                        );
        
        g = arc->g;
        b = arc->b;
        g_ff = arc->g/pow(arc->tr, 2.);
        g_ft = (-arc->g*arc->cc + arc->b*arc->dd)/(pow(arc->cc, 2) + pow(arc->dd, 2));
        b_ft = (-arc->b*arc->cc - arc->g*arc->dd)/(pow(arc->cc, 2) + pow(arc->dd, 2));
        
        g_tt = arc->g;
        g_tf = (-arc->g*arc->cc - arc->b*arc->dd)/(pow(arc->cc, 2) + pow(arc->dd, 2));
        b_tf = (-arc->b*arc->cc + arc->g*arc->dd)/(pow(arc->cc, 2) + pow(arc->dd, 2));
        
        b_ff = (arc->ch/2 + arc->b)/pow(arc->tr, 2.);
        b_tt = (arc->ch/2 + arc->b);
        ch = arc->ch;
        S_max = arc->limit;
        th_min = arc->tbound.min;
        th_max = arc->tbound.max;
        
        if(arc->status == 1){
            arc->connect();
            add_arc(arc);
        }
        else {
            delete arc;
        }

        getline(file, word,'\n');
        file >> word;
    }
    file.close();
    return 0;
}

//template<typename type> var_gen<type>::var_gen(const string& name):var<type>(name){
//};
