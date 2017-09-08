//
//  PowerNet.cpp
//
//
//  Created by Guagnlei on 03/06/2017.
//

#include <gravity/PowerNet.h>
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


static int max_line_len;
static char* line = nullptr;

PowerNet::PowerNet(){
    bMVA=0;
}

// Read a grid 
int PowerNet::readgrid(const char* fname) {
    //FILE *fp = fopen(fn,"r");
    //if(fp == NULL)
    //{
    //    fprintf(stderr,"can’t open input file %s\n",fn);
    //    exit(1);
    //}
    //
    //this->Net::readFile(fn);
    //
    //fclose(fp);
    double pl = 0, ql = 0, gs = 0, bs = 0, kvb = 0, vmin = 0, vmax = 0, vs = 0;
    int id = 0;
    cout << "Loading file " << fname << endl;
    ifstream file(fname, std::ifstream::in);
    if(!file.is_open()){
        cout << "Could not open file\n";
        return -1;
    }

    _clone = new Net(); // just for chordal extension
    string word;
    while (word.compare("function")){
        file >> word;
    }
    
    file.ignore(6);
    file >> word;
//    getline(file, word);
    _name = word;
//    cout << _name << endl;
    while (word.compare("mpc.baseMVA")){
        file >> word;
    }
    file.ignore(3);
    getline(file, word,';');
    bMVA = atoi(word.c_str());
//    cout << "BaseMVA = " << bMVA << endl;
    
    /* Nodes data */
    while (word.compare("mpc.bus")){
        file >> word;
    }
    getline(file, word);
    Node* node = NULL;
    Node* node_clone = NULL;
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
        vs = atof(word.c_str());
        file >> ws >> word >> ws >> word;
        kvb = atof(word.c_str());
        file >> ws >> word >> ws >> word;
        vmax = atof(word.c_str());
        getline(file, word,';');
        vmin = atof(word.c_str());
        node = new Node(name, pl, ql, gs, bs, vmin, vmax, kvb, 1);
        node_clone = new Node(name, pl, ql, gs, bs, vmin, vmax, kvb, 1);
        node->vs = vs;
        add_node(node);
        _clone->add_node(node_clone);
//        node->print();
        file >> word;
    }
    file.seekg (0, file.beg);
    /* Generator data */
    while (word.compare("mpc.gen")){
        file >> word;
    }
    double qmin = 0, qmax = 0, pmin = 0, pmax = 0, ps = 0, qs = 0;
    int status = 0;
    getline(file, word);
    file >> word;
    std::vector<bool> gen_status;
    while(word.compare("];")){
        name = word.c_str();
        node = get_node(name);
        file >> word;
        ps = atof(word.c_str())/bMVA;
        file >> word;
        qs = atof(word.c_str())/bMVA;
        file >> word;
        qmax = atof(word.c_str())/bMVA;
        file >> word;
        qmin = atof(word.c_str())/bMVA;
        file >> ws >> word >> ws >> word >> ws >> word;
        status = atof(word.c_str());
        file >> word;
        pmax = atof(word.c_str())/bMVA;
        file >> word;
        pmin = atof(word.c_str())/bMVA;
        getline(file, word,'\n');
        gen_status.push_back(status==1);
        if(status==1){
            node->_has_gen = true;
            Gen* g = new Gen(node, to_string(node->_gen.size()), pmin, pmax, qmin, qmax);
            g->ps = ps;
            g->qs = qs;
            gens.push_back(g);
            node->_gen.push_back(g);
//            g->print();
        }
//        getline(file, word);
        file >> word;
    }
    file.seekg (0, file.beg);
    /* Generator costs */
    while (word.compare("mpc.gencost")){
        file >> word;
    }
    double c0 = 0, c1 = 0,c2 = 0;
    getline(file, word);
//    cout<<"Number of generators = " << gens.size() << endl;
    int gen_counter = 0;
    for (int i = 0; i < gen_status.size(); ++i) {
        file >> ws >> word >> ws >> word >> ws >> word >> ws >> word >> ws >> word;
        c2 = atof(word.c_str());
        file >> word;
        c1 = atof(word.c_str());
        file >> word;
        c0 = atof(word.c_str());
        if (gen_status[i]) {
          gens[gen_counter]->set_costs(c0, c1, c2);
//        gens[gen_counter]->print();
          gen_counter++;
        }
        getline(file, word);
    }
    file.seekg (0, file.beg);
    /* Lines data */
    m_theta_lb = 0;
    m_theta_ub = 0;
    while (word.compare("mpc.branch")){
        file >> word;
    }
    getline(file, word);
    double res = 0;
    Arc* arc = NULL;
    Arc* arc_clone = NULL;
    string src, dest;
    file >> word;
    while(word.compare("];")){
        src = word;
        file >> dest;
        id = (int)arcs.size() + 1;
        //arc = new Arc(src+dest);
        //arc_clone = new Arc(src+dest);
        arc = new Arc(to_string(id));
        arc_clone = new Arc(to_string(id));
        arc->id = id;
        arc_clone->id = id;
        arc->src = get_node(src);
        arc->dest= get_node(dest);
        arc_clone->src = _clone->get_node(src);
        arc_clone->dest = _clone->get_node(dest);
        file >> word;
        arc->r = atof(word.c_str());
        file >> word;
        arc->x = atof(word.c_str());
        res = pow(arc->r,2) + pow(arc->x,2);
        if (res==0) {
            cerr << " line with r = x = 0" << endl;
            exit(-1);
        }
        arc->g = arc->r/res;
        arc->b = -arc->x/res;
        
        file >> word;
        arc->ch = atof(word.c_str());
        file >> word;
        arc->limit = atof(word.c_str())/bMVA;
        file >> ws >> word >> ws >> word >> ws >> word;
        if(atof(word.c_str()) == 0)
            arc->tr = 1.0;
        else
            arc->tr = atof(word.c_str());
        file >> word;
        arc->as = atof(word.c_str())*M_PI/180;
        file >> word;
        arc->cc = arc->tr*cos(arc->as);
        arc->dd = arc->tr*sin(arc->as);
        arc->status = atof(word.c_str());
        arc_clone->status = arc->status;
        file >> word;
        arc->tbound.min = atof(word.c_str())*M_PI/180;
        arc_clone->tbound.min = arc->tbound.min;
        m_theta_lb += arc->tbound.min;
        file >> word;
        arc->tbound.max = atof(word.c_str())*M_PI/180;
        arc_clone->tbound.max = arc->tbound.max;
        m_theta_ub += arc->tbound.max;
        arc->smax = max(pow(arc->src->vbound.max,2)*(arc->g*arc->g+arc->b*arc->b)*(pow(arc->src->vbound.max,2) + pow(arc->dest->vbound.max,2)), pow(arc->dest->vbound.max,2)*(arc->g*arc->g+arc->b*arc->b)*(pow(arc->dest->vbound.max,2) + pow(arc->src->vbound.max,2)));
        if(arc->status==1){
            arc->connect();
            if(!add_arc(arc)){// not a parallel line
                arc_clone->connect();
                _clone->add_arc(arc_clone);
            }
            else {
                delete arc_clone;
            }
        }
        else {
            delete arc_clone;
            delete arc;
        }
//        arc->print();
        getline(file, word,'\n');
        file >> word;
    }
    file.close();
//    for (auto n:nodes) {
//        n->print();
//        cout << "node" << n->ID << ": fill_in = " << n->fill_in << endl;
//    }
    get_tree_decomp_bags();
    return 0;
}
