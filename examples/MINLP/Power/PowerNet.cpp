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
    cphi.set_name("cphi");
    sphi.set_name("sphi");
    cos_d.set_name("cos_d");
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
        if (arc->tbound.min==0 && arc->tbound.max==0) {
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
        DebugOn("\nNew iteration");
        for(auto b_it = _bags.begin(); b_it != _bags.end();) {
            std::vector<Node*> b = *b_it;
            if(b.size() == 3) {
                cout << "\nBag: " << b[0]->_name << ", " << b[1]->_name << ", " << b[2]->_name;
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
                    DebugOn("\nFixing arc a12 (" << a12->_src->_name << ", " << a12->_dest->_name << "), adding bag #" << id_sorted);
                    fixed++;
                }
                if (a13->_free) {
                    a13->_free = false;
                    DebugOn("\nFixing arc a13 (" << a13->_src->_name << ", " << a13->_dest->_name << "), adding bag #" << id_sorted);
                    fixed++;
                }
                if (a32->_free) {
                    a32->_free = false;
                    DebugOn("\nFixing arc a32 (" << a32->_src->_name << ", " << a32->_dest->_name << "), adding bag #" << id_sorted);
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
                                bags_sorted.push_back(bag);
                                id_sorted++;
                                DebugOn("\nFixing arc in a larger bag (" << a->_src->_name << ", " << a->_dest->_name << ")");
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
            if(b.size() > 2) bags_sorted.push_back(b);
            _bags.erase(b_it);
//            id_sorted++;
    }
    _bags = bags_sorted;

    for(auto& a: arcs) {
        if(a->_imaginary) a->_free = true;
    }


    for(auto& k: _bus_pairs._keys){
        _bus_pairs_chord._keys.push_back(k);
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
                sin_max_ = sin(m_theta_lb);
                sin_min_ = sin(m_theta_ub);
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
}

