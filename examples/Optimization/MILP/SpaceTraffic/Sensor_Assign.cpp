//
//  main.cpp
//  bilevel_sensor
//
//  Created by Svetlana Riabova on 6/8/22.
//

#include <iostream>
#include "Sensor_Assign.hpp"

int main(int argc, const char * argv[]) {
    myModel m = myModel();
    m.readData(argc, argv);
    m.InitBilevel();
    //int s = 0;
    return 0;
}
//Agent::Agent(int name) { n = name; }

void myModel::readData(int argc, const char * argv[]){
    N = 4; M = 7; K = 2;
    int degree = 100;

    if(argc>=2){
            string fname = argv[1];
            auto dims = graph.read_pairwise_list(fname);
            N = dims.first;
            M = dims.second;
        }
        else {
            graph.generate_bipartite_random(N,M,degree);
        }

    /* Graph nodes are indexed in {0,...,n+m-1})*/
    assert(N + M == graph.nodes.size());/* Make sure we have the right number of nodes */
    
    //init ownr...
    random_device rd; // obtain a random number from hardware
    mt19937 gen(rd()); // seed the generator
    uniform_int_distribution<> distr(0, K-1); // define the range

    for(int i = 0; i < N; i++)
        owner.push_back(distr(gen)); // generate numbers

    /* Init agents and weights */
    /*int wUB = 100; //UB on weights
    uniform_int_distribution<> distr1(0, wUB); // define the range
    agents.resize(K);
    for (int k = 0; k < K; k++) {
        agents[k].n = k;
        for (auto a: graph.arcs) {
            agents[k].w.push_back(distr1(gen));
        }
    }

    for (int i = 0; i < N; i++) {
        agents[owner[i]].own.push_back(i);
        for (int k = 0; k < owner[i]; k++) {
            agents[k].oths.push_back(i);
        }
        for (int k = owner[i] + 1; k < K; k++) {
            agents[k].oths.push_back(i);
        }
    }
    
    for (int k = 0; k < K; k++) {
        for (int n : agents[k].own) {
            indices tmp(graph.get_node(to_string(n))->get_out());
            agents[k].own_arcs = tmp;
        }
        for (int n : agents[k].oths) {
            indices tmp(graph.get_node(to_string(n))->get_out());
            agents[k].oths_arcs = tmp;
        }
    }*/
    
    /* Indexing sets */
    sensors = range (0, N - 1);
    objects = range (N, N + M - 1);
    arcs.add(graph.arcs);
    /*for (int i = 0; i < N; i++) {
        own_sens.add(to_string(i) + "," + to_string(owner[i]));
        for (int k = 0; k < owner[i]; k++) {
            bought_sens.add(to_string(i) + "," + to_string(k));
        }
        for (int k = owner[i] + 1; k < K; k++) {
            bought_sens.add(to_string(i) + "," + to_string(k));
        }
        for (Arc* a : graph.get_node(to_string(i))->get_out()) {
            own_arcs.add(a->_src->_name + "," + a->_dest->_name +  "," +  to_string(owner[i]));
            for (int k = 0; k < owner[i]; k++) {
                bought_arcs.add(a->_src->_name + "," + a->_dest->_name + "," + to_string(k));
            }
            for (int k = owner[i] + 1; k < K; k++) {
                bought_arcs.add(a->_src->_name + "," + a->_dest->_name + "," + to_string(k));
            }
        }
    }*/
    
    /*for (int k = 0; k < K; k++) {
        for (int i = 0; i < N; i++) {
            if (owner[i] == k) {
                own_sens.add(to_string(i) + "," + to_string(k));
                for (Arc* a : graph.get_node(to_string(i))->get_out()) {
                    own_arcs.add(a->_src->_name + "," + a->_dest->_name +  "," +  to_string(k));
                }
            }
            else {
                bought_sens.add(to_string(i) + "," + to_string(k));
                for (Arc* a : graph.get_node(to_string(i))->get_out()) {
                    bought_arcs.add(a->_src->_name + "," + a->_dest->_name +  "," +  to_string(k));
                }
            }
        }
    }*/
    
    for (int k = 0; k < K; k++) {
        for (int i = 0; i < N; i++) {
            if (owner[i] == k) {
                own_sens.add(to_string(i) + "," + to_string(k));
            }
            else {
                bought_sens.add(to_string(i) + "," + to_string(k));
            }
        }
        for (int j = 0; j < M; j++) {
            for (Arc* a: graph.get_node(to_string(j))->get_in()) {
                if (owner[stoi(a->_src->_name)] == k) {
                    own_arcs.add(a->_src->_name + "," + a->_dest->_name +  "," +  to_string(k));
                    for (Arc* b: graph.get_node(to_string(j))->get_in()) {
                        if ((owner[stoi(b->_src->_name)] == k) && (stoi(b->_src->_name) != stoi(a->_src->_name))) {
                            own_rplc.add(a->_src->_name + "," + b->_src->_name + "," + a->_dest->_name +  "," +  to_string(k));
                            /*for (Arc* c: graph.get_node(to_string(j))->get_in()) {
                                if (owner[stoi(c->_src->_name)] != k) {
                                    own_oths_rplc1.add(a->_src->_name + "," + b->_src->_name + "," + a->_dest->_name +  "," +  to_string(k));
                                    own_oths_rplc2.add(a->_src->_name + "," + c->_src->_name + "," + a->_dest->_name +  "," +  to_string(k));
                                }
                            }*/
                        }
                        else if (stoi(b->_src->_name) != stoi(a->_src->_name)) {
                            oths_rplc.add(a->_src->_name + "," + b->_src->_name + "," + a->_dest->_name +  "," +  to_string(k));
                        }
                    }
                }
                else {
                    bought_arcs.add(a->_src->_name + "," + a->_dest->_name +  "," +  to_string(k));
                }
            }
        }
    }
    
    jk = indices(objects, range(0,K-1));
    

}

void myModel::InitPrimal() {

}
void myModel::InitDual() {

}
void myModel::InitBilevel() {

    Model<> model("BilevelSensor");

    /*Variables*/

    var<double> p("p", pos_);
    model.add(p.in(sensors));
    var<double> y("y", pos_);
    model.add(y.in(sensors));
    var<int> s("s", 0, 1);
    model.add(s.in(own_arcs));
    var<int> sn("sn", 0, 1);
    model.add(sn.in(sensors));
    var<int> z0("z0");
    model.add(z0.in(arcs));
    var<int> z("z", 0, 1);
    model.add(z.in(bought_arcs));

    /*Objective*/
    param<double> w_own("w_own");
    w_own.in(own_arcs);
    w_own.initialize_normal(2, 1);
    param<double> w_bought("w_bought");
    w_bought.in(bought_arcs);
    w_bought.initialize_normal(2, 1);
    param<double> w0("w0");
    w0.in(arcs);
    w0.initialize_normal(2, 1);
    
    func<> obj;
    obj += product(w_own + w0.in_ignore_ith(2, 1, own_arcs), s);
    obj += product(p, sn);
    obj += product((w_bought - p), z);
    obj += product(w0.in_ignore_ith(2, 1, bought_arcs), z);
    
    model.max(obj);
    
    /*Constraints*/
    //Upper level
    
    Constraint<> ub("Unique_Bought_Obsrvn");
    ub = sum(z.in_matrix(1, 2)) - sn;
    model.add(ub.in(sensors) == 0);

    Constraint<> luo("Leader_Unique_Object");
    luo = sum(z0.in_matrix(0, 1));
    model.add(luo.in(objects) <= 1);
    
    //Lower level
        //----Primal Feasibility----

    Constraint<> fua("Unique_Own_Assignment");
    fua = sum(s.in_matrix(1, 1)) + sn;
    model.add(fua.in(own_sens) == 1);
    
    Constraint<> fuab("Unique_Bought_Assignment");
    fuab = sum(z.in_matrix(1, 1));
    model.add(fuab.in(bought_sens) <= 1);
    
    indices z_ids = z.get_matrix_ids(0, 1);
    Constraint<> fub("Follower_Unique_Object");
    fub = sum(s.sum_over(z_ids, 0)) + sum(z.in(z_ids));
    model.add(fub.in(z_ids.ignore_ith(0, 1)) <= 1);
    
        //----Fair price----
    Constraint<> fp("FairPrice");
    fp = p.in_ignore_ith(1, 2, bought_arcs) - (product(w_bought, z) + y.in_ignore_ith(1, 2, bought_arcs))/2;
    model.add(fp.in(bought_arcs) >= 0);
    fp.print();
    /*indices s_ids = s.get_matrix_ids(0, 1);
    s_ids.print();
    Constraint<> sl1("Seller lb1");
    sl1 = y.in_ignore_ith(1, 2, s_ids) - w_own * (1 - sum(s.in(s_ids)) - sum(z.sum_over(s_ids, 0)));
    model.add(sl1.in(s_ids) >= 0);*/
    
    Constraint<> sl2("Seller lb2");
    sl2 = y.in_ignore_ith(1, 3, own_rplc) - (w_own.in_ignore_ith(1, 1, own_rplc) - w_own.in_ignore_ith(0, 1, own_rplc)) * s.in_ignore_ith(0, 1, own_rplc);
    model.add(sl2.in(own_rplc) >= 0);
    
    Constraint<> sl3("Seller lb3");
    sl3 = y.in_ignore_ith(1, 3, oths_rplc) - (w_own.in_ignore_ith(1, 1, oths_rplc) - (1 - p.in_ignore_ith(1, 3, oths_rplc))) * z.in_ignore_ith(0, 1, oths_rplc);
    model.add(sl3.in(oths_rplc) >= 0);
    
    model.print();

}
