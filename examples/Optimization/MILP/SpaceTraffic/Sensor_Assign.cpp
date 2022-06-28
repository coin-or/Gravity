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
    int wUB = 100; //UB on weights
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
    }
    
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
    var<int> s("s", 0, 1);
    model.add(s.in(own_arcs));
    var<int> sn("sn", 0, 1);
    model.add(sn.in(sensors));
    var<int> z0("z0");
    model.add(z0.in(arcs));
    var<int> z("z", 0, 1);
    model.add(z.in(bought_arcs));
    var<double> u("u", pos_);
    model.add(u.in(own_sens));
    var<double> up("up");
    model.add(up.in(bought_sens));
    var<double> q("q", pos_);
    model.add(q.in(own_arcs));
    var<double> qp("qp", pos_);
    model.add(qp.in(bought_arcs));
    var<double> r("r", pos_);
    model.add(r.in(jk));

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
    w_own.print();
    w_bought.print();
    jk.print();
    
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
    
    
    /*Constraint<> fub("Follower_Unique_Object");
    fub = sum(s.in_matrix(0, 1)) + sum(z.in_matrix(0, 1));
    model.add(fub.in(jk) <= 1);*/
    
        //----Dual Feasibility----
    
    Constraint<> d1("DualConstr1");
    d1 = u.in_ignore_ith(1, 1, own_arcs) + q + r.in_ignore_ith(0, 1, own_arcs) - w_own;
    model.add(d1.in(own_arcs) >= 0);
    
    Constraint<> d2("DualConstr2");
    d2 = u - p;
    model.add(d2.in(own_sens) >= 0);
    
    Constraint<> d3("DualConstr3");
    d3 = up.in_ignore_ith(1, 1, bought_arcs) + qp + r.in_ignore_ith(0, 1, bought_arcs) - w_bought + p.in_ignore_ith(1, 2, bought_arcs);
    model.add(d3.in(bought_arcs) >= 0);
    
        //----Strong Duality----
    
    indices agents = range(0, K - 1);
    Constraint<> sd("Strong_Duality");
    sd = sum(u.in_matrix(0, 1)) + sum(up.in_matrix(0, 1)) + sum(q.in_matrix(0, 2)) + sum(qp.in_matrix(0, 2)) + sum(r.in_matrix(0, 1)) - sum(product(w_own.in_matrix(0, 2), s.in_matrix(0, 2))) - sum(product(p, s.in_matrix(0, 2))) - sum(product(w_bought - p.in_ignore_ith(1, 2, bought_arcs), z.in_matrix(0, 2)));
    model.add(sd.in(agents) == 0);
    
    model.print();

}
