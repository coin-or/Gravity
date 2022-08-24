//
//  Sensor_Assign2.cpp
//  sensor_assign
//
//  Created by Svetlana Riabova on 8/18/22.
//

#include <iostream>
#include "Sensor_Assign2.hpp"
#include <chrono>
using namespace std::chrono;

int main(int argc, const char * argv[]) {
    
    myModel m = myModel();
    vector<param<double>> par = m.readData(argc, argv, 1, 2);
    auto start = high_resolution_clock::now();
    m.InitBilevel(par[0], par[1], par[2], 0.001);
    //m.GreedyStart(par[0], par[1], par[2]);
    auto stop = high_resolution_clock::now();
    auto duration1 = duration_cast<seconds>(stop - start);
    cout << "Init + greedy time: " << duration1.count() << endl;
    m.mSolve();
    auto stop2 = high_resolution_clock::now();
    auto duration2 = duration_cast<seconds>(stop2 - stop);
    cout << m.N << " " << m.M << " " << m.K << " " << duration1.count() + duration2.count() << endl;

    return 0;
}

vector<param<double>> myModel::readData(int argc, const char * argv[], int n1, int n2){
    N = 4; M = 7; K = 3;
    int degree = 100;

    if(argc>=2){
            string fname = argv[n1];
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
    /*random_device rd; // obtain a random number from hardware
    mt19937 gen(rd()); // seed the generator
    uniform_int_distribution<> distr(0, K-1); // define the range

    for(int i = 0; i < N; i++)
        owner.push_back(distr(gen)); // generate numbers*/
    
    /* Indexing sets */
    //sensors = range (0, N - 1);
    //objects = range (N, N + M - 1);
    arcs.add(graph.arcs);
    for (int i = 0; i < N; i++) {
        sensors.add("sensor" + to_string(i));
    }
    for (int i = N; i < N + M; i++) {
        objects.add("object" + to_string(i));
    }
    
    DebugOn("Graph has " << graph.nodes.size() << " nodes" << endl);
    DebugOn("Graph has " << graph.arcs.size() << " arcs" << endl);
    //Parameters
    vector<param<double>> par (3);
    param<double> w0("w0");
    param<double> w("w");
    //param<double> w_own("w_own");
    //param<double> w_bought("w_bought");
    if (argc >= 3) {
        string fname = argv[n2];
        fstream file;
        file.open(fname);
        //FILE *fp = fopen(fname.c_str(),"r");
        /*if(f == NULL)
        {
                cout << "Can’t open input file " << fname;
                exit(1);
        }*/
        string tmp;
        string tmp1;
        string tmp2;
        file >> tmp1;
        K = stoi(tmp1);
        for (int i = 0; i < N; i++) {
            file >> tmp2;
            owner.push_back(stoi(tmp2));
        }
        for (int k = 0; k < K; k++) {
            agents.add("agent" + to_string(k));
            for (int i = 0; i < N; i++) {
                if (owner[i] == k) {
                    own_sens.add("sensor" + to_string(i) + ",agent" + to_string(k));
                    for (Arc* a: graph.get_node("sensor" + to_string(i))->get_out()) {
                        own_arcs.add("sensor" + to_string(i) + "," + a->_dest->_name +  ",agent" +  to_string(k));
                        agents_arcs.add("sensor" + to_string(i) + "," + a->_dest->_name +  ",agent" +  to_string(k));
                        for (Arc* b: graph.get_node(a->_dest->_name)->get_in()) {
                            if (owner[stoi(b->_src->_name.substr(6, b->_src->_name.find(",")))] == k) {
                                if (stoi(b->_src->_name.substr(6, b->_src->_name.find(","))) != i) {
                                own_rplc.add(a->_src->_name + "," + b->_src->_name + "," + a->_dest->_name +  ",agent" +  to_string(k));
                                }
                            }
                            else {
                                oths_rplc.add("sensor" + to_string(i) + "," + b->_src->_name + "," + a->_dest->_name +  ",agent" +  to_string(k));
                            }
                        }
                    }
                }
                else {
                    bought_sens.add("sensor" + to_string(i) + "," + to_string(k));
                    for (Arc* a: graph.get_node("sensor" + to_string(i))->get_out()) {
                        bought_arcs.add(a->_src->_name + "," + a->_dest->_name +  ",agent" +  to_string(k));
                        agents_arcs.add(a->_src->_name + "," + a->_dest->_name +  ",agent" +  to_string(k));
                    }
                }
            }
        }
        operations = indices(sensors, agents);
        
        w0.in(arcs);
        w.in(agents_arcs);
        for (Arc* a: graph.arcs) {
            for (int k = 0; k < K; k++) {
                file >> tmp;
                string tmpstr = a->_src->_name + "," + a->_dest->_name + ",agent" + to_string(k);
                w(a->_src->_name + "," + a->_dest->_name + ",agent" + to_string(k)) = stod(tmp);
            }
        }
    }
    else {
        w0.initialize_normal(2, 1);
        w.initialize_normal(2, 1);
    }
    par[0] = w0;
    par[1] = w;
    par[2] = w0;
    
    return par;
}

void myModel::InitBilevel(param<double> &w0, param<double> &w, param<double> &w_bought, double eps) {
    
    e = eps;
    
    w.reset_range();
    w_bought.reset_range();
    /*Variables*/
    var<double> p("p", 0, w._range->second);
    model.add(p.in(sensors));
    var<double> y("y", 0, w._range->second);
    model.add(y.in(sensors));
    
    var<int> s("s", 0, 1);
    model.add(s.in(arcs));
    var<int> z("z", 0, 1);
    model.add(z.in(operations));
    
    func<> f = w.in(own_arcs) + w0.in_ignore_ith(2, 1, own_arcs);
    f.eval_all();
    param<> w_own0("w_own0");
    w_own0.in(own_arcs);
    w_own0.copy_vals(f);
    
    func<> f2 = w.in(bought_arcs) + w0.in_ignore_ith(2, 1, bought_arcs);
    f2.eval_all();
    param<> w_bought0("w_bought0");
    w_bought0.in(bought_arcs);
    w_bought0.copy_vals(f2);
    
    /*Objective*/
    func<> obj;
    obj += product(w.in(own_arcs), s.in_ignore_ith(2, 1, own_arcs) * z.in_ignore_ith(1, 1, own_arcs));
    obj += product(w.in(bought_arcs), s.in_ignore_ith(2, 1, bought_arcs) * z.in_ignore_ith(1, 1, bought_arcs));
    obj -= product(1,p.in_ignore_ith(1, 2, own_arcs) * z.in_ignore_ith(1, 1, own_arcs));
    obj -= product(1,p.in_ignore_ith(1, 2, bought_arcs) * z.in_ignore_ith(1, 1, bought_arcs));
    obj -= product(e,p.in_ignore_ith(1, 1, operations) * z);
    obj += sum(p);
    model.max(obj);
    
    /*Constraints*/
    
    Constraint<> ua("Unique_Assignment");
    ua = sum(s.in_matrix(1, 1));
    model.add(ua.in(sensors) == 1);

    Constraint<> uo("Unique_Operation");
    uo = sum(z.in_matrix(1, 1));
    model.add(uo.in(sensors) == 1);
    
//    Constraint<> uobj("Unique_Object");
//    uobj = sum(s.in_ignore_ith(2, 1, agents_arcs).in_matrix(0, 1)) + sum(z.in_ignore_ith(1, 1, agents_arcs).in_matrix(0, 1)) - N;
//    model.add(uobj.in(agents_arcs.ignore_ith(0, 1)) <= 1);
    
    indices uobj_ids("uobj_ids"), s_ids("s_ids"), z_ids("z_ids");
    s_ids = arcs;
    z_ids = operations;
    int row_id = 0;
    for (int j = 0; j < M; j++) {
        for (int k = 0; k < K; k++) {
            uobj_ids.add("Unique_Object:object" + to_string(j) + ",agent" + to_string(k));
            for (Arc* a: graph.get_node("object" + to_string(j))->get_in()) {
                s_ids.add_in_row(row_id, a->_src->_name + ",object" + to_string(j));
                z_ids.add_in_row(row_id, a->_src->_name + ",agent" + to_string(k));
            }
            row_id++;
        }
    }
    
    Constraint<> uobj("Unique_Object");
    uobj = s.in(s_ids) + z.in(z_ids) - N;
    model.add(uobj.in(uobj_ids) <= 1);
    
    Constraint<> fulb("Followers_Utility_lb");
    fulb = w.in(bought_arcs) - p.in_ignore_ith(1, 2, bought_arcs) * z.in_ignore_ith(1, 1, bought_arcs);
    model.add(fulb.in(bought_arcs) >= 0);
    
        //----Fair price----
    Constraint<> fp("FairPrice");
    fp = p.in_ignore_ith(1, 2, bought_arcs) - (w.in(bought_arcs) * s.in_ignore_ith(2, 1, bought_arcs) * z.in_ignore_ith(1, 1, bought_arcs) + y.in_ignore_ith(1, 2, bought_arcs))/2;
    model.add(fp.in(bought_arcs) >= 0);
    
    Constraint<> sl1("Seller lb1");
    sl1 = y.in_ignore_ith(1, 2, agents_arcs) + w * s.in_ignore_ith(2, 1, agents_arcs) * z.in_ignore_ith(1, 1, agents_arcs) - w;
    model.add(sl1.in(agents_arcs) >= 0);
    
    Constraint<> sl2("Seller lb2");
    sl2 = y.in_ignore_ith(1, 3, own_rplc) - (w.in_ignore_ith(1, 1, own_rplc) - w.in_ignore_ith(0, 1, own_rplc)) * s.in_ignore_ith(0, 1, own_rplc);
    model.add(sl2.in(own_rplc) >= 0);
    
    indices tmp_ids = oths_rplc.ignore_ith(1, 1).ignore_ith(2, 1);
    indices tmp_ids1 = oths_rplc.ignore_ith(0, 1).ignore_ith(2, 1);
    Constraint<> sl3("Seller lb3");
    sl3 = y.in_ignore_ith(1, 3, oths_rplc) - (w.in_ignore_ith(1, 1, oths_rplc) - w.in_ignore_ith(0, 1, oths_rplc)) * s.in(tmp_ids) * z.in_ignore_ith(1, 2, oths_rplc) - p.in(tmp_ids1) * z.in_ignore_ith(1, 2, oths_rplc);
    model.add(sl3.in(oths_rplc) >= 0);
    
    //For comparison
    /*Constraint<> no_colab("nc");
    no_colab = sn;
    model.add(no_colab.in(sensors) == 0);*/
    
//   model.print_symbolic();
    model.print();
//    model.print_constraints_stats(1e-4);
    
}

void myModel::mSolve() {
    solver<> sol(model, gurobi);
    sol.run();
    //model.print_solution();
}
