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
    m.InitBilevel(par[0], par[1], 0.001);
    m.readGreedySol("/Users/svetlanariabova/Projects/Sensor/Data/sol_tmp/sol1000.dat"); //comment if not using greedy sol
    auto stop = high_resolution_clock::now();
    auto duration1 = duration_cast<seconds>(stop - start);
    cout << "Init + greedy time: " << duration1.count() << endl;
    m.mSolve();
    auto stop2 = high_resolution_clock::now();
    auto duration2 = duration_cast<seconds>(stop2 - stop);
    cout << m.N << " " << m.M << " " << m.K << " " << duration1.count() + duration2.count() << endl; //prints num sensors; num objetcs; num agents; total time after reading input (init + greedy + solver)

    return 0;
}

vector<param<double>> myModel::readData(int argc, const char * argv[], int n1, int n2){

    if(argc>=2){
        /*read graph from file*/
        string fname = argv[n1];
        auto dims = graph.read_pairwise_list(fname);
        N = dims.first;
        M = dims.second;
    }
    else {
        /*generate graph*/
        N = 4; M = 7; K = 3;
        int degree = 100;
        graph.generate_bipartite_random(N,M,degree);
    }

    /* Graph nodes are indexed in {0,...,n+m-1})*/
    assert(N + M == graph.nodes.size());/* Make sure we have the right number of nodes */
    
    /* Indexing sets */
    arcs.add(graph.arcs);
    for (int i = 0; i < N; i++) {
        sensors.add("sensor" + to_string(i));
    }
    for (int i = N; i < N + M; i++) {
        objects.add("object" + to_string(i));
    }
    own_arcs = indices("own_arcs");
    bought_arcs = indices("bought_arcs");
    
    DebugOn("Graph has " << graph.nodes.size() << " nodes" << endl);
    DebugOn("Graph has " << graph.arcs.size() << " arcs" << endl);
    
    //Parameters
    vector<param<double>> par (3);
    param<double> w0("w0");
    param<double> w("w");
    
    if (argc >= 3) {
        string fname = argv[n2];
        fstream file;
        file.open(fname);
        string tmp;
        string tmp1;
        string tmp2;
        
        /*num agents; owner*/
        file >> tmp1;
        K = stoi(tmp1);
        for (int i = 0; i < N; i++) {
            file >> tmp2;
            owner.push_back(stoi(tmp2));
        }
        
        /*define index sets for weights, then read weights (need owner to sepsrate own and bought weights)*/
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
                    bought_sens.add("sensor" + to_string(i) + ",agent" + to_string(k));
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
        file.close();
    }
    else {
        /*generate data; probably won't work because index sets are define inside the reading from file case*/
        random_device rd; // obtain a random number from hardware
        mt19937 gen(rd()); // seed the generator
        uniform_int_distribution<> distr(0, K-1); // define the range
        for(int i = 0; i < N; i++)
            owner.push_back(distr(gen)); // generate numbers
        w0.initialize_normal(2, 1);
        w.initialize_normal(2, 1);
    }
    par[0] = w0;
    par[1] = w;
    
    return par;
}

void myModel::InitBilevel(param<double> &w0, param<double> &w, double eps) {
    
    e = eps; //for e * p_sn in obj
    
    /*reset_range to define ub on prices*/
    w.reset_range();

    /*---------Variables---------*/
    /*price + aux for price*/
    var<double> p("p", 0, w._range->second);
    model.add(p.in(sensors));
    var<double> y("y", 0, w._range->second);
    model.add(y.in(sensors));
    
    /*assignment*/
    var<int> s("s", 0, 1);
    model.add(s.in(arcs));
    var<int> z("z", 0, 1);
    model.add(z.in(operations));
    
    /*weights for obj (own + leader, bought + leader); leader receives utility when an object is observed, no matter by whom*/
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
    
    /*--------Objective---------*/
    func<> obj;
    obj += product(w_own0, s.in_ignore_ith(2, 1, own_arcs) * z.in_ignore_ith(1, 1, own_arcs)); //use own sens
    obj += sum(p); //sell sens pt.1
    obj -= product(1, p.in_ignore_ith(1, 2, own_arcs) * z.in_ignore_ith(1, 1, own_arcs)); //sell sens pt.2
    obj += product(w_bought0, s.in_ignore_ith(2, 1, bought_arcs) * z.in_ignore_ith(1, 1, bought_arcs)); //buy sens pt.1
    obj -= product(1, p.in_ignore_ith(1, 1, bought_sens) * z.in(bought_sens)); //buy sens pt.2
    obj -= product(e, p.in_ignore_ith(1, 1, bought_sens) * z.in(bought_sens)); //regularization term; sets prices to their lb from Fair_price constraints (can't use equality there)
    model.max(obj);
    
    /*--------Constraints--------*/
    
    Constraint<> ua("Unique_Assignment");
    ua = sum(s.in_matrix(1, 1)); //sensor must observe exactly one object
    model.add(ua.in(sensors) == 1);

    Constraint<> uo("Unique_Operation");
    uo = sum(z.in_matrix(1, 1)); //sensor operated by exactly one agent
    model.add(uo.in(sensors) == 1);
    
//    Constraint<> uobj("Unique_Object");
//    uobj = sum(s.in_ignore_ith(2, 1, agents_arcs).in_matrix(0, 1)) + sum(z.in_ignore_ith(1, 1, agents_arcs).in_matrix(0, 1)) - N;
//    model.add(uobj.in(agents_arcs.ignore_ith(0, 1)) <= 1);
    
    /*matching index sets for s and z in Unique_Object*/
    indices uobj_ids("uobj_ids"), s_ids("s_ids"), z_ids("z_ids");
    s_ids = arcs;
    z_ids = operations;
    int row_id = 0;
    for (int j = 0; j < M; j++) {
        for (int k = 0; k < K; k++) {
            uobj_ids.add("object" + to_string(j) + ",agent" + to_string(k));
            for (Arc* a: graph.get_node("object" + to_string(j))->get_in()) {
                s_ids.add_in_row(row_id, a->_src->_name + ",object" + to_string(j));
                z_ids.add_in_row(row_id, a->_src->_name + ",agent" + to_string(k));
            }
            row_id++;
        }
    }
    
    Constraint<> uobj("Unique_Object");
    uobj = s.in(s_ids) + z.in(z_ids) - N; //agent observes each object no more than once (sum_i(s(i, j) * z(i, k)) <= 1 or sum_i(s(i, j) + z(i, k) - 1) <= 1)
    model.add(uobj.in(uobj_ids) <= 1);
    
    Constraint<> fulb("Followers_Utility_lb"); //agents don't pay more than they get
    fulb = w.in(bought_arcs) * s.in_ignore_ith(2, 1, bought_arcs) * z.in_ignore_ith(1, 1, bought_arcs) + w._range->second * (1 - s.in_ignore_ith(2, 1, bought_arcs) * z.in_ignore_ith(1, 1, bought_arcs)) - p.in_ignore_ith(1, 2, bought_arcs) * z.in_ignore_ith(1, 1, bought_arcs);
    model.add(fulb.in(bought_arcs) >= 0);
    
        //----Fair price----
    Constraint<> fp("FairPrice"); //set price to the mid-point btw buyer and seller
    fp = p.in_ignore_ith(1, 2, bought_arcs) - (w.in(bought_arcs) * s.in_ignore_ith(2, 1, bought_arcs) * z.in_ignore_ith(1, 1, bought_arcs) + y.in_ignore_ith(1, 2, bought_arcs))/2;
    model.add(fp.in(bought_arcs) >= 0);
    
    /*Price proxies (marginal utility of using own sensor)*/
    indices c_lb1("c_lb1"), y_lb1("y_lb1"), w_own_lb1("w_own_lb1"), z_lb1("z_lb1"), s_lb1("s_lb1");
    y_lb1 = sensors;
    z_lb1 = bought_arcs;
    s_lb1 = own_arcs;
    w_own_lb1 = own_arcs;
    row_id = 0;
    for (int i = 0; i < N; i++) {
        for (Arc* b: graph.get_node("sensor" + to_string(i))->get_out()) {
            string j = b->_dest->_name;
            c_lb1.add("sensor" + to_string(i) + "," + b->_dest->_name);
            y_lb1.add_ref("sensor" + to_string(i));
            w_own_lb1.add_ref("sensor" + to_string(i) + "," + b->_dest->_name + ",agent" + to_string(owner[i]));
            for (Arc* a: graph.get_node(j)->get_in()) {
                s_lb1.add_in_row(row_id, a->_src->_name + "," + j);
                z_lb1.add_in_row(row_id, a->_src->_name + ",agent" + to_string(owner[i]));
            }
            row_id++;
        }
    }
    
    Constraint<> sl1("Seller lb1"); //utility of adding this observation (if possible)
    sl1 = y.in(y_lb1) + w.in(w_own_lb1) * s.in(s_lb1) * z.in(z_lb1) - w.in(w_own_lb1);
    model.add(sl1.in(c_lb1) >= 0);
    
    /*Probably better to replace this with .add_in_row()*/
    indices tmp_ids2 = own_rplc.ignore_ith(0, 1).ignore_ith(3, 1);
    Constraint<> sl2("Seller lb2"); //utility of replacing observation done by another sensor of their own by this one
    sl2 = y.in_ignore_ith(1, 3, own_rplc) - (w.in_ignore_ith(1, 1, own_rplc) - w.in_ignore_ith(0, 1, own_rplc)) * s.in_ignore_ith(2, 3, tmp_ids2);
    model.add(sl2.in(own_rplc) >= 0);

    indices tmp_ids = oths_rplc.ignore_ith(0, 1).ignore_ith(2, 1);
    indices tmp_ids1 = oths_rplc.ignore_ith(0, 1);
    Constraint<> sl3("Seller lb3"); //utility of replacing observation done by a bought sensor by this one
    sl3 = y.in_ignore_ith(1, 3, oths_rplc) - (w.in_ignore_ith(1, 1, oths_rplc) - w.in_ignore_ith(0, 1, oths_rplc)) * s.in_ignore_ith(2, 3, tmp_ids) * z.in_ignore_ith(1, 1, tmp_ids1) - p.in_ignore_ith(1, 2, tmp_ids1) * z.in_ignore_ith(1, 1, tmp_ids1);
    model.add(sl3.in(oths_rplc) >= 0);
    
    //For comparison: case with no collaboration (only using own sensors)
    /*Constraint<> no_colab("nc");
    no_colab = sn;
    model.add(no_colab.in(sensors) == 0);*/
    
}

void myModel::mSolve() {
    model.print_constraints_stats(1e-4);
    solver<> sol(model, gurobi);
//    model.set_name("Greedy_read");
//    model.write_solution();
    sol.run();
//    model.set_name("Sensor_assign2");
//    model.write_solution();
    //model.print_solution();
}

void myModel::readGreedySol(string fname) {
    
    auto p = model.get_var<double>("p");
    auto y = model.get_var<double>("y");
    auto s = model.get_var<int>("s");
    auto z = model.get_var<int>("z");
    
    fstream file;
    file.open(fname);
    string id;
    string tmpval;
    for (int i = 0; i < N; i++) {
        file >> tmpval;
        p("sensor" + to_string(i)).set_val(stod(tmpval));
        file >> tmpval;
        y("sensor" + to_string(i)).set_val(stod(tmpval));
    }
    for (int i = 0; i < N; i++) {
        for (Arc* a: graph.get_node("sensor" + to_string(i))->get_out()) {
            file >> id >> tmpval;
            s(id).set_val(stoi(tmpval));
        }
    }
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < K; k++) {
            file >> id >> tmpval;
            z(id).set_val(stoi(tmpval));
        }
    }
    file.close();
}
