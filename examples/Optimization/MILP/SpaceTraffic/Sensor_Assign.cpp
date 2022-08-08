//
//  main.cpp
//  bilevel_sensor
//
//  Created by Svetlana Riabova on 6/8/22.
//

#include <iostream>
#include "Sensor_Assign.hpp"
#include <chrono>
using namespace std::chrono;

int main(int argc, const char * argv[]) {
    /*cout << "Sensors Objects time" << endl;
    for (int i = 1; i < 7; i++) {
        myModel m = myModel();
        auto start = high_resolution_clock::now();
        vector<param<double>> par = m.readData(argc, argv, 2*i-1, 2*i);
        m.InitBilevel(par[0], par[1], par[2]);
        m.mSolve();
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<seconds>(stop - start);
        cout << m.sensors.size() << " " << m.M << " " << duration.count() << endl;
    }*/
    myModel m = myModel();
    //auto start = high_resolution_clock::now();
    vector<param<double>> par = m.readData(argc, argv, 1, 2);
    m.InitBilevel(par[0], par[1], par[2]);
    //m.GreedyStart(par[0], par[1], par[2]);
    m.mSolve();
    auto y = m.model.get_var<double>("y");
    auto psn = m.model.get_var<double>("p_sn");
    ofstream yFile;
    yFile.open("/Users/svetlanariabova/Projects/Sensor/Data/Stats/lb_stats100000.csv");
    ofstream pFile;
    pFile.open("/Users/svetlanariabova/Projects/Sensor/Data/Stats/p_stats100000.csv");
    ofstream ubFile;
    ubFile.open("/Users/svetlanariabova/Projects/Sensor/Data/Stats/ub_stats100000.csv");
    vector<double> y_vals = *y.get_vals();
    vector<double> p_vals = *psn.get_vals();
    double w_max;
    for (int i = 0; i < m.N; i++) {
        yFile << y_vals[i] << endl;
        pFile << p_vals[i] << endl;
        w_max = 0;
        for (Arc* a: m.graph.get_node("sensor" + to_string(i))->get_out()) {
            for (int k = 0; k < m.owner[i]; k++) {
                if (par[2].eval("sensor" + to_string(i) + "," + a->_dest->_name + "," +"agent" + to_string(k)) > w_max) {
                    w_max = par[2].eval("sensor" + to_string(i) + "," + a->_dest->_name + "," +"agent" + to_string(k));
                }
            }
            for (int k = m.owner[i] + 1; k < m.K; k++) {
                if (par[2].eval("sensor" + to_string(i) + "," + a->_dest->_name + "," +"agent" + to_string(k)) > w_max) {
                    w_max = par[2].eval("sensor" + to_string(i) + "," + a->_dest->_name + "," +"agent" + to_string(k));
                }
            }
        }
        ubFile << w_max << endl;
    }
    yFile.close();
    pFile.close();
    ubFile.close();
    //auto stop = high_resolution_clock::now();
    //auto duration = duration_cast<seconds>(stop - start);
    //cout << m.sensors.size() << " " << m.M << " " << duration.count() << endl;
    //cout << duration.count() << endl;
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
    random_device rd; // obtain a random number from hardware
    mt19937 gen(rd()); // seed the generator
    uniform_int_distribution<> distr(0, K-1); // define the range

    for(int i = 0; i < N; i++)
        owner.push_back(distr(gen)); // generate numbers
    
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
    
    jk = indices(objects, range(0,K-1));
    DebugOn("Graph has " << graph.nodes.size() << " nodes" << endl);
    DebugOn("Graph has " << graph.arcs.size() << " arcs" << endl);
    //Parameters
    vector<param<double>> par (3);
    param<double> w0("w0");
    param<double> w_own("w_own");
    param<double> w_bought("w_bought");
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
            owner[i] = stoi(tmp2);
        }
        for (int k = 0; k < K; k++) {
            for (int i = 0; i < N; i++) {
                if (owner[i] == k) {
                    own_sens.add("sensor" + to_string(i) + ",agent" + to_string(k));
                    for (Arc* a: graph.get_node("sensor" + to_string(i))->get_out()) {
                        own_arcs.add("sensor" + to_string(i) + "," + a->_dest->_name +  ",agent" +  to_string(k));
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
                    }
                }
            }
        }
        w0.in(arcs);
        w_own.in(own_arcs);
        w_bought.in(bought_arcs);
        for (Arc* a: graph.arcs) {
            for (int k = 0; k < owner[stoi(a->_src->_name.substr(6, a->_src->_name.find(",")))]; k++) {
                file >> tmp;
                string tmpstr = a->_src->_name + "," + a->_dest->_name + ",agent" + to_string(k);
                w_bought(a->_src->_name + "," + a->_dest->_name + ",agent" + to_string(k)) = stoi(tmp);
            }
            file >> tmp;
            w_own(a->_src->_name + "," + a->_dest->_name + ",agent" + to_string(owner[stoi(a->_src->_name.substr(6, a->_src->_name.find(",")))])) = stod(tmp);
            for (int k = owner[stoi(a->_src->_name.substr(6, a->_src->_name.find(",")))] + 1; k < K; k++) {
                file >> tmp;
                w_bought(a->_src->_name + "," + a->_dest->_name + ",agent" + to_string(k)) = stod(tmp);
            }
        }
//        w0.initialize_normal(2, 1);
    }
    else {
        w0.initialize_normal(2, 1);
        w_own.initialize_normal(2, 1);
        w_bought.initialize_normal(2, 1);
    }
    par[0] = w0;
    par[1] = w_own;
    par[2] = w_bought;
    
    return par;
}

void myModel::InitBilevel(param<double> w0, param<double> w_own, param<double> w_bought) {
    
    double e = 0.001;

    w_own.reset_range();
    w_bought.reset_range();
    /*Variables*/
    var<double> p("p", 0, std::max(w_own._range->second,w_bought._range->second));
    model.add(p.in(sensors));
    var<double> y("y", 0, std::max(w_own._range->second,w_bought._range->second));
    model.add(y.in(sensors));
    var<int> s("s", 0, 1);
    model.add(s.in(own_arcs));
    var<int> sn("sn", 0, 1);
    model.add(sn.in(sensors));
    var<double> p_sn("p_sn", 0, std::max(w_own._range->second,w_bought._range->second));
    model.add(p_sn.in(sensors));
    
    var<double> p_z("p_z", 0, std::max(w_own._range->second,w_bought._range->second));
    model.add(p_z.in(bought_arcs));
    
    var<int> z0("z0", 0, 1);
    model.add(z0.in(arcs));
    var<int> z("z", 0, 1);
    model.add(z.in(bought_arcs));
    
    
    Constraint<> p_sn_def("p_sn_def");
    p_sn_def = p*sn - p_sn;
    model.add(p_sn_def.in(sensors)==0);
    
    
    Constraint<> p_z_def("p_z_def");
    p_z_def = p.in_ignore_ith(1, 2, bought_arcs)*z - p_z;
    model.add(p_z_def.in(bought_arcs)==0);
    
    func<> f = w_own + w0.in_ignore_ith(2, 1, own_arcs);
    f.eval_all();
    param<> w_own0("w_own0");
    w_own0.in(own_arcs);
    w_own0.copy_vals(f);
    
    func<> f2 = w_bought + w0.in_ignore_ith(2, 1, bought_arcs);
    f2.eval_all();
    param<> w_bought0("w_bought0");
    w_bought0.in(bought_arcs);
    w_bought0.copy_vals(f2);
    
    /*Objective*/
    func<> obj;
    obj += product(w_own0, s);
    obj += sum(p_sn);
    obj += product(w_bought0, z);
    obj -= sum(p_z);
    obj -= e * sum(p_z);
//    model.add(obj.in(range(0,0)) == 0);
    model.max(obj);
    
    /*Constraints*/
    //Upper level
//    m.addConstrs((quicksum(z[i, j, k] for k in range(ownr[i]) for j in range(M)) +
//                     quicksum(z[i, j, k] for k in range(ownr[i] + 1, K) for j in range(M)) + quicksum(z0[i, j] for j in range(M))
//                     == sn[i] for i in range(S)), name="unique_bought_obsrvn")
    
    indices z0_ids("z0_ids"), z_ids("z_ids");
    z0_ids = arcs;
    z_ids = bought_arcs;
    for (int i = 0; i<N; i++) {
        for (Arc* b: graph.get_node("sensor" + to_string(i))->get_out()) {
            string j = b->_dest->_name;
            z0_ids.add_in_row(i, "sensor" + to_string(i) + "," + b->_dest->_name);
            for (int k = 0; k < K; k++) {
                if (k != owner[i]) {
                    z_ids.add_in_row(i,"sensor" + to_string(i) + "," + b->_dest->_name + ",agent" + to_string(k));
                }
            }
        }
    }
    
    Constraint<> ub("Unique_Bought_Obsrvn");
    ub = z0.in(z0_ids) + z.in(z_ids) - sn;
    model.add(ub.in(sensors) == 0);

    Constraint<> luo("Leader_Unique_Object");
    luo = sum(z0.in_matrix(0, 1));
    model.add(luo.in(objects) <= 1);
    
    Constraint<> luub("Leader_Utility_ub");
    luub = p.in_ignore_ith(1, 1, arcs) * z0.in(arcs) - w0.in(arcs);
    model.add(luub.in(arcs) <= 0);
    
    //Lower level
        //----Primal Feasibility----

    Constraint<> fua("Unique_Own_Assignment");
    fua = sum(s.in_matrix(1, 1)) + sn.in_ignore_ith(1, 1, own_sens);
    model.add(fua.in(own_sens) <= 1);
    
    Constraint<> fuab("Unique_Bought_Assignment");
    fuab = sum(z.in_matrix(1, 1)) -  sn.in_ignore_ith(1, 1, bought_sens);
    model.add(fuab.in(bought_sens) <= 0);
    
    
    indices c_fub("c_fub"), z_fub("z_fub"), s_fub("s_fub");
    z_fub = *z._indices;
    s_fub = *s._indices;
    size_t row_id = 0;
    bool no_s = true, no_z = true;
    for (int j = 0; j<M; j++) {
        for (int k = 0; k < K; k++) {
            c_fub.insert("object" + to_string(j)+ ",agent" + to_string(k));
            no_s = true, no_z = true;
            for (Arc* a: graph.get_node("object" + to_string(j))->get_in()) {
                int i = stoi(a->_src->_name.substr(6, a->_src->_name.find(",")));
                if (k != owner[i]) {
                    no_z = false;
                    z_fub.add_in_row(row_id, a->_src->_name + "," + a->_dest->_name + ",agent" + to_string(k));
                }
                else{
                    no_s = false;
                    s_fub.add_in_row(row_id, a->_src->_name + "," + a->_dest->_name + ",agent" + to_string(k));
                }
            }
            if(no_z)
                z_fub.add_empty_row();
            if(no_s)
                s_fub.add_empty_row();
            row_id++;
        }
    }
    
    Constraint<> fub("Follower_Unique_Object");
    fub = s.in(s_fub) + z.in(z_fub);
    model.add(fub.in(c_fub) <= 1);
    
    Constraint<> fuub("Followers_Utility_ub");
    fuub = p.in_ignore_ith(1, 2, bought_arcs) * z - w_bought;
    model.add(fuub.in(bought_arcs) <= 0);
    
        //----Fair price----
    Constraint<> fp("FairPrice");
    fp = p.in_ignore_ith(1, 2, bought_arcs) - (w_bought*z + y.in_ignore_ith(1, 2, bought_arcs))/2;
    model.add(fp.in(bought_arcs) >= 0);

    //indices s_ids = s.get_matrix_ids(0, 1);
    //s_ids.print();
    
    indices c_lb1("c_lb1"), y_lb1("y_lb1"), w_own_lb1("w_own_lb1"), w_own_z_lb1("w_own_z_lb1"), w_own_s_lb1("w_own_s_lb1"), z_lb1("z_lb1"), s_lb1("s_lb1");
    y_lb1 = sensors;
    z_lb1 = bought_arcs;
    s_lb1 = own_arcs;
    w_own_lb1 = own_arcs;
    w_own_s_lb1 = own_arcs;
    w_own_z_lb1 = own_arcs;
    row_id = 0;
    for (int i = 0; i < N; i++) {
        for (Arc* b: graph.get_node("sensor" + to_string(i))->get_out()) {
            no_z = true;
            string j = b->_dest->_name;
            c_lb1.add("Seller lb1:" + to_string(i) + "," + b->_dest->_name + ",agent" + to_string(owner[i]));
            y_lb1.add_ref("sensor" + to_string(i));
            w_own_lb1.add_ref("sensor" + to_string(i) + "," + b->_dest->_name + ",agent" + to_string(owner[i]));
            for (Arc* a: graph.get_node(j)->get_in()) {
                if (owner[stoi(a->_src->_name.substr(6, a->_src->_name.find(",")))] == owner[i]) {
                    w_own_s_lb1.add_in_row(row_id, "sensor" + to_string(i) + "," + b->_dest->_name + ",agent" + to_string(owner[i]));
                    s_lb1.add_in_row(row_id,a->_src->_name + "," + j +  ",agent" + to_string(owner[i]));
                }
                else if (owner[stoi(a->_src->_name.substr(6, a->_src->_name.find(",")))] != owner[i]) {
                    w_own_z_lb1.add_in_row(row_id, "sensor" + to_string(i) + "," + b->_dest->_name + ",agent" + to_string(owner[i]));
                    z_lb1.add_in_row(row_id, a->_src->_name + "," + j +  ",agent" + to_string(owner[i]));
                    no_z = false;
                }
            }
            if(no_z){
                z_lb1.add_empty_row();
            }
            row_id++;
        }
    }
    
    Constraint<> sl1("Seller lb1");
    sl1 = y.in(y_lb1) + w_own.in(w_own_s_lb1)*s.in(s_lb1) + w_own.in(w_own_z_lb1)*z.in(z_lb1) - w_own.in(w_own_lb1);
    model.add(sl1.in(c_lb1) >= 0);
    
    //w_own.print();
    //w_bought.print();
    
//    own_rplc.print();
    Constraint<> sl2("Seller lb2");
    sl2 = y.in_ignore_ith(1, 3, own_rplc) - (w_own.in_ignore_ith(1, 1, own_rplc) - w_own.in_ignore_ith(0, 1, own_rplc)) * s.in_ignore_ith(0, 1, own_rplc);
    model.add(sl2.in(own_rplc) >= 0);
    
    Constraint<> sl3("Seller lb3");
    sl3 = y.in_ignore_ith(1, 3, oths_rplc) - (w_own.in_ignore_ith(1, 1, oths_rplc) - w_bought.in_ignore_ith(0, 1, oths_rplc))*z.in_ignore_ith(0, 1, oths_rplc) - p_z.in_ignore_ith(0, 1, oths_rplc);
    model.add(sl3.in(oths_rplc) >= 0);
    
    //For comparison
    /*Constraint<> no_colab("nc");
    no_colab = sn;
    model.add(no_colab.in(sensors) == 0);*/
    
    
//    model.print_symbolic();
 //   model.print();
//    model.write("before.txt");
//    model.replace_integers();
//    model.restructure();
//    model.write(3);
    //model.print();
//    model.replace_integers();
//    auto R = model.relax();
//    R->print();
//    solver<> sol(model, gurobi);
//    model.restructure();
//
    //int time_limit = 300;//seconds
 //   sol.run();//(1e-5, time_limit);
//    model.print_solution();
//    model.print_constraints_stats(1e-4);
//    model.print();
    
    //return &model;
}

void myModel::mSolve() {
    solver<> sol(model, gurobi);
    sol.run();
}

void myModel::GreedyStart(param<double> w0, param<double> w_own, param<double> w_bought) {

    param<double> wt0 = w0.deep_copy();
    param<double> wt_own = w_own.deep_copy();
    param<double> wt_bought = w_bought.deep_copy();
    auto s = model.get_var<int>("s");
    auto sn = model.get_var<int>("sn");
    auto z0 = model.get_var<int>("z0");
    auto z = model.get_var<int>("z");
    auto p = model.get_var<double>("p");
    auto y = model.get_var<double>("y");

    
    string idx1;
    string idx2;
    string idx3;
    double m1;
    double m2;
    double m3;
    int ownr;
    int sensor;
    int object;
    int agent;

    while(parSum(wt0) + parSum(wt_own) + parSum(wt_bought) > 0) {
        idx1 = findMax(wt0);
        idx2 = findMax(wt_own);
        idx3 = findMax(wt_bought);
        m1 = wt0.eval(idx1);
        m2 = wt_own.eval(idx2);
        m3 = wt_bought.eval(idx3);
        if (m1 >= m2) {
            if (m1 >= m3) {
                //assign leader
                assignLeader(idx1, wt0, wt_own, wt_bought);
            }
            else {
                //assign bought
                assignBought(idx3, wt0, wt_own, wt_bought);
            }
        }
        else if (m2 >= m3) {
            //assign own
            assignOwn(idx2, wt0, wt_own, wt_bought);
        }
        else {
            //assign bought
            assignBought(idx3, wt0, wt_own, wt_bought);
        }
    }
    
    //fix prices
    double y1 = 0;
    double y2 = 0;
    double y3 = 0;
    double tmp_wt = 0;
    int tmp_add = 0;
    vector<int> srs_k;
    vector<int> srs_oths;
    string j;
    int t;
    
    for (int i = 0; i < N; i++) {
        for (Arc* b: graph.get_node("sensor" + to_string(i))->get_out()) {
            j = b->_dest->_name;
            srs_k.clear();
            for (Arc* a: graph.get_node(j)->get_in()) {
                t = stoi(a->_src->_name.substr(6, a->_src->_name.length()));
                if (owner[t] == owner[i]) { srs_k.push_back(t); }
                else { srs_oths.push_back(t); }
            }
            for (int t: srs_k) {
                if (s("sensor" + to_string(t) + "," + j + "," + "agent" + to_string(owner[i])).getvalue() > 0) {
                    tmp_add = 1;
                }
                if (tmp_add > 0) { break; }
            }
            for (int r: srs_k) {
                if (z("sensor" + to_string(r) + "," + j + "," + "agent" + to_string(owner[i])).getvalue() > 0) {
                    tmp_add = 1;
                }
                if (tmp_add > 0) { break; }
            }
            if (y1 < w_own.eval("sensor" + to_string(i) + "," + j + "," + "agent" + to_string(owner[i]))) {
                y1 = w_own.eval("sensor" + to_string(i) + "," + j + "," + "agent" + to_string(owner[i]));
            }
        }
        if (y1 > y2) {
            if (y1 > y3) { p("sensor" + to_string(i)).set_val(y1); y("sensor" + to_string(i)).set_val(y1); }
            else { p("sensor" + to_string(i)).set_val(y3); y("sensor" + to_string(i)).set_val(y3); }
        }
        else {
            if (y2 > y3) { p("sensor" + to_string(i)).set_val(y2); y("sensor" + to_string(i)).set_val(y2); }
            else { p("sensor" + to_string(i)).set_val(y3); y("sensor" + to_string(i)).set_val(y3); }
        }
    }

}


void myModel::assignLeader(string &idx, param<double> wt0, param<double> wt_own, param<double> wt_bought) {
    auto sn = model.get_var<int>("sn");
    auto z0 = model.get_var<int>("z0");
    int ownr;
    int sensor;
    int object;
    string j;
    
    ownr = owner[stoi(idx.substr(6, idx.find(",")))];
    sensor = stoi(idx.substr(6, idx.find(",")));
    cout << "z0:" << idx << endl;
    object = stoi(idx.substr(idx.find("t") + 1, idx.find(",")));
    z0(idx.substr(0, nthOccurrence(idx, ",", 2))).set_val(1);
    //sensor not used for other objs leader + owner
    for (Arc* b: graph.get_node("sensor" + to_string(sensor))->get_out()) {
        j = b->_dest->_name;
        wt_own("sensor" + to_string(sensor) + "," + "object" + j + "," + "agent" + to_string(ownr)).set_val(0);
        wt0("sensor" + to_string(sensor) + "," + "object" + j).set_val(0);
    }
    //object not observed twice
    for (int i = 0; i < sensor; i ++) {
        wt0("sensor" + to_string(i) + "," + "object" + to_string(object)).set_val(0);
    }
    for (int i = sensor + 1; i < N; i ++) {
        wt0("sensor" + to_string(i) + "," + "object" + to_string(object)).set_val(0);
    }
    //sensor not used by other agents
    for (int k = 0; k < ownr; k++) {
        for (Arc* b: graph.get_node("sensor" + to_string(sensor))->get_out()) {
            j = b->_dest->_name;
            wt_bought("sensor" + to_string(sensor) + "," + "object" + j + "," + "agent" + to_string(k)).set_val(0);
        }
    }
    for (int k = ownr + 1; k < K; k++) {
        for (Arc* b: graph.get_node("sensor" + to_string(sensor))->get_out()) {
            j = b->_dest->_name;
            wt_bought("sensor" + to_string(sensor) + "," + "object" + j + "," + "agent" + to_string(k)).set_val(0);
        }
    }
    //sold
    sn("sensor" + to_string(sensor)).set_val(1);
}

void myModel::assignOwn(string &idx, param<double> wt0, param<double> wt_own, param<double> wt_bought) {
    auto s = model.get_var<int>("s");
    int ownr;
    int sensor;
    string object;
    string j;
    string i;
    
    ownr = owner[stoi(idx.substr(6, idx.find(",")))];
    sensor = stoi(idx.substr(6, idx.find(",")));
    object = idx.substr(idx.find(",") + 1, idx.find(","));
    cout << "s:" << idx << endl;
    s(idx).set_val(1);
    //sensor not used by leader
    for (Arc* b: graph.get_node("sensor" + to_string(sensor))->get_out()) {
        j = b->_dest->_name;
        wt0("sensor" + to_string(sensor) + "," + j).set_val(0);
        //sensor not used by other agents
        for (int k = 0; k < owner[sensor]; k++) {
            wt_bought("sensor" + to_string(sensor) + "," + j + "," + "agent" + to_string(k)).set_val(0);
        }
        for (int k = owner[sensor] + 1; k < K; k++) {
            wt_bought("sensor" + to_string(sensor) + "," + j + "," + "agent" + to_string(k)).set_val(0);
        }
    }
    //object not observed twice
    for (Arc* a: graph.get_node(object)->get_in()) {
        i = a->_src->_name;
        wt_own(i + "," + object + "," + "agent" + to_string(ownr)).set_val(0);
    }
}

void myModel::assignBought(string &idx, param<double> wt0, param<double> wt_own, param<double> wt_bought) {
    auto sn = model.get_var<int>("sn");
    auto z = model.get_var<int>("z");
    int ownr;
    int sensor;
    int object;
    int agent;
    string j;
    string i;
    
    ownr = owner[stoi(idx.substr(6, idx.find(",")))];
    sensor = stoi(idx.substr(6, idx.find(",")));
    object = stoi(idx.substr(idx.find("t") + 1, idx.find(",")));
    agent = stoi(idx.substr(nthOccurrence(idx, "t", 2) + 1, idx.length()));
    cout << "z:" << idx << endl;
    z(idx).set_val(1);
    //sensor not used for other objs leader + owner
    for (Arc* b: graph.get_node("sensor" + to_string(sensor))->get_out()) {
        j = b->_dest->_name;
        wt_own("sensor" + to_string(sensor) + "," + j + "," + "agent" + to_string(ownr)).set_val(0);
        wt0("sensor" + to_string(sensor) + "," + j).set_val(0);
    }
    //object not observed twice
    for (Arc* a: graph.get_node("object" + to_string(object))->get_in()) {
        i = a->_src->_name;
        if (owner[stoi(i.substr(6, i.length()))] == agent) { wt_own(i + "," + "object" + to_string(object) + "," + "agent" + to_string(agent)).set_val(0); }
        else { wt_bought(i + "," + "object" + to_string(object) + "," + "agent" + to_string(agent)).set_val(0); }
    }
    //sensor not used by other agents
    for (Arc* b: graph.get_node("sensor" + to_string(sensor))->get_out()) {
        j = b->_dest->_name;
        for (int k = 0; k < std::min(agent, owner[sensor]); k++) {
            wt_bought("sensor" + to_string(sensor) + "," + j + "," + "agent" + to_string(k)).set_val(0);
        }
        for (int k = std::min(agent, owner[sensor]) + 1; k < std::max(agent, owner[sensor]); k++) {
            wt_bought("sensor" + to_string(sensor) + "," + j + "," + "agent" + to_string(k)).set_val(0);
        }
        for (int k = std::max(agent, owner[sensor]) + 1; k < K; k++) {
            wt_bought("sensor" + to_string(sensor) + "," + j + "," + "agent" + to_string(k)).set_val(0);
        }
    }
    //sold
    sn("sensor" + to_string(sensor)).set_val(1);
}

double myModel::parSum(param<double> w) {
    double s = 0;
    for (auto& n : *w.get_vals()) {
        s += n;
    }
    return s;
}

string myModel::findMax(param<double> w) {
    string max_idx = (*w.get_keys())[0];
    double max_el = 0;
    for (auto i: *w.get_keys()) {
        if (w.eval(i) > max_el) {
            max_el = w.eval(i);
            max_idx = i;
        }
    }
    return max_idx;
}

int myModel::nthOccurrence(const std::string& str, const std::string& findMe, int nth)
{
    size_t  pos = -1;
    int     cnt = 0;

    while( cnt != nth )
    {
        pos+=1;
        pos = str.find(findMe, pos);
        if ( pos == std::string::npos )
            return -1;
        cnt++;
    }
    return pos;
}
