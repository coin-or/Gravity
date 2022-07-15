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
    vector<param<double>> par = m.readData(argc, argv);
    m.InitBilevel(par[0], par[1], par[2]);
    //int s = 0;
    return 0;
}

vector<param<double>> myModel::readData(int argc, const char * argv[]){
    N = 4; M = 7; K = 3;
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
    
    //Parameters
    vector<param<double>> par (3);
    param<double> w0("w0");
    param<double> w_own("w_own");
    param<double> w_bought("w_bought");
    if (argc >= 3) {
        string fname = argv[2];
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
                    own_sens.add("sensor" + to_string(i) + "," + to_string(k));
                    for (Arc* a: graph.get_node("sensor" + to_string(i))->get_out()) {
                        own_arcs.add("sensor" + to_string(i) + "," + a->_dest->_name +  "," +  to_string(k));
                        for (Arc* b: graph.get_node(a->_dest->_name)->get_in()) {
                            if ((owner[stoi(b->_src->_name.substr(6, b->_src->_name.find(",")))] == k) && (stoi(b->_src->_name.substr(6, b->_src->_name.find(","))) != i)) {
                                own_rplc.add(a->_src->_name + "," + b->_src->_name + "," + a->_dest->_name +  "," +  to_string(k));
                            }
                            else if (stoi(b->_src->_name.substr(6, b->_src->_name.find(","))) != i) {
                                oths_rplc.add("sensor" + to_string(i) + "," + b->_src->_name + "," + a->_dest->_name +  "," +  to_string(k));
                            }
                        }
                    }
                }
                else {
                    bought_sens.add("sensor" + to_string(i) + "," + to_string(k));
                    for (Arc* a: graph.get_node("sensor" + to_string(i))->get_out()) {
                        bought_arcs.add(a->_src->_name + "," + a->_dest->_name +  "," +  to_string(k));
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
                w_bought(a->_src->_name + "," + a->_dest->_name + "," + to_string(k)) = stoi(tmp);
            }
            file >> tmp;
            w_own(a->_src->_name + "," + a->_dest->_name + "," + to_string(owner[stoi(a->_src->_name.substr(6, a->_src->_name.find(",")))])) = stod(tmp);
            for (int k = owner[stoi(a->_src->_name.substr(6, a->_src->_name.find(",")))] + 1; k < K; k++) {
                file >> tmp;
                w_bought(a->_src->_name + "," + a->_dest->_name + "," + to_string(k)) = stod(tmp);
            }
        }
        w0.initialize_normal(2, 1);
    }
    /*else {
        w0.initialize_normal(2, 1);
        w_own.initialize_normal(2, 1);
        w_bought.initialize_normal(2, 1);
    }*/
    par[0] = w0;
    par[1] = w_own;
    par[2] = w_bought;
    
    return par;
}

void myModel::InitBilevel(param<double> w0, param<double> w_own, param<double> w_bought) {

    Model<> model("BilevelSensor");
    
    double e = 0.001;

    /*Variables*/
    var<double> p("p", pos_);
    model.add(p.in(sensors));
    var<double> y("y", pos_);
    model.add(y.in(sensors));
    var<int> s("s", 0, 1);
    model.add(s.in(own_arcs));
    var<int> sn("sn", 0, 1);
    model.add(sn.in(sensors));
    var<int> z0("z0", 0, 1);
    model.add(z0.in(arcs));
    var<int> z("z", 0, 1);
    model.add(z.in(bought_arcs));

    var<> obj_var("obj_var", 0,10000);
    model.add(obj_var.in(range(0,0)));
    
    param<> ones("ones");
    ones.in(sensors);
    ones = 1;
    
    
    /*Objective*/
    Constraint<> obj("obj_def");
    obj += product(w_own + w0.in_ignore_ith(2, 1, own_arcs), s);
    obj += ones.tr()*p*sn;
    obj += product(w_bought + w0.in_ignore_ith(2, 1, bought_arcs), z);
    obj -= ones.in_ignore_ith(1, 2, bought_arcs).tr()*p.in_ignore_ith(1, 2, bought_arcs)*z;
    obj -= e * ones.in_ignore_ith(1, 2, bought_arcs).tr()*p.in_ignore_ith(1, 2, bought_arcs)*z;
    obj -= obj_var;
    model.add(obj.in(range(0,0)) == 0);
    
    model.max(obj_var);
    
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
    fp = p.in_ignore_ith(1, 2, bought_arcs) - (w_bought*z + y.in_ignore_ith(1, 2, bought_arcs))/2;
    model.add(fp.in(bought_arcs) >= 0);

    //indices s_ids = s.get_matrix_ids(0, 1);
    //s_ids.print();
    /*Constraint<> sl1("Seller lb1");
    sl1 = y.in_ignore_ith(1, 2, own_arcs) - w_own * (1 - sum(s.sum_over(own_rplc.ignore_ith(0, 1), 0)) - sum(z.sum_over(own_rplc.ignore_ith(0, 1), 0)));
    model.add(sl1.in(own_arcs) >= 0);*/
    
    //w_own.print();
    //w_bought.print();
    for (int i = 0; i < N; i++) {
        for (Arc* b: graph.get_node("sensor" + to_string(i))->get_out()) {
            string j = b->_dest->_name;
            Constraint<> sl1("Seller lb1:" + to_string(i) + "," + b->_dest->_name + "," + to_string(owner[i]));
            sl1 += y("sensor" + to_string(i)) - w_own("sensor" + to_string(i) + "," + b->_dest->_name + "," + to_string(owner[i]));
            for (Arc* a: graph.get_node(j)->get_in()) {
                if ((owner[stoi(a->_src->_name.substr(6, a->_src->_name.find(",")))] == owner[i]) && (stoi(a->_src->_name.substr(6, a->_src->_name.find(","))) != i)) {
                    sl1 += w_own("sensor" + to_string(i) + "," + b->_dest->_name + "," + to_string(owner[i])) * s(a->_src->_name + "," + j +  "," + to_string(owner[i]));
                }
                else if (owner[stoi(a->_src->_name.substr(6, a->_src->_name.find(",")))] != owner[i]) {
                    sl1 += w_own("sensor" + to_string(i) + "," + b->_dest->_name + "," + to_string(owner[i])) * z(a->_src->_name + "," + j +  "," + to_string(owner[i]));
                }
            }
            model.add(sl1.in(range(0,0)) >= 0);
        }
    }
    
    
    Constraint<> sl2("Seller lb2");
    sl2 = y.in_ignore_ith(1, 3, own_rplc) - (w_own.in_ignore_ith(1, 1, own_rplc) - w_own.in_ignore_ith(0, 1, own_rplc)) * s.in_ignore_ith(0, 1, own_rplc);
    model.add(sl2.in(own_rplc) >= 0);
    
    Constraint<> sl3("Seller lb3");
    sl3 = y.in_ignore_ith(1, 3, oths_rplc) - (w_own.in_ignore_ith(1, 1, oths_rplc) - 1)*z.in_ignore_ith(0, 1, oths_rplc) - p.in_ignore_ith(1, 3, oths_rplc) * z.in_ignore_ith(0, 1, oths_rplc);
    model.add(sl3.in(oths_rplc) >= 0);
    
    
    model.print();
    solver<> sol(model, ipopt);
    model.restructure();
//
    sol.run();
    //
//    model.print();

}
