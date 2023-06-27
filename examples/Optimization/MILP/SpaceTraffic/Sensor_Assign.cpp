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

#ifdef USE_H5CPP
#include <h5cpp/hdf5.hpp>
using namespace hdf5;
#endif


int main(int argc, const char * argv[]) {
    bool run_MIP = true;//false;
    string fname = "/Users/svetlanariabova/Projects/Sensor/Data/hdf5/new_weights_0.h5";//string(prj_dir)+"/data_sets/sensor/test.h5";
//    if(argc<2){
//        DebugOn("Please enter the input file path as first argument, using the default file test.h5 found under Gravity/datasets/sensor\n");
//    }
//    else{
//        fname = argv[1];
//    }
    if(argc>3)
        run_MIP = true;
    myModel m = myModel();
    auto par = m.readHD5(fname);
//    vector<param<double>> par = m.readData(argc, argv, 1, 2);
    auto start = high_resolution_clock::now();
    m.InitBilevel2(par[0], par[1], par[2], 0.001);
//    m.GreedyStart(par[0], par[1], par[2]); //comment if no greedy start not needed
//    cout << "Writing solution file." << endl;
//    m.writeGreedySol(); //writing greedy sol to a file to load it to sensor_assign2
    auto stop = high_resolution_clock::now();
    auto duration1 = duration_cast<seconds>(stop - start);
    cout << "Init + greedy time: " << duration1.count() << " seconds." << endl;
//    bool run_MIP = true;
    m.mSolve(run_MIP);
    //m.writeGreedySol(); //writing opt sol to a file to load it to sensor_assign2
    auto stop2 = high_resolution_clock::now();
    auto duration2 = duration_cast<seconds>(stop2 - stop);
//    cout << "Done." << endl;
//    cout << m.N << " " << m.M << " " << m.K << " " << duration1.count() + duration2.count() << endl; //prints num sensors; num objetcs; num agents; total time after reading input (init + greedy + solver)
    return 0;
}

vector<param<double>> myModel::readHD5(const string& fname){

#ifdef USE_H5CPP
    auto hd5file = file::open(fname);
    auto RootGroup = hd5file.root();
    /* Read agent names */
    auto agents_group = RootGroup.get_group("agents");
    auto agents_names = agents_group.get_dataset("agent_name");
    agents = indices("agents");
    dataspace::Simple Dataspace(agents_names.dataspace());
    auto Dimensions = Dataspace.current_dimensions();
    auto nb_agents = Dimensions[0];
//    budget.resize(nb_agents);
//    auto budget_set = agents_group.get_dataset("budget");
//    budget_set.read(budget);
    auto MaxDimensions = Dataspace.maximum_dimensions();
    std::vector<string> AllElements(Dataspace.size());
    agents_names.read(AllElements);
    for (auto Value : AllElements) {
        agents.insert("Agent_"+to_string(agents.size()));
    }
    /* Read objects */
    auto objects_group = RootGroup.get_group("objects");
    auto object_names_set = objects_group.get_dataset("object_name");
    auto priority_set = objects_group.get_dataset("priority");
    dataspace::Simple object_names_space(object_names_set.dataspace());
    objects = indices("objects");
    auto nb_objects = object_names_space.size();
    std::vector<string> object_names(nb_objects);
    object_names_set.read(object_names);
    for (auto name : object_names) {
        objects.insert("Object_"+to_string(objects.size()));
    }
    auto pr_size = priority_set.dataspace().size();
    assert(pr_size==nb_objects*nb_agents);
    vector<double> priority(pr_size);
    priority_set.read(priority);
    /* Read sensors */
    auto sensors_group = RootGroup.get_group("sensors");
    auto sensors_names_set = sensors_group.get_dataset("sensor_name");
    auto agent_own_set = sensors_group.get_dataset("agent_id");
    dataspace::Simple sensors_names_space(sensors_names_set.dataspace());
    dataspace::Simple agent_own_space(agent_own_set.dataspace());
    sensors = indices("sensors");
    own_sens = indices("own_sens");

    auto nb_sensors = sensors_names_space.size();
    owner.resize(nb_sensors);/*< vector storing ownership of each sensor */
    vector<string> sensors_names(nb_sensors);
    sensors_names_set.read(sensors_names);
    for (auto name : sensors_names) {
        sensors.insert("Sensor_"+to_string(sensors.size()));
    }
    vector<int> agent_own(nb_sensors);
    agent_own_set.read(agent_own);
    int s_id = 0;
    for (auto agent_id : agent_own) {
        owner[s_id] = agent_id;
        owner_map[sensors.get_key(s_id)] = agent_id;
        own_sens.insert(sensors.get_key(s_id++)+","+agents.get_key(agent_id));
    }
//    agents.print();
//    objects.print();
//    sensors.print();
//    own_sens.print();
    /* Read arcs */
    auto arcs_group = RootGroup.get_group("arcs");
    auto object_id_set = arcs_group.get_dataset("object_id");
    auto sensor_id_set = arcs_group.get_dataset("sensor_id");
    auto quality_set = arcs_group.get_dataset("quality");
    auto nb_arcs = object_id_set.dataspace().size();
    vector<int> object_ids(nb_arcs), sensor_ids(nb_arcs);
    vector<double> quality(nb_arcs);
    object_id_set.read(object_ids);
    sensor_id_set.read(sensor_ids);
    quality_set.read(quality);
    Node* node = NULL;
    Arc* arc = NULL;
    string src, dest;

    for (int i = 0; i <nb_sensors; i++){
        node = new Node(sensors.get_key(i),graph.nodes.size());
        node->owner_id = owner[i];
        graph.add_node(node);
    }
    for (int i = 0; i <nb_objects; i++){
        node = new Node(objects.get_key(i),graph.nodes.size());
        graph.add_node(node);
    }

    for (int i = 0; i < nb_arcs; i++) {
        src = sensors.get_key(sensor_ids[i]);
        dest = objects.get_key(object_ids[i]);
        arc = new Arc(src + "," + dest);
        arc->_id = i;
        arc->_src = graph.get_node(src);
        arc->_dest= graph.get_node(dest);
        arc->weight = quality[i];
        graph.add_arc(arc);
        arc->connect(false);
    }
//    graph.print();
    DebugOn("Done reading h5 file and building graph" << endl);
    DebugOn("Graph has " << graph.nodes.size() << " nodes" << endl);
    DebugOn("Graph has " << graph.arcs.size() << " arcs" << endl);
    DebugOn("Number of objects =  " << nb_objects << endl);
    DebugOn("Number of sensors =  " << nb_sensors << endl);
    DebugOn("Number of agents =  " << nb_agents << endl);
    N = nb_sensors;
    M = nb_objects;
    K = nb_agents;

    /* Graph nodes are indexed in {0,...,n+m-1})*/
    assert(N + M == graph.nodes.size());/* Make sure we have the right number of nodes */
    arcs = indices("arcs");
    own_arcs = indices("own_arcs");
    bought_arcs = indices("bought_arcs");

    /* Parameters */
    vector<param<double>> par (3);
    param<double> w0("w0");
    param<double> w_own("w_own");
    param<double> w_bought("w_bought");
    arcs.add(graph.arcs);
    w0.in(arcs);
    w_own.in(own_arcs);
    w_bought.in(bought_arcs);
    DebugOn("Building coefficient matrix..." << endl);
    string sensor_name, agent_name, object_name;
    /*define index sets for weights, then read weights (need owner to sepsrate own and bought weights)*/
    for (int k = 0; k < K; k++) {
        agent_name = agents.get_key(k);
        for (int i = 0; i < N; i++) {
            sensor_name = sensors.get_key(i);
            auto sensor_node = graph.get_node(sensor_name);
            if (owner[i] == k) {
                for (const Arc* a: sensor_node->get_out()) {
                    own_arcs.add(sensor_name + "," + a->_dest->_name + "," + agent_name);
                    w_own.add_val(a->weight*priority[a->_dest->_id+nb_objects*k]);
//                    for (const Arc* b: a->_dest->get_in()) {
//                        if (b->_src->owner_id == k) {
//                            if (b->_src->_id != i) {
//                                own_rplc.add(a->_src->_name + "," + b->_src->_name + "," + a->_dest->_name + "," + agent_name);
//                            }
//                        }
//                        else {
//                            oths_rplc.add(sensor_name + "," + b->_src->_name + "," + a->_dest->_name + "," + agent_name);
//                        }
//                    }
                }
            }
            else {
                bought_sens.add(sensor_name + "," + agent_name);
                for (const Arc* a: sensor_node->get_out()) {
                    bought_arcs.add(a->_src->_name + "," + a->_dest->_name +  "," + agent_name);
                    w_bought.add_val(a->weight*priority[a->_dest->_id+nb_objects*k]);
                }
            }
        }
    }


    /*pass weights to init*/
    par[0] = w0;
    par[1] = w_own;
    par[2] = w_bought;

    return par;
#else
    cerr << "Can't read Hd5 as a solver: this version of Gravity was compiled without H5CPP. Rerun cmake with -DH5CPP=ON." << endl;
    exit(1);
#endif
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
        graph.generate_bipartite_random(N, M, degree);
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
        objects2.add("object" + to_string(i - N));
    }
    //arcs = indices("arcs");
    own_arcs = indices("own_arcs");
    bought_arcs = indices("bought_arcs");
    
    DebugOn("Graph has " << graph.nodes.size() << " nodes" << endl);
    DebugOn("Graph has " << graph.arcs.size() << " arcs" << endl);
    
    //Parameters
    vector<param<double>> par (3);
    param<double> w0("w0");
    param<double> w_own("w_own");
    param<double> w_bought("w_bought");
    
    if (argc >= 3) {
        /*read data from file*/
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
        for (int k = 0; k < K; k++) {
            agents.add("agent" + to_string(k));
        }
        
        /*define index sets for weights, then read weights (need owner to sepsrate own and bought weights)*/
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
                    bought_sens.add("sensor" + to_string(i) + ",agent" + to_string(k));
                    for (Arc* a: graph.get_node("sensor" + to_string(i))->get_out()) {
                        bought_arcs.add(a->_src->_name + "," + a->_dest->_name +  ",agent" +  to_string(k));
                    }
                }
            }
        }
//        own_sens.print();
        w0.in(arcs);
        w_own.in(own_arcs);
        w_bought.in(bought_arcs);
        for (Arc* a: graph.arcs) {
            for (int k = 0; k < owner[stoi(a->_src->_name.substr(6, a->_src->_name.find(",")))]; k++) {
                file >> tmp;
                w_bought(a->_src->_name + "," + a->_dest->_name + ",agent" + to_string(k)) = stod(tmp);
            }
            file >> tmp;
            w_own(a->_src->_name + "," + a->_dest->_name + ",agent" + to_string(owner[stoi(a->_src->_name.substr(6, a->_src->_name.find(",")))])) = stod(tmp);
            for (int k = owner[stoi(a->_src->_name.substr(6, a->_src->_name.find(",")))] + 1; k < K; k++) {
                file >> tmp;
                w_bought(a->_src->_name + "," + a->_dest->_name + ",agent" + to_string(k)) = stod(tmp);
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
        w_own.initialize_normal(2, 1);
        w_bought.initialize_normal(2, 1);
    }
    
    /*pass weights to init*/
    par[0] = w0;
    par[1] = w_own;
    par[2] = w_bought;
    
    return par;
}

void myModel::InitBilevel(param<double> &w0, param<double> &w_own, param<double> &w_bought, double eps) {
    
    e = eps; //for e * p_sn in obj
    
    /*reset_range to define ub on prices*/
    w_own.reset_range();
    w_bought.reset_range();
    w0.reset_range();
    
    /*---------Variables----------*/
    /*price + aux for price*/
    var<double> p("p", 0, std::max(w_own._range->second,w_bought._range->second));
    model.add(p.in(sensors));
    var<double> y("y", 0, std::max(w_own._range->second,w_bought._range->second));
    model.add(y.in(sensors));
    
    /*assignment*/
    var<int> s("s", 0, 1); //use own sens
    model.add(s.in(own_arcs));
    var<int> sn("sn", 0, 1); //sensor sold
    model.add(sn.in(sensors));
    var<int> z0("z0", 0, 1); //leader buys sens (observation)
    model.add(z0.in(arcs));
    var<int> z("z", 0, 1); //agebts buy sens (observation)
    model.add(z.in(bought_arcs));
    
    /*multiplication vars*/
    var<double> p_sn("p_sn", 0, std::max(w_own._range->second,w_bought._range->second));
    model.add(p_sn.in(sensors));
    var<double> p_z("p_z", 0, std::max(w_own._range->second,w_bought._range->second));
    model.add(p_z.in(bought_arcs));
    /* experiimental */
    var<double> p_z0("p_z0", 0, std::max(w_own._range->second,w0._range->second));
    model.add(p_z0.in(arcs));
    
//    /*tmp var to see difference in price terms*/
//    var<double> p_diff("p_diff", 0, std::max(w_own._range->second,w_bought._range->second));
//    model.add(p_diff);
    
    /*multiplication vars def*/
    Constraint<> p_sn_def("p_sn_def");
    p_sn_def = p*sn - p_sn;
    model.add(p_sn_def.in(sensors)==0);
    
    Constraint<> p_z_def("p_z_def");
    p_z_def = p.in_ignore_ith(1, 2, bought_arcs)*z - p_z;
    model.add(p_z_def.in(bought_arcs)==0);
    
    /* p*z0 for cancelling out price terms in the obj (experimental) */
    Constraint<> p_z0_def("p_z0_def");
    p_z0_def = p.in_ignore_ith(1, 1, arcs) * z0 - p_z0;
    model.add(p_z0_def.in(arcs) == 0);
    
//    /*p-diff def*/
//    Constraint<> p_diff_def("p_diff_def");
//    p_diff_def = p_diff - sum(p_sn) + sum(p_z) + sum(p_z0);
//    model.add(p_diff_def == 0);
    
    /*weights for obj (own + leader, bought + leader); leader receives utility when an object is observed, no matter by whom*/
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
    
    /*--------Objective---------*/
    func<> obj;
    obj += product(w_own0, s); //use own sens
    obj += sum(p_sn); //sell sens
    obj += product(w_bought0, z); //buy sens pt.1
    obj += product(w0, z0);
    obj -= sum(p_z); //buy sens pt.2
    obj -= sum(p_z0);
    //obj += p_diff;
    obj -= e * sum(p_sn); //regularization term; sets prices to their lb from Fair_price constraints (can't use equality there)
    model.max(obj);
    
    /*--------Constraints---------*/
    
    /*matching index sets for z0 and z in Unique_Bought_Obsrvn*/
    indices z0_ids("z0_ids"), z_ids("z_ids"), ub_idx("ub_idx"), ub_idx2("ub_idx2");
    z0_ids = arcs;
    z_ids = bought_arcs;
    ub_idx = sensors;
    ub_idx2 = sensors;

    for (int i = 0; i<N; i++) {
        for (Arc* b: graph.get_node("Sensor_" + to_string(i))->get_out()) {
            string j = b->_dest->_name;
            z0_ids.add_in_row(i, "Sensor_" + to_string(i) + "," + b->_dest->_name);
            for (int k = 0; k < K; k++) {
                if (k != owner[i]) {
                    z_ids.add_in_row(i,"Sensor_" + to_string(i) + "," + b->_dest->_name + ",Agent_" + to_string(k));
                }
            }
        }
        if (graph.get_node("Sensor_" + to_string(i))->get_out().size() > 0) {
            ub_idx.add_in_row(i, "Sensor_" + to_string(i));
        }
        else {
            ub_idx2.add_in_row(i, "Sensor_" + to_string(i));
        }
    }
    
    Constraint<> ub("Unique_Bought_Obsrvn"); //if sensor is sold, it is sold to exactly one agent; not sold - no agent uses it except owner
    ub = z0.in(z0_ids) + z.in(z_ids) - sn.in(ub_idx);
    model.add(ub.in(ub_idx) == 0);
    
    Constraint<> ub2("Unique_Bought_Obsrvn2"); //if sensor has no arcs, it cannot be sold
    ub2 = sn.in(ub_idx2);
    model.add(ub2.in(ub_idx2) == 0);

    
    indices z0_luo("z0_luo"), luo_idx("luo_idx");
    z0_luo = objects;
    luo_idx = objects;
    for (int i = 0; i < M; i++) {
        for (Arc* a: graph.get_node("Object_" + to_string(i))->get_in()) {
            z0_luo.add_in_row(i, a->_src->_name + "," + a->_dest->_name);
        }
        if (graph.get_node("Object_" + to_string(i))->get_in().size() > 0) {
            luo_idx.add_in_row(i, "Object_" + to_string(i));
        }
    }
    
    Constraint<> luo("Leader_Unique_Object"); //leader observes each obj no more than once
    luo = z0.in(z0_luo);//sum(z0.in_matrix(0, 1));
    model.add(luo.in(luo_idx) <= 1);
    
    Constraint<> lulb("Leader_Utility_lb"); //leader doesn't pay more than they get
    lulb = p.in_ignore_ith(1, 1, arcs) * z0.in(arcs) - w0.in(arcs);
    model.add(lulb.in(arcs) <= 0);

    indices fua_idx("fua_idx"), s_fua("s_fua");
    s_fua = own_arcs;
    fua_idx = sensors;
    for (int i = 0; i < N; i++) {
        for (Arc* a: graph.get_node("Sensor_" + to_string(i))->get_out()) {
            s_fua.add_in_row(i, a->_src->_name + "," + a->_dest->_name + ",Agent_" + to_string(owner[i]));
        }
        if (graph.get_node("Sensor_" + to_string(i))->get_out().size() > 0) {
            fua_idx.add_in_row(i, "Sensor_" + to_string(i));
        }
    }
    
    Constraint<> fua("Unique_Own_Assignment"); //sensor cannot do 2 or more observations
    fua = s.in(s_fua) + sn.in(fua_idx);
    model.add(fua.in(fua_idx) <= 1);
    
    /*matching indices for Follower_Unique_Object*/
    indices c_fub("c_fub"), z_fub("z_fub"), s_fub("s_fub");
    z_fub = *z._indices;
    s_fub = *s._indices;
    size_t row_id = 0;
    bool no_s = true, no_z = true;
    for (int j = 0; j < M; j++) {
        for (int k = 0; k < K; k++) {
            c_fub.insert("Object_" + to_string(j)+ ",Agent_" + to_string(k));
            no_s = true, no_z = true;
            for (Arc* a: graph.get_node("Object_" + to_string(j))->get_in()) {
                int i = stoi(a->_src->_name.substr(7, a->_src->_name.find(",")));
                if (k != owner[i]) {
                    no_z = false;
                    z_fub.add_in_row(row_id, a->_src->_name + "," + a->_dest->_name + ",Agent_" + to_string(k));
                }
                else{
                    no_s = false;
                    s_fub.add_in_row(row_id, a->_src->_name + "," + a->_dest->_name + ",Agent_" + to_string(k));
                }
            }
            if(no_z)
                z_fub.add_empty_row();
            if(no_s)
                s_fub.add_empty_row();
            row_id++;
        }
    }
    
    Constraint<> fub("Follower_Unique_Object"); //agents observe each object no more than once
    fub = s.in(s_fub) + z.in(z_fub);
    model.add(fub.in(c_fub) <= 1);
    
    Constraint<> fulb("Followers_Utility_lb"); //agents don't pay more than they get
    fulb = w_bought - p_z;
    model.add(fulb.in(bought_arcs) >= 0);
    
//        //----Fair price----
//    Constraint<> fp("FairPrice"); //set price to the mid-point btw buyer and seller
//    fp = p.in_ignore_ith(1, 2, bought_arcs) - (w_bought.in(bought_arcs) * z.in(bought_arcs) + y.in_ignore_ith(1, 2, bought_arcs))/2;
//    model.add(fp.in(bought_arcs) >= 0);
//
//    indices c_lb1("c_lb1"), y_lb1("y_lb1"), w_own_lb1("w_own_lb1"), w_own_z_lb1("w_own_z_lb1"), w_own_s_lb1("w_own_s_lb1"), z_lb1("z_lb1"), s_lb1("s_lb1");
//    y_lb1 = sensors;
//    z_lb1 = bought_arcs;
//    s_lb1 = own_arcs;
//    w_own_lb1 = own_arcs;
//    w_own_s_lb1 = own_arcs;
//    w_own_z_lb1 = own_arcs;
//    row_id = 0;
//    for (int i = 0; i < N; i++) {
//        for (Arc* b: graph.get_node("Sensor_" + to_string(i))->get_out()) {
//            no_z = true;
//            string j = b->_dest->_name;
//            c_lb1.add("Seller lb1:" + to_string(i) + "," + b->_dest->_name + ",Agent_" + to_string(owner[i]));
//            y_lb1.add_ref("Sensor_" + to_string(i));
//            w_own_lb1.add_ref("Sensor_" + to_string(i) + "," + b->_dest->_name + ",Agent_" + to_string(owner[i]));
//            for (Arc* a: graph.get_node(j)->get_in()) {
//                if (owner[stoi(a->_src->_name.substr(7, a->_src->_name.find(",")))] == owner[i]) {
//                    w_own_s_lb1.add_in_row(row_id, "Sensor_" + to_string(i) + "," + b->_dest->_name + ",Agent_" + to_string(owner[i]));
//                    s_lb1.add_in_row(row_id,a->_src->_name + "," + j +  ",Agent_" + to_string(owner[i]));
//                }
//                else if (owner[stoi(a->_src->_name.substr(6, a->_src->_name.find(",")))] != owner[i]) {
//                    w_own_z_lb1.add_in_row(row_id, "Sensor_" + to_string(i) + "," + b->_dest->_name + ",Agent_" + to_string(owner[i]));
//                    z_lb1.add_in_row(row_id, a->_src->_name + "," + j +  ",Agent_" + to_string(owner[i]));
//                    no_z = false;
//                }
//            }
//            if(no_z){
//                z_lb1.add_empty_row();
//            }
//            row_id++;
//        }
//    }
//
//    /*Price proxies (marginal utility of using own sensor)*/
//    Constraint<> sl1("Seller lb1"); //utility of adding this observation (if possible)
//    sl1 = y.in(y_lb1) + w_own.in(w_own_s_lb1)*s.in(s_lb1) + w_own.in(w_own_z_lb1)*z.in(z_lb1) - w_own.in(w_own_lb1);
//    model.add(sl1.in(c_lb1) >= 0);
//
//    Constraint<> sl2("Seller lb2"); //utility of replacing observation done by another sensor of their own by this one
//    sl2 = y.in_ignore_ith(1, 3, own_rplc) - (w_own.in_ignore_ith(1, 1, own_rplc) - w_own.in_ignore_ith(0, 1, own_rplc)) * s.in_ignore_ith(0, 1, own_rplc);
//    model.add(sl2.in(own_rplc) >= 0);
//
//    Constraint<> sl3("Seller lb3"); //utility of replacing observation done by a bought sensor by this one
//    sl3 = y.in_ignore_ith(1, 3, oths_rplc) - (w_own.in_ignore_ith(1, 1, oths_rplc) - w_bought.in_ignore_ith(0, 1, oths_rplc))*z.in_ignore_ith(0, 1, oths_rplc) - p_z.in_ignore_ith(0, 1, oths_rplc);
//    model.add(sl3.in(oths_rplc) >= 0);
    
    //For comparison: case with no collaboration (only using own sensors)
    /*Constraint<> no_colab("nc");
    no_colab = sn;
    model.add(no_colab.in(sensors) == 0);*/
}


void myModel::mSolve(bool run_mip) {

    model.set_name("Sensor_assign");
//    model.print();
//    model.print_constraints_stats(1e-6);
//    model.write();
//    model.write_solution();
    if(run_mip){
        solver<> sol(model, gurobi);
        sol.run();
//        model.write_solution();
    }
//#ifdef USE_H5CPP
////    auto hd5file = file::open(string(prj_dir)+"/data_sets/sensor/horizon10deg_0.hd5");
////    auto RootGroup = hd5file.root();
////    auto Dataset = RootGroup.get_dataset("masn");
////    dataspace::Simple Dataspace(Dataset.dataspace());
////    auto Dimensions = Dataspace.current_dimensions();
////    auto MaxDimensions = Dataspace.maximum_dimensions();
////    std::cout << "Dataset dimensions\n";
////    std::cout << "   Current | Max\n";
////    for (int i = 0; i < Dimensions.size(); i++) {
////        std::cout << "i:" << i << "      " << Dimensions[i] << " | "
////        << MaxDimensions[i] << "\n";
////    }
////
////    auto CreationProperties = Dataset.creation_list();
////    auto ChunkDims = CreationProperties.chunk();
////    std::cout << "\nChunk size\n";
////    for (int i = 0; i < ChunkDims.size(); i++) {
////        std::cout << "i:" << i << "     " << ChunkDims[i] << "\n";
////    }
////
////    std::cout << "\nData type\n";
////    auto Int32Type = datatype::create<std::int32_t>();
////    auto UInt32Type = datatype::create<std::uint32_t>();
////    auto FloatType = datatype::create<float>();
////    auto DataTypeClass = Dataset.datatype().get_class();
////    auto CurrentType = Dataset.datatype();
////    std::cout << "Is:        " << DataTypeClass << std::endl;
////    std::cout << "Is  int32: " << (Int32Type == CurrentType) << std::endl;
////    std::cout << "Is uint32: " << (UInt32Type == CurrentType) << std::endl;
////    std::cout << "Is  float: " << (FloatType == CurrentType) << std::endl;
////
////    std::cout << "\nAll elements\n";
////    std::vector<int> AllElements(Dataspace.size());
////    Dataset.read(AllElements);
////    for (auto Value : AllElements) {
////        std::cout << Value << " ";
////    }
////    std::cout << "\n\nRow access\n";
////    std::vector<int> RowData(static_cast<size_t>(Dimensions[1]));
////    for (size_t i = 0; i < Dimensions[0]; i++) {
////        dataspace::Hyperslab RowSelection{{i, 0}, {1, 3}};
////        Dataset.read(RowData, RowSelection);
////        std::cout << "i: " << i << " | ";
////        for (auto Value : RowData) {
////            std::cout << Value << " ";
////        }
////        std::cout << "\n";
////    }
////    std::cout << "\nElement access\n     j:0  j:1 j:2\n";
////    for (size_t i = 0; i < Dimensions[0]; i++) {
////        std::cout << "i:" << i << "    ";
////        for (size_t j = 0; j < Dimensions[1]; j++) {
////            int Value;
////            dataspace::Hyperslab ElementSelection{{i, j}, {1, 1}};
////            Dataset.read(Value, ElementSelection);
////            std::cout << Value << "    ";
////        }
////        std::cout << "\n";
////    }
//    // create a file
//    file::File f = file::create("sol.h5",file::AccessFlags::Truncate);
//
//    // create a group
//    node::Group root_group = f.root();
//    node::Group my_group = root_group.create_group("prices");
//
//    auto p = model.get_var<double>("p");
//    auto sn = model.get_var<int>("sn");
//    auto p_vals = *p.get_vals();
//    for (auto i = 0; i<p_vals.size(); i++) {
//        if(sn.eval(i)!=1)
//            p_vals[i] = 0;
//    }
//    // create a dataset
//    node::Dataset dataset = my_group.create_dataset("p",
//                                                    datatype::create<vector<double>>(),
//                                                    dataspace::create(p_vals));
//
//    // write to dataset
//    dataset.write(p_vals);
//    node::Group binary_group = root_group.create_group("assignment");
//
//    vector<string> binaries;
//    auto z = model.get_var<int>("z");
//    auto z_vals = *z.get_vals();
//    for (int i = 0; i<z_vals.size(); i++) {
//        if(z_vals[i]==1){
//            binaries.push_back(z._indices->get_key(i));
//        }
//
//    }
//    auto s = model.get_var<int>("s");
//    auto s_vals = *s.get_vals();
//    for (int i = 0; i<s_vals.size(); i++) {
//        if(s_vals[i]==1){
//            binaries.push_back(s._indices->get_key(i));
//        }
//
//    }
//
//    // create a dataset
//    dataset = binary_group.create_dataset("binaries: (sensor id, object id, agent id)",
//                                                    datatype::create<vector<string>>(),
//                                                    dataspace::create(binaries));
//
//    // write to dataset
//    dataset.write(binaries);
//#endif
}

class TCompare
{
public:
    bool operator() (const tuple<double, int, string>& w1, const tuple<double, int, string>& w2)
    {
        return get<0>(w1)>get<0>(w2);
    }
};

void myModel::GreedyStart(const param<double> &w0, const param<double> &w_own, const param<double> &w_bought) {

    DebugOn("Running Greedy Algorithm..."<<endl);
    /*copying params and putting them in a queue to change them in greedy*/
    param<double> wt0 = w0.deep_copy();
    param<double> wt_own = w_own.deep_copy();
    param<double> wt_bought = w_bought.deep_copy();

    priority_queue<tuple<double, int, string>, vector<tuple<double, int, string>>, TCompare> w_q; //tuples compared by first element
    double wt_sum = 0; //while loop stopping criterion

    for (int i = 0; i < w0.get_dim(); i++) {
        w_q.push(make_tuple(w0.eval(i), i, "w0"));
        wt_sum += w0.eval(i);
    }
    for (int i = 0; i < w_own.get_dim(); i++) {
        w_q.push(make_tuple(w_own.eval(i), i, "w_own"));
        wt_sum += w_own.eval(i);
    }
    for (int i = 0; i < w_bought.get_dim(); i++) {
        w_q.push(make_tuple(w_bought.eval(i), i, "w_bought"));
        wt_sum += w_bought.eval(i);
    }
    DebugOn("Done building priority queue " << endl);
    /*copying budget to change it in greedy*/
//    map<string,double> bdgt_map;
//    for (int i = 0; i<budget.size(); i++) {
//        bdgt_map[agents.get_key(i)] = budget[i];
//    }
    /*tmp; when no hdf5 and no bdgt*/
    map<string,double> bdgt_map;
    for (int k = 0; k<K; k++) {
        bdgt_map["Agent_" + to_string(k)] = 1000;
    }
    
    /*getting vars from model to set values*/
    auto s = model.get_var<int>("s");
    auto sn = model.get_var<int>("sn");
    auto z0 = model.get_var<int>("z0");
    auto z = model.get_var<int>("z");
    auto p = model.get_var<double>("p");
    auto p_sn = model.get_var<double>("p_sn");
    auto p_z = model.get_var<double>("p_z");
    auto y = model.get_var<double>("y");
    
    /*index and value of max weights*/
    string idx_str;
    string sensor_name;
    string object_name;
    string agent_name;
    tuple<double, int, string> max_wt;
    double w;
    int idx;

    //double obj = 0; //used to eval greedy objective
    map<string,double> p_ub; //price upper bound (for fair price)
    for (const auto &sensor_name:*sensors._keys) {
        p_ub[sensor_name] = std::max(w_own._range->second,w_bought._range->second);
    }

    double psum_wt0 = parSum(wt0), psum_wt_own = parSum(wt_own), psum_wt_bought = parSum(wt_bought);
    while(!w_q.empty()) {
//        DebugOn("parSum(wt0) = " << psum_wt0 << endl);
//        DebugOn("parSum(wt_own) = " << psum_wt_own << endl);
//        DebugOn("parSum(wt_bought) = " << psum_wt_bought << endl);

        /*finding max weights (idx) in 3 weight sets*/
        max_wt = w_q.top();
        w_q.pop();

        w = get<0>(max_wt);
        idx = get<1>(max_wt);

        if (get<2>(max_wt) == "w0") {
            if (wt0.eval(idx) > 0) {
                idx_str = w0.get_key(idx);
                sensor_name = idx_str.substr(0, idx_str.find(","));
                auto sub_idx = idx_str.substr(idx_str.find(",")+1);
                string object_name = sub_idx.substr(0, sub_idx.find(","));
                string agent_name  = sub_idx.substr(sub_idx.find(",")+1);
                int ownr = owner_map[sensor_name];
                string owner_name = agents.get_key(ownr);
                /*budget check might be too strict as it checks for (w_bought + w_own)/2, the actual price is likely to be lower (but cannot tell at this point)*/
                if (bdgt_map.at(agents.get_key(0)) >= (w + w_own.eval(sensor_name + "," + object_name + "," + owner_name))/2) {
                    /*assign to central authority*/
                    DebugOn("Sensor " + sensor_name + " assigned to Central Authority!\n");
                    assignLeader(idx_str, wt0, wt_own, wt_bought);
                    p_ub.at(sensor_name) = w; //price upper bound (for fair price)
                    bdgt_map.at(agents.get_key(0)) -= (w + wt_own.eval(sensor_name + "," + object_name + "," + owner_name))/2;
                }
                else{
                    wt0(sensor_name).set_val(0);
                }
            }
        }
        else if (get<2>(max_wt) == "w_own") {
            if (wt_own.eval(idx) > 0) {
                /*assign own*/
                idx_str = wt_own.get_key(idx);
                assignOwn(idx_str, wt0, wt_own, wt_bought);
            }
        }
        else {
            if (wt_bought.eval(idx) > 0) {
                idx_str = wt_bought.get_key(idx);
                sensor_name = idx_str.substr(0, idx_str.find(","));
                auto sub_idx = idx_str.substr(idx_str.find(",")+1);
                string object_name = sub_idx.substr(0, sub_idx.find(","));
                string agent_name  = sub_idx.substr(sub_idx.find(",")+1);
                //int ownr = owner_map[sensor_name];
                int ownr = owner[stoi(sensor_name.substr(sensor_name.find("_") + 1))];
                string owner_name = agents.get_key(ownr);
                if (bdgt_map.at(agent_name) >= (w + wt_own.eval(sensor_name + "," + object_name + "," + owner_name))/2) {
                    /*assign bought*/
                    assignBought(idx_str, wt0, wt_own, wt_bought, ownr, owner_name);
                    p_ub.at(sensor_name) = w; //price upper bound (for fair price)
                    bdgt_map.at(agent_name) -= (w + wt_own.eval(sensor_name + "," + object_name + "," + owner_name))/2;
                }
                else{
                    wt_bought(idx_str).set_val(0);
                }
            }
        }
//        psum_wt0 = parSum(wt0);
//        psum_wt_own = parSum(wt_own);
//        psum_wt_bought = parSum(wt_bought);
    }
   
    DebugOff("parSum(wt0) = " << psum_wt0 << endl);
    DebugOff("parSum(wt_own) = " << psum_wt_own << endl);
    DebugOff("parSum(wt_bought) = " << psum_wt_bought << endl);
    DebugOn("----------------------------------------" << endl);
    DebugOn("Computing fair prices..." << endl);

    /*set prices*/
    double y1 = 0; //Seller_lb1
    double y2 = 0; //Seller_lb2
    double y3 = 0; //Seller_lb3
    int s_sum = 0; //sum over own sens that this one could replace (for lb2)
    int z_sum = 0; //sum over bought sens one could replace (for lb3)
    string t; //var for replaced sensor
    for (int i = 0; i < N; i++) {
        auto sensor_name = sensors.get_key(i);
        if(sn.eval(i)!=1){
            DebugOn(sensor_name << " was not sold\n");
            continue;
        }
        DebugOn("Computing fair price for " << sensor_name << endl);
        auto sensor_node = graph.get_node(sensor_name);
        auto agent_name = agents.get_key(owner[i]);
        for (const Arc* a : sensor_node->get_out()) {
            for (const Arc* b : a->_dest->get_in()) {
                t = b->_src->_name;
                int o_id = owner[stoi(t.substr(t.find("_") + 1))]; //b->_src->owner_id;
                if (o_id == owner[i]) {
                    /*Seller_lb2*/
                    s_sum += s.eval(t + "," + a->_dest->_name + "," + agent_name);
                    y2 = std::max(y2, (w_own.eval(sensor_name + "," + a->_dest->_name + "," + agent_name) - w_own.eval(t + "," + a->_dest->_name + "," + agent_name)) * s.eval(t + "," + a->_dest->_name + "," + agent_name));
                }
                else {
                    /*Seller_lb3*/
                    string id_tmp = t + "," + a->_dest->_name + "," + agent_name;
                    z_sum += z.eval(id_tmp);
                    y3 = std::max(y3, std::min(p_ub.at(sensor_name), (w_own.eval(sensor_name + "," + a->_dest->_name + "," + agent_name) - w_bought.eval(t + "," + a->_dest->_name + "," + agent_name) + p_ub.at(t)) * z.eval(t + "," + a->_dest->_name + "," + agent_name)));
                }
            }
            if (y1 < w_own.eval(sensor_name + "," + a->_dest->_name + "," + agent_name) * (1 - s_sum - z_sum)) {
                /*Seller_lb1*/
                y1 = w_own.eval(sensor_name + "," + a->_dest->_name + "," + agent_name) * (1 - s_sum - z_sum);
            }
            s_sum = 0;
            z_sum = 0;
        }
        /*Fair price*/
        y.set_val(sensor_name,std::max(std::max(y1, y2), y3));
        p.set_val(sensor_name,(p_ub.at(sensor_name) + std::max(std::max(y1, y2), y3))/2.);
        /*update p_sn and p_z*/
//        p_sn.set_val(sensor_name,sn.eval(sensor_name) * (p_ub.at(sensor_name) + std::max(std::max(y1, y2), y3))/2.);
//        for (auto a : sensor_node->get_out()) {
//            for (int k = 0; k < K; k++) {
//                if(k==owner[i])
//                    continue;
//                string agent_k_name = agents.get_key(k);
//                p_z.set_val(sensor_name + "," + a->_dest->_name + "," + agent_k_name, z.eval(sensor_name + "," + a->_dest->_name + "," + agent_k_name) * (p_ub.at(sensor_name) + std::max(std::max(y1, y2), y3))/2.);
//            }
//        }
        //obj -= e * p_sn.eval("Sensor_" + to_string(i));
        y1 = 0;
        y2 = 0;
        y3 = 0;
    }
    //cout << "Greedy objective: " << obj << endl;
}

void myModel::writeGreedySol() {
    
    /*getting vars to eval them*/
    auto p = model.get_var<double>("p");
    auto y = model.get_var<double>("y");
    auto sn = model.get_var<int>("sn");
    auto z0 = model.get_var<int>("z0");
    auto s = model.get_var<int>("s");
    auto z = model.get_var<int>("z");
    
    //format: p y; id s; id z
    
    ofstream solFile;
    //solFile.open("sol.dat");
    solFile.open("/Users/svetlanariabova/Projects/Sensor/Data/sol_tmp/sol16.dat"); //hardcoded file name; might need to change later
    for (int i = 0; i < N; i++) {
        /*greedy solution might have vars at -inf when sensor is not used; replacing it with 0 in else statement*/
        if (p.eval("Sensor_" + to_string(i)) >= 0) {
            solFile << p.eval("Sensor_" + to_string(i)) << " " << y.eval("Sensor_" + to_string(i)) << " " << sn.eval("Sensor_" + to_string(i)) << endl;
        }
        else { solFile << 0 << " " << 0 << endl; }
    }
    
    bool no_use = true; //used to write a row if an arc is not used
    for (int i = 0; i < N; i++) {
        for (Arc* a: graph.get_node("Sensor_" + to_string(i))->get_out()) {
            if (s.eval("Sensor_" + to_string(i) + "," + a->_dest->_name + ",Agent_" + to_string(owner[i])) > 0.5) {
                solFile << "Sensor_" + to_string(i) + "," + a->_dest->_name << " " << 1 << endl; //arc used by sensor owner
            }
            else {
                for (int k = 0; k < owner[i]; k++) {
                    if (z.eval("Sensor_" + to_string(i) + "," + a->_dest->_name + ",Agent_" + to_string(k)) > 0.5) {
                        solFile << "Sensor_" + to_string(i) + "," + a->_dest->_name << " " << 1 << endl; //arc bought
                        no_use = false;
                        break;
                    }
                }
                for (int k = owner[i] + 1; k < K; k++) {
                    if (z.eval("Sensor_" + to_string(i) + "," + a->_dest->_name + ",Agent_" + to_string(k)) > 0.5) {
                        solFile << "Sensor_" + to_string(i) + "," + a->_dest->_name << " " << 1 << endl; //arc bought
                        no_use = false;
                        break;
                    }
                }
                if (no_use) { solFile << "Sensor_" + to_string(i) + "," + a->_dest->_name << " " << 0 << endl; } //arc not used
                no_use = true;
            }
        }
    }
    
    bool no_operate = true; //used to write a row when z(i, k) = 0
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < owner[i]; k++) {
            for (Arc* a: graph.get_node("Sensor_" + to_string(i))->get_out()) {
                if (z.eval("Sensor_" + to_string(i) + "," + a->_dest->_name + ",Agent_" + to_string(k)) > 0.5) {
                    solFile << "Sensor_" + to_string(i) + ",Agent_" + to_string(k) << " " << 1 << endl; //sensor operated by an agent who bought it
                    no_operate = false;
                    break;
                }
            }
            if (no_operate) { solFile << "Sensor_" + to_string(i) + ",Agent_" + to_string(k) << " " << 0 << endl; } //sensor not operated by that agent
            no_operate = true;
        }
        for (Arc* a: graph.get_node("Sensor_" + to_string(i))->get_out()) {
            if (s.eval("Sensor_" + to_string(i) + "," + a->_dest->_name + ",Agent_" + to_string(owner[i])) > 0.5) {
                solFile << "Sensor_" + to_string(i) + ",Agent_" + to_string(owner[i]) << " " << 1 << endl; //sensor operated by owner
                no_operate = false;
                break;
            }
        }
        if (no_operate) { solFile << "Sensor_" + to_string(i) + ",Agent_" + to_string(owner[i]) << " " << 0 << endl; } //sensor not operated by owner
        no_operate = true;
        
        for (int k = owner[i] + 1; k < K; k++) {
            for (Arc* a: graph.get_node("Sensor_" + to_string(i))->get_out()) {
                if (z.eval("Sensor_" + to_string(i) + "," + a->_dest->_name + ",Agent_" + to_string(k)) > 0.5) {
                    solFile << "Sensor_" + to_string(i) + ",Agent_" + to_string(k) << " " << 1 << endl; //sensor operated by an agent who bought it
                    no_operate = false;
                    break;
                }
            }
            if (no_operate) { solFile << "Sensor_" + to_string(i) + ",Agent_" + to_string(k) << " " << 0 << endl; } //sensor not operated by that agent
            no_operate = true;
        }
        for (Arc* a: graph.get_node("Sensor_" + to_string(i))->get_out()) {
            if (z0.eval("Sensor_" + to_string(i) + "," + a->_dest->_name) > 0.5) {
                solFile << "Sensor_" + to_string(i) << " " << 1 << endl; //sensor operated by CA
                no_operate = false;
                break;
            }
        }
        if (no_operate) { solFile << "Sensor_" + to_string(i) << " " << 0 << endl; } //sensor not operated by CA
        no_operate = true;
    }
    solFile.close();
}


void myModel::assignLeader(string &idx, param<double> &wt0, param<double> &wt_own, param<double> &wt_bought) {
    auto sn = model.get_var<int>("sn");
    auto z0 = model.get_var<int>("z0");
    int ownr;
    int sensor;
    int object;
    string j;
    
    ownr = owner[stoi(idx.substr(7, idx.find(",")))];
    sensor = stoi(idx.substr(7, idx.find(",")));
    object = stoi(idx.substr(idx.find("t") + 1, idx.find(",")));
    
    //assign to leader
    z0(idx.substr(0, nthOccurrence(idx, ",", 2))).set_val(1);

    //sensor not used for other objs leader + owner
    for (Arc* b: graph.get_node("Sensor_" + to_string(sensor))->get_out()) {
        j = b->_dest->_name;
        wt_own("Sensor_" + to_string(sensor) + "," + "Object_" + j + "," + "agent" + to_string(ownr)).set_val(0);
        wt0("Sensor_" + to_string(sensor) + "," + "Object_" + j).set_val(0);
    }
    //object not observed twice
    for (int i = 0; i < sensor; i ++) {
        wt0("Sensor_" + to_string(i) + "," + "Object_" + to_string(object)).set_val(0);
    }
    for (int i = sensor + 1; i < N; i ++) {
        wt0("Sensor_" + to_string(i) + "," + "Object_" + to_string(object)).set_val(0);
    }
    //sensor not used by other agents
    for (int k = 0; k < ownr; k++) {
        for (Arc* b: graph.get_node("Sensor_" + to_string(sensor))->get_out()) {
            j = b->_dest->_name;
            wt_bought("Sensor_" + to_string(sensor) + "," + "Object_" + j + "," + "agent" + to_string(k)).set_val(0);
        }
    }
    for (int k = ownr + 1; k < K; k++) {
        for (Arc* b: graph.get_node("Sensor_" + to_string(sensor))->get_out()) {
            j = b->_dest->_name;
            wt_bought("Sensor_" + to_string(sensor) + "," + "Object_" + j + "," + "agent" + to_string(k)).set_val(0);
        }
    }
    //sold
    sn("Sensor_" + to_string(sensor)).set_val(1);

}

void myModel::assignOwn(string &idx, param<double> &wt0, param<double> &wt_own, param<double> &wt_bought) {
    auto s = model.get_var<int>("s");
//    int ownr;
//    int sensor;
//    string object;
    string j;
    string i;
    
    string sensor_name = idx.substr(0, idx.find(","));
    auto sub_idx = idx.substr(idx.find(",")+1);
    string object_name = sub_idx.substr(0, sub_idx.find(","));
    string agent_name  = sub_idx.substr(sub_idx.find(",")+1);
    int ownr = owner_map[sensor_name];
    string owner_name = agents.get_key(ownr);

    
    //assign to owner
    s.set_val(idx,1);
    auto sensor_node = graph.get_node(sensor_name);
    if(!sensor_node){
        cerr << "Cannot find node with name " << sensor_name << endl;
        exit(-1);
    }
    //sensor not used by leader
    for (Arc* b: sensor_node->get_out()) {
        j = b->_dest->_name;
        wt0.set_val((sensor_name + "," + j),0);
        //sensor not used twice
        wt_own.set_val((sensor_name + "," + j + "," + owner_name),0);
        //sensor not used by other agents
        for (int k = 0; k < K; k++) {
            if(k==ownr)
                continue;
            string agent_k_name = agents.get_key(k);
            wt_bought.set_val((sensor_name + "," + j + "," + agent_k_name),0);
        }
    }
    auto object_node = graph.get_node(object_name);
    if(!object_node){
        cerr << "In assignOwn, cannot find node with name " << object_name << endl;
        cerr << "idx = " << idx << endl;
        exit(-1);
    }
//    object not observed twice
    for (Arc* a: object_node->get_in()) {
        i = a->_src->_name;
        if (a->_src->owner_id == ownr) { wt_own.set_val((i + "," + object_name + "," + owner_name),0); }
        else { wt_bought.set_val((i + "," + object_name + "," + owner_name),0); }
    }
}

void myModel::assignBought(string &idx, param<double> &wt0, param<double> &wt_own, param<double> &wt_bought, int ownr, string owner_name) {
    auto sn = model.get_var<int>("sn");
    auto z = model.get_var<int>("z");
    //int ownr;
    string i,j;
    string sensor_name = idx.substr(0, idx.find(","));
    auto sub_idx = idx.substr(idx.find(",")+1);
    string object_name = sub_idx.substr(0, sub_idx.find(","));
    string agent_name  = sub_idx.substr(sub_idx.find(",")+1);
    //ownr = owner_map[sensor_name];
    //string owner_name = agents.get_key(ownr);
    //assign to agent
    z.set_val(idx,1);

    auto sensor_node = graph.get_node(sensor_name);
    if(!sensor_node){
        cerr << "Cannot find node with name " << sensor_name << endl;
        exit(-1);
    }
        
    //sensor not used for other objs leader + owner + buyer
    for (Arc* b: sensor_node->get_out()) {
        j = b->_dest->_name;
        wt_own.set_val((sensor_name + "," + j + "," + owner_name),0);
        wt0.set_val((sensor_name + "," + j),0);
//        if(agent_name!=owner_name)
        wt_bought.set_val((sensor_name + "," + j + "," + agent_name),0);
    }
//    object not observed twice
    auto object_node = graph.get_node(object_name);
    if(!object_node){
        cerr << "In assignBought, cannot find node with name " << object_name << endl;
        cerr << "idx = " << idx << endl;
        exit(-1);
    }
    for (Arc* a: object_node->get_in()) {
        i = a->_src->_name;
        if (agents.get_key(a->_src->owner_id) == agent_name) { wt_own.set_val(i + "," + object_name + "," + agent_name,0); }
        else { wt_bought.set_val(i + "," + object_name + "," + agent_name,0); }
    }
//    sensor not used by other agents
    for (Arc* b: sensor_node->get_out()) {
        j = b->_dest->_name;
        for (int k = 0; k < K; k++) {
            string it_agent_name = agents.get_key(k);
            if(k==ownr || it_agent_name==agent_name)
                continue;
            wt_bought.set_val((sensor_name + "," + j + "," + it_agent_name),0);
        }
    }
    //sold
    sn.set_val((sensor_name),1);
}

double myModel::parSum(const param<double>& w) {
    double s = 0;
    for (auto& n : *w.get_vals()) {
        s += n;
    }
    return s;
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

void myModel::InitBilevel2(param<double> &w0, param<double> &w_own, param<double> &w_bought, double eps) {
    
    e = eps; //for e * p_sn in obj
    
    /*reset_range to define ub on prices*/
    w_own.reset_range();
    w_bought.reset_range();
    w0.reset_range();
    
    /*---------Variables----------*/
    /*price + aux for price*/
    var<double> p("p", 0, std::max(w_own._range->second,w_bought._range->second));
    model.add(p.in(sensors));
    var<double> y("y", 0, std::max(w_own._range->second,w_bought._range->second));
    model.add(y.in(sensors));
    
    /*assignment*/
    var<int> s("s", 0, 1); //use own sens
    model.add(s.in(own_arcs));
    var<int> sn("sn", 0, 1); //sensor sold
    model.add(sn.in(sensors));
    var<int> z0("z0", 0, 1); //leader buys sens (observation)
    model.add(z0.in(arcs));
    var<int> z("z", 0, 1); //agebts buy sens (observation)
    model.add(z.in(bought_arcs));
    
    /*duals*/
    var<double> a("a", 0, 100); //CHANGE UB!!!
    model.add(a.in(own_sens));
//    var<double> b("b", 0, 100);
//    model.add(b.in(bought_sens));
    var<double> g("g", 0, 1);
    model.add(g.in(objects));
    
    /*multiplication vars*/
    var<double> p_sn("p_sn", 0, std::max(w_own._range->second,w_bought._range->second));
    model.add(p_sn.in(sensors));
    var<double> p_z("p_z", 0, std::max(w_own._range->second,w_bought._range->second));
    model.add(p_z.in(bought_arcs));
    /* experiimental */
    var<double> p_z0("p_z0", 0, std::max(w_own._range->second,w0._range->second));
    model.add(p_z0.in(arcs));
    
    /*multiplication vars def*/
    Constraint<> p_sn_def("p_sn_def");
    p_sn_def = p*sn - p_sn;
    model.add(p_sn_def.in(sensors)==0);
    
    Constraint<> p_z_def("p_z_def");
    p_z_def = p.in_ignore_ith(1, 2, bought_arcs)*z - p_z;
    model.add(p_z_def.in(bought_arcs)==0);
    
    /* p*z0 for cancelling out price terms in the obj (experimental) */
    Constraint<> p_z0_def("p_z0_def");
    p_z0_def = p.in_ignore_ith(1, 1, arcs) * z0 - p_z0;
    model.add(p_z0_def.in(arcs) == 0);
    
    
    /*weights for obj (own + leader, bought + leader); leader receives utility when an object is observed, no matter by whom*/
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
    
    /*--------Objective---------*/
    func<> obj;
    obj += product(w_own0, s); //use own sens
    obj += sum(p_sn); //sell sens
    obj += product(w_bought0, z); //buy sens pt.1
    obj += product(w0, z0);
    obj -= sum(p_z); //buy sens pt.2
    obj -= e * sum(p_sn); //regularization term; sets prices to their lb from Fair_price constraints (can't use equality there)
    model.max(obj);
    
    /*--------Constraints---------*/
    
    /*matching index sets for z0 and z in Unique_Bought_Obsrvn*/
    indices z0_ids("z0_ids"), z_ids("z_ids"), ub_idx("ub_idx"), ub_idx2("ub_idx2");
    z0_ids = arcs;
    z_ids = bought_arcs;
    ub_idx = sensors;
    ub_idx2 = sensors;

    for (int i = 0; i<N; i++) {
        for (Arc* b: graph.get_node("Sensor_" + to_string(i))->get_out()) {
            string j = b->_dest->_name;
            z0_ids.add_in_row(i, "Sensor_" + to_string(i) + "," + b->_dest->_name);
            for (int k = 0; k < K; k++) {
                if (k != owner[i]) {
                    z_ids.add_in_row(i,"Sensor_" + to_string(i) + "," + b->_dest->_name + ",Agent_" + to_string(k));
                }
            }
        }
        if (graph.get_node("Sensor_" + to_string(i))->get_out().size() > 0) {
            ub_idx.add_in_row(i, "Sensor_" + to_string(i));
        }
        else {
            ub_idx2.add_in_row(i, "Sensor_" + to_string(i));
        }
    }
    
    Constraint<> ub("Unique_Bought_Obsrvn"); //if sensor is sold, it is sold to exactly one agent; not sold - no agent uses it except owner
    ub = z0.in(z0_ids) + z.in(z_ids) - sn.in(ub_idx);
    model.add(ub.in(ub_idx) == 0);
    
    Constraint<> ub2("Unique_Bought_Obsrvn2"); //if sensor has no arcs, it cannot be sold
    ub2 = sn.in(ub_idx2);
    model.add(ub2.in(ub_idx2) == 0);
    
    indices z0_luo("z0_luo"), luo_idx("luo_idx");
    z0_luo = objects;
    luo_idx = objects;
    for (int i = 0; i < M; i++) {
        for (Arc* a: graph.get_node("Object_" + to_string(i))->get_in()) {
            z0_luo.add_in_row(i, a->_src->_name + "," + a->_dest->_name);
        }
        if (graph.get_node("Object_" + to_string(i))->get_in().size() > 0) {
            luo_idx.add_in_row(i, "Object_" + to_string(i));
        }
    }
    
    Constraint<> luo("Leader_Unique_Object"); //leader observes each obj no more than once
    luo = z0.in(z0_luo);//sum(z0.in_matrix(0, 1));
    model.add(luo.in(luo_idx) <= 1);
    
    Constraint<> lulb("Leader_Utility_lb"); //leader doesn't pay more than they get
    lulb = p.in_ignore_ith(1, 1, arcs) * z0.in(arcs) - w0.in(arcs);
    model.add(lulb.in(arcs) <= 0);
    
    /*duality*/
    
    indices s_sd("s_sd"), z_sd("z_sd"), sn_sd("sn_sd"), a_sd("a_sd"), b_sd("b_sd"), g_sd("g_sd"), sd_idx("sd_idx"), s_sd2("s_sd2"), z_sd2("z_sd2"), sn_sd2("sn_sd2"), a_sd2("a_sd2"), b_sd2("b_sd2"), g_sd2("g_sd2"), sd_idx2("sd_idx2");
    s_sd = own_arcs;
    z_sd = bought_arcs;
    sn_sd = sensors;
    a_sd = own_sens;
    b_sd = bought_sens;
    g_sd = objects;
    sd_idx = agents;
    sd_idx2 = agents;
    string idx, sensor_name, sub_idx, object_name, agent_name;
    
    for (int k = 0; k < K; k++) {
        int own_arc_count = 0;
        for (int i = 0; i < own_arcs.size(); i++) {
            idx = own_arcs.get_key(i);
            sensor_name = idx.substr(0, idx.find(","));
            auto sub_idx = idx.substr(idx.find(",")+1);
            string object_name = sub_idx.substr(0, sub_idx.find(","));
            string agent_name  = sub_idx.substr(sub_idx.find(",")+1);
            if (agent_name == "Agent_" + to_string(k)) {
                s_sd.add_in_row(k, idx);
                own_arc_count += 1;
            }
        }
        if (own_arc_count > 0) {
            sd_idx.add_in_row(k, "Agent_" + to_string(k));
            for (int i = 0; i < bought_arcs.size(); i++) {
                idx = bought_arcs.get_key(i);
                sensor_name = idx.substr(0, idx.find(","));
                auto sub_idx = idx.substr(idx.find(",")+1);
                string object_name = sub_idx.substr(0, sub_idx.find(","));
                string agent_name  = sub_idx.substr(sub_idx.find(",")+1);
                if (agent_name == "Agent_" + to_string(k)) {
                    z_sd.add_in_row(k, idx);
                }
            }
            for (int i = 0; i < own_sens.size(); i++) {
                idx = own_sens.get_key(i);
                sensor_name = idx.substr(0, idx.find(","));
                auto sub_idx = idx.substr(idx.find(",")+1);
                string agent_name = sub_idx.substr(0, sub_idx.find(","));
                if (agent_name == "Agent_" + to_string(k)) {
                    sn_sd.add_in_row(k, sensor_name);
                    a_sd.add_in_row(k, idx);
                }
            }
        }
        else {
            sd_idx2.add_in_row(k, "Agent_" + to_string(k));
            for (int i = 0; i < bought_arcs.size(); i++) {
                idx = bought_arcs.get_key(i);
                sensor_name = idx.substr(0, idx.find(","));
                auto sub_idx = idx.substr(idx.find(",")+1);
                string object_name = sub_idx.substr(0, sub_idx.find(","));
                string agent_name  = sub_idx.substr(sub_idx.find(",")+1);
                if (agent_name == "Agent_" + to_string(k)) {
                    z_sd2.add_in_row(k, idx);
                }
            }
            for (int i = 0; i < own_sens.size(); i++) {
                idx = own_sens.get_key(i);
                sensor_name = idx.substr(0, idx.find(","));
                auto sub_idx = idx.substr(idx.find(",")+1);
                string agent_name = sub_idx.substr(0, sub_idx.find(","));
                if (agent_name == "Agent_" + to_string(k)) {
                    sn_sd2.add_in_row(k, sensor_name);
                    a_sd2.add_in_row(k, idx);
                }
            }
        }
        
    }
    for (int i = 0; i < M; i++) {
        if (graph.get_node("Object_" + to_string(i))->get_in().size() > 0) {
            g_sd.add_in_row(i, "Object_" + to_string(i));
        }
    }
    
    Constraint<> sd1("Strong_Duality_wOwnArcs"); //primal obj = dual obj
    sd1 = product(w_own.in(s_sd), s.in(s_sd)) + p_sn.in(sn_sd) + product(w_bought.in(z_sd), z.in(z_sd)) - p_z.in(z_ids) - a.in(a_sd) - sum(g.in(g_sd));
    model.add(sd1.in(sd_idx) == 0);

    Constraint<> sd2("Strong_Duality_woOwnArcs"); //primal obj = dual obj
    sd2 = product(w_bought.in(z_sd2), z.in(z_sd2)) - p_z.in(z_sd2) - a.in(a_sd2) - sum(g.in(g_sd2));
    model.add(sd2.in(sd_idx2) == 0);

    indices fua_idx("fua_idx"), s_fua("s_fua");
    s_fua = own_arcs;
    fua_idx = sensors;
    for (int i = 0; i < N; i++) {
        for (Arc* a: graph.get_node("Sensor_" + to_string(i))->get_out()) {
            s_fua.add_in_row(i, a->_src->_name + "," + a->_dest->_name + ",Agent_" + to_string(owner[i]));
        }
        if (graph.get_node("Sensor_" + to_string(i))->get_out().size() > 0) {
            fua_idx.add_in_row(i, "Sensor_" + to_string(i));
        }
    }
    
    Constraint<> fua("Unique_Own_Assignment"); //sensor cannot do 2 or more observations
    fua = s.in(s_fua) + sn.in(fua_idx);
    model.add(fua.in(fua_idx) <= 1);
    
    /*matching indices for Follower_Unique_Object*/
    indices c_fub("c_fub"), z_fub("z_fub"), s_fub("s_fub");
    z_fub = *z._indices;
    s_fub = *s._indices;
    size_t row_id = 0;
    bool no_s = true, no_z = true;
    for (int j = 0; j < M; j++) {
        for (int k = 0; k < K; k++) {
            c_fub.insert("Object_" + to_string(j)+ ",Agent_" + to_string(k));
            no_s = true, no_z = true;
            for (Arc* a: graph.get_node("Object_" + to_string(j))->get_in()) {
                int i = stoi(a->_src->_name.substr(7, a->_src->_name.find(",")));
                if (k != owner[i]) {
                    no_z = false;
                    z_fub.add_in_row(row_id, a->_src->_name + "," + a->_dest->_name + ",Agent_" + to_string(k));
                }
                else{
                    no_s = false;
                    s_fub.add_in_row(row_id, a->_src->_name + "," + a->_dest->_name + ",Agent_" + to_string(k));
                }
            }
            if(no_z)
                z_fub.add_empty_row();
            if(no_s)
                s_fub.add_empty_row();
            row_id++;
        }
    }
    
    Constraint<> fub("Follower_Unique_Object"); //agents observe each object no more than once
    fub = s.in(s_fub) + z.in(z_fub);
    model.add(fub.in(c_fub) <= 1);
    
    indices g_df1("g_df1");
    g_df1 = objects;
    for (int i = 0; i < own_arcs.size(); i++) {
        idx = own_arcs.get_key(i);
        sensor_name = idx.substr(0, idx.find(","));
        auto sub_idx = idx.substr(idx.find(",")+1);
        string object_name = sub_idx.substr(0, sub_idx.find(","));
        string agent_name  = sub_idx.substr(sub_idx.find(",")+1);
        g_df1.add_in_row(i, object_name);
    }
    
    Constraint<> df1("Dual_feas1");
    df1 = a.in_ignore_ith(1, 1, own_arcs) + g.in(g_df1) - w_own.in(own_arcs);
    model.add(df1.in(own_arcs) >= 0);
    
    Constraint<> df2("Dual_feas2");
    df2 = a.in(own_sens) - p.in_ignore_ith(1, 1, own_sens);
    model.add(df2.in(own_sens) >= 0);
    
    indices g_df3("g_df3");
    g_df3 = objects;
    for (int i = 0; i < bought_arcs.size(); i++) {
        idx = bought_arcs.get_key(i);
        sensor_name = idx.substr(0, idx.find(","));
        auto sub_idx = idx.substr(idx.find(",")+1);
        string object_name = sub_idx.substr(0, sub_idx.find(","));
        string agent_name  = sub_idx.substr(sub_idx.find(",")+1);
        g_df3.add_in_row(i, object_name);
    }
    
    Constraint<> df3("Dual_feas3");
    df3 = g.in(g_df3) - w_bought.in(bought_arcs) + p.in_ignore_ith(1, 2, bought_arcs);
    model.add(df3.in(bought_arcs) >= 0);
    
   
//    Constraint<> fulb("Followers_Utility_lb"); //agents don't pay more than they get
//    fulb = w_bought - p_z;
//    model.add(fulb.in(bought_arcs) >= 0);
    
}
