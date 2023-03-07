//
//  misdp.cpp
//  misdp
//
//  Created by Smitha on 7/7/22.
//

#include <stdio.h>
#include "read_misdp.h"
#include <gravity/solver.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <fstream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;


vector<vector<int>> allPossibleSubset(int n)
{
    vector<vector<int>> res;
    int count = pow(2, n);
    // The outer for loop will run 2^n times to print all subset .
    // Here variable i will act as a binary counter
    for (int i = 0; i < count; i++) {
        vector<int> v;
        // The inner for loop will run n times ,
        // As the maximum number of elements a set can have is n
        // This loop will generate a subset
        for (int j = 0; j < n; j++) {
            // This if condition will check if jth bit in
            // binary representation of  i  is set or not
            // if the value of (i & (1 << j)) is not 0 ,
            // include j in the current subset
            // otherwise exclude j
            if ((i & (1 << j)) != 0)
                v.push_back(j);
        }
        if(v.size()>=3)
            res.push_back(v);
    }
    return res;
}


/*returns matrix of dimension n-1, without row rowno and col colno*/
vector<vector<double>> get_minor(vector<vector<double>>& X, int rowno, int colno, int n){
    if(n==1)
        throw invalid_argument ("n=1");
    vector<vector<double>> Y(n-1);
    int count=0;
    for(auto i=0; i<n;i++){
        for(auto j=0; j<n && i!=rowno;j++){
            if(j!=colno)
                Y[count].push_back(X[i][j]);
        }
        if(i!=rowno){
            count++;
        }
    }
    return Y;
}

double determinant(vector<vector<double>> X, int n){
    double det;
    if(n==1){
        det=X[0][0];
    }
    else{
        det=0;
        for(auto i=0;i<n;i++){
            auto Y=get_minor(X, 0, i, n);
            det+=std::pow(-1, i)*X[0][i]*determinant(Y,n-1);
        }
    }
    return det;
}
using namespace gravity;
using namespace std;

int main(int argc, char * argv[]){
    //string fname=string(prj_dir)+"/data_sets/MISDP/2x3_3bars.cbf";
string fname=string(prj_dir)+"/data_sets/MISDP/2x7_3bars.cbf";
    //string fname=string(prj_dir)+"/data_sets/MISDP/2x4_2scen_3bars.cbf";
   // string fname=string(prj_dir)+"/data_sets/MISDP/8_1.json";
    //string fname=string(prj_dir)+"/data_sets/MISDP/coloncancer_1_100_5.cbf";
    //string fname=string(prj_dir)+"/data_sets/MISDP/2g_4_164_k3_5_6.cbf";
    
//    if(argc>=2){
//        fname=argv[1];
//    }
//
//    std::ifstream i(fname);
//    json e;
//    i >> e;
//    vector<pair<pair<int,int>,double>> edges_to_augment;
//    int num_nodes = e["num_nodes"];
//    DebugOn("Number of nodes = " << num_nodes << endl);
//    edges_to_augment.insert(end(edges_to_augment), begin(e["edges_to_augment"]), end(e["edges_to_augment"]));
//    DebugOn("Number of edges = " << edges_to_augment.size() << endl);
//
//    /* Graph */
//    Net graph;
//    for (int i = 0; i < num_nodes; i++) {
//        auto n = new Node(to_string(i));
//        graph.add_node(n);
//    }
//    int index = 0;
//    for (const auto &e:edges_to_augment) {
//        string src = to_string(e.first.first-1);
//        string dest = to_string(e.first.second-1);
//        double weight = e.second;
//        string name = src+","+dest;
//        auto arc = new Arc(name);
//        arc->_id = index++;
//        arc->_src = graph.get_node(src);
//        arc->_dest = graph.get_node(dest);
//        arc->_weight = weight;
//        arc->connect();
//        graph.add_arc(arc);
//    }
//
//    /* Indices */
//    indices E = indices(graph.arcs);
//    indices N = indices(graph.nodes);
//    indices S("S"), E_S("E_S");
//    S = range(1,std::pow(2,num_nodes) - (num_nodes*(num_nodes-1))/2 - num_nodes - 1);/* excluding empty S, |S| = 1, and |S| = 2*/
//    param<> rhs("rhs");
//    rhs.in(S);
//    E_S = E;
//    int row_id = 0, count = 0;
//    vector<pair<int,std::vector<string>>> _bag_names;
//    auto subsets = allPossibleSubset(num_nodes);
//    for (auto const &subset: subsets) {
//        pair<int,vector<string>> bn;
//        bn.first=count++;
//        for(auto n:subset){
//            bn.second.push_back(to_string(n));
//            DebugOn(n <<"\t");
//        }
//        DebugOn(endl);
//        _bag_names.push_back(bn);
//        E_S.add_empty_row();
//        rhs.set_val(row_id, subset.size()-1);
//        for (int i = 0; i<subset.size()-1; i++) {
//            for (int j = i+1; j<subset.size(); j++) {
//                E_S.add_in_row(row_id, to_string(subset[i])+","+to_string(subset[j]));
//            }
//        }
//        row_id++;
//    }
//    indices delta_i("delta_i");
//    delta_i = E;
//    for (int i = 0; i<num_nodes; i++) {
//        auto node = graph.nodes.at(i);
//        delta_i.add_empty_row();
//        for (const Arc* a: node->get_out()) {
//            delta_i.add_in_row(i, a->_name);
//        }
//        for (const Arc* a: node->get_in()) {
//            delta_i.add_in_row(i, a->_name);
//        }
//    }
//
//    indices delta_i_excl0("delta_i_excl0"), delta_i_excln("delta_i_excln"), delta_0("delta_0");
//    delta_i_excl0 = E;
//    delta_i_excln = E;
//    delta_0 = E;
//    auto node_0 = graph.nodes.at(0);
//    for (const Arc* a: node_0->get_out()) {
//        delta_0.add_in_row(0, a->_name);
//    }
//    for (const Arc* a: node_0->get_in()) {
//        delta_0.add_in_row(0, a->_name);
//    }
//    for (int i = 0; i<num_nodes-1; i++) {
//        auto node_0 = graph.nodes.at(i);
//        auto node_1 = graph.nodes.at(i+1);
//        delta_i_excl0.add_empty_row();
//        delta_i_excln.add_empty_row();
//        for (const Arc* a: node_0->get_out()) {
//            delta_i_excln.add_in_row(i, a->_name);
//        }
//        for (const Arc* a: node_0->get_in()) {
//            delta_i_excln.add_in_row(i, a->_name);
//        }
//        for (const Arc* a: node_1->get_out()) {
//            delta_i_excl0.add_in_row(i, a->_name);
//        }
//        for (const Arc* a: node_1->get_in()) {
//            delta_i_excl0.add_in_row(i, a->_name);
//        }
//    }
//
//    /* Parameters*/
//    param<> w("w");
//    w.in(E);
//    index = 0;
//    for(const Arc* a: graph.arcs){
//        w.set_val(index++, a->_weight);
//    }
//
//    bool project = false;
//    /* Model */
//    Model<> Laplacian("Laplacian");
////    Laplacian._bag_names = _bag_names;
////    DebugOn("Adding " << _bag_names.size() << " 3x3 bags\n");
//    /* Variables */    var<> gamma("𝛾");
//    if(!project)
//        Laplacian.add(gamma.in(R(1)));
//    var<int> x("x", 0, 1);
//    Laplacian.add(x.in(E));
//    var<> Wii("Wii", pos_);
//    Laplacian.add(Wii.in(N));
//    var<> Wij("Wij");
//    Laplacian.add(Wij.in(E));
//
//
//    /* Constraints */
//    if(project){
//        indices E_excl0("E_excl0"), E_excln("E_excln");
//        for(const Arc* a: graph.arcs){
//            if(a!=graph.arcs.at(0))
//                E_excl0.add(a->_name);
//            if(a!=graph.arcs.back())
//                E_excln.add(a->_name);
//        }
//        Constraint<> Wij_def("Wij_def");
//        Wij_def = Wij.in(E_excln) + w.in(E_excln)*x.in(E_excln) - (Wij.in(E_excl0) + w.in(E_excl0)*x.in(E_excl0));
//        Laplacian.add(Wij_def.in(E_excl0) == 0);
//
//        indices N_excl0("N_excl0"), N_excln("N_excln");
//        for(const Node* n: graph.nodes){
//            if(n!=graph.nodes.at(0))
//                N_excl0.add(n->_name);
//            if(n!=graph.nodes.back())
//                N_excln.add(n->_name);
//        }
//
//        Constraint<> Wii_def("Wii_def");
//        Wii_def = w.in(delta_i_excln)*x.in(delta_i_excln) - Wii.in(N_excln) - (w.in(delta_i_excl0)*x.in(delta_i_excl0) - Wii.in(N_excl0));
//        Laplacian.add(Wii_def.in(N_excl0) == 0);
//
//        Constraint<> Wii_W_ij_def("Wii_W_ij_def");
//        Wii_W_ij_def = (num_nodes/(num_nodes-1.))*(w.in(delta_0)*x.in(delta_0) - Wii(N._keys->front())) - num_nodes*(Wij(E._keys->front()) + w(E._keys->front())*x(E._keys->front()));
//        Laplacian.add(Wii_W_ij_def == 0);
//    }
//    else{
//        Constraint<> Wij_def("Wij_def");
//        Wij_def = Wij + w*x - (1./num_nodes)*gamma;
//        Laplacian.add(Wij_def.in(E) == 0);
//
//        Constraint<> Wii_def("Wii_def");
//        Wii_def = Wii - (w.in(delta_i)*x.in(delta_i) - (num_nodes-1.)/num_nodes*gamma);
//        Laplacian.add(Wii_def.in(N) == 0);
//    }
//
//    Constraint<> NodeCut("NodeCut");
//    NodeCut = x.in(delta_i);
//    Laplacian.add(NodeCut.in(N) >= 1);
//
//    Constraint<> Spanning_tree("Spanning_tree");
//    Spanning_tree = sum(x);
//    Laplacian.add(Spanning_tree == num_nodes - 1);
//
//    Constraint<> Subtour("Subtour");
//    Subtour = x.in(E_S) - rhs;
//    Laplacian.add(Subtour.in(S) <= 0);
//
//    Constraint<> Minor2("Minor2");
//    Minor2 = Wij*Wij - Wii.from(E)*Wii.to(E);
//    Laplacian.add(Minor2.in(E) <= 0);
//
//    Laplacian.make_PSD(Wii,Wij);
////    Laplacian.max(gamma);
//
////    double coef_ij = num_nodes/(2.*E.size());
////    Laplacian.max(coef*(sum(Wij) + product(w,x)));
//    if(project){
//        double coef_i = 1./((num_nodes - 1));
//        Laplacian.max(coef_i*(2.*product(w,x) - sum(Wii)));// + coef_ij*(sum(Wij) + product(w,x)));
//        //    Laplacian.max(coef_i*(w.in(delta_0)*x.in(delta_0) - Wii(N._keys->front())));// + coef_ij*(sum(Wij) + product(w,x)));
//    }
//    else{
//        Laplacian.max(gamma);
//        /* Map linking matrix entries to decision vars */
//        vector<vector<pair<pair<size_t,size_t>, double>>> Wij_gamma_map(E.size());
//        vector<vector<pair<pair<size_t,size_t>, double>>> Wii_gamma_map(N.size());
//        vector<vector<pair<pair<size_t,size_t>, double>>> Wij_x_map(E.size());
//        vector<vector<pair<pair<size_t,size_t>, double>>> Wii_x_map(N.size());
//
//
//
//        for (int i = 0; i<num_nodes; i++) {
//            auto node = graph.nodes.at(i);
//            Wii_gamma_map[i].push_back({{gamma.get_vec_id(),0},-(num_nodes-1.)/num_nodes});
//            for (const Arc* a: node->get_out()) {
//                Wii_x_map[i].push_back({{x.get_vec_id(),a->_id},a->_weight});
//            }
//            for (const Arc* a: node->get_in()) {
//                Wii_x_map[i].push_back({{x.get_vec_id(),a->_id},a->_weight});
//            }
//        }
//
//        index = 0;
//        for(const Arc* a: graph.arcs){
//            Wij_gamma_map[index].push_back({{gamma.get_vec_id(),0},1./num_nodes});
//            Wij_x_map[index].push_back({{x.get_vec_id(),a->_id},-a->_weight});
//            index++;
//        }
//        bool orig_cuts = true;/* Generate the SDP cuts using original variables (get rid of Wii and Wij) */
//        if(!project && orig_cuts){
//            Laplacian.Xii_x_map = Wii_gamma_map;
//            Laplacian.Xij_x_map = Wij_gamma_map;
//            Laplacian.Xii_y_map = Wii_x_map;
//            Laplacian.Xij_y_map = Wij_x_map;
//        }
//    }
//
//    Laplacian.print();
////    return 0;
//    solver<> sv(Laplacian,cplex);
//    sv.run();
//    Laplacian.print_solution();
//    return 0;
    
    
    
    bool root_refine = false;
    if(argc>=2){
        fname=argv[1];
    }
    if(argc>=3){
        root_refine = true;
    }
    auto m=make_shared<Model<double>>("misdp_test");
    
    auto g=CBF_read(fname.c_str(), m);
    m->print();
    double ub_solver_tol=1e-6, lb_solver_tol=1e-6, range_tol=1e-4, opt_rel_tol=1e-2, opt_abs_tol=1e6;
    unsigned max_iter=1e3;
    SolverType ub_solver_type = ipopt, lb_solver_type = gurobi;
    
    double max_time = 300;
    DebugOn("Instance "<<fname<<endl);
    
//    m->reset();
    bool upper_bound_heur=false;

    if(upper_bound_heur){
        vector<double> sol(m->_nb_vars);
        auto m2=m->copy();
        auto m1=m->copy();
        var<double> Xa=m2->get_var<double>("X");
        var<double> Xija=m2->get_var<double>("Xij");
        // m2->min(sum(Xija)-sum(Xa));
        auto o=25*(*m2->_obj)+sum(Xija)-sum(Xa);
        m2->min(o);
        Constraint<> diag_cut("diag_cut");
        diag_cut=sum(Xija)-sum(Xa);
        m2->add(diag_cut<=0);
        //m2->print();
        solver<> s(m2,ipopt);
        s.run();
        m2->print_solution();
        auto y2d=m2->get_var<double>("y");
        auto y1=m1->get_var<int>("y");
        
        for(auto k:*y2d.get_keys()){
            if(y2d.eval(k)<0.5){
                y1.set_lb(k,0);
                y1.set_ub(k,0);
                
            }
            else{
                y1.set_lb(k,1);
                y1.set_ub(k,1);
            }
        }
        solver<> s1(m1,gurobi);
        s1.run();
        //check_PSD(m1);
        if(m1->_status==0){
            m1->get_solution(sol);
        }
    }
    solver<> sc(m,gurobi);
    bool relax = false;
    auto ts=get_wall_time();
    if(root_refine){
        //auto m2=m->copy();
        sc.run(relax = true);
       // m->print_solution();
//      upper_bound_heur=true;
//        if(upper_bound_heur){
//            auto y2d=m->get_var<int>("y");
//            auto y2=m2->get_var<int>("y");
//
//            for(auto k:*y2d.get_keys()){
//                if(y2d.eval(k)<0.5){
//                    y2.set_lb(k,0);
//                    y2.set_ub(k,0);
//
//                }
//                else{
//                    y2.set_lb(k,1);
//                    y2.set_ub(k,1);
//                }
//            }
//            solver<> s2(m2,gurobi);
//            s2.run();
//            //check_PSD(m1);
//            if(m2->_status==0){
//                vector<double> sol(m->_nb_vars);
//                m2->get_solution(sol);
//                m->set_solution(sol);
//            }
//
//        }
    }
    sc.run();
    auto tf=get_wall_time();
    m->print_solution();
    m->print_constraints_stats(1e-9);
    
    auto eig_value=m->check_PSD();
    
    string out_file_name=fname;
    auto pos=out_file_name.find_last_of("/");
    out_file_name=out_file_name.substr(pos+1);
    pos=out_file_name.find_first_of(".");
    auto out_file_name1=out_file_name.substr(0,pos);
    out_file_name=string(prj_dir)+"/results_misdp/"+out_file_name1+".txt";
    ofstream fout(out_file_name.c_str());
    fout<<out_file_name1<<"&"<<m->get_obj_val()<<"&"<<tf-ts<<"&"<<eig_value<<"&"<<m->_rel_obj_val<<"&"<<m->num_cuts[0]<<"&"<<m->num_cuts[1]<<"&"<<m->num_cuts[2]<<"&"<<m->num_cuts[3]<<"&"<<m->num_cuts[4]<<"&"<<m->num_cuts[5]<<"\n";
    fout.close();
    cout<<out_file_name1<<" obj "<<m->get_obj_val()<<" time "<<tf-ts<<" smallest eig "<<eig_value<<" lower bound  "<<m->_rel_obj_val<<endl;
    cout<<"soc incumbent "<<m->num_cuts[0]<<" eig bags incumbent "<<m->num_cuts[1]<<" eig full incumbent "<<m->num_cuts[2]<<" soc mipnode "<<m->num_cuts[3]<<" eig bags mipnode "<<m->num_cuts[4]<<" eig full mipnode "<<m->num_cuts[5]<<"\n";
}
