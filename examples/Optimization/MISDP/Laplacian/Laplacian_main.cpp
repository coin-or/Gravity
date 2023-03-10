//
//  Laplacian_main.cpp
//
//  Created by Hassan Hijazi on March 9 2023.
//

#include <stdio.h>
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
        if(v.size()==2 || v.size()==3)
            res.push_back(v);
    }
    return res;
}

using namespace gravity;
using namespace std;

int main(int argc, char * argv[]){
    string fname=string(prj_dir)+"/data_sets/MISDP/8_1.json";
    
    if(argc>=2){
        fname=argv[1];
    }

    std::ifstream i(fname);
    json e;
    i >> e;
    vector<pair<pair<int,int>,double>> edges_to_augment;
    int num_nodes = e["num_nodes"];
    DebugOn("Number of nodes = " << num_nodes << endl);
    edges_to_augment.insert(end(edges_to_augment), begin(e["edges_to_augment"]), end(e["edges_to_augment"]));
    DebugOn("Number of edges = " << edges_to_augment.size() << endl);

    /* Graph */
    Net graph;
    for (int i = 0; i < num_nodes; i++) {
        auto n = new Node(to_string(i));
        graph.add_node(n);
    }
    int index = 0;
    for (const auto &e:edges_to_augment) {
        string src = to_string(e.first.first-1);
        string dest = to_string(e.first.second-1);
        double weight = e.second;
        string name = src+","+dest;
        auto arc = new Arc(name);
        arc->_id = index++;
        arc->_src = graph.get_node(src);
        arc->_dest = graph.get_node(dest);
        arc->_weight = weight;
        arc->connect();
        graph.add_arc(arc);
    }

    /* Indices */
    indices E = indices(graph.arcs);
    indices N = indices(graph.nodes);
    indices S("S"), E_S("E_S");
    S = range(1,std::pow(2,num_nodes) - (num_nodes*(num_nodes-1))/2 - num_nodes - 1);/* excluding empty S, |S| = 1, and |S| = 2*/
    param<> rhs("rhs");
    rhs.in(S);
    E_S = E;
    int row_id = 0, count = 0;
    vector<pair<int,std::vector<string>>> _bag_names;
    auto subsets = allPossibleSubset(num_nodes);
    for (auto const &subset: subsets) {
        pair<int,vector<string>> bn;
        bn.first=count++;
        for(auto n:subset){
            bn.second.push_back(to_string(n));
            DebugOn(n <<"\t");
        }
        DebugOn(endl);
        _bag_names.push_back(bn);
        E_S.add_empty_row();
        rhs.set_val(row_id, subset.size()-1);
        for (int i = 0; i<subset.size()-1; i++) {
            for (int j = i+1; j<subset.size(); j++) {
                E_S.add_in_row(row_id, to_string(subset[i])+","+to_string(subset[j]));
            }
        }
        row_id++;
    }
    indices delta_i("delta_i");
    delta_i = E;
    for (int i = 0; i<num_nodes; i++) {
        auto node = graph.nodes.at(i);
        delta_i.add_empty_row();
        for (const Arc* a: node->get_out()) {
            delta_i.add_in_row(i, a->_name);
        }
        for (const Arc* a: node->get_in()) {
            delta_i.add_in_row(i, a->_name);
        }
    }

    indices delta_i_excl0("delta_i_excl0"), delta_i_excln("delta_i_excln"), delta_0("delta_0");
    delta_i_excl0 = E;
    delta_i_excln = E;
    delta_0 = E;
    auto node_0 = graph.nodes.at(0);
    for (const Arc* a: node_0->get_out()) {
        delta_0.add_in_row(0, a->_name);
    }
    for (const Arc* a: node_0->get_in()) {
        delta_0.add_in_row(0, a->_name);
    }
    for (int i = 0; i<num_nodes-1; i++) {
        auto node_0 = graph.nodes.at(i);
        auto node_1 = graph.nodes.at(i+1);
        delta_i_excl0.add_empty_row();
        delta_i_excln.add_empty_row();
        for (const Arc* a: node_0->get_out()) {
            delta_i_excln.add_in_row(i, a->_name);
        }
        for (const Arc* a: node_0->get_in()) {
            delta_i_excln.add_in_row(i, a->_name);
        }
        for (const Arc* a: node_1->get_out()) {
            delta_i_excl0.add_in_row(i, a->_name);
        }
        for (const Arc* a: node_1->get_in()) {
            delta_i_excl0.add_in_row(i, a->_name);
        }
    }

    /* Parameters*/
    param<> w("w");
    w.in(E);
    index = 0;
    for(const Arc* a: graph.arcs){
        w.set_val(index++, a->_weight);
    }

    bool project = false;
    /* Model */
    Model<> Laplacian("Laplacian");
//    Laplacian._bag_names = _bag_names;
    DebugOn("Adding " << Laplacian._bag_names.size() << " bags\n");
    /* Variables */    var<> gamma("𝛾");
    if(!project)
        Laplacian.add(gamma.in(R(1)));
    var<int> x("x", 0, 1);
    Laplacian.add(x.in(E));
    var<> Wii("Wii", pos_);
    Laplacian.add(Wii.in(N));
    var<> Wij("Wij");
    Laplacian.add(Wij.in(E));


    /* Constraints */
    if(project){
        indices E_excl0("E_excl0"), E_excln("E_excln");
        for(const Arc* a: graph.arcs){
            if(a!=graph.arcs.at(0))
                E_excl0.add(a->_name);
            if(a!=graph.arcs.back())
                E_excln.add(a->_name);
        }
        Constraint<> Wij_def("Wij_def");
        Wij_def = Wij.in(E_excln) + w.in(E_excln)*x.in(E_excln) - (Wij.in(E_excl0) + w.in(E_excl0)*x.in(E_excl0));
        Laplacian.add(Wij_def.in(E_excl0) == 0);

        indices N_excl0("N_excl0"), N_excln("N_excln");
        for(const Node* n: graph.nodes){
            if(n!=graph.nodes.at(0))
                N_excl0.add(n->_name);
            if(n!=graph.nodes.back())
                N_excln.add(n->_name);
        }

        Constraint<> Wii_def("Wii_def");
        Wii_def = w.in(delta_i_excln)*x.in(delta_i_excln) - Wii.in(N_excln) - (w.in(delta_i_excl0)*x.in(delta_i_excl0) - Wii.in(N_excl0));
        Laplacian.add(Wii_def.in(N_excl0) == 0);

        Constraint<> Wii_W_ij_def("Wii_W_ij_def");
        Wii_W_ij_def = (num_nodes/(num_nodes-1.))*(w.in(delta_0)*x.in(delta_0) - Wii(N._keys->front())) - num_nodes*(Wij(E._keys->front()) + w(E._keys->front())*x(E._keys->front()));
        Laplacian.add(Wii_W_ij_def == 0);
    }
    else{
        Constraint<> Wij_def("Wij_def");
        Wij_def = Wij + w*x - (1./num_nodes)*gamma;
        Laplacian.add(Wij_def.in(E) == 0);

        Constraint<> Wii_def("Wii_def");
        Wii_def = Wii - (w.in(delta_i)*x.in(delta_i) - (num_nodes-1.)/num_nodes*gamma);
        Laplacian.add(Wii_def.in(N) == 0);
    }

    Constraint<> NodeCut("NodeCut");
    NodeCut = x.in(delta_i);
    Laplacian.add(NodeCut.in(N) >= 1);

    Constraint<> Spanning_tree("Spanning_tree");
    Spanning_tree = sum(x);
    Laplacian.add(Spanning_tree == num_nodes - 1);

//    Constraint<> Subtour("Subtour");
//    Subtour = x.in(E_S) - rhs;
//    Laplacian.add(Subtour.in(S) <= 0);

//    Constraint<> Minor2("Minor2");
//    Minor2 = Wij*Wij - Wii.from(E)*Wii.to(E);
//    Laplacian.add(Minor2.in(E) <= 0);

    Laplacian.make_PSD(Wii,Wij);
//    Laplacian.max(gamma);

//    double coef_ij = num_nodes/(2.*E.size());
//    Laplacian.max(coef*(sum(Wij) + product(w,x)));
    if(project){
        double coef_i = 1./((num_nodes - 1));
        Laplacian.max(coef_i*(2.*product(w,x) - sum(Wii)));// + coef_ij*(sum(Wij) + product(w,x)));
        //    Laplacian.max(coef_i*(w.in(delta_0)*x.in(delta_0) - Wii(N._keys->front())));// + coef_ij*(sum(Wij) + product(w,x)));
    }
    else{
        Laplacian.max(gamma);
        /* Map linking matrix entries to decision vars */
        vector<vector<pair<pair<size_t,size_t>, double>>> Wij_gamma_map(E.size());
        vector<vector<pair<pair<size_t,size_t>, double>>> Wii_gamma_map(N.size());
        vector<vector<pair<pair<size_t,size_t>, double>>> Wij_x_map(E.size());
        vector<vector<pair<pair<size_t,size_t>, double>>> Wii_x_map(N.size());



        for (int i = 0; i<num_nodes; i++) {
            auto node = graph.nodes.at(i);
            Wii_gamma_map[i].push_back({{gamma.get_vec_id(),0},-(num_nodes-1.)/num_nodes});
            for (const Arc* a: node->get_out()) {
                Wii_x_map[i].push_back({{x.get_vec_id(),a->_id},a->_weight});
            }
            for (const Arc* a: node->get_in()) {
                Wii_x_map[i].push_back({{x.get_vec_id(),a->_id},a->_weight});
            }
        }

        index = 0;
        for(const Arc* a: graph.arcs){
            Wij_gamma_map[index].push_back({{gamma.get_vec_id(),0},1./num_nodes});
            Wij_x_map[index].push_back({{x.get_vec_id(),a->_id},-a->_weight});
            index++;
        }
        bool orig_cuts = true;/* Generate the SDP cuts using original variables (get rid of Wii and Wij) */
        if(!project && orig_cuts){
            Laplacian.Xii_x_map = Wii_gamma_map;
            Laplacian.Xij_x_map = Wij_gamma_map;
            Laplacian.Xii_y_map = Wii_x_map;
            Laplacian.Xij_y_map = Wij_x_map;
        }
    }

    Laplacian.print();
//    return 0;
    solver<> sv(Laplacian,cplex);
    sv.run();
    Laplacian.print_solution();
    return 0;
}
