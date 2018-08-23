//
//  Net.cpp
//
//
//  Created by Guanglei Wang on 03/06/2017.
//

#include <gravity/Net.h>
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
//#include <armadillo>
#ifdef USE_BOOST
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <deque>
#include <iterator>
#endif

using namespace std;
using namespace gravity;

static int max_line_len;
static char* line = nullptr;

Net::Net() {}

/* returns true if an arc is already present between the given nodes */
bool Net::duplicate(std::string n1, std::string n2, int id1) {
    int id2 = get_arc(n1, n2)->_id;
    
    if (id2 < id1)
        return true;
    else
        return false;
}


Net* Net::clone() {
    Net* copy_net = new Net();
    Node* node = NULL;
    
    for (int i=0; i<nodes.size(); i++) {
        node = this->nodes[i];
        copy_net->add_node(node->clone());
    }
    
    Arc* arc = NULL;
    for (int i=0; i < arcs.size(); i++) {
        
        /* ignores if the arc is a paralel line to an already existing arc */
        //if (duplicate(arcs[i]->_src->_name, arcs[i]->_dest->_name, arcs[i]->_id)) {
//            continue;
//        }
        arc = arcs[i]->clone();
        
        /* Update the source and destination to the new nodes in copy_net */
        arc->_src = copy_net->get_node(arc->_src->_name);
        arc->_dest = copy_net->get_node(arc->_dest->_name);
        
        /* Add the new arc to the list of arcs */
        copy_net->add_arc(arc);
        
        /* Connects it to its source and destination */
        arc->connect();
    }
    return copy_net;
}

Net* Net::clone_undirected() {
    Net* copy_net = new Net();
    Node* node = NULL;
    
    for (int i=0; i<nodes.size(); i++) {
        node = this->nodes[i];
        copy_net->add_node(node->clone());
    }
    
    Arc* arc = NULL;
    for (int i=0; i < arcs.size(); i++) {
        if (copy_net->get_arc(arcs[i]->_src->_name, arcs[i]->_dest->_name)!=nullptr) {
            continue;
        }
        arc = arcs[i]->clone();
        
        /* Update the source and destination to the new nodes in copy_net */
        arc->_src = copy_net->get_node(arc->_src->_name);
        arc->_dest = copy_net->get_node(arc->_dest->_name);
        
        /* Add the undirected arc to the list of arcs */
        copy_net->add_undirected_arc(arc);
        
        /* Connects it to its source and destination */
        arc->connect();
    }
    return copy_net;
}

const bool bag_compare(const vector<Node*> & a,const vector<Node*>& b) {
    return a.size() > b.size();
}


const bool node_compare(const Node* n1, const Node* n2) {
    return n1->fill_in > n2->fill_in;
}

void Net::add_node(Node* node) {
    node->_id = (int) nodes.size();
    
    if (!nodeID.insert(pair<string,Node*>(node->_name, node)).second) {
        cerr << "ERROR: adding the same node twice!";
    }
    
    nodes.push_back(node);
}

Node* Net::get_node(string name){
    return nodeID.find(name)->second;
}

/* returns the undirected arc formed by node n1 and n2 */
Arc* Net::get_arc(Node* n1, Node* n2) {
    string src, dest, key, inv_key;
    src = n1->_name;
    dest = n2->_name;
    key.clear();
    inv_key.clear();
    key.append(src);
    inv_key.append(dest);
    key.append(",");
    inv_key.append(",");
    key.append(dest);
    inv_key.append(src);
    map<string, set<Arc*>*>::iterator it= arcID.find(key);
    if (it != arcID.end()) {
        for (auto a: *it->second) {
            //   if (!a->parallel) {
            return a;
            // }
        }
    }
    it = arcID.find(inv_key);
    if (it != arcID.end()) {
        for (auto a: *it->second) {
            //   if (!a->parallel) {
            return a;
            // }
        }
    }
    return nullptr;
}

/* returns the Id of the arc formed by nodes names n1 and n2 */
// this assumes that the graph is undirected.
Arc* Net::get_arc(std::string src, std::string dest) {
    std::string key, inv_key;
    key.clear();
    inv_key.clear();
    key.append(src);
    inv_key.append(dest);
    key.append(",");
    inv_key.append(",");
    key.append(dest);
    inv_key.append(src);
    map<string, set<Arc*>*>::iterator it= arcID.find(key);
    if (it != arcID.end()) {
        for (auto a: *it->second) {
            return a;
        }
    }
    
    it = arcID.find(inv_key);
    if (it != arcID.end()) {
        for (auto a: *it->second) {
            return a;
        }
    }
    
    return nullptr;
}

Arc* Net::get_directed_arc(std::string src, std::string dest) {
    std::string key;
    key.clear();
    key.append(src);
    key.append(",");
    key.append(dest);
    map<string, set<Arc*>*>::iterator it= arcID.find(key);
    if (it != arcID.end()) {
        for (auto a: *it->second) {
            return a;
        }
    }

    return nullptr;
}

bool Net::add_arc(Arc* a) {
    bool parallel = false;
    set<Arc*>* s = NULL;
    string src, dest, key;
    src = a->_src->_name;
    dest = a->_dest->_name;
    if (src == dest){
        throw invalid_argument ("It is now allowed to make a node self connected in gravity. \n");
        
    }
    
    
    key.clear();
    key.append(src);
    key.append(",");
    key.append(dest);
    
    if(arcID.find(key)==arcID.end()) {
        s = new set<Arc*>;
        s->insert(a);
        arcID.insert(pair<string, set<Arc*>*>(key,s));
    }
    else {
        if(arcID.find(key)!=arcID.end())
            s = arcID[key];
        s->insert(a);
        DebugOff("\nWARNING: adding another Directed line between same nodes! \n Node ID: " << src << " and Node ID: " << dest << endl);
        a->_parallel = true;
        parallel = true;
    }
    arcs.push_back(a);
    return parallel;
}
// undirected
void Net::add_undirected_arc(Arc* a) {
    bool parallel = false;
    set<Arc*>* s = NULL;
    string src, dest, key, key_inv;
    src = a->_src->_name;
    dest = a->_dest->_name;

    if (src == dest){
        throw invalid_argument ("It is now allowed to make a node self connected in gravity. \n");
    
    }
    
    key.clear();
    key.append(src);
    key.append(",");
    key.append(dest);
    
    key_inv.clear();
    key_inv.append(dest);
    key_inv.append(",");
    key_inv.append(src);
    
    if(arcID.find(key)==arcID.end()&& arcID.find(key_inv)==arcID.end()) {
        s = new set<Arc*>;
        s->insert(a);
        arcID.insert(pair<string, set<Arc*>*>(key,s));
        arcs.push_back(a);
    }
}


/** remove the arc by
1. removing it from the arcs list, but without changing the arc id.
2. removing it from the map container.
*/

void Net::remove_arc(Arc* a) {
    arcs.erase(arcs.begin()+(a->_id));
    //arcs[a->_id] = nullptr;
    arcID.erase(a->_src->_name+","+a->_dest->_name);
}

// Reading files
char* Net::readline(FILE *input)
{
    size_t len;
    // line, max_line_len have been declared
    if(std::fgets(line,max_line_len,input)==NULL) return NULL;

    while(strrchr(line,'\n') == NULL)
    {
        max_line_len *= 2;
        line = (char *)realloc(line,max_line_len);
        len = strlen(line);
        if(fgets(line+len,max_line_len-len,input) == NULL)
            break;
    }
    return line;
}

void Net::exit_input_error(int line_num) {
    fprintf(stderr,"Wrong input format at line %d\n", line_num);
    exit(1);
}

// Reading graphs with rudy format
void Net::readrudy(const char* fname) {
    int Num_nodes=0;
    int Num_edges=0;

    ifstream infile(fname);

    string sLine;

    if (infile.good())
    {
        getline(infile, sLine);
        istringstream iss(sLine);
        iss >> Num_nodes;
        iss >> Num_edges;
    }
    else{
        fprintf(stderr,"can’t open input file %s\n",fname);
        exit(1);
    }
    // get nodes

    string name;
    Node* node = nullptr;

    for (int i= 1; i< Num_nodes + 1; i++) {
        name = to_string(i);
        node = new Node(name,i-1);
        add_node(node);
    }


    // get arcs
    Arc* arc = NULL;

    // note that src, dest are names of nodes.
    string src, dest;
    double weight;
    while(getline(infile,sLine,'\n'))
    {
        istringstream iss(sLine);
        iss >> src >> dest >> weight;

        name = (int)arcs.size()+1;
        arc = new Arc(name);

        arc->_id = (int)arcs.size();

        arc->_src = get_node(src);
        arc->_dest= get_node(dest);
        arc->_weight=weight;
        add_arc(arc);
        arc->connect();
    }
    infile.close();
}



/** construct a graph by reading an adjacency matrix */
void Net::read_adjacency_matrix(const char* fname) {
    FILE *fp = fopen(fname,"r");
    if(fp == NULL)
    {
        fprintf(stderr,"can’t open input file %s\n",fname);
        exit(1);
    }

    max_line_len = 1024;
    line = new char[max_line_len];

    vector<vector<int>> matrix;
    int temp;
    while(readline(fp)!=NULL)
    {
        vector<int> row;
        stringstream linestream(line);
        while (linestream>>temp)
            row.push_back(temp);
        matrix.push_back(row);
    }
    
    int n=0;
    n =matrix.size();

    string name;
    int id = 0;

    Node* node = NULL;
    for (int i= 0; i<n; i++) {
        name = to_string(i);
        node = new Node(name,i);
        add_node(node);
    }

    Arc* arc = NULL;
    string src, dest;
    unsigned index = 0;
    for (int i = 0; i <(n); i++)
        for (int j=i+1; j<n; j++) {
            if (matrix[i][j] > 0)
            {
                src = to_string(i);
                dest = to_string(j);
                id = index;
                arc = new Arc(src + "," + dest);
                arc->_id = id;
                arc->_src = get_node(src);
                arc->_dest= get_node(dest);
                add_arc(arc);
                arc->connect();
            }
            index++;
        }
    delete[] line;
    fclose(fp);
}

void Net::get_complement(const char* fname) {
    FILE *fp = fopen(fname,"r");
    if(fp == NULL)
    {
        fprintf(stderr,"can’t open input file %s\n",fname);
        exit(1);
    }

    max_line_len = 1024;
    line = new char[max_line_len];

    vector<vector<int>> matrix;
    int temp;
    while(readline(fp)!=NULL)
    {
        vector<int> row;
        stringstream linestream(line);
        while (linestream>>temp)
            row.push_back(temp);
        matrix.push_back(row);
    }
    int n=0;
    n =matrix.size();

    string name;
    int id = 0;

    Node* node = NULL;
    for (int i= 0; i<n; i++) {
        name = to_string(i);
        node = new Node(name,i);
        add_node(node);
    }

    Arc* arc = NULL;
    string src, dest;

    for (int i = 0; i <(n-1); i++)
        for (int j=i+1; j<n; j++) {
            if (matrix[i][j] == 0)
            {
                src = to_string(i);
                dest = to_string(j);
                id = (int)arcs.size();
                arc = new Arc(to_string(id));
                arc->_id = id;
                arc->_src = get_node(src);
                arc->_dest= get_node(dest);
                add_arc(arc);
                arc->connect();
            }
        }
    delete[] line;
    fclose(fp);
}



/*  @brief Remove node and all incident arcs from the network
 @note Does not remove the incident arcs from the list of arcs in the network!
 @return the id of the node removed
 */
string Net::remove_end_node() {
    Node* n = nodes.back();
    Node * nn = nullptr;
    string n_id = n->_name;
    for (auto a: n->branches) {
        nn = a->neighbour(n);
        nn->removeArc(a);
        for (auto aa: nn->branches) {
            if (!aa->neighbour(nn)->is_connected(n)) {
                nn->fill_in--;
                assert(nn->fill_in >=0);
            }
        }
    }
    nodes.pop_back();
    return n_id;
}

// use greedy fill-in algorithm.
void Net::get_tree_decomp_bags(bool print_bags) {
    Node* n = nullptr;
    Node* u = nullptr;
    Node* nn = nullptr;
    Arc* arc = nullptr;
    
    string name="";
    Net* graph_clone = clone_undirected(); //
    int nb = 0;
    
    /** cliques with less than 1 nodes are useless for us.*/
    while (graph_clone->nodes.size()> 2) {
        sort(graph_clone->nodes.begin(), graph_clone->nodes.end(),node_compare);
        
        // last element has the minimum fill-in.
        n = graph_clone->nodes.back();
        if(!n->_active) {
            graph_clone->remove_end_node();
            continue;
        }
        Debug(n->_name << endl);
        Debug(graph_clone->nodes.size() << endl);
        vector<Node*> bag_copy;
        vector<Node*> bag;
        DebugOff("new bag = { ");
        for (auto nn: n->get_neighbours()) {
            if(!nn->_active) continue;
            bag_copy.push_back(nn);
            bag.push_back(get_node(nn->_name)); // Note it takes original node.
            DebugOff(nn->_name << ", ");
        }
        DebugOff(n->_name << "}\n");
        graph_clone->remove_end_node();
        bag_copy.push_back(n);
        bag.push_back(get_node(n->_name)); // node in this graph
        sort(bag_copy.begin(), bag_copy.end(), [](const Node* a, const Node* b) -> bool{return a->_id < b->_id;});
        sort(bag.begin(), bag.end(), [](const Node* a, const Node* b) -> bool{return a->_id < b->_id;});
        
        // update clone_graph and construct chordal extension.
        for (int i = 0; i < bag_copy.size(); i++) {
            u = bag_copy.at(i);
            for (int j = i+1; j<bag_copy.size(); j++) {
                nn = bag_copy.at(j);
                if (u->is_connected(nn)) {
                    if(get_arc(u,nn) && !get_arc(u,nn)->_active) {
                        Arc* off_arc = get_arc(u,nn);
                        off_arc->_imaginary = true;
                        off_arc->_free = true;
                    }
                    continue;
                }
                name = to_string((int) graph_clone->arcs.size()+1);
                arc = new Arc(name);
                
                arc->_id = arcs.size();
                arc->_src = u;
                arc->_dest = nn;
                arc->_imaginary = true;
                arc->_free = true;
                arc->connect();
                graph_clone->add_undirected_arc(arc);
            }
        }
        
        if (true) {
            //            DebugOn("bag_copy = {");
            //            for (int i=0; i<bag_copy.size();     i++) {
            //                cout << bag_copy.at(i)->_name << " ";
            //            }
            //            DebugOn("}" << endl);
            DebugOff("bag = {");
//            for (int i=0; i<bag.size();     i++) {
//                cout << bag.at(i)->_name << " ";
//            }
            DebugOff("}" << endl);
        }
        //        _bags_copy.push_back(bag_copy);
        _bags.push_back(bag); // bag original
        
        if (bag_copy.size()==3) {
            nb++;
        }
//        else if(bag_copy.size()>3){
//            DebugOn("Decomposing bigger bag into 3d bags\n");
//            
//            for (auto i = 0; i<bag_copy.size()-2; i++) {
//                for (auto j = i+1; j<bag_copy.size()-1; j++) {
//                    for (auto k = j+1; k<bag_copy.size(); k++) {
//                        vector<Node*> new_bag;
//                        new_bag.push_back(bag[i]);
//                        new_bag.push_back(bag[j]);
//                        new_bag.push_back(bag[k]);
//                        DebugOn("new bag = {");
//                        for (int i=0; i<new_bag.size();     i++) {
//                            cout << new_bag.at(i)->_name << " ";
//                        }
//                        DebugOn("}" << endl);
//                        _bags.push_back(new_bag);
//                    }
//                }
//            }
//        }
        delete n;
    }
    //    sort(_bags.begin(), _bags.end(), bag_compare);
    
    
    Debug("\n Number of 3D bags = " << nb << endl);
//    DebugOn("\n Total number of bags = " << _bags.size() << endl);
    
    delete graph_clone;
    
}

/** Return the vector of arcs ignoring parallel lines **/
std::vector<gravity::index_pair*> Net::get_bus_pairs(){
    return _bus_pairs._keys;
}

/** Return the vector of arcs of the chordal completion ignoring parallel lines **/
std::vector<gravity::index_pair*> Net::get_bus_pairs_chord(){
    return _bus_pairs_chord._keys;
}

Net* Net::get_chordal_extension() {
    Node* n = nullptr;
    Node* u = nullptr;
    Node* nn = nullptr;
    Arc* arc = nullptr;
    Arc* arc_chordal = nullptr;

    Node* u_chordal = nullptr;
    Node* nn_chordal = nullptr;

    string name="";
    string name_chordal="";
    Net* chordal_extension = clone();
    Net* graph_clone = clone_undirected();
    int nb = 0;

    /** cliques with less than 1 nodes are useless for us.*/
    while (graph_clone->nodes.size() > 1) {
        sort(graph_clone->nodes.begin(), graph_clone->nodes.end(),node_compare);
        // last element has the minimum fill-in.
        n = graph_clone->nodes.back();         
        Debug(n->_name << endl);
        Debug(_clone->nodes.size() << endl);
        vector<Node*> bag_copy;
        vector<Node*> bag;
        Debug("new bag_copy = { ");

        for (auto nn: n->get_neighbours()) {
            bag_copy.push_back(nn);
            bag.push_back(get_node(nn->_name));
            Debug(nn->_name << ", ");
        }

        graph_clone->remove_end_node();
        bag_copy.push_back(n);
        bag.push_back(get_node(n->_name)); // node in this graph
        sort(bag_copy.begin(), bag_copy.end(),[](const Node* a, const Node* b) -> bool{return a->_id < b->_id;});
        sort(bag.begin(), bag.end(),[](const Node* a, const Node* b) -> bool{return a->_id < b->_id;});

        // update graph_graph and construct chordal extension.
        for (int i = 0; i < bag_copy.size() - 1; i++) {
            u = bag_copy.at(i);
            u_chordal = chordal_extension->get_node(u->_name);
            for (int j = i+1; j<bag_copy.size(); j++) {
                nn = bag_copy.at(j);
                nn_chordal=chordal_extension->get_node(nn->_name);
                if (u->is_connected(nn)) {
                    continue;
                }
                name = to_string((int) graph_clone->arcs.size()+1);
                name_chordal = to_string((int)chordal_extension->arcs.size()+1);

                arc = new Arc(name);
                arc_chordal = new Arc(name_chordal);

                arc->_id = arcs.size();
                arc->_src = u;
                arc->_dest = nn;
                arc->connect();
                graph_clone->add_undirected_arc(arc);

                arc_chordal->_id = chordal_extension->arcs.size();
                arc_chordal->_src = u_chordal;
                arc_chordal->_dest = nn_chordal;
                arc_chordal->connect();
                chordal_extension->add_undirected_arc(arc_chordal);
            }
        }
        _bags_copy.push_back(bag_copy);
        _bags.push_back(bag);
        if (bag_copy.size()==3) {
            nb++;
        }
        delete n;
    }
    // sort the bags by its size (descending order)
    sort(_bags.begin(), _bags.end(), bag_compare);
    printf("With greedy fill-in algirithm, the chordal graph added  %lu edges \n", (chordal_extension->arcs.size() - arcs.size()));

    delete graph_clone;
    return chordal_extension;
}

// get cliques from the tree decomposition
// Two methods
// first one: check the inclusion relationship
// second one: use the RIP property of the tree decomposition, thus just need to check every leaf..
// One need to execute either get_tree_decomposition or get_chordal_extension first, then run get_clique_tree.

// use _bags instead of bag_copy
void Net::get_cliquebags (bool print) {
    for (unsigned i = 0; i < _bags.size(); i++) {
        for (unsigned j = i+1; j < _bags.size();) {
            if (std::includes(_bags[i].begin(),_bags[i].end(),
                              _bags[j].begin(), _bags[j].end()))
            {
                _bags.erase(_bags.begin()+j);
            }
            else
                j++;
        }
    }
    cout << "Number of maximal cliques of the chordal extension = " << _bags.size() << endl <<endl;
}

/* Destructors */
Net::~Net() {
    if (!nodes.empty()) {
        for (vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); it++) {
            delete (*it);
        }
        nodes.clear();
    }
    if(!arcs.empty()) {
        for (vector<Arc*>::iterator it = arcs.begin(); it != arcs.end(); it++) {
            if(*it)
                delete (*it);
        }
        arcs.clear();
    }

    if(!cycle_basis.empty()) {
        for (vector<Path*>::iterator it = cycle_basis.begin(); it != cycle_basis.end(); it++) {
            delete (*it);
        }
        arcs.clear();
    }
    for (pair<string,set<Arc*>*> it:arcID) {
        delete it.second;
    }
}

Net* Net::get_clique_tree(){
    Net* cliquetree = new Net();
    Node* node = nullptr;
    Arc*  a = nullptr;
    string name;
    get_cliquebags(true);
#ifdef USE_BOOST
    /** Note that we also need the edge information of the clique tree **/
    /** boost graph library or implement the expanded version of MCS algorithm by Blair and Peyton */
    typedef boost::adjacency_list <boost::vecS,
    boost::vecS,
    boost::undirectedS,
    boost::no_property,
    boost::property < boost::edge_weight_t, int >> Graph;
    typedef boost::graph_traits <Graph>::edge_descriptor Edge;
    //typedef boost::graph_traits <Graph>::vertex_descriptor Vertex;
    
    // BUILD THE INTERSECTION GRAPH OF THE CLIQUES
    typedef std::pair<int, int> E;
    std::vector<E> edges;
    std::vector<int> weights;
    int nb_cliques = this->_bags.size();
    for (int i = 0; i < nb_cliques; i++) {
        DebugOn("bag " << i << " has " << this->_bags[i].size() << " nodes." <<endl);
        sort(this->_bags[i].begin(), this->_bags[i].end());
        for (int j = i +1; j < nb_cliques; j++) {
            vector<Node*> v3;
            sort(this->_bags[j].begin(), this->_bags[j].end());
            set_intersection(this->_bags[i].begin(), this->_bags[i].end(), this->_bags[j].begin(), this->_bags[j].end(), back_inserter(v3));
            if (v3.size() > 0) {
                edges.push_back(E(i, j));
                weights.push_back(-v3.size());
            }
        }
    }
    //size_t num_edges = edges.size();
    
#if defined(BOOST_MSVC) && BOOST_MSVC <= 1300
    Graph g(num_nodes);
    boost::property_map<Graph, edge_weight_t>::type weightmap = get(edge_weight, g);
    for (std::size_t j = 0; j < num_edges; ++j) {
        Edge e;
        bool inserted;
        boost::tie(e, inserted) = boost::add_edge(edges[j].first, edges[j].second, g);
        boost::weightmap[e] = weights[j];
    }
#else
    Graph g(edges.begin(), edges.end(), weights.begin(), nb_cliques);
#endif
    boost::property_map < Graph, boost::edge_weight_t >::type weight = get(boost::edge_weight, g);
    std::vector < Edge > spanning_tree;
    boost::kruskal_minimum_spanning_tree(g, std::back_inserter(spanning_tree));
    
    DebugOn("Print the total " << spanning_tree.size() << " edges in the clique tree:" << endl);
    
    //////////CLIQUE TREE /////////////////////////////
    for (int i = 0; i < nb_cliques; i++) {
        node= new Node(to_string(i), i);
        cliquetree->add_node(node);
    }
    
    for (std::vector < Edge >::iterator ei = spanning_tree.begin();
         ei != spanning_tree.end(); ++ei) {
        int u = source(*ei, g);
        int v = target(*ei, g);
        DebugOn(u << " <--> " << v
                << " with weight of " << -weight[*ei]
                << endl);
        name = (int) cliquetree->arcs.size();
        a = new Arc(name);
        a->_id = cliquetree->arcs.size();
        
        // intersection
        vector<Node*> v3;
        sort(this->_bags[u].begin(), this->_bags[u].end());
        sort(this->_bags[v].begin(), this->_bags[v].end());
        set_intersection(this->_bags[u].begin(), this->_bags[u].end(),
                         this->_bags[v].begin(), this->_bags[v].end(),
                         back_inserter(v3));
        a->_src = cliquetree->get_node(to_string(u));
        a->_dest = cliquetree->get_node(to_string(v));
        a->_weight = -weight[*ei];
        a->_intersection = v3;
        cliquetree->add_arc(a);
        a->connect();
        
        for (int i = 0; i < v3.size(); i++){
                auto  node = v3.at(i);
            for (int j = i+1; j < v3.size(); j++){
                auto arc = get_arc(node, v3.at(j)); 
                if (arc != nullptr){
                    a->_intersection_clique.push_back(new index_pair(index_(arc->_src->_name), index_(arc->_dest->_name), arc->_active)); 
                }
             //   else
               //     a->_intersection_clique.push_back(new index_pair(index_(node->_name), index_(v3.at(j)->_name), true));
            }
        }
    }
#endif
    return cliquetree;
}


//void Net::chol_decompose(bool print){
//    arma::mat adjacency_matrix = arma::zeros(nodes.size(), nodes.size());
//    for (auto &arc: arcs){
//        adjacency_matrix(arc->_src->_id, arc->_dest->_id) = 1;
//        adjacency_matrix(arc->_dest->_id, arc->_src->_id) = 1;
//
//    }
//    arma::mat pertubation = arma::zeros(nodes.size(),nodes.size());
//    // if fails,  make epsilon larger.
//    unsigned epsilon = 2;
//    arma::mat adjacency_matrix_psd = adjacency_matrix + epsilon*pertubation.eye();
//    arma::mat chordal_sparsity = arma::zeros(nodes.size(),nodes.size());
//    chordal_sparsity = arma::chol(adjacency_matrix_psd);
//    if (print){
//        cout << "chordal sparsity matrix: \n " << chordal_sparsity << endl;
//    }
//    unsigned num_nonzeros = 0;
//    for (int i = 0; i < nodes.size(); i++)
//        for (int j = i+1; j< nodes.size(); j++){
//            if (chordal_sparsity(i, j) != 0){
//                num_nonzeros +=1;
//            }
//        }
//    printf("With cholesky decomposition, the chordal graph added  %lu edges \n", (num_nonzeros - arcs.size()));
//}


std::vector<gravity::index_pair*> Net::get_bus_pairs_all(){
    vector<gravity::index_pair*> res;
    string ni, nj;
    for(int i = 0; i < nodes.size()-1; i++) {
        for(int j = i+1; j < nodes.size(); j++) {
            ni = nodes[i]->_name;
            nj = nodes[j]->_name;
            res.push_back(new index_pair(index_(ni), index_(nj), 1));
        }
    }
    return res;
}
