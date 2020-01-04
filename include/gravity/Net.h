//  Net.h
//  Network
//  Guanglei WANG
//  Copyright (c) All rights reserved
//

#ifndef Cycle_Basis_PF_Net_h
#define Cycle_Basis_PF_Net_h
#include <map>
#include <set>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <memory>
#include <assert.h>
//#include <gravity/types.h>
#include <gravity/Node.h>
#include <gravity/Arc.h>
#include <gravity/Path.h>

class Net{

public:
    std::string _name;

    /** Set of nodes */
    std::vector<Node*> nodes;

    /** Set of existing + potential arcs */
    std::vector<Arc*> arcs;

    /** Set of contingency arcs */
    std::vector<Arc*> conting_arcs;

    /** Set of existing arcs */
    std::vector<Arc*> _exist_arcs;

    /** Set of bus pairs */
    gravity::node_pairs _bus_pairs;

    /** Set of bus pairs in the chordal completion */
    gravity::node_pairs _bus_pairs_chord;

    /** Mapping the directed arcs to their source-_destination by their names, i.e, (name_src, name_dest)*/
    std::map<std::string, std::map<string,Arc*>*> arcID;

    /** Mapping the line name to the line pointer */
    std::map<std::string, Arc*> arcMap;


    /** Mapping the node name to its position in the vector, key = node name */
    std::map<std::string, Node*> nodeID;

    /** Vector of cycles forming a cycle basis */
    std::vector<Path*> cycle_basis;
    
    /** Horton subnetwork */
    Net* horton_net = nullptr;
    
    /** Is a tree */
    bool _tree = false;
    /** Indices */
    gravity::indices bus_pairs = gravity::indices("bus_pairs"), bus_pairs_chord = gravity::indices("bus_pairs_chordal");

    bool duplicate(std::string name1, std::string name2, int id1);

    // bags are sorted in an ascending order of ids.
    std::vector<pair<string,std::vector<Node*>>> _bags; // each node is from this nodes.
    std::vector<pair<string,std::vector<Node*>>> _bags_copy; // each node is a copy of the original node (not by reference).

    /** Clone the graph exactly **/
    Net* clone();

    /** Clone to get a copied (Undirected) graph: no parallel lines, at most one
     * arc between two nodes*/
    // not that clone_undirected is often used for implementing graph algorithms.
    Net* clone_undirected();

    Net();
    ~Net();

    /** Modifiers */
    void add_node(Node* n);
    /* Add node at specified position in nodes */
    void add_node(Node* n, size_t pos);
    bool add_arc(Arc* a);
    void add_undirected_arc(Arc* a);

    /* Accessors */
    Node* get_node(std::string name) const;

    /** returns the arc formed by node n1 and n2 */
    Arc* get_arc(Node* n1, Node* n2);

    /** returns the arc formed by node names n1 and n2 */
    Arc* get_arc(std::string n1, std::string n2);

    bool has_arc(std::string n1, std::string n2) {
        return get_arc(n1,n2)!=nullptr;
    }

    bool has_directed_arc(Node* n1, Node* n2){
        Arc* a = get_arc(n1, n2);
            if (n1->_id == a->_src->_id && n2->_id==a->_dest->_id) {
                return true;
            }
        return false;
    }



    char* readline(FILE *input);
    void exit_input_error(int line_num);
    void read_adjacency_matrix(const string& fname);
    void readrudy(const char* fname);
    void get_complement(const char* fname);//  (i, j): i<j \E

    /**  @brief Remove node and all incident arcs from the network
     @note Does not remove the incident arcs from the list of arcs in the network!
     @return the id of the node removed
     */
    std::string remove_end_node();

    /** Erase Horton network */
    void clear_horton_net();

    /**  @brief Remove node and all incident arcs from the network
     @note Does not remove the incident arcs from the list of arcs in the network!
     @return the id of the node removed
     */
    void remove_arc(Arc* a);

    /** Compute the tree decomposition bags **/
    void get_tree_decomp_bags();
    
    std::vector<std::vector<Node*>> decompose_bags_4d(bool print_bags=false);
    std::vector<pair<string,vector<Node*>>> decompose_bags_3d(bool print_bags=false);

    void print() const;
    
    /** Sort nodes in decreasing degree */
    void orderNodes(Net* net);
    
    /** Sort arcs in decreasing weight */
    void orderArcs(Net* net);
    
    void add_horton_nodes(Net* net);
    
    void add_horton_branches(Net* net);
    
    /* Algorithms */
    
    /** Reset all distances to max value */
    void resetDistance();
    
    /** sets the cycle member fo all Nodes to be false */
    void reset_nodeCycle();
    
    /** sets all nodes to unexplored */
    void reset_nodeExplored();
    
    /** Computes and returns the shortest path between src and dest in net */
    Path* Dijkstra(Node* src, Node* dest, Net* net);
    
    /** Computes a spanning tree of minimal weight in horton's network, then adds the corresponding cycles to the original network by connecting to src
     @note Implements Kruskal’s greedy algorithm
     */
    void minimal_spanning_tree(Node* src, Net* net);
    
    /** Combines the src node with the path p to form a cycle */
    void combine(Node* src, Path* p);
    
    vector<Path*> get_cycle_basis();
    
    void Fast_Horton(Net *net);

    gravity::indices get_pairs_chord(const vector<pair<string,vector<Node*>>>& bags);
    
    
    /** get algorithmic graph */
    void get_algorithmic_graph(); // a cloned graph without in-active, parallel lines.

    /** Return a chordal extension graph with tree decomposition **/
    Net* get_chordal_extension();

    /** Compute the vector of bus pairs, ignoring parallel lines **/
    gravity::indices get_bus_pairs();
    
    
     /** Compute the vector of reference bus pairs, ignoring parallel lines **/
     gravity::indices get_ref_bus_pairs();
    

    /** Compute the tree decomposition bags **/
    void  get_cliquebags(bool print=false); // remove bags that are not maximal cliques.
    Net* get_clique_tree();

    /** Linear algebra based methods based on Armadillo*/
    void chol_decompose(bool print=false);

    Arc *get_directed_arc(std::string src, std::string dest);


    vector<gravity::index_pair *> get_bus_pairs_all();
};
#endif
