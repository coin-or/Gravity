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
#include <assert.h>
#include <gravity/types.h>
#include <gravity/Node.h>
#include <gravity/Arc.h>
#include <gravity/Path.h>

class Net{
    
public:
    std::string _name;
    
    /** Set of nodes */
    std::vector<Node*> nodes;
    
    /** Set of arcs */
    std::vector<Arc*> arcs;
    
    /** Set of bus pairs */
    gravity::node_pairs _bus_pairs;
    
    /** Mapping the arcs to their source-_destination by their names, i.e, (name_src, name_dest)*/
    
    std::map<std::string, std::set<Arc*>*> arcID;
    
    /** Mapping the node name to its position in the vector, key = node name */
    std::map<std::string, Node*> nodeID;
    
    /** Vector of cycles forming a cycle basis */
    std::vector<Path*> cycle_basis;
    
    bool duplicate(std::string name1, std::string name2, int id1);
    
    /** Tree decomposition bags */
    std::vector<std::vector<Node*>> _bags;
    
    /** Cloning */
    Net* clone();
    
    Net();
    ~Net();
   
    /** Modifiers */
    void add_node(Node* n);
    bool add_arc(Arc* a);
    
    /* Accessors */
    Node* get_node(std::string name);

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
    void read_adjacency_matrix(const char* fname);
    void readrudy(const char* fname);
    void get_complement(const char* fname);//  (i, j): i<j \E
    
    /**  @brief Remove node and all incident arcs from the network
     @note Does not remove the incident arcs from the list of arcs in the network!
     @return the id of the node removed
     */
    std::string remove_end_node();
    
    
    /**  @brief Remove node and all incident arcs from the network
     @note Does not remove the incident arcs from the list of arcs in the network!
     @return the id of the node removed
     */
    void remove_arc(Arc* a);
    
    /** Compute the tree decomposition bags **/
    void get_tree_decomp_bags(bool print_bags = false);
    
    
    
    /** Return a chordal extension graph with tree decomposition **/
    Net* get_chordal_extension();
    
    /** Compute the vector of bus pairs, ignoring parallel lines **/
    std::vector<gravity::index_pair*> get_bus_pairs();

    
    
    /** Compute the tree decomposition bags **/
    void get_clique_tree (bool print_clique = false);
};
#endif
