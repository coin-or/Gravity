//
//  Node.h
//  Cycle_Basis_PF
//
//  Created by Sumiran on 17/06/2014.
//  Copyright (c) 2014 NICTA. All rights reserved.
//

#ifndef Cycle_Basis_PF_Node_h
#define Cycle_Basis_PF_Node_h
#include <vector>
#include <string>

//class Path;
class Arc;
//class Net;

class Node{
    
public:
    
    std::string _name;
    int ID;
    std::vector<Arc*> branches;
    
    /* the number of edges needed to make the subgraph formed by adjacent nodes a clique */
    int fill_in;

    
    // constructions
    Node();
    Node(int id);
    Node(std::string name, int id);
    ~Node();
    Node* clone();
    
    void addArc(Arc* a);
    int removeArc(Arc* a);
    void update_fill_in(Node* n);
    
    /*
     @brief Returns true if n is an adjacent node.
     */
    bool is_connected(Node* n);
    
};

#endif
