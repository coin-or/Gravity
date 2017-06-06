//
//  Arc.h
//  Cycle_Basis_PF
//
//  Created by Sumiran on 18/06/2014.
//  Copyright (c) 2014 NICTA. All rights reserved.
//

#ifndef Cycle_Basis_PF_Arc_h
#define Cycle_Basis_PF_Arc_h
#include <Gravity/Node.h>
#include "assert.h"
#include "string"

class Arc{
public:
    int id;
    std::string _name;
    Node* src;
    Node* dest;
  //  double weight;
    
    /* @brief Returns the neighbour of n if n is a node of the arc, null otherwise */
    Node* neighbour(Node* n);
    
//  bool in_cycle;
//  Path* horton_path;
    
    
 /* @brief Returns the neighbour of n if n is a node of the arc, null otherwise */
 //   Node* neighbour(Node* n);
    
    Arc(std::string name);
    ~Arc();
    Arc(Node* s, Node* d);
    Arc* clone();
    
    /* Connects the current arc to its source and destination, adding itself to the list of branches in these nodes */
    void connect();
    
    void print();
};

#endif
