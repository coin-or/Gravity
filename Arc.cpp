//
//  Arc.cpp
//  Cycle_Basis_PF
//
//  Created by Sumiran on 18/06/2014.
//  Copyright (c) 2014 NICTA. All rights reserved.
//

#include "Arc.h"
#include <iostream>

using namespace std;

Arc::~Arc(){}

Arc::Arc(string name):src(NULL), dest(NULL){}


Arc::Arc(Node* s, Node* d){
    src = s;
    dest = d;
 //   in_cycle = false;
  //  parallel = false;
  //  connect();
}

Arc* Arc::clone(){
    Arc* copy = new Arc(_name);
    copy->src = src;
    copy->dest = dest;
 //   copy->in_cycle = in_cycle;
    return copy;
}

/* @brief Returns the neighbour of n if n is a node of the arc, null otherwise */
Node* Arc::neighbour(Node* n){
    Node* neigh = NULL;
    if (src == n)
        neigh = dest;
    if (dest == n)
        neigh = src;
    return neigh;
}


/* Connects the current arc to its source and destination, adding itself to the list of branches in these nodes */
void Arc::connect(){
    src->update_fill_in(dest);
    dest->update_fill_in(src);
    Node* common = nullptr;
    for (auto a:src->branches) {
        common = a->neighbour(src);
        if (common->is_connected(dest)) {
            common->fill_in--;
            assert(common->fill_in >=0);
        }
    }
    src->addArc(this);
    dest->addArc(this);
}

void Arc::print(){
    std::cout << "(" << src->ID << ", " << dest->ID << ")" <<std::endl;

}
