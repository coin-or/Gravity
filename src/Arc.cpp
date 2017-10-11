//
//  Arc.cpp
//  Cycle_Basis_PF
//
//  Created by Sumiran on 18/06/2014.
//  Copyright (c) 2014 NICTA. All rights reserved.
//

#include <gravity/Arc.h>
#include <iostream>

using namespace std;

Arc::Arc(){ _name=  ""; _src = nullptr; _dest = nullptr;}

Arc::~Arc(){}

Arc::Arc(string name):_src(NULL), _dest(NULL){
    _name = name;
}


Arc::Arc(Node* s, Node* d){
    _src = s;
    _dest = d;
    _weight = 1;
 //   in_cycle = false;
  //  parallel = false;
  //  connect();
}

Arc::Arc(Node*s, Node* d, double w){
    _src = s;
    _dest = d;
    _weight = w;
}

Arc* Arc::clone(){
    Arc* copy = new Arc(_name);
    copy->_src = _src;
    copy->_dest = _dest;
    copy->_weight = _weight;
    copy->_id = _id;
    return copy;
}

/* @brief Returns the neighbour of n if n is a node of the arc, null otherwise */
Node* Arc::neighbour(Node* n){
    Node* neigh = NULL;
    if (_src == n)
        neigh = _dest;
    if (_dest == n)
        neigh = _src;
    return neigh;
}


/* Connects the current arc to its source and _destination, adding itself to the list of branches in these nodes */
void Arc::connect(){
    _src->update_fill_in(_dest);// update the fill-ins
    _dest->update_fill_in(_src);
    Node* common = nullptr;
    // just for source. 
    for (auto a:_src->branches) {
        common = a->neighbour(_src);
        if (common->is_connected(_dest)) {
            common->fill_in--;
            assert(common->fill_in >=0);
        }
    }
    _src->addArc(this);
    _dest->addArc(this);
}

void Arc::print(){
    std::cout << "(" << _src->_id << ", " << _dest->_id << ")" <<std::endl;

}
