//
//  Node.cpp
//  Cycle_Basis_PF
//  adapt from power tools 
//  Created by Sumiran on 17/06/2014.
//  Copyright (c) 2014 NICTA. All rights reserved.
//

#include <Gravity/Node.h>
#include <Gravity/Arc.h>
#include <iostream>
#include <limits.h>

using namespace std;

Node::Node(){};

Node::Node(int id):_name(""),ID(id),fill_in(0){};
Node::Node(string name, int id):_name(name),ID(id),fill_in(0){};
Node::~Node(){};

Node* Node::clone(){
    Node* copy = new Node();
    copy->ID = ID;
    return copy;
};

/*
 @brief Adds a to the list of incident arcs
 */
void Node::addArc(Arc* a){
    branches.push_back(a);
}


/*
 @brief Find and remove incident arc from list of branches
 @return 0 if a was found and removed, -1 oterwise
 */
int Node::removeArc(Arc* a){
    vector<Arc*>::iterator it = branches.begin();
    while (it != branches.end()) {
        if((*it) == a){
            it = branches.erase(it);
            return 0;
        }
        it++;
    }
    return -1;
}



bool Node::is_connected(Node* n){
    for (auto a:branches) {
        if (n->ID==a->neighbour(this)->ID) {
            return true;
        }
    }
    for (auto a:n->branches) {
        if (ID==a->neighbour(n)->ID) {
            return true;
        }
    }
    return false;
}

void Node::update_fill_in(Node* n){
    Node * nn = nullptr;
    for(auto a:branches){
        nn = a->neighbour(this);
        if (nn->ID==n->ID) {
            continue;
        }
        if (!n->is_connected(nn)) {
            fill_in++;
        }
        
    }
}




