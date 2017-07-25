//
//  Path.cpp
//  Cycle_Basis_PF
//
//  Created by Sumiran on 17/06/2014.
//  Copyright (c) 2014 NICTA. All rights reserved.
//

#include <gravity/Path.h>
#include <fstream>


/* @brief Returns true if the pair (n1,n2) is a source-destination pair for this path */
bool Path::source_dest(Node* n1, Node* n2){
    if(nodes.front()->ID==n1->ID && nodes.back()->ID==n2->ID)
        return true;
    
    if(nodes.front()->ID==n2->ID && nodes.back()->ID==n1->ID)
        return true;
    return false;
}

/* Returns the length of the path */
int Path::length(){
    return (int)nodes.size();
}

/* Returns true if the path is a cycle */
bool Path::cycle(){
    return (nodes.front()==nodes.back());
}

Path* Path::clone(){
    Path* newp = new Path();
    for (std::list<Node*>::iterator i = nodes.begin(); i != nodes.end(); i++) {
        newp->nodes.push_back((*i));
    }
    return newp;
}
