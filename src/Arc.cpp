//
//  Arc.cpp
//  Cycle_Basis_PF
//
//

#include <gravity/Arc.h>
#include <iostream>

using namespace std;

Arc::Arc() {
    _name=  "";
    _src = nullptr;
    _dest = nullptr;
}

Arc::~Arc() {}

Arc::Arc(string name):_src(NULL), _dest(NULL) {
    _name = name;
}


Arc::Arc(Node* s, Node* d) {
    _src = s;
    _dest = d;
    _weight = 1;
    _name = s->_name + "," + d->_name;
//   in_cycle = false;
    //  parallel = false;
    //  connect();
}

Arc::Arc(Node*s, Node* d, double w) {
    _src = s;
    _dest = d;
    _weight = w;
}

Arc* Arc::clone() {
    Arc* copy = new Arc(_name);
    copy->_src = _src;
    copy->_dest = _dest;
    copy->_weight = _weight;
    copy->_id = _id;
    copy->_active = _active;
    return copy;
}

/* @brief Returns the neighbour of n if n is a node of the arc, null otherwise */
Node* Arc::neighbour(Node* n) {
    Node* neigh = NULL;
    if (_src == n)
        neigh = _dest;
    if (_dest == n)
        neigh = _src;
    return neigh;
}


/* Connects the current arc to its source and _destination, adding itself to the list of branches in these nodes */
void Arc::connect() {
    if (_src == _dest){
        //throw invalid_argument ("It is now allowed to make a node self connected in gravity. \n");
        std::cout << "It is now allowed to make a node self connected in gravity" << endl;
    }
    _src->update_fill_in(_dest);// update the fill-ins
    _dest->update_fill_in(_src);
    Node* common = nullptr;
    // just for source.
    for (auto a:_src->branches) {
        common = a->neighbour(_src);
        if (common->is_connected(_dest)) {
            common->fill_in--;
//            assert(common->fill_in >=0);
        }
    }
    _src->addArc(this);
    _dest->addArc(this);
}


void Arc::print() {
    std::cout << "(" << _src->_id << ", " << _dest->_id << ")" <<std::endl;
}
