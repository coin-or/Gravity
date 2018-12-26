////
////  Node.cpp
////  Cycle_Basis_PF
////  adapt from power tools 
////  Created by Hassan on 17/06/2014.
//
////
//
//#include <gravity/Node.h>
//#include <gravity/Arc.h>
//#include <iostream>
//#include <limits.h>
//
//using namespace std;
//
//Node::Node(){};
//
//Node::Node(string name, int id):_name(name),_id(id),fill_in(0){};
//
//Node::~Node(){};
//
//Node* Node::clone(){
//    Node* copy = new Node();
//    copy->_id = _id;
//    copy->_name = _name;
//    copy->fill_in = 0; // fill_in is associated with the topology. so it should be 0.
//    return copy;
//};
//
///*
// @brief Adds a to the list of incident arcs
// */
//
//void Node::addArc(Arc* a){
//    branches.push_back(a);
//}
//
//
///*
// @brief Find and remove incident arc from list of branches
// @return 0 if a was found and removed, -1 oterwise
// */
//int Node::removeArc(Arc* a){
//    vector<Arc*>::iterator it = branches.begin();
//    while (it != branches.end()) {
//        if((*it) == a){            
//            it = branches.erase(it);
//            return 0;
//        }
//        it++;
//    }
//    return -1;
//}
//
//bool Node::is_connected(Node* n){
//    for (auto a:branches) {
//        if (n->_id==a->neighbour(this)->_id) {
//            return true;
//        }
//    }
//    for (auto a:n->branches) {
//        if (_id==a->neighbour(n)->_id) {
//            return true;
//        }
//    }
//    return false;
//}
//
//void Node::update_fill_in(Node* n){
//    Node * nn = nullptr;
//    
//    for(auto a:branches){
//        nn = a->neighbour(this); //this node.
//        // if nn is null
//        if (nn->_id==n->_id) {
//            continue; //self connect
//        }
//        if (!n->is_connected(nn)) {
//            fill_in++; // if this node is connected to node.
//        }
//        
//    }
//}
//
//
//std::vector<Arc*> Node::get_out(){
//    vector<Arc*> res;
//    for (auto a:branches) {
//        if(a->_src->_id==_id){
//            res.push_back(a);
//        }
//    }
//    return res;
//}
//
//std::vector<Arc*> Node::get_in(){
//    vector<Arc*> res;
//    for (auto a:branches) {
//        if(a->_dest->_id==_id){
//            res.push_back(a);
//        }
//    }
//    return res;
//}
//
//std::set<Node*> Node::get_neighbours(){
//    set<Node*> res;
//    for (auto a:branches) {
//        //if(a->_dest->_id=_id && std::find(res.begin(),res.end(), a->_src)== res.end()){
//       if(a->_dest->_id== _id){
//            res.insert(a->_src);
//        }
//        
//        //if(a->_src->_id==_id && std::find(res.begin(),res.end(), a->_dest)== res.end() ){
//        //if(a->_src->_id==_id && std::find(res.begin(),res.end(), a->_dest)== res.end() ){
//            if(a->_src->_id==_id ){
//            res.insert(a->_dest);
//        }
//    }
//    // uniqueness.
//    
//    return res;
//}
//
