//
//  PowerNet.cpp
//
//
//  Created by Guanglei Wang on 03/06/2017.
//

#include "PoolNet.h"
#include <algorithm>
#include <map>
#define _USE_MATH_DEFINES
#include <cmath>
#include <list>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <queue>
#include <time.h>
//#include <xlnt/xlnt.hpp>
#include <gravity/solver.h>

using namespace std;

PoolNet::PoolNet() {
    //in inputs_pools
    x_min.set_name("x_min");
    x_max.set_name("x_max");
    //in pools_outputs
    y_min.set_name("y_min");
    y_max.set_name("y_max");
    //in inputs_outputs
    z_min.set_name("z_min");
    z_max.set_name("z_max");
    //in inputs
    cost_ip.set_name("cost_ip");
    cost_io.set_name("cost_io");
    cost_po.set_name("cost_po");
    avail_min.set_name("avail_min");
    avail_max.set_name("avail_max");
    inqual.set_name("inqual");
    //in outputs
    rev.set_name("rev");
    dem_min.set_name("dem_min");
    dem_max.set_name("dem_max");
    outqual_min.set_name("outqual_min");
    outqual_max.set_name("outqual_max");
    //in pools
    pool_cap.set_name("pool_cap");
 
}

PoolNet::~PoolNet() {
    
    for (Node* n:nodes) {
        delete n;
    }
    nodes.clear();
    for (Arc* a:arcs) {
        delete a;
    }
    arcs.clear();
}

indices PoolNet::out_arcs_per_node() const{
    auto ids = indices(arcs);
    ids._type = matrix_;
    ids.set_name("out_arcs_per_node");
    ids._ids = make_shared<vector<vector<size_t>>>();
    ids._ids->resize(nodes.size());
    string key;
    size_t inst = 0;
    for (auto n: nodes) {
        if (n->_active) {
            for (auto a:n->get_out()) {
                if (!a->_active) {
                    continue;
                }
                key = a->_name;
                ids._ids->at(inst).push_back(a->_id);
            }
            
            
            inst++;
        }
    }
    return ids;
}

indices PoolNet::out_arcs_to_pool_per_input() const{
    indices ids = indices("out_arcs_to_pool_per_input");
    int row_id = 0;
    for(const string& in_pool_id: *this->Inputs._keys){
        auto in_pool = nodeID.at(in_pool_id);
        ids.add_empty_row();
        for (const Arc* out: in_pool->get_out()) {
            if(!out->_free)
                ids.add_in_row(row_id, out->_name);
        }
        row_id++;
    }
    return ids;
}

indices PoolNet::out_arcs_to_output_per_input() const{
    indices ids = indices("out_arcs_to_output_per_input");
    int row_id = 0;
    for(const string& in_pool_id: *this->Inputs._keys){
        auto in_pool = nodeID.at(in_pool_id);
        ids.add_empty_row();
        for (const Arc* out: in_pool->get_out()) {
            if(out->_free)
                ids.add_in_row(row_id, out->_name);
        }
        row_id++;
    }
    return ids;
}


indices PoolNet::in_arcs_per_pool_attr() const{
    indices ids = indices("in_arcs_per_pool_attrs");
    int row_id = 0;
    for(const string& key_id: *this->pools_attr._keys){
        ids.add_empty_row();
        auto pos = nthOccurrence(key_id, ",", 1);
        auto pool_id = key_id.substr(0,pos);
        auto pool = nodeID.at(pool_id);
        for (const Arc* out: pool->get_in()) {
            ids.add_in_row(row_id, out->_name);
        }
        row_id++;
    }
    return ids;
}


indices PoolNet::in_arcs_per_pool() const{
    indices ids = indices("in_arcs_per_pool");
    int row_id = 0;
    for(const string& key_id: *this->Pools._keys){
        auto pool = nodeID.at(key_id);
        ids.add_empty_row();
        for (const Arc* out: pool->get_in()) {
            ids.add_in_row(row_id, out->_name);
        }
        row_id++;
    }
    return ids;
}

indices PoolNet::in_arcs_attr_per_pool() const{
    indices ids = indices("in_arcs_attr_per_pool");
    int row_id = 0;
    for(const string& key_id: *this->pools_attr._keys){
        auto pos = nthOccurrence(key_id, ",", 1);
        auto attr = key_id.substr(pos+1);
        auto pool_id = key_id.substr(0,pos);
        auto pool = nodeID.at(pool_id);
        ids.add_empty_row();
        for (const Arc* out: pool->get_in()) {
            auto input_id = out->_src->_name;
            ids.add_in_row(row_id, input_id+","+attr);
        }
        row_id++;
    }
    return ids;
}


indices PoolNet::out_arcs_per_pool_attr() const{
    indices ids = indices("out_arcs_per_pool_attr");
    int row_id = 0;
    for(const string& key_id: *this->pools_attr._keys){
        auto pos = nthOccurrence(key_id, ",", 1);
        auto pool_id = key_id.substr(0,pos);
        auto pool = nodeID.at(pool_id);
        ids.add_empty_row();
        for (const Arc* out: pool->get_out()) {
            ids.add_in_row(row_id, out->_name);
        }
        row_id++;
    }
    return ids;
}


indices PoolNet::out_arcs_per_pool() const{
    indices ids = indices("out_arcs_per_pool");
    int row_id = 0;
    for(const string& key_id: *this->Pools._keys){
        auto pool = nodeID.at(key_id);
        ids.add_empty_row();
        for (const Arc* out: pool->get_out()) {
            ids.add_in_row(row_id, out->_name);
        }
        row_id++;
    }
    return ids;
}


indices PoolNet::in_arcs_from_pool_per_output_attr() const{
    indices ids = indices("in_arcs_from_pool_per_output_attr");
    int row_id = 0;
    for(const string& key_id: *this->outputs_attr._keys){
        auto pos = nthOccurrence(key_id, ",", 1);
        auto output_id = key_id.substr(0,pos);
        auto output = nodeID.at(output_id);
        ids.add_empty_row();
        for (const Arc* out: output->get_in()) {
            if(!out->_free)
                ids.add_in_row(row_id, out->_name);
        }
        row_id++;
    }
    return ids;
}

indices PoolNet::in_arcs_from_pool_per_output() const{
    indices ids = indices("in_arcs_from_pool_per_output");
    int row_id = 0;
    for(const string& key_id: *this->Outputs._keys){
        auto output = nodeID.at(key_id);
        ids.add_empty_row();
        for (const Arc* out: output->get_in()) {
            if(!out->_free)
                ids.add_in_row(row_id, out->_name);
        }
        row_id++;
    }
    return ids;
}

indices PoolNet::in_arcs_from_input_per_output() const{
    indices ids = indices("in_arcs_from_input_per_output");
    int row_id = 0;
    for(const string& out_id: *this->Outputs._keys){
        auto out_node = nodeID.at(out_id);
        ids.add_empty_row();
        for (const Arc* out: out_node->get_in()) {
            if(out->_free)
                ids.add_in_row(row_id, out->_name);
        }
        row_id++;
    }
    return ids;
}

indices PoolNet::in_arcs_from_input_per_output_attr() const{
    indices ids = indices("in_arcs_from_input_per_output_attr");
    int row_id = 0;
    for(const string& outat_id: *this->outputs_attr._keys){
        auto pos = nthOccurrence(outat_id, ",", 1);
        auto output_id = outat_id.substr(0,pos);
        auto out_node = nodeID.at(output_id);
        ids.add_empty_row();
        for (const Arc* out: out_node->get_in()) {
            if(out->_free)
                ids.add_in_row(row_id, out->_name);
        }
        row_id++;
    }
    return ids;
}

indices PoolNet::in_arcs_per_node(indices arcset) const{
    auto ids = indices(arcs);
    ids._type = matrix_;
    ids.set_name("in_arcs_per_node_in_"+arcset.get_name());
    ids._ids = make_shared<vector<vector<size_t>>>();
    ids._ids->resize(nodes.size());
    string key;
    size_t inst = 0;
    for (auto n: nodes) {
        for (auto a:*(arcset._keys)) {
             if(arcMap.find(a)!=arcMap.end()){
                auto arca= arcMap.find(a)->second;
                if(arca->_dest->_name==n->_name){
                    ids._ids->at(inst).push_back(arca->_id);
                }
            }
        }
        inst++;
    }
    return ids;
}



indices PoolNet::in_arcs_per_node() const{
    auto ids = indices(arcs);
    ids._type = matrix_;
    ids.set_name("in_arcs_per_node");
    ids._ids = make_shared<vector<vector<size_t>>>();
    ids._ids->resize(nodes.size());
    string key;
    size_t inst = 0;
    for (auto n: nodes) {
        if (n->_active) {
            for (auto a:n->get_in()) {
                if (!a->_active) {
                    continue;
                }
                key = a->_name;
                auto it1 = ids._keys_map->find(key);
                if (it1 == ids._keys_map->end()){
                    throw invalid_argument("In function in_arcs_per_node(), unknown key.");
                }
                ids._ids->at(inst).push_back(it1->second);
            }
            inst++;
        }
    }
    return ids;
}


void PoolNet::readgrid() {
    
    string fname=string(prj_dir)+"/data_sets/Pooling/Adhya1_gms.txt";
    string word="", tempwrd, numwrd;
    
    int flag;
    double val;
    
    int N_node,N_input,N_output,N_pool, N_attr;
    
    ifstream file(fname.c_str(), std::ifstream::in);
    
    if(!file.is_open()) {
        throw invalid_argument("Could not open file " + fname);
    }
    
    
    while(word.find("Declare")==string::npos){
    getline(file, word);
    }

    getline(file, word);
    
    tempwrd=word.substr(word.find_first_of("*")+1);
    numwrd=tempwrd.substr(0, tempwrd.find_first_of(" "));

   
    N_node = atoi(numwrd.c_str());
    
    getline(file, word);
    
    tempwrd=word.substr(word.find_first_of("*")+1);
    numwrd=tempwrd.substr(0, tempwrd.find_first_of(" "));
    
    
    N_input = atoi(numwrd.c_str());

    getline(file, word);
    
    tempwrd=word.substr(word.find_first_of("/")+2);
    numwrd=tempwrd.substr(0, tempwrd.find_first_of("*"));
    
    
    N_output = N_node+1-atoi(numwrd.c_str());
    
    
    getline(file, word);
    
    tempwrd=word.substr(word.find_first_of("*")+1);
    numwrd=tempwrd.substr(0, tempwrd.find_first_of(" "));
    
    
    N_attr = atoi(numwrd.c_str());
    Attr=indices("Attr");
    Attr = range(1,N_attr);
    while(word.find("table c(i,j)")==string::npos){
        getline(file, word);
    }


    Inputs=indices("Inputs");
    Outputs=indices("Outputs");
    Pools=indices("Pools");
    
    
    N_pool=N_node-N_input-N_output;
    for (auto i=1;i<=N_input;i++)
    {
        Node* node = NULL;
        node = new Node(to_string(i));
        add_node(node);
        Inputs.add(to_string(i));
    }
    
    
    
    
    for (auto i=N_input+1;i<=N_input+N_pool;i++)
    {
        Node* node = NULL;
        node = new Node(to_string(i));
        add_node(node);
        Pools.add(to_string(i));
    }
    
    pools_attr = indices(Pools,Attr);
    
    for (auto i=N_input+N_pool+1;i<=N_input+N_pool+N_output;i++)
    {
        Node* node = NULL;
        node = new Node(to_string(i));
        add_node(node);
        Outputs.add(to_string(i));
    }
    
    unsigned index = 0;
    string src, dest;
    
    while(word.find("table a(i,j)")==string::npos){
        getline(file, word);
    }
     getline(file, word);
    
    inputs_pools=indices("inputs_pools");
    inputs_outputs=indices("inputs_outputs");
    pools_outputs=indices("pools_outputs");
    
    for(auto i=1;i<=N_input;i++){
        file>>flag;
        for(auto j=N_input+1;j<=N_input+N_pool;j++){
            file>>flag;
            if(flag==1){
                Arc* arc = NULL;
                
                auto src = get_node(to_string(i));
                auto dest= get_node(to_string(j));
                
                arc = new Arc(to_string(i) + "," + to_string(j));
                arc->_id = index++;
                arc->_src = src;
                arc->_dest= dest;
                this->add_arc(arc);
                arc->connect();
                inputs_pools.add(to_string(i) + "," + to_string(j));
            }
        }
        
        for(auto j=N_input+N_pool+1;j<=N_input+N_pool+N_output;j++){
            file>>flag;
            if(flag==1){
                Arc* arc = NULL;
                
                auto src = get_node(to_string(i));
                auto dest= get_node(to_string(j));
                
                arc = new Arc(to_string(i) + "," + to_string(j));
                arc->_id = index++;
                arc->_src = src;
                arc->_dest= dest;
                arc->_free=true;
                this->add_arc(arc);
                arc->connect();
                inputs_outputs.add(to_string(i) + "," + to_string(j));
                
            }
        }
    }
    
    for(auto i=N_input+1;i<=N_input+N_pool;i++){
        file>>flag;
        for(auto j=1;j<=N_pool;j++){
            file>>flag;
                    }
        
        for(auto j=N_input+N_pool+1;j<=N_input+N_pool+N_output;j++){
            file>>flag;
            if(flag==1){
                Arc* arc = NULL;
                
                auto src = get_node(to_string(i));
                auto dest= get_node(to_string(j));
                
                arc = new Arc(to_string(i) + "," + to_string(j));
                arc->_id = index++;
                arc->_src = src;
                arc->_dest= dest;
                this->add_arc(arc);
                arc->connect();
                pools_outputs.add(to_string(i) + "," + to_string(j));
                
            }
        }
    }
    
    while(word.find("table q(i,k)")==string::npos){
        getline(file, word);
    }
    getline(file, word);
 
    inputs_attr=indices("inputs_attr");
    outputs_attr=indices("outputs_attr");
    for(auto i=1;i<=N_input;i++){
        file>>flag;
        for(auto j=1;j<=N_attr;j++){
            file>>val;
            inqual.add_val(to_string(i)+","+to_string(j), val);
            inputs_attr.add(to_string(i) + "," + to_string(j));
        }
    }
    for(auto i=N_input+N_pool+1;i<=N_input+N_pool+N_output;i++){
        file>>flag;
        for(auto j=1;j<=N_attr;j++){
            file>>val;
            outqual_max.add_val(to_string(i)+","+to_string(j), val);
            outqual_min.add_val(to_string(i)+","+to_string(j), 0);
            outputs_attr.add(to_string(i) + "," + to_string(j));
        }
    }
    
 
    while(word.find("capacity lower bound")==string::npos){
        getline(file, word);
    }
    file>>word;
      file>>word;
      file>>word;
    
    for(auto i=1;i<=N_input;i++){
        file>>flag;
       
            file>>val;
            avail_min.add_val(to_string(i), val);
       
    }
    
    for(auto i=N_input+N_pool+1;i<=N_input+N_pool+N_output;i++){
        file>>flag;
        file>>val;
        dem_min.add_val(to_string(i), val);
    }
    while(word.find("capacity upper bound")==string::npos){
        getline(file, word);
    }
    file>>word;
    file>>word;
    file>>word;
    
    for(auto i=1;i<=N_input;i++){
        file>>flag;
        
        file>>val;
        avail_max.add_val(to_string(i), val);
        
    }
    for(auto i=N_input+1;i<=N_input+N_pool;i++){
        file>>flag;
        file>>val;
        pool_cap.add_val(to_string(i), val);
    }
    for(auto i=N_input+N_pool+1;i<=N_input+N_pool+N_output;i++){
        file>>flag;
        file>>val;
        dem_max.add_val(to_string(i), val);
    }
    
    for(auto key: *(inputs_pools._keys)){
        x_min.add_val(key, 0);
        x_max.add_val(key, 100);
     
        
    }
    for(auto key: *(pools_outputs._keys)){
        y_min.add_val(key, 0);
        y_max.add_val(key, 100);
        
    }
    for(auto key: *(inputs_outputs._keys)){
        z_min.add_val(key, 0);
        z_max.add_val(key, 100);
        
    }
    
    
    
    for(auto key: *(Outputs._keys)){
        rev.add_val(key, 0);
        
    }
    file.clear();
    file.seekg(0, ios::beg);
    while(word.find("table c(i,j)")==string::npos){
        getline(file, word);
    }
    getline(file, word);
    for(auto i=1;i<=N_input;i++){
        file>>flag;
        for(auto j=N_input+1;j<=N_input+N_pool;j++){
            file>>val;
    auto key=to_string(i)+","+to_string(j);
            cost_ip.add_val(key, val);
      }

        for(auto j=N_input+N_pool+1;j<=N_input+N_pool+N_output;j++){
           file>>val;
            auto key=to_string(i)+","+to_string(j);
            cost_io.add_val(key, val);
        }
    }
    for(auto i=N_input+1;i<=N_input+N_pool;i++){
        file>>flag;
        for(auto j=1;j<=N_pool;j++){
            file>>val;
        }
                for(auto j=N_input+N_pool+1;j<=N_input+N_pool+N_output;j++){
               file>>val;
        auto key=to_string(i)+","+to_string(j);
        cost_po.add_val(key, val);
    }
    }
    file.close();
    
}



