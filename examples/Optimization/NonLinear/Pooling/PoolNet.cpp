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
    inqual.set_name("inqual");
    //in pools_outputs
    y_min.set_name("y_min");
    y_max.set_name("y_max");
    //in inputs_outputs
    z_min.set_name("z_min");
    z_max.set_name("z_max");
    //in inputs
    cost.set_name("cost");
    avail_min.set_name("avail_min");
    avail_max.set_name("avail_max");
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

indices PoolNet::out_arcs_per_node(indices arcset) const{
    auto ids = indices(arcs);
    ids._type = matrix_;
    ids.set_name("out_arcs_per_node_in_"+arcset.get_name());
    ids._ids = make_shared<vector<vector<size_t>>>();
    ids._ids->resize(nodes.size());
    string key;
    size_t inst = 0;
    for (auto n: nodes) {
        for (auto a:*(arcset._keys)) {
                if(arcMap.find(a)!=arcMap.end()){
                auto arca= arcMap.find(a)->second;
                if(arca->_src->_name==n->_name){
                    ids._ids->at(inst).push_back(arca->_id);
                }
            }
        }
        inst++;
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
    
    string name;
    
    int N_input=5;
    int N_output=3;
    int N_pool=4;
    
    Inputs=indices("Inputs");
    Outputs=indices("Outputs");
    Pools=indices("Pools");

    
    for (auto i=0;i<N_input;i++)
    {
        Node* node = NULL;
        node = new Node(to_string(i));
        add_node(node);
        Inputs.add(to_string(i));
    }
    
    
    
    
    for (auto i=N_input;i<N_input+N_pool;i++)
    {
        Node* node = NULL;
        node = new Node(to_string(i));
        add_node(node);
        Pools.add(to_string(i));
    }
    
    for (auto i=N_input+N_pool;i<N_input+N_pool+N_output;i++)
    {
        Node* node = NULL;
        node = new Node(to_string(i));
        add_node(node);
        Outputs.add(to_string(i));
    }
    unsigned index = 0;
    string src, dest;
    
    
    inputs_pools=indices("inputs_pools");
    for (auto i=0;i<N_input;i++)
    {
        for (auto j=N_input;j<N_input+N_pool;j++)
        {
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
    
    for (auto i=N_input;i<N_input+N_pool;i++)
    {
        for (auto j=N_input+N_pool;j<N_input+N_pool+N_output;j++)
        {
            Arc* arc = NULL;
            
            
            
            arc = new Arc(to_string(i) + "," + to_string(j));
            arc->_id = index++;
            arc->_src = get_node(to_string(i));
            arc->_dest= get_node(to_string(j));
            this->add_arc(arc);
            arc->connect();
            pools_outputs.add(to_string(i) + "," +to_string(j));
        }
        
    }
    
    for (auto i=0;i<N_input;i++)
    {
        for (auto j=N_input+N_pool;j<N_input+N_pool+N_output;j++)
        {
            Arc* arc = NULL;
            
            
            
            arc = new Arc(to_string(i) + "," + to_string(j));
            arc->_id = index++;
            arc->_src = get_node(to_string(i));
            arc->_dest= get_node(to_string(j));
            this->add_arc(arc);
            arc->connect();
            inputs_outputs.add(to_string(i) + "," + to_string(j));
        }
        
    }
    
    for(auto key: *(inputs_pools._keys)){
        x_min.add_val(key, 0);
        x_max.add_val(key, 10);
        inqual.add_val(key, 0.5);
        
    }
    for(auto key: *(pools_outputs._keys)){
        y_min.add_val(key, 0);
        y_max.add_val(key, 10);
        
    }
    for(auto key: *(inputs_outputs._keys)){
        z_min.add_val(key, 0);
        z_max.add_val(key, 10);
        
    }
    

    
    for(auto key: *(Inputs._keys)){
        cost.add_val(key, 10);
        avail_min.add_val(key, 0);
        avail_max.add_val(key, 10);
        inqual.add_val(key, 0.5);
        
    }
    
    for(auto key: *(Outputs._keys)){
        rev.add_val(key, 10);
        dem_min.add_val(key, 1);
        dem_max.add_val(key, 100);
        outqual_min.add_val(key, 0.3);
        outqual_max.add_val(key, 0.7);
    }
    for(auto key: *(Pools._keys)){
        pool_cap.add_val(key, 50);
    }
    
}




