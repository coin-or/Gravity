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
    dem_min.set_name("dem_min");
    dem_max.set_name("dem_max");
    outqual_min.set_name("outqual_min");
    outqual_max.set_name("outqual_max");
    //in pools
    pool_cap.set_name("pool_cap");
    sumyk.set_name("sumyk");
    sumyk.in(range(0,0));
    sumyk.set_val(0);
    
 
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

indices PoolNet::inputs_pools_outputs() const{
    indices in_po_out("in_po_out");

    for(const string& ip: *this->inputs_pools._keys){
        auto i_node=ip.substr(0, ip.find_first_of(","));
        auto p_node=ip.substr(ip.find_first_of(",")+1);
        for(const string& po: *this->pools_outputs._keys){
            auto p_node1=po.substr(0, po.find_first_of(","));
            auto o_node=po.substr(po.find_first_of(",")+1);
            if(p_node==p_node1){
                in_po_out.add(ip+","+o_node);
            }
    }
    
    
}
    return in_po_out;
}



void PoolNet::readgrid(string fname) {
    
   // string fname=string(prj_dir)+"/data_sets/Pooling/Adhya1_gms.txt";
    string word="", tempwrd, numwrd1, numwrd;
    
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
    numwrd=tempwrd.substr(0, tempwrd.find_first_of("/"));

   
    N_node = atoi(numwrd.c_str());
    
    getline(file, word);
    
    tempwrd=word.substr(word.find_first_of("*")+1);
    numwrd=tempwrd.substr(0, tempwrd.find_first_of("/"));
    
    
    N_input = atoi(numwrd.c_str());

    getline(file, word);
    
    numwrd1=word.substr(word.find_first_of("/")+1,word.find_first_of("*"));
    tempwrd=word.substr(word.find_first_of("*")+1);
    numwrd=tempwrd.substr(0, tempwrd.find_first_of("/"));
    
    
    
    N_output = atoi(numwrd.c_str())-atoi(numwrd1.c_str())+1;
    
    
    getline(file, word);
    
    tempwrd=word.substr(word.find_first_of("*")+1);
    numwrd=tempwrd.substr(0, tempwrd.find_first_of("/"));
    
    
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
shared_ptr<Model<>> build_pool_qform(PoolNet& poolnet)
{
int output = 0;


//do bounds on x,y,z using preprocessign in paper!
//This is p-q-formulaiton of pooling problem!



auto SPP= make_shared<Model<>>("Std-Pooling-Prob-PQ");
indices Inputs=poolnet.Inputs;
indices Pools=poolnet.Pools;
indices Outputs=poolnet.Outputs;
indices Attr=poolnet.Attr;

//indices Nodes=pool.nodes;

indices inputs_pools=poolnet.inputs_pools;
indices pools_outputs=poolnet.pools_outputs;
indices inputs_outputs=poolnet.inputs_outputs;
indices inputs_attr=poolnet.inputs_attr;
indices outputs_attr=poolnet.outputs_attr;
indices pool_attr = indices(Pools,Attr);

auto out_arcs_to_pool_per_input = poolnet.out_arcs_to_pool_per_input();
auto out_arcs_to_output_per_input = poolnet.out_arcs_to_output_per_input();
auto in_arcs_per_pool = poolnet.in_arcs_per_pool();
auto in_arcs_per_pool_attr = poolnet.in_arcs_per_pool_attr();
auto in_arcs_attr_per_pool = poolnet.in_arcs_attr_per_pool();
auto out_arcs_per_pool = poolnet.out_arcs_per_pool();
auto out_arcs_per_pool_attr = poolnet.out_arcs_per_pool_attr();
auto in_arcs_from_pool_per_output = poolnet.in_arcs_from_pool_per_output();
auto in_arcs_from_pool_per_output_attr=poolnet.in_arcs_from_pool_per_output_attr();
auto in_arcs_from_input_per_output = poolnet.in_arcs_from_input_per_output();
auto in_arcs_from_input_per_output_attr = poolnet.in_arcs_from_input_per_output_attr();




auto y_min=poolnet.y_min.in(pools_outputs);
auto y_max=poolnet.y_max.in(pools_outputs);

auto z_min=poolnet.z_min.in(inputs_outputs);
auto z_max=poolnet.z_max.in(inputs_outputs);

// auto cost=poolnet.cost.in(Inputs);
auto avail_min=poolnet.avail_min.in(Inputs);
auto avail_max=poolnet.avail_max.in(Inputs);
auto p_in=poolnet.inqual.in(inputs_attr);

auto rev=poolnet.rev.in(Outputs);
auto dem_min=poolnet.dem_min.in(Outputs);
auto dem_max=poolnet.dem_max.in(Outputs);
auto p_out_min=poolnet.outqual_min.in(outputs_attr);
auto p_out_max=poolnet.outqual_max.in(outputs_attr);

auto pool_cap=poolnet.pool_cap.in(Pools);

auto cost_ip=poolnet.cost_ip.in(inputs_pools);
auto cost_io=poolnet.cost_io.in(inputs_outputs);
auto cost_po=poolnet.cost_po.in(pools_outputs);

var<> q("q", 0, 1);
var<> cq("cq", 0, 100);

var<> y("y", y_min, y_max), z("z", z_min, z_max);
var<> p_pool("p_pool", 0, 5);
SPP->add(q.in(inputs_pools));
SPP->add(y.in(pools_outputs));
SPP->add(z.in(inputs_outputs));
SPP->add(p_pool.in(pool_attr));
SPP->add(cq.in(range(0,1)));

q.initialize_all(0.5);


int row_id = 0;
indices ypo_per_input_matrix = indices("ypo_per_input_matrix");
for (const string& input_key:*Inputs._keys) {
    ypo_per_input_matrix.add_empty_row();
    for (auto &po:*pools_outputs._keys) {
        auto pos = nthOccurrence(po, ",", 1);
        auto pool = po.substr(0,pos);
        auto in_po=input_key+","+pool;
        if(inputs_pools._keys_map->find(in_po) !=inputs_pools._keys_map->end()){
            ypo_per_input_matrix.add_in_row(row_id, po);
        }
    }
    row_id++;
}

row_id = 0;
indices q_per_ypo_per_input_matrix = indices("q_per_ypo_per_input_matrix");
for (const string& input_key:*Inputs._keys) {
    q_per_ypo_per_input_matrix.add_empty_row();
    for (auto &po:*pools_outputs._keys) {
        auto pos = nthOccurrence(po, ",", 1);
        auto pool = po.substr(0,pos);
        auto in_po=input_key+","+pool;
        if(inputs_pools._keys_map->find(in_po) !=inputs_pools._keys_map->end()){
            q_per_ypo_per_input_matrix.add_in_row(row_id, in_po);
        }
    }
    row_id++;
}


Constraint<> avail_lb("avail_lb");
avail_lb=q.in(q_per_ypo_per_input_matrix)*y.in(ypo_per_input_matrix)+sum(z, out_arcs_to_output_per_input)-avail_min;
SPP->add(avail_lb.in(Inputs)>=0);
// SPP->print();



Constraint<> avail_ub("avail_ub");
avail_ub=q.in(q_per_ypo_per_input_matrix)*y.in(ypo_per_input_matrix)+sum(z, out_arcs_to_output_per_input)-avail_max;
SPP->add(avail_ub.in(Inputs)<=0);
//SPP->print();



Constraint<> pool_capacity("pool_capacity");
pool_capacity=sum(y, out_arcs_per_pool)-pool_cap;
SPP->add(pool_capacity.in(Pools)<=0);



Constraint<> demand_lb("demand_lb");
demand_lb=sum(y, in_arcs_from_pool_per_output)+sum(z,in_arcs_from_input_per_output)-dem_min;
SPP->add(demand_lb.in(Outputs)>=0);

Constraint<> demand_ub("demand_ub");
demand_ub=sum(y, in_arcs_from_pool_per_output)+sum(z,in_arcs_from_input_per_output)-dem_max;
SPP->add(demand_ub.in(Outputs)<=0);

Constraint<> mass_balance("mass_balance");
mass_balance=sum(q, in_arcs_per_pool)-1;
SPP->add(mass_balance.in(Pools)==0);



row_id = 0;
indices pool_matrix = indices("pool_matrix");
for (const string& pool_key:*pool_attr._keys) {
    auto pos = nthOccurrence(pool_key, ",", 1);
    auto pool = pool_key.substr(0,pos);
    pool_matrix.add_empty_row();
    for (auto &out:*Outputs._keys) {
        auto po_out=pool+","+out;
        if(pools_outputs._keys_map->find(po_out) !=pools_outputs._keys_map->end()){
            pool_matrix.add_in_row(row_id, pool_key);
        }
    }
    row_id++;
}





Constraint<> quality_balance("quality_balance");//TODO debug transpose version
quality_balance=p_in.in(in_arcs_attr_per_pool)*q.in(in_arcs_per_pool_attr) - p_pool;// - p_pool* sum(y, out_arcs_per_pool)
SPP->add(quality_balance.in(pool_attr)==0);



row_id = 0;
indices input_attr_per_output_attr_matrix = indices("input_attr_per_output_attr_matrix");
for (const string& outat_key:*outputs_attr._keys) {
    auto pos = nthOccurrence(outat_key, ",", 1);
    auto out=outat_key.substr(0, pos);
    auto attr = outat_key.substr(pos+1);
    input_attr_per_output_attr_matrix.add_empty_row();
    for (auto &input:*Inputs._keys) {
        
        auto key=input+","+attr;
        input_attr_per_output_attr_matrix.add_in_row(row_id, key);
        
    }
    row_id++;
}



row_id = 0;
indices pool_attr_per_output_attr_matrix = indices("pool_attr_per_output_attr_matrix");
for (const string& outat_key:*outputs_attr._keys) {
    auto pos = nthOccurrence(outat_key, ",", 1);
    auto out=outat_key.substr(0, pos);
    auto attr = outat_key.substr(pos+1);
    pool_attr_per_output_attr_matrix.add_empty_row();
    for (auto &pool:*Pools._keys) {
        auto po_out= pool+","+out;
        if(pools_outputs._keys_map->find(po_out) !=pools_outputs._keys_map->end()){
            auto key=pool+","+attr;
            pool_attr_per_output_attr_matrix.add_in_row(row_id, key);
        }
    }
    row_id++;
}



row_id = 0;
indices output_attr_per_ypo_matrix = indices("output_attr_per_ypo_matrix");
for (const string& outat_key:*outputs_attr._keys) {
    auto pos = nthOccurrence(outat_key, ",", 1);
    auto out = outat_key.substr(0,pos);
    output_attr_per_ypo_matrix.add_empty_row();
    for (auto &pool:*Pools._keys) {
        auto po_out=pool+","+out;
        if(pools_outputs._keys_map->find(po_out) !=pools_outputs._keys_map->end()){
            output_attr_per_ypo_matrix.add_in_row(row_id, outat_key);
        }
    }
    row_id++;
}

row_id = 0;
indices output_attr_per_zio_matrix = indices("output_attr_per_zio_matrix");
for (const string& outat_key:*outputs_attr._keys) {
    auto pos = nthOccurrence(outat_key, ",", 1);
    auto out = outat_key.substr(0,pos);
    output_attr_per_zio_matrix.add_empty_row();
    for (auto &input:*Inputs._keys) {
        auto io_out=input+","+out;
        if(inputs_outputs._keys_map->find(io_out) !=inputs_outputs._keys_map->end()){
            output_attr_per_zio_matrix.add_in_row(row_id, outat_key);
        }
    }
    row_id++;
}



//    Constraint<> product_quality_ub("product_quality_ub");
//    product_quality_ub=y.in(in_arcs_from_pool_per_output_attr)*p_pool.in(pool_attr_per_output_attr_matrix)+z.in(in_arcs_from_input_per_output_attr)*p_in.in(input_attr_per_output_attr_matrix)-p_out_max.in(output_attr_per_ypo_matrix)*y.in(in_arcs_from_pool_per_output_attr)-p_out_max.in(output_attr_per_zio_matrix)*z.in(in_arcs_from_input_per_output_attr);
//    SPP->add(product_quality_ub.in(outputs_attr)<=0);
//func<> a= (p_in.in(in_arcs_attr_per_pool)*q.in(in_arcs_per_pool_attr));

Constraint<> product_quality_ub("product_quality_ub");
product_quality_ub=y.in(in_arcs_from_pool_per_output_attr)*p_pool.in(pool_attr_per_output_attr_matrix)+z.in(in_arcs_from_input_per_output_attr)*p_in.in(input_attr_per_output_attr_matrix)-p_out_max.in(output_attr_per_ypo_matrix)*y.in(in_arcs_from_pool_per_output_attr)-p_out_max.in(output_attr_per_zio_matrix)*z.in(in_arcs_from_input_per_output_attr);
SPP->add(product_quality_ub.in(outputs_attr)<=0);










//Constraint<> sumy("sumy");
//sumy=sum(y);
// SPP->add_lazy(sumy>=11);

//row_id = 0;
//indices pool_per_output = indices("pool_per_output");
//for (auto &out:*Outputs._keys) {
//    pool_per_output.add_empty_row();
//    for (const string& pool_out:*pools_outputs._keys) {
//        auto pos = nthOccurrence(pool_out, ",", 1);
//        auto out1 = pool_out.substr(pos+1);
//
//        if(out1==out){
//            pool_per_output.add_in_row(row_id, pool_out);
//        }
//    }
//    row_id++;
//}

//Constraint<> costq("costq");
//costq=(cost_ip.in(in_arcs_per_pool)*q.in(in_arcs_per_pool))-cq;
//SPP->add(costq==0);

//SPP->min((cost_ip.in(q_per_ypo_per_input_matrix).tr()*q).tr()*y);
    SPP->min(q.tr()*cost_ip.in(q_per_ypo_per_input_matrix)*y);
//auto obj=(cost_ip.in(in_arcs_per_pool)*q.in(in_arcs_per_pool));
//auto obj1=(cost_ip.in(in_arcs_per_pool)*q.in(in_arcs_per_pool));

//    func<> a=(product(y.in(pool_per_output), obj1));
//    a.print();

//auto obj= y.in(pool_per_output)*(cq);
// obj1.eval_all();
// obj.print();
//
//   auto obj= product(cost_ip, )+product(cost_io, z)+product(cost_po, y);
//   SPP->min(obj);

    return SPP;

}
shared_ptr<Model<>> build_pool_pform(PoolNet& poolnet,  SolverType solv_type)
{
    //This is p-formulaiton of pooling problem, yet to explore, p-q and q!
    
    
    auto SPP= make_shared<Model<>>("Std-Pooling-Prob-P");
   
    indices Inputs=poolnet.Inputs;
    indices Pools=poolnet.Pools;
    indices Outputs=poolnet.Outputs;
    indices Attr=poolnet.Attr;
    
    //indices Nodes=pool.nodes;
    
    indices inputs_pools=poolnet.inputs_pools;
    indices pools_outputs=poolnet.pools_outputs;
    indices inputs_outputs=poolnet.inputs_outputs;
    indices inputs_attr=poolnet.inputs_attr;
    indices outputs_attr=poolnet.outputs_attr;
    indices pool_attr = indices(Pools,Attr);
    
    auto out_arcs_to_pool_per_input = poolnet.out_arcs_to_pool_per_input();
    auto out_arcs_to_output_per_input = poolnet.out_arcs_to_output_per_input();
    auto in_arcs_per_pool = poolnet.in_arcs_per_pool();
    auto in_arcs_per_pool_attr = poolnet.in_arcs_per_pool_attr();
    auto in_arcs_attr_per_pool = poolnet.in_arcs_attr_per_pool();
    auto out_arcs_per_pool = poolnet.out_arcs_per_pool();
    auto out_arcs_per_pool_attr = poolnet.out_arcs_per_pool_attr();
    auto in_arcs_from_pool_per_output = poolnet.in_arcs_from_pool_per_output();
    auto in_arcs_from_pool_per_output_attr=poolnet.in_arcs_from_pool_per_output_attr();
    auto in_arcs_from_input_per_output = poolnet.in_arcs_from_input_per_output();
    auto in_arcs_from_input_per_output_attr = poolnet.in_arcs_from_input_per_output_attr();
    
    auto x_min=poolnet.x_min.in(inputs_pools);
    auto x_max=poolnet.x_max.in(inputs_pools);
    
    
    auto y_min=poolnet.y_min.in(pools_outputs);
    auto y_max=poolnet.y_max.in(pools_outputs);
    
    auto z_min=poolnet.z_min.in(inputs_outputs);
    auto z_max=poolnet.z_max.in(inputs_outputs);
    
    // auto cost=poolnet.cost.in(Inputs);
    auto avail_min=poolnet.avail_min.in(Inputs);
    auto avail_max=poolnet.avail_max.in(Inputs);
    auto p_in=poolnet.inqual.in(inputs_attr);
    
    auto rev=poolnet.rev.in(Outputs);
    auto dem_min=poolnet.dem_min.in(Outputs);
    auto dem_max=poolnet.dem_max.in(Outputs);
    auto p_out_min=poolnet.outqual_min.in(outputs_attr);
    auto p_out_max=poolnet.outqual_max.in(outputs_attr);
    
    auto pool_cap=poolnet.pool_cap.in(Pools);
    
    auto cost_ip=poolnet.cost_ip.in(inputs_pools);
    auto cost_io=poolnet.cost_io.in(inputs_outputs);
    auto cost_po=poolnet.cost_po.in(pools_outputs);
    auto sumyk=poolnet.sumyk;
    
    var<> x("x", x_min, x_max);

    
    
    var<> y("y", y_min, y_max), z("z", z_min, z_max);
    var<> p_pool("p_pool", 0, 5);
    SPP->add(x.in(inputs_pools));
    SPP->add(y.in(pools_outputs));
    SPP->add(z.in(inputs_outputs));
    SPP->add(p_pool.in(pool_attr));
//    SPP->add(sumyk);
//    sumyk.set_lb(0);
    x.initialize_all(2.0);
    y.initialize_all(2.0);
    
    Constraint<> avail_lb("avail_lb");
    avail_lb=sum(x, out_arcs_to_pool_per_input)+sum(z, out_arcs_to_output_per_input)-avail_min;
    //    avail_lb=sum(x, out_arcs_to_pool_per_input)-avail_min;
    SPP->add(avail_lb.in(Inputs)>=0);
   // SPP->print();
    
    
    Constraint<> avail_ub("avail_ub");
    avail_ub=sum(x, out_arcs_to_pool_per_input)+sum(z, out_arcs_to_output_per_input)-avail_max;
    SPP->add(avail_ub.in(Inputs)<=0);
    
    
    Constraint<> pool_capacity("pool_capacity");
    pool_capacity=sum(x, in_arcs_per_pool)-pool_cap;
    SPP->add(pool_capacity.in(Pools)<=0);
    
    
    Constraint<> demand_lb("demand_lb");
    demand_lb=sum(y, in_arcs_from_pool_per_output)+sum(z,in_arcs_from_input_per_output)-dem_min;
    SPP->add(demand_lb.in(Outputs)>=0);
    
    Constraint<> demand_ub("demand_ub");
    demand_ub=sum(y, in_arcs_from_pool_per_output)+sum(z,in_arcs_from_input_per_output)-dem_max;
    SPP->add(demand_ub.in(Outputs)<=0);
    
    Constraint<> mass_balance("mass_balance");
    mass_balance=sum(x, in_arcs_per_pool)-sum(y, out_arcs_per_pool);
    SPP->add(mass_balance.in(Pools)==0);
    
    int row_id = 0;
    indices pool_matrix = indices("pool_matrix");
    for (const string& pool_key:*pool_attr._keys) {
        auto pos = nthOccurrence(pool_key, ",", 1);
        auto pool = pool_key.substr(0,pos);
        pool_matrix.add_empty_row();
        for (auto &out:*Outputs._keys) {
            auto po_out=pool+","+out;
            if(pools_outputs._keys_map->find(po_out) !=pools_outputs._keys_map->end()){
                pool_matrix.add_in_row(row_id, pool_key);
            }
        }
        row_id++;
    }
    
    if(solv_type==gurobi){
    
        
    Constraint<> quality_balance_le("quality_balance_le");//TODO debug transpose version
    quality_balance_le=p_in.in(in_arcs_attr_per_pool)*x.in(in_arcs_per_pool_attr) - p_pool.in(pool_matrix)*y.in(out_arcs_per_pool_attr);// - p_pool* sum(y, out_arcs_per_pool)
    SPP->add(quality_balance_le.in(pool_attr)<=0);
    
    Constraint<> quality_balance_ge("quality_balance_ge");//TODO debug transpose version
    quality_balance_ge=p_in.in(in_arcs_attr_per_pool)*x.in(in_arcs_per_pool_attr) - p_pool.in(pool_matrix)*y.in(out_arcs_per_pool_attr);// - p_pool* sum(y, out_arcs_per_pool)
    SPP->add(quality_balance_ge.in(pool_attr)>=0);
    //SPP->print();
    }
    else{
        Constraint<> quality_balance_le("quality_balance_le");//TODO debug transpose version
        quality_balance_le=p_in.in(in_arcs_attr_per_pool)*x.in(in_arcs_per_pool_attr) - p_pool.in(pool_matrix)*y.in(out_arcs_per_pool_attr);// - p_pool* sum(y, out_arcs_per_pool)
        SPP->add(quality_balance_le.in(pool_attr)<=0);
        
        Constraint<> quality_balance_ge("quality_balance_ge");//TODO debug transpose version
        quality_balance_ge=p_in.in(in_arcs_attr_per_pool)*x.in(in_arcs_per_pool_attr) - p_pool.in(pool_matrix)*y.in(out_arcs_per_pool_attr);// - p_pool* sum(y, out_arcs_per_pool)
        SPP->add(quality_balance_ge.in(pool_attr)>=0);
        
    }
        
    
    
    
    
    row_id = 0;
    indices pool_attr_per_output_attr_matrix = indices("pool_attr_per_output_attr_matrix");
    for (const string& outat_key:*outputs_attr._keys) {
        auto pos = nthOccurrence(outat_key, ",", 1);
        auto out=outat_key.substr(0, pos);
        auto attr = outat_key.substr(pos+1);
        pool_attr_per_output_attr_matrix.add_empty_row();
        for (auto &pool:*Pools._keys) {
            auto po_out= pool+","+out;
            if(pools_outputs._keys_map->find(po_out) !=pools_outputs._keys_map->end()){
                auto key=pool+","+attr;
                pool_attr_per_output_attr_matrix.add_in_row(row_id, key);
            }
        }
        row_id++;
    }
    
    row_id = 0;
    indices input_attr_per_output_attr_matrix = indices("input_attr_per_output_attr_matrix");
    for (const string& outat_key:*outputs_attr._keys) {
        auto pos = nthOccurrence(outat_key, ",", 1);
        auto out=outat_key.substr(0, pos);
        auto attr = outat_key.substr(pos+1);
        input_attr_per_output_attr_matrix.add_empty_row();
        for (auto &input:*Inputs._keys) {
            auto io_out= input+","+out;
            if(inputs_outputs._keys_map->find(io_out) !=inputs_outputs._keys_map->end()){
                auto key=input+","+attr;
                input_attr_per_output_attr_matrix.add_in_row(row_id, key);
            }
        }
        row_id++;
    }
    
    row_id = 0;
    indices output_attr_per_ypo_matrix = indices("output_attr_per_ypo_matrix");
    for (const string& outat_key:*outputs_attr._keys) {
        auto pos = nthOccurrence(outat_key, ",", 1);
        auto out = outat_key.substr(0,pos);
        output_attr_per_ypo_matrix.add_empty_row();
        for (auto &pool:*Pools._keys) {
            auto po_out=pool+","+out;
            if(pools_outputs._keys_map->find(po_out) !=pools_outputs._keys_map->end()){
                output_attr_per_ypo_matrix.add_in_row(row_id, outat_key);
            }
        }
        row_id++;
    }
    
    row_id = 0;
    indices output_attr_per_zio_matrix = indices("output_attr_per_zio_matrix");
    for (const string& outat_key:*outputs_attr._keys) {
        auto pos = nthOccurrence(outat_key, ",", 1);
        auto out = outat_key.substr(0,pos);
        output_attr_per_zio_matrix.add_empty_row();
        for (auto &input:*Inputs._keys) {
            auto io_out=input+","+out;
            if(inputs_outputs._keys_map->find(io_out) !=inputs_outputs._keys_map->end()){
                output_attr_per_zio_matrix.add_in_row(row_id, outat_key);
            }
        }
        row_id++;
    }
    
    
    //
    //    Constraint<> product_quality_lb("product_quality_lb");//TODO debug transpose version and propagate matrix indexing to function
    //    product_quality_lb=y.in(in_arcs_from_pool_per_output)*p_pool+z.in(in_arcs_from_input_per_output)*(p_in-p_out_min.in(outinput_matrix)).in(outinput_matrix) - p_out_min.in(outpool_matrix)*y.in(in_arcs_from_pool_per_output);//-p_out_min.in(outinput_matrix)*z.in(in_arcs_from_input_per_output);
    //    SPP->add(product_quality_lb.in(Outputs)>=0);
    //
    //    Constraint<> product_quality_ub("product_quality_ub");//TODO debug transpose version and propagate matrix indexing to function
    //    product_quality_ub=y.in(in_arcs_from_pool_per_output)*p_pool+z.in(in_arcs_from_input_per_output)*(p_in-p_out_max.in(outinput_matrix)).in(outinput_matrix) - p_out_max.in(outpool_matrix)*y.in(in_arcs_from_pool_per_output);//-p_out_min.in(outinput_matrix)*z.in(in_arcs_from_input_per_output);
    //    SPP->add(product_quality_ub.in(Outputs)<=0);
    
    
//    Constraint<> product_quality_ub_le("product_quality_ub_le");//TODO debug transpose version and propagate matrix indexing to function
//    product_quality_ub_le=y.in(in_arcs_from_pool_per_output_attr)*p_pool.in(pool_attr_per_output_attr_matrix)+z.in(in_arcs_from_input_per_output_attr)*p_in.in(input_attr_per_output_attr_matrix)-p_out_max.in(output_attr_per_ypo_matrix)*y.in(in_arcs_from_pool_per_output_attr)-p_out_max.in(output_attr_per_zio_matrix)*z.in(in_arcs_from_input_per_output_attr);
//    SPP->add(product_quality_ub_le.in(outputs_attr)<=0);
    
    Constraint<> product_quality_ub("product_quality_ub"); product_quality_ub=y.in(in_arcs_from_pool_per_output_attr)*p_pool.in(pool_attr_per_output_attr_matrix)+z.in(in_arcs_from_input_per_output_attr)*p_in.in(input_attr_per_output_attr_matrix)-p_out_max.in(output_attr_per_ypo_matrix)*y.in(in_arcs_from_pool_per_output_attr)-p_out_max.in(output_attr_per_zio_matrix)*z.in(in_arcs_from_input_per_output_attr);
    SPP->add(product_quality_ub.in(outputs_attr)<=0);
    
    Constraint<> sumy_con("sumy_con");
    sumy_con=sum(y)-sumyk;
    SPP->add(sumy_con.in(range(0,0))>=0);
    
    auto obj= product(cost_ip, x)+product(cost_io, z)+product(cost_po, y);
    SPP->min(obj);
    
    return SPP;
}
shared_ptr<Model<>> build_pool_pqform(PoolNet& poolnet,  SolverType solv_type)
{
    //This is pq-formulaiton of pooling problem as written in DOI 10.1007/s10898-012-9875-6!
    
    
    auto SPP= make_shared<Model<>>("Std-Pooling-Prob-P");
    
    indices Inputs=poolnet.Inputs;
    indices Pools=poolnet.Pools;
    indices Outputs=poolnet.Outputs;
    indices Attr=poolnet.Attr;
    
    //indices Nodes=pool.nodes;
    
    
    indices inputs_pools=poolnet.inputs_pools;
    indices pools_outputs=poolnet.pools_outputs;
    indices inputs_outputs=poolnet.inputs_outputs;
    indices inputs_attr=poolnet.inputs_attr;
    indices outputs_attr=poolnet.outputs_attr;
    indices pool_attr = indices(Pools,Attr);
    
    indices inputs_pools_outputs=poolnet.inputs_pools_outputs();
    auto out_arcs_to_pool_per_input = poolnet.out_arcs_to_pool_per_input();
    auto out_arcs_to_output_per_input = poolnet.out_arcs_to_output_per_input();
    auto in_arcs_per_pool = poolnet.in_arcs_per_pool();
    auto in_arcs_per_pool_attr = poolnet.in_arcs_per_pool_attr();
    auto in_arcs_attr_per_pool = poolnet.in_arcs_attr_per_pool();
    auto out_arcs_per_pool = poolnet.out_arcs_per_pool();
    auto out_arcs_per_pool_attr = poolnet.out_arcs_per_pool_attr();
    auto in_arcs_from_pool_per_output = poolnet.in_arcs_from_pool_per_output();
    auto in_arcs_from_pool_per_output_attr=poolnet.in_arcs_from_pool_per_output_attr();
    auto in_arcs_from_input_per_output = poolnet.in_arcs_from_input_per_output();
    auto in_arcs_from_input_per_output_attr = poolnet.in_arcs_from_input_per_output_attr();
    

    
    // auto cost=poolnet.cost.in(Inputs);
    auto avail_min=poolnet.avail_min.in(Inputs);
    auto avail_max=poolnet.avail_max.in(Inputs);
    auto p_in=poolnet.inqual.in(inputs_attr);
    
    auto rev=poolnet.rev.in(Outputs);
    auto dem_min=poolnet.dem_min.in(Outputs);
    auto dem_max=poolnet.dem_max.in(Outputs);
    auto p_out_min=poolnet.outqual_min.in(outputs_attr);
    auto p_out_max=poolnet.outqual_max.in(outputs_attr);
    
    auto pool_cap=poolnet.pool_cap.in(Pools);
    
    auto cost_ip=poolnet.cost_ip.in(inputs_pools);
    auto cost_io=poolnet.cost_io.in(inputs_outputs);
    auto cost_po=poolnet.cost_po.in(pools_outputs);
    auto sumyk=poolnet.sumyk;
    
//    auto x_min=poolnet.x_min.in(inputs_pools_outputs);
//    auto x_max=poolnet.x_max.in(inputs_pools_outputs);
//
    
//    auto y_min=poolnet.y_min.in(pools_outputs);
//    auto y_max=poolnet.y_max.in(pools_outputs);
//
    auto z_min=poolnet.z_min.in(inputs_outputs);
    auto z_max=poolnet.z_max.in(inputs_outputs);
    
    var<> x("x",0, 100), y("y", 0, 100);
    var<> q("q", 0, 1), z("z", z_min, z_max);
  
    SPP->add(x.in(inputs_pools_outputs));
    SPP->add(q.in(inputs_pools));
    SPP->add(z.in(inputs_outputs));
    SPP->add(y.in(pools_outputs));

    //    SPP->add(sumyk);
    //    sumyk.set_lb(0);
    x.initialize_all(2.0);
    q.initialize_all(0.5);
    
        int row_id = 0;
        indices input_x_matrix = indices("input_x_matrix");
        for (const string& input_key:*Inputs._keys) {
          
            input_x_matrix.add_empty_row();
            for (auto &ipo:*inputs_pools_outputs._keys) {
                
                auto pos = nthOccurrence(ipo, ",", 1);
                auto input1 = ipo.substr(0,pos);
                
                if(input_key==input1){
                    input_x_matrix.add_in_row(row_id, ipo);
                }
            }
            row_id++;
        }
    

    
    
    Constraint<> avail_ub("avail_ub");
  //  avail_lb=sum(x, out_arcs_to_pool_per_input)+sum(z, out_arcs_to_output_per_input)-avail_min;
        avail_ub=sum(x, input_x_matrix)+sum(z, out_arcs_to_output_per_input)-avail_max;;
    SPP->add(avail_ub.in(Inputs)<=0);
    
    row_id = 0;
    indices pool_x_matrix = indices("pool_x_matrix");
    for (const string& pool_key:*Pools._keys) {
        
        pool_x_matrix.add_empty_row();
        for (auto &ipo:*inputs_pools_outputs._keys) {
            
            auto pos = nthOccurrence(ipo, ",", 1);
            auto pos1 = nthOccurrence(ipo, ",", 2);
            auto pool1 = ipo.substr(pos+1, pos1-(pos+1));
            
            if(pool_key==pool1){
                pool_x_matrix.add_in_row(row_id, ipo);
            }
        }
        row_id++;
    }
    
        Constraint<> pool_capacity("pool_capacity");
        pool_capacity=sum(x, pool_x_matrix)-pool_cap;
        SPP->add(pool_capacity.in(Pools)<=0);

    
    
    row_id = 0;
    indices output_x_matrix = indices("pool_x_matrix");
    for (const string& out_key:*Outputs._keys) {
        
        output_x_matrix.add_empty_row();
        for (auto &ipo:*inputs_pools_outputs._keys) {
            
            auto pos = nthOccurrence(ipo, ",", 2);
            auto pos1 = nthOccurrence(ipo, ",", 3);
            auto out1 = ipo.substr(pos+1, pos1-(pos+1));
            
            if(out_key==out1){
                output_x_matrix.add_in_row(row_id, ipo);
            }
        }
        row_id++;
    }
    
    
    Constraint<> demand_ub("demand_ub");
    demand_ub=sum(x, output_x_matrix)+sum(z,in_arcs_from_input_per_output)-dem_max;
    SPP->add(demand_ub.in(Outputs)<=0);

    SPP->print();
    
    row_id = 0;
    indices outattr_x_matrix = indices("outattr_x_matrix");
    for (const string& outat_key:*outputs_attr._keys) {
        
        outattr_x_matrix.add_empty_row();
        auto out_key = outat_key.substr(0, outat_key.find_first_of(","));
        for (auto &ipo:*inputs_pools_outputs._keys) {
            
            auto pos = nthOccurrence(ipo, ",", 2);
            auto pos1 = nthOccurrence(ipo, ",", 3);
            auto out1 = ipo.substr(pos+1, pos1-(pos+1));
            
            if(out_key==out1){
                outattr_x_matrix.add_in_row(row_id, ipo);
            }
        }
        row_id++;
    }
    
    row_id = 0;
    indices outattr_pin_matrix = indices("outattr_pin_matrix");
    for (const string& outat_key:*outputs_attr._keys) {
        
        outattr_pin_matrix.add_empty_row();
        auto out_key = outat_key.substr(0, outat_key.find_first_of(","));
        auto at_key = outat_key.substr(outat_key.find_first_of(",")+1);
        for (auto &ipo:*inputs_pools_outputs._keys) {
            auto pos1 = nthOccurrence(ipo, ",", 1);
            auto pos2 = nthOccurrence(ipo, ",", 2);
            auto pos3 = nthOccurrence(ipo, ",", 3);
            auto out1 = ipo.substr(pos2+1, pos3-(pos2+1));
            auto in=ipo.substr(0, pos1);
            if(out_key==out1){
                outattr_pin_matrix.add_in_row(row_id, in+","+at_key);
            }
        }
        row_id++;
    }
    
    row_id = 0;
    indices outattr_pout_matrix = indices("outattr_pout_matrix");
    for (const string& outat_key:*outputs_attr._keys) {
        
        outattr_pout_matrix.add_empty_row();
        auto out_key = outat_key.substr(0, outat_key.find_first_of(","));
        auto at_key = outat_key.substr(outat_key.find_first_of(",")+1);
        for (auto &ipo:*inputs_pools_outputs._keys) {
            auto pos1 = nthOccurrence(ipo, ",", 1);
            auto pos2 = nthOccurrence(ipo, ",", 2);
            auto pos3 = nthOccurrence(ipo, ",", 3);
            auto out1 = ipo.substr(pos2+1, pos3-(pos2+1));
            if(out_key==out1){
                outattr_pout_matrix.add_in_row(row_id, outat_key);
            }
        }
        row_id++;
    }
    
    row_id = 0;
    indices outattrz_pin_matrix = indices("outattrz_pin_matrix");
    for (const string& outat_key:*outputs_attr._keys) {
        
        outattrz_pin_matrix.add_empty_row();
        auto out_key = outat_key.substr(0, outat_key.find_first_of(","));
        auto at_key = outat_key.substr(outat_key.find_first_of(",")+1);
        for (auto &ipo:*inputs_outputs._keys) {
            auto pos1 = nthOccurrence(ipo, ",", 1);
            auto out1 = ipo.substr(pos1+1);
            auto in=ipo.substr(0, pos1);
            if(out_key==out1){
               outattrz_pin_matrix.add_in_row(row_id, in+","+at_key);
            }
        }
        row_id++;
    }
    
    row_id = 0;
    indices outattrz_pout_matrix = indices("outattrz_pout_matrix");
    for (const string& outat_key:*outputs_attr._keys) {
        
        outattrz_pout_matrix.add_empty_row();
        auto out_key = outat_key.substr(0, outat_key.find_first_of(","));
        auto at_key = outat_key.substr(outat_key.find_first_of(",")+1);
        for (auto &ipo:*inputs_outputs._keys) {
            auto pos1 = nthOccurrence(ipo, ",", 1);
            auto out1=ipo.substr(pos1+1);
            if(out_key==out1){
                outattrz_pout_matrix.add_in_row(row_id, outat_key);
            }
        }
        row_id++;
    }

     Constraint<> quality_balance_le("quality_balance_le");//TODO debug transpose version
    quality_balance_le=(p_in.in(outattr_pin_matrix)-p_out_max.in(outattr_pout_matrix)).in(outattr_x_matrix)*x.in(outattr_x_matrix)+(p_in.in(outattrz_pin_matrix)-p_out_max.in(outattrz_pout_matrix)).in(in_arcs_from_input_per_output_attr)*z.in(in_arcs_from_input_per_output_attr);
    SPP->add(quality_balance_le.in(outputs_attr)<=0);
            SPP->print();
    
    row_id = 0;
    indices pool_q_matrix = indices("pool_q_matrix");
    for (const string& pool_key:*Pools._keys) {
        
        pool_q_matrix.add_empty_row();
        for (auto &ip:*inputs_pools._keys) {
            
            auto pos = nthOccurrence(ip, ",", 1);
            auto pool1 = ip.substr(pos+1);
            
            if(pool_key==pool1){
                pool_q_matrix.add_in_row(row_id, ip);
            }
        }
        row_id++;
    }
    
    Constraint<> simplex("simplex");//TODO debug transpose version
    simplex=q.in(pool_q_matrix)-1;
    SPP->add(simplex.in(Pools)==0);
    SPP->print();
    
    row_id = 0;
    indices pooloutput_x_matrix = indices("pooloutput_x_matrix");
    for (const string& poolout_key:*pools_outputs._keys) {
        
        pooloutput_x_matrix.add_empty_row();
        for (auto &ipo:*inputs_pools_outputs._keys) {
            
            auto pos = nthOccurrence(ipo, ",", 1);
            auto poolout1 = ipo.substr(pos+1);
            
            if(poolout_key==poolout1){
                pooloutput_x_matrix.add_in_row(row_id, ipo);
            }
        }
        row_id++;
    }
    
    
    Constraint<> PQ("PQ");//TODO debug transpose version
    PQ=x.in(pooloutput_x_matrix)-y;
    SPP->add(PQ.in(pools_outputs)==0);
    SPP->print();
    
    row_id = 0;
    indices inputpool_x_matrix = indices("inputpool_x_matrix");
    for (const string& inputpool_key:*inputs_pools._keys) {
        
        inputpool_x_matrix.add_empty_row();
        for (auto &ipo:*inputs_pools_outputs._keys) {
            
            auto pos = nthOccurrence(ipo, ",", 2);
            auto inputpool1 = ipo.substr(0, pos);
            
            if(inputpool_key==inputpool1){
                inputpool_x_matrix.add_in_row(row_id, ipo);
            }
        }
        row_id++;
    }
    row_id = 0;
    indices inputpool_q_matrix = indices("inputpool_q_matrix");
    indices inputpool_poolcap_matrix = indices("inputpool_poolcap_matrix");
    for (const string& inputpool_key:*inputs_pools._keys) {
        
        inputpool_q_matrix.add_empty_row();
        inputpool_poolcap_matrix.add_empty_row();
        inputpool_q_matrix.add_in_row(row_id, inputpool_key);
        auto pos=nthOccurrence(inputpool_key, ",", 1);
        auto pool_key=inputpool_key.substr(pos+1);
        inputpool_poolcap_matrix.add_in_row(row_id, pool_key);
        row_id++;
    }
        
        
    
    
    
    
    Constraint<> PQ1("PQ1");//TODO debug transpose version
    PQ1=x.in(inputpool_x_matrix)-q.in(inputpool_q_matrix)*pool_cap.in(inputpool_poolcap_matrix);
    SPP->add(PQ1.in(inputs_pools)<=0);
    SPP->print();
    
    row_id = 0;
    indices inpoolout_y_matrix = indices("inpoolout_y_matrix");
    indices inpoolout_q_matrix = indices("inpoolout_q_matrix");
    for (const string& inpoout_key:*inputs_pools_outputs._keys) {
            auto pos = nthOccurrence(inpoout_key, ",", 1);
            auto poout=inpoout_key.substr(pos+1);
            auto pos1 = nthOccurrence(inpoout_key, ",", 2);
            auto inpo=inpoout_key.substr(0,pos1);
            inpoolout_y_matrix.add_in_row(row_id, poout);
            inpoolout_q_matrix.add_in_row(row_id, inpo);
            row_id++;
    }
    
        Constraint<> mass_balance("mass_balance");
        mass_balance=x.in(inputs_pools_outputs)-q.in(inpoolout_q_matrix)*y.in(inpoolout_y_matrix);
        SPP->add(mass_balance.in(inputs_pools_outputs)==0);
    
    SPP->print();
    
    
//
//    Constraint<> mass_balance("mass_balance");
//    mass_balance=sum(x, in_arcs_per_pool)-sum(y, out_arcs_per_pool);
//    SPP->add(mass_balance.in(Pools)==0);
//    
//    int row_id = 0;
//    indices pool_matrix = indices("pool_matrix");
//    for (const string& pool_key:*pool_attr._keys) {
//        auto pos = nthOccurrence(pool_key, ",", 1);
//        auto pool = pool_key.substr(0,pos);
//        pool_matrix.add_empty_row();
//        for (auto &out:*Outputs._keys) {
//            auto po_out=pool+","+out;
//            if(pools_outputs._keys_map->find(po_out) !=pools_outputs._keys_map->end()){
//                pool_matrix.add_in_row(row_id, pool_key);
//            }
//        }
//        row_id++;
//    }
//    
//    if(solv_type==gurobi){
//        
//        
//        Constraint<> quality_balance_le("quality_balance_le");//TODO debug transpose version
//        quality_balance_le=p_in.in(in_arcs_attr_per_pool)*x.in(in_arcs_per_pool_attr) - p_pool.in(pool_matrix)*y.in(out_arcs_per_pool_attr);// - p_pool* sum(y, out_arcs_per_pool)
//        SPP->add(quality_balance_le.in(pool_attr)<=0);
//        
//        Constraint<> quality_balance_ge("quality_balance_ge");//TODO debug transpose version
//        quality_balance_ge=p_in.in(in_arcs_attr_per_pool)*x.in(in_arcs_per_pool_attr) - p_pool.in(pool_matrix)*y.in(out_arcs_per_pool_attr);// - p_pool* sum(y, out_arcs_per_pool)
//        SPP->add(quality_balance_ge.in(pool_attr)>=0);
//        //SPP->print();
//    }
//    else{
//        Constraint<> quality_balance_le("quality_balance_le");//TODO debug transpose version
//        quality_balance_le=p_in.in(in_arcs_attr_per_pool)*x.in(in_arcs_per_pool_attr) - p_pool.in(pool_matrix)*y.in(out_arcs_per_pool_attr);// - p_pool* sum(y, out_arcs_per_pool)
//        SPP->add(quality_balance_le.in(pool_attr)<=0);
//        
//        Constraint<> quality_balance_ge("quality_balance_ge");//TODO debug transpose version
//        quality_balance_ge=p_in.in(in_arcs_attr_per_pool)*x.in(in_arcs_per_pool_attr) - p_pool.in(pool_matrix)*y.in(out_arcs_per_pool_attr);// - p_pool* sum(y, out_arcs_per_pool)
//        SPP->add(quality_balance_ge.in(pool_attr)>=0);
//        
//    }
//    
//    
//    
//    
//    
//    row_id = 0;
//    indices pool_attr_per_output_attr_matrix = indices("pool_attr_per_output_attr_matrix");
//    for (const string& outat_key:*outputs_attr._keys) {
//        auto pos = nthOccurrence(outat_key, ",", 1);
//        auto out=outat_key.substr(0, pos);
//        auto attr = outat_key.substr(pos+1);
//        pool_attr_per_output_attr_matrix.add_empty_row();
//        for (auto &pool:*Pools._keys) {
//            auto po_out= pool+","+out;
//            if(pools_outputs._keys_map->find(po_out) !=pools_outputs._keys_map->end()){
//                auto key=pool+","+attr;
//                pool_attr_per_output_attr_matrix.add_in_row(row_id, key);
//            }
//        }
//        row_id++;
//    }
//    
//    row_id = 0;
//    indices input_attr_per_output_attr_matrix = indices("input_attr_per_output_attr_matrix");
//    for (const string& outat_key:*outputs_attr._keys) {
//        auto pos = nthOccurrence(outat_key, ",", 1);
//        auto out=outat_key.substr(0, pos);
//        auto attr = outat_key.substr(pos+1);
//        input_attr_per_output_attr_matrix.add_empty_row();
//        for (auto &input:*Inputs._keys) {
//            auto io_out= input+","+out;
//            if(inputs_outputs._keys_map->find(io_out) !=inputs_outputs._keys_map->end()){
//                auto key=input+","+attr;
//                input_attr_per_output_attr_matrix.add_in_row(row_id, key);
//            }
//        }
//        row_id++;
//    }
//    
//    row_id = 0;
//    indices output_attr_per_ypo_matrix = indices("output_attr_per_ypo_matrix");
//    for (const string& outat_key:*outputs_attr._keys) {
//        auto pos = nthOccurrence(outat_key, ",", 1);
//        auto out = outat_key.substr(0,pos);
//        output_attr_per_ypo_matrix.add_empty_row();
//        for (auto &pool:*Pools._keys) {
//            auto po_out=pool+","+out;
//            if(pools_outputs._keys_map->find(po_out) !=pools_outputs._keys_map->end()){
//                output_attr_per_ypo_matrix.add_in_row(row_id, outat_key);
//            }
//        }
//        row_id++;
//    }
//    
//    row_id = 0;
//    indices output_attr_per_zio_matrix = indices("output_attr_per_zio_matrix");
//    for (const string& outat_key:*outputs_attr._keys) {
//        auto pos = nthOccurrence(outat_key, ",", 1);
//        auto out = outat_key.substr(0,pos);
//        output_attr_per_zio_matrix.add_empty_row();
//        for (auto &input:*Inputs._keys) {
//            auto io_out=input+","+out;
//            if(inputs_outputs._keys_map->find(io_out) !=inputs_outputs._keys_map->end()){
//                output_attr_per_zio_matrix.add_in_row(row_id, outat_key);
//            }
//        }
//        row_id++;
//    }
//    
//    
//    //
//    //    Constraint<> product_quality_lb("product_quality_lb");//TODO debug transpose version and propagate matrix indexing to function
//    //    product_quality_lb=y.in(in_arcs_from_pool_per_output)*p_pool+z.in(in_arcs_from_input_per_output)*(p_in-p_out_min.in(outinput_matrix)).in(outinput_matrix) - p_out_min.in(outpool_matrix)*y.in(in_arcs_from_pool_per_output);//-p_out_min.in(outinput_matrix)*z.in(in_arcs_from_input_per_output);
//    //    SPP->add(product_quality_lb.in(Outputs)>=0);
//    //
//    //    Constraint<> product_quality_ub("product_quality_ub");//TODO debug transpose version and propagate matrix indexing to function
//    //    product_quality_ub=y.in(in_arcs_from_pool_per_output)*p_pool+z.in(in_arcs_from_input_per_output)*(p_in-p_out_max.in(outinput_matrix)).in(outinput_matrix) - p_out_max.in(outpool_matrix)*y.in(in_arcs_from_pool_per_output);//-p_out_min.in(outinput_matrix)*z.in(in_arcs_from_input_per_output);
//    //    SPP->add(product_quality_ub.in(Outputs)<=0);
//    
//    
//    //    Constraint<> product_quality_ub_le("product_quality_ub_le");//TODO debug transpose version and propagate matrix indexing to function
//    //    product_quality_ub_le=y.in(in_arcs_from_pool_per_output_attr)*p_pool.in(pool_attr_per_output_attr_matrix)+z.in(in_arcs_from_input_per_output_attr)*p_in.in(input_attr_per_output_attr_matrix)-p_out_max.in(output_attr_per_ypo_matrix)*y.in(in_arcs_from_pool_per_output_attr)-p_out_max.in(output_attr_per_zio_matrix)*z.in(in_arcs_from_input_per_output_attr);
//    //    SPP->add(product_quality_ub_le.in(outputs_attr)<=0);
//    
//    Constraint<> product_quality_ub("product_quality_ub"); product_quality_ub=y.in(in_arcs_from_pool_per_output_attr)*p_pool.in(pool_attr_per_output_attr_matrix)+z.in(in_arcs_from_input_per_output_attr)*p_in.in(input_attr_per_output_attr_matrix)-p_out_max.in(output_attr_per_ypo_matrix)*y.in(in_arcs_from_pool_per_output_attr)-p_out_max.in(output_attr_per_zio_matrix)*z.in(in_arcs_from_input_per_output_attr);
//    SPP->add(product_quality_ub.in(outputs_attr)<=0);
//    
//    Constraint<> sumy_con("sumy_con");
//    sumy_con=sum(y)-sumyk;
//    SPP->add(sumy_con.in(range(0,0))>=0);
//    
//    auto obj= product(cost_ip, x)+product(cost_io, z)+product(cost_po, y);
//    SPP->min(obj);
//    
    return SPP;
}

