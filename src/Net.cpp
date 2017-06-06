//
//  Net.cpp
//
//
//  Created by Guagnlei on 03/06/2017.
//

#include "Net.h"
#include "Node.h"
#include "Path.h"
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

using namespace std;

Net::Net(){
    //horton_net = nullptr;
    _clone = nullptr;
    _bags = nullptr;
}

/* returns true if an arc is already present between the given nodes */
bool Net::duplicate(int n1, int n2, int id1){
    int id2 = get_arc(n1, n2)->id;
    if (id2 < id1)
        return true;
    else return false;
}


Net* Net::clone(){
    Net *copy_net = new Net();
    Node* node = NULL;
    for (int i=0; i<nodes.size(); i++) {
        node = this->nodes[i];
        copy_net->add_node(node->clone());
    }
    Arc* arc = NULL;
    for (int i=0; i<arcs.size(); i++){
        /* ignores if the arc is a paralel line to an already existing line */
        if (duplicate(arcs[i]->src->ID, arcs[i]->dest->ID, arcs[i]->id)) {
            continue;
        }
        arc = arcs[i]->clone();
        /* Update the source and destination to the new nodes in copy_net */
        arc->src = copy_net->get_node(arc->src->_name);
        arc->dest = copy_net->get_node(arc->dest->_name);
        /* Add the new arc to the list of arcs IS THIS REALLY USEFULL ? */
        copy_net->add_arc(arc);
        /* Connects it to its source and destination */
        arc->connect();
    }
    return copy_net;
}


const bool node_compare(const Node* n1, const Node* n2) {
    return n1->fill_in > n2->fill_in;
}

void Net::add_node(Node* node){
    node->ID = (int)nodes.size();
    
  //  nodeID.insert(pair<string,Node*>(_name,node)).

    
    if (nodeID.empty())
        {cout << "empty nodeID" << endl;}
    else
        {cout << nodeID.size()<< endl;}

    
    if(!nodeID.insert(pair<string,Node*>(node->_name, node)).second){
        cerr << "ERROR: adding the same node twice!";
    }
    nodes.push_back(node);
}

Node* Net::get_node(string id){
    return nodeID.find(id)->second;
}


/* returns the arc formed by node ids n1 and n2 */
Arc* Net::get_arc(Node* n1, Node* n2){
    string src, dest, key, inv_key;
    src = n1->_name;
    dest = n2->_name;
    key.clear();
    inv_key.clear();
    key.append(src);
    inv_key.append(dest);
    key.append(",");
    inv_key.append(",");
    key.append(dest);
    inv_key.append(src);
    map<string, set<Arc*>*>::iterator it= lineID.find(key);
    if (it != lineID.end()) {
        for (auto a: *it->second){
         //   if (!a->parallel) {
                return a;
           // }
        }
    }
    it = lineID.find(inv_key);
    if (it != lineID.end()) {
        for (auto a: *it->second){
         //   if (!a->parallel) {
                return a;
           // }
        }
    }
    return nullptr;
}


/* returns the Id of the arc formed by node ids n1 and n2 */
Arc* Net::get_arc(int n1, int n2){
    string src, dest, key, inv_key;
    src = to_string(n1);
    dest = to_string(n2);
    key.clear();
    inv_key.clear();
    key.append(src);
    inv_key.append(dest);
    key.append(",");
    inv_key.append(",");
    key.append(dest);
    inv_key.append(src);
    map<string, set<Arc*>*>::iterator it= lineID.find(key);
    if (it != lineID.end()) {
        for (auto a: *it->second){
          //  if (!a->parallel) {
                return a;
            //}
        }
    }
    it = lineID.find(inv_key);
    if (it != lineID.end()) {
        for (auto a: *it->second){
            //if (!a->parallel) {
                return a;
            //}
        }
    }
    return nullptr;
}



bool Net::add_arc(Arc* a){
    bool parallel = false;
    set<Arc*>* s = NULL;
    string src, dest, key, inv_key;
    src = a->src->_name;
    dest = a->dest->_name;
    key.clear();
    inv_key.clear();
    key.append(src);
    inv_key.append(dest);
    key.append(",");
    inv_key.append(",");
    key.append(dest);
    inv_key.append(src);
    if(lineID.find(key)==lineID.end() && lineID.find(inv_key)==lineID.end()){
        s = new set<Arc*>;
        s->insert(a);
        lineID.insert(pair<string, set<Arc*>*>(key,s));
    }
    else {
        if(lineID.find(key)!=lineID.end())
            s = lineID[key];
        if(lineID.find(inv_key)!=lineID.end())
            s = lineID[inv_key];
        s->insert(a);
        cout << "\nWARNING: adding another line between same nodes!\n";
      //  a->parallel = true;
        parallel = true;
    }
    arcs.push_back(a);
    return parallel;
}


void Net::remove_arc(Arc* a){
    //arcs.erase(arcs.at(a->id));
    arcs[a->id] = nullptr;
    lineID.erase(a->src->_name+","+a->dest->_name);
}

// Reading files
char* Net::readline(FILE *input)
{
    size_t len;
    int max_line_len=1024;
    char* line = new char[max_line_len];
    if(std::fgets(line,max_line_len,input)==NULL) return NULL;
    
    while(strrchr(line,'\n') == NULL)
    {
        max_line_len *= 2;
        line = (char *)realloc(line,max_line_len);
        len = strlen(line);
        if(fgets(line+len,max_line_len-len,input) == NULL)
            break;
    }
    return line;
}

void Net::exit_input_error(int line_num){
    fprintf(stderr,"Wrong input format at line %d\n", line_num);
    exit(1);
}

// readFile: just read a matrix, nothing new!
void Net::readFile(string fn){
    auto fname = fn.c_str();
    FILE *fp = fopen(fname,"r");
    if(fp == NULL)
    {
        fprintf(stderr,"can’t open input file %s\n",fname);
        exit(1);
    }
    
    size_t max_line_len = 1024;
    char* line = new char[max_line_len];
    
    vector<vector<int>> matrix;
    int temp;
    while((line = readline(fp))!=NULL)
    {
       vector<int> row;
       stringstream linestream(line);
       while (linestream>>temp)
           row.push_back(temp);
        matrix.push_back(row);
       // cout <<matrix.size()<< endl;
    }
    rewind(fp);
    delete[] line;
    fclose(fp);
}

// populates the graph
void Net::topology(string fn){
    auto fname = fn.c_str();
    FILE *fp = fopen(fname,"r");
    if(fp == NULL)
    {
        fprintf(stderr,"can’t open input file %s\n",fname);
        exit(1);
    }
    
    size_t max_line_len = 1024;
    char* line = new char[max_line_len];
    
    vector<vector<int>> matrix;
    int temp;
    while((line = readline(fp))!=NULL)
    {
        vector<int> row;
        stringstream linestream(line);
        while (linestream>>temp)
            row.push_back(temp);
        matrix.push_back(row);
        
        cout <<matrix.size()<< endl;
    }
    rewind(fp);
    delete[] line;
    fclose(fp);
}

/*  @brief Remove node and all incident arcs from the network
 @note Does not remove the incident arcs from the list of arcs in the network!
 @return the id of the node removed
 */
string Net::remove_end_node(){
    Node* n = nodes.back();
    Node * nn = nullptr;
    string n_id = n->_name;
    for (auto a: n->branches) {
        nn = a->neighbour(n);
        nn->removeArc(a);
        for (auto aa: nn->branches){
            if (!aa->neighbour(nn)->is_connected(n)) {
                nn->fill_in--;
                assert(nn->fill_in >=0);
            }
        }
    }
    nodes.pop_back();
    return n_id;
}


void Net::get_tree_decomp_bags(){
    _bags = new vector<vector<Node*>*>();
    vector<Node*>* bag = nullptr;
    Node* n = nullptr;
    Node* u = nullptr;
    Node* nn = nullptr;
    Arc* arc = nullptr;
    int id = 0;
    int nb3 = 0;
    while (_clone->nodes.size()>2) {
        sort(_clone->nodes.begin(), _clone->nodes.end(),node_compare);
        n = _clone->nodes.back();
        cout << n->_name << endl;
        cout<<_clone->nodes.size() << endl;
        bag = new vector<Node*>();
              //  cout << "new bag = { ";
        for (auto a: n->branches) {
            nn = a->neighbour(n);
            bag->push_back(nn);
            cout << nn->_name << ", ";
        }
        _clone->remove_end_node();
        for (int i = 0; i<bag->size(); i++) {
            u = bag->at(i);
            for (int j = i+1; j<bag->size(); j++) {
                nn = bag->at(j);
                if (u->is_connected(nn)) {
                    continue;
                }
                id = (int)_clone->arcs.size() + 1;
                //arc = new Arc(u->_name+nn->_name);
                arc = new Arc(to_string(id));
                arc->id = id;
                arc->src = u;
                arc->dest = nn;
                arc->connect();
                _clone->add_arc(arc);
            }
        }
        bag->push_back(n);
        //        cout << n->_name << "}\n";
        _bags->push_back(bag);
        if (bag->size()==3) {
            nb3++;
        }
    }
       cout << "\nNumber of 3D bags = " << nb3 << endl;
}


/* Destructors */
Net::~Net(){
    if (!nodes.empty()) {
        for (vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); it++){
            delete (*it);
        }
        nodes.clear();
    }
    if(!arcs.empty()){
        for (vector<Arc*>::iterator it = arcs.begin(); it != arcs.end(); it++){
            if(*it)
                delete (*it);
        }
        arcs.clear();
    }
   
    if(!cycle_basis.empty()){
        for (vector<Path*>::iterator it = cycle_basis.begin(); it != cycle_basis.end(); it++){
            delete (*it);
        }
        arcs.clear();
    }
    for (pair<string,set<Arc*>*> it:lineID) {
        delete it.second;
    }
}

int Net::test(){

    string name;
    int id = 0;
    /*
    cout << "Loading file " << fname << endl;
    ifstream file(fname.c_str());
    if(!file.is_open()){
        cout << "Could not open file\n";
        return -1;
    }
     */
    _clone = new Net();
    //string word;
  
    //file.ignore(6);
    //file >> word;
    //    getline(file, word);
    //_name = word;
    //    cout << _name << endl;
    //while (word.compare("mpc.baseMVA")){
    //    file >> word;
    //}
    //file.ignore(3);
    //getline(file, word,';');
   // bMVA = atoi(word.c_str());
    //    cout << "BaseMVA = " << bMVA << endl;
    
    /* Nodes data */
    /*
    while (word.compare("mpc.bus")){
        file >> word;
    }
    getline(file, word);
     */
    Node* node = NULL;
    Node* node_clone = NULL;
    //file >> word;
  //  while(word.compare("];")){
     //   name = word.c_str();
      //  id = atoi(name.c_str());
    /*
        file >> ws >> word >> ws >> word;
        pl = atof(word.c_str())/bMVA;
        file >> word;
        ql = atof(word.c_str())/bMVA;
        file >> word;
        gs = atof(word.c_str())/bMVA;
        file >> word;
        bs = atof(word.c_str())/bMVA;
        file >> ws >> word >> ws >> word;
        vs = atof(word.c_str());
        file >> ws >> word >> ws >> word;
        kvb = atof(word.c_str());
        file >> ws >> word >> ws >> word;
        vmax = atof(word.c_str());
        getline(file, word,';');
        vmin = atof(word.c_str());
    */
       // node = new Node(name, pl, ql, gs, bs, vmin, vmax, kvb, 1);
        //node_clone = new Node(name, pl, ql, gs, bs, vmin, vmax, kvb, 1);
    for (int i= 0; i<3; i++){
        name = to_string(i);
        id = i;
        cout << "name " << name << " ID: " << i << endl;
        node = new Node(name,id);
        node_clone = new Node(name,i);
        add_node(node);
        _clone->add_node(node_clone);
    }
      //  node->vs = vs;
       // add_node(node);
       // _clone->add_node(node_clone);
        //        node->print();
     //   file >> word;
   // }
   // file.seekg (0, file.beg);
    /* Generator data */
    /*
    while (word.compare("mpc.gen")){
        file >> word;
    }
    double qmin = 0, qmax = 0, pmin = 0, pmax = 0, ps = 0, qs = 0;
    int status = 0;
    getline(file, word);
    file >> word;
    std::vector<bool> gen_status;
    while(word.compare("];")){
        name = word.c_str();
        node = get_node(name);
        file >> word;
        ps = atof(word.c_str())/bMVA;
        file >> word;
        qs = atof(word.c_str())/bMVA;
        file >> word;
        qmax = atof(word.c_str())/bMVA;
        file >> word;
        qmin = atof(word.c_str())/bMVA;
        file >> ws >> word >> ws >> word >> ws >> word;
        status = atof(word.c_str());
        file >> word;
        pmax = atof(word.c_str())/bMVA;
        file >> word;
        pmin = atof(word.c_str())/bMVA;
        getline(file, word,'\n');
        gen_status.push_back(status==1);
        if(status==1){
            node->_has_gen = true;
            Gen* g = new Gen(node, to_string(node->_gen.size()), pmin, pmax, qmin, qmax);
            g->ps = ps;
            g->qs = qs;
            gens.push_back(g);
            node->_gen.push_back(g);
            //            g->print();
        }
        //        getline(file, word);
        file >> word;
    }
    file.seekg (0, file.beg);
     */
    /* Generator costs */
    /*
    while (word.compare("mpc.gencost")){
        file >> word;
    }
    double c0 = 0, c1 = 0,c2 = 0;
    getline(file, word);
    //  cout<<"Number of generators = " << gens.size() << endl;
    int gen_counter = 0;
    for (int i = 0; i < gen_status.size(); ++i) {
        file >> ws >> word >> ws >> word >> ws >> word >> ws >> word >> ws >> word;
        c2 = atof(word.c_str());
        file >> word;
        c1 = atof(word.c_str());
        file >> word;
        c0 = atof(word.c_str());
        if (gen_status[i]) {
            gens[gen_counter]->set_costs(c0, c1, c2);
            //        gens[gen_counter]->print();
            gen_counter++;
        }
        getline(file, word);
    }
    file.seekg (0, file.beg);
    /* Lines data */
    //m_theta_lb = 0;
    //m_theta_ub = 0;
    
    /*
    while (word.compare("mpc.branch")){
        file >> word;
    }
     
    getline(file, word);
    double res = 0;
     */
    Arc* arc = NULL;
    Arc* arc_clone = NULL;
    string src, dest;
   // file >> word;
     
   // while(word.compare("];")){
   //     src = word;
   //     file >> dest;
     
        id = (int)arcs.size() + 1;
        //arc = new Arc(src+dest);
        //arc_clone = new Arc(src+dest);
        arc = new Arc(to_string(id));
        arc_clone = new Arc(to_string(id));
        arc->id = id;
        arc_clone->id = id;
        arc->src = get_node(src);
        arc->dest= get_node(dest);
        arc_clone->src = _clone->get_node(src);
        arc_clone->dest = _clone->get_node(dest);
       // file >> word;
       // arc->r = atof(word.c_str());
       // file >> word;
       // arc->x = atof(word.c_str());
       // res = pow(arc->r,2) + pow(arc->x,2);
       // if (res==0) {
       //     cerr << " line with r = x = 0" << endl;
       //     exit(-1);
       // }
       // arc->g = arc->r/res;
       // arc->b = -arc->x/res;
        
        //file >> word;
        //arc->ch = atof(word.c_str());
        //file >> word;
        //arc->limit = atof(word.c_str())/bMVA;
        //file >> ws >> word >> ws >> word >> ws >> word;
        //if(atof(word.c_str()) == 0)
        //    arc->tr = 1.0;
        //else
        //    arc->tr = atof(word.c_str());
        //file >> word;
        //arc->as = atof(word.c_str())*M_PI/180;
        //file >> word;
        //arc->cc = arc->tr*cos(arc->as);
        //arc->dd = arc->tr*sin(arc->as);
        //arc->status = atof(word.c_str());
        //arc_clone->status = arc->status;
    /*
        file >> word;
        arc->tbound.min = atof(word.c_str())*M_PI/180;
        arc_clone->tbound.min = arc->tbound.min;
        m_theta_lb += arc->tbound.min;
        file >> word;
        arc->tbound.max = atof(word.c_str())*M_PI/180;
        arc_clone->tbound.max = arc->tbound.max;
        m_theta_ub += arc->tbound.max;
        arc->smax = max(pow(arc->src->vbound.max,2)*(arc->g*arc->g+arc->b*arc->b)*(pow(arc->src->vbound.max,2) + pow(arc->dest->vbound.max,2)), pow(arc->dest->vbound.max,2)*(arc->g*arc->g+arc->b*arc->b)*(pow(arc->dest->vbound.max,2) + pow(arc->src->vbound.max,2)));
     */
       // if(arc->status==1){
            arc->connect();
            if(!add_arc(arc)){// not a parallel line
                arc_clone->connect();
                _clone->add_arc(arc_clone);
            }
            else {
                delete arc_clone;
            }
       // }
        //else {
         //   delete arc_clone;
         //   delete arc;
        //}
        //        arc->print();
       // getline(file, word,'\n');
        //file >> word;
    //}
     
    //file.close();
    //    for (auto n:nodes) {
    //        n->print();
    //        cout << "node" << n->ID << ": fill_in = " << n->fill_in << endl;
    //    }
    get_tree_decomp_bags();
    return 0;
}
