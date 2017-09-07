//
//  PowerNet.cpp
//
//
//  Created by Guagnlei on 03/06/2017.
//

#include <gravity/PowerNet.h>
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
//#define USEDEBUG
#ifdef USEDEBUG
#define Debug(x) cout << x
#else
#define Debug(x)
#endif
#define DebugOn(x) cout << x
#define DebugOff(x)

using namespace std;


static int max_line_len;
static char* line = nullptr;

PowerNet::PowerNet(){
    bMVA=0;
//    horton_net = nullptr;
    _clone = nullptr;
 //   _bags = nullptr;
}

// readFile: just read a matrix, nothing new!
void PowerNet::readFile(string fn) {
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

    stringstream linestream(line);

    while(readline(fp)!=NULL)
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

// read rudy
void PowerNet::readrudy(string fn) {
    auto fname = fn.c_str();
    int Num_nodes=0;
    int Num_edges=0;


    ifstream infile(fname);

    string sLine;

    if (infile.good())
    {
        getline(infile, sLine);
        istringstream iss(sLine);
        iss >> Num_nodes;
        iss >> Num_edges;
    }
    else{
        fprintf(stderr,"can’t open input file %s\n",fname);
        exit(1);
    }
    // get nodes

    string name;
    _clone = new PowerNet();
    _chordalextension = new PowerNet();
    Node* node = nullptr;
    Node* node_clone = nullptr;
    Node* node_chordal = nullptr;

    for (int i= 1; i<Num_nodes+1; i++) {
        name = to_string(i);
        node = new Node(name,i-1);
        node_clone = new Node(name,i-1);
        node_chordal = new Node(name,i-1);
        add_node(node);
        _clone->add_node(node_clone);
        _chordalextension->add_node(node_chordal);
        // cout << "size of chordal extension is: " << _chordalextension->nodes.size() << endl;
    }


    // get arcs
    Arc* arc = NULL;
    Arc* arc_clone = NULL;
    Arc* arc_chordal = NULL;

    // note that src, dest are names of nodes.
    string src, dest;
    double weight;
    while(getline(infile,sLine,'\n'))
    {
        istringstream iss(sLine);
        iss >> src >> dest >> weight;
        //cout << src  << ", " << dest << ", " << weight << endl;

        name = (int)arcs.size()+1; //

        arc = new Arc(name);
        arc_clone = new Arc(name);
        arc_chordal = new Arc(name);

        arc->id = (int)arcs.size();
        arc_clone->id = (int)_clone->arcs.size();
        arc_chordal->id = (int)_chordalextension->arcs.size();

        arc->src = get_node(src);
        arc->dest= get_node(dest);
        arc->weight=weight;
        add_arc(arc);
        arc->connect();


        arc_clone->src = _clone->get_node(src);
        arc_clone->dest = _clone->get_node(dest);
        arc_clone->weight =weight;
        _clone->add_arc(arc_clone);
        arc_clone->connect();

        arc_chordal->src = _chordalextension->get_node(src);
        arc_chordal->dest = _chordalextension->get_node(dest);
        arc_chordal->weight =weight;
        _chordalextension->add_arc(arc_chordal);
        arc_chordal->connect();
    }
    infile.close();
}




/* Destructors */
PowerNet::~PowerNet() {
    if (!nodes.empty()) {
        for (vector<Node*>::iterator it = nodes.begin(); it != nodes.end(); it++) {
            delete (*it);
        }
        nodes.clear();
    }
    if(!arcs.empty()) {
        for (vector<Arc*>::iterator it = arcs.begin(); it != arcs.end(); it++) {
            if(*it)
                delete (*it);
        }
        arcs.clear();
    }

    if(!cycle_basis.empty()) {
        for (vector<Path*>::iterator it = cycle_basis.begin(); it != cycle_basis.end(); it++) {
            delete (*it);
        }
        arcs.clear();
    }
    for (pair<string,set<Arc*>*> it:lineID) {
        delete it.second;
    }
    delete _clone;
}

