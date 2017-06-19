//
//  MinkMinkmodel.cpp
//  Gravity
//
//  Created by Guanglei Wang on 19/6/17.
//
//

#include "Minkmodel.hpp"
using namespace std;

Minkmodel::Minkmodel(){};
Minkmodel::~Minkmodel(){};

Minkmodel::Minkmodel(ModelType type, Net* graph):_type(type), _graph(graph){};

void Minkmodel::build(){
    switch (_type) {
        case MIP:
            add_vars_origin();
            add_vars_lifted();
            add_obj();
            break;
        case SDP:
            add_vars_lifted();
            add_obj();
            break;
        default:
            break;
    }
}
void Minkmodel::reset(){};

void Minkmodel::add_vars_origin(){
    var<double> zij("z",0,1);
    
    _Minkmodel.add_var(zij^(_graph->nodes.size()*(_graph->nodes.size()-1)/2));

}

