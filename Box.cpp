//
//  Box.cpp
//  Gravity
//
//  Created by Guanglei Wang on 23/2/18.
//
//

#include "Box.hpp"
#include <assert.h>
box::box() {};
box::~box() {};

using namespace std;
vector<double> box::get_ub() const{
    return ub_;
};

vector<double> box::get_lb() const{
    return lb_;
};
unsigned box::get_dim() const{
    return dim_;
};
void box::set_ub(const vector<double> ub) {
     ub_ = ub;
};

void box::set_lb(const vector<double> lb) {
    lb_ = lb;
};
void  box::set_dim(const unsigned dim) {
     dim_ =dim;
};

void box::set_vertices(const double l, const double u, const unsigned dim){
    if (dim == 1){
        V_.at(0).push_back(l);
        V_.at(1).push_back(u);
    }
    else if (dim < 1)
        std::cerr << "Dim should be as least 1!!" << std::endl;
    else{
        set_vertices(l, u, dim -1);
        unsigned n = pow(2, dim-1);
        for (unsigned i = 0; i < n; i++){
            V_.at(n+i) = V_.at(i);
            V_.at(i).push_back(l);
            V_.at(n+i).push_back(u);
        }
    }
};

void box::set_vertices(const vector<double> l, const vector<double> u, const unsigned dim){
    assert(l.size() == dim_);
    assert(u.size() == dim_);
    assert(dim_== dim);
    V_.resize(dim);
    if (dim == 1){
        V_.at(0).push_back(l[0]);
        V_.at(1).push_back(u[0]);
    }
    else if (dim < 1)
        std::cerr << "Dim should be as least 1!!" << std::endl;
    else{
        vector<double> lb_sub(l.begin()+1, l.end());
        vector<double> ub_sub(u.begin()+1, u.end());
        set_vertices(lb_sub, ub_sub, dim-1);
        unsigned n = pow(2, dim-1);
        for (unsigned i = 0; i < n; i++){
            V_.at(n+i) = V_.at(i);
            V_.at(i).push_back(l[i]);
            V_.at(n+i).push_back(u[i]);
        }
    }
};

