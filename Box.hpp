//
//  Box.hpp
//  Gravity
//
//  Created by Guanglei Wang on 23/2/18.
//
//

#ifndef Box_hpp
#define Box_hpp

#include <stdio.h>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;
class box{
private:
    vector<double> lb_, ub_;
    unsigned dim_;
    vector<vector<double> > V_;
    
public:
    box();
    ~box();
    box(double l, double u, unsigned dim):dim_(dim){
        V_.resize(dim);
        lb_.resize(dim);
        ub_.resize(dim);
        for (auto i=0; i < dim; i++){
            lb_[i] = l;
            ub_[i] = u;
        }
    };
    box(const vector<double> l, const vector<double> u, unsigned dim):dim_(dim){
        V_.resize(dim);
        lb_.resize(dim);
        ub_.resize(dim);
        for (auto i=0; i < dim; i++){
            lb_[i] = l[i];
            ub_[i] = u[i];
        }
    };
    unsigned get_dim() const;
    vector<double> get_lb() const;
    vector<double> get_ub() const;
    void set_dim(unsigned);
    void set_lb(const vector<double>);
    void set_ub(const vector<double>);
    void set_vertices(const double, const double, const unsigned);
    void set_vertices(const vector<double>, const vector<double>, const unsigned);

};
#endif /* Box_hpp */
