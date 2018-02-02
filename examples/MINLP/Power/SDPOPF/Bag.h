//
// Created by kbestuzheva on 12/18/17.
//

#ifndef GRAVITY_BAG_H
#define GRAVITY_BAG_H

#include <stdio.h>
#include <gravity/Net.h>
#include <gravity/model.h>
#include <gravity/solver.h>
#include <stdlib.h>
#include "../PowerNet.h"
#include "armadillo"

class Bag{
public:
    int _id;
    bool _is_psd;
    PowerNet* _grid;
    std::vector<Node*> _nodes;
//    bool _all_lines = true;

    vector <gravity::index_> _indices;
    param<double> _wmin;
    param<double> _wmax;

//    gravity::node_pairs _bus_pairs;

    gravity::Model _model;

    param<double> _wstarp;
    param<double> _W_star;

    Bag() {};
    Bag(int id, const PowerNet& grid, vector<Node*> nodes);
    ~Bag();

    bool is_PSD();

    /* find nearest feasible point */
    param<double> nfp();

    bool add_lines();

    param<double> fill_wstar();

    void update_PSD();
};

#endif //GRAVITY_BAG_H
