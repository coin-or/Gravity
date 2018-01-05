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
#include "/nh/nest/u/kbestuzheva/usr/include/armadillo"

class Bag{
public:
    int _id;
    PowerNet *_grid;
    std::vector<Node*> _nodes;
    bool _all_lines = true;
    bool _to_add = false;

    vector <gravity::index_*> _indices;
    param<double> _wmin;
    param<double> _wmax;

//    gravity::node_pairs _bus_pairs;

    gravity::Model _model;

    param<double> _wstarp;

    Bag() {};
    Bag(int id, PowerNet* grid, vector<Node*> nodes);
    ~Bag();

    bool is_PSD();

    /* find nearest feasible point */
    param<double> nfp();

    bool add_lines();
};

#endif //GRAVITY_BAG_H
