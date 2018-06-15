//
//  Partition.hpp
//  Gravity
//
//  Created by Guanglei Wang on 12/2/18.
//
//

#ifndef Partition_hpp
#define Partition_hpp

#include <stdio.h>
#include <cstring>
#include <fstream>
#include "../PowerNet.h"
#include <gravity/solver.h>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include <armadillo>
class Partition{
public: 
    Net G_part;
    vector<vector<Node*>> bag_bus;
    std::vector<std::vector<Node*>> bag_bus_union_out;
    vector<vector<Node*>> bag_bus_out;
    vector<vector<Gen*>> bag_gens;
    vector<vector<Line*>> bag_arcs_disjoint;
    vector<vector<Line*>> bag_arcs_neighbour;
    vector<vector<Line*>> bag_arcs_in;
    vector<vector<Line*>> bag_arcs_out;
    vector<vector<Line*>> bag_arcs_union;
    vector<vector<Line*>> bag_arcs_union_in;
    vector<vector<Line*>> bag_arcs_union_out;
    vector<vector<gravity::index_pair*>> bag_bus_pairs_disjoint; // bus_pairs in each bag.
    vector<vector<gravity::index_pair*>> bag_bus_pairs_neighbour;
    vector<vector<gravity::index_pair*>> bag_bus_pairs_neighbour_directed;
    vector<vector<gravity::index_pair*>> bag_bus_pairs_union;
    vector<vector<gravity::index_pair*>> bag_bus_pairs_union_directed;
    vector<gravity::index_pair*> inter_pairs;
    Partition();
    Partition(int Num_parts);
    ~Partition();
    void get_ncut(const PowerNet&, const unsigned&);
};

#endif /* Partition_hpp */
