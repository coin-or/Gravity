//
//  Net.h
//
//  Created by Guanglei Wang on 16/06/2014.
//

#ifndef Net_h
#define Net_h

#include <map>
#include <math.h>
#include "Bus.h"
#include "Gen.h"
#include "Line.h"
#include <gravity/Net.h>
#include <gravity/model.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <assert.h>

using namespace gravity;



// A PowerNet extends a standard network by incorporating additional parameters related to power systems.
class PoolNet: public Net {
    
public:
    
    param<double> x_min, x_max, inqual;//in inputs_pools
    param<double> y_min, y_max;//in pools_outputs
    param<double> z_min, z_max,outqual_min, outqual_max;//in inputs_outputs
    param<double> cost, avail_max, avail_min; //in inputs
    param<double> rev, dem_max, dem_min; //in outputs
    param<double> pool_cap; //in pools

    
    /** Set of all diesel generators */
     indices Inputs;
    
     indices Pools;
    
     indices Outputs;
    
    indices Attr;                         
    
    indices inputs_pools;
    
    indices pools_outputs;
    
    indices pools_attr;
    
    indices inputs_outputs;
    
    indices inputs_attr;
    
    indices outputs_attr;

    /** Constructors */
    PoolNet();
    ~PoolNet();
    
    /** Power grid data parser from Matpower*/
    void readgrid();
    
        void readgrid1();
    
     /** Accessors */
    
        void update_net();
    

    
    /** get set indexed by bus pairs in the chordal extension */
   // gravity::indices get_bus_pairs_chord(const vector<pair<string,vector<Node*>>>& bags);
    
    
    
    /** Power Model<>s */
   
    indices out_arcs_to_pool_per_input() const;

    indices out_arcs_to_output_per_input() const;

    indices in_arcs_attr_per_pool() const;
    
    indices in_arcs_per_pool() const;

    indices out_arcs_per_pool() const;

    indices in_arcs_from_pool_per_output() const;

    indices in_arcs_from_input_per_output() const;
    
    indices out_arcs_per_node() const;
    indices in_arcs_per_node() const;
  
    indices out_arcs_per_node(indices arcs) const;
    indices in_arcs_per_node(indices arcs) const;
    
    void fill_wbnds();
};

#endif
