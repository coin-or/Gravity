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
    
    param<double> x_min,x_max;//in inputs_pools
    param<double> y_min, y_max;//in pools_outputs
    param<double> z_min, z_max;//in inputs_outputs
    param<double> A_L, A_U; //in inputs
    param<double> D_L, D_U; //in outputs
    param<double> S; //in pools
    param<double> c_tx,c_tz, c_ty;
    param<double> C; //in inputs_attr
    param<double> P_U; //in outputs_attr
    param<double> sumyk;
    param<double> Wij_min, Wij_max, Wii_min, Wii_max, q_min, q_max;

    
    /** Set of all diesel generators */
     indices Inputs;
    
     indices Pools;
    
     indices Outputs;
    
    indices Attr;                         
    
    indices inputs_pools;
    
    indices pools_outputs;
    
    
    indices inputs_outputs;
    
    
    indices inputs_pools_outputs;
    
    indices inputs_attr=indices("inputs_attr");
    indices outputs_attr=indices("outputs_attr");
    indices pools_attr=indices("pools_attr");

    


    /** Constructors */
    PoolNet();
    ~PoolNet();
    
    /** Power grid data parser from Matpower*/
    void readgrid(string fname);
    
        void readgrid1();
    
    

    
    /** get set indexed by bus pairs in the chordal extension */
   // gravity::indices get_bus_pairs_chord(const vector<pair<string,vector<Node*>>>& bags);
    
    
    
    /** Power Model<>s */
   
    indices out_arcs_to_pool_per_input() const;

    indices out_arcs_to_output_per_input() const;

    indices in_arcs_attr_per_pool() const;
    
    indices in_arcs_per_pool_attr() const;

    indices in_arcs_per_pool() const;
    
    indices inputs_pools_outputs_fill() const;

    indices out_arcs_per_pool() const;
    
    indices out_arcs_per_pool_attr() const;

    indices in_arcs_from_pool_per_output_attr() const;
    
     indices in_arcs_from_pool_per_output() const;

    indices in_arcs_from_input_per_output() const;
    
     indices in_arcs_from_input_per_output_attr() const;
    
    indices out_arcs_per_node() const;
    indices in_arcs_per_node() const;
  
    indices out_arcs_per_node(indices arcs) const;
    indices in_arcs_per_node(indices arcs) const;
    
    indices pool_x_matrix_fill() const;
    indices input_x_matrix_fill() const;
    indices output_x_matrix_fill() const;
    indices outattrz_pin_matrix_fill() const;
    indices outattrz_pout_matrix_fill() const;
    vector<indices> outattr_x_p_matrix_fill() const;
    vector<indices> outattrz_p_matrix_fill() const;
    indices pool_q_matrix_fill() const;
    indices pooloutput_x_matrix_fill() const;
    vector<indices> inputpool_x_q_S_matrix_fill() const;
    vector<indices> inpoolout_y_q_matrix_fill() const;
    vector<indices> inpoolout_x_c_matrix_fill() const;
    vector<indices> TxplusTy_fill() const;
 
  
    
    
    
    void fill_wbnds();
};
shared_ptr<Model<>> build_pool_qform(PoolNet& poolnet);
shared_ptr<Model<>> build_pool_pform(PoolNet& poolnet,  SolverType solv_type);
shared_ptr<Model<>> build_pool_pqform(PoolNet& poolnet,  SolverType solv_type);
vector<param<>> fill_wijbounds(PoolNet& poolnet, vector<indices>& vec_pairs);

#endif
