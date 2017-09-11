//
//  Bus.h
//  Cycle_Basis_PF
//
//  Modifed by Guanglei Wang on 08/07/2017.
//  Created by Sumiran on 16/06/2014.
//  Copyright (c) 2014 NICTA. All rights reserved.
//

#ifndef   Bus_h
#define   Bus_h
#include <gravity/Complex.h>
#include <gravity/Gen.h>
#include <gravity/Conductor.h>
#include <gravity/Node.h>
#include <stdio.h>

// A bus is a node.
// Additionally, it has some physical features.
class Line;

class Bus: public Node{
    
public:
    /** @brief Bus base kvolts */
    double _kvb;
    
    /** @brief Bus status, in/out of the network */
    bool _active;
    
    /** @brief Indicates if bus has integrated generation */
    bool _has_gen;

    
    /** @brief voltage magnitude bounds */
    Bound vbound;
    
    /** @brief Corresponding set of generators if installed  */
    std::vector<Gen*> _gen;
    
    /** @brief Corresponding set of conductors
     @note By default, the network is in one phase, thus exactly one conductor is added to the bus when created. In a three phase network, this array will contain three conductors.
     */
    std::vector<Conductor*> _cond;
    
    /** @brief Set of lines connected to this bus */
    std::map<int, Line*> _lines;
    
    /** Complex voltage */
   // Complex _V_;
    
    /*
    * Phase angle variable 
    var<> theta;
    
    * Voltage magnitude variable 
    var<>  v;
    
    * Voltage magnitude squared variable 
    var<>  w;
    
    * Real Voltage 
    var<>  vr;
    * Imaginary Voltage 
    var<>  vi;

    * Active Power Load violation 
    var<>  plv;

    * Reactive Power Load violation 
    var<>  qlv;

    * Voltage bound violation 
    var<>  vbv;
   */
    
    /** @brief snapshot value of voltage magnitude */
    double vs;

    
    // constructor
    Bus();
    
    /** @brief Initialiser with Bus id and conductor properties
     @note By default, the network is in one phase, thus exactly one conductor is added to the bus when created. In a three phase network, this array will contain three conductors.
     */
    Bus(std::string name, double pl, double ql, double bs, double gs, double v_min, double v_max, double kvb,  int phase);
    
   // void init_complex(bool polar);
    
    void init_lifted_complex();
    
    /** @brief Connect a line to the current bus */
    void connect_line(Line* l);
    
    /** @brief Returns the number of active lines connected to this bus */
    int get_degree();

    /** @brief Returns the active power load at this bus */
    double pl();

    /** @brief Returns the real part of the bus shunt */
    double gs();

    /** @brief Returns the real part of the bus shunt */
    double bs();

    
    /** @brief Returns the reactive power load at this bus */
    double ql();

    /*
     @brief Returns the vector of outgoing active arcs
     */
    std::vector<Line*> get_out();
    
    /*
     @brief Returns the vector of incoming active arcs
     */
    std::vector<Line*> get_in();

    /** @brief Returns the lower bound on the voltage magnitude at this bus */
    double vmin();

    /** @brief Returns the upper bound on the voltage magnitude at this bus */
    double vmax();

    /** @brief Connect a generator to the current bus */
    void connect_gen(Gen* g);
    
    /** @brief Create and connect a generator to the current bus */
    void connect_gen(double pmin, double p_max, double q_min, double q_max);
    
    /** @brief Memory release */
    ~Bus();
    void print();
};

#endif
