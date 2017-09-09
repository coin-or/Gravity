    //
    //  Gen.h
    //  PowerTools++
    //
    //  Created by Hassan Hijazi on 30/01/2015.
    //  Copyright (c) 2015 NICTA. All rights reserved.
    //

#ifndef Gen_h
#define Gen_h

#include <stdio.h>
#include <gravity/var.h>
#include <gravity/Bound.h>
#include <gravity/GenCost.h>
#include <gravity/Complex.h>
#include <iostream>
#include <stdio.h>

class Bus;

class Gen{
public:
    /** @brief Identifier
     @note One bus may have multiple generators
     */
    std::string _name;
    
    /** @brief Gen status, in/out of the network */
    bool _active;
    
    /** @brief Corresponding Bus */
    Bus* _bus;

    /** @brief Active Power Generation Bounds */
    Bound pbound;

    /** @brief Reactive Power Generation Bounds */
    Bound qbound;

    /** @brief generation cost function/coefficients */
    GenCost* _cost;
    
    /** @brief snapshot value of active power generation */
    double ps;
    /** @brief snapshot value of reactive power generation */
    double qs;
    
    /** @brief Complex power generation variable */
    Complex<double> _S_;
    /** @brief Active power generation variable */
    var<> pg;
    /** @brief Reactive power generation variable */
    var<> qg;
    
    
    /** @brief Initialiser with linked Bus and generator properties
     */
    Gen(Bus* b, std::string name, double p_min, double p_max, double q_min, double q_max);
    
    
    /** @brief Destructor
     */
    ~Gen();
    
    void init_complex();
    
    /** @brief Prints generator infos */
    void print();
    
    /** @brief Set cost coefficients fot the current generator */
    void set_costs(double c0, double c1, double c2);
    
};
#endif /* defined(__PowerTools____Gen__) */
