//
//  Conductor.h
//  PowerTools++
//
//  Created by Hassan on 3/02/2015.
//  Copyright (c) 2015 NICTA. All rights reserved.
//

#ifndef Conductor_h
#define Conductor_h

#include <stdio.h>

class Bus;

class Conductor {
    
public:
    /** @brief Conductor active power load */
    double _pl;
    
    /** @brief Conductor reactive power load */
    double _ql;
    
    /** @brief Conductor shunt b value */
    double _bs;
    
    /** @brief Conductor shunt g value */
    double _gs;
    
    /** @brief Conductor phase */
    int _phase;
    
    /** @brief Corresponding bus */
    Bus* bus;
    /**
     @brief Initialiser
     @note Designated initialiser
     */
    Conductor();
    /**
     @brief Initialiser with linked bus and conductor properties
     */
    Conductor(Bus* b, double pl, double ql, double gs, double bs, int phase);
    
    /** @brief Memory release */
    ~Conductor();
    
};
#endif /* defined(__PowerTools____Conductor__) */
