//
//  Gen.cpp
//  PowerTools++
//
//  Created by Hassan Hijazi on 30/01/2015.
//  Copyright (c) 2015 NICTA. All rights reserved.
//

#include <gravity/Gen.h>

using namespace std;

/** @brief Initialiser with linked Bus and generator properties
 @note Designated initialiser
 */
Gen::Gen(Bus* bus, string name, double p_min, double p_max, double q_min, double q_max):_bus(bus), _name(name), _cost(new GenCost()), _active(true), ps(0), qs(0){
    pbound.min = p_min;
    pbound.max = p_max;
    qbound.min = q_min;
    qbound.max = q_max;
    ID = -1;
};

Gen::~Gen(){
    delete _cost;
}

//void Gen::init_complex(){
//    _S_ = Complex("S", &pg, &qg);
//    _S_._name.append(_name);    
//}

/** @brief Prints generator infos */
void Gen::print(){

    cout << "Gen Id: "<< std::atoi(_name.c_str());
    if(!_active)
        cout << " | WARNING: Inactive generator";
    printf(" | Cost coefficients: (c0=%.02f,c1=%.02f,c2=%.02f)", _cost->c0, _cost->c1, _cost->c2);
    printf(" | Active Power Bounds: (%.02f,%.02f)", pbound.min, pbound.max);
    printf(" | Reactive Power Bounds: (%.02f,%.02f)\n", qbound.min, qbound.max);

};

/** @brief Set cost coefficients fot the current generator */
void Gen::set_costs(double c0, double c1, double c2){
    _cost->c0 = c0;
    _cost->c1 = c1;
    _cost->c2 = c2;
};

