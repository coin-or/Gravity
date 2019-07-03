//
//  Gen.cpp
//  Gravity++
//
//  Created by Hassan Hijazi on 30/01/2015.

//

#include "Gen.h"

using namespace std;

/** @brief Initialiser with linked Bus and generator properties
 @note Designated initialiser
 */
Gen::Gen(Bus* bus, string name, double p_min, double p_max, double q_min, double q_max):_bus(bus), _cost(new GenCost()), _ps(0), _qs(0){
    _active = true;
    _name = name;
    _pbound.min = p_min;
    _pbound.max = p_max;
    _qbound.min = q_min;
    _qbound.max = q_max;
    _id = -1;
};

Gen::Gen(const string& name):_cost(new GenCost()){
    _name = name;
    _active = true;
}

Gen::~Gen(){
    delete _cost;
}

//void Gen::init_complex(){
//    _S_ = Complex("S", &pg, &qg);
//    _S_._name.append(_name);
//}

/** @brief Prints generator infos */
void Gen::print(){

//    cout << "Gen Id: "<< std::atoi(_name.c_str());
    cout << _name;
    if(!_active)
        cout << " | WARNING: Inactive generator";
    printf(" | Cost coefficients: (c0=%.02f,c1=%.02f,c2=%.02f)", _cost->c0, _cost->c1, _cost->c2);
    printf(" | Active Power Bounds: (%.02f,%.02f)", _pbound.min, _pbound.max);
    printf(" | Reactive Power Bounds: (%.02f,%.02f)", _qbound.min, _qbound.max);
    printf(" | Participation Factor: %.02f%%\n", _part_factor);

};

/** @brief Set cost coefficients fot the current generator */
void Gen::set_costs(double c0, double c1, double c2){
    _cost->c0 = c0;
    _cost->c1 = c1;
    _cost->c2 = c2;
};

