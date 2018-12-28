//    //
//    //  Gen.h
//    //  PowerTools++
//    //
//    //  Created by Hassan Hijazi on 30/01/2015.
//    //  Copyright (c) 2015 NICTA. All rights reserved.
//    //
//
//#ifndef Gen_h
//#define Gen_h
//
//#include <stdio.h>
//#include <gravity/var.h>
//#include "Bound.h"
//#include "GenCost.h"
//#include <iostream>
//#include <stdio.h>
//#include <gravity/Auxiliary.h>
//
//class Bus;
//
//class Gen: public gravity::aux{
//public:
//    /** @brief Unit name
//     @note One bus may have multiple generators
//     */    
//    string _unit_name;
//    
//    string _type_name="Gen";
//
//    
//    
//    int _gen_type = -1;
//    
//    /** Position of a generator in a container. **/
//    int _id;
//    
//    /** @brief Corresponding Bus */
//    Bus* _bus;
//
//    /** @brief Active Power Generation Bounds */
//    Bound _pbound;
//
//    /** @brief Reactive Power Generation Bounds */
//    Bound _qbound;
//    
//    /** @brief Participation Factor */
//    double _part_factor = 0;
//
//    /** @brief generation cost function/coefficients */
//    GenCost* _cost;
//    
//    /** @brief snapshot value of active power generation */
//    double _ps;
//    /** @brief snapshot value of reactive power generation */
//    double _qs;
//    
//    /** @brief Complex power generation variable */
//    //Complex<double> _S_;
//    /** @brief Active power generation variable */
//    //var<> pg;
//    /** @brief Reactive power generation variable */
//    //var<> qg;
//    Bus* _src = nullptr;
//    Bus* _dest = nullptr;
//    
//    /** @brief Initialiser with linked Bus and generator properties
//     */
//    Gen(Bus* b, std::string name, double p_min, double p_max, double q_min, double q_max);
//    
//    Gen(const string& name);
//    
//    /** @brief Destructor
//     */
//    ~Gen();
//    
//    //void init_complex();
//    
//    /** @brief Prints generator infos */
//    void print();
//    
//    /** @brief Set cost coefficients fot the current generator */
//    void set_costs(double c0, double c1, double c2);
//    
//};
//#endif /* defined(__PowerTools____Gen__) */
