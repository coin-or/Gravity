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
#include "Gen.h"
#include "Conductor.h"
#include <gravity/Node.h>
#include <gravity/func.h>
#include <stdio.h>

// A bus is a node.
// Additionally, it has some physical features.
class BatteryInverter: public gravity::aux{
public:
    
    
    string _type_name="BattInv";
    int _bat_type = -1;
    /** @brief Max Apparent Power. */
    double _max_s = 0;
    /** @brief Life time. */
    unsigned _lifetime = 0;
    /** @brief Capital Cost. */
    double _capcost = 0;
    /** @brief x coefficients for efficiency curve */
    vector<double> _x_eff;
    /** @brief y coefficients for efficiency curve. */
    vector<double> _y_eff;
    
    BatteryInverter(){};
    BatteryInverter(string name, double max_s, unsigned lifetime, double capcost, vector<double> x, vector<double> y): _max_s(max_s), _lifetime(lifetime), _capcost(capcost){
        _active = true;
        _name = name;
        auto n = x.size();
        _x_eff.resize(n);
        _y_eff.resize(n);
        for (int i = 0; i<n; i++) {
            _x_eff[i] = x[i];
            _y_eff[i] = y[i];
        }
    };
    
    void print() const{
        cout << "name = " << _name << ", max apparent power = " <<_max_s << ", lifetime = " <<_lifetime << ", capital cost = " << _capcost << endl;
    };
};

    
    class WindGen: public gravity::aux{
    public:
        
        
        string _type_name="WindGen";
        
        /** @brief Max Apparent Power. */
        double _max_s = 0;
        /** @brief Life time. */
        unsigned _lifetime = 0;
        /** @brief Capital Cost. */
        double _capcost = 0;
        /** @brief x coefficients for efficiency curve */
        double _x_eff[5];
        /** @brief y coefficients for efficiency curve. */
        double _y_eff[5];
        
        WindGen(){};
        WindGen(string name, double max_s=0, unsigned lifetime=100, double capcost=0): _max_s(max_s), _lifetime(lifetime), _capcost(capcost){
            _active = true;
            _name = name;
        };
    
        void print() const{
            cout << "name = " << _name << ", max apparent power = " <<_max_s << ", lifetime = " <<_lifetime << ", capital cost = " << _capcost << endl;
        };
    };


class PV: public gravity::aux{
public:
    
    
    string _type_name="PV";
    
    /** @brief Min Capacity. */
    double _min_cap = 0;
    
    /** @brief Max Capacity. */
    double _max_cap = 0;
    
    
    /** @brief Life time. */
    unsigned _lifetime = 0;
    
    /** @brief Capital Cost. */
    double _capcost = 0;
    
    /** @brief Variable Cost (Depends on capacity). */
    double _varcost = 0;
    
    /** @brief x coefficients for efficiency curve */
    double _x_eff[5];
    /** @brief y coefficients for efficiency curve. */
    double _y_eff[5];
    
    PV(){};
    PV(string name, double min_cap=0, double max_cap=0, unsigned lifetime=100, double capcost=0, double varcost=0): _min_cap(min_cap), _max_cap(max_cap), _lifetime(lifetime), _capcost(capcost), _varcost(varcost){
        _active = true;
        _name = name;
    };
    
    void print() const{
        cout << "PV name = " << _name << ", min capcacity = " << _min_cap << ", max capcacity = " << _max_cap << ", lifetime = " <<_lifetime << ", capital cost = " << _capcost << ", variable cost = " << _varcost <<endl;
    };
};

class DieselGen: public gravity::aux{
public:
    
    string _type_name="DieselGen";
    
    
    /** @brief Max Real Power. */
    double _max_p = 0;
    /** @brief Max Apparent Power. */
    double _max_s = 0;
    /** @brief Life time. */
    unsigned _lifetime = 0;
    /** @brief Capital Cost. */
    double _capcost = 0;
    /** @brief Costs appearing in the quadratic objective. */
    double _c0 = 0, _c1=0, _c2=0;
    /** @brief Type. */
    unsigned _type = 0;
    /** @brief Efficiency. */
    double _eff = 0;
    /** @brief Max ramping up/down factors. */
    double _max_ramp_down = 0, _max_ramp_up = 0;
    /** @brief Max ramping up/down factors. */
    unsigned _min_down_time = 0, _min_up_time = 0;
    DieselGen(){};
    DieselGen(string name, double max_p, double max_s, unsigned lifetime, double capcost, double c0, double c1, double c2, unsigned type, double eff, double max_ramp_down, double max_ramp_up, double min_down_time, double min_up_time): _max_p(max_p), _max_s(max_s), _lifetime(lifetime), _capcost(capcost), _c0(c0), _c1(c1), _c2(c2), _type(type), _eff(eff), _max_ramp_down(max_ramp_down), _max_ramp_up(max_ramp_up),_min_down_time(min_down_time), _min_up_time(min_up_time){_name = name; _active = true;};
    void print() const{
        cout << "name = " << _name << ", max real power = " << _max_p << ", " << "max apparent power = " <<_max_s << "min down time = " <<_min_down_time << "min up time = " <<_min_up_time << ", " << "lifetime = " <<_lifetime << ", " <<"capital cost = " << _capcost << ", " << "obj costs = (" <<_c0 << ", " << _c1 << ", " << _c2 << "), "<< "type = " <<_type << ", efficiency = " <<_eff << ", max ramp down = " <<_max_ramp_down << ", max ramp up = " <<_max_ramp_up << endl;
    }
};



class DieselData {
public:
    /** @brief Minimal number of potential discrete investment, e.g., diesel generators. */
    unsigned min_invest = 0;
    /** @brief Maximal number of potential discrete investment, e.g., diesel generators. */
    unsigned max_invest = 0;
    /** @brief Number of existing discrete investment options, e.g., diesel generators. */
    unsigned existing = 0;
    /** @brief Age of existing discrete investment options, e.g., diesel generators. */
    unsigned existing_age = 0;
    DieselData(){};
    DieselData(unsigned min, unsigned max, unsigned exist, unsigned exist_age):min_invest(min), max_invest(max), existing(exist), existing_age(exist_age){};
};

class BatteryData {
public:
    /** @brief Minimal number of potential battery investment. */
    unsigned min_invest = 0;
    /** @brief Maximal number of potential battery investment. */
    unsigned max_invest = 0;
    /** @brief Number of existing batteries. */
    unsigned existing = 0;
    /** @brief Age of existing batteries. */
    unsigned existing_age = 0;
    BatteryData(){};
    BatteryData(unsigned min, unsigned max, unsigned exist, unsigned exist_age):min_invest(min), max_invest(max), existing(exist), existing_age(exist_age){};
};

class Line;

class Bus: public Node{
    
public:
    /** @brief Bus base kvolts */
    double _kvb;

    double w;
    
    
    /** @brief Indicates bus type */
    unsigned _type = 0;
    
    /** @brief Indicates if bus has integrated generation */
    bool _has_gen = false;
    
    /** @brief Indicates if bus can have continuous investment, e.g., batteries installed */
    bool _batt_invest = false;
    
    /** @brief Indicates the minimum capacity of batteries that need to be installed on current bus*/
    double _min_batt_cap = 0;
    
    /** @brief Indicates the maximum capacity of batteries that need to be installed on current bus*/
    double _max_batt_cap = 0;
    
    /** @brief Indicates the minimum capacity of PV that need to be installed on current bus*/
    double _min_PV_cap = 0;
    
    /** @brief Indicates the maximum capacity of PV that need to be installed on current bus*/
    double _max_PV_cap = 0;
    
    /** @brief Indicates the maximum area of PV that can be installed on current bus*/
    double _max_PV_area = 0;
    
    /** @brief Indicates the capacity of existing batteries on current bus */
    double _existing_batt_cap = 0;
    
    /** @brief Indicates the capacity of existing PV on current bus */
    double _existing_PV_cap = 0;
    
    /** @brief Indicates the age of existing batteries on current bus */
    unsigned _existing_batt_age = 0;
    
    /** @brief Indicates the age of existing PV on current bus */
    unsigned _existing_PV_age = 0;
    
    /** @brief Indicates if bus can have discrete investment, e.g., diesel generators installed */
    bool _diesel_invest = false;
    
    /** @brief Enable load demand at this node */
    bool _enable_load = true;
    
    /** @brief Enable PV generation at this node */
    bool _enable_PV = false;
    
    
    /** @brief Data on discrete investment per technology type , e.g., diesel generators min/max/existing numbers */
    map<unsigned,DieselData> _diesel_data;
    
    /** @brief Data on battery investment per technology type , min/max/existing numbers */
    map<unsigned,BatteryData> _battery_data;
    
    
    /** @brief Number of wind generators */
    size_t _wind_gens = 0;
    
    /** @brief voltage magnitude bounds */
    Bound vbound;
    
    /** @brief Corresponding set of generators (potential + installed)  */
    std::vector<Gen*> _gen;
    
    /** @brief Corresponding set of generators (potential)  */
    std::vector<Gen*> _pot_gen;
    
    /** @brief Corresponding set of PV generators (potential + installed)  */
    std::vector<PV*> _pv;
    
    /** @brief Corresponding set of PV generators (potential)  */
    std::vector<PV*> _pot_pv;
    
    /** @brief Corresponding set of battery inverters (potential + installed)  */
    std::vector<BatteryInverter*> _bat;
    
    /** @brief Corresponding set of battery inverters (potential)  */
    std::vector<BatteryInverter*> _pot_bat;
    
    /** @brief Corresponding set of wind gens (potential + installed)  */
    std::vector<WindGen*> _wind;
    
    /** @brief Corresponding set of conductors
     @note By default, the network is in one phase, thus exactly one conductor is added to the bus when created. In a three phase network, this array will contain three conductors.
     */
    std::vector<Conductor*> _cond;
    
    /** @brief Set of lines connected to this bus */
    std::map<int, Line*> _lines;
    
    /** @brief Function representing the sum of real power injections at this bus */
    gravity::func_ p_injec;
    
    /** @brief Function representing the sum of reactive power injections at this bus */
    gravity::func_ q_injec;
    
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
    Bus(std::string name, double pl=0, double ql=0, double bs=0, double gs=0, double v_min=0, double v_max=0, double kvb=0,  int phase=0);
    
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
    
    /*
     @brief Returns the vector of connected gens
     */
    vector<gravity::aux*> get_gens();

    /*
     @brief Returns the vector of connected batteries
     */
    vector<gravity::aux*> get_bats();
    
    /*
     @brief Returns the vector of potential gens
     */
    vector<gravity::aux*> get_pot_gens();
    
    /*
     @brief Returns the vector of potential batteries
     */
    vector<gravity::aux*> get_pot_bats();
    
    /*
     @brief Returns the vector of connected wind gens
     */
    vector<gravity::aux*> get_wind();
    
    /*
     @brief Returns the vector of connected PV gens
     */
    vector<gravity::aux*> get_pv();
    
    /** @brief Returns the lower bound on the voltage magnitude at this bus */
    double vmin(void);

    /** @brief Returns the upper bound on the voltage magnitude at this bus */
    double vmax(void);

    /** @brief Connect a generator to the current bus */
    void connect_gen(Gen* g);
    
    /** @brief Create and connect a generator to the current bus */
    void connect_gen(double pmin, double p_max, double q_min, double q_max);
    
    /** @brief Memory release */
    ~Bus();
    void print();
};

#endif
