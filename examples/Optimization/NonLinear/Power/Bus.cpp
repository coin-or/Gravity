//

#include "Bus.h"
#include "Line.h"
#include <iostream>

using namespace std;

Bus::Bus(){}

Bus::Bus(string name, double pl, double ql, double gs, double bs, double v_min, double v_max, double kvb, int phase): Node(name), _kvb(kvb), _has_gen(false), vs(1){
    Conductor* c = new Conductor(this, pl, ql, gs, bs, phase);
    _cond.push_back(c);
    vbound.min = v_min;
    vbound.max = v_max;
};

Bus::~Bus(){
    for (auto cond:_cond) {
        delete cond;
    }
    for (auto b: _bat) {
        delete b;
    }
    _bat.clear();
    for (auto w: _wind) {
        delete w;
    }
    _wind.clear();
    for (auto g: _pv) {
        delete g;
    }
    _pv.clear();
}

//void Bus::init_complex(bool polar){
//    if (polar) {
//        _V_ = Complex("V", &vr, &vi, &theta, &v);
//    }
//    else {
//        _V_ = Complex("V", &vr, &vi);
//    }
//    _V_._name.append(_name);
//}

//void Bus::init_lifted_complex(){
//    _V_ = Complex("V", &vr, &vi, &w);
//    _V_._name.append(_name);
//    _V_.lift();
//}
/** @brief Returns the active power load at this bus */
double Bus::pl(){
    return _cond[0]->_pl;
};

/** @brief Returns the reactive power load at this bus */
double Bus::ql(){
    return _cond[0]->_ql;
};

/** @brief Returns the real part of the bus shunt */
double Bus::gs(){
    return _cond[0]->_gs;
};

/** @brief Returns the real part of the bus shunt */
double Bus::bs(){
    return _cond[0]->_bs;
};

std::vector<Line*> Bus::get_out(){
    vector<Line*> res;
    for (auto a: branches) {
        Line* la = (Line*) a;
        if (la->status == 0) {
            continue;
        }
        if(la->_src->_id==_id){
            res.push_back(la);
        }
    }
    return res;
}

std::vector<Line*> Bus::get_in(){
    vector<Line*> res;
    for (auto a:branches) {
        Line* la = (Line*) a;
        if (la->status==0) {
            continue;
        }
        if(la->_dest->_id==_id){
            res.push_back(la);
        }
    }
    return res;
}

vector<gravity::aux*> Bus::get_aux(const string& aux_type){
    if(aux_type=="gens"){
        return get_gens();
    }
    if(aux_type=="pot_gens"){
        return get_pot_gens();
    }
    if(aux_type=="bats"){
        return get_bats();
    }
    if(aux_type=="pot_bats"){
        return get_pot_bats();
    }
    if(aux_type=="wind"){
        return get_wind();
    }
    if(aux_type=="pv"){
        return get_pv();
    }
    throw invalid_argument("In function vector<gravity::aux*> Bus::get_aux(const string& aux_type), unsupported aux_type");
}

vector<gravity::aux*> Bus::get_gens(){
    vector<gravity::aux*> res;
    for (auto g:_gen) {
        res.push_back(g);
    }
    return res;
}

vector<gravity::aux*> Bus::get_pot_gens(){
    vector<gravity::aux*> res;
    for (auto g:_pot_gen) {
        res.push_back(g);
    }
    return res;
}

vector<gravity::aux*> Bus::get_bats(){
    vector<gravity::aux*> res;
    for (auto b:_bat) {
        res.push_back(b);
    }
    return res;
}

vector<gravity::aux*> Bus::get_pot_bats(){
    vector<gravity::aux*> res;
    for (auto b:_pot_bat) {
        res.push_back(b);
    }
    return res;
}

vector<gravity::aux*> Bus::get_wind(){
    vector<gravity::aux*> res;
    for (auto w:_wind) {
        res.push_back(w);
    }
    return res;
}


vector<gravity::aux*> Bus::get_pv(){
    vector<gravity::aux*> res;
    for (auto w:_pv) {
        res.push_back(w);
    }
    return res;
}

/** @brief Returns the lower bound on the voltage magnitude at this bus */
double Bus::vmin(){
    return vbound.min;
};

/** @brief Returns the upper bound on the voltage magnitude at this bus */
double Bus::vmax(){
    return vbound.max;
}


void Bus::print(){
    printf("\nBus Id: %s | load = %.02f | shunt = (%.02f,%.02f) | vbounds = (%.02f,%.02f)\n", _name.c_str(), _cond[0]->_pl, _cond[0]->_gs, _cond[0]->_bs, vbound.min, vbound.max);
        //v.print();
        //theta.print();
    
    if (_wind_gens>0) {
        printf("\tHas %zu wind generators.\n", _wind_gens);
    }
//    if (_diesel_invest) {
//        for (auto &data: _diesel_data) {
//            if (data.second.existing>0) {
//                printf("\tHas %u installed diesel technology of type %u: (age=%u)\n", data.second.existing,data.first, data.second.existing_age);
//            }
//        }
//    }
//    if (_diesel_invest) {
//        printf("\tHas Diesel Generation Investment Option:\n");
//        for (auto &data: _diesel_data) {
//            if (data.second.max_invest>0) {
//                printf("\t\tHas potential diesel technology of type %u: (min=%u, max=%u)\n", data.first, data.second.min_invest,data.second.max_invest);
//            }
//        }
//    }
//    for (auto &data: _battery_data) {
//        if (data.second.existing>0) {
//            printf("\tHas %u installed battery technology of type %u: (age=%u)\n", data.second.existing,data.first, data.second.existing_age);
//        }
//    }
//    if (_batt_invest) {
//        printf("\tHas Battery Investment Option.\n");
//        printf("\tMinimum Battery Investment Option = %g.\n", _min_batt_cap);
//        printf("\tMaximum Battery Investment Option = %g.\n", _max_batt_cap);
//        printf("\tBattery Investment Option:\n");
//        for (auto &data: _battery_data) {
//            if (data.second.max_invest>0) {
//                printf("\t\tHas potential battery inverter technology of type %u: (min=%u, max=%u)\n", data.first, data.second.min_invest,data.second.max_invest);
//            }
//        }
//    }
//    if(_existing_batt_cap > 0){
//        printf("\tHas Existing Battery inverter with Capacity %g anf age %u.\n", _existing_batt_cap, _existing_batt_age);
//    }
    
    printf("\tList of connected lines: ");
    for (auto a:branches) {
        cout << "(" << a->_src->_name<< "," << a->_dest->_name << ") ";
    }
    if(!_bat.empty()){
        printf("\tList of existing + potential batteries:\n");
        for(auto b:_bat) {
            b->print();
        }
    }
    if(!_pot_bat.empty()){
        printf("\tList of potential batteries:\n");
        for(auto b:_pot_bat) {
            b->print();
        }
    }
    if(!_gen.empty()){
        printf("\tList of installed + potential diesnel generators:\n");
        for(Gen * g:_gen) {
            g->print();
        }
    }
    if(!_pot_gen.empty()){
        printf("\tList of potential diesel generators:\n");
        for (auto g:_pot_gen) {
            g->print();
        }
    }
    cout << ";\n";
}
