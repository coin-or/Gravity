//

#include <gravity/Bus.h>
#include <iostream>

using namespace std;

Bus::Bus(){}

Bus::Bus(string name, double pl, double ql, double gs, double bs, double v_min, double v_max, double kvb, int phase): Node(name), _kvb(kvb), _active(true), _has_gen(false), vs(1){
    Conductor* c = new Conductor(this, pl, ql, gs, bs, phase);
    _cond.push_back(c);
    vbound.min = v_min;
    vbound.max = v_max;
};

Bus::~Bus(){
    for (Conductor* cond:_cond) {
        delete cond;
    }
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
    if(_has_gen){
        printf("    List of installed generators:\n");
        for(Gen * g:_gen) {
            g->print();
        }
    }
    printf("    List of connected lines:\n");
    for (auto it:_lines) {
        cout << " " << it.first;
    }
    cout << ";\n";
}
