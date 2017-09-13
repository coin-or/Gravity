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
#include <gravity/var.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <assert.h>


// A powerNet is a network.
// Each Node is a Bus.
// Each Arc is a line.
class PowerNet: public Net {
public:
    double bMVA;
    double m_theta_ub;
    double m_theta_lb;

    /** Set of generators */
    std::vector<Gen*> gens;

    /** Constructors */
    PowerNet();
    ~PowerNet();

    /** Power grid data parser */
    int readgrid(const char* fname);
};

//template<typename type = double>
//class var_gen: public var<type> {
//public:
//    var_gen();
//    var_gen(const string& name);
//    ~var_gen(){};
//
//    // member function
//    var<type> in(const vector<Gen*>& gens) {
//        var<type> res(this->get_name());
//        res.set_id(this->get_id());
//        res.set_vec_id(this->get_vec_id());
//        res.set_intype(this->get_intype());
//        res._range = this->_range;
//        //res._val = this->_val;
//        res._lb = this->_lb;
//        res._ub = this->_ub;
//        string key;
//        for(auto it = gens.begin(); it!= gens.end(); it++) {
//            key = (*it)->_name;
//            auto pp = param_::_indices->insert(make_pair<>(key,param_::_indices->size()));
//            if(pp.second) { //new index inserted
//                res._indices->insert(make_pair<>(key,param_::_indices->size()-1));
//                res._ids->push_back(param_::_indices->size()-1);
//            }
//            else {
//                res._indices->insert(make_pair<>(key,pp.first->second));
//                res._ids->push_back(pp.first->second);
//            }
//            res._dim++;
//        }
//        res._name += ".in_gens";
//        res._is_indexed = true;
//        return res;
//    }
//};
#endif
