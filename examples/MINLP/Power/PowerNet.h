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
#include <iostream>
#include <fstream>
#include <vector>
#include <assert.h>

using namespace gravity;

// A PowerNet extends a standard network by incorporating additional parameters related to power systems.
class PowerNet: public Net {

public:
    
    double bMVA; /**< Base MVA */
    double m_theta_lb = 0, m_theta_ub = 0; /**< BigM values for phase angles */
    
    param<double> pg_min, pg_max, qg_min, qg_max, pg_s, qg_s; /**< Upper and lower bounds on generation along woth nominal values (default set points)*/
    param<double> c0, c1, c2; /**< Generation costs */
    param<double> th_min, th_max; /**< Upper and lower bounds on phase angles */
    param<double> v_min, v_max, v_s; /**< Voltage bounds and nominal values (default set points) */
    param<double> pl, ql; /**< Load vectors */
    param<double> gs, bs; /**< Transformer params */
    param<double> g, b, ch, S_max; /**< Power lines parameters, resp., impedance, line charging, and thermal limits */
    param<double> rtr, itr; /**< Transformers phase shifters parameters, a->cc = a->tr*cos(a->as) and a->dd = a->tr*sin(a->as);*/

    
    /** Set of generators */
    std::vector<Gen*> gens;

    /** Constructors */
    PowerNet();
    ~PowerNet();

    /** Power grid data parser */
    int readgrid(const char* fname);
};

#endif
