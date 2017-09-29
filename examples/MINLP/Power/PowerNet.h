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
    
    string ref_bus;
    double bMVA; /**< Base MVA */
    double m_theta_lb = 0, m_theta_ub = 0; /**< BigM values for phase angles */
    
    param<double> pg_min, pg_max, qg_min, qg_max, pg_s, qg_s; /**< Upper and lower bounds on generation along woth nominal values (default set points)*/
    param<double> c0, c1, c2; /**< Generation costs */
    param<double> th_min, th_max, tan_th_min, tan_th_max; /**< Upper and lower bounds on phase angles. tan is for the tangent of the angles */
    param<double> v_min, v_max, v_s; /**< Voltage bounds and nominal values (default set points) */
    param<double> w_min, w_max; /**< Voltage bounds in lifted W space */
    param<double> pl, ql; /**< Load vectors */
    param<double> gs, bs; /**< Transformer params */
    param<double> tbound_max_tan, tbound_min_tan;  /** tan (th_min), tan(th_max) **/
    param<double> g, b, ch, S_max, wr_min, wr_max, wi_min, wi_max; /**< Power lines parameters, resp., impedance, line charging, and thermal limits. w params are for lifted variavles in W space */
    param<double> g_ff, g_ft, g_tt, g_tf, b_ff, b_ft, b_tf, b_tt; /**< Transformers phase shifters parameters, e.g., g_ft = (-a->b*a->cc - a->g*a->dd)/(pow(a->cc,2)+pow(a->dd,2)) where a->cc = a->tr*cos(a->as) and a->dd = a->tr*sin(a->as);*/

    
    /** Set of generators */
    std::vector<Gen*> gens;

    /** Constructors */
    PowerNet();
    ~PowerNet();
    

    /** Power grid data parser */
    int readgrid(const char* fname);
    
    /** Accessors */
    string get_ref_bus();
};

#endif
