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

typedef enum { ACPOL, ACRECT, QC, QC_SDP, OTS, DF, SOCP, SDP, DC, QC_OTS_L, QC_OTS_N, QC_OTS_O, SOCP_OTS, GRB_TEST } PowerModelType;
// A PowerNet extends a standard network by incorporating additional parameters related to power systems.
class PowerNet: public Net {

public:
    string ref_bus;
    double bMVA; /**< Base MVA */
    double bV; /**< Base Voltage */
    double m_theta_lb = 0, m_theta_ub = 0; /**< BigM values for phase angles */
    size_t nb_nodes = 0, nb_branches = 0, nb_gens = 0;
    
    param<Real> pg_min, pg_max, qg_min, qg_max, pg_s, qg_s; /**< Upper and lower bounds on generation along with nominal values (default set points)*/
    param<Real> c0, c1, c2; /**< Generation costs */
    param<Real> th_min, th_max, tan_th_min, tan_th_max, cphi, sphi, cos_d; /**< Upper and lower bounds on phase angles. tan is for the tangent of the angles */
    param<Real> v_min, v_max, v_s; /**< Voltage bounds and nominal values (default set points) */
    param<Real> w_min, w_max; /**< Voltage bounds in lifted W space */
    param<Real> pl, ql; /**< Load vectors */
    param<Real> gs, bs; /**< Transformer params */
    param<Real> tbound_max_tan, tbound_min_tan;  /** tan (th_min), tan(th_max) **/
    param<Real> g, b, ch, tr, as, S_max, wr_min, wr_max, wi_min, wi_max; /**< Power lines parameters, resp., impedance, line charging, and thermal limits. w params are for lifted variavles in W space */
    param<Real> g_ff, g_ft, g_tt, g_tf, b_ff, b_ft, b_tf, b_tt, Y, Y_t, Y_charge, Y_charge_t; /**< Transformers phase shifters parameters, e.g., g_ft = (-a->b*a->cc - a->g*a->dd)/(pow(a->cc,2)+pow(a->dd,2)) where a->cc = a->tr*cos(a->as) and a->dd = a->tr*sin(a->as);*/
    
    /** Set of generators */
    std::vector<Gen*> gens;

    /** Constructors */
    PowerNet();
    ~PowerNet();
    
    /** Power grid data parser */
    int readgrid(const char* fname);
    
    /** Accessors */
    string get_ref_bus();
    unsigned get_nb_active_gens() const;
    unsigned get_nb_active_bus_pairs() const;
    unsigned get_nb_active_arcs() const;
    unsigned get_nb_active_nodes() const;
    void time_expand(unsigned T); /* < Time expansion of the grid parameters */
    void update_net();
};
#endif
