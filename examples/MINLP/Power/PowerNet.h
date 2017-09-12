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
#include "Line.h"
#include <gravity/Net.h>
#include <gravity/model.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <assert.h>


// A powerNet is a network.
// Each Node is a Bus.
// Each Arc is a line.
class PowerNet: public Net{
    
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
#endif
