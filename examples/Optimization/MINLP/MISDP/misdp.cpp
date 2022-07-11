//
//  misdp.cpp
//  misdp
//
//  Created by Smitha on 7/7/22.
//

#include <stdio.h>
#include "read_misdp.h"
#include <gravity/solver.h>

using namespace gravity;
using namespace std;

int main(){
string fname="/Users/smitha/FrameworkMISDP-instances-CBF/TrussTopology/2x4_2scen_3bars.cbf";
auto m=make_shared<Model<double>>("misdp_test");
CBF_read(fname.c_str(), m);
    m->print();
}
