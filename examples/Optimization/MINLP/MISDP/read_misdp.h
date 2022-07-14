//
//  read_misdp.hpp
//  misdp
//
//  Created by Smitha on 7/7/22.
//
#ifndef read_misdp_h
#define read_misdp_h
#include <gravity/solver.h>
#include <stdio.h>

using namespace gravity;
using namespace std;

Net CBF_read(const char *file, shared_ptr<Model<double>>& m);

#endif /* read_misdp_hpp */
