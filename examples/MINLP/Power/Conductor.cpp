//
//  Conductor.cpp
//  PowerTools++
//
//  Created by Hassan on 3/02/2015.
//  Copyright (c) 2015 NICTA. All rights reserved.
//

#include "Conductor.h"

Conductor::Conductor(Bus* b, double pl, double ql, double gs, double bs, int phase):bus(b), _pl(pl), _ql(ql), _bs(bs), _gs(gs), _phase(phase){
}

Conductor::~Conductor(){};
