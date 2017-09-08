//
//  Node.cpp
//  Cycle_Basis_PF
//
//  Created by Sumiran on 17/06/2014.
//  Copyright (c) 2014 NICTA. All rights reserved.
//

#include <gravity/PowerNode.h>
#include <iostream>

using namespace std;

PowerNode::PowerNode(){}

PowerNode::PowerNode(string name, int id, double pl, double ql, double gs, double bs, double v_min, double v_max, double kvb,  int phase):Bus(name,pl,ql,gs,bs,v_min,v_max,kvb,phase), Node(name,id){}

PowerNode::PowerNode(string name, double pl, double ql, double gs, double bs, double v_min, double v_max, double kvb,  int phase):PowerNode(name, -1, pl, ql, gs, bs, v_min, v_max, kvb, phase){}

PowerNode::~PowerNode(){}
