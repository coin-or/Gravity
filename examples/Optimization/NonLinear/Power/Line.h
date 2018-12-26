////
////  Line.h
////  Cycle_Basis_PF
////
////  Created by Sumiran on 16/06/2014.
////  Modifed by Guanglei on 07/07/2017 for gravity.
////  Copyright (c) 2014 NICTA. All rights reserved.
////
//
//#ifndef Line_h
//#define Line_h
//#include <gravity/Arc.h>
//#include "Bound.h"
//
//// A line is an arc.
//// It has additional physical properties.
//
//class Line: public Arc {
//public:
//    double limit = 0;
//    double ch = 0;
//    double tr = 0;
//    double as = 0;
//    double r = 0;
//    double x = 0;
//    double g = 0;
//    double b = 0;
//    double cc = 0;
//    double dd = 0;
//    double smax = 0;
//    double cost = 0;
//    int status = 1; // on /off, failure or not.
//    unsigned b_type = 0;
//    double wr = 1.0;
//    double wi = 0.0;
//    Bound tbound;
//    //Complex _Si_;
//    //Complex _Sj_;
//    //Complex _W_;
////    var<>   wr;
////    var<>   wi;
////    var<>   pi;
////    var<>   qi;
////    var<>   pj;
////    var<>   qj;
////    var<>   vcs;    /** vi*vj*cs */
////    var<>   vsn;    /** vi*vj*sn */
////    var<>   vv;     /** vi*vj */
////    var<>   cs;     /** cos(thetai - thetaj - as) */
////    var<>   sn;     /** sin(thetai - thetaj - as) */
////    var<>   ci;     /** square of the current magnitude */
////    var<>   w_line_ij;
////    var<>   w_line_ji;
////    var<>   delta;
////    var<bool> on;
//    Line();
//    Line(const string& name);
//    ~Line();
// //   void init_complex();
//    void print();
//};
//
//#endif
