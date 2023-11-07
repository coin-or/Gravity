//
//  LABS.cpp
//
//  Created by Hassan Hijazi on 11/07/2023.
//

#include <stdio.h>
#include <gravity/solver.h>

/*
    The Low Autocorrelation Binary Sequences (LABS) problem
    see https://www.csplib.org/Problems/prob005/ for more details
 */

using namespace gravity;
using namespace std;

int main(int argc, char * argv[]){
    /* binary string dimension*/
    int n = 3;
    if(argc>=2){
        n=stoi(argv[1]);
    }
    Model<> m("LABS_"+to_string(n));
    var<> s("s", -1, 1);
    var<> z("z", -1, 1);
    var<int> y("y", 0, 1);
    var<> cs("cs", pos_);
    var<> c("c");
    indices s_ids = range(0,n-1);
    indices c_ids = range(1,n-1);
    indices z_ids("z_ids");
    for (int i = 0; i<n-1; i++) {
        for (int k = i+1; k<n; k++) {
            z_ids.add(to_string(i)+","+to_string(k));
        }
    }
    m.add(s.in(s_ids), c.in(c_ids), z.in(z_ids));
    m.add(y.in(s_ids));
    
//    Constraint<> y_on_off("y_on_off");
//    y_on_off = s - 2*y + 1;
//    m.add(y_on_off.in(s_ids) == 0);
    
    Constraint<> z_def("z_def");
    z_def = z - (2*y.from(z_ids) - 1)*(2*y.to(z_ids) - 1);
    m.add(z_def.in(z_ids) == 0);
    
    indices z_sum("z_sum");
    for (int k = 1; k<=n-1; k++) {
        for (int i = 0; i<n-k; i++) {
            z_sum.add_in_row(k-1, to_string(i)+","+to_string(i+k));
        }
    }
    
    Constraint<> c_def("c_def");
    c_def = c - z.in(z_sum);
    m.add(c_def.in(c_ids) == 0);
    
//    Constraint<> cs_def("cs_def");
//    cs_def = cs - c*c;
//    m.add(cs_def.in(c_ids) >= 0);
    
//    m.min(sum(cs));
    param<> ones("1");
    ones.in(c_ids);
    ones = 1;
    m.min(ones.tr()*c*c);
    
//    m.print();
    solver<> g_sol(m,gurobi);
    g_sol.run();
    m.print_solution();
    
}
