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
    Model<> M_obj("M_obj_LABS_"+to_string(n));
    Model<> M("LABS_"+to_string(n));
    var<> s("s", -1, 1);
    var<> z("z", -1, 1);
    var<int> y("y", 0, 1);
    var<> cs("cs", pos_);
    var<> c("c");
    indices s_ids = range(0,n-1);

    bool unconstrained = false;
    if(unconstrained){
        M_obj.add(y.in(s_ids));
        func<> obj;        
        for (int k = 1; k<=n-1; k++) {
            func<> cterm;
            for (int i = 0; i<n-k; i++) {
                cterm += (2*y(to_string(i))-1)*(2*y(to_string(i+k))-1);
            }
            obj += pow(cterm,2);
        }
        M_obj.min(obj);
        M_obj.print();
        solver<> g_sol(M_obj,ipopt);
        g_sol.run();
        M_obj.print_solution();
        M_obj.round_solution();
        M_obj.print_solution();
    }
    else{
        indices c_ids = range(1,n-1);
        indices z_ids("z_ids");
        string idx;
        for (int i = 0; i<n-1; i++) {
            for (int k = i+1; k<n; k++) {
                idx = to_string(i)+","+to_string(k);
                z_ids.add(idx);
            }
        }
        M.add(s.in(s_ids), c.in(c_ids), z.in(z_ids));
        M.add(y.in(s_ids));
        
        
    //    Constraint<> y_on_off("y_on_off");
    //    y_on_off = s - 2*y + 1;
    //    m.add(y_on_off.in(s_ids) == 0);
        
        Constraint<> z_def("z_def");
        z_def = z - (2*y.from(z_ids) - 1)*(2*y.to(z_ids) - 1);
        M.add(z_def.in(z_ids) == 0);
        
        indices z_sum("z_sum");
        for (int k = 1; k<=n-1; k++) {
            for (int i = 0; i<n-k; i++) {
                z_sum.add_in_row(k-1, to_string(i)+","+to_string(i+k));
            }
        }
        
        Constraint<> c_def("c_def");
        c_def = c - z.in(z_sum);
        M.add(c_def.in(c_ids) == 0);
        
    //    Constraint<> cs_def("cs_def");
    //    cs_def = cs - c*c;
    //    m.add(cs_def.in(c_ids) >= 0);
        
    //    m.min(sum(cs));
        param<> ones("1");
        ones.in(c_ids);
        ones = 1;
        M.min(ones.tr()*c*c);
        
        solver<> mip_solver(M,gurobi);
//        solver<> nlp_solver(M,ipopt);
        mip_solver.run();
//        nlp_solver.run();
//        M.print_solution();
//        M.round_solution();
        M.print_solution();
    }
    
}
