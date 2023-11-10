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
    DebugOn("Optimizing for N = " << n << endl);
    bool new_algorithm = false;
    if(new_algorithm){
        vector<Model<>> models;
        for (int i = 0; i<n-1; i++) {
            Model<> M("LABS_"+to_string(n)+"_"+to_string(i));
            
            models.push_back(M);
        }
        return 0;
    }
    
    
    
    Model<> M_obj("M_obj_LABS_"+to_string(n));
    Model<> M("LABS_"+to_string(n));
    var<int> s("s", -1, 1);
    var<> z("z", -1, 1);
    var<int> y("y", 0, 1);
    var<> cs("cs", pos_);
    var<> c("c");
    indices s_ids = range(0,n-1);
    indices c_ids = range(1,n-1);
    bool unconstrained = true;
    if(unconstrained){
        s.exclude_zero();
        M_obj.add(s.in(s_ids));
        M_obj.add(c.in(c_ids));
        func<int> obj;
        for (int k = 1; k<=n-1; k++) {
            func<int> cterm;
            for (int i = 0; i<n-k; i++) {
                cterm += (s(to_string(i)))*(s(to_string(i+k)));
            }
            cterm.eval_all();
            Constraint<> C_sq_k("C_sq_"+to_string(k));
            C_sq_k = c(k) - pow(cterm,2);
            M_obj.add(C_sq_k>=0);
//            obj += pow(cterm,2);
        }
//        obj.print();
        
        M_obj.min(sum(c));
        M_obj.print();
        s.initialize_binary();
        solver<> g_sol(M_obj,ipopt);
        g_sol.run();
        M_obj.print_solution();
        M_obj.round_solution();
        M_obj.print_solution();
    }
    else{
        indices z_ids("z_ids");
        string idx;
        for (int i = 0; i<n-1; i++) {
            for (int k = i+1; k<n; k++) {
                idx = to_string(i)+","+to_string(k);
                z_ids.add(idx);
            }
        }
        M.add(c.in(c_ids), z.in(z_ids));
        M.add(s.in(s_ids), y.in(s_ids));
//        s.set_lb("0",1);
//        y.set_lb("0",1);
        
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
