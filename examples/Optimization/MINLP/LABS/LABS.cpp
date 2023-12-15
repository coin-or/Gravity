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
    
    double total_time_start = get_wall_time();
    
    /* binary string dimension*/
    int n = 3;
    bool new_algorithm = false;
    bool nlp = false;
    string alg = "mip";
    if(argc>=2){
        n=stoi(argv[1]);
    }
    if(argc>=3){
        alg = argv[2];
    }

    if(alg=="new")
        new_algorithm =true;
    if(alg=="nlp")
        nlp =true;
    DebugOn("Optimizing for N = " << n << endl);
    
    if(new_algorithm){
        vector<Model<>> models;
        for (int i = 0; i<n; i++) {
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
    var<int> cs("cs", pos_);
    var<> c("c");
    indices s_ids = range(0,n-1);
    indices c_ids = range(1,n-1);
    int opt_obj = 0;
    if(nlp){
        s.exclude_zero();
        s.in(s_ids);
//        M_obj.add(s.in(s_ids));
        M_obj.add(y.in(s_ids));
//        s.set_lb("0",1);
        y.set_lb("0",1);

//        M_obj.add(c.in(c_ids));
        indices pairs("pairs"), quad_terms("quad_terms"), multi_terms("multi_terms"), multi_quad_terms("multi_quad_terms");
        func<> obj;
        for (int k = 1; k<=n-1; k++) {
            func<> cterm;
            for (int i = 0; i<n-k; i++) {
                cterm += (s(to_string(i)))*(s(to_string(i+k)));
            }
            obj += pow(cterm,2);
        }
        obj.print();
        int nb_quad = obj._qterms->size();
        DebugOn("Number of quadratic terms = " << nb_quad << endl);
        int nb_mult = obj._pterms->size();
        DebugOn("Number of multilinear terms = " << nb_mult << endl);

        for(auto p: *obj._qterms)
        {
            auto idx1 = p.second._p->first->_indices->_ids->front().at(0);
            auto idx2 = p.second._p->second->_indices->_ids->front().at(0);
            pairs.add(to_string(idx1)+","+to_string(idx2));
        }
        quad_terms = pairs.subset();
        multi_quad_terms = pairs.subset();
        string pair_idx;
        int nb_row = 0;
        for(auto p: *obj._pterms)
        {
            string mult_idx;
            set<int> tabu;
            for(auto it = p.second._l->begin(); it != p.second._l->end(); it++) {
                bool existing_pair = false;
                auto idx1 = it->first->_indices->_ids->front().at(0);
                mult_idx += to_string(idx1);
                if(tabu.count(idx1)>0)
                    continue;
                for(auto it2 = next(it); it2 != p.second._l->end(); it2++) {
                    auto idx2 = it2->first->_indices->_ids->front().at(0);
                    if(tabu.count(idx2)>0)
                        continue;
                    pair_idx = to_string(idx1)+","+to_string(idx2);
                    if(pairs.has_key(pair_idx)){
                        existing_pair = true;
                        multi_quad_terms.add_in_row(nb_row, pair_idx);
                        tabu.insert(idx1);
                        tabu.insert(idx2);
                        break;
                    }
                }
                if(existing_pair){
                    continue;
                }
                else{
                    pairs.add(pair_idx);
                }
            }
            multi_terms.add(mult_idx);
            nb_row++;
        }
        indices multi_p1("multi_p1"),multi_p2("multi_p2");
        multi_p1 = pairs.subset();
        multi_p2 = pairs.subset();
        for(auto i = 0; i < multi_quad_terms.get_nb_rows(); i++){
            multi_p1._ids->at(0).push_back(multi_quad_terms._ids->at(i).at(0));
            multi_p2._ids->at(0).push_back(multi_quad_terms._ids->at(i).at(1));
        }
        var<int> p("p", -1, 1);
        var<int> pp("pp", -1, 1);
        
        M_obj.add(p.in(pairs));
        M_obj.add(pp.in(multi_terms));
        
        Constraint<> p_def("p_def");
        p_def = p - (2*y.from(pairs) - 1)*(2*y.to(pairs) - 1);
        M_obj.add(p_def.in(pairs) == 0);
        
        Constraint<> pp_def("pp_def");
        pp_def = pp - p.in(multi_p1)*p.in(multi_p2);
        M_obj.add(pp_def.in(multi_terms) == 0);
        
//        var<> obj_var("obj_var", pos_);
//        M_obj.add(obj_var.in(R(1)));
//        
//        Constraint<> obj_def("obj_def");
//        obj_def = obj_var - (4*sum(pp)+ 2*sum(p) + obj.eval_cst(0));
//        M_obj.add(obj_def == 0);
//        
//        M_obj.min(obj_var);
        M_obj.min(4*sum(pp)+ 2*sum(p) + obj.eval_cst(0));
//        M_obj.print();
//        auto g = M_obj.get_interaction_graph();
//        g.print();
//        g.get_tree_decomp_bags();
//        auto ConvM = M_obj.relax();
//        ConvM->print();
//        s.initialize_binary();
//        s.initialize_all(1);
        solver<> g_sol(M_obj,gurobi);
        g_sol.run();
        M_obj.print_solution();
        M_obj.round_solution();
        M_obj.print_solution();
        opt_obj = round(M_obj.get_obj_val());
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
        M.add(c.in(c_ids));
//        M.add(cs.in(c_ids));
//        M.add(z.in(z_ids));
//        M.add(s.in(s_ids));
        M.add(y.in(s_ids));
//        s.set_lb("0",1);
        y.set_lb("0",1);
        
    //    Constraint<> y_on_off("y_on_off");
    //    y_on_off = s - 2*y + 1;
    //    m.add(y_on_off.in(s_ids) == 0);
        
//        Constraint<> z_def("z_def");
//        z_def = z - (2*y.from(z_ids) - 1)*(2*y.to(z_ids) - 1);
//        M.add(z_def.in(z_ids) == 0);
        
        indices z_sum("z_sum");
        indices y_sum_fr("y_sum_fr"), y_sum_to("y_sum_to");
        param<> rhs("rhs");
        rhs.in(c_ids);
        for (int k = 1; k<=n-1; k++) {
            for (int i = 0; i<n-k; i++) {
                z_sum.add_in_row(k-1, to_string(i)+","+to_string(i+k));
                y_sum_fr.add_in_row(k-1, to_string(i));
                y_sum_to.add_in_row(k-1, to_string(i+k));
            }
            rhs.set_val(k-1, n-1-k);
        }
        
        Constraint<> c_def("c_def");
//        c_def = c - z.in(z_sum);
        c_def = c - (2*y.in(y_sum_fr) - 1)*(2*y.in(y_sum_to) - 1) - rhs;
        M.add(c_def.in(c_ids) == 0);
        for(auto i = 0; i < c_ids.size(); i++){
            c.set_lb(c_ids._keys->at(i), -1.*y_sum_fr._ids->at(i).size());
            c.set_ub(c_ids._keys->at(i), y_sum_fr._ids->at(i).size());
        }
//        Constraint<> c_fix("c_fix");
//        c_fix = c[1];
//        M.add(c_fix == 1);
        
//        Constraint<> cs_def("cs_def");
//        cs_def = cs - c*c;
//        M.add(cs_def.in(c_ids) >= 0);
        
//        Constraint<> cs_abs_l("cs_abs_l");
//        cs_abs_l = cs - c;
//        M.add(cs_abs_l.in(c_ids) >= 0);
//        
//        Constraint<> cs_abs_u("cs_abs_u");
//        cs_abs_u = cs + c;
//        M.add(cs_abs_u.in(c_ids) >= 0);
//        
//        M.min(sum(cs));
        param<> ones("1");
        ones.in(c_ids);
        ones = 1;
        auto f = ones.tr()*c*c;
        M.min(ones.tr()*c*c);
//        M.min(sum(cs));
//        M.print();
        solver<> mip_solver(M,gurobi);
//        solver<> nlp_solver(M,ipopt);
        mip_solver.run(1e-6, 7200);
//        f.eval_all();
//        DebugOn("Obj = " << f._val->at(0) << endl);
//        nlp_solver.run();
//        M.print_solution();
//        M.round_solution();
        opt_obj = round(M.get_obj_val());
        M.print_solution();
        
    }
    double total_time_end = get_wall_time();
    auto total_time = total_time_end - total_time_start;
    
    DebugOn("Optimal objective for n = " << n << " is " << opt_obj << endl);
    DebugOn("Wall clock time = " << total_time << endl);
}
