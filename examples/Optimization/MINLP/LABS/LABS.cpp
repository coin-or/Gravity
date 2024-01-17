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
    bool nlp = false, compact = false, ones = false;
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
    if(alg=="compact")
        compact =true;
    if(alg=="ones")
        ones =true;
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
    var<int> z("z", -1, 1);
    var<int> y("y", 0, 1);
    var<int> cs("cs", pos_);
    var<> c("c", -1.*n, n);
    indices s_ids = range(0,n-1);
    indices c_ids = range(1,n-1);
    int opt_obj = 0;
    if(nlp){
        s.exclude_zero();
        s.in(s_ids);
        M_obj.add(s.in(s_ids));
        indices pairs("pairs"), pairs_fr("pairs_fr"), pairs_to("pairs_to"), quad_terms("quad_terms"), multi_terms("multi_terms"), multi_lin_terms("multi_lin_terms"), multi_quad_terms("multi_quad_terms");
        func<> obj;
        for (int k = 1; k<=n-1; k++) {
            func<> cterm;
            for (int i = 0; i<n-k; i++) {
                cterm += (s(to_string(i)))*(s(to_string(i+k)));
            }
            obj += pow(cterm,2);
        }
        
        
        M_obj.min(obj);
        M_obj.print();
        s.initialize_binary();
        solver<> g_sol(M_obj,ipopt);
        g_sol.run();
        M.round_solution();
        opt_obj = round(M_obj.get_obj_val());
        M_obj.print_solution();
    }
    else if(compact){
        s.exclude_zero();
        s.in(s_ids);
//        M_obj.add(s.in(s_ids));
        M_obj.add(y.in(s_ids));
//        s.set_lb("0",1);
        y.set_lb("0",1);

//        M_obj.add(c.in(c_ids));
        indices pairs("pairs"), pairs_fr("pairs_fr"), pairs_to("pairs_to"), quad_terms("quad_terms"), multi_terms("multi_terms"), multi_lin_terms("multi_lin_terms"), multi_quad_terms("multi_quad_terms");
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
        pairs_fr = s_ids.subset();
        pairs_to = s_ids.subset();
        for(auto p: *obj._qterms)
        {
            auto idx1 = p.second._p->first->_indices->_ids->front().at(0);
            auto idx2 = p.second._p->second->_indices->_ids->front().at(0);
            pairs.add(to_string(idx1)+","+to_string(idx2));
            pairs_fr.add_ref(to_string(idx1));
            pairs_to.add_ref(to_string(idx2));
        }
        quad_terms = pairs.subset();
        multi_quad_terms = pairs.subset();
        multi_lin_terms = s_ids.subset();
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
                multi_lin_terms.add_in_row(nb_row, to_string(idx1));
                if(next(it)!= p.second._l->end())
                    mult_idx += ",";
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
//        var<> p("p", -1, 1);
//        var<> pp("pp", -1, 1);
//        
////        M_obj.add(nz.in(multi_terms));
//        M_obj.add(p.in(pairs));
//        M_obj.add(pp.in(multi_terms));
        
//        Constraint<> p_def("p_def");
//        p_def = p - (2*y.from(pairs) - 1)*(2*y.to(pairs) - 1);
//        M_obj.add(p_def.in(pairs) == 0);
        
        var<int> y2("y2", 0, 1);
        M_obj.add(y2.in(pairs));
        
        var<int> yp("yp", 0, 1);
        M_obj.add(yp.in(pairs));
        
        Constraint<> y2_def("y2_def");
        y2_def = y2 - (2*yp.in(pairs)  - y.from(pairs) - y.to(pairs) + 1);
        M_obj.add(y2_def.in(pairs) == 0);
        
        
        Constraint<> yp_def("yp_def");
        yp_def = yp - y.from(pairs)*y.to(pairs);
        M_obj.add(yp_def.in(pairs) == 0);
        
        
//        Constraint<> pp_def("pp_def");
//        pp_def = pp - (2*y2.in(multi_p1) - 1)*(2*y2.in(multi_p2) - 1);
//        M_obj.add(pp_def.in(multi_terms) == 0);
        
        
        
//        var<> obj_var("obj_var", pos_);
//        M_obj.add(obj_var.in(R(1)));
        
        param<> eight("8");
        eight.in(pairs);
        eight = 8;
        
        param<> four_fr("4_fr");
        four_fr.in(pairs_fr);
        four_fr = 4;
        
        param<> four_to("4_to");
        four_to.in(pairs_to);
        four_to = 4;
        
        param<> eight_p1("8_p1");
        eight_p1.in(multi_p1);
        eight_p1 = 8;
        
        param<> eight_p2("8_p2");
        eight_p2.in(multi_p2);
        eight_p2 = 8;
        
        
        param<> sixteen("16");
        sixteen.in(multi_terms);
        sixteen = 16;
        
        M_obj.min((sixteen.tr()*y2.in(multi_p1)*y2.in(multi_p2) - eight_p1.tr()*y2.in(multi_p1) - eight_p2.tr()*y2.in(multi_p2) + 4*nb_mult) + (eight.tr()*yp.in(pairs) - four_fr.tr()*y.in(pairs_fr) - four_to.tr()*y.in(pairs_to) + 2*nb_quad) + obj.eval_cst(0));//4*sum(pp)+
//        M_obj.add(obj_def.in(range(1,1)) == 0);
//        
//        M_obj.min(obj_var);

//        M_obj.print();
//        auto g = M_obj.get_interaction_graph();
//        g.print();
//        g.get_tree_decomp_bags();
//        auto ConvM = M_obj.relax();
//        ConvM->print();
//        s.initialize_binary();
        y.initialize_all(1);
        y2.initialize_all(1);
        yp.initialize_all(1);
        solver<> g_sol(M_obj,ipopt);
        g_sol.run(5,1e-6,10000);
        M_obj.round_solution();
        opt_obj = round(M_obj.get_obj_val());
        M_obj.print_solution();
//        M_obj.print_solution();
        
    }
    else if(ones){
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
        M.add(s.in(s_ids));
//        s.set_lb("0",1);
//        y.set_lb("0",1);
        
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
        
        var<> sp("sp", -1, 1);
        M.add(sp.in(z_ids));
        
        Constraint<> sp_def_ub1("sp_def_ub1");
        sp_def_ub1 = sp - s.from(z_ids);
        M.add(sp_def_ub1.in(z_ids) <= 0);
        
        Constraint<> sp_def_ub2("sp_def_ub2");
        sp_def_ub2 = sp - s.to(z_ids);
        M.add(sp_def_ub2.in(z_ids) <= 0);
        
        Constraint<> sp_def_lb("sp_def_lb");
        sp_def_lb = sp - (s.from(z_ids) + s.to(z_ids) - 1);
        M.add(sp_def_lb.in(z_ids) >= 0);
        
        Constraint<> c_def("c_def");
//        c_def = c - z.in(z_sum);
//        c_def = c - (2*y.in(y_sum_fr) - 1)*(2*y.in(y_sum_to) - 1) - rhs;
        c_def = c - 2*sp.in(z_sum)  - s.in(y_sum_fr) - s.in(y_sum_to);
        M.add(c_def.in(c_ids) == 0);
//        for(auto i = 0; i < c_ids.size(); i++){
//            c.set_lb(c_ids._keys->at(i), -1.*y_sum_fr._ids->at(i).size());
//            c.set_ub(c_ids._keys->at(i), y_sum_fr._ids->at(i).size());
//        }
//
//        Constraint<> c_fix("c_fix");
//        c_fix = c[1];
//        M.add(c_fix == -1);
        
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
        M.print();
        solver<> mip_solver(M,gurobi);
//        solver<> nlp_solver(M,ipopt);
        mip_solver.run(1e-6, 3600);
//        f.eval_all();
//        DebugOn("Obj = " << f._val->at(0) << endl);
//        nlp_solver.run();
//        M.print_solution();
//        M.round_solution();
        opt_obj = round(M.get_obj_val());
        M.print_solution();
        
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
        
        var<> yp("yp", 0, 1);
        M.add(yp.in(z_ids));
        
        Constraint<> yp_def_ub1("yp_def_ub1");
        yp_def_ub1 = yp - y.from(z_ids);
        M.add(yp_def_ub1.in(z_ids) <= 0);
        
        Constraint<> yp_def_ub2("yp_def_ub2");
        yp_def_ub2 = yp - y.to(z_ids);
        M.add(yp_def_ub2.in(z_ids) <= 0);
        
        Constraint<> yp_def_lb("yp_def_lb");
        yp_def_lb = yp - (y.from(z_ids) + y.to(z_ids) - 1);
        M.add(yp_def_lb.in(z_ids) >= 0);
        
        Constraint<> c_def("c_def");
//        c_def = c - z.in(z_sum);
//        c_def = c - (2*y.in(y_sum_fr) - 1)*(2*y.in(y_sum_to) - 1) - rhs;
        c_def = c - 4*yp.in(z_sum)  + 2*y.in(y_sum_fr) + 2*y.in(y_sum_to) - 1 - rhs;
        M.add(c_def.in(c_ids) == 0);
//        for(auto i = 0; i < c_ids.size(); i++){
//            c.set_lb(c_ids._keys->at(i), -1.*y_sum_fr._ids->at(i).size());
//            c.set_ub(c_ids._keys->at(i), y_sum_fr._ids->at(i).size());
//        }
//        
//        Constraint<> c_fix("c_fix");
//        c_fix = c[1];
//        M.add(c_fix == -1);
        
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
        M.print();
        solver<> mip_solver(M,gurobi);
//        solver<> nlp_solver(M,ipopt);
        mip_solver.run(1e-6, 3600);
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
