//
//  model.cpp
//  Gravity
//
//  Created by Hijazi, Hassan.
//
//
//
#include <gravity/model.h>
#include <gravity/solver.h>
#include <math.h> //for setting the rounding direction

using namespace std;
namespace gravity {


/** Lift and linearize the nonlinear constraint c, return the linearized form and add linking constraints to the model.
 @param[in] c: constraint to linearize
 @param[in] partition_model: formulation used for partitionning the nonconvex parts of the constraint
 @return the linearized constraint
 @note This function will add constraints linking the lifted variables to the original ones, if a variable's partition is greater than 1, it will also add the disjunctive constraints corresponding to the partitionning of the variables.
 **/
template <typename type>
template<class T, typename enable_if<is_same<T, Cpx>::value>::type*>
Constraint<type> Model<type>::lift(Constraint<type>& c, string model_type){
    if(c.is_constant() || c.is_linear()){
        return c;
    }
    if(c.is_nonlinear() || c.is_polynomial()){
        throw invalid_argument("lift can only be called on quadratic constraints for now");
    }
    /* Lambda models are taken from Padberg's paper as they are described in type II and type III */
    if((model_type != "on/off") && (model_type != "lambda_II") && (model_type != "lambda_III")){
        throw invalid_argument("model_type can only be one of the following: 'on/off', 'lambda_II', 'lambda_III' ");
    }
    Constraint<Cpx> lifted(c._name+"_lifted");
    if (!c.get_cst()->is_zero()) {
        if (c.get_cst()->is_number()) {
            auto f_cst = static_pointer_cast<constant<Cpx>>(c.get_cst());
            lifted.add_cst(*f_cst);
        }
        else if (c.get_cst()->is_param()) {
            auto f_cst = static_pointer_cast<param<Cpx>>(c.get_cst());
            lifted.add_cst(*f_cst);
        }
        else {
            auto f_cst = static_pointer_cast<func<Cpx>>(c.get_cst());
            lifted.add_cst(*f_cst);
        }
        if (lifted._cst->is_function()) {
            lifted.embed(*static_pointer_cast<func<Cpx>>(lifted._cst));
        }
    }
    for (auto &pair:*c._lterms) {
        auto term = pair.second;
        if (term._coef->is_function()) {
            auto coef = *static_pointer_cast<func<Cpx>>(term._coef);
            term._coef = func<Cpx>(coef).copy();
        }
        else if(term._coef->is_param()) {
            auto coef = *static_pointer_cast<param<Cpx>>(term._coef);
            term._coef = param<Cpx>(coef).copy();
        }
        else if(term._coef->is_number()) {
            auto coef = *static_pointer_cast<constant<Cpx>>(term._coef);
            term._coef = constant<Cpx>(coef).copy();//TODO if T2==type no need to cast
        }
        lifted.insert(term);
    }
    bool lift_sign; /* create lift_sign for correct lower/upper bounding of the variables */
    for (auto &pair:*c._qterms) {
        auto term = pair.second;
        lterm lt;
        lt._sign = term._sign;
        if (term._coef->is_function()) {
            auto coef = *static_pointer_cast<func<Cpx>>(term._coef);
            lt._coef = func<Cpx>(coef).copy();
        }
        else if(term._coef->is_param()) {
            auto coef = *static_pointer_cast<param<Cpx>>(term._coef);
            lt._coef = param<Cpx>(coef).copy();
        }
        else if(term._coef->is_number()) {
            auto coef = *static_pointer_cast<constant<Cpx>>(term._coef);
            lt._coef = constant<Cpx>(coef).copy();
            lift_sign = (term._sign ^ coef.is_negative()); //TODO: update prod_sign in other cases of coef type. Don't know how to do!
        }
        
        if (c.func<Cpx>::is_concave()) //reverse the sign if the constraint is concave
        {
            DebugOn("Changing the sign of the lifted variable." << endl);
            lift_sign = !lift_sign;
        }
        else{
            DebugOn("Keeping the sign of the lifted variable same." << endl);
        }
        
        //arrange the variables so that if they have the same base name, use them ordered in name
        auto o1 = *static_pointer_cast<var<Cpx>>(term._p->first);
        auto o2 = *static_pointer_cast<var<Cpx>>(term._p->second);
        if((o1 != o2) && (o1.get_name(true,true) == o2.get_name(true,true)) && (o1._name > o2._name) ){
            o2 = *static_pointer_cast<var<Cpx>>(term._p->first);
            o1 = *static_pointer_cast<var<Cpx>>(term._p->second);
            DebugOff("O1 name "<< o1._name << endl);
            DebugOff("O2 name "<< o2._name << endl);
        }
        
        string name;
        indices ids;
        if(o1==o2){
            name = "Lift("+o1.get_name(true,true)+"^2)";
            ids = *o1._indices;
        }
        else {
            name = "Lift("+o1.get_name(true,true)+o2.get_name(true,true)+")";
            ids = combine(*o1._indices,*o2._indices);
        }
        auto unique_ids = ids.get_unique_keys(); /* In case of an indexed variable, keep the unique keys only */
        auto o1_ids = *o1._indices;
        auto o2_ids = *o2._indices;
        if(unique_ids.size()!=ids.size()){/* If some keys are repeated, remove them from the refs of o1 and o2 */
            auto keep_refs = ids.get_unique_refs();
            o1_ids.filter_refs(keep_refs);
            o2_ids.filter_refs(keep_refs);
        }
        param<Cpx> lb("lb"), ub("ub");
        lb.in(unique_ids);ub.in(unique_ids);
        auto it = _vars_name.find(name);
        auto name1 = o1.get_name(true,true);
        auto name2 = o2.get_name(true,true);
        if(it==_vars_name.end()){
            /* create the lifted variable with proper lower and upper bounds */
            var<Cpx> vlift(name, lb, ub);
            vlift._lift = true;
            add(vlift.in(unique_ids));
            lt._p = make_shared<var<Cpx>>(vlift.in(ids));
        }
        else {
            auto vlift = static_pointer_cast<var<Cpx>>(it->second);
            auto added = vlift->add_bounds(lb,ub);
            lt._p = make_shared<var<Cpx>>(vlift->in(ids));
            if(!added.empty()){
                assert(o1._indices->size()==o2._indices->size());
                if(added.size()!=o1._indices->size()){/* If some keys are repeated, remove them from the refs of o1 and o2 */
                    auto keep_refs = ids.diff_refs(added);
                    o1_ids.filter_refs(keep_refs);
                    o2_ids.filter_refs(keep_refs);
                }
                reindex_vars();
                // If some keys are repeated in individual indices, remove them from the refs of o1 and o2
                auto o1_ids_uq = o1_ids;
                auto o2_ids_uq = o2_ids;
                auto keep_refs1 = o1_ids_uq.get_unique_refs();
                auto keep_refs2 = o2_ids_uq.get_unique_refs();
                o1_ids_uq.filter_refs(keep_refs1);
                o2_ids_uq.filter_refs(keep_refs2);
                reindex_vars();
            }
        }
        lifted.insert(lt);
    }
    for (auto &pair:*c._pterms) {
        auto term = pair.second;
        lterm lt;
        lt._sign = term._sign;
        if (term._coef->is_function()) {
            auto coef = *static_pointer_cast<func<Cpx>>(term._coef);
            lt._coef = func<Cpx>(coef).copy();
        }
        else if(term._coef->is_param()) {
            auto coef = *static_pointer_cast<param<Cpx>>(term._coef);
            lt._coef = param<Cpx>(coef).copy();
        }
        else if(term._coef->is_number()) {
            auto coef = *static_pointer_cast<constant<Cpx>>(term._coef);
            lt._coef = constant<Cpx>(coef).copy();
        }
        func<Cpx> prod = 1;
        string prod_name = "Lift(";
        auto list = pair.second._l;
        for (auto &ppi: *list) {
            auto p = ppi.first;
            auto orig_var = *static_pointer_cast<var<Cpx>>(p);
            if(ppi.second>1){
                prod_name += orig_var.get_name(true,true)+"("+orig_var._indices->get_name()+")^"+to_string(ppi.second);
                //TODO Lift univarite power function
            }
            else{
                prod_name += orig_var.get_name(true,true)+"("+orig_var._indices->get_name()+")";
            }
            prod *= pow(orig_var,ppi.second);
        }
        prod_name += ")";
        
        auto ids = *c._indices;
        param<Cpx> lb("lb"), ub("ub");
        lb.in(ids);ub.in(ids);
        lb.set_val(prod._range->first);
        ub.set_val(prod._range->second);
        var<Cpx> vlift(prod_name, lb, ub);
        auto it = _vars_name.find(prod_name);
        if(it==_vars_name.end()){
            vlift._lift=true;
            add(vlift.in(ids));
            lt._p = make_shared<var<Cpx>>(vlift);
        }
        else {
            vlift = *static_pointer_cast<var<Cpx>>(it->second);
            lt._p = make_shared<var<Cpx>>(vlift);
        }
        lifted.insert(lt);
    }
    lifted._range = c._range;
    lifted._all_convexity = linear_;
    lifted._all_sign = c._all_sign;
    lifted._ftype = lin_;
    lifted._ctype = c._ctype;
    lifted._indices = c._indices;
    lifted._dim[0] = c._dim[0];
    lifted._dim[1] = c._dim[1];
    return lifted;
}

template <typename type>
template<typename T,typename std::enable_if<is_arithmetic<T>::value>::type*>
Constraint<type> Model<type>::lift(Constraint<type>& c, string model_type){
    if(c.is_constant() || c.is_linear()){
        return c;
    }
    if(c.is_nonlinear() || c.is_polynomial()){
        throw invalid_argument("lift can only be called on quadratic constraints for now");
    }
    /* Lambda models are taken from Padberg's paper as they are described in type II and type III */
    if((model_type != "on/off") && (model_type != "lambda_II") && (model_type != "lambda_III")){
        throw invalid_argument("model_type can only be one of the following: 'on/off', 'lambda_II', 'lambda_III' ");
    }
    Constraint<type> lifted(c._name+"_lifted");
    if (!c.get_cst()->is_zero()) {
        if (c.get_cst()->is_number()) {
            auto f_cst = static_pointer_cast<constant<type>>(c.get_cst());
            lifted.add_cst(*f_cst);
        }
        else if (c.get_cst()->is_param()) {
            auto f_cst = static_pointer_cast<param<type>>(c.get_cst());
            lifted.add_cst(*f_cst);
        }
        else {
            auto f_cst = static_pointer_cast<func<type>>(c.get_cst());
            lifted.add_cst(*f_cst);
        }
        if (lifted._cst->is_function()) {
            lifted.embed(*static_pointer_cast<func<type>>(lifted._cst));
        }
    }
    for (auto &pair:*c._lterms) {
        auto term = pair.second;
        if (term._coef->is_function()) {
            auto coef = *static_pointer_cast<func<type>>(term._coef);
            term._coef = func<type>(coef).copy();
        }
        else if(term._coef->is_param()) {
            auto coef = *static_pointer_cast<param<type>>(term._coef);
            term._coef = param<type>(coef).copy();
        }
        else if(term._coef->is_number()) {
            auto coef = *static_pointer_cast<constant<type>>(term._coef);
            term._coef = constant<type>(coef).copy();//TODO if T2==type no need to cast
        }
        lifted.insert(term);
    }
    bool lift_sign; /* create lift_sign for correct lower/upper bounding of the variables */
    for (auto &pair:*c._qterms) {
        auto term = pair.second;
        lterm lt;
        lt._sign = term._sign;
        if (term._coef->is_function()) {
            auto coef = *static_pointer_cast<func<type>>(term._coef);
            lt._coef = func<type>(coef).copy();
        }
        else if(term._coef->is_param()) {
            auto coef = *static_pointer_cast<param<type>>(term._coef);
            lt._coef = param<type>(coef).copy();
        }
        else if(term._coef->is_number()) {
            auto coef = *static_pointer_cast<constant<type>>(term._coef);
            lt._coef = constant<type>(coef).copy();
            lift_sign = (term._sign ^ coef.is_negative()); //TODO: update prod_sign in other cases of coef type. Don't know how to do!
        }
        
        if (c.func<type>::is_concave()) //reverse the sign if the constraint is concave
        {
            DebugOn("Changing the sign of the lifted variable." << endl);
            lift_sign = !lift_sign;
        }
        else{
            DebugOn("Keeping the sign of the lifted variable same." << endl);
        }
        
        //arrange the variables so that if they have the same base name, use them ordered in name
        auto o1 = *static_pointer_cast<var<type>>(term._p->first);
        auto o2 = *static_pointer_cast<var<type>>(term._p->second);
        if((o1 != o2) && (o1.get_name(true,true) == o2.get_name(true,true)) && (o1._name > o2._name) ){
            o2 = *static_pointer_cast<var<type>>(term._p->first);
            o1 = *static_pointer_cast<var<type>>(term._p->second);
            DebugOff("O1 name "<< o1._name << endl);
            DebugOff("O2 name "<< o2._name << endl);
        }
        
        string name;
        indices ids;
        if(o1==o2){
            name = "Lift("+o1.get_name(true,true)+"^2)";
            ids = *o1._indices;
        }
        else {
            name = "Lift("+o1.get_name(true,true)+o2.get_name(true,true)+")";
            ids = combine(*o1._indices,*o2._indices);
        }
        auto unique_ids = ids.get_unique_keys(); /* In case of an indexed variable, keep the unique keys only */
        auto o1_ids = *o1._indices;
        auto o2_ids = *o2._indices;
        if(unique_ids.size()!=ids.size()){/* If some keys are repeated, remove them from the refs of o1 and o2 */
            auto keep_refs = ids.get_unique_refs();
            o1_ids.filter_refs(keep_refs);
            o2_ids.filter_refs(keep_refs);
        }
        
        // collect the number of partitions of each variable
        int num_partns1 = *o1._num_partns;
        int num_partns2 = *o2._num_partns;
        
        param<type> lb("lb"), ub("ub");
        lb.in(unique_ids);ub.in(unique_ids);
        
        //calculate the tightest valid bounds
        if(o1==o2) //if variables are same, calculate the bounds more efficiently
        {
            for (int i=0; i<unique_ids.size(); i++) {
                //calculate all the possibilities and assign the worst case
                size_t id1;
                if(o1_ids._ids == nullptr){
                    id1 = i;
                }
                else id1 = o1_ids._ids->at(0).at(i);
                auto key1 = o1_ids._keys->at(id1);
                
                auto prod_b1 = o1.get_lb(key1)*o1.get_lb(key1);
                auto prod_b2 = o1.get_lb(key1)*o1.get_ub(key1);
                auto prod_b3 = o1.get_ub(key1)*o1.get_ub(key1);
                
                
                
                lb.set_val(key1, std::max(std::min(std::min(prod_b1,prod_b2), prod_b3), (type)0 ));
                ub.set_val(key1, std::max(std::max(prod_b1,prod_b2),prod_b3));
            }
        }
        else //if variables are different, need to check all four combinations
        {
            for (int i=0; i<unique_ids.size(); i++) {
                //calculate all the possibilities and assign the worst case
                size_t id1;
                size_t id2;
                if(o1_ids._ids == nullptr){
                    id1 = i;
                }
                else id1 = o1_ids._ids->at(0).at(i);
                if(o2_ids._ids == nullptr){
                    id2 = i;
                }
                else id2 = o2_ids._ids->at(0).at(i);
                auto key1 = o1_ids._keys->at(id1);
                auto key2 = o2_ids._keys->at(id2);
                
                auto prod_b1 = o1.get_lb(key1)*o2.get_lb(key2);
                auto prod_b2 = o1.get_lb(key1)*o2.get_ub(key2);
                auto prod_b3 = o1.get_ub(key1)*o2.get_lb(key2);
                auto prod_b4 = o1.get_ub(key1)*o2.get_ub(key2);
                
                lb.set_val(key1+","+key2, std::min(std::min(prod_b1,prod_b2),std::min(prod_b3,prod_b4)));
                ub.set_val(key1+","+key2, std::max(std::max(prod_b1,prod_b2),std::max(prod_b3,prod_b4)));
            }
        }
        
        auto it = _vars_name.find(name);
        
        auto name1 = o1.get_name(true,true);
        auto name2 = o2.get_name(true,true);
        
        if(it==_vars_name.end()){
            
            
            
            //create the lifted variable with proper lower and upper bounds
            var<type> vlift(name, lb, ub);
            vlift._lift = true;
            add(vlift.in(unique_ids));
            lt._p = make_shared<var<type>>(vlift.in(ids));
            
            //check the sign of the lift and the correspoinding bounding functions
            if(c.check_soc() && c.is_eq()){
                if(lift_sign){
                    vlift._lift_ub = true;
                    vlift._lift_lb = false;
                }
                else{
                    vlift._lift_ub = false;
                    vlift._lift_lb = true;
                }
            }
            else{
                vlift._lift_ub = true;
                vlift._lift_lb = true;
            }
            
            
            if((num_partns1 > 1) || (num_partns2 > 1)) {
#ifdef PARTITION
                // If some keys are repeated in individual indices, remove them from the refs of o1 and o2
                auto o1_ids_uq = o1_ids;
                auto o2_ids_uq = o2_ids;
                auto keep_refs1 = o1_ids_uq.get_unique_refs();
                auto keep_refs2 = o2_ids_uq.get_unique_refs();
                o1_ids_uq.filter_refs(keep_refs1);
                o2_ids_uq.filter_refs(keep_refs2);
                reindex_vars();
                if (o1 == o2) //if the variables are same add 1d partition
                {
                    
                    DebugOn("<<<<<<<<<< THIS IS NOT SEEN LIFT -> SINGLE <<<<<<<<<<<" << endl);
                    //create the binary variables for the partitions
                    var<int> on(name1+"_binary",0,1);
                    
                    //create the proper indices and add the binary variables to the model
                    indices partns("partns");
                    for (int i = 0; i < num_partns1 ; ++i)
                    {
                        partns.add(name1+ "{" +to_string(i+1) + "}");
                    }
                    auto inst_partition = indices(unique_ids,partns);
                    add(on.in(inst_partition));
                    
                    //collect the number of entries in each of the index set
                    auto nb_entries_v1 = o1_ids.get_nb_entries();
                    auto nb_entries = unique_ids.get_nb_entries();
                    auto total_entries = inst_partition.get_nb_entries();
                    
                    //add the binary assignment constraint
                    Constraint<> onSum(pair.first + "_binarySum");
                    onSum += sum(on.in_matrix(nb_entries,total_entries-nb_entries));
                    add(onSum.in(unique_ids) == 1);
                    
                    //if the model type is selected as on/off, call on_off formulation for activation of individual constraints
                    if(model_type == "on/off"){
                        add_on_off_McCormick_refined(pair.first, vlift.in(unique_ids), o1.in(o1_ids), o2.in(o2_ids), on);
                    }
                    
                    else{ //means it is one of the lambda formulations
                        
                        //difference is this has one more partition index
                        indices partns_lambda("partns_lambda");
                        for (int i = 0; i < num_partns1+1 ; ++i)
                        {
                            partns_lambda.add(name1+ "{" +to_string(i+1) + "}");
                        }
                        auto inst_partition_lambda = indices(unique_ids,partns_lambda);
                        
                        // Convex combination variables
                        var<> lambda(name1+"_lambda",pos_);
                        add(lambda.in(inst_partition_lambda));
                        
                        /** Parameters */
                        // Bounds on variable v1 & v2
                        param<> bounds(name1+"_bounds");
                        bounds.in(inst_partition_lambda);
                        
                        // Function values on the extreme points
                        param<> EP(name1+name2+"_grid_values");
                        EP.in(inst_partition_lambda);
                        
                        size_t nb_ins = vlift.in(unique_ids).get_nb_inst();
                        auto o1_global_lb = o1.get_lb();
                        auto increment = (o1.get_ub() - o1_global_lb)/num_partns1;
                        
                        // fill bounds and function values
                        for (int i=0 ; i<num_partns1+1; ++i) {
                            auto bound_partn = o1_global_lb + increment*i;
                            bound_partn.eval_all();
                            for (size_t inst = 0; inst< nb_ins; inst++){
                                auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"}";
                                bounds.set_val(cur_idx,bound_partn.eval(inst));
                                EP.set_val(cur_idx,(bound_partn.eval(inst)*bound_partn.eval(inst)));
                            }
                        }
                        
                        // Lambda coefficient matrix when linking with partition variables
                        param<> lambda_coef(name1+"_lambda_linking_coefficients");
                        
                        // Partition coefficient matrix when linking with lambda variables
                        param<> on_coef(name1+"_partition_linking_coefficients");
                        
                        // create constraint indices
                        indices const_idx("const_idx");
                        
                        if(model_type == "lambda_II"){
                            
                            //fill constraint indices
                            for (int i=0; i<num_partns1+1; ++i){
                                const_idx.add(to_string(i+1));
                            }
                            
                            // Lambda coefficient matrix when linking with partition variables
                            lambda_coef.in(indices(inst_partition_lambda, const_idx));
                            
                            // Partition coefficient matrix when linking with lambda variables
                            on_coef.in(indices(inst_partition, const_idx));
                            
                            // fill lambda_coef
                            for (size_t inst = 0; inst< nb_ins; inst++){
                                auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                for (int i=0 ; i<num_partns1+1; ++i) {
                                    string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"}," + to_string(i+1);
                                    lambda_coef.set_val(cur_idx,1);
                                }
                            }
                            
                            // fill on_coef
                            for (size_t inst = 0; inst< nb_ins; inst++){
                                auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"}," + to_string(1);
                                on_coef.set_val(cur_idx,1);
                                for (int i=1 ; i<num_partns1; ++i) {
                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i)+"}," + to_string(i+1);
                                    on_coef.set_val(cur_idx,1);
                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"}," + to_string(i+1);
                                    on_coef.set_val(cur_idx,1);
                                }
                                cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1)+"}," + to_string(num_partns1+1);
                                on_coef.set_val(cur_idx,1);
                            }
                        }
                        
                        else /*means model_type == "lambda_III" */{
                            
                            //fill constraint indices
                            for (int i=0; i<(num_partns1-2)*2+2; ++i){
                                const_idx.add(to_string(i+1));
                            }
                            
                            // Lambda coefficient matrix when linking with partition variables
                            lambda_coef.in(indices(inst_partition_lambda, const_idx));
                            
                            // Partition coefficient matrix when linking with lambda variables
                            on_coef.in(indices(inst_partition, const_idx));
                            
                            // fill lambda_coef
                            for (size_t inst = 0; inst< nb_ins; inst++){
                                auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"}," + to_string(1);
                                lambda_coef.set_val(cur_idx,1);
                                for (int i=1 ; i<(num_partns1-2)*2+1; i=i+2) {
                                    for (int j=(i-1)/2 + 2; j<num_partns1+1; ++j) {
                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"}," + to_string(i+1) ;
                                        lambda_coef.set_val(cur_idx,1);
                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"}," + to_string(i+2);
                                        lambda_coef.set_val(cur_idx,-1);
                                    }
                                }
                                cur_idx =  cur_var_idx+","+name1+"{"+to_string(num_partns1+1)+"}," + to_string((num_partns1-2)*2+2);
                                lambda_coef.set_val(cur_idx,1);
                            }
                            
                            // fill on_coef
                            for (size_t inst = 0; inst< nb_ins; inst++){
                                auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"}," + to_string(1);
                                on_coef.set_val(cur_idx,1);
                                
                                for (int i=1; i<num_partns1; ++i) {
                                    cur_idx =  cur_var_idx+","+name1+"{"+to_string(i+1)+"}," + to_string(2);
                                    on_coef.set_val(cur_idx, 1);
                                }
                                
                                for (int i=2 ; i<(num_partns1-2)*2+2; i=i+2) {
                                    for (int j=i/2+1; j<num_partns1; ++j) {
                                        cur_idx =  cur_var_idx +","+name1+"{"+to_string(j+1)+"}," +  to_string(i+1);
                                        on_coef.set_val(cur_idx,-1);
                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"}," + to_string(i+2) ;
                                        on_coef.set_val(cur_idx,1);
                                    }
                                }
                            }
                            
                        }
                        
                        
                        /** Constraints */
                        if (vlift._lift_ub){
                            // Representation of the quadratic term with secant
                            Constraint<> quad_ub(pair.first+"_quad_ub");
                            quad_ub = EP.in_matrix(nb_entries,total_entries-nb_entries)*lambda.in_matrix(nb_entries,total_entries-nb_entries) - vlift.in(unique_ids);
                            add(quad_ub.in(unique_ids) >= 0); /*using it as the upper bound to be valid*/
                        }
                        if (vlift._lift_lb){
                            Constraint<> quad_lb(pair.first+"_quad_lb");
                            quad_lb = pow(o1.from_ith(0,unique_ids),2) - vlift.in(unique_ids);
                            quad_lb._relaxed = true;
                            add(quad_lb.in(unique_ids) <= 0); /*using it as the lower bound to be valid*/
                        }
                        
                        // Representation of o1 with convex combination
                        Constraint<> o1_rep(pair.first+"_o1_rep");
                        o1_rep = bounds.from_ith(0,inst_partition_lambda).in_matrix(nb_entries, 1) * lambda.in_matrix(nb_entries,total_entries-nb_entries) - o1.in(o1_ids);
                        add(o1_rep.in(unique_ids) == 0);
                        
                        // Linking partition variables with lambda for both lambda formulations
                        if(model_type == "lambda_II"){
                            Constraint<> on_link_lambda(pair.first+"_on_link_lambda_II");
                            on_link_lambda = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef.in_matrix(nb_entries,total_entries-nb_entries) - on.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,on_coef.get_matrix_ids(nb_entries,total_entries-nb_entries)) * on_coef.in_matrix(nb_entries,total_entries-nb_entries);
                            add(on_link_lambda.in(indices(unique_ids,const_idx)) <= 0);
                            
                        }
                        else{ //means model_type == "lambda_III"
                            Constraint<> on_link_lambda(pair.first+"_on_link_lambda_III");
                            on_link_lambda = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef.in_matrix(nb_entries,total_entries-nb_entries) -
                            on.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,on_coef.get_matrix_ids(nb_entries,total_entries-nb_entries)) * on_coef.in_matrix(nb_entries,total_entries-nb_entries);
                            add(on_link_lambda.in(indices(unique_ids,const_idx)) <= 0);
                            
                        }
                        
                        
                        // sum over lambda
                        Constraint<> lambdaSum(pair.first+"_lambdaSum");
                        lambdaSum = sum(lambda.in_matrix(nb_entries,total_entries-nb_entries));
                        add(lambdaSum.in(unique_ids) == 1);
                    }
                    
                }
                else{ //else add 2d partition
                    
                    auto binvar_ptr1 = _vars_name.find(name1+"_binary");
                    auto binvar_ptr2 = _vars_name.find(name2+"_binary");
                    
                    if(binvar_ptr1 !=_vars_name.end() && binvar_ptr2 !=_vars_name.end()){ //means v1 has been partitioned before
                        
                        //if the variables are same core name (means they are same symbolic variable with different indexing)
                        if(name1 == name2){
                            DebugOn("<<<<<<<<<< THIS IS NOT SEEN LIFT -> DOUBLE -> SEEN BOTH BINARIES -> SAME CORE VARS <<<<<<<<<<<" << endl);
                            
                            //create the proper index set for partitions
                            indices partns1("partns1");
                            for (int i = 0; i < num_partns1 ; ++i)
                            {
                                partns1.add(name1+ "{" +to_string(i+1) + "}");
                            }
                            
                            indices partns2("partns2");
                            for (int i = 0; i < num_partns2 ; ++i)
                            {
                                partns2.add(name2+ "{" +to_string(i+1) + "}");
                            }
                            
                            //cast the variable pointer for further use
                            auto binvar1 = static_pointer_cast<var<int>>(binvar_ptr1->second);
                            
                            //define the lower and upper bounds
                            param<int> lb1("lb1"), ub1("ub1");
                            lb1.in(union_ids(o1_ids_uq, o2_ids_uq),partns1);
                            ub1.in(union_ids(o1_ids_uq, o2_ids_uq),partns1);
                            lb1.set_val(0), ub1.set_val(1);
                            auto added1 = binvar1->add_bounds(lb1,ub1);
                            reindex_vars();
                            
                            //collect the number of entries in each index set
                            auto nb_entries_v1 = o1_ids.get_nb_entries();
                            auto nb_entries_v2 = o2_ids.get_nb_entries();
                            auto nb_entries = unique_ids.get_nb_entries();
                            
                            //if there are new indices for the previously defined variable add the corresponding constraint for partitions
                            if(!added1.empty()){
                                Constraint<> onSum1(o1._name+"_binarySum");
                                onSum1 = sum(binvar1->in(added1).in_matrix(nb_entries_v1,1));
                                auto vset1 = added1.from_ith(0,nb_entries_v1);
                                vset1.filter_refs(vset1.get_unique_refs());
                                add(onSum1.in(vset1) == 1);
                            }
                            
                            //if the on/off formulation is chosen for activating constraint
                            if(model_type == "on/off"){
                                var<int> on(name1+name2+"_binary",0,1);
                                
                                indices partns("partns");
                                partns = indices(partns1,partns2);
                                auto inst_partition = indices(unique_ids,partns);
                                add(on.in(inst_partition));
                                auto total_entries = inst_partition.get_nb_entries();
                                
                                Constraint<> onLink1(pair.first+"_binaryLink1");
                                onLink1 = binvar1->from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) - on;
                                add(onLink1.in(inst_partition) >= 0);
                                
                                Constraint<> onLink2(pair.first+"_binaryLink2");
                                onLink2 = binvar1->in_ignore_ith(nb_entries_v2,1,inst_partition.ignore_ith(0,nb_entries_v1)) - on;
                                add(onLink2.in(inst_partition) >= 0);
                                
                                Constraint<> onLink3(pair.first+"_binaryLink3");
                                onLink3 = binvar1->from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) + binvar1->in_ignore_ith(nb_entries_v2,1,inst_partition.ignore_ith(0,nb_entries_v1)) - 1 - on;
                                add(onLink3.in(inst_partition) <= 0);
                                
                                Constraint<> onSumComb(pair.first+"_binarySum");
                                onSumComb = sum(on.in_matrix(nb_entries,total_entries-nb_entries));
                                add(onSumComb.in(unique_ids) == 1);
                                
                                add_on_off_McCormick_refined(pair.first, vlift.in(unique_ids), o1.in(o1_ids), o2.in(o2_ids), on);
                            }
                            
                            else{ //means it is one of the lambda formulations
                                
                                //difference is this has one more partition index
                                indices partns1_lambda("partns1_lambda");
                                for (int i = 0; i < num_partns1+1; ++i)
                                {
                                    partns1_lambda.add(name1+ "{" +to_string(i+1) + "}");
                                }
                                
                                indices partns2_lambda("partns2_lambda");
                                for (int i = 0; i < num_partns2+1; ++i)
                                {
                                    partns2_lambda.add(name2+ "{" +to_string(i+1) + "}");
                                }
                                
                                indices partns_lambda("partns_lambda");
                                partns_lambda = indices(partns1_lambda,partns2_lambda);
                                auto inst_partition_lambda = indices(unique_ids,partns_lambda);
                                auto inst_partition_bounds1 = indices(unique_ids,partns1_lambda);
                                auto inst_partition_bounds2 = indices(unique_ids,partns2_lambda);
                                
                                // Convex combination variables
                                var<> lambda(name1+name2+"_lambda",pos_);
                                add(lambda.in(inst_partition_lambda));
                                
                                /** Parameters */
                                // Bounds on variable v1 & v2
                                param<> bounds1(name1+"_bounds1");
                                bounds1.in(inst_partition_bounds1);
                                
                                param<> bounds2(name2+"_bounds2");
                                bounds2.in(inst_partition_bounds2);
                                
                                // Function values on the extreme points
                                param<> EP(name1+name2+"_grid_values");
                                EP.in(inst_partition_lambda);
                                auto total_entries = inst_partition_lambda.get_nb_entries();
                                
                                size_t nb_ins = vlift.in(unique_ids).get_nb_inst();
                                auto o1_global_lb = o1.get_lb();
                                auto increment1 = (o1.get_ub() - o1_global_lb)/num_partns1;
                                
                                auto o2_global_lb = o2.get_lb();
                                auto increment2 = (o2.get_ub() - o2_global_lb)/num_partns2;
                                
                                // fill bounds and function values
                                for (int i=0 ; i<num_partns1+1; ++i) {
                                    auto bound_partn1 = o1_global_lb + increment1*i;
                                    bound_partn1.eval_all();
                                    auto bound_partn2 = o2_global_lb + increment2*i;
                                    bound_partn2.eval_all();
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                        auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"}";
                                        bounds1.set_val(cur_idx,bound_partn1.eval(inst));
                                        bounds2.set_val(cur_idx,bound_partn2.eval(inst));
                                        for(int j=0; j<num_partns2+1; ++j){
                                            auto bound_partn2_temp = o2_global_lb + increment2*j;
                                            bound_partn2_temp.eval_all();
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"}";
                                            EP.set_val(cur_idx,(bound_partn1.eval(inst)*bound_partn2_temp.eval(inst)));
                                        }
                                    }
                                }
                                
                                // Lambda coefficient matrix when linking with partition variables
                                param<> lambda_coef1(name1+"_lambda_linking_coefficients1");
                                param<> lambda_coef2(name2+"_lambda_linking_coefficients2");
                                
                                // Partition coefficient matrix when linking with lambda variables
                                param<> on_coef1(name1+"_partition_linking_coefficients1");
                                param<> on_coef2(name2+"_partition_linking_coefficients2");
                                
                                // create constraint indices
                                indices const_idx1("const_idx1");
                                indices const_idx2("const_idx2");
                                
                                if(model_type == "lambda_II"){
                                    
                                    //fill constraint indices
                                    for (int i=0; i<num_partns1+1; ++i){
                                        const_idx1.add(to_string(i+1));
                                    }
                                    
                                    //fill constraint indices
                                    for (int i=0; i<num_partns2+1; ++i){
                                        const_idx2.add(to_string(i+1));
                                    }
                                    
                                    // Lambda coefficient matrix when linking with partition variables
                                    lambda_coef1.in(indices(inst_partition_lambda, const_idx1));
                                    lambda_coef2.in(indices(inst_partition_lambda, const_idx2));
                                    
                                    // Lambda coefficient matrix when linking with partition variables
                                    on_coef1.in(indices(unique_ids, partns1, const_idx1));
                                    on_coef2.in(indices(unique_ids, partns2, const_idx2));
                                    
                                    // fill lambda_coef1 and lambda_coef2
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                        auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                        for (int i=0 ; i<num_partns1+1; ++i) {
                                            for (int j=0 ; j<num_partns2+1; ++j) {
                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                lambda_coef1.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(j+1);
                                                lambda_coef2.set_val(cur_idx,1);
                                            }
                                        }
                                    }
                                    
                                    // fill on_coef1 and on_coef2
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                        auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                        on_coef1.set_val(cur_idx,1);
                                        cur_idx = cur_var_idx+","+name2+"{"+to_string(1)+"},"+to_string(1);
                                        on_coef2.set_val(cur_idx,1);
                                        for (int i=1 ; i<num_partns1; ++i) {
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(i)+"},"+to_string(i+1);
                                            on_coef1.set_val(cur_idx,1);
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(i+1);
                                            on_coef1.set_val(cur_idx,1);
                                        }
                                        for (int i=1 ; i<num_partns2; ++i) {
                                            cur_idx = cur_var_idx+","+name2+"{"+to_string(i)+"},"+to_string(i+1);
                                            on_coef2.set_val(cur_idx,1);
                                            cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"},"+to_string(i+1);
                                            on_coef2.set_val(cur_idx,1);
                                        }
                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1)+"},"+to_string(num_partns1+1);
                                        on_coef1.set_val(cur_idx,1);
                                        cur_idx = cur_var_idx+","+name2+"{"+to_string(num_partns2)+"},"+to_string(num_partns2+1);
                                        on_coef2.set_val(cur_idx,1);
                                    }
                                }
                                
                                
                                else /*means model_type == "lambda_III" */{
                                    
                                    //fill constraint indices
                                    for (int i=0; i<(num_partns1-2)*2+2; ++i){
                                        const_idx1.add(to_string(i+1));
                                    }
                                    
                                    //fill constraint indices
                                    for (int i=0; i<(num_partns2-2)*2+2; ++i){
                                        const_idx2.add(to_string(i+1));
                                    }
                                    
                                    // Lambda coefficient matrix when linking with partition variables
                                    lambda_coef1.in(indices(inst_partition_lambda, const_idx1));
                                    lambda_coef2.in(indices(inst_partition_lambda, const_idx2));
                                    
                                    // Partition coefficient matrix when linking with lambda variables
                                    on_coef1.in(indices(unique_ids, partns1, const_idx1));
                                    on_coef2.in(indices(unique_ids, partns2, const_idx2));
                                    
                                    // fill lambda_coef1 and lambda_coef2
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                        auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                        for (int j=0; j<num_partns2+1; ++j) {
                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(1);
                                            lambda_coef1.set_val(cur_idx,1);
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string((num_partns1-2)*2+2);
                                            lambda_coef1.set_val(cur_idx,1);
                                        }
                                        
                                        for (int i=1 ; i<(num_partns1-2)*2+1; i=i+2) {
                                            for (int j=(i-1)/2 + 2; j<num_partns1+1; ++j) {
                                                for(int k=0; k<num_partns2+1; ++k){
                                                    string cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+name2+"{"+to_string(k+1)+"},"+to_string(i+1);
                                                    lambda_coef1.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+name2+"{"+to_string(k+1)+"},"+to_string(i+2);
                                                    lambda_coef1.set_val(cur_idx,-1);
                                                }
                                            }
                                        }
                                        
                                        for (int i=0; i<num_partns1+1; ++i) {
                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(1)+"},"+to_string(1);
                                            lambda_coef2.set_val(cur_idx,1);
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(num_partns2+1)+"},"+to_string((num_partns2-2)*2+2);
                                            lambda_coef2.set_val(cur_idx,1);
                                        }
                                        
                                        for (int i=1 ; i<(num_partns2-2)*2+1; i=i+2) {
                                            for (int j=(i-1)/2 + 2; j<num_partns2+1; ++j) {
                                                for(int k=0; k<num_partns1+1; ++k){
                                                    string cur_idx = cur_var_idx+","+name1+"{"+to_string(k+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                    lambda_coef2.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(k+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                    lambda_coef2.set_val(cur_idx,-1);
                                                }
                                            }
                                        }
                                    }
                                    
                                    
                                    
                                    // fill on_coef1 and on_coef2
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                        auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                        on_coef1.set_val(cur_idx,1);
                                        
                                        for (int i=1; i<num_partns1; ++i) {
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(2);
                                            on_coef1.set_val(cur_idx, 1);
                                        }
                                        
                                        for (int i=2 ; i<(num_partns1-2)*2+2; i=i+2) {
                                            for (int j=i/2+1; j<num_partns1; ++j) {
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                on_coef1.set_val(cur_idx,-1);
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                on_coef1.set_val(cur_idx,1);
                                            }
                                        }
                                        
                                        cur_idx = cur_var_idx+","+name2+"{"+to_string(1)+"},"+to_string(1);
                                        on_coef2.set_val(cur_idx,1);
                                        for (int i=1; i<num_partns2; ++i) {
                                            cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"},"+to_string(2);
                                            on_coef2.set_val(cur_idx, 1);
                                        }
                                        
                                        for (int i=2 ; i<(num_partns2-2)*2+2; i=i+2) {
                                            for (int j=i/2+1; j<num_partns2; ++j) {
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                on_coef2.set_val(cur_idx,-1);
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                on_coef2.set_val(cur_idx,1);
                                            }
                                        }
                                    }
                                }
                                
                                
                                /** Constraints */
                                // Representation of the bilinear term with convex combination
                                Constraint<> bln_rep(pair.first+"_bln_rep");
                                bln_rep = EP.in_matrix(nb_entries,total_entries-nb_entries)*lambda.in_matrix(nb_entries,total_entries-nb_entries) - vlift.in(unique_ids);
                                add(bln_rep.in(unique_ids) == 0);
                                
                                // Representation of o1 with convex combination
                                Constraint<> o1_rep(pair.first+"_o1_rep");
                                
                                o1_rep = bounds1.from_ith(0,inst_partition_lambda).in_matrix(nb_entries, 1) * lambda.in_matrix(nb_entries,total_entries-nb_entries) - o1.in(o1_ids);
                                
                                add(o1_rep.in(unique_ids) == 0);
                                
                                // Representation of o2 with convex combination
                                Constraint<> o2_rep(pair.first+"_o2_rep");
                                
                                o2_rep = bounds2.in_ignore_ith(nb_entries, 1, inst_partition_lambda).in_matrix(nb_entries,1) * lambda.in_matrix(nb_entries,total_entries-nb_entries) - o2.in(o2_ids);
                                add(o2_rep.in(unique_ids) == 0);
                                
                                // Linking partition variables1 with lambda
                                if(model_type == "lambda_II"){
                                    Constraint<> on_link_lambda1(pair.first+"_on_link_lambda1_II");
                                    
                                    on_link_lambda1 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef1.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef1.in_matrix(nb_entries,total_entries-nb_entries) - binvar1->in_ignore_ith(nb_entries_v1,nb_entries_v2,on_coef1.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef1.in_matrix(nb_entries,1);
                                    add(on_link_lambda1.in(indices(unique_ids,const_idx1)) <= 0);
                                    
                                    Constraint<> on_link_lambda2(pair.first+"_on_link_lambda2_II");
                                    on_link_lambda2 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef2.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef2.in_matrix(nb_entries,total_entries-nb_entries) - binvar1->in_ignore_ith(0,nb_entries_v1,on_coef2.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef2.in_matrix(nb_entries,1);
                                    add(on_link_lambda2.in(indices(unique_ids,const_idx2)) <= 0);
                                    
                                }
                                else{
                                    Constraint<> on_link_lambda1(pair.first+"_on_link_lambda1_III");
                                    
                                    on_link_lambda1 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef1.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef1.in_matrix(nb_entries,total_entries-nb_entries) - binvar1->in_ignore_ith(nb_entries_v1,nb_entries_v2,on_coef1.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef1.in_matrix(nb_entries,1);
                                    add(on_link_lambda1.in(indices(unique_ids,const_idx1)) <= 0);
                                    
                                    Constraint<> on_link_lambda2(pair.first+"_on_link_lambda2_III");
                                    on_link_lambda2 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef2.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef2.in_matrix(nb_entries,total_entries-nb_entries) - binvar1->in_ignore_ith(0,nb_entries_v1,on_coef2.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef2.in_matrix(nb_entries,1);
                                    add(on_link_lambda2.in(indices(unique_ids,const_idx2)) <= 0);
                                }
                                // sum over lambda
                                Constraint<> lambdaSum(pair.first+"_lambdaSum");
                                lambdaSum = sum(lambda.in_matrix(nb_entries,total_entries-nb_entries));
                                add(lambdaSum.in(unique_ids) == 1);
                            }
                        }
                        
                        else{
                            DebugOn("<<<<<<<<<< THIS IS NOT SEEN LIFT -> DOUBLE -> SEEN BOTH BINARIES -> DIFF CORE VARS <<<<<<<<<<<" << endl);
                            
                            indices partns1("partns1");
                            for (int i = 0; i < num_partns1 ; ++i)
                            {
                                partns1.add(name1+ "{" +to_string(i+1) + "}");
                            }
                            indices partns2("partns2");
                            for (int i = 0; i < num_partns2 ; ++i)
                            {
                                partns2.add(name2+ "{" +to_string(i+1) + "}");
                            }
                            
                            auto binvar1 = static_pointer_cast<var<int>>(binvar_ptr1->second);
                            param<int> lb1("lb1"), ub1("ub1");
                            lb1.in(o1_ids_uq,partns1);
                            ub1.in(o1_ids_uq,partns1);
                            lb1.set_val(0), ub1.set_val(1);
                            auto added1 = binvar1->add_bounds(lb1,ub1);
                            reindex_vars();
                            
                            auto binvar2 = static_pointer_cast<var<int>>(binvar_ptr2->second);
                            param<int> lb2("lb2"), ub2("ub2");
                            lb2.in(o2_ids_uq,partns2);
                            ub2.in(o2_ids_uq,partns2);
                            lb2.set_val(0), ub2.set_val(1);
                            auto added2 = binvar2->add_bounds(lb2,ub2);
                            reindex_vars();
                            
                            auto nb_entries_v1 = o1_ids.get_nb_entries();
                            auto nb_entries_v2 = o2_ids.get_nb_entries();
                            auto nb_entries = unique_ids.get_nb_entries();
                            
                            if(!added1.empty()){
                                Constraint<> onSum1(o1._name+"_binarySum");
                                onSum1 = sum(binvar1->in(added1).in_matrix(nb_entries_v1,1));
                                auto vset1 = added1.from_ith(0,nb_entries_v1);
                                vset1.filter_refs(vset1.get_unique_refs());
                                add(onSum1.in(vset1) == 1);
                            }
                            
                            if(!added2.empty()){
                                Constraint<> onSum2(o2._name+"_binarySum");
                                onSum2 = sum(binvar2->in(added2).in_matrix(nb_entries_v2,1));
                                auto vset2 = added2.from_ith(0,nb_entries_v2);
                                vset2.filter_refs(vset2.get_unique_refs());
                                add(onSum2.in(vset2) == 1);
                            }
                            
                            if(model_type == "on/off"){ //if on/off is chosen
                                
                                var<int> on(name1+name2+"_binary",0,1);
                                
                                indices partns("partns");
                                partns = indices(partns1,partns2);
                                auto inst_partition = indices(unique_ids,partns);
                                add(on.in(inst_partition));
                                auto total_entries = inst_partition.get_nb_entries();
                                
                                Constraint<> onLink1(pair.first+"_binaryLink1");
                                onLink1 = binvar1->from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) - on;
                                add(onLink1.in(inst_partition) >= 0);
                                
                                Constraint<> onLink2(pair.first+"_binaryLink2");
                                onLink2 = binvar2->in_ignore_ith(nb_entries_v2,1,inst_partition.ignore_ith(0,nb_entries_v1)) - on;
                                add(onLink2.in(inst_partition) >= 0);
                                
                                Constraint<> onLink3(pair.first+"_binaryLink3");
                                onLink3 = binvar1->from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) + binvar2->in_ignore_ith(nb_entries_v2,1,inst_partition.ignore_ith(0,nb_entries_v1)) - 1 - on;
                                add(onLink3.in(inst_partition) <= 0);
                                
                                Constraint<> onSumComb(pair.first+"_binarySum");
                                onSumComb = sum(on.in_matrix(nb_entries,total_entries-nb_entries));
                                add(onSumComb.in(unique_ids) == 1);
                                
                                add_on_off_McCormick_refined(pair.first, vlift.in(unique_ids), o1.in(o1_ids), o2.in(o2_ids), on);
                            }
                            
                            else{ //means it is one of the lambda formulations
                                
                                //difference is this has one more partition index
                                indices partns1_lambda("partns1_lambda");
                                for (int i = 0; i < num_partns1+1; ++i)
                                {
                                    partns1_lambda.add(name1+ "{" +to_string(i+1) + "}");
                                }
                                
                                indices partns2_lambda("partns2_lambda");
                                for (int i = 0; i < num_partns2+1; ++i)
                                {
                                    partns2_lambda.add(name2+ "{" +to_string(i+1) + "}");
                                }
                                
                                indices partns_lambda("partns_lambda");
                                partns_lambda = indices(partns1_lambda,partns2_lambda);
                                auto inst_partition_lambda = indices(unique_ids,partns_lambda);
                                auto inst_partition_bounds1 = indices(unique_ids,partns1_lambda);
                                auto inst_partition_bounds2 = indices(unique_ids,partns2_lambda);
                                
                                // Convex combination variables
                                var<> lambda(name1+name2+"_lambda",pos_);
                                add(lambda.in(inst_partition_lambda));
                                
                                /** Parameters */
                                // Bounds on variable v1 & v2
                                param<> bounds1(name1+"_bounds1");
                                bounds1.in(inst_partition_bounds1);
                                
                                param<> bounds2(name2+"_bounds2");
                                bounds2.in(inst_partition_bounds2);
                                
                                // Function values on the extreme points
                                param<> EP(name1+name2+"_grid_values");
                                EP.in(inst_partition_lambda);
                                auto total_entries = inst_partition_lambda.get_nb_entries();
                                
                                size_t nb_ins = vlift.in(unique_ids).get_nb_inst();
                                auto o1_global_lb = o1.get_lb();
                                auto increment1 = (o1.get_ub() - o1_global_lb)/num_partns1;
                                
                                auto o2_global_lb = o2.get_lb();
                                auto increment2 = (o2.get_ub() - o2_global_lb)/num_partns2;
                                
                                // fill bounds1 and function values
                                for (int i=0 ; i<num_partns1+1; ++i) {
                                    auto bound_partn1 = o1_global_lb + increment1*i;
                                    bound_partn1.eval_all();
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                        auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"}";
                                        bounds1.set_val(cur_idx,bound_partn1.eval(inst));
                                        for(int j=0; j<num_partns2+1; ++j){
                                            auto bound_partn2 = o2_global_lb + increment2*j;
                                            bound_partn2.eval_all();
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"}";
                                            EP.set_val(cur_idx,(bound_partn1.eval(inst)*bound_partn2.eval(inst)));
                                        }
                                    }
                                }
                                // fill bounds2
                                for (int i=0 ; i<num_partns2+1; ++i) {
                                    auto bound_partn2 = o2_global_lb + increment2*i;
                                    bound_partn2.eval_all();
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                        auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                        string cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"}";
                                        bounds2.set_val(cur_idx,bound_partn2.eval(inst));
                                    }
                                }
                                
                                // Lambda coefficient matrix when linking with partition variables
                                param<> lambda_coef1(name1+"_lambda_linking_coefficients1");
                                param<> lambda_coef2(name2+"_lambda_linking_coefficients2");
                                
                                // Partition coefficient matrix when linking with lambda variables
                                param<> on_coef1(name1+"_partition_linking_coefficients1");
                                param<> on_coef2(name2+"_partition_linking_coefficients2");
                                
                                // create constraint indices
                                indices const_idx1("const_idx1");
                                indices const_idx2("const_idx2");
                                
                                if(model_type == "lambda_II"){
                                    
                                    //fill constraint indices
                                    for (int i=0; i<num_partns1+1; ++i){
                                        const_idx1.add(to_string(i+1));
                                    }
                                    
                                    //fill constraint indices
                                    for (int i=0; i<num_partns2+1; ++i){
                                        const_idx2.add(to_string(i+1));
                                    }
                                    
                                    // Lambda coefficient matrix when linking with partition variables
                                    if(num_partns1 > 1) lambda_coef1.in(indices(inst_partition_lambda, const_idx1));
                                    if(num_partns2 > 1) lambda_coef2.in(indices(inst_partition_lambda, const_idx2));
                                    
                                    // Lambda coefficient matrix when linking with partition variables
                                    if(num_partns1 > 1) on_coef1.in(indices(unique_ids, partns1, const_idx1));
                                    if(num_partns2 > 1) on_coef2.in(indices(unique_ids, partns2, const_idx2));
                                    
                                    // fill lambda_coef1 and lambda_coef2
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                        auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                        for (int i=0 ; i<num_partns1+1; ++i) {
                                            for (int j=0 ; j<num_partns2+1; ++j) {
                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                if(num_partns1 > 1) lambda_coef1.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(j+1);
                                                if(num_partns2 > 1) lambda_coef2.set_val(cur_idx,1);
                                            }
                                        }
                                    }
                                    
                                    // fill on_coef1 and on_coef2
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                        auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                        if(num_partns1 > 1) on_coef1.set_val(cur_idx,1);
                                        cur_idx = cur_var_idx+","+name2+"{"+to_string(1)+"},"+to_string(1);
                                        if(num_partns2 > 1) on_coef2.set_val(cur_idx,1);
                                        if(num_partns1 > 1){
                                            for (int i=1 ; i<num_partns1; ++i) {
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(i)+"},"+to_string(i+1);
                                                on_coef1.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(i+1);
                                                on_coef1.set_val(cur_idx,1);
                                            }
                                        }
                                        if(num_partns2 > 1){
                                            for (int i=1 ; i<num_partns2; ++i) {
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(i)+"},"+to_string(i+1);
                                                on_coef2.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"},"+to_string(i+1);
                                                on_coef2.set_val(cur_idx,1);
                                            }
                                        }
                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1)+"},"+to_string(num_partns1+1);
                                        if(num_partns1 > 1) on_coef1.set_val(cur_idx,1);
                                        cur_idx = cur_var_idx+","+name2+"{"+to_string(num_partns2)+"},"+to_string(num_partns2+1);
                                        if(num_partns2 > 1) on_coef2.set_val(cur_idx,1);
                                    }
                                }
                                
                                
                                else /*means model_type == "lambda_III" */{
                                    
                                    //fill constraint indices
                                    for (int i=0; i<(num_partns1-2)*2+2; ++i){
                                        const_idx1.add(to_string(i+1));
                                    }
                                    
                                    //fill constraint indices
                                    for (int i=0; i<(num_partns2-2)*2+2; ++i){
                                        const_idx2.add(to_string(i+1));
                                    }
                                    
                                    // Lambda coefficient matrix when linking with partition variables
                                    if(num_partns1 > 1) lambda_coef1.in(indices(inst_partition_lambda, const_idx1));
                                    if(num_partns2 > 1) lambda_coef2.in(indices(inst_partition_lambda, const_idx2));
                                    
                                    // Partition coefficient matrix when linking with lambda variables
                                    if(num_partns1 > 1) on_coef1.in(indices(unique_ids, partns1, const_idx1));
                                    if(num_partns2 > 1) on_coef2.in(indices(unique_ids, partns2, const_idx2));
                                    
                                    // fill lambda_coef1 and lambda_coef2
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                        auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                        if(num_partns1 > 1) {
                                            for (int j=0; j<num_partns2+1; ++j) {
                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(1);
                                                lambda_coef1.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string((num_partns1-2)*2+2);
                                                lambda_coef1.set_val(cur_idx,1);
                                            }
                                            
                                            for (int i=1 ; i<(num_partns1-2)*2+1; i=i+2) {
                                                for (int j=(i-1)/2 + 2; j<num_partns1+1; ++j) {
                                                    for(int k=0; k<num_partns2+1; ++k){
                                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+name2+"{"+to_string(k+1)+"},"+to_string(i+1);
                                                        lambda_coef1.set_val(cur_idx,1);
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+name2+"{"+to_string(k+1)+"},"+to_string(i+2);
                                                        lambda_coef1.set_val(cur_idx,-1);
                                                    }
                                                }
                                            }
                                        }
                                        
                                        if(num_partns2 > 1){
                                            for (int i=0; i<num_partns1+1; ++i) {
                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(1)+"},"+to_string(1);
                                                lambda_coef2.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(num_partns2+1)+"},"+to_string((num_partns2-2)*2+2);
                                                lambda_coef2.set_val(cur_idx,1);
                                            }
                                            
                                            for (int i=1 ; i<(num_partns2-2)*2+1; i=i+2) {
                                                for (int j=(i-1)/2 + 2; j<num_partns2+1; ++j) {
                                                    for(int k=0; k<num_partns1+1; ++k){
                                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(k+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                        lambda_coef2.set_val(cur_idx,1);
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(k+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                        lambda_coef2.set_val(cur_idx,-1);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    
                                    
                                    
                                    // fill on_coef1 and on_coef2
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                        auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                        if(num_partns1 > 1) {
                                            on_coef1.set_val(cur_idx,1);
                                            
                                            
                                            for (int i=1; i<num_partns1; ++i) {
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(2);
                                                on_coef1.set_val(cur_idx, 1);
                                            }
                                            
                                            for (int i=2 ; i<(num_partns1-2)*2+2; i=i+2) {
                                                for (int j=i/2+1; j<num_partns1; ++j) {
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                    on_coef1.set_val(cur_idx,-1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                    on_coef1.set_val(cur_idx,1);
                                                }
                                            }
                                        }
                                        
                                        if(num_partns2 > 1) {
                                            cur_idx = cur_var_idx+","+name2+"{"+to_string(1)+"},"+to_string(1);
                                            on_coef2.set_val(cur_idx,1);
                                            for (int i=1; i<num_partns2; ++i) {
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"},"+to_string(2);
                                                on_coef2.set_val(cur_idx, 1);
                                            }
                                            
                                            for (int i=2 ; i<(num_partns2-2)*2+2; i=i+2) {
                                                for (int j=i/2+1; j<num_partns2; ++j) {
                                                    cur_idx = cur_var_idx+","+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                    on_coef2.set_val(cur_idx,-1);
                                                    cur_idx = cur_var_idx+","+name2+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                    on_coef2.set_val(cur_idx,1);
                                                }
                                            }
                                        }
                                    }
                                }
                                
                                
                                /** Constraints */
                                // Representation of the bilinear term with convex combination
                                Constraint<> bln_rep(pair.first+"_bln_rep");
                                bln_rep = EP.in_matrix(nb_entries,total_entries-nb_entries)*lambda.in_matrix(nb_entries,total_entries-nb_entries) - vlift.in(unique_ids);
                                add(bln_rep.in(unique_ids) == 0);
                                
                                // Representation of o1 with convex combination
                                Constraint<> o1_rep(pair.first+"_o1_rep");
                                o1_rep = bounds1.from_ith(0,inst_partition_lambda).in_matrix(nb_entries, 1) * lambda.in_matrix(nb_entries,total_entries-nb_entries) - o1.in(o1_ids);
                                add(o1_rep.in(unique_ids) == 0);
                                
                                // Representation of o2 with convex combination
                                Constraint<> o2_rep(pair.first+"_o2_rep");
                                o2_rep = bounds2.in_ignore_ith(nb_entries, 1, inst_partition_lambda).in_matrix(nb_entries,1) * lambda.in_matrix(nb_entries,total_entries-nb_entries) - o2.in(o2_ids);
                                add(o2_rep.in(unique_ids) == 0);
                                
                                // Linking partition variables1 with lambda
                                if(model_type == "lambda_II"){
                                    if(num_partns1 > 1) {
                                        Constraint<> on_link_lambda1(pair.first+"_on_link_lambda1_II");
                                        on_link_lambda1 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef1.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef1.in_matrix(nb_entries,total_entries-nb_entries) - binvar1->in_ignore_ith(nb_entries_v1,nb_entries_v2,on_coef1.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef1.in_matrix(nb_entries,1);
                                        add(on_link_lambda1.in(indices(unique_ids,const_idx1)) <= 0);
                                    }
                                    if(num_partns2 > 1) {
                                        Constraint<> on_link_lambda2(pair.first+"_on_link_lambda2_II");
                                        on_link_lambda2 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef2.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef2.in_matrix(nb_entries,total_entries-nb_entries) - binvar2->in_ignore_ith(0,nb_entries_v1,on_coef2.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef2.in_matrix(nb_entries,1);
                                        add(on_link_lambda2.in(indices(unique_ids,const_idx2)) <= 0);
                                        
                                    }
                                }
                                else{
                                    if(num_partns1 > 1) {
                                        Constraint<> on_link_lambda1(pair.first+"_on_link_lambda1_III");
                                        on_link_lambda1 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef1.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef1.in_matrix(nb_entries,total_entries-nb_entries) - binvar1->in_ignore_ith(nb_entries_v1,nb_entries_v2,on_coef1.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef1.in_matrix(nb_entries,1);
                                        add(on_link_lambda1.in(indices(unique_ids,const_idx1)) <= 0);
                                        
                                        
                                    }
                                    if(num_partns2 > 1) {
                                        Constraint<> on_link_lambda2(pair.first+"_on_link_lambda2_III");
                                        on_link_lambda2 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef2.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef2.in_matrix(nb_entries,total_entries-nb_entries) - binvar2->in_ignore_ith(0,nb_entries_v1,on_coef2.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef2.in_matrix(nb_entries,1);
                                        add(on_link_lambda2.in(indices(unique_ids,const_idx2)) <= 0);
                                    }
                                }
                                // sum over lambda
                                Constraint<> lambdaSum(pair.first+"_lambdaSum");
                                lambdaSum = sum(lambda.in_matrix(nb_entries,total_entries-nb_entries));
                                add(lambdaSum.in(unique_ids) == 1);
                            }
                        }
                    }
                    else if(binvar_ptr1 !=_vars_name.end() && binvar_ptr2 ==_vars_name.end()){ //if first variable v1 has been partitioned before
                        DebugOn("<<<<<<<<<< THIS IS NOT SEEN LIFT -> DOUBLE -> SEEN FIRST BINARY -> DIFF CORE VARS <<<<<<<<<<<" << endl);
                        
                        indices partns1("partns1");
                        for (int i = 0; i < num_partns1 ; ++i)
                        {
                            partns1.add(name1+ "{" +to_string(i+1) + "}");
                        }
                        
                        var<int> on2(name2+"_binary",0,1);
                        indices partns2("partns2");
                        for (int i = 0; i < num_partns2 ; ++i)
                        {
                            partns2.add(name2+ "{" + to_string(i+1) + "}");
                        }
                        add(on2.in(o2_ids_uq,partns2));
                        
                        auto binvar1 = static_pointer_cast<var<int>>(binvar_ptr2->second);
                        param<int> lb1("lb1"), ub1("ub1");
                        lb1.in(o1_ids_uq,partns1);
                        ub1.in(o1_ids_uq,partns1);
                        lb1.set_val(0), ub1.set_val(1);
                        auto added1 = binvar1->add_bounds(lb1,ub1);
                        reindex_vars();
                        
                        auto nb_entries_v1 = o1_ids.get_nb_entries();
                        auto nb_entries_v2 = o2_ids.get_nb_entries();
                        auto nb_entries = unique_ids.get_nb_entries();
                        
                        if(!added1.empty()){
                            Constraint<> onSum1(o1._name+"_binarySum");
                            onSum1 = sum(binvar1->in(added1).in_matrix(nb_entries_v1,1));
                            auto vset1 = added1.from_ith(0,nb_entries_v1);
                            vset1.filter_refs(vset1.get_unique_refs());
                            add(onSum1.in(vset1) == 1);
                        }
                        
                        Constraint<> onSum2(o2._name+"_binarySum");
                        onSum2 = sum(on2.in_matrix(nb_entries_v2,1));
                        add(onSum2.in(o2_ids_uq) == 1);
                        
                        if(model_type == "on/off"){//if on/off is chosen
                            var<int> on(name1+name2+"_binary",0,1);
                            
                            indices partns("partns");
                            partns = indices(partns1,partns2);
                            auto inst_partition = indices(unique_ids,partns);
                            add(on.in(inst_partition));
                            auto total_entries = inst_partition.get_nb_entries();
                            
                            Constraint<> onLink1(pair.first+"_binaryLink1");
                            onLink1 = binvar1->from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) - on;
                            add(onLink1.in(inst_partition) >= 0);
                            
                            Constraint<> onLink2(pair.first+"_binaryLink2");
                            onLink2 = on2.in_ignore_ith(nb_entries_v2,1,inst_partition.ignore_ith(0,nb_entries_v1)) - on;
                            add(onLink2.in(inst_partition) >= 0);
                            
                            Constraint<> onLink3(pair.first+"_binaryLink3");
                            onLink3 = binvar1->from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) + on2.in_ignore_ith(nb_entries_v2,1,inst_partition.ignore_ith(0,nb_entries_v1)) - 1 - on;
                            add(onLink3.in(inst_partition) <= 0);
                            
                            Constraint<> onSumComb(pair.first+"_binarySum");
                            onSumComb = sum(on.in_matrix(nb_entries,total_entries-nb_entries));
                            add(onSumComb.in(unique_ids) == 1);
                            
                            add_on_off_McCormick_refined(pair.first, vlift.in(unique_ids), o1.in(o1_ids), o2.in(o2_ids), on);
                        }
                        
                        
                        else{ //means it is one of the lambda formulations
                            
                            //difference is this has one more partition index
                            indices partns1_lambda("partns1_lambda");
                            for (int i = 0; i < num_partns1+1; ++i)
                            {
                                partns1_lambda.add(name1+ "{" +to_string(i+1) + "}");
                            }
                            
                            indices partns2_lambda("partns2_lambda");
                            for (int i = 0; i < num_partns2+1; ++i)
                            {
                                partns2_lambda.add(name2+ "{" +to_string(i+1) + "}");
                            }
                            
                            indices partns_lambda("partns_lambda");
                            partns_lambda = indices(partns1_lambda,partns2_lambda);
                            auto inst_partition_lambda = indices(unique_ids,partns_lambda);
                            auto inst_partition_bounds1 = indices(unique_ids,partns1_lambda);
                            auto inst_partition_bounds2 = indices(unique_ids,partns2_lambda);
                            
                            // Convex combination variables
                            var<> lambda(name1+name2+"_lambda",pos_);
                            add(lambda.in(inst_partition_lambda));
                            
                            /** Parameters */
                            // Bounds on variable v1 & v2
                            param<> bounds1(name1+"_bounds1");
                            bounds1.in(inst_partition_bounds1);
                            
                            param<> bounds2(name2+"_bounds2");
                            bounds2.in(inst_partition_bounds2);
                            
                            // Function values on the extreme points
                            param<> EP(name1+name2+"_grid_values");
                            EP.in(inst_partition_lambda);
                            auto total_entries = inst_partition_lambda.get_nb_entries();
                            
                            size_t nb_ins = vlift.in(unique_ids).get_nb_inst();
                            auto o1_global_lb = o1.get_lb();
                            auto increment1 = (o1.get_ub() - o1_global_lb)/num_partns1;
                            
                            auto o2_global_lb = o2.get_lb();
                            auto increment2 = (o2.get_ub() - o2_global_lb)/num_partns2;
                            
                            // fill bounds1 and function values
                            for (int i=0 ; i<num_partns1+1; ++i) {
                                auto bound_partn1 = o1_global_lb + increment1*i;
                                bound_partn1.eval_all();
                                for (size_t inst = 0; inst< nb_ins; inst++){
                                    auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                    auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                    string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"}";
                                    bounds1.set_val(cur_idx,bound_partn1.eval(inst));
                                    for(int j=0; j<num_partns2+1; ++j){
                                        auto bound_partn2 = o2_global_lb + increment2*j;
                                        bound_partn2.eval_all();
                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"}";
                                        EP.set_val(cur_idx,(bound_partn1.eval(inst)*bound_partn2.eval(inst)));
                                    }
                                }
                            }
                            // fill bounds2
                            for (int i=0 ; i<num_partns2+1; ++i) {
                                auto bound_partn2 = o2_global_lb + increment2*i;
                                bound_partn2.eval_all();
                                for (size_t inst = 0; inst< nb_ins; inst++){
                                    auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                    auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                    string cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"}";
                                    bounds2.set_val(cur_idx,bound_partn2.eval(inst));
                                }
                            }
                            
                            // Lambda coefficient matrix when linking with partition variables
                            param<> lambda_coef1(name1+"_lambda_linking_coefficients1");
                            param<> lambda_coef2(name2+"_lambda_linking_coefficients2");
                            
                            // Partition coefficient matrix when linking with lambda variables
                            param<> on_coef1(name1+"_partition_linking_coefficients1");
                            param<> on_coef2(name2+"_partition_linking_coefficients2");
                            
                            // create constraint indices
                            indices const_idx1("const_idx1");
                            indices const_idx2("const_idx2");
                            
                            if(model_type == "lambda_II"){
                                
                                //fill constraint indices
                                for (int i=0; i<num_partns1+1; ++i){
                                    const_idx1.add(to_string(i+1));
                                }
                                
                                //fill constraint indices
                                for (int i=0; i<num_partns2+1; ++i){
                                    const_idx2.add(to_string(i+1));
                                }
                                
                                // Lambda coefficient matrix when linking with partition variables
                                if(num_partns1 > 1) lambda_coef1.in(indices(inst_partition_lambda, const_idx1));
                                if(num_partns2 > 1) lambda_coef2.in(indices(inst_partition_lambda, const_idx2));
                                
                                // Lambda coefficient matrix when linking with partition variables
                                if(num_partns1 > 1)  on_coef1.in(indices(unique_ids, partns1, const_idx1));
                                if(num_partns2 > 1)  on_coef2.in(indices(unique_ids, partns2, const_idx2));
                                
                                // fill lambda_coef1 and lambda_coef2
                                for (size_t inst = 0; inst< nb_ins; inst++){
                                    auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                    auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                    for (int i=0 ; i<num_partns1+1; ++i) {
                                        for (int j=0 ; j<num_partns2+1; ++j) {
                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                            if(num_partns1 > 1)  lambda_coef1.set_val(cur_idx,1);
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(j+1);
                                            if(num_partns2 > 1)  lambda_coef2.set_val(cur_idx,1);
                                        }
                                    }
                                }
                                
                                // fill on_coef1 and on_coef2
                                for (size_t inst = 0; inst< nb_ins; inst++){
                                    auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                    auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                    string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                    if(num_partns1 > 1)  on_coef1.set_val(cur_idx,1);
                                    cur_idx = cur_var_idx+","+name2+"{"+to_string(1)+"},"+to_string(1);
                                    if(num_partns2 > 1)  on_coef2.set_val(cur_idx,1);
                                    
                                    if(num_partns1 > 1) {
                                        for (int i=1 ; i<num_partns1; ++i) {
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(i)+"},"+to_string(i+1);
                                            on_coef1.set_val(cur_idx,1);
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(i+1);
                                            on_coef1.set_val(cur_idx,1);
                                        }
                                    }
                                    if(num_partns2 > 1) {
                                        for (int i=1 ; i<num_partns2; ++i) {
                                            cur_idx = cur_var_idx+","+name2+"{"+to_string(i)+"},"+to_string(i+1);
                                            on_coef2.set_val(cur_idx,1);
                                            cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"},"+to_string(i+1);
                                            on_coef2.set_val(cur_idx,1);
                                        }
                                    }
                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1)+"},"+to_string(num_partns1+1);
                                    if(num_partns1 > 1)  on_coef1.set_val(cur_idx,1);
                                    cur_idx = cur_var_idx+","+name2+"{"+to_string(num_partns2)+"},"+to_string(num_partns2+1);
                                    if(num_partns1 > 2)  on_coef2.set_val(cur_idx,1);
                                }
                            }
                            
                            
                            else /*means model_type == "lambda_III" */{
                                
                                //fill constraint indices
                                for (int i=0; i<(num_partns1-2)*2+2; ++i){
                                    const_idx1.add(to_string(i+1));
                                }
                                
                                //fill constraint indices
                                for (int i=0; i<(num_partns2-2)*2+2; ++i){
                                    const_idx2.add(to_string(i+1));
                                }
                                
                                // Lambda coefficient matrix when linking with partition variables
                                if(num_partns1 > 1)  lambda_coef1.in(indices(inst_partition_lambda, const_idx1));
                                if(num_partns2 > 1)  lambda_coef2.in(indices(inst_partition_lambda, const_idx2));
                                
                                // Partition coefficient matrix when linking with lambda variables
                                if(num_partns1 > 1)  on_coef1.in(indices(unique_ids, partns1, const_idx1));
                                if(num_partns2 > 1)  on_coef2.in(indices(unique_ids, partns2, const_idx2));
                                
                                // fill lambda_coef1 and lambda_coef2
                                for (size_t inst = 0; inst< nb_ins; inst++){
                                    auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                    auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                    if(num_partns1 > 1){
                                        for (int j=0; j<num_partns2+1; ++j) {
                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(1);
                                            lambda_coef1.set_val(cur_idx,1);
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string((num_partns1-2)*2+2);
                                            lambda_coef1.set_val(cur_idx,1);
                                        }
                                        
                                        for (int i=1 ; i<(num_partns1-2)*2+1; i=i+2) {
                                            for (int j=(i-1)/2 + 2; j<num_partns1+1; ++j) {
                                                for(int k=0; k<num_partns2+1; ++k){
                                                    string cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+name2+"{"+to_string(k+1)+"},"+to_string(i+1);
                                                    lambda_coef1.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+name2+"{"+to_string(k+1)+"},"+to_string(i+2);
                                                    lambda_coef1.set_val(cur_idx,-1);
                                                }
                                            }
                                        }
                                    }
                                    if(num_partns2 > 1) {
                                        for (int i=0; i<num_partns1+1; ++i) {
                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(1)+"},"+to_string(1);
                                            lambda_coef2.set_val(cur_idx,1);
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(num_partns2+1)+"},"+to_string((num_partns2-2)*2+2);
                                            lambda_coef2.set_val(cur_idx,1);
                                        }
                                        
                                        for (int i=1 ; i<(num_partns2-2)*2+1; i=i+2) {
                                            for (int j=(i-1)/2 + 2; j<num_partns2+1; ++j) {
                                                for(int k=0; k<num_partns1+1; ++k){
                                                    string cur_idx = cur_var_idx+","+name1+"{"+to_string(k+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                    lambda_coef2.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(k+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                    lambda_coef2.set_val(cur_idx,-1);
                                                }
                                            }
                                        }
                                    }
                                }
                                
                                
                                
                                // fill on_coef1 and on_coef2
                                for (size_t inst = 0; inst< nb_ins; inst++){
                                    auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                    auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                    string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                    if(num_partns1 > 1) {
                                        on_coef1.set_val(cur_idx,1);
                                        
                                        for (int i=1; i<num_partns1; ++i) {
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(2);
                                            on_coef1.set_val(cur_idx, 1);
                                        }
                                        
                                        for (int i=2 ; i<(num_partns1-2)*2+2; i=i+2) {
                                            for (int j=i/2+1; j<num_partns1; ++j) {
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                on_coef1.set_val(cur_idx,-1);
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                on_coef1.set_val(cur_idx,1);
                                            }
                                        }
                                    }
                                    
                                    if(num_partns2 > 1) {
                                        cur_idx = cur_var_idx+","+name2+"{"+to_string(1)+"},"+to_string(1);
                                        on_coef2.set_val(cur_idx,1);
                                        for (int i=1; i<num_partns2; ++i) {
                                            cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"},"+to_string(2);
                                            on_coef2.set_val(cur_idx, 1);
                                        }
                                        
                                        for (int i=2 ; i<(num_partns2-2)*2+2; i=i+2) {
                                            for (int j=i/2+1; j<num_partns2; ++j) {
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                on_coef2.set_val(cur_idx,-1);
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                on_coef2.set_val(cur_idx,1);
                                            }
                                        }
                                    }
                                }
                            }
                            
                            
                            /** Constraints */
                            // Representation of the bilinear term with convex combination
                            Constraint<> bln_rep(pair.first+"_bln_rep");
                            bln_rep = EP.in_matrix(nb_entries,total_entries-nb_entries)*lambda.in_matrix(nb_entries,total_entries-nb_entries) - vlift.in(unique_ids);
                            add(bln_rep.in(unique_ids) == 0);
                            
                            // Representation of o1 with convex combination
                            Constraint<> o1_rep(pair.first+"_o1_rep");
                            
                            o1_rep = bounds1.from_ith(0,inst_partition_lambda).in_matrix(nb_entries, 1) * lambda.in_matrix(nb_entries,total_entries-nb_entries) - o1.in(o1_ids);
                            add(o1_rep.in(unique_ids) == 0);
                            
                            // Representation of o2 with convex combination
                            Constraint<> o2_rep(pair.first+"_o2_rep");
                            
                            o2_rep = bounds2.in_ignore_ith(nb_entries, 1, inst_partition_lambda).in_matrix(nb_entries,1) * lambda.in_matrix(nb_entries,total_entries-nb_entries) - o2.in(o2_ids);
                            add(o2_rep.in(unique_ids) == 0);
                            
                            // Linking partition variables1 with lambda
                            if(model_type == "lambda_II"){
                                if(num_partns1 > 1) {
                                    Constraint<> on_link_lambda1(pair.first+"_on_link_lambda1_II");
                                    on_link_lambda1 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef1.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef1.in_matrix(nb_entries,total_entries-nb_entries) - binvar1->in_ignore_ith(nb_entries_v1,nb_entries_v2,on_coef1.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef1.in_matrix(nb_entries,1);
                                    add(on_link_lambda1.in(indices(unique_ids,const_idx1)) <= 0);
                                }
                                if(num_partns2 > 1) {
                                    Constraint<> on_link_lambda2(pair.first+"_on_link_lambda2_II");
                                    on_link_lambda2 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef2.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef2.in_matrix(nb_entries,total_entries-nb_entries) - on2.in_ignore_ith(0,nb_entries_v1,on_coef2.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef2.in_matrix(nb_entries,1);
                                    add(on_link_lambda2.in(indices(unique_ids,const_idx2)) <= 0);
                                    
                                }
                            }
                            else{
                                if(num_partns1 > 1) {
                                    Constraint<> on_link_lambda1(pair.first+"_on_link_lambda1_III");
                                    on_link_lambda1 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef1.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef1.in_matrix(nb_entries,total_entries-nb_entries) - binvar1->in_ignore_ith(nb_entries_v1,nb_entries_v2,on_coef1.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef1.in_matrix(nb_entries,1);
                                    add(on_link_lambda1.in(indices(unique_ids,const_idx1)) <= 0);
                                }
                                if(num_partns2 > 1) {
                                    Constraint<> on_link_lambda2(pair.first+"_on_link_lambda2_III");
                                    on_link_lambda2 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef2.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef2.in_matrix(nb_entries,total_entries-nb_entries) - on2.in_ignore_ith(0,nb_entries_v1,on_coef2.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef2.in_matrix(nb_entries,1);
                                    add(on_link_lambda2.in(indices(unique_ids,const_idx2)) <= 0);
                                }
                            }
                            // sum over lambda
                            Constraint<> lambdaSum(pair.first+"_lambdaSum");
                            lambdaSum = sum(lambda.in_matrix(nb_entries,total_entries-nb_entries));
                            add(lambdaSum.in(unique_ids) == 1);
                        }
                        
                    }
                    else if(binvar_ptr1 ==_vars_name.end() && binvar_ptr2 !=_vars_name.end()){ //means v2 has been partitioned before)
                        DebugOn("<<<<<<<<<< THIS IS NOT SEEN LIFT -> DOUBLE -> SEEN SECOND BINARY -> DIFF CORE VARS <<<<<<<<<<<" << endl);
                        
                        var<int> on1(name1+"_binary",0,1);
                        indices partns1("partns1");
                        for (int i = 0; i < num_partns1 ; ++i)
                        {
                            partns1.add(name1+ "{" +to_string(i+1) + "}");
                        }
                        add(on1.in(o1_ids_uq,partns1));
                        
                        indices partns2("partns2");
                        for (int i = 0; i < num_partns2 ; ++i)
                        {
                            partns2.add(name2+ "{" + to_string(i+1) + "}");
                        }
                        
                        auto binvar2 = static_pointer_cast<var<int>>(binvar_ptr2->second);
                        param<int> lb2("lb2"), ub2("ub2");
                        lb2.in(o2_ids_uq,partns2);
                        ub2.in(o2_ids_uq,partns2);
                        lb2.set_val(0), ub2.set_val(1);
                        auto added2 = binvar2->add_bounds(lb2,ub2);
                        reindex_vars();
                        
                        auto nb_entries_v1 = o1_ids.get_nb_entries();
                        auto nb_entries_v2 = o2_ids.get_nb_entries();
                        auto nb_entries = unique_ids.get_nb_entries();
                        
                        if(!added2.empty()){
                            Constraint<> onSum2(o2._name+"_binarySum");
                            onSum2 = sum(binvar2->in(added2).in_matrix(nb_entries_v2,1));
                            auto vset2 = added2.from_ith(0,nb_entries_v2);
                            vset2.filter_refs(vset2.get_unique_refs());
                            add(onSum2.in(vset2) == 1);
                        }
                        
                        Constraint<> onSum1(o1._name+"_binarySum");
                        onSum1 = sum(on1.in_matrix(nb_entries_v1,1));
                        add(onSum1.in(o1_ids_uq) == 1);
                        
                        if(model_type == "on/off"){//if on/off is chosen
                            var<int> on(name1+name2+"_binary",0,1);
                            
                            indices partns("partns");
                            partns = indices(partns1,partns2);
                            auto inst_partition = indices(unique_ids,partns);
                            add(on.in(inst_partition));
                            auto total_entries = inst_partition.get_nb_entries();
                            
                            Constraint<> onLink1(pair.first+"_binaryLink1");
                            onLink1 = on1.from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) - on;
                            add(onLink1.in(inst_partition) >= 0);
                            
                            Constraint<> onLink2(pair.first+"_binaryLink2");
                            onLink2 = binvar2->in_ignore_ith(nb_entries_v2,1,inst_partition.ignore_ith(0,nb_entries_v1)) - on;
                            add(onLink2.in(inst_partition) >= 0);
                            
                            Constraint<> onLink3(pair.first+"_binaryLink3");
                            onLink3 = on1.from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) + binvar2->in_ignore_ith(nb_entries_v2,1,inst_partition.ignore_ith(0,nb_entries_v1)) - 1 - on;
                            add(onLink3.in(inst_partition) <= 0);
                            
                            Constraint<> onSumComb(pair.first+"_binarySum");
                            onSumComb = sum(on.in_matrix(nb_entries,total_entries-nb_entries));
                            add(onSumComb.in(unique_ids) == 1);
                            
                            add_on_off_McCormick_refined(pair.first, vlift.in(unique_ids), o1.in(o1_ids), o2.in(o2_ids), on);
                        }
                        
                        
                        else{ //means it is one of the lambda formulations
                            
                            //difference is this has one more partition index
                            indices partns1_lambda("partns1_lambda");
                            for (int i = 0; i < num_partns1+1; ++i)
                            {
                                partns1_lambda.add(name1+ "{" +to_string(i+1) + "}");
                            }
                            
                            indices partns2_lambda("partns2_lambda");
                            for (int i = 0; i < num_partns2+1; ++i)
                            {
                                partns2_lambda.add(name2+ "{" +to_string(i+1) + "}");
                            }
                            
                            indices partns_lambda("partns_lambda");
                            partns_lambda = indices(partns1_lambda,partns2_lambda);
                            auto inst_partition_lambda = indices(unique_ids,partns_lambda);
                            auto inst_partition_bounds1 = indices(unique_ids,partns1_lambda);
                            auto inst_partition_bounds2 = indices(unique_ids,partns2_lambda);
                            
                            // Convex combination variables
                            var<> lambda(name1+name2+"_lambda",pos_);
                            add(lambda.in(inst_partition_lambda));
                            
                            /** Parameters */
                            // Bounds on variable v1 & v2
                            param<> bounds1(name1+"_bounds1");
                            bounds1.in(inst_partition_bounds1);
                            
                            param<> bounds2(name2+"_bounds2");
                            bounds2.in(inst_partition_bounds2);
                            
                            // Function values on the extreme points
                            param<> EP(name1+name2+"_grid_values");
                            EP.in(inst_partition_lambda);
                            auto total_entries = inst_partition_lambda.get_nb_entries();
                            
                            size_t nb_ins = vlift.in(unique_ids).get_nb_inst();
                            auto o1_global_lb = o1.get_lb();
                            auto increment1 = (o1.get_ub() - o1_global_lb)/num_partns1;
                            
                            auto o2_global_lb = o2.get_lb();
                            auto increment2 = (o2.get_ub() - o2_global_lb)/num_partns2;
                            
                            // fill bounds1 and function values
                            for (int i=0 ; i<num_partns1+1; ++i) {
                                auto bound_partn1 = o1_global_lb + increment1*i;
                                bound_partn1.eval_all();
                                for (size_t inst = 0; inst< nb_ins; inst++){
                                    auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                    auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                    string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"}";
                                    bounds1.set_val(cur_idx,bound_partn1.eval(inst));
                                    for(int j=0; j<num_partns2+1; ++j){
                                        auto bound_partn2 = o2_global_lb + increment2*j;
                                        bound_partn2.eval_all();
                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"}";
                                        EP.set_val(cur_idx,(bound_partn1.eval(inst)*bound_partn2.eval(inst)));
                                    }
                                }
                            }
                            // fill bounds2
                            for (int i=0 ; i<num_partns2+1; ++i) {
                                auto bound_partn2 = o2_global_lb + increment2*i;
                                bound_partn2.eval_all();
                                for (size_t inst = 0; inst< nb_ins; inst++){
                                    auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                    auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                    string cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"}";
                                    bounds2.set_val(cur_idx,bound_partn2.eval(inst));
                                }
                            }
                            
                            // Lambda coefficient matrix when linking with partition variables
                            param<> lambda_coef1(name1+"_lambda_linking_coefficients1");
                            param<> lambda_coef2(name2+"_lambda_linking_coefficients2");
                            
                            // Partition coefficient matrix when linking with lambda variables
                            param<> on_coef1(name1+"_partition_linking_coefficients1");
                            param<> on_coef2(name2+"_partition_linking_coefficients2");
                            
                            // create constraint indices
                            indices const_idx1("const_idx1");
                            indices const_idx2("const_idx2");
                            
                            if(model_type == "lambda_II"){
                                
                                //fill constraint indices
                                for (int i=0; i<num_partns1+1; ++i){
                                    const_idx1.add(to_string(i+1));
                                }
                                
                                //fill constraint indices
                                for (int i=0; i<num_partns2+1; ++i){
                                    const_idx2.add(to_string(i+1));
                                }
                                
                                // Lambda coefficient matrix when linking with partition variables
                                if(num_partns1 > 1) lambda_coef1.in(indices(inst_partition_lambda, const_idx1));
                                if(num_partns2 > 1) lambda_coef2.in(indices(inst_partition_lambda, const_idx2));
                                
                                // Lambda coefficient matrix when linking with partition variables
                                if(num_partns1 > 1)  on_coef1.in(indices(unique_ids, partns1, const_idx1));
                                if(num_partns2 > 1)  on_coef2.in(indices(unique_ids, partns2, const_idx2));
                                
                                // fill lambda_coef1 and lambda_coef2
                                for (size_t inst = 0; inst< nb_ins; inst++){
                                    auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                    auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                    for (int i=0 ; i<num_partns1+1; ++i) {
                                        for (int j=0 ; j<num_partns2+1; ++j) {
                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                            if(num_partns1 > 1)  lambda_coef1.set_val(cur_idx,1);
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(j+1);
                                            if(num_partns2 > 1)  lambda_coef2.set_val(cur_idx,1);
                                        }
                                    }
                                }
                                
                                // fill on_coef1 and on_coef2
                                for (size_t inst = 0; inst< nb_ins; inst++){
                                    auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                    auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                    string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                    if(num_partns1 > 1)  on_coef1.set_val(cur_idx,1);
                                    cur_idx = cur_var_idx+","+name2+"{"+to_string(1)+"},"+to_string(1);
                                    if(num_partns2 > 1)  on_coef2.set_val(cur_idx,1);
                                    
                                    if(num_partns1 > 1) {
                                        for (int i=1 ; i<num_partns1; ++i) {
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(i)+"},"+to_string(i+1);
                                            on_coef1.set_val(cur_idx,1);
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(i+1);
                                            on_coef1.set_val(cur_idx,1);
                                        }
                                    }
                                    if(num_partns2 > 1) {
                                        for (int i=1 ; i<num_partns2; ++i) {
                                            cur_idx = cur_var_idx+","+name2+"{"+to_string(i)+"},"+to_string(i+1);
                                            on_coef2.set_val(cur_idx,1);
                                            cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"},"+to_string(i+1);
                                            on_coef2.set_val(cur_idx,1);
                                        }
                                    }
                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1)+"},"+to_string(num_partns1+1);
                                    if(num_partns1 > 1)  on_coef1.set_val(cur_idx,1);
                                    cur_idx = cur_var_idx+","+name2+"{"+to_string(num_partns2)+"},"+to_string(num_partns2+1);
                                    if(num_partns1 > 2)  on_coef2.set_val(cur_idx,1);
                                }
                            }
                            
                            
                            else /*means model_type == "lambda_III" */{
                                
                                //fill constraint indices
                                for (int i=0; i<(num_partns1-2)*2+2; ++i){
                                    const_idx1.add(to_string(i+1));
                                }
                                
                                //fill constraint indices
                                for (int i=0; i<(num_partns2-2)*2+2; ++i){
                                    const_idx2.add(to_string(i+1));
                                }
                                
                                // Lambda coefficient matrix when linking with partition variables
                                if(num_partns1 > 1)  lambda_coef1.in(indices(inst_partition_lambda, const_idx1));
                                if(num_partns2 > 1)  lambda_coef2.in(indices(inst_partition_lambda, const_idx2));
                                
                                // Partition coefficient matrix when linking with lambda variables
                                if(num_partns1 > 1)  on_coef1.in(indices(unique_ids, partns1, const_idx1));
                                if(num_partns2 > 1)  on_coef2.in(indices(unique_ids, partns2, const_idx2));
                                
                                // fill lambda_coef1 and lambda_coef2
                                for (size_t inst = 0; inst< nb_ins; inst++){
                                    auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                    auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                    if(num_partns1 > 1){
                                        for (int j=0; j<num_partns2+1; ++j) {
                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(1);
                                            lambda_coef1.set_val(cur_idx,1);
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string((num_partns1-2)*2+2);
                                            lambda_coef1.set_val(cur_idx,1);
                                        }
                                        
                                        for (int i=1 ; i<(num_partns1-2)*2+1; i=i+2) {
                                            for (int j=(i-1)/2 + 2; j<num_partns1+1; ++j) {
                                                for(int k=0; k<num_partns2+1; ++k){
                                                    string cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+name2+"{"+to_string(k+1)+"},"+to_string(i+1);
                                                    lambda_coef1.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+name2+"{"+to_string(k+1)+"},"+to_string(i+2);
                                                    lambda_coef1.set_val(cur_idx,-1);
                                                }
                                            }
                                        }
                                    }
                                    if(num_partns2 > 1) {
                                        for (int i=0; i<num_partns1+1; ++i) {
                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(1)+"},"+to_string(1);
                                            lambda_coef2.set_val(cur_idx,1);
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(num_partns2+1)+"},"+to_string((num_partns2-2)*2+2);
                                            lambda_coef2.set_val(cur_idx,1);
                                        }
                                        
                                        for (int i=1 ; i<(num_partns2-2)*2+1; i=i+2) {
                                            for (int j=(i-1)/2 + 2; j<num_partns2+1; ++j) {
                                                for(int k=0; k<num_partns1+1; ++k){
                                                    string cur_idx = cur_var_idx+","+name1+"{"+to_string(k+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                    lambda_coef2.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(k+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                    lambda_coef2.set_val(cur_idx,-1);
                                                }
                                            }
                                        }
                                    }
                                }
                                
                                
                                
                                // fill on_coef1 and on_coef2
                                for (size_t inst = 0; inst< nb_ins; inst++){
                                    auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                    auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                    string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                    if(num_partns1 > 1) {
                                        on_coef1.set_val(cur_idx,1);
                                        
                                        for (int i=1; i<num_partns1; ++i) {
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(2);
                                            on_coef1.set_val(cur_idx, 1);
                                        }
                                        
                                        for (int i=2 ; i<(num_partns1-2)*2+2; i=i+2) {
                                            for (int j=i/2+1; j<num_partns1; ++j) {
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                on_coef1.set_val(cur_idx,-1);
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                on_coef1.set_val(cur_idx,1);
                                            }
                                        }
                                    }
                                    
                                    if(num_partns2 > 1) {
                                        cur_idx = cur_var_idx+","+name2+"{"+to_string(1)+"},"+to_string(1);
                                        on_coef2.set_val(cur_idx,1);
                                        for (int i=1; i<num_partns2; ++i) {
                                            cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"},"+to_string(2);
                                            on_coef2.set_val(cur_idx, 1);
                                        }
                                        
                                        for (int i=2 ; i<(num_partns2-2)*2+2; i=i+2) {
                                            for (int j=i/2+1; j<num_partns2; ++j) {
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                on_coef2.set_val(cur_idx,-1);
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                on_coef2.set_val(cur_idx,1);
                                            }
                                        }
                                    }
                                }
                            }
                            
                            
                            /** Constraints */
                            // Representation of the bilinear term with convex combination
                            Constraint<> bln_rep(pair.first+"_bln_rep");
                            bln_rep = EP.in_matrix(nb_entries,total_entries-nb_entries)*lambda.in_matrix(nb_entries,total_entries-nb_entries) - vlift.in(unique_ids);
                            add(bln_rep.in(unique_ids) == 0);
                            
                            // Representation of o1 with convex combination
                            Constraint<> o1_rep(pair.first+"_o1_rep");
                            
                            o1_rep = bounds1.from_ith(0,inst_partition_lambda).in_matrix(nb_entries, 1) * lambda.in_matrix(nb_entries,total_entries-nb_entries) - o1.in(o1_ids);
                            add(o1_rep.in(unique_ids) == 0);
                            
                            // Representation of o2 with convex combination
                            Constraint<> o2_rep(pair.first+"_o2_rep");
                            
                            o2_rep = bounds2.in_ignore_ith(nb_entries, 1, inst_partition_lambda).in_matrix(nb_entries,1) * lambda.in_matrix(nb_entries,total_entries-nb_entries) - o2.in(o2_ids);
                            add(o2_rep.in(unique_ids) == 0);
                            
                            // Linking partition variables1 with lambda
                            if(model_type == "lambda_II"){
                                if(num_partns1 > 1) {
                                    Constraint<> on_link_lambda1(pair.first+"_on_link_lambda1_II");
                                    on_link_lambda1 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef1.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef1.in_matrix(nb_entries,total_entries-nb_entries) - on1.in_ignore_ith(nb_entries_v1,nb_entries_v2,on_coef1.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef1.in_matrix(nb_entries,1);
                                    add(on_link_lambda1.in(indices(unique_ids,const_idx1)) <= 0);
                                }
                                if(num_partns2 > 1) {
                                    Constraint<> on_link_lambda2(pair.first+"_on_link_lambda2_II");
                                    on_link_lambda2 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef2.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef2.in_matrix(nb_entries,total_entries-nb_entries) - binvar2->in_ignore_ith(0,nb_entries_v1,on_coef2.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef2.in_matrix(nb_entries,1);
                                    add(on_link_lambda2.in(indices(unique_ids,const_idx2)) <= 0);
                                    
                                }
                            }
                            else{
                                if(num_partns1 > 1) {
                                    Constraint<> on_link_lambda1(pair.first+"_on_link_lambda1_III");
                                    on_link_lambda1 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef1.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef1.in_matrix(nb_entries,total_entries-nb_entries) - on1.in_ignore_ith(nb_entries_v1,nb_entries_v2,on_coef1.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef1.in_matrix(nb_entries,1);
                                    add(on_link_lambda1.in(indices(unique_ids,const_idx1)) <= 0);
                                }
                                if(num_partns2 > 1) {
                                    Constraint<> on_link_lambda2(pair.first+"_on_link_lambda2_III");
                                    on_link_lambda2 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef2.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef2.in_matrix(nb_entries,total_entries-nb_entries) - binvar2->in_ignore_ith(0,nb_entries_v1,on_coef2.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef2.in_matrix(nb_entries,1);
                                    add(on_link_lambda2.in(indices(unique_ids,const_idx2)) <= 0);
                                }
                            }
                            // sum over lambda
                            Constraint<> lambdaSum(pair.first+"_lambdaSum");
                            lambdaSum = sum(lambda.in_matrix(nb_entries,total_entries-nb_entries));
                            add(lambdaSum.in(unique_ids) == 1);
                        }
                        
                    }
                    else{ //means both variables v1 and v2 haven't been partitioned
                        if(name1 == name2){
                            DebugOn("<<<<<<<<<< THIS IS NOT SEEN LIFT -> DOUBLE -> NOT SEEN BINARIES -> SAME CORE VARS <<<<<<<<<<<" << endl);
                            
                            var<int> on1(name1+"_binary",0,1);
                            indices partns1("partns1");
                            for (int i = 0; i < num_partns1 ; ++i)
                            {
                                partns1.add(name1+ "{" +to_string(i+1) + "}");
                            }
                            add(on1.in(union_ids(o1_ids_uq, o2_ids_uq),partns1));
                            indices partns2("partns2");
                            for (int i = 0; i < num_partns2 ; ++i)
                            {
                                partns2.add(name2+ "{" +to_string(i+1) + "}");
                            }
                            
                            auto nb_entries_v1 = o1_ids.get_nb_entries();
                            auto nb_entries_v2 = o2_ids.get_nb_entries();
                            auto nb_entries = unique_ids.get_nb_entries();
                            
                            Constraint<> onSum1(o1._name+"_binarySum");
                            onSum1 = sum(on1.in_matrix(nb_entries_v1,1));
                            add(onSum1.in(union_ids(o1_ids_uq,o2_ids_uq)) == 1);
                            
                            if(model_type == "on/off")
                            {
                                var<int> on(name1+name2+"_binary",0,1);
                                
                                indices partns("partns");
                                partns = indices(partns1,partns2);
                                auto inst_partition = indices(unique_ids,partns);
                                add(on.in(inst_partition));
                                auto total_entries = inst_partition.get_nb_entries();
                                
                                Constraint<> onLink1(pair.first+"_binaryLink1");
                                onLink1 = on1.from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) - on;
                                add(onLink1.in(inst_partition) >= 0);
                                
                                Constraint<> onLink2(pair.first+"_binaryLink2");
                                onLink2 = on1.in_ignore_ith(nb_entries_v1,1,inst_partition.ignore_ith(0,nb_entries_v1)) - on;
                                add(onLink2.in(inst_partition) >= 0);
                                
                                Constraint<> onLink3(pair.first+"_binaryLink3");
                                onLink3 = on1.from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) + on1.in_ignore_ith(nb_entries_v1,1,inst_partition.ignore_ith(0,nb_entries_v1)) - 1 - on;
                                add(onLink3.in(inst_partition) <= 0);
                                
                                Constraint<> onSumComb(pair.first+"_binarySum");
                                onSumComb = sum(on.in_matrix(nb_entries,total_entries-nb_entries));
                                add(onSumComb.in(unique_ids) == 1);
                                
                                add_on_off_McCormick_refined(pair.first, vlift.in(unique_ids), o1.in(o1_ids), o2.in(o2_ids), on);
                            }
                            
                            else{ //means it is one of the lambda formulations
                                
                                //difference is this has one more partition index
                                indices partns1_lambda("partns1_lambda");
                                for (int i = 0; i < num_partns1+1; ++i)
                                {
                                    partns1_lambda.add(name1+ "{" +to_string(i+1) + "}");
                                }
                                
                                indices partns2_lambda("partns2_lambda");
                                for (int i = 0; i < num_partns2+1; ++i)
                                {
                                    partns2_lambda.add(name2+ "{" +to_string(i+1) + "}");
                                }
                                
                                indices partns_lambda("partns_lambda");
                                partns_lambda = indices(partns1_lambda,partns2_lambda);
                                auto inst_partition_lambda = indices(unique_ids,partns_lambda);
                                auto inst_partition_bounds1 = indices(unique_ids,partns1_lambda);
                                auto inst_partition_bounds2 = indices(unique_ids,partns2_lambda);
                                
                                // Convex combination variables
                                var<> lambda(name1+name2+"_lambda",pos_);
                                add(lambda.in(inst_partition_lambda));
                                
                                /** Parameters */
                                // Bounds on variable v1 & v2
                                param<> bounds1(name1+"_bounds1");
                                bounds1.in(inst_partition_bounds1);
                                
                                param<> bounds2(name2+"_bounds2");
                                bounds2.in(inst_partition_bounds2);
                                
                                // Function values on the extreme points
                                param<> EP(name1+name2+"_grid_values");
                                EP.in(inst_partition_lambda);
                                auto total_entries = inst_partition_lambda.get_nb_entries();
                                
                                size_t nb_ins = vlift.in(unique_ids).get_nb_inst();
                                
                                auto o1_global_lb = o1.get_lb();
                                auto increment1 = (o1.get_ub() - o1_global_lb)/num_partns1;
                                
                                auto o2_global_lb = o2.get_lb();
                                auto increment2 = (o2.get_ub() - o2_global_lb)/num_partns2;
                                
                                // fill bounds and function values
                                for (int i=0 ; i<num_partns1+1; ++i) {
                                    auto bound_partn1 = o1_global_lb + increment1*i;
                                    bound_partn1.eval_all();
                                    auto bound_partn2 = o2_global_lb + increment2*i;
                                    bound_partn2.eval_all();
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                        auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"}";
                                        bounds1.set_val(cur_idx,bound_partn1.eval(inst));
                                        bounds2.set_val(cur_idx,bound_partn2.eval(inst));
                                        for(int j=0; j<num_partns2+1; ++j){
                                            auto bound_partn2_temp = o2_global_lb + increment2*j;
                                            bound_partn2_temp.eval_all();
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"}";
                                            EP.set_val(cur_idx,(bound_partn1.eval(inst)*bound_partn2_temp.eval(inst)));
                                        }
                                    }
                                }
                                
                                // Lambda coefficient matrix when linking with partition variables
                                param<> lambda_coef1(name1+"_lambda_linking_coefficients1");
                                param<> lambda_coef2(name2+"_lambda_linking_coefficients2");
                                
                                // Partition coefficient matrix when linking with lambda variables
                                param<> on_coef1(name1+"_partition_linking_coefficients1");
                                param<> on_coef2(name2+"_partition_linking_coefficients2");
                                
                                // create constraint indices
                                indices const_idx1("const_idx1");
                                indices const_idx2("const_idx2");
                                
                                if(model_type == "lambda_II"){
                                    
                                    //fill constraint indices
                                    for (int i=0; i<num_partns1+1; ++i){
                                        const_idx1.add(to_string(i+1));
                                    }
                                    
                                    //fill constraint indices
                                    for (int i=0; i<num_partns2+1; ++i){
                                        const_idx2.add(to_string(i+1));
                                    }
                                    
                                    // Lambda coefficient matrix when linking with partition variables
                                    lambda_coef1.in(indices(inst_partition_lambda, const_idx1));
                                    lambda_coef2.in(indices(inst_partition_lambda, const_idx2));
                                    
                                    // Lambda coefficient matrix when linking with partition variables
                                    on_coef1.in(indices(unique_ids, partns1, const_idx1));
                                    on_coef2.in(indices(unique_ids, partns2, const_idx2));
                                    
                                    // fill lambda_coef1 and lambda_coef2
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                        auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                        for (int i=0 ; i<num_partns1+1; ++i) {
                                            for (int j=0 ; j<num_partns2+1; ++j) {
                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                lambda_coef1.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(j+1);
                                                lambda_coef2.set_val(cur_idx,1);
                                            }
                                        }
                                    }
                                    
                                    // fill on_coef1 and on_coef2
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                        auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                        on_coef1.set_val(cur_idx,1);
                                        cur_idx = cur_var_idx+","+name2+"{"+to_string(1)+"},"+to_string(1);
                                        on_coef2.set_val(cur_idx,1);
                                        for (int i=1 ; i<num_partns1; ++i) {
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(i)+"},"+to_string(i+1);
                                            on_coef1.set_val(cur_idx,1);
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(i+1);
                                            on_coef1.set_val(cur_idx,1);
                                        }
                                        for (int i=1 ; i<num_partns2; ++i) {
                                            cur_idx = cur_var_idx+","+name2+"{"+to_string(i)+"},"+to_string(i+1);
                                            on_coef2.set_val(cur_idx,1);
                                            cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"},"+to_string(i+1);
                                            on_coef2.set_val(cur_idx,1);
                                        }
                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1)+"},"+to_string(num_partns1+1);
                                        on_coef1.set_val(cur_idx,1);
                                        cur_idx = cur_var_idx+","+name2+"{"+to_string(num_partns2)+"},"+to_string(num_partns2+1);
                                        on_coef2.set_val(cur_idx,1);
                                    }
                                }
                                
                                
                                else /*means model_type == "lambda_III" */{
                                    
                                    //fill constraint indices
                                    for (int i=0; i<(num_partns1-2)*2+2; ++i){
                                        const_idx1.add(to_string(i+1));
                                    }
                                    
                                    //fill constraint indices
                                    for (int i=0; i<(num_partns2-2)*2+2; ++i){
                                        const_idx2.add(to_string(i+1));
                                    }
                                    
                                    // Lambda coefficient matrix when linking with partition variables
                                    lambda_coef1.in(indices(inst_partition_lambda, const_idx1));
                                    lambda_coef2.in(indices(inst_partition_lambda, const_idx2));
                                    
                                    // Partition coefficient matrix when linking with lambda variables
                                    on_coef1.in(indices(unique_ids, partns1, const_idx1));
                                    on_coef2.in(indices(unique_ids, partns2, const_idx2));
                                    
                                    // fill lambda_coef1 and lambda_coef2
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                        auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                        for (int j=0; j<num_partns2+1; ++j) {
                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(1);
                                            lambda_coef1.set_val(cur_idx,1);
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string((num_partns1-2)*2+2);
                                            lambda_coef1.set_val(cur_idx,1);
                                        }
                                        
                                        for (int i=1 ; i<(num_partns1-2)*2+1; i=i+2) {
                                            for (int j=(i-1)/2 + 2; j<num_partns1+1; ++j) {
                                                for(int k=0; k<num_partns2+1; ++k){
                                                    string cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+name2+"{"+to_string(k+1)+"},"+to_string(i+1);
                                                    lambda_coef1.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+name2+"{"+to_string(k+1)+"},"+to_string(i+2);
                                                    lambda_coef1.set_val(cur_idx,-1);
                                                }
                                            }
                                        }
                                        
                                        for (int i=0; i<num_partns1+1; ++i) {
                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(1)+"},"+to_string(1);
                                            lambda_coef2.set_val(cur_idx,1);
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(num_partns2+1)+"},"+to_string((num_partns2-2)*2+2);
                                            lambda_coef2.set_val(cur_idx,1);
                                        }
                                        
                                        for (int i=1 ; i<(num_partns2-2)*2+1; i=i+2) {
                                            for (int j=(i-1)/2 + 2; j<num_partns2+1; ++j) {
                                                for(int k=0; k<num_partns1+1; ++k){
                                                    string cur_idx = cur_var_idx+","+name1+"{"+to_string(k+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                    lambda_coef2.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(k+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                    lambda_coef2.set_val(cur_idx,-1);
                                                }
                                            }
                                        }
                                    }
                                    
                                    
                                    
                                    // fill on_coef1 and on_coef2
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                        auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                        on_coef1.set_val(cur_idx,1);
                                        
                                        for (int i=1; i<num_partns1; ++i) {
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(2);
                                            on_coef1.set_val(cur_idx, 1);
                                        }
                                        
                                        for (int i=2 ; i<(num_partns1-2)*2+2; i=i+2) {
                                            for (int j=i/2+1; j<num_partns1; ++j) {
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                on_coef1.set_val(cur_idx,-1);
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                on_coef1.set_val(cur_idx,1);
                                            }
                                        }
                                        
                                        cur_idx = cur_var_idx+","+name2+"{"+to_string(1)+"},"+to_string(1);
                                        on_coef2.set_val(cur_idx,1);
                                        for (int i=1; i<num_partns2; ++i) {
                                            cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"},"+to_string(2);
                                            on_coef2.set_val(cur_idx, 1);
                                        }
                                        
                                        for (int i=2 ; i<(num_partns2-2)*2+2; i=i+2) {
                                            for (int j=i/2+1; j<num_partns2; ++j) {
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                on_coef2.set_val(cur_idx,-1);
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                on_coef2.set_val(cur_idx,1);
                                            }
                                        }
                                    }
                                }
                                
                                
                                /** Constraints */
                                // Representation of the bilinear term with convex combination
                                Constraint<> bln_rep(pair.first+"_bln_rep");
                                bln_rep = EP.in_matrix(nb_entries,total_entries-nb_entries)*lambda.in_matrix(nb_entries,total_entries-nb_entries) - vlift.in(unique_ids);
                                add(bln_rep.in(unique_ids) == 0);
                                
                                // Representation of o1 with convex combination
                                Constraint<> o1_rep(pair.first+"_o1_rep");
                                o1_rep = bounds1.from_ith(0,inst_partition_lambda).in_matrix(nb_entries, 1) * lambda.in_matrix(nb_entries,total_entries-nb_entries) - o1.in(o1_ids);
                                add(o1_rep.in(unique_ids) == 0);
                                
                                // Representation of o2 with convex combination
                                Constraint<> o2_rep(pair.first+"_o2_rep");
                                o2_rep = bounds2.in_ignore_ith(nb_entries, 1, inst_partition_lambda).in_matrix(nb_entries,1) * lambda.in_matrix(nb_entries,total_entries-nb_entries) - o2.in(o2_ids);
                                add(o2_rep.in(unique_ids) == 0);
                                
                                // Linking partition variables1 with lambda
                                if(model_type == "lambda_II"){
                                    Constraint<> on_link_lambda1(pair.first+"_on_link_lambda1_II");
                                    on_link_lambda1 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef1.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef1.in_matrix(nb_entries,total_entries-nb_entries) - on1.in_ignore_ith(nb_entries_v1,nb_entries_v2,on_coef1.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef1.in_matrix(nb_entries,1);
                                    add(on_link_lambda1.in(indices(unique_ids,const_idx1)) <= 0);
                                    
                                    Constraint<> on_link_lambda2(pair.first+"_on_link_lambda2_II");
                                    on_link_lambda2 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef2.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef2.in_matrix(nb_entries,total_entries-nb_entries) - on1.in_ignore_ith(0,nb_entries_v1,on_coef2.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef2.in_matrix(nb_entries,1);
                                    add(on_link_lambda2.in(indices(unique_ids,const_idx2)) <= 0);
                                }
                                else{
                                    Constraint<> on_link_lambda1(pair.first+"_on_link_lambda1_III");
                                    on_link_lambda1 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef1.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef1.in_matrix(nb_entries,total_entries-nb_entries) - on1.in_ignore_ith(nb_entries_v1,nb_entries_v2,on_coef1.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef1.in_matrix(nb_entries,1);
                                    add(on_link_lambda1.in(indices(unique_ids,const_idx1)) <= 0);
                                    
                                    Constraint<> on_link_lambda2(pair.first+"_on_link_lambda2_III");
                                    on_link_lambda2 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef2.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef2.in_matrix(nb_entries,total_entries-nb_entries) - on1.in_ignore_ith(0,nb_entries_v1,on_coef2.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef2.in_matrix(nb_entries,1);
                                    add(on_link_lambda2.in(indices(unique_ids,const_idx2)) <= 0);
                                }
                                // sum over lambda
                                Constraint<> lambdaSum(pair.first+"_lambdaSum");
                                lambdaSum = sum(lambda.in_matrix(nb_entries,total_entries-nb_entries));
                                add(lambdaSum.in(unique_ids) == 1);
                            }
                        }
                        else{
                            DebugOn("<<<<<<<<<< THIS IS NOT SEEN LIFT -> DOUBLE -> NOT SEEN BINARIES -> DIFF CORE VARS <<<<<<<<<<<" << endl);
                            
                            var<int> on1(name1+"_binary",0,1);
                            indices partns1("partns1");
                            for (int i = 0; i < num_partns1 ; ++i)
                            {
                                partns1.add(name1+ "{" + to_string(i+1) + "}");
                            }
                            add(on1.in(o1_ids_uq,partns1));
                            
                            var<int> on2(name2+"_binary",0,1);
                            indices partns2("partns2");
                            for (int i = 0; i < num_partns2 ; ++i)
                            {
                                partns2.add(name2+ "{" + to_string(i+1) + "}");
                            }
                            add(on2.in(o2_ids_uq,partns2));
                            
                            auto nb_entries_v1 = o1_ids.get_nb_entries();
                            auto nb_entries_v2 = o2_ids.get_nb_entries();
                            auto nb_entries = unique_ids.get_nb_entries();
                            
                            Constraint<> onSum1(o1._name+"_binarySum");
                            onSum1 = sum(on1.in_matrix(nb_entries_v1,1));
                            add(onSum1.in(o1_ids_uq) == 1);
                            
                            Constraint<> onSum2(o2._name+"_binarySum");
                            onSum2 = sum(on2.in_matrix(nb_entries_v2,1));
                            add(onSum2.in(o2_ids_uq) == 1);
                            
                            if(model_type == "on/off"){//if on/off is chosen
                                var<int> on(name1+name2+"_binary",0,1);
                                
                                indices partns("partns");
                                partns = indices(partns1,partns2);
                                auto inst_partition = indices(unique_ids,partns);
                                add(on.in(inst_partition));
                                auto total_entries = inst_partition.get_nb_entries();
                                
                                Constraint<> onLink1(pair.first+"_binaryLink1");
                                onLink1 = on1.from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) - on;
                                add(onLink1.in(inst_partition) >= 0);
                                
                                Constraint<> onLink2(pair.first+"_binaryLink2");
                                onLink2 = on2.in_ignore_ith(nb_entries_v2,1,inst_partition.ignore_ith(0,nb_entries_v1)) - on;
                                add(onLink2.in(inst_partition) >= 0);
                                
                                Constraint<> onLink3(pair.first+"_binaryLink3");
                                onLink3 = on1.from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) + on2.in_ignore_ith(nb_entries_v2,1,inst_partition.ignore_ith(0,nb_entries_v1)) - 1 - on;
                                add(onLink3.in(inst_partition) <= 0);
                                
                                Constraint<> onSumComb(pair.first+"_binarySum");
                                onSumComb = sum(on.in_matrix(nb_entries,total_entries-nb_entries));
                                add(onSumComb.in(unique_ids) == 1);
                                
                                add_on_off_McCormick_refined(pair.first, vlift.in(unique_ids), o1.in(o1_ids), o2.in(o2_ids), on);
                            }
                            
                            else{//means it is one of the lambda formulations
                                
                                //difference is this has one more partition index
                                indices partns1_lambda("partns1_lambda");
                                for (int i = 0; i < num_partns1+1; ++i)
                                {
                                    partns1_lambda.add(name1+ "{" +to_string(i+1) + "}");
                                }
                                
                                indices partns2_lambda("partns2_lambda");
                                for (int i = 0; i < num_partns2+1; ++i)
                                {
                                    partns2_lambda.add(name2+ "{" +to_string(i+1) + "}");
                                }
                                
                                indices partns_lambda("partns_lambda");
                                partns_lambda = indices(partns1_lambda,partns2_lambda);
                                auto inst_partition_lambda = indices(unique_ids,partns_lambda);
                                auto inst_partition_bounds1 = indices(unique_ids,partns1_lambda);
                                auto inst_partition_bounds2 = indices(unique_ids,partns2_lambda);
                                
                                // Convex combination variables
                                var<> lambda(name1+name2+"_lambda",pos_);
                                add(lambda.in(inst_partition_lambda));
                                
                                /** Parameters */
                                // Bounds on variable v1 & v2
                                param<> bounds1(name1+"_bounds1");
                                bounds1.in(inst_partition_bounds1);
                                
                                param<> bounds2(name2+"_bounds2");
                                bounds2.in(inst_partition_bounds2);
                                
                                // Function values on the extreme points
                                param<> EP(name1+name2+"_grid_values");
                                EP.in(inst_partition_lambda);
                                auto total_entries = inst_partition_lambda.get_nb_entries();
                                
                                size_t nb_ins = vlift.in(unique_ids).get_nb_inst();
                                auto o1_global_lb = o1.get_lb();
                                auto increment1 = (o1.get_ub() - o1_global_lb)/num_partns1;
                                
                                auto o2_global_lb = o2.get_lb();
                                auto increment2 = (o2.get_ub() - o2_global_lb)/num_partns2;
                                
                                // fill bounds1 and function values
                                for (int i=0 ; i<num_partns1+1; ++i) {
                                    auto bound_partn1 = o1_global_lb + increment1*i;
                                    bound_partn1.eval_all();
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                        auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"}";
                                        bounds1.set_val(cur_idx,bound_partn1.eval(inst));
                                        for(int j=0; j<num_partns2+1; ++j){
                                            auto bound_partn2 = o2_global_lb + increment2*j;
                                            bound_partn2.eval_all();
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"}";
                                            EP.set_val(cur_idx,(bound_partn1.eval(inst)*bound_partn2.eval(inst)));
                                        }
                                    }
                                }
                                // fill bounds2
                                for (int i=0 ; i<num_partns2+1; ++i) {
                                    auto bound_partn2 = o2_global_lb + increment2*i;
                                    bound_partn2.eval_all();
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                        auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                        string cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"}";
                                        bounds2.set_val(cur_idx,bound_partn2.eval(inst));
                                    }
                                }
                                
                                // Lambda coefficient matrix when linking with partition variables
                                param<> lambda_coef1(name1+"_lambda_linking_coefficients1");
                                param<> lambda_coef2(name2+"_lambda_linking_coefficients2");
                                
                                // Partition coefficient matrix when linking with lambda variables
                                param<> on_coef1(name1+"_partition_linking_coefficients1");
                                param<> on_coef2(name2+"_partition_linking_coefficients2");
                                
                                // create constraint indices
                                indices const_idx1("const_idx1");
                                indices const_idx2("const_idx2");
                                
                                if(model_type == "lambda_II"){
                                    
                                    //fill constraint indices
                                    for (int i=0; i<num_partns1+1; ++i){
                                        const_idx1.add(to_string(i+1));
                                    }
                                    
                                    //fill constraint indices
                                    for (int i=0; i<num_partns2+1; ++i){
                                        const_idx2.add(to_string(i+1));
                                    }
                                    
                                    // Lambda coefficient matrix when linking with partition variables
                                    if(num_partns1 > 1) lambda_coef1.in(indices(inst_partition_lambda, const_idx1));
                                    if(num_partns2 > 1) lambda_coef2.in(indices(inst_partition_lambda, const_idx2));
                                    
                                    // Lambda coefficient matrix when linking with partition variables
                                    if(num_partns1 > 1) on_coef1.in(indices(unique_ids, partns1, const_idx1));
                                    if(num_partns2 > 1) on_coef2.in(indices(unique_ids, partns2, const_idx2));
                                    
                                    // fill lambda_coef1 and lambda_coef2
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                        auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                        for (int i=0 ; i<num_partns1+1; ++i) {
                                            for (int j=0 ; j<num_partns2+1; ++j) {
                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                if(num_partns1 > 1) lambda_coef1.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(j+1);
                                                if(num_partns2 > 1) lambda_coef2.set_val(cur_idx,1);
                                            }
                                        }
                                    }
                                    
                                    // fill on_coef1 and on_coef2
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                        auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                        if(num_partns1 > 1) on_coef1.set_val(cur_idx,1);
                                        cur_idx = cur_var_idx+","+name2+"{"+to_string(1)+"},"+to_string(1);
                                        if(num_partns2 > 1) on_coef2.set_val(cur_idx,1);
                                        if(num_partns1 > 1) {
                                            for (int i=1 ; i<num_partns1; ++i) {
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(i)+"},"+to_string(i+1);
                                                on_coef1.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(i+1);
                                                on_coef1.set_val(cur_idx,1);
                                            }
                                        }
                                        if(num_partns2 > 1) {
                                            for (int i=1 ; i<num_partns2; ++i) {
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(i)+"},"+to_string(i+1);
                                                on_coef2.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"},"+to_string(i+1);
                                                on_coef2.set_val(cur_idx,1);
                                            }
                                        }
                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1)+"},"+to_string(num_partns1+1);
                                        if(num_partns1 > 1) on_coef1.set_val(cur_idx,1);
                                        cur_idx = cur_var_idx+","+name2+"{"+to_string(num_partns2)+"},"+to_string(num_partns2+1);
                                        if(num_partns2 > 1) on_coef2.set_val(cur_idx,1);
                                    }
                                }
                                
                                
                                else /*means model_type == "lambda_III" */{
                                    
                                    //fill constraint indices
                                    for (int i=0; i<(num_partns1-2)*2+2; ++i){
                                        const_idx1.add(to_string(i+1));
                                    }
                                    
                                    //fill constraint indices
                                    for (int i=0; i<(num_partns2-2)*2+2; ++i){
                                        const_idx2.add(to_string(i+1));
                                    }
                                    
                                    // Lambda coefficient matrix when linking with partition variables
                                    if(num_partns1 > 1) lambda_coef1.in(indices(inst_partition_lambda, const_idx1));
                                    if(num_partns2 > 1) lambda_coef2.in(indices(inst_partition_lambda, const_idx2));
                                    
                                    // Partition coefficient matrix when linking with lambda variables
                                    if(num_partns1 > 1) on_coef1.in(indices(unique_ids, partns1, const_idx1));
                                    if(num_partns2 > 1) on_coef2.in(indices(unique_ids, partns2, const_idx2));
                                    
                                    // fill lambda_coef1 and lambda_coef2
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                        auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                        if(num_partns1 > 1) {
                                            for (int j=0; j<num_partns2+1; ++j) {
                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(1);
                                                lambda_coef1.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string((num_partns1-2)*2+2);
                                                lambda_coef1.set_val(cur_idx,1);
                                            }
                                            
                                            for (int i=1 ; i<(num_partns1-2)*2+1; i=i+2) {
                                                for (int j=(i-1)/2 + 2; j<num_partns1+1; ++j) {
                                                    for(int k=0; k<num_partns2+1; ++k){
                                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+name2+"{"+to_string(k+1)+"},"+to_string(i+1);
                                                        lambda_coef1.set_val(cur_idx,1);
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+name2+"{"+to_string(k+1)+"},"+to_string(i+2);
                                                        lambda_coef1.set_val(cur_idx,-1);
                                                    }
                                                }
                                            }
                                        }
                                        if(num_partns2 > 1) {
                                            for (int i=0; i<num_partns1+1; ++i) {
                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(1)+"},"+to_string(1);
                                                lambda_coef2.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(num_partns2+1)+"},"+to_string((num_partns2-2)*2+2);
                                                lambda_coef2.set_val(cur_idx,1);
                                            }
                                            for (int i=1 ; i<(num_partns2-2)*2+1; i=i+2) {
                                                for (int j=(i-1)/2 + 2; j<num_partns2+1; ++j) {
                                                    for(int k=0; k<num_partns1+1; ++k){
                                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(k+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                        lambda_coef2.set_val(cur_idx,1);
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(k+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                        lambda_coef2.set_val(cur_idx,-1);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    
                                    
                                    // fill on_coef1 and on_coef2
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift.in(unique_ids).get_id_inst(inst);
                                        auto cur_var_idx = unique_ids._keys->at(cur_var_id);
                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                        if(num_partns1 > 1) {
                                            on_coef1.set_val(cur_idx,1);
                                            
                                            for (int i=1; i<num_partns1; ++i) {
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(2);
                                                on_coef1.set_val(cur_idx, 1);
                                            }
                                            
                                            for (int i=2 ; i<(num_partns1-2)*2+2; i=i+2) {
                                                for (int j=i/2+1; j<num_partns1; ++j) {
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                    on_coef1.set_val(cur_idx,-1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                    on_coef1.set_val(cur_idx,1);
                                                }
                                            }
                                        }
                                        if(num_partns2 > 1) {
                                            cur_idx = cur_var_idx+","+name2+"{"+to_string(1)+"},"+to_string(1);
                                            on_coef2.set_val(cur_idx,1);
                                            for (int i=1; i<num_partns2; ++i) {
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"},"+to_string(2);
                                                on_coef2.set_val(cur_idx, 1);
                                            }
                                            
                                            for (int i=2 ; i<(num_partns2-2)*2+2; i=i+2) {
                                                for (int j=i/2+1; j<num_partns2; ++j) {
                                                    cur_idx = cur_var_idx+","+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                    on_coef2.set_val(cur_idx,-1);
                                                    cur_idx = cur_var_idx+","+name2+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                    on_coef2.set_val(cur_idx,1);
                                                }
                                            }
                                        }
                                    }
                                }
                                
                                
                                /** Constraints */
                                // Representation of the bilinear term with convex combination
                                Constraint<> bln_rep(pair.first+"_bln_rep");
                                bln_rep = EP.in_matrix(nb_entries,total_entries-nb_entries)*lambda.in_matrix(nb_entries,total_entries-nb_entries) - vlift.in(unique_ids);
                                add(bln_rep.in(unique_ids) == 0);
                                
                                // Representation of o1 with convex combination
                                Constraint<> o1_rep(pair.first+"_o1_rep");
                                o1_rep = bounds1.from_ith(0,inst_partition_lambda).in_matrix(nb_entries, 1) * lambda.in_matrix(nb_entries,total_entries-nb_entries) - o1.in(o1_ids);
                                add(o1_rep.in(unique_ids) == 0);
                                
                                // Representation of o2 with convex combination
                                Constraint<> o2_rep(pair.first+"_o2_rep");
                                o2_rep = bounds2.in_ignore_ith(nb_entries, 1, inst_partition_lambda).in_matrix(nb_entries,1) * lambda.in_matrix(nb_entries,total_entries-nb_entries) - o2.in(o2_ids);
                                add(o2_rep.in(unique_ids) == 0);
                                
                                // Linking partition variables1 with lambda
                                if(model_type == "lambda_II"){
                                    if(num_partns1 > 1) {
                                        Constraint<> on_link_lambda1(pair.first+"_on_link_lambda1_II");
                                        on_link_lambda1 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef1.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef1.in_matrix(nb_entries,total_entries-nb_entries) - on1.in_ignore_ith(nb_entries_v1,nb_entries_v2,on_coef1.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef1.in_matrix(nb_entries,1);
                                        add(on_link_lambda1.in(indices(unique_ids,const_idx1)) <= 0);
                                    }
                                    if(num_partns2 > 1) {
                                        Constraint<> on_link_lambda2(pair.first+"_on_link_lambda2_II");
                                        on_link_lambda2 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef2.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef2.in_matrix(nb_entries,total_entries-nb_entries) - on2.in_ignore_ith(0,nb_entries_v1,on_coef2.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef2.in_matrix(nb_entries,1);
                                        add(on_link_lambda2.in(indices(unique_ids,const_idx2)) <= 0);
                                    }
                                }
                                else{
                                    if(num_partns1 > 1) {
                                        Constraint<> on_link_lambda1(pair.first+"_on_link_lambda1_III");
                                        on_link_lambda1 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef1.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef1.in_matrix(nb_entries,total_entries-nb_entries) - on1.in_ignore_ith(nb_entries_v1,nb_entries_v2,on_coef1.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef1.in_matrix(nb_entries,1);
                                        add(on_link_lambda1.in(indices(unique_ids,const_idx1)) <= 0);
                                    }
                                    if(num_partns2 > 1) {
                                        Constraint<> on_link_lambda2(pair.first+"_on_link_lambda2_III");
                                        on_link_lambda2 = lambda.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef2.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef2.in_matrix(nb_entries,total_entries-nb_entries) - on2.in_ignore_ith(0,nb_entries_v1,on_coef2.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef2.in_matrix(nb_entries,1);
                                        add(on_link_lambda2.in(indices(unique_ids,const_idx2)) <= 0);
                                    }
                                }
                                // sum over lambda
                                Constraint<> lambdaSum(pair.first+"_lambdaSum");
                                lambdaSum = sum(lambda.in_matrix(nb_entries,total_entries-nb_entries));
                                add(lambdaSum.in(unique_ids) == 1);
                            }
                        }
                    }
                }
#endif
            }
            else {
                add_McCormick(pair.first, vlift.in(unique_ids), o1.in(o1_ids), o2.in(o2_ids));
            }
        }
        else {
            auto vlift = static_pointer_cast<var<type>>(it->second);
            auto added = vlift->add_bounds(lb,ub);
            lt._p = make_shared<var<type>>(vlift->in(ids));
            if(!added.empty()){
                assert(o1._indices->size()==o2._indices->size());
                if(added.size()!=o1._indices->size()){/* If some keys are repeated, remove them from the refs of o1 and o2 */
                    auto keep_refs = ids.diff_refs(added);
                    o1_ids.filter_refs(keep_refs);
                    o2_ids.filter_refs(keep_refs);
                }
                reindex_vars();
                // If some keys are repeated in individual indices, remove them from the refs of o1 and o2
                auto o1_ids_uq = o1_ids;
                auto o2_ids_uq = o2_ids;
                auto keep_refs1 = o1_ids_uq.get_unique_refs();
                auto keep_refs2 = o2_ids_uq.get_unique_refs();
                o1_ids_uq.filter_refs(keep_refs1);
                o2_ids_uq.filter_refs(keep_refs2);
                reindex_vars();
                
                //check the sign of the lift and the correspoinding boudning functions
                if(c.check_soc() && c.is_eq()){
                    if(lift_sign){
                        vlift->_lift_ub = true;
                        vlift->_lift_lb = false;
                    }
                    else{
                        vlift->_lift_ub = false;
                        vlift->_lift_lb = true;
                    }
                }
                else{
                    vlift->_lift_ub = true;
                    vlift->_lift_lb = true;
                }
                if((num_partns1 > 1) || (num_partns2 > 1)) {
#ifdef PARTITION
                    auto binvar_ptr1 = _vars_name.find(name1+"_binary");
                    auto binvar_ptr2 = _vars_name.find(name2+"_binary");
                    auto binvar_ptr3 = _vars_name.find(name1+name2+"_binary");
                    
                    if (o1 == o2) //if the variables are same add 1d partition
                    {
                        if(binvar_ptr1 !=_vars_name.end()){ //means v1 has been partitioned before
                            DebugOn("<<<<<<<<<< THIS IS SEEN LIFT -> SINGLE -> SEEN BINARY <<<<<<<<<<<" << endl);
                            
                            auto binvar1 = static_pointer_cast<var<int>>(binvar_ptr1->second);
                            
                            indices partns("partns");
                            for (int i = 0; i < num_partns1 ; ++i)
                            {
                                partns.add(name1+  "{" + to_string(i+1) + "}");
                            }
                            auto inst_partition = indices(added,partns);
                            
                            param<int> lb1("lb1"), ub1("ub1");
                            lb1.in(added,partns);
                            ub1.in(added,partns);
                            lb1.set_val(0), ub1.set_val(1);
                            auto added1 = binvar1->add_bounds(lb1,ub1);
                            reindex_vars();
                            
                            auto nb_entries_v1 = o1_ids.get_nb_entries();
                            auto nb_entries = added.get_nb_entries();
                            auto total_entries = inst_partition.get_nb_entries();
                            
                            Constraint<> onSumComb(pair.first+"_binarySum");
                            onSumComb = sum((binvar1->in(added1)).in_matrix(nb_entries,total_entries-nb_entries));
                            add(onSumComb.in(added) == 1);
                            
                            if(model_type == "on/off"){//if on/off is chosen
                                add_on_off_McCormick_refined(pair.first, vlift->in(added), o1.in(o1_ids), o2.in(o2_ids), binvar1->in(added1));
                            }
                            
                            else{ //means it is one of the lambda formulations
                                
                                //difference is this has one more partition index
                                indices partns_lambda("partns_lambda");
                                for (int i = 0; i < num_partns1+1 ; ++i)
                                {
                                    partns_lambda.add(name1+ "{" +to_string(i+1) + "}");
                                }
                                auto inst_partition_lambda = indices(added,partns_lambda);
                                
                                // Convex combination variables
                                auto lambda_ptr = _vars_name.find(name1+"_lambda");
                                auto lambda = static_pointer_cast<var<double>>(lambda_ptr->second);
                                param<double> lb_lambda("lb_lambda"), ub_lambda("ub_lambda");
                                lb_lambda.in(added,partns_lambda);
                                ub_lambda.in(added,partns_lambda);
                                lb_lambda.set_val(0), ub_lambda.set_val(1);
                                auto added_lambda = lambda->add_bounds(lb_lambda,ub_lambda);
                                reindex_vars();
                                
                                /** Parameters */
                                // Bounds on variable v1 & v2
                                param<> bounds(name1+"_bounds");
                                bounds.in(inst_partition_lambda);
                                
                                // Function values on the extreme points
                                param<> EP(name1+name2+"_grid_values");
                                EP.in(inst_partition_lambda);
                                
                                size_t nb_ins = vlift->in(added).get_nb_inst();
                                auto o1_global_lb = o1.get_lb();
                                auto increment = (o1.get_ub() - o1_global_lb)/num_partns1;
                                
                                // fill bounds and function values
                                for (int i=0 ; i<num_partns1+1; ++i) {
                                    auto bound_partn = o1_global_lb + increment*i;
                                    bound_partn.eval_all();
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift->get_id_inst(inst);
                                        auto cur_var_idx = added._keys->at(cur_var_id);
                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"}";
                                        bounds.set_val(cur_idx,bound_partn.eval(inst));
                                        EP.set_val(cur_idx,(bound_partn.eval(inst)*bound_partn.eval(inst)));
                                    }
                                }
                                
                                // Lambda coefficient matrix when linking with partition variables
                                param<> lambda_coef(name1+"_lambda_linking_coefficients");
                                
                                // Partition coefficient matrix when linking with lambda variables
                                param<> on_coef(name1+"_partition_linking_coefficients");
                                
                                // create constraint indices
                                indices const_idx("const_idx");
                                
                                if(model_type == "lambda_II"){
                                    
                                    //fill constraint indices
                                    for (int i=0; i<num_partns1+1; ++i){
                                        const_idx.add(to_string(i+1));
                                    }
                                    
                                    // Lambda coefficient matrix when linking with partition variables
                                    lambda_coef.in(indices(inst_partition_lambda, const_idx));
                                    
                                    // Partition coefficient matrix when linking with lambda variables
                                    on_coef.in(indices(inst_partition, const_idx));
                                    
                                    // fill lambda_coef
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift->get_id_inst(inst);
                                        auto cur_var_idx = added._keys->at(cur_var_id);
                                        for (int i=0 ; i<num_partns1+1; ++i) {
                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(i+1);
                                            lambda_coef.set_val(cur_idx,1);
                                        }
                                    }
                                    
                                    // fill on_coef
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift->get_id_inst(inst);
                                        auto cur_var_idx = added._keys->at(cur_var_id);
                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                        on_coef.set_val(cur_idx,1);
                                        for (int i=1 ; i<num_partns1; ++i) {
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(i)+"},"+to_string(i+1);
                                            on_coef.set_val(cur_idx,1);
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(i+1);
                                            on_coef.set_val(cur_idx,1);
                                        }
                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1)+"},"+to_string(num_partns1+1);
                                        on_coef.set_val(cur_idx,1);
                                    }
                                }
                                
                                else /*means model_type == "lambda_III" */{
                                    
                                    //fill constraint indices
                                    for (int i=0; i<(num_partns1-2)*2+2; ++i){
                                        const_idx.add(to_string(i+1));
                                    }
                                    
                                    // Lambda coefficient matrix when linking with partition variables
                                    lambda_coef.in(indices(inst_partition_lambda, const_idx));
                                    
                                    // Partition coefficient matrix when linking with lambda variables
                                    on_coef.in(indices(inst_partition, const_idx));
                                    
                                    // fill lambda_coef
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift->get_id_inst(inst);
                                        auto cur_var_idx = added._keys->at(cur_var_id);
                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                        lambda_coef.set_val(cur_idx,1);
                                        for (int i=1 ; i<(num_partns1-2)*2+1; i=i+2) {
                                            for (int j=(i-1)/2 + 2; j<num_partns1+1; ++j) {
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                lambda_coef.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                lambda_coef.set_val(cur_idx,-1);
                                            }
                                        }
                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1+1)+"},"+to_string((num_partns1-2)*2+2);
                                        lambda_coef.set_val(cur_idx,1);
                                    }
                                    
                                    // fill on_coef
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift->get_id_inst(inst);
                                        auto cur_var_idx = added._keys->at(cur_var_id);
                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                        on_coef.set_val(cur_idx,1);
                                        
                                        for (int i=1; i<num_partns1; ++i) {
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(2);
                                            on_coef.set_val(cur_idx, 1);
                                        }
                                        
                                        for (int i=2 ; i<(num_partns1-2)*2+2; i=i+2) {
                                            for (int j=i/2+1; j<num_partns1; ++j) {
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                on_coef.set_val(cur_idx,-1);
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                on_coef.set_val(cur_idx,1);
                                            }
                                        }
                                    }
                                    
                                }
                                
                                
                                /** Constraints */
                                
                                if (vlift->_lift_ub){
                                    // Representation of the quadratic term with secant
                                    Constraint<> quad_ub(pair.first+"_quad_ub");
                                    quad_ub = EP.in_matrix(nb_entries,total_entries-nb_entries)*lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries) - vlift->in(added);
                                    add(quad_ub.in(added) >= 0); /*using it as the upper bound to be valid*/
                                }
                                
                                if (vlift->_lift_lb){
                                    Constraint<> quad_lb(pair.first+"_quad_lb");
                                    quad_lb = o1.from_ith(0,added)*o2.from_ith(nb_entries_v1,added) - vlift->in(added);
                                    quad_lb._relaxed = true;
                                    add(quad_lb.in(added) <= 0); /*using it as the lower bound to be valid*/
                                }
                                
                                // Representation of o1 with convex combination
                                Constraint<> o1_rep(pair.first+"_o1_rep");
                                o1_rep = bounds.from_ith(0,inst_partition_lambda).in_matrix(nb_entries, 1) * lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries) - o1.in(o1_ids);
                                add(o1_rep.in(added) == 0);
                                
                                // Linking partition variables with lambda
                                if(model_type == "lambda_II"){
                                    Constraint<> on_link_lambda(pair.first+"_on_link_lambda_II");
                                    on_link_lambda = lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef.in_matrix(nb_entries,total_entries-nb_entries) - binvar1->in(added1).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,on_coef.get_matrix_ids(nb_entries,total_entries-nb_entries)) * on_coef.in_matrix(nb_entries,total_entries-nb_entries);
                                    add(on_link_lambda.in(indices(added,const_idx)) <= 0);
                                }
                                else{
                                    Constraint<> on_link_lambda(pair.first+"_on_link_lambda_III");
                                    on_link_lambda = lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef.in_matrix(nb_entries,total_entries-nb_entries) - binvar1->in(added1).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,on_coef.get_matrix_ids(nb_entries,total_entries-nb_entries)) * on_coef.in_matrix(nb_entries,total_entries-nb_entries);
                                    add(on_link_lambda.in(indices(added,const_idx)) <= 0);
                                }
                                
                                
                                // sum over lambda
                                Constraint<> lambdaSum(pair.first+"_lambdaSum");
                                lambdaSum = sum(lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries));
                                add(lambdaSum.in(added) == 1);
                            }
                        }
                        else{ //means v1 has not been partitioned before
                            DebugOn("<<<<<<<<<< THIS IS SEEN LIFT -> DOUBLE -> NOT SEEN BINARY <<<<<<<<<<<" << endl);
                            
                            var<int> on(name1+"_binary",0,1);
                            indices partns("partns");
                            for (int i = 0; i < num_partns1 ; ++i)
                            {
                                partns.add(name1+  "{" + to_string(i+1) + "}");
                            }
                            auto inst_partition = indices(added,partns);
                            add(on.in(inst_partition));
                            
                            auto nb_entries_v1 = o1_ids.get_nb_entries();
                            auto nb_entries = added.get_nb_entries();
                            auto total_entries = inst_partition.get_nb_entries();
                            
                            Constraint<> onSumComb(pair.first+"_binarySum");
                            onSumComb = sum(on.in_matrix(nb_entries,total_entries-nb_entries));
                            add(onSumComb.in(added) == 1);
                            
                            if(model_type == "on/off"){//if on/off is chosen
                                add_on_off_McCormick_refined(pair.first, vlift->in(added), o1.in(o1_ids), o2.in(o2_ids), on);
                            }
                            
                            else{ //means it is one of the lambda formulations
                                
                                //difference is this has one more partition index
                                indices partns_lambda("partns_lambda");
                                for (int i = 0; i < num_partns1+1 ; ++i)
                                {
                                    partns_lambda.add(name1+ "{" +to_string(i+1) + "}");
                                }
                                auto inst_partition_lambda = indices(added,partns_lambda);
                                
                                // Convex combination variables
                                auto lambda_ptr = _vars_name.find(name1+"_lambda");
                                auto lambda = static_pointer_cast<var<double>>(lambda_ptr->second);
                                param<double> lb_lambda("lb_lambda"), ub_lambda("ub_lambda");
                                lb_lambda.in(added,partns_lambda);
                                ub_lambda.in(added,partns_lambda);
                                lb_lambda.set_val(0), ub_lambda.set_val(1);
                                auto added_lambda = lambda->add_bounds(lb_lambda,ub_lambda);
                                reindex_vars();
                                
                                /** Parameters */
                                // Bounds on variable v1 & v2
                                param<> bounds(name1+"_bounds");
                                bounds.in(inst_partition_lambda);
                                
                                // Function values on the extreme points
                                param<> EP(name1+name2+"_grid_values");
                                EP.in(inst_partition_lambda);
                                
                                size_t nb_ins = vlift->in(added).get_nb_inst();
                                auto o1_global_lb = o1.get_lb();
                                auto increment = (o1.get_ub() - o1_global_lb)/num_partns1;
                                
                                // fill bounds and function values
                                for (int i=0 ; i<num_partns1+1; ++i) {
                                    auto bound_partn = o1_global_lb + increment*i;
                                    bound_partn.eval_all();
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift->get_id_inst(inst);
                                        auto cur_var_idx = added._keys->at(cur_var_id);
                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"}";
                                        bounds.set_val(cur_idx,bound_partn.eval(inst));
                                        EP.set_val(cur_idx,(bound_partn.eval(inst)*bound_partn.eval(inst)));
                                    }
                                }
                                
                                // Lambda coefficient matrix when linking with partition variables
                                param<> lambda_coef(name1+"_lambda_linking_coefficients");
                                
                                // Partition coefficient matrix when linking with lambda variables
                                param<> on_coef(name1+"_partition_linking_coefficients");
                                
                                // create constraint indices
                                indices const_idx("const_idx");
                                
                                if(model_type == "lambda_II"){
                                    
                                    //fill constraint indices
                                    for (int i=0; i<num_partns1+1; ++i){
                                        const_idx.add(to_string(i+1));
                                    }
                                    
                                    // Lambda coefficient matrix when linking with partition variables
                                    lambda_coef.in(indices(inst_partition_lambda, const_idx));
                                    
                                    // Partition coefficient matrix when linking with lambda variables
                                    on_coef.in(indices(inst_partition, const_idx));
                                    
                                    // fill lambda_coef
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift->get_id_inst(inst);
                                        auto cur_var_idx = added._keys->at(cur_var_id);
                                        for (int i=0 ; i<num_partns1+1; ++i) {
                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(i+1);
                                            lambda_coef.set_val(cur_idx,1);
                                        }
                                    }
                                    
                                    // fill on_coef
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift->get_id_inst(inst);
                                        auto cur_var_idx = added._keys->at(cur_var_id);
                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                        on_coef.set_val(cur_idx,1);
                                        for (int i=1 ; i<num_partns1; ++i) {
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(i)+"},"+to_string(i+1);
                                            on_coef.set_val(cur_idx,1);
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(i+1);
                                            on_coef.set_val(cur_idx,1);
                                        }
                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1)+"},"+to_string(num_partns1+1);
                                        on_coef.set_val(cur_idx,1);
                                    }
                                }
                                
                                else /*means model_type == "lambda_III" */{
                                    
                                    //fill constraint indices
                                    for (int i=0; i<(num_partns1-2)*2+2; ++i){
                                        const_idx.add(to_string(i+1));
                                    }
                                    
                                    // Lambda coefficient matrix when linking with partition variables
                                    lambda_coef.in(indices(inst_partition_lambda, const_idx));
                                    
                                    // Partition coefficient matrix when linking with lambda variables
                                    on_coef.in(indices(inst_partition, const_idx));
                                    
                                    // fill lambda_coef
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift->get_id_inst(inst);
                                        auto cur_var_idx = added._keys->at(cur_var_id);
                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                        lambda_coef.set_val(cur_idx,1);
                                        for (int i=1 ; i<(num_partns1-2)*2+1; i=i+2) {
                                            for (int j=(i-1)/2 + 2; j<num_partns1+1; ++j) {
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                lambda_coef.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                lambda_coef.set_val(cur_idx,-1);
                                            }
                                        }
                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1+1)+"},"+to_string((num_partns1-2)*2+2);
                                        lambda_coef.set_val(cur_idx,1);
                                    }
                                    
                                    // fill on_coef
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift->get_id_inst(inst);
                                        auto cur_var_idx = added._keys->at(cur_var_id);
                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                        on_coef.set_val(cur_idx,1);
                                        
                                        for (int i=1; i<num_partns1; ++i) {
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(2);
                                            on_coef.set_val(cur_idx, 1);
                                        }
                                        
                                        for (int i=2 ; i<(num_partns1-2)*2+2; i=i+2) {
                                            for (int j=i/2+1; j<num_partns1; ++j) {
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                on_coef.set_val(cur_idx,-1);
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                on_coef.set_val(cur_idx,1);
                                            }
                                        }
                                    }
                                    
                                }
                                
                                
                                /** Constraints */
                                
                                if (vlift->_lift_ub){
                                    // Representation of the quadratic term with secant
                                    Constraint<> quad_ub(pair.first+"_quad_ub");
                                    quad_ub = EP.in_matrix(nb_entries,total_entries-nb_entries)*lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries) - vlift->in(added);
                                    add(quad_ub.in(added) >= 0); /*using it as the upper bound to be valid*/
                                }
                                
                                if (vlift->_lift_lb){
                                    Constraint<> quad_lb(pair.first+"_quad_lb");
                                    quad_lb = o1.from_ith(0,added)*o2.from_ith(nb_entries_v1,added) - vlift->in(added);
                                    quad_lb._relaxed = true;
                                    add(quad_lb.in(added) <= 0); /*using it as the lower bound to be valid*/
                                }
                                
                                // Representation of o1 with convex combination
                                Constraint<> o1_rep(pair.first+"_o1_rep");
                                o1_rep = bounds.from_ith(0,inst_partition_lambda).in_matrix(nb_entries, 1) * lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries) - o1.in(o1_ids);
                                add(o1_rep.in(added) == 0);
                                
                                // Linking partition variables with lambda
                                if(model_type == "lambda_II"){
                                    Constraint<> on_link_lambda(pair.first+"_on_link_lambda_II");
                                    on_link_lambda = lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef.in_matrix(nb_entries,total_entries-nb_entries) - on.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,on_coef.get_matrix_ids(nb_entries,total_entries-nb_entries)) * on_coef.in_matrix(nb_entries,total_entries-nb_entries);
                                    add(on_link_lambda.in(indices(added,const_idx)) <= 0);
                                }
                                else{
                                    Constraint<> on_link_lambda(pair.first+"_on_link_lambda_III");
                                    on_link_lambda = lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef.in_matrix(nb_entries,total_entries-nb_entries) - on.in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,on_coef.get_matrix_ids(nb_entries,total_entries-nb_entries)) * on_coef.in_matrix(nb_entries,total_entries-nb_entries);
                                    add(on_link_lambda.in(indices(added,const_idx)) <= 0);
                                }
                                
                                
                                // sum over lambda
                                Constraint<> lambdaSum(pair.first+"_lambdaSum");
                                lambdaSum = sum(lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries));
                                add(lambdaSum.in(added) == 1);
                            }
                        }
                    }
                    else{ //else add 2d partition
                        if(binvar_ptr1 !=_vars_name.end() && binvar_ptr2 !=_vars_name.end()){ //means v1 and v2 has been partitioned before
                            if(name1 == name2){
                                DebugOn("<<<<<<<<<< THIS IS SEEN LIFT -> DOUBLE -> SEEN BOTH BINARIES -> SAME CORE VARS <<<<<<<<<<<" << endl);
                                
                                indices partns1("partns1");
                                for (int i = 0; i < num_partns1 ; ++i)
                                {
                                    partns1.add(name1+ "{" + to_string(i+1) + "}");
                                }
                                
                                auto binvar1 = static_pointer_cast<var<int>>(binvar_ptr1->second);
                                param<int> lb1("lb1"), ub1("ub1");
                                lb1.in(union_ids(o1_ids_uq,o2_ids_uq),partns1);
                                ub1.in(union_ids(o1_ids_uq,o2_ids_uq),partns1);
                                lb1.set_val(0), ub1.set_val(1);
                                auto added1 = binvar1->add_bounds(lb1,ub1);
                                reindex_vars();
                                
                                indices partns2("partns2");
                                for (int i = 0; i < num_partns2 ; ++i)
                                {
                                    partns2.add(name2+ "{" + to_string(i+1) + "}");
                                }
                                
                                auto nb_entries_v1 = o1_ids.get_nb_entries();
                                auto nb_entries_v2 = o2_ids.get_nb_entries();
                                auto nb_entries = added.get_nb_entries();
                                
                                if(!added1.empty()){
                                    Constraint<> onSum1(o1._name+"_binarySum");
                                    onSum1 = sum(binvar1->in(added1).in_matrix(nb_entries_v1,1));
                                    auto vset1 = added1.from_ith(0,nb_entries_v1);
                                    vset1.filter_refs(vset1.get_unique_refs());
                                    add(onSum1.in(vset1) == 1);
                                }
                                
                                if(model_type == "on/off"){//if on/off is chosen
                                    
                                    indices partns("partns");
                                    partns = indices(partns1,partns2);
                                    auto inst_partition = indices(added,partns);
                                    auto total_entries = inst_partition.get_nb_entries();
                                    
                                    if(binvar_ptr3 !=_vars_name.end()){ //means the combined binary variable has been used before
                                        auto binvar3 = static_pointer_cast<var<int>>(binvar_ptr3->second);
                                        param<int> lb3("lb3"), ub3("ub3");
                                        lb3.in(added,partns);
                                        ub3.in(added,partns);
                                        lb3.set_val(0), ub3.set_val(1);
                                        auto added3 = binvar3->add_bounds(lb3,ub3);
                                        reindex_vars();
                                        
                                        Constraint<> onLink1(pair.first+"_binaryLink1");
                                        onLink1 = binvar1->from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v1)) - binvar3->in(inst_partition);
                                        add(onLink1.in(inst_partition) >= 0);
                                        
                                        Constraint<> onLink2(pair.first+"_binaryLink2");
                                        onLink2 = binvar1->in_ignore_ith(nb_entries_v1,1,inst_partition.ignore_ith(0,nb_entries_v1)) - binvar3->in(inst_partition);
                                        add(onLink2.in(inst_partition) >= 0);
                                        
                                        Constraint<> onLink3(pair.first+"_binaryLink3");
                                        onLink3 = binvar1->from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v1)) + binvar1->in_ignore_ith(nb_entries_v1,1,inst_partition.ignore_ith(0,nb_entries_v1)) - 1 - binvar3->in(inst_partition);
                                        add(onLink3.in(inst_partition) <= 0);
                                        
                                        Constraint<> onSumComb(pair.first+"_binarySum");
                                        onSumComb = sum((binvar3->in(added3)).in_matrix(nb_entries,total_entries-nb_entries));
                                        add(onSumComb.in(added) == 1);
                                        
                                        add_on_off_McCormick_refined(pair.first, vlift->in(added), o1.in(o1_ids), o2.in(o2_ids), binvar3->in(added3));
                                    }
                                    else{ // means the combined binary variable has not been used before
                                        
                                        var<int> on(name1+name2+"_binary",0,1);
                                        add(on.in(inst_partition));
                                        
                                        Constraint<> onLink1(pair.first+"_binaryLink1");
                                        onLink1 = binvar1->from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v1)) - on;
                                        add(onLink1.in(inst_partition) >= 0);
                                        
                                        Constraint<> onLink2(pair.first+"_binaryLink2");
                                        onLink2 = binvar1->in_ignore_ith(nb_entries_v1,1,inst_partition.ignore_ith(0,nb_entries_v1)) - on;
                                        add(onLink2.in(inst_partition) >= 0);
                                        
                                        Constraint<> onLink3(pair.first+"_binaryLink3");
                                        onLink3 = binvar1->from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v1)) + binvar1->in_ignore_ith(nb_entries_v1,1,inst_partition.ignore_ith(0,nb_entries_v1)) - 1 - on;
                                        add(onLink3.in(inst_partition) <= 0);
                                        
                                        Constraint<> onSumComb(pair.first+"_binarySum");
                                        onSumComb = sum(on.in_matrix(nb_entries,total_entries-nb_entries));
                                        add(onSumComb.in(added) == 1);
                                        
                                        add_on_off_McCormick_refined(pair.first, vlift->in(added), o1.in(o1_ids), o2.in(o2_ids), on);
                                    }
                                }
                                
                                else{//means it is one of the lambda formulations
                                    
                                    //difference is this has one more partition index
                                    indices partns1_lambda("partns1_lambda");
                                    for (int i = 0; i < num_partns1+1; ++i)
                                    {
                                        partns1_lambda.add(name1+ "{" +to_string(i+1) + "}");
                                    }
                                    
                                    indices partns2_lambda("partns2_lambda");
                                    for (int i = 0; i < num_partns2+1; ++i)
                                    {
                                        partns2_lambda.add(name2+ "{" +to_string(i+1) + "}");
                                    }
                                    
                                    indices partns_lambda("partns_lambda");
                                    partns_lambda = indices(partns1_lambda,partns2_lambda);
                                    auto inst_partition_lambda = indices(added,partns_lambda);
                                    auto inst_partition_bounds1 = indices(added,partns1_lambda);
                                    auto inst_partition_bounds2 = indices(added,partns2_lambda);
                                    
                                    // Convex combination variables
                                    auto lambda_ptr = _vars_name.find(name1+name2+"_lambda");
                                    auto lambda = static_pointer_cast<var<double>>(lambda_ptr->second);
                                    param<double> lb_lambda("lb_lambda"), ub_lambda("ub_lambda");
                                    lb_lambda.in(added,partns_lambda);
                                    ub_lambda.in(added,partns_lambda);
                                    lb_lambda.set_val(0), ub_lambda.set_val(1);
                                    auto added_lambda = lambda->add_bounds(lb_lambda,ub_lambda);
                                    reindex_vars();
                                    
                                    /** Parameters */
                                    // Bounds on variable v1 & v2
                                    param<> bounds1(name1+"_bounds1");
                                    bounds1.in(inst_partition_bounds1);
                                    
                                    param<> bounds2(name2+"_bounds2");
                                    bounds2.in(inst_partition_bounds2);
                                    
                                    // Function values on the extreme points
                                    param<> EP(name1+name2+"_grid_values");
                                    EP.in(inst_partition_lambda);
                                    auto total_entries = inst_partition_lambda.get_nb_entries();
                                    
                                    size_t nb_ins = vlift->in(added).get_nb_inst();
                                    auto o1_global_lb = o1.get_lb();
                                    auto increment1 = (o1.get_ub() - o1_global_lb)/num_partns1;
                                    
                                    auto o2_global_lb = o2.get_lb();
                                    auto increment2 = (o2.get_ub() - o2_global_lb)/num_partns2;
                                    
                                    // fill bounds and function values
                                    for (int i=0 ; i<num_partns1+1; ++i) {
                                        auto bound_partn1 = o1_global_lb + increment1*i;
                                        bound_partn1.eval_all();
                                        auto bound_partn2 = o2_global_lb + increment2*i;
                                        bound_partn2.eval_all();
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift->get_id_inst(inst);
                                            auto cur_var_idx = added._keys->at(cur_var_id);
                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"}";
                                            bounds1.set_val(cur_idx,bound_partn1.eval(inst));
                                            bounds2.set_val(cur_idx,bound_partn2.eval(inst));
                                            for(int j=0; j<num_partns2+1; ++j){
                                                auto bound_partn2_temp = o2_global_lb + increment2*j;
                                                bound_partn2_temp.eval_all();
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"}";
                                                EP.set_val(cur_idx,(bound_partn1.eval(inst)*bound_partn2_temp.eval(inst)));
                                            }
                                        }
                                    }
                                    
                                    // Lambda coefficient matrix when linking with partition variables
                                    param<> lambda_coef1(name1+"_lambda_linking_coefficients1");
                                    param<> lambda_coef2(name2+"_lambda_linking_coefficients2");
                                    
                                    // Partition coefficient matrix when linking with lambda variables
                                    param<> on_coef1(name1+"_partition_linking_coefficients1");
                                    param<> on_coef2(name2+"_partition_linking_coefficients2");
                                    
                                    // create constraint indices
                                    indices const_idx1("const_idx1");
                                    indices const_idx2("const_idx2");
                                    
                                    if(model_type == "lambda_II"){
                                        
                                        //fill constraint indices
                                        for (int i=0; i<num_partns1+1; ++i){
                                            const_idx1.add(to_string(i+1));
                                        }
                                        
                                        //fill constraint indices
                                        for (int i=0; i<num_partns2+1; ++i){
                                            const_idx2.add(to_string(i+1));
                                        }
                                        
                                        // Lambda coefficient matrix when linking with partition variables
                                        lambda_coef1.in(indices(inst_partition_lambda, const_idx1));
                                        lambda_coef2.in(indices(inst_partition_lambda, const_idx2));
                                        
                                        // Lambda coefficient matrix when linking with partition variables
                                        on_coef1.in(indices(added, partns1, const_idx1));
                                        on_coef2.in(indices(added, partns2, const_idx2));
                                        
                                        // fill lambda_coef1 and lambda_coef2
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift->get_id_inst(inst);
                                            auto cur_var_idx = added._keys->at(cur_var_id);
                                            for (int i=0 ; i<num_partns1+1; ++i) {
                                                for (int j=0 ; j<num_partns2+1; ++j) {
                                                    string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                    lambda_coef1.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(j+1);
                                                    lambda_coef2.set_val(cur_idx,1);
                                                }
                                            }
                                        }
                                        
                                        // fill on_coef1 and on_coef2
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift->get_id_inst(inst);
                                            auto cur_var_idx = added._keys->at(cur_var_id);
                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                            on_coef1.set_val(cur_idx,1);
                                            cur_idx = cur_var_idx+","+name2+"{"+to_string(1)+"},"+to_string(1);
                                            on_coef2.set_val(cur_idx,1);
                                            for (int i=1 ; i<num_partns1; ++i) {
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(i)+"},"+to_string(i+1);
                                                on_coef1.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(i+1);
                                                on_coef1.set_val(cur_idx,1);
                                            }
                                            for (int i=1 ; i<num_partns2; ++i) {
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(i)+"},"+to_string(i+1);
                                                on_coef2.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"},"+to_string(i+1);
                                                on_coef2.set_val(cur_idx,1);
                                            }
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1)+"},"+to_string(num_partns1+1);
                                            on_coef1.set_val(cur_idx,1);
                                            cur_idx = cur_var_idx+","+name2+"{"+to_string(num_partns2)+"},"+to_string(num_partns2+1);
                                            on_coef2.set_val(cur_idx,1);
                                        }
                                    }
                                    
                                    
                                    else /*means model_type == "lambda_III" */{
                                        
                                        //fill constraint indices
                                        for (int i=0; i<(num_partns1-2)*2+2; ++i){
                                            const_idx1.add(to_string(i+1));
                                        }
                                        
                                        //fill constraint indices
                                        for (int i=0; i<(num_partns2-2)*2+2; ++i){
                                            const_idx2.add(to_string(i+1));
                                        }
                                        
                                        // Lambda coefficient matrix when linking with partition variables
                                        lambda_coef1.in(indices(inst_partition_lambda, const_idx1));
                                        lambda_coef2.in(indices(inst_partition_lambda, const_idx2));
                                        
                                        // Partition coefficient matrix when linking with lambda variables
                                        on_coef1.in(indices(added, partns1, const_idx1));
                                        on_coef2.in(indices(added, partns2, const_idx2));
                                        
                                        // fill lambda_coef1 and lambda_coef2
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift->get_id_inst(inst);
                                            auto cur_var_idx = added._keys->at(cur_var_id);
                                            for (int j=0; j<num_partns2+1; ++j) {
                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(1);
                                                lambda_coef1.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string((num_partns1-2)*2+2);
                                                lambda_coef1.set_val(cur_idx,1);
                                            }
                                            
                                            for (int i=1 ; i<(num_partns1-2)*2+1; i=i+2) {
                                                for (int j=(i-1)/2 + 2; j<num_partns1+1; ++j) {
                                                    for(int k=0; k<num_partns2+1; ++k){
                                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+name2+"{"+to_string(k+1)+"},"+to_string(i+1);
                                                        lambda_coef1.set_val(cur_idx,1);
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+name2+"{"+to_string(k+1)+"},"+to_string(i+2);
                                                        lambda_coef1.set_val(cur_idx,-1);
                                                    }
                                                }
                                            }
                                            
                                            for (int i=0; i<num_partns1+1; ++i) {
                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(1)+"},"+to_string(1);
                                                lambda_coef2.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(num_partns2+1)+"},"+to_string((num_partns2-2)*2+2);
                                                lambda_coef2.set_val(cur_idx,1);
                                            }
                                            
                                            for (int i=1 ; i<(num_partns2-2)*2+1; i=i+2) {
                                                for (int j=(i-1)/2 + 2; j<num_partns2+1; ++j) {
                                                    for(int k=0; k<num_partns1+1; ++k){
                                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(k+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                        lambda_coef2.set_val(cur_idx,1);
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(k+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                        lambda_coef2.set_val(cur_idx,-1);
                                                    }
                                                }
                                            }
                                        }
                                        
                                        
                                        
                                        // fill on_coef1 and on_coef2
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift->get_id_inst(inst);
                                            auto cur_var_idx = added._keys->at(cur_var_id);
                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                            on_coef1.set_val(cur_idx,1);
                                            
                                            for (int i=1; i<num_partns1; ++i) {
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(2);
                                                on_coef1.set_val(cur_idx, 1);
                                            }
                                            
                                            for (int i=2 ; i<(num_partns1-2)*2+2; i=i+2) {
                                                for (int j=i/2+1; j<num_partns1; ++j) {
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                    on_coef1.set_val(cur_idx,-1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                    on_coef1.set_val(cur_idx,1);
                                                }
                                            }
                                            
                                            cur_idx = cur_var_idx+","+name2+"{"+to_string(1)+"},"+to_string(1);
                                            on_coef2.set_val(cur_idx,1);
                                            for (int i=1; i<num_partns2; ++i) {
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"},"+to_string(2);
                                                on_coef2.set_val(cur_idx, 1);
                                            }
                                            
                                            for (int i=2 ; i<(num_partns2-2)*2+2; i=i+2) {
                                                for (int j=i/2+1; j<num_partns2; ++j) {
                                                    cur_idx = cur_var_idx+","+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                    on_coef2.set_val(cur_idx,-1);
                                                    cur_idx = cur_var_idx+","+name2+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                    on_coef2.set_val(cur_idx,1);
                                                }
                                            }
                                        }
                                    }
                                    
                                    
                                    /** Constraints */
                                    // Representation of the bilinear term with convex combination
                                    Constraint<> bln_rep(pair.first+"_bln_rep");
                                    bln_rep = EP.in_matrix(nb_entries,total_entries-nb_entries)*lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries) - vlift->in(added);
                                    add(bln_rep.in(added) == 0);
                                    
                                    // Representation of o1 with convex combination
                                    Constraint<> o1_rep(pair.first+"_o1_rep");
                                    o1_rep = bounds1.from_ith(0,inst_partition_lambda).in_matrix(nb_entries, 1) * lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries) - o1.in(o1_ids);
                                    add(o1_rep.in(added) == 0);
                                    
                                    // Representation of o2 with convex combination
                                    Constraint<> o2_rep(pair.first+"_o2_rep");
                                    o2_rep = bounds2.in_ignore_ith(nb_entries, 1, inst_partition_lambda).in_matrix(nb_entries,1) * lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries) - o2.in(o2_ids);
                                    add(o2_rep.in(added) == 0);
                                    
                                    // Linking partition variables1 with lambda
                                    if(model_type == "lambda_II"){
                                        Constraint<> on_link_lambda1(pair.first+"_on_link_lambda1_II");
                                        on_link_lambda1 = lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef1.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef1.in_matrix(nb_entries,total_entries-nb_entries) - binvar1->in_ignore_ith(nb_entries_v1,nb_entries_v2,on_coef1.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef1.in_matrix(nb_entries,1);
                                        add(on_link_lambda1.in(indices(added,const_idx1)) <= 0);
                                        
                                        Constraint<> on_link_lambda2(pair.first+"_on_link_lambda2_II");
                                        on_link_lambda2 = lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef2.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef2.in_matrix(nb_entries,total_entries-nb_entries) - binvar1->in_ignore_ith(0,nb_entries_v1,on_coef2.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef2.in_matrix(nb_entries,1);
                                        add(on_link_lambda2.in(indices(added,const_idx2)) <= 0);
                                    }
                                    else{
                                        Constraint<> on_link_lambda1(pair.first+"_on_link_lambda1_III");
                                        on_link_lambda1 = lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef1.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef1.in_matrix(nb_entries,total_entries-nb_entries) - binvar1->in_ignore_ith(nb_entries_v1,nb_entries_v2,on_coef1.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef1.in_matrix(nb_entries,1);
                                        add(on_link_lambda1.in(indices(added,const_idx1)) <= 0);
                                        
                                        Constraint<> on_link_lambda2(pair.first+"_on_link_lambda2_III");
                                        on_link_lambda2 = lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef2.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef2.in_matrix(nb_entries,total_entries-nb_entries) - binvar1->in_ignore_ith(0,nb_entries_v1,on_coef2.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef2.in_matrix(nb_entries,1);
                                        add(on_link_lambda2.in(indices(added,const_idx2)) <= 0);
                                    }
                                    // sum over lambda
                                    Constraint<> lambdaSum(pair.first+"_lambdaSum");
                                    lambdaSum = sum(lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries));
                                    add(lambdaSum.in(added) == 1);
                                }
                                
                            }
                            else{
                                DebugOn("<<<<<<<<<< THIS IS SEEN LIFT -> DOUBLE -> SEEN BOTH BINARIES -> DIFF CORE VARS <<<<<<<<<<<" << endl);
                                
                                auto binvar1 = static_pointer_cast<var<int>>(binvar_ptr1->second);
                                indices partns1("partns1");
                                for (int i = 0; i < num_partns1 ; ++i)
                                {
                                    partns1.add(name1+ "{" +to_string(i+1) + "}");
                                }
                                param<int> lb1("lb1"), ub1("ub1");
                                lb1.in(o1_ids_uq,partns1);
                                ub1.in(o1_ids_uq,partns1);
                                lb1.set_val(0), ub1.set_val(1);
                                auto added1 = binvar1->add_bounds(lb1,ub1);
                                reindex_vars();
                                
                                auto binvar2 = static_pointer_cast<var<int>>(binvar_ptr2->second);
                                indices partns2("partns2");
                                for (int i = 0; i < num_partns2 ; ++i)
                                {
                                    partns2.add(name2+ "{" + to_string(i+1) + "}");
                                }
                                param<int> lb2("lb2"), ub2("ub2");
                                lb2.in(o2_ids_uq,partns2);
                                ub2.in(o2_ids_uq,partns2);
                                lb2.set_val(0), ub2.set_val(1);
                                auto added2 = binvar2->add_bounds(lb2,ub2);
                                reindex_vars();
                                
                                auto nb_entries_v1 = o1_ids.get_nb_entries();
                                auto nb_entries_v2 = o2_ids.get_nb_entries();
                                auto nb_entries = added.get_nb_entries();
                                
                                if(!added1.empty()){
                                    Constraint<> onSum1(o1._name+"_binarySum");
                                    onSum1 = sum(binvar1->in(added1).in_matrix(nb_entries_v1,1));
                                    auto vset1 = added1.from_ith(0,nb_entries_v1);
                                    vset1.filter_refs(vset1.get_unique_refs());
                                    add(onSum1.in(vset1) == 1);
                                }
                                
                                if(!added2.empty()){
                                    Constraint<> onSum2(o2._name+"_binarySum");
                                    onSum2 = sum(binvar2->in(added2).in_matrix(nb_entries_v2,1));
                                    auto vset2 = added2.from_ith(0,nb_entries_v2);
                                    vset2.filter_refs(vset2.get_unique_refs());
                                    add(onSum2.in(vset2) == 1);
                                }
                                
                                if(model_type == "on/off"){//if on/off is chosen
                                    
                                    indices partns("partns");
                                    partns = indices(partns1,partns2);
                                    auto inst_partition = indices(added,partns);
                                    auto total_entries = inst_partition.get_nb_entries();
                                    
                                    if(binvar_ptr3 !=_vars_name.end()){ //means the combined binary variable has been used before
                                        auto binvar3 = static_pointer_cast<var<int>>(binvar_ptr3->second);
                                        param<int> lb3("lb3"), ub3("ub3");
                                        lb3.in(added,partns);
                                        ub3.in(added,partns);
                                        lb3.set_val(0), ub3.set_val(1);
                                        auto added3 = binvar3->add_bounds(lb3,ub3);
                                        reindex_vars();
                                        
                                        Constraint<> onLink1(pair.first+"_binaryLink1");
                                        onLink1 = binvar1->from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) - binvar3->in(inst_partition);
                                        add(onLink1.in(inst_partition) >= 0);
                                        
                                        Constraint<> onLink2(pair.first+"_binaryLink2");
                                        onLink2 = binvar2->in_ignore_ith(nb_entries_v2,1,inst_partition.ignore_ith(0,nb_entries_v1)) - binvar3->in(inst_partition);
                                        add(onLink2.in(inst_partition) >= 0);
                                        
                                        Constraint<> onLink3(pair.first+"_binaryLink3");
                                        onLink3 = binvar1->from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) + binvar2->in_ignore_ith(nb_entries_v2,1,inst_partition.ignore_ith(0,nb_entries_v1)) - 1 - binvar3->in(inst_partition);
                                        add(onLink3.in(inst_partition) <= 0);
                                        
                                        Constraint<> onSumComb(pair.first+"_binarySum");
                                        onSumComb = sum((binvar3->in(added3)).in_matrix(nb_entries,total_entries-nb_entries));
                                        add(onSumComb.in(added) == 1);
                                        
                                        add_on_off_McCormick_refined(pair.first, vlift->in(added), o1.in(o1_ids), o2.in(o2_ids), binvar3->in(added3));
                                    }
                                    else{ //means the combined binary variable has not been used before
                                        var<int> on(name1+name2+"_binary",0,1);
                                        add(on.in(inst_partition));
                                        
                                        Constraint<> onLink1(pair.first+"_binaryLink1");
                                        onLink1 = binvar1->from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) - on;
                                        add(onLink1.in(inst_partition) >= 0);
                                        
                                        Constraint<> onLink2(pair.first+"_binaryLink2");
                                        onLink2 = binvar2->in_ignore_ith(nb_entries_v2,1,inst_partition.ignore_ith(0,nb_entries_v1)) - on;
                                        add(onLink2.in(inst_partition) >= 0);
                                        
                                        Constraint<> onLink3(pair.first+"_binaryLink3");
                                        onLink3 = binvar1->from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) + binvar2->in_ignore_ith(nb_entries_v2,1,inst_partition.ignore_ith(0,nb_entries_v1)) - 1 - on;
                                        add(onLink3.in(inst_partition) <= 0);
                                        
                                        Constraint<> onSumComb(pair.first+"_binarySum");
                                        onSumComb = sum(on.in_matrix(nb_entries,total_entries-nb_entries));
                                        add(onSumComb.in(added) == 1);
                                        
                                        add_on_off_McCormick_refined(pair.first, vlift->in(added), o1.in(o1_ids), o2.in(o2_ids), on);
                                    }
                                }
                                
                                else{//means it is one of the lambda formulations
                                    
                                    //difference is this has one more partition index
                                    indices partns1_lambda("partns1_lambda");
                                    for (int i = 0; i < num_partns1+1; ++i)
                                    {
                                        partns1_lambda.add(name1+ "{" +to_string(i+1) + "}");
                                    }
                                    
                                    indices partns2_lambda("partns2_lambda");
                                    for (int i = 0; i < num_partns2+1; ++i)
                                    {
                                        partns2_lambda.add(name2+ "{" +to_string(i+1) + "}");
                                    }
                                    
                                    indices partns_lambda("partns_lambda");
                                    partns_lambda = indices(partns1_lambda,partns2_lambda);
                                    auto inst_partition_lambda = indices(added,partns_lambda);
                                    auto inst_partition_bounds1 = indices(added,partns1_lambda);
                                    auto inst_partition_bounds2 = indices(added,partns2_lambda);
                                    
                                    // Convex combination variables
                                    auto lambda_ptr = _vars_name.find(name1+name2+"_lambda");
                                    auto lambda = static_pointer_cast<var<double>>(lambda_ptr->second);
                                    param<double> lb_lambda("lb_lambda"), ub_lambda("ub_lambda");
                                    lb_lambda.in(added,partns_lambda);
                                    ub_lambda.in(added,partns_lambda);
                                    lb_lambda.set_val(0), ub_lambda.set_val(1);
                                    auto added_lambda = lambda->add_bounds(lb_lambda,ub_lambda);
                                    reindex_vars();
                                    
                                    /** Parameters */
                                    // Bounds on variable v1 & v2
                                    param<> bounds1(name1+"_bounds1");
                                    bounds1.in(inst_partition_bounds1);
                                    
                                    param<> bounds2(name2+"_bounds2");
                                    bounds2.in(inst_partition_bounds2);
                                    
                                    // Function values on the extreme points
                                    param<> EP(name1+name2+"_grid_values");
                                    EP.in(inst_partition_lambda);
                                    auto total_entries = inst_partition_lambda.get_nb_entries();
                                    
                                    size_t nb_ins = vlift->in(added).get_nb_inst();
                                    auto o1_global_lb = o1.get_lb();
                                    auto increment1 = (o1.get_ub() - o1_global_lb)/num_partns1;
                                    
                                    auto o2_global_lb = o2.get_lb();
                                    auto increment2 = (o2.get_ub() - o2_global_lb)/num_partns2;
                                    
                                    // fill bounds1 and function values
                                    for (int i=0 ; i<num_partns1+1; ++i) {
                                        auto bound_partn1 = o1_global_lb + increment1*i;
                                        bound_partn1.eval_all();
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift->get_id_inst(inst);
                                            auto cur_var_idx = added._keys->at(cur_var_id);
                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"}";
                                            bounds1.set_val(cur_idx,bound_partn1.eval(inst));
                                            for(int j=0; j<num_partns2+1; ++j){
                                                auto bound_partn2 = o2_global_lb + increment2*j;
                                                bound_partn2.eval_all();
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"}";
                                                EP.set_val(cur_idx,(bound_partn1.eval(inst)*bound_partn2.eval(inst)));
                                            }
                                        }
                                    }
                                    // fill bounds2
                                    for (int i=0 ; i<num_partns2+1; ++i) {
                                        auto bound_partn2 = o2_global_lb + increment2*i;
                                        bound_partn2.eval_all();
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift->get_id_inst(inst);
                                            auto cur_var_idx = added._keys->at(cur_var_id);
                                            string cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"}";
                                            bounds2.set_val(cur_idx,bound_partn2.eval(inst));
                                        }
                                    }
                                    
                                    // Lambda coefficient matrix when linking with partition variables
                                    param<> lambda_coef1(name1+"_lambda_linking_coefficients1");
                                    param<> lambda_coef2(name2+"_lambda_linking_coefficients2");
                                    
                                    // Partition coefficient matrix when linking with lambda variables
                                    param<> on_coef1(name1+"_partition_linking_coefficients1");
                                    param<> on_coef2(name2+"_partition_linking_coefficients2");
                                    
                                    // create constraint indices
                                    indices const_idx1("const_idx1");
                                    indices const_idx2("const_idx2");
                                    
                                    if(model_type == "lambda_II"){
                                        
                                        //fill constraint indices
                                        for (int i=0; i<num_partns1+1; ++i){
                                            const_idx1.add(to_string(i+1));
                                        }
                                        
                                        //fill constraint indices
                                        for (int i=0; i<num_partns2+1; ++i){
                                            const_idx2.add(to_string(i+1));
                                        }
                                        
                                        // Lambda coefficient matrix when linking with partition variables
                                        if(num_partns1 > 1) lambda_coef1.in(indices(inst_partition_lambda, const_idx1));
                                        if(num_partns2 > 1) lambda_coef2.in(indices(inst_partition_lambda, const_idx2));
                                        
                                        // Lambda coefficient matrix when linking with partition variables
                                        if(num_partns1 > 1) on_coef1.in(indices(added, partns1, const_idx1));
                                        if(num_partns2 > 1) on_coef2.in(indices(added, partns2, const_idx2));
                                        
                                        // fill lambda_coef1 and lambda_coef2
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift->get_id_inst(inst);
                                            auto cur_var_idx = added._keys->at(cur_var_id);
                                            for (int i=0 ; i<num_partns1+1; ++i) {
                                                for (int j=0 ; j<num_partns2+1; ++j) {
                                                    string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                    if(num_partns1 > 1) lambda_coef1.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(j+1);
                                                    if(num_partns2 > 1) lambda_coef2.set_val(cur_idx,1);
                                                }
                                            }
                                        }
                                        
                                        // fill on_coef1 and on_coef2
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift->get_id_inst(inst);
                                            auto cur_var_idx = added._keys->at(cur_var_id);
                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                            if(num_partns1 > 1) on_coef1.set_val(cur_idx,1);
                                            cur_idx = cur_var_idx+","+name2+"{"+to_string(1)+"},"+to_string(1);
                                            if(num_partns2 > 1) on_coef2.set_val(cur_idx,1);
                                            if(num_partns1 > 1) {
                                                for (int i=1 ; i<num_partns1; ++i) {
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i)+"},"+to_string(i+1);
                                                    on_coef1.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(i+1);
                                                    on_coef1.set_val(cur_idx,1);
                                                }
                                            }
                                            if(num_partns2 > 1) {
                                                for (int i=1 ; i<num_partns2; ++i) {
                                                    cur_idx = cur_var_idx+","+name2+"{"+to_string(i)+"},"+to_string(i+1);
                                                    on_coef2.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"},"+to_string(i+1);
                                                    on_coef2.set_val(cur_idx,1);
                                                }
                                            }
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1)+"},"+to_string(num_partns1+1);
                                            if(num_partns1 > 1) on_coef1.set_val(cur_idx,1);
                                            cur_idx = cur_var_idx+","+name2+"{"+to_string(num_partns2)+"},"+to_string(num_partns2+1);
                                            if(num_partns2 > 1) on_coef2.set_val(cur_idx,1);
                                        }
                                    }
                                    
                                    
                                    else /*means model_type == "lambda_III" */{
                                        
                                        //fill constraint indices
                                        for (int i=0; i<(num_partns1-2)*2+2; ++i){
                                            const_idx1.add(to_string(i+1));
                                        }
                                        
                                        //fill constraint indices
                                        for (int i=0; i<(num_partns2-2)*2+2; ++i){
                                            const_idx2.add(to_string(i+1));
                                        }
                                        
                                        // Lambda coefficient matrix when linking with partition variables
                                        if(num_partns1 > 1) lambda_coef1.in(indices(inst_partition_lambda, const_idx1));
                                        if(num_partns2 > 1) lambda_coef2.in(indices(inst_partition_lambda, const_idx2));
                                        
                                        // Partition coefficient matrix when linking with lambda variables
                                        if(num_partns1 > 1) on_coef1.in(indices(added, partns1, const_idx1));
                                        if(num_partns2 > 1) on_coef2.in(indices(added, partns2, const_idx2));
                                        
                                        // fill lambda_coef1 and lambda_coef2
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift->get_id_inst(inst);
                                            auto cur_var_idx = added._keys->at(cur_var_id);
                                            if(num_partns1 > 1) {
                                                for (int j=0; j<num_partns2+1; ++j) {
                                                    string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(1);
                                                    lambda_coef1.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string((num_partns1-2)*2+2);
                                                    lambda_coef1.set_val(cur_idx,1);
                                                }
                                                
                                                for (int i=1 ; i<(num_partns1-2)*2+1; i=i+2) {
                                                    for (int j=(i-1)/2 + 2; j<num_partns1+1; ++j) {
                                                        for(int k=0; k<num_partns2+1; ++k){
                                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+name2+"{"+to_string(k+1)+"},"+to_string(i+1);
                                                            lambda_coef1.set_val(cur_idx,1);
                                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+name2+"{"+to_string(k+1)+"},"+to_string(i+2);
                                                            lambda_coef1.set_val(cur_idx,-1);
                                                        }
                                                    }
                                                }
                                            }
                                            if(num_partns2 > 1) {
                                                for (int i=0; i<num_partns1+1; ++i) {
                                                    string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(1)+"},"+to_string(1);
                                                    lambda_coef2.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(num_partns2+1)+"},"+to_string((num_partns2-2)*2+2);
                                                    lambda_coef2.set_val(cur_idx,1);
                                                }
                                                
                                                for (int i=1 ; i<(num_partns2-2)*2+1; i=i+2) {
                                                    for (int j=(i-1)/2 + 2; j<num_partns2+1; ++j) {
                                                        for(int k=0; k<num_partns1+1; ++k){
                                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(k+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                            lambda_coef2.set_val(cur_idx,1);
                                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(k+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                            lambda_coef2.set_val(cur_idx,-1);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        
                                        
                                        
                                        // fill on_coef1 and on_coef2
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift->get_id_inst(inst);
                                            auto cur_var_idx = added._keys->at(cur_var_id);
                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                            if(num_partns1 > 1) {
                                                on_coef1.set_val(cur_idx,1);
                                                
                                                for (int i=1; i<num_partns1; ++i) {
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(2);
                                                    on_coef1.set_val(cur_idx, 1);
                                                }
                                                
                                                for (int i=2 ; i<(num_partns1-2)*2+2; i=i+2) {
                                                    for (int j=i/2+1; j<num_partns1; ++j) {
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                        on_coef1.set_val(cur_idx,-1);
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                        on_coef1.set_val(cur_idx,1);
                                                    }
                                                }
                                            }
                                            if(num_partns2 > 1) {
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(1)+"},"+to_string(1);
                                                on_coef2.set_val(cur_idx,1);
                                                for (int i=1; i<num_partns2; ++i) {
                                                    cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"},"+to_string(2);
                                                    on_coef2.set_val(cur_idx, 1);
                                                }
                                                
                                                for (int i=2 ; i<(num_partns2-2)*2+2; i=i+2) {
                                                    for (int j=i/2+1; j<num_partns2; ++j) {
                                                        cur_idx = cur_var_idx+","+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                        on_coef2.set_val(cur_idx,-1);
                                                        cur_idx = cur_var_idx+","+name2+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                        on_coef2.set_val(cur_idx,1);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    
                                    
                                    /** Constraints */
                                    // Representation of the bilinear term with convex combination
                                    Constraint<> bln_rep(pair.first+"_bln_rep");
                                    bln_rep = EP.in_matrix(nb_entries,total_entries-nb_entries)*lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries) - vlift->in(added);
                                    add(bln_rep.in(added) == 0);
                                    
                                    // Representation of o1 with convex combination
                                    Constraint<> o1_rep(pair.first+"_o1_rep");
                                    o1_rep = bounds1.from_ith(0,inst_partition_lambda).in_matrix(nb_entries, 1) * lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries) - o1.in(o1_ids);
                                    add(o1_rep.in(added) == 0);
                                    
                                    // Representation of o2 with convex combination
                                    Constraint<> o2_rep(pair.first+"_o2_rep");
                                    o2_rep = bounds2.in_ignore_ith(nb_entries, 1, inst_partition_lambda).in_matrix(nb_entries,1) * lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries) - o2.in(o2_ids);
                                    add(o2_rep.in(added) == 0);
                                    
                                    // Linking partition variables1 with lambda
                                    if(model_type == "lambda_II"){
                                        if(num_partns1 > 1) {
                                            Constraint<> on_link_lambda1(pair.first+"_on_link_lambda1_II");
                                            on_link_lambda1 = lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef1.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef1.in_matrix(nb_entries,total_entries-nb_entries) - binvar1->in_ignore_ith(nb_entries_v1,nb_entries_v2,on_coef1.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef1.in_matrix(nb_entries,1);
                                            add(on_link_lambda1.in(indices(added,const_idx1)) <= 0);
                                        }
                                        if(num_partns2 > 1) {
                                            Constraint<> on_link_lambda2(pair.first+"_on_link_lambda2_II");
                                            on_link_lambda2 = lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef2.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef2.in_matrix(nb_entries,total_entries-nb_entries) - binvar2->in_ignore_ith(0,nb_entries_v1,on_coef2.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef2.in_matrix(nb_entries,1);
                                            add(on_link_lambda2.in(indices(added,const_idx2)) <= 0);
                                        }
                                    }
                                    else{
                                        if(num_partns1 > 1) {
                                            Constraint<> on_link_lambda1(pair.first+"_on_link_lambda1_III");
                                            on_link_lambda1 = lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef1.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef1.in_matrix(nb_entries,total_entries-nb_entries) - binvar1->in_ignore_ith(nb_entries_v1,nb_entries_v2,on_coef1.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef1.in_matrix(nb_entries,1);
                                            add(on_link_lambda1.in(indices(added,const_idx1)) <= 0);
                                        }
                                        if(num_partns2 > 1) {
                                            Constraint<> on_link_lambda2(pair.first+"_on_link_lambda2_III");
                                            on_link_lambda2 = lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef2.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef2.in_matrix(nb_entries,total_entries-nb_entries) - binvar2->in_ignore_ith(0,nb_entries_v1,on_coef2.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef2.in_matrix(nb_entries,1);
                                            add(on_link_lambda2.in(indices(added,const_idx2)) <= 0);
                                        }
                                    }
                                    // sum over lambda
                                    Constraint<> lambdaSum(pair.first+"_lambdaSum");
                                    lambdaSum = sum(lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries));
                                    add(lambdaSum.in(added) == 1);
                                }
                            }
                        }
                        else if(binvar_ptr1 !=_vars_name.end() && binvar_ptr2 ==_vars_name.end()){ //means v2 has not been partitioned before
                            {
                                DebugOn("<<<<<<<<<< THIS IS SEEN LIFT -> DOUBLE -> SEEN FIRST BINARY -> DIFF CORE VARS <<<<<<<<<<<" << endl);
                                
                                auto binvar1 = static_pointer_cast<var<int>>(binvar_ptr1->second);
                                indices partns1("partns1");
                                for (int i = 0; i < num_partns1 ; ++i)
                                {
                                    partns1.add(name1+ "{" +to_string(i+1) + "}");
                                }
                                param<int> lb1("lb1"), ub1("ub1");
                                lb1.in(o1_ids_uq,partns1);
                                ub1.in(o1_ids_uq,partns1);
                                lb1.set_val(0), ub1.set_val(1);
                                auto added1 = binvar1->add_bounds(lb1,ub1);
                                reindex_vars();
                                
                                var<int> on2(name2+"_binary",0,1);
                                indices partns2("partns2");
                                for (int i = 0; i < num_partns2 ; ++i)
                                {
                                    partns2.add(name2+ "{" + to_string(i+1) + "}");
                                }
                                add(on2.in(o2_ids_uq,partns2));
                                
                                auto nb_entries_v1 = o1_ids.get_nb_entries();
                                auto nb_entries_v2 = o2_ids.get_nb_entries();
                                auto nb_entries = unique_ids.get_nb_entries();
                                
                                if(!added1.empty()){
                                    Constraint<> onSum1(o1._name+"_binarySum");
                                    onSum1 = sum(binvar1->in(added1).in_matrix(nb_entries_v1,1));
                                    auto vset1 = added1.from_ith(0,nb_entries_v1);
                                    vset1.filter_refs(vset1.get_unique_refs());
                                    add(onSum1.in(vset1) == 1);
                                }
                                
                                Constraint<> onSum2(o2._name+"_binarySum");
                                onSum2 = sum(on2.in_matrix(nb_entries_v2,1));
                                add(onSum2.in(o2_ids_uq) == 1);
                                
                                if(model_type == "on/off"){//if on/off is chosen
                                    
                                    indices partns("partns");
                                    partns = indices(partns1,partns2);
                                    auto inst_partition = indices(added,partns);
                                    auto total_entries = inst_partition.get_nb_entries();
                                    
                                    if(binvar_ptr3 !=_vars_name.end()){ //means the combined binary variable has been used before
                                        
                                        auto binvar3 = static_pointer_cast<var<int>>(binvar_ptr3->second);
                                        param<int> lb3("lb3"), ub3("ub3");
                                        lb3.in(added,partns);
                                        ub3.in(added,partns);
                                        lb3.set_val(0), ub3.set_val(1);
                                        auto added3 = binvar3->add_bounds(lb3,ub3);
                                        reindex_vars();
                                        
                                        Constraint<> onLink1(pair.first+"_binaryLink1");
                                        onLink1 = binvar1->from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) - binvar3->in(inst_partition);
                                        add(onLink1.in(inst_partition) >= 0);
                                        
                                        Constraint<> onLink2(pair.first+"_binaryLink2");
                                        onLink2 = on2.in_ignore_ith(nb_entries_v2,1,inst_partition.ignore_ith(0,nb_entries_v1)) - binvar3->in(inst_partition);
                                        add(onLink2.in(inst_partition) >= 0);
                                        
                                        Constraint<> onLink3(pair.first+"_binaryLink3");
                                        onLink3 = binvar1->from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) + on2.in_ignore_ith(nb_entries_v2,1,inst_partition.ignore_ith(0,nb_entries_v1)) - 1 - binvar3->in(inst_partition);
                                        add(onLink3.in(inst_partition) <= 0);
                                        
                                        Constraint<> onSumComb(pair.first+"_binarySum");
                                        onSumComb = sum((binvar3->in(added3)).in_matrix(nb_entries,total_entries-nb_entries));
                                        add(onSumComb.in(added) == 1);
                                        
                                        add_on_off_McCormick_refined(pair.first, vlift->in(added), o1.in(o1_ids), o2.in(o2_ids), binvar3->in(added3));
                                    }
                                    else{ //means the combined binary variable has not been used before
                                        var<int> on(name1+name2+"_binary",0,1);
                                        add(on.in(inst_partition));
                                        
                                        Constraint<> onLink1(pair.first+"_binaryLink1");
                                        onLink1 = binvar1->from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) - on;
                                        add(onLink1.in(inst_partition) >= 0);
                                        
                                        Constraint<> onLink2(pair.first+"_binaryLink2");
                                        onLink2 = on2.in_ignore_ith(nb_entries_v2,1,inst_partition.ignore_ith(0,nb_entries_v1)) - on;
                                        add(onLink2.in(inst_partition) >= 0);
                                        
                                        Constraint<> onLink3(pair.first+"_binaryLink3");
                                        onLink3 = binvar1->from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) + on2.in_ignore_ith(nb_entries_v2,1,inst_partition.ignore_ith(0,nb_entries_v1)) - 1 - on;
                                        add(onLink3.in(inst_partition) <= 0);
                                        
                                        Constraint<> onSumComb(pair.first+"_binarySum");
                                        onSumComb = sum(on.in_matrix(nb_entries,total_entries-nb_entries));
                                        add(onSumComb.in(added) == 1);
                                        
                                        add_on_off_McCormick_refined(pair.first, vlift->in(added), o1.in(o1_ids), o2.in(o2_ids), on);
                                    }
                                }
                                
                                else{//means it is one of the lambda formulations
                                    
                                    //difference is this has one more partition index
                                    indices partns1_lambda("partns1_lambda");
                                    for (int i = 0; i < num_partns1+1; ++i)
                                    {
                                        partns1_lambda.add(name1+ "{" +to_string(i+1) + "}");
                                    }
                                    
                                    indices partns2_lambda("partns2_lambda");
                                    for (int i = 0; i < num_partns2+1; ++i)
                                    {
                                        partns2_lambda.add(name2+ "{" +to_string(i+1) + "}");
                                    }
                                    
                                    indices partns_lambda("partns_lambda");
                                    partns_lambda = indices(partns1_lambda,partns2_lambda);
                                    auto inst_partition_lambda = indices(added,partns_lambda);
                                    auto inst_partition_bounds1 = indices(added,partns1_lambda);
                                    auto inst_partition_bounds2 = indices(added,partns2_lambda);
                                    
                                    // Convex combination variables
                                    auto lambda_ptr = _vars_name.find(name1+name2+"_lambda");
                                    auto lambda = static_pointer_cast<var<double>>(lambda_ptr->second);
                                    param<double> lb_lambda("lb_lambda"), ub_lambda("ub_lambda");
                                    lb_lambda.in(added,partns_lambda);
                                    ub_lambda.in(added,partns_lambda);
                                    lb_lambda.set_val(0), ub_lambda.set_val(1);
                                    auto added_lambda = lambda->add_bounds(lb_lambda,ub_lambda);
                                    reindex_vars();
                                    
                                    /** Parameters */
                                    // Bounds on variable v1 & v2
                                    param<> bounds1(name1+"_bounds1");
                                    bounds1.in(inst_partition_bounds1);
                                    
                                    param<> bounds2(name2+"_bounds2");
                                    bounds2.in(inst_partition_bounds2);
                                    
                                    // Function values on the extreme points
                                    param<> EP(name1+name2+"_grid_values");
                                    EP.in(inst_partition_lambda);
                                    auto total_entries = inst_partition_lambda.get_nb_entries();
                                    
                                    size_t nb_ins = vlift->in(added).get_nb_inst();
                                    auto o1_global_lb = o1.get_lb();
                                    auto increment1 = (o1.get_ub() - o1_global_lb)/num_partns1;
                                    
                                    auto o2_global_lb = o2.get_lb();
                                    auto increment2 = (o2.get_ub() - o2_global_lb)/num_partns2;
                                    
                                    // fill bounds1 and function values
                                    for (int i=0 ; i<num_partns1+1; ++i) {
                                        auto bound_partn1 = o1_global_lb + increment1*i;
                                        bound_partn1.eval_all();
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift->get_id_inst(inst);
                                            auto cur_var_idx = added._keys->at(cur_var_id);
                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"}";
                                            bounds1.set_val(cur_idx,bound_partn1.eval(inst));
                                            for(int j=0; j<num_partns2+1; ++j){
                                                auto bound_partn2 = o2_global_lb + increment2*j;
                                                bound_partn2.eval_all();
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"}";
                                                EP.set_val(cur_idx,(bound_partn1.eval(inst)*bound_partn2.eval(inst)));
                                            }
                                        }
                                    }
                                    // fill bounds2
                                    for (int i=0 ; i<num_partns2+1; ++i) {
                                        auto bound_partn2 = o2_global_lb + increment2*i;
                                        bound_partn2.eval_all();
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift->get_id_inst(inst);
                                            auto cur_var_idx = added._keys->at(cur_var_id);
                                            string cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"}";
                                            bounds2.set_val(cur_idx,bound_partn2.eval(inst));
                                        }
                                    }
                                    
                                    // Lambda coefficient matrix when linking with partition variables
                                    param<> lambda_coef1(name1+"_lambda_linking_coefficients1");
                                    param<> lambda_coef2(name2+"_lambda_linking_coefficients2");
                                    
                                    // Partition coefficient matrix when linking with lambda variables
                                    param<> on_coef1(name1+"_partition_linking_coefficients1");
                                    param<> on_coef2(name2+"_partition_linking_coefficients2");
                                    
                                    // create constraint indices
                                    indices const_idx1("const_idx1");
                                    indices const_idx2("const_idx2");
                                    
                                    if(model_type == "lambda_II"){
                                        
                                        //fill constraint indices
                                        for (int i=0; i<num_partns1+1; ++i){
                                            const_idx1.add(to_string(i+1));
                                        }
                                        
                                        //fill constraint indices
                                        for (int i=0; i<num_partns2+1; ++i){
                                            const_idx2.add(to_string(i+1));
                                        }
                                        
                                        // Lambda coefficient matrix when linking with partition variables
                                        if(num_partns1 > 1) lambda_coef1.in(indices(inst_partition_lambda, const_idx1));
                                        if(num_partns2 > 1) lambda_coef2.in(indices(inst_partition_lambda, const_idx2));
                                        
                                        // Lambda coefficient matrix when linking with partition variables
                                        if(num_partns1 > 1) on_coef1.in(indices(added, partns1, const_idx1));
                                        if(num_partns2 > 1) on_coef2.in(indices(added, partns2, const_idx2));
                                        
                                        // fill lambda_coef1 and lambda_coef2
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift->get_id_inst(inst);
                                            auto cur_var_idx = added._keys->at(cur_var_id);
                                            for (int i=0 ; i<num_partns1+1; ++i) {
                                                for (int j=0 ; j<num_partns2+1; ++j) {
                                                    string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                    if(num_partns1 > 1) lambda_coef1.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(j+1);
                                                    if(num_partns2 > 1) lambda_coef2.set_val(cur_idx,1);
                                                }
                                            }
                                        }
                                        
                                        // fill on_coef1 and on_coef2
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift->get_id_inst(inst);
                                            auto cur_var_idx = added._keys->at(cur_var_id);
                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                            if(num_partns1 > 1) on_coef1.set_val(cur_idx,1);
                                            cur_idx = cur_var_idx+","+name2+"{"+to_string(1)+"},"+to_string(1);
                                            if(num_partns2 > 1) on_coef2.set_val(cur_idx,1);
                                            if(num_partns1 > 1) {
                                                for (int i=1 ; i<num_partns1; ++i) {
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i)+"},"+to_string(i+1);
                                                    on_coef1.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(i+1);
                                                    on_coef1.set_val(cur_idx,1);
                                                }
                                            }
                                            if(num_partns2 > 1) {
                                                for (int i=1 ; i<num_partns2; ++i) {
                                                    cur_idx = cur_var_idx+","+name2+"{"+to_string(i)+"},"+to_string(i+1);
                                                    on_coef2.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"},"+to_string(i+1);
                                                    on_coef2.set_val(cur_idx,1);
                                                }
                                            }
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1)+"},"+to_string(num_partns1+1);
                                            if(num_partns1 > 1) on_coef1.set_val(cur_idx,1);
                                            cur_idx = cur_var_idx+","+name2+"{"+to_string(num_partns2)+"},"+to_string(num_partns2+1);
                                            if(num_partns2 > 1) on_coef2.set_val(cur_idx,1);
                                        }
                                    }
                                    
                                    
                                    else /*means model_type == "lambda_III" */{
                                        
                                        //fill constraint indices
                                        for (int i=0; i<(num_partns1-2)*2+2; ++i){
                                            const_idx1.add(to_string(i+1));
                                        }
                                        
                                        //fill constraint indices
                                        for (int i=0; i<(num_partns2-2)*2+2; ++i){
                                            const_idx2.add(to_string(i+1));
                                        }
                                        
                                        // Lambda coefficient matrix when linking with partition variables
                                        if(num_partns1 > 1) lambda_coef1.in(indices(inst_partition_lambda, const_idx1));
                                        if(num_partns2 > 1) lambda_coef2.in(indices(inst_partition_lambda, const_idx2));
                                        
                                        // Partition coefficient matrix when linking with lambda variables
                                        if(num_partns1 > 1) on_coef1.in(indices(added, partns1, const_idx1));
                                        if(num_partns2 > 1) on_coef2.in(indices(added, partns2, const_idx2));
                                        
                                        // fill lambda_coef1 and lambda_coef2
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift->get_id_inst(inst);
                                            auto cur_var_idx = added._keys->at(cur_var_id);
                                            if(num_partns1 > 1) {
                                                for (int j=0; j<num_partns2+1; ++j) {
                                                    string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(1);
                                                    lambda_coef1.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string((num_partns1-2)*2+2);
                                                    lambda_coef1.set_val(cur_idx,1);
                                                }
                                                
                                                for (int i=1 ; i<(num_partns1-2)*2+1; i=i+2) {
                                                    for (int j=(i-1)/2 + 2; j<num_partns1+1; ++j) {
                                                        for(int k=0; k<num_partns2+1; ++k){
                                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+name2+"{"+to_string(k+1)+"},"+to_string(i+1);
                                                            lambda_coef1.set_val(cur_idx,1);
                                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+name2+"{"+to_string(k+1)+"},"+to_string(i+2);
                                                            lambda_coef1.set_val(cur_idx,-1);
                                                        }
                                                    }
                                                }
                                            }
                                            if(num_partns2 > 1) {
                                                for (int i=0; i<num_partns1+1; ++i) {
                                                    string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(1)+"},"+to_string(1);
                                                    lambda_coef2.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(num_partns2+1)+"},"+to_string((num_partns2-2)*2+2);
                                                    lambda_coef2.set_val(cur_idx,1);
                                                }
                                                
                                                for (int i=1 ; i<(num_partns2-2)*2+1; i=i+2) {
                                                    for (int j=(i-1)/2 + 2; j<num_partns2+1; ++j) {
                                                        for(int k=0; k<num_partns1+1; ++k){
                                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(k+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                            lambda_coef2.set_val(cur_idx,1);
                                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(k+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                            lambda_coef2.set_val(cur_idx,-1);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        
                                        
                                        
                                        // fill on_coef1 and on_coef2
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift->get_id_inst(inst);
                                            auto cur_var_idx = added._keys->at(cur_var_id);
                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                            if(num_partns1 > 1) {
                                                on_coef1.set_val(cur_idx,1);
                                                
                                                for (int i=1; i<num_partns1; ++i) {
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(2);
                                                    on_coef1.set_val(cur_idx, 1);
                                                }
                                                
                                                for (int i=2 ; i<(num_partns1-2)*2+2; i=i+2) {
                                                    for (int j=i/2+1; j<num_partns1; ++j) {
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                        on_coef1.set_val(cur_idx,-1);
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                        on_coef1.set_val(cur_idx,1);
                                                    }
                                                }
                                            }
                                            if(num_partns2 > 1) {
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(1)+"},"+to_string(1);
                                                on_coef2.set_val(cur_idx,1);
                                                for (int i=1; i<num_partns2; ++i) {
                                                    cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"},"+to_string(2);
                                                    on_coef2.set_val(cur_idx, 1);
                                                }
                                                
                                                for (int i=2 ; i<(num_partns2-2)*2+2; i=i+2) {
                                                    for (int j=i/2+1; j<num_partns2; ++j) {
                                                        cur_idx = cur_var_idx+","+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                        on_coef2.set_val(cur_idx,-1);
                                                        cur_idx = cur_var_idx+","+name2+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                        on_coef2.set_val(cur_idx,1);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    
                                    
                                    /** Constraints */
                                    // Representation of the bilinear term with convex combination
                                    Constraint<> bln_rep(pair.first+"_bln_rep");
                                    bln_rep = EP.in_matrix(nb_entries,total_entries-nb_entries)*lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries) - vlift->in(added);
                                    add(bln_rep.in(added) == 0);
                                    
                                    // Representation of o1 with convex combination
                                    Constraint<> o1_rep(pair.first+"_o1_rep");
                                    o1_rep = bounds1.from_ith(0,inst_partition_lambda).in_matrix(nb_entries, 1) * lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries) - o1.in(o1_ids);
                                    add(o1_rep.in(added) == 0);
                                    
                                    // Representation of o2 with convex combination
                                    Constraint<> o2_rep(pair.first+"_o2_rep");
                                    o2_rep = bounds2.in_ignore_ith(nb_entries, 1, inst_partition_lambda).in_matrix(nb_entries,1) * lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries) - o2.in(o2_ids);
                                    add(o2_rep.in(added) == 0);
                                    
                                    // Linking partition variables1 with lambda
                                    if(model_type == "lambda_II"){
                                        if(num_partns1 > 1) {
                                            Constraint<> on_link_lambda1(pair.first+"_on_link_lambda1_II");
                                            on_link_lambda1 = lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef1.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef1.in_matrix(nb_entries,total_entries-nb_entries) - binvar1->in_ignore_ith(nb_entries_v1,nb_entries_v2,on_coef1.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef1.in_matrix(nb_entries,1);
                                            add(on_link_lambda1.in(indices(added,const_idx1)) <= 0);
                                        }
                                        if(num_partns2 > 1) {
                                            Constraint<> on_link_lambda2(pair.first+"_on_link_lambda2_II");
                                            on_link_lambda2 = lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef2.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef2.in_matrix(nb_entries,total_entries-nb_entries) - on2.in_ignore_ith(0,nb_entries_v1,on_coef2.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef2.in_matrix(nb_entries,1);
                                            add(on_link_lambda2.in(indices(added,const_idx2)) <= 0);
                                        }
                                    }
                                    else{
                                        if(num_partns1 > 1) {
                                            Constraint<> on_link_lambda1(pair.first+"_on_link_lambda1_III");
                                            on_link_lambda1 = lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef1.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef1.in_matrix(nb_entries,total_entries-nb_entries) - binvar1->in_ignore_ith(nb_entries_v1,nb_entries_v2,on_coef1.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef1.in_matrix(nb_entries,1);
                                            add(on_link_lambda1.in(indices(added,const_idx1)) <= 0);
                                        }
                                        if(num_partns2 > 1) {
                                            Constraint<> on_link_lambda2(pair.first+"_on_link_lambda2_III");
                                            on_link_lambda2 = lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef2.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef2.in_matrix(nb_entries,total_entries-nb_entries) - on2.in_ignore_ith(0,nb_entries_v1,on_coef2.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef2.in_matrix(nb_entries,1);
                                            add(on_link_lambda2.in(indices(added,const_idx2)) <= 0);
                                        }
                                    }
                                    // sum over lambda
                                    Constraint<> lambdaSum(pair.first+"_lambdaSum");
                                    lambdaSum = sum(lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries));
                                    add(lambdaSum.in(added) == 1);
                                }
                            }
                        }
                        else if(binvar_ptr1 ==_vars_name.end() && binvar_ptr2 !=_vars_name.end()){ //means v1 has not been partitioned before
                             DebugOn("<<<<<<<<<< THIS IS SEEN LIFT -> DOUBLE -> SEEN SECOND BINARY -> DIFF CORE VARS <<<<<<<<<<<" << endl);
                            
                            var<int> on1(name1+"_binary",0,1);
                            indices partns1("partns1");
                            for (int i = 0; i < num_partns1 ; ++i)
                            {
                                partns1.add(name1+ "{" + to_string(i+1) + "}");
                            }
                            add(on1.in(o1_ids_uq,partns1));
                            
                            auto binvar2 = static_pointer_cast<var<int>>(binvar_ptr2->second);
                            indices partns2("partns2");
                            for (int i = 0; i < num_partns2 ; ++i)
                            {
                                partns2.add(name2+ "{" + to_string(i+1) + "}");
                            }
                            param<int> lb2("lb2"), ub2("ub2");
                            lb2.in(o2_ids_uq,partns2);
                            ub2.in(o2_ids_uq,partns2);
                            lb2.set_val(0), ub2.set_val(1);
                            auto added2 = binvar2->add_bounds(lb2,ub2);
                            reindex_vars();
                            
                            auto nb_entries_v1 = o1_ids.get_nb_entries();
                            auto nb_entries_v2 = o2_ids.get_nb_entries();
                            auto nb_entries = added.get_nb_entries();
                            
                            Constraint<> onSum1(o1._name+"_binarySum");
                            onSum1 = sum(on1.in_matrix(nb_entries_v1,1));
                            add(onSum1.in(o1_ids_uq) == 1);
                            
                            if(!added2.empty()){
                                Constraint<> onSum2(o2._name+"_binarySum");
                                onSum2 = sum(binvar2->in(added2).in_matrix(nb_entries_v2,1));
                                auto vset2 = added2.from_ith(0,nb_entries_v2);
                                vset2.filter_refs(vset2.get_unique_refs());
                                add(onSum2.in(vset2) == 1);
                            }
                            
                            if(model_type == "on/off"){//if on/off is chosen
                                
                                indices partns("partns");
                                partns = indices(partns1,partns2);
                                auto inst_partition = indices(added,partns);
                                auto total_entries = inst_partition.get_nb_entries();
                                
                                if(binvar_ptr3 !=_vars_name.end()){ //means the combined binary variable has been used before
                                    auto binvar3 = static_pointer_cast<var<int>>(binvar_ptr3->second);
                                    param<int> lb3("lb3"), ub3("ub3");
                                    lb3.in(added,partns);
                                    ub3.in(added,partns);
                                    lb3.set_val(0), ub3.set_val(1);
                                    auto added3 = binvar3->add_bounds(lb3,ub3);
                                    reindex_vars();
                                    
                                    Constraint<> onLink1(pair.first+"_binaryLink1");
                                    onLink1 = on1.from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) - binvar3->in(inst_partition);
                                    add(onLink1.in(inst_partition) >= 0);
                                    
                                    Constraint<> onLink2(pair.first+"_binaryLink2");
                                    onLink2 = binvar2->in_ignore_ith(nb_entries_v2,1,inst_partition.ignore_ith(0,nb_entries_v1)) - binvar3->in(inst_partition);
                                    add(onLink2.in(inst_partition) >= 0);
                                    
                                    Constraint<> onLink3(pair.first+"_binaryLink3");
                                    onLink3 = on1.from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) + binvar2->in_ignore_ith(nb_entries_v2,1,inst_partition.ignore_ith(0,nb_entries_v1)) - 1 - binvar3->in(inst_partition);
                                    add(onLink3.in(inst_partition) <= 0);
                                    
                                    Constraint<> onSumComb(pair.first+"_binarySum");
                                    onSumComb = sum((binvar3->in(added3)).in_matrix(nb_entries,total_entries-nb_entries));
                                    add(onSumComb.in(added) == 1);
                                    
                                    add_on_off_McCormick_refined(pair.first, vlift->in(added), o1.in(o1_ids), o2.in(o2_ids), binvar3->in(added3));
                                }
                                
                                else{ //means the combined binary variable has not been used before
                                    var<int> on(name1+name2+"_binary",0,1);
                                    add(on.in(inst_partition));
                                    
                                    Constraint<> onLink1(pair.first+"_binaryLink1");
                                    onLink1 = on1.from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) - on;
                                    add(onLink1.in(inst_partition) >= 0);
                                    
                                    Constraint<> onLink2(pair.first+"_binaryLink2");
                                    onLink2 = binvar2->in_ignore_ith(nb_entries_v2,1,inst_partition.ignore_ith(0,nb_entries_v1)) - on;
                                    add(onLink2.in(inst_partition) >= 0);
                                    
                                    Constraint<> onLink3(pair.first+"_binaryLink3");
                                    onLink3 = on1.from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) + binvar2->in_ignore_ith(nb_entries_v2,1,inst_partition.ignore_ith(0,nb_entries_v1)) - 1 - on;
                                    add(onLink3.in(inst_partition) <= 0);
                                    
                                    Constraint<> onSumComb(pair.first+"_binarySum");
                                    onSumComb = sum(on.in_matrix(nb_entries,total_entries-nb_entries));
                                    add(onSumComb.in(added) == 1);
                                    
                                    add_on_off_McCormick_refined(pair.first, vlift->in(added), o1.in(o1_ids), o2.in(o2_ids), on);
                                }
                                
                            }
                            
                            else{//means it is one of the lambda formulations
                                
                                //difference is this has one more partition index
                                indices partns1_lambda("partns1_lambda");
                                for (int i = 0; i < num_partns1+1; ++i)
                                {
                                    partns1_lambda.add(name1+ "{" +to_string(i+1) + "}");
                                }
                                
                                indices partns2_lambda("partns2_lambda");
                                for (int i = 0; i < num_partns2+1; ++i)
                                {
                                    partns2_lambda.add(name2+ "{" +to_string(i+1) + "}");
                                }
                                
                                indices partns_lambda("partns_lambda");
                                partns_lambda = indices(partns1_lambda,partns2_lambda);
                                auto inst_partition_lambda = indices(added,partns_lambda);
                                auto inst_partition_bounds1 = indices(added,partns1_lambda);
                                auto inst_partition_bounds2 = indices(added,partns2_lambda);
                                
                                // Convex combination variables
                                auto lambda_ptr = _vars_name.find(name1+name2+"_lambda");
                                auto lambda = static_pointer_cast<var<double>>(lambda_ptr->second);
                                param<double> lb_lambda("lb_lambda"), ub_lambda("ub_lambda");
                                lb_lambda.in(added,partns_lambda);
                                ub_lambda.in(added,partns_lambda);
                                lb_lambda.set_val(0), ub_lambda.set_val(1);
                                auto added_lambda = lambda->add_bounds(lb_lambda,ub_lambda);
                                reindex_vars();
                                
                                /** Parameters */
                                // Bounds on variable v1 & v2
                                param<> bounds1(name1+"_bounds1");
                                bounds1.in(inst_partition_bounds1);
                                
                                param<> bounds2(name2+"_bounds2");
                                bounds2.in(inst_partition_bounds2);
                                
                                // Function values on the extreme points
                                param<> EP(name1+name2+"_grid_values");
                                EP.in(inst_partition_lambda);
                                auto total_entries = inst_partition_lambda.get_nb_entries();
                                
                                size_t nb_ins = vlift->in(added).get_nb_inst();
                                auto o1_global_lb = o1.get_lb();
                                auto increment1 = (o1.get_ub() - o1_global_lb)/num_partns1;
                                
                                auto o2_global_lb = o2.get_lb();
                                auto increment2 = (o2.get_ub() - o2_global_lb)/num_partns2;
                                
                                // fill bounds1 and function values
                                for (int i=0 ; i<num_partns1+1; ++i) {
                                    auto bound_partn1 = o1_global_lb + increment1*i;
                                    bound_partn1.eval_all();
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift->get_id_inst(inst);
                                        auto cur_var_idx = added._keys->at(cur_var_id);
                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"}";
                                        bounds1.set_val(cur_idx,bound_partn1.eval(inst));
                                        for(int j=0; j<num_partns2+1; ++j){
                                            auto bound_partn2 = o2_global_lb + increment2*j;
                                            bound_partn2.eval_all();
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"}";
                                            EP.set_val(cur_idx,(bound_partn1.eval(inst)*bound_partn2.eval(inst)));
                                        }
                                    }
                                }
                                // fill bounds2
                                for (int i=0 ; i<num_partns2+1; ++i) {
                                    auto bound_partn2 = o2_global_lb + increment2*i;
                                    bound_partn2.eval_all();
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift->get_id_inst(inst);
                                        auto cur_var_idx = added._keys->at(cur_var_id);
                                        string cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"}";
                                        bounds2.set_val(cur_idx,bound_partn2.eval(inst));
                                    }
                                }
                                
                                // Lambda coefficient matrix when linking with partition variables
                                param<> lambda_coef1(name1+"_lambda_linking_coefficients1");
                                param<> lambda_coef2(name2+"_lambda_linking_coefficients2");
                                
                                // Partition coefficient matrix when linking with lambda variables
                                param<> on_coef1(name1+"_partition_linking_coefficients1");
                                param<> on_coef2(name2+"_partition_linking_coefficients2");
                                
                                // create constraint indices
                                indices const_idx1("const_idx1");
                                indices const_idx2("const_idx2");
                                
                                if(model_type == "lambda_II"){
                                    
                                    //fill constraint indices
                                    for (int i=0; i<num_partns1+1; ++i){
                                        const_idx1.add(to_string(i+1));
                                    }
                                    
                                    //fill constraint indices
                                    for (int i=0; i<num_partns2+1; ++i){
                                        const_idx2.add(to_string(i+1));
                                    }
                                    
                                    // Lambda coefficient matrix when linking with partition variables
                                    if(num_partns1 > 1) lambda_coef1.in(indices(inst_partition_lambda, const_idx1));
                                    if(num_partns2 > 1) lambda_coef2.in(indices(inst_partition_lambda, const_idx2));
                                    
                                    // Lambda coefficient matrix when linking with partition variables
                                    if(num_partns1 > 1) on_coef1.in(indices(added, partns1, const_idx1));
                                    if(num_partns2 > 1) on_coef2.in(indices(added, partns2, const_idx2));
                                    
                                    // fill lambda_coef1 and lambda_coef2
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift->get_id_inst(inst);
                                        auto cur_var_idx = added._keys->at(cur_var_id);
                                        for (int i=0 ; i<num_partns1+1; ++i) {
                                            for (int j=0 ; j<num_partns2+1; ++j) {
                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                if(num_partns1 > 1) lambda_coef1.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(j+1);
                                                if(num_partns2 > 1) lambda_coef2.set_val(cur_idx,1);
                                            }
                                        }
                                    }
                                    
                                    // fill on_coef1 and on_coef2
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift->get_id_inst(inst);
                                        auto cur_var_idx = added._keys->at(cur_var_id);
                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                        if(num_partns1 > 1) on_coef1.set_val(cur_idx,1);
                                        cur_idx = cur_var_idx+","+name2+"{"+to_string(1)+"},"+to_string(1);
                                        if(num_partns2 > 1) on_coef2.set_val(cur_idx,1);
                                        if(num_partns1 > 1) {
                                            for (int i=1 ; i<num_partns1; ++i) {
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(i)+"},"+to_string(i+1);
                                                on_coef1.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(i+1);
                                                on_coef1.set_val(cur_idx,1);
                                            }
                                        }
                                        if(num_partns2 > 1) {
                                            for (int i=1 ; i<num_partns2; ++i) {
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(i)+"},"+to_string(i+1);
                                                on_coef2.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"},"+to_string(i+1);
                                                on_coef2.set_val(cur_idx,1);
                                            }
                                        }
                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1)+"},"+to_string(num_partns1+1);
                                        if(num_partns1 > 1) on_coef1.set_val(cur_idx,1);
                                        cur_idx = cur_var_idx+","+name2+"{"+to_string(num_partns2)+"},"+to_string(num_partns2+1);
                                        if(num_partns2 > 1) on_coef2.set_val(cur_idx,1);
                                    }
                                }
                                
                                
                                else /*means model_type == "lambda_III" */{
                                    
                                    //fill constraint indices
                                    for (int i=0; i<(num_partns1-2)*2+2; ++i){
                                        const_idx1.add(to_string(i+1));
                                    }
                                    
                                    //fill constraint indices
                                    for (int i=0; i<(num_partns2-2)*2+2; ++i){
                                        const_idx2.add(to_string(i+1));
                                    }
                                    
                                    // Lambda coefficient matrix when linking with partition variables
                                    if(num_partns1 > 1) lambda_coef1.in(indices(inst_partition_lambda, const_idx1));
                                    if(num_partns2 > 1) lambda_coef2.in(indices(inst_partition_lambda, const_idx2));
                                    
                                    // Partition coefficient matrix when linking with lambda variables
                                    if(num_partns1 > 1) on_coef1.in(indices(added, partns1, const_idx1));
                                    if(num_partns2 > 1) on_coef2.in(indices(added, partns2, const_idx2));
                                    
                                    // fill lambda_coef1 and lambda_coef2
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift->get_id_inst(inst);
                                        auto cur_var_idx = added._keys->at(cur_var_id);
                                        if(num_partns1 > 1) {
                                            for (int j=0; j<num_partns2+1; ++j) {
                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(1);
                                                lambda_coef1.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string((num_partns1-2)*2+2);
                                                lambda_coef1.set_val(cur_idx,1);
                                            }
                                            
                                            for (int i=1 ; i<(num_partns1-2)*2+1; i=i+2) {
                                                for (int j=(i-1)/2 + 2; j<num_partns1+1; ++j) {
                                                    for(int k=0; k<num_partns2+1; ++k){
                                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+name2+"{"+to_string(k+1)+"},"+to_string(i+1);
                                                        lambda_coef1.set_val(cur_idx,1);
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+name2+"{"+to_string(k+1)+"},"+to_string(i+2);
                                                        lambda_coef1.set_val(cur_idx,-1);
                                                    }
                                                }
                                            }
                                        }
                                        if(num_partns2 > 1) {
                                            for (int i=0; i<num_partns1+1; ++i) {
                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(1)+"},"+to_string(1);
                                                lambda_coef2.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(num_partns2+1)+"},"+to_string((num_partns2-2)*2+2);
                                                lambda_coef2.set_val(cur_idx,1);
                                            }
                                            
                                            for (int i=1 ; i<(num_partns2-2)*2+1; i=i+2) {
                                                for (int j=(i-1)/2 + 2; j<num_partns2+1; ++j) {
                                                    for(int k=0; k<num_partns1+1; ++k){
                                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(k+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                        lambda_coef2.set_val(cur_idx,1);
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(k+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                        lambda_coef2.set_val(cur_idx,-1);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    
                                    
                                    
                                    // fill on_coef1 and on_coef2
                                    for (size_t inst = 0; inst< nb_ins; inst++){
                                        auto cur_var_id = vlift->get_id_inst(inst);
                                        auto cur_var_idx = added._keys->at(cur_var_id);
                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                        if(num_partns1 > 1) {
                                            on_coef1.set_val(cur_idx,1);
                                            
                                            for (int i=1; i<num_partns1; ++i) {
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(2);
                                                on_coef1.set_val(cur_idx, 1);
                                            }
                                            
                                            for (int i=2 ; i<(num_partns1-2)*2+2; i=i+2) {
                                                for (int j=i/2+1; j<num_partns1; ++j) {
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                    on_coef1.set_val(cur_idx,-1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                    on_coef1.set_val(cur_idx,1);
                                                }
                                            }
                                        }
                                        if(num_partns2 > 1) {
                                            cur_idx = cur_var_idx+","+name2+"{"+to_string(1)+"},"+to_string(1);
                                            on_coef2.set_val(cur_idx,1);
                                            for (int i=1; i<num_partns2; ++i) {
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"},"+to_string(2);
                                                on_coef2.set_val(cur_idx, 1);
                                            }
                                            
                                            for (int i=2 ; i<(num_partns2-2)*2+2; i=i+2) {
                                                for (int j=i/2+1; j<num_partns2; ++j) {
                                                    cur_idx = cur_var_idx+","+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                    on_coef2.set_val(cur_idx,-1);
                                                    cur_idx = cur_var_idx+","+name2+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                    on_coef2.set_val(cur_idx,1);
                                                }
                                            }
                                        }
                                    }
                                }
                                
                                
                                /** Constraints */
                                // Representation of the bilinear term with convex combination
                                Constraint<> bln_rep(pair.first+"_bln_rep");
                                bln_rep = EP.in_matrix(nb_entries,total_entries-nb_entries)*lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries) - vlift->in(added);
                                add(bln_rep.in(added) == 0);
                                
                                // Representation of o1 with convex combination
                                Constraint<> o1_rep(pair.first+"_o1_rep");
                                o1_rep = bounds1.from_ith(0,inst_partition_lambda).in_matrix(nb_entries, 1) * lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries) - o1.in(o1_ids);
                                add(o1_rep.in(added) == 0);
                                
                                // Representation of o2 with convex combination
                                Constraint<> o2_rep(pair.first+"_o2_rep");
                                o2_rep = bounds2.in_ignore_ith(nb_entries, 1, inst_partition_lambda).in_matrix(nb_entries,1) * lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries) - o2.in(o2_ids);
                                add(o2_rep.in(added) == 0);
                                
                                // Linking partition variables1 with lambda
                                if(model_type == "lambda_II"){
                                    if(num_partns1 > 1) {
                                        Constraint<> on_link_lambda1(pair.first+"_on_link_lambda1_II");
                                        on_link_lambda1 = lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef1.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef1.in_matrix(nb_entries,total_entries-nb_entries) - on1.in_ignore_ith(nb_entries_v1,nb_entries_v2,on_coef1.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef1.in_matrix(nb_entries,1);
                                        add(on_link_lambda1.in(indices(added,const_idx1)) <= 0);
                                    }
                                    if(num_partns2 > 1) {
                                        Constraint<> on_link_lambda2(pair.first+"_on_link_lambda2_II");
                                        on_link_lambda2 = lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef2.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef2.in_matrix(nb_entries,total_entries-nb_entries) - binvar2->in_ignore_ith(0,nb_entries_v1,on_coef2.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef2.in_matrix(nb_entries,1);
                                        add(on_link_lambda2.in(indices(added,const_idx2)) <= 0);
                                    }
                                }
                                else{
                                    if(num_partns1 > 1) {
                                        Constraint<> on_link_lambda1(pair.first+"_on_link_lambda1_III");
                                        on_link_lambda1 = lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef1.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef1.in_matrix(nb_entries,total_entries-nb_entries) - on1.in_ignore_ith(nb_entries_v1,nb_entries_v2,on_coef1.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef1.in_matrix(nb_entries,1);
                                        add(on_link_lambda1.in(indices(added,const_idx1)) <= 0);
                                    }
                                    if(num_partns2 > 1) {
                                        Constraint<> on_link_lambda2(pair.first+"_on_link_lambda2_III");
                                        on_link_lambda2 = lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef2.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef2.in_matrix(nb_entries,total_entries-nb_entries) - binvar2->in_ignore_ith(0,nb_entries_v1,on_coef2.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef2.in_matrix(nb_entries,1);
                                        add(on_link_lambda2.in(indices(added,const_idx2)) <= 0);
                                    }
                                }
                                // sum over lambda
                                Constraint<> lambdaSum(pair.first+"_lambdaSum");
                                lambdaSum = sum(lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries));
                                add(lambdaSum.in(added) == 1);
                            }
                        }
                        else{ //means v1 and v2 have not been partitioned before
                            if(name1 == name2){
                                 DebugOn("<<<<<<<<<< THIS IS SEEN LIFT -> DOUBLE -> NOT SEEN BINARIES -> SAME CORE VARS <<<<<<<<<<<" << endl);
                                
                                var<int> on1(name1+"_binary",0,1);
                                indices partns1("partns1");
                                for (int i = 0; i < num_partns1 ; ++i)
                                {
                                    partns1.add(name1+ "{" +to_string(i+1) + "}");
                                }
                                add(on1.in(union_ids(o1_ids_uq, o2_ids_uq),partns1));
                                indices partns2("partns2");
                                for (int i = 0; i < num_partns2 ; ++i)
                                {
                                    partns2.add(name2+ "{" +to_string(i+1) + "}");
                                }
                                
                                auto nb_entries_v1 = o1_ids.get_nb_entries();
                                auto nb_entries_v2 = o2_ids.get_nb_entries();
                                auto nb_entries = unique_ids.get_nb_entries();
                                
                                Constraint<> onSum1(o1._name+"_binarySum");
                                onSum1 = sum(on1.in_matrix(nb_entries_v1,1));
                                add(onSum1.in(union_ids(o1_ids_uq,o2_ids_uq)) == 1);
                                
                                if(model_type == "on/off"){//if on/off is chosen
                                    
                                    indices partns("partns");
                                    partns = indices(partns1,partns2);
                                    auto inst_partition = indices(added,partns);
                                    auto total_entries = inst_partition.get_nb_entries();
                                    
                                    if(binvar_ptr3 !=_vars_name.end()){ //means the combined binary variable has been used before
                                        auto binvar3 = static_pointer_cast<var<int>>(binvar_ptr3->second);
                                        param<int> lb3("lb3"), ub3("ub3");
                                        lb3.in(added,partns);
                                        ub3.in(added,partns);
                                        lb3.set_val(0), ub3.set_val(1);
                                        auto added3 = binvar3->add_bounds(lb3,ub3);
                                        reindex_vars();
                                        
                                        Constraint<> onLink1(pair.first+"_binaryLink1");
                                        onLink1 = on1.from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v1)) - binvar3->in(inst_partition);
                                        add(onLink1.in(inst_partition) >= 0);
                                        
                                        Constraint<> onLink2(pair.first+"_binaryLink2");
                                        onLink2 = on1.in_ignore_ith(nb_entries_v1,1,inst_partition.ignore_ith(0,nb_entries_v1)) - binvar3->in(inst_partition);
                                        add(onLink2.in(inst_partition) >= 0);
                                        
                                        Constraint<> onLink3(pair.first+"_binaryLink3");
                                        onLink3 = on1.from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v1)) + on1.in_ignore_ith(nb_entries_v1,1,inst_partition.ignore_ith(0,nb_entries_v1)) - 1 - binvar3->in(inst_partition);
                                        add(onLink3.in(inst_partition) <= 0);
                                        
                                        Constraint<> onSumComb(pair.first+"_binarySum");
                                        onSumComb = sum((binvar3->in(added3)).in_matrix(nb_entries,total_entries-nb_entries));
                                        add(onSumComb.in(added) == 1);
                                        
                                        add_on_off_McCormick_refined(pair.first, vlift->in(added), o1.in(o1_ids), o2.in(o2_ids), binvar3->in(added3));
                                    }
                                    
                                    else{ //means the combined binary variable has not been used before
                                        
                                        var<int> on(name1+name2+"_binary",0,1);
                                        add(on.in(inst_partition));
                                        
                                        Constraint<> onLink1(pair.first+"_binaryLink1");
                                        onLink1 = on1.from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v1)) - on;
                                        add(onLink1.in(inst_partition) >= 0);
                                        
                                        Constraint<> onLink2(pair.first+"_binaryLink2");
                                        onLink2 = on1.in_ignore_ith(nb_entries_v1,1,inst_partition.ignore_ith(0,nb_entries_v1)) - on;
                                        add(onLink2.in(inst_partition) >= 0);
                                        
                                        Constraint<> onLink3(pair.first+"_binaryLink3");
                                        onLink3 = on1.from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v1)) + on1.in_ignore_ith(nb_entries_v1,1,inst_partition.ignore_ith(0,nb_entries_v1)) - 1 - on;
                                        add(onLink3.in(inst_partition) <= 0);
                                        
                                        Constraint<> onSumComb(pair.first+"_binarySum");
                                        onSumComb = sum(on.in_matrix(nb_entries,total_entries-nb_entries));
                                        add(onSumComb.in(added) == 1);
                                        
                                        add_on_off_McCormick_refined(pair.first, vlift->in(added), o1.in(o1_ids), o2.in(o2_ids), on);
                                    }
                                }
                                
                                else{//means it is one of the lambda formulations
                                    
                                    //difference is this has one more partition index
                                    indices partns1_lambda("partns1_lambda");
                                    for (int i = 0; i < num_partns1+1; ++i)
                                    {
                                        partns1_lambda.add(name1+ "{" +to_string(i+1) + "}");
                                    }
                                    
                                    indices partns2_lambda("partns2_lambda");
                                    for (int i = 0; i < num_partns2+1; ++i)
                                    {
                                        partns2_lambda.add(name2+ "{" +to_string(i+1) + "}");
                                    }
                                    
                                    indices partns_lambda("partns_lambda");
                                    partns_lambda = indices(partns1_lambda,partns2_lambda);
                                    auto inst_partition_lambda = indices(added,partns_lambda);
                                    auto inst_partition_bounds1 = indices(added,partns1_lambda);
                                    auto inst_partition_bounds2 = indices(added,partns2_lambda);
                                    
                                    // Convex combination variables
                                    auto lambda_ptr = _vars_name.find(name1+name2+"_lambda");
                                    auto lambda = static_pointer_cast<var<double>>(lambda_ptr->second);
                                    param<double> lb_lambda("lb_lambda"), ub_lambda("ub_lambda");
                                    lb_lambda.in(added,partns_lambda);
                                    ub_lambda.in(added,partns_lambda);
                                    lb_lambda.set_val(0), ub_lambda.set_val(1);
                                    auto added_lambda = lambda->add_bounds(lb_lambda,ub_lambda);
                                    reindex_vars();
                                    
                                    /** Parameters */
                                    // Bounds on variable v1 & v2
                                    param<> bounds1(name1+"_bounds1");
                                    bounds1.in(inst_partition_bounds1);
                                    
                                    param<> bounds2(name2+"_bounds2");
                                    bounds2.in(inst_partition_bounds2);
                                    
                                    // Function values on the extreme points
                                    param<> EP(name1+name2+"_grid_values");
                                    EP.in(inst_partition_lambda);
                                    auto total_entries = inst_partition_lambda.get_nb_entries();
                                    
                                    size_t nb_ins = vlift->in(added).get_nb_inst();
                                    auto o1_global_lb = o1.get_lb();
                                    auto increment1 = (o1.get_ub() - o1_global_lb)/num_partns1;
                                    
                                    auto o2_global_lb = o2.get_lb();
                                    auto increment2 = (o2.get_ub() - o2_global_lb)/num_partns2;
                                    
                                    // fill bounds and function values
                                    for (int i=0 ; i<num_partns1+1; ++i) {
                                        auto bound_partn1 = o1_global_lb + increment1*i;
                                        bound_partn1.eval_all();
                                        auto bound_partn2 = o2_global_lb + increment2*i;
                                        bound_partn2.eval_all();
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift->get_id_inst(inst);
                                            auto cur_var_idx = added._keys->at(cur_var_id);
                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"}";
                                            bounds1.set_val(cur_idx,bound_partn1.eval(inst));
                                            bounds2.set_val(cur_idx,bound_partn2.eval(inst));
                                            for(int j=0; j<num_partns2+1; ++j){
                                                auto bound_partn2_temp = o2_global_lb + increment2*j;
                                                bound_partn2_temp.eval_all();
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"}";
                                                EP.set_val(cur_idx,(bound_partn1.eval(inst)*bound_partn2_temp.eval(inst)));
                                            }
                                        }
                                    }
                                    
                                    // Lambda coefficient matrix when linking with partition variables
                                    param<> lambda_coef1(name1+"_lambda_linking_coefficients1");
                                    param<> lambda_coef2(name2+"_lambda_linking_coefficients2");
                                    
                                    // Partition coefficient matrix when linking with lambda variables
                                    param<> on_coef1(name1+"_partition_linking_coefficients1");
                                    param<> on_coef2(name2+"_partition_linking_coefficients2");
                                    
                                    // create constraint indices
                                    indices const_idx1("const_idx1");
                                    indices const_idx2("const_idx2");
                                    
                                    if(model_type == "lambda_II"){
                                        
                                        //fill constraint indices
                                        for (int i=0; i<num_partns1+1; ++i){
                                            const_idx1.add(to_string(i+1));
                                        }
                                        
                                        //fill constraint indices
                                        for (int i=0; i<num_partns2+1; ++i){
                                            const_idx2.add(to_string(i+1));
                                        }
                                        
                                        // Lambda coefficient matrix when linking with partition variables
                                        lambda_coef1.in(indices(inst_partition_lambda, const_idx1));
                                        lambda_coef2.in(indices(inst_partition_lambda, const_idx2));
                                        
                                        // Lambda coefficient matrix when linking with partition variables
                                        on_coef1.in(indices(added, partns1, const_idx1));
                                        on_coef2.in(indices(added, partns2, const_idx2));
                                        
                                        // fill lambda_coef1 and lambda_coef2
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift->get_id_inst(inst);
                                            auto cur_var_idx = added._keys->at(cur_var_id);
                                            for (int i=0 ; i<num_partns1+1; ++i) {
                                                for (int j=0 ; j<num_partns2+1; ++j) {
                                                    string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                    lambda_coef1.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(j+1);
                                                    lambda_coef2.set_val(cur_idx,1);
                                                }
                                            }
                                        }
                                        
                                        // fill on_coef1 and on_coef2
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift->get_id_inst(inst);
                                            auto cur_var_idx = added._keys->at(cur_var_id);
                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                            on_coef1.set_val(cur_idx,1);
                                            cur_idx = cur_var_idx+","+name2+"{"+to_string(1)+"},"+to_string(1);
                                            on_coef2.set_val(cur_idx,1);
                                            for (int i=1 ; i<num_partns1; ++i) {
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(i)+"},"+to_string(i+1);
                                                on_coef1.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(i+1);
                                                on_coef1.set_val(cur_idx,1);
                                            }
                                            for (int i=1 ; i<num_partns2; ++i) {
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(i)+"},"+to_string(i+1);
                                                on_coef2.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"},"+to_string(i+1);
                                                on_coef2.set_val(cur_idx,1);
                                            }
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1)+"},"+to_string(num_partns1+1);
                                            on_coef1.set_val(cur_idx,1);
                                            cur_idx = cur_var_idx+","+name2+"{"+to_string(num_partns2)+"},"+to_string(num_partns2+1);
                                            on_coef2.set_val(cur_idx,1);
                                        }
                                    }
                                    
                                    
                                    else /*means model_type == "lambda_III" */{
                                        
                                        //fill constraint indices
                                        for (int i=0; i<(num_partns1-2)*2+2; ++i){
                                            const_idx1.add(to_string(i+1));
                                        }
                                        
                                        //fill constraint indices
                                        for (int i=0; i<(num_partns2-2)*2+2; ++i){
                                            const_idx2.add(to_string(i+1));
                                        }
                                        
                                        // Lambda coefficient matrix when linking with partition variables
                                        lambda_coef1.in(indices(inst_partition_lambda, const_idx1));
                                        lambda_coef2.in(indices(inst_partition_lambda, const_idx2));
                                        
                                        // Partition coefficient matrix when linking with lambda variables
                                        on_coef1.in(indices(added, partns1, const_idx1));
                                        on_coef2.in(indices(added, partns2, const_idx2));
                                        
                                        // fill lambda_coef1 and lambda_coef2
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift->get_id_inst(inst);
                                            auto cur_var_idx = added._keys->at(cur_var_id);
                                            for (int j=0; j<num_partns2+1; ++j) {
                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(1);
                                                lambda_coef1.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string((num_partns1-2)*2+2);
                                                lambda_coef1.set_val(cur_idx,1);
                                            }
                                            
                                            for (int i=1 ; i<(num_partns1-2)*2+1; i=i+2) {
                                                for (int j=(i-1)/2 + 2; j<num_partns1+1; ++j) {
                                                    for(int k=0; k<num_partns2+1; ++k){
                                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+name2+"{"+to_string(k+1)+"},"+to_string(i+1);
                                                        lambda_coef1.set_val(cur_idx,1);
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+name2+"{"+to_string(k+1)+"},"+to_string(i+2);
                                                        lambda_coef1.set_val(cur_idx,-1);
                                                    }
                                                }
                                            }
                                            
                                            for (int i=0; i<num_partns1+1; ++i) {
                                                string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(1)+"},"+to_string(1);
                                                lambda_coef2.set_val(cur_idx,1);
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(num_partns2+1)+"},"+to_string((num_partns2-2)*2+2);
                                                lambda_coef2.set_val(cur_idx,1);
                                            }
                                            
                                            for (int i=1 ; i<(num_partns2-2)*2+1; i=i+2) {
                                                for (int j=(i-1)/2 + 2; j<num_partns2+1; ++j) {
                                                    for(int k=0; k<num_partns1+1; ++k){
                                                        string cur_idx = cur_var_idx+","+name1+"{"+to_string(k+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                        lambda_coef2.set_val(cur_idx,1);
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(k+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                        lambda_coef2.set_val(cur_idx,-1);
                                                    }
                                                }
                                            }
                                        }
                                        
                                        
                                        
                                        // fill on_coef1 and on_coef2
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift->get_id_inst(inst);
                                            auto cur_var_idx = added._keys->at(cur_var_id);
                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                            on_coef1.set_val(cur_idx,1);
                                            
                                            for (int i=1; i<num_partns1; ++i) {
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(2);
                                                on_coef1.set_val(cur_idx, 1);
                                            }
                                            
                                            for (int i=2 ; i<(num_partns1-2)*2+2; i=i+2) {
                                                for (int j=i/2+1; j<num_partns1; ++j) {
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                    on_coef1.set_val(cur_idx,-1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                    on_coef1.set_val(cur_idx,1);
                                                }
                                            }
                                            
                                            cur_idx = cur_var_idx+","+name2+"{"+to_string(1)+"},"+to_string(1);
                                            on_coef2.set_val(cur_idx,1);
                                            for (int i=1; i<num_partns2; ++i) {
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"},"+to_string(2);
                                                on_coef2.set_val(cur_idx, 1);
                                            }
                                            
                                            for (int i=2 ; i<(num_partns2-2)*2+2; i=i+2) {
                                                for (int j=i/2+1; j<num_partns2; ++j) {
                                                    cur_idx = cur_var_idx+","+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                    on_coef2.set_val(cur_idx,-1);
                                                    cur_idx = cur_var_idx+","+name2+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                    on_coef2.set_val(cur_idx,1);
                                                }
                                            }
                                        }
                                    }
                                    
                                    
                                    /** Constraints */
                                    // Representation of the bilinear term with convex combination
                                    Constraint<> bln_rep(pair.first+"_bln_rep");
                                    bln_rep = EP.in_matrix(nb_entries,total_entries-nb_entries)*lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries) - vlift->in(added);
                                    add(bln_rep.in(added) == 0);
                                    
                                    // Representation of o1 with convex combination
                                    Constraint<> o1_rep(pair.first+"_o1_rep");
                                    o1_rep = bounds1.from_ith(0,inst_partition_lambda).in_matrix(nb_entries, 1) * lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries) - o1.in(o1_ids);
                                    add(o1_rep.in(added) == 0);
                                    
                                    // Representation of o2 with convex combination
                                    Constraint<> o2_rep(pair.first+"_o2_rep");
                                    o2_rep = bounds2.in_ignore_ith(nb_entries, 1, inst_partition_lambda).in_matrix(nb_entries,1) * lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries) - o2.in(o2_ids);
                                    add(o2_rep.in(added) == 0);
                                    
                                    // Linking partition variables1 with lambda
                                    if(model_type == "lambda_II"){
                                        Constraint<> on_link_lambda1(pair.first+"_on_link_lambda1_II");
                                        on_link_lambda1 = lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef1.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef1.in_matrix(nb_entries,total_entries-nb_entries) - on1.in_ignore_ith(nb_entries_v1,nb_entries_v2,on_coef1.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef1.in_matrix(nb_entries,1);
                                        add(on_link_lambda1.in(indices(added,const_idx1)) <= 0);
                                        
                                        Constraint<> on_link_lambda2(pair.first+"_on_link_lambda2_II");
                                        on_link_lambda2 = lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef2.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef2.in_matrix(nb_entries,total_entries-nb_entries) - on1.in_ignore_ith(0,nb_entries_v1,on_coef2.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef2.in_matrix(nb_entries,1);
                                        add(on_link_lambda2.in(indices(added,const_idx2)) <= 0);
                                    }
                                    else{
                                        Constraint<> on_link_lambda1(pair.first+"_on_link_lambda1_III");
                                        on_link_lambda1 = lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef1.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef1.in_matrix(nb_entries,total_entries-nb_entries) - on1.in_ignore_ith(nb_entries_v1,nb_entries_v2,on_coef1.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef1.in_matrix(nb_entries,1);
                                        add(on_link_lambda1.in(indices(added,const_idx1)) <= 0);
                                        
                                        Constraint<> on_link_lambda2(pair.first+"_on_link_lambda2_III");
                                        on_link_lambda2 = lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef2.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef2.in_matrix(nb_entries,total_entries-nb_entries) - on1.in_ignore_ith(0,nb_entries_v1,on_coef2.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef2.in_matrix(nb_entries,1);
                                        add(on_link_lambda2.in(indices(added,const_idx2)) <= 0);
                                    }
                                    // sum over lambda
                                    Constraint<> lambdaSum(pair.first+"_lambdaSum");
                                    lambdaSum = sum(lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries));
                                    add(lambdaSum.in(added) == 1);
                                }
                                
                            }
                            else{
                                DebugOn("<<<<<<<<<< THIS IS SEEN LIFT -> DOUBLE -> NOT SEEN BINARIES -> DIFF CORE VARS <<<<<<<<<<<" << endl);
                                
                                var<int> on1(name1+"_binary",0,1);
                                indices partns1("partns1");
                                for (int i = 0; i < num_partns1 ; ++i)
                                {
                                    partns1.add(name1+ "{" + to_string(i+1) + "}");
                                }
                                add(on1.in(o1_ids_uq,partns1));
                                
                                var<int> on2(name2+"_binary",0,1);
                                indices partns2("partns2");
                                for (int i = 0; i < num_partns2 ; ++i)
                                {
                                    partns2.add(name2+ "{" + to_string(i+1) + "}");
                                }
                                add(on2.in(o2_ids_uq,partns2));
                                
                                auto nb_entries_v1 = o1_ids.get_nb_entries();
                                auto nb_entries_v2 = o2_ids.get_nb_entries();
                                auto nb_entries = unique_ids.get_nb_entries();
                                
                                Constraint<> onSum1(o1._name+"_binarySum");
                                onSum1 = sum(on1.in_matrix(nb_entries_v1,1));
                                add(onSum1.in(o1_ids_uq) == 1);
                                
                                Constraint<> onSum2(o2._name+"_binarySum");
                                onSum2 = sum(on2.in_matrix(nb_entries_v2,1));
                                add(onSum2.in(o2_ids_uq) == 1);
                                
                                if(model_type == "on/off"){//if on/off is chosen
                                    
                                    indices partns("partns");
                                    partns = indices(partns1,partns2);
                                    auto inst_partition = indices(added,partns);
                                    auto total_entries = inst_partition.get_nb_entries();
                                    
                                    if(binvar_ptr3 !=_vars_name.end()){ //means the combined binary variable has been used before
                                        auto binvar3 = static_pointer_cast<var<int>>(binvar_ptr3->second);
                                        param<int> lb3("lb3"), ub3("ub3");
                                        lb3.in(added,partns);
                                        ub3.in(added,partns);
                                        lb3.set_val(0), ub3.set_val(1);
                                        auto added3 = binvar3->add_bounds(lb3,ub3);
                                        reindex_vars();
                                        
                                        Constraint<> onLink1(pair.first+"_binaryLink1");
                                        onLink1 = on1.from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) - binvar3->in(inst_partition);
                                        add(onLink1.in(inst_partition) >= 0);
                                        
                                        Constraint<> onLink2(pair.first+"_binaryLink2");
                                        onLink2 = on2.in_ignore_ith(nb_entries_v2,1,inst_partition.ignore_ith(0,nb_entries_v1)) - binvar3->in(inst_partition);
                                        add(onLink2.in(inst_partition) >= 0);
                                        
                                        Constraint<> onLink3(pair.first+"_binaryLink3");
                                        onLink3 = on1.from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) + on2.in_ignore_ith(nb_entries_v2,1,inst_partition.ignore_ith(0,nb_entries_v1)) - 1 - binvar3->in(inst_partition);
                                        add(onLink3.in(inst_partition) <= 0);
                                        
                                        Constraint<> onSumComb(pair.first+"_binarySum");
                                        onSumComb = sum((binvar3->in(added3)).in_matrix(nb_entries,total_entries-nb_entries));
                                        add(onSumComb.in(added) == 1);
                                        
                                        add_on_off_McCormick_refined(pair.first, vlift->in(added), o1.in(o1_ids), o2.in(o2_ids), binvar3->in(added3));
                                    }
                                    else{ //means the combined binary variable has not been used before
                                        var<int> on(name1+name2+"_binary",0,1);
                                        add(on.in(inst_partition));
                                        
                                        Constraint<> onLink1(pair.first+"_binaryLink1");
                                        onLink1 = on1.from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) - on;
                                        add(onLink1.in(inst_partition) >= 0);
                                        
                                        Constraint<> onLink2(pair.first+"_binaryLink2");
                                        onLink2 = on2.in_ignore_ith(nb_entries_v2,1,inst_partition.ignore_ith(0,nb_entries_v1)) - on;
                                        add(onLink2.in(inst_partition) >= 0);
                                        
                                        Constraint<> onLink3(pair.first+"_binaryLink3");
                                        onLink3 = on1.from_ith(0,inst_partition.ignore_ith(nb_entries_v1, nb_entries_v2)) + on2.in_ignore_ith(nb_entries_v2,1,inst_partition.ignore_ith(0,nb_entries_v1)) - 1 - on;
                                        add(onLink3.in(inst_partition) <= 0);
                                        
                                        Constraint<> onSumComb(pair.first+"_binarySum");
                                        onSumComb = sum(on.in_matrix(nb_entries,total_entries-nb_entries));
                                        add(onSumComb.in(added) == 1);
                                        
                                        add_on_off_McCormick_refined(pair.first, vlift->in(added), o1.in(o1_ids), o2.in(o2_ids), on);
                                    }
                                }
                                
                                else{//means it is one of the lambda formulations
                                    
                                    //difference is this has one more partition index
                                    indices partns1_lambda("partns1_lambda");
                                    for (int i = 0; i < num_partns1+1; ++i)
                                    {
                                        partns1_lambda.add(name1+ "{" +to_string(i+1) + "}");
                                    }
                                    
                                    indices partns2_lambda("partns2_lambda");
                                    for (int i = 0; i < num_partns2+1; ++i)
                                    {
                                        partns2_lambda.add(name2+ "{" +to_string(i+1) + "}");
                                    }
                                    
                                    indices partns_lambda("partns_lambda");
                                    partns_lambda = indices(partns1_lambda,partns2_lambda);
                                    auto inst_partition_lambda = indices(added,partns_lambda);
                                    auto inst_partition_bounds1 = indices(added,partns1_lambda);
                                    auto inst_partition_bounds2 = indices(added,partns2_lambda);
                                    
                                    // Convex combination variables
                                    auto lambda_ptr = _vars_name.find(name1+name2+"_lambda");
                                    auto lambda = static_pointer_cast<var<double>>(lambda_ptr->second);
                                    param<double> lb_lambda("lb_lambda"), ub_lambda("ub_lambda");
                                    lb_lambda.in(added,partns_lambda);
                                    ub_lambda.in(added,partns_lambda);
                                    lb_lambda.set_val(0), ub_lambda.set_val(1);
                                    auto added_lambda = lambda->add_bounds(lb_lambda,ub_lambda);
                                    reindex_vars();
                                    
                                    /** Parameters */
                                    // Bounds on variable v1 & v2
                                    param<> bounds1(name1+"_bounds1");
                                    bounds1.in(inst_partition_bounds1);
                                    
                                    param<> bounds2(name2+"_bounds2");
                                    bounds2.in(inst_partition_bounds2);
                                    
                                    // Function values on the extreme points
                                    param<> EP(name1+name2+"_grid_values");
                                    EP.in(inst_partition_lambda);
                                    auto total_entries = inst_partition_lambda.get_nb_entries();
                                    
                                    size_t nb_ins = vlift->in(added).get_nb_inst();
                                    auto o1_global_lb = o1.get_lb();
                                    auto increment1 = (o1.get_ub() - o1_global_lb)/num_partns1;
                                    
                                    auto o2_global_lb = o2.get_lb();
                                    auto increment2 = (o2.get_ub() - o2_global_lb)/num_partns2;
                                    
                                    // fill bounds1 and function values
                                    for (int i=0 ; i<num_partns1+1; ++i) {
                                        auto bound_partn1 = o1_global_lb + increment1*i;
                                        bound_partn1.eval_all();
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift->get_id_inst(inst);
                                            auto cur_var_idx = added._keys->at(cur_var_id);
                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"}";
                                            bounds1.set_val(cur_idx,bound_partn1.eval(inst));
                                            for(int j=0; j<num_partns2+1; ++j){
                                                auto bound_partn2 = o2_global_lb + increment2*j;
                                                bound_partn2.eval_all();
                                                cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"}";
                                                EP.set_val(cur_idx,(bound_partn1.eval(inst)*bound_partn2.eval(inst)));
                                            }
                                        }
                                    }
                                    // fill bounds2
                                    for (int i=0 ; i<num_partns2+1; ++i) {
                                        auto bound_partn2 = o2_global_lb + increment2*i;
                                        bound_partn2.eval_all();
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift->get_id_inst(inst);
                                            auto cur_var_idx = added._keys->at(cur_var_id);
                                            string cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"}";
                                            bounds2.set_val(cur_idx,bound_partn2.eval(inst));
                                        }
                                    }
                                    
                                    // Lambda coefficient matrix when linking with partition variables
                                    param<> lambda_coef1(name1+"_lambda_linking_coefficients1");
                                    param<> lambda_coef2(name2+"_lambda_linking_coefficients2");
                                    
                                    // Partition coefficient matrix when linking with lambda variables
                                    param<> on_coef1(name1+"_partition_linking_coefficients1");
                                    param<> on_coef2(name2+"_partition_linking_coefficients2");
                                    
                                    // create constraint indices
                                    indices const_idx1("const_idx1");
                                    indices const_idx2("const_idx2");
                                    
                                    if(model_type == "lambda_II"){
                                        
                                        //fill constraint indices
                                        for (int i=0; i<num_partns1+1; ++i){
                                            const_idx1.add(to_string(i+1));
                                        }
                                        
                                        //fill constraint indices
                                        for (int i=0; i<num_partns2+1; ++i){
                                            const_idx2.add(to_string(i+1));
                                        }
                                        
                                        // Lambda coefficient matrix when linking with partition variables
                                        if(num_partns1 > 1) lambda_coef1.in(indices(inst_partition_lambda, const_idx1));
                                        if(num_partns2 > 1) lambda_coef2.in(indices(inst_partition_lambda, const_idx2));
                                        
                                        // Lambda coefficient matrix when linking with partition variables
                                        if(num_partns1 > 1) on_coef1.in(indices(added, partns1, const_idx1));
                                        if(num_partns2 > 1) on_coef2.in(indices(added, partns2, const_idx2));
                                        
                                        // fill lambda_coef1 and lambda_coef2
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift->get_id_inst(inst);
                                            auto cur_var_idx = added._keys->at(cur_var_id);
                                            for (int i=0 ; i<num_partns1+1; ++i) {
                                                for (int j=0 ; j<num_partns2+1; ++j) {
                                                    string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                    if(num_partns1 > 1) lambda_coef1.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(j+1);
                                                    if(num_partns2 > 1) lambda_coef2.set_val(cur_idx,1);
                                                }
                                            }
                                        }
                                        
                                        // fill on_coef1 and on_coef2
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift->get_id_inst(inst);
                                            auto cur_var_idx = added._keys->at(cur_var_id);
                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                            if(num_partns1 > 1) on_coef1.set_val(cur_idx,1);
                                            cur_idx = cur_var_idx+","+name2+"{"+to_string(1)+"},"+to_string(1);
                                            if(num_partns2 > 1) on_coef2.set_val(cur_idx,1);
                                            if(num_partns1 > 1) {
                                                for (int i=1 ; i<num_partns1; ++i) {
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i)+"},"+to_string(i+1);
                                                    on_coef1.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(i+1);
                                                    on_coef1.set_val(cur_idx,1);
                                                }
                                            }
                                            if(num_partns2 > 1) {
                                                for (int i=1 ; i<num_partns2; ++i) {
                                                    cur_idx = cur_var_idx+","+name2+"{"+to_string(i)+"},"+to_string(i+1);
                                                    on_coef2.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"},"+to_string(i+1);
                                                    on_coef2.set_val(cur_idx,1);
                                                }
                                            }
                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1)+"},"+to_string(num_partns1+1);
                                            if(num_partns1 > 1) on_coef1.set_val(cur_idx,1);
                                            cur_idx = cur_var_idx+","+name2+"{"+to_string(num_partns2)+"},"+to_string(num_partns2+1);
                                            if(num_partns2 > 1) on_coef2.set_val(cur_idx,1);
                                        }
                                    }
                                    
                                    
                                    else /*means model_type == "lambda_III" */{
                                        
                                        //fill constraint indices
                                        for (int i=0; i<(num_partns1-2)*2+2; ++i){
                                            const_idx1.add(to_string(i+1));
                                        }
                                        
                                        //fill constraint indices
                                        for (int i=0; i<(num_partns2-2)*2+2; ++i){
                                            const_idx2.add(to_string(i+1));
                                        }
                                        
                                        // Lambda coefficient matrix when linking with partition variables
                                        if(num_partns1 > 1) lambda_coef1.in(indices(inst_partition_lambda, const_idx1));
                                        if(num_partns2 > 1) lambda_coef2.in(indices(inst_partition_lambda, const_idx2));
                                        
                                        // Partition coefficient matrix when linking with lambda variables
                                        if(num_partns1 > 1) on_coef1.in(indices(added, partns1, const_idx1));
                                        if(num_partns2 > 1) on_coef2.in(indices(added, partns2, const_idx2));
                                        
                                        // fill lambda_coef1 and lambda_coef2
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift->get_id_inst(inst);
                                            auto cur_var_idx = added._keys->at(cur_var_id);
                                            if(num_partns1 > 1) {
                                                for (int j=0; j<num_partns2+1; ++j) {
                                                    string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(1);
                                                    lambda_coef1.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(num_partns1+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string((num_partns1-2)*2+2);
                                                    lambda_coef1.set_val(cur_idx,1);
                                                }
                                                
                                                for (int i=1 ; i<(num_partns1-2)*2+1; i=i+2) {
                                                    for (int j=(i-1)/2 + 2; j<num_partns1+1; ++j) {
                                                        for(int k=0; k<num_partns2+1; ++k){
                                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+name2+"{"+to_string(k+1)+"},"+to_string(i+1);
                                                            lambda_coef1.set_val(cur_idx,1);
                                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+name2+"{"+to_string(k+1)+"},"+to_string(i+2);
                                                            lambda_coef1.set_val(cur_idx,-1);
                                                        }
                                                    }
                                                }
                                            }
                                            if(num_partns2 > 1) {
                                                for (int i=0; i<num_partns1+1; ++i) {
                                                    string cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(1)+"},"+to_string(1);
                                                    lambda_coef2.set_val(cur_idx,1);
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+name2+"{"+to_string(num_partns2+1)+"},"+to_string((num_partns2-2)*2+2);
                                                    lambda_coef2.set_val(cur_idx,1);
                                                }
                                                
                                                for (int i=1 ; i<(num_partns2-2)*2+1; i=i+2) {
                                                    for (int j=(i-1)/2 + 2; j<num_partns2+1; ++j) {
                                                        for(int k=0; k<num_partns1+1; ++k){
                                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(k+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                            lambda_coef2.set_val(cur_idx,1);
                                                            cur_idx = cur_var_idx+","+name1+"{"+to_string(k+1)+"},"+name2+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                            lambda_coef2.set_val(cur_idx,-1);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        
                                        
                                        
                                        // fill on_coef1 and on_coef2
                                        for (size_t inst = 0; inst< nb_ins; inst++){
                                            auto cur_var_id = vlift->get_id_inst(inst);
                                            auto cur_var_idx = added._keys->at(cur_var_id);
                                            string cur_idx = cur_var_idx+","+name1+"{"+to_string(1)+"},"+to_string(1);
                                            if(num_partns1 > 1) {
                                                on_coef1.set_val(cur_idx,1);
                                                
                                                for (int i=1; i<num_partns1; ++i) {
                                                    cur_idx = cur_var_idx+","+name1+"{"+to_string(i+1)+"},"+to_string(2);
                                                    on_coef1.set_val(cur_idx, 1);
                                                }
                                                
                                                for (int i=2 ; i<(num_partns1-2)*2+2; i=i+2) {
                                                    for (int j=i/2+1; j<num_partns1; ++j) {
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                        on_coef1.set_val(cur_idx,-1);
                                                        cur_idx = cur_var_idx+","+name1+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                        on_coef1.set_val(cur_idx,1);
                                                    }
                                                }
                                            }
                                            if(num_partns2 > 1) {
                                                cur_idx = cur_var_idx+","+name2+"{"+to_string(1)+"},"+to_string(1);
                                                on_coef2.set_val(cur_idx,1);
                                                for (int i=1; i<num_partns2; ++i) {
                                                    cur_idx = cur_var_idx+","+name2+"{"+to_string(i+1)+"},"+to_string(2);
                                                    on_coef2.set_val(cur_idx, 1);
                                                }
                                                
                                                for (int i=2 ; i<(num_partns2-2)*2+2; i=i+2) {
                                                    for (int j=i/2+1; j<num_partns2; ++j) {
                                                        cur_idx = cur_var_idx+","+name2+"{"+to_string(j+1)+"},"+to_string(i+1);
                                                        on_coef2.set_val(cur_idx,-1);
                                                        cur_idx = cur_var_idx+","+name2+"{"+to_string(j+1)+"},"+to_string(i+2);
                                                        on_coef2.set_val(cur_idx,1);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    
                                    
                                    /** Constraints */
                                    // Representation of the bilinear term with convex combination
                                    Constraint<> bln_rep(pair.first+"_bln_rep");
                                    bln_rep = EP.in_matrix(nb_entries,total_entries-nb_entries)*lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries) - vlift->in(added);
                                    add(bln_rep.in(added) == 0);
                                    
                                    // Representation of o1 with convex combination
                                    Constraint<> o1_rep(pair.first+"_o1_rep");
                                    o1_rep = bounds1.from_ith(0,inst_partition_lambda).in_matrix(nb_entries, 1) * lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries) - o1.in(o1_ids);
                                    add(o1_rep.in(added) == 0);
                                    
                                    // Representation of o2 with convex combination
                                    Constraint<> o2_rep(pair.first+"_o2_rep");
                                    o2_rep = bounds2.in_ignore_ith(nb_entries, 1, inst_partition_lambda).in_matrix(nb_entries,1) * lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries) - o2.in(o2_ids);
                                    add(o2_rep.in(added) == 0);
                                    
                                    // Linking partition variables1 with lambda
                                    if(model_type == "lambda_II"){
                                        if(num_partns1 > 1) {
                                            Constraint<> on_link_lambda1(pair.first+"_on_link_lambda1_II");
                                            on_link_lambda1 = lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef1.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef1.in_matrix(nb_entries,total_entries-nb_entries) - on1.in_ignore_ith(nb_entries_v1,nb_entries_v2,on_coef1.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef1.in_matrix(nb_entries,1);
                                            add(on_link_lambda1.in(indices(added,const_idx1)) <= 0);
                                        }
                                        if(num_partns2 > 1) {
                                            Constraint<> on_link_lambda2(pair.first+"_on_link_lambda2_II");
                                            on_link_lambda2 = lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef2.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef2.in_matrix(nb_entries,total_entries-nb_entries) - on2.in_ignore_ith(0,nb_entries_v1,on_coef2.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef2.in_matrix(nb_entries,1);
                                            add(on_link_lambda2.in(indices(added,const_idx2)) <= 0);
                                        }
                                    }
                                    else{
                                        if(num_partns1 > 1) {
                                            Constraint<> on_link_lambda1(pair.first+"_on_link_lambda1_III");
                                            on_link_lambda1 = lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef1.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef1.in_matrix(nb_entries,total_entries-nb_entries) - on1.in_ignore_ith(nb_entries_v1,nb_entries_v2,on_coef1.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef1.in_matrix(nb_entries,1);
                                            add(on_link_lambda1.in(indices(added,const_idx1)) <= 0);
                                        }
                                        if(num_partns2 > 1) {
                                            Constraint<> on_link_lambda2(pair.first+"_on_link_lambda2_III");
                                            on_link_lambda2 = lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries).from_ith(0,lambda_coef2.get_matrix_ids(nb_entries,total_entries-nb_entries))*lambda_coef2.in_matrix(nb_entries,total_entries-nb_entries) - on2.in_ignore_ith(0,nb_entries_v1,on_coef2.get_matrix_ids(nb_entries,1).from_ith(0,nb_entries+1)) * on_coef2.in_matrix(nb_entries,1);
                                            add(on_link_lambda2.in(indices(added,const_idx2)) <= 0);
                                        }
                                    }
                                    // sum over lambda
                                    Constraint<> lambdaSum(pair.first+"_lambdaSum");
                                    lambdaSum = sum(lambda->in(added_lambda).in_matrix(nb_entries,total_entries-nb_entries));
                                    add(lambdaSum.in(added) == 1);
                                }
                            }
                        }
                    }
#endif
                }
                else {
                    add_McCormick(pair.first, vlift->in(added), o1.in(o1_ids), o2.in(o2_ids));
                }
            }
        }
        
        lifted.insert(lt);
    }
    for (auto &pair:*c._pterms) {
        auto term = pair.second;
        lterm lt;
        lt._sign = term._sign;
        if (term._coef->is_function()) {
            auto coef = *static_pointer_cast<func<type>>(term._coef);
            lt._coef = func<type>(coef).copy();
        }
        else if(term._coef->is_param()) {
            auto coef = *static_pointer_cast<param<type>>(term._coef);
            lt._coef = param<type>(coef).copy();
        }
        else if(term._coef->is_number()) {
            auto coef = *static_pointer_cast<constant<type>>(term._coef);
            lt._coef = constant<type>(coef).copy();
        }
        func<type> prod = 1;
        string prod_name = "Lift(";
        auto list = pair.second._l;
        for (auto &ppi: *list) {
            auto p = ppi.first;
            auto orig_var = *static_pointer_cast<var<type>>(p);
            if(ppi.second>1){
                prod_name += orig_var.get_name(true,true)+"("+orig_var._indices->get_name()+")^"+to_string(ppi.second);
                //TODO Lift univarite power function
            }
            else{
                prod_name += orig_var.get_name(true,true)+"("+orig_var._indices->get_name()+")";
            }
            prod *= pow(orig_var,ppi.second);
        }
        prod_name += ")";
        
        auto ids = *c._indices;
        param<type> lb("lb"), ub("ub");
        lb.in(ids);ub.in(ids);
        lb.set_val(prod._range->first);
        ub.set_val(prod._range->second);
        var<type> vlift(prod_name, lb, ub);
        auto it = _vars_name.find(prod_name);
        if(it==_vars_name.end()){
            vlift._lift=true;
            add(vlift.in(ids));
            lt._p = make_shared<var<type>>(vlift);
        }
        else {
            vlift = *static_pointer_cast<var<type>>(it->second);
            lt._p = make_shared<var<type>>(vlift);
        }
        lifted.insert(lt);
    }

    lifted._range = c._range;
    lifted._all_convexity = linear_;
    lifted._all_sign = c._all_sign;
    lifted._ftype = lin_;
    lifted._ctype = c._ctype;
    lifted._indices = c._indices;
    lifted._dim[0] = c._dim[0];
    lifted._dim[1] = c._dim[1];
    return lifted;
}


    template <typename type>
    template<typename T,
    typename std::enable_if<is_same<T,double>::value>::type*>
    std::tuple<bool,int,double,double,double,bool> Model<type>::run_obbt(double max_time, unsigned max_iter, const pair<bool,double>& upper_bound, unsigned precision, shared_ptr<Model<type>> ub_model, shared_ptr<Model<type>> nonlin_model, bool nonlin) {
        std::tuple<bool,int,double, double, double, bool> res;

#ifdef USE_MPI
    int worker_id, nb_workers;
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
    auto err_size = MPI_Comm_size(MPI_COMM_WORLD, &nb_workers);
#endif
    int nb_threads = 12;
    int nb_total_threads = nb_threads; /** Used when MPI is ON to multipply with the number of workers */
#ifdef USE_MPI
    nb_total_threads *= nb_workers;
#endif
        
        vector<shared_ptr<Model<>>> batch_models;
        map<string, bool> fixed_point;
        map<string, double> interval_original, interval_new, ub_original, lb_original;
        string var_key,var_key_k,key_lb,key_ub, key_lb_k, key_ub_k;
        string vname;
        string mname, mkname, vkname, keyk, dirk;
        string var_vp_key, vp_key_lb, vp_key_ub;
        string dir_array[2]={"LB", "UB"};
        var<> vark, vk, v, var_ub;
        
        double boundk1, objk, left, right, mid, temp, tempa, vk_ub, vk_lb, vk_xb;
        
        bool infeasible=false;
        bool terminate=false;
        bool break_flag=false, time_limit = false, lifted_var=false, close=false;
        bool xb_true=true;
        double sum=0, avg=0, num_var=0.0;
        const double rel_tol=1e-2, abs_tol=1e6, fixed_tol_abs=1e-3, fixed_tol_rel=1e-3, zero_tol=1e-6, range_tol=1e-2, zero_val=1e-6;
        int gap_count_int=1, iter=0;
        double ub_vp, lb_vp, ub_vp_new, lb_vp_new;
        SolverType solv_type = ipopt;
        const double tol = 1e-8;
          int output = 0;
        
        
        
        double solver_time =0, solver_time_end, gapnl,gap, solver_time_start = get_wall_time();
       
        shared_ptr<map<string,size_t>> p_map;
//        this->reindex();
        solver<> SDPLB2(*this,solv_type);

        SDPLB2.run(output = 0, 1e-9, "ma27");
        double lower_bound_init=-999, lower_bound;

//        this->print();

        //Check if gap is already not zero at root node
        if(this->_status==0)
        {
            lower_bound_init=this->get_obj_val()*upper_bound.second;
            gapnl=(upper_bound.second-lower_bound_init)/(upper_bound.second)*100;
            DebugOn("Initial gap "<<gapnl<<endl);
            double lower_bound=lower_bound_init;
            

            
          if ((upper_bound.second-lower_bound)>=abs_tol || (upper_bound.second-lower_bound)/(upper_bound.second+zero_tol)>=rel_tol)
        {
            terminate=false;
            for(auto &it:this->_vars_name)
            {
                string vname=it.first;
                v=this->get_var<double>(vname);
                auto v_keys=v.get_keys();
                auto v_key_map=v.get_keys_map();
                for(auto &key: *v_keys)
                {
                    var_key = vname+"|"+ key;
                    key_lb= var_key +"|LB";
                    key_ub= var_key +"|UB";
                    //Do not do OBBT on lifted variables
                    if(v._lift){
                        fixed_point[key_lb]=true;
                        fixed_point[key_ub]=true;
                        DebugOff("Skipping OBBT for "<<vname<<"\t"<<key<<endl);
                    }
                    else{
                        fixed_point[key_lb]=false;
                        fixed_point[key_ub]=false;
                    }
                    auto key_pos=v_key_map->at(key);
                    
                    if(v._off[key_pos]==true)
                    {
                        fixed_point[key_lb]=true;
                        fixed_point[key_ub]=true;
                        DebugOn("Skipping OBBT for "<<vname<<"\t"<<key<<endl);
                    }
                    
                    interval_original[var_key]=v.get_ub(key)-v.get_lb(key);
                    ub_original[var_key]=v.get_ub(key);
                    lb_original[var_key]=v.get_lb(key);
                    interval_new[var_key]=v.get_ub(key)-v.get_lb(key);
                    
                }
                
            }
        solver_time= get_wall_time()-solver_time_start;
        while(solver_time<=max_time && !terminate && iter<=max_iter)
        {
            iter++;
            terminate=true;
            for (auto it=this->_vars_name.begin(); it!=this->_vars_name.end(); it++)
            {
                vname=it->first;
                v = this->get_var<double>(vname);
                auto v_keys=v.get_keys();
                for(auto it_key=v.get_keys()->begin(); it_key!=v.get_keys()->end(); it_key++)
                {
                    
                    auto key = *it_key;
                    solver_time_end=get_wall_time();
                    solver_time= solver_time_end-solver_time_start;
                    if(solver_time>=max_time)
                        
                    {
                        break_flag=true;
                        time_limit = true;
                        break;
                    }
                    var_key=vname+"|"+ key;
                    key_lb= var_key +"|LB";
                    key_ub= var_key +"|UB";
                    interval_new[var_key]=v.get_ub(key)-v.get_lb(key);
                    if(std::abs(v.get_ub(key)-v.get_lb(key))<=range_tol)
                    {
                        fixed_point[key_lb]=true;
                        fixed_point[key_ub]=true;
                        
                    }
                    //Either if not fixed point, or if at the last key of the last variable
                    if(fixed_point[key_lb]==false || fixed_point[key_ub]==false || (next(it)==this->_vars_name.end() && next(it_key)==v.get_keys()->end()))
                    {
                        //Loop on directions, upper bound and lower bound
                        for(auto &dir: dir_array)
                        {
                           
                            auto modelk = this->copy();
                            
                            mname=vname+"|"+key+"|"+dir;
                            modelk->set_name(mname);
                            
                            vark=modelk->template get_var<T>(vname);
                            vark.initialize_midpoint();
                            if(dir=="LB")
                            {
                                modelk->min(vark(key));
                            }
                            else
                            {
                                modelk->max(vark(key));
                                
                            }
                            modelk->reindex();
                            if(fixed_point[mname]==false){
                                batch_models.push_back(modelk);
                            }
                            //When batch models has reached size of nb_threads or when at the last key of last avriable
                            if (batch_models.size()==nb_total_threads || (next(it)==this->_vars_name.end() && next(it_key)==v.get_keys()->end() && dir=="UB"))
                            {
                                double batch_time_start = get_wall_time();
#ifdef USE_MPI
                                run_MPI(batch_models,solv_type,tol,nb_threads,"ma27",2000,2000, false,true);

#else
                              run_parallel(batch_models,solv_type,tol,nb_threads, "ma27", 2000);
                              //  run_parallel(batch_models,cplex,1e-6,nb_threads);
#endif
                                double batch_time_end = get_wall_time();
                                auto batch_time = batch_time_end - batch_time_start;
                                DebugOn("Done running batch models, solve time = " << to_string(batch_time) << endl);
                                for (auto model:batch_models)
                                {
                                    //Update bounds only of the model status is solved to optimal                                }
                                    if(model->_status==0)
                                    {
                                        mkname=model->get_name();
                                        std::size_t pos = mkname.find("|");
                                        vkname.assign(mkname, 0, pos);
                                        mkname=mkname.substr(pos+1);
                                        pos=mkname.find("|");
                                        keyk.assign(mkname, 0, pos);
                                        dirk=mkname.substr(pos+1);
                                        vk=this->get_var<T>(vkname);
                                        var_key_k=vkname+"|"+keyk;
                                        objk=model->get_obj_val();
                                    
                                        if(dirk=="LB")
                                        {
                                            boundk1=vk.get_lb(keyk);
                                            //Uncertainty in objk=obk+-solver_tolerance, here we choose lowest possible value in uncertainty interval
                                            objk=std::max(objk-range_tol, boundk1);
                                        }
                                        else
                                        {
                                            boundk1=vk.get_ub(keyk);
                                            //Uncertainty in objk=obk+-solver_tolerance, here we choose highest possible value in uncertainty interval
                                            objk=std::min(objk+range_tol, boundk1);
                                            
                                        }
                                        if(std::abs(boundk1-objk) <= fixed_tol_abs || std::abs((boundk1-objk)/(boundk1+zero_tol))<=fixed_tol_rel)
                                        {//do not close intervals to OBBT before finishing at least one full iteration over all variables
                                            if(iter>1)
                                                fixed_point[model->get_name()]=true;
                                            
                                        }
                                        else
                                        {
                                            if(dirk=="LB"){
                                                vk.set_lb(keyk, objk);/* IN MPI this needs to be broadcasted back to the other workers */
                                            }
                                            else{
                                                vk.set_ub(keyk, objk);
                                            }
                                            //If crossover in bounds,just exchange them
                                            if(vk.get_ub(keyk)<vk.get_lb(keyk))
                                            {
                                                fixed_point[var_key_k+"|LB"]=true;
                                                fixed_point[var_key_k+"|UB"]=true;
                                                temp=vk.get_ub(keyk);
                                                tempa=vk.get_lb(keyk);
                                                vk.set_ub(keyk, tempa);
                                                vk.set_lb(keyk, temp);
                                                
                                            }
                                            else if(!vk._lift){
                                                fixed_point[model->get_name()]=false;
                                                terminate=false;
                                            }
//                                            for (auto &vp: model->_vars) {
//                                                auto nb_inst = vp.second->get_dim();
//                                                auto ldv= vp.second->_l_dual;
//                                                auto udv=vp.second->_u_dual;
//                                                auto vpkeys=vp.second->get_keys();
//                                                bool in_orig_model=false;
//                                                if(ub_model->_vars_name.find(vk._name)!=ub_model->_vars_name.end())
//                                                {
//                                                     var_ub=ub_model->get_var<T>(vk._name);
//                                                    in_orig_model=true;
//                                                }
//
//                                                for(auto vpiter=0;vpiter<nb_inst;vpiter++)
//                                                {
//                                                    var_vp_key = vp.second->_name+"|"+ (*vpkeys)[vpiter];
//                                                    vp_key_lb= var_vp_key +"|LB";
//                                                    vp_key_ub= var_vp_key +"|UB";
//
//                                                    if(fixed_point[vp_key_lb]==false || fixed_point[vp_key_ub]==false){
//
//                                                    if(udv[vpiter] >= 1E-3)
//                                                    {
//                                                        //DebugOn("Ith variable "<<vk._name<<endl);
//                                                        //DebugOn("Ith variable index "<<keyk<<endl);
//
//
//
//                                                       var<> var_vp=this->get_var<T>(vp.second->_name);
//
//                                                        lb_vp=vp.second->get_double_lb(vpiter);
//                                                        ub_vp=vp.second->get_double_ub(vpiter);
//
//                                                        vk_ub=vk.get_ub(keyk);
//                                                        vk_lb=vk.get_lb(keyk);
//                                                        if(in_orig_model)
//                                                        {
//                                                            vk_xb=var_ub.eval(keyk);
//                                                            if(dirk=="UB"){
//                                                                vk_lb=vk_xb;
//                                                            }
//                                                            if(dirk=="LB"){
//                                                                vk_ub=vk_xb;
//                                                            }
//
//                                                        }
//
//                                                        lb_vp_new=ub_vp-(vk_ub-vk_lb)/udv[vpiter];
//
//
//
//
//                                                        if((lb_vp_new-lb_vp)>=range_tol){
//                                                           //DebugOn("Lower Bound update"<<endl);
////                                                            DebugOn("Ith variable "<<vk._name<<endl);
////                                                            DebugOn("Ith variable index "<<keyk<<endl);
////                                                            DebugOn("Jth variable "<<var_vp._name<<endl);
////                                                            DebugOn("Jth variable index "<<(*vpkeys)[vpiter]<<endl);
////                                                            DebugOn(lb_vp<<"\t"<< lb_vp_new<<endl);
//                                                        var_vp.set_lb((*vpkeys)[vpiter], lb_vp_new);
//                                                        }
//
//                                                    }
//                                                    if(ldv[vpiter] >= 1E-3){
//                                                        lb_vp=vp.second->get_double_lb(vpiter);
//                                                        ub_vp=vp.second->get_double_ub(vpiter);
//
//                                                        vk_ub=vk.get_ub(keyk);
//                                                        vk_lb=vk.get_lb(keyk);
//                                                        if(in_orig_model)
//                                                        {
//                                                            vk_xb=var_ub.eval(keyk);
//                                                            if(dirk=="UB"){
//                                                                vk_lb=vk_xb;
//                                                            }
//                                                            if(dirk=="LB"){
//                                                                vk_ub=vk_xb;
//                                                            }
//
//                                                        }
//
//                                                        ub_vp_new=lb_vp+(vk.get_ub(keyk)-vk.get_lb(keyk))/ldv[vpiter];
//                                                        var<> var_vp=this->get_var<T>(vp.second->_name);
//                                                        if((ub_vp-ub_vp_new)>=range_tol){
//                                                          // DebugOn("Upper Bound update"<<endl);
////                                                            DebugOn("Ith variable "<<vk._name<<endl);
////                                                            DebugOn("Ith variable index "<<keyk<<endl);
////                                                            DebugOn("Jth variable "<<var_vp._name<<endl);
////                                                            DebugOn("Jth variable index "<<(*vpkeys)[vpiter]<<endl);
////                                                             DebugOn(ub_vp<<"\t"<< ub_vp_new<<endl);
//                                                           var_vp.set_ub((*vpkeys)[vpiter], ub_vp_new);
//                                                        }
//                                                    }
//
//                                                    }
//
//                                                }
//                                            }

                                        }
                                        //If interval becomes smaller than range_tol, reset bounds so that interval=range_tol
                                        if(std::abs(vk.get_ub(keyk)-vk.get_lb(keyk))<range_tol)
                                        {
                                            //If original interval is itself smaller than range_tol, do not have to reset interval
                                            if(interval_original[var_key_k]>=range_tol)
                                            {
                                                DebugOn("Entered reset");
                                                //Mid is the midpoint of interval
                                                mid=(vk.get_ub(keyk)+vk.get_lb(keyk))/2.0;
                                                left=mid-range_tol/2.0;
                                                right=mid+range_tol/2.0;
                                                //If resized interval does not cross original bounds, reset
                                                if(right<=ub_original[var_key_k] && left>=lb_original[var_key_k])
                                                {
                                                    vk.set_ub(keyk, right);
                                                    vk.set_lb(keyk, left);
                                                }
                                                //If resized interval crosses original upperbound, set the new bound to upperbound, and lower bound is expanded to upperbound-range_tolerance
                                                else if(right>ub_original[var_key_k])
                                                {
                                                    
                                                    vk.set_ub(keyk, ub_original[var_key_k]);
                                                    vk.set_lb(keyk, ub_original[var_key_k]-range_tol);
                                                }
                                                //If resized interval crosses original lowerbound, set the new bound to lowerbound, and upper bound is expanded to lowerbound+range_tolerance
                                                else if(left<lb_original[var_key_k])
                                                {
                                                    vk.set_lb(keyk, lb_original[var_key_k]);
                                                    vk.set_ub(keyk, lb_original[var_key_k]+range_tol);
                                                    
                                                }
                                                //In the resized interval both original lower and upper bounds can not be crosses, because original interval is greater
                                                //than range_tol
                                                
                                            }
                                        }
                                        auto lift_name="Lift("+vkname+"^2)";
                                        if(this->_vars_name.find(lift_name)!=this->_vars_name.end()){
                                            auto liftv=this->get_var<T>(lift_name);
                                            auto lbvk=vk.get_lb(keyk);
                                            auto ubvk=vk.get_ub(keyk);
                                            auto temp_a=std::max(std::max(0.0, lbvk), ubvk*(-1));
                                            auto lift_lb=std::min(temp_a*temp_a, ubvk*ubvk);
                                            liftv.set_lb(keyk, lift_lb);
                                            liftv.set_ub(keyk, std::max(lbvk*lbvk, ubvk*ubvk));
                                            
                                        }
                                    }
                                    else
                                    {
                                              //    model->print();
                                        solver<> SDPLB_model(*model,solv_type);
                                        SDPLB_model.run(output = 0, tol, "ma27");
                                        DebugOn("OBBT step has failed in iteration\t"<<iter<<endl);
                                        
                                        
                                        //                                        fixed_point[pk]=true;
                                    }
                                }
                                batch_models.clear();
                            }
                        }
                    }
                }
            }
            
            //Check if OBBT has converged, can check every gap_count_int intervals
            if(iter%gap_count_int==0)
            {    solver_time= get_wall_time()-solver_time_start;

                                   //this->print();

                this->reset_constrs();
                this->reindex();
                solver<> SDPLB1(*this,solv_type);
               // SDPLB1.run(output = 5, tol, "ma57");
                SDPLB1.run(output = 0, tol);
                if(this->_status==0)
                {
                    lower_bound=this->get_obj_val()*upper_bound.second;
                    auto gap = 100*(upper_bound.second - lower_bound)/upper_bound.second;
                    DebugOn("Gap "<<gap<<" at iteration "<<iter<<" and solver time "<<solver_time<<endl);
                }
                
              
                if (upper_bound.second- lower_bound<=abs_tol && (upper_bound.second- lower_bound)/(upper_bound.second+zero_tol)<=rel_tol)
                {
                    
                    //                        this->print();
                    //                        this->print_solution();
                    //                         this->print_nonzero_constraints(tol);
                    
                    DebugOn("Gap closed at iter "<< iter<<endl);
                    DebugOn("Initial Gap Nonlinear = " << to_string(gapnl) << "%."<<endl);
                    gap = 100*(upper_bound.second - lower_bound)/upper_bound.second;
                    DebugOn("Final Gap = " << to_string(gap) << "%."<<endl);
                    DebugOn("Upper bound = " << to_string(upper_bound.second) << "."<<endl);
                    DebugOn("Lower bound = " << to_string(lower_bound) << "."<<endl);
                    DebugOn("Time\t"<<solver_time<<endl);
                    close=true;
                    terminate=true;
                    
                    
                }
//                else
//                {
//                    if(!nonlin) //If linearized model is used for obbt
//                    {
//                        
//                    }
//                }
            }
            
            if(break_flag==true)
            {
                DebugOn("Maximum Time Exceeded\t"<<max_time<<endl);
                DebugOn("Iterations\t"<<iter<<endl);
                
                break;
            }
            solver_time= get_wall_time()-solver_time_start;
            DebugOn("Solved Fixed Point iteration " << iter << endl);
        }
          
        vector<double> interval_gap;
        
        for(auto &it:this->_vars_name)
        {
            string vname=it.first;
            v=this->get_var<double>(vname);
            auto v_keys=v.get_keys();
            bool in_orig_model=false;
            if(ub_model->_vars_name.find(vname)!=ub_model->_vars_name.end())
            {
                var_ub=ub_model->template get_var<T>(vname);
                in_orig_model=true;
            }
            
            for(auto &key: *v_keys)
            { num_var++;
                var_key=vname+"|"+ key;
                interval_gap.push_back((interval_original[var_key]-interval_new[var_key])/(interval_original[var_key]+zero_tol)*100.0);
                sum+=interval_gap.back();
                if( in_orig_model)
                {
                    var_ub.uneval();
                    if((var_ub.eval(key)-v.get_lb(key)) <0.000 || (var_ub.eval(key)-v.get_ub(key))>0.000){
                        xb_true=false;
                         DebugOn("xb false Variable " <<vname<< " key "<< key<< " UB_value " <<var_ub.eval(key) <<"OBBT, lb, ub "<< v.get_lb(key)<<" "<< v.get_ub(key)<<endl);
                    }
                }
                DebugOff(var_key<<" " << interval_gap.back()<< " LB flag = " << fixed_point.at(var_key+"|LB") << endl);
                DebugOff(var_key<<" " << interval_gap.back()<< " UB flag = " << fixed_point.at(var_key+"|UB") << endl);
            }
            
        }
        avg=sum/num_var;
        
        DebugOff("Average interval reduction\t"<<avg<<endl);
        
        if(!close)
        {
            
            this->reset_constrs();
            
                      
            solver<> SDPLB1(*this,solv_type);
            
            SDPLB1.run(output = 0, tol, "ma27");
        }
        
    }
    //    avg=sum/num_var;
    //
    //    DebugOn("Average interval reduction\t"<<avg<<endl);
    
    if(!close)
    {
#ifdef USE_MPI
        if(worker_id==0){
#endif
            this->reset_constrs();
            solver<> SDPLB1(*this,solv_type);
            
            
            SDPLB1.run(output = 0, tol);
            //                this->print_constraints_stats(tol);
            //                bool print_only_relaxed;
            //                this->print_nonzero_constraints(tol,print_only_relaxed=true);
            
            //        SDP->print_solution();
            
            //        SDP->print();
            
            if(this->_status==0)
            {
                
                DebugOff("\nResults: " << " " << to_string(this->get_obj_val()*upper_bound.second) << " " <<endl);  
                DebugOff("Solution Print"<<endl);
                //                    this->print();
                //                SDP->print_solution();
                //                    this->print_constraints_stats(tol);
                
                
                DebugOn("Initial Gap Nonlinear = " << to_string(gapnl) << "%."<<endl);
                lower_bound=this->get_obj_val()*upper_bound.second;
                gap = 100*(upper_bound.second - lower_bound)/upper_bound.second;
                DebugOn("Final Gap = " << to_string(gap) << "%."<<endl);
                DebugOn("Upper bound = " << to_string(upper_bound.second) << "."<<endl);
                DebugOn("Lower bound = " << to_string(lower_bound) << "."<<endl);
                DebugOn("Time\t"<<solver_time<<endl);
                
            }
            else
            {
                DebugOn("Initial Gap = " << to_string(gapnl) << "%."<<endl);
                DebugOn("Lower bounding problem status = " << this->_status <<endl);
                DebugOn("Lower bounding problem not solved to optimality, cannot compute final gap"<<endl);
            }
            if(time_limit){
                DebugOn("Reached Time limit!"<<endl);
            }
            else {
                DebugOn("Terminate\t"<<terminate<<endl);
            }
            
            
            DebugOn("Time\t"<<solver_time<<endl);
            DebugOn("Iterations\t"<<iter<<endl);
#ifdef USE_MPI
        }
        #endif

    }
        }
        else
        {
         DebugOn("Lower bounding problem not solved to optimality, cannot compute initial gap"<<endl);
        }
         std::get<0>(res) = terminate;
         std::get<1>(res) = iter;
         std::get<2>(res) = solver_time;
         std::get<3>(res) = lower_bound_init;
         std::get<4>(res) = avg;
         std::get<5>(res) = xb_true;
        //this->print();
        return res;
}



//    string result_name=string(prj_dir)+"/results_obbt/"+grid._name+".txt";
//    ofstream fout(result_name.c_str(), ios_base::app);
//#ifdef USE_MPI
//    if(worker_id==0){
//#endif
//        fout<<grid._name<<"\t"<<std::fixed<<std::setprecision(5)<<gapnl<<"\t"<<std::setprecision(5)<<upper_bound<<"\t"<<std::setprecision(5)<<lower_bound<<"\t"<<std::setprecision(5)<<gap<<"\t"<<terminate<<"\t"<<iter<<"\t"<<std::setprecision(5)<<solver_time<<endl;
//#ifdef USE_MPI
//        }
//#endif
//

    
    
    template std::tuple<bool,int,double,double,double,bool> gravity::Model<double>::run_obbt<double, (void*)0>(double, unsigned int, const pair<bool,double>&, unsigned int, shared_ptr<Model<double>>, shared_ptr<Model<double>>, bool);

    template Constraint<Cpx> Model<Cpx>::lift(Constraint<Cpx>& c, string model_type);



//    template void Model<double>::run_obbt(double max_time, unsigned max_iter);
//    template func<double> constant<double>::get_real() const;
//    template class Model<double>;
//    template class Model<Cpx>;

}


