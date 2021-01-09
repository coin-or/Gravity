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


const bool var_compare(const pair<string,shared_ptr<param_>>& v1, const pair<string,shared_ptr<param_>>& v2) {
    return v1.second->get_nb_rows() > v2.second->get_nb_rows();
}

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
            name = "Lift("+o1.get_name(true,true)+"|"+o2.get_name(true,true)+")";
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
                    auto keep_refs = ids.get_diff_refs(added);
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
            DebugOff("Changing the sign of the lifted variable." << endl);
            lift_sign = !lift_sign;
        }
        else{
            DebugOff("Keeping the sign of the lifted variable same." << endl);
        }
        
        //arrange the variables so that if they have the same base name, use them ordered in name
        auto o1_ptr = static_pointer_cast<var<type>>(term._p->first);
        auto o2_ptr = static_pointer_cast<var<type>>(term._p->second);
        auto o1 = *o1_ptr;
        auto o2 = *o2_ptr;
        if((o1 != o2) && (o1._indices->_keys->at(0) > o2._indices->_keys->at(0)) ){
            //        if((o1 != o2) && (o1.get_name(false,false) > o2.get_name(false,false)) ){
            o2_ptr = static_pointer_cast<var<type>>(term._p->first);
            o1_ptr = static_pointer_cast<var<type>>(term._p->second);
            o1 = *o1_ptr;
            o2 = *o2_ptr;
            DebugOff("O1 name "<< o1._name << endl);
            DebugOff("O2 name "<< o2._name << endl);
        }
        
        string name;
        indices ids;
        if(o1==o2){
            name = "Lift("+o1.get_name(true,true)+"^2)";
            ids = *o1._indices;
            ids.set_name(o1._name);
        }
        else {
            name = "Lift("+o1.get_name(true,true)+"|"+o2.get_name(true,true)+")";
            ids = combine(*o1._indices,*o2._indices);
            ids.set_name(o1._name+"|"+o2._name);
        }
        auto o1_ids = *o1._indices;
        auto o2_ids = *o2._indices;
        auto flat_ids = ids;
        if(ids.is_matrix_indexed()){/* Flatten matrix indexed sets */
            flat_ids.flatten();
            o1_ids.flatten();
            o2_ids.flatten();
        }
        auto unique_ids = flat_ids.get_unique_keys(); /* In case of an indexed variable, keep the unique keys only */
        if(unique_ids.size()!=flat_ids.size()){/* If some keys are repeated, remove them from the refs of o1 and o2 */
            auto keep_refs = flat_ids.get_unique_refs();
            o1_ids.filter_refs(keep_refs);
            o2_ids.filter_refs(keep_refs);
        }
        
        // collect the number of partitions of each variable
        int num_partns1 = *o1._num_partns;
        int num_partns2 = *o2._num_partns;
        
        func<double> lb, ub;
        
        //calculate the tightest valid bounds
        if(o1==o2) //if variables are same, calculate the bounds more efficiently
        {
            auto it = _vars_name.find(name);
            if(it!=_vars_name.end()){
                auto vlift = static_pointer_cast<var<type>>(it->second);
                vlift->_lb->merge_vars(*vlift->_ub);
                auto new_ids = unique_ids.get_diff_refs(*vlift->_indices);
                unique_ids.filter_refs(new_ids);
                o1_ids.filter_refs(new_ids);
                o2_ids.filter_refs(new_ids);
                auto o1_lb = vlift->get_square_lb();
                auto o1_ub = vlift->get_square_ub();
                o1_lb->_indices->add_refs(o1_ids);
                o1_ub->_indices->add_refs(o1_ids);
                vlift->_original_vars[0]->_indices->add_refs(o1_ids);
                vlift->_original_vars[0]->_lb->index_in(*vlift->_original_vars[0]->_indices);
                vlift->_original_vars[0]->_ub->index_in(*vlift->_original_vars[0]->_indices);
                vlift->_original_vars[1]->_indices->add_refs(o2_ids);
                vlift->_original_vars[1]->_lb->index_in(*vlift->_original_vars[1]->_indices);
                vlift->_original_vars[1]->_ub->index_in(*vlift->_original_vars[1]->_indices);
                
            }
            func<double> prod_b1 = (o1.get_lb()*o1.get_lb()).in(unique_ids);
            func<double> prod_b2 = (o1.get_lb()*o1.get_ub()).in(unique_ids);
            func<double> prod_b3 = (o1.get_ub()*o1.get_ub()).in(unique_ids);
            
            lb = gravity::max(gravity::min(gravity::min(prod_b1,prod_b2), prod_b3).in(unique_ids), func<type>());
            ub = gravity::max(gravity::max(prod_b1,prod_b2).in(unique_ids),prod_b3);
            
            
            //            lb = max(min(min(o1.get_lb()*o1.get_lb(),o1.get_lb()*o1.get_ub()), o1.get_ub()*o1.get_ub()), 0);
            //            ub = max(max(o1.get_lb()*o1.get_lb(),o1.get_lb()*o1.get_ub()), o1.get_ub()*o1.get_ub());
            //            for (int i=0; i<unique_ids.size(); i++) {
            //                //calculate all the possibilities and assign the worst case
            //                size_t id1;
            //                if(o1_ids._ids == nullptr){
            //                    id1 = i;
            //                }
            //                else id1 = o1_ids._ids->at(0).at(i);
            //                auto key1 = o1_ids._keys->at(id1);
            //
            //                auto prod_b1 = o1.get_lb(key1)*o1.get_lb(key1);
            //                auto prod_b2 = o1.get_lb(key1)*o1.get_ub(key1);
            //                auto prod_b3 = o1.get_ub(key1)*o1.get_ub(key1);
            //
            //
            //
            //                lb.set_val(key1, std::max(std::min(std::min(prod_b1,prod_b2), prod_b3), (type)0 ));
            //                ub.set_val(key1, std::max(std::max(prod_b1,prod_b2),prod_b3));
            //            }
        }
        else /* if variables are different */
        {
            auto it = _vars_name.find(name);
            if(it!=_vars_name.end()){
                auto vlift = static_pointer_cast<var<type>>(it->second);
                vlift->_lb->merge_vars(*vlift->_ub);
                auto new_ids = unique_ids.get_diff_refs(*vlift->_indices);
                unique_ids.filter_refs(new_ids);
                o1_ids.filter_refs(new_ids);
                o2_ids.filter_refs(new_ids);
                auto o1_lb = vlift->get_bilinear_lb1();
                auto o1_ub = vlift->get_bilinear_ub1();
                o1_lb->_indices->add_refs(o1_ids);
                o1_ub->_indices->add_refs(o1_ids);
                auto o2_lb = vlift->get_bilinear_lb2();
                auto o2_ub = vlift->get_bilinear_ub2();
                o2_lb->_indices->add_refs(o2_ids);
                o2_ub->_indices->add_refs(o2_ids);
                vlift->_original_vars[0]->_indices->add_refs(o1_ids);
                vlift->_original_vars[1]->_indices->add_refs(o2_ids);
                vlift->_original_vars[0]->_lb->index_in(*vlift->_original_vars[0]->_indices);
                vlift->_original_vars[0]->_ub->index_in(*vlift->_original_vars[0]->_indices);
                vlift->_original_vars[1]->_lb->index_in(*vlift->_original_vars[1]->_indices);
                vlift->_original_vars[1]->_ub->index_in(*vlift->_original_vars[1]->_indices);
            }
            
            auto lb1 = o1.get_lb().in(o1_ids);
            auto ub1 = o1.get_ub().in(o1_ids);
            auto lb2 = o2.get_lb().in(o2_ids);
            auto ub2 = o2.get_ub().in(o2_ids);
            func<double> prod_b1 = (lb1*lb2).in(unique_ids);
            func<double> prod_b2 = (lb1*ub2).in(unique_ids);
            func<double> prod_b3 = (ub1*lb2).in(unique_ids);
            func<double> prod_b4 = (ub1*ub2).in(unique_ids);
            
            lb = gravity::min(gravity::min(prod_b1,prod_b2).in(unique_ids),gravity::min(prod_b3,prod_b4).in(unique_ids));
            ub = gravity::max(gravity::max(prod_b1,prod_b2).in(unique_ids),gravity::max(prod_b3,prod_b4).in(unique_ids));
            //            for (int i=0; i<unique_ids.size(); i++) {
            //                //calculate all the possibilities and assign the worst case
            //                size_t id1;
            //                size_t id2;
            //                if(o1_ids._ids == nullptr){
            //                    id1 = i;
            //                }
            //                else id1 = o1_ids._ids->at(0).at(i);
            //                if(o2_ids._ids == nullptr){
            //                    id2 = i;
            //                }
            //                else id2 = o2_ids._ids->at(0).at(i);
            //                auto key1 = o1_ids._keys->at(id1);
            //                auto key2 = o2_ids._keys->at(id2);
            //
            //                auto prod_b1 = o1.get_lb(key1)*o2.get_lb(key2);
            //                auto prod_b2 = o1.get_lb(key1)*o2.get_ub(key2);
            //                auto prod_b3 = o1.get_ub(key1)*o2.get_lb(key2);
            //                auto prod_b4 = o1.get_ub(key1)*o2.get_ub(key2);
            //
            //                lb.set_val(key1+","+key2, std::min(std::min(prod_b1,prod_b2),std::min(prod_b3,prod_b4)));
            //                ub.set_val(key1+","+key2, std::max(std::max(prod_b1,prod_b2),std::max(prod_b3,prod_b4)));
            //            }
        }
        //        lb.update_var_indices(un)
        lb.index_in(unique_ids);
        ub.index_in(unique_ids);
        lb.eval_all();
        ub.eval_all();
        //        auto common_refs = o1._indices->get_common_refs(unique_ids);
        //        lb.update_rows(common_refs);
        //        ub.update_rows(common_refs);
        auto it = _vars_name.find(name);
        
        auto name1 = o1.get_name(true,true);
        auto name2 = o2.get_name(true,true);
        
        if(it==_vars_name.end()){
            
            
            
            //create the lifted variable with proper lower and upper bounds
            var<type> vlift(name, lb, ub);
            vlift._lift = true;
            vlift._original_vars.push_back(make_shared<var<type>>(*o1_ptr));
            vlift._original_vars.push_back(make_shared<var<type>>(*o2_ptr));
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
            param<T> lb_p("lb_p");
            lb_p = lb;
            param<T> ub_p("ub_p");
            ub_p = ub;
            auto added = vlift->add_bounds(lb_p,ub_p);
            lt._p = make_shared<var<type>>(vlift->in(ids));
            if(!added.empty()){
                vlift->_lb->update_vars();
                vlift->_ub->update_vars();
                assert(o1._indices->size()==o2._indices->size());
                if(added.size()!=o1_ids.size()){/* If some keys are repeated, remove them from the refs of o1 and o2 */
                    auto keep_refs = flat_ids.get_diff_refs(added);
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
template<typename T>
void Model<type>::populate_original_interval(shared_ptr<Model<type>>& obbt_model, map<string, bool>& fixed_point, map<string, double>& ub_original,map<string, double>& lb_original,map<string, double>& interval_original,map<string, double>& interval_new, int& count_skip, int& count_var, double range_tol){
    var<> v, var_ub;
    std::string var_key, key_lb, key_ub;
    double vi_ub_val;
    bool in_orig_model=false;
    for(auto &it:obbt_model->_vars)
    {
        string vname=(*it.second)._name;
        v=obbt_model->template get_var<double>(vname);
        auto v_keys=v.get_keys();
        auto v_key_map=v.get_keys_map();
        in_orig_model=false;
        if(this->_vars_name.find(vname)!=this->_vars_name.end()){
                var_ub=this->template get_var<T>(vname);
                in_orig_model=true;
        }
        for(auto &key: *v_keys)
        {
            var_key = vname+"|"+ key;
            key_lb= var_key +"|LB";
            key_ub= var_key +"|UB";
            /* Do not do OBBT on lifted variables */
            if(v._lift){
                fixed_point[key_lb]=true;
                fixed_point[key_ub]=true;
                count_skip+=2;
                DebugOff("Skipping OBBT for "<<vname<<"\t"<<key<<endl);
            }
            else{
                fixed_point[key_lb]=false;
                fixed_point[key_ub]=false;
                if(in_orig_model){
                    vi_ub_val=var_ub.eval(key);
                    if(vi_ub_val-v.get_lb(key)<=range_tol){
                        fixed_point[var_key+"|LB"]=true;
                        DebugOff("vi_ub "<<vi_ub_val<<" vlb "<<v.get_lb(key)<<" "<<key_lb<<endl);
                    }
                    if(v.get_ub(key)-vi_ub_val<=range_tol){
                        fixed_point[var_key+"|UB"]=true;
                        DebugOff("vi_ub "<<vi_ub_val<<" vlb "<<v.get_ub(key)<<" "<<key_ub<<endl);
                    }
                }
            }
            auto key_pos=v_key_map->at(key);
            
            if(v._off[key_pos]==true)
            {
                fixed_point[key_lb]=true;
                fixed_point[key_ub]=true;
                count_skip+=2;
                DebugOff("Off var: "<<vname<<"\t"<<key<<endl);
            }
            count_var++;
            interval_original[var_key]=v.get_ub(key)-v.get_lb(key);
            ub_original[var_key]=v.get_ub(key);
            lb_original[var_key]=v.get_lb(key);
            interval_new[var_key]=v.get_ub(key)-v.get_lb(key);
            
        }
        
    }
}
template <typename type>
template<typename T>
double Model<type>::populate_final_interval_gap(const shared_ptr<Model<type>>& obbt_model, const map<string, double>& interval_original, map<string, double>& interval_new, double& sum, bool& xb_true, const double zero_tol, int count_var){
    var<> v, var_ub;
    std::string var_key;
    double avg=0;
    for(auto &it:obbt_model->_vars_name)
    {
        string vname=it.first;
        v=obbt_model->template get_var<double>(vname);
        auto v_keys=v.get_keys();
        bool in_orig_model=false;
        if(this->_vars_name.find(vname)!=this->_vars_name.end())
        {
            var_ub=this->template get_var<T>(vname);
            in_orig_model=true;
        }
        for(auto &key: *v_keys)
        {
            var_key=vname+"|"+ key;
            interval_new[var_key]=v.get_ub(key)-v.get_lb(key);
            sum+=((interval_original.at(var_key)-interval_new.at(var_key))/(interval_original.at(var_key)+zero_tol)*100.0);
            if( in_orig_model)
            {
                var_ub.uneval();
                if((var_ub.eval(key)-v.get_lb(key)) < - 1e-6 || (var_ub.eval(key)-v.get_ub(key))>1e-6){
                    xb_true=false;
                    DebugOn("xb false Variable " <<vname<< " key "<< key<< " UB_value " <<var_ub.eval(key) <<"OBBT, lb, ub "<< v.get_lb(key)<<" "<< v.get_ub(key)<<endl);
                }
            }
        }
        
    }
    avg=sum/count_var;
    DebugOff("Average interval reduction "<<avg<<endl);
    return avg;
}
template <typename type>
template<typename T>
void Model<type>::create_batch_models(vector<shared_ptr<Model<type>>>& batch_models, int nb_threads, double ub_scale_value){
    for(auto i=0;i<nb_threads;i++){
        auto modelk = this->copy();
        param<> ub("ub");
        ub = ub_scale_value;
        auto obj = *modelk->_obj;
        if(modelk->_cons_name.count("obj|ub")==0){
            Constraint<type> obj_ub("obj|ub");
            obj_ub = obj - ub;
            modelk->add(obj_ub<=0);
        }
        batch_models.push_back(modelk);
        batch_models.at(i)->set_name(to_string(i));
    }
}
template <typename type>
template<typename T>
void Model<type>::batch_models_obj_lb_constr(vector<shared_ptr<Model<type>>>& batch_models, int nb_threads, double lower_bound_lin, double lower_bound_old, double lower_bound_nonlin_init, double upper_bound, double ub_scale_value){
    double lb, lb_old;
    lb=std::max(lower_bound_lin, lower_bound_nonlin_init)/upper_bound*ub_scale_value;
    lb_old=std::max(lower_bound_old, lower_bound_nonlin_init)/upper_bound*ub_scale_value;
    for(auto& modelk:batch_models){
        if(modelk->_cons_name.count("obj|lb")==0){
            auto obj = *modelk->_obj;
            param<> lb_p("lb_p");
            lb_p=lb;
            Constraint<type> obj_lb("obj|lb");
            obj_lb = obj - lb;
            modelk->add(obj_lb>=0);
        }
        else{
            auto con=modelk->get_constraint("obj|lb");
            if(lb>lb_old){
                modelk->remove("obj|lb");
                Constraint<> a(*con);
                modelk->add(a>=(lb-lb_old));
            }
        }
    }
}

//Check if OBBT has converged, can check every gap_count_int intervals
template <typename type>
template<typename T>
void Model<type>::compute_iter_gap(double& gap, double& active_tol, bool& terminate, bool linearize, int iter, shared_ptr<Model<type>>& obbt_model, const Model<type>& interior_model, SolverType lb_solver_type, int nb_root_refine, const double upper_bound, double& lower_bound, const double ub_scale_value, double lb_solver_tol, double& active_root_tol, int& oacuts, const double abs_tol, const double rel_tol, const double zero_tol, string lin_solver, int max_iter, int max_time, vector<double>& vrbasis, std::map<string,double>& crbasis, bool initialize_primal){
    gap=-999;
#ifdef USE_MPI
    int worker_id, nb_workers;
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
    auto err_size = MPI_Comm_size(MPI_COMM_WORLD, &nb_workers);
#endif
    int output=0;
    bool close=false;
    obbt_model->reset_lifted_vars_bounds();
    obbt_model->reset_constrs();
    if(!linearize){
        solver<> LB_solver(obbt_model,lb_solver_type);
        LB_solver.run(output = 0, lb_solver_tol, max_iter, max_time);
        if(obbt_model->_status==0)
        {
            lower_bound=obbt_model->get_obj_val()*upper_bound/ub_scale_value;
	     if (std::abs(upper_bound- lower_bound)<=abs_tol && ((upper_bound- lower_bound))/(std::abs(upper_bound)+zero_tol)<=rel_tol)
            {
                close=true;
                //obbt_model->print();
                //obbt_model->print_solution();
            }
        }
    }
    else{
        if(active_tol>lb_solver_tol){
                active_tol*=0.1;
        }
        if(active_root_tol>lb_solver_tol){
                active_root_tol*=0.1;
        }
        close=this->root_refine(interior_model, obbt_model, lb_solver_type, nb_root_refine, upper_bound, lower_bound, ub_scale_value, lb_solver_tol, active_root_tol, oacuts,  abs_tol, rel_tol, zero_tol, "ma27", max_iter, max_time, vrbasis, crbasis, initialize_primal);
    }
    DebugOff("lower bound "<<lower_bound<<endl);
    if(obbt_model->_status==0)
    {
        gap = 100*(upper_bound - lower_bound)/std::abs(upper_bound);
        if(close)
	{
            terminate=true;
        }
    }
    else{
        DebugOn("Failed to solve lower bounding problem"<<endl);
        lower_bound=numeric_limits<double>::min();
        terminate=true;
    }
#ifdef USE_MPI
    if(worker_id==0)
#endif
        DebugOn("Gap "<<gap<<" at iteration "<<iter<<endl);
}
template <typename type>
template<typename T,
typename std::enable_if<is_same<T,double>::value>::type*>
std::tuple<bool,int,double,double,double,double,double,double,int,int,int> Model<type>::run_obbt(shared_ptr<Model<T>> relaxed_model, double max_time, unsigned max_iter, double rel_tol, double abs_tol, unsigned nb_threads, SolverType ub_solver_type, SolverType lb_solver_type, double ub_solver_tol, double lb_solver_tol, double range_tol, bool linearize, bool scale_objective, int nb_refine,  int nb_root_refine,  double viol_obbt_init, double viol_root_init, bool initialize_primal) {
#ifdef USE_MPI
    int worker_id, nb_workers;
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
    auto err_size = MPI_Comm_size(MPI_COMM_WORLD, &nb_workers);
#endif
    std::tuple<bool,int,double,double,double,double,double,double,int,int,int> res;
    int total_iter=0, global_iter=1,fail=0;
    int output, oacuts_init=0, oacuts=0;
    double total_time =0, time_start = get_wall_time(), time_end = 0, lower_bound_nonlin_init = numeric_limits<double>::min();
#ifdef USE_MPI
    if(worker_id==0){
        time_start = get_wall_time();
    }
#else
    time_start = get_wall_time();
#endif
    double gap_new=-999, gap=0, ub_scale_value;
    const double gap_tol=rel_tol;
    solver<> UB_solver(*this,ub_solver_type);
    UB_solver.run(output = 0, ub_solver_tol, 2000, 600);
    ub_scale_value=this->get_obj_val();
    solver<> LBnonlin_solver(relaxed_model,ub_solver_type);
    if(scale_objective){
        auto obj = *relaxed_model->_obj/ub_scale_value;
        relaxed_model->min(obj);
        relaxed_model->reset();
        ub_scale_value=1.0;
    }
    LBnonlin_solver.run(output = 0 , lb_solver_tol, 2000, 600);
    if(relaxed_model->_status==0)
    {
        lower_bound_nonlin_init = relaxed_model->get_obj_val()*this->get_obj_val()/ub_scale_value;;
        DebugOff("Initial lower bound = "<<lower_bound_nonlin_init<<endl);
    }
    shared_ptr<Model<>> obbt_model=relaxed_model->copy();
    obbt_model->_status=relaxed_model->_status;
    Model<> interior_model;
    if(linearize){
        auto lin_model=relaxed_model->buildOA();
        interior_model=lin_model->add_outer_app_solution(*relaxed_model);
        obbt_model=lin_model;
        oacuts_init=lin_model->_nb_cons;
        oacuts=lin_model->_nb_cons;
    }
    int run_obbt_iter=1;
    auto status = run_obbt_one_iteration(relaxed_model, max_time, max_iter, rel_tol, abs_tol, nb_threads, ub_solver_type, lb_solver_type, ub_solver_tol, lb_solver_tol, range_tol, linearize, obbt_model, interior_model, oacuts, oacuts_init, run_obbt_iter, ub_scale_value, time_start, nb_refine, nb_root_refine, viol_obbt_init, viol_root_init, initialize_primal);
    double upper_bound = get<5>(status);
    
    total_iter += get<1>(status);
    auto lower_bound=get<6>(status);
    gap = (upper_bound - lower_bound)/std::abs(upper_bound)*100;
    fail+=get<10>(status);
    while((gap > rel_tol*100.0 || (upper_bound-lower_bound)>abs_tol) && obbt_model->_status==0){
        //MPI_Barrier(MPI_COMM_WORLD);
        auto time=get_wall_time()-time_start;
        if(time>=max_time){
            DebugOff("Maximum time exceeded"<<endl);
            break;
        }
        if(total_iter>= max_iter){
            DebugOff("Maximum iterations exceeded"<<endl);
            break;
        }
        if(!linearize && get<1>(status)==1 && max_iter>1)
            break;
        if(linearize && get<1>(status)==1 && run_obbt_iter>=4)
            break;
        oacuts=get<8>(status);
        run_obbt_iter++;
        status = run_obbt_one_iteration(relaxed_model, max_time, max_iter, rel_tol, abs_tol, nb_threads, ub_solver_type, lb_solver_type, ub_solver_tol, lb_solver_tol, range_tol, linearize, obbt_model, interior_model, oacuts, oacuts_init, run_obbt_iter, ub_scale_value, time_start, nb_refine, nb_root_refine, viol_obbt_init, viol_root_init, initialize_primal);
        lower_bound=get<6>(status);
        gap_new = (upper_bound - lower_bound)/std::abs(upper_bound)*100;
        total_iter += get<1>(status);
        if(get<1>(status)>0)
            global_iter++;
        gap=gap_new;
        fail+=get<10>(status);
    }
#ifdef USE_MPI
    if(worker_id==0){
        time_end = get_wall_time();
    }
#else
    time_end = get_wall_time();
#endif
    total_time = time_end - time_start;
    //obbt_model->print_constraints_stats(1e-8);
    get<0>(res)=get<0>(status);
    get<1>(res)=total_iter;
    get<2>(res)=total_time;
    get<3>(res)=lower_bound_nonlin_init;
    get<4>(res)=get<4>(status);
    get<5>(res)=get<5>(status);
    get<6>(res)=get<6>(status);
    get<7>(res)=get<7>(status);
    get<8>(res)=get<8>(status);
    get<9>(res)=oacuts_init;
    get<10>(res)=fail;
    upper_bound=get<5>(status);
    DebugOff("Total wall-clock time spent in OBBT = " << total_time << endl);
    DebugOff("Total number of OBBT iterations = " << total_iter << endl);
    DebugOff("Number of global iterations = " << global_iter << endl);
    auto gapnl=(upper_bound-lower_bound_nonlin_init)/std::abs(upper_bound)*100;
    DebugOff("Initial gap = "<<gapnl<<"%"<<endl);
    if(obbt_model->_status==0){
        auto lower_bound_final=get<6>(status);
        auto gap_final = 100*(upper_bound - lower_bound_final)/std::abs(upper_bound);
        DebugOff("Final gap = " << to_string(gap_final) << "%."<<endl);
    }
    return res;
}
/** function to run one global iteration of the obbt algorithm
 @param[in] relaxed model: original nonlinear nonconvex model
 @param[in] max_time: max_time to run the run_obbt algorithm
 @param[in] max_iter: max_iter to run the run_obbt algorithm
 @param[in] rel_tol,abs_tol: To determine when the run_obbt is converged
 @param[in] nb_threads: threads per machine
 @param[in] ub_solver_type, lb_solver_type: solver types for upper bound problem and all lower and obbt problems
 @param[in] ub_solver_tol, lb_solver_tol: tolerances for upper bound and lower bound problems
 @param[in] range_tol: a variable interval is closed if <= range_tol
 @param[in] linearize: true if linear obbt algorithm is used
 @param[in] obbt_model: model which is used to compute lower bound. Linear in the case of linearize, else is initially identical to relaxed_model
 @param[in] interior_model: model to compute interior point of relaxed model
 @param[in] oacuts, oacuts_init: initial oa cuts at start of run_obbt
 @param[in] run_obbt_iter: global iteration of run_obbt
 @param[in] ub_scale_value: if scaling is true in run_obbt, ub_scale_value=1, else equals upper bound
 @param[in] solver_time_start: time at which run_obbt started
 @param[in] nb_refine: number of refinement steps, used for linearizs only
 @return returns tuple of terminate status, no of iter of run_obbt_one_iter, solver_time,lower_bound_nonlin_init,lower_bound_init(linear), upper_bound,lower_bound, avergae interval reduction,oacuts (on obbt_model) at end of run_obbt_one_iter,oacuts_init(no. cuts on obbt model after root refine in run_obbt_one_iter, total number of failed obbt subproblems in run_obbt_one_iter
 */
template <typename type>
template<typename T,
typename std::enable_if<is_same<T,double>::value>::type*>
std::tuple<bool,int,double,double,double,double,double,double,int,int,int> Model<type>::run_obbt_one_iteration(shared_ptr<Model<T>> relaxed_model, double max_time, unsigned max_iter, double rel_tol, double abs_tol, unsigned nb_threads, SolverType ub_solver_type, SolverType lb_solver_type, double ub_solver_tol, double lb_solver_tol, double range_tol, bool linearize, shared_ptr<Model<T>> obbt_model, Model<T> & interior_model, int oacuts, int oacuts_init, int run_obbt_iter, double ub_scale_value, double solver_time_start, int nb_refine,  int nb_root_refine,  double viol_obbt_init, double viol_root_init, bool initialize_primal){
    std::tuple<bool,int,double, double, double, double, double, double, int, int,int> res;
#ifdef USE_MPI
    int worker_id, nb_workers, nb_workers_;
    auto err_rank = MPI_Comm_rank(MPI_COMM_WORLD, &worker_id);
    auto err_size = MPI_Comm_size(MPI_COMM_WORLD, &nb_workers);
#endif
    int nb_total_threads = nb_threads; /** Used when MPI is ON to multiply with the number of workers */
#ifdef USE_MPI
    nb_total_threads *= nb_workers;
#endif
    vector<shared_ptr<Model<>>> batch_models;
    vector<string> objective_models;
    vector<double> sol_obj, vrbasis;
    vrbasis.push_back(0);
    std::map<string,double> crbasis;
    vector<int> sol_status;
    vector<vector<double>> vbasis;
    vector<std::map<string,double>> cbasis;
    vbasis.resize(nb_threads);
    cbasis.resize(nb_threads);
    map<string, bool> fixed_point;
    map<string, double> interval_original, interval_new, ub_original, lb_original;
    map<int, double> map_lb,map_ub;
    string vname, var_key, mname, cut_type="allvar";
    string dir_array[2]={"LB", "UB"};
    var<> v;
    bool close=false, terminate=false, xb_true=true, alg_batch_reset=true;
    const double fixed_tol_abs=1e-3, fixed_tol_rel=1e-3, zero_tol=1e-6, obbt_subproblem_tol=1e-6;
    int iter=0, fail=0, count_var=0, count_skip=0, nb_init_refine=nb_refine;
    double solver_time =0, gapnl,gap, gaplin=-999, sum=0, avg=0, active_root_tol=lb_solver_tol, active_tol=1e-6;
    double lower_bound_nonlin_init = numeric_limits<double>::min(), lower_bound_init = numeric_limits<double>::min(), upper_bound = 0, lower_bound = numeric_limits<double>::min(), lower_bound_old;
    map<string,int> old_map;
    if(this->_status==0){
        upper_bound=this->get_obj_val();
        //if(relaxed_model->_status==0){
            lower_bound_nonlin_init=relaxed_model->get_obj_val()*upper_bound/ub_scale_value;
            lower_bound_init = lower_bound_nonlin_init;
            lower_bound = lower_bound_nonlin_init;
            gapnl=(upper_bound-lower_bound_nonlin_init)/std::abs(upper_bound)*100;
            /* Check if gap is already not zero at root node */
            if ((upper_bound-lower_bound_nonlin_init)>=abs_tol || (upper_bound-lower_bound_nonlin_init)/(std::abs(upper_bound)+zero_tol)>=rel_tol){
                if(linearize){
                    /*Set values of active tol and nb_init_refine*/
                    set_activetol_initrefine(active_tol, active_root_tol, viol_obbt_init, viol_root_init, nb_init_refine, nb_root_refine, lb_solver_tol, run_obbt_iter);
                    /*Root refine obbt_model*/
                    close=relaxed_model->root_refine(interior_model, obbt_model, lb_solver_type, nb_init_refine, upper_bound, lower_bound_init, ub_scale_value, lb_solver_tol, active_root_tol, oacuts,  abs_tol, rel_tol, zero_tol, "ma27", 10000, 2000, vrbasis, crbasis, initialize_primal);
                    oacuts_init=oacuts;
                    gaplin=(upper_bound-lower_bound_init)/std::abs(upper_bound)*100;
                    lower_bound_old=lower_bound_init;
                }
                //if(obbt_model->_status==0){
                    /*Initialize fixed point, interval original and new, bounds original*/
                    this->populate_original_interval(obbt_model, fixed_point, ub_original,lb_original, interval_original,interval_new,  count_skip, count_var,range_tol);
                    solver_time= get_wall_time()-solver_time_start;
                    /*Create nb_threads copy of obbt_models*/
                    obbt_model->create_batch_models(batch_models, nb_threads, ub_scale_value);
                    if(linearize){
                        if(initialize_primal && lb_solver_type==gurobi){
                            initialize_basis_vectors(lb_solver_type, vbasis,cbasis,vrbasis,crbasis,nb_threads);
                        }
                    }
                    /*Run obbt algorithm until termiante is true, iter and time less than max iter and max time*/
                    while(solver_time<=max_time && !terminate && iter<max_iter){
                        iter++;
                        terminate=true;
                        for (auto it=obbt_model->_vars.begin(); it!=obbt_model->_vars.end(); it++){
                            vname=it->second->_name;
                            v = obbt_model->template get_var<double>(vname);
                            auto v_keys=v.get_keys();
                            for(auto it_key=v.get_keys()->begin(); it_key!=v.get_keys()->end(); it_key++){
                                auto key = *it_key;
                                var_key=vname+"|"+ key;
                                /* Add to batch if not reached fixed point, or if we're at the last key of the last variable */
                                if(fixed_point[var_key +"|LB"]==false || fixed_point[var_key +"|UB"]==false || (next(it)==obbt_model->_vars.end() && next(it_key)==v.get_keys()->end())){
                                    /* Loop on Min/Max, upper bound and lower bound */
                                    for(auto &dir: dir_array){
                                        mname=var_key+"|"+dir;
                                        if(fixed_point[mname]==false){
                                            objective_models.push_back(mname);
                                        }
                                        /* When batch models has reached size of nb_threads or when at the last key of last variable */
                                        if (objective_models.size()==nb_total_threads || (next(it)==obbt_model->_vars.end() && next(it_key)==v.get_keys()->end() && dir=="UB")){
#ifdef USE_MPI
                                            nb_workers_= run_MPI_new(objective_models, sol_obj, sol_status, batch_models, relaxed_model, interior_model, cut_type, active_tol, lb_solver_type, obbt_subproblem_tol, nb_threads, "ma27", 10000, 600, linearize, nb_refine, old_map, vbasis, cbasis, initialize_primal);
#else
                                            auto viol= run_parallel_new(objective_models, sol_obj, sol_status, batch_models, relaxed_model, interior_model, cut_type, active_tol, lb_solver_type, obbt_subproblem_tol, nb_threads, "ma27", 10000, 600, linearize, nb_refine, vbasis, cbasis, initialize_primal);
#endif
                                            auto b=this->obbt_batch_update_bounds( objective_models,  sol_obj, sol_status,  batch_models,obbt_model,  fixed_point,  interval_original,  ub_original,  lb_original, terminate,  fail, range_tol, fixed_tol_abs, fixed_tol_rel,  zero_tol, iter);
                                            this->generate_lagrange_bounds(objective_models, batch_models, obbt_model, fixed_point, range_tol, zero_tol, map_lb, map_ub);
#ifdef USE_MPI
                                            send_lagrange_bounds(nb_workers_, map_lb, map_ub);
#endif
                                            this->obbt_update_lagrange_bounds(batch_models, obbt_model,   fixed_point, interval_original, ub_original, lb_original, terminate, fail, range_tol, fixed_tol_abs, fixed_tol_rel, zero_tol, iter, map_lb, map_ub);
                                            sol_status.clear();
                                            sol_obj.clear();
                                            objective_models.clear();
                                            map_lb.clear();
                                            map_ub.clear();
                                            solver_time=get_wall_time()-solver_time_start;     
#ifdef USE_MPI
                                            if(worker_id==0)
#endif
                                                DebugOff("batch: "<<solver_time<<endl); 
                                        }
                                        solver_time=get_wall_time()-solver_time_start;
                                        if(solver_time>=max_time){
                                            break;
                                        }
                                    }
                                }
                                if(solver_time>=max_time){
                                    break;
                                }
                            }
                            if(solver_time>=max_time){
                                break;
                            }
                        }
                        /*Compute gap at the end of iter, adjusts active tol and root refine if linearize*/
                        relaxed_model->compute_iter_gap(gap, active_tol, terminate, linearize,iter, obbt_model, interior_model, lb_solver_type, nb_root_refine, upper_bound, lower_bound, ub_scale_value, lb_solver_tol, active_root_tol, oacuts, abs_tol, rel_tol, zero_tol, "ma27", 10000, 2000, vrbasis, crbasis, initialize_primal);
                        if(linearize && !terminate){
                                batch_models.clear();
                                obbt_model->create_batch_models(batch_models, nb_threads, ub_scale_value);
                                if(initialize_primal && lb_solver_type==gurobi){
                                     initialize_basis_vectors(lb_solver_type, vbasis,cbasis,vrbasis,crbasis,nb_threads);
                            }
                        }
                        solver_time= get_wall_time()-solver_time_start;
#ifdef USE_MPI
                        if(worker_id==0)
#endif
                            DebugOn("Gap time "<<solver_time<<endl);
                    }
                    /*average interval reduction,final interval, sanity check on bounds*/
                    avg=this->populate_final_interval_gap(obbt_model,  interval_original, interval_new, sum,  xb_true, zero_tol, count_var);
//                }
//                else{
//                    DebugOn("Initial lower bounding problem not solved to optimality, cannot compute initial gap"<<endl);
//                    lower_bound=numeric_limits<double>::min();
//                }
            }
//        }
//        else{
//            DebugOn("Lower bounding problem not solved to optimality, cannot compute initial gap"<<endl);
//            lower_bound=numeric_limits<double>::min();
//        }
    }
    else{
        DebugOn("Upper bounding problem not solved to optimality, cannot compute gap"<<endl);
        upper_bound=numeric_limits<double>::min();
    }
    std::get<0>(res) = terminate;
    std::get<1>(res) = iter;
    std::get<2>(res) = solver_time;
    std::get<3>(res) = lower_bound_nonlin_init;
    std::get<4>(res) = lower_bound_init;
    std::get<5>(res) = upper_bound;
    std::get<6>(res) = lower_bound;
    std::get<7>(res) = avg;
    std::get<8>(res) = oacuts;
    std::get<9>(res) = oacuts_init;
    std::get<10>(res)=fail;
    return res;
}
template std::tuple<bool,int,double,double,double,double,double,double,int,int,int> gravity::Model<double>::run_obbt<double, (void*)0>(shared_ptr<Model<double>> relaxed_model, double max_time, unsigned max_iter, double rel_tol, double abs_tol, unsigned nb_threads, SolverType ub_solver_type, SolverType lb_solver_type, double ub_solver_tol, double lb_solver_tol, double range_tol, bool linearize, bool scale_objective, int nb_refine,  int nb_root_refine,  double viol_obbt_init, double viol_root_init, bool initialize_primal);

template std::tuple<bool,int,double,double,double,double,double,double,int,int,int> gravity::Model<double>::run_obbt_one_iteration<double, (void*)0>(shared_ptr<Model<double>> relaxed_model, double max_time, unsigned max_iter, double rel_tol, double abs_tol, unsigned nb_threads, SolverType ub_solver_type, SolverType lb_solver_type, double ub_solver_tol, double lb_solver_tol, double range_tol, bool linearize, shared_ptr<Model<double>> obbt_model, Model<double> & interior_model, int oacuts, int oacuts_init, int run_obbt_iter, double ub_value, double solver_time_start, int nb_refine, int nb_root_refine, double viol_obbt_init, double viol_root_init, bool initialize_primal);

template Constraint<Cpx> Model<Cpx>::lift(Constraint<Cpx>& c, string model_type);
template Constraint<> Model<>::lift(Constraint<>& c, string model_type);
template void Model<double>::populate_original_interval(shared_ptr<Model<double>>& obbt_model, map<string, bool>& fixed_point, map<string, double>& ub_original,map<string, double>& lb_original,map<string, double>& interval_original,map<string, double>& interval_new, int& count_skip, int& count_var, double range_tol);
//template void Model<double>::populate_original_interval(map<string, bool>& fixed_point, map<string, double>& ub_original,map<string, double>& lb_original,map<string, double>& interval_original,map<string, double>& interval_new, int& count_skip, int& count_var);
template double Model<double>::populate_final_interval_gap(const shared_ptr<Model<double>>& obbt_model, const map<string, double>& interval_original, map<string, double>& interval_new, double& sum, bool& xb_true, const double zero_tol, int count_var);
template void Model<double>::create_batch_models(vector<shared_ptr<Model<double>>>& batch_models, int nb_threads, double ub_scale_value);
template void Model<double>::compute_iter_gap(double& gap, double& active_tol, bool& terminate, bool linearize, int iter, shared_ptr<Model<double>>& obbt_model, const Model<double>& interior_model, SolverType lb_solver_type, int nb_root_refine, const double upper_bound, double& lower_bound, const double ub_scale_value, double lb_solver_tol, double& active_root_tol, int& oacuts, const double abs_tol, const double rel_tol, const double zero_tol, string lin_solver, int max_iter, int max_time, vector<double>& vrbasis, std::map<string,double>& crbasis, bool initialize_primal);
template void Model<double>::batch_models_obj_lb_constr(vector<shared_ptr<Model<double>>>& batch_models, int nb_threads, double lower_bound_lin, double lower_bound_old, double lower_bound_nonlin_init, double upper_bound, double ub_scale_value);


//    template func<double> constant<double>::get_real() const;
//    template class Model<double>;
//    template class Model<Cpx>;

}





