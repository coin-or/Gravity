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
    Constraint<type> Model<type>::lift(Constraint<type>& c, string model_type, bool add_McCormick){
        if(c.is_constant() || c.is_linear()){
            return c;
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
                name = "Lift("+o1.get_name(true,true)+";"+o2.get_name(true,true)+")";
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
//                vlift._lift = true;
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
//                vlift._lift=true;
                add(vlift.in(ids));
                lt._p = make_shared<var<Cpx>>(vlift);
            }
            else {
                vlift = *static_pointer_cast<var<Cpx>>(it->second);
                lt._p = make_shared<var<Cpx>>(vlift);
            }
        }
        lifted.insert(lt);
    }


template <typename type>
template<typename T,typename std::enable_if<is_arithmetic<T>::value>::type*>
Constraint<type> Model<type>::lift(Constraint<type>& c, string model_type, bool add_McCormick_Constraints){
    if(c.is_constant() || c.is_linear()){
        return c;
    }
    if(c.is_nonlinear() || c.is_polynomial()){
        throw invalid_argument("lift can only be called on quadratic constraints for now");
    }
    bool soc_eq = c.check_soc() && c.is_eq();
    
    Constraint<type> lifted(c._name+"_lifted");
    /* Add constant part */
    if (!c.get_cst()->is_zero()) {
        if (c.get_cst()->is_number()) {
            auto f_cst = static_pointer_cast<constant<type>>(c.get_cst());
            lifted += *f_cst;
        }
        else if (c.get_cst()->is_param()) {
            auto f_cst = static_pointer_cast<param<type>>(c.get_cst());
            lifted += *f_cst;
        }
        else if (c.get_cst()->is_function()) {
            auto f_cst = static_pointer_cast<func<type>>(c.get_cst());
            lifted += *f_cst;
            if(lifted._cst->is_function())
                lifted.embed(*static_pointer_cast<func<type>>(lifted._cst));
        }
        else{
            throw invalid_argument("In Model<type>::lift(), constant with unrecognized type!");
        }
    }
    /* Linear part */
    for (auto &pair:*c._lterms) {
        lifted.insert(pair.second);
    }
    /* Quadratic part */
    for (auto &pair:*c._qterms) {
        bool lift_sign; /* if C is an equality SOC, only add one side relaxation since the other side is convex */
        auto qterm = pair.second;
        lift_sign = (qterm._sign ^ qterm._coef->is_negative());
        if (c.func<type>::is_concave())
        {
            DebugOff("Changing the sign of the lifted variable." << endl);
            lift_sign = !lift_sign;
        }
        else{
            DebugOff("Keeping the sign of the lifted variable same." << endl);
        }
        lterm lt;
        lt._sign = qterm._sign;
        auto x1_ptr = static_pointer_cast<var<type>>(qterm._p->first);
        auto x2_ptr = static_pointer_cast<var<type>>(qterm._p->second);
        /* Use the same order for lifting */
        if(x1_ptr->_indices->_keys->at(0) > x2_ptr->_indices->_keys->at(0)){
            x2_ptr = static_pointer_cast<var<type>>(qterm._p->first);
            x1_ptr = static_pointer_cast<var<type>>(qterm._p->second);
            DebugOff("x1 name "<< x1._name << endl);
            DebugOff("x2 name "<< x2._name << endl);
        }
        auto x1 = *x1_ptr;
        auto x2 = *x2_ptr;
        string lifted_name;
        indices ids;
        bool is_square = (x1_ptr==x2_ptr || x1==x2);
        if(is_square){
            lifted_name = "Lift("+x1.get_name(true,true)+"^2)";
            ids = *x1._indices;
            ids.set_name(x1._name);
        }
        else {
            lifted_name = "Lift("+x1.get_name(true,true)+";"+x2.get_name(true,true)+")";
            ids = combine(*x1._indices,*x2._indices);
            ids.set_name(x1._name+";"+x2._name);
        }
        auto x1_ids = *x1._indices;
        auto x2_ids = *x2._indices;
        auto flat_ids = ids;
        if(ids.is_matrix_indexed()){/* Flatten matrix indexed sets */
            flat_ids.flatten();
            x1_ids.flatten();
            x2_ids.flatten();
        }
        auto unique_ids = flat_ids.get_unique_keys(); /* In case of an indexed variable, keep the unique keys only */
        if(unique_ids.size()!=flat_ids.size()){/* If some keys are repeated, remove them from the refs of o1 and o2 */
            auto keep_refs = flat_ids.get_unique_refs();
            x1_ids.filter_refs(keep_refs);
            x2_ids.filter_refs(keep_refs);
        }
        auto name1 = x1.get_name(true,true);
        auto name2 = x2.get_name(true,true);
        func<double> lb, ub;/* lower and upper bounds for lifted variable */
        
        auto x1_unindexed = this->get_var<type>(x1_ptr->get_name(true,true));
        auto x2_unindexed = this->get_var<type>(x2_ptr->get_name(true,true));
        
        
        auto it = _vars_name.find(lifted_name);
        if(it==_vars_name.end()){/* New lifted variable */
            /* Create the lifted variable with correct lower and upper bounds */
            /* Get the unindexed variables from the model */
            /* Get the */
            auto lb1 = x1_unindexed.get_lb().in(x1_ids);auto ub1 = x1_unindexed.get_ub().in(x1_ids);
            if(is_square){
                func<double> prod_b1 = pow(lb1,2);
                func<double> prod_b2 = pow(ub1,2);
                func<double> prod_b3 = lb1*ub1;
                lb = gravity::max(gravity::min(gravity::min(prod_b1,prod_b2).in(unique_ids), prod_b3).in(unique_ids), func<type>());
                ub = gravity::max(prod_b1,prod_b2);/* max(lb^2,ub^2) */
            }
            else {
                auto lb2 = x2_unindexed.get_lb().in(x2_ids);auto ub2 = x2_unindexed.get_ub().in(x2_ids);
                func<double> prod_b1 = lb1*lb2;
                func<double> prod_b2 = lb1*ub2;
                func<double> prod_b3 = ub1*lb2;
                func<double> prod_b4 = ub1*ub2;
                lb = gravity::min(gravity::min(prod_b1,prod_b2).in(unique_ids),gravity::min(prod_b3,prod_b4).in(unique_ids));
                ub = gravity::max(gravity::max(prod_b1,prod_b2).in(unique_ids),gravity::max(prod_b3,prod_b4).in(unique_ids));
            }
            var<type> vlift(lifted_name, lb, ub);
            vlift._lift = true;
            vlift._original_vars.push_back(make_shared<var<type>>(x1.in(x1_ids)));
            vlift._original_vars.push_back(make_shared<var<type>>(x2.in(x2_ids)));
            add(vlift.in(unique_ids));
            if(soc_eq){
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
            lt._p = make_shared<var<type>>(vlift.in(ids));
            if(add_McCormick_Constraints){
                add_McCormick(pair.first, vlift.in(unique_ids), x1.in(x1_ids), x2.in(x2_ids));
            }
        }
        else {/* Lifted variable already added to model*/
            auto vlift = static_pointer_cast<var<type>>(it->second);
            auto new_ids = unique_ids.get_diff_refs(*vlift->_indices);/* Get the new indices */
            unique_ids.filter_refs(new_ids);/* Only keep new indices */
            vlift->_lb->merge_vars(*vlift->_ub);/* Make sure the parameters are shared among both functions */
            if(is_square){
                x1_ids.filter_refs(new_ids);
                auto x_lb = vlift->get_square_lb();
                auto x_ub = vlift->get_square_ub();
                x_lb->_indices->add_refs(x1_ids);
                x_ub->_indices->add_refs(x1_ids);
                x_lb->update_dim();
                x_ub->update_dim();
                vlift->_original_vars[0]->_indices->add_refs(x1_ids);
                vlift->_original_vars[0]->_lb->index_in(*vlift->_original_vars[0]->_indices);
                vlift->_original_vars[0]->_ub->index_in(*vlift->_original_vars[0]->_indices);
            }
            else{
                x1_ids.filter_refs(new_ids);
                x2_ids.filter_refs(new_ids);
                auto x1_lb = vlift->get_bilinear_lb1();
                auto x1_ub = vlift->get_bilinear_ub1();
                x1_lb->_indices->add_refs(x1_ids);
                x1_ub->_indices->add_refs(x1_ids);
                x1_lb->update_dim();
                x1_ub->update_dim();
                auto x2_lb = vlift->get_bilinear_lb2();
                auto x2_ub = vlift->get_bilinear_ub2();
                x2_lb->_indices->add_refs(x2_ids);
                x2_ub->_indices->add_refs(x2_ids);
                x2_lb->update_dim();
                x2_ub->update_dim();
                vlift->_original_vars[0]->_indices->add_refs(x1_ids);
                vlift->_original_vars[1]->_indices->add_refs(x2_ids);
                vlift->_original_vars[0]->_lb->index_in(*vlift->_original_vars[0]->_indices);
                vlift->_original_vars[0]->_ub->index_in(*vlift->_original_vars[0]->_indices);
                vlift->_original_vars[1]->_lb->index_in(*vlift->_original_vars[1]->_indices);
                vlift->_original_vars[1]->_ub->index_in(*vlift->_original_vars[1]->_indices);
            }
            auto added = vlift->_indices->add(unique_ids);
            if(!added.empty()){
                vlift->get_lb()._indices->add_refs(unique_ids);
                vlift->get_ub()._indices->add_refs(unique_ids);
                vlift->_lb->index_in(*vlift->_indices);
                vlift->_lb->uneval();
                vlift->_lb->eval_all();
                vlift->_ub->index_in(*vlift->_indices);
                vlift->_ub->uneval();
                vlift->_ub->eval_all();
                /* Create the lifted variable with correct lower and upper bounds */
                /* Get the unindexed variables from the model */

//                auto lb1 = x1_unindexed.get_lb().in(*vlift->_original_vars[0]->_indices);auto ub1 = x1_unindexed.get_ub().in(*vlift->_original_vars[0]->_indices);
//                if(is_square){
//                    func<double> prod_b1 = pow(lb1,2);
//                    func<double> prod_b2 = pow(ub1,2);
//                    func<double> prod_b3 = lb1*ub1;
//                    lb = gravity::max(gravity::min(gravity::min(prod_b1,prod_b2).in(*vlift->_indices), prod_b3).in(*vlift->_indices), func<type>());
//                    ub = gravity::max(prod_b1,prod_b2);/* max(lb^2,ub^2) */
//                }
//                else {
//                    auto lb2 = x2_unindexed.get_lb().in(*vlift->_original_vars[1]->_indices);auto ub2 = x2_unindexed.get_ub().in(*vlift->_original_vars[1]->_indices);
//                    func<double> prod_b1 = lb1*lb2;
//                    func<double> prod_b2 = lb1*ub2;
//                    func<double> prod_b3 = ub1*lb2;
//                    func<double> prod_b4 = ub1*ub2;
//                    lb = gravity::min(gravity::min(prod_b1,prod_b2).in(*vlift->_indices),gravity::min(prod_b3,prod_b4).in(*vlift->_indices));
//                    ub = gravity::max(gravity::max(prod_b1,prod_b2).in(*vlift->_indices),gravity::max(prod_b3,prod_b4).in(*vlift->_indices));
//                }
//                lb.uneval();ub.uneval();
//                lb.eval_all();ub.eval_all();
//                auto lb_param = vlift->_lb->_params->begin()->second.first;
//                auto ub_param = vlift->_ub->_params->begin()->second.first;
//                static_pointer_cast<param<T>>(lb_param)->_val = lb._val;
//                static_pointer_cast<param<T>>(ub_param)->_val = ub._val;
//                vlift->_lb->_val = lb._val;
//                vlift->_ub->_val = ub._val;
//                vlift->_lb->uneval();vlift->_ub->uneval();
//                vlift->_lb->eval_all();vlift->_ub->eval_all();
                vlift->update_dim();
                reindex_vars();
            }
            lt._p = make_shared<var<type>>(vlift->in(ids));
            if(add_McCormick_Constraints){
                add_McCormick(pair.first, vlift->in(unique_ids), x1.in(x1_ids), x2.in(x2_ids));
            }
        }
        lt._coef = qterm._coef;
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
int Model<type>::readNL(const string& fname){    
    mp::Problem p;
    mp::ReadNLFile(fname, p);
    auto nb_vars = p.num_vars();
    auto nb_cstr = p.num_algebraic_cons();
    DebugOn("The number of variables is " << nb_vars << endl);
    DebugOn("The number of constraints is " << nb_cstr << endl);
    indices C("C"), I("I"), LinConstr("LinConstr"), QuadConstr("QuadConstr"), NonLinConstr("NonLinConstr");
    int nb_cont = 0;
    int nb_int = 0;
    int nb_other = 0;
    vector<int> C_ids,I_ids;
    for (const auto v: p.vars()) {
        if(v.type()==mp::var::CONTINUOUS){
            nb_cont++;
            C.insert(to_string(v.index()));
            C_ids.push_back(v.index());
        }
        else if(v.type()==mp::var::INTEGER){
            I.insert(to_string(v.index()));
            I_ids.push_back(v.index());
            nb_int++;
        }
        else{
            throw invalid_argument("Unrecognised variable type, can only be continuous or integer");
        }
    }
    DebugOn("Number of continuous variables = " << nb_cont << endl);
    DebugOn("Number of integer variables = " << nb_int << endl);
    
    param<> x_ub("x-ub"), x_lb("x-lb");
    param<int> y_ub("y-ub"), y_lb("y-lb");
    x_ub.in(C);x_lb.in(C);
    y_ub.in(I);y_lb.in(I);
    for (int i = 0; i<C.size(); i++) {
        x_lb.set_val(i, std::max(-1e12,p.var(C_ids[i]).lb()));
        x_ub.set_val(i, std::min(1e12,p.var(C_ids[i]).ub()));
    }
    for (int i = 0; i<I.size(); i++) {
        y_lb.set_val(i, std::max(-1e12,p.var(I_ids[i]).lb()));
        y_ub.set_val(i, std::min(1e12,p.var(I_ids[i]).ub()));
    }
    var<> x("x", x_lb, x_ub);
    var<int> y("y", y_lb, y_ub);

    param<> rhs("rhs");
    int nb_lin = 0;
    int nb_nonlin = 0;
    int index = 0;
    _name = fname;

    if(!C.empty())
        add(x.in(C));
    if(!I.empty()){
        add(y.in(I));
        replace_integers();
    }

    MPConverter converter(*this);
    map<int,vector<int>> constr_sparsity;
    vector<int> C_lin, C_nonlin, C_quad;
    
    int num_objs = p.num_objs();
    if(num_objs>=1){
        if(num_objs>=2){
            DebugOn("Gravity currently supports only one objective, will only add the first one");
        }
        mp::Problem::Objective obj = p.obj(0);
        mp::obj::Type main_obj_type = p.obj(0).type();
        func<> objective;
        auto lexpr = obj.linear_expr();
        for (const auto term: lexpr){
            auto coef = term.coef();
            auto var_id = term.var_index();
            if(coef!=0)
                objective += coef*converter.get_cont_int_var(var_id);
        }
        auto nl_expr = obj.nonlinear_expr();
        if (nl_expr){
            auto expr = converter.Visit(nl_expr);
            objective += expr;
        }
        auto sense = main_obj_type == mp::obj::MIN;
        if(sense){
            min(objective);
        }
        else {
            max(objective);
        }
    }
    
    for (const auto con: p.algebraic_cons()) {
        auto lexpr = con.linear_expr();
        auto nl_expr = con.nonlinear_expr();
        if (nl_expr){
            auto expr = converter.Visit(nl_expr);
            nb_nonlin++;
            C_nonlin.push_back(index);
            NonLinConstr.insert(to_string(index));
            for (const auto term: lexpr){
                auto coef = term.coef();
                auto var_id = term.var_index();
                if(coef!=0)
                    expr += coef*converter.get_cont_int_var(var_id);
            }
            
            auto c_lb = con.lb();
            auto c_ub = con.ub();
            if(c_lb==c_ub){
                Constraint<> c("NL_C_eq_"+to_string(index));
                c += expr;
                add(c.in(range(1,1)) == c_lb);
            }
            else {
                if(c_lb>numeric_limits<double>::lowest()){
                    Constraint<> c("NL_C_geq_"+to_string(index));
                    c += expr;
                    add(c.in(range(1,1)) >= c_lb);
                }
                if(c_ub<numeric_limits<double>::max()){
                    Constraint<> c("NL_C_leq_"+to_string(index));
                    c += expr;
                    add(c.in(range(1,1)) <= c_ub);
                }
            }
        }
        else{
            int nb_terms = 0;
            func<> expr;
            for (const auto term: lexpr){
                auto coef = term.coef();
                auto var_id = term.var_index();
                if(coef!=0){
                    expr += coef*converter.get_cont_int_var(var_id);
                    nb_terms++;
                }
            }
            auto c_lb = con.lb();
            auto c_ub = con.ub();
            if(nb_terms==1 && c_lb!=c_ub){/* this is just a bound constraint */
                const auto term = *lexpr.begin();
                auto coef = term.coef();
                auto var_id = term.var_index();
                auto v = converter.get_cont_int_var(var_id);
                if(coef != 0 && c_lb>numeric_limits<double>::lowest()){
                    if(v._is_relaxed){/* an integer */
                        v._lb->_val->at(v.get_id_inst()) = c_lb/coef;
                    }
                    else {
                        v._lb->_val->at(v.get_id_inst()) = c_lb/coef;
                    }
                }
                if(coef != 0 && c_ub<numeric_limits<double>::max()){
                    if(v._is_relaxed){/* an integer */
                        v._ub->_val->at(v.get_id_inst()) = c_ub/coef;
                    }
                    else {
                        v._ub->_val->at(v.get_id_inst()) = c_ub/coef;
                    }
                }
            }
            else{
                constr_sparsity[nb_terms].push_back(index);
                if(c_lb==c_ub){
                    Constraint<> c("Lin_C_eq_"+to_string(index));
                    c += expr;
                    add(c.in(range(1,1)) == c_lb);
                }
                else {
                    if(c_lb>numeric_limits<double>::lowest()){
                        Constraint<> c("Lin_C_geq_"+to_string(index));
                        c += expr;
                        add(c.in(range(1,1)) >= c_lb);
                    }
                    if(c_ub<numeric_limits<double>::max()){
                        Constraint<> c("Lin_C_leq_"+to_string(index));
                        c += expr;
                        add(c.in(range(1,1)) <= c_ub);
                    }
                }
                LinConstr.insert(to_string(index));
                C_lin.push_back(index);
                nb_lin++;
            }
        }
        index++;
    }
    DebugOn("Number of linear constraints = " << nb_lin << endl);
    DebugOn("Number of non linear constraints = " << nb_nonlin << endl);
    DebugOn("Number of sparsity degrees for linear constraints = " << constr_sparsity.size() << endl);

    
    return 0;
}
    
template <typename type>
template<typename T>
void Model<type>::populate_original_interval(map<string, bool>& fixed_point, map<string, double>& ub_original,map<string, double>& lb_original,map<string, double>& interval_original,map<string, double>& interval_new, int& count_skip, int& count_var){
    var<> v;
    std::string var_key, key_lb, key_ub;
    for(auto &it:_vars)
    {
        string vname=(*it.second)._name;
        v=this->template get_var<double>(vname);
        auto v_keys=v.get_keys();
        auto v_key_map=v.get_keys_map();
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
            obj_ub = (obj - ub);
           // obj_ub = (obj - ub)*1000/ub_scale_value;
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
                close= true;
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
    int worker_id, nb_workers;
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
    string vname, var_key, mname, cut_type="allvar";
    string dir_array[2]={"LB", "UB"};
    var<> v;
    bool close=false, terminate=false, xb_true=true, alg_batch_reset=true;
    const double fixed_tol_abs=range_tol, fixed_tol_rel=range_tol, zero_tol=1e-6, obbt_subproblem_tol=1e-6;
    int iter=0, fail=0, count_var=0, count_skip=0, nb_init_refine=nb_refine;
    double solver_time =0, gapnl,gap, gaplin=-999, sum=0, avg=0, active_root_tol=lb_solver_tol, active_tol=1e-6;
    double lower_bound_nonlin_init = numeric_limits<double>::min(), lower_bound_init = numeric_limits<double>::min(), upper_bound = 0, lower_bound = numeric_limits<double>::min(), lower_bound_old;
    map<string,int> old_map;
    if(this->_status==0){
        upper_bound=this->get_obj_val();
        if(relaxed_model->_status==0){
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
                if(obbt_model->_status==0){
                    /*Initialize fixed point, interval original and new, bounds original*/
                    obbt_model->populate_original_interval(fixed_point, ub_original,lb_original, interval_original,interval_new,  count_skip, count_var);
                    solver_time= get_wall_time()-solver_time_start;
                    /*Create nb_threads copy of obbt_models*/
                    obbt_model->create_batch_models(batch_models, nb_threads, ub_scale_value);
                    if(linearize){
                        if(initialize_primal && lb_solver_type==gurobi){
                            initialize_basis_vectors(lb_solver_type, vbasis,cbasis,vrbasis,crbasis,nb_threads);
                        }
                        if(!alg_batch_reset){
                        obbt_model->batch_models_obj_lb_constr(batch_models, nb_threads, lower_bound, lower_bound_old,lower_bound_nonlin_init, upper_bound, ub_scale_value);
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
                                            auto viol= run_MPI_new(objective_models, sol_obj, sol_status, batch_models, relaxed_model, interior_model, cut_type, active_tol, lb_solver_type, obbt_subproblem_tol, nb_threads, "ma27", 10000, 600, linearize, nb_refine, old_map, vbasis, cbasis, initialize_primal);
#else
                                            auto viol= run_parallel_new(objective_models, sol_obj, sol_status, batch_models, relaxed_model, interior_model, cut_type, active_tol, lb_solver_type, obbt_subproblem_tol, nb_threads, "ma27", 10000, 600, linearize, nb_refine, vbasis, cbasis, initialize_primal);
#endif
                                            auto b=obbt_model->obbt_update_bounds( objective_models, sol_obj,  sol_status, batch_models,  fixed_point, interval_original, interval_new, ub_original, lb_original, terminate, fail, range_tol, fixed_tol_abs, fixed_tol_rel, zero_tol, iter);
                                            sol_status.clear();
                                            sol_obj.clear();
                                            objective_models.clear();
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
                            if(!alg_batch_reset){
                            obbt_model->batch_models_obj_lb_constr(batch_models, nb_threads, lower_bound,lower_bound_old,lower_bound_nonlin_init, upper_bound, ub_scale_value);
                            lower_bound_old=lower_bound;
                            }
			    else{
                                batch_models.clear();
                                obbt_model->create_batch_models(batch_models, nb_threads, ub_scale_value);
                                if(initialize_primal && lb_solver_type==gurobi){
                                     initialize_basis_vectors(lb_solver_type, vbasis,cbasis,vrbasis,crbasis,nb_threads);
                               }
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
                }
                else{
                    DebugOn("Initial lower bounding problem not solved to optimality, cannot compute initial gap"<<endl);
                    lower_bound=numeric_limits<double>::min();
                }
            }
        }
        else{
            DebugOn("Lower bounding problem not solved to optimality, cannot compute initial gap"<<endl);
            lower_bound=numeric_limits<double>::min();
        }
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
template void Model<double>::populate_original_interval(map<string, bool>& fixed_point, map<string, double>& ub_original,map<string, double>& lb_original,map<string, double>& interval_original,map<string, double>& interval_new, int& count_skip, int& count_var);
template double Model<double>::populate_final_interval_gap(const shared_ptr<Model<double>>& obbt_model, const map<string, double>& interval_original, map<string, double>& interval_new, double& sum, bool& xb_true, const double zero_tol, int count_var);
template void Model<double>::create_batch_models(vector<shared_ptr<Model<double>>>& batch_models, int nb_threads, double ub_scale_value);
template void Model<double>::compute_iter_gap(double& gap, double& active_tol, bool& terminate, bool linearize, int iter, shared_ptr<Model<double>>& obbt_model, const Model<double>& interior_model, SolverType lb_solver_type, int nb_root_refine, const double upper_bound, double& lower_bound, const double ub_scale_value, double lb_solver_tol, double& active_root_tol, int& oacuts, const double abs_tol, const double rel_tol, const double zero_tol, string lin_solver, int max_iter, int max_time, vector<double>& vrbasis, std::map<string,double>& crbasis, bool initialize_primal);
template void Model<double>::batch_models_obj_lb_constr(vector<shared_ptr<Model<double>>>& batch_models, int nb_threads, double lower_bound_lin, double lower_bound_old, double lower_bound_nonlin_init, double upper_bound, double ub_scale_value);


}





