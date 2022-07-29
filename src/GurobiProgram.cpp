#include <gravity/GurobiProgram.h>
#include <thread>         // std::this_thread::sleep_for
#include <chrono>         // std::chrono::seconds

class cuts: public GRBCallback
{
public:
    vector<GRBVar> vars;
    int nb_vars;
    Model<>* m;
    double *x;
    vector<double> cont_x, int_x;
    double objc=999;
    int soc_viol, soc_found,soc_added,det_viol, det_found, det_added;
    int soc_viol_user=0, soc_found_user=0,soc_added_user=0,det_viol_user=0, det_found_user=0, det_added_user=0;
    Model<> interior;
    cuts(const vector<GRBVar>& _grb_vars, int xn, Model<>* mod, Model<>& mod_int, int& soc_violn, int& soc_foundn, int& soc_addedn, int& det_violn, int& det_foundn, int& det_addedn) {
        vars = _grb_vars;
        nb_vars = xn;
        x = new double[nb_vars];
        m = mod;
        cont_x.resize(nb_vars);
        int_x.resize(nb_vars);
        soc_viol=soc_violn;
        soc_found=soc_foundn;
        soc_added=soc_addedn;
        det_viol=det_violn;
        det_found=det_foundn;
        det_added=det_addedn;
        interior=mod_int;
    }
    ~cuts(){
        DebugOff("soc_viol "<<soc_viol_user<<endl);
        DebugOff("soc_found "<<soc_found_user<<endl);
        DebugOff("soc_added "<<soc_added_user<<endl);
        DebugOff("det_viol "<<det_viol_user<<endl);
        DebugOff("det_found "<<det_found_user<<endl);
        DebugOff("det_added "<<det_added_user<<endl);
     
        delete [] x;
    }
protected:
    void callback() {
        try {
            bool incumbent=true;
            bool mipnode=true;
            if(incumbent){
                if (where == GRB_CB_MIPSOL) {
                    int i,j;
                    x=getSolution(vars.data(),nb_vars);
                    for(i=0;i<nb_vars;i++){
                        int_x[i] = x[i];
                    }
                    m->set_solution(int_x);
                    auto res=m->cutting_planes_solution(interior, 1e-9, soc_viol, soc_found, soc_added, det_viol, det_found, det_added);
                    if(res.size()>=1){
                        for(i=0;i<res.size();i++){
                            GRBLinExpr expr = 0;
                            for(j=0;j<res[i].size()-1;j+=2){
                                int c=res[i][j];
                                expr += res[i][j+1]*vars[c];
                            }
                            expr+=res[i][j];
                            addLazy(expr, GRB_LESS_EQUAL, 0);
                        }
                    }
                }
            }
            if(mipnode){
                if (where == GRB_CB_MIPNODE){
                    int stat=getIntInfo(GRB_CB_MIPNODE_STATUS);
                    if(stat==2){
                        DebugOff(getIntInfo(GRB_CB_MIPNODE_STATUS)<<endl);
                        int nct=getDoubleInfo(GRB_CB_MIPNODE_NODCNT);
                        if(nct%100==0){
                            double obj=getDoubleInfo(GRB_CB_MIPNODE_OBJBST);
                            double obj1=getDoubleInfo(GRB_CB_MIPNODE_OBJBND);
                            DebugOff(obj<<"\t"<<obj1<<"\t"<<endl);
                            int i,j;
                            m->get_solution(int_x);
                            x=getNodeRel(vars.data(),nb_vars);
                            for(i=0;i<nb_vars;i++){
                                cont_x[i] = x[i];
                            }
                            m->set_solution(cont_x);
                            auto res=m->cutting_planes_solution(interior, 1e-9,soc_viol_user, soc_found_user,soc_added_user,det_viol_user, det_found_user, det_added_user);
                            if(res.size()>=1){
                                for(i=0;i<res.size();i++){
                                    GRBLinExpr expr = 0;
                                    for(j=0;j<res[i].size()-1;j+=2){
                                        int c=res[i][j];
                                        expr += res[i][j+1]*vars[c];
                                    }
                                    expr+=res[i][j];
                                    addCut(expr, GRB_LESS_EQUAL, 0);
                                }
                            }
                        }
                        m->set_solution(int_x);
                    }
                    
                }
            }
        } catch (GRBException e) {
            cout << "Error number: " << e.getErrorCode() << endl;
            cout << e.getMessage() << endl;
        } catch (...) {
            cout << "Error during callback" << endl;
        }
    }
};

//GurobiProgram::GurobiProgram(){
////    model = m;
//    grb_env = new GRBEnv();
////    grb_env->set(GRB_IntParam_Presolve,0);
//    //grb_env->set(GRB_DoubleParam_NodeLimit,1);
//    grb_env->set(GRB_DoubleParam_TimeLimit,7200);
////    grb_env->set(GRB_DoubleParam_MIPGap,0.01);
//    grb_env->set(GRB_IntParam_Threads,1);
////       grb_env->set(GRB_IntParam_Presolve,0);
////      grb_env->set(GRB_IntParam_NumericFocus,3);
////     grb_env->set(GRB_IntParam_NonConvex,2);
//    //grb_env->set(GRB_DoubleParam_FeasibilityTol, 1E-6);
//     grb_env->set(GRB_DoubleParam_OptimalityTol, 1E-6);
//
//    grb_env->set(GRB_IntParam_OutputFlag,1);
////    grb_mod = new GRBModel(*grb_env);
//    grb_mod = NULL;
//}


GurobiProgram::GurobiProgram(Model<>* m) {
    bool found_token = false;
    while (!found_token) {
        try{
            grb_env = new GRBEnv();
        //    grb_env->set(GRB_IntParam_Presolve,0);
            //grb_env->set(GRB_DoubleParam_NodeLimit,1);
            //grb_env->set(GRB_DoubleParam_TimeLimit,7200);
            //    grb_env->set(GRB_DoubleParam_MIPGap,0.01);
            //grb_env->set(GRB_IntParam_Threads,1);
        //    grb_env->set(GRB_IntParam_Presolve,0);
        //    grb_env->set(GRB_IntParam_NumericFocus,3);
            grb_env->set(GRB_IntParam_NonConvex,2);
         //   grb_env->set(GRB_DoubleParam_FeasibilityTol, 1E-6);
         //   grb_env->set(GRB_DoubleParam_OptimalityTol, 1E-6);
            
            //grb_env->set(GRB_IntParam_OutputFlag,0);
            grb_mod = new GRBModel(*grb_env);
            //grb_mod->set(GRB_IntParam_Method, 0);
            found_token = true;
        }
        catch(GRBException e) {
//            cerr << "\nWas not able to create Gurobi environment or model, Error code = " << e.getErrorCode() << endl;
//            cerr << e.getMessage() << endl;
//            exit(-1);
            found_token = false;
            this_thread::sleep_for (chrono::seconds(1));
        }
    }
    //    grb_env->set(GRB_IntParam_OutputFlag,2);
    _model = m;
    m->fill_in_maps();
    m->compute_funcs();
    //DebugOn("found token "<<endl);
}

//GurobiProgram::GurobiProgram(const shared_ptr<Model<>>& m) {
//    grb_env = new GRBEnv();
//    //grb_env->set(GRB_IntParam_Presolve,0);
//    //grb_env->set(GRB_DoubleParam_NodeLimit,1);
//    grb_env->set(GRB_DoubleParam_TimeLimit,7200);
//    //    grb_env->set(GRB_DoubleParam_MIPGap,0.01);
//        grb_env->set(GRB_IntParam_Threads,1);
////    grb_env->set(GRB_IntParam_Presolve,0);
////    grb_env->set(GRB_IntParam_NumericFocus,3);
////    grb_env->set(GRB_IntParam_NonConvex,2);
//    //grb_env->set(GRB_DoubleParam_FeasibilityTol, 1E-6);
//    grb_env->set(GRB_DoubleParam_OptimalityTol, 1E-6);
//
//    grb_env->set(GRB_IntParam_OutputFlag,1);
//    grb_mod = new GRBModel(*grb_env);
////    grb_env->set(GRB_IntParam_ƒFlag,2);
////    _model = m;
//    m->fill_in_maps();
//    m->compute_funcs();
//}


GurobiProgram::~GurobiProgram() {
//    for (auto p : _grb_vars) delete p.second;
    if (grb_mod) {
        delete grb_mod;
        grb_mod = nullptr;
        DebugOff("deleted model"<<endl);
    }
    if (grb_env){
        delete grb_env;
        grb_env = nullptr;
        DebugOff("deleted environment"<<endl);
    }
   
}

void GurobiProgram::reset_model(){
    if (grb_mod) delete grb_mod;
    _grb_vars.clear();
//    grb_env->set(GRB_IntParam_OutputFlag,_output);
    grb_mod = new GRBModel(*grb_env);
}

bool GurobiProgram::solve(int output, bool relax, double tol, double mipgap, bool gurobi_crossover){
    //cout << "\n Presolve = " << grb_env->get(GRB_IntParam_Presolve) << endl;
//    print_constraints();
 //   if (relax) relax_model();
//    relax_model();
    //grb_mod->set(GRB_DoubleParam_BarConvTol, 1e-8);
      //grb_mod->set(GRB_IntParam_ScaleFlag, 1);
      grb_mod->set(GRB_DoubleParam_FeasibilityTol, tol);
      grb_mod->set(GRB_DoubleParam_OptimalityTol, tol);
    //grb_mod->set(GRB_DoubleParam_MIPGap, mipgap);
      grb_mod->set(GRB_IntParam_Threads, 1);
      grb_mod->set(GRB_DoubleParam_TimeLimit, 3600.0);
    ///grb_mod->set(GRB_IntParam_NumericFocus, 1);
//    if(grb_first_run){
//        grb_mod->set(GRB_IntParam_Method, 0);
//    }
//    else{
//        grb_mod->set(GRB_IntParam_Method, 1);
//    }
    
//    if(!gurobi_crossover){
//        grb_mod->set(GRB_IntParam_Crossover, 0);
//    }
    grb_mod->set(GRB_IntParam_OutputFlag, 1);
    grb_mod->update();
        /*int n=grb_mod->get(GRB_IntAttr_NumVars);
    auto lin=_model->buildOA();
        auto interior=lin->add_outer_app_solution(*_model);
        int soc_viol=0,soc_found=0,soc_added=0,det_viol=0,det_found=0,det_added=0;
        DebugOn("solved interior "<<endl);
        cuts cb(_grb_vars, n, _model, interior, soc_viol,soc_found,soc_added,det_viol,det_found,det_added);
        grb_mod->setCallback(&cb);*/
//    warm_start(); // No need to reset variables if Gurobi model has not changed.
    grb_mod->write("gurobiprint.lp");
    try{
        grb_mod->optimize();
    }
    catch(GRBException e) {
        cerr << "\nWas not able to optimize Gurobi model, Error code = " << e.getErrorCode() << endl;
        cerr << e.getMessage() << endl;
        exit(-1);
    }
    //cout<<"Status "<<grb_mod->get(GRB_IntAttr_Status)<<endl;

    // cout<<"BoundVio "<<grb_mod->get(GRB_DoubleAttr_BoundVio)<<endl;
    // cout<<"BoundSVio "<<grb_mod->get(GRB_DoubleAttr_BoundSVio)<<endl;
    // cout<<"ConstrVio "<<grb_mod->get(GRB_DoubleAttr_ConstrVio)<<endl;
    // cout<<"ConstrSVio "<<grb_mod->get(GRB_DoubleAttr_ConstrSVio)<<endl;
    // cout<<"ConstrResidual "<<grb_mod->get(GRB_DoubleAttr_ConstrResidual)<<endl;
    // cout<<"ConstrSResidual "<<grb_mod->get(GRB_DoubleAttr_ConstrSResidual)<<endl;
 
    //grb_mod->write("~/mod.mps");
    grb_first_run=false;
    if (grb_mod->get(GRB_IntAttr_Status) != 2) {
        cerr << "\nModel has not been solved to optimality, error code = " << grb_mod->get(GRB_IntAttr_Status) << endl;
        return false;
    }
    update_solution();
//    GRBVar* gvars = grb_mod->getVars();
//    for(int i = 0; i < grb_mod->get(GRB_IntAttr_NumVars); ++i) {
////        cout << gvars[i].get(GRB_StringAttr_VarName) << "  " << gvars[i].get(GRB_DoubleAttr_X) << endl;
//        if (gvars[i].get(GRB_CharAttr_VType)==GRB_BINARY) {
//            cout << gvars[i].get(GRB_StringAttr_VarName) << "  ";
//            cout << gvars[i].get(GRB_DoubleAttr_X);
//            cout << "\n";
//        }
//    }
    _model->_obj->set_val(grb_mod->get(GRB_DoubleAttr_ObjVal));
//    cout << "\n***** Optimal Objective = " << _model->get_obj_val() << " *****\n";
    if (grb_mod->get(GRB_IntAttr_IsMIP)) {
        cout.setf(ios::fixed);
        cout.precision(3);
        cout << "Results: " << grb_mod->get(GRB_DoubleAttr_ObjVal) << " & ";
        cout.precision(4);
        cout << (grb_mod->get(GRB_DoubleAttr_MIPGap))*100 << "% & ";
        cout.precision(0);
        cout << grb_mod->get(GRB_DoubleAttr_NodeCount) << " & ";
        cout.precision(2);
        cout << grb_mod->get(GRB_DoubleAttr_Runtime) << " & " << endl;
    }
//    delete[] gvars;
    return true;
}

void GurobiProgram::prepare_model(){
    _model->fill_in_maps();
    _model->compute_funcs();
//    DebugOn("going to fill in vmap"<<endl);
    fill_in_grb_vmap();
//    DebugOn("going to fill in cmap"<<endl);
    create_grb_constraints();
//    DebugOn("going to fill in omap"<<endl);
    set_grb_objective();

//    grb_mod->write("gurobiprint.lp");

//    print_constraints();
}
void GurobiProgram::initialize_basis(const std::vector<double>& vbasis, const std::vector<double>& cbasis){
    grb_mod->update();
    int nv=grb_mod->get(GRB_IntAttr_NumVars);
    int nc=grb_mod->get(GRB_IntAttr_NumConstrs);
    int count=0;
    for(auto &gv:_grb_vars){
        gv.set(GRB_IntAttr_VBasis, vbasis[count++]);
    }

    GRBConstr* gcons= grb_mod->getConstrs();
  
//DebugOn("nc "<<nc<<endl);
//DebugOn("cbasis "<<cbasis.size()<<endl);
    for(auto i=0;i<cbasis.size();i++){
        gcons[i].set(GRB_IntAttr_CBasis, cbasis[i]);
    }
  for(auto i=cbasis.size();i<nc;i++){
        gcons[i].set(GRB_IntAttr_CBasis, 0);
    }
delete[] gcons;
gcons=nullptr;
}

void GurobiProgram::initialize_pstart(const std::vector<double>& pstart, const std::map<std::string, double>& dstart){
    grb_mod->update();
    int nv=grb_mod->get(GRB_IntAttr_NumVars);
    int nc=grb_mod->get(GRB_IntAttr_NumConstrs);
    if(pstart.size()==nv){
            int count=0;
    for(auto &gv:_grb_vars){
        gv.set(GRB_DoubleAttr_PStart, pstart[count++]);
    }

    GRBConstr* gcons= grb_mod->getConstrs();
  
//DebugOn("nc "<<nc<<endl);
//DebugOn("cbasis "<<cbasis.size()<<endl);
    for(auto i=0;i<nc;i++){
        auto cname=gcons[i].get(GRB_StringAttr_ConstrName);
        if(dstart.find(cname)!=dstart.end()){
            gcons[i].set(GRB_DoubleAttr_DStart, dstart.at(cname));
        }
        else{
            gcons[i].set(GRB_DoubleAttr_DStart,0);
        }
    }

delete[] gcons;
    gcons=nullptr;
    }
    //DebugOn(grb_mod->get(GRB_IntAttr_NumConstrs));
}


void GurobiProgram::get_basis(std::vector<double>& vbasis, std::vector<double>& cbasis){
    int nv=grb_mod->get(GRB_IntAttr_NumVars);
    int nc=grb_mod->get(GRB_IntAttr_NumConstrs);
    vbasis.resize(nv);
    cbasis.resize(nc);
    GRBVar* gvars= grb_mod->getVars();
    for(auto i=0;i<nv;i++){
        vbasis[i]=gvars[i].get(GRB_IntAttr_VBasis);
    }

    GRBConstr* gcons= grb_mod->getConstrs();

    for(auto i=0;i<nc;i++){
        cbasis[i]=gcons[i].get(GRB_IntAttr_CBasis);
    }
delete[] gcons;
delete[] gvars;
    gcons=nullptr;
    gvars=nullptr;
}

void GurobiProgram::get_pstart(std::vector<double>& pstart, std::map<std::string, double>& dstart){
    int nv=grb_mod->get(GRB_IntAttr_NumVars);
    int nc=grb_mod->get(GRB_IntAttr_NumConstrs);
    pstart.resize(nv);
    dstart.clear();
    GRBVar* gvars= grb_mod->getVars();
    for(auto i=0;i<nv;i++){
        pstart[i]=gvars[i].get(GRB_DoubleAttr_X);
    }

    GRBConstr* gcons= grb_mod->getConstrs();

    for(auto i=0;i<nc;i++){
        auto cname=gcons[i].get(GRB_StringAttr_ConstrName);
        dstart[cname]=gcons[i].get(GRB_DoubleAttr_Pi);
    }
    delete[] gcons;
    delete[] gvars;
    gcons=nullptr;
    gvars=nullptr;
}

void GurobiProgram::update_model(){
    _model->fill_in_maps();
    _model->compute_funcs();
    fill_in_grb_vmap();
    create_grb_constraints();
    set_grb_objective();
    
}

void GurobiProgram::update_solution(){
    size_t vid, vid_inst;
    GRBVar gvar;
    param_* v;
//    for (auto i = 0; i < _grb_vars.size(); i++) {
//        gvar = _grb_vars.at(i);
//        auto dim = _model->_vars[i]->get_dim();
//        for (auto j = 0; j < _model->_vars[i]->get_dim(); j++) {
//            poly_set_val(j, gvar.get(GRB_DoubleAttr_X), _model->_vars[i]);
//        }
//    }
    for(auto& v_p: _model->_vars)
    {
        v = v_p.second.get();
        auto idx = v->get_id();
        auto dim = v->_dim[0];
        for (auto i = 0; i < dim; i++) {
            auto vid = idx + v->get_id_inst(i);
            gvar = _grb_vars.at(vid);
            v->set_double_val(i,gvar.get(GRB_DoubleAttr_X));
        }
    }
}

void GurobiProgram::warm_start(){
    GRBVar gvar;
    param_* v;
    double value;
    for(auto& v_p: _model->_vars)
    {
        v = v_p.second.get();
        auto idx = v->get_id();
        auto dim = v->_dim[0];
        for (auto i = 0; i < dim; i++) {
            auto vid = idx + i;
            v->get_double_val(i,value);
            _grb_vars.at(vid).set(GRB_DoubleAttr_PStart, value);
        }
    }
}

void GurobiProgram::relax_model(){
    GRBVar* gvars = grb_mod->getVars();
    for(int i = 0; i < grb_mod->get(GRB_IntAttr_NumVars); ++i) {
        if (gvars[i].get(GRB_CharAttr_VType) == 'B') gvars[i].set(GRB_CharAttr_VType,'C');
    }
}

void GurobiProgram::fill_in_grb_vmap(){
    param_* v;
    auto size_init=_grb_vars.size();
    _grb_vars.resize(_model->get_nb_vars());
    bool add;
    for(auto& v_p: _model->_vars)
    {
        v = v_p.second.get();
        if (!v->_new) {
            continue;
        }
        if(v->_new && size_init>v->get_id())
        {
            add=false;
        }
        else{
            add=true;
        }
        //DebugOn("to add v"<<endl);
        v->_new = false;
        auto idx = v->get_id();
        switch (v->get_intype()) {
            case float_: {
                auto real_var = (var<float>*)v;
                for (int i = 0; i < real_var->_dim[0]; i++) {
                    auto vid = idx + i;
                    if(add){
                    _grb_vars.at(vid) = (GRBVar(grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_CONTINUOUS, v->get_name(true,true)+"("+v->_indices->_keys->at(i)+")")));
                }
                else{
                    _grb_vars.at(vid).set(GRB_DoubleAttr_LB, real_var->get_lb(i));
                    _grb_vars.at(vid).set(GRB_DoubleAttr_UB, real_var->get_ub(i));
                }
                }
                break;
            }
            case long_:{
                auto real_var = (var<long double>*)v;
                for (int i = 0; i < real_var->_dim[0]; i++) {
                    auto vid = idx + i;
                    if(add){
                    _grb_vars.at(vid) = (GRBVar(grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_CONTINUOUS, v->get_name(true,true)+"("+v->_indices->_keys->at(i)+")")));
                }
                else{
                    _grb_vars.at(vid).set(GRB_DoubleAttr_LB, real_var->get_lb(i));
                    _grb_vars.at(vid).set(GRB_DoubleAttr_UB, real_var->get_ub(i));
                }
                }
                break;
            }
            case double_:{
                auto real_var = (var<double>*)v;
                for (size_t i = 0; i < real_var->_dim[0]; i++) {
                    auto vid = idx + i;
                    if(add){
                    _grb_vars.at(vid) = (GRBVar(grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_CONTINUOUS, v->get_name(true,true)+"("+v->_indices->_keys->at(i)+")")));
                    //DebugOn("added var"<<endl);
                }
                    else{
                    _grb_vars.at(vid).set(GRB_DoubleAttr_LB, real_var->get_lb(i));
                    _grb_vars.at(vid).set(GRB_DoubleAttr_UB, real_var->get_ub(i));
                        DebugOff("updated bounds"<<endl);
                }
                }
                break;
            }
            case integer_:{
                auto real_var = (var<int>*)v;
                for (int i = 0; i < real_var->_dim[0]; i++) {
                    auto vid = idx + i;
                    if(add){
                    _grb_vars.at(vid) = (GRBVar(grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_INTEGER, v->get_name(true,true)+"("+v->_indices->_keys->at(i)+")")));
                }
                else{
                    _grb_vars.at(vid).set(GRB_DoubleAttr_LB, real_var->get_lb(i));
                    _grb_vars.at(vid).set(GRB_DoubleAttr_UB, real_var->get_ub(i));
                }

            }
                break;
            }
            case short_:{
                auto real_var = (var<short>*)v;
                for (int i = 0; i < real_var->_dim[0]; i++) {
                    auto vid = idx + i;
                    if(add){
                    _grb_vars.at(vid) = (GRBVar(grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_INTEGER, v->get_name(true,true)+"("+v->_indices->_keys->at(i)+")")));
                }
                else{
                    _grb_vars.at(vid).set(GRB_DoubleAttr_LB, real_var->get_lb(i));
                    _grb_vars.at(vid).set(GRB_DoubleAttr_UB, real_var->get_ub(i));
                }
                }
                break;
            }
            case binary_:{
                auto real_var = (var<bool>*)v;
                for (int i = 0; i < real_var->_dim[0]; i++) {
                    auto vid = idx + i;
                    if(add){
                    _grb_vars.at(vid) = (GRBVar(grb_mod->addVar(real_var->get_lb(i), real_var->get_ub(i), 0.0, GRB_BINARY, v->get_name(true,true)+"("+v->_indices->_keys->at(i)+")")));
                   }
                   else{
                    _grb_vars.at(vid).set(GRB_DoubleAttr_LB, real_var->get_lb(i));
                    _grb_vars.at(vid).set(GRB_DoubleAttr_UB, real_var->get_ub(i));
                }
                }
                break;
            }
            default:
                break;
        }
        
    }
//    for(auto& v_p: _model->_vars)
//    {
//        v = v_p.second;
//        auto real_var = (var<double>*)v;
//        for (int i = 0; i < real_var->_dim[0]; i++) {
//            auto vid = v->get_id() + v->get_id_inst(i);
//            DebugOn("VID = "<< vid <<" : " << _grb_vars.at(vid).get(GRB_StringAttr_VarName) << " in [" << _grb_vars.at(vid).get(GRB_DoubleAttr_LB) << "," << _grb_vars.at(vid).get(GRB_DoubleAttr_UB) << "]\n" );
//        }
//    }
}

void GurobiProgram::create_grb_constraints(bool use_gravity_name){
    char sense;
    size_t idx = 0, idx_inst = 0, idx1 = 0, idx2 = 0, idx_inst1 = 0, idx_inst2 = 0, nb_inst = 0, inst = 0, c_idx = 0;
    GRBLinExpr lterm, linlhs;
    GRBQuadExpr quadlhs;
    GRBVar gvar1, gvar2;
    double coeff;
    for(auto& p: _model->_cons){
        auto c = p.second;
        if (!c->_new && c->_all_satisfied) {
            continue;
        }
        c->_new = false;
        //DebugOn("to add constraint"<<endl);
        if (c->is_nonlinear()) {
            DebugOn("nonlinear"<<endl);
            throw invalid_argument("Gurobi cannot handle nonlinear constraints that are not convex quadratic.\n");
        }
        switch(c->get_ctype()) {
            case geq:
                sense = GRB_GREATER_EQUAL;
                break;
            case leq:
                sense = GRB_LESS_EQUAL;
                break;
            case eq:
                sense = GRB_EQUAL;
                break;
            default:
                break;
        }
        nb_inst = c->get_nb_inst();
        inst = 0;
        if (c->is_linear()) {
            for (size_t i = 0; i< nb_inst; i++){
                if (c->_violated[i]) {
                    linlhs = 0;
                    for (auto& it1: c->get_lterms()) {
                        lterm = 0;
                        if (it1.second._p->_is_vector || it1.second._p->is_matrix_indexed() || it1.second._coef->is_matrix()) {
                            auto dim =it1.second._p->get_dim(i);
                            for (int j = 0; j<dim; j++) {
                                coeff = c->eval(it1.second._coef,i,j);
                                    gvar1 = _grb_vars[it1.second._p->get_id() + it1.second._p->get_id_inst(i,j)];
                                    lterm += coeff*gvar1;
                            }
                        }
                        else {
                            coeff = c->eval(it1.second._coef,i);
                                gvar1 = _grb_vars[it1.second._p->get_id() + it1.second._p->get_id_inst(i)];
                                lterm += coeff*gvar1;
                        }
                        if (!it1.second._sign) {
                            lterm *= -1;
                        }
                        linlhs += lterm;
                    }
                    linlhs += c->eval(c->get_cst(), i);
                    if(use_gravity_name){
                        if(c->_indices)
                            grb_mod->addConstr(linlhs,sense,0,c->get_name()+"("+c->_indices->_keys->at(i)+")");
                        else
                            grb_mod->addConstr(linlhs,sense,0,c->get_name());
                    }
                    else {
                        grb_mod->addConstr(linlhs,sense,0,"C("+to_string(c_idx++)+")");
                    }
                }
            }
        }
        else {
            for (size_t i = 0; i< nb_inst; i++){
                if (c->_violated[i]) {
                    quadlhs = 0;
                    for (auto& it1: c->get_lterms()) {
                        lterm = 0;
                        if (it1.second._coef->_is_transposed || it1.second._coef->is_matrix() || it1.second._p->is_matrix_indexed()) {
                            auto dim = it1.second._p->get_dim(i);
                            for (size_t j = 0; j<dim; j++) {
                                coeff = c->eval(it1.second._coef,i,j);
                                gvar1 = _grb_vars[it1.second._p->get_id() + it1.second._p->get_id_inst(i,j)];
                                lterm += coeff*gvar1;
                            }
                        }
                        else {
                            coeff = c->eval(it1.second._coef,i);
                            gvar1 = _grb_vars[it1.second._p->get_id() + it1.second._p->get_id_inst(i)];
                            lterm += coeff*gvar1;
                        }
                        if (!it1.second._sign) {
                            lterm *= -1;
                        }
                        quadlhs += lterm;
                    }
                    for (auto& it1: c->get_qterms()) {
                        if (it1.second._coef_p1_tr) { // qterm = (coef*p1)^T*p2
                            for (auto i = 0; i<it1.second._p->first->get_dim(); i++) {
                                for (auto j = 0; j<it1.second._p->first->get_dim(); j++) {
                                    coeff = _model->_obj->eval(it1.second._coef,i,j);
                                    gvar1 = _grb_vars[it1.second._p->first->get_id() + it1.second._p->first->get_id_inst(j)];
                                    gvar2 = _grb_vars[it1.second._p->second->get_id() + it1.second._p->second->get_id_inst(i)];
                                    if (!it1.second._sign) {
                                        quadlhs -= coeff*gvar1*gvar2;
                                    }
                                    else {
                                        quadlhs += coeff*gvar1*gvar2;
                                    }
                                }
                            }
                        }
                        else if (it1.second._coef->_is_transposed || it1.second._coef->is_matrix_indexed() || it1.second._p->first->is_matrix_indexed()) {
                            auto dim =it1.second._p->first->get_dim(i);
                            for (int j = 0; j<dim; j++) {
                                coeff = c->eval(it1.second._coef,i,j);
                                gvar1 = _grb_vars[it1.second._p->first->get_id() + it1.second._p->first->get_id_inst(i,j)];
                                gvar2 = _grb_vars[it1.second._p->second->get_id() + it1.second._p->second->get_id_inst(i,j)];
                                if (!it1.second._sign) {
                                    quadlhs += -1*coeff*gvar1*gvar2;
                                }
                                else {
                                    quadlhs += coeff*gvar1*gvar2;
                                }
                            }
                        }
                        else {
                            gvar1 = _grb_vars[it1.second._p->first->get_id() + it1.second._p->first->get_id_inst(i)];
                            gvar2 = _grb_vars[it1.second._p->second->get_id() + it1.second._p->second->get_id_inst(i)];
                            coeff = c->eval(it1.second._coef,i);
                            if (!it1.second._sign) {
                                quadlhs += -1*coeff*gvar1*gvar2;
                            }
                            else {
                                quadlhs += coeff*gvar1*gvar2;
                            }
                        }
                    }
                    quadlhs += c->eval(c->get_cst(), i);
                    if(use_gravity_name){
                        if(c->_indices)
                            grb_mod->addQConstr(quadlhs,sense,0,c->get_name()+"("+c->_indices->_keys->at(i)+")");
                        else
                            grb_mod->addQConstr(quadlhs,sense,0,c->get_name());
                    }
                    else{
                        grb_mod->addQConstr(quadlhs,sense,0,"C("+to_string(c_idx++)+")");
                    }
                }
            }
        }
    }
}

void GurobiProgram::set_grb_objective(){
//    size_t idx = 0;
    GRBLinExpr lterm;
    GRBQuadExpr qobj;
    GRBVar gvar1, gvar2;
    int objt;
    double coeff;
    if (!_model->_obj->_new) {
        return;
    }
    _model->_obj->_new = false;
    if (_model->_objt == minimize) objt = GRB_MINIMIZE;
    else objt = GRB_MAXIMIZE;
        qobj = 0;
        for (auto& it1: _model->_obj->get_lterms()) {
            lterm = 0;
            if (it1.second._coef->_is_transposed || it1.second._coef->is_matrix() || it1.second._p->is_matrix_indexed()) {
                auto dim = it1.second._p->get_dim(0);
                for (size_t j = 0; j<dim; j++) {
                    coeff = _model->_obj->eval(it1.second._coef,0,j);
                    gvar1 = _grb_vars[it1.second._p->get_id() + it1.second._p->get_id_inst(0,j)];
                    lterm += coeff*gvar1;
                }
            }
            else {
                coeff = _model->_obj->eval(it1.second._coef);
                gvar1 = _grb_vars[it1.second._p->get_id() + it1.second._p->get_id_inst()];
                lterm += coeff*gvar1;
            }
            if (!it1.second._sign) {
                lterm *= -1;
            }
            qobj += lterm;
        }
    for (auto& it1: _model->_obj->get_qterms()) {
        if (it1.second._coef_p1_tr) { // qterm = (coef*p1)^T*p2
//            assert(it_qterm.second._p->first->_dim[1]==1 && it_qterm.second._coef->_dim[0]==it_qterm.second._p->second->_dim[0]);
            for (auto i = 0; i<it1.second._p->first->get_dim(); i++) {
                for (auto j = 0; j<it1.second._p->first->get_dim(); j++) {
                    coeff = _model->_obj->eval(it1.second._coef,i,j);
                    gvar1 = _grb_vars[it1.second._p->first->get_id() + it1.second._p->first->get_id_inst(j)];
                    gvar2 = _grb_vars[it1.second._p->second->get_id() + it1.second._p->second->get_id_inst(i)];
                    if (!it1.second._sign) {
                        qobj -= coeff*gvar1*gvar2;
                    }
                    else {
                        qobj += coeff*gvar1*gvar2;
                    }
                }
            }
        }
        else if (it1.second._coef->_is_transposed) {
            auto dim =it1.second._p->first->get_dim();
            for (int j = 0; j<dim; j++) {
                coeff = _model->_obj->eval(it1.second._coef,j);
                gvar1 = _grb_vars[it1.second._p->first->get_id() + it1.second._p->first->get_id_inst(j)];
                gvar2 = _grb_vars[it1.second._p->second->get_id() + it1.second._p->second->get_id_inst(j)];
                if (!it1.second._sign) {
                    qobj += -1*coeff*gvar1*gvar2;
                }
                else {
                    qobj += coeff*gvar1*gvar2;
                }
            }
        }
        else {
            coeff = _model->_obj->eval(it1.second._coef);
            if (!it1.second._sign) {
                qobj += -1*coeff*gvar1*gvar2;
            }
            else {
                qobj += coeff*gvar1*gvar2;
            }
        }
    }
    qobj += _model->_obj->eval(_model->_obj->get_cst());
    grb_mod->setObjective(qobj,objt);
//    grb_mod->update();
}

void GurobiProgram::print_constraints(){
    GRBConstr* gconstrs = grb_mod->getConstrs();
    for(int i = 0; i < grb_mod->get(GRB_IntAttr_NumConstrs); ++i) {
        if(gconstrs[i].get(GRB_CharAttr_Sense)!='=') {
            cout << gconstrs[i].get(GRB_StringAttr_ConstrName) << "  ";
//            cout << gconstrs[i].get(GRB_DoubleAttr_Slack);
            cout << "\n";
        }
    }
    GRBQConstr* gqconstrs = grb_mod->getQConstrs();
    for(int i = 0; i < grb_mod->get(GRB_IntAttr_NumQConstrs); ++i) {
        if(gqconstrs[i].get(GRB_CharAttr_QCSense)!='=') {
            cout << gqconstrs[i].get(GRB_StringAttr_QCName) << "  ";
//            cout << gqconstrs[i].get(GRB_DoubleAttr_Slack);
//            cout << gqconstrs[i].get(GRB_DoubleAttr_QCSlack);
            cout << "\n";
        }
    }
}

