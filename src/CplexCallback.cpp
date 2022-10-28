//
//  CplexCallback.cpp
//  Gravity
//
//  Created by Mertcan Yetkin on 20/06/2019.
//
//

#include <ilcplex/ilocplex.h>
#include <gravity/CplexCallback.h>

using namespace gravity;

//CplexCallback::CplexCallback(const IloNumVarArray &_callback_vars):callback_vars(_callback_vars){}

//void CplexCallback::dynamicCuts(const IloCplex::Callback::Context &context) const{
//    n = _cplex_vars.size();
//    int n = context.getEnv()
//    double *x = new double[n];
//    vector<double> vec_x;
//    // double obj=getDoubleInfo(GRB_CB_MIPSOL_OBJ);
//    int i,j;
//    x=getNodeRel(vars.data(),n);
//    //x=getSolution(vars.data(),n);
//    for(i=0;i<n;i++){
//        vec_x.push_back(x[i]);
//    }
//    m->set_solution(vec_x);
//
//    auto res= m->cutting_planes_soc(1e-9, soc_found, soc_added);
//    if(res.size()>=1){
//        for(auto i=0;i<res.size();i++){
//            GRBLinExpr expr = 0;
//            int j;
//            DebugOff("soc cut at ");
//            for(j=0;j<res[i].size()-1;j+=2){
//                int c=res[i][j];
//                if(std::abs(c)>=0)
//                    expr += res[i][j+1]*vars[c];
//                DebugOff(to_string_with_precision(vec_x[c],10)<<" "<<to_string_with_precision(c,10)<<" "<<to_string_with_precision(res[i][j+1],10)<<" ");
//            }
//            DebugOff(endl);
//            addLazy(expr, GRB_LESS_EQUAL, 0);
////                                addCut(expr, GRB_LESS_EQUAL, 0);
//
//            //vec_expi.push_back(expr);
//        }
//    }
//    if(res.size()==0 || !hierarc){
//        m->set_solution(vec_x);
//        auto res1=m->cuts_eigen_bags(1e-9);
//        if(res1.size()>=1){
//            for(auto i=0;i<res1.size();i++){
//                GRBLinExpr expr = 0;
//                int j;
//                DebugOff("eig cut at");
//                for(j=0;j<res1[i].size()-1;j+=2){
//                    int c=res1[i][j];
//                    expr += res1[i][j+1]*vars[c];
//                    DebugOff(to_string_with_precision(vec_x[c],10)<<" "<<to_string_with_precision(c,10)<<" "<<to_string_with_precision(res1[i][j+1],10)<<" ");
//                }
//                DebugOff(to_string_with_precision(res1[i][j],10));
//                DebugOff(endl);
//                if(std::abs(res1[i][j])>=1e-12)
//                    expr += res1[i][j];
//                addLazy(expr, GRB_LESS_EQUAL, 0);
////                                    addCut(expr, GRB_LESS_EQUAL, 0);
//                //vec_expi.push_back(expr);
//            }
//        }
//
//
//        if(res1.size()==0 || !hierarc){
//            m->set_solution(vec_x);
//            auto res2=m->cuts_eigen_full(1e-9);
//            if(res2.size()>=1){
//                for(auto i=0;i<res2.size();i++){
//                    GRBLinExpr expr = 0;
//                    int j;
//                    for(j=0;j<res2[i].size()-1;j+=2){
//                        int c=res2[i][j];
//                        expr += res2[i][j+1]*vars[c];
//                    }
//                    if(std::abs(res2[i][j])>=1e-12)
//                        expr += res2[i][j];
//                    addLazy(expr, GRB_LESS_EQUAL, 0);
//                    //                                    addCut(expr, GRB_LESS_EQUAL, 0);
//                    //vec_expi.push_back(expr);
//                }
//            }
//        }
//    }
//}

//void CplexCallback::lazyCuts (const IloCplex::Callback::Context &context) const{}



CplexCallback::~CplexCallback(){
    IloInt numWorkers = workers.size();
    for (IloInt w = 0; w < numWorkers; w++) {
       delete workers[w];
    }
    workers.clear();
    /* Delete all the model copies */
    for(auto i = 1; i<_models.size(); i++){
        delete _models[i];
    }
    _models.clear();
 }
