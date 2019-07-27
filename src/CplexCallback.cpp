////
////  CplexCallback.cpp
////  Gravity
////
////  Created by Mertcan Yetkin on 20/06/2019.
////
////
//
//#include <ilcplex/ilocplex.h>
//#include <gravity/CplexCallback.h>
//
//using namespace gravity;
//
////CplexCallback::CplexCallback(const IloNumVarArray &_callback_vars):callback_vars(_callback_vars){}
//
//void CplexCallback::dynamicCuts(const IloCplex::Callback::Context &context) const{}
//
//void CplexCallback::lazyCuts (const IloCplex::Callback::Context &context) const{}
//
//void CplexCallback::invoke (const IloCplex::Callback::Context &context){}
//
//
///* Implementation of the invoke function */
//void invoke (const IloCplex::Callback::Context &context){};
////{
////    if ( context.inRelaxation() ) {
////        if ( cuts.getSize() > 0 ) {
////            cutsFromTable(context);
////        }
////        else {
////            separateDisagregatedCuts(context);
////        }
////    }
////
////    if ( context.inCandidate() )
////        lazyCapacity (context);
//
////}
//
//CplexCallback::~CplexCallback(){}
