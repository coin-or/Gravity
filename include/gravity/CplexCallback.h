//
//  CplexProgram.h
//  Gravity
//
//  Created by Mertcan Yetkin on 20/06/2019.
//
//

#ifndef CplexCallback_h
#define CplexCallback_h

#include <stdio.h>

#ifndef __Gravity____CplexCallback
#define __Gravity____CplexCallback

#ifdef USE_CPLEX
#include <ilcplex/ilocplex.h>
#endif
#include <gravity/model.h>
using namespace gravity;   


class CplexCallback: public IloCplex::Callback::Function {
private:
    /* Variables to be considered in the callback. */
    IloNumVarArray callback_vars;
    
public:
    
    CplexCallback(){};
    
    /* Constructor with data */
    CplexCallback(const IloNumVarArray &_callback_vars):callback_vars(_callback_vars){}
    
    
    
    void dynamicCuts (const IloCplex::Callback::Context &context) const;
    //        IloInt const nbLocations = opened.getSize();
    //        IloInt const nbClients = supply.getSize();
    //
    //        // For each j and c check whether in the current solution (obtained by
    //        // calls to getValue()) we have supply[c][j] > opened[j]. If so, then we have
    //        // found a violated constraint and add it as a cut.
    //        for (IloInt j = 0; j < nbLocations; ++j) {
    //            for (IloInt c = 0; c < nbClients; ++c) {
    //                IloNum const s = context.getRelaxationPoint(supply[c][j]);
    //                IloNum const o = context.getRelaxationPoint(opened[j]);
    //                if ( s > o + EPS) {
    //                    cout << "Adding: " << supply[c][j].getName() << " <= "
    //                    << opened[j].getName() << " [" << s << " > " << o << "]" << endl;
    //                    context.addUserCut( supply[c][j] - opened[j] <= 0,
    //                                       IloCplex::UseCutPurge, IloFalse);
    //                }
    //            }
    //        }
    //    }
    
    void lazyCuts(const IloCplex::Callback::Context &context) const;
    //        IloInt const nbLocations = opened.getSize();
    //        IloInt const nbClients = supply.getSize();
    //        if ( !context.isCandidatePoint() )
    //            throw IloCplex::Exception(-1, "Unbounded solution");
    //        for (IloInt j = 0; j < nbLocations; ++j) {
    //            IloNum isUsed = context.getCandidatePoint(opened[j]);
    //            IloNum served = 0.0; // Number of clients currently served from j
    //            for (IloInt c = 0; c < nbClients; ++c)
    //                served += context.getCandidatePoint(supply[c][j]);
    //            if ( served > (nbClients - 1.0) * isUsed + EPS ) {
    //                IloNumExpr sum = IloExpr(context.getEnv());
    //                for (IloInt c = 0; c < nbClients; ++c)
    //                    sum += supply[c][j];
    //                sum -= (nbClients - 1) * opened[j];
    //                cout << "Adding lazy capacity constraint " << sum << " <= 0" << endl;
    //                context.rejectCandidate(sum <= 0.0);
    //                sum.end();
    //            }
    //        }
    //    }
    
    // This is the function that we have to implement and that CPLEX will call
    // during the solution process at the places that we asked for.
    virtual void invoke (const IloCplex::Callback::Context &context);
    
    /// Destructor
    virtual ~CplexCallback();
    
};

#endif /* defined(__Gravity____CplexCallback) */
#endif /* CplexCallback_hpp */
