////
////
////  Gravity
////
////  Created by Hassan Hijazi on August 29 2018.
////
////
//
//
//#include <iostream>
//#include <string>
//#include <fstream>
//#include <gravity/solver.h>
//#include <stdio.h>
//#include <stdlib.h>
//#include <optionParser.hpp>
//
//
//using namespace std;
//using namespace gravity;
//
//
//void run_qtdesign(const int n, const int d, const param<>& M, bool relax){
//#ifdef USE_QPP
//    double tol = 1e-6;
//    bool ibm = false;
//    bool cont = true;
//    bool add_rx = true;
//    bool rev_cnot = false;
//    int output = 0;
//    if (ibm && cont) {
//        cerr << "Choose IBM or Continuous, or none.";
//        exit(-1);
//    }
//
//    /** Parameters **/
//
//    /* Qubit number */
//    const int m = pow(2,n);
//    DebugOn("Number of Qubits = " << n << endl);
//    /* Circuit depth */
//    DebugOn("Depth = " << d << endl);
//    if (d<2) {
//        cerr << "Depth should be greater or equal than 2!" << endl;
//        return;
//    }
//    /* T and T conjugate transpose gate matrices, one for each row (qubit) */
//    param<> T[n], Tt[n];
//
//    /* S and St conjugate transpose gate matrices, one for each row (qubit) */
//    param<> S[n], St[n];
//
//    /* Hadamard gate matrices */
//    param<> H[n];
//
//    /* Rx gate matrices */
//    param<> Rx[n];
//
//    /* Cnot gate matrices */
//    param<> Cnot[n*(n-1)/2];
//    param<> CnotR[n*(n-1)/2];
//    /* SWAP gate matrices */
//    param<> Swap[n*(n-1)/2];
//
//
//
//
//    /** Initialization **/
//    int idx = 0;
//    for (auto i = 0; i<n; i++) { // Qubit (row) iterator
//        T[i].set_name("T_"+to_string(i+1));
//        T[i].set_size(2*m, 2*m);/* Representing Complex number as 2x2 real matrices */
//        T[i].QuantumT(i,n);/* Fill nonzero values corresponding to a T gate */
//        Tt[i].set_name("Tt_"+to_string(i+1));
//        Tt[i].set_size(2*m, 2*m);
//        Tt[i].QuantumT(i,n,true);/* Fill nonzero values corresponding to a T conjugate transpose gate */
//        Rx[i].set_name("Rx_"+to_string(i+1));
//        Rx[i].set_size(2*m, 2*m);
//        Rx[i].QuantumRx(i,n);/* Fill nonzero values corresponding to an H gate */
//        S[i].set_name("S_"+to_string(i+1));
//        S[i].set_size(2*m, 2*m);
//        S[i].QuantumS(i,n);/* Fill nonzero values corresponding to an H gate */
//        St[i].set_name("St_"+to_string(i+1));
//        St[i].set_size(2*m, 2*m);
//        St[i].QuantumSt(i,n);/* Fill nonzero values corresponding to an H gate */
//        H[i].set_name("H_"+to_string(i+1));
//        H[i].set_size(2*m, 2*m);
//        H[i].QuantumH(i,n);/* Fill nonzero values corresponding to an H gate */
//        for (auto k = i+1; k<n; k++) {
//            Cnot[idx].set_name("Cnot_"+to_string(i+1)+to_string(k+1));
//            Cnot[idx].set_size(2*m, 2*m);
//            Cnot[idx].QuantumCnot(i,k,n);/* Fill nonzero values corresponding to a CNOT gate */
//            CnotR[idx].set_name("CnotR_"+to_string(k+1)+to_string(i+1));
//            CnotR[idx].set_size(2*m, 2*m);
//            CnotR[idx].QuantumCnot(k,i,n);/* Fill nonzero values corresponding to a CNOT gate */
//            Swap[idx].set_name("Swap"+to_string(i+1)+to_string(k+1));
//            Swap[idx].set_size(2*m, 2*m);
//            Swap[idx].QuantumSwap(i,k,n);/* Fill nonzero values corresponding to a CNOT gate */
//
//            idx++;
//        }
//    }
//
//
//    /** Indices **/
//    auto d_ids = indices(1,d);/*< Depth indices */
//    auto q_ids = indices(1,n);/*< Qubit indices */
//    auto pairs = indices(ordered_pairs(2, n+1)); /*< Indices pairs (control-taget) for Cnot gates */
//    auto rpairs = indices(ordered_pairs(2, n+1, true)); /*< Indices pairs (control-taget) for Cnot gates */
//    auto m2 = indices(indices(1,2*m), indices(1,2*m)); /*< Indices for square mxm matrices */
//
//
//    /** Variables **/
//    /** Binaries **/
////    var<bool> zT("zT"); /*< Binary variables for T gates */
////    var<bool> zTt("zTt"); /*< Binary variables for T conjugate gates */
////    var<bool> zH("zH"); /*< Binary variables for H gates */
////    var<bool> zS("zS"); /*< Binary variables for S gates */
////    var<bool> zSt("zSt"); /*< Binary variables for S gates */
////    var<bool> zRx("zRx"); /*< Binary variables for Rx gates */
////    var<bool> zCnot("zCnot"); /*< Binary variables for CNOT gates */
////    var<bool> zCnotR("zCnotR"); /*< Binary variables for CNOT gates */
////    var<bool> zSwap("zSwap"); /*< Binary variables for SWAP gates */
//
//    /** Gate Matrices **/
//    var<> G[d]; /*< Base Gate Matrices */
//    param<> Gl[d], Gu[d]; /*< Lower and Upper Bounds on Base Gate Matrices */
//    bool zeros[2*m][2*m]; /*< Zeros in Gate Matrices */
//    var<> L[d-2]; /*< Lifted Product Matrices */
//    param<> Ll[d], Lu[d]; /*< Lower and Upper Bounds on Lifted Matrices */
////    var<> lambda_r("lambda_r", -1, 1);
////    var<> lambda_i("lambda_i", -1, 1);
////    var<> slack("slack",-1,1);
//
//
//
//    /** Model **/
//
//    /** Adding Variables **/
//    Model Qdesign("Qdesign");
//    /** Binaries **/
////    if(ibm){
////        if(add_rx) {
////                Qdesign.add(zS.in(d_ids,q_ids), zSt.in(d_ids,q_ids), zT.in(d_ids,q_ids),
////                        zTt.in(d_ids,q_ids),zRx.in(d_ids,q_ids), zCnot.in(d_ids,pairs));
////        }
////        else {
////            Qdesign.add(zS.in(d_ids,q_ids), zSt.in(d_ids,q_ids), zT.in(d_ids,q_ids), zTt.in(d_ids,q_ids),zCnot.in(d_ids,pairs));
////        }
////    }
////    else if(cont){
////        Qdesign.add(lambda_r.in(d_ids), lambda_i.in(d_ids));
////        if(add_rx){
////            if (rev_cnot) {
////                Qdesign.add(zS.in(d_ids,q_ids), zRx.in(d_ids,q_ids) ,zCnot.in(d_ids,pairs), zCnotR.in(d_ids,rpairs));
////            }
////            else {
////                Qdesign.add(zS.in(d_ids,q_ids), zRx.in(d_ids,q_ids) ,zCnot.in(d_ids,pairs));
////            }
////        }
////        else {
////            if (rev_cnot) {
////                Qdesign.add(zS.in(d_ids,q_ids), zCnot.in(d_ids,pairs), zCnotR.in(d_ids,rpairs));
////            }
////            else {
////                Qdesign.add(zS.in(d_ids,q_ids), zCnot.in(d_ids,pairs));
////            }
////        }
////    }
////    else {
////        Qdesign.add(zT.in(d_ids,q_ids), zTt.in(d_ids,q_ids), zH.in(d_ids,q_ids), zCnot.in(d_ids,pairs), zCnotR.in(d_ids,rpairs));
////    }
////    Qdesign.add(slack.in(R(2*m,2*m)));
//    DebugOff("Total number of variables after adding the binaries = " << Qdesign.get_nb_vars() << endl);
//    /** Gate Matrices **/
//    for (int depth = 0; depth < d; depth++) {
//        Gl[depth] = param<>("Gl_"+to_string(depth+1));
//        Gu[depth] = param<>("Gu_"+to_string(depth+1));
//        if (relax) {
//            G[depth] = var<>("G_"+to_string(depth+1), Gl[depth].in(m2), Gu[depth].in(m2));
//        }
//        else{
//            G[depth] = var<>("G_"+to_string(depth+1), -1, 1);
//        }
//
//        Qdesign.add(G[depth].in(m2));
//    }
//    if (d>2) {
//        for (int depth = 0; depth < d-2; depth++) {
//            Ll[depth] = param<>("Ll_"+to_string(depth+1));
//            Lu[depth] = param<>("Lu_"+to_string(depth+1));
//            if (relax) {
//                L[depth] = var<>("L_"+to_string(depth+1), Ll[depth].in(m2), Lu[depth].in(m2));
//            }
//            else {
//                L[depth] = var<>("L_"+to_string(depth+1), -1, 1);
//            }
//            Qdesign.add(L[depth].in(R(2*m,2*m)));
//        }
//    }
//
//    DebugOff("Total number of variables after adding the G vars= " << Qdesign.get_nb_vars() << endl);
//
//
//    /** Adding Constraints **/
//    /** Building indexed variables as matrices (for a faster access in constraints) **/
//    var<bool> zT_[d][n], zTt_[d][n], zH_[d][n], zS_[d][n], zSt_[d][n], zRx_[d][n], zCnot_[d][n][n], zCnotR_[d][n][n], zSwap_[d][n][n];
//    var<> G_[d][2*m][2*m], L_[d-2][2*m][2*m];
//    param<> Gl_[d][2*m][2*m], Gu_[d][2*m][2*m], Ll_[d-2][2*m][2*m], Lu_[d-2][2*m][2*m];
//    var<> slack_[2*m][2*m];
//    var<> lambda_r_[d], lambda_i_[d];
//    for (int i = 0; i < 2*m; i++) {
//        for (int j = 0; j < 2*m; j++) {
//            slack_[i][j] = slack(i+1,j+1);
//        }
//    }
//    for (size_t depth = 0; depth < d; depth++) {
//        lambda_r_[depth] = lambda_r(to_string(depth+1));
//        lambda_i_[depth] = lambda_i(to_string(depth+1));
//        for (size_t qubit1 = 0; qubit1 < n; qubit1++) {
//            zT_[depth][qubit1] = zT(depth+1,qubit1+1);
//            zTt_[depth][qubit1] = zTt(depth+1,qubit1+1);
//            zH_[depth][qubit1] = zH(depth+1,qubit1+1);
//            zS_[depth][qubit1] = zS(depth+1,qubit1+1);
//            zSt_[depth][qubit1] = zSt(depth+1,qubit1+1);
//            zRx_[depth][qubit1] = zRx(depth+1,qubit1+1);
//            for (size_t qubit2 = qubit1+1; qubit2 < n; qubit2++) {
//                zCnot_[depth][qubit1][qubit2] = zCnot(to_string(depth+1)+","+to_string(qubit1+1)+","+to_string(qubit2+1));
//                zCnotR_[depth][qubit1][qubit2] = zCnotR(to_string(depth+1)+","+to_string(qubit2+1)+","+to_string(qubit1+1));
//                zSwap_[depth][qubit1][qubit2] = zSwap(to_string(depth+1)+","+to_string(qubit1+1)+","+to_string(qubit2+1));
//
//            }
//        }
//        for (int i = 0; i < 2*m; i++) {
//            for (int j = 0; j < 2*m; j++) {
//                G_[depth][i][j] = G[depth](i+1,j+1);
//                Gl_[depth][i][j] = Gl[depth](i+1,j+1);
//                Gu_[depth][i][j] = Gu[depth](i+1,j+1);
//            }
//        }
//    }
//    if (d>2) {
//        for (size_t depth = 0; depth < d-2; depth++) {
//            for (int i = 0; i < 2*m; i++) {
//                for (int j = 0; j < 2*m; j++) {
//                    L_[depth][i][j] = L[depth](i+1,j+1);
//                    Ll_[depth][i][j] = Ll[depth](i+1,j+1);
//                    Lu_[depth][i][j] = Lu[depth](i+1,j+1);
//                }
//            }
//        }
//    }
//
//    /** Objective **/
////    Qdesign.min(sum(zRx) + sum(zH) + sum(zT) + sum(zTt) + sum(zCnot));
////    Qdesign.min(sum(zRx) + sum(zS) + sum(zCnot));
//
//    func_ obj;
//    if(!relax){
//        for (int i = 0; i < 2*m; i++) {
//            for (int j = 0; j < 2*m; j++) {
//                obj += slack_[i][j]*slack_[i][j];
//            }
//        }
//    }
//    if(relax){
//        obj += sum(zCnot);
//    }
//    Qdesign.min(obj);
//
//
//    /** Base gate matrix definition **/
//    for (int depth = 0; depth < d; depth++) {
//        if(cont){
//            for (int qubit = 0; qubit < n; qubit++) {
//                Constraint Complex("Complex"+to_string(depth+1)+to_string(qubit+1));
//                Complex += lambda_i_[depth]*lambda_i_[depth] + lambda_r_[depth]*lambda_r_[depth];
//                if(relax){
//                    Qdesign.add(Complex==1);
//                }
//                else{
//                    Qdesign.add(Complex>=1);
//                }
//            }
//        }
//
//
//        for (int i = 0; i < 2*m; i++) {
//            for (int j = 0; j < 2*m; j++) {
//                double lb = 1, ub = -1;
//                Constraint BaseGate("BaseGate_"+to_string(depth+1)+"_"+to_string(i+1)+to_string(j+1));
//                if (i%2==0) {
//                    BaseGate = G_[depth][i][j];
//                }
//                else {
//                    if (j%2==0) {//Imaginary
//                        BaseGate = -1*G_[depth][i-1][j+1];
//                    }
//                    else {
//                        BaseGate = G_[depth][i-1][j-1];
//                    }
//                }
//                int idx = 0;
//                for (int qubit = 0; qubit < n; qubit++) {
//                    if(ibm){
//                        if(add_rx){
//                            BaseGate -= zS_[depth][qubit]*S[qubit].eval(i,j)+ zSt_[depth][qubit]*St[qubit].eval(i,j) + zT_[depth][qubit]*T[qubit].eval(i,j)+ zTt_[depth][qubit]*Tt[qubit].eval(i,j) + zRx_[depth][qubit]*Rx[qubit].eval(i,j);
//                        }
//                        else {
//                            BaseGate -= zS_[depth][qubit]*S[qubit].eval(i,j)+ zSt_[depth][qubit]*St[qubit].eval(i,j) + zT_[depth][qubit]*T[qubit].eval(i,j)+ zTt_[depth][qubit]*Tt[qubit].eval(i,j);
//                        }
//                    }
//                    else if (cont){
//                        if ((i%2!=0 && j%2==0) || (i%2==0 && j%2!=0)) {//imaginary part
//                            if(add_rx){
//                                BaseGate -= S[qubit].eval(i,j)*lambda_i_[depth]*zS_[depth][qubit] + zRx_[depth][qubit]*Rx[qubit].eval(i,j);
//                                if(i%2==0){
//                                    lb = min(lb,S[qubit].eval(i,j)*-1);
//                                    lb = min(lb,Rx[qubit].eval(i,j));
//                                    ub = max(ub,S[qubit].eval(i,j));
//                                    ub = max(ub,Rx[qubit].eval(i,j));
//                                }
//                            }
//                            else {
//                                BaseGate -= S[qubit].eval(i,j)*lambda_i_[depth]*zS_[depth][qubit];
//                                if(i%2==0){
//                                    lb = min(lb,S[qubit].eval(i,j)*-1);
//                                    ub = max(ub,S[qubit].eval(i,j));
//                                }
//                            }
//                        }
//                        else {
//                            if(add_rx){
//                                if (i%2==0) {
//                                    if(S[qubit].eval(i,j+1)!=0){
//                                        BaseGate -= lambda_r_[depth]*zS_[depth][qubit] + zRx_[depth][qubit]*Rx[qubit].eval(i,j);
//                                        lb = min(lb,-1.);
//                                        lb = min(lb,Rx[qubit].eval(i,j));
//                                        ub = max(ub,1.);
//                                        ub = max(ub,Rx[qubit].eval(i,j));
//                                    }
//                                    else {
//                                        BaseGate -= zS_[depth][qubit]*S[qubit].eval(i,j) + zRx_[depth][qubit]*Rx[qubit].eval(i,j);
//                                        lb = min(lb,S[qubit].eval(i,j));
//                                        lb = min(lb,Rx[qubit].eval(i,j));
//                                        ub = max(ub,S[qubit].eval(i,j));
//                                        ub = max(ub,Rx[qubit].eval(i,j));
//                                    }
//
//                                }
//                                else {
//                                    if(S[qubit].eval(i,j-1)!=0){
//                                        BaseGate -= lambda_r_[depth]*zS_[depth][qubit] + zRx_[depth][qubit]*Rx[qubit].eval(i,j);
//                                    }
//                                    else {
//                                        BaseGate -= zS_[depth][qubit]*S[qubit].eval(i,j) + zRx_[depth][qubit]*Rx[qubit].eval(i,j);
//                                    }
//                                }
//                            }
//                            else {
//                                if (i%2==0) {
//                                    if(S[qubit].eval(i,j+1)==-1){
//                                        BaseGate -= lambda_r_[depth]*zS_[depth][qubit];
//                                        lb = min(lb,-1.);
//                                        ub = max(ub,1.);
//                                    }
//                                    else {
//                                        BaseGate -= zS_[depth][qubit]*S[qubit].eval(i,j);
//                                        lb = min(lb,S[qubit].eval(i,j));
//                                        ub = max(ub,S[qubit].eval(i,j));
//                                    }
//                                }
//                                else {
//                                    if(S[qubit].eval(i,j-1)==1){
//                                        BaseGate -= lambda_r_[depth]*zS_[depth][qubit];
//                                    }
//                                    else {
//                                        BaseGate -= zS_[depth][qubit]*S[qubit].eval(i,j);
//                                    }
//                                }
//                            }
//                        }
//                    }
//                    else{
//                        BaseGate -= zT_[depth][qubit]*T[qubit].eval(i,j) + zTt_[depth][qubit]*Tt[qubit].eval(i,j) + zH_[depth][qubit]*H[qubit].eval(i,j);
//                        lb = min(lb,T[qubit].eval(i,j));
//                        lb = min(lb,Tt[qubit].eval(i,j));
//                        lb = min(lb,H[qubit].eval(i,j));
//                        ub = max(ub,T[qubit].eval(i,j));
//                        ub = max(ub,Tt[qubit].eval(i,j));
//                        ub = max(ub,H[qubit].eval(i,j));
//                    }
//
//                    for (int qubit2 = qubit+1; qubit2 < n; qubit2++) {
//                        BaseGate -= zCnot_[depth][qubit][qubit2]*Cnot[idx].eval(i,j);
//                        lb = min(lb, Cnot[idx].eval(i,j));
//                        ub = max(ub, Cnot[idx].eval(i,j));
//                        if (rev_cnot) {
//                            BaseGate -= zCnotR_[depth][qubit][qubit2]*CnotR[idx].eval(i,j);
//                            lb = min(lb, CnotR[idx].eval(i,j));
//                            ub = max(ub, CnotR[idx].eval(i,j));
//                        }
//                        idx++;
//                    }
//                }
//                if (BaseGate.get_nb_vars()==1) { /* This implies that the corresponding entry is zero */
//                    zeros[i][j]=true;
//                }
//                else {
//                    zeros[i][j]=false;
//                    if(i%2==0){
//                        Gl_[depth][i][j] = lb;
//                        Gu_[depth][i][j] = ub;
//                        if (j%2==0) {
//                            Gl_[depth][i+1][j+1] = lb;
//                            Gu_[depth][i+1][j+1] = ub;
//                        }
//                        else {
//                            Gl_[depth][i+1][j-1] = -ub;
//                            Gu_[depth][i+1][j-1] = -lb;
//                        }
//                        Qdesign.add(BaseGate==0);
//                    }
////                    BaseGate.print();
//                }
//            }
//        }
//    }
//
//    /** Lifted + Target gate constraint **/
//    /* No need to use lifted matrices if the depth = 2 */
//    func_ G1, G2;
//    if (d==2) {
//        for (int i = 0; i < 2*m; i+=2) {
//            for (int j = 0; j < 2*m; j++) {
//                Constraint TargetGate("TargetGate_"+to_string(i+1)+to_string(j+1));
//                for (int k = 0; k < 2*m; k++) {
//                    if (!zeros[i][k] && !zeros[k][j]) {
//                        if (k%2==0) {
//                            G2 = G_[1][k][j];
//                        }
//                        else {
//                            if (j%2==0) {//Imaginary
//                                G2 = -1*G_[1][k-1][j+1];
//                            }
//                            else {
//                                G2 = G_[1][k-1][j-1];
//                            }
//                        }
//                        TargetGate += G_[0][i][k]*G2;
//                    }
//                }
//                TargetGate -=  slack_[i][j];
//                Qdesign.add(TargetGate==M.eval(i,j));
//            }
//        }
//    }
//    else {
//        /** Lifted gate constraint **/
//        /** Depth 0 **/
//        for (int i = 0; i < 2*m; i+=2) {
//            for (int j = 0; j < 2*m; j++) {
//                double lb = 1, ub = -1;
//                Constraint LiftedGate("LiftedGate_0_"+to_string(i+1)+to_string(j+1));
//                LiftedGate += L_[0][i][j];
//                for (int k = 0; k < 2*m; k++) {
//                    if (!zeros[i][k] && !zeros[k][j]) {
//                        lb = min(lb, (Gl_[0][i][k].eval()*Gl_[1][k][j].eval()));
//                        lb = min(lb, (Gl_[0][i][k].eval()*Gu_[1][k][j].eval()));
//                        lb = min(lb, (Gu_[0][i][k].eval()*Gl_[1][k][j].eval()));
//                        lb = min(lb, (Gu_[0][i][k].eval()*Gu_[1][k][j].eval()));
//                        ub = max(ub, (Gl_[0][i][k].eval()*Gl_[1][k][j].eval()));
//                        ub = max(ub, (Gl_[0][i][k].eval()*Gu_[1][k][j].eval()));
//                        ub = max(ub, (Gu_[0][i][k].eval()*Gl_[1][k][j].eval()));
//                        ub = max(ub, (Gu_[0][i][k].eval()*Gu_[1][k][j].eval()));
//                        if (k%2==0) {
//                            G2 = G_[1][k][j];
//                        }
//                        else {
//                            if (j%2==0) {//Imaginary
//                                G2 = -1*G_[1][k-1][j+1];
//                            }
//                            else {
//                                G2 = G_[1][k-1][j-1];
//                            }
//                        }
//                        LiftedGate -= G_[0][i][k]*G2;
//                        Ll_[0][i][j] = Ll_[0][i][j].eval() + lb;
//                        Lu_[0][i][j] = Lu_[0][i][j].eval() + ub;
//                        if(Ll_[0][i][j].eval()<-1){
//                            Ll_[0][i][j] = -1;
//                        }
//                        if(Lu_[0][i][j].eval()>1){
//                            Lu_[0][i][j] = 1;
//                        }
//                        if (j%2==0) {//Real part
//                            Ll_[0][i+1][j+1] = Ll_[0][i][j].eval();
//                            Lu_[0][i+1][j+1] = Lu_[0][i][j].eval();
//                        }
//                        else{
//                            Ll_[0][i+1][j-1] = -Lu_[0][i][j].eval();
//                            Lu_[0][i+1][j-1] = -Ll_[0][i][j].eval();
//                        }
//                    }
//                }
//                Qdesign.add(LiftedGate==0);
//            }
//        }
//        for (int depth = 1; depth < d-2; depth++) {
//            for (int i = 0; i < 2*m; i+=2) {
//                for (int j = 0; j < 2*m; j++) {
//                    double lb = 1, ub = -1;
//                    Constraint LiftedGate("LiftedGate"+to_string(depth)+"_"+to_string(i+1)+to_string(j+1));
//                    LiftedGate += L_[depth][i][j];
//                    for (int k = 0; k < 2*m; k++) {
//                        if (!zeros[k][j]) {
//                            lb = min(lb, (Ll_[depth-1][i][k].eval()*Gl_[depth+1][k][j].eval()));
//                            lb = min(lb, (Ll_[depth-1][i][k].eval()*Gu_[depth+1][k][j].eval()));
//                            lb = min(lb, (Lu_[depth-1][i][k].eval()*Gl_[depth+1][k][j].eval()));
//                            lb = min(lb, (Lu_[depth-1][i][k].eval()*Gu_[depth+1][k][j].eval()));
//                            ub = max(ub, (Ll_[depth-1][i][k].eval()*Gl_[depth+1][k][j].eval()));
//                            ub = max(ub, (Ll_[depth-1][i][k].eval()*Gu_[depth+1][k][j].eval()));
//                            ub = max(ub, (Lu_[depth-1][i][k].eval()*Gl_[depth+1][k][j].eval()));
//                            ub = max(ub, (Lu_[depth-1][i][k].eval()*Gu_[depth+1][k][j].eval()));
//                            if (k%2==0) {
//                                G2 = G_[depth+1][k][j];
//                            }
//                            else {
//                                if (j%2==0) {//Imaginary
//                                    G2 = -1*G_[depth+1][k-1][j+1];
//                                }
//                                else {
//                                    G2 = G_[depth+1][k-1][j-1];
//                                }
//                            }
//                            LiftedGate -= L_[depth-1][i][k]*G2;
//                            Ll_[depth][i][j] = Ll_[depth][i][j].eval() + lb;
//                            Lu_[depth][i][j] = Lu_[depth][i][j].eval() + ub;
//                            if(Ll_[depth][i][j].eval()<-1){
//                                Ll_[depth][i][j] = -1;
//                            }
//                            if(Lu_[depth][i][j].eval()>1){
//                                Lu_[depth][i][j] = 1;
//                            }
//                        }
//                    }
//                    Qdesign.add(LiftedGate==0);
//                }
//            }
//        }
//
//        /** Target gate constraint **/
//        for (int i = 0; i < 2*m; i+=2) {
//            for (int j = 0; j < 2*m; j++) {
//                Constraint TargetGate("TargetGate_"+to_string(i+1)+to_string(j+1));
//                for (int k = 0; k < 2*m; k++) {
//                    if (!zeros[k][j]) {
//                        if (k%2==0) {
//                            G2 = G_[d-1][k][j];
//                        }
//                        else {
//                            if (j%2==0) {//Imaginary
//                                G2 = -1*G_[d-1][k-1][j+1];
//                            }
//                            else {
//                                G2 = G_[d-1][k-1][j-1];
//                            }
//                        }
//                        TargetGate += L_[d-3][i][k]*G2;
//                    }
//                }
//                if(!relax){
//                    TargetGate -=  slack_[i][j];
//                }
//                Qdesign.add(TargetGate==M.eval(i,j));
//            }
//        }
//    }
//
//
//    Constraint TOFF_Tt("TOFF_Tt");
//    TOFF_Tt += sum(zCnot);
//    Qdesign.add(TOFF_Tt<=1);
//
//    /** One gate per depth constraints **/
//    for(int depth = 0; depth < d ; depth++){
//        Constraint OneGate("OneGate_"+to_string(depth+1));
//        for(int qubit = 0; qubit < n ; qubit++){
//            if(ibm){
//                if(add_rx){
//                    OneGate += zS_[depth][qubit] + zSt_[depth][qubit] + zT_[depth][qubit] + zTt_[depth][qubit] + zRx_[depth][qubit];
//                }
//                else {
//                    OneGate += zS_[depth][qubit] + zSt_[depth][qubit] + zT_[depth][qubit] + zTt_[depth][qubit];
//                }
//            }
//            else if (cont){
//                if(add_rx){
//                    OneGate += zS_[depth][qubit] + zRx_[depth][qubit];
//                }
//                else {
//                    OneGate += zS_[depth][qubit];
//                }
//            }
//            else {
//                OneGate += zT_[depth][qubit] + zTt_[depth][qubit] + zH_[depth][qubit];
//            }
//            for (int qubit2 = qubit+1; qubit2 < n; qubit2++) {
//                OneGate += zCnot_[depth][qubit][qubit2];
//                if (rev_cnot) {
//                    OneGate += zCnotR_[depth][qubit][qubit2];
//                }
//            }
//        }
//        Qdesign.add(OneGate==1);
//    }
//    Qdesign.print();
//    auto solver_time_start = get_wall_time();
//    if(relax){
//        auto Relaxation = Qdesign.build_McCormick();
////        Relaxation->print_expanded();
//        solver slvr(*Relaxation,cplex);
//        slvr.run(output, relax = false, tol=1e-6, 0.01, "ma27");
//        Qdesign.print_solution(false);
//    }
//    else {
//        solver slvr(Qdesign,bonmin);
//        slvr.run(output=0, relax = false, tol=1e-6, 0.01, "ma27");
//        Qdesign.print_solution(false);
//    }
//    auto solver_time_end = get_wall_time();
//    auto solve_time = solver_time_end - solver_time_start;
//
//    auto total_time_start = get_wall_time();
//    auto total_time_end = get_wall_time();
//    auto total_time = total_time_end - total_time_start;
//    DebugOn("Solver Time = " << to_string(solve_time) << " seconds" << endl);
//    DebugOn("Total Computing Time = " << to_string(total_time) << " seconds" << endl);
//    DebugOn("Optimal Objective = " << to_string(Qdesign._obj_val) << endl);
//#else
//    cerr << "Error: this version of Gravity "
//    "was compiled without QPP support. Please rerun cmake with -DENABLE_QPP=1." << endl;
//#endif
//}
//int main (int argc, char * argv[])
//{
//    int output = 0;
//    //  Start Timers
//    string path = argv[0];
//
//
//    string mehrotra = "no", log_level="5", depth_str="2", nb_qubits_str="2", relax = "n";
//
//    /** create a OptionParser with options */
//    op::OptionParser opt;
//    opt.add_option("h", "help", "shows option help"); // no default value means boolean options, which default value is false
//    opt.add_option("l", "log", "Log level (def. 5)", log_level );
//    opt.add_option("r", "relax", "Convex Relaxation y/n (def. n)", relax );
//    opt.add_option("d", "depth", "Circuit depth", depth_str );
//    opt.add_option("n", "nb_qubits", "Number of Qubits", nb_qubits_str );
//
//    /** parse the options and verify that all went well. If not, errors and help will be shown */
//    bool correct_parsing = opt.parse_options(argc, argv);
//
//    if(!correct_parsing){
//        return EXIT_FAILURE;
//    }
//
//    output = op::str2int(opt["l"]);
//    bool has_help = op::str2bool(opt["h"]);
//    /** show help */
//    if(has_help) {
//        opt.show_help();
//        exit(0);
//    }
//    const int n = op::str2int(opt["n"]);
//    const int d = op::str2int(opt["d"]);
//    const int m = pow(2,n);
//
//#ifdef USE_QPP
//    /* Target Matrix */
//    param<> M("M");
//    M.set_size(2*m, 2*m);
//    if(true){
//        using namespace qpp;
//        auto Mp = gt.expandout(gt.CTRL(gt.Z, {0}, {1}, n), 0, 1, m);
////        auto Mp = gt.expandout(gt.Fd(m), 0, 1, m);
////        auto Mp = gt.expandout(gt.TOF, 0, 1, m);
//        DebugOn("Target matrix = " << endl);
//        DebugOn(disp(Mp) << "\n");
//        M.set_vals(Mp.sparseView());
//        DebugOff(M.to_str(true));
//        if (opt["r"]=="y") {
//            run_qtdesign(n, d, M, true);
//        }
//        else {
//            run_qtdesign(n, d, M, false);
//        }
//    }
//
//#endif
//    return 0;
//}
