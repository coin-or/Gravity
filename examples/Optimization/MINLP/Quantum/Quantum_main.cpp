//
//
//  Gravity
//
//  Created by Hassan Hijazi on August 29 2018.
//
//


#include <iostream>
#include <string>
#include <fstream>
#include <gravity/solver.h>
#include <stdio.h>
#include <stdlib.h>
#include <optionParser.hpp>


using namespace std;
using namespace gravity;



param<> run_reverse(const int n, const int d, const param<>& M){
    param<> Mres("Mres");
#ifdef USE_QPP
    double tol = 1e-6;
    bool ibm = false;
    bool cont = true;
    bool add_rx = true;
    bool rev_cnot = false;
    int output = 0;
    bool relax = false;
    if (ibm && cont) {
        cerr << "Choose IBM or Continuous, or none.";
        exit(-1);
    }
    
    /** Parameters **/
    
    /* Qubit number */
    const int m = pow(2,n);
    DebugOn("Number of Qubits = " << n << endl);
    /* Circuit depth */
    DebugOn("Depth = " << d << endl);
    if (d<2) {
        cerr << "Depth should be greater or equal than 2!" << endl;
        exit(-1);
    }

    /* S gate matrices, one for each row (qubit) */
    param<> S[n];
    
    /* Rx gate matrices */
    param<> Rx[n];
    
    /* Cnot gate matrices */
    param<> Cnot[n*(n-1)/2];
    param<> CnotR[n*(n-1)/2];
    
    
    /** Initialization **/
    int idx = 0;
    for (auto i = 0; i<n; i++) { // Qubit (row) iterator
        Rx[i].set_name("Rx_"+to_string(i+1));
        Rx[i].set_size(2*m, 2*m, 0);
        Rx[i].QuantumRx(i,n);/* Fill nonzero values corresponding to an H gate */
        S[i].set_name("S_"+to_string(i+1));
        S[i].set_size(2*m, 2*m, 0);
        S[i].QuantumS(i,n);/* Fill nonzero values corresponding to an H gate */
        for (auto k = i+1; k<n; k++) {
            Cnot[idx].set_name("Cnot_"+to_string(i+1)+to_string(k+1));
            Cnot[idx].set_size(2*m, 2*m, 0);
            Cnot[idx].QuantumCnot(i,k,n);/* Fill nonzero values corresponding to a CNOT gate */
            CnotR[idx].set_name("CnotR_"+to_string(k+1)+to_string(i+1));
            CnotR[idx].set_size(2*m, 2*m, 0);
            CnotR[idx].QuantumCnot(k,i,n);/* Fill nonzero values corresponding to a CNOT gate */
            idx++;
        }
    }
    
    
    /** Indices **/
    auto d_ids = indices(1,d);/*< Depth indices */
    auto q_ids = indices(1,n);/*< Qubit indices */
    auto pairs = indices(ordered_pairs(2, n+1)); /*< Indices pairs (control-taget) for Cnot gates */
    auto rpairs = indices(ordered_pairs(2, n+1, true)); /*< Indices pairs (control-taget) for Cnot gates */
    auto m2 = indices(indices(1,2*m), indices(1,2*m)); /*< Indices for square mxm matrices */
    
    
    /** Variables **/
    /** Binaries **/
    var<bool> zS("zS"); /*< Binary variables for S gates */
    var<bool> zRx("zRx"); /*< Binary variables for Rx gates */
    var<bool> zCnot("zCnot"); /*< Binary variables for CNOT gates */
    var<bool> zCnotR("zCnotR"); /*< Binary variables for CNOT gates */

    
    /** Gate Matrices **/
    var<> G[d]; /*< Base Gate Matrices */
    bool zeros[2*m][2*m]; /*< Zeros in Gate Matrices */
    var<> L[d-2]; /*< Lifted Product Matrices */
    var<> lambda_r("lambda_r", -1, 1);
    var<> lambda_i("lambda_i", -1, 1);
    var<> slack("slack", -10, 10);
    
    
    
    /** Model **/
    
    /** Adding Variables **/
    Model Qdesign("Qdesign");
    /** Binaries **/
    
    Qdesign.add(lambda_r.in(d_ids), lambda_i.in(d_ids));
    //        Qdesign.add(lambda_i.in(d_ids));
    if (rev_cnot) {
        Qdesign.add(zS.in(d_ids,q_ids), zRx.in(d_ids,q_ids) ,zCnot.in(d_ids,pairs), zCnotR.in(d_ids,rpairs));
    }
    else {
        Qdesign.add(zS.in(d_ids,q_ids), zRx.in(d_ids,q_ids) ,zCnot.in(d_ids,pairs));
    }

    Qdesign.add(slack.in(m2));
    DebugOff("Total number of variables after adding the binaries = " << Qdesign.get_nb_vars() << endl);
    /** Gate Matrices **/
    for (int depth = 0; depth < d; depth++) {
        G[depth] = var<>("G_"+to_string(depth+1), -1,1);//TODO automate bounds on G
        Qdesign.add(G[depth].in(m2));
    }
    if (d>2) {
        for (int depth = 0; depth < d-2; depth++) {
            L[depth] = var<>("L_"+to_string(depth+1), -10,10);
            Qdesign.add(L[depth].in(m2));
        }
    }
    
    DebugOff("Total number of variables after adding the G vars= " << Qdesign.get_nb_vars() << endl);
    
    
    /** Adding Constraints **/
    /** Building indexed variables as matrices (for a faster access in constraints) **/
    var<bool> zS_[d][n], zRx_[d][n], zCnot_[d][n][n], zCnotR_[d][n][n];
    var<> G_[d][2*m][2*m], L_[d-2][2*m][2*m];
    var<> slack_[2*m][2*m];
    var<> lambda_r_[d], lambda_i_[d];
    for (int i = 0; i < 2*m; i++) {
        for (int j = 0; j < 2*m; j++) {
            slack_[i][j] = slack(i+1,j+1);
        }
    }
    for (size_t depth = 0; depth < d; depth++) {
        lambda_r_[depth] = lambda_r(depth+1);
        lambda_i_[depth] = lambda_i(depth+1);
        for (size_t qubit1 = 0; qubit1 < n; qubit1++) {
            zS_[depth][qubit1] = zS(depth+1,qubit1+1);
            zRx_[depth][qubit1] = zRx(depth+1,qubit1+1);
            for (size_t qubit2 = qubit1+1; qubit2 < n; qubit2++) {
                zCnot_[depth][qubit1][qubit2] = zCnot(depth+1,qubit1+1,qubit2+1);
                zCnotR_[depth][qubit1][qubit2] = zCnotR(depth+1,qubit2+1,qubit1+1);
            }
        }
        for (int i = 0; i < 2*m; i++) {
            for (int j = 0; j < 2*m; j++) {
                G_[depth][i][j] = G[depth](i+1,j+1);
            }
        }
    }
    if (d>2) {
        for (size_t depth = 0; depth < d-2; depth++) {
            for (int i = 0; i < 2*m; i++) {
                for (int j = 0; j < 2*m; j++) {
                    L_[depth][i][j] = L[depth](i+1,j+1);
                }
            }
        }
    }
    
    /** Objective **/
    func_ obj;
    for (int i = 0; i < 2*m; i++) {
        for (int j = 0; j < 2*m; j++) {
            if (i==j) {
                obj += (1-slack_[i][j])*(1-slack_[i][j]);
            }
            else {
                obj += slack_[i][j]*slack_[i][j];
            }
        }
    }
    
    Qdesign.min(obj);
    
    
    /** Base gate matrix definition **/
    for (int depth = 0; depth < d; depth++) {
        Constraint Complex("Complex"+to_string(depth+1));
        Complex += lambda_i_[depth]*lambda_i_[depth] + lambda_r_[depth]*lambda_r_[depth];
        Qdesign.add(Complex>=1);
        
        for (int i = 0; i < 2*m; i++) {
            for (int j = 0; j < 2*m; j++) {
                Constraint BaseGate("BaseGate_"+to_string(depth+1)+"_"+to_string(i+1)+to_string(j+1));
                if (i%2==0) {
                    BaseGate = G_[depth][i][j];
                }
                else {
                    if (j%2==0) {//Imaginary
                        BaseGate = -1*G_[depth][i-1][j+1];
                    }
                    else {
                        BaseGate = G_[depth][i-1][j-1];
                    }
                }
                int idx = 0;
                for (int qubit = 0; qubit < n; qubit++) {
                    if ((i%2!=0 && j%2==0) || (i%2==0 && j%2!=0)) {//imaginary part
                            BaseGate -= S[qubit].eval(i,j)*lambda_i_[depth]*zS_[depth][qubit];
                    }
                    else {
                        if (i%2==0) {
                            if(S[qubit].eval(i,j+1)!=0){
                                BaseGate -= lambda_r_[depth]*zS_[depth][qubit] + zRx_[depth][qubit]*Rx[qubit].eval(i,j);
                            }
                            else {
                                BaseGate -= zS_[depth][qubit]*S[qubit].eval(i,j) + zRx_[depth][qubit]*Rx[qubit].eval(i,j);
                            }
                        }
                        else {
                            if(S[qubit].eval(i,j-1)!=0){
                                BaseGate -= lambda_r_[depth]*zS_[depth][qubit] + zRx_[depth][qubit]*Rx[qubit].eval(i,j);
                            }
                            else {
                                BaseGate -= zS_[depth][qubit]*S[qubit].eval(i,j) + zRx_[depth][qubit]*Rx[qubit].eval(i,j);
                            }
                        }
                    }
                    for (int qubit2 = qubit+1; qubit2 < n; qubit2++) {
                        BaseGate -= zCnot_[depth][qubit][qubit2]*Cnot[idx].eval(i,j);
                        if (rev_cnot) {
                            BaseGate -= zCnotR_[depth][qubit][qubit2]*CnotR[idx].eval(i,j);
                        }
                        idx++;
                    }
                }
                if (BaseGate.get_nb_vars()==1) { /* This implies that the corresponding entry is zero */
                    zeros[i][j]=true;
                }
                else {
                    zeros[i][j]=false;
                    Qdesign.add(BaseGate==0);
                }
            }
        }
    }
    
    /** Lifted + Target gate constraint **/
    /* No need to use lifted matrices if the depth = 2 */
    if (d==2) {
        for (int i = 0; i < 2*m; i+=2) {
            for (int j = 0; j < 2*m; j++) {
                Constraint TargetGate("TargetGate_"+to_string(i+1)+to_string(j+1));
                for (int k = 0; k < 2*m; k++) {
                    if (!zeros[i][k] && !zeros[k][j]) {
                        TargetGate += G_[0][i][k]*G_[1][k][j];
                    }
                }
                    Qdesign.add(TargetGate==M.eval(i,j));
            }
        }
    }
    else {
        /** Lifted gate constraint **/
        /** Depth 0 **/
        for (int i = 0; i < 2*m; i+=2) {
            for (int j = 0; j < 2*m; j++) {
                Constraint LiftedGate("LiftedGate_0_"+to_string(i+1)+to_string(j+1));
                LiftedGate += L_[0][i][j];
                for (int k = 0; k < 2*m; k++) {
                    if (!zeros[i][k] && !zeros[k][j]) {
                        LiftedGate += G_[0][i][k]*G_[1][k][j];
                    }
                }
                Qdesign.add(LiftedGate==0);
            }
        }
        for (int depth = 1; depth < d-2; depth++) {
            for (int i = 0; i < 2*m; i+=2) {
                for (int j = 0; j < 2*m; j++) {
                    Constraint LiftedGate("LiftedGate"+to_string(depth)+"_"+to_string(i+1)+to_string(j+1));
                    LiftedGate += L_[depth][i][j];
                    for (int k = 0; k < 2*m; k++) {
                        if (!zeros[k][j]) {
                            LiftedGate -= L_[depth-1][i][k]*G_[depth+1][k][j];
                        }
                    }
                    Qdesign.add(LiftedGate==0);
                }
            }
        }
        
        /** Target gate constraint **/
        for (int i = 0; i < 2*m; i+=2) {
            for (int j = 0; j < 2*m; j++) {
                Constraint TargetGate("TargetGate_"+to_string(i+1)+to_string(j+1));
                for (int k = 0; k < 2*m; k++) {
                    if (!zeros[k][j]) {
                        TargetGate += L_[d-3][i][k]*G_[d-1][k][j];
                    }
                }
//                if (TargetGate.get_nb_vars()>1) {
                    Qdesign.add(TargetGate==M.eval(i,j));
//                }
                
            }
        }
    }
    
//    for (size_t depth = 0; depth < d-1; depth++) {
//        for (size_t qubit1 = 0; qubit1 < n; qubit1++) {
////            Constraint NonConseqS("NonConseq_S_" +to_string(depth)+"_"+to_string(qubit1+1));
////            NonConseqS += zS_[depth][qubit1] + zS_[depth+1][qubit1];
////            Qdesign.add(NonConseqS <= 1);
//            for (size_t qubit2 = qubit1+1; qubit2 < n; qubit2++) {
//                Constraint NonConseqCnot("NonConseq_Cnot_" +to_string(depth)+"_"+to_string(qubit1+1)+"_"+to_string(qubit2+1));
//                NonConseqCnot += zCnot_[depth][qubit1][qubit2] + zCnot_[depth+1][qubit1][qubit2];
//                Qdesign.add(NonConseqCnot <= 1);
//            }
//        }
//    }
//    Constraint Test("Test");
//    Test += sum(zCnot);
//    Qdesign.add(Test>=1);
    
    
    /** One gate per depth constraints **/
    for(int depth = 0; depth < d ; depth++){
        Constraint OneGate("OneGate_"+to_string(depth+1));
        for(int qubit = 0; qubit < n ; qubit++){
            OneGate += zS_[depth][qubit] + zRx_[depth][qubit];
            for (int qubit2 = qubit+1; qubit2 < n; qubit2++) {
                OneGate += zCnot_[depth][qubit][qubit2];
                if (rev_cnot) {
                    OneGate += zCnotR_[depth][qubit][qubit2];
                }
            }
        }
        Qdesign.add(OneGate==1);
    }
    Qdesign.print_expanded();
    solver slvr(Qdesign,bonmin);
    auto solver_time_start = get_wall_time();
    slvr.run(output, relax = false, tol=1e-6, 0.01, "ma27");
    Qdesign.print_solution();
    slack.gravity::param<double>::print(true);
    lambda_r.gravity::param<double>::print(true);
    lambda_i.gravity::param<double>::print(true);
    auto solver_time_end = get_wall_time();
    auto solve_time = solver_time_end - solver_time_start;
    
    auto total_time_start = get_wall_time();
    auto total_time_end = get_wall_time();
    auto total_time = total_time_end - total_time_start;
    DebugOn("Solver Time = " << to_string(solve_time) << " seconds" << endl);
    DebugOn("Total Computing Time = " << to_string(total_time) << " seconds" << endl);
    Mres.set_vals(slack);
#else
    cerr << "Error: this version of Gravity "
    "was compiled without QPP support. Please rerun cmake with -DENABLE_QPP=1." << endl;
#endif
    return Mres;
}

void run_qtdesign_fixed(const int n, const int d, const param<>& M){
#ifdef USE_QPP
    double tol = 1e-6;
    bool ibm = false;
    bool cont = true;
    bool add_rx = true;
    bool rev_cnot = false;
    int output = 0;
    bool relax = false;
    if (ibm && cont) {
        cerr << "Choose IBM or Continuous, or none.";
        exit(-1);
    }
    
    /** Parameters **/
    
    /* Qubit number */
    const int m = pow(2,n);
    DebugOn("Number of Qubits = " << n << endl);
    /* Circuit depth */
    DebugOn("Depth = " << d << endl);
    if (d<2) {
        cerr << "Depth should be greater or equal than 2!" << endl;
        return;
    }
    
    /* U1 gate matrices, one for each row (qubit) */
    param<> U1[n];
    
    /* Rx gate matrices */
    param<> Rx[n];
    
    /* Cnot gate matrices */
    param<> Cnot[n*(n-1)/2];
    param<> CnotR[n*(n-1)/2];

    
//    /* Target Matrix */
//    param<> M("M");
//    M.set_size(2*m, 2*m, 0);
//    {
//        using namespace qpp;
//        double theta = pi/2;
//        auto M1 = gt.expandout(gt.Rn(theta, {0,0,1}), 0, n);
//        auto Mt1 = M1 * complex<double>(cos(theta/2),sin(theta/2));
//        auto Mt2 = gt.expandout(gt.Rn(pi/2, {1,0,0}), 0, n);
//        auto Mt3 = gt.expandout(gt.CTRL(gt.X, {1}, {0}, n), 0, 1, m);
//        theta = 3*pi/4;
//        auto M4 = gt.expandout(gt.Rn(theta, {0,0,1}), 0, n);
//        auto Mt4 = M4 * complex<double>(cos(theta/2),sin(theta/2));
//        theta = -pi/4;
//        auto M5 = gt.expandout(gt.Rn(theta, {0,0,1}), 1, n);
//        auto Mt5 = M5 * complex<double>(cos(theta/2),sin(theta/2));
//        auto Mt6 = gt.expandout(gt.CTRL(gt.X, {0}, {1}, n), 0, 1, m);
//        auto Mt7 = gt.expandout(gt.Rn(pi/2, {1,0,0}), 0, n);
//        auto Mt8 = gt.expandout(gt.CTRL(gt.X, {1}, {0}, n), 0, 1, m);
//        theta = pi/2;
//        auto M9 = gt.expandout(gt.Rn(theta, {0,0,1}), 0, n);
//        auto Mt9 = M9 * complex<double>(cos(theta/2),sin(theta/2));
//        theta = 7*pi/4;
//        auto M10 = gt.expandout(gt.Rn(theta, {0,0,1}), 0, n);
//        auto Mt10 = M10 * complex<double>(cos(theta/2),sin(theta/2));
//        auto Mp = Mt1*Mt2*Mt3*Mt4*Mt5*Mt6*Mt7*Mt8*Mt9*Mt10;
//        //        DebugOn("M5 matrix at position " << to_string(1) <<" = " << endl);
//        //        DebugOn(disp(M5) << "\n");
//        
////        auto Mp = gt.expandout(gt.Fd(m), 0, 1, m);
////        Mp *= sqrt(2);
//        DebugOn("Target matrix = " << endl);
//        DebugOn(disp(Mp) << "\n");
//        M.set_vals(Mp.sparseView());
//        DebugOff(M.to_str(true));
//    }
    /** Initialization **/
    int idx = 0;
    for (auto i = 0; i<n; i++) { // Qubit (row) iterator
        Rx[i].set_name("Rx_"+to_string(i+1));
        Rx[i].set_size(2*m, 2*m, 0);
        Rx[i].QuantumRx(i,n);/* Fill nonzero values corresponding to an H gate */
        U1[i].set_name("U1_"+to_string(i+1));
        U1[i].set_size(2*m, 2*m, 0);
        U1[i].QuantumS(i,n);/* Fill nonzero values corresponding to an H gate */
        for (auto k = i+1; k<n; k++) {
            Cnot[idx].set_name("Cnot_"+to_string(i+1)+to_string(k+1));
            Cnot[idx].set_size(2*m, 2*m, 0);
            Cnot[idx].QuantumCnot(i,k,n);/* Fill nonzero values corresponding to a CNOT gate */
            CnotR[idx].set_name("CnotR_"+to_string(k+1)+to_string(i+1));
            CnotR[idx].set_size(2*m, 2*m, 0);
            CnotR[idx].QuantumCnot(k,i,n);/* Fill nonzero values corresponding to a CNOT gate */
            idx++;
        }
    }
    
    
    /** Indices **/
    auto d_ids = indices(1,d);/*< Depth indices */
    auto q_ids = indices(1,n);/*< Qubit indices */
    auto m2 = indices(indices(1,2*m), indices(1,2*m)); /*< Indices for square mxm matrices */
    
    
    /** Variables **/
    /** Binaries **/
    var<bool> zU1("zU1"); /*< Binary variables for S gates */
    var<bool> zRx("zRx"); /*< Binary variables for Rx gates */
    
    
    /** Gate Matrices **/
    var<> G[d]; /*< Base Gate Matrices */
    bool zeros[2*m][2*m]; /*< Zeros in Gate Matrices */
    var<> L[d-2]; /*< Lifted Product Matrices */

    var<> lambda_r("lambda_r", -1, 1);
    var<> lambda_i("lambda_i", -1, 1);
    var<> slack("slack",-10,10);
    
    
    
    /** Model **/
    
    /** Adding Variables **/
    Model Qdesign("Qdesign");
    Qdesign.add(lambda_i.in(d_ids), lambda_r.in(d_ids));
    Qdesign.add(zU1.in(d_ids,q_ids), zRx.in(d_ids,q_ids));
    Qdesign.add(slack.in(m2));
    DebugOff("Total number of variables after adding the binaries = " << Qdesign.get_nb_vars() << endl);
    /** Gate Matrices **/
    for (int depth = 0; depth < d; depth++) {
        G[depth] = var<>("G_"+to_string(depth+1), -1,1);//TODO automate bounds on G
        if (d==2 || d==5 || d==7) {
            continue;
        }
        Qdesign.add(G[depth].in(m2));
    }
    for (int depth = 0; depth < d-2; depth++) {
        L[depth] = var<>("L_"+to_string(depth+1), -10,10);
        Qdesign.add(L[depth].in(m2));
    }
    
    DebugOff("Total number of variables after adding the G vars= " << Qdesign.get_nb_vars() << endl);
    
    
    /** Adding Constraints **/
    /** Building indexed variables as matrices (for a faster access in constraints) **/
    var<bool> zU1_[d][n], zRx_[d][n];
    var<> G_[d][2*m][2*m], L_[d-2][2*m][2*m];
    var<> slack_[2*m][2*m];
    var<> lambda_r_[d], lambda_i_[d];
    for (int i = 0; i < 2*m; i++) {
        for (int j = 0; j < 2*m; j++) {
            slack_[i][j] = slack(i+1,j+1);
        }
    }
    for (size_t depth = 0; depth < d; depth++) {
        lambda_r_[depth] = lambda_r(depth+1);
        lambda_i_[depth] = lambda_i(depth+1);
        for (size_t qubit1 = 0; qubit1 < n; qubit1++) {
            zU1_[depth][qubit1] = zU1(depth+1,qubit1+1);
            zRx_[depth][qubit1] = zRx(depth+1,qubit1+1);
        }
        for (int i = 0; i < 2*m; i++) {
            for (int j = 0; j < 2*m; j++) {
                G_[depth][i][j] = G[depth](i+1,j+1);
            }
        }
    }
    for (size_t depth = 0; depth < d-2; depth++) {
        for (int i = 0; i < 2*m; i++) {
            for (int j = 0; j < 2*m; j++) {
                L_[depth][i][j] = L[depth](i+1,j+1);
            }
        }
    }
    
    /** Objective **/
    func_ obj;
    for (int i = 0; i < 2*m; i++) {
        for (int j = 0; j < 2*m; j++) {
            obj += slack_[i][j]*slack_[i][j];
        }
    }
    Qdesign.min(obj);
    
    
    /** Base gate matrix definition **/
    for (int depth = 0; depth < d; depth++) {
        if (depth==2 || depth==5 || depth==7) {
            continue;
        }

        for (int i = 0; i < 2*m; i++) {
            for (int j = 0; j < 2*m; j++) {
                Constraint BaseGate("BaseGate_"+to_string(depth+1)+"_"+to_string(i+1)+to_string(j+1));
                if (i%2==0) {
                    BaseGate = G_[depth][i][j];
                }
                else {
                    if (j%2==0) {//Imaginary
                        BaseGate = -1*G_[depth][i-1][j+1];
                    }
                    else {
                        BaseGate = G_[depth][i-1][j-1];
                    }
                }
                for (int qubit = 0; qubit < n; qubit++) {
                    if ((i%2!=0 && j%2==0) || (i%2==0 && j%2!=0)) {//imaginary part
                        BaseGate -= U1[qubit].eval(i,j)*lambda_i_[depth]*zU1_[depth][qubit] + zRx_[depth][qubit]*Rx[qubit].eval(i,j);
                    }
                    else {
                        if (i%2==0) {
                            if(U1[qubit].eval(i,j+1)!=0){
                                BaseGate -= lambda_r_[depth]*zU1_[depth][qubit] + zRx_[depth][qubit]*Rx[qubit].eval(i,j);
                            }
                            else {
                                BaseGate -= zU1_[depth][qubit]*U1[qubit].eval(i,j) + zRx_[depth][qubit]*Rx[qubit].eval(i,j);
                            }
                        }
                        else {
                            if(U1[qubit].eval(i,j-1)!=0){
                                BaseGate -= lambda_r_[depth]*zU1_[depth][qubit] + zRx_[depth][qubit]*Rx[qubit].eval(i,j);
                            }
                            else {
                                BaseGate -= zU1_[depth][qubit]*U1[qubit].eval(i,j) + zRx_[depth][qubit]*Rx[qubit].eval(i,j);
                            }
                        }
                    }
                }
                if (BaseGate.get_nb_vars()==1) { /* This implies that the corresponding entry is zero */
                    zeros[i][j]=true;
                }
                else {
                    zeros[i][j]=false;
                    Qdesign.add(BaseGate==0);
                    //                    BaseGate.print();
                }
            }
        }
    }
    
    /** Lifted + Target gate constraint **/
    /* No need to use lifted matrices if the depth = 2 */
    
    /** Lifted gate constraint **/
    /** Depth 0 **/
    for (int i = 0; i < 2*m; i+=2) {
        for (int j = 0; j < 2*m; j++) {
            Constraint LiftedGate("LiftedGate_0_"+to_string(i+1)+to_string(j+1));
            LiftedGate += L_[0][i][j];
            for (int k = 0; k < 2*m; k++) {
                if (!zeros[i][k] && !zeros[k][j]) {
                    LiftedGate -= G_[0][i][k]*G_[1][k][j];
                }
            }
            Qdesign.add(LiftedGate==0);
        }
    }
    for (int depth = 1; depth < d-2; depth++) {
        if (depth==1 || depth==4 || depth==6) {
            continue;
        }
        for (int i = 0; i < 2*m; i+=2) {
            for (int j = 0; j < 2*m; j++) {
                Constraint LiftedGate("LiftedGate"+to_string(depth)+"_"+to_string(i+1)+to_string(j+1));
                LiftedGate += L_[depth][i][j];
                for (int k = 0; k < 2*m; k++) {
                    if (!zeros[k][j]) {
                        LiftedGate -= L_[depth-1][i][k]*G_[depth+1][k][j];
                    }
                }
                Qdesign.add(LiftedGate==0);
            }
        }
    }
    
    //d = 1
    for (int i = 0; i < 2*m; i+=2) {
        for (int j = 0; j < 2*m; j++) {
            Constraint LiftedGate("LiftedGate"+to_string(1)+"_"+to_string(i+1)+to_string(j+1));
            LiftedGate += L_[1][i][j];
            for (int k = 0; k < 2*m; k++) {
                    LiftedGate -= L_[0][i][k]*CnotR[0].eval(k,j);
            }
            Qdesign.add(LiftedGate==0);
        }
    }
    
    //d = 4
    for (int i = 0; i < 2*m; i+=2) {
        for (int j = 0; j < 2*m; j++) {
            Constraint LiftedGate("LiftedGate"+to_string(4)+"_"+to_string(i+1)+to_string(j+1));
            LiftedGate += L_[4][i][j];
            for (int k = 0; k < 2*m; k++) {
                LiftedGate -= L_[3][i][k]*Cnot[0].eval(k,j);
            }
            Qdesign.add(LiftedGate==0);
        }
    }
    
    //d = 6
    for (int i = 0; i < 2*m; i+=2) {
        for (int j = 0; j < 2*m; j++) {
            Constraint LiftedGate("LiftedGate"+to_string(6)+"_"+to_string(i+1)+to_string(j+1));
            LiftedGate += L_[6][i][j];
            for (int k = 0; k < 2*m; k++) {
                LiftedGate -= L_[5][i][k]*CnotR[0].eval(k,j);
            }
            Qdesign.add(LiftedGate==0);
        }
    }
    
    /** Target gate constraint **/
    for (int i = 0; i < 2*m; i+=2) {
        for (int j = 0; j < 2*m; j++) {
            Constraint TargetGate("TargetGate_"+to_string(i+1)+to_string(j+1));
            for (int k = 0; k < 2*m; k++) {
                if (!zeros[k][j]) {
                    TargetGate += L_[d-3][i][k]*G_[d-1][k][j];
                }
            }
            TargetGate -=  slack_[i][j];
//            if (TargetGate.get_nb_vars()>1) {
                Qdesign.add(TargetGate==M.eval(i,j));
                //                    Qdesign.add(TargetGate==0);
//            }
            
        }
    }
    

//    Constraint D2("D2");
//    D2 += zU1_[0][0];
//    Qdesign.add(D2==1);
//    
//    Constraint D5("D5");
//    D5 += zRx_[1][0];
//    Qdesign.add(D5==1);
//
//    Constraint D7("D7");
//    D7 += zU1_[3][0];
//    Qdesign.add(D7==1);
//    
//    Constraint D8("D8");
//    D8 += zU1_[4][1];
//    Qdesign.add(D8==1);
//    
//    Constraint D9("D9");
//    D9 += zRx_[6][0];
//    Qdesign.add(D9==1);
//    
//    Constraint D10("D10");
//    D10 += zU1_[8][0];
//    Qdesign.add(D10==1);
//
//    Constraint D11("D11");
//    D11 += zU1_[9][1];
//    Qdesign.add(D11==1);
    
    Constraint D12("D12");
    D12 += zU1_[2][1] + zU1_[2][0] + zRx_[2][1] + zRx_[2][0];
    Qdesign.add(D12==0);
    
    Constraint D13("D13");
    D13 += zU1_[5][1] + zU1_[5][0] + zRx_[5][1] + zRx_[5][0];
    Qdesign.add(D13==0);
    
    Constraint D14("D14");
    D14 += zU1_[7][1] + zU1_[7][0] + zRx_[7][1] + zRx_[7][0];
    Qdesign.add(D14==0);

    
    /** One gate per depth constraints **/
    for(int depth = 0; depth < d ; depth++){
        if (depth==2 || depth==5 || depth==7) {
            continue;
        }
        Constraint OneGate("OneGate_"+to_string(depth+1));
        for(int qubit = 0; qubit < n ; qubit++){
            OneGate += zU1_[depth][qubit] + zRx_[depth][qubit];
        }
        Qdesign.add(OneGate==1);
    }
    //    Qdesign.print_expanded();
    solver slvr(Qdesign,bonmin);
    auto solver_time_start = get_wall_time();
    slvr.run(5, false, 1e-6, 0.01, "ma27");
    Qdesign.print_solution();
    slack.gravity::param<double>::print(true);
    lambda_r.gravity::param<double>::print(true);
    lambda_i.gravity::param<double>::print(true);
    auto solver_time_end = get_wall_time();
    auto solve_time = solver_time_end - solver_time_start;
    
    auto total_time_start = get_wall_time();
    auto total_time_end = get_wall_time();
    auto total_time = total_time_end - total_time_start;
    DebugOn("Solver Time = " << to_string(solve_time) << " seconds" << endl);
    DebugOn("Total Computing Time = " << to_string(total_time) << " seconds" << endl);
#else
    cerr << "Error: this version of Gravity "
    "was compiled without QPP support. Please rerun cmake with -DENABLE_QPP=1." << endl;
#endif
}

void run_qtdesign_polar(const int n, const int d, const param<>& M, bool relax){
#ifdef USE_QPP
    double tol = 1e-6;
    bool ibm = false;
    bool cont = true;
    bool add_rx = true;
    bool rev_cnot = false;
    int output = 0;
    if (ibm && cont) {
        cerr << "Choose IBM or Continuous, or none.";
        exit(-1);
    }

    /** Parameters **/
    
    /* Qubit number */
    const int m = pow(2,n);
    DebugOn("Number of Qubits = " << n << endl);
    /* Circuit depth */
    DebugOn("Depth = " << d << endl);
    if (d<2) {
        cerr << "Depth should be greater or equal than 2!" << endl;
        return;
    }
    /* T and T conjugate transpose gate matrices, one for each row (qubit) */
    param<> T[n], Tt[n];
    
    /* S and St conjugate transpose gate matrices, one for each row (qubit) */
    param<> S[n], St[n];
    
    /* Hadamard gate matrices */
    param<> H[n];
    
    /* Rx gate matrices */
    param<> Rx[n];
    
    /* Cnot gate matrices */
    param<> Cnot[n*(n-1)/2];
    param<> CnotR[n*(n-1)/2];
    /* SWAP gate matrices */
    param<> Swap[n*(n-1)/2];
    
    
    /** Initialization **/
    int idx = 0;
    for (auto i = 0; i<n; i++) { // Qubit (row) iterator
        T[i].set_name("T_"+to_string(i+1));
        T[i].set_size(2*m, 2*m, 0);/* Representing Complex number as 2x2 real matrices */
        T[i].QuantumT(i,n);/* Fill nonzero values corresponding to a T gate */
        Tt[i].set_name("Tt_"+to_string(i+1));
        Tt[i].set_size(2*m, 2*m, 0);
        Tt[i].QuantumT(i,n,true);/* Fill nonzero values corresponding to a T conjugate transpose gate */
        Rx[i].set_name("Rx_"+to_string(i+1));
        Rx[i].set_size(2*m, 2*m, 0);
        Rx[i].QuantumRx(i,n);/* Fill nonzero values corresponding to an H gate */
        S[i].set_name("S_"+to_string(i+1));
        S[i].set_size(2*m, 2*m, 0);
        S[i].QuantumS(i,n);/* Fill nonzero values corresponding to an H gate */
        St[i].set_name("St_"+to_string(i+1));
        St[i].set_size(2*m, 2*m, 0);
        St[i].QuantumSt(i,n);/* Fill nonzero values corresponding to an H gate */
        H[i].set_name("H_"+to_string(i+1));
        H[i].set_size(2*m, 2*m, 0);
        H[i].QuantumH(i,n);/* Fill nonzero values corresponding to an H gate */
        for (auto k = i+1; k<n; k++) {
            Cnot[idx].set_name("Cnot_"+to_string(i+1)+to_string(k+1));
            Cnot[idx].set_size(2*m, 2*m, 0);
            Cnot[idx].QuantumCnot(i,k,n);/* Fill nonzero values corresponding to a CNOT gate */
            CnotR[idx].set_name("CnotR_"+to_string(k+1)+to_string(i+1));
            CnotR[idx].set_size(2*m, 2*m, 0);
            CnotR[idx].QuantumCnot(k,i,n);/* Fill nonzero values corresponding to a CNOT gate */
            Swap[idx].set_name("Swap"+to_string(i+1)+to_string(k+1));
            Swap[idx].set_size(2*m, 2*m, 0);
            Swap[idx].QuantumSwap(i,k,n);/* Fill nonzero values corresponding to a CNOT gate */
            
            idx++;
        }
    }
    
    
    /** Indices **/
    auto d_ids = indices(1,d);/*< Depth indices */
    auto q_ids = indices(1,n);/*< Qubit indices */
    auto pairs = indices(ordered_pairs(2, n+1)); /*< Indices pairs (control-taget) for Cnot gates */
    auto rpairs = indices(ordered_pairs(2, n+1, true)); /*< Indices pairs (control-taget) for Cnot gates */
    auto m2 = indices(indices(1,2*m), indices(1,2*m)); /*< Indices for square mxm matrices */
    
    
    /** Variables **/
    /** Binaries **/
    var<bool> zT("zT"); /*< Binary variables for T gates */
    var<bool> zTt("zTt"); /*< Binary variables for T conjugate gates */
    var<bool> zH("zH"); /*< Binary variables for H gates */
    var<bool> zS("zS"); /*< Binary variables for S gates */
    var<bool> zSt("zSt"); /*< Binary variables for S gates */
    var<bool> zRx("zRx"); /*< Binary variables for Rx gates */
    var<bool> zCnot("zCnot"); /*< Binary variables for CNOT gates */
    var<bool> zCnotR("zCnotR"); /*< Binary variables for CNOT gates */
    var<bool> zSwap("zSwap"); /*< Binary variables for SWAP gates */
    
    
    /** Gate Matrices **/
    var<> G[d]; /*< Base Gate Matrices */
    bool zeros[2*m][2*m]; /*< Zeros in Gate Matrices */
    var<> L[d-2]; /*< Lifted Product Matrices */
    double pi = 4.*atan(1.);
    var<> theta("theta", -2.*pi, 2.*pi);
    theta.initialize_all(pi/2.);
    var<> slack("slack", -10, 10);
    
    
    
    /** Model **/
    
    /** Adding Variables **/
    Model Qdesign("Qdesign");
    /** Binaries **/
    
    //    Qdesign.add(zT.in(d_ids,q_ids), zTt.in(d_ids,q_ids), zH.in(d_ids,q_ids), zCnot.in(d_ids,pairs));
    //    Qdesign.add(zT.in(d_ids,q_ids), zTt.in(d_ids,q_ids), zH.in(d_ids,q_ids),zRx.in(d_ids,q_ids), zCnot.in(d_ids,pairs));
    if(ibm){
        //        Qdesign.add(lambda.in(d_ids));
        if(add_rx) {
            Qdesign.add(zS.in(d_ids,q_ids), zSt.in(d_ids,q_ids), zT.in(d_ids,q_ids),
                        zTt.in(d_ids,q_ids),zRx.in(d_ids,q_ids), zCnot.in(d_ids,pairs));
        }
        else {
            Qdesign.add(zS.in(d_ids,q_ids), zSt.in(d_ids,q_ids), zT.in(d_ids,q_ids), zTt.in(d_ids,q_ids),zCnot.in(d_ids,pairs));
        }
        //
    }
    else if(cont){
        Qdesign.add(theta.in(d_ids));
        //        Qdesign.add(lambda_i.in(d_ids));
        if(add_rx){
            if (rev_cnot) {
                Qdesign.add(zS.in(d_ids,q_ids), zRx.in(d_ids,q_ids) ,zCnot.in(d_ids,pairs), zCnotR.in(d_ids,rpairs));
            }
            else {
                Qdesign.add(zS.in(d_ids,q_ids), zRx.in(d_ids,q_ids) ,zCnot.in(d_ids,pairs));
            }
        }
        else {
            if (rev_cnot) {
                Qdesign.add(zS.in(d_ids,q_ids), zCnot.in(d_ids,pairs), zCnotR.in(d_ids,rpairs));
            }
            else {
                Qdesign.add(zS.in(d_ids,q_ids), zCnot.in(d_ids,pairs));
            }
        }
    }
    else {
        Qdesign.add(zS.in(d_ids,q_ids),zT.in(d_ids,q_ids), zTt.in(d_ids,q_ids), zH.in(d_ids,q_ids), zCnot.in(d_ids,pairs), zCnotR.in(d_ids,rpairs), zSwap.in(d_ids,pairs));
    }
    Qdesign.add(slack.in(m2));
    DebugOff("Total number of variables after adding the binaries = " << Qdesign.get_nb_vars() << endl);
    /** Gate Matrices **/
    for (int depth = 0; depth < d; depth++) {
        G[depth] = var<>("G_"+to_string(depth+1), -1,1);//TODO automate bounds on G
        Qdesign.add(G[depth].in(m2));
    }
    if (d>2) {
        for (int depth = 0; depth < d-2; depth++) {
            L[depth] = var<>("L_"+to_string(depth+1), -10,10);
            Qdesign.add(L[depth].in(m2));
        }
    }
    
    DebugOff("Total number of variables after adding the G vars= " << Qdesign.get_nb_vars() << endl);
    
    
    /** Adding Constraints **/
    /** Building indexed variables as matrices (for a faster access in constraints) **/
    var<bool> zT_[d][n], zTt_[d][n], zH_[d][n], zS_[d][n], zSt_[d][n], zRx_[d][n], zCnot_[d][n][n], zCnotR_[d][n][n], zSwap_[d][n][n];
    //    var<> zT_[d][n], zTt_[d][n], zH_[d][n], zS_[d][n], zSt_[d][n], zRx_[d][n], zCnot_[d][n][n], zCnotR_[d][n][n], zSwap_[d][n][n];
    var<> G_[d][2*m][2*m], L_[d-2][2*m][2*m];
    var<> slack_[2*m][2*m];
    var<> theta_[d];
    for (int i = 0; i < 2*m; i++) {
        for (int j = 0; j < 2*m; j++) {
            slack_[i][j] = slack(i+1,j+1);
        }
    }
    for (size_t depth = 0; depth < d; depth++) {
        theta_[depth] = theta(depth+1);
        for (size_t qubit1 = 0; qubit1 < n; qubit1++) {
            zT_[depth][qubit1] = zT(depth+1,qubit1+1);
            zTt_[depth][qubit1] = zTt(depth+1,qubit1+1);
            zH_[depth][qubit1] = zH(depth+1,qubit1+1);
            zS_[depth][qubit1] = zS(depth+1,qubit1+1);
            zSt_[depth][qubit1] = zSt(depth+1,qubit1+1);
            zRx_[depth][qubit1] = zRx(depth+1,qubit1+1);
            for (size_t qubit2 = qubit1+1; qubit2 < n; qubit2++) {
                zCnot_[depth][qubit1][qubit2] = zCnot(depth+1,qubit1+1,qubit2+1);
                zCnotR_[depth][qubit1][qubit2] = zCnotR(depth+1,qubit2+1,qubit1+1);
                zSwap_[depth][qubit1][qubit2] = zSwap(depth+1,qubit1+1,qubit2+1);
                
            }
        }
        for (int i = 0; i < 2*m; i++) {
            for (int j = 0; j < 2*m; j++) {
                G_[depth][i][j] = G[depth](i+1,j+1);
            }
        }
    }
    if (d>2) {
        for (size_t depth = 0; depth < d-2; depth++) {
            for (int i = 0; i < 2*m; i++) {
                for (int j = 0; j < 2*m; j++) {
                    L_[depth][i][j] = L[depth](i+1,j+1);
                }
            }
        }
    }
    
    /** Objective **/
    //    Qdesign.min(sum(zRx) + sum(zH) + sum(zT) + sum(zTt) + sum(zCnot));
    //    Qdesign.min(sum(zRx) + sum(zS) + sum(zCnot));
    
    func_ obj;
    for (int i = 0; i < 2*m; i++) {
        for (int j = 0; j < 2*m; j++) {
            obj += slack_[i][j]*slack_[i][j];
        }
    }
    //    obj -=zSt_[0][0] + zSt_[1][1] + zCnot_[2][0][1] + zS_[3][1] + zCnot_[4][0][1];
    //    obj += sum(zCnot);
    //    func_ obj;
    //    for (size_t depth = 0; depth < d; depth++) {
    //        for (size_t qubit1 = 0; qubit1 < n; qubit1++) {
    //            obj -= zT_[depth][qubit1]*zT_[depth][qubit1] + zTt_[depth][qubit1]*zTt_[depth][qubit1] + zH_[depth][qubit1]*zH_[depth][qubit1]+ zRx_[depth][qubit1]*zRx_[depth][qubit1];
    //            for (size_t qubit2 = qubit1+1; qubit2 < n; qubit2++) {
    //                obj -= zCnot_[depth][qubit1][qubit2]*zCnot_[depth][qubit1][qubit2];
    ////                obj -= zCnot_[depth][qubit1][qubit2]*zCnot_[depth][qubit1][qubit2];
    //            }
    //        }
    //    }
    Qdesign.min(obj);
    
    
    /** Base gate matrix definition **/
    for (int depth = 0; depth < d; depth++) {
        
        for (int i = 0; i < 2*m; i++) {
            for (int j = 0; j < 2*m; j++) {
                Constraint BaseGate("BaseGate_"+to_string(depth+1)+"_"+to_string(i+1)+to_string(j+1));
                if (i%2==0) {
                    BaseGate = G_[depth][i][j];
                }
                else {
                    if (j%2==0) {//Imaginary
                        BaseGate = -1*G_[depth][i-1][j+1];
                    }
                    else {
                        BaseGate = G_[depth][i-1][j-1];
                    }
                }
                int idx = 0;
                for (int qubit = 0; qubit < n; qubit++) {
                    if(ibm){
                        //                        if ((i%2!=0 && j%2==0) || (i%2==0 && j%2!=0)) {
                        //                            BaseGate -= lambda_[depth]*zS_[depth][qubit]*S[qubit].eval(i,j) + zRx_[depth][qubit]*Rx[qubit].eval(i,j);
                        //                        }
                        //                        else {
                        //                            BaseGate -= zS_[depth][qubit]*S[qubit].eval(i,j) + zRx_[depth][qubit]*Rx[qubit].eval(i,j);
                        if(add_rx){
                            BaseGate -= zS_[depth][qubit]*S[qubit].eval(i,j)+ zSt_[depth][qubit]*St[qubit].eval(i,j) + zT_[depth][qubit]*T[qubit].eval(i,j)+ zTt_[depth][qubit]*Tt[qubit].eval(i,j) + zRx_[depth][qubit]*Rx[qubit].eval(i,j);
                        }
                        else {
                            BaseGate -= zS_[depth][qubit]*S[qubit].eval(i,j)+ zSt_[depth][qubit]*St[qubit].eval(i,j) + zT_[depth][qubit]*T[qubit].eval(i,j)+ zTt_[depth][qubit]*Tt[qubit].eval(i,j);
                        }
                        //                        BaseGate -= zS_[depth][qubit]*S[qubit].eval(i,j)+ zSt_[depth][qubit]*St[qubit].eval(i,j) + zT_[depth][qubit]*T[qubit].eval(i,j)+ zTt_[depth][qubit]*Tt[qubit].eval(i,j);
                        //                        }
                    }
                    else if (cont){
                        if ((i%2!=0 && j%2==0) || (i%2==0 && j%2!=0)) {//imaginary part
                            if(add_rx){
                                BaseGate -= S[qubit].eval(i,j)*sin(theta_[depth])*zS_[depth][qubit] + zRx_[depth][qubit]*Rx[qubit].eval(i,j);
                            }
                            else {
                                BaseGate -= S[qubit].eval(i,j)*sin(theta_[depth])*zS_[depth][qubit];
                            }
                        }
                        else {
                            if(add_rx){
                                if (i%2==0) {
                                    if(S[qubit].eval(i,j+1)!=0){
                                        BaseGate -= cos(theta_[depth])*zS_[depth][qubit] + zRx_[depth][qubit]*Rx[qubit].eval(i,j);
                                    }
                                    else {
                                        BaseGate -= zS_[depth][qubit]*S[qubit].eval(i,j) + zRx_[depth][qubit]*Rx[qubit].eval(i,j);
                                    }
                                }
                                else {
                                    if(S[qubit].eval(i,j-1)!=0){
                                        BaseGate -= cos(theta_[depth])*zS_[depth][qubit] + zRx_[depth][qubit]*Rx[qubit].eval(i,j);
                                    }
                                    else {
                                        BaseGate -= zS_[depth][qubit]*S[qubit].eval(i,j) + zRx_[depth][qubit]*Rx[qubit].eval(i,j);
                                    }
                                }
                            }
                            else {
                                if (i%2==0) {
                                    if(S[qubit].eval(i,j+1)==-1){
                                        BaseGate -= cos(theta_[depth])*zS_[depth][qubit];
                                    }
                                    else {
                                        BaseGate -= zS_[depth][qubit]*S[qubit].eval(i,j);
                                    }
                                }
                                else {
                                    if(S[qubit].eval(i,j-1)==1){
                                        BaseGate -= cos(theta_[depth])*zS_[depth][qubit];
                                    }
                                    else {
                                        BaseGate -= zS_[depth][qubit]*S[qubit].eval(i,j);
                                    }
                                }
                            }
                        }
                    }
                    else{
                        BaseGate -= zS_[depth][qubit]*S[qubit].eval(i,j) + zT_[depth][qubit]*T[qubit].eval(i,j) + zTt_[depth][qubit]*Tt[qubit].eval(i,j) + zH_[depth][qubit]*H[qubit].eval(i,j);
                    }
                    
                    for (int qubit2 = qubit+1; qubit2 < n; qubit2++) {
                        BaseGate -= zCnot_[depth][qubit][qubit2]*Cnot[idx].eval(i,j);
                        if (rev_cnot) {
                            BaseGate -= zCnotR_[depth][qubit][qubit2]*CnotR[idx].eval(i,j);
                        }
                        //                        BaseGate -= zSwap_[depth][qubit][qubit2]*Swap[idx].eval(i,j);
                        idx++;
                    }
                }
                if (BaseGate.get_nb_vars()==1) { /* This implies that the corresponding entry is zero */
                    zeros[i][j]=true;
                }
                else {
                    zeros[i][j]=false;
                    Qdesign.add(BaseGate==0);
                    //                    BaseGate.print();
                }
            }
        }
    }
    
    /** Lifted + Target gate constraint **/
    /* No need to use lifted matrices if the depth = 2 */
    if (d==2) {
        for (int i = 0; i < 2*m; i+=2) {
            for (int j = 0; j < 2*m; j++) {
                Constraint TargetGate("TargetGate_"+to_string(i+1)+to_string(j+1));
                for (int k = 0; k < 2*m; k++) {
                    if (!zeros[i][k] && !zeros[k][j]) {
                        TargetGate += G_[0][i][k]*G_[1][k][j];
                    }
                }
                TargetGate -=  slack_[i][j];
//                if (TargetGate.get_nb_vars()>1) {
                    Qdesign.add(TargetGate==M.eval(i,j));
//                }
            }
        }
    }
    else {
        /** Lifted gate constraint **/
        /** Depth 0 **/
        for (int i = 0; i < 2*m; i+=2) {
            for (int j = 0; j < 2*m; j++) {
                Constraint LiftedGate("LiftedGate_0_"+to_string(i+1)+to_string(j+1));
                LiftedGate += L_[0][i][j];
                for (int k = 0; k < 2*m; k++) {
                    if (!zeros[i][k] && !zeros[k][j]) {
                        LiftedGate -= G_[0][i][k]*G_[1][k][j];
                    }
                }
                Qdesign.add(LiftedGate==0);
            }
        }
        for (int depth = 1; depth < d-2; depth++) {
            for (int i = 0; i < 2*m; i+=2) {
                for (int j = 0; j < 2*m; j++) {
                    Constraint LiftedGate("LiftedGate"+to_string(depth)+"_"+to_string(i+1)+to_string(j+1));
                    LiftedGate += L_[depth][i][j];
                    for (int k = 0; k < 2*m; k++) {
                        if (!zeros[k][j]) {
                            LiftedGate -= L_[depth-1][i][k]*G_[depth+1][k][j];
                        }
                    }
                    Qdesign.add(LiftedGate==0);
                }
            }
        }
        
        /** Target gate constraint **/
        for (int i = 0; i < 2*m; i+=2) {
            for (int j = 0; j < 2*m; j++) {
                Constraint TargetGate("TargetGate_"+to_string(i+1)+to_string(j+1));
                for (int k = 0; k < 2*m; k++) {
                    if (!zeros[k][j]) {
                        TargetGate += L_[d-3][i][k]*G_[d-1][k][j];
                    }
                }
                TargetGate -=  slack_[i][j];
//                if (TargetGate.get_nb_vars()>1) {
                    Qdesign.add(TargetGate==M.eval(i,j));
                    //                    Qdesign.add(TargetGate==0);
//                }
                
            }
        }
    }
    
    //    for (size_t depth = 0; depth < d-1; depth++) {
    //        for (size_t qubit1 = 0; qubit1 < n; qubit1++) {
    //            Constraint NonConseqT("NonConseq_T_" +to_string(depth)+"_"+to_string(qubit1+1));
    //            NonConseqT += zT_[depth][qubit1] + zTt_[depth+1][qubit1];
    //            Qdesign.add(NonConseqT <= 1);
    //            Constraint NonConseqTt("NonConseq_Tt_" +to_string(depth)+"_"+to_string(qubit1+1));
    //            NonConseqTt += zTt_[depth][qubit1] + zTt_[depth+1][qubit1];
    //            Qdesign.add(NonConseqTt <= 1);
    //            Constraint NonConseqH("NonConseq_H_" +to_string(depth)+"_"+to_string(qubit1+1));
    //            NonConseqH += zH_[depth][qubit1] + zH_[depth+1][qubit1];
    //            Qdesign.add(NonConseqH <= 1);
    //            for (size_t qubit2 = qubit1+1; qubit2 < n; qubit2++) {
    //                Constraint NonConseqCnot("NonConseq_Cnot_" +to_string(depth)+"_"+to_string(qubit1+1)+"_"+to_string(qubit2+1));
    //                NonConseqCnot += zCnot_[depth][qubit1][qubit2] + zCnot_[depth+1][qubit1][qubit2];
    //                Qdesign.add(NonConseqCnot <= 1);
    ////                Constraint Symmetry1("Symmetry1" +to_string(depth)+"_"+to_string(qubit1+1)+"_"+to_string(qubit2+1));
    ////                Symmetry1 += zH_[depth][qubit1] - zH_[depth+1][qubit2];
    ////                Qdesign.add(Symmetry1 >= 0);
    ////                Constraint Symmetry2("Symmetry2" +to_string(depth)+"_"+to_string(qubit1+1)+"_"+to_string(qubit2+1));
    ////                Symmetry2 += zT_[depth][qubit1] - zT_[depth+1][qubit2];
    ////                Qdesign.add(Symmetry2 >= 0);
    ////                Constraint Symmetry3("Symmetry3" +to_string(depth)+"_"+to_string(qubit1+1)+"_"+to_string(qubit2+1));
    ////                Symmetry3 += zTt_[depth][qubit1] - zTt_[depth+1][qubit2];
    ////                Qdesign.add(Symmetry3 >= 0);
    ////                Constraint Symmetry4("Symmetry4" +to_string(depth)+"_"+to_string(qubit1+1)+"_"+to_string(qubit2+1));
    ////                Symmetry4 += zTt_[depth][qubit1] - zH_[depth+1][qubit2];
    ////                Qdesign.add(Symmetry4 >= 0);
    ////                for (size_t qubit3 = 0; qubit3 < n; qubit3++) {
    ////                    if (qubit3!=qubit1 && qubit3!=qubit2) {
    ////                        Constraint Symmetry5("Symmetry5" +to_string(depth)+"_"+to_string(qubit1+1)+"_"+to_string(qubit2+1)+"_"+to_string(qubit3+1));
    ////                        Symmetry5 += zH_[depth][qubit3] - zCnot_[depth][qubit1][qubit2];
    ////                        Qdesign.add(Symmetry5 >= 0);
    ////                    }
    ////                    if (qubit3<qubit2) {
    ////                        Constraint Symmetry6("Symmetry6" +to_string(depth)+"_"+to_string(qubit1+1)+"_"+to_string(qubit2+1)+"_"+to_string(qubit3+1));
    ////                        Symmetry6 += zT_[depth][qubit3] - zCnot_[depth+1][qubit1][qubit2];
    ////                        Qdesign.add(Symmetry6 <= 0);
    ////                        Constraint Symmetry7("Symmetry7" +to_string(depth)+"_"+to_string(qubit1+1)+"_"+to_string(qubit2+1)+"_"+to_string(qubit3+1));
    ////                        Symmetry7 += zTt_[depth][qubit3] - zCnot_[depth][qubit1][qubit2];
    ////                        Qdesign.add(Symmetry7 >= 0);
    ////                    }
    ////                }
    //
    //            }
    //        }
    //    }
    //    Constraint TargetGate("TargetGate");
    //    TargetGate += G[d-2]*G[d-1] - M;
    //    Qdesign.add(TargetGate==0);
    //    TargetGate.print();
    
    //    Constraint TOFF_H("TOFF_H");
    //    TOFF_H += sum(zCnot);
    //    Qdesign.add(TOFF_H>=1);
    ////
    //    Constraint TOFF_T("TOFF_T");
    //    TOFF_T += sum(zT)+sum(zTt);
    //    Qdesign.add(TOFF_T==0);
    
    //    Constraint Config1("Config1");
    //    Config1 += zCnotR_[2][0][1];
    //    Qdesign.add(Config1==1);
    //
    //    Constraint Config2("Config2");
    //    Config2 += zCnot_[5][0][1];
    //    Qdesign.add(Config2==1);
    //
    //    Constraint Config3("Config3");
    //    Config3 += zCnotR_[7][0][1];
    //    Qdesign.add(Config3==1);
    
    
    //
    //    Constraint TOFF_C1("TOFF_C1");
    //    TOFF_C1 += zRx_[1][1];
    //    Qdesign.add(TOFF_C1==1);
    //////
    //    Constraint TOFF_C2("TOFF_C2");
    //    TOFF_C2 += zRx_[5][1];
    //    Qdesign.add(TOFF_C2==1);
    ////
    //    Constraint TOFF_C3("TOFF_C3");
    //    TOFF_C3 += zS_[0][1];
    //    Qdesign.add(TOFF_C3==1);
    //////
    //    Constraint TOFF_C4("TOFF_C4");
    //    TOFF_C4 += zS_[2][1];
    //    Qdesign.add(TOFF_C4==1);
    ////
        Constraint TOFF_C5("TOFF_C5");
        TOFF_C5 += theta_[2];
        Qdesign.add(TOFF_C5==-pi/2);
//    //    //
        Constraint TOFF_C6("TOFF_C6");
        TOFF_C6 += theta_[0];
        Qdesign.add(TOFF_C6==pi/2);
//    ////
//        Constraint TOFF_C7("TOFF_C7");
//        TOFF_C7 += sum(zRx);
//        Qdesign.add(TOFF_C7==0);
    //
//        Constraint TOFF_C8("TOFF_C8");
//        TOFF_C8 += zCnot_[1][0][1];
//        Qdesign.add(TOFF_C8==1);
//    
    //
    
    //    Constraint TOFF_Tt("TOFF_Tt");
    //    TOFF_Tt += sum(zCnot);
    //    Qdesign.add(TOFF_Tt>=6);
    
    /** One gate per depth constraints **/
    for(int depth = 0; depth < d ; depth++){
        Constraint OneGate("OneGate_"+to_string(depth+1));
        for(int qubit = 0; qubit < n ; qubit++){
            if(ibm){
                //                OneGate += zS_[depth][qubit] + zRx_[depth][qubit];
                if(add_rx){
                    OneGate += zS_[depth][qubit] + zSt_[depth][qubit] + zT_[depth][qubit] + zTt_[depth][qubit] + zRx_[depth][qubit];
                }
                else {
                    OneGate += zS_[depth][qubit] + zSt_[depth][qubit] + zT_[depth][qubit] + zTt_[depth][qubit];
                }
                //                OneGate += zS_[depth][qubit] + zSt_[depth][qubit] + zT_[depth][qubit] + zTt_[depth][qubit];
            }
            else if (cont){
                if(add_rx){
                    OneGate += zS_[depth][qubit] + zRx_[depth][qubit];
                }
                else {
                    OneGate += zS_[depth][qubit];
                }
            }
            else {
                OneGate += zS_[depth][qubit] + zT_[depth][qubit] + zTt_[depth][qubit] + zH_[depth][qubit];
            }
            for (int qubit2 = qubit+1; qubit2 < n; qubit2++) {
                OneGate += zCnot_[depth][qubit][qubit2];
                if (rev_cnot) {
                    OneGate += zCnotR_[depth][qubit][qubit2];
                }
                //                OneGate += zSwap_[depth][qubit][qubit2];
            }
        }
        Qdesign.add(OneGate==1);
    }
    //    Qdesign.print_expanded();
    auto solver_time_start = get_wall_time();
    if(relax){
        auto Relaxation = Qdesign.build_McCormick();
        Relaxation->print_expanded();
        solver slvr(*Relaxation,cplex);
        slvr.run(output, relax = false, tol=1e-6, 0.01, "ma27");
        Qdesign.print_solution(false);
    }
    else {
        solver slvr(Qdesign,bonmin);
        slvr.run(output, relax = false, tol=1e-6, 0.01, "ma27");
        Qdesign.print_solution(false);
    }
    auto solver_time_end = get_wall_time();
    auto solve_time = solver_time_end - solver_time_start;
    
    auto total_time_start = get_wall_time();
    auto total_time_end = get_wall_time();
    auto total_time = total_time_end - total_time_start;
    DebugOn("Solver Time = " << to_string(solve_time) << " seconds" << endl);
    DebugOn("Total Computing Time = " << to_string(total_time) << " seconds" << endl);
#else
    cerr << "Error: this version of Gravity "
    "was compiled without QPP support. Please rerun cmake with -DENABLE_QPP=1." << endl;
#endif
}

void run_qtdesign(const int n, const int d, const param<>& M, bool relax){
#ifdef USE_QPP
    double tol = 1e-6;
    bool ibm = false;
    bool cont = true;
    bool add_rx = true;
    bool rev_cnot = false;
    int output = 0;
    if (ibm && cont) {
        cerr << "Choose IBM or Continuous, or none.";
        exit(-1);
    }

    /** Parameters **/
    
    /* Qubit number */
    const int m = pow(2,n);
    DebugOn("Number of Qubits = " << n << endl);
    /* Circuit depth */
    DebugOn("Depth = " << d << endl);
    if (d<2) {
        cerr << "Depth should be greater or equal than 2!" << endl;
        return;
    }
    /* T and T conjugate transpose gate matrices, one for each row (qubit) */
    param<> T[n], Tt[n];
    
    /* S and St conjugate transpose gate matrices, one for each row (qubit) */
    param<> S[n], St[n];
    
    /* Hadamard gate matrices */
    param<> H[n];
    
    /* Rx gate matrices */
    param<> Rx[n];
    
    /* Cnot gate matrices */
    param<> Cnot[n*(n-1)/2];
    param<> CnotR[n*(n-1)/2];
    /* SWAP gate matrices */
    param<> Swap[n*(n-1)/2];

    
    
    
    /** Initialization **/
    int idx = 0;
    for (auto i = 0; i<n; i++) { // Qubit (row) iterator
        T[i].set_name("T_"+to_string(i+1));
        T[i].set_size(2*m, 2*m, 0);/* Representing Complex number as 2x2 real matrices */
        T[i].QuantumT(i,n);/* Fill nonzero values corresponding to a T gate */
        Tt[i].set_name("Tt_"+to_string(i+1));
        Tt[i].set_size(2*m, 2*m, 0);
        Tt[i].QuantumT(i,n,true);/* Fill nonzero values corresponding to a T conjugate transpose gate */
        Rx[i].set_name("Rx_"+to_string(i+1));
        Rx[i].set_size(2*m, 2*m, 0);
        Rx[i].QuantumRx(i,n);/* Fill nonzero values corresponding to an H gate */
        S[i].set_name("S_"+to_string(i+1));
        S[i].set_size(2*m, 2*m, 0);
        S[i].QuantumS(i,n);/* Fill nonzero values corresponding to an H gate */
        St[i].set_name("St_"+to_string(i+1));
        St[i].set_size(2*m, 2*m, 0);
        St[i].QuantumSt(i,n);/* Fill nonzero values corresponding to an H gate */
        H[i].set_name("H_"+to_string(i+1));
        H[i].set_size(2*m, 2*m, 0);
        H[i].QuantumH(i,n);/* Fill nonzero values corresponding to an H gate */
        for (auto k = i+1; k<n; k++) {
            Cnot[idx].set_name("Cnot_"+to_string(i+1)+to_string(k+1));
            Cnot[idx].set_size(2*m, 2*m, 0);
            Cnot[idx].QuantumCnot(i,k,n);/* Fill nonzero values corresponding to a CNOT gate */
            CnotR[idx].set_name("CnotR_"+to_string(k+1)+to_string(i+1));
            CnotR[idx].set_size(2*m, 2*m, 0);
            CnotR[idx].QuantumCnot(k,i,n);/* Fill nonzero values corresponding to a CNOT gate */
            Swap[idx].set_name("Swap"+to_string(i+1)+to_string(k+1));
            Swap[idx].set_size(2*m, 2*m, 0);
            Swap[idx].QuantumSwap(i,k,n);/* Fill nonzero values corresponding to a CNOT gate */

            idx++;
        }
    }
    
    
    /** Indices **/
    auto d_ids = indices(1,d);/*< Depth indices */
    auto q_ids = indices(1,n);/*< Qubit indices */
    auto pairs = indices(ordered_pairs(2, n+1)); /*< Indices pairs (control-taget) for Cnot gates */
    auto rpairs = indices(ordered_pairs(2, n+1, true)); /*< Indices pairs (control-taget) for Cnot gates */
    auto m2 = indices(indices(1,2*m), indices(1,2*m)); /*< Indices for square mxm matrices */
    
    
    /** Variables **/
    /** Binaries **/
    var<bool> zT("zT"); /*< Binary variables for T gates */
    var<bool> zTt("zTt"); /*< Binary variables for T conjugate gates */
    var<bool> zH("zH"); /*< Binary variables for H gates */
    var<bool> zS("zS"); /*< Binary variables for S gates */
    var<bool> zSt("zSt"); /*< Binary variables for S gates */
    var<bool> zRx("zRx"); /*< Binary variables for Rx gates */
    var<bool> zCnot("zCnot"); /*< Binary variables for CNOT gates */
    var<bool> zCnotR("zCnotR"); /*< Binary variables for CNOT gates */
    var<bool> zSwap("zSwap"); /*< Binary variables for SWAP gates */
    
//    var<> zT("zT", 0,1); /*< Binary variables for T gates */
//    var<> zTt("zTt", 0,1); /*< Binary variables for T conjugate gates */
//    var<> zH("zH", 0,1); /*< Binary variables for H gates */
//    var<> zS("zS", 0,1); /*< Binary variables for S gates */
//    var<> zSt("zSt", 0,1); /*< Binary variables for S gates */
//    var<> zRx("zRx", 0,1); /*< Binary variables for Rx gates */
//    var<> zCnot("zCnot", 0,1); /*< Binary variables for CNOT gates */
//    var<> zCnotR("zCnotR", 0,1); /*< Binary variables for CNOT gates */
//    var<> zSwap("zSwap", 0,1); /*< Binary variables for SWAP gates */

    /** Gate Matrices **/
    var<> G[d]; /*< Base Gate Matrices */
    param<> Gl[d], Gu[d]; /*< Lower and Upper Bounds on Base Gate Matrices */
    bool zeros[2*m][2*m]; /*< Zeros in Gate Matrices */
    var<> L[d-2]; /*< Lifted Product Matrices */
    param<> Ll[d], Lu[d]; /*< Lower and Upper Bounds on Lifted Matrices */
    var<> lambda_r("lambda_r", -1, 1);
//    lambda_r.initialize_all(0.1);
//    lambda_r.initialize_uniform();
    var<> lambda_i("lambda_i", -1, 1);
    var<> slack("slack",-1,1);
    
    
    
    /** Model **/
    
    /** Adding Variables **/
    Model Qdesign("Qdesign");
    /** Binaries **/

    //    Qdesign.add(zT.in(d_ids,q_ids), zTt.in(d_ids,q_ids), zH.in(d_ids,q_ids), zCnot.in(d_ids,pairs));
//    Qdesign.add(zT.in(d_ids,q_ids), zTt.in(d_ids,q_ids), zH.in(d_ids,q_ids),zRx.in(d_ids,q_ids), zCnot.in(d_ids,pairs));
    if(ibm){
//        Qdesign.add(lambda.in(d_ids));
        if(add_rx) {
                Qdesign.add(zS.in(d_ids,q_ids), zSt.in(d_ids,q_ids), zT.in(d_ids,q_ids),
                        zTt.in(d_ids,q_ids),zRx.in(d_ids,q_ids), zCnot.in(d_ids,pairs));
        }
        else {
            Qdesign.add(zS.in(d_ids,q_ids), zSt.in(d_ids,q_ids), zT.in(d_ids,q_ids), zTt.in(d_ids,q_ids),zCnot.in(d_ids,pairs));
        }
//
    }
    else if(cont){
//        Qdesign.add(lambda_r.in(d_ids,q_ids), lambda_i.in(d_ids,q_ids));
        Qdesign.add(lambda_r.in(d_ids), lambda_i.in(d_ids));
//        Qdesign.add(lambda_i.in(d_ids));
        if(add_rx){
            if (rev_cnot) {
                Qdesign.add(zS.in(d_ids,q_ids), zRx.in(d_ids,q_ids) ,zCnot.in(d_ids,pairs), zCnotR.in(d_ids,rpairs));
            }
            else {
                Qdesign.add(zS.in(d_ids,q_ids), zRx.in(d_ids,q_ids) ,zCnot.in(d_ids,pairs));
            }
        }
        else {
            if (rev_cnot) {
                Qdesign.add(zS.in(d_ids,q_ids), zCnot.in(d_ids,pairs), zCnotR.in(d_ids,rpairs));
            }
            else {
                Qdesign.add(zS.in(d_ids,q_ids), zCnot.in(d_ids,pairs));
            }
        }
    }
    else {
        Qdesign.add(zT.in(d_ids,q_ids), zTt.in(d_ids,q_ids), zH.in(d_ids,q_ids), zCnot.in(d_ids,pairs), zCnotR.in(d_ids,rpairs));
    }
    Qdesign.add(slack.in(m2));
    DebugOff("Total number of variables after adding the binaries = " << Qdesign.get_nb_vars() << endl);
    /** Gate Matrices **/
    for (int depth = 0; depth < d; depth++) {
        Gl[depth] = param<>("Gl_"+to_string(depth+1));
        Gu[depth] = param<>("Gu_"+to_string(depth+1));
        if (relax) {
            G[depth] = var<>("G_"+to_string(depth+1), Gl[depth].in(m2), Gu[depth].in(m2));
        }
        else{
            G[depth] = var<>("G_"+to_string(depth+1), -1, 1);
        }
        
        Qdesign.add(G[depth].in(m2));
    }
    if (d>2) {
        for (int depth = 0; depth < d-2; depth++) {
            Ll[depth] = param<>("Ll_"+to_string(depth+1));
            Lu[depth] = param<>("Lu_"+to_string(depth+1));
            if (relax) {
                L[depth] = var<>("L_"+to_string(depth+1), Ll[depth].in(m2), Lu[depth].in(m2));
            }
            else {
                L[depth] = var<>("L_"+to_string(depth+1), -1, 1);
            }
            Qdesign.add(L[depth].in(m2));
        }
    }

    DebugOff("Total number of variables after adding the G vars= " << Qdesign.get_nb_vars() << endl);
    
    
    /** Adding Constraints **/
    /** Building indexed variables as matrices (for a faster access in constraints) **/
    var<bool> zT_[d][n], zTt_[d][n], zH_[d][n], zS_[d][n], zSt_[d][n], zRx_[d][n], zCnot_[d][n][n], zCnotR_[d][n][n], zSwap_[d][n][n];
//    var<> zT_[d][n], zTt_[d][n], zH_[d][n], zS_[d][n], zSt_[d][n], zRx_[d][n], zCnot_[d][n][n], zCnotR_[d][n][n], zSwap_[d][n][n];
    var<> G_[d][2*m][2*m], L_[d-2][2*m][2*m];
    param<> Gl_[d][2*m][2*m], Gu_[d][2*m][2*m], Ll_[d-2][2*m][2*m], Lu_[d-2][2*m][2*m];
    var<> slack_[2*m][2*m];
    var<> lambda_r_[d], lambda_i_[d];
//    var<> lambda_r_[d][n], lambda_i_[d][n];
    for (int i = 0; i < 2*m; i++) {
        for (int j = 0; j < 2*m; j++) {
            slack_[i][j] = slack(i+1,j+1);
        }
    }
    for (size_t depth = 0; depth < d; depth++) {
        lambda_r_[depth] = lambda_r(depth+1);
        lambda_i_[depth] = lambda_i(depth+1);
        for (size_t qubit1 = 0; qubit1 < n; qubit1++) {
//            lambda_r_[depth][qubit1] = lambda_r(depth+1,qubit1+1);
//            lambda_i_[depth][qubit1] = lambda_i(depth+1,qubit1+1);
            zT_[depth][qubit1] = zT(depth+1,qubit1+1);
            zTt_[depth][qubit1] = zTt(depth+1,qubit1+1);
            zH_[depth][qubit1] = zH(depth+1,qubit1+1);
            zS_[depth][qubit1] = zS(depth+1,qubit1+1);
            zSt_[depth][qubit1] = zSt(depth+1,qubit1+1);
            zRx_[depth][qubit1] = zRx(depth+1,qubit1+1);
            for (size_t qubit2 = qubit1+1; qubit2 < n; qubit2++) {
                zCnot_[depth][qubit1][qubit2] = zCnot(depth+1,qubit1+1,qubit2+1);
                zCnotR_[depth][qubit1][qubit2] = zCnotR(depth+1,qubit2+1,qubit1+1);
                zSwap_[depth][qubit1][qubit2] = zSwap(depth+1,qubit1+1,qubit2+1);

            }
        }
        for (int i = 0; i < 2*m; i++) {
            for (int j = 0; j < 2*m; j++) {
                G_[depth][i][j] = G[depth](i+1,j+1);
                Gl_[depth][i][j] = Gl[depth](i+1,j+1);
                Gu_[depth][i][j] = Gu[depth](i+1,j+1);
            }
        }
    }
    if (d>2) {
        for (size_t depth = 0; depth < d-2; depth++) {
            for (int i = 0; i < 2*m; i++) {
                for (int j = 0; j < 2*m; j++) {
                    L_[depth][i][j] = L[depth](i+1,j+1);
                    Ll_[depth][i][j] = Ll[depth](i+1,j+1);
                    Lu_[depth][i][j] = Lu[depth](i+1,j+1);
                }
            }
        }
    }
    
    /** Objective **/
//    Qdesign.min(sum(zRx) + sum(zH) + sum(zT) + sum(zTt) + sum(zCnot));
//    Qdesign.min(sum(zRx) + sum(zS) + sum(zCnot));
    
    func_ obj;
    if(!relax){
        for (int i = 0; i < 2*m; i++) {
            for (int j = 0; j < 2*m; j++) {
                obj += slack_[i][j]*slack_[i][j];
            }
        }
    }
//    obj -=zSt_[0][0] + zSt_[1][1] + zCnot_[2][0][1] + zS_[3][1] + zCnot_[4][0][1];
    if(relax){
        obj += sum(zCnot);
    }
//    func_ obj;
//    for (size_t depth = 0; depth < d; depth++) {
//        for (size_t qubit1 = 0; qubit1 < n; qubit1++) {
//            obj -= zT_[depth][qubit1]*zT_[depth][qubit1] + zTt_[depth][qubit1]*zTt_[depth][qubit1] + zH_[depth][qubit1]*zH_[depth][qubit1]+ zRx_[depth][qubit1]*zRx_[depth][qubit1];
//            for (size_t qubit2 = qubit1+1; qubit2 < n; qubit2++) {
//                obj -= zCnot_[depth][qubit1][qubit2]*zCnot_[depth][qubit1][qubit2];
////                obj -= zCnot_[depth][qubit1][qubit2]*zCnot_[depth][qubit1][qubit2];
//            }
//        }
//    }
    Qdesign.min(obj);
//    Qdesign.min(1);

    
    /** Base gate matrix definition **/
    for (int depth = 0; depth < d; depth++) {
//        Constraint Complex("Complex"+to_string(depth+1));
//        Complex += lambda_i_[depth]*lambda_i_[depth] + lambda_r_[depth]*lambda_r_[depth];
//        Qdesign.add(Complex>=1);
        if(cont){
            for (int qubit = 0; qubit < n; qubit++) {
                Constraint Complex("Complex"+to_string(depth+1)+to_string(qubit+1));
//                Complex += lambda_i_[depth][qubit]*lambda_i_[depth][qubit] + lambda_r_[depth][qubit]*lambda_r_[depth][qubit];
                Complex += lambda_i_[depth]*lambda_i_[depth] + lambda_r_[depth]*lambda_r_[depth];
                if(relax){
                    Qdesign.add(Complex==1);
                }
                else{
                    Qdesign.add(Complex>=1);
                }
            }
        }
        
        
        for (int i = 0; i < 2*m; i++) {
            for (int j = 0; j < 2*m; j++) {
                double lb = 1, ub = -1;
                Constraint BaseGate("BaseGate_"+to_string(depth+1)+"_"+to_string(i+1)+to_string(j+1));
                if (i%2==0) {
                    BaseGate = G_[depth][i][j];
                }
                else {
                    if (j%2==0) {//Imaginary
                        BaseGate = -1*G_[depth][i-1][j+1];
                    }
                    else {
                        BaseGate = G_[depth][i-1][j-1];
                    }
                }
                int idx = 0;
                for (int qubit = 0; qubit < n; qubit++) {
                    if(ibm){
//                        if ((i%2!=0 && j%2==0) || (i%2==0 && j%2!=0)) {
//                            BaseGate -= lambda_[depth]*zS_[depth][qubit]*S[qubit].eval(i,j) + zRx_[depth][qubit]*Rx[qubit].eval(i,j);
//                        }
//                        else {
//                            BaseGate -= zS_[depth][qubit]*S[qubit].eval(i,j) + zRx_[depth][qubit]*Rx[qubit].eval(i,j);
                        if(add_rx){
                            BaseGate -= zS_[depth][qubit]*S[qubit].eval(i,j)+ zSt_[depth][qubit]*St[qubit].eval(i,j) + zT_[depth][qubit]*T[qubit].eval(i,j)+ zTt_[depth][qubit]*Tt[qubit].eval(i,j) + zRx_[depth][qubit]*Rx[qubit].eval(i,j);
                        }
                        else {
                            BaseGate -= zS_[depth][qubit]*S[qubit].eval(i,j)+ zSt_[depth][qubit]*St[qubit].eval(i,j) + zT_[depth][qubit]*T[qubit].eval(i,j)+ zTt_[depth][qubit]*Tt[qubit].eval(i,j);
                        }
//                        BaseGate -= zS_[depth][qubit]*S[qubit].eval(i,j)+ zSt_[depth][qubit]*St[qubit].eval(i,j) + zT_[depth][qubit]*T[qubit].eval(i,j)+ zTt_[depth][qubit]*Tt[qubit].eval(i,j);
//                        }
                    }
                    else if (cont){
                        if ((i%2!=0 && j%2==0) || (i%2==0 && j%2!=0)) {//imaginary part
                            if(add_rx){
//                                BaseGate -= S[qubit].eval(i,j)*lambda_i_[depth][qubit]*zS_[depth][qubit] + zRx_[depth][qubit]*Rx[qubit].eval(i,j);
                                BaseGate -= S[qubit].eval(i,j)*lambda_i_[depth]*zS_[depth][qubit] + zRx_[depth][qubit]*Rx[qubit].eval(i,j);
                                if(i%2==0){
                                    lb = min(lb,S[qubit].eval(i,j)*-1);
                                    lb = min(lb,Rx[qubit].eval(i,j));
                                    ub = max(ub,S[qubit].eval(i,j));
                                    ub = max(ub,Rx[qubit].eval(i,j));
                                }
                            }
                            else {
//                                BaseGate -= S[qubit].eval(i,j)*lambda_i_[depth][qubit]*zS_[depth][qubit];
                                BaseGate -= S[qubit].eval(i,j)*lambda_i_[depth]*zS_[depth][qubit];
                                if(i%2==0){
                                    lb = min(lb,S[qubit].eval(i,j)*-1);
                                    ub = max(ub,S[qubit].eval(i,j));
                                }
                            }
                        }
                        else {
                            if(add_rx){
                                if (i%2==0) {
                                    if(S[qubit].eval(i,j+1)!=0){
//                                       BaseGate -= lambda_r_[depth][qubit]*zS_[depth][qubit] + zRx_[depth][qubit]*Rx[qubit].eval(i,j);
                                        BaseGate -= lambda_r_[depth]*zS_[depth][qubit] + zRx_[depth][qubit]*Rx[qubit].eval(i,j);
                                        lb = min(lb,-1.);
                                        lb = min(lb,Rx[qubit].eval(i,j));
                                        ub = max(ub,1.);
                                        ub = max(ub,Rx[qubit].eval(i,j));
                                    }
                                    else {
                                        BaseGate -= zS_[depth][qubit]*S[qubit].eval(i,j) + zRx_[depth][qubit]*Rx[qubit].eval(i,j);
                                        lb = min(lb,S[qubit].eval(i,j));
                                        lb = min(lb,Rx[qubit].eval(i,j));
                                        ub = max(ub,S[qubit].eval(i,j));
                                        ub = max(ub,Rx[qubit].eval(i,j));
                                    }
                                    
                                }
                                else {
                                    if(S[qubit].eval(i,j-1)!=0){
//                                        BaseGate -= lambda_r_[depth][qubit]*zS_[depth][qubit] + zRx_[depth][qubit]*Rx[qubit].eval(i,j);
                                        BaseGate -= lambda_r_[depth]*zS_[depth][qubit] + zRx_[depth][qubit]*Rx[qubit].eval(i,j);
                                    }
                                    else {
                                        BaseGate -= zS_[depth][qubit]*S[qubit].eval(i,j) + zRx_[depth][qubit]*Rx[qubit].eval(i,j);
                                    }
                                }
                            }
                            else {
                                if (i%2==0) {
                                    if(S[qubit].eval(i,j+1)==-1){
//                                        BaseGate -= lambda_r_[depth][qubit]*zS_[depth][qubit];
                                        BaseGate -= lambda_r_[depth]*zS_[depth][qubit];
                                        lb = min(lb,-1.);
                                        ub = max(ub,1.);
                                    }
                                    else {
                                        BaseGate -= zS_[depth][qubit]*S[qubit].eval(i,j);
                                        lb = min(lb,S[qubit].eval(i,j));
                                        ub = max(ub,S[qubit].eval(i,j));
                                    }
                                }
                                else {
                                    if(S[qubit].eval(i,j-1)==1){
//                                        BaseGate -= lambda_r_[depth][qubit]*zS_[depth][qubit];
                                        BaseGate -= lambda_r_[depth]*zS_[depth][qubit];
                                    }
                                    else {
                                        BaseGate -= zS_[depth][qubit]*S[qubit].eval(i,j);
                                    }
                                }
                            }
                        }
                    }
                    else{
                        BaseGate -= zT_[depth][qubit]*T[qubit].eval(i,j) + zTt_[depth][qubit]*Tt[qubit].eval(i,j) + zH_[depth][qubit]*H[qubit].eval(i,j);
                        lb = min(lb,T[qubit].eval(i,j));
                        lb = min(lb,Tt[qubit].eval(i,j));
                        lb = min(lb,H[qubit].eval(i,j));
                        ub = max(ub,T[qubit].eval(i,j));
                        ub = max(ub,Tt[qubit].eval(i,j));
                        ub = max(ub,H[qubit].eval(i,j));
                    }
                    
                    for (int qubit2 = qubit+1; qubit2 < n; qubit2++) {
                        BaseGate -= zCnot_[depth][qubit][qubit2]*Cnot[idx].eval(i,j);
                        lb = min(lb, Cnot[idx].eval(i,j));
                        ub = max(ub, Cnot[idx].eval(i,j));
                        if (rev_cnot) {
                            BaseGate -= zCnotR_[depth][qubit][qubit2]*CnotR[idx].eval(i,j);
                            lb = min(lb, CnotR[idx].eval(i,j));
                            ub = max(ub, CnotR[idx].eval(i,j));
                        }
//                        BaseGate -= zSwap_[depth][qubit][qubit2]*Swap[idx].eval(i,j);
                        idx++;
                    }
                }
                if (BaseGate.get_nb_vars()==1) { /* This implies that the corresponding entry is zero */
                    zeros[i][j]=true;
//                    if(i%2==0){
//                        Gl_[depth][i][j].set_val(0);
//                        Gu_[depth][i][j].set_val(0);
//                    }
                }
                else {
                    zeros[i][j]=false;
                    if(i%2==0){
                        Gl_[depth][i][j] = lb;
                        Gu_[depth][i][j] = ub;
                        if (j%2==0) {
                            Gl_[depth][i+1][j+1] = lb;
                            Gu_[depth][i+1][j+1] = ub;
                        }
                        else {
                            Gl_[depth][i+1][j-1] = -ub;
                            Gu_[depth][i+1][j-1] = -lb;
                        }
                        Qdesign.add(BaseGate==0);
                    }
//                    BaseGate.print();
                }
            }
        }
    }
    
    /** Lifted + Target gate constraint **/
    /* No need to use lifted matrices if the depth = 2 */
    func_ G1, G2;
    if (d==2) {
        for (int i = 0; i < 2*m; i+=2) {
            for (int j = 0; j < 2*m; j++) {
                Constraint TargetGate("TargetGate_"+to_string(i+1)+to_string(j+1));
                for (int k = 0; k < 2*m; k++) {
                    if (!zeros[i][k] && !zeros[k][j]) {
                        if (k%2==0) {
                            G2 = G_[1][k][j];
                        }
                        else {
                            if (j%2==0) {//Imaginary
                                G2 = -1*G_[1][k-1][j+1];
                            }
                            else {
                                G2 = G_[1][k-1][j-1];
                            }
                        }
                        TargetGate += G_[0][i][k]*G2;
                    }
                }
                TargetGate -=  slack_[i][j];
//                if (TargetGate.get_nb_vars()>1) {
                    Qdesign.add(TargetGate==M.eval(i,j));
//                }
            }
        }
    }
    else {
        /** Lifted gate constraint **/
        /** Depth 0 **/
        for (int i = 0; i < 2*m; i+=2) {
            for (int j = 0; j < 2*m; j++) {
                double lb = 1, ub = -1;
                Constraint LiftedGate("LiftedGate_0_"+to_string(i+1)+to_string(j+1));
                LiftedGate += L_[0][i][j];
                for (int k = 0; k < 2*m; k++) {
                    if (!zeros[i][k] && !zeros[k][j]) {
                        lb = min(lb, (Gl_[0][i][k].eval()*Gl_[1][k][j].eval()));
                        lb = min(lb, (Gl_[0][i][k].eval()*Gu_[1][k][j].eval()));
                        lb = min(lb, (Gu_[0][i][k].eval()*Gl_[1][k][j].eval()));
                        lb = min(lb, (Gu_[0][i][k].eval()*Gu_[1][k][j].eval()));
                        ub = max(ub, (Gl_[0][i][k].eval()*Gl_[1][k][j].eval()));
                        ub = max(ub, (Gl_[0][i][k].eval()*Gu_[1][k][j].eval()));
                        ub = max(ub, (Gu_[0][i][k].eval()*Gl_[1][k][j].eval()));
                        ub = max(ub, (Gu_[0][i][k].eval()*Gu_[1][k][j].eval()));
                        if (k%2==0) {
                            G2 = G_[1][k][j];
                        }
                        else {
                            if (j%2==0) {//Imaginary
                                G2 = -1*G_[1][k-1][j+1];
                            }
                            else {
                                G2 = G_[1][k-1][j-1];
                            }
                        }
                        LiftedGate -= G_[0][i][k]*G2;
                        Ll_[0][i][j] = Ll_[0][i][j].eval() + lb;
                        Lu_[0][i][j] = Lu_[0][i][j].eval() + ub;
                        if(Ll_[0][i][j].eval()<-1){
                            Ll_[0][i][j] = -1;
                        }
                        if(Lu_[0][i][j].eval()>1){
                            Lu_[0][i][j] = 1;
                        }
                        if (j%2==0) {//Real part
                            Ll_[0][i+1][j+1] = Ll_[0][i][j].eval();
                            Lu_[0][i+1][j+1] = Lu_[0][i][j].eval();
                        }
                        else{
                            Ll_[0][i+1][j-1] = -Lu_[0][i][j].eval();
                            Lu_[0][i+1][j-1] = -Ll_[0][i][j].eval();
                        }
                    }
//                    else {
//                        Ll_[0][i][j].set_val(0);
//                        Lu_[0][i][j].set_val(0);
//
//                        Ll[0].set_val(i,j,0);
//                        Lu[0].set_val(i,j,0);
//                        if (j%2==0) {
//                            Ll[0].set_val(i+1,j+1,0);
//                            Lu[0].set_val(i+1,j+1,0);
//                        }
//                        else{
//                            Ll[0].set_val(i+1,j-1,0);
//                            Lu[0].set_val(i+1,j-1,0);
//                        }
//                    }
                }
                Qdesign.add(LiftedGate==0);
            }
        }
        for (int depth = 1; depth < d-2; depth++) {
            for (int i = 0; i < 2*m; i+=2) {
                for (int j = 0; j < 2*m; j++) {
                    double lb = 1, ub = -1;
                    Constraint LiftedGate("LiftedGate"+to_string(depth)+"_"+to_string(i+1)+to_string(j+1));
                    LiftedGate += L_[depth][i][j];
                    for (int k = 0; k < 2*m; k++) {
                        if (!zeros[k][j]) {
                            lb = min(lb, (Ll_[depth-1][i][k].eval()*Gl_[depth+1][k][j].eval()));
                            lb = min(lb, (Ll_[depth-1][i][k].eval()*Gu_[depth+1][k][j].eval()));
                            lb = min(lb, (Lu_[depth-1][i][k].eval()*Gl_[depth+1][k][j].eval()));
                            lb = min(lb, (Lu_[depth-1][i][k].eval()*Gu_[depth+1][k][j].eval()));
                            ub = max(ub, (Ll_[depth-1][i][k].eval()*Gl_[depth+1][k][j].eval()));
                            ub = max(ub, (Ll_[depth-1][i][k].eval()*Gu_[depth+1][k][j].eval()));
                            ub = max(ub, (Lu_[depth-1][i][k].eval()*Gl_[depth+1][k][j].eval()));
                            ub = max(ub, (Lu_[depth-1][i][k].eval()*Gu_[depth+1][k][j].eval()));
                            if (k%2==0) {
                                G2 = G_[depth+1][k][j];
                            }
                            else {
                                if (j%2==0) {//Imaginary
                                    G2 = -1*G_[depth+1][k-1][j+1];
                                }
                                else {
                                    G2 = G_[depth+1][k-1][j-1];
                                }
                            }
                            LiftedGate -= L_[depth-1][i][k]*G2;
                            Ll_[depth][i][j] = Ll_[depth][i][j].eval() + lb;
                            Lu_[depth][i][j] = Lu_[depth][i][j].eval() + ub;
                            if(Ll_[depth][i][j].eval()<-1){
                                Ll_[depth][i][j] = -1;
                            }
                            if(Lu_[depth][i][j].eval()>1){
                                Lu_[depth][i][j] = 1;
                            }
//                            if (j%2!=0) {
//                                Ll_[depth][i+1][j+1].set_val(-Lu_[depth][i][j].eval());
//                                Lu_[depth][i+1][j+1].set_val(-Ll_[depth][i][j].eval());
//                            }
//                            else{
//                                Ll_[depth][i+1][j-1].set_val(Ll_[depth][i][j].eval());
//                                Lu_[depth][i+1][j-1].set_val(Lu_[depth][i][j].eval());
//                            }
                        }
//                        else {
//                            Ll[depth].set_val(i,j,0);
//                            Lu[depth].set_val(i,j,0);
//                            if (j%2==0) {
//                                Ll[depth].set_val(i+1,j+1,0);
//                                Lu[depth].set_val(i+1,j+1,0);
//                            }
//                            else{
//                                Ll[depth].set_val(i+1,j-1,0);
//                                Lu[depth].set_val(i+1,j-1,0);
//                            }
//                        }
                    }
                    Qdesign.add(LiftedGate==0);
                }
            }
        }

        /** Target gate constraint **/
        for (int i = 0; i < 2*m; i+=2) {
            for (int j = 0; j < 2*m; j++) {
                Constraint TargetGate("TargetGate_"+to_string(i+1)+to_string(j+1));
                for (int k = 0; k < 2*m; k++) {
                    if (!zeros[k][j]) {
                        if (k%2==0) {
                            G2 = G_[d-1][k][j];
                        }
                        else {
                            if (j%2==0) {//Imaginary
                                G2 = -1*G_[d-1][k-1][j+1];
                            }
                            else {
                                G2 = G_[d-1][k-1][j-1];
                            }
                        }
                        TargetGate += L_[d-3][i][k]*G2;
                    }
                }
                if(!relax){
                    TargetGate -=  slack_[i][j];
                }
//                auto key = to_string(i+1)+","+to_string(j+1);
//                if (TargetGate.get_nb_vars()>1) {
                    Qdesign.add(TargetGate==M.eval(i,j));
//                    Qdesign.add(TargetGate==0);
//                }
                
            }
        }
    }
    
//    for (size_t depth = 0; depth < d-1; depth++) {
//        for (size_t qubit1 = 0; qubit1 < n; qubit1++) {
//            Constraint NonConseqT("NonConseq_T_" +to_string(depth)+"_"+to_string(qubit1+1));
//            NonConseqT += zT_[depth][qubit1] + zTt_[depth+1][qubit1];
//            Qdesign.add(NonConseqT <= 1);
//            Constraint NonConseqTt("NonConseq_Tt_" +to_string(depth)+"_"+to_string(qubit1+1));
//            NonConseqTt += zTt_[depth][qubit1] + zTt_[depth+1][qubit1];
//            Qdesign.add(NonConseqTt <= 1);
//            Constraint NonConseqH("NonConseq_H_" +to_string(depth)+"_"+to_string(qubit1+1));
//            NonConseqH += zH_[depth][qubit1] + zH_[depth+1][qubit1];
//            Qdesign.add(NonConseqH <= 1);
//            for (size_t qubit2 = qubit1+1; qubit2 < n; qubit2++) {
//                Constraint NonConseqCnot("NonConseq_Cnot_" +to_string(depth)+"_"+to_string(qubit1+1)+"_"+to_string(qubit2+1));
//                NonConseqCnot += zCnot_[depth][qubit1][qubit2] + zCnot_[depth+1][qubit1][qubit2];
//                Qdesign.add(NonConseqCnot <= 1);
////                Constraint Symmetry1("Symmetry1" +to_string(depth)+"_"+to_string(qubit1+1)+"_"+to_string(qubit2+1));
////                Symmetry1 += zH_[depth][qubit1] - zH_[depth+1][qubit2];
////                Qdesign.add(Symmetry1 >= 0);
////                Constraint Symmetry2("Symmetry2" +to_string(depth)+"_"+to_string(qubit1+1)+"_"+to_string(qubit2+1));
////                Symmetry2 += zT_[depth][qubit1] - zT_[depth+1][qubit2];
////                Qdesign.add(Symmetry2 >= 0);
////                Constraint Symmetry3("Symmetry3" +to_string(depth)+"_"+to_string(qubit1+1)+"_"+to_string(qubit2+1));
////                Symmetry3 += zTt_[depth][qubit1] - zTt_[depth+1][qubit2];
////                Qdesign.add(Symmetry3 >= 0);
////                Constraint Symmetry4("Symmetry4" +to_string(depth)+"_"+to_string(qubit1+1)+"_"+to_string(qubit2+1));
////                Symmetry4 += zTt_[depth][qubit1] - zH_[depth+1][qubit2];
////                Qdesign.add(Symmetry4 >= 0);
////                for (size_t qubit3 = 0; qubit3 < n; qubit3++) {
////                    if (qubit3!=qubit1 && qubit3!=qubit2) {
////                        Constraint Symmetry5("Symmetry5" +to_string(depth)+"_"+to_string(qubit1+1)+"_"+to_string(qubit2+1)+"_"+to_string(qubit3+1));
////                        Symmetry5 += zH_[depth][qubit3] - zCnot_[depth][qubit1][qubit2];
////                        Qdesign.add(Symmetry5 >= 0);
////                    }
////                    if (qubit3<qubit2) {
////                        Constraint Symmetry6("Symmetry6" +to_string(depth)+"_"+to_string(qubit1+1)+"_"+to_string(qubit2+1)+"_"+to_string(qubit3+1));
////                        Symmetry6 += zT_[depth][qubit3] - zCnot_[depth+1][qubit1][qubit2];
////                        Qdesign.add(Symmetry6 <= 0);
////                        Constraint Symmetry7("Symmetry7" +to_string(depth)+"_"+to_string(qubit1+1)+"_"+to_string(qubit2+1)+"_"+to_string(qubit3+1));
////                        Symmetry7 += zTt_[depth][qubit3] - zCnot_[depth][qubit1][qubit2];
////                        Qdesign.add(Symmetry7 >= 0);
////                    }
////                }
//
//            }
//        }
//    }
//    Constraint TargetGate("TargetGate");
//    TargetGate += G[d-2]*G[d-1] - M;
//    Qdesign.add(TargetGate==0);
//    TargetGate.print();
    
//    Constraint TOFF_H("TOFF_H");
//    TOFF_H += sum(zCnot);
//    Qdesign.add(TOFF_H>=1);
////
//    Constraint TOFF_T("TOFF_T");
//    TOFF_T += sum(zT)+sum(zTt);
//    Qdesign.add(TOFF_T==0);

//    Constraint Config1("Config1");
//    Config1 += zCnot_[2][0][1];
//    Qdesign.add(Config1==1);
//    
//    Constraint Config2("Config2");
//    Config2 += zT_[5][1];
//    Qdesign.add(Config2==1);
////
//    Constraint Config3("Config3");
//    Config3 += zH_[6][1];
//    Qdesign.add(Config3==1);
//
//    
////    
//    Constraint TOFF_C1("TOFF_C1");
//    TOFF_C1 += zT_[5][1];
//    Qdesign.add(TOFF_C1==1);
//////////
//    Constraint TOFF_C2("TOFF_C2");
//    TOFF_C2 += zT_[4][0];
//    Qdesign.add(TOFF_C2==1);
////////
//    Constraint TOFF_C3("TOFF_C3");
//    TOFF_C3 += zTt_[2][0];
//    Qdesign.add(TOFF_C3==1);
//////////
//    Constraint TOFF_C4("TOFF_C4");
//    TOFF_C4 += zH_[0][0];
//    Qdesign.add(TOFF_C4==1);
//////
//    Constraint TOFF_C5("TOFF_C5");
//    TOFF_C5 += zTt_[2][1];
//    Qdesign.add(TOFF_C5==1);
//    //
//    Constraint TOFF_C6("TOFF_C6");
//    TOFF_C6 += lambda_r_[1];
//    Qdesign.add(TOFF_C6==-1);
////
//    Constraint TOFF_C7("TOFF_C7");
//    TOFF_C7 += sum(zS);
//    Qdesign.add(TOFF_C7>=7);
//
//    Constraint TOFF_C5("TOFF_C5");
//    TOFF_C5 += zCnotR_[1][0][1];
//    Qdesign.add(TOFF_C5==1);
//    
//    Constraint TOFF_C6("TOFF_C6");
//    TOFF_C6 += zCnotR_[3][0][1];
//    Qdesign.add(TOFF_C6==1);
//
//
    
//    Constraint TOFF_Tt2("TOFF_Tt2");
//    TOFF_Tt2 += sum(zCnotR);
//    Qdesign.add(TOFF_Tt2==2);

    
    Constraint TOFF_Tt("TOFF_Tt");
    TOFF_Tt += sum(zCnot);
    Qdesign.add(TOFF_Tt<=1);

    /** One gate per depth constraints **/
    for(int depth = 0; depth < d ; depth++){
        Constraint OneGate("OneGate_"+to_string(depth+1));
        for(int qubit = 0; qubit < n ; qubit++){
            if(ibm){
//                OneGate += zS_[depth][qubit] + zRx_[depth][qubit];
                if(add_rx){
                    OneGate += zS_[depth][qubit] + zSt_[depth][qubit] + zT_[depth][qubit] + zTt_[depth][qubit] + zRx_[depth][qubit];
                }
                else {
                    OneGate += zS_[depth][qubit] + zSt_[depth][qubit] + zT_[depth][qubit] + zTt_[depth][qubit];
                }
//                OneGate += zS_[depth][qubit] + zSt_[depth][qubit] + zT_[depth][qubit] + zTt_[depth][qubit];
            }
            else if (cont){
                if(add_rx){
                    OneGate += zS_[depth][qubit] + zRx_[depth][qubit];
                }
                else {
                    OneGate += zS_[depth][qubit];
                }
            }
            else {
                OneGate += zT_[depth][qubit] + zTt_[depth][qubit] + zH_[depth][qubit];
            }
            for (int qubit2 = qubit+1; qubit2 < n; qubit2++) {
                OneGate += zCnot_[depth][qubit][qubit2];
                if (rev_cnot) {
                    OneGate += zCnotR_[depth][qubit][qubit2];
                }
//                OneGate += zSwap_[depth][qubit][qubit2];
            }
        }
        Qdesign.add(OneGate==1);
    }
    Qdesign.print_expanded();
    auto solver_time_start = get_wall_time();
    if(relax){
        auto Relaxation = Qdesign.build_McCormick();
//        Relaxation->print_expanded();
        solver slvr(*Relaxation,cplex);
        slvr.run(output, relax = false, tol=1e-6, 0.01, "ma27");
        Qdesign.print_solution(false);
    }
    else {
        solver slvr(Qdesign,bonmin);
        slvr.run(output=0, relax = false, tol=1e-6, 0.01, "ma27");
        Qdesign.print_solution(false);
    }
    auto solver_time_end = get_wall_time();
    auto solve_time = solver_time_end - solver_time_start;
    
    auto total_time_start = get_wall_time();
    auto total_time_end = get_wall_time();
    auto total_time = total_time_end - total_time_start;
    DebugOn("Solver Time = " << to_string(solve_time) << " seconds" << endl);
    DebugOn("Total Computing Time = " << to_string(total_time) << " seconds" << endl);
    DebugOn("Optimal Objective = " << to_string(Qdesign._obj_val) << endl);
#else
    cerr << "Error: this version of Gravity "
    "was compiled without QPP support. Please rerun cmake with -DENABLE_QPP=1." << endl;
#endif
}
int main (int argc, char * argv[])
{
    int output = 0;
    //  Start Timers
    string path = argv[0];
    
    
    string mehrotra = "no", log_level="5", depth_str="2", nb_qubits_str="2", relax = "n";
    
    /** create a OptionParser with options */
    op::OptionParser opt;
    opt.add_option("h", "help", "shows option help"); // no default value means boolean options, which default value is false
    opt.add_option("l", "log", "Log level (def. 5)", log_level );
    opt.add_option("r", "relax", "Convex Relaxation y/n (def. n)", relax );
    opt.add_option("d", "depth", "Circuit depth", depth_str );
    opt.add_option("n", "nb_qubits", "Number of Qubits", nb_qubits_str );
    
    /** parse the options and verify that all went well. If not, errors and help will be shown */
    bool correct_parsing = opt.parse_options(argc, argv);
    
    if(!correct_parsing){
        return EXIT_FAILURE;
    }
    
    output = op::str2int(opt["l"]);
    bool has_help = op::str2bool(opt["h"]);
    /** show help */
    if(has_help) {
        opt.show_help();
        exit(0);
    }
    const int n = op::str2int(opt["n"]);
    const int d = op::str2int(opt["d"]);
    const int m = pow(2,n);
    
#ifdef USE_QPP
    /* Target Matrix */
    param<> M("M");
    M.set_size(2*m, 2*m);
    if(true){
        using namespace qpp;
        //        auto M1 = gt.expandout(gt.T, 0, n);
        //        auto M2 = gt.expandout(gt.CTRL(gt.X, {0}, {1}, n), 0, 1, pow(2,n));
        ////        auto M2 = gt.expandout(gt.T, 1, n);
        //        auto M3 = gt.expandout(gt.CTRL(gt.X, {1}, {2}, n), 0, 1, pow(2,n));
        //        auto M4 = gt.expandout(gt.T, 2, n);
        //        auto M5 = gt.expandout(gt.T, 1, n);
        //        auto M1 = gt.expandout(gt.T, 0, n);
        //        auto M1 = gt.expandout(gt.H, 2, n);
        //        auto M2 = gt.expandout(gt.CTRL(gt.X, {1}, {2}, n), 0, 1, m);
        //        auto M3 = gt.expandout(adjoint(gt.T), 2, n);
        //        auto M4 = gt.expandout(gt.CTRL(gt.X, {0}, {2}, n), 0, 1, m);
        //        auto M5 = gt.expandout(gt.T, 2, n);
        //        auto M6 = gt.expandout(gt.CTRL(gt.X, {1}, {2}, n), 0, 1, m);
        //        auto M7 = gt.expandout(adjoint(gt.T), 2, n);
        //        auto M8 = gt.expandout(gt.CTRL(gt.X, {0}, {2}, n), 0, 1, m);
        //        auto M9 = gt.expandout(gt.T, 2, n);
        //        auto M10 = gt.expandout(gt.H, 2, n);
//                auto M11 = gt.expandout(gt.H, 0, n);
//        DebugOn("H matrix at position " << to_string(0) <<" = " << endl);
//        DebugOn(disp(M11) << "\n");

        //        auto M12 = gt.expandout(gt.CTRL(gt.X, {0}, {1}, n), 0, 1, m);
        //        auto M13 = gt.expandout(adjoint(gt.T), 1, n);
        //        auto M14 = gt.expandout(gt.CTRL(gt.X, {0}, {1}, n), 0, 1, m);
        //        auto M15 = gt.expandout(gt.S, 1, n);
        //        auto M16 = gt.expandout(gt.T, 0, n);
//        auto theta = -pi/2;
//        auto M0 = gt.expandout(gt.Rn(theta, {0,0,1}), 0, n);
//        auto Mt0 = M0 * complex<double>(cos(theta/2),sin(theta/2));
//        DebugOn("Mt0 matrix at position " << to_string(0) <<" = " << endl);
//        DebugOn(disp(Mt0) << "\n");
                auto theta = pi;
                auto M1 = gt.expandout(gt.Rn(theta, {0,0,1}), 1, n);
                auto Mt1 = M1 * complex<double>(cos(theta/2),sin(theta/2));
//        DebugOn("Mt1 matrix at position " << to_string(0) <<" = " << endl);
//                DebugOn(disp(Mt1) << "\n");
        auto Mt2 = gt.expandout(gt.Rn(pi/2, {1,0,0}), 1, n);
//        auto Mp = Mt0*Rx;
//        
        //        auto Mt1 = M1 * complex<double>(cos(theta/2),sin(theta/2));
        
                theta = pi/2;
                auto M3 = gt.expandout(gt.Rn(theta, {0,0,1}), 1, n);
                auto Mt3 = M3 * complex<double>(cos(theta/2),sin(theta/2));
                auto Mt4 = gt.expandout(gt.CTRL(gt.X, {0}, {1}, n), 0, 1, m);
                auto Mt6 = gt.expandout(gt.Rn(pi/2, {1,0,0}), 1, n);
                theta = pi/2;
                auto M5 = gt.expandout(gt.Rn(theta, {0,0,1}), 1, n);
                auto Mt5 = M5 * complex<double>(cos(theta/2),sin(theta/2));
        auto Mp = Mt1*Mt2*Mt3*Mt4*Mt5*Mt6;
        ////        DebugOn("M5 matrix at position " << to_string(1) <<" = " << endl);
        ////        DebugOn(disp(M5) << "\n");
        //        auto Mt6 = gt.expandout(gt.Rn(pi/2, {1,0,0}), 1, n);
        
        
//                auto M1 = gt.expandout(gt.H, 0, n);
//                auto M2 = gt.expandout(gt.CTRL(gt.X, {1}, {0}, n), 0, 1, m);
//                auto M3 = gt.expandout(adjoint(gt.T), 0, n);
//                auto M4 = gt.expandout(gt.CTRL(gt.X, {1}, {0}, n), 0, 1, m);
//                auto M5 = gt.expandout(gt.T, 0, n);
//                auto M6 = gt.expandout(gt.T, 1, n);
//                auto M7 = gt.expandout(gt.H, 1, n);
//                auto Mp = M1*M2*M3*M4*M5*M6*M7;
        //        auto M5 = gt.expandout(gt.T, 2, n);
        //        auto M6 = gt.expandout(gt.CTRL(gt.X, {1}, {2}, n), 0, 1, m);
        //        auto M7 = gt.expandout(adjoint(gt.T), 2, n);
        //        auto M8 = gt.expandout(gt.CTRL(gt.X, {0}, {2}, n), 0, 1, m);
        //        auto M9 = gt.expandout(gt.T, 2, n);
        //        auto M10 = gt.expandout(gt.H, 2, n);
        //        auto M11 = gt.expandout(adjoint(gt.T), 1, n);
        //        auto M12 = gt.expandout(gt.CTRL(gt.X, {0}, {1}, n), 0, 1, m);
        //        auto M13 = gt.expandout(adjoint(gt.T), 1, n);
        //        auto M14 = gt.expandout(gt.T, 0, n);
        //        auto M15 = gt.expandout(gt.CTRL(gt.X, {0}, {1}, n), 0, 1, m);
        
        
//                auto Mp = gt.expandout(gt.H, 0, n);
        //        auto M2 = gt.expandout(gt.CTRL(gt.X, {1}, {2}, n), 0, 1, m);
        //        auto M3 = gt.expandout(adjoint(gt.T), 2, n);
        //        auto M4 = gt.expandout(gt.CTRL(gt.X, {0}, {2}, n), 0, 1, m);
        //        auto M5 = gt.expandout(gt.T, 2, n);
        //        auto M6 = gt.expandout(gt.CTRL(gt.X, {1}, {2}, n), 0, 1, m);
        //        auto M7 = gt.expandout(adjoint(gt.T), 2, n);
        //        auto M8 = gt.expandout(gt.CTRL(gt.X, {0}, {2}, n), 0, 1, m);
        //        auto M9 = gt.expandout(gt.T, 2, n);
        //        auto M11 = gt.expandout(gt.T, 1, n);
        //        auto M10 = gt.expandout(gt.H, 2, n);
        //        auto M12 = gt.expandout(gt.CTRL(gt.H, {1}, {2}, n), 0, 1, m);
        //        auto M13 = gt.expandout(gt.CTRL(gt.H, {2}, {1}, n), 0, 1, m);
        //        auto M14 = gt.expandout(gt.CTRL(gt.H, {1}, {2}, n), 0, 1, m);
        //        auto M15 = gt.expandout(gt.CTRL(gt.H, {0}, {2}, n), 0, 1, m);
        //        auto M16 = gt.expandout(gt.T, 0, n);
        //        auto M17 = gt.expandout(adjoint(gt.T), 2, n);
        //        auto M18 = gt.expandout(gt.CTRL(gt.H, {0}, {2}, n), 0, 1, m);
        //        auto M19 = gt.expandout(gt.CTRL(gt.H, {1}, {2}, n), 0, 1, m);
        //        auto M20 = gt.expandout(gt.CTRL(gt.H, {2}, {1}, n), 0, 1, m);
        //        auto M21 = gt.expandout(gt.CTRL(gt.H, {1}, {2}, n), 0, 1, m);
        
//                auto Mp = gt.expandout(gt.H, 0, n);
        //        auto M1 = gt.expandout(adjoint(gt.S), 0, n);
        //        auto M2 = gt.expandout(adjoint(gt.S), 1, n);
//                auto Mp = gt.expandout(gt.CTRL(gt.X, {1}, {0}, n), 0, 1, m);
        //        auto M4 = gt.expandout((gt.S), 1, n);
        //        auto M5 = gt.expandout(gt.CTRL(gt.X, {0}, {1}, n), 0, 1, m);
        
        //        auto Mp = M1*M2*M3*M4*M5*M6;
        //        auto Mp = Mt1*Mt2*Mt3*Mt4*Mt5*Mt6;
        //        auto Mp = gt.expandout(gt.Id(), 0, n);
        //        auto Mp = gt.expandout(gt.CZ, 0, 1, m);
        //            auto Mp = gt.expandout(gt.Fd(m), 0, 1, m);
        //        Mp *= 2/sqrt(2);
//        double theta = pi/2;
//        auto M1 = gt.expandout(gt.Rn(theta, {0,0,1}), 0, n);
//        auto Mt1 = M1 * complex<double>(cos(theta/2),sin(theta/2));
//        auto Mt2 = gt.expandout(gt.Rn(pi/2, {1,0,0}), 0, n);
//        auto Mt3 = gt.expandout(gt.CTRL(gt.X, {1}, {0}, n), 0, 1, m);
//        theta = 3*pi/4;
//        auto M4 = gt.expandout(gt.Rn(theta, {0,0,1}), 0, n);
//        auto Mt4 = M4 * complex<double>(cos(theta/2),sin(theta/2));
//        theta = -pi/4;
//        auto M5 = gt.expandout(gt.Rn(theta, {0,0,1}), 1, n);
//        auto Mt5 = M5 * complex<double>(cos(theta/2),sin(theta/2));
//        auto Mt6 = gt.expandout(gt.CTRL(gt.X, {0}, {1}, n), 0, 1, m);
//        auto Mt7 = gt.expandout(gt.Rn(pi/2, {1,0,0}), 0, n);
//        auto Mt8 = gt.expandout(gt.CTRL(gt.X, {1}, {0}, n), 0, 1, m);
//        theta = pi/2;
//        auto M9 = gt.expandout(gt.Rn(theta, {0,0,1}), 0, n);
//        auto Mt9 = M9 * complex<double>(cos(theta/2),sin(theta/2));
//        theta = 7*pi/4;
//        auto M10 = gt.expandout(gt.Rn(theta, {0,0,1}), 1, n);
//        auto Mt10 = M10 * complex<double>(cos(theta/2),sin(theta/2));
//        auto Mp = Mt1*Mt2*Mt3*Mt4*Mt5*Mt6*Mt7*Mt8*Mt9*Mt10;
//        auto Mp = M1*M2*M3*M4;
//        auto Mp = gt.expandout(gt.Fd(m), 0, 1, m);
//        auto Mp = gt.expandout(gt.CTRL(gt.Z, {0}, {1}, n), 0, 1, m);
        //        auto Mp = gt.expandout(gt.SWAP, 0, 1, m);
        //        auto M1 = gt.expandout(gt.Rn(pi/2, {1,0,0}), 0, n);
        //        auto Mp = gt.expandout(gt.T, 1, n);
        //        auto Mp = M1*M2*M3*M4*M5*M6*M7*M8*M9*M10*M11*M12*M13*M14*M15*M16*M17*M18*M19*M20*M21;
        //        auto Mp = M1*M2*M3*M4*M5*M6*M7*M8*M9*M10*M11*M12*M13*M14*M15*M16;
//        auto Mp = gt.expandout(gt.TOF, 0, 1, m);
//        auto Mp = gt.expandout(gt.FRED, 0, 1, m);
        DebugOn("Target matrix = " << endl);
        DebugOn(disp(Mp) << "\n");
        M.set_vals(Mp.sparseView());
        DebugOff(M.to_str(true));
        exit(1);
//        run_relaxation(n, d, M);
        if (opt["r"]=="y") {
            run_qtdesign(n, d, M, true);
        }
        else {
            run_qtdesign(n, d, M, false);
        }
    }
    
    
//    Model mod("test");
//    var<> x("x",-10,10);
//    var<> y("y",-10,10);
//    mod.add(x);
//    mod.add(y);
//    Constraint C1("C1");
//    C1 += y - x*x;
//    mod.add(C1>=0);
//    
//    Constraint C2("C2");
//    C2 += y - x*x*(x-2) + 1e-5;
//    mod.add(C1<=0);
//    
//    mod.min(x);
//    
//    solver slvr(mod,ipopt);
//    slvr.run(output, false, 1e-6, 0.01, "ma27");
//    mod.print_solution(false);
    
//    auto Mres= run_reverse(n,d/2+2,M);
//    run_qtdesign(n, d/2, Mres);
//    run_qtdesign_fixed(2,10);
#endif
    return 0;
}
