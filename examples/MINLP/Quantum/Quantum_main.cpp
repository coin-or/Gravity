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


int main (int argc, char * argv[])
{
#ifdef USE_QPP
    //  Start Timers
    string path = argv[0];
    int output = 0;
    bool relax = false;
    
    double tol = 1e-6;
    string mehrotra = "no", log_level="5", depth_str="2", nb_qubits_str="2";
    
    /** create a OptionParser with options */
    op::OptionParser opt;
    opt.add_option("h", "help", "shows option help"); // no default value means boolean options, which default value is false
    opt.add_option("l", "log", "Log level (def. 5)", log_level );
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
    
    /** Parameters **/
    
    /* Qubit number */
    const unsigned n = op::str2int(opt["n"]);
    const unsigned m = pow(2,n);
    DebugOn("Number of Qubits = " << n << endl);
    /* Circuit depth */
    const unsigned d = op::str2int(opt["d"]);;
    DebugOn("Depth = " << d << endl);
    if (d<2) {
        cerr << "Depth should be greater or equal than 2!" << endl;
        return -1;
    }
    /* T and T conjugate transpose gate matrices, one for each row (qubit) */
    param<> T[n], Tt[n];
    
    /* Hadamard gate matrices */
    param<> H[n];
    
    /* Cnot gate matrices */
    param<> Cnot[n*(n-1)/2];
    
    
    /* Target Matrix */
    param<> M("M");
    M.set_size(2*m, 2*m, 0);
    {
        using namespace qpp;
//        auto M1 = gt.expandout(gt.T, 0, n);
//        auto M2 = gt.expandout(gt.CTRL(gt.X, {0}, {1}, n), 0, 1, pow(2,n));
////        auto M2 = gt.expandout(gt.T, 1, n);
//        auto M3 = gt.expandout(gt.CTRL(gt.X, {1}, {2}, n), 0, 1, pow(2,n));
//        auto M4 = gt.expandout(gt.T, 2, n);
//        auto M5 = gt.expandout(gt.T, 1, n);
//        auto M1 = gt.expandout(gt.T, 0, n);
        auto M1 = gt.expandout(gt.CTRL(gt.X, {0}, {1}, n), 0, 1, pow(2,n));
        auto M2 = gt.expandout(gt.T, 1, n);
        auto M3 = gt.expandout(gt.T, 0, n);
        auto M4 = gt.expandout(gt.CTRL(gt.X, {1}, {2}, n), 0, 1, pow(2,n));
        auto M5 = gt.expandout(gt.T, 2, n);
        
//        auto M2 = gt.expandout(gt.CNOT, 0, 1, m);
//        auto M2 = gt.expandout(gt.CTRL(gt.X, {0}, {1}, n), 0, 1, pow(2,n));
//        auto M3 = gt.expandout(gt.H, 2, n);
        auto Mp = M1*M2*M3*M4*M5;
//        auto Mp = gt.expandout(gt.Id(), 0, n);
//        auto Mp = gt.expandout(gt.TOF, 0, 1, m);
        DebugOn("Target matrix = " << endl);
        DebugOn(disp(Mp) << "\n");
        M.set_vals(Mp.sparseView());
        DebugOff(M.to_str(true));
    }
    /** Initialization **/
    unsigned idx = 0;
    for (auto i = 0; i<n; i++) { // Qubit (row) iterator
        T[i].set_name("T_"+to_string(i+1));
        T[i].set_size(2*m, 2*m, 0);/* Representing Complex number as 2x2 real matrices */
        T[i].QuantumT(i,n);/* Fill nonzero values corresponding to a T gate */
        Tt[i].set_name("Tt_"+to_string(i+1));
        Tt[i].set_size(2*m, 2*m, 0);
        Tt[i].QuantumT(i,n,true);/* Fill nonzero values corresponding to a T conjugate transpose gate */
        H[i].set_name("H_"+to_string(i+1));
        H[i].set_size(2*m, 2*m, 0);
        H[i].QuantumH(i,n);/* Fill nonzero values corresponding to an H gate */
        for (auto k = i+1; k<n; k++) {
            Cnot[idx].set_name("Cnot_"+to_string(i+1)+to_string(k+1));
            Cnot[idx].set_size(2*m, 2*m, 0);
            Cnot[idx].QuantumCnot(i,k,n);/* Fill nonzero values corresponding to a CNOT gate */
            idx++;
        }
    }
    
    
    /** Indices **/
    auto d_ids = indices(1,d);/*< Depth indices */
    auto q_ids = indices(1,n);/*< Qubit indices */
    auto pairs = indices(ordered_pairs(2, n+1)); /*< Indices pairs (control-taget) for Cnot gates */
    auto m2 = indices(indices(1,2*m), indices(1,2*m)); /*< Indices for square mxm matrices */
    
    
    /** Variables **/
    /** Binaries **/
    var<bool> zT("zT"); /*< Binary variables for T gates */
    var<bool> zTt("zTt"); /*< Binary variables for T conjugate gates */
    var<bool> zH("zH"); /*< Binary variables for H gates */
    var<bool> zCnot("zCnot"); /*< Binary variables for CNOT gates */
    /** Gate Matrices **/
    var<> G[d]; /*< Base Gate Matrices */
    bool zeros[2*m][2*m]; /*< Zeros in Gate Matrices */
    var<> L[d-2]; /*< Lifted Product Matrices */
    
    
    
    /** Model **/
    
    /** Adding Variables **/
    Model Qdesign("Qdesign");
    /** Binaries **/
    Qdesign.add(zT.in(d_ids,q_ids), zTt.in(d_ids,q_ids), zH.in(d_ids,q_ids), zCnot.in(d_ids,pairs));
    DebugOff("Total number of variables after adding the binaries = " << Qdesign.get_nb_vars() << endl);
    /** Gate Matrices **/
    for (unsigned depth = 0; depth < d; depth++) {
        G[depth] = var<>("G_"+to_string(depth+1), -1/sqrt(2),1);
        Qdesign.add(G[depth].in(m2));
    }
    if (d>2) {
        for (unsigned depth = 0; depth < d-2; depth++) {
            L[depth] = var<>("L_"+to_string(depth+1), -3,3);
            Qdesign.add(L[depth].in(m2));
        }
    }

    DebugOff("Total number of variables after adding the G vars= " << Qdesign.get_nb_vars() << endl);
    
    
    /** Adding Constraints **/
    /** Building indexed variables as matrices (for a faster access in constraints) **/
    var<bool> zT_[d][n], zTt_[d][n], zH_[d][n], zCnot_[d][n][n];
    var<> G_[d][2*m][2*m], L_[d-2][2*m][2*m];
    for (size_t depth = 0; depth < d; depth++) {
        for (size_t qubit1 = 0; qubit1 < n; qubit1++) {
            zT_[depth][qubit1] = zT(depth+1,qubit1+1);
            zTt_[depth][qubit1] = zTt(depth+1,qubit1+1);
            zH_[depth][qubit1] = zH(depth+1,qubit1+1);
            for (size_t qubit2 = qubit1+1; qubit2 < n; qubit2++) {
                zCnot_[depth][qubit1][qubit2] = zCnot(depth+1,qubit1+1,qubit2+1);
            }
        }
        for (unsigned i = 0; i < 2*m; i++) {
            for (unsigned j = 0; j < 2*m; j++) {
                G_[depth][i][j] = G[depth](i+1,j+1);
            }
        }
    }
    if (d>2) {
        for (size_t depth = 0; depth < d-2; depth++) {
            for (unsigned i = 0; i < 2*m; i++) {
                for (unsigned j = 0; j < 2*m; j++) {
                    L_[depth][i][j] = L[depth](i+1,j+1);
                }
            }
        }
    }
    
    /** Objective **/
    Qdesign.min(sum(zH) + sum(zT) + sum(zTt) + sum(zCnot));
    
//    func_ obj;
//    for (size_t depth = 0; depth < d; depth++) {
//        for (size_t qubit1 = 0; qubit1 < n; qubit1++) {
//            obj -= zT_[depth][qubit1]*zT_[depth][qubit1] + zTt_[depth][qubit1]*zTt_[depth][qubit1] + zH_[depth][qubit1]*zH_[depth][qubit1];
//            for (size_t qubit2 = qubit1+1; qubit2 < n; qubit2++) {
//                obj += 1e3*zCnot_[depth][qubit1][qubit2];
////                obj -= zCnot_[depth][qubit1][qubit2]*zCnot_[depth][qubit1][qubit2];
//            }
//        }
//    }
//    Qdesign.min(obj);

    
    /** Base gate matrix definition **/
    for (unsigned depth = 0; depth < d; depth++) {
        for (unsigned i = 0; i < 2*m; i++) {
            for (unsigned j = 0; j < 2*m; j++) {
                Constraint BaseGate("BaseGate_"+to_string(depth+1)+"_"+to_string(i+1)+to_string(j+1));
                BaseGate = G_[depth][i][j];
                unsigned idx = 0;
                for (unsigned qubit = 0; qubit < n; qubit++) {
                    BaseGate -= zT_[depth][qubit]*T[qubit].eval(i,j) + zTt_[depth][qubit]*Tt[qubit].eval(i,j) + zH_[depth][qubit]*H[qubit].eval(i,j);
                    for (unsigned qubit2 = qubit+1; qubit2 < n; qubit2++) {
                        BaseGate -= zCnot_[depth][qubit][qubit2]*Cnot[idx++].eval(i,j);
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
        for (unsigned i = 0; i < 2*m; i++) {
            for (unsigned j = 0; j < 2*m; j++) {
                Constraint TargetGate("TargetGate_"+to_string(i+1)+to_string(j+1));
                for (unsigned k = 0; k < 2*m; k++) {
                    if (!zeros[i][k] && !zeros[k][j]) {
                        TargetGate += G_[d-2][i][k]*G_[d-1][k][j];
                    }
                }
                if (TargetGate.get_nb_vars()>0) {
                    Qdesign.add(TargetGate==M.eval(i,j));
                }
            }
        }
    }
    else {
        /** Lifted gate constraint **/
        /** Depth 0 **/
        for (unsigned i = 0; i < 2*m; i++) {
            for (unsigned j = 0; j < 2*m; j++) {
                Constraint LiftedGate("LiftedGate_0_"+to_string(i+1)+to_string(j+1));
                LiftedGate += L_[0][i][j];
                for (unsigned k = 0; k < 2*m; k++) {
                    if (!zeros[i][k] && !zeros[k][j]) {
                        LiftedGate -= G_[0][i][k]*G_[1][k][j];
                    }
                }
                Qdesign.add(LiftedGate==0);
            }
        }
        for (unsigned depth = 1; depth < d-2; depth++) {
            for (unsigned i = 0; i < 2*m; i++) {
                for (unsigned j = 0; j < 2*m; j++) {
                    Constraint LiftedGate("LiftedGate"+to_string(depth)+"_"+to_string(i+1)+to_string(j+1));
                    LiftedGate += L_[depth][i][j];
                    for (unsigned k = 0; k < 2*m; k++) {
                        if (!zeros[k][j]) {
                            LiftedGate -= L_[depth-1][i][k]*G_[depth+1][k][j];
                        }
                    }
                    Qdesign.add(LiftedGate==0);
                }
            }
        }

        /** Target gate constraint **/
        for (unsigned i = 0; i < 2*m; i+=2) {
            for (unsigned j = 0; j < 2*m; j++) {
                Constraint TargetGate("TargetGate_"+to_string(i+1)+to_string(j+1));
                for (unsigned k = 0; k < 2*m; k++) {
                    if (!zeros[k][j]) {
                        TargetGate += L_[d-3][i][k]*G_[d-1][k][j];
                    }
                }
                if (TargetGate.get_nb_vars()>0) {
                    Qdesign.add(TargetGate==M.eval(i,j));
                }
                
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
//                Constraint Symmetry1("Symmetry1" +to_string(depth)+"_"+to_string(qubit1+1)+"_"+to_string(qubit2+1));
//                Symmetry1 += zH_[depth][qubit1] - zH_[depth+1][qubit2];
//                Qdesign.add(Symmetry1 >= 0);
////                Constraint Symmetry2("Symmetry2" +to_string(depth)+"_"+to_string(qubit1+1)+"_"+to_string(qubit2+1));
////                Symmetry2 += zT_[depth][qubit1] - zT_[depth+1][qubit2];
////                Qdesign.add(Symmetry2 >= 0);
////                Constraint Symmetry3("Symmetry3" +to_string(depth)+"_"+to_string(qubit1+1)+"_"+to_string(qubit2+1));
////                Symmetry3 += zTt_[depth][qubit1] - zTt_[depth+1][qubit2];
////                Qdesign.add(Symmetry3 >= 0);
////                Constraint Symmetry4("Symmetry4" +to_string(depth)+"_"+to_string(qubit1+1)+"_"+to_string(qubit2+1));
////                Symmetry4 += zTt_[depth][qubit1] - zH_[depth+1][qubit2];
////                Qdesign.add(Symmetry4 >= 0);
//                for (size_t qubit3 = 0; qubit3 < n; qubit3++) {
////                    if (qubit3!=qubit1 && qubit3!=qubit2) {
////                        Constraint Symmetry5("Symmetry5" +to_string(depth)+"_"+to_string(qubit1+1)+"_"+to_string(qubit2+1)+"_"+to_string(qubit3+1));
////                        Symmetry5 += zH_[depth][qubit3] - zCnot_[depth][qubit1][qubit2];
////                        Qdesign.add(Symmetry5 >= 0);
////                    }
//                    if (qubit3<qubit2) {
//                        Constraint Symmetry6("Symmetry6" +to_string(depth)+"_"+to_string(qubit1+1)+"_"+to_string(qubit2+1)+"_"+to_string(qubit3+1));
//                        Symmetry6 += zT_[depth][qubit3] - zCnot_[depth+1][qubit1][qubit2];
//                        Qdesign.add(Symmetry6 <= 0);
////                        Constraint Symmetry7("Symmetry7" +to_string(depth)+"_"+to_string(qubit1+1)+"_"+to_string(qubit2+1)+"_"+to_string(qubit3+1));
////                        Symmetry7 += zTt_[depth][qubit3] - zCnot_[depth][qubit1][qubit2];
////                        Qdesign.add(Symmetry7 >= 0);
//                    }
//                }
//
//            }
//        }
//    }
//    Constraint TargetGate("TargetGate");
//    TargetGate += G[d-2]*G[d-1] - M;
//    Qdesign.add(TargetGate==0);
//    TargetGate.print();
    
    /** One gate per depth constraints **/
    for(unsigned depth = 0; depth < d ; depth++){
        Constraint OneGate("OneGate_"+to_string(depth+1));
        for(unsigned qubit = 0; qubit < n ; qubit++){
            OneGate += zT_[depth][qubit] + zTt_[depth][qubit] + zH_[depth][qubit];
            for (unsigned qubit2 = qubit+1; qubit2 < n; qubit2++) {
                OneGate += zCnot_[depth][qubit][qubit2];
            }
        }
        Qdesign.add(OneGate==1);
//        OneGate.print();
    }
    Qdesign.print_expanded();
    
    solver slvr(Qdesign,ipopt);
    auto solver_time_start = get_wall_time();
    slvr.run(output, relax = false, tol=1e-6, 0.01, "ma27");
    Qdesign.print_solution();
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
    
    return 0;
}
