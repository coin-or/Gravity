//
//  SOCOPF.cpp
//  Gravity
//
//  Created by Hassan Hijazi on 21 January 2018.
//
//
#include <stdio.h>
#include <PowerNet.h>
#include <gravity/solver.h>
#include <optionParser.hpp>

using namespace std;
using namespace gravity;

int main (int argc, char * argv[])
{
    int output = 0;
    double tol = 1e-6;
    double solver_time_end, total_time_end, solve_time, total_time;
    string mehrotra = "no", log_level="0";
    
    // Specify the additional constraints
    bool current = true;
    bool current_partition_on_off_automated = true;
    
    //    Specify the use of partitioning scheme without current
    string model_type = "Model_II"; //the default relaxation model is Model_II
    
    //    Switch the data file to another instance
//    string fname = string(prj_dir)+"/data_sets/Power/nesta_case5_pjm.m";
    string fname = string(prj_dir)+"/data_sets/Power/nesta_case9_bgm__nco_tree.m";
//    string fname = string(prj_dir)+"/data_sets/Power/nesta_case39_1_bgm__nco.m";
    
    string path = argv[0];
    string solver_str="ipopt";
    
    /** Create a OptionParser with options */
    op::OptionParser opt;
    opt.add_option("h", "help", "shows option help"); // no default value means boolean options, which default value is false
    opt.add_option("f", "file", "Input file name", fname);
    opt.add_option("l", "log", "Log level (def. 0)", log_level );
    
    /** Parse the options and verify that all went well. If not, errors and help will be shown */
    bool correct_parsing = opt.parse_options(argc, argv);
    
    if(!correct_parsing){
        return EXIT_FAILURE;
    }
    
    output = op::str2int(opt["l"]);
    
    fname = opt["f"];
    bool has_help = op::str2bool(opt["h"]);
    if(has_help) {
        opt.show_help();
        exit(0);
    }
    
    double total_time_start = get_wall_time();
    PowerNet grid;
    grid.readgrid(fname);
    grid.get_tree_decomp_bags();
    
    /* Grid Stats */
    auto bus_pairs = grid.get_bus_pairs();
    auto nb_gen = grid.get_nb_active_gens();
    auto nb_lines = grid.get_nb_active_arcs();
    auto nb_buses = grid.get_nb_active_nodes();
    DebugOn("nb active gens = " << nb_gen << endl);
    DebugOn("nb active lines = " << nb_lines << endl);
    DebugOn("nb active buses = " << nb_buses << endl);
    
    /** Sets */
    auto nodes = indices(grid.nodes);
    auto arcs = indices(grid.arcs);
    auto gens = indices(grid.gens);
    auto gen_nodes = grid.gens_per_node();
    auto out_arcs = grid.out_arcs_per_node();
    auto in_arcs = grid.in_arcs_per_node();
    
    /* Grid Parameters */
    auto pg_min = grid.pg_min.in(gens);
    auto pg_max = grid.pg_max.in(gens);
    auto qg_min = grid.qg_min.in(gens);
    auto qg_max = grid.qg_max.in(gens);
    auto c1 = grid.c1.in(gens);
    auto c2 = grid.c2.in(gens);
    auto c0 = grid.c0.in(gens);
    auto lij_min=grid.lij_min.in(arcs);
    auto lij_max=grid.lij_max.in(arcs);
    auto lji_min=grid.lji_min.in(arcs);
    auto lji_max=grid.lji_max.in(arcs);
    auto cc=grid.cc.in(arcs);
    auto dd=grid.dd.in(arcs);
    auto ch_half=grid.ch_half.in(arcs);
    auto b = grid.b.in(arcs);
    auto g = grid.g.in(arcs);
    auto tr = grid.tr.in(arcs);
    
    
    
    /** MODEL DECLARATION */
    Model<> SOCP("SCOPF Model"); //model for the relaxation
    Model<> ACOPF("ACOPF Model"); //exact model
    
    /** Variables */
    /* power generation variables */
    var<> Pg("Pg", pg_min, pg_max);
    var<> Qg ("Qg", qg_min, qg_max);
    SOCP.add(Pg.in(gens),Qg.in(gens));
    ACOPF.add(Pg.in(gens),Qg.in(gens));
    
    /* power flow variables */
    var<> Pf_from("Pf_from", -1*grid.S_max,grid.S_max);
    var<> Qf_from("Qf_from", -1*grid.S_max,grid.S_max);
    var<> Pf_to("Pf_to", -1*grid.S_max,grid.S_max);
    var<> Qf_to("Qf_to", -1*grid.S_max,grid.S_max);
    SOCP.add(Pf_from.in(arcs), Qf_from.in(arcs), Pf_to.in(arcs), Qf_to.in(arcs));
    ACOPF.add(Pf_from.in(arcs), Qf_from.in(arcs), Pf_to.in(arcs), Qf_to.in(arcs));
    
    /* Real part of Wij = ViVj */
    var<>  R_Wij("R_Wij", grid.wr_min, grid.wr_max);
    /* Imaginary part of Wij = ViVj */
    var<>  Im_Wij("Im_Wij", grid.wi_min, grid.wi_max);
    /* Magnitude of Wii = Vi^2 */
    var<>  Wii("Wii", grid.w_min, grid.w_max);
    SOCP.add(Wii.in(nodes));
    SOCP.add(R_Wij.in(bus_pairs));
    SOCP.add(Im_Wij.in(bus_pairs));
    ACOPF.add(Wii.in(nodes));
    ACOPF.add(R_Wij.in(bus_pairs));
    ACOPF.add(Im_Wij.in(bus_pairs));
    
    var<> lij("lij", lij_min,lij_max);
    var<> lji("lji", lji_min,lji_max);
    if(current){
        SOCP.add(lij.in(arcs),lji.in(arcs));
    }
    
    /* Initialize variables */
    R_Wij.initialize_all(1.0);
    Wii.initialize_all(1.001);
    
    /************** Add the lifted variables *****************/
    if(current){
        
        param<Cpx> T("T"), Y("Y"), Ych("Ych");
        var<Cpx> L_from("L_from"), Wij("Wij");
        T.real_imag(cc.in(arcs), dd.in(arcs));
        Y.real_imag(g.in(arcs), b.in(arcs));
        Ych.set_imag(ch_half.in(arcs));
        
        
        L_from.set_real(lij.in(arcs));
        Wij.real_imag(R_Wij.in_pairs(arcs), Im_Wij.in_pairs(arcs));
        var<Cpx> Sij("Sij"), Sji("Sji");
        Sij.real_imag(Pf_from.in(arcs), Qf_from.in(arcs));
        Sji.real_imag(Pf_to.in(arcs), Qf_to.in(arcs));
        
        
        Constraint<Cpx> I_from("I_from");
        I_from=(Y+Ych)*(conj(Y)+conj(Ych))*Wii.from(arcs)-T*Y*(conj(Y)+conj(Ych))*conj(Wij)-conj(T)*conj(Y)*(Y+Ych)*Wij+pow(tr,2)*Y*conj(Y)*Wii.to(arcs);
        SOCP.add_real(I_from.in(arcs)==pow(tr,2)*L_from);
        
        var<Cpx> L_to("L_to");
        L_to.set_real(lji.in(arcs));
        
        Constraint<Cpx> I_to("I_to");
        I_to=pow(tr,2)*(Y+Ych)*(conj(Y)+conj(Ych))*Wii.to(arcs)-conj(T)*Y*(conj(Y)+conj(Ych))*Wij-T*conj(Y)*(Y+Ych)*conj(Wij)+Y*conj(Y)*Wii.from(arcs);
        SOCP.add_real(I_to.in(arcs)==pow(tr,2)*L_to);
        
        Constraint<> I_from_Pf("I_from_Pf");
        I_from_Pf=lij*Wii.from(arcs)-pow(tr,2)*(pow(Pf_from,2) + pow(Qf_from,2));
        //        SOCP.add(I_from_Pf.in(arcs)==0,true);
        SOCP.add(I_from_Pf.in(arcs)>=0);
        //        SOCP.get_constraint("I_from_Pf")->_relaxed = true;
        
        Constraint<> I_to_Pf("I_to_Pf");
        I_to_Pf=lji.in(arcs)*Wii.to(arcs)-(pow(Pf_to.in(arcs),2) + pow(Qf_to.in(arcs), 2));
        SOCP.add(I_to_Pf.in(arcs)>=0);
        SOCP.get_constraint("I_to_Pf")->_relaxed = true;
        
        /* Second-order cone */
        Constraint<> SOC("SOC");
        SOC = pow(R_Wij.in(bus_pairs), 2) + pow(Im_Wij.in(bus_pairs), 2) - Wii.from(bus_pairs)*Wii.to(bus_pairs);
        SOCP.add(SOC.in(bus_pairs) <= 0);
        
    }
    
    
    
    
    /**  Objective */
    auto obj = product(c1,Pg) + product(c2,pow(Pg,2)) + sum(c0);
    SOCP.min(obj);
    ACOPF.min(obj);
    
    /** Constraints */
    
    
    /* Equality of Second-order cone (for upperbound) */
    Constraint<> Equality_SOC("Equality_SOC");
    Equality_SOC = pow(R_Wij, 2) + pow(Im_Wij, 2) - Wii.from(bus_pairs)*Wii.to(bus_pairs);
    ACOPF.add(Equality_SOC.in(bus_pairs) == 0);
    
//    SOCP.add(Equality_SOC.in(bus_pairs) == 0, true);
    
    /* Flow conservation */
    Constraint<> KCL_P("KCL_P");
    KCL_P  = sum(Pf_from, out_arcs) + sum(Pf_to, in_arcs) + grid.pl - sum(Pg, gen_nodes) + grid.gs*Wii;
    SOCP.add(KCL_P.in(nodes) == 0);
    ACOPF.add(KCL_P.in(nodes) == 0);
    
    Constraint<> KCL_Q("KCL_Q");
    KCL_Q  = sum(Qf_from, out_arcs) + sum(Qf_to, in_arcs) + grid.ql - sum(Qg, gen_nodes) - grid.bs*Wii;
    SOCP.add(KCL_Q.in(nodes) == 0);
    ACOPF.add(KCL_Q.in(nodes) == 0);
    
    /* AC Power Flow */
    Constraint<> Flow_P_From("Flow_P_From");
    Flow_P_From = Pf_from - (grid.g_ff*Wii.from(arcs) + grid.g_ft*R_Wij + grid.b_ft*Im_Wij);
    SOCP.add(Flow_P_From.in(arcs) == 0);
    ACOPF.add(Flow_P_From.in(arcs) == 0);
    
    Constraint<> Flow_P_To("Flow_P_To");
    Flow_P_To = Pf_to - (grid.g_tt*Wii.to(arcs) + grid.g_tf*R_Wij - grid.b_tf*Im_Wij);
    SOCP.add(Flow_P_To.in(arcs) == 0);
    ACOPF.add(Flow_P_To.in(arcs) == 0);
    
    Constraint<> Flow_Q_From("Flow_Q_From");
    Flow_Q_From = Qf_from - (grid.g_ft*Im_Wij - grid.b_ff*Wii.from(arcs) - grid.b_ft*R_Wij);
    SOCP.add(Flow_Q_From.in(arcs) == 0);
    ACOPF.add(Flow_Q_From.in(arcs) == 0);
    
    Constraint<> Flow_Q_To("Flow_Q_To");
    Flow_Q_To = Qf_to + grid.b_tt*Wii.to(arcs) + grid.b_tf*R_Wij + grid.g_tf*Im_Wij;
    SOCP.add(Flow_Q_To.in(arcs) == 0);
    ACOPF.add(Flow_Q_To.in(arcs) == 0);
    
    /* Phase Angle Bounds constraints */
    Constraint<> PAD_UB("PAD_UB");
    PAD_UB = Im_Wij;
    PAD_UB <= grid.tan_th_max*R_Wij;
    SOCP.add(PAD_UB.in(bus_pairs));
    ACOPF.add(PAD_UB.in(bus_pairs));
    
    Constraint<> PAD_LB("PAD_LB");
    PAD_LB =  Im_Wij;
    PAD_LB >= grid.tan_th_min*R_Wij;
    SOCP.add(PAD_LB.in(bus_pairs));
    ACOPF.add(PAD_LB.in(bus_pairs));
    
    /* Thermal Limit Constraints */
    Constraint<> Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from = pow(Pf_from, 2) + pow(Qf_from, 2);
    Thermal_Limit_from <= pow(grid.S_max,2);
    SOCP.add(Thermal_Limit_from.in(arcs));
    ACOPF.add(Thermal_Limit_from.in(arcs));
    
    Constraint<> Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to = pow(Pf_to, 2) + pow(Qf_to, 2);
    Thermal_Limit_to <= pow(grid.S_max,2);
    SOCP.add(Thermal_Limit_to.in(arcs));
    ACOPF.add(Thermal_Limit_to.in(arcs));
    
    /* Lifted Nonlinear Cuts */
    Constraint<> LNC1("LNC1");
    LNC1 += (grid.v_min.from(bus_pairs)+grid.v_max.from(bus_pairs))*(grid.v_min.to(bus_pairs)+grid.v_max.to(bus_pairs))*(grid.sphi*Im_Wij + grid.cphi*R_Wij);
    LNC1 -= grid.v_max.to(bus_pairs)*grid.cos_d*(grid.v_min.to(bus_pairs)+grid.v_max.to(bus_pairs))*Wii.from(bus_pairs);
    LNC1 -= grid.v_max.from(bus_pairs)*grid.cos_d*(grid.v_min.from(bus_pairs)+grid.v_max.from(bus_pairs))*Wii.to(bus_pairs);
    LNC1 -= grid.v_max.from(bus_pairs)*grid.v_max.to(bus_pairs)*grid.cos_d*(grid.v_min.from(bus_pairs)*grid.v_min.to(bus_pairs) - grid.v_max.from(bus_pairs)*grid.v_max.to(bus_pairs));
    SOCP.add(LNC1.in(bus_pairs) >= 0);
    ACOPF.add(LNC1.in(bus_pairs) >= 0);
    
    Constraint<> LNC2("LNC2");
    LNC2 += (grid.v_min.from(bus_pairs)+grid.v_max.from(bus_pairs))*(grid.v_min.to(bus_pairs)+grid.v_max.to(bus_pairs))*(grid.sphi*Im_Wij + grid.cphi*R_Wij);
    LNC2 -= grid.v_min.to(bus_pairs)*grid.cos_d*(grid.v_min.to(bus_pairs)+grid.v_max.to(bus_pairs))*Wii.from(bus_pairs);
    LNC2 -= grid.v_min.from(bus_pairs)*grid.cos_d*(grid.v_min.from(bus_pairs)+grid.v_max.from(bus_pairs))*Wii.to(bus_pairs);
    LNC2 += grid.v_min.from(bus_pairs)*grid.v_min.to(bus_pairs)*grid.cos_d*(grid.v_min.from(bus_pairs)*grid.v_min.to(bus_pairs) - grid.v_max.from(bus_pairs)*grid.v_max.to(bus_pairs));
    SOCP.add(LNC2.in(bus_pairs) >= 0);
    ACOPF.add(LNC2.in(bus_pairs) >= 0);
    
    
    /***************** CALLING OBBT BEFORE CALLING RUN **********************/
    //    SOCP.reset_constrs();
    double max_time = 100000;
    int max_iter = 5;
    int precision = 4;
    double upperbound = grid.solve_acopf(ACRECT);
    SOCP.run_obbt(max_time,max_iter,{true,upperbound},precision);
    auto original_SOC = grid.build_SCOPF();
    solver<> SOCOPF_ORIG(original_SOC, ipopt);
    SOCOPF_ORIG.run(output, tol = 1e-6);
    double original_LB = original_SOC->get_obj_val();
    double gap = 100*(upperbound - original_LB)/upperbound;
    DebugOn("Initial Gap = " << to_string(gap) << "%."<<endl);
    gap = 100*(upperbound - SOCP.get_obj_val())/upperbound;
    DebugOn("Gap after OBBT = " << to_string(gap) << "%."<<endl);
    
    auto nonzero_idx = SOCP.sorted_nonzero_constraint_indices(tol, true, "I_to_Pf");
    nonzero_idx.print();
    
    if(current){
        //THESE ARE ALREADY INCLUDED BEFORE OBBT
        
        //        param<Cpx> T("T"), Y("Y"), Ych("Ych");
        //        var<Cpx> L_from("L_from"), Wij("Wij");
        //        T.real_imag(cc.in(arcs), dd.in(arcs));
        //        Y.real_imag(g.in(arcs), b.in(arcs));
        //        Ych.set_imag(ch_half.in(arcs));
        //
        //
        //        L_from.set_real(lij.in(arcs));
        //        Wij.real_imag(R_Wij.in_pairs(arcs), Im_Wij.in_pairs(arcs));
        //        var<Cpx> Sij("Sij"), Sji("Sji");
        //        Sij.real_imag(Pf_from.in(arcs), Qf_from.in(arcs));
        //        Sji.real_imag(Pf_to.in(arcs), Qf_to.in(arcs));
        //
        //
        //        Constraint<Cpx> I_from("I_from");
        //        I_from=(Y+Ych)*(conj(Y)+conj(Ych))*Wii.from(arcs)-T*Y*(conj(Y)+conj(Ych))*conj(Wij)-conj(T)*conj(Y)*(Y+Ych)*Wij+pow(tr,2)*Y*conj(Y)*Wii.to(arcs);
        //        SOCP.add_real(I_from.in(arcs)==pow(tr,2)*L_from);
        //
        //        var<Cpx> L_to("L_to");
        //        L_to.set_real(lji.in(arcs));
        //
        //        Constraint<Cpx> I_to("I_to");
        //        I_to=pow(tr,2)*(Y+Ych)*(conj(Y)+conj(Ych))*Wii.to(arcs)-conj(T)*Y*(conj(Y)+conj(Ych))*Wij-T*conj(Y)*(Y+Ych)*conj(Wij)+Y*conj(Y)*Wii.from(arcs);
        //        SOCP.add_real(I_to.in(arcs)==pow(tr,2)*L_to);
        //
        //        Constraint<> I_from_Pf("I_from_Pf");
        //        I_from_Pf=lij*Wii.from(arcs)-pow(tr,2)*(pow(Pf_from,2) + pow(Qf_from,2));
        //        SOCP.add(I_from_Pf.in(arcs)>=0);
        //        SOCP.get_constraint("I_from_Pf")->_relaxed = true;
        
        
        indices nonzero_arcs("nonzero_arcs");
        nonzero_arcs.add("40,25,37", "4,2,30" ,"13,6,31", "19,10,32", "36,22,35", "9,5,6", "28,16,24", "18,10,13", "38,23,36", "14,7,8", "1,1,39", "16,9,39", "17,10,11", "11,6,7", "29,17,18");
        
        
        indices arcs1("arcs1");
        arcs1.add("0,1,4", "3,3,6");
        
        indices arcs2("arcs2");
        arcs2.add("2,5,6");
        
        indices arcs3("arcs3");
        arcs3.add("5,7,8", "6,2,8", "7,8,9");
        
        indices bus_pairs1("bus_pairs1");
        bus_pairs1.add("1,4", "4,5", "5,6");
        
        indices bus_pairs2("bus_pairs2");
        bus_pairs2.add("8,9", "5,6","3,6"); //"5,6", "3,6", "6,7", "7,8", "2,8"
        
        
        
        if (current_partition_on_off_automated){
            
            // ********************* THIS PART IS FOR LIFT & PARTITION *********************
            /* Set the number of partitions (default is 1)*/
            Pf_to._num_partns = 3;
            Qf_to._num_partns = 3;
            Wii._num_partns = 2;
            lji._num_partns = 2;
            
            Constraint<> I_to_Pf_EQ("I_to_Pf_EQ");
            I_to_Pf_EQ = lji.in(arcs)*Wii.to(arcs)-(pow(Pf_to.in(arcs),2) + pow(Qf_to.in(arcs), 2));
            auto I_to_Pf_EQ_Standard = SOCP.get_standard_SOC(I_to_Pf_EQ);
            SOCP.add(I_to_Pf_EQ_Standard.in(arcs)==0, true, "lambda_II");

            
            // ********************* THIS PART IS FOR SOC_PARTITION FUNCTION *********************
//            Constraint<> I_to_Pf_temp("I_to_Pf_temp");
//            I_to_Pf_temp = lji.in(arcs)*Wii.to(arcs)-(pow(Pf_to.in(arcs),2) + pow(Qf_to.in(arcs), 2));
//            I_to_Pf_temp.in(arcs) >= 0;
//
//            //trial use SOC_partition
//            SOCP.SOC_partition(I_to_Pf_temp,20,20,true);
////            SOCP.SOC_partition(I_to_Pf_temp,12,12,false);
            
            
        }
    }
    
    
    
//    solver<> SOCPOPF(SOCP,ipopt);
//    double solver_time_start = get_wall_time();
//
//    SOCPOPF.run(output = 5, tol = 1e-7, "ma27");
//    SOCP.print_solution();
//    SOCP.print_constraints_stats(tol);
//    SOCP.print_nonzero_constraints(tol,true);
//    auto lower_bound = SOCP.get_obj_val()*upperbound;
//
//    gap = 100*(upperbound - lower_bound)/upperbound;
//    solver_time_end = get_wall_time();
//    total_time_end = get_wall_time();
//    solve_time = solver_time_end - solver_time_start;
//    total_time = total_time_end - total_time_start;
//    string out = "\nDATA_OPF, " + grid._name + ", " + to_string(nb_buses) + ", " + to_string(nb_lines) +", " + to_string(lower_bound) + ", " + to_string(-numeric_limits<double>::infinity()) + ", " + to_string(solve_time) + ", LocalOptimal, " + to_string(total_time);
//    DebugOn(out <<endl);
//    DebugOn("Final Gap = " << to_string(gap) << "%."<<endl);
//    DebugOn("Upper bound = " << to_string(upperbound) << "."<<endl);
//    DebugOn("Lower bound = " << to_string(lower_bound) << "."<<endl);
    
    /***************** OUTER APPROXIMATION BEFORE RUN *****************/
//    auto SOCPI=SOCP.build_model_interior();
//    solver<> SOCPINT(SOCPI,ipopt);
//    SOCPINT.run(output = 5, tol = 1e-7, "ma27");
//
//    vector<double> xint(SOCPI._nb_vars);
//
//    SOCPI.get_solution(xint);
    
//    DebugOn("THIS IS THE ORIGINAL FORMULATION" << endl);
//    SOCP.print();
//    auto SOCPOA = SOCP.buildOA(4);
//    DebugOn("THIS IS THE OA FORMULATION" << endl);
//    SOCPOA->print();
    
    /***************** OUTER APPROXIMATION DONE *****************/
    /***************** IF YOU WANT TO OMIT OUTER APPROXIMATION CHANGE THE MODEL IN THE SOLVER TO SOCP *****************/
    /* Solver selection */
    solver<> SOCOPF_CPX(SOCP, cplex);
    auto solver_time_start = get_wall_time();
    
    /** use the following line if you want to relax the integer variables **/
    //        SOCOPF_CPX.run(true);
    SOCOPF_CPX.run(output,tol = 1e-6);
    solver_time_end = get_wall_time();
    total_time_end = get_wall_time();
    solve_time = solver_time_end - solver_time_start;
    total_time = total_time_end - total_time_start;
    
    /* Solver selection */
    //    solver<> ACOPF_IPOPT(ACOPF, ipopt);
    //    ACOPF_IPOPT.run(output, tol = 1e-6);
    
    
    
    auto out = "DATA_OPF, " + grid._name + ", # of Buses:" + to_string(nb_buses) + ", # of Lines:" + to_string(nb_lines) +", Objective:" + to_string_with_precision(SOCP.get_obj_val(),10) + ", Upper bound:" + to_string(upperbound) + ", Solve time:" + to_string(solve_time) + ", Total time: " + to_string(total_time);
    DebugOn(out <<endl);
    
    //    double gap = 100*(ACOPF.get_obj_val() - SOCP.get_obj_val())/ACOPF.get_obj_val();
    
    gap = 100*(upperbound - SOCP.get_obj_val())/upperbound;
    DebugOn("Final Gap = " << to_string(gap) << "%."<<endl);
    
    auto nonzero_idx2 = SOCP.sorted_nonzero_constraint_indices(tol, true, "I_to_Pf");
    nonzero_idx2.print();
    
    
    return 0;
    
}

