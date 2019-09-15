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
    string fname = string(prj_dir)+"/data_sets/Power/nesta_case5_pjm.m";
  //  string fname = string(prj_dir)+"/data_sets/Power/nesta_case9_bgm__nco_tree.m";
//    string fname = string(prj_dir)+"/data_sets/Power/nesta_case39_1_bgm__nco.m";
//    string fname="/Users/smitha/Desktop/nesta-0.7.0/opf/nco/nesta_case9_tree.m";
    
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
    
    double upperbound = grid.solve_acopf(ACRECT);
    
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
    auto pg_max_sq = grid.pg_max_sq.in(gens);
    auto pg_min_sq = grid.pg_min_sq.in(gens);
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
    auto pl = grid.pl.in(nodes);
    auto ql = grid.ql.in(nodes);
    auto gs = grid.gs.in(nodes);
    auto bs = grid.bs.in(nodes);
    auto g_ft = grid.g_ft.in(arcs);
    auto g_ff = grid.g_ff.in(arcs);
    auto g_tt = grid.g_tt.in(arcs);
    auto g_tf = grid.g_tf.in(arcs);
    auto b_ft = grid.b_ft.in(arcs);
    auto b_ff = grid.b_ff.in(arcs);
    auto b_tf = grid.b_tf.in(arcs);
    auto b_tt = grid.b_tt.in(arcs);
    auto S_max = grid.S_max.in(arcs);
    auto v_max = grid.v_max.in(nodes);
    auto v_min = grid.v_min.in(nodes);
    auto w_max = grid.w_max.in(nodes);
    auto w_min = grid.w_min.in(nodes);
    auto tan_th_min = grid.tan_th_min.in(bus_pairs);
    auto tan_th_max = grid.tan_th_max.in(bus_pairs);
    
    
    /** MODEL DECLARATION */
   // Model<> SOCP("SCOPF Model"); //model for the relaxation
    auto SOCP = make_shared<Model<>>("SOCP Model");
    Model<> ACOPF("ACOPF Model"); //exact model
    
    /** Variables */
    /* power generation variables */
    var<> Pg("Pg", pg_min, pg_max);
    var<> Qg ("Qg", qg_min, qg_max);
    var<> etag("etag", pg_min_sq, pg_max_sq);
    SOCP->add(Pg.in(gens),Qg.in(gens));
    SOCP->add(etag.in(gens));
    //ACOPF.add(Pg.in(gens),Qg.in(gens));
    
    /* power flow variables */
    var<> Pf_from("Pf_from", -1*grid.S_max,grid.S_max);
    var<> Qf_from("Qf_from", -1*grid.S_max,grid.S_max);
    var<> Pf_to("Pf_to", -1*grid.S_max,grid.S_max);
    var<> Qf_to("Qf_to", -1*grid.S_max,grid.S_max);
    SOCP->add(Pf_from.in(arcs), Qf_from.in(arcs), Pf_to.in(arcs), Qf_to.in(arcs));
    //ACOPF.add(Pf_from.in(arcs), Qf_from.in(arcs), Pf_to.in(arcs), Qf_to.in(arcs));
    
    /* Real part of Wij = ViVj */
    var<>  R_Wij("R_Wij", grid.wr_min, grid.wr_max);
    /* Imaginary part of Wij = ViVj */
    var<>  Im_Wij("Im_Wij", grid.wi_min, grid.wi_max);
    /* Magnitude of Wii = Vi^2 */
    var<>  Wii("Wii", grid.w_min, grid.w_max);
    SOCP->add(Wii.in(nodes));
    SOCP->add(R_Wij.in(bus_pairs));
    SOCP->add(Im_Wij.in(bus_pairs));
    //ACOPF.add(Wii.in(nodes));
    //ACOPF.add(R_Wij.in(bus_pairs));
    //ACOPF.add(Im_Wij.in(bus_pairs));
    
    var<> lij("lij", lij_min,lij_max);
    var<> lji("lji", lji_min,lji_max);
    if(current){
        SOCP->add(lij.in(arcs));
        SOCP->add(lji.in(arcs));
    }
    
    /* Initialize variables */
    R_Wij.initialize_all(1.0);
    //lij.initialize_all(1.0);
    Wii.initialize_all(1.001);
    
    Constraint<> obj_cost("obj_cost");
    obj_cost=etag-pow(Pg,2);
    SOCP->add(obj_cost.in(gens)>=0);
    
    /************** Add the lifted variables *****************/
    if(current){
        
        param<Cpx> T("T"), Y("Y"), Ych("Ych");
        var<Cpx> L_from("L_from"), Wij("Wij");
        T.real_imag(cc.in(arcs), dd.in(arcs));
        Y.real_imag(g.in(arcs), b.in(arcs));
        Ych.set_imag(ch_half.in(arcs));
        
        
        L_from.set_real(lij.in(arcs));
        Wij.real_imag(R_Wij.in_pairs(arcs), Im_Wij.in_pairs(arcs));
//        var<Cpx> Sij("Sij"), Sji("Sji");
//        Sij.real_imag(Pf_from.in(arcs), Qf_from.in(arcs));
//        Sji.real_imag(Pf_to.in(arcs), Qf_to.in(arcs));
        
        
        Constraint<Cpx> I_from("I_from");
        I_from=(Y+Ych)*(conj(Y)+conj(Ych))*Wii.from(arcs)-T*Y*(conj(Y)+conj(Ych))*conj(Wij)-conj(T)*conj(Y)*(Y+Ych)*Wij+pow(tr,2)*Y*conj(Y)*Wii.to(arcs)-pow(tr,2)*L_from;
        SOCP->add_real(I_from.in(arcs)==0);
        
        var<Cpx> L_to("L_to");
        L_to.set_real(lji.in(arcs));

        Constraint<Cpx> I_to("I_to");
        I_to=pow(tr,2)*(Y+Ych)*(conj(Y)+conj(Ych))*Wii.to(arcs)-conj(T)*Y*(conj(Y)+conj(Ych))*Wij-T*conj(Y)*(Y+Ych)*conj(Wij)+Y*conj(Y)*Wii.from(arcs)-pow(tr,2)*L_to;
        SOCP->add_real(I_to.in(arcs)==0);
//
        Constraint<> I_from_Pf("I_from_Pf");
        I_from_Pf=lij*Wii.from(arcs)-pow(tr,2)*(pow(Pf_from,2) + pow(Qf_from,2));
        SOCP->add(I_from_Pf.in(arcs)==0, true);
       
            //   SOCP->get_constraint("I_from_Pf")->_relaxed = true;
        
        Constraint<> I_to_Pf("I_to_Pf");
        I_to_Pf=lji.in(arcs)*Wii.to(arcs)-(pow(Pf_to.in(arcs),2) + pow(Qf_to.in(arcs), 2));
        SOCP->add(I_to_Pf.in(arcs)>=0);
    //    SOCP->get_constraint("I_to_Pf")->_relaxed = true;
        
        /* Second-order cone */
        Constraint<> SOC("SOC");
        SOC = pow(R_Wij, 2) + pow(Im_Wij, 2) - Wii.from(bus_pairs)*Wii.to(bus_pairs);
        SOCP->add(SOC.in(bus_pairs) == 0, true);
        
    }
    
    
    
    
    /**  Objective */
    auto obj = product(c1,Pg) + product(c2,etag) + sum(c0);
    SOCP->min(obj);
    //ACOPF.min(obj);
    
    /** Constraints */
    
    
    /* Equality of Second-order cone (for upperbound) */
    Constraint<> Equality_SOC("Equality_SOC");
    Equality_SOC = pow(R_Wij, 2) + pow(Im_Wij, 2) - Wii.from(bus_pairs)*Wii.to(bus_pairs);
    ////ACOPF.add(Equality_SOC.in(bus_pairs) == 0);
    
    /* Flow conservation */
    Constraint<> KCL_P("KCL_P");
    KCL_P  = sum(Pf_from, out_arcs) + sum(Pf_to, in_arcs) + pl - sum(Pg, gen_nodes) + gs*Wii;
    SOCP->add(KCL_P.in(nodes) == 0);
    //ACOPF.add(KCL_P.in(nodes) == 0);
    
    Constraint<> KCL_Q("KCL_Q");
    KCL_Q  = sum(Qf_from, out_arcs) + sum(Qf_to, in_arcs) + ql - sum(Qg, gen_nodes) - bs*Wii;
    SOCP->add(KCL_Q.in(nodes) == 0);
    //ACOPF.add(KCL_Q.in(nodes) == 0);
    
    /* AC Power Flow */
    Constraint<> Flow_P_From("Flow_P_From");
    Flow_P_From = Pf_from - (g_ff*Wii.from(arcs) + g_ft*R_Wij.in_pairs(arcs) + b_ft*Im_Wij.in_pairs(arcs));
    SOCP->add(Flow_P_From.in(arcs) == 0);
    //ACOPF.add(Flow_P_From.in(arcs) == 0);
    
    Constraint<> Flow_P_To("Flow_P_To");
    Flow_P_To = Pf_to - (g_tt*Wii.to(arcs) + g_tf*R_Wij.in_pairs(arcs) - b_tf*Im_Wij.in_pairs(arcs));
    SOCP->add(Flow_P_To.in(arcs) == 0);
    //ACOPF.add(Flow_P_To.in(arcs) == 0);
    
    Constraint<> Flow_Q_From("Flow_Q_From");
    Flow_Q_From = Qf_from - (g_ft*Im_Wij.in_pairs(arcs) - b_ff*Wii.from(arcs) - b_ft*R_Wij.in_pairs(arcs));
    SOCP->add(Flow_Q_From.in(arcs) == 0);
    //ACOPF.add(Flow_Q_From.in(arcs) == 0);
    
    Constraint<> Flow_Q_To("Flow_Q_To");
    Flow_Q_To = Qf_to + b_tt*Wii.to(arcs) + b_tf*R_Wij.in_pairs(arcs) + g_tf*Im_Wij.in_pairs(arcs);
    SOCP->add(Flow_Q_To.in(arcs) == 0);
    //ACOPF.add(Flow_Q_To.in(arcs) == 0);

    
    
    /* Phase Angle Bounds constraints */
    Constraint<> PAD_UB("PAD_UB");
    PAD_UB = Im_Wij.in(bus_pairs);
    PAD_UB <= tan_th_max*R_Wij.in(bus_pairs);
        SOCP->add(PAD_UB.in(bus_pairs));
        //ACOPF.add(PAD_UB.in(bus_pairs));
    
    Constraint<> PAD_LB("PAD_LB");
    PAD_LB =  Im_Wij.in(bus_pairs);
    PAD_LB >= tan_th_min*R_Wij.in(bus_pairs);
    SOCP->add(PAD_LB.in(bus_pairs));
    //ACOPF.add(PAD_LB.in(bus_pairs));

    
    
    /* Thermal Limit Constraints */
    Constraint<> Thermal_Limit_from("Thermal_Limit_from");
    Thermal_Limit_from = pow(Pf_from, 2) + pow(Qf_from, 2);
    Thermal_Limit_from <= pow(S_max,2);
    //SDP.add(Thermal_Limit_from.in(arcs));
    SOCP->add(Thermal_Limit_from.in(arcs), true);
    //ACOPF.add(Thermal_Limit_from.in(arcs));
    
    
    
    Constraint<> Thermal_Limit_to("Thermal_Limit_to");
    Thermal_Limit_to = pow(Pf_to, 2) + pow(Qf_to, 2);
    Thermal_Limit_to <= pow(S_max,2);
    SOCP->add(Thermal_Limit_to.in(arcs), true);
    //ACOPF.add(Thermal_Limit_to.in(arcs));

    
    func<> theta_L = atan(min(Im_Wij.get_lb().in(bus_pairs)/R_Wij.get_ub().in(bus_pairs),Im_Wij.get_lb().in(bus_pairs)/R_Wij.get_lb().in(bus_pairs)));
    func<> theta_U = atan(max(Im_Wij.get_ub().in(bus_pairs)/R_Wij.get_lb().in(bus_pairs),Im_Wij.get_ub().in(bus_pairs)/R_Wij.get_ub().in(bus_pairs)));
    func<> phi=(theta_U.in(bus_pairs)+theta_L.in(bus_pairs))/2.0;
    func<> del=(theta_U.in(bus_pairs)-theta_L.in(bus_pairs))/2.0;
    
    
    Constraint<> LNC1("LNC1");
    LNC1 += (sqrt(Wii.get_lb().from(bus_pairs))+sqrt(Wii.get_ub().from(bus_pairs)))*(sqrt(Wii.get_lb().to(bus_pairs))+sqrt(Wii.get_ub().to(bus_pairs)))*(Im_Wij.in(bus_pairs)*sin(phi.in(bus_pairs)) + R_Wij.in(bus_pairs)*cos(phi.in(bus_pairs)));
    
    LNC1 -=sqrt(Wii.get_ub().to(bus_pairs))*cos(del.in(bus_pairs))*(sqrt(Wii.get_lb().to(bus_pairs))+sqrt(Wii.get_ub().to(bus_pairs)))*Wii.from(bus_pairs);
    
    LNC1 -=sqrt(Wii.get_ub().from(bus_pairs))*cos(del.in(bus_pairs))*(sqrt(Wii.get_lb().from(bus_pairs))+sqrt(Wii.get_ub().from(bus_pairs)))*Wii.to(bus_pairs);
    
    LNC1-=sqrt(Wii.get_ub().from(bus_pairs))*sqrt(Wii.get_ub().to(bus_pairs))*cos(del)*(sqrt(Wii.get_lb().from(bus_pairs))*
                                                                                        sqrt(Wii.get_lb().to(bus_pairs)) - sqrt(Wii.get_ub().from(bus_pairs))*sqrt(Wii.get_ub().to(bus_pairs)));
    SOCP->add(LNC1.in(bus_pairs) >= 0);
    //ACOPF.add(LNC1.in(bus_pairs) >= 0);
    
    
    Constraint<> LNC2("LNC2");
    LNC2 += (sqrt(Wii.get_lb().from(bus_pairs))+sqrt(Wii.get_ub().from(bus_pairs)))*(sqrt(Wii.get_lb().to(bus_pairs))+sqrt(Wii.get_ub().to(bus_pairs)))*(sin(phi.in(bus_pairs))*Im_Wij.in(bus_pairs) + cos(phi.in(bus_pairs))*R_Wij.in(bus_pairs));
    LNC2 -=sqrt(Wii.get_lb().to(bus_pairs))*cos(del.in(bus_pairs))*(sqrt(Wii.get_lb().to(bus_pairs))+sqrt(Wii.get_ub().to(bus_pairs)))*Wii.from(bus_pairs);
    LNC2 -=sqrt(Wii.get_lb().from(bus_pairs))*cos(del.in(bus_pairs))*(sqrt(Wii.get_lb().from(bus_pairs))+sqrt(Wii.get_ub().from(bus_pairs)))*Wii.to(bus_pairs);
    LNC2 -=sqrt(Wii.get_lb().from(bus_pairs))*sqrt(Wii.get_lb().to(bus_pairs))*cos(del.in(bus_pairs))*(sqrt(Wii.get_ub().from(bus_pairs))*
                                                                                                       sqrt(Wii.get_ub().to(bus_pairs))-sqrt(Wii.get_lb().from(bus_pairs))*sqrt(Wii.get_lb().to(bus_pairs)));
    SOCP->add(LNC2.in(bus_pairs) >= 0);
    //ACOPF.add(LNC2.in(bus_pairs) >= 0);
    
    Constraint<> obj_UB("obj_UB");
    obj_UB=(product(c1,Pg) + product(c2,etag) + sum(c0))-upperbound;
    SOCP->add(obj_UB.in(range(0,0))<=0);
    
    /***************** CALLING OBBT BEFORE CALLING RUN **********************/
    //    SOCP.reset_constrs();
//    SOCP->print();
//    solver<> SOCOPF(SOCP, ipopt);
//    SOCOPF.run(output=5, tol = 1e-6);
    
  
    
    double max_time = 100000;
    int max_iter = 5;
    int precision = 4;
    solver<> SOCPOPF(SOCP, ipopt);
    SOCPOPF.run(output, tol = 1e-6);
    
    SOCP->run_obbt(max_time,max_iter,{true,upperbound},precision);
    auto original_SOC = grid.build_SCOPF();
    solver<> SOCOPF_ORIG(original_SOC, ipopt);
    SOCOPF_ORIG.run(output, tol = 1e-6);
    double original_LB = original_SOC->get_obj_val();
    double gap = 100*(upperbound - original_LB)/upperbound;
    DebugOn("Initial Gap = " << to_string(gap) << "%."<<endl);
    if(SOCP->_status==0||SOCP->_status==1)
    {
    gap = 100*(upperbound - SOCP->get_obj_val())/upperbound;
    DebugOn("Gap after OBBT = " << to_string(gap) << "%."<<endl);
    }
    auto nonzero_idx = SOCP->sorted_nonzero_constraint_indices(tol, true, "I_to_Pf");
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
        
        
        
        if (true){
            
            // ********************* THIS PART IS FOR LIFT & PARTITION *********************
            /* Set the number of partitions (default is 1)*/
            Pf_to._num_partns = 20;
            Qf_to._num_partns = 20;
            Wii._num_partns = 10;
            lji._num_partns = 10;
            
            Constraint<> I_to_Pf_EQ("I_to_Pf_EQ");
            I_to_Pf_EQ = lji.in(arcs)*Wii.to(arcs)-(pow(Pf_to.in(arcs),2) + pow(Qf_to.in(arcs), 2));
            auto I_to_Pf_EQ_Standard = SOCP->get_standard_SOC(I_to_Pf_EQ);
            SOCP->add(I_to_Pf_EQ_Standard.in(arcs)==0, true, "lambda_II");

            
            // ********************* THIS PART IS FOR SOC_PARTITION FUNCTION *********************
            Constraint<> I_to_Pf_temp("I_to_Pf_temp");
            I_to_Pf_temp = lji.in(arcs)*Wii.to(arcs)-(pow(Pf_to.in(arcs),2) + pow(Qf_to.in(arcs), 2));
            I_to_Pf_temp.in(arcs) >= 0;

            //trial use SOC_partition
//            SOCP.SOC_partition(I_to_Pf_temp,20,20,true);
            SOCP->SOC_partition(I_to_Pf_temp,12,12,false);
            
            
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
//    auto SOCPI=SOCP->build_model_interior();
//    SOCPI->print();
//    solver<> SOCPINT(SOCPI,ipopt);
//    SOCPINT.run(output = 5, tol = 1e-7, "ma27");
//
//    vector<double> xint(SOCPI->_nb_vars);
//
//    SOCPI->get_solution(xint);
    
    
    DebugOn("THIS IS THE ORIGINAL FORMULATION" << endl);
    SOCP->print();
    auto SOCPOA = SOCP->buildOA(10);
    DebugOn("THIS IS THE OA FORMULATION" << endl);
    SOCPOA->print();
    
    /***************** OUTER APPROXIMATION DONE *****************/
    /***************** IF YOU WANT TO OMIT OUTER APPROXIMATION CHANGE THE MODEL IN THE SOLVER TO SOCP *****************/
    /* Solver selection */
    solver<> SOCOPF_CPX(SOCPOA, ipopt);
    auto solver_time_start = get_wall_time();
    
    /** use the following line if you want to relax the integer variables **/
    //        SOCOPF_CPX.run(true);
    SOCOPF_CPX.run(output,tol = 1e-6);
    gap = 100*(upperbound - SOCPOA->get_obj_val())/upperbound;
    DebugOn("Gap after OA = " << to_string(gap) << "%."<<endl);
    solver_time_end = get_wall_time();
    total_time_end = get_wall_time();
    solve_time = solver_time_end - solver_time_start;
    total_time = total_time_end - total_time_start;
    
    /* Solver selection */
    //    solver<> //ACOPF_IPOPT(//ACOPF, ipopt);
    //    //ACOPF_IPOPT.run(output, tol = 1e-6);
    
    
    
    auto out = "DATA_OPF, " + grid._name + ", # of Buses:" + to_string(nb_buses) + ", # of Lines:" + to_string(nb_lines) +", Objective:" + to_string_with_precision(SOCP->get_obj_val(),10) + ", Upper bound:" + to_string(upperbound) + ", Solve time:" + to_string(solve_time) + ", Total time: " + to_string(total_time);
    DebugOn(out <<endl);
    
    //    double gap = 100*(//ACOPF.get_obj_val() - SOCP.get_obj_val())///ACOPF.get_obj_val();
    
    gap = 100*(upperbound - SOCP->get_obj_val())/upperbound;
    DebugOn("Final Gap = " << to_string(gap) << "%."<<endl);
    
    auto nonzero_idx2 = SOCP->sorted_nonzero_constraint_indices(tol, true, "I_to_Pf");
    nonzero_idx2.print();
    
    
    return 0;
    
}

