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
    
    //    Specify the use of partitioning scheme
    bool do_partition = true;
    bool do_Model_III = false;
    string model_type = "Model_II"; //the default relaxation model is Model_II
    
    //    Switch the data file to another instance
    string fname = string(prj_dir)+"/data_sets/Power/nesta_case9_bgm__nco_tree.m";
    //    string fname = string(prj_dir)+"/data_sets/Power/nesta_case5_pjm.m";
    
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
    
    /* Initialize variables */
    R_Wij.initialize_all(1.0);
    Wii.initialize_all(1.001);
    
    /************** Add the lifted variables *****************/
//     If the partition scheme is on
    if (do_partition){
        var<> WijWji("WijWji");
        SOCP.add(WijWji.in(bus_pairs));
        
        var<> R_WijWij("R_WijWij",pos_);
        SOCP.add(R_WijWij.in(bus_pairs));
        
        var<> Im_WijWij("Im_WijWij",pos_);
        SOCP.add(Im_WijWij.in(bus_pairs));
        
        if (do_Model_III){
            model_type = "Model_III";
        }
        
        /* Second-order cone constraints */
        Constraint<> SOC("SOC");
        SOC = R_WijWij + Im_WijWij - WijWji;
        SOCP.add(SOC.in(bus_pairs) == 0);
        
        //collect the bounds on variables
        
        //bounds for Wii
        auto LB_Wii = grid.w_min;
        auto UB_Wii = grid.w_max;
    
        auto LB_R_Wij = grid.wr_min;
        auto UB_R_Wij = grid.wr_max;
        
        auto LB_Im_Wij = grid.wi_min;
        auto UB_Im_Wij = grid.wi_max;
        
        
        // define the number of partitions for variables
        int num_partitions1 = 7; //number of partitions for Wii(from)
        int num_partitions2 = 7; //number of partitions for Wii(to)
        int num_partitions3 = 7; //number of partitions for R_Wi
        
        /************** THIS SHOULD BE AN EVEN NUMBER FOR BETTER ACCURACY ***************/
        int num_partitions4 = 7; //number of partitions for Im_Wij
        
        // allocate the partition bound arrays
        vector<double> p1(num_partitions1+1); //bounds for Wii(from)
        vector<double> p2(num_partitions2+1); //bounds for Wii(to)
        vector<double> p3(num_partitions3+1); //bounds for R_Wij
        vector<double> p4(num_partitions4+1); //bounds for Im_Wij
        
        //parsing related items
        string delimiter = ","; //delimiter for correcly seperating the keys
        int fromIDX; //from index of the key
        int toIDX; //to index of the key
        string myString; //temporary string
        size_t pos; //position of the delimiter
        size_t delimiter_lenght = delimiter.length();
        
        
        for (int i=0; i<bus_pairs.size(); ++i) {
//        for (int i=0; i<1; ++i) {
            myString = bus_pairs._keys->at(i);
            pos = bus_pairs._keys->at(i).find(delimiter);
            
            fromIDX = stoi(myString.substr(0, pos)) - 1 ;
            toIDX = stoi(myString.substr(pos+delimiter_lenght)) - 1;
            
           
            //fill the partition bounds for the variables
            for (int j=0; j<num_partitions1+1; ++j) {
                p1[j] = (num_partitions1-j)*LB_Wii.eval(fromIDX)+(j)*UB_Wii.eval(fromIDX);
            }
            transform(p1.begin(), p1.end(), p1.begin(), bind(divides<double>(), placeholders::_1, num_partitions1));
          
            for (int j=0; j<num_partitions2+1; ++j) {
                p2[j] = (num_partitions2-j)*LB_Wii.eval(toIDX)+(j)*UB_Wii.eval(toIDX);
            }
            transform(p2.begin(), p2.end(), p2.begin(), bind(divides<double>(), placeholders::_1, num_partitions2));
            
            for (int j=0; j<num_partitions3+1; ++j) {
                p3[j] = (num_partitions3-j)*LB_R_Wij.eval(i)+(j)*UB_R_Wij.eval(i);
            }
            transform(p3.begin(), p3.end(), p3.begin(), bind(divides<double>(), placeholders::_1, num_partitions3));
           
            for (int j=0; j<num_partitions4+1; ++j) {
                p4[j] = (num_partitions4-j)*LB_Im_Wij.eval(i)+(j)*UB_Im_Wij.eval(i);
            }
            transform(p4.begin(), p4.end(), p4.begin(), bind(divides<double>(), placeholders::_1, num_partitions4));
            
            //add the partitions&relaxation on the variables
            SOCP.partition("WijWji" + to_string(i), model_type, WijWji(i), Wii(fromIDX), Wii(toIDX),p1,p2);
            SOCP.partition("R_WijWij" + to_string(i), model_type, R_WijWij(i), R_Wij(i), R_Wij(i),p3,p3);
            SOCP.partition("Im_WijWij" + to_string(i), model_type, Im_WijWij(i), Im_Wij(i), Im_Wij(i),p4,p4);
        }
    }
    
    
    /**  Objective */
    auto obj = product(c1,Pg) + product(c2,pow(Pg,2)) + sum(c0);
    SOCP.min(obj);
    ACOPF.min(obj);
    
    /** Constraints */
    
    if (!do_partition){
    /* Second-order cone constraints */
    Constraint<> SOC("SOC");
    SOC = pow(R_Wij, 2) + pow(Im_Wij, 2) - Wii.from(bus_pairs)*Wii.to(bus_pairs);
    SOCP.add(SOC.in(bus_pairs) <= 0);
    }
    
    /* Equality of Second-order cone (for upperbound) */
    Constraint<> Equality_SOC("Equality_SOC");
    Equality_SOC = pow(R_Wij, 2) + pow(Im_Wij, 2) - Wii.from(bus_pairs)*Wii.to(bus_pairs);
    ACOPF.add(Equality_SOC.in(bus_pairs) == 0);
    
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
    
  
    SOCP.print();
    
    /* Solver selection */
        solver<> SOCOPF_CPX(SOCP, cplex);
        auto solver_time_start = get_wall_time();
        SOCOPF_CPX.run(output, tol = 1e-6);
        solver_time_end = get_wall_time();
        total_time_end = get_wall_time();
        solve_time = solver_time_end - solver_time_start;
        total_time = total_time_end - total_time_start;
    
//        SOCP.print_solution();
    
    /* Solver selection */
    solver<> ACOPF_IPOPT(ACOPF, ipopt);
    ACOPF_IPOPT.run(output, tol = 1e-6);

//        ACOPF.print_solution();

    
    string out = "DATA_OPF, " + grid._name + ", # of Buses:" + to_string(nb_buses) + ", # of Lines:" + to_string(nb_lines) +", Objective:" + to_string_with_precision(SOCP.get_obj_val(),10) + ", Upper bound:" + to_string(ACOPF.get_obj_val()) + ", Solve time:" + to_string(solve_time) + ", Total time: " + to_string(total_time);
    DebugOn(out <<endl);
    
    double gap = 100*(ACOPF.get_obj_val() - SOCP.get_obj_val())/ACOPF.get_obj_val();
    DebugOn("Final Gap = " << to_string(gap) << "%."<<endl);
    
    return 0;
    
    
}
