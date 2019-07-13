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
    bool current_partition_lambda = false;
    bool current_partition_on_off = true;
    
    //    Specify the use of partitioning scheme without current
    bool do_partition = false;
    bool do_Model_III = false;
    string model_type = "Model_II"; //the default relaxation model is Model_II
    
    //    Switch the data file to another instance
    string fname = string(prj_dir)+"/data_sets/Power/pglib_opf_case3_lmbd__api.m";
    //   string fname = string(prj_dir)+"/data_sets/Power/nesta_case9_bgm__nco_tree.m";
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
        Constraint<> SOC("SOC_original");
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
        int num_partitions1 = 2; //number of partitions for Wii(from)
        int num_partitions2 = 2; //number of partitions for Wii(to)
        int num_partitions3 = 10; //number of partitions for R_Wi
        
        /************** THIS SHOULD BE AN EVEN NUMBER FOR BETTER ACCURACY ***************/
        int num_partitions4 = 10; //number of partitions for Im_Wij
        
        // allocate the partition bound arrays
        vector<double> p1(num_partitions1+1); //bounds for Wii(from)
        vector<double> p2(num_partitions2+1); //bounds for Wii(to)
        vector<double> p3(num_partitions3+1); //bounds for R_Wij
        vector<double> p4(num_partitions4+1); //bounds for Im_Wij
        
        //parsing related items
        string delimiter = ","; //delimiter for correcly seperating the keys
        string fromIDX; //from index of the key
        string toIDX; //to index of the key
        string myString; //temporary string
        size_t pos; //position of the delimiter
        size_t delimiter_lenght = delimiter.length();
        
        vector<int> constraint_idx = {0,3,6};
        
        //        for (int i=0; i<bus_pairs.size(); ++i) {
        for (int j=0; j<constraint_idx.size(); ++j) {
            int i = constraint_idx[j];
            myString = bus_pairs._keys->at(i);
            pos = bus_pairs._keys->at(i).find(delimiter);
            
            fromIDX = myString.substr(0, pos);
            toIDX = myString.substr(pos+delimiter_lenght);
            
            
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
        SOCP.add(I_from_Pf.in(arcs)>=0);
//        SOCP.get_constraint("I_from_Pf")->_relaxed = true;
        //        SOCP.add(I_from_Pf.in(arcs)==0, true);
        
        Constraint<> I_to_Pf("I_to_Pf");
        I_to_Pf=lji*Wii.to(arcs)-(pow(Pf_to,2) + pow(Qf_to, 2));
        SOCP.add(I_to_Pf.in(arcs)>=0);
        SOCP.get_constraint("I_to_Pf")->_relaxed = true;
        //        SOCP.add(I_to_Pf.in(arcs)==0, true);
        
        
        if (current_partition_lambda){
            if (do_Model_III){
                model_type = "Model_III";
            }
            
            var<> Pf_to_squared("Pf_to_squared",pos_);
            SOCP.add(Pf_to_squared.in(arcs));
            
            var<> Qf_to_squared("Qf_to_squared",pos_);
            SOCP.add(Qf_to_squared.in(arcs));
            
            var<> ljiWii_to("ljiWii_to");
            SOCP.add(ljiWii_to.in(arcs));
            
            Constraint<> I_to_Pf_EQ("I_to_Pf_EQ");
            I_to_Pf_EQ=ljiWii_to-(Pf_to_squared + Qf_to_squared);
            SOCP.add(I_to_Pf_EQ.in(arcs)==0);
            
            //collect the bounds on variables
            
            //bounds for Pf_to and Qf_to
            auto UB_Pf_to = grid.S_max; // LB is the negative of UB
            auto UB_Qf_to = grid.S_max; // LB is the negative of UB
            
            //bounds for Wii
            auto LB_Wii = grid.w_min;
            auto UB_Wii = grid.w_max;
            
            // define the number of partitions for variables
            /************** THESE SHOULD BE AN EVEN NUMBER FOR BETTER ACCURACY ***************/
            int num_partitions1 = 20; //number of partitions for Pf_to
            int num_partitions2 = 20; //number of partitions for Qf_to
            
            int num_partitions3 = 4; //number of partitions for Wii(to)
            int num_partitions4 = 4; //number of partitions for lji
            
            // allocate the partition bound arrays
            vector<double> p1(num_partitions1+1); //bounds for Pf_to
            vector<double> p2(num_partitions2+1); //bounds for Qf_to
            vector<double> p3(num_partitions3+1); //bounds for Wii(to)
            vector<double> p4(num_partitions4+1); //bounds for lji
            
            //parsing related items
            string delimiter = ","; //delimiter for correcly seperating the keys
            string toIDX; //to index of the key
            string myString; //temporary string
            size_t pos; //position of the delimiter
            size_t delimiter_lenght = delimiter.length();
            
            vector<int> constraint_idx = {0,1,2,3,4,5,6,7};
            
            //            for (int i=0; i<arcs.size(); ++i) {
            for (int k=0; k<constraint_idx.size(); ++k) {
                int i = constraint_idx[k];
                myString = bus_pairs._keys->at(i);
                pos = bus_pairs._keys->at(i).find(delimiter);
                toIDX = myString.substr(pos+delimiter_lenght);
                
                //fill the partition bounds for the variables
                for (int j=0; j<num_partitions1+1; ++j) {
                    p1[j] = (num_partitions1-j)*(-UB_Pf_to.eval(i))+(j)*UB_Pf_to.eval(i);
                }
                transform(p1.begin(), p1.end(), p1.begin(), bind(divides<double>(), placeholders::_1, num_partitions1));
                
                for (int j=0; j<num_partitions2+1; ++j) {
                    p2[j] = (num_partitions2-j)*(-UB_Qf_to.eval(i))+(j)*UB_Qf_to.eval(i);
                }
                transform(p2.begin(), p2.end(), p2.begin(), bind(divides<double>(), placeholders::_1, num_partitions2));
                
                for (int j=0; j<num_partitions3+1; ++j) {
                    p3[j] = (num_partitions3-j)*LB_Wii.eval(toIDX)+(j)*UB_Wii.eval(toIDX);
                }
                transform(p3.begin(), p3.end(), p3.begin(), bind(divides<double>(), placeholders::_1, num_partitions3));
                
                for (int j=0; j<num_partitions4+1; ++j) {
                    p4[j] = (num_partitions4-j)*lji_min.eval(i)+(j)*lji_max.eval(i);
                }
                transform(p4.begin(), p4.end(), p4.begin(), bind(divides<double>(), placeholders::_1, num_partitions4));
                
                //add the partitions&relaxation on the variables
                SOCP.partition("Pf_to_squared" + to_string(i), model_type, Pf_to_squared(i), Pf_to(i), Pf_to(i),p1,p1);
                SOCP.partition("Qf_to_squared" + to_string(i), model_type, Qf_to_squared(i), Qf_to(i), Qf_to(i),p2,p2);
                SOCP.partition("ljiWii_to" + to_string(i), model_type, ljiWii_to(i),  Wii(toIDX), lji(i),p3,p4);
            }
        }
        
        if (current_partition_on_off){
            
            var<> Pf_to_squared("Pf_to_squared", 0, grid.S_max*grid.S_max);
            SOCP.add(Pf_to_squared.in(arcs));
            Pf_to_squared._lift = true;
            
            var<> Qf_to_squared("Qf_to_squared", 0, grid.S_max*grid.S_max);
            SOCP.add(Qf_to_squared.in(arcs));
            Qf_to_squared._lift = true;
            
            /*need to provide bounds for the variables,
             have a scheme to provide bounds for the bilinear case*/
            auto Wii_to = Wii.to(arcs);
            auto id_set = indices("Wii_to,Arcs");
            id_set = combine(*Wii_to._indices, *lji._indices);
            var<> ljiWii_to("ljiWii_to",0,lji_max*grid.w_max.to(arcs));
            SOCP.add(ljiWii_to.in(id_set));
            ljiWii_to._lift = true;
            
            Constraint<> I_to_Pf_EQ("I_to_Pf_EQ");
            I_to_Pf_EQ=ljiWii_to-(Pf_to_squared + Qf_to_squared);
            SOCP.add(I_to_Pf_EQ.in(arcs)==0);
            
            //collect the bounds on variables
            
            //bounds for Pf_to and Qf_to
            auto UB_Pf_to = grid.S_max; // LB is the negative of UB
            auto UB_Qf_to = grid.S_max; // LB is the negative of UB
            
            //bounds for Wii
            auto LB_Wii = grid.w_min;
            auto UB_Wii = grid.w_max;
            
            // define the number of partitions for variables
            /************** THESE SHOULD BE AN EVEN NUMBER FOR BETTER ACCURACY ***************/
            int num_partitions1 = 50; //number of partitions for Pf_to
            int num_partitions2 = 50; //number of partitions for Qf_to
            
            int num_partitions3 = 8; //number of partitions for Wii(to)
            int num_partitions4 = 8; //number of partitions for lji
            
            
            /* create an index set for all z and unify them maybe later */
            var<int> z1("z1",0,1);
            indices partns1("partns1");
            partns1 = indices(range(1,num_partitions1));
            auto var_indices1 = *Pf_to._indices;
            auto inst_partition1 = indices(var_indices1,partns1);
            SOCP.add(z1.in(inst_partition1));
            
            var<int> z2("z2",0,1);
            indices partns2("partns2");
            partns2 = indices(range(1,num_partitions2));
            auto inst_partition2 = indices(var_indices1,partns2);
            SOCP.add(z2.in(inst_partition2));
            
            var<int> z3("z3",0,1);
            indices partns3("partns3");
            partns3 = indices(range(1,num_partitions3),range(1,num_partitions4));
            auto inst_partition3 = indices(var_indices1,partns3);
            SOCP.add(z3.in(inst_partition3));
            
            //parsing related items
            string delimiter = ","; //delimiter for correcly seperating the keys
            string toIDX; //to index of the key
            string myString; //temporary string
            size_t pos; //position of the delimiter
            size_t delimiter_lenght = delimiter.length();
            
            vector<int> constraint_idx = {0,1,2};
            
            //            for (int i=0; i<arcs.size(); ++i) {
            for (int k=0; k<constraint_idx.size(); ++k) {
                int i = constraint_idx[k];
                myString = bus_pairs._keys->at(i);
                pos = bus_pairs._keys->at(i).find(delimiter);
                toIDX = myString.substr(pos+delimiter_lenght);
                
                
//                add the partitions&relaxation on the variables
                auto cur_key = var_indices1._keys->at(i);
                auto Pf_to1 = Pf_to(cur_key);
                auto Pf_to_s1 = Pf_to_squared(cur_key);
                indices var_indices_temp("var_indices_temp");
                var_indices_temp.add({cur_key});
                auto z1temp = z1.in(var_indices_temp,partns1);
                
                Constraint<> z1Sum("z1Sum"+to_string(i));
                z1Sum = sum(z1temp);
                SOCP.add(z1Sum==1);
                
                auto Qf_to1 = Qf_to(cur_key);
                auto Qf_to_s1 = Qf_to_squared(cur_key);
                auto z2temp = z2.in(var_indices_temp,partns2);
                
                Constraint<> z2Sum("z2Sum"+to_string(i));
                z2Sum = sum(z2temp);
                SOCP.add(z2Sum==1);
                
                auto ljiWii_to1 = ljiWii_to(toIDX+","+cur_key);
                auto Wii_to1 = Wii(toIDX);
                auto lji1 = lji(cur_key);
                auto z3temp = z3.in(var_indices_temp,partns3);
                
                Constraint<> z3Sum("z3Sum"+to_string(i));
                z3Sum = sum(z3temp);
                SOCP.add(z3Sum==1);
                
                SOCP.add_on_off_McCormick_new("Pf_to_squared"+to_string(i), Pf_to_s1, Pf_to1, Pf_to1, z1temp, num_partitions1,num_partitions1);
                SOCP.add_on_off_McCormick_new("Qf_to_squared" + to_string(i), Qf_to_s1, Qf_to1, Qf_to1,  z2temp, num_partitions2, num_partitions2);
              SOCP.add_on_off_McCormick_new("ljiWii_to" + to_string(i), ljiWii_to1,  Wii_to1, lji1, z3temp,  num_partitions3, num_partitions4);
            }
        }
    }
    
    
    
    
    
    /**  Objective */
    auto obj = product(c1,Pg) + product(c2,pow(Pg,2)) + sum(c0);
    SOCP.min(obj);
    ACOPF.min(obj);
    
    /** Constraints */
    
    /* Equality of Second-order cone (for upperbound) */
    Constraint<> SOC("SOC");
    SOC = pow(R_Wij, 2) + pow(Im_Wij, 2) - Wii.from(bus_pairs)*Wii.to(bus_pairs);
    SOCP.add(SOC.in(bus_pairs) <= 0);
    SOCP.get_constraint("SOC")->_relaxed = true;
    
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
    
    
    /* Solver selection */
    solver<> SOCOPF_CPX(SOCP, cplex);
    auto solver_time_start = get_wall_time();
    
    /** use the following line if you want to relax the integer variables **/
//    SOCOPF_CPX.run(true);
    SOCOPF_CPX.run(output,tol = 1e-6);
    solver_time_end = get_wall_time();
    total_time_end = get_wall_time();
    solve_time = solver_time_end - solver_time_start;
    total_time = total_time_end - total_time_start;
    
    //        SOCP.print_solution();
    
    /* Solver selection */
    //    solver<> ACOPF_IPOPT(ACOPF, ipopt);
    //    ACOPF_IPOPT.run(output, tol = 1e-6);
    
    double upperbound = grid.solve_acopf(ACRECT);
    
    //    solver<> SOCOPF_CPX2(SOCP, cplex);
    //    SOCOPF_CPX2.run(output, tol = 1e-6);
    //        ACOPF.print_solution();
    
    
    string out = "DATA_OPF, " + grid._name + ", # of Buses:" + to_string(nb_buses) + ", # of Lines:" + to_string(nb_lines) +", Objective:" + to_string_with_precision(SOCP.get_obj_val(),10) + ", Upper bound:" + to_string(upperbound) + ", Solve time:" + to_string(solve_time) + ", Total time: " + to_string(total_time);
    DebugOn(out <<endl);
    
    //    double gap = 100*(ACOPF.get_obj_val() - SOCP.get_obj_val())/ACOPF.get_obj_val();
    
    
    
    SOCP.print();
    SOCP.print_solution();

    auto v = SOCP.sorted_nonzero_constraints(tol,true,true);
    double gap = 100*(upperbound - SOCP.get_obj_val())/upperbound;
    DebugOn("Final Gap = " << to_string(gap) << "%."<<endl);
//    for (int i = 0; i < v.size(); i++)
//        cout << get<0>(v[i])<< " "
//        << get<1>(v[i]) << " "
//        << get<2>(v[i]) << "\n";
    
    
    
    
    
    /************ Collecting the indices of variables that appear in the nonzero constraints ************/
    //    indices myIdx;
    //    string inst;
    //
    //    for (int i = 0; i < v.size(); i++){
    //        inst = SOCP._cons[get<1>(v[i])]->_vars->begin()->second.first->_indices->_keys->at(get<2>(v[i]));
    //        cout << inst << endl;
    //        cout << SOCP._cons[get<1>(v[i])]->_indices->_keys->at(get<2>(v[i])) << endl;
    //        cout << SOCP._cons[get<1>(v[i])]->_vars->begin()->second.first->_indices->_keys->at(get<2>(v[i])) << endl;
    ////        inst = SOCP._cons[get<1>(v[i])]->_indices->_keys->at(get<2>(v[i]));
    //        myIdx.add(inst);
    //        //   alternatively, you can do
    //        //   SOCP._cons[get<1>(v[i])]->_vars->begin()->second.first->_indices->_keys->at(get<2>(v[i]));
    //    }
    //
    //
    //    /* THERE IS A BIG BUG IN HERE */
    //    /* The constraint indices are correct but the variable indices are wrong and basically the first 3 constraints are produces in here. (1,4), (4,5), (5,6) are the first three constraints. */
    //    Model<> asd("asd");
    //    asd.add(R_Wij);
    //    asd.add(Im_Wij);
    //    asd.add(Wii);
    //    Constraint<> SOC("SOC");
    //    SOC = pow(R_Wij, 2) + pow(Im_Wij, 2) - Wii.from(bus_pairs)*Wii.to(bus_pairs);
    //    asd.add(SOC.in(myIdx) <= 0);
    //
    //    asd.print();
    
    /* Think about the extra-partition and representation of the same W_ii in the bilinear product, we need to link lambdas somehow*/
    
    /* Implement on-off constraints for the squared term in a symbolic way */
    
//    var<> xtemp("xtemp");
//    xtemp.in(R(5));
//    Constraint<> xtempcons("xtempcons");
//    param<> temp("temp");
//    temp.add_val(1);
//    temp.add_val(2);
//    temp.add_val(2);
//    temp.add_val(4);
//    temp.add_val(5);
//    Constraint<> xtempcons2("xtempcons2");
//    xtempcons2 = xtemp;
//    xtempcons = xtemp;
//    SOCP.add(xtempcons >= temp);
//    SOCP.add(xtempcons2 <= temp);
//
//
//    xtempcons.get_cst()->print();
//    xtempcons2.get_cst()->print();
    
 
    
    return 0;
    
}

