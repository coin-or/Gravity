//
//  SOC_partition_3d.cpp
//  Gravity
//
//  Created by Mertcan Yetkin on 8/5/19.
//

#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <PowerNet.h>
#include <gravity/solver.h>
#include <optionParser.hpp>

using namespace std;
using namespace gravity;

int main (int argc, char * argv[])
{
    int output = 0;
    double tol = 1e-6;
    double solver_time_end, solve_time;
    
    //partition parameters declaration
    int var_partns_quad = 10;
    int var_partns_bln = 5;
    int SOC_partns1 = 100; //for the first SOC constraint
    int SOC_partns2 = 10; // for the second SOC constraint
    
    /** MODEL DECLARATION */
    Model<> P_Vars("P_Vars"); //model for the partition on the variables
    Model<> P_SOC("P_SOC"); //model for the partition on the SOC itself
    
    indices varIdx("varIdx");
    varIdx.add("1","2");
    
    /** Variables */
    var<> W1("W1",-5,5);
    P_Vars.add(W1.in(range(1,1)));
    P_SOC.add(W1.in(range(1,1)));
    
    var<> W2("W2",-5,5);
    P_Vars.add(W2.in(range(1,1)));
    P_SOC.add(W2.in(range(1,1)));
    
    var<> W3("W3", 3,5);
    P_Vars.add(W3.in(range(1,1)));
    P_SOC.add(W3.in(range(1,1)));
    
    var<> W4("W4",3,5);
    P_Vars.add(W4.in(range(2,2)));
    P_SOC.add(W4.in(range(2,2)));
    
    var<> t("t",0,5);
    P_Vars.add(t.in(range(1,1)));
    P_SOC.add(t.in(range(1,1)));
    
    
    //sets
    indices myIdx1("myIdx1");
    for (int i=0; i<SOC_partns1; ++i) {
//        myIdx1.add(to_string(i+1)+",W1{1},W2{1},t{1}");
        myIdx1.add(to_string(i+1));
    }
    
    indices myIdx2("myIdx2");
    for (int i=0; i<SOC_partns2; ++i) {
//        myIdx2.add(to_string(i+1)+",t{1},W3{1},W4{1}");
        myIdx2.add(to_string(i+1));
    }
    
    //indicator variables for the SOC partitions
    var<int> on1("on1",0,1); //for the first partitioning of the SOC constraint
    P_SOC.add(on1.in(myIdx1));
    
    var<int> on2("on2",0,1); //for the second partitioning of the SOC constraint
    P_SOC.add(on2.in(myIdx2));
    
    //declare the parameters and their size
    //for the first SOC planes
    param<> W1_coef("W1_coef");
    W1_coef.in(myIdx1);
    
    param<> W2_coef("W2_coef");
    W2_coef.in(myIdx1);
    
    param<> t_coef1("t_coef1");
    t_coef1.in(myIdx1);
    
    //for the second SOC planes
    param<> W3_minus_W4_coef("W3_minus_W4_coef");
    W3_minus_W4_coef.in(myIdx2);
    
    param<> W3_plus_W4_coef("W3_plus_W4_coef");
    W3_plus_W4_coef.in(myIdx2);
    
    param<> t_coef2("t_coef2");
    t_coef2.in(myIdx2);
    
    //fill the values of parameters
    //normal coefficients for SOC cone
    for (int i=0; i<SOC_partns1; ++i) {
        auto t_val = std::cos(2*M_PI*i/SOC_partns1)*std::sin(2*M_PI*(i+1)/SOC_partns1) - std::cos(2*M_PI*(i+1)/SOC_partns1)*std::sin(2*M_PI*i/SOC_partns1);
        auto W1_val = std::sin(2*M_PI*i/SOC_partns1)-std::sin(2*M_PI*(i+1)/SOC_partns1);
        auto W2_val = std::cos(2*M_PI*(i+1)/SOC_partns1)-std::cos(2*M_PI*i/SOC_partns1);
        if (t_val < 0)
        {
            t_val = -t_val;
            W1_val = -W1_val;
            W2_val = -W2_val;
        }
//        t_coef1.set_val(to_string(i+1)+",W1{1},W2{1},t{1}", t_val);
//        W1_coef.set_val(to_string(i+1)+",W1{1},W2{1},t{1}", W1_val);
//        W2_coef.set_val(to_string(i+1)+",W1{1},W2{1},t{1}", W2_val);
        t_coef1.set_val(to_string(i+1), t_val);
        W1_coef.set_val(to_string(i+1), W1_val);
        W2_coef.set_val(to_string(i+1), W2_val);
        
    }
    //coefficients for the rotated SOC cone (t^2 <= W3*W4 =======> 4t^2 + (W3-W4)^2 <= (W3+W4)^2)
    for (int j=0; j<SOC_partns2; ++j) {
        auto W3_plus_W4_val = std::cos(2*M_PI*j/SOC_partns2)*std::sin(2*M_PI*(j+1)/SOC_partns2) - std::cos(2*M_PI*(j+1)/SOC_partns2)*std::sin(2*M_PI*j/SOC_partns2);
        auto t_val = 2*(std::sin(2*M_PI*j/SOC_partns2)-std::sin(2*M_PI*(j+1)/SOC_partns2));
        auto W3_minus_W4_val = std::cos(2*M_PI*(j+1)/SOC_partns2)-std::cos(2*M_PI*j/SOC_partns2);
        if (W3_plus_W4_val < 0)
        {
            W3_plus_W4_val = -W3_plus_W4_val;
            t_val = -t_val;
            W3_minus_W4_val = -W3_minus_W4_val;
        }
//        W3_plus_W4_coef.set_val(to_string(j+1)+",t{1},W3{1},W4{1}", W3_plus_W4_val);
//        t_coef2.set_val(to_string(j+1)+",t{1},W3{1},W4{1}", t_val);
//        W3_minus_W4_coef.set_val(to_string(j+1)+",t{1},W3{1},W4{1}", W3_minus_W4_val);
        W3_plus_W4_coef.set_val(to_string(j+1), W3_plus_W4_val);
        t_coef2.set_val(to_string(j+1), t_val);
        W3_minus_W4_coef.set_val(to_string(j+1), W3_minus_W4_val);
        
    }
    
    /**  Objective */
    auto obj = pow(W1,2) + pow(W2,2);
    //    auto obj = pow(t,2);
    P_Vars.min(obj);
    P_SOC.min(obj);
    
    /** Constraints */
    
    // set the partition numbers for the variables
    W1._num_partns = var_partns_quad;
    W2._num_partns = var_partns_quad;
    W3._num_partns = var_partns_bln;
    W4._num_partns = var_partns_bln;
    
    //SOC constraint
    Constraint<> SOC("SOC");
    SOC = pow(W1,2) + pow(W2,2) - W3*W4;
    P_Vars.add(SOC.in(range(1,1)) == 0, true);
    P_SOC.add(SOC.in(range(1,1)) <= 0);
    
    //SOC intermediaries
    Constraint<> SOC1("SOC1");
    SOC1 = pow(W1,2) + pow(W2,2) - pow(t,2);
    SOC1 <= 0;
    P_SOC.add(SOC1.in(range(1,1)));
    
    Constraint<> SOC2("SOC2");
    SOC2 = pow(2*t,2) + pow(W3-W4,2) - pow(W3+W4,2);
    SOC2 <= 0;
    SOC2.print();
    P_SOC.add(SOC2.in(range(1,1)));
    
//    Constraint<> SOC3("SOC3");
//    SOC3 = pow(t,2) - W3*W4;
//    SOC3 == 0;
//    P_SOC.add(SOC3.in(range(1,1)));
    
    
    //reset the partition numbers
    
    // set the partition numbers for the variables
    W1._in_SOC_partn = true;
    W2._in_SOC_partn = true;
    W3._in_SOC_partn = true;
    W4._in_SOC_partn = true;
    
    //SOC hyperplanes
    Constraint<> SOC_Planes1("SOC_Planes1");
    SOC_Planes1 = W1_coef * W1 + W2_coef * W2 + t_coef1*t;
    SOC_Planes1.in(myIdx1) <= 0;
    P_SOC.add_on_off_multivariate_refined(SOC_Planes1, on1);
    
    Constraint<> SOC_Planes2("SOC_Planes2");
    SOC_Planes2 = t_coef2 * t + W3_minus_W4_coef * (W3-W4) + W3_plus_W4_coef*(W3+W4);
    SOC_Planes2.in(myIdx2) <= 0;
    P_SOC.add_on_off_multivariate_refined(SOC_Planes2, on2);
    
    //Indicator variables for SOC hyperplanes
    Constraint<> onSum1("onSum1");
    onSum1 = sum(on1);
    P_SOC.add(onSum1 == 1);
    
    Constraint<> onSum2("onSum2");
    onSum2 = sum(on2);
    P_SOC.add(onSum2 == 1);
    
    //Create solvers and run
    solver<> P_Vars_CPX(P_Vars, cplex);
    auto solver_time_start = get_wall_time();
    P_Vars_CPX.run(output,tol = 1e-6);
    solver_time_end = get_wall_time();
    solve_time = solver_time_end - solver_time_start;
    
    //    P_Vars.print_solution();
    
    DebugOn("For P_Var: Solve time: " << to_string(solve_time) << endl);
    auto nonzero_idx1 = P_Vars.sorted_nonzero_constraint_indices(tol, true, "SOC_convex");
    
    solver<> P_SOC_CPX(P_SOC, cplex);
    solver_time_start = get_wall_time();
    P_SOC_CPX.run(output,tol = 1e-6);
    solver_time_end = get_wall_time();
    solve_time = solver_time_end - solver_time_start;
    
    //    P_SOC.print_solution();
    
    DebugOn("For P_SOC: Solve time: " << to_string(solve_time) << endl);
    auto nonzero_idx2 = P_SOC.sorted_nonzero_constraint_indices(tol, true, "SOC");
    
    auto SOC1_bool_soc = SOC1.check_soc();
    auto SOC1_bool_rotated_soc = SOC1.check_rotated_soc();
    
    auto SOC2_bool_soc = SOC2.check_soc();
    auto SOC2_bool_rotated_soc = SOC2.check_rotated_soc();
    
//    auto SOC3_bool_soc = SOC3.check_soc();
//    auto SOC3_bool_rotated_soc = SOC3.check_rotated_soc();
    
    DebugOn("SOC1 stats: soc-" << SOC1_bool_soc << ", rotated-" << SOC1_bool_rotated_soc << endl);
    DebugOn("SOC2 stats: soc-" << SOC2_bool_soc << ", rotated-" << SOC2_bool_rotated_soc << endl);
//    DebugOn("SOC3 stats: soc-" << SOC3_bool_soc << ", rotated-" << SOC3_bool_rotated_soc << endl);
    
//    P_SOC.on_off_SORC_partition(SOC, 10);
    
    return 0;
}
