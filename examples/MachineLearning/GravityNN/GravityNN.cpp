#include <iostream>
#include <fstream>
#include <onnx.pb.h>
#include <vector>
#include "Layers.hpp"
#include <gravity/solver.h>
#include "NeuralNet.hpp"

using namespace gravity;

int main(int argc, char * argv[]){
    string fname = string(prj_dir)+"/data_sets/VNN/tll_bound.onnx";
    if(argc >= 2) {
        fname = argv[1];
    }

    NeuralNet nn(fname);

    // Global indices
    indices hidden_states("hidden_states"), y_ids("y_ids");/*< x_ids for continuous vars, y_ids for binary vars */

    // Params
    param<> B("B"), C("C"), W("W");
    param<> Min("Min"), Max("Max");
    indices B_ids("B_ids"), C_ids("C_ids"), W_ids("W_ids"), Min_ids("Min_ids"), Max_ids("Max_ids");
    B.in(B_ids);
    C.in(C_ids);
    W.in(W_ids);
    Min.in(Min_ids);
    Max.in(Max_ids);

    // Gemm indices
    IndexSet Gemms({"Constr", "In", "Out", "B", "C"});
    Gemms["Out"] = hidden_states;
    Gemms["In"]  = hidden_states;
    Gemms["B"]   = B_ids;
    Gemms["C"]   = C_ids;

    // Conv indices
    IndexSet Convs({"Constr", "In", "Out", "W", "B"});
    Convs["Out"] = hidden_states;
    Convs["In"]  = hidden_states;
    Convs["W"]   = W_ids;
    Convs["B"]   = B_ids;

    // ReLU indices
    IndexSet ReLUs({"Constr", "In", "Out"});
    ReLUs["Out"] = hidden_states;
    ReLUs["In"]  = hidden_states;

    // NoOp indices
    IndexSet NoOps({"Constr", "In", "Out"});
    NoOps["Out"] = hidden_states;
    NoOps["In"]  = hidden_states;

    // Add indices
    IndexSet Adds({"Constr", "Out", "A", "B"});
    Adds["Out"] = hidden_states;
    Adds["A"]   = hidden_states;
    Adds["B"]   = hidden_states;

    // Sub indices
    IndexSet Subs({"Constr", "Out", "A", "B"});
    Subs["Out"] = hidden_states;
    Subs["A"]   = hidden_states;
    Subs["B"]   = hidden_states;

    // Cos indices
    IndexSet Coss({"Constr", "Out", "In"});
    Coss["Out"] = hidden_states;
    Coss["In"]  = hidden_states;

    // Cos indices
    IndexSet Sins({"Constr", "Out", "In"});
    Sins["Out"] = hidden_states;
    Sins["In"]  = hidden_states;

    // Neg indices
    IndexSet Negs({"Constr", "Out", "In"});
    Negs["Out"] = hidden_states;
    Negs["In"]  = hidden_states;

    // Pow indices
    IndexSet Pows({"Constr", "Out", "In"});
    Pows["Out"] = hidden_states;
    Pows["In"]  = hidden_states;

    // Mul indices
    IndexSet Muls({"Constr", "Out", "A", "B"});
    Muls["Out"] = hidden_states;
    Muls["A"]   = hidden_states;
    Muls["B"]   = hidden_states;

    // Div indices
    IndexSet Divs({"Constr", "Out", "A", "B"});
    Divs["Out"] = hidden_states;
    Divs["A"]   = hidden_states;
    Divs["B"]   = hidden_states;

    nn.index_hidden_states(hidden_states, y_ids);

    size_t gemm_row_id = 0, conv_row_id = 0;
    for (auto l: nn.layers) {
        switch (l->operator_type) {
            case _relu: {
                auto relu = dynamic_cast<Relu*>(l);
                relu->build_constraint(ReLUs, {});
                break;
            }
            case _gemm:{
                auto gemm = dynamic_cast<GEMM*>(l);
                gemm->build_constraint(Gemms, {&B, &C});
                break;
            }
            case _conv:{
                auto conv = dynamic_cast<Conv*>(l);
                conv->build_constraint(Convs, {&W, &B});
                break;
            }
            case _noop:{
                auto noop = dynamic_cast<NoOp*>(l);
                noop->build_constraint(NoOps, {});
                break;
            }
            case _add:{
                auto add = dynamic_cast<Add*>(l);
                add->build_constraint(Adds, {});
                break;
            }
            case _sub:{
                auto sub = dynamic_cast<Sub*>(l);
                sub->build_constraint(Subs, {});
                break;
            }
            case _cos:{
                auto cos = dynamic_cast<Cos*>(l);
                cos->build_constraint(Coss, {});
                break;
            }
            case _sin:{
                auto sin = dynamic_cast<Sin*>(l);
                sin->build_constraint(Sins, {});
                break;
            }
            case _neg:{
                auto neg = dynamic_cast<Neg*>(l);
                neg->build_constraint(Negs, {});
                break;
            }
            case _pow:{
                auto neg = dynamic_cast<Pow*>(l);
                neg->build_constraint(Pows, {});
                break;
            }
            case _mul:{
                auto mul = dynamic_cast<Mul*>(l);
                mul->build_constraint(Muls, {});
                break;
            }
            case _div:{
                auto div = dynamic_cast<Div*>(l);
                div->build_constraint(Divs, {});
                break;
            }
            default:{
                break;
            }
        }
    }

    Model<> NN("NN_"+fname.substr(fname.find_last_of("/")));
    param<> x_lb("x_lb"), x_ub("x_ub");
    x_lb.in(hidden_states);
    x_ub.in(hidden_states);
    x_lb = std::numeric_limits<double>::lowest();
    x_ub = std::numeric_limits<double>::max();

    nn.set_bounds(x_lb, x_ub);

    var<> x("x", x_lb, x_ub);
    var<int> y("y", 0, 1);

    NN.add(x.in(hidden_states));
    NN.add(y.in(y_ids));
    nn.initialize_state(x, y);

    ///* Objective function */
    //NN.min(
    //    x(nn.layers.back()->outputs[0]->strkey(4))
    //   -x(nn.layers.back()->outputs[0]->strkey(9))
    //);
    NN.max(
      x(nn.layers.back()->outputs[0]->strkey(0))
    );

    /* Constraints */
    Constraint<> ReLU("ReLU");
    ReLU = x.in(ReLUs["Out"]) - x.in(ReLUs["In"]);
    NN.add(ReLU.in(ReLUs["Constr"]) >= 0);

    Constraint<> ReLU_on("ReLU_on");
    ReLU_on = x.in(ReLUs["Out"]) - x.in(ReLUs["In"]);
    NN.add_on_off(ReLU_on.in(ReLUs["Constr"]) <= 0, y.in(ReLUs["Constr"]), true);

    Constraint<> ReLU_y_off("ReLU_y_off");
    ReLU_y_off = x.in(ReLUs["Out"]);
    NN.add_on_off(ReLU_y_off.in(ReLUs["Constr"]) <= 0, y.in(ReLUs["Constr"]), false);

    Constraint<> Gemm("Gemm");
    Gemm = x.in(Gemms["Out"]) - (x.in(Gemms["In"])*B.in(Gemms["B"]) + C.in(Gemms["C"]));
    NN.add(Gemm.in(Gemms["Constr"]) == 0);

    Constraint<> Conv_("Conv");
    Conv_ = x.in(Convs["Out"]) - (x.in(Convs["In"])*W.in(Convs["W"]) + B.in(Convs["B"]));
    NN.add(Conv_.in(Convs["Constr"]) == 0);

    Constraint<> NoOp("NoOp");
    NoOp = x.in(NoOps["Out"]) - x.in(NoOps["In"]);
    NN.add(NoOp.in(NoOps["Constr"]) == 0);

    Constraint<> Add_("Add");
    Add_ = x.in(Adds["Out"]) - (x.in(Adds["A"]) + x.in(Adds["B"]));
    NN.add(Add_.in(Adds["Constr"]) == 0);

    Constraint<> Sub_("Sub");
    Sub_ = x.in(Subs["Out"]) - (x.in(Subs["A"]) - x.in(Subs["B"]));
    NN.add(Sub_.in(Subs["Constr"]) == 0);

    Constraint<> Cos_("Cos");
    Cos_ = x.in(Coss["Out"]) - cos(x.in(Coss["In"]));
    NN.add(Cos_.in(Coss["Constr"]) == 0);

    Constraint<> Sin_("Sin");
    Sin_ = x.in(Sins["Out"]) - sin(x.in(Sins["In"]));
    NN.add(Sin_.in(Sins["Constr"]) == 0);

    Constraint<> Neg_("Neg");
    Neg_ = x.in(Negs["Out"]) + x.in(Negs["In"]);
    NN.add(Neg_.in(Negs["Constr"]) == 0);

    Constraint<> Pow_("Pow");
    Pow_ = x.in(Pows["Out"]) - pow(x.in(Pows["In"]), 2.0);
    NN.add(Pow_.in(Pows["Constr"]) == 0);

    Constraint<> Mul_("Mul");
    Mul_ = x.in(Muls["Out"]) - (x.in(Muls["A"]) * x.in(Muls["B"]));
    NN.add(Mul_.in(Muls["Constr"]) == 0);

    Constraint<> Div_("Div");
    Div_ = x.in(Divs["Out"])*x.in(Divs["B"]) - x.in(Divs["A"]);
    NN.add(Div_.in(Divs["Constr"]) == 0);

    // NN.print();
    NN.write();

    solver<> S(NN,gurobi);
    S.run();

    // NN.print_solution();

    auto sol = std::vector<double>();
    NN.get_solution(sol);
    sol.resize(nn.input_numel);

    std::cout << "Solution: " << print_vector(sol) << std::endl;
    std::cout << "Obj. value: " << std::setprecision(8) << NN.get_obj_val() << std::endl;

    return 0;
}
