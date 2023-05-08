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
    indices B_ids("B_ids"), C_ids("C_ids"), W_ids("W_ids");
    B.in(B_ids);
    C.in(C_ids);
    W.in(W_ids);

    // Gemm indices
    indices Gemms("Gemms"), Gemms_out("Gemms_out"), Gemms_in("Gemms_in"), B_Gemm("B_Gemm"), C_Gemm("C_Gemm");
    Gemms_out = hidden_states;
    Gemms_in  = hidden_states;
    B_Gemm    = B_ids;
    C_Gemm    = C_ids;

    // Conv indices
    indices Convs("Convs"); // Constraint indices, 1 per output
    indices Convs_out("Convs_out"), Convs_in("Convs_in"); // Input/output indices, 1 per input/output
    indices W_Conv("W_Conv"), B_Conv("B_Conv"); // Parameter indices, 1 per parameter (weight / bias)
    Convs_out = hidden_states;
    Convs_in  = hidden_states;
    W_Conv    = W_ids;
    B_Conv    = B_ids;

    // ReLU indices
    indices ReLUs("ReLUs"), ReLUs_out("ReLUs"), ReLUs_in("ReLUs_in");
    ReLUs_out = hidden_states;
    ReLUs_in  = hidden_states;

    // NoOp indices
    indices NoOps("NoOps"), NoOps_out("NoOps_out"), NoOps_in("NoOps_in");
    NoOps_out = hidden_states;
    NoOps_in  = hidden_states;

    nn.index_hidden_states(hidden_states, y_ids);

    size_t gemm_row_id = 0, conv_row_id = 0;
    for (auto l: nn.layers) {
        switch (l->operator_type) {
            case _relu: {
                auto relu = dynamic_cast<Relu*>(l);
                relu->build_constraints(ReLUs, ReLUs_in, ReLUs_out);
                break;
            }
            case _gemm:{
                auto gemm = dynamic_cast<GEMM*>(l);
                gemm->build_constraints(Gemms, Gemms_in, Gemms_out, B_Gemm, C_Gemm, B, C, gemm_row_id);
                break;
            }
            case _conv:{
                auto conv = dynamic_cast<Conv*>(l);
                conv->build_constraints(Convs, Convs_in, Convs_out, W_Conv, B_Conv, W, B, conv_row_id);
                break;
            }
            case _noop:{
                auto noop = dynamic_cast<NoOp*>(l);
                noop->build_constraints(NoOps, NoOps_in, NoOps_out);
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

    /* Objective function */
    // NN.min(
    //    x(nn.layers.back()->name+"_out,4")
    //    -x(nn.layers.back()->name+"_out,9")
    // );
    NN.max(x(nn.layers.back()->outputs[0]->strkey(0)));

    /* Constraints */
    Constraint<> ReLU("ReLU");
    ReLU = x.in(ReLUs_out) - x.in(ReLUs_in);
    NN.add(ReLU.in(ReLUs) >= 0);

    Constraint<> ReLU_on("ReLU_on");
    ReLU_on = x.in(ReLUs_out) - x.in(ReLUs_in);
    NN.add_on_off(ReLU_on.in(ReLUs) <= 0, y.in(ReLUs), true);

    Constraint<> ReLU_y_off("ReLU_y_off");
    ReLU_y_off = x.in(ReLUs_out);
    NN.add_on_off(ReLU_y_off.in(ReLUs) <= 0, y.in(ReLUs), false);

    Constraint<> Gemm("Gemm");
    Gemm = x.in(Gemms_out) - (x.in(Gemms_in)*B.in(B_Gemm) + C.in(C_Gemm));
    NN.add(Gemm.in(Gemms) == 0);

    Constraint<> Conv_("Conv");
    Conv_ = x.in(Convs_out) - (x.in(Convs_in)*W.in(W_Conv) + B.in(B_Conv));
    NN.add(Conv_.in(Convs) == 0);

    Constraint<> NoOp("NoOp");
    NoOp = x.in(NoOps_out) - x.in(NoOps_in);
    NN.add(NoOp.in(NoOps) == 0);

    NN.print();
    NN.write();

    solver<> S(NN,gurobi);
    S.run();

    NN.print_solution();

    auto sol = std::vector<double>();
    NN.get_solution(sol);
    sol.resize(nn.input_numel);

    std::cout << "Solution: " << print_vector(sol) << std::endl;
    std::cout << "Obj. value: " << std::setprecision(8) << NN.get_obj_val() << std::endl;

    return 0;
}
