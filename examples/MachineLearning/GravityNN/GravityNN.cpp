#include <iostream>
#include <fstream>
#include <onnx.pb.h>
#include <vector>
#include "Layers.hpp"
#include <gravity/solver.h>

using namespace gravity;

std::tuple<std::vector<Layer*>, size_t, Tensors> build_graph(std::string fname) {
    std::fstream input(fname, std::ios::in | std::ios::binary);
    onnx::ModelProto model;
    bool isSuccess = model.ParseFromIstream(&input);
    onnx::GraphProto graph = model.graph();

    std::vector<Layer*> layers;
    Tensors tensors;
    // Tensors with data
    for (const auto& initializer : graph.initializer()) {
        tensors[initializer.name()] = Tensor(initializer);
    }

    // Tensor with shape/metadata only
    for (const auto& vinfo : graph.value_info()) {
        tensors[vinfo.name()] = Tensor(vinfo);
    }
    for (const auto& input : graph.input()) {
        tensors[input.name()] = Tensor(input);
    }
    for (const auto& output : graph.output()) {
        tensors[output.name()] = Tensor(output);
    }

    for (const auto& node : graph.node()) {
        if (node.op_type() == "Gemm") {
            layers.push_back(new GEMM(node, tensors));
        } else if (node.op_type() == "Relu") {
            layers.push_back(new Relu(node, tensors));
        } else if (node.op_type() == "Conv") {
            layers.push_back(new Conv(node, tensors));
        } else {
            throw std::runtime_error("Unsupported operator " + node.op_type());
        }
    }

    for (auto i = 0; i < layers.size()-1; i++) {
        layers[i]->is_pre_activation = layers[i+1]->is_activation_func;
    }

    for (const auto& layer : layers) {
        layer->print();
    }

    std::vector<size_t> input_dims;
    if (graph.input_size() > 1) {
        throw std::runtime_error("Network has more than one input. Not supported.");
    } else {
        input_dims = tensors[graph.input(0).name()].shape;
    }
    size_t input_numel = vecprod(input_dims);

    return {layers, input_numel, tensors};
}

int main(int argc, char * argv[]){
    string fname = string(prj_dir)+"/data_sets/VNN/tll_bound.onnx";
    if(argc>=2){
        fname=argv[1];
    }

    auto tmp = build_graph(fname);
    auto layers = std::get<0>(tmp);
    auto input_numel = std::get<1>(tmp);
    auto tensors = std::get<2>(tmp);

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
    C_Gemm    = Gemms;

    // Conv indices
    indices Convs("Convs"); // Constraint indices, 1 per output
    indices Convs_out("Convs_out"), Convs_in("Convs_in"); // Input/output indices, 1 per input/output
    indices W_Conv("W_Conv"), B_Conv("B_Conv"); // Parameter indices, 1 per parameter (weight / bias)
    Convs_out = hidden_states;
    Convs_in  = hidden_states;
    W_Conv    = W_ids;
    B_Conv    = Convs;

    // ReLU indices
    indices ReLUs("ReLUs"), ReLUs_out("ReLUs"), ReLUs_in("ReLUs_in");
    ReLUs_out = hidden_states;
    ReLUs_in  = hidden_states;

    for (auto i = 0; i < input_numel; i++) {
        hidden_states.add("input,"+to_string(i));
    }

    for(auto i = 0; i < layers.size(); i++){
        auto l = layers[i];
        for(auto j = 0; j < l->outputs.at(0).numel; j++){
            std::string key = l->name+"_out,"+to_string(j);
            hidden_states.add(key);
            if(l->operator_type==_relu){
                y_ids.add(l->name+","+to_string(j));
            }
        }
    }

    size_t gemm_row_id = 0, conv_row_id = 0;
    for(auto i = 0; i < layers.size(); i++){
        auto l = layers[i];
        // Input indexing
        std::string input_key = (i == 0) ? "input," : layers[i-1]->name+"_out,";
        std::string output_key = l->name+"_out,";
        switch (l->operator_type) {
            case _relu: {
                for(auto j = 0; j < l->inputs.at(0).numel;j++){
                    ReLUs.add(l->name + "," + to_string(j));
                    ReLUs_out.add_ref(output_key+to_string(j));
                    ReLUs_in.add_ref(input_key+to_string(j));
                }
                break;
            }
            case _gemm:{
                auto gemm = reinterpret_cast<GEMM*>(l);
                gemm->add_parameters({&B, &C});

                // Output indexing
                for (auto j = 0; j < l->outputs.at(0).numel; j++) {
                    Gemms.add(gemm->name + "," + to_string(j));
                    Gemms_out.add_ref(output_key+to_string(j));
                }

                // Expression
                for(auto j = 0; j < l->outputs[0].shape[1];j++){
                    for(auto k = 0; k < l->inputs[0].shape[1]; k++){
                        std::string weight_id = gemm->name+","+to_string(k)+","+to_string(j);
                        std::string input_id = input_key+to_string(k);

                        B_Gemm.add_in_row(gemm_row_id, weight_id);
                        Gemms_in.add_in_row(gemm_row_id, input_id);
                    }
                    gemm_row_id++;
                }
                break;
            }
            case _conv:{
                auto conv = reinterpret_cast<Conv*>(l);
                conv->add_parameters({&W, &B});

                // Output indexing
                for (auto j = 0; j < l->outputs.at(0).numel; j++) {
                    Convs.add(conv->name + "," + to_string(j));
                }

                for (int oh = 0; oh < conv->out_h; oh++) {
                    for (int ow = 0; ow < conv->out_w; ow++) {
                        for (int oc = 0; oc < conv->out_c; oc++) {
                            Convs_out.add_ref(output_key+to_string(conv->outputs.at(0).flatten(0, oc, oh, ow)));
                            for (int kh = 0; kh < conv->kern_h; kh++) {
                                for (int kw = 0; kw < conv->kern_w; kw++) {
                                    for (int kc = 0; kc < conv->kern_c; kc++) {
                                        int h_ind = (conv->strides[0]*oh + conv->dilations[0]*kh - conv->pads[0]);
                                        int w_ind = (conv->strides[1]*ow + conv->dilations[1]*kw - conv->pads[3]);
                                        if ((h_ind < conv->inp_h) && (h_ind >= 0) && (w_ind < conv->inp_w) && (w_ind >= 0)) {
                                            std::string w_idx = to_string(oc)+","+to_string(kc)+","+to_string(kh)+","+to_string(kw);
                                            std::string weight_id = conv->name + "," + w_idx;
                                            std::string input_id = input_key+to_string(conv->inputs.at(0).flatten(0, kc, h_ind, w_ind));

                                            W_Conv.add_in_row(conv_row_id, weight_id);
                                            Convs_in.add_in_row(conv_row_id, input_id);
                                        }
                                    }
                                }
                            }

                            conv_row_id++;
                        }
                    }
                }
                break;
            }
            default:{
                break;
            }
        }
    }


    Model<> NN("NN_"+fname.substr(fname.find_last_of("/")));
    param<> x_lb("x_lb"), x_ub("x_ub");
    x_lb.in(hidden_states);x_ub.in(hidden_states);
    x_lb = std::numeric_limits<double>::lowest();
    x_ub = std::numeric_limits<double>::max();

    if (tensors.count("Input0_lower") != 0) {
        auto lower = tensors.at("Input0_lower");
        auto upper = tensors.at("Input0_upper");
        for(auto j = 0; j < input_numel;j++){
            x_lb.set_val("input,"+to_string(j), lower(j));
            x_ub.set_val("input,"+to_string(j), upper(j));
        }
    }

    for(auto i = 0; i < ReLUs_out.size(); i++){
       auto key = ReLUs_out._keys->at(ReLUs_out.get_id_inst(i));
       x_lb.set_val(key, 0);
    }

    for(auto i = 0; i < layers.size(); i++){
       auto l = layers[i];
       for(auto j = 0; j < l->lowers[0].numel;j++){
           x_lb.set_val(l->name+"_out,"+to_string(j), l->lowers[0](j));
           x_ub.set_val(l->name+"_out,"+to_string(j), l->uppers[0](j));
       }
    }

    var<> x("x", x_lb, x_ub);
    var<int> y("y", 0, 1);
    NN.add(x.in(hidden_states));
    NN.add(y.in(y_ids));

    /* Objective function */
    NN.max(x(layers.back()->name+"_out,0"));

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

    // NN.print();
    NN.write();

    solver<> S(NN,gurobi);
    S.run();

    NN.print_solution();

    auto sol = std::vector<double>();
    NN.get_solution(sol);
    sol.resize(input_numel);

    std::cout << "Solution: " << print_vector(sol) << std::endl;

    for(auto l:layers){
        delete l;
    }

    return 0;
}
