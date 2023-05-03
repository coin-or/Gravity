#include <iostream>
#include <fstream>
#include <onnx.pb.h>
#include <vector>
#include "Layers.hpp"
#include <gravity/solver.h>

using namespace gravity;

std::pair<std::vector<Layer*>, std::vector<size_t>> build_graph(std::string fname) {
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
        }
        else if (node.op_type() == "MatMul") {
            layers.push_back(new MatMul(node, tensors));
        }
        else if (node.op_type() == "Add") {
            layers.push_back(new Add(node, tensors));
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

    return std::make_pair(layers, input_dims);
}

int main(int argc, char * argv[]){
    string fname = string(prj_dir)+"/data_sets/VNN/tll_bound.onnx";
    if(argc>=2){
        fname=argv[1];
    }

    auto tmp = build_graph(fname);
    auto layers = tmp.first;
    auto input_dims = tmp.second;

    /* INDEX SETS */
    /* Indexing variables */
    indices x_ids("x_ids"), y_ids("y_ids");/*< x_ids for continuous vars, y_ids for binary vars */

    /* Indexing params */
    indices B_ids("B_ids"), C_ids("C_ids");/*< x_ids for continuous vars, y_ids for binary vars */

    param<> B("B"), C("C");
    B.in(B_ids);
    C.in(C_ids);
    /* Creating variables' indexing sets */
    for(auto j = 0; j < input_dims[1];j++){
        x_ids.add("input"+to_string(j));
    }
    for(auto i = 0; i < layers.size(); i++){
        auto l = layers[i];
        if (l->is_pre_activation) { /*< The next layer is an activation layer (e.g. ReLU), no need to introduce variables for this layer */
            continue;
        }
        for(auto j = 0; j < l->outputs.at(0).numel; j++){
            x_ids.add(l->name+","+to_string(j));
            if(l->operator_type==_relu){
                y_ids.add(l->name+","+to_string(j));
            }
        }
    }
    /* Indexing constraints */
    indices ReLUs("ReLUs"), x_ReLUs("x_ReLUs"), B_ReLUs("B_ReLUs"), C_ReLUs("C_ReLUs"), Gemms("Gemms"), x_Gemms_in("x_Gemms_in"), x_Gemms_out("x_Gemms_out"), B_Gemm("B_Gemm");
    B_Gemm = B_ids;
    B_ReLUs = B_ids;
    C_ReLUs = C_ids;
    x_ReLUs = x_ids;
    x_Gemms_in = x_ids;
    x_Gemms_out = x_ids;
    size_t relu_row_id = 0, gemm_row_id = 0;
    string key;
    for(auto i = 0; i < layers.size(); i++){
        auto l = layers[i];
        switch (l->operator_type) {
            case _relu: {
                for(auto j = 0; j < l->outputs.at(0).numel;j++){
                    ReLUs.add(l->name+","+to_string(j));
                }
                break;
            }
            case _gemm:{
                bool add_Gemm_constraint = !l->is_pre_activation;
                auto gemm = (GEMM*)l;
                gemm->B.add_params(B, gemm->name);
                if (gemm->has_optional_C) {
                    gemm->C.add_params(C, gemm->name);
                }

                for(auto j = 0; j < l->outputs[0].shape[1];j++){
                    if(add_Gemm_constraint){
                        Gemms.add(l->name+","+to_string(j));
                    }
                    else{
                        C_ReLUs.add_ref(l->name+","+to_string(j));
                    }
                    for(auto k = 0; k < l->inputs[0].shape[1]; k++){
                        key = gemm->name+","+to_string(k)+","+to_string(j);
                        if(add_Gemm_constraint){
                            B_Gemm.add_in_row(gemm_row_id, key);
                            if(i==0){
                                x_Gemms_in.add_in_row(gemm_row_id, "input"+to_string(k));
                            }
                            else{
                                x_Gemms_in.add_in_row(gemm_row_id, layers[i-1]->name+","+to_string(k));
                            }
                        }
                    }
                    if(add_Gemm_constraint){
                        gemm_row_id++;
                    }
                }
                break;
            }
            default:
                break;
        }
    }

    for(auto i = 0; i<layers.size(); i++){
        auto l = layers[i];
        if (l->is_pre_activation) {/*< The next layer is an activation layer (e.g. ReLU), no need to introduce variables for this layer */
            for(auto j = 0; j < l->outputs[0].shape[1];j++){
                for(auto k = 0; k < l->inputs[0].shape[1];k++){
                    if(i==0)
                        x_ReLUs.add_in_row(relu_row_id, "input"+to_string(k));
                    else
                        x_ReLUs.add_in_row(relu_row_id, layers[i-1]->name+","+to_string(k));
                    B_ReLUs.add_in_row(relu_row_id, l->name+","+to_string(k)+","+to_string(j));
                }
                relu_row_id++;

            }
        }
    }

    Model<> NN("NN_"+fname.substr(fname.find_last_of("/")));
    param<> x_lb("x_lb"), x_ub("x_ub");
    x_lb.in(x_ids);x_ub.in(x_ids);
    x_lb = -1000;
    x_ub = 1000;
    for(auto j = 0; j < input_dims[1];j++){
        x_lb.set_val(j, -2);
        x_ub.set_val(j, 2);
    }

    for(auto const &key: *ReLUs._keys){
        x_lb.set_val(key, 0);
    }

    for(auto i = 0; i < layers.size(); i++){
        auto l = layers[i];
        if(l->is_pre_activation)
            continue;

        for(auto j = 0; j < l->lowers[0].numel;j++){
            x_lb.set_val(l->name+","+to_string(j), l->lowers[0](j));
            x_ub.set_val(l->name+","+to_string(j), l->uppers[0](j));
        }
    }

    var<> x("x", x_lb, x_ub);
    var<int> y("y", 0, 1);
    NN.add(x.in(x_ids));
    NN.add(y.in(y_ids));

    /* Objective function */
    NN.max(x(layers.back()->name+",0"));

    /* Constraints */
    Constraint<> ReLU("ReLU");
    ReLU = x.in(ReLUs) - (x.in(x_ReLUs)*B.in(B_ReLUs) + C.in(C_ReLUs));
    NN.add(ReLU.in(ReLUs) >= 0);

    Constraint<> ReLU_on("ReLU_on");
    ReLU_on = x.in(ReLUs) - (x.in(x_ReLUs)*B.in(B_ReLUs) + C.in(C_ReLUs));
    NN.add_on_off(ReLU_on.in(ReLUs) <= 0, y.in(ReLUs), true);

    Constraint<> ReLU_y_off("ReLU_y_off");
    ReLU_y_off = x.in(ReLUs);
    NN.add_on_off(ReLU_y_off.in(ReLUs) <= 0, y.in(ReLUs), false);

    Constraint<> Gemm("Gemm");
    Gemm = x.in(Gemms) - (x.in(x_Gemms_in)*B.in(B_Gemm) + C.in(Gemms));
    NN.add(Gemm.in(Gemms) == 0);

    NN.print();
    NN.write();
    
    solver<> S(NN,gurobi);
    S.run();

    NN.print_solution();
    for(auto l:layers){
        delete l;
    }

    return 0;
}