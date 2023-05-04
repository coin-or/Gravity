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


    /* INDEX SETS */
    /* Indexing variables */
    indices Layer_ids("Layer_ids"), y_ids("y_ids");/*< x_ids for continuous vars, y_ids for binary vars */

    /* Indexing params */
    indices B_ids("B_ids"), C_ids("C_ids");/*< x_ids for continuous vars, y_ids for binary vars */
    indices ReLUs("ReLUs"), ReLUs_in("ReLUs_in"), ReLUs_out("ReLUs_out"), Layers_in("Layers_in"), Layers_out("Layers_out"), Gemms("Gemms"), Gemms_in("Gemms_in"), Gemms_out("Gemms_out"), B_Gemm("B_Gemm"), C_Gemm("C_Gemm");
    param<> B("B"), C("C");
    B.in(B_ids);
    C.in(C_ids);
    Gemms_in = Layer_ids;
    Gemms_out = Layer_ids;
    ReLUs_in = Layer_ids;
    ReLUs_out = Layer_ids;
    C_Gemm = Gemms;
    Layers_in = Layer_ids;
    Layers_out = Layer_ids;

    for(auto i = 0; i < layers.size(); i++){
        auto l = layers[i];
        for(auto j = 0; j < l->inputs.at(0).numel; j++){
            Layer_ids.add(l->name+"_in,"+to_string(j));
            if(i>0)
                Layers_out.add_ref(l->name+"_in,"+to_string(j));
            if(l->operator_type==_relu){
                y_ids.add(l->name+","+to_string(j));
            }
        }
        for(auto j = 0; j < l->outputs.at(0).numel; j++){
            Layer_ids.add(l->name+"_out,"+to_string(j));
            if(i< layers.size()-1)
                Layers_in.add_ref(l->name+"_out,"+to_string(j));
        }
    }
    /* Indexing constraints */
    B_Gemm = B_ids;
    Gemms_in = Layer_ids;
    size_t relu_row_id = 0, gemm_row_id = 0;
    for(auto i = 0; i < layers.size(); i++){
        auto l = layers[i];
        switch (l->operator_type) {
            case _relu: {
                for(auto j = 0; j < l->outputs.at(0).numel;j++){
                    ReLUs.add(l->name+","+to_string(j));
                    ReLUs_in.add_ref(l->name+"_in,"+to_string(j));
                    ReLUs_out.add_ref(l->name+"_out,"+to_string(j));
                }
                break;
            }
            case _gemm:{
                auto gemm = reinterpret_cast<GEMM*>(l);
                gemm->B.add_params(B, gemm->name);
                if (gemm->has_optional_C) {
                    gemm->C.add_params(C, gemm->name);
                }


                for(auto j = 0; j < l->outputs[0].shape[1];j++){
                    Gemms.add(l->name+","+to_string(j));
                    Gemms_out.add_ref(l->name+"_out,"+to_string(j));
                    for(auto k = 0; k < l->inputs[0].shape[1]; k++){
                        B_Gemm.add_in_row(gemm_row_id, gemm->name+","+to_string(k)+","+to_string(j));
                        Gemms_in.add_in_row(gemm_row_id, gemm->name+"_in,"+to_string(k));
                    }
                    gemm_row_id++;
                }
                break;
            }
            default:
                break;
        }
    }


    Model<> NN("NN_"+fname.substr(fname.find_last_of("/")));
    param<> x_lb("x_lb"), x_ub("x_ub");
    x_lb.in(Layer_ids);x_ub.in(Layer_ids);
    x_lb = -1000;
    x_ub = 1000;

    if (tensors.count("Input0_lower") != 0) {
        auto lower = tensors.at("Input0_lower");
        auto upper = tensors.at("Input0_upper");
        for(auto j = 0; j < input_numel;j++){
            x_lb.set_val(layers[0]->name+"_in,"+to_string(j), lower(j));
            x_ub.set_val(layers[0]->name+"_in,"+to_string(j), upper(j));
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
    NN.add(x.in(Layer_ids));
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

    Constraint<> Linking_Layers("Linking_Layers");
    Linking_Layers = x.in(Layers_in) - x.in(Layers_out);
    NN.add(Linking_Layers.in(range(1,Layers_in.size())) == 0);

    NN.print();
    NN.write();

    solver<> S(NN,gurobi);
    S.run();

    // NN.print_solution();

    auto sol = std::vector<double>();
    NN.get_solution(sol);
    sol.resize(input_numel);

    std::cout << "Solution: " << print_vector(sol) << std::endl;


    for(auto l:layers){
        delete l;
    }

    return 0;
}
