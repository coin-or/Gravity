#include <iostream>
#include <fstream>
#include <onnx.pb.h>
#include <vector>
#include "Layers.hpp"
#include <gravity/model.h>
#include <gravity/func.h>

std::vector<Layer*> build_graph(const onnx::GraphProto& graph) {
    std::vector<Layer*> layers;
    Initializers initializers;

    // Parse initializers
    std::cout << "Initializers: " << std::endl;
    for (const auto& initializer : graph.initializer()) {
        initializers[initializer.name()] = initializer;
    }

    for (const auto& node : graph.node()) {
        if (node.op_type() == "Gemm") {
            layers.push_back(new GEMM(node, initializers));
        } else if (node.op_type() == "Relu") {
            layers.push_back(new ReLU(node, initializers));
        }
    }

    for (const auto& layer : layers) {
        layer->print();
    }

//    gravity::param<float> input("input");
//    input = {1, 2, 3, 4};
//    HiddenStates hidden_states;
//    hidden_states[layers[0]->inputs[0]] = input;
//
//    for (const auto& layer : layers) {
//        layer->forward(hidden_states);
//    }
//
//    auto output = hidden_states[(*layers.end())->outputs[0]];
//    output.print();
    return layers;
}




using namespace gravity;
int main() {
    string fname = string(prj_dir)+"/data_sets/VNN/simple.onnx";
    std::fstream input(fname, std::ios::in | std::ios::binary);
    onnx::ModelProto model;
    bool isSuccess = model.ParseFromIstream(&input);
    onnx::GraphProto graph = model.graph();
    
    auto layers = build_graph(graph);
    auto input_dims = get_input_dim(graph.input());/* Getting input layer dim */
    indices ids("ids");
    for(auto j = 0; j < input_dims[1];j++){
        ids.add("input"+to_string(j));
    }
    auto nb_layers = layers.size();
    vector<size_t> dims;
    for(auto i = 0; i<nb_layers; i++){
        auto l = layers[i];
        if(i+1<nb_layers && layers[i+1]->var_dims.empty())/*< The next layer is an activation layer (e.g. ReLU), no need to introduce variables for this layer */
            continue;
        indices M("M_"+l->name);
        if(l->var_dims.empty()){/*< Activation function such as ReLU, use the output of the prvious layer to index the RelU output vars */
            dims = layers[i-1]->var_dims;
        }
        else {
            dims = l->var_dims;
        }
        for(auto j = 0; j < dims[1];j++){
            ids.add(l->name+","+to_string(j));
        }
    }
    
    
    Model<> NN("NN_"+fname.substr(fname.find_last_of("/")));
    var<> x("x");
    NN.add(x.in(ids));
    NN.print();
    for(auto l:layers){
        delete l;
    }
    return 0;
}
