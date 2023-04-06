#include <string>
#include <gravity/solver.h>
#include "GravityNN.hpp"
#include "onnx.hpp"
#include "layers.hpp"

int main(int argc, char** argv) {
    string fname = string(prj_dir)+"/data_sets/VNN/mip_lay-16_0.mps";
    Model<> model("new_mip");
    model.readMPS(fname);
//    model = model_from_file(fname);
    solver<> MIP_Solver(model, cplex);
    MIP_Solver.run();
//    model.print();
    model.print_symbolic();
    model.print_solution();
    model.print_int_solution();

    // ReLU upper bound: ReLUX_Y
    // ReLU Bool: aReLUX_Y

    // for (auto v : model._vars_name) {
        // std::cout << v.first << std::endl;
    // }
    
    
	return 0;
}


int readONNX(string fname){
    std::fstream input(fname, std::ios::in | std::ios::binary);
    onnx::ModelProto model;
    bool isSuccess = model.ParseFromIstream(&input);
    onnx::GraphProto graph = model.graph();

    std::map<std::string, onnx::TensorProto> initializers;

    for (const onnx::TensorProto& initializer: graph.initializer()) {
        initializers[initializer.name()] = initializer;
    }

    std::vector<Layer*> layers;
    std::map<std::string, Layer*> layers_by_output;

    NeuralNetwork net;

    int layer_idx = 0;
    for (auto& node: graph.node()) {
        Layer* layer = new Layer(layer_idx, node, initializers);
        layers.push_back(layer);
        for (auto& output: layer->outputs()) {
            layers_by_output[output] = layer;
        }
        layer_idx++;

        net.add_node(layer);
    }

    int arc_idx = 0;
    for (auto& layer: layers) {
        for (auto& input: layer->inputs()) {
            if (layers_by_output.count(input) == 0) {
                continue;
            }

            Layer* src = layers_by_output[input];
            Arc* arc = new Arc(src, layer);
            arc->_id = arc_idx;
            arc_idx++;

            net.add_arc(arc);
        }
    }
}
