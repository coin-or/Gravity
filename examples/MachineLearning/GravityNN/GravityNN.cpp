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
        else if (node.op_type() == "MatMul") {
            layers.push_back(new MatMul(node, initializers));
        }
        else if (node.op_type() == "Add") {
            layers.push_back(new Add(node, initializers));
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
    string fname = string(prj_dir)+"/data_sets/VNN/MatMul_Add.onnx";
    std::fstream input(fname, std::ios::in | std::ios::binary);
    onnx::ModelProto model;
    bool isSuccess = model.ParseFromIstream(&input);
    onnx::GraphProto graph = model.graph();
    
    auto layers = build_graph(graph);
    auto input_dims = get_input_dim(graph.input());/* Getting input layer dim */
    /* INDEX SETS */
    /* Indexing variables */
    indices x_ids("x_ids"), y_ids("y_ids");/*< x_ids for continuous vars, y_ids for binary vars */
    
    /* Indexing params */
    indices B_ids("B_ids"), C_ids("C_ids");/*< x_ids for continuous vars, y_ids for binary vars */

    /* Creating variables' indexing sets */
    for(auto j = 0; j < input_dims[1];j++){
        x_ids.add("input"+to_string(j));
    }
    auto nb_layers = layers.size();
    vector<size_t> dims;
    bool is_activation_func = false;
    for(auto i = 0; i<nb_layers; i++){
        auto l = layers[i];
        if(i+1<nb_layers && l->var_dims.size()==2 && layers[i+1]->operator_type==_relu){/*< The next layer is an activation layer (e.g. ReLU), no need to introduce variables for this layer */
            continue;
        }
        if(l->is_activation_func){/*< Activation functions (e.g., ReLU) use the output of the previous layer to index the output vars */
            dims = layers[i-1]->var_dims;
        }
        else {
            dims = l->var_dims;
        }
        if(dims.size()==2){
            for(auto j = 0; j < dims[1];j++){
                x_ids.add(l->name+","+to_string(j));
                if(l->operator_type==_relu){
                    y_ids.add(l->name+","+to_string(j));
                }
            }
        }
        else{
            for(auto j = 0; j < dims[0];j++){
                x_ids.add(l->name+","+to_string(j));
            }
        }
            
    }
    /* Indexing constraints */
    indices ReLUs("ReLUs"), x_ReLUs("x_ReLUs"), Gemms("Gemms"), x_Gemms("x_Gemms"), Adds("Adds"), x_Adds("x_Adds"), MatMuls("MatMuls"), x_MatMuls("MatMuls");
    param<> B("B");
    x_ReLUs = x_ids;
    x_Gemms = x_ids;
    x_MatMuls = x_ids;
    x_Adds = x_ids;
    size_t idx = 0;
    for(auto i = 0; i<nb_layers; i++){
        auto l = layers[i];
        if(i+1<nb_layers && l->var_dims.size()==2 && layers[i+1]->is_activation_func){/*< The next layer is an activation layer (e.g. ReLU), no need to introduce variables for this layer */
            for(auto j = 0; j < l->var_dims[1];j++){
                if(layers[i+1]->operator_type==_relu){
                    for(auto k = 0; k < l->var_dims[0];k++){
                            x_ReLUs.add_in_row(idx, l->name+","+to_string(k));
                    }
                    idx++;
                }
            }
        }
        switch (l->operator_type) {
            case _relu:
                /* Creating ReLU indexing sets */
                if(layers[i-1]->var_dims.size()==2){
                    for(auto j = 0; j < layers[i-1]->var_dims[1];j++){
                        ReLUs.add(l->name+","+to_string(j));
                    }
                }
                else{
                    for(auto j = 0; j < layers[i-1]->var_dims[0];j++){
                        ReLUs.add(l->name+","+to_string(j));
                    }
                }
                break;
                
            case _gemm:
                for(auto j = 0; j < l->var_dims[1];j++){
                    Gemms.add(l->name+","+to_string(j));
                }
                
            case _add:
                for(auto j = 0; j < l->var_dims[0];j++){
                    auto Add_layer = (Add*)l;
                    Adds.add(l->name+","+to_string(j));
                    B.add_val(l->name+","+to_string(j), Add_layer->B.data[j]);
                }

            default:
                break;
        }
        
        
    }
        
    Model<> NN("NN_"+fname.substr(fname.find_last_of("/")));
    var<> x("x", -1, 1);
    var<int> y("y", 0, 1);
    NN.add(x.in(x_ids));
    NN.add(y.in(y_ids));
   
    Constraint<> ReLU("ReLU");
    ReLU = x.in(x_ReLUs) - x.in(ReLUs);// make this one disjunctive in y
    NN.add(ReLU.in(ReLUs) == 0);
    
    Constraint<> Gemm("Gemm");
    Gemm = x.in(x_Gemms);
    NN.add(Gemm.in(Gemms) == 0);
    
    Constraint<> Add("Add");
    Add = x.in(Adds) + B.in(Adds);
    NN.add(Add.in(Adds) == 0);

    NN.print();
    for(auto l:layers){
        delete l;
    }
    

    return 0;
}
