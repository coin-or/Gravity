#include <iostream>
#include <fstream>
#include <onnx.pb.h>
#include <vector>
#include "Layers.hpp"

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


#include <gravity/solver.h>


using namespace gravity;
int main (int argc, char * argv[]){
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
    indices ReLUs("ReLUs"), x_ReLUs("x_ReLUs"), B_ReLUs("B_ReLUs"), C_ReLUs("C_ReLUs"), Gemms("Gemms"), x_Gemms_in("x_Gemms_in"), x_Gemms_out("x_Gemms_out"), Adds("Adds"), x_Adds("x_Adds"), MatMuls("MatMuls"), x_MatMuls("MatMuls"), B_Gemm("B_Gemm");
    B_Gemm = B_ids;
    B_ReLUs = B_ids;
    C_ReLUs = C_ids;
    x_ReLUs = x_ids;
    x_Gemms_in = x_ids;
    x_Gemms_out = x_ids;
    x_MatMuls = x_ids;
    x_Adds = x_ids;
    size_t relu_row_id = 0, gemm_row_id = 0;
    string key;
    for(auto i = 0; i<nb_layers; i++){
        auto l = layers[i];
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

            case _gemm:{
                bool add_Gemm_constraint = true;
                if(i+1<nb_layers && layers[i+1]->operator_type==_relu)
                    add_Gemm_constraint = false;
                auto gemm = (GEMM*)l;
                size_t B_idx = 0, C_idx = 0;
                for(auto j = 0; j < l->var_dims[0];j++){
                    for(auto k = 0; k < l->var_dims[1]; k++){
                        key = gemm->name+","+to_string(j)+","+to_string(k);
                        B.add_val(key, gemm->B(B_idx++));
                    }
                }
//                B.print_vals(3);
                for(auto j = 0; j < l->var_dims[1];j++){
                    if(gemm->has_optional_C){
                        C.add_val(l->name+","+to_string(j), gemm->C(C_idx++));
                    }
                    if(add_Gemm_constraint){
                        Gemms.add(l->name+","+to_string(j));
                    }
                    else{
                        C_ReLUs.add_ref(l->name+","+to_string(j));
                    }
                    for(auto k = 0; k < l->var_dims[0]; k++){
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

            }
                break;
            case _add:
                for(auto j = 0; j < l->var_dims[0];j++){
                    auto Add_layer = (Add*)l;
                    Adds.add(l->name+","+to_string(j));
                    B.add_val(l->name+","+to_string(j), Add_layer->B(j));
                }
                break;
            default:
                break;
        }
    }

    for(auto i = 0; i<nb_layers; i++){
        auto l = layers[i];
        if(i+1<nb_layers && l->var_dims.size()==2 && layers[i+1]->operator_type==_relu){/*< The next layer is an activation layer (e.g. ReLU), no need to introduce variables for this layer */
            for(auto j = 0; j < l->var_dims[1];j++){
                for(auto k = 0; k < l->var_dims[0];k++){
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

    for(auto i = 0; i<nb_layers; i++){
        auto l = layers[i];
        if(i+1<nb_layers && layers[i+1]->operator_type==_relu)
            continue;
        auto layer_name = l->name;
        size_t dim = 0;
        if(l->lowers.size()>0){
            if(l->var_dims.size()==2){
                dim = l->var_dims[1];
            }
            else{
                dim = l->var_dims[0];
            }
            for(auto j = 0; j < dim;j++){
                x_lb.set_val(layer_name+","+to_string(j), l->lowers[0](j));
                x_ub.set_val(layer_name+","+to_string(j), l->uppers[0](j));
            }
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

//    Constraint<> ReLU_on("ReLU_on");
//    ReLU_on = x.in(ReLUs) - (x.in(x_ReLUs)*B.in(B_ReLUs) + C.in(C_ReLUs)) - 1000*(1-y.in(ReLUs));
//    NN.add(ReLU_on.in(ReLUs) <= 0);
    
    Constraint<> ReLU_on("ReLU_on");
    ReLU_on = x.in(ReLUs) - (x.in(x_ReLUs)*B.in(B_ReLUs) + C.in(C_ReLUs));
    NN.add_on_off(ReLU_on.in(ReLUs) <= 0, y.in(ReLUs), true);
//
    
//    Constraint<> ReLU_off("ReLU_off");
//    ReLU_off = (B.in(B_ReLUs)*x.in(x_ReLUs) + C.in(C_ReLUs)) - 1000*y.in(ReLUs);
//    NN.add(ReLU_off.in(ReLUs) <= 0);
//

    Constraint<> ReLU_y_off("ReLU_y_off");
    ReLU_y_off = x.in(ReLUs);
    NN.add_on_off(ReLU_y_off.in(ReLUs) <= 0, y.in(ReLUs), false);


    Constraint<> Gemm("Gemm");
    Gemm = x.in(Gemms) - (x.in(x_Gemms_in)*B.in(B_Gemm) + C.in(Gemms));
    NN.add(Gemm.in(Gemms) == 0);

    Constraint<> Add("Add");
    Add = x.in(Adds) + B.in(Adds);
    NN.add(Add.in(Adds) == 0);

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
