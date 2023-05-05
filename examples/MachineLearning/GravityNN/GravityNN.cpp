#include <iostream>
#include <fstream>
#include <onnx.pb.h>
#include <vector>
#include "Layers.hpp"
//#include "smtlib2yices.h"

std::tuple<std::vector<Layer*>, std::vector<size_t>, Tensors> build_graph(std::string fname) {
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
        } else if (node.op_type() == "MatMul") {
            layers.push_back(new MatMul(node, tensors));
        } else if (node.op_type() == "Add") {
            layers.push_back(new Add(node, tensors));
        } else if (node.op_type() == "Conv") {
            layers.push_back(new Conv(node, tensors));
        } else if (node.op_type() == "Flatten") {
            layers.push_back(new Flatten(node, tensors));
        } else {
            throw std::runtime_error("Unsupported operator " + node.op_type());
        }
    }

//    for (const auto& layer : layers) {
//        layer->print();
//    }

    std::vector<size_t> input_dims;
    if (graph.input_size() > 1) {
        throw std::runtime_error("Network has more than one input. Not supported.");
    } else {
        input_dims = tensors[graph.input(0).name()].shape;
    }
    return {layers,input_dims,tensors};
}


#include <gravity/solver.h>


using namespace gravity;
int main (int argc, char * argv[]){
    string fname_onnx = string(prj_dir)+"/data_sets/VNN/mnist_x6_bounded.onnx";
    if(argc>=2){
        fname_onnx=argv[1];
    }
    string fname_vnnlib = string(prj_dir)+"/data_sets/VNN/mnist_x6_bounded.vnnlib";
    if(argc>=2){
        fname_onnx=argv[1];
    }
//    smtlib2_yices_parser *yp = smtlib2_yices_parser_new();
//    smtlib2_abstract_parser_parse((smtlib2_abstract_parser *)yp, stdin);
//    smtlib2_yices_parser_delete(yp);

    auto tmp = build_graph(fname_onnx);
    auto layers = std::get<0>(tmp);
    auto input_dims = std::get<1>(tmp);
    auto tensors = std::get<2>(tmp);

    /* INDEX SETS */
    /* Indexing variables */
    indices x_ids("x_ids"), y_ids("y_ids");/*< x_ids for continuous vars, y_ids for binary vars */

    /* Indexing params */
    indices B_ids("B_ids"), C_ids("C_ids"), W_ids("W_ids");/*< x_ids for continuous vars, y_ids for binary vars */

    param<> B("B"), C("C"), W("W");
    B.in(B_ids);
    C.in(C_ids);
    W.in(W_ids);
    /* Creating variables' indexing sets */
    size_t n_input = accumulate(input_dims.begin(), input_dims.end(), 1, multiplies<size_t>());
    for(auto j = 0; j < n_input;j++){
        x_ids.add("input"+to_string(j));
    }
    auto nb_layers = layers.size();
    vector<size_t> dims;
    bool is_activation_func = false;
    for(auto i = 0; i<nb_layers; i++){
        auto l = layers[i];
        if(i+1<nb_layers && layers[i+1]->operator_type==_relu){/*< The next layer is an activation layer (e.g. ReLU), no need to introduce variables for this layer */
            continue;
        }
        for(auto j = 0; j < l->outputs[0].numel;j++){
            x_ids.add(l->name+","+to_string(j));
            if(l->operator_type==_relu){
                y_ids.add(l->name+","+to_string(j));
            }
        }
    }
    /* Indexing constraints */
    indices ReLUs("ReLUs"), Conv_ReLUs("Conv_ReLUs"), Gemm_ReLUs("Gemm_ReLUs"), x_Gemm_ReLUs("x_Gemm_ReLUs"), x_Conv_ReLUs("x_Conv_ReLUs"), B_Gemm_ReLUs("B_Gemm_ReLUs"), C_Gemm_ReLUs("C_Gemm_ReLUs"), B_Conv_ReLUs("B_Conv_ReLUs"), Gemms("Gemms"), x_Gemms_in("x_Gemms_in"), x_Convs_in("x_Convs_in"), x_Gemms_out("x_Gemms_out"), Adds("Adds"), x_Adds("x_Adds"), MatMuls("MatMuls"), x_MatMuls("MatMuls"), B_Gemm("B_Gemm"), W_Conv("W_Conv"), Convs("Convs");
    Conv_ReLUs = ReLUs;
    Gemm_ReLUs = ReLUs;
    B_Gemm = B_ids;
    B_Gemm_ReLUs = B_ids;
    B_Conv_ReLUs = B_ids;
    C_Gemm_ReLUs = C_ids;
    x_Gemm_ReLUs = x_ids;
    x_Conv_ReLUs = x_ids;
    x_Gemms_in = x_ids;
    x_Gemms_out = x_ids;
    x_MatMuls = x_ids;
    x_Adds = x_ids;
    size_t relu_row_id = 0, gemm_row_id = 0, conv_row_id = 0;
    string key;
    for(auto i = 0; i<nb_layers; i++){
        auto l = layers[i];
        switch (l->operator_type) {
            case _relu:
                /* Creating ReLU indexing sets */
                if(layers[i-1]->operator_type==_gemm){
                    for(auto j = 0; j < layers[i-1]->outputs[0].numel;j++){
                        ReLUs.add(l->name+","+to_string(j));
                        Gemm_ReLUs.add_ref(l->name+","+to_string(j));
                    }
                }
                else if(layers[i-1]->operator_type==_conv){
                    for(auto j = 0; j < layers[i-1]->outputs[0].numel;j++){
                        ReLUs.add(l->name+","+to_string(j));
                        Conv_ReLUs.add_ref(l->name+","+to_string(j));
                    }
                }
                else{
                    for(auto j = 0; j < layers[i-1]->var_dims[0];j++){
                        ReLUs.add(l->name+","+to_string(j));
                    }
                }
                break;
            case _conv:{
                bool add_Conv_constraint = true;
                if(i+1<nb_layers && layers[i+1]->operator_type==_relu)
                    add_Conv_constraint = false;
                auto conv = (Conv*)l;
                size_t W_idx = 0, B_idx = 0;
                for(auto h = 0; h < l->inputs[1].shape[0]; h++) {
                    for(auto i = 0; i < l->inputs[1].shape[1]; i++) {
                        for(auto j = 0; j < l->inputs[1].shape[2]; j++) {
                            for(auto k = 0; k < l->inputs[1].shape[3]; k++) {
                                key = conv->name+","+to_string(h)+","+to_string(i)+","+to_string(j)+","+to_string(k);
                                W.add_val(key, conv->inputs[1](W_idx++));
                            }
                        }
                    }
                }
                for(auto h = 0; h < l->inputs[2].shape[0]; h++) {
                    B.add_val(l->name+","+to_string(h), conv->inputs[2](B_idx++));
                    if(add_Conv_constraint){
                        Convs.add(l->name+","+to_string(h));
                    }
                }
                
                if(!add_Conv_constraint){
                    for(auto h = 0; h < l->outputs[0].numel; h++) {
                            B_Conv_ReLUs.add_ref(l->name+","+to_string(h%l->inputs[2].shape[0]));
                    }
                }
                for(auto k = 0; k < l->inputs[0].numel; k++){
                    key = conv->name+","+to_string(k);
                    if(add_Conv_constraint){
                        if(i==0){
                            x_Convs_in.add_in_row(conv_row_id, "input"+to_string(k));
                        }
                        else{
                            x_Convs_in.add_in_row(conv_row_id, layers[i-1]->name+","+to_string(k));
                        }
                    }
                }
                    if(add_Conv_constraint){
                        conv_row_id++;
                    }
//                }
//                B.print_vals(6);
                break;
            }
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
                        C_Gemm_ReLUs.add_ref(l->name+","+to_string(j));
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
        if(i+1<nb_layers && l->operator_type==_gemm && layers[i+1]->operator_type==_relu){/*< The next layer is an activation layer (e.g. ReLU), no need to introduce variables for this layer */
            for(auto j = 0; j < l->var_dims[1];j++){
                for(auto k = 0; k < l->var_dims[0];k++){
                    if(i==0)
                        x_Gemm_ReLUs.add_in_row(relu_row_id, "input"+to_string(k));
                    else
                        x_Gemm_ReLUs.add_in_row(relu_row_id, layers[i-1]->name+","+to_string(k));
                    B_Gemm_ReLUs.add_in_row(relu_row_id, l->name+","+to_string(k)+","+to_string(j));
                }
                relu_row_id++;

            }
        }
        // The equation convolves the filter with the specified input region
        // Iterate over the filter
        if(i+1<nb_layers && l->operator_type==_conv && layers[i+1]->operator_type==_relu){/*< The next layer is an activation layer (e.g. ReLU), no need to introduce variables for this layer */
            auto conv = (Conv*)l;
            string lname = l->name+",";
            if(i==0)
                lname = "input";
            for(auto i = 0; i < l->outputs[0].shape[3]; i++) {
                for(auto j = 0; j < l->outputs[0].shape[2]; j++) {
                    for(auto k = 0; k < l->outputs[0].shape[1]; k++) {
                        
                        for(auto di = 0; di < l->inputs[1].shape[3]; di++) {
                            for(auto dj = 0; dj < l->inputs[1].shape[2]; dj++) {
                                for(auto dk = 0; dk < l->inputs[1].shape[1]; dk++) {
                                    auto w_ind = (conv->strides[0]*i+di - conv->pads[0]);
                                    auto h_ind = (conv->strides[1]*j+dj - conv->pads[3]);
                                    if(h_ind < conv->inputs[0].shape[2] && h_ind >= 0 && w_ind < conv->inputs[0].shape[3] && w_ind >= 0){
                                        x_Conv_ReLUs.add_in_row(conv_row_id, lname+to_string(dk*conv->inputs[0].shape[2]*conv->inputs[0].shape[3] + w_ind*conv->inputs[0].shape[3] + h_ind));
                                        W_Conv.add_in_row(conv_row_id, conv->name+","+to_string(k)+","+to_string(dk)+","+to_string(di)+","+to_string(dj));
                                    }
                                }
                            }
                        }
                        conv_row_id++;
                    }
                }
            }
        }
    }
    Model<> NN("NN_"+fname_onnx.substr(fname_onnx.find_last_of("/")+1));
    param<> x_lb("x_lb"), x_ub("x_ub");
    param<> s_lb("s_lb"), s_ub("s_ub");
    param<int> y_lb("y_lb"), y_ub("y_ub");
    x_lb.in(x_ids);x_ub.in(x_ids);
    y_lb.in(y_ids);y_ub.in(y_ids);
    s_lb.in(y_ids);s_ub.in(y_ids);
    x_lb = numeric_limits<double>::lowest();
    x_ub = numeric_limits<double>::max();
    s_lb = 0;
    s_ub = numeric_limits<double>::max();
    y_lb = 0;
    y_ub = 1;

    
    for(auto const &key: *ReLUs._keys){
        x_lb.set_val(key, 0);
    }
    size_t nb_Relus_fixed = 0;
    for(auto i = 0; i<nb_layers; i++){
        auto l = layers[i];
        if(i+1<nb_layers && layers[i+1]->operator_type==_relu)
            continue;
        auto layer_name = l->name;
        if(l->operator_type!=_flatten & l->lowers.size()>0){
            for(auto j = 0; j < l->outputs[0].numel;j++){
                x_lb.set_val(layer_name+","+to_string(j), l->lowers[0](j));
                x_ub.set_val(layer_name+","+to_string(j), l->uppers[0](j));
                if(l->operator_type==_relu){
                    if(l->uppers[0](j)==0){
                        y_ub.set_val(layer_name+","+to_string(j), 0);
                        x_ub.set_val(layer_name+","+to_string(j), 0);
                        nb_Relus_fixed++;
                    }
                    if(l->lowers[0](j)>0){
                        y_lb.set_val(layer_name+","+to_string(j), 1);
                        s_ub.set_val(layer_name+","+to_string(j), 0);
                        nb_Relus_fixed++;
                    }
                }
            }
        }
    }
    DebugOn("Fixed " << nb_Relus_fixed << " ReLUs\n");
    if (tensors.count("Input0_lower") != 0) {
        auto lower = tensors.at("Input0_lower");
        auto upper = tensors.at("Input0_upper");
        for(auto j = 0; j < input_dims[1];j++){
            x_lb.set_val("input"+to_string(j), lower(j));
            x_ub.set_val("input"+to_string(j), upper(j));
        }
    }
//    param<> x_lb("x_lb"), x_ub("x_ub");
//    x_lb.in(x_ids);x_ub.in(x_ids);
//    x_lb = -1000;
//    x_ub = 1000;
//    for(auto j = 0; j < input_dims[1];j++){
//        x_lb.set_val(j, -2);
//        x_ub.set_val(j, 2);
//    }
//
//    for(auto const &key: *ReLUs._keys){
//        x_lb.set_val(key, 0);
//    }
//
//    for(auto i = 0; i<nb_layers; i++){
//        auto l = layers[i];
//        if(i+1<nb_layers && layers[i+1]->operator_type==_relu)
//            continue;
//        auto layer_name = l->name;
//        size_t dim = 0;
//        if(l->lowers.size()>0){
//            if(l->var_dims.size()==2){
//                dim = l->var_dims[1];
//            }
//            else{
//                dim = l->var_dims[0];
//            }
//            for(auto j = 0; j < dim;j++){
//                x_lb.set_val(layer_name+","+to_string(j), l->lowers[0](j));
//                x_ub.set_val(layer_name+","+to_string(j), l->uppers[0](j));
//            }
//        }
//    }
    var<> x("x", x_lb, x_ub);
    var<int> y("y", y_lb, y_ub);
    var<> s("s", s_lb, s_ub);
    NN.add(x.in(x_ids));
    NN.add(y.in(y_ids));

    /* Objective function */
    NN.min(x(layers.back()->name+",4") - x(layers.back()->name+",1"));
//    NN.max(x(layers.back()->name+",0"));

    
//    Constraint<> Conv_ReLU("Conv_ReLU");
//    Conv_ReLU = x.in(Conv_ReLUs) -  (x.in(x_Conv_ReLUs)*W.in(W_Conv) + B.in(B_Conv_ReLUs));
//    NN.add(Conv_ReLU.in(Conv_ReLUs) >= 0);
    
    /* Constraints */
    Constraint<> Gemm_ReLU("Gemm_ReLU");
    Gemm_ReLU = x.in(Gemm_ReLUs) - (x.in(x_Gemm_ReLUs)*B.in(B_Gemm_ReLUs) + C.in(C_Gemm_ReLUs));
    NN.add(Gemm_ReLU.in(Gemm_ReLUs) >= 0);

//    Constraint<> ReLU_on("ReLU_on");
//    ReLU_on = x.in(ReLUs) - (x.in(x_ReLUs)*B.in(B_ReLUs) + C.in(C_ReLUs)) - 1000*(1-y.in(ReLUs));
//    NN.add(ReLU_on.in(ReLUs) <= 0);
    
    Constraint<> ReLU_on("ReLU_on");
    ReLU_on = x.in(Gemm_ReLUs) - (x.in(x_Gemm_ReLUs)*B.in(B_Gemm_ReLUs) + C.in(C_Gemm_ReLUs));
    NN.add_on_off(ReLU_on.in(Gemm_ReLUs) <= 0, y.in(Gemm_ReLUs), true);
//
    
//    Constraint<> ReLU_off("ReLU_off");
//    ReLU_off = (B.in(B_ReLUs)*x.in(x_ReLUs) + C.in(C_ReLUs)) - 1000*y.in(ReLUs);
//    NN.add(ReLU_off.in(ReLUs) <= 0);
//

    Constraint<> ReLU_y_off("ReLU_y_off");
    ReLU_y_off = x.in(Gemm_ReLUs);
    NN.add_on_off(ReLU_y_off.in(Gemm_ReLUs) <= 0, y.in(Gemm_ReLUs), false);
    
//    Constraint<> ReLU_off("ReLU_off");
//    ReLU_off = x.in(Gemm_ReLUs) - 1e-4;
//    NN.add_on_off(ReLU_off.in(Gemm_ReLUs) >= 0, y.in(Gemm_ReLUs), true);


    Constraint<> Gemm("Gemm");
    Gemm = x.in(Gemms) - (x.in(x_Gemms_in)*B.in(B_Gemm) + C.in(Gemms));
    NN.add(Gemm.in(Gemms) == 0);

//    Constraint<> Conv("Conv");
//    Conv = x.in(Convs) - (x.in(x_Convs_in)*W.in(W_Conv) + B.in(Convs));
//    NN.add(Conv.in(Convs) == 0);

    
//    Constraint<> Add("Add");
//    Add = x.in(Adds) + B.in(Adds);
//    NN.add(Add.in(Adds) == 0);

//    NN.print();
    NN.write();
    

//    solver<> S(NN,gurobi);
//    S.run(1e-4, 1800);
//    NN.print_solution();
    bool build_NLP=true;
    if(build_NLP){
        Model<> NLP("NLP_"+fname_onnx.substr(fname_onnx.find_last_of("/")+1));
        
        NLP.add(x.in(x_ids));
//        NLP.add(y.in(y_ids));
        NLP.add(s.in(y_ids));

        /* Objective function */
        NLP.min(x(layers.back()->name+",4") - x(layers.back()->name+",0"));
    //    NLP.max(x(layers.back()->name+",0"));

        
    //    Constraint<> Conv_ReLU("Conv_ReLU");
    //    Conv_ReLU = x.in(Conv_ReLUs) -  (x.in(x_Conv_ReLUs)*W.in(W_Conv) + B.in(B_Conv_ReLUs));
    //    NLP.add(Conv_ReLU.in(Conv_ReLUs) >= 0);
        
        /* Constraints */
        Constraint<> Gemm_ReLU("Gemm_ReLU");
        Gemm_ReLU = x.in(Gemm_ReLUs) - s.in(Gemm_ReLUs) - (x.in(x_Gemm_ReLUs)*B.in(B_Gemm_ReLUs) + C.in(C_Gemm_ReLUs));
        NLP.add(Gemm_ReLU.in(Gemm_ReLUs) == 0);
        
        Constraint<> Complement("Complement");
        Complement = x.in(ReLUs)*s.in(ReLUs);
        NLP.add(Complement.in(ReLUs) <= 0);

    //    Constraint<> ReLU_on("ReLU_on");
    //    ReLU_on = x.in(ReLUs) - (x.in(x_ReLUs)*B.in(B_ReLUs) + C.in(C_ReLUs)) - 1000*(1-y.in(ReLUs));
    //    NLP.add(ReLU_on.in(ReLUs) <= 0);
        
//        Constraint<> ReLU_on("ReLU_on");
//        ReLU_on = x.in(Gemm_ReLUs) - (x.in(x_Gemm_ReLUs)*B.in(B_Gemm_ReLUs) + C.in(C_Gemm_ReLUs));
//        NLP.add_on_off(ReLU_on.in(Gemm_ReLUs) <= 0, y.in(Gemm_ReLUs), true);
    //
        
    //    Constraint<> ReLU_off("ReLU_off");
    //    ReLU_off = (B.in(B_ReLUs)*x.in(x_ReLUs) + C.in(C_ReLUs)) - 1000*y.in(ReLUs);
    //    NLP.add(ReLU_off.in(ReLUs) <= 0);
    //

//        Constraint<> ReLU_y_off("ReLU_y_off");
//        ReLU_y_off = x.in(Gemm_ReLUs);
//        NLP.add_on_off(ReLU_y_off.in(Gemm_ReLUs) <= 0, y.in(Gemm_ReLUs), false);


        Constraint<> Gemm("Gemm");
        Gemm = x.in(Gemms) - (x.in(x_Gemms_in)*B.in(B_Gemm) + C.in(Gemms));
        NLP.add(Gemm.in(Gemms) == 0);
        NLP.write();
        solver<> S(NLP,ipopt);
        S.run(1e-4, 1800);
        NLP.print_solution();
        size_t nb_z_comp = 0;
        for (auto i = 0; i<y_ids.size(); i++) {
            auto key = y_ids._keys->at(i);
            if(s.eval(i)<=1e-6 && std::abs(x.eval(key)) <= 1e-6)
                nb_z_comp++;
        }
        DebugOn("Number of zero complementarities = " << nb_z_comp << endl);
    }
    
    for(auto l:layers){
        delete l;
    }


    return 0;
}
