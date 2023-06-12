#include <iostream>
#include <fstream>
#include <onnx.pb.h>
#include <vector>
#include "Layers.hpp"
//#include "smtlib2yices.h"
using namespace gravity;

shared_ptr<Model<>> build_NN_MIP(const indices& ReLUs, const indices& x_ids, const indices& z_ids, const indices& B_ids,const indices& C_ids, const vector<Layer*>& layers, const param<>& x_lb, const param<>& x_ub, const param<int>& z_lb, const param<int>& z_ub, param<>& B, param<>& C, int relu_id_min, int relu_id_max, int nb_relus){
    auto nb_layers = layers.size();
    int min_i = 0, max_i = nb_layers, relu_id = 0;
    for(auto i = 0; i<nb_layers; i++){
        auto l = layers[i];
        if(l->operator_type==_relu){
            if(relu_id<relu_id_min)
                min_i = i;
            if(relu_id>=relu_id_max){
                max_i = i;
                break;
            }
            relu_id++;
        }
    }
    if(relu_id_max == nb_relus-1)
        max_i = layers.size();
    
    /* Indexing constraints */
    indices Conv_ReLUs("Conv_ReLUs"), Gemm_ReLUs("Gemm_ReLUs"), x_Gemm_ReLUs("x_Gemm_ReLUs"), x_Conv_ReLUs("x_Conv_ReLUs"), B_Gemm_ReLUs("B_Gemm_ReLUs"), C_Gemm_ReLUs("C_Gemm_ReLUs"), B_Conv_ReLUs("B_Conv_ReLUs"), Gemms("Gemms"), x_Gemms_in("x_Gemms_in"), x_Convs_in("x_Convs_in"), x_Gemms_out("x_Gemms_out"), x_Adds("x_Adds"), MatMuls("MatMuls"), x_MatMuls("MatMuls"), B_Gemm("B_Gemm"), W_Conv("W_Conv"), Convs("Convs");
    Conv_ReLUs = ReLUs;
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
    string key;
    size_t relu_row_id = 0, gemm_row_id = 0, conv_row_id = 0;
        for(auto i = min_i; i<max_i; i++){
            auto l = layers[i];
            switch (l->operator_type) {
                case _relu:
                    /* Creating ReLU indexing sets */
//                    if(layers[i-1]->operator_type==_gemm){
//                        for(auto j = 0; j < layers[i-1]->outputs[0].numel;j++){
//                            Gemm_ReLUs.add(l->name+","+to_string(j));
//                        }
//                    }
//                    else if(layers[i-1]->operator_type==_conv){
//                        for(auto j = 0; j < layers[i-1]->outputs[0].numel;j++){
//                            Conv_ReLUs.add_ref(l->name+","+to_string(j));
//                        }
//                    }
                    break;
                case _conv:{
                    bool add_Conv_constraint = true;
                    if(i+1<nb_layers && layers[i+1]->operator_type==_relu)
                        add_Conv_constraint = false;
                    auto conv = (Conv*)l;
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
                    break;
                }
                case _gemm:{
                    bool add_Gemm_constraint = true;
                    if(i+1<nb_layers && layers[i+1]->operator_type==_relu)
                        add_Gemm_constraint = false;
                    auto gemm = (GEMM*)l;
                    if(add_Gemm_constraint){
                        for(auto j = 0; j < l->var_dims[1];j++){
                            Gemms.add(l->name+","+to_string(j));
                            for(auto k = 0; k < l->var_dims[0]; k++){
                                key = gemm->name+","+to_string(k)+","+to_string(j);
                                B_Gemm.add_in_row(gemm_row_id, key);
                                if(i==0){
                                    x_Gemms_in.add_in_row(gemm_row_id, "input"+to_string(k));
                                }
                                else{
                                    x_Gemms_in.add_in_row(gemm_row_id, layers[i-1]->name+","+to_string(k));
                                }
                            }
                            gemm_row_id++;
                        }
                    }

                }
                    break;
                default:
                    break;
            }
        }

    for(auto i = min_i; i<max_i; i++){
        auto l = layers[i];
        if(i+1<nb_layers && l->operator_type==_gemm && layers[i+1]->operator_type==_relu){/*< The next layer is an activation layer (e.g. ReLU), no need to introduce variables for this layer */
            for(auto j = 0; j < l->var_dims[1];j++){
                Gemm_ReLUs.add(layers[i+1]->name+","+to_string(j));
                for(auto k = 0; k < l->var_dims[0];k++){
                    if(i==0)
                        x_Gemm_ReLUs.add_in_row(relu_row_id, "input"+to_string(k));
                    else
                        x_Gemm_ReLUs.add_in_row(relu_row_id, layers[i-1]->name+","+to_string(k));
                    B_Gemm_ReLUs.add_in_row(relu_row_id, l->name+","+to_string(k)+","+to_string(j));
                }
                C_Gemm_ReLUs.add_ref(l->name+","+to_string(j));
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
    
    
    
    shared_ptr<Model<>> NN = make_shared<Model<>>(("NN_"+to_string(relu_id_min)+"_"+to_string(relu_id_max)));
    

    var<> x("x", x_lb, x_ub);
    var<> z("z", z_lb, z_ub);
    z._is_relaxed = true;
    NN->add(x.in(x_ids));
    NN->add(z.in(z_ids));

    /* Objective function */
    NN->min(x(layers.back()->name+",4") - x(layers.back()->name+",1"));
//    NN.max(x(layers.back()->name+",0"));

    
//    Constraint<> Conv_ReLU("Conv_ReLU");
//    Conv_ReLU = x.in(Conv_ReLUs) -  (x.in(x_Conv_ReLUs)*W.in(W_Conv) + B.in(B_Conv_ReLUs));
//    NN.add(Conv_ReLU.in(Conv_ReLUs) >= 0);
    
    /* Constraints */
    Constraint<> Gemm_ReLU("Gemm_ReLU");
    Gemm_ReLU = x.in(Gemm_ReLUs) - (x.in(x_Gemm_ReLUs)*B.in(B_Gemm_ReLUs) + C.in(C_Gemm_ReLUs));
    NN->add(Gemm_ReLU.in(Gemm_ReLUs) >= 0);

//    Constraint<> ReLU_on("ReLU_on");
//    ReLU_on = x.in(ReLUs) - (x.in(x_ReLUs)*B.in(B_ReLUs) + C.in(C_ReLUs)) - 1000*(1-y.in(ReLUs));
//    NN.add(ReLU_on.in(ReLUs) <= 0);
    
    Constraint<> ReLU_on("ReLU_on");
    ReLU_on = x.in(Gemm_ReLUs) - (x.in(x_Gemm_ReLUs)*B.in(B_Gemm_ReLUs) + C.in(C_Gemm_ReLUs));
    NN->add_on_off(ReLU_on.in(Gemm_ReLUs) <= 0, z.in(Gemm_ReLUs), true);
    
//    Constraint<> ReLU_off("ReLU_off");
//    ReLU_off = (B.in(B_ReLUs)*x.in(x_ReLUs) + C.in(C_ReLUs)) - 1000*y.in(ReLUs);
//    NN.add(ReLU_off.in(ReLUs) <= 0);
//

    Constraint<> ReLU_y_off("ReLU_y_off");
    ReLU_y_off = x.in(Gemm_ReLUs);
    NN->add_on_off(ReLU_y_off.in(Gemm_ReLUs) <= 0, z.in(Gemm_ReLUs), false);
    
//    Constraint<> ReLU_off("ReLU_off");
//    ReLU_off = x.in(Gemm_ReLUs) - 1e-4;
//    NN.add_on_off(ReLU_off.in(Gemm_ReLUs) >= 0, y.in(Gemm_ReLUs), true);


    Constraint<> Gemm("Gemm");
    Gemm = x.in(Gemms) - (x.in(x_Gemms_in)*B.in(B_Gemm) + C.in(Gemms));
    NN->add(Gemm.in(Gemms) == 0);
    return NN;
}

std::tuple<std::vector<Layer*>, std::vector<size_t>, Tensors> build_graph(std::string fname) {
    
	std::fstream input(fname, std::ios::in | std::ios::binary);
    std::stringstream buffer;
    buffer << input.rdbuf();
    std::string onnx_str = buffer.str();

    onnx::ModelProto model;
    bool isSuccess = model.ParseFromString(onnx_str);
    if (!isSuccess) {
        throw std::runtime_error("Failed to parse onnx file.");
    }
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
    indices x_ids("x_ids"), z_ids("z_ids");/*< x_ids for continuous vars, z_ids for binary vars */

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
    int nb_relus = 0;
    vector<size_t> dims;
    bool is_activation_func = false;
    for(auto i = 0; i<nb_layers; i++){
        auto l = layers[i];
        if(i+1<nb_layers && layers[i+1]->operator_type==_relu){/*< The next layer is an activation layer (e.g. ReLU), no need to introduce variables for this layer */
            nb_relus++;
            continue;
        }
        for(auto j = 0; j < l->outputs[0].numel;j++){
            x_ids.add(l->name+","+to_string(j));
            if(l->operator_type==_relu){
                z_ids.add(l->name+","+to_string(j));
            }
        }
    }
    
    param<> x_lb("x_lb"), x_ub("x_ub");
    param<> s_lb("s_lb"), s_ub("s_ub");
    param<int> z_lb("z_lb"), z_ub("z_ub");
    x_lb.in(x_ids);x_ub.in(x_ids);
    z_lb.in(z_ids);z_ub.in(z_ids);
    s_lb.in(z_ids);s_ub.in(z_ids);
    x_lb = numeric_limits<double>::lowest();
    x_ub = numeric_limits<double>::max();
    s_lb = 0;
    s_ub = numeric_limits<double>::max();
    z_lb = 0;
    z_ub = 1;
    indices ReLUs("ReLUs"), Adds("Adds");
    string key;
    for(auto i = 0; i<nb_layers; i++){
        auto l = layers[i];
        switch (l->operator_type) {
            case _relu:
                for(auto j = 0; j < l->outputs[0].numel;j++){
                    ReLUs.add(l->name+","+to_string(j));
                }
                break;
            case _conv:{
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
                }
                break;
            }
            case _gemm:{
                auto gemm = (GEMM*)l;
                size_t B_idx = 0, C_idx = 0;
                for(auto j = 0; j < l->var_dims[0];j++){
                    for(auto k = 0; k < l->var_dims[1]; k++){
                        key = gemm->name+","+to_string(j)+","+to_string(k);
                        B.add_val(key, gemm->B(B_idx++));
                    }
                }
                for(auto j = 0; j < l->var_dims[1];j++){
                    if(gemm->has_optional_C){
                        C.add_val(l->name+","+to_string(j), gemm->C(C_idx++));
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
                        z_ub.set_val(layer_name+","+to_string(j), 0);
                        x_ub.set_val(layer_name+","+to_string(j), 0);
                        nb_Relus_fixed++;
                    }
                    if(l->lowers[0](j)>0){
                        z_lb.set_val(layer_name+","+to_string(j), 1);
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
    
    
//    int relu_id_min = 0, relu_id_max = 1;
//    auto NN = build_NN_MIP(ReLUs, x_ids, y_ids, B_ids, C_ids, layers, x_lb, x_ub, y_lb, y_ub, B, C, relu_id_min, relu_id_max, nb_relus);
//
////    NN->print();
//    NN->write();
//
//    DebugOn("Total number of relu layers = " << nb_relus << endl);
//    int horizon_len = 3;
//    DebugOn("Creating rolling horizon of size " << horizon_len << endl);
//    /* Split nb_layers into horizon_len parts */
//    vector<size_t> limits = bounds(horizon_len, nb_relus);
////    NN->replace_integers();
//    double max_time = 600;
//    unsigned max_iter=1e3;
//    double rel_tol=1e-2;
//    double abs_tol=1e6;
//    unsigned nb_threads = 2;
//    double ub_solver_tol=1e-4;
//    double lb_solver_tol=1e-4;
//    double range_tol=0;
//    SolverType ub_solver_type = gurobi;
//    SolverType lb_solver_type = gurobi;
//    auto x_NN = NN->get_ptr_var<double>("x");
//    for(auto i = 0; i < x_ids.size(); i++){
//        key = x_ids._keys->at(i);
//        if(key.find("Relu1")==string::npos)
//            x_NN->_off[i] = true;
//    }
//    auto y_NN = NN->get_ptr_var<double>("y");
//    for(auto i = 0; i < y_ids.size(); i++){
//        y_NN->_off[i] = true;
//    }
//    NN->run_obbt(NN, max_time, max_iter, rel_tol, abs_tol, nb_threads, ub_solver_type, lb_solver_type, ub_solver_tol, lb_solver_tol, range_tol);
//    DebugOn("Done running OBBT on NN1\n");
////    solver<> S(NN,gurobi);
////    S.run(1e-4, 30);
////    NN->print();
//    auto x_NN_lb = x_NN->get_lb();
//    auto x_NN_ub = x_NN->get_ub();
//    x_lb.copy_vals(x_NN_lb);
//    x_ub.copy_vals(x_NN_ub);
////    x_ub.print_vals(6);
////    x_ub("Relu1,0").print_vals(6);
//
//    auto NN2 = build_NN_MIP(ReLUs, x_ids, y_ids, B_ids, C_ids, layers, x_lb, x_ub, y_lb, y_ub, B, C, 1, 2, nb_relus);
//    NN2->write();
//    auto x_NN2 = NN2->get_ptr_var<double>("x");
//    for(auto i = 0; i < x_ids.size(); i++){
//        key = x_ids._keys->at(i);
//        if(key.find("Relu2")==string::npos)
//            x_NN2->_off[i] = true;
//    }
//    auto y_NN2 = NN2->get_ptr_var<double>("y");
//    for(auto i = 0; i < y_ids.size(); i++){
//        y_NN2->_off[i] = true;
//    }
//
//    DebugOn("running OBBT on NN2\n");
//    NN2->run_obbt(NN2, max_time, max_iter, rel_tol, abs_tol, nb_threads, ub_solver_type, lb_solver_type, ub_solver_tol, lb_solver_tol, range_tol);
//    DebugOn("Done running OBBT on NN2\n");
//    auto x_NN2_lb = x_NN2->get_lb();
//    auto x_NN2_ub = x_NN2->get_ub();
//
//    x_lb.copy_vals(x_NN2_lb);
//    x_ub.copy_vals(x_NN2_ub);
//
//    auto NN3 = build_NN_MIP(ReLUs, x_ids, y_ids, B_ids, C_ids, layers, x_lb, x_ub, y_lb, y_ub, B, C, 2, 3, nb_relus);
//    NN3->write();
//    auto x_NN3 = NN3->get_ptr_var<double>("x");
//    for(auto i = 0; i < x_ids.size(); i++){
//        key = x_ids._keys->at(i);
//        if(key.find("Relu3")==string::npos && key.find("Gemm4")==string::npos)
//            x_NN3->_off[i] = true;
//    }
//    auto y_NN3 = NN3->get_ptr_var<double>("y");
//    for(auto i = 0; i < y_ids.size(); i++){
//        y_NN3->_off[i] = true;
//    }
//
//    DebugOn("running OBBT on NN3\n");
//    NN3->run_obbt(NN3, max_time, max_iter, rel_tol, abs_tol, nb_threads, ub_solver_type, lb_solver_type, ub_solver_tol, lb_solver_tol, range_tol);
//    DebugOn("Done running OBBT on NN3\n");
//    auto x_NN3_lb = x_NN3->get_lb();
//    auto x_NN3_ub = x_NN3->get_ub();
//
//    x_lb.copy_vals(x_NN3_lb);
//    x_ub.copy_vals(x_NN3_ub);
//    x_ub.print_vals(6);
//    x_ub("Relu1,0").print_vals(6);
    auto NN_full = build_NN_MIP(ReLUs, x_ids, z_ids, B_ids, C_ids, layers, x_lb, x_ub, z_lb, z_ub, B, C, 0, nb_relus-1, nb_relus);
    NN_full->write();
//
//
//    auto x_star = NN->get_ptr_var<double>("x");
//    auto x = NN2->get_ptr_var<double>("x");
//    for (auto i = 0; i < 256; i++){
//        auto key = "Relu1,"+to_string(i);
//        x->set_lb(key, x_star->eval(key));
//        x->set_ub(key, x_star->eval(key));
//    }
//    NN2->write();
    solver<> S2(NN_full,gurobi);
    S2.run(1e-4, 180);
    
    NN_full->print_solution();
    bool build_NLP=false;
    if(build_NLP){
        
        /* Indexing constraints */
        indices Conv_ReLUs("Conv_ReLUs"), Gemm_ReLUs("Gemm_ReLUs"), x_Gemm_ReLUs("x_Gemm_ReLUs"), x_Conv_ReLUs("x_Conv_ReLUs"), B_Gemm_ReLUs("B_Gemm_ReLUs"), C_Gemm_ReLUs("C_Gemm_ReLUs"), B_Conv_ReLUs("B_Conv_ReLUs"), Gemms("Gemms"), x_Gemms_in("x_Gemms_in"), x_Convs_in("x_Convs_in"), x_Gemms_out("x_Gemms_out"), x_Adds("x_Adds"), MatMuls("MatMuls"), x_MatMuls("MatMuls"), B_Gemm("B_Gemm"), W_Conv("W_Conv"), Convs("Convs");
        Conv_ReLUs = ReLUs;
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
        string key;
        size_t relu_row_id = 0, gemm_row_id = 0, conv_row_id = 0;
            for(auto i = 0; i<nb_layers; i++){
                auto l = layers[i];
                switch (l->operator_type) {
                    case _relu:
                        /* Creating ReLU indexing sets */
                        if(layers[i-1]->operator_type==_gemm){
                            for(auto j = 0; j < layers[i-1]->outputs[0].numel;j++){
                                Gemm_ReLUs.add(l->name+","+to_string(j));
                            }
                        }
                        else if(layers[i-1]->operator_type==_conv){
                            for(auto j = 0; j < layers[i-1]->outputs[0].numel;j++){
                                Conv_ReLUs.add_ref(l->name+","+to_string(j));
                            }
                        }
                        break;
                    case _conv:{
                        bool add_Conv_constraint = true;
                        if(i+1<nb_layers && layers[i+1]->operator_type==_relu)
                            add_Conv_constraint = false;
                        auto conv = (Conv*)l;
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
                        break;
                    }
                    case _gemm:{
                        bool add_Gemm_constraint = true;
                        if(i+1<nb_layers && layers[i+1]->operator_type==_relu)
                            add_Gemm_constraint = false;
                        auto gemm = (GEMM*)l;
                        for(auto j = 0; j < l->var_dims[1];j++){
                            if(add_Gemm_constraint){
                                Gemms.add(l->name+","+to_string(j));
                                for(auto k = 0; k < l->var_dims[0]; k++){
                                    key = gemm->name+","+to_string(k)+","+to_string(j);
                                    B_Gemm.add_in_row(gemm_row_id, key);
                                    if(i==0){
                                        x_Gemms_in.add_in_row(gemm_row_id, "input"+to_string(k));
                                    }
                                    else{
                                        x_Gemms_in.add_in_row(gemm_row_id, layers[i-1]->name+","+to_string(k));
                                    }
                                }
                                gemm_row_id++;
                            }
                        }

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
                    C_Gemm_ReLUs.add_ref(l->name+","+to_string(j));
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
        
        Model<> NLP("NLP_"+fname_onnx.substr(fname_onnx.find_last_of("/")+1));
        var<> x("x", x_lb, x_ub);
//        var<int> y("y", y_lb, y_ub);
        var<> s("s", s_lb, s_ub);
        
        NLP.add(x.in(x_ids));
//        NLP.add(y.in(y_ids));
        NLP.add(s.in(z_ids));

        /* Objective function */
        NLP.min(x(layers.back()->name+",4") - x(layers.back()->name+",0"));

        /* Constraints */
        Constraint<> Gemm_ReLU("Gemm_ReLU");
        Gemm_ReLU = x.in(Gemm_ReLUs) - s.in(Gemm_ReLUs) - (x.in(x_Gemm_ReLUs)*B.in(B_Gemm_ReLUs) + C.in(C_Gemm_ReLUs));
        NLP.add(Gemm_ReLU.in(Gemm_ReLUs) == 0);

        Constraint<> Complement("Complement");
        Complement = x.in(ReLUs)*s.in(ReLUs);
        NLP.add(Complement.in(ReLUs) <= 0);

        Constraint<> Gemm("Gemm");
        Gemm = x.in(Gemms) - (x.in(x_Gemms_in)*B.in(B_Gemm) + C.in(Gemms));
        NLP.add(Gemm.in(Gemms) == 0);
        NLP.write();
        solver<> S(NLP,ipopt);

        bool run_NLP = true;
        while(run_NLP){
            S.run(5, 1e-6, "ma57", 100);
//            NLP.print_solution();



            param<> x_star("x_star");
            x_star.in(x_ids);
            x_star.copy_vals(x);
            param<> s_star("s_star");
            s_star.in(z_ids);
            s_star.copy_vals(s);
            
            indices I_s("I_s"), I_y("I_y"), I_free("I_free");
            indices x_Gemm_ReLUs_s("x_Gemm_ReLUs_s"), x_Gemm_ReLUs_y("x_Gemm_ReLUs_y"), x_Gemm_ReLUs_free("x_Gemm_ReLUs_free");
            x_Gemm_ReLUs_s = x_ids;
            x_Gemm_ReLUs_y = x_ids;
            x_Gemm_ReLUs_free = x_ids;
            indices B_Gemm_ReLUs_s("B_Gemm_ReLUs_s"), B_Gemm_ReLUs_y("B_Gemm_ReLUs_y"), B_Gemm_ReLUs_free("B_Gemm_ReLUs_free");
            B_Gemm_ReLUs_s = B_ids;
            B_Gemm_ReLUs_y = B_ids;
            B_Gemm_ReLUs_free = B_ids;
            size_t nb_z_comp = 0;
            for (auto i = 0; i<z_ids.size(); i++) {
                auto key = z_ids._keys->at(i);
                if(s.eval(i)<=1e-3 && std::abs(x.eval(key)) <= 1e-3){
                    nb_z_comp++;
                    DebugOn("zero complementarities at " << key << endl);
                    I_free.add(key);
                    x_Gemm_ReLUs_free.add_empty_row();
                    B_Gemm_ReLUs_free.add_empty_row();
                    x_Gemm_ReLUs_free._ids->back() = (x_Gemm_ReLUs._ids->at(i));
                    B_Gemm_ReLUs_free._ids->back() = (B_Gemm_ReLUs._ids->at(i));
                }
                else{
                    if(s.eval(i)<=1e-3){
                        I_s.add(key);
                        x_Gemm_ReLUs_s.add_empty_row();
                        B_Gemm_ReLUs_s.add_empty_row();
                        x_Gemm_ReLUs_s._ids->back() = (x_Gemm_ReLUs._ids->at(i));
                        B_Gemm_ReLUs_s._ids->back() = (B_Gemm_ReLUs._ids->at(i));
                    }
                    else{
                        I_y.add(key);
                        x_Gemm_ReLUs_y.add_empty_row();
                        B_Gemm_ReLUs_y.add_empty_row();
                        x_Gemm_ReLUs_y._ids->back() = (x_Gemm_ReLUs._ids->at(i));
                        B_Gemm_ReLUs_y._ids->back() = (B_Gemm_ReLUs._ids->at(i));
                    }
                }
            }
            DebugOn("I_s size = " << I_s.size() << endl);
            DebugOn("I_y size = " << I_y.size() << endl);
            DebugOn("I_free size = " << I_free.size() << endl);
            DebugOn("z_ids size = " << z_ids.size() << endl);
            
            DebugOn("Number of zero complementarities = " << nb_z_comp << endl);
            bool solve_mip = true;
            if(solve_mip){
                Model<> MIP("MIP_"+fname_onnx.substr(fname_onnx.find_last_of("/")+1));
                var<> d_x("d_x", x_lb - x_star, x_ub - x_star), d_s("d_s", s_lb - s_star, s_ub - s_star);
                var<int> z("z", 0, 1);
                MIP.add(d_x.in(x_ids), d_s.in(z_ids));
                MIP.add(z.in(I_free));
//                for(auto const &key: *I_free._keys){
//                    d_x.set_lb(key, 0);
//                    d_s.set_lb(key, 0);
//                }

                /* Objective function */
                MIP.min(d_x(layers.back()->name+",4") - d_x(layers.back()->name+",0"));

                /* Constraints */
                
                Constraint<> Gemm_ReLU("Gemm_ReLU");
                Gemm_ReLU = d_x.in(Gemm_ReLUs) - d_s.in(Gemm_ReLUs) - (d_x.in(x_Gemm_ReLUs)*B.in(B_Gemm_ReLUs));
                MIP.add(Gemm_ReLU.in(Gemm_ReLUs) == 0);
                
                Constraint<> Gemm_ReLU_s("Gemm_ReLU_s");
                Gemm_ReLU_s = d_s.in(I_s);
                MIP.add(Gemm_ReLU_s.in(I_s) == 0);
                
                Constraint<> Gemm_ReLU_y("Gemm_ReLU_y");
                Gemm_ReLU_y = d_x.in(I_y);
                MIP.add(Gemm_ReLU_y.in(I_y) == 0);

                
//                Constraint<> Gemm_ReLU_s("Gemm_ReLU_s");
//                Gemm_ReLU_s = d_x.in(I_s) - (d_x.in(x_Gemm_ReLUs_s)*B.in(B_Gemm_ReLUs_s));
//                MIP.add(Gemm_ReLU_s.in(I_s) == 0);
//
//                Constraint<> Gemm_ReLU_y("Gemm_ReLU_y");
//                Gemm_ReLU_y = -1*d_s.in(I_y) - (d_x.in(x_Gemm_ReLUs_y)*B.in(B_Gemm_ReLUs_y));
//                MIP.add(Gemm_ReLU_y.in(I_y) == 0);
//
//
//                Constraint<> Gemm_ReLU_free("Gemm_ReLU_free");
//                Gemm_ReLU_free = d_x.in(I_free) - d_s.in(I_free) - (d_x.in(x_Gemm_ReLUs_free)*B.in(B_Gemm_ReLUs_free));
//                MIP.add(Gemm_ReLU_free.in(I_free) == 0);

                Constraint<> Complement_x("Complement_x");
                Complement_x = d_x.in(I_free);
                MIP.add_on_off(Complement_x.in(I_free) <= 0, z, true);/*< if z = 1 d_x is zero */

                Constraint<> Complement_s("Complement_s");
                Complement_s = d_s.in(I_free);
                MIP.add_on_off(Complement_s.in(I_free) <= 0, z, false);/*< if z = 0 d_s is zero */



                Constraint<> Gemm("Gemm");
                Gemm = d_x.in(Gemms) - (d_x.in(x_Gemms_in)*B.in(B_Gemm));
                MIP.add(Gemm.in(Gemms) == 0);

                MIP.write();
                solver<> S(MIP,gurobi);
                S.run();
//                MIP.print_solution();
//                MIP.print_constraints_stats(1e-6);
                if(std::abs(MIP.get_obj_val())<= 1e-8) {
                    run_NLP = false;
                }
                else{
//                    for(auto const &key: *I_s._keys){
//                        d_s.param<>::set_val(key, 0);
//                    }
//                    for(auto const &key: *I_y._keys){
//                        d_x.param<>::set_val(key, 0);
//                    }
//                    d_x.print_vals(6);
//                    d_s.print_vals(6);
                    
                    
                    auto f = x + d_x;
                    f.eval_all();
                    auto f_s = s + d_s;
                    f_s.eval_all();
                    x.copy_vals(f);
                    s.copy_vals(f_s);
                    
                    x.print_vals(6);
                    NLP.reset_constrs();
                    NLP.print_constraints_stats(1e-6);
//                    break;
                }
            }
        }

    }
    
    for(auto l:layers){
        delete l;
    }


    return 0;
}
