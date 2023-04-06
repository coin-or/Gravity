#include <gravity/Node.h>
#include <gravity/Net.h>
#include <onnx.pb.h>


class Layer: public Node {
public:
    Layer(
        int idx,
        onnx::NodeProto node,
        std::map<std::string, onnx::TensorProto> global_initializers
    ): Node(node.name(), idx) {
        this->_type_name = node.op_type();

        for (auto input: node.input()) {
            if (global_initializers.count(input) == 0) {
                continue;
            }
            this->initializers[input] = global_initializers[input];
        }
    }

    std::vector<std::string> inputs() {
        std::vector<std::string> inputs;
        for (auto& input: this->node.input()) {
            if (this->initializers.count(input) == 0) {
                inputs.push_back(input);
            }
        }
        return inputs;
    }

    std::vector<std::string> outputs() {
        return std::vector<std::string>{this->node.output().begin(), this->node.output().end()};
    }

    onnx::NodeProto node;
    std::map<std::string, onnx::TensorProto> initializers;
};


class NeuralNetwork: public Net {
public:
	NeuralNetwork(): Net() {}
    
    void readONNX(string fname){
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
        

        int layer_idx = 0;
        for (auto& node: graph.node()) {
            Layer* layer = new Layer(layer_idx, node, initializers);
            layers.push_back(layer);
            for (auto& output: layer->outputs()) {
                layers_by_output[output] = layer;
            }
            layer_idx++;

            add_node(layer);
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

                add_arc(arc);
            }
        }
    }
};


