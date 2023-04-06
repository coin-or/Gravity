#include <gravity/Node.h>
#include <gravity/Net.h>
#include <onnx.pb.h>

class NeuralNetwork: public Net {
public:
	NeuralNetwork(): Net() {}
};

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
