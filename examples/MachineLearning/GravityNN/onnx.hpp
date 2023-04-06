#pragma once
#include <onnx.pb.h>
#include <cassert>
#include <fstream>

class Tensor {
public:
	Tensor(onnx::TensorProto tensor) {
		this->dtype = tensor.data_type();

		if (this->dtype == onnx::TensorProto_DataType_FLOAT) {
			for (auto f : tensor.float_data()) {
				this->data.push_back(f);
			}
			this->dtype_str = "float";
		} else if (this->dtype == onnx::TensorProto_DataType_INT64) {
			for (auto i : tensor.int64_data()) {
				this->data.push_back(i);
			}
			this->dtype_str = "int64";
		} else {
			throw std::runtime_error("Unsupported data type.");
		}

		this->name = tensor.name();

		for (auto dim : tensor.dims()) {
			this->shape.push_back(dim);
		}

	}

	std::string name;
	std::vector<float> data;
	std::vector<size_t> shape;

	size_t dtype;
	std::string dtype_str;
};

class IO {
public:
	IO(onnx::ValueInfoProto value_info_proto, bool is_input) {
		this->name = value_info_proto.name();
		this->is_input = is_input;

		const onnx::TypeProto& type_proto = value_info_proto.type();
		const onnx::TypeProto::Tensor& tensor_type = type_proto.tensor_type();
		const onnx::TensorShapeProto& shape = tensor_type.shape();
		for(int j = 0; j < shape.dim_size(); j++) {
			const onnx::TensorShapeProto::Dimension& dim = shape.dim(j);
			this->shape.push_back(dim.dim_value());
		}
	}

	std::string name;
	std::vector<size_t> shape;

	bool is_input;
};
