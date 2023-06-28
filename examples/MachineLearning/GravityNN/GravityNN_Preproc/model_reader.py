import gzip
import tempfile

import onnx2torch
import onnxruntime as ort
import onnxsim
import torch

import onnx
from onnx import version_converter


def calculate_io_shapes(proto: onnx.ModelProto):
    assert len(proto.graph.input) == 1, f"Expected 1 input, got {len(proto.graph.input)}"
    assert len(proto.graph.output) == 1, f"Expected 1 output, got {len(proto.graph.output)}"

    input_shape  = [dim.dim_value for dim in proto.graph.input[0].type.tensor_type.shape.dim]
    output_shape = [dim.dim_value for dim in proto.graph.output[0].type.tensor_type.shape.dim]

    # check for dynamic shapes
    dynamic = [s == 0 for s in input_shape]
    if any(dynamic):
        assert any(dynamic[1:]) == False
    input_shape[0] = 1 if dynamic[0] else input_shape[0]
    output_shape[0] = 1 if dynamic[0] else output_shape[0]

    return input_shape, output_shape

def load_onnx(model_path: str):
    if model_path.endswith(".gz"):
        with gzip.open(model_path, "rb") as f:
            proto = onnx.load(f)
    else:
        proto = onnx.load(model_path)

    proto = remove_initializer_from_input(proto)

    input_shape = proto.graph.input[0].type.tensor_type.shape.dim
    input_shape = [x.dim_value for x in input_shape]
    input_shape[0] = 1 if input_shape[0] == 0 else input_shape[0]

    for i, v in enumerate(input_shape):
        proto.graph.input[0].type.tensor_type.shape.dim[i].dim_value = v

    proto = ort_optimize(proto)

    if proto.opset_import[0].version < 12:
        print(f"Converting model from {proto.opset_import[0].version} to opset 12")
        proto = version_converter.convert_version(proto, 12)

    proto, check = onnxsim.simplify(proto)
    if not check:
        raise RuntimeError("Simplification failed")

    proto = ort_optimize(proto)
    return proto

def load_model(model_path: str):
    proto = load_onnx(model_path)
    input_shape, output_shape = calculate_io_shapes(proto)

    model = onnx2torch.convert(proto)

    model.eval()
    for p in model.parameters():
        p.requires_grad = False

    inp = torch.randn(input_shape)

    return model, input_shape, output_shape

def remove_initializer_from_input(model: onnx.ModelProto):
    if model.ir_version < 4:
        # tqdm.tqdm.write("Model with ir_version below 4 requires initilizers in graph input. Updating IR")
        model = onnx.helper.make_model(model.graph)

    inputs = model.graph.input
    name_to_input = {}
    for inp in inputs:
        name_to_input[inp.name] = inp

    for initializer in model.graph.initializer:
        if initializer.name in name_to_input:
            inputs.remove(name_to_input[initializer.name])
    return model

def ort_optimize(model: onnx.ModelProto):
    # Taken from onnxsim, TODO: Credit them
    model_bytes = model.SerializeToString()

    tmp_file = tempfile.NamedTemporaryFile()
    sess_options = ort.SessionOptions()
    sess_options.graph_optimization_level = ort.GraphOptimizationLevel.ORT_ENABLE_BASIC
    sess_options.optimized_model_filepath = tmp_file.name
    _ = ort.InferenceSession(model_bytes, sess_options, providers=['CPUExecutionProvider'])

    model = onnx.load(tmp_file.name)
    return model