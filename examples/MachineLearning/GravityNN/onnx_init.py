from collections import defaultdict
from copy import deepcopy
from pathlib import Path

import onnx.helper as helper
import onnxruntime as ort
import pandas as pd
import torch
from fire import Fire

import onnx
from GravityNN_Preproc.model_reader import load_onnx
from GravityNN_Preproc.vnnlib_parser import read_vnnlib


def tensor_to_onnx(tensor: torch.Tensor, name: str):
    npten = tensor.detach().cpu().numpy()
    return onnx.helper.make_tensor(
        name=name,
        data_type=helper.np_dtype_to_tensor_dtype(npten.dtype),
        dims=npten.shape,
        vals=npten.flatten(),
    )

def get_forward_initialization(model: onnx.ModelProto, input_: torch.Tensor):
    fmodel = deepcopy(model)
    # Add all intermediate outputs to graph outputs
    out_names = []
    for node in fmodel.graph.node:
        for out in node.output:
            out_names.append(helper.make_tensor_value_info(out, onnx.TensorProto.FLOAT, None))

    fmodel.graph.output.extend(out_names)
    
    inp = {fmodel.graph.input[0].name: input_.numpy()}
    sess = ort.InferenceSession(fmodel.SerializeToString())
    out = sess.run(None, inp)

    out = {name.name: ten for name, ten in zip(fmodel.graph.output, out)}

    # Convert to tensors
    inits = []
    for name, ten in out.items():
        inits.append(tensor_to_onnx(torch.from_numpy(ten), f"{name}_forward"))
    inits.append(tensor_to_onnx(input_, f"{fmodel.graph.input[0].name}_forward"))

    init_model = deepcopy(model)
    # add tensors to graph initializer
    init_model.graph.initializer.extend(inits)
    return init_model

def add_input_bounds(model: onnx.ModelProto, lb: torch.Tensor, ub: torch.Tensor):
    input_name = model.graph.input[0].name
    lb_proto = tensor_to_onnx(lb, f"{input_name}_lower")
    ub_proto = tensor_to_onnx(ub, f"{input_name}_upper")
    model.graph.initializer.extend([lb_proto, ub_proto])
    return model

def run_shape_inference(model: onnx.ModelProto):
    model = onnx.shape_inference.infer_shapes(model)
    return model

def add_specification(model: onnx.ModelProto, spec_matrix: torch.Tensor, values: torch.Tensor):
    spec_matrix_proto = tensor_to_onnx(spec_matrix, "obj_spec_matrix")
    values_proto = tensor_to_onnx(values, "obj_spec_values")
    model.graph.initializer.extend([spec_matrix_proto, values_proto])
    return model

def get_input_shape(model: onnx.ModelProto):
    input_shape = model.graph.input[0].type.tensor_type.shape.dim
    input_shape = [x.dim_value for x in input_shape]
    return input_shape

def rename_nodes(model: onnx.ModelProto):
    # Rename nodes to f"{OpType}_{num_of_op}"
    opcounts = defaultdict(int)
    for node in model.graph.node:
        node.name = f"{node.op_type}{opcounts[node.op_type]}"
        opcounts[node.op_type] += 1

    tensor_remap = {}
    # rename every node output to f"{node.name}_{output_index}"
    for node in model.graph.node:
        for i, output in enumerate(node.output):
            tensor_remap[output] = f"{node.name}/{i}"
            node.output[i] = f"{node.name}/{i}"
    
    # rename every input to nodes based on tensor_remap
    for node in model.graph.node:
        for i, input_ in enumerate(node.input):
            if input_ in tensor_remap:
                node.input[i] = tensor_remap[input_]
    
    # Rename outputs
    for i, output in enumerate(model.graph.output):
        if output.name in tensor_remap:
            model.graph.output[i].name = tensor_remap[output.name]

    return model


def main(root_dir: str, instance_idx: int, out_path: str, bound: bool = True):
    """
        Usage:
            python onnx_init.py PATH_TO_INSTANCE_ROOT INSTANCE_IDX ONNX_OUTPUT_NAME [BOUND_WITH_LIRPA_IBP]

            example:
            python onnx_init.py ../vnncomp2023_benchmarks/benchmarks/collins_rul_cnn/ 10 ./onnx/CollinsRulCNN.onnx False
    """

    # load instance.csv from root_dir
    instances     = pd.read_csv(str(Path(root_dir) / "instances.csv"), names=["onnx", "vnn", "time"])
    onnx_path     = str(Path(root_dir) / instances.iloc[instance_idx]['onnx'])
    instance_path = str(Path(root_dir) / instances.iloc[instance_idx]['vnn'])

    proto = load_onnx(onnx_path)

    lb, ub, specs, values = read_vnnlib(instance_path)

    if not torch.all(lb <= ub):
        raise ValueError("Lower bound is not less than upper bound")

    proto = rename_nodes(proto)
    proto = run_shape_inference(proto)

    input_shape = get_input_shape(proto)
    lb = lb.reshape(input_shape)
    ub = ub.reshape(input_shape)
    inp = (ub + lb) / 2

    proto = get_forward_initialization(proto, inp)
    proto = add_input_bounds(proto, lb, ub)
    proto = add_specification(proto, specs, values)

    onnx.save(proto, out_path)

if __name__ == "__main__":
    Fire(main)