import torch
from vnnlib.compat import CompatTransformer
from vnnlib.parser import parse_file


def read_vnnlib(instance):
    ast_node = parse_file(instance, strict=False)
    bounds, obj = CompatTransformer("X", "Y").transform(ast_node)[0]

    lb, ub = torch.tensor(bounds, dtype=torch.float32).chunk(2, dim=1)
    lb = lb.t()
    ub = ub.t()

    specs =  torch.vstack([torch.from_numpy(o[0]).float() for o in obj])
    values = torch.vstack([torch.from_numpy(o[1]).float() for o in obj])

    return lb, ub, specs, values