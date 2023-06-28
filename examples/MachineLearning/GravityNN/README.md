# Using onnx_init

First set up a conda environment with python version 3.10
```bash
conda create -n vnn2023 python=3.10
conda activate vnn2023
```

Then install requirements:
```bash
pip install -r requirements.txt
```

Now we can create onnx models. Here is an example:
```bash
python onnx_init.py ~/Documents/vnncomp2023_benchmarks/benchmarks/tllverifybench/ 0 ../../../data_sets/VNN/tll_new.onnx
```

This will initialize an ONNX model with instance 0 in the TLL benchmark and write it out to the file `data_sets/VNN/tll_new.onnx`