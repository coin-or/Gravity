# To Do

## Gravity Operators
- [ ] AveragePool
- [ ] ArgMax
- [ ] Cast
- [ ] ConvTranspose
- [ ] Equal
- [ ] Gather
- [ ] OneHot
- [ ] Pad
- [ ] !!! New operators in 2023 benchmarks !!!

## Tasks
- [ ] Generalize rolling horizon
- [ ] IBP in Gravity
- [ ] Smart ReLU
- [ ] Envelopes for nonlinear functions
  - [ ] Sin
  - [ ] Cos
  - [ ] Exp
- [ ] GRBTune eventually

## Stretch
- [ ] Zonotopes
- [ ] Support ViT (involves supporting higher rank tensors in various operators)


# Benchmarks
- [X] acasxu
- [X] tllverifybench
- [X] collins_rul_cnn
- [X] dist_shift | BROKEN BECAUSE OF NONLINEAR
  - Concat, Gemm, Relu, Reshape, Sigmoid
- [X] ml4acopf | BROKEN BECAUSE OF NONLINEAR
  - Transpose, Sub, Pow, Gemm, Neg, Slice, Relu, MatMul, Cos, Sin, Add, Mul, Concat, Gather, Sigmoid
- [X] nn4sys
  - Sigmoid, Sub, Pow, Gemm, Concat, Slice, Flatten, Relu, Mul, Div, MatMul, Reshape, ReduceSum, Conv, Gather, Add, Split
- [X] vit
- [X] yolo
- [X] traffic_signs_recognition
  - Binary convs
- [X] vggnet16
  - [ ] Simply too big
- [ ] cctsdb_yolo
  - Argmax + A bunch of other unsupported stuff
- [ ] cgan
  - ConvTranspose, Resize, Tanh
- [ ] collins_yolo_robustness
  - Resize + MaxPool 
- [ ] metaroom
  - BROKEN BECAUSE OF PROJECTION OP