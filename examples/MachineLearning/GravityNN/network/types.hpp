#pragma once

#include <string>

typedef enum { 
    _add,
    _average_pool,
    _batchnorm,
    _clip,
    _concat,
    _conv,
    _cos,
    _div,
    _exp,
    _flatten,
    _gather,
    _gemm,
    _input,
    _leaky_relu,
    _matmul,
    _mul,
    _neg,
    _pow,
    _reduce_sum,
    _relu,
    _reshape,
    _sigmoid,
    _signum,
    _sin,
    _slice,
    _softmax,
    _split,
    _squeeze,
    _sub,
    _transpose,
    _unsupported,
} OType; /* Operator Type */
typedef enum { LOWER, UPPER} Side;

class Bound {
public:
    Bound(std::string layer_name, std::string neuron_name, double value, Side side) {
        this->layer_name = layer_name;
        this->neuron_name = neuron_name;
        this->value = value;
        this->side = side;
        this->old_value = value;
    }
    
    Bound(std::string layer_name, std::string neuron_name, double value, double old_value, Side side) {
        this->layer_name = layer_name;
        this->neuron_name = neuron_name;
        this->value = value;
        this->side = side;
        this->old_value = old_value;
    }

    std::string layer_name;
    std::string neuron_name;
    double value;
    double old_value;
    Side side;
};
