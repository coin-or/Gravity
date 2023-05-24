#pragma once

#include <string>

typedef enum { _gemm, _relu, _conv, _input, _noop, _add, _sub, _cos, _sin, _neg, _pow, _mul, _div, _clip} OType; /* Operator Type */
typedef enum { LOWER, UPPER} Side;

class Bound {
public:
    Bound(std::string layer_name, std::string neuron_name, float value, Side side) {
        this->layer_name = layer_name;
        this->neuron_name = neuron_name;
        this->value = value;
        this->side = side;
    }

    std::string layer_name;
    std::string neuron_name;
    float value;
    Side side;
};
