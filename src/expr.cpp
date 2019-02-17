////
////  expr.cpp
////  Gravity
////
////  Created by Hijazi, Hassan on 24/10/16.
////
////
#include <cmath>
#include <gravity/expr.h>
#include <gravity/func.h>

using namespace std;
string operator_str(gravity::OperatorType ot){
    switch (ot) {
        case gravity::log_:
            return "log";
        case gravity::exp_:
            return "exp";
        case gravity::cos_:
            return "cos";
        case gravity::sin_:
            return "sin";
        case gravity::tan_:
            return "tan";
        case gravity::sqrt_:
            return "sqrt";
        case gravity::relu_:
            return "ReLU";
        case gravity::unit_step_:
            return "UnitStep";
        default:
            break;
    }
    throw invalid_argument("Unsupported unitary operator");
}
