#pragma once

#include <network/Layers/Linear.hpp>
#include <network/Layers/Unary.hpp>
#include <network/Layers/LayerBase.hpp>

class FusedGEMMRelu: public Layer {
public:
    FusedGEMMRelu(GEMM* gemm, Relu* relu): Layer(), gemm(gemm), relu(relu) {
        this->inputs = gemm->inputs;
        this->outputs = relu->outputs;

        // So GEMM output indexes into ReLU output
        this->gemm->Y = this->relu->Y;

        this->operator_type = _fused_gemm_relu;

        this->name = "Fused/" + this->gemm->name + "/" + this->relu->name;

        this->range_lower = 0.0;
    }

    void index_hidden_states(indices& hidden_states, indices& y_ids) override {
        for (auto i = 0; i < this->relu->Y->numel; i++) {
            hidden_states.add(this->relu->Y->strkey(i));
            y_ids.add(this->relu->Y->strkey(i));
        }
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return this->gemm->get_indices();
    }

    void add_parameters(gravity::param<>& w) const override {
        this->gemm->add_parameters(w);
    }

    void index_constraint(IndexSet& inds) override {
        this->gemm->index_constraint(inds);
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {
        gravity::func<> expr = x.in(inds["In"])*w.in(inds["B"]) + w.in(inds["C"]);

        /* Constraints */
        Constraint<> ReLU("FusedGEMMReLU");
        ReLU = x.in(inds["Out"]) - expr;
        NN.add(ReLU.in(inds["Constr"]) >= 0);

        // Using == 0 instead of <= 0 on these constraints seems better
        Constraint<> ReLU_on("FusedGEMMReLU_on");
        ReLU_on = x.in(inds["Out"]) - expr;
        NN.add_on_off(ReLU_on.in(inds["Constr"]) == 0, y.in(inds["Constr"]), true);

        Constraint<> ReLU_y_off("FusedGEMMReLU_off");
        ReLU_y_off = x.in(inds["Out"]);
        NN.add_on_off(ReLU_y_off.in(inds["Constr"]) == 0, y.in(inds["Constr"]), false);
    }

    GEMM* gemm;
    Relu* relu;
};