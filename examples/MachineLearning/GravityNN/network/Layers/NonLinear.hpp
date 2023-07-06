#pragma once

#include <network/Layers/LayerBase.hpp>

using namespace gravity;

class Cos : public Layer {
public:
    Cos(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _cos;
        this->X = &tensors[node.input(0)];
        this->Y = &tensors[node.output(0)];

        this->range_lower = -1.0;
        this->range_upper =  1.0;
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"Out", "In"}, {}, {"z"}};
    }
    
    void index_hidden_states(indices& hidden_states, indices& y_ids) override {
        for (auto i = 0; i < this->Y->numel; i++) {
            auto key = this->Y->strkey(i);
            hidden_states.add(key);
            y_ids.add(key);
        }
    }

    void index_constraint(IndexSet& inds) override {
        for(auto j = 0; j < this->X->numel;j++){
            inds["Constr"].add(this->Y->strkey(j));
            inds["In"].add_ref(this->X->strkey(j));
            inds["Out"].add_ref(this->Y->strkey(j));
            inds["z"].add_ref(this->Y->strkey(j));
        }
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {
        gravity::param<double> l(this->lname() + "l"), u(this->lname() + "u"), f_l(this->lname() + "f_l"), f_u(this->lname() + "f_u"), df_l(this->lname() + "df_l"), df_u(this->lname() + "df_u");
        gravity::param<double> ifln(this->lname() + "ifln"), f_ifln(this->lname() + "f_ifln"), df_ifln(this->lname() + "df_ifln");
        gravity::param<double> m_lhs(this->lname() + "m_lhs"), m_rhs(this->lname() + "m_rhs");
        gravity::param<double> z_l(this->lname() + "z_l"), z_u(this->lname() + "z_u");

        for (int i = 0; i < this->X->numel; ++i) {
            l.add_val(std::max(-0.5*M_PI, this->X->lb[i]));
            u.add_val(std::min(1.5*M_PI, this->X->ub[i]));
            f_l.add_val(cos(l.eval(i)));
            f_u.add_val(cos(u.eval(i)));
            df_l.add_val(-sin(l.eval(i)));
            df_u.add_val(-sin(u.eval(i)));
            ifln.add_val(0.5*M_PI);
            f_ifln.add_val(0.0);
            df_ifln.add_val(-1.0);
            z_l.add_val(0.0);
            z_u.add_val(1.0);
        }
        for (int i = 0; i < this->X->numel; ++i) {
            if (u.eval(i) <= ifln.eval(i)) {
                ifln.set_val(i, u.eval(i));
                f_ifln.set_val(i, f_u.eval(i));
                df_ifln.set_val(i, df_u.eval(i));
                z_u.set_val(i, 0.0);
            }
            if (l.eval(i) >= ifln.eval(i)) {
                ifln.set_val(i, l.eval(i));
                f_ifln.set_val(i, f_l.eval(i));
                df_ifln.set_val(i, df_l.eval(i));
                z_l.set_val(i, 1.0);
            }
            if (z_l.eval(i) > z_u.eval(i)) {
                z_l.set_val(i, z_u.eval(i));
            }
            
           double d1 = (l.eval(i) - ifln.eval(i));
           double d2 = (u.eval(i) - ifln.eval(i));
           if (d1 != 0) {
               m_lhs.add_val((f_l.eval(i) - f_ifln.eval(i)) / d1);
           } else {
               m_lhs.add_val(0.0);
           }
           if (d2 != 0) {
               m_rhs.add_val((f_u.eval(i) - f_ifln.eval(i)) / d2);
           } else {
               m_rhs.add_val(0.0);
           }
        }
        
        // LHS of the inflection point
        Constraint<> Cos_LHS_Upper_1(this->lname() + "_Cos_LHS_Upper_1");
        Cos_LHS_Upper_1 = x.in(inds["Out"]) - (df_l * (x.in(inds["In"]) - l) + f_l);
        NN.add_on_off(Cos_LHS_Upper_1.in(inds["Constr"]) <= 0, y.in(inds["Constr"]), false);
        
        Constraint<> Cos_LHS_Upper_2(this->lname() + "_Cos_LHS_Upper_2");
        Cos_LHS_Upper_2 = x.in(inds["Out"]) - (df_ifln * (x.in(inds["In"]) - ifln) + f_ifln);
        NN.add_on_off(Cos_LHS_Upper_2.in(inds["Constr"]) <= 0, y.in(inds["Constr"]), false);
        
        Constraint<> Cos_LHS_Lower(this->lname() + "_Cos_LHS_Lower");
        Cos_LHS_Lower = x.in(inds["Out"]) - (m_lhs * (x.in(inds["In"]) - l) + f_l);
        NN.add_on_off(Cos_LHS_Lower.in(inds["Constr"]) >= 0, y.in(inds["Constr"]), false);
        
        // RHS of the inflection point
        Constraint<> Cos_RHS_Lower_1(this->lname() + "_Cos_RHS_Lower_1");
        Cos_RHS_Lower_1 = x.in(inds["Out"]) - (df_u * (x.in(inds["In"]) - u) + f_u);
        NN.add_on_off(Cos_RHS_Lower_1.in(inds["Constr"]) >= 0, y.in(inds["Constr"]), true);
        
        Constraint<> Cos_RHS_Lower_2(this->lname() + "_Cos_RHS_Lower_2");
        Cos_RHS_Lower_2 = x.in(inds["Out"]) - (df_ifln * (x.in(inds["In"]) - ifln) + f_ifln);
        NN.add_on_off(Cos_RHS_Lower_2.in(inds["Constr"]) >= 0, y.in(inds["Constr"]), true);
        
        Constraint<> Cos_RHS_Upper(this->lname() + "_Cos_RHS_Upper");
        Cos_RHS_Upper = x.in(inds["Out"]) - (m_rhs * (x.in(inds["In"]) - ifln) + f_ifln);
        NN.add_on_off(Cos_RHS_Upper.in(inds["Constr"]) <= 0, y.in(inds["Constr"]), true);
        
        // enforce bounds on z
        Constraint<> Cos_Z_Lower(this->lname() + "_Cos_Z_Lower");
        Cos_Z_Lower = y.in(inds["z"]) - z_l;
        NN.add(Cos_Z_Lower.in(inds["Constr"]) >= 0);
        
        Constraint<> Cos_Z_Upper(this->lname() + "_Cos_Z_Upper");
        Cos_Z_Upper = y.in(inds["z"]) - z_u;
        NN.add(Cos_Z_Upper.in(inds["Constr"]) <= 0);
    }

    Tensor *X, *Y; // Input and output
};

class Sin : public Layer {
public:
    Sin(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _sin;
        this->X = &tensors[node.input(0)];
        this->Y = &tensors[node.output(0)];

        this->range_lower = -1.0;
        this->range_upper =  1.0;
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"Out", "In"}, {}, {"z"}};
    }
    
    void index_hidden_states(indices& hidden_states, indices& y_ids) override {
        for (auto i = 0; i < this->Y->numel; i++) {
            auto key = this->Y->strkey(i);
            hidden_states.add(key);
            y_ids.add(key);
        }
    }

    void index_constraint(IndexSet& inds) override {
        for(auto j = 0; j < this->X->numel;j++){
            inds["Constr"].add(this->Y->strkey(j));
            inds["In"].add_ref(this->X->strkey(j));
            inds["Out"].add_ref(this->Y->strkey(j));
            inds["z"].add_ref(this->Y->strkey(j));
        }
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {
        gravity::param<double> l(this->lname() + "l"), u(this->lname() + "u"), f_l(this->lname() + "f_l"), f_u(this->lname() + "f_u"), df_l(this->lname() + "df_l"), df_u(this->lname() + "df_u");
        gravity::param<double> ifln(this->lname() + "ifln"), f_ifln(this->lname() + "f_ifln"), df_ifln(this->lname() + "df_ifln");
        gravity::param<double> m_lhs(this->lname() + "m_lhs"), m_rhs(this->lname() + "m_rhs");
        gravity::param<double> z_l(this->lname() + "z_l"), z_u(this->lname() + "z_u");

        for (int i = 0; i < this->X->numel; ++i) {
            l.add_val(std::max(-M_PI, this->X->lb[i]));
            u.add_val(std::min(M_PI, this->X->ub[i]));
            f_l.add_val(sin(l.eval(i)));
            f_u.add_val(sin(u.eval(i)));
            df_l.add_val(cos(l.eval(i)));
            df_u.add_val(cos(u.eval(i)));
            ifln.add_val(0.0);
            f_ifln.add_val(0.0);
            df_ifln.add_val(1.0);
            z_l.add_val(0.0);
            z_u.add_val(1.0);
        }
        for (int i = 0; i < this->X->numel; ++i) {
            if (u.eval(i) <= ifln.eval(i)) {
                ifln.set_val(i, u.eval(i));
                f_ifln.set_val(i, f_u.eval(i));
                df_ifln.set_val(i, df_u.eval(i));
                z_u.set_val(i, 0.0);
            }
            if (l.eval(i) >= ifln.eval(i)) {
                ifln.set_val(i, l.eval(i));
                f_ifln.set_val(i, f_l.eval(i));
                df_ifln.set_val(i, df_l.eval(i));
                z_l.set_val(i, 1.0);
            }
            if (z_l.eval(i) > z_u.eval(i)) {
                z_l.set_val(i, z_u.eval(i));
            }
            double d1 = (l.eval(i) - ifln.eval(i));
            double d2 = (u.eval(i) - ifln.eval(i));
            if (d1 != 0) {
                m_lhs.add_val((f_l.eval(i) - f_ifln.eval(i)) / d1);
            } else {
                m_lhs.add_val(0.0);
            }
            if (d2 != 0) {
                m_rhs.add_val((f_u.eval(i) - f_ifln.eval(i)) / d2);
            } else {
                m_rhs.add_val(0.0);
            }
        }
        // LHS of the inflection point
        Constraint<> Sin_LHS_Lower_1(this->lname() + "_Sin_LHS_Lower_1");
        Sin_LHS_Lower_1 = x.in(inds["Out"]) - (df_l * (x.in(inds["In"]) - l) + f_l);
        NN.add_on_off(Sin_LHS_Lower_1.in(inds["Constr"]) >= 0, y.in(inds["Constr"]), false);
        
        Constraint<> Sin_LHS_Lower_2(this->lname() + "_Sin_LHS_Lower_2");
        Sin_LHS_Lower_2 = x.in(inds["Out"]) - (df_ifln * (x.in(inds["In"]) - ifln) + f_ifln);
        NN.add_on_off(Sin_LHS_Lower_2.in(inds["Constr"]) >= 0, y.in(inds["Constr"]), false);
        
        Constraint<> Sin_LHS_Upper(this->lname() + "_Sin_LHS_Upper");
        Sin_LHS_Upper = x.in(inds["Out"]) - (m_lhs * (x.in(inds["In"]) - l) + f_l);
        NN.add_on_off(Sin_LHS_Upper.in(inds["Constr"]) <= 0, y.in(inds["Constr"]), false);
        
        // RHS of the inflection point
        Constraint<> Sin_RHS_Upper_1(this->lname() + "_Sin_RHS_Upper_1");
        Sin_RHS_Upper_1 = x.in(inds["Out"]) - (df_u * (x.in(inds["In"]) - u) + f_u);
        NN.add_on_off(Sin_RHS_Upper_1.in(inds["Constr"]) <= 0, y.in(inds["Constr"]), true);
        
        Constraint<> Sin_RHS_Upper_2(this->lname() + "_Sin_RHS_Upper_2");
        Sin_RHS_Upper_2 = x.in(inds["Out"]) - (df_ifln * (x.in(inds["In"]) - ifln) + f_ifln);
        NN.add_on_off(Sin_RHS_Upper_2.in(inds["Constr"]) <= 0, y.in(inds["Constr"]), true);
        
        Constraint<> Sin_RHS_Lower(this->lname() + "_Sin_RHS_Lower");
        Sin_RHS_Lower = x.in(inds["Out"]) - (m_rhs * (x.in(inds["In"]) - ifln) + f_ifln);
        NN.add_on_off(Sin_RHS_Lower.in(inds["Constr"]) >= 0, y.in(inds["Constr"]), true);
        
        // enforce bounds on z
        Constraint<> Sin_Z_Lower(this->lname() + "_Sin_Z_Lower");
        Sin_Z_Lower = y.in(inds["z"]) - z_l;
        NN.add(Sin_Z_Lower.in(inds["Constr"]) >= 0);
        
        Constraint<> Sin_Z_Upper(this->lname() + "_Sin_Z_Upper");
        Sin_Z_Upper = y.in(inds["z"]) - z_u;
        NN.add(Sin_Z_Upper.in(inds["Constr"]) <= 0);
    }

    Tensor *X, *Y; // Input and output
};

class Pow : public Layer {
public:
    Pow(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _pow;
        this->X = &tensors[node.input(0)];
        this->exp = &tensors[node.input(1)];
        this->Y = &tensors[node.output(0)];

        if (this->exp->numel != 1) {
            throw std::runtime_error("Pow: exponent must be a scalar.");
        }

        if (this->exp->operator()(0) != 2.0) {
            throw std::runtime_error("Pow: exponent must be 2.");
        }

        this->range_lower = 0.0;
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"Out", "In"}, {}};
    }

    void index_constraint(IndexSet& inds) override {
        for(auto j = 0; j < this->X->numel;j++){
            inds["Constr"].add(this->Y->strkey(j));
            inds["In"].add_ref(this->X->strkey(j));
            inds["Out"].add_ref(this->Y->strkey(j));
        }
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {
        Constraint<> Pow_(this->lname() + "_Pow");
        Pow_ = x.in(inds["Out"]) - pow(x.in(inds["In"]), 2.0);
        NN.add(Pow_.in(inds["Constr"]) == 0);
    }


    Tensor *X, *exp, *Y; // Input and output
};

class Exp: public Layer {
public:
    Exp(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _exp;
        this->X = &tensors[node.input(0)];
        this->Y = &tensors[node.output(0)];

        this->range_lower = 0.0;
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"Out", "In"}, {}};
    }

    void index_constraint(IndexSet& inds) override {
        for(auto j = 0; j < this->X->numel;j++){
            inds["Constr"].add(this->Y->strkey(j));
            inds["In"].add_ref(this->X->strkey(j));
            inds["Out"].add_ref(this->Y->strkey(j));
        }
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {
        Constraint<> Exp_(this->lname() + "_Exp");
        Exp_ = x.in(inds["Out"]) - gravity::exp(x.in(inds["In"]));
        NN.add(Exp_.in(inds["Constr"]) == 0);
    }

    Tensor *X, *Y; // Input and output
};

class Sigmoid: public Layer {
public:
    Sigmoid(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _sigmoid;
        this->X = &tensors[node.input(0)];
        this->Y = &tensors[node.output(0)];

        this->range_lower = 0.0;
        this->range_upper = 1.0;
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"Out", "In"}, {}, {"z"}};
    }

    void index_hidden_states(indices& hidden_states, indices& y_ids) override {
        for (auto i = 0; i < this->Y->numel; i++) {
            auto key = this->Y->strkey(i);
            hidden_states.add(key);
            y_ids.add(key);
        }
    }
    
    void index_constraint(IndexSet& inds) override {
        for(auto j = 0; j < this->X->numel;j++){
            inds["Constr"].add(this->Y->strkey(j));
            inds["In"].add_ref(this->X->strkey(j));
            inds["Out"].add_ref(this->Y->strkey(j));
            inds["z"].add_ref(this->Y->strkey(j));
        }
    }
    
    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {
        gravity::param<double> l(this->lname() + "l"), u(this->lname() + "u"), f_l(this->lname() + "f_l"), f_u(this->lname() + "f_u"), df_l(this->lname() + "df_l"), df_u(this->lname() + "df_u");
        gravity::param<double> ifln(this->lname() + "ifln"), f_ifln(this->lname() + "f_ifln"), df_ifln(this->lname() + "df_ifln");
        gravity::param<double> m_lhs(this->lname() + "m_lhs"), m_rhs(this->lname() + "m_rhs");
        gravity::param<double> z_l(this->lname() + "z_l"), z_u(this->lname() + "z_u");

        for (int i = 0; i < this->X->numel; ++i) {
            l.add_val(std::max(HMIN, this->X->lb[i]));
            u.add_val(std::min(HMAX, this->X->ub[i]));
            f_l.add_val(1.0 / (1.0 + exp(-l.eval(i))));
            f_u.add_val(1.0 / (1.0 + exp(-u.eval(i))));
            df_l.add_val(f_l.eval(i) * (1 - f_l.eval(i)));
            df_u.add_val(f_u.eval(i) * (1 - f_u.eval(i)));
            ifln.add_val(0.0);
            f_ifln.add_val(0.5);
            df_ifln.add_val(0.25);
            z_l.add_val(0.0);
            z_u.add_val(1.0);
        }
        for (int i = 0; i < this->X->numel; ++i) {
            if (u.eval(i) <= ifln.eval(i)) {
                ifln.set_val(i, u.eval(i));
                f_ifln.set_val(i, f_u.eval(i));
                df_ifln.set_val(i, df_u.eval(i));
                z_u.set_val(i, 0.0);
            }
            if (l.eval(i) >= ifln.eval(i)) {
                ifln.set_val(i, l.eval(i));
                f_ifln.set_val(i, f_l.eval(i));
                df_ifln.set_val(i, df_l.eval(i));
                z_l.set_val(i, 1.0);
            }
            if (z_l.eval(i) > z_u.eval(i)) {
                z_l.set_val(i, z_u.eval(i));
            }
            double d1 = (l.eval(i) - ifln.eval(i));
            double d2 = (u.eval(i) - ifln.eval(i));
            if (d1 != 0) {
                m_lhs.add_val((f_l.eval(i) - f_ifln.eval(i)) / d1);
            } else {
                m_lhs.add_val(0.0);
            }
            if (d2 != 0) {
                m_rhs.add_val((f_u.eval(i) - f_ifln.eval(i)) / d2);
            } else {
                m_rhs.add_val(0.0);
            }
        }
        // Negative input
        Constraint<> Sigmoid_LHS_Lower_1(this->lname() + "_Sigmoid_LHS_Lower_1");
        Sigmoid_LHS_Lower_1 = x.in(inds["Out"]) - (df_ifln * (x.in(inds["In"]) - ifln) + f_ifln); // y >= 0.25*x + 0.5
        NN.add_on_off(Sigmoid_LHS_Lower_1.in(inds["Constr"]) >= 0, y.in(inds["z"]), false);
        
        Constraint<> Sigmoid_LHS_Lower_2(this->lname() + "_Sigmoid_LHS_Lower_2");
        Sigmoid_LHS_Lower_2 = x.in(inds["Out"]) - (df_l * (x.in(inds["In"]) - l) + f_l); // y >= df_l*x + f_l - df_l * l
        NN.add_on_off(Sigmoid_LHS_Lower_2.in(inds["Constr"]) >= 0, y.in(inds["z"]), false);
        
        Constraint<> Sigmoid_LHS_Upper(this->lname() + "_Sigmoid_LHS_Upper");
        Sigmoid_LHS_Upper = x.in(inds["Out"]) - (m_lhs * (x.in(inds["In"]) - ifln) + f_ifln); // y <= (f_l - 0.5)/l * x + 0.5
        NN.add_on_off(Sigmoid_LHS_Upper.in(inds["Constr"]) <= 0, y.in(inds["z"]), false);
        
        // Positive input
        Constraint<> Sigmoid_RHS_Upper_1(this->lname() + "_Sigmoid_RHS_Upper_1");
        Sigmoid_RHS_Upper_1 = x.in(inds["Out"]) - (df_ifln * (x.in(inds["In"]) - ifln) + f_ifln); // y <= 0.25*x + 0.5
        NN.add_on_off(Sigmoid_RHS_Upper_1.in(inds["Constr"]) <= 0, y.in(inds["z"]), true);
        
        Constraint<> Sigmoid_RHS_Upper_2(this->lname() + "_Sigmoid_RHS_Upper_2");
        Sigmoid_RHS_Upper_2 = x.in(inds["Out"]) - (df_u * (x.in(inds["In"]) - u) + f_u); // y <= df_u*x + f_u - df_u * u
        NN.add_on_off(Sigmoid_RHS_Upper_2.in(inds["Constr"]) <= 0, y.in(inds["z"]), true);
        
        Constraint<> Sigmoid_RHS_Lower(this->lname() + "_Sigmoid_RHS_Lower");
        Sigmoid_RHS_Lower = x.in(inds["Out"]) - (m_rhs * (x.in(inds["In"]) - ifln) + f_ifln); // y >= (f_u - 0.5) / u * x + 0.5
        NN.add_on_off(Sigmoid_RHS_Lower.in(inds["Constr"]) >= 0, y.in(inds["z"]), true);
        
        // enforce bounds on z
        Constraint<> Sigmoid_Z_Lower(this->lname() + "_Sigmoid_Z_Lower");
        Sigmoid_Z_Lower = y.in(inds["z"]) - z_l;
        NN.add(Sigmoid_Z_Lower.in(inds["Constr"]) >= 0);
        
        Constraint<> Sigmoid_Z_Upper(this->lname() + "_Sigmoid_Z_Upper");
        Sigmoid_Z_Upper = y.in(inds["z"]) - z_u;
        NN.add(Sigmoid_Z_Upper.in(inds["Constr"]) <= 0);
    }

    Tensor *X, *Y; // Input and output
};

class Softmax: public Layer {
public:
    Softmax(const onnx::NodeProto& node, Tensors& tensors): Layer(node, tensors) {
        operator_type = _softmax;
        this->X = &tensors[node.input(0)];
        this->Y = &tensors[node.output(0)];

        if (const auto* axis = find_attribute("axis", node)) {
            int64_t i = axis->i();
            if (i < 0) {
                i += this->X->ndims;
            }
            this->axis = i;
        } else {
            throw std::runtime_error("Softmax: axis attribute not found");
        }

        this->range_lower = 0.0;
        this->range_upper = 1.0;
    }

    std::vector<std::vector<std::string>> get_indices() const override {
        return {{"Out", "In", "ExpAux", "ExpSum", "SumAux", "SumProd"}, {}};
    }

    void index_hidden_states(indices& hidden_states, indices& y_ids) override {
        auto outer_shape = this->X->shape;
        outer_shape.at(this->axis) = 1;

        for (auto outer_ind: ShapeIter(outer_shape)) {
            hidden_states.add(this->Y->strkey(outer_ind) + "_sum_aux");
            for (size_t axind = 0; axind < this->X->shape.at(this->axis); axind++) {
                auto inner_ind = outer_ind;
                inner_ind.at(this->axis) = axind;

                hidden_states.add(this->Y->strkey(inner_ind));
                hidden_states.add(this->Y->strkey(inner_ind) + "_exp_aux");
            }
        }
    }

    void index_constraint(IndexSet& inds) override {
        auto outer_shape = this->X->shape;
        outer_shape.at(this->axis) = 1;

        for (auto outer_ind: ShapeIter(outer_shape)) {
            inds["ConstrB"].add(this->Y->strkey(outer_ind) + "_sum");
            inds["SumAux"].add_ref(this->Y->strkey(outer_ind) + "_sum_aux");

            for (size_t axind = 0; axind < this->X->shape.at(this->axis); axind++) {
                auto inner_ind = outer_ind;
                inner_ind.at(this->axis) = axind;

                inds["Constr"].add(this->Y->strkey(inner_ind));
                inds["Out"].add_ref(this->Y->strkey(inner_ind));

                inds["In"].add_ref(this->X->strkey(inner_ind));
                inds["ExpAux"].add_ref(this->Y->strkey(inner_ind) + "_exp_aux");

                inds["ExpSum"].add_in_row(inds.row_id, this->Y->strkey(inner_ind) + "_exp_aux");
                inds["SumProd"].add_ref(this->Y->strkey(outer_ind) + "_sum_aux");
            }
            inds.row_id++;
        }
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {
        // Exp constraint
        Constraint<> Exp_(this->lname() + "_Softmax_Exp");
        Exp_ = x.in(inds["ExpAux"]) - gravity::exp(x.in(inds["In"]));
        NN.add(Exp_.in(inds["Constr"]) == 0);

        // Sum constraint
        Constraint<> Sum_(this->lname() + "_Softmax_Sum");
        Sum_ = x.in(inds["SumAux"]) - x.in(inds["ExpSum"]);
        Sum_.print();
        NN.add(Sum_.in(inds["ConstrB"]) == 0);

        // Out constraint
        Constraint<> Out_(this->lname() + "_Softmax_Out");
        Out_ = x.in(inds["Out"])*x.in(inds["SumProd"]) - x.in(inds["ExpAux"]);
        NN.add(Out_.in(inds["Constr"]) == 0);
    }

    void set_bounds(gravity::param<>& x_lb, gravity::param<>& x_ub) override {
        auto outer_shape = this->X->shape;
        outer_shape.at(this->axis) = 1;

        for (auto outer_ind: ShapeIter(outer_shape)) {
            x_lb.set_val(this->Y->strkey(outer_ind) + "_sum_aux", 0.0);
            x_ub.set_val(this->Y->strkey(outer_ind) + "_sum_aux", HMAX);
            for (size_t axind = 0; axind < this->X->shape.at(this->axis); axind++) {
                auto inner_ind = outer_ind;
                inner_ind.at(this->axis) = axind;

                x_lb.set_val(this->Y->strkey(inner_ind), 0.0);
                x_ub.set_val(this->Y->strkey(inner_ind), 1.0);

                x_lb.set_val(this->Y->strkey(inner_ind) + "_exp_aux", 0.0);
                x_ub.set_val(this->Y->strkey(inner_ind) + "_exp_aux", HMAX);
            }
        }
    }

    Tensor *X, *Y; // Input and output
    size_t axis;
};
