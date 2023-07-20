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
        gravity::param<double> l(this->lname() + "l"), u(this->lname() + "u"), f_l(this->lname() + "f_l"), f_u(this->lname() + "f_u"), df_l(this->lname() + "df_l"), df_u(this->lname() + "df_u");
        gravity::param<double> k(this->lname() + "k");
        gravity::param<double> m(this->lname() + "m"), f_m(this->lname() + "f_m"), df_m(this->lname() + "df_m");
        
        for (int i = 0; i < this->X->numel; ++i) {
            l.add_val(std::max(HMIN, this->X->lb[i]));
            u.add_val(std::min(HMAX, this->X->ub[i]));
            m.add_val((l.eval(i) + u.eval(i))/2);
            f_l.add_val(exp(l.eval(i)));
            f_u.add_val(exp(u.eval(i)));
            f_m.add_val(exp(m.eval(i)));
            df_l.add_val(exp(l.eval(i)));
            df_u.add_val(exp(u.eval(i)));
            df_m.add_val(exp(m.eval(i)));
            k.add_val((f_u.eval(i) - f_l.eval(i))/(u.eval(i)-l.eval(i)));
        }
        
        // upper bound
        Constraint<> Exp_Upper(this->lname() + "_Exp_Upper");
        Exp_Upper = x.in(inds["Out"]) - (k * (x.in(inds["In"]) - l) + f_l);
        NN.add(Exp_Upper.in(inds["Constr"]) <= 0);
        
        // lower bound at l
        Constraint<> Exp_OA1(this->lname() + "_OA1");
        Exp_OA1 = x.in(inds["Out"]) - (df_l * (x.in(inds["In"]) - l) + f_l);
        NN.add(Exp_OA1.in(inds["Constr"]) >= 0);
        
        // lower bound at u
        Constraint<> Exp_OA2(this->lname() + "_OA2");
        Exp_OA2 = x.in(inds["Out"]) - (df_u * (x.in(inds["In"]) - u) + f_u);
        NN.add(Exp_OA2.in(inds["Constr"]) >= 0);

        // lower bound at (l+u)/2
        Constraint<> Exp_OA3(this->lname() + "_OA3");
        Exp_OA3 = x.in(inds["Out"]) - (df_m * (x.in(inds["In"]) - m) + f_m);
        NN.add(Exp_OA3.in(inds["Constr"]) >= 0);
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
private:
    template <typename T>
    std::vector<T> exclude_index(const std::vector<T>& vec, int idx) {
        std::vector<T> result;
        for (int i = 0; i < vec.size(); ++i) {
            if (i != idx) {
                result.push_back(vec[i]);
            }
        }
        return result;
    }
    
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
        return {{"Out", "In", "ExpAux", "ExpSum", "SumAux", "SumProd", "OutRow"}, {}};
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
                inds["OutRow"].add_in_row(inds.row_id, this->Y->strkey(inner_ind));

                inds["In"].add_ref(this->X->strkey(inner_ind));
                inds["ExpAux"].add_ref(this->Y->strkey(inner_ind) + "_exp_aux");

                inds["ExpSum"].add_in_row(inds.row_id, this->Y->strkey(inner_ind) + "_exp_aux");
                inds["SumProd"].add_ref(this->Y->strkey(outer_ind) + "_sum_aux");
            }
            inds.row_id++;
        }
    }

    void add_constraints(gravity::Model<>& NN, IndexSet& inds, gravity::param<>& w, gravity::var<>& x, gravity::var<int>& y) override {
        int d = this->X->shape.at(this->axis);
        auto outer_shape = this->X->shape;
        outer_shape.at(this->axis) = 1;
        
        int num_outer_iterations = 1; // The number of iterations in the outer loop, initialized to 1
        for (int dim : outer_shape) {
            num_outer_iterations *= dim;
        }

        int stride = 1; // The stride length when accessing the desired axis
        for (int i = this->X->shape.size() - 1; i > this->axis; --i) {
            stride *= this->X->shape[i];
        }

        for (int i = 0; i < num_outer_iterations; ++i) {
            std::vector<double> l(d, 0.0), u(d, 0.0);
            for (int j = 0; j < d; ++j) {
                int flat_index = (i * d + j) * stride;
                l[j] = this->X->lb[flat_index];
                u[j] = this->X->ub[flat_index];
            }
            for (int j = 0; j < d; ++j) {
                std::vector<double> l_excluded = exclude_index(l, j);
                std::vector<double> u_excluded = exclude_index(u, j);

                std::vector<double> diffs_l(d-1, 0.0), diffs_u(d-1, 0.0);
                std::transform(l_excluded.begin(), l_excluded.end(), diffs_l.begin(), [&](double val){return val - u[j];});
                std::transform(u_excluded.begin(), u_excluded.end(), diffs_u.begin(), [&](double val){return val - l[j];});
                
                // Tangent points for bounding exponentials from below
                std::vector<double> diffs_t(diffs_u.size());
                std::transform(diffs_u.begin(), diffs_u.end(), diffs_l.begin(), diffs_t.begin(), [](double u_val, double l_val){
                    return std::min(std::log((std::exp(u_val) - std::exp(l_val)) / (u_val - l_val)), l_val + 1);
                });
                
                // Lower and upper bounds on denominator of softmax
                double den_l = 1 + std::inner_product(diffs_t.begin(), diffs_t.end(), diffs_l.begin(), 0.0, std::plus<>(), [](double t_val, double l_val){return std::exp(t_val) * (l_val - t_val + 1);});
                double den_u = 1 + std::accumulate(diffs_u.begin(), diffs_u.end(), 0.0, [](double sum, double u_val){return sum + std::exp(u_val);});
                
                // Tangent point for bounding reciprocal from below
                double den_t = std::max(std::sqrt(den_l * den_u), den_u / 2);
                
                // Coefficients of linear upper bound
                std::vector<double> a_lin_u(d, 0.0), a_lin_l(d, 0.0);

                for (int i = 0; i < d; ++i) {
                    if (i != j) {
                        a_lin_u[i] = -std::exp(diffs_t[i - (i > j)]) / (den_l * den_u);
                        a_lin_l[i] = -(std::exp(diffs_u[i - (i > j)]) - std::exp(diffs_l[i - (i > j)])) / (diffs_u[i - (i > j)] - diffs_l[i - (i > j)]) / std::pow(den_t, 2);
                    }
                }
                a_lin_u[j] = -std::accumulate(a_lin_u.begin(), a_lin_u.end(), 0.0);
                a_lin_l[j] = -std::accumulate(a_lin_l.begin(), a_lin_l.end(), 0.0);

                double b_lin_u = 1 / den_l + 1 / den_u - (1 + std::inner_product(diffs_t.begin(), diffs_t.end(), diffs_l.begin(), 0.0, std::plus<>(), [](double t_val, double l_val){return std::exp(t_val) * (1 - t_val);})) / (den_l * den_u);
                
                double b_lin_l = std::inner_product(diffs_u.begin(), diffs_u.end(), diffs_l.begin(), 0.0, std::plus<>(), [](double u_val, double l_val){
                    return (u_val * std::exp(l_val) - l_val * std::exp(u_val)) / (u_val - l_val);
                });
                b_lin_l = 1 / den_t * (2 - 1 / den_t * (1 + b_lin_l));
                
                //TODO: add linear bounds for output j
            }
        }

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
        
        // Sum to 1 constraint
        Constraint<> SumTo1_(this->lname() + "_Softmax_SumTo1");
        SumTo1_ = x.in(inds["OutRow"]) - 1;
        SumTo1_.print();
        NN.add(SumTo1_.in(inds["ConstrB"]) == 0);
        
        // Ouput non-negative
        Constraint<> OutNonNeg_(this->lname() + "_Softmax_OutNonNeg");
        OutNonNeg_ = x.in(inds["Out"]);
        NN.add(OutNonNeg_.in(inds["Constr"]) >= 0);
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
