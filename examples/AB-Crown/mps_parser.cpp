#include <spdlog/spdlog.h>
#include <string>
#include <gurobi_c++.h>
#include <gravity/model.h>
#include <gravity/func.h>
#include <map>

using namespace gravity;

template<typename T>
std::pair<T, T> check_inf_bounds(T lb, T ub) {
	// If ub == 1e100, then it is infinity. If lb == -1e100, then it is -infinity
	ub = (ub ==  1e100) ? numeric_limits<T>::max() : ub;
	lb = (lb == -1e100) ? numeric_limits<T>::lowest() : lb;
	return std::make_pair(lb, ub);
}

gravity::func<> from_lin_expr(GRBLinExpr& row, Model<>& gravity_model) {
	gravity::func<> expr;
	for (int j = 0; j < row.size(); j++) {
		double coeff = row.getCoeff(j);
		if (coeff == 0.0) {
			continue;
		}
		
		GRBVar var = row.getVar(j);

		char vtype = var.get(GRB_CharAttr_VType);
		std::string var_name = var.get(GRB_StringAttr_VarName);

		if (vtype == GRB_CONTINUOUS) {
			expr += coeff * gravity_model.get_var<double>(var_name);
		} else if (vtype == GRB_INTEGER) {
			expr += coeff * gravity_model.get_var<int>(var_name);
		} else if (vtype == GRB_BINARY) {
			expr += coeff * gravity_model.get_var<bool>(var_name);
		} else {
			spdlog::error("Unknown variable type: {}", vtype);
			exit(1);
		}
	}

	return expr;
}

void add_vars(GRBModel& model, Model<>& gravity_model) {
	for (int i = 0; i < model.get(GRB_IntAttr_NumVars); i++) {
		GRBVar gurobi_var = model.getVar(i);
		double lb = gurobi_var.get(GRB_DoubleAttr_LB);
		double ub = gurobi_var.get(GRB_DoubleAttr_UB);
		const char type = gurobi_var.get(GRB_CharAttr_VType);
		const std::string name = gurobi_var.get(GRB_StringAttr_VarName);

		if (type == GRB_CONTINUOUS) {
			auto bounds = check_inf_bounds<double>(static_cast<double>(lb), static_cast<double>(ub));
			var<double> gravity_var(name, bounds.first, bounds.second);
			gravity_model.add(gravity_var);
		} else if (type == GRB_INTEGER) {
			auto bounds = check_inf_bounds<int>(static_cast<int>(lb), static_cast<int>(ub));
			var<int> gravity_var(name, bounds.first, bounds.second);
			gravity_model.add(gravity_var);
		} else if (type == GRB_BINARY) {
			auto bounds = check_inf_bounds<bool>(static_cast<bool>(lb), static_cast<bool>(ub));
			var<bool> gravity_var(name, bounds.first, bounds.second);
			gravity_model.add(gravity_var);
		} else {
			spdlog::error("Unknown variable type: {}", type);
			exit(1);
		}
	}
}

void add_objective(GRBModel& gurobi_model, Model<>& gravity_model) {
	GRBQuadExpr quad_obj = gurobi_model.getObjective();
	GRBLinExpr lin_obj = quad_obj.getLinExpr();

	gravity::func<> gravityObj = from_lin_expr(lin_obj, gravity_model);

	// Get model objsense
	int objsense = gurobi_model.get(GRB_IntAttr_ModelSense);
	if (objsense == GRB_MAXIMIZE) {
		gravity_model.max(gravityObj);
	} else if (objsense == GRB_MINIMIZE) {
		gravity_model.min(gravityObj);
	} else {
		spdlog::error("Unknown objsense: {}", objsense);
		exit(1);
	}

}

void add_constraints(GRBModel& gurobi_model, Model<>& gravity_model) {
	int n_constraints = gurobi_model.get(GRB_IntAttr_NumConstrs);
	int n_vars = gurobi_model.get(GRB_IntAttr_NumVars);

	for (int i = 0; i < n_constraints; i++) {
		GRBConstr constr = gurobi_model.getConstr(i);
		std::string constr_name = constr.get(GRB_StringAttr_ConstrName);
		char sense = constr.get(GRB_CharAttr_Sense);
		double rhs = constr.get(GRB_DoubleAttr_RHS);

		GRBLinExpr row = gurobi_model.getRow(constr);
		Constraint<double> expr(constr_name);
		expr = from_lin_expr(row, gravity_model);

		if (sense == GRB_LESS_EQUAL) {
			gravity_model.add(expr <= rhs);
		} else if (sense == GRB_GREATER_EQUAL) {
			gravity_model.add(expr >= rhs);
		} else if (sense == GRB_EQUAL) {
			gravity_model.add(expr == rhs);
		} else {
			spdlog::error("Unknown constraint sense: {}", sense);
			exit(1);
		}
	}
}

int main(int argc, char** argv) {
	GRBEnv env = GRBEnv();
	GRBModel gurobi_model = GRBModel(env, "/mnt/trail_test/haydnj/vnn/Gravity/simple.mps");
	Model<> gravity_model("GravityModel");

	int n_vars = gurobi_model.get(GRB_IntAttr_NumVars);
	int n_constraints = gurobi_model.get(GRB_IntAttr_NumConstrs);
	int n_int_vars = gurobi_model.get(GRB_IntAttr_NumIntVars);
	int n_bin_vars = gurobi_model.get(GRB_IntAttr_NumBinVars);
	spdlog::info("NumVars: {}, NumConstraints: {}", n_vars, n_constraints);
	spdlog::info("NumIntVars: {}, NumBinVars: {}", n_int_vars, n_bin_vars);

	add_vars(gurobi_model, gravity_model);

	spdlog::info("NVars: {}", gravity_model._vars.size());

	// Make sure all variables were added
	if (gravity_model._vars.size() != n_vars) {
		spdlog::error("Not all variables were added. Only Continuous, Integer, and Binary variables are supported.");
		exit(1);
	}

	spdlog::info("Done adding variables.");
	add_constraints(gurobi_model, gravity_model);

	// Build objective
	add_objective(gurobi_model, gravity_model);

	gravity_model.print();
	return 0;
}