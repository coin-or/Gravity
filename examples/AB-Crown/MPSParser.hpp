#include <string>
#include <gurobi_c++.h>
#include <gravity/model.h>
#include <gravity/func.h>
#include <map>

using namespace gravity;

// GRB_DoubleAttr_MaxBound?
std::pair<bool, bool> build_bool_bounds(double lb, double ub) {
	return std::make_pair(static_cast<bool>(lb), static_cast<bool>(ub));
}

std::pair<int, int> build_int_bounds(double lb, double ub) {
	int n_ub = (ub ==  GRB_MAXINT) ? std::numeric_limits<int>::max() : ub;
	int n_lb = (lb == -GRB_MAXINT) ? std::numeric_limits<int>::lowest() : lb;

	return std::make_pair(n_ub, n_lb);
}

std::pair<double, double> build_double_bounds(double lb, double ub) {
	ub = (ub ==  GRB_INFINITY) ? std::numeric_limits<double>::max() : ub;
	lb = (lb == -GRB_INFINITY) ? std::numeric_limits<double>::lowest() : lb;

	return std::make_pair(static_cast<double>(lb), static_cast<double>(ub));
}

gravity::func<> func_from_lin_expr(GRBLinExpr& row, Model<>& gravity_model) {
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
			expr += coeff * gravity_model.get_var<int>(var_name);
		} else {
			Warning("Unknown variable type: " << vtype << std::endl);
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
			auto bounds = build_double_bounds(lb, ub);
			var<double> gravity_var(name, lb, ub);
			gravity_model.add(gravity_var.in(R(1)));
		} else if (type == GRB_INTEGER) {
			auto bounds = build_int_bounds(lb, ub);
			var<int> gravity_var(name, (int)lb, (int)ub);
			gravity_model.add(gravity_var.in(R(1)));
		} else if (type == GRB_BINARY) {
			auto bounds = build_bool_bounds(lb, ub);
			var<int> gravity_var(name, (int)lb, (int)ub);
			gravity_model.add(gravity_var.in(R(1)));
		} else {
			Warning("Unknown variable type: " << type << std::endl);
			exit(1);
		}
	}
}

void add_objective(GRBModel& gurobi_model, Model<>& gravity_model) {
	GRBQuadExpr quad_obj = gurobi_model.getObjective();
	GRBLinExpr lin_obj = quad_obj.getLinExpr();

	gravity::func<> gravityObj = func_from_lin_expr(lin_obj, gravity_model);

	// Get model objsense
	int objsense = gurobi_model.get(GRB_IntAttr_ModelSense);
	if (objsense == GRB_MAXIMIZE) {
		gravity_model.max(gravityObj);
	} else if (objsense == GRB_MINIMIZE) {
		gravity_model.min(gravityObj);
	} else {
		Warning("Unknown objective sense: " << objsense << std::endl);
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
		expr = func_from_lin_expr(row, gravity_model);

		if (sense == GRB_LESS_EQUAL) {
			gravity_model.add(expr <= rhs);
		} else if (sense == GRB_GREATER_EQUAL) {
			gravity_model.add(expr >= rhs);
		} else if (sense == GRB_EQUAL) {
			gravity_model.add(expr == rhs);
		} else {
			Warning("Unknown constraint sense: " << sense << std::endl);
			exit(1);
		}
	}
}

Model<> model_from_file(std::string path) {
	GRBEnv env = GRBEnv();
	GRBModel gurobi_model = GRBModel(env, path);
	Model<> gravity_model("GravityModel");

	int n_vars = gurobi_model.get(GRB_IntAttr_NumVars);
	int n_constraints = gurobi_model.get(GRB_IntAttr_NumConstrs);
	int n_int_vars = gurobi_model.get(GRB_IntAttr_NumIntVars);
	int n_bin_vars = gurobi_model.get(GRB_IntAttr_NumBinVars);

	add_vars(gurobi_model, gravity_model);

	// Make sure all variables were added
	if (gravity_model._vars.size() != n_vars) {
		Warning("Not all variables were added. Only Continuous, Integer, and Binary variables are supported." << std::endl);
		exit(1);
	}

	 add_constraints(gurobi_model, gravity_model);

	// Build objective
	 add_objective(gurobi_model, gravity_model);

	return gravity_model;
}
