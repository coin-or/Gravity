
#include <fstream>
#include <map>
#include <set>
#include <spdlog/spdlog.h>
#include <sstream>
#include <string>
#include <gurobi_c++.h>
#include <gravity/model.h>
#include <gravity/func.h>
#include <map>

using namespace gravity;

int main(int argc, char** argv) {
	std::map<std::string, var<double>> cont_vars;
	std::map<std::string, var<int>> int_vars;
	std::map<std::string, var<bool>> bin_vars;
	std::map<std::string, Constraint<double>> constraints;

	try {
		GRBEnv env = GRBEnv();
    	GRBModel m = GRBModel(env, "/mnt/trail_test/haydnj/vnn/Gravity/simple.mps");

		Model<> gravityModel("GravityModel");

		int n_vars = m.get(GRB_IntAttr_NumVars);
		int n_constraints = m.get(GRB_IntAttr_NumConstrs);
		int n_int_vars = m.get(GRB_IntAttr_NumIntVars);
		int n_bin_vars = m.get(GRB_IntAttr_NumBinVars);
		spdlog::info("NumVars: {}, NumConstraints: {}", n_vars, n_constraints);
		spdlog::info("NumIntVars: {}, NumBinVars: {}", n_int_vars, n_bin_vars);

		// Print variable names
		for (int i = 0; i < n_vars; i++) {
			spdlog::info("Getting var {}/{}", i, n_vars);

			GRBVar gurobi_var = m.getVar(i);
			const double lb = gurobi_var.get(GRB_DoubleAttr_LB);
			const double ub = gurobi_var.get(GRB_DoubleAttr_UB);
			const char type = gurobi_var.get(GRB_CharAttr_VType);
			const std::string name = gurobi_var.get(GRB_StringAttr_VarName);

			// spdlog::info("\tVarName: {}, LB: {}, UB: {}, Type: {}", name, lb, ub, type);
			if (type == GRB_CONTINUOUS) {
				cont_vars.insert({name, {name, lb, ub}});
				gravityModel.add(cont_vars[name]);
			} else if (type == GRB_INTEGER) {
				int_vars.insert({name, {name, (int)lb, (int)ub}});
				gravityModel.add(cont_vars[name]);
			} else if (type == GRB_BINARY) {
				bin_vars.insert({name, {name, (bool)lb, (bool)ub}});
				gravityModel.add(bin_vars[name]);
			} else {
				spdlog::error("Unknown variable type: {}", type);
				exit(1);
			}
		}

		// Make sure all variables were added
		if (cont_vars.size() + int_vars.size() + bin_vars.size() != n_vars) {
			spdlog::error("Not all variables were added.");
			exit(1);
		}

		spdlog::info("Done adding variables.");

		// Print Gravity variables
		spdlog::info("~~~~~~~~ ContVars ~~~~~~~~");
		for (auto& v : cont_vars) {
			spdlog::info("{}: [{}, {}]", v.first, v.second.get_lb(0), v.second.get_ub(0));
		}

		spdlog::info("~~~~~~~~ IntVars ~~~~~~~~");
		for (auto& v : int_vars) {
			spdlog::info("{}: [{}, {}]", v.first, v.second.get_lb(0), v.second.get_ub(0));
		}

		spdlog::info("~~~~~~~~ BinVars ~~~~~~~~");
		for (auto& v : bin_vars) {
			spdlog::info("{}: [{}, {}]", v.first, v.second.get_lb(0), v.second.get_ub(0));
		}

		// Print constraint names
		for (int i = 0; i < n_constraints; i++) {
			GRBConstr constr = m.getConstr(i);
			std::string name = constr.get(GRB_StringAttr_ConstrName);
			char sense = constr.get(GRB_CharAttr_Sense);
			double rhs = constr.get(GRB_DoubleAttr_RHS);

			Constraint<double> gravityConstraint(name);

			for (int j = 0; j < n_vars; j++) {
				GRBVar var = m.getVar(j);
				double coeff = m.getCoeff(constr, var);
				if (coeff == 0.0) {
					continue;
				}

				std::string var_name = var.get(GRB_StringAttr_VarName);
				if (cont_vars.find(var_name) != cont_vars.end()) {
					gravityConstraint += coeff * cont_vars[var_name];
				} else if (int_vars.find(var_name) != int_vars.end()) {
					gravityConstraint += coeff * int_vars[var_name];
				} else if (bin_vars.find(var_name) != bin_vars.end()) {
					gravityConstraint += coeff * bin_vars[var_name];
				} else {
					spdlog::error("Unknown variable name: {}", var_name);
					exit(1);
				}
			}

			if (sense == GRB_LESS_EQUAL) {
				gravityModel.add(gravityConstraint <= rhs);
			} else if (sense == GRB_GREATER_EQUAL) {
				gravityModel.add(gravityConstraint >= rhs);
			} else if (sense == GRB_EQUAL) {
				gravityModel.add(gravityConstraint == rhs);
			} else {
				spdlog::error("Unknown constraint sense: {}", sense);
				exit(1);
			}
		}

		// Build objective
		GRBQuadExpr quad_obj = m.getObjective();
		GRBLinExpr lin_obj = quad_obj.getLinExpr();
		unsigned int n_terms = lin_obj.size();

		gravity::func<> gravityObj;

		for (int i = 0; i < n_terms; i++) {
			GRBVar var = lin_obj.getVar(i);
			double coeff = lin_obj.getCoeff(i);
			std::string var_name = var.get(GRB_StringAttr_VarName);

			if (cont_vars.find(var_name) != cont_vars.end()) {
				gravityObj += coeff * cont_vars[var_name];
			} else if (int_vars.find(var_name) != int_vars.end()) {
				gravityObj += coeff * int_vars[var_name];
			} else if (bin_vars.find(var_name) != bin_vars.end()) {
				gravityObj += coeff * bin_vars[var_name];
			} else {
				spdlog::error("Unknown variable name: {}", var_name);
				exit(1);
			}
		}

		// Get model objsense
		int objsense = m.get(GRB_IntAttr_ModelSense);
		if (objsense == GRB_MAXIMIZE) {
			gravityModel.max(gravityObj);
		} else if (objsense == GRB_MINIMIZE) {
			gravityModel.min(gravityObj);
		} else {
			spdlog::error("Unknown objsense: {}", objsense);
			exit(1);
		}

		// Add objective
		gravityModel.print();
	} catch (GRBException e) {
		spdlog::error("Error code = {}", e.getErrorCode());
		spdlog::error("{}", e.getMessage());
	}
	return 0;
}