
#include <fstream>
#include <spdlog/spdlog.h>
#include <sstream>
#include <string>
#include <map>

class MPSData {
public:
	std::string name;
	std::string objsense = "MIN";

	std::vector<std::tuple<std::string, std::string>> rows;
	// vector of tuple of 6 strings for multi-objective costs
	std::vector<std::tuple<std::string, std::string, std::string, std::string, std::string, std::string>> multi_objective_rows;

	// Column -> (Row, Coeff, IsInteger)
	std::map<std::string, std::vector<std::tuple<std::string, float, bool>>> columns;

	// RHS
	std::vector<std::tuple<std::string, float>> rhs;

	MPSData() {}
};

void read_name(std::ifstream& infile, MPSData& mps_data) {
	std::string line, _field;
	std::getline(infile, line);

	if (line.find("NAME") == std::string::npos) {
		spdlog::error("NAME not found");
		exit(1);
	}

	std::istringstream iss{line};
	iss >> _field >> mps_data.name;
}

void read_obj_sense(std::ifstream& infile, MPSData& mps_data) {
	std::string line, _field;
	std::getline(infile, line);

	if (line.find("OBJSENSE") == std::string::npos) {
		spdlog::error("OBJSENSE not found");
		exit(1);
	}

	std::istringstream iss{line};
	iss >> _field >> mps_data.objsense;
}

void read_rows(std::ifstream& infile, MPSData& mps_data) {
	std::string line;
	std::string a, b, c, d, e, f;

	// skip rows header
	std::getline(infile, line);

	while (std::getline(infile, line))
	{
		if (
			(line.find("COLUMNS") == 0)  ||
			(line.find("LAZYCONS") == 0) ||
			(line.find("USERCUTS") == 0)
		) {
			break;
		}
		// a: constraint type
		// b: constraint name
	    std::istringstream iss(line);
	   	iss >> a >> b >> c >> d >> e >> f;

		// if f is not empty, then we are reading a multi-objective cost
		if (!f.empty()) {
			mps_data.multi_objective_rows.push_back(std::make_tuple(a, b, c, d, e, f));
		} else {
			mps_data.rows.push_back(std::make_tuple(a, b));
		}
	}

	// Rewind to previous line
	infile.seekg(-(line.size() + 1), std::ios_base::cur);
}

void read_columns(std::ifstream& infile, MPSData& mps_data) {
	std::string line;
	std::string col_name, row;
	float coeff;

	bool integrality = false;
	while (std::getline(infile, line) && line.find("RHS") != 0)
	{
		// Check if we are entering or exiting the integrality section
		if (line.find("'MARKER'") != std::string::npos) {
			if (line.find("'INTORG'") != std::string::npos) {
				spdlog::info("Entering integrality section");
				integrality = true;
				continue;
			} else if (line.find("'INTEND'") != std::string::npos) {
				spdlog::info("Exiting integrality section");
				integrality = false;
				continue;
			}
		}

	    std::istringstream iss(line);
	   	iss >> col_name;

		spdlog::info("Column: {}", col_name);

		// 0, 1, or 2 coefficients
		while (iss >> row >> coeff) {
			mps_data.columns[col_name].push_back(std::make_tuple(row, coeff, integrality));
			spdlog::info("\tRow: {}, coeff: {}, Integrality: {}", row, coeff, integrality);
		}
	}

	// Rewind to previous line
	infile.seekg(-(line.size() + 1), std::ios_base::cur);
}

void read_rhs(std::ifstream& infile, MPSData& mps_data) {
	std::string line;
	std::string _rhs_name, row;
	float coeff;

	// skip rhs header
	std::getline(infile, line);
	while (std::getline(infile, line))
	{
		if (
			(line.find("BOUNDS") == 0) ||
			(line.find("ENDATA") == 0)
		) {
			break;
		}

	    std::istringstream iss(line);
	   	iss >> _rhs_name;

		spdlog::info("RHS: {}", _rhs_name);
		// 1 or 2 rhs values
		while (iss >> row >> coeff) {
			mps_data.rhs.push_back(std::make_tuple(row, coeff));
			spdlog::info("\tRow: {}, coeff: {}", row, coeff);
		}
	}

	// Rewind to previous line
	infile.seekg(-(line.size() + 1), std::ios_base::cur);
}

void read_bounds(std::ifstream& infile, MPSData& mps_data) {
	std::string line;
	std::string bound_type, _bound_name, var_name, bound_val;

	// skip bounds header
	std::getline(infile, line);
	while (std::getline(infile, line))
	{
		if (
			(line.find("QUADOBJ") == 0)    ||
			(line.find("QCMATRIX") == 0)   ||
			(line.find("PWLOBJ") == 0)     ||
			(line.find("SOS") == 0)        ||
			(line.find("INDICATORS") == 0) ||
			(line.find("GENCONS") == 0)    ||
			(line.find("SCENARIOS") == 0)  ||
			(line.find("ENDATA") == 0)
		) {
			break;
		}

	    std::istringstream iss(line);
		iss >> bound_type >> _bound_name >> var_name >> bound_val;
		spdlog::info("Bound: {}, {}, {}, {}", bound_type, _bound_name, var_name, bound_val);
	}

	// Rewind to previous line
	infile.seekg(-(line.size() + 1), std::ios_base::cur);
}

std::string read_rewind(std::ifstream& infile) {
	std::string line;
	std::getline(infile, line);
	infile.seekg(-(line.size() + 1), std::ios_base::cur);
	return line;
}

int main(int argc, char** argv) {
	MPSData mps_data;

	std::string line;
	std::ifstream infile;
	infile.open("/home/haydnj/Documents/Gravity/out.mps");
	if (!infile.good()) {
		spdlog::error("File badd :(");
		return -1;
	}

	// Begin reading
	while (true)
	{
		line = read_rewind(infile);
		spdlog::info("Line: {}", line);

		if (line.find("NAME") == 0) {
			read_name(infile, mps_data);
			spdlog::info("Found NAME: {}", mps_data.name);
		} else if (line.find("OBJSENSE") == 0) {
			read_obj_sense(infile, mps_data);
			spdlog::info("Found OBJSENSE: {}", mps_data.objsense);
		} else if (line.find("ROWS") == 0) {
			read_rows(infile, mps_data);
		} else if (line.find("COLUMNS") == 0) {
			read_columns(infile, mps_data);
		} else if (line.find("RHS") == 0) {
			read_rhs(infile, mps_data);
		} else if (line.find("BOUNDS") == 0) {
			read_bounds(infile, mps_data);
		} else if (line.find("ENDATA") == 0) {
			break;
		} else {
			spdlog::error("Unsupported section: {}", line);
			return -1;
		}
	}

	return 0;
}