
#include <fstream>
#include <spdlog/spdlog.h>
#include <sstream>
#include <string>
#include <map>
#include <set>

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

bool read_line(std::ifstream& infile, std::string& line) {
	// Read a line and return true if it is not a comment
	// Comments are lines that start with a '*' or a '$'
	while (std::getline(infile, line)) {
		if (line[0] != '*' && line[0] != '$') {
			return true;
		}
	}
	return false;
}

void read_name(std::ifstream& infile, MPSData& mps_data) {
	spdlog::info("Reading NAME");

	std::string line, _field;
	read_line(infile, line);

	if (line.find("NAME") == std::string::npos) {
		spdlog::error("NAME not found");
		exit(1);
	}

	std::istringstream iss{line};
	iss >> _field >> mps_data.name;
}

void read_obj_sense(std::ifstream& infile, MPSData& mps_data) {
	spdlog::info("Reading OBJSENSE");

	std::string line, _field;
	read_line(infile, line);

	if (line.find("OBJSENSE") != std::string::npos) {
		std::istringstream iss{line};
		iss >> _field >> mps_data.objsense;
	}
}

void read_rows(std::ifstream& infile, MPSData& mps_data) {
	spdlog::info("Reading ROWS");

	std::string line;
	std::string ctype, cname, priority, weight, reltol, abstol;

	// skip rows header
	read_line(infile, line);

	while (read_line(infile, line) && !line.empty() && line[0] == ' ')
	{
	    std::istringstream iss(line);
	   	iss >> ctype >> cname >> priority >> weight >> reltol >> abstol;
		spdlog::debug("ctype: {}, cname: {}", ctype, cname);

		// if f is not empty, then we are reading a multi-objective cost
		if (!abstol.empty()) {
			mps_data.multi_objective_rows.push_back(
				std::make_tuple(ctype, cname, priority, weight, reltol, abstol)
			);
		} else {
			mps_data.rows.push_back(std::make_tuple(ctype, cname));
		}
	}

	// Rewind to previous line
	infile.seekg(-(line.size() + 1), std::ios_base::cur);
}

void read_columns(std::ifstream& infile, MPSData& mps_data) {
	spdlog::info("Reading COLUMNS");

	std::string line;
	std::string col, row;
	float coeff;

	bool integrality = false;

	// skip columns header
	read_line(infile, line);
	while (read_line(infile, line) && !line.empty() && line[0] == ' ')
	{
		// Check if we are entering or exiting the integrality section
		if (line.find("'MARKER'") != std::string::npos) {
			if (line.find("'INTORG'") != std::string::npos) {
				integrality = true;
				continue;
			} else if (line.find("'INTEND'") != std::string::npos) {
				integrality = false;
				continue;
			}
		}

	    std::istringstream iss(line);
	   	iss >> col;

		// 0, 1, or 2 coefficients
		while (iss >> row >> coeff) {
			mps_data.columns[col].push_back(std::make_tuple(row, coeff, integrality));
		}
	}

	// Rewind to previous line
	infile.seekg(-(line.size() + 1), std::ios_base::cur);
}

void read_rhs(std::ifstream& infile, MPSData& mps_data) {
	spdlog::info("Reading RHS");

	std::string line;
	std::string _rhs_name, row;
	float coeff;

	// skip rhs header
	read_line(infile, line);
	while (read_line(infile, line) && !line.empty() && line[0] == ' ')
	{
	    std::istringstream iss(line);
	   	iss >> _rhs_name;

		// 1 or 2 rhs values
		while (iss >> row >> coeff) {
			mps_data.rhs.push_back(std::make_tuple(row, coeff));
		}
	}

	// Rewind to previous line
	infile.seekg(-(line.size() + 1), std::ios_base::cur);
}

void read_bounds(std::ifstream& infile, MPSData& mps_data) {
	spdlog::info("Reading BOUNDS");

	std::string line;
	std::string bound_type, _bound_name, var_name, bound_val;

	// skip bounds header
	read_line(infile, line);
	while (read_line(infile, line) && !line.empty() && line[0] == ' ')
	{
	    std::istringstream iss(line);
		iss >> bound_type >> _bound_name >> var_name >> bound_val;
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