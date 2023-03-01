
#include <fstream>
#include <spdlog/spdlog.h>
#include <sstream>
#include <string>
#include <map>
#include <set>

const std::set<std::string> ROW_TYPES = {
	"N", // free row, first instance is objective function
	"L", // less than or equal to
	"G", // greater than or equal to
	"E", // equal to
};

const std::set<std::string> BOUND_TYPES = {
	"LO", // lower bound
	"UP", // upper bound
	"FX", // fixed at specified value
	"FR", // free
	"MI", // infinite lower bound
	"PL", // infinite upper bound
	"BV", // binary variable
	"LI", // lower bound on integer variable
	"UI", // upper bound on integer variable
	"SC", // upper bound for semi-continuous variable
	"SI", // upper bound for semi-integer variable
};

class MPSReader {
public:
	MPSReader(std::string file_path) {
		this->file_path = file_path;
		this->infile.open(file_path);

		if (!this->infile.good()) {
			spdlog::error("Could not open file {}", file_path);
			exit(1);
		}
	}

	void parse() {
		std::string line;

		// Begin reading
		while (true)
		{
			line = this->read_rewind();

			if (line.find("NAME") == 0) {
				this->read_name();
				spdlog::info("Found NAME: {}", this->name);
			} else if (line.find("OBJSENSE") == 0) {
				this->read_obj_sense();
				spdlog::info("Found OBJSENSE: {}", this->objsense);
			} else if (line.find("ROWS") == 0) {
				read_rows();
			} else if (line.find("COLUMNS") == 0) {
				read_columns();
			} else if (line.find("RHS") == 0) {
				read_rhs();
			} else if (line.find("BOUNDS") == 0) {
				read_bounds();
			} else if (line.find("ENDATA") == 0) {
				break;
			} else {
				spdlog::error("Unsupported section: {}", line);
				exit(1);
			}
		}
	}

	// Name of the problem
	std::string name;
	// Optimization type (MIN or MAX)
	std::string objsense = "MIN";

	// Constraint type, Name
	std::vector<std::tuple<std::string, std::string>> rows;

	// Column -> (Row, Coeff, IsInteger)
	std::map<std::string, std::vector<std::tuple<std::string, double, bool>>> columns;

	// Row, RHS value
	std::vector<std::tuple<std::string, double>> rhs;

	std::vector<std::tuple<std::string, std::string, double>> bounds;

	// Name of objective function
	std::string obj_fn;
	// Offset for objective function
	double obj_fn_offset = 0.0;

private:
	void read_name();
	void read_obj_sense();
	void read_rows();
	void read_columns();
	void read_rhs();
	void read_bounds();

	bool read_line(std::string &line);
	std::string read_rewind();
	void rewind_line(std::string &line);

	std::string file_path;
	std::ifstream infile;
};

std::string MPSReader::read_rewind() {
	std::string line;
	std::getline(this->infile, line);
	this->rewind_line(line);
	return line;
}

void MPSReader::rewind_line(std::string &line) {
	infile.seekg(-(line.size() + 1), std::ios_base::cur);
}

bool MPSReader::read_line(std::string& line) {
	while (std::getline(this->infile, line)) {
		// MPS Comments
		if (line[0] != '*' && line[0] != '$') {
			return true;
		}
	}
	return false;
}

void MPSReader::read_name() {
	spdlog::info("Reading NAME");

	std::string line, _field;
	this->read_line(line);

	if (line.find("NAME") == std::string::npos) {
		spdlog::error("NAME not found");
		exit(1);
	}

	std::istringstream iss{line};
	iss >> _field >> this->name;
}

void MPSReader::read_obj_sense() {
	spdlog::info("Reading OBJSENSE");

	std::string line, _field;
	this->read_line(line);

	std::istringstream iss{line};
	iss >> _field >> this->objsense;
}

void MPSReader::read_rows() {
	spdlog::info("Reading ROWS");

	std::string line;
	bool found_obj = false;

	// skip rows header
	this->read_line(line);

	while (this->read_line(line) && !line.empty() && line[0] == ' ')
	{
		std::string ctype, cname;
		std::string priority, weight, reltol, abstol;

	    std::istringstream iss(line);
	   	iss >> ctype >> cname >> priority >> weight >> reltol >> abstol;

		if (ROW_TYPES.count(ctype) == 0) {
			spdlog::error("Unsupported row type: {}", ctype);
			exit(1);
		}

		// if f is not empty, then we are reading a multi-objective cost
		if (!abstol.empty()) {
			spdlog::error("Multi-objective costs not supported");
			exit(1);
		}
		if (ctype == "N" && !found_obj) {
			this->obj_fn = cname;
			found_obj = true;
		}

		this->rows.push_back(std::make_tuple(ctype, cname));
	}
	this->rewind_line(line);

	if (!found_obj) {
		spdlog::error("Could not find objective function. First free row interpreted as objective function");
		exit(1);
	}
}

void MPSReader::read_columns() {
	spdlog::info("Reading COLUMNS");

	std::string line;
	bool is_integer = false;

	// skip columns header
	this->read_line(line);
	while (read_line(line) && !line.empty() && line[0] == ' ')
	{
		// Check if we are entering or exiting the integrality section
		if (line.find("'MARKER'") != std::string::npos) {
			if (line.find("'INTORG'") != std::string::npos) {
				is_integer = true;
				continue;
			} else if (line.find("'INTEND'") != std::string::npos) {
				is_integer = false;
				continue;
			} else {
				spdlog::error("Malformed MARKER line: ", line);
				exit(1);
			}
		}

		std::string col, row;
		double coeff;

	    std::istringstream iss(line);
	   	iss >> col;

		// 0, 1, or 2 coefficients
		while (iss >> row >> coeff) {
			this->columns[col].push_back(std::make_tuple(row, coeff, is_integer));
		}
	}

	this->rewind_line(line);
}

void MPSReader::read_rhs() {
	spdlog::info("Reading RHS");

	std::string line;

	// skip rhs header
	this->read_line(line);
	while (read_line(line) && !line.empty() && line[0] == ' ')
	{
	    std::istringstream iss(line);
		std::string _rhs_name, row;
		double coeff;

	   	iss >> _rhs_name;

		// 1 or 2 rhs values
		while (iss >> row >> coeff) {
			if (row == this->obj_fn) {
				this->obj_fn_offset = coeff;
				continue;
			}

			this->rhs.push_back(std::make_tuple(row, coeff));
		}
	}

	// Rewind to previous line
	this->rewind_line(line);
}

void MPSReader::read_bounds() {
	spdlog::info("Reading BOUNDS");

	std::string line;
	std::string bound_type, _bound_name, var_name, bound_val;

	// skip bounds header
	this->read_line(line);
	while (this->read_line(line) && !line.empty() && line[0] == ' ')
	{
		spdlog::info("Line: {}", line);
	    std::istringstream iss(line);
		iss >> bound_type >> _bound_name >> var_name >> bound_val;
		spdlog::info("Bound type: {}, Bound name: {}, Var name: {}, Bound val: {}", bound_type, _bound_name, var_name, bound_val);

		if (BOUND_TYPES.count(bound_type) == 0) {
			spdlog::error("Unknown bound type: {}", bound_type);
			exit(1);
		}

		this->bounds.push_back(std::make_tuple(bound_type, var_name, std::stod(bound_val)));
	}

	// Rewind to previous line
	this->infile.seekg(-(line.size() + 1), std::ios_base::cur);
}

int main(int argc, char** argv) {
	MPSReader mps_data("/home/haydnj/Documents/Gravity/out.mps");
	mps_data.parse();

	return 0;
}