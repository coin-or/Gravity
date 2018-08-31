#include "optionParser.hpp"

using namespace op;


string op::bool2str(bool value){
 return (value)? "true": "false";
}

bool op::str2bool(string value){
 return (value == "true" || value == "TRUE" );
}

int op::str2int(string value){
 return atoi(value.c_str());
}

float op::str2float(string value){
 return (float)atof(value.c_str());
}

double op::str2double(string value){
 return atof(value.c_str());
}


op::Option::Option(string short_name, string long_name, string description, string default_value){
	this->short_name = "-" + short_name;
	this->long_name = "--" + long_name;
    this->description = description;
    this->default_value = default_value;
	var = default_value;
	_is_boolean = false;
}

op::Option::Option(string short_name, string long_name, string description, bool default_value){
	this->short_name = "-" + short_name;
	this->long_name = "--" + long_name;
    this->description = description;
    this->default_value = bool2str(default_value);
	var = bool2str(default_value);
	_is_boolean = true;
}



bool op::OptionParser::has_option(string name){
    for( int i=0; i<options.size(); i++ )
        
		if(options[i].get_short_name() == "-" + name || options[i].get_long_name() == "--" + name){
			// an argument name has been recognised
			return true;				
		}
		// else option_recognised stays to false
	return false;
}

std::string op::OptionParser::repeat_str( const std::string &word, int times ) {
   std::string result;
   result.reserve(times*word.length()); // avoid repeated reallocation
   for ( int a = 0 ; a < times ; a++ ) 
      result += word ;
   return result;
}


int op::OptionParser::max(int a, int b){
	return (a<b)? b: a;
}


op::OptionParser::OptionParser(int max_size_desc_line){
	max_size_short = max_size_long = max_size_desc = 0;
	description_max_line_size = max_size_desc_line;
}

op::OptionParser::~OptionParser(){
}


void op::OptionParser::show_help(){
	string sn = "Short_name";
	string ln = "Long_name";
	string desc = "Description";
	// fisrt calculate the max size of each printed column :
	int max_size_s = max(max_size_short, sn.size());
	int max_size_l = max(max_size_long, ln.size());
	int max_size_d = max(max_size_desc, desc.size());
	//max_size_d = min(max_size_d, description_max_line_size); // a description line has a limit

	// then print nicely the help
	string separator = repeat_str("=", max_size_s + max_size_l + max_size_d + 26);

	cout << endl << "Help :" << endl;
	cout << separator << endl;


	string sn_spaces = repeat_str( " " , max_size_s - sn.size() +4 ); // garanty 2 spaces at least
	string ln_spaces = repeat_str( " " , max_size_l - ln.size() +4 ); // garanty 2 spaces at least
	string desc_spaces = repeat_str( " " , max_size_d - desc.size() +4 ); // garanty 2 spaces at least

	// print the title line
	cout << sn << sn_spaces << ln << ln_spaces << desc << desc_spaces << "Default_value" << endl;
	cout << separator << endl;

	// print the options
	for (vector<Option>::iterator it=options.begin(); it!=options.end(); it++){
		sn = it->get_short_name();
		ln = it->get_long_name();
		desc = it->get_description();

		sn_spaces = repeat_str( " " , max_size_s - sn.size() +4 ); // garanty 2 spaces at least
		ln_spaces = repeat_str( " " , max_size_l - ln.size() +4 ); // garanty 2 spaces at least
		desc_spaces = repeat_str( " " , max_size_d - desc.size() +4 ); // garanty 2 spaces at least

		cout << sn << sn_spaces << ln << ln_spaces << desc << desc_spaces << it->get_default() << endl;
	}

	cout << separator << endl;
	cout << endl;
}

void op::OptionParser::add_option(string short_name, string long_name, string description, string default_value){
	if(has_option(short_name) || has_option(long_name))
		throw string("Error when adding option " + long_name + ". Option already exists !");
	
	Option opt(short_name, long_name, description, default_value);
	options.push_back(opt);
	// update max sizes for a nice help printing
	max_size_short = max(max_size_short, short_name.size());
	max_size_long = max(max_size_long, long_name.size());
	max_size_desc = max(max_size_desc, description.size());
}

// to add a boolean option
void op::OptionParser::add_option(string short_name, string long_name, string description){
	if(has_option(short_name) || has_option(long_name))
		throw string("Error when adding option " + long_name + ". Option already exists !");

	Option opt(short_name, long_name, description, false);
	options.push_back(opt);
	// update max sizes for a nice help printing
	max_size_short = max(max_size_short, short_name.size());
	max_size_long = max(max_size_long, long_name.size());
	max_size_desc = max(max_size_desc, description.size());
}



bool op::OptionParser::parse_options(int& argc, char**& argv){
	// for each argument, check all options
	for(int i = 1; i<argc; i++){ // begin at 1 because of programm name
		string arg(argv[i]);
		bool option_recognised = false;
		for( vector<Option>::iterator it=options.begin(); it!=options.end(); it++ ){
			if(it->get_short_name() == arg || it->get_long_name() == arg){
				// an argument name has been recognised
				option_recognised = true;

				if( it->is_boolean() ){
					it->set_var(true);
				}
				else{
                    i++; // increase i
					// check if there is an argument after
					if(i==argc || has_option(argv[i]) ){
						cout << "Non boolean option " << arg << " needs an argument !" << endl;
						return false;
					}
                    string var(argv[i]);
					it->set_var(var);
				}
			}
			// else option_recognised stays to false
		}
		// if option is not recognised
		if(!option_recognised){
			cout << "Error : Unknown option " << arg << endl;
			show_help();
			return false;
		}
	}
	return true;
}



string op::OptionParser::operator[](string name) const{
	bool option_recognised = false;
	string var;

    for( int i=0; i<options.size(); i++ )
        
		if( (options[i].get_short_name()).substr(1) == name || (options[i].get_long_name()).substr(2) == name){
			// an argument name has been recognised
			option_recognised = true;
			var = options[i].get_var();
		}
		// else option_recognised stays to false

	// if option is not recognised
	if(!option_recognised){
		throw string( "Unknown option " + name );
	}
	
	return var;
}




