#ifndef OPTION_PARSER_H
#define OPTION_PARSER_H

#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <vector>

#include <typeinfo>

using namespace std;



namespace op{
   
/* convert a boolean to a string */
string bool2str(bool value);

/* convert a string to a boolean */
bool str2bool(string value);

/* convert a string to an int */
int str2int(string value);

/* convert a string to an int */
float str2float(string value);

/* convert a string to an int */
double str2double(string value);


/* CLASS OPTION  ===================================

	An option is defined by a short and a long name
	The - and -- in the option names are added automatically

*/
class Option{
protected:
	string var, default_value;
	string short_name, long_name, description;
	bool _is_boolean;

public:
	Option(string short_name, string long_name, string description, string default_value);

	/* create a boolean option */
	Option(string short_name, string long_name, string description, bool default_value);

	~Option(){}
    
    string get_short_name()const{return short_name;}
    string get_long_name()const{return long_name;}
    string get_description()const{return description;}
    string get_default() const{ return default_value;}
    string get_var() const{ return var;}
    void set_var(string& new_var){ var = new_var;}
    void set_var(bool new_val){ var = bool2str(new_val);}
    bool is_boolean(){return _is_boolean;}

};







/* CLASS OPTIONPARSER  ===================================

	This class implements a very easy to use c++ option parser
	Options can be added and then their value accessed with [] operator
	All options values are string, the user must convert them himself

*/
class OptionParser{

protected:
    vector<Option> options;
    int max_size_short; // max size of short names of options, usefull for a beautifull help print
    int max_size_long;  // max size of long names of options, usefull for a beautifull help print
    int max_size_desc;  // max size of desciption of options, usefull for a beautifull help print
    int description_max_line_size; // maximum size of a descritption line, if reached, a line is skipped


    bool has_option(string name);
    int get_description_max_line_size() const{ return description_max_line_size;}
    void set_description_max_line_size(int size){ description_max_line_size = size;}


	// returns a string of n spaces
    std::string repeat_str( const std::string &word, int times );

    // returns max
    int max(int a, int b);


public:
	OptionParser(int max_size_desc_line = 100);
	~OptionParser();

	/* Shows help with aligned columns */	
	void show_help();

	void add_option(string short_name, string long_name, string description, string default_value);

	/* to add a boolean option */	
	void add_option(string short_name, string long_name, string description);

	/* This function changes default option values to those which are read on the command line */	
	bool parse_options(int& argc, char**& argv);

	/* Get the option values as string */	
	string operator[](string name) const;


};


} // end namespace op



#endif