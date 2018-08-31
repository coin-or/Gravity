// C++ OPTION PARSER TOOL
// creator : Nerzadler

#include <iostream>
#include "optionParser.hpp"

using namespace std;


int main(int argc, char* argv[]){

	// create a OptionParser with options
	op::OptionParser opt;
	opt.add_option("h", "help", "shows option help"); // no default value means boolean options, which default value is false
	opt.add_option("w", "window_size", "window's size", "256" );
	opt.add_option("r", "rate", "learning rate", "0.01" );
	opt.add_option("m", "mode", "learning mode", "random" );

	// parse the options and verify that all went well. If not, errors and help will be shown
	bool correct_parsing = opt.parse_options(argc, argv);
    
    if(!correct_parsing){
        return EXIT_FAILURE;
    }

    cout << opt["rate"] << endl;
    cout << opt["w"] << endl;
    cout << opt["h"] << endl;
    
    // you have to convert yourself arguments if needed. This gives you full control over your options
    int window_size     =  op::str2int(opt["w"]);
    float learning_rate =  op::str2float(opt["r"]);
    bool has_help       =  op::str2bool(opt["h"]);
    string mode         =  opt["m"];

    // show help
    if(has_help) opt.show_help();
    
	return 0;
}