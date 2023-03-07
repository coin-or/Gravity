//
//  DCOPF.cpp
//  Gravity
//
//  Created by Hassan Hijazi on 19 Jan 18.
//
//
#include <stdio.h>
#include <PowerNet.h>
#include <gravity/solver.h>
#ifdef USE_OPT_PARSER
#include <optionParser.hpp>
#endif
#include <gravity/rapidcsv.h>
using namespace std;
using namespace gravity;

int main (int argc, char * argv[])
{
    string fname = string(prj_dir)+"/data_sets/Roster/Enrollment_Details_23.csv";
    
    
#ifdef USE_OPT_PARSER
    /** Create a OptionParser with options */
    auto options = readOptions(argc, argv);
    options.add_option("f", "file", "Input file name", fname);
    
    /** Parse the options and verify that all went well. If not, errors and help will be shown */
    bool correct_parsing = options.parse_options(argc, argv);
    
    if(!correct_parsing){
        return EXIT_FAILURE;
    }
    
    fname = options["f"];
    bool has_help = op::str2bool(options["h"]);
    if(has_help) {
        options.show_help();
        exit(0);
    }
#else
    if(argc>=2){
        fname=argv[1];
    }
#endif
    rapidcsv::Document doc(fname);
    map<string, int> LA_div_sizes, WR_div_sizes;
    cout << "There are " << doc.GetRowCount() << " registered players." << std::endl;
    vector<int> ages = doc.GetColumn<int>("Player Age");
    vector<string> emails_vec = doc.GetColumn<string>("User Email");
    vector<string> loc_full = doc.GetColumn<string>(doc.GetColumnCount()-2);
    vector<string> loc_season = doc.GetColumn<string>(doc.GetColumnCount()-3);
    set<string> emails;
    for(auto i = 0; i<ages.size(); i++){
        int age = ages[i];
        emails.insert(emails_vec[i]);
        if(loc_full[i]=="Los Alamos" || loc_season[i]=="Los Alamos")
            LA_div_sizes["U"+to_string(age+1)]++;
        else if(loc_full[i]=="White Rock" || loc_season[i]=="White Rock")
            WR_div_sizes["U"+to_string(age+1)]++;
        else
            cout << "WARNING: preferred location unspecified, not sure where to place player! Row number: " << i << endl;
    }
    DebugOn("There are " << LA_div_sizes.size() << " divisions in LA " << endl);
    for(const auto &p: LA_div_sizes){
        DebugOn(p.first << " : " << p.second << endl);
    }
    DebugOn("There are " << WR_div_sizes.size() << " divisions in WR " << endl);
    for(const auto &p: WR_div_sizes){
        DebugOn(p.first << " : " << p.second << endl);
    }
    bool print_emails = true;
    if(print_emails){
        DebugOn("Printing email lists in batches of 100 recipients or less:" <<endl);
        int idx = 1;
        for(const auto &em: emails){
            if(idx%100==0)
                DebugOn("\nNEW BATCH\n");
            DebugOn(em<<";");
            idx++;
        }
    }
    DebugOn("\nDone parsing enrollment file" <<endl);
    return 0;
}
