//
//  Roster_main.cpp
//
//  This program can read in csv files with player info and create balanced teams
//  The program takes the path to the CSV file as a first argument.
//  If you wish to specify the min and max team sizes, you can add them as second (for min) and third arguments (max).
//  For example, you can type in: ./roster ../../data_sets/Roster/Enrollment_Details_Example.csv 10 20
//  The default min/max is set to 8/16.
//
//  Gravity
//
//  Created by Hassan Hijazi on 8 March 2023.
//
//
#include <stdio.h>
#include <PowerNet.h>
#include <gravity/solver.h>
#ifdef USE_OPT_PARSER
#include <optionParser.hpp>
#endif
#include "Roster.hpp"
using namespace std;
using namespace gravity;

int main (int argc, char * argv[])
{
    string fname = string(prj_dir)+"/data_sets/Roster/Enrollment_Details_Example.csv";
    int min_team_size = 8, max_team_size = 16;
    cout << "Welcome to the rostering program!\nThe program takes the path to the CSV file as a first argument.\nIf you wish to specify the min and max team sizes, you can add them as second (for min) and third arguments (max).\n For example, you can type in: ./roster ../../data_sets/Roster/Enrollment_Details_Example.csv 10 20\n The default min/max is set to 8/16.";
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
    if(argc>=4){
        min_team_size=stoi(argv[2]);
        max_team_size=stoi(argv[3]);
    }
#endif
    rapidcsv::Document doc(fname);
    map<int, vector<Team>> LA_divs_sorted, WR_divs_sorted;
    auto nb_players = doc.GetRowCount();
    cout << "There are " <<  nb_players << " registered players." << std::endl;
    vector<int> ages = doc.GetColumn<int>("Player Age");
    vector<string> fnames = doc.GetColumn<string>("Player First Name");
    vector<string> lnames = doc.GetColumn<string>("Player Last Name");
    vector<string> genders = doc.GetColumn<string>("Player Gender");
    vector<string> emails_vec = doc.GetColumn<string>("User Email");
    vector<string> loc_full = doc.GetColumn<string>(doc.GetColumnCount()-2);
    vector<string> loc_season = doc.GetColumn<string>(doc.GetColumnCount()-3);
    set<string> emails;
    for(auto i = 0; i<nb_players; i++){
        int age = ages[i];
        Player p(fnames[i], lnames[i], ages[i], emails_vec[i], genders[i]);
        emails.insert(emails_vec[i]);
        if(loc_full[i]=="Los Alamos" || loc_season[i]=="Los Alamos"){
            int div_num = age+1;
            string team_name = "LA_U"+to_string(div_num);
            if(LA_divs_sorted.count(div_num)==0)
                LA_divs_sorted[div_num].push_back(Team(team_name));
            LA_divs_sorted[div_num].back().players.push_back(p);
        }
        else if(loc_full[i]=="White Rock" || loc_season[i]=="White Rock"){
            int div_num = age+1;
            string team_name = "WR_U"+to_string(div_num);
            if(WR_divs_sorted.count(div_num)==0)
                WR_divs_sorted[div_num].push_back(Team(team_name));
            WR_divs_sorted[div_num].back().players.push_back(p);
        }
        else
            cout << "WARNING: preferred location unspecified, not sure where to place player! Row number: " << i << endl;
    }
    DebugOn("There are " << LA_divs_sorted.size() << " divisions in LA " << endl);
    DebugOn("There are " << WR_divs_sorted.size() << " divisions in WR " << endl);
    

    map<string, vector<Team>> LA_divs, WR_divs;
    for(auto const &div: LA_divs_sorted){
        LA_divs["U"+to_string(div.first)] = div.second;
    }
    
    for(auto const &div: WR_divs_sorted){
        WR_divs["U"+to_string(div.first)] = div.second;
    }
    
    /* First, combining small teams */
    
    bool min_size_ok = false;
    while(!min_size_ok){
        map<string, vector<Team>> new_LA_divs = LA_divs;
        for(auto const &div: LA_divs){
            DebugOn("LA Division " << div.first);
            Team t =  div.second.back();
            DebugOn(" has " << t.players.size() << " players" << endl);
            if(t.players.size()<min_team_size){
                DebugOn("Merging teams to satisfy min size constraint" << endl);
                string div_name = div.first;
                int div_num = stoi(div_name.substr(div_name.find_first_of("U")+1));
                string div_up_name = "U"+to_string(div_num+1);
                auto it = LA_divs.find(div_up_name);
                if(it!=LA_divs.end()){
                    Team new_team = t.combine_with(it->second.back());
                    new_team.team_name = "LA_"+div_name+"_"+div_up_name;
                    new_LA_divs[div_name+"_"+div_up_name].push_back(new_team);
                    new_LA_divs.erase(div_name);
                    new_LA_divs.erase(div_up_name);
                }
                else{
                    string div_down_name = "U"+to_string(div_num-1);
                    auto it = LA_divs.find(div_down_name);
                    if(it==LA_divs.end()){
                        goto endloop;
                    }
                    Team new_team = t.combine_with(it->second.back());
                    new_team.team_name = "LA_"+div_down_name+"_"+div_name;
                    new_LA_divs[div_down_name+"_"+div_name].push_back(new_team);
                    new_LA_divs.erase(div_name);
                    new_LA_divs.erase(div_down_name);
                }
                break;
            }
        }
        LA_divs = new_LA_divs;
        min_size_ok = true;
        for(auto const &div: LA_divs){
            Team t =  div.second.back();
            if(t.players.size()<min_team_size){
                min_size_ok = false;
            }
        }
    }
    min_size_ok = false;
    while(!min_size_ok){
        map<string, vector<Team>> new_WR_divs = WR_divs;
        for(auto const &div: WR_divs){
            DebugOn("WR Division " << div.first);
            Team t =  div.second.back();
            DebugOn(" has " << t.players.size() << " players" << endl);
            if(t.players.size()<min_team_size){
                DebugOn("Merging teams to satisfy min size constraint" << endl);
                string div_name = div.first;
                int div_num = stoi(div_name.substr(div_name.find_first_of("U")+1));
                string div_up_name = "U"+to_string(div_num+1);
                auto it = WR_divs.find(div_up_name);
                if(it!=WR_divs.end()){
                    Team new_team = t.combine_with(it->second.back());
                    new_team.team_name = "WR_"+div_name+"_"+div_up_name;
                    new_WR_divs[div_name+"_"+div_up_name].push_back(new_team);
                    new_WR_divs.erase(div_name);
                    new_WR_divs.erase(div_up_name);
                }
                else{
                    string div_down_name = "U"+to_string(div_num-1);
                    auto it = WR_divs.find(div_down_name);
                    if(it==WR_divs.end()){
                        goto endloop;
                    }
                    Team new_team = t.combine_with(it->second.back());
                    new_team.team_name = "WR_"+div_down_name+"_"+div_name;
                    new_WR_divs[div_down_name+"_"+div_name].push_back(new_team);
                    new_WR_divs.erase(div_name);
                    new_WR_divs.erase(div_down_name);

                }
                break;
            }
        }
        WR_divs = new_WR_divs;
        min_size_ok = true;
        for(auto const &div: WR_divs){
            Team t =  div.second.back();
            if(t.players.size()<min_team_size){
                min_size_ok = false;
            }
        }
    }
    endloop:
    for(auto const &div: LA_divs){
        DebugOn("LA Division " << div.first);
        for(Team t: div.second){
            DebugOn(" has " << t.players.size() << " players" << endl);
            if(t.players.size()>max_team_size){
                DebugOn("Splitting team to satisfy max size constraint" << endl);
                int nb_factors = ceil(t.players.size()/(double)max_team_size);
                vector<Team> split_t = t.split_into(nb_factors);
                LA_divs[div.first] = split_t;
            }
        }
    }
    for(auto const &div: WR_divs){
        DebugOn("WR Division " << div.first);
        for(Team t: div.second){
            DebugOn(" has " << t.players.size() << " players" << endl);
            if(t.players.size()>max_team_size){
                DebugOn("Splitting team to satisfy max size constraint" << endl);
                int nb_factors = ceil(t.players.size()/(double)max_team_size);
                vector<Team> split_t = t.split_into(nb_factors);
                WR_divs[div.first] = split_t;
            }
        }
    }
    for(auto const &div: LA_divs){
        DebugOn("\n\nLA Division " << div.first << endl);
        for(auto &t: div.second){
            DebugOn("\n\tLA Team " << t.team_name << endl);
            DebugOn("\t " << t.players.size() << " players" << endl);
            for(auto const &p: t.players){
                DebugOn("\t\t" << p.first_name << " " << p.last_name  << endl);
            }
            DebugOn("Email list for team:\n");
            set<string> unique_emails;
            for(auto const &p: t.players){
                if(unique_emails.insert(p.email).second){
                    DebugOn(p.email << ";");
                }
            }
            DebugOn(endl);
            t.save_team();
            DebugOn("Done writing csv file for team\n");
        }
    }
    for(auto const &div: WR_divs){
        DebugOn("\n\nWR Division " << div.first << endl);
        for(auto const &t: div.second){
            DebugOn("\n\tWR Team " << t.team_name << endl);
            DebugOn("\t " << t.players.size() << " players" << endl);
            for(auto const &p: t.players){
                DebugOn("\t\t" << p.first_name << " " << p.last_name  << endl);
            }
            DebugOn("Email list for team:\n");
            set<string> unique_emails;
            for(auto const &p: t.players){
                if(unique_emails.insert(p.email).second){
                    DebugOn(p.email << ";");
                }
            }
            DebugOn(endl);
            t.save_team();
            DebugOn("Done writing csv file for team\n");
        }
    }
    bool print_emails = true;
    if(print_emails){
        DebugOn("\n\nPrinting full email lists in batches of 100 recipients or less:" <<endl);
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
