//
//  Roster.hpp
//
//  Created by Hassan on March 07 2023.

#ifndef Roster_hpp
#define Roster_hpp

#include <stdio.h>
#include <math.h>
#include <gravity/rapidcsv.h>

class Player {

public:
    /* First name */
    string first_name;
    /* Last name */
    string last_name;
    /* Age */
    int age;
    /* email */
    string email;
    /* Gender */
    string gender;
    
    Player(const string& fname, const string& lname, int _age, const string& _email, const string& _gender):first_name(fname), last_name(lname), age(_age), email(_email), gender(_gender){};
    
    
};

class Team {

public:
    /* Team name */
    string team_name;
    /* Head coach name */
    string coach_name;
    /* Head coach email */
    string coach_email;
    /* Assistant coach name */
    string assist_coach_name;
    /* Assistant coach name */
    string assist_coach_email;
    /* List of players */
    vector<Player> players;
    
    Team(const string& _team_name):team_name(_team_name){};
    
    vector<Team> split_into(int nb_teams){
        int nb_players = std::floor(players.size()/nb_teams);
        vector<Team> res;
        string tname = team_name;
        string city_name = tname.substr(0, tname.find_first_of("_"));
        size_t idx = 0;
        for(int i = 0; i<nb_teams; i++){
            Team t(city_name+to_string(i+1)+tname.substr(tname.find_first_of("_")));
            if(i==nb_teams-1)
                nb_players = players.size() - nb_players*(nb_teams-1);
            for(int n = 0; n < nb_players; n++){
                t.players.push_back(players[idx++]);
            }
            res.push_back(t);
        }
        assert(idx==players.size());
        return res;
    }
    
    Team combine_with(const Team& t) const{
        string tname  = t.team_name;
        Team res(team_name+tname.substr(tname.find_first_of("_")));
        for(int n = 0; n < players.size(); n++){
            res.players.push_back(players[n]);
        }
        for(int n = 0; n < t.players.size(); n++){
            res.players.push_back(t.players[n]);
        }        
        return res;
    }
    
    void save_team() const{
        string csv = "First Name,Last Name,Age,Email\n";
        for(int n = 0; n < players.size(); n++){
            csv += players[n].first_name+","+players[n].last_name+","+to_string(players[n].age)+","+players[n].email+"\n";
        }
        string path = team_name+".csv";
        std::ofstream outfile;
        outfile.open(path, std::ifstream::out | std::ifstream::binary);
        outfile << csv;
        outfile.close();
    }

};
#endif
