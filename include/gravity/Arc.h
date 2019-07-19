//
//  Arc.h
//  Cycle_Basis_PF
//
//  Created by Hassan on 18/06/2014.

//

#ifndef Cycle_Basis_PF_Arc_h
#define Cycle_Basis_PF_Arc_h
#include <gravity/Node.h>
//#include <gravity/types.h>
#include <gravity/Path.h>
#include <assert.h>
#include <string>
#include <vector>


class Arc{
public:
    int _id;
    std::string _name;
    std::string _type_name="Arc";
    Node* _src;
    Node* _dest;
    double _weight;
    double _len = 0;
    bool _is_transformer = false;
    bool _active = true;
    bool _expansion = false;
    bool _parallel = false;
    set<int> _phases;
    bool _imaginary = false;
    int _free = false;
    bool in_cycle = false;
    Path* horton_path = nullptr;
    double weight = 0;
    
    std::vector<Node*> _intersection; // intersection of node _src and node _dest
    std::vector<gravity::index_pair*> _intersection_clique; // useful for clique tree

    /* @brief Returns the neighbour of n if n is a node of the arc, null otherwise */
    Node* neighbour(Node* n);
    
//  bool in_cycle;
//  Path* horton_path;
    
    
 /* @brief Returns the neighbour of n if n is a node of the arc, null otherwise */
 //   Node* neighbour(Node* n);
    
    Arc();
    Arc(std::string name);
    virtual ~Arc();
    Arc(Node* s, Node* d);
    Arc(Node* s, Node* d, double weight);
    Arc* clone();
    
    bool has_phase(const string& ph) const{
        for(auto ph_i: _phases){
            if("ph"+to_string(ph_i)==ph){
                return true;
            }
        }
        return false;
    }
    
    void set_phases(string phases){
        _phases = gravity::get_phases(phases);
    }
    
    /* Connects the current arc to its source and destination, adding itself to the list of branches in these nodes */
    void connect();
    
    void print();
    
    std::vector<gravity::index_pair*> get_intersection_clique();
};




#endif
