//
//  Auxiliary.h
//  Gravity
//
//
//

#ifndef Auxiliary_h
#define Auxiliary_h
//#include <gravity/utils.h>

using namespace std;

namespace gravity {

    set<int> get_phases(string phases);
    
    /** Backbone class for auxiliary objects that can be attached to nodes, e.g., generators. */
    class aux{
    public:
        bool _active;
        string _name;
        set<int> _phases;
        virtual ~aux(){};
        void set_phases(string phases){
            _phases = get_phases(phases);
        }
        bool has_phase(const string& ph) const{
            for(auto ph_i: _phases){
                if("ph"+to_string(ph_i)==ph){
                    return true;
                }
            }
            return false;
        }
    };
};
#endif /* Auxiliary_h */
