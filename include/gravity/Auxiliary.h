//
//  Auxiliary.h
//  Gravity
//
//  Created by Hassan Hijazi on 1/3/18.
//
//

#ifndef Auxiliary_h
#define Auxiliary_h

using namespace std;

namespace gravity {
    
    /** Backbone class for auxiliary objects that can be attached to nodes, e.g., generators. */
    class aux{
    public:
        bool _active;
        string _name;
        virtual ~aux(){};
    };
};
#endif /* Auxiliary_h */
