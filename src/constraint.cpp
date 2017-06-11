//
//  constraint.cpp
//  Gravity
//
//  Created by Hijazi, Hassan (Data61, Canberra City) on 6/5/17.
//
//

#include <Gravity/constraint.h>

/** Constructor */
//@{
Constraint::Constraint():Constraint("noname"){};
Constraint::Constraint(string name):Constraint(name, leq){};
Constraint::Constraint(string name, ConstraintType ctype){
    _name = name;
    _ctype = ctype;
    _rhs = 0;
};

Constraint::Constraint(const Constraint& c):func_(c){
    _name = c._name;
    _ctype = c._ctype;
    _rhs = c._rhs;
    _id = c._id;    
}

//@}

/* Destructor */
Constraint::~Constraint(){};



Constraint& Constraint::operator<=(double rhs) {
    _ctype = leq;
    _rhs = rhs;
    return *this;
}

Constraint& Constraint::operator=(double rhs) {
    _ctype = eq;
    _rhs = rhs;
    return *this;
}

Constraint& Constraint::operator>=(double rhs) {
    _ctype = geq;
    _rhs = rhs;
    return *this;
}


/* Accessors */

double Constraint::get_rhs() const{
    return _rhs;
};


string Constraint::get_name() const{
    return _name;
};


int Constraint::get_type() const{
    return _ctype;
};


/* Modifiers */

Constraint& Constraint::operator=(const func_& f){
    this->func_::operator=(f);
    return *this;
}

/* Output */
void Constraint::print() const{
    cout << _name << " : ";
    
    this->func_::print(false);
    switch (_ctype) {
        case leq:
            cout << " <= ";
            break;
        case geq:
            cout << " >=  ";
            break;
        case eq:
            cout << " = ";
            break;
        default:
            break;
    }
    //    printf("%f;\n", _rhs);
    cout << _rhs << ";\n";    
};
