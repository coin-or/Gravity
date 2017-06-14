//
//  constraint.cpp
//  Gravity
//
//  Created by Hijazi, Hassan (Data61, Canberra City) on 6/5/17.
//
//
#include <math.h>
#include <Gravity/constraint.h>
#define EPS 0.000001

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

bool Constraint::is_active() const{
    return fabs(_dual) >  EPS;
}

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
