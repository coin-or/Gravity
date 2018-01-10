//
//  constraint.cpp
//  Gravity
//
//  Created by Hijazi, Hassan (Data61, Canberra City) on 6/5/17.
//
//
#include <math.h>
#include <gravity/constraint.h>


using namespace gravity;
/** Constructor */
//@{
Constraint::Constraint():Constraint("noname"){};
Constraint::Constraint(string name):Constraint(name, leq){};
Constraint::Constraint(string name, ConstraintType ctype){
    _name = name;
    _ctype = ctype;
    _rhs = 0;
    _is_constraint = true;
};

Constraint::Constraint(const Constraint& c):func_(c){
    _name = c._name;
    _ctype = c._ctype;
    _rhs = c._rhs;
    _id = c._id;
    _is_constraint = true;
}

//@}

/* Destructor */
Constraint::~Constraint(){};





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

bool Constraint::is_active(unsigned inst) const{
    return fabs(_dual[inst]) >  EPS;
}

/* Operators */

Constraint& Constraint::operator=(const Constraint& c){
    _ctype = c._ctype;
    _rhs = c._rhs;
    _id = c._id;
    _name = c._name;
    _dual = c._dual;
    this->func_::operator=(c);
    return *this;
}

Constraint& Constraint::operator=(Constraint&& c){
    _ctype = c._ctype;
    _rhs = c._rhs;
    _id = c._id;
    _name = c._name;
    _dual = c._dual;
    this->func_::operator=(move(c));
    return *this;
}

Constraint& Constraint::operator<=(double rhs) {
    _ctype = leq;
    _rhs = rhs;
    return *this;
}

Constraint& Constraint::operator==(double rhs) {
    _ctype = eq;
    _rhs = rhs;
    return *this;
}

Constraint& Constraint::operator>=(double rhs) {
    _ctype = geq;
    _rhs = rhs;
    return *this;
}

Constraint& Constraint::operator <=(const func_& f){
    _ctype = leq;
    (*this) -= f;
    return *this;
};
Constraint& Constraint::operator >=(const func_& f){
    _ctype = geq;
    (*this) -= f;
    return *this;
};

Constraint& Constraint::operator==(const func_& f){
    _ctype = eq;
    (*this) -= f;    
    return *this;
}

Constraint& Constraint::operator=(const func_& f){
    this->func_::operator=(f);
    return *this;
}


size_t Constraint::get_id_inst(size_t ind) const{
    if (_nb_instances==1) {
        return 0;
    }
    return ind;
}


/* Output */

void Constraint::print_expanded(){
    eval();
    auto nb_inst = get_nb_instances();
    for (unsigned inst = 0; inst<nb_inst; inst++) {
        print(inst);
    }
}

void Constraint::print(unsigned inst){
    cout << _name;
    if (_nb_instances>1) {
        cout << "[" << to_string(inst) << "]";
    }
    if (is_convex()) {
        cout << " (Convex Constraint) : ";
    }
    else if (is_concave()){
        cout << " (Concave Constraint) : ";
    }
    else {
        cout << " : ";
    }        
    this->func_::print(inst);
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

bool Constraint::is_convex() const{
    return ((_all_convexity==convex_ &&_ctype==leq) || (_all_convexity==concave_ &&_ctype==geq));
}

bool Constraint::is_concave() const{
    return ((_all_convexity==convex_ &&_ctype==geq) || (_all_convexity==concave_ &&_ctype==leq));
}

void Constraint::print(){
    cout << _name;
    if (is_convex()) {
        cout << " (Convex Constraint) : ";
    }
    else if (is_concave()){
        cout << " (Concave Constraint) : ";
    }
    else {
        cout << " : ";
    }
    this->func_::print(false,false);
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
