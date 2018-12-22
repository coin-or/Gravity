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
    _all_lazy = make_shared<bool>(false);
    _dim[0] = 1;
};

Constraint::Constraint(const Constraint& c):func_(c){
    _name = c._name;
    _ctype = c._ctype;
    _rhs = c._rhs;
    _id = c._id;
    _is_constraint = true;
    _all_lazy = c._all_lazy;
}

//@}

/* Destructor */
Constraint::~Constraint(){};





/* Accessors */

size_t Constraint::get_nb_instances() const{
    size_t nb = 0;
    if (!*_all_lazy) {
        return _dim[0];
    }
    if(_lazy.size()==0){
        return 0;
    }
    for (size_t i = 0; i<_dim[0]; i++) {
        if (!_lazy[i]) {
            nb++;
        }
    }
    return nb;
}

double Constraint::get_rhs() const{
    return _rhs;
};


string Constraint::get_name() const{
    return _name;
};


int Constraint::get_type() const{
    return _ctype;
};

bool Constraint::is_active(size_t inst, double tol) const{
    return fabs(get_val(inst) - _rhs) < tol;
//    return fabs(_dual[inst]) >  EPS;
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
    if (_dim[0]==1) {
        return 0;
    }
    return ind;
}


/* Output */

void Constraint::print(){
    auto nb_inst = _dim[0];
    allocate_mem();
    for (size_t inst = 0; inst<nb_inst; inst++) {
        if (*_all_lazy && _lazy[inst]) {
            continue;
        }
        eval(inst);
        print(inst);
    }
}

void Constraint::print(size_t inst){
    cout << _name;
    if (_dim[0]>1) {
        if (!_indices || _indices->empty()) {
            cout << "[" << inst << "]";
        }
        else {
            cout << "[" << _indices->_keys->at(inst) << "]";
        }
    }
    if (is_complex()) {
        cout << " (Complex) : ";
    }
    else if (is_linear()) {
        cout << " (Linear) : ";
    }
    else if (is_convex() || is_rotated_soc() || is_soc()) {
        cout << " (Convex) : ";
    }
    else if (is_concave()){
        cout << " (Concave) : ";
    }
    else {
        cout << " (Unknown) : ";
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
    return (_all_convexity==linear_ || (_all_convexity==convex_ &&_ctype==leq) || (_all_convexity==concave_ &&_ctype==geq));
}

bool Constraint::is_concave() const{
    return (_all_convexity==linear_ || (_all_convexity==convex_ &&_ctype==geq) || (_all_convexity==concave_ &&_ctype==leq));
}

bool Constraint::is_ineq() const{
    return (_ctype==leq || _ctype==geq);
}

void Constraint::print_symbolic(){
    cout << _name;
    if (is_complex()) {
        cout << " (Complex) : ";
    }
    else if (is_linear()) {
        cout << " (Linear) : ";
    }
    else if (is_convex()) {
        cout << " (Convex) : ";
    }
    else if (is_concave()){
        cout << " (Concave) : ";
    }
    else {
        cout << " (Unknown) : ";
    }
    this->func_::print_symbolic(false,false);
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
