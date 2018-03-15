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
Constraint::Constraint():Constraint("noname") {
    _indices = make_shared<map<string,unsigned>>();
    _rev_indices = make_shared<vector<string>>();
    _ids = make_shared<vector<vector<unsigned>>>();
    _ids->resize(1);
};
Constraint::Constraint(string name):Constraint(name, leq) {
    _indices = make_shared<map<string,unsigned>>();
    _rev_indices = make_shared<vector<string>>();
    _ids = make_shared<vector<vector<unsigned>>>();
    _ids->resize(1);
};
Constraint::Constraint(string name, ConstraintType ctype) {
    _name = name;
    _ctype = ctype;
    _rhs = 0;
    _is_constraint = true;
    _all_lazy = make_shared<bool>(false);
    _indices = make_shared<map<string,unsigned>>();
    _rev_indices = make_shared<vector<string>>();
    _ids = make_shared<vector<vector<unsigned>>>();
    _ids->resize(1);
};

Constraint::Constraint(const Constraint& c):func_(c) {
    _name = c._name;
    _ctype = c._ctype;
    _rhs = c._rhs;
    _id = c._id;
    _ids = c._ids;
    _is_constraint = true;
    _all_lazy = make_shared<bool>(false);
    _indices = c._indices;
    _rev_indices = c._rev_indices;
    _ids->resize(1);
}
//@}
/* Destructor */
Constraint::~Constraint() {};

/* Accessors */
size_t Constraint::get_nb_instances() const {
    size_t nb = 0;
    if (!*_all_lazy) {
        return _nb_instances;
    }
    for (unsigned i = 0; i<_nb_instances; i++) {
        if (!_lazy[i]) {
            nb++;
        }
    }
    return nb;
}

double Constraint::get_rhs() const {
    return _rhs;
};


string Constraint::get_name() const {
    return _name;
};


int Constraint::get_type() const {
    return _ctype;
};

bool Constraint::is_active(unsigned inst) const {
    return fabs(get_val(inst) - _rhs) < EPS;
//    return fabs(_dual[inst]) >  EPS;
}

/* Operators */

Constraint& Constraint::operator=(const Constraint& c) {
    _ctype = c._ctype;
    _rhs = c._rhs;
    _id = c._id;
    _name = c._name;
    _dual = c._dual;
    _rev_indices = c._rev_indices;
    this->func_::operator=(c);
    return *this;
}

Constraint& Constraint::operator=(Constraint&& c) {
    _ctype = c._ctype;
    _rhs = c._rhs;
    _id = c._id;
    _name = c._name;
    _dual = c._dual;
    _rev_indices = move(c._rev_indices);
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

Constraint& Constraint::operator <=(const func_& f) {
    _ctype = leq;
    (*this) -= f;
    return *this;
};
Constraint& Constraint::operator >=(const func_& f) {
    _ctype = geq;
    (*this) -= f;
    return *this;
};

Constraint& Constraint::operator==(const func_& f) {
    _ctype = eq;
    (*this) -= f;
    return *this;
}

Constraint& Constraint::operator=(const func_& f) {
    this->func_::operator=(f);
    return *this;
}


size_t Constraint::get_id_inst(size_t ind) const {
    if (_nb_instances==1) {
        return 0;
    }
    return ind;
}


/* Output */

void Constraint::print_expanded() {
    auto nb_inst = get_nb_instances();
    for (unsigned inst = 0; inst<nb_inst; inst++) {
        eval(inst);
        print(inst);
    }
}

void Constraint::print(unsigned inst) {
    cout << _name;
    if (_nb_instances>1) {
        cout << "[" << to_string(inst) << "]";
    }
    if (is_convex()) {
        cout << " (Convex Constraint) : ";
    }
    else if (is_concave()) {
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

bool Constraint::is_convex() const {
    return ((_all_convexity==convex_ &&_ctype==leq) || (_all_convexity==concave_ &&_ctype==geq));
}

bool Constraint::is_concave() const {
    return ((_all_convexity==convex_ &&_ctype==geq) || (_all_convexity==concave_ &&_ctype==leq));
}

bool Constraint::is_ineq() const {
    return (_ctype==leq || _ctype==geq);
}

void Constraint::print() {
    cout << _name;
    if (is_convex()) {
        cout << " (Convex Constraint) : ";
    }
    else if (is_concave()) {
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
