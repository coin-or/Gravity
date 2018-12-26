////
////  constraint.cpp
////  Gravity
////
////  Created by Hijazi, Hassan (Data61, Canberra City) on 6/5/17.
////
////
//#include <math.h>
//#include <gravity/constraint.h>
//
//
//using namespace gravity;
///** Constructor */
////@{
//Constraint_::Constraint_():Constraint_("noname"){};
//Constraint_::Constraint_(string name):Constraint_(name, leq){};
//Constraint_::Constraint_(string name, ConstraintType ctype){
//    _name = name;
//    _ctype = ctype;
//    _is_constraint = true;
//    _all_lazy = make_shared<bool>(false);
//    _dim[0] = 1;
//};
//
//Constraint_::Constraint_(const Constraint_& c):func_(c){
//    _name = c._name;
//    _ctype = c._ctype;
//    _id = c._id;
//    _is_constraint = true;
//    _all_lazy = c._all_lazy;
//}
//
////@}
//
///* Destructor */
//Constraint_::~Constraint_(){};
//
//
//
//
//
///* Accessors */
//
//size_t Constraint_::get_nb_instances() const{
//    size_t nb = 0;
//    if (!*_all_lazy) {
//        return _dim[0];
//    }
//    if(_lazy.size()==0){
//        return 0;
//    }
//    for (size_t i = 0; i<_dim[0]; i++) {
//        if (!_lazy[i]) {
//            nb++;
//        }
//    }
//    return nb;
//}
//
//
//
//string Constraint_::get_name() const{
//    return _name;
//};
//
//
//int Constraint_::get_type() const{
//    return _ctype;
//};
//
//
///* Operators */
//
//Constraint_& Constraint_::operator=(const Constraint_& c){
//    _ctype = c._ctype;
//    _id = c._id;
//    _name = c._name;
//    _dual = c._dual;
//    this->func_::operator=(c);
//    return *this;
//}
//
//Constraint_& Constraint_::operator=(Constraint_&& c){
//    _ctype = c._ctype;
//    _id = c._id;
//    _name = c._name;
//    _dual = c._dual;
//    this->func_::operator=(move(c));
//    return *this;
//}
//
//Constraint_& Constraint_::operator<=(double rhs) {
//    _ctype = leq;
//    return *this;
//}
//
//Constraint_& Constraint_::operator==(double rhs) {
//    _ctype = eq;
//    return *this;
//}
//
//Constraint_& Constraint_::operator>=(double rhs) {
//    _ctype = geq;
//    return *this;
//}
//
//Constraint_& Constraint_::operator <=(const func_& f){
//    _ctype = leq;
//    (*this) -= f;
//    return *this;
//};
//Constraint_& Constraint_::operator >=(const func_& f){
//    _ctype = geq;
//    (*this) -= f;
//    return *this;
//};
//
//Constraint_& Constraint_::operator==(const func_& f){
//    _ctype = eq;
//    (*this) -= f;    
//    return *this;
//}
//
//Constraint_& Constraint_::operator=(const func_& f){
//    this->func_::operator=(f);
//    return *this;
//}
//
//
//size_t Constraint_::get_id_inst(size_t ind) const{
//    if (_dim[0]==1) {
//        return 0;
//    }
//    return ind;
//}
//
//
///* Output */
//
//void Constraint_::print(size_t inst){
//    cout << _name;
//    if (_dim[0]>1) {
//        if (!_indices || _indices->empty()) {
//            cout << "[" << inst << "]";
//        }
//        else {
//            cout << "[" << _indices->_keys->at(inst) << "]";
//        }
//    }
//    if (is_complex()) {
//        cout << " (Complex) : ";
//    }
//    else if (is_linear()) {
//        cout << " (Linear) : ";
//    }
//    else if (is_convex() || is_rotated_soc() || is_soc()) {
//        cout << " (Convex) : ";
//    }
//    else if (is_concave()){
//        cout << " (Concave) : ";
//    }
//    else {
//        cout << " (Unknown) : ";
//    }
//    this->func_::print(inst);
//    switch (_ctype) {
//        case leq:
//            cout << " <= ";
//            break;
//        case geq:
//            cout << " >=  ";
//            break;
//        case eq:
//            cout << " = ";
//            break;
//        default:
//            break;
//    }
//    //    printf("%f;\n", _rhs);
//    cout << 0 << ";\n";
//};
//
//bool Constraint_::is_convex() const{
//    return (_all_convexity==linear_ || (_all_convexity==convex_ &&_ctype==leq) || (_all_convexity==concave_ &&_ctype==geq));
//}
//
//bool Constraint_::is_concave() const{
//    return (_all_convexity==linear_ || (_all_convexity==convex_ &&_ctype==geq) || (_all_convexity==concave_ &&_ctype==leq));
//}
//
//bool Constraint_::is_ineq() const{
//    return (_ctype==leq || _ctype==geq);
//}
//
//void Constraint_::print_symbolic(){
//    cout << _name;
//    if (is_complex()) {
//        cout << " (Complex) : ";
//    }
//    else if (is_linear()) {
//        cout << " (Linear) : ";
//    }
//    else if (is_convex()) {
//        cout << " (Convex) : ";
//    }
//    else if (is_concave()){
//        cout << " (Concave) : ";
//    }
//    else {
//        cout << " (Unknown) : ";
//    }
//    this->func_::print_symbolic(false,false);
//    switch (_ctype) {
//        case leq:
//            cout << " <= ";
//            break;
//        case geq:
//            cout << " >=  ";
//            break;
//        case eq:
//            cout << " = ";
//            break;
//        default:
//            break;
//    }
//    //    printf("%f;\n", _rhs);
//    cout << 0 << ";\n";
//};
