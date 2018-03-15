//
//  constraint.hpp
//  Gravity
//
//  Created by Hijazi, Hassan (Data61, Canberra City) on 6/5/17.
//
//

#ifndef constraint_hpp
#define constraint_hpp

#include <stdio.h>
#include <gravity/func.h>

namespace gravity {
class Constraint :public func_ {

protected:
    string                      _name = "no_name";
    shared_ptr<map<string,unsigned>>       _indices = nullptr; /*<< A map storing all the indices this parameter has, the key is represented by a string, while the entry indicates the right position in the values and bounds vectors */
    shared_ptr<vector<string>>  _rev_indices = nullptr; /*<< A vector storing all the indices this parameter has */

public:
    unsigned                    _jac_cstr_idx; /* First index of the corresponding non-zero values in the Jacobian */
    unsigned                    _id = -1;
    ConstraintType              _ctype = leq; /**< Constraint type: leq, geq or eq */
    double                      _rhs = 0;
    vector<double>              _dual ; /**< Lagrange multipliers at a KKT point */
    bool                        _all_active = true;
    vector<bool>                _active;
    shared_ptr<bool>            _all_lazy;
    vector<bool>                _lazy;
    bool                        _all_satisfied = true;
    vector<bool>                _violated;

    /** Constructor */
    //@{
    Constraint();
    Constraint(const Constraint& c);
    Constraint(std::string name);
    Constraint(std::string name, ConstraintType ctype);
    //@}
    /* Destructor */
    ~Constraint();
    /* Boolean Requests */

    /* Operators */
    Constraint& operator=(const Constraint& c);
    Constraint& operator=(Constraint&& c);

    Constraint& operator <=(double rhs);
    Constraint& operator >=(double rhs);
    Constraint& operator ==(double rhs);
    Constraint& operator <=(const func_& rhs);
    Constraint& operator >=(const func_& rhs);
    Constraint& operator ==(const func_& rhs);
    Constraint& operator =(const func_& rhs);

    /* Accessors */
    size_t get_nb_instances() const;
    string get_name() const;
    int get_type() const;
    double get_rhs() const;
    bool is_active(unsigned inst = 0) const;
    bool is_convex() const;
    bool is_concave() const;
    bool is_ineq() const;

    size_t get_id_inst(size_t ind) const;


    /* Modifiers */

    void make_lazy() {
        *_all_lazy = true;
        _lazy.resize(_nb_instances,true);
    }

    Constraint& in(const node_pairs& np) {
        this->func_::in(np);
        return *this;
    };

    template<typename Tobj> Constraint& in(const vector<Tobj*>& vec) {
        this->func_::in(vec);
        //_rev_indices->resize(_val->size());
        //unsigned index = 0;
        //for(auto it = vec.begin(); it!= vec.end(); it++) {
        //    _rev_indices->at(index) = (*it)->_name;
        //    index++;
        //}
        return *this;
    };

    template<typename Tobj> Constraint& in(const vector<Tobj>& vec) {
        this->func_::in(vec);
        //_rev_indices->resize(_val->size());
        //unsigned index = 0;
        //for (auto it = vec.begin(); it !=vec.end(); it++){
        //    _rev_indices->at(index) = (*it)._name;
        //    index++;
        //}
        return *this;
    };

    template<typename Tobj> Constraint& in_at(const vector<Tobj>& vec, const unsigned t) {
        this->func_::in_at(vec, t);
        return *this;
    };

    template<typename Tobj> Constraint& in(const vector<Tobj*>& vec, const unsigned T) {
        this->func_::in(vec, T);
        string key;
        for (unsigned t = 0; t < T; t++) {
            for (auto it = vec.begin(); it != vec.end(); it++) {
                if(!(*it)->_active) {
                    continue;
                }
                key = (*it)->_name + "," + to_string(t);
                auto index = this->_indices->size();
                auto pp = this->_indices->insert(make_pair<>(key, index));
                if(pp.second) { //new index inserted
                    this->_val->resize(max(_val->size(), index+1));
                    this->_dim[0] = max(_dim[0],_val->size());
                    this->_rev_indices->resize(_val->size());
                    this->_rev_indices->at(index) = key;
                    this->func_::_ids->at(0).push_back(index);
                }
                else {
                    this->func_::_ids->at(0).push_back(pp.first->second);
                }
            }
        }
        return *this;
    };

    template<typename Tobj> Constraint& in_at(const vector<Tobj*>& vec, const unsigned T) {
        this->func_::in_at(vec, T);
        return *this;
    };

    template<typename... Args>
    Constraint operator()(string t1, Args&&... args) {
        Constraint res(this->_name);
        res._id = this->_id;
        res._range = this->_range;
        res._val = this->_val;
        res._is_vector = this->_is_vector;
        res._is_matrix = this->_is_matrix;
        res._is_transposed = _is_transposed;
        res._rev_indices = this->_rev_indices;
        res._indices = this->_indices;
        res._dual = this->_dual;
        list<string> indices;
        //indices = {forward<size_t>(args)...};
        indices = {forward<Args>(args)...};
        indices.push_front(t1);
        string key;
        auto it = indices.begin();
        for (int i= 0; i < indices.size(); i++) {
            key += *it;
            if (i< indices.size()-1) {
                key += ",";
            }
            it++;
        }
        if (indices.size()==2) {
            _is_matrix = true;
        }
        auto index = _indices->size();
        auto pp = _indices->insert(make_pair<>(key,index));
        if(pp.second) { //new index inserted
            _val->resize(max(_val->size(),index+1));
            _dim[0] = max(_dim[0],_val->size());
            _rev_indices->resize(_val->size());
            _rev_indices->at(index) = key;
            res._ids->at(0).push_back(_indices->size()-1);
        }
        else {
            res._ids->at(0).push_back(pp.first->second);
        }
        res._dim[0]=1;
        res._name += "["+key+"]";
        return res;
    };

    /* Output */
    void print_expanded();
    void print(unsigned);
    void print();


};
}
#endif /* constraint_hpp */
