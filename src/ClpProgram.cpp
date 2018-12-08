//
//  ClpProgram.cpp
//
//
//  Created by Guanglei Wang on 2018/12/6.
//


#include <gravity/ClpProgram.h>

ClpProgram::ClpProgram(Model* m){
  _model = m;
  //fill derivatives and nonzeros
  m->fill_in_maps();
  m->compute_funcs();
}

/**
 *@Description: build the model
 * either build by adding rows or
 * by adding columns by not
 * simultanenouly.
 */

void ClpProgram::prepare_model(){
  fill_in_clp_vars();
  create_clp_constraints();
  set_clp_objective();
}

bool ClpProgram::solve(SolveType type){
  ClpSimplex solver(*_clp);
  if (type == SOLVE_PRIMAL){
    solver.primal();
    return true;
  }
  solver.dual();
  return true;
}

void ClpProgram::update_model(){
  _model->fill_in_maps();
  _model->reset_funcs();
  _model->compute_funcs();
  //fill_in_clp_vars();
  create_clp_constraints();
  set_clp_objective();
}

void ClpProgram::set_clp_objective(){
  int objt;
  double coeff;
  int sign = 1;
  if (!_model->_obj._new) {
    return;
  }
  _model->_obj._new = false;
  if (_model->_objt == minimize)
    objt = 1;
  else
    objt = -1;
  
  for (auto& it1: _model->_obj.get_lterms()) {
    if (!it1.second._sign) {
      sign = -1;
    }
    sign = sign*objt;
    if (it1.second._coef->_is_transposed) {
      auto dim = it1.second._p->_dim[0];
      auto idx = it1.second._p->get_id();
      for (int j = 0; j<dim; j++) {
        coeff = sign*t_eval(it1.second._coef, j);
        _clp->setObjectiveCoefficient(idx + it1.second._p->get_id_inst(j), coeff);
      }
    }
    else {
      coeff = sign*t_eval(it1.second._coef);
      _clp->setObjectiveCoefficient(it1.second._p->get_id() + it1.second._p->get_id_inst(), coeff);
    }
  }
}

ClpProgram::~ClpProgram(){
  delete _clp;
}


/**
 * @description: filling clp vars
 * as all vars should be continuous vars
 * all others are discarded.
 */
void ClpProgram::fill_in_clp_vars(){
  _clp->resize(0, _model->get_nb_vars());
  param_* v;
  for(auto& v_p: _model->_vars)
  {
    v = v_p.second;
    if (!v->_new) {
      continue;//Variable already added to the program
    }
    v->_new = false;
    auto idx = v->get_id();
    if( idx == -1) {
      throw invalid_argument("Variable needs to be added to model first: use add_var(v) function:" + v->get_name());
    }
    switch (v->get_intype()) {
        case float_: {
          auto real_var = (var<float>*)v;
          for (int i = 0; i < real_var->get_dim(); i++) {
            auto vid = idx + v->get_id_inst(i);
            _clp->setColLower(vid, real_var->get_lb(i));
            _clp->setColUpper(vid, real_var->get_ub(i));
          }
          break;
        }
        case double_: {
          auto real_var = (var<double>*)v;
          for (int i = 0; i < real_var->get_dim(); i++) {
            auto vid = idx + v->get_id_inst(i);
            _clp->setColLower(vid, real_var->get_lb(i));
            _clp->setColUpper(vid, real_var->get_ub(i));
          }
          break;
        }
      default:
        break;
    }
  }
}


void ClpProgram::create_clp_constraints(){
  CoinBuild _build;
  size_t idx = 0, nb_inst = 0, inst = 0;
  Constraint* c;
  bool sign = 1;
  for(auto& p: _model->_cons) {
    c = p.second.get();
    if (!c->_new) {
      continue;//Constraint already added to the program
    }
    c->_new = false;
    if (c->is_nonlinear()) {
      throw invalid_argument("Clp cannot handle nonlinear constraints.\n");
    }
    if (c->is_quadratic()){
      throw invalid_argument("Clp cannot handle nonlinear constraints even it is convex quadratic. \n");
    }
    nb_inst = c->_dim[0];
    inst = 0;
    for (size_t i = 0; i< nb_inst; i++) {
      vector<int> row2Index;
      vector<double> row2Coeff;
      for (auto& it_lterm: c->get_lterms()) {
        idx = it_lterm.second._p->get_vec_id();
        if (!it_lterm.second._sign) {
          sign = -1;
        }
        if (it_lterm.second._coef->_is_transposed || it_lterm.second._coef->is_matrix()) {
          auto dim = it_lterm.second._p->get_dim(i); // dimension of this term
          for (int j = 0; j< dim; j++) {
            row2Coeff.push_back(sign*t_eval(it_lterm.second._coef,i,j));
            row2Index.push_back(it_lterm.second._p->get_id() +it_lterm.second._p->get_id_inst(i,j));
          }
        }
        else {
          row2Coeff.push_back(sign*t_eval(it_lterm.second._coef,i));
          row2Index.push_back(it_lterm.second._p->get_id() + it_lterm.second._p->get_id_inst(i));
        }
      }
      int row2Ind[row2Index.size()];
      double row2Val[row2Index.size()];
      for (size_t l=0;l < row2Index.size(); l++){
        row2Ind[l] = row2Index.at(l);
        row2Val[l] = row2Coeff.at(l);
      }
      if(c->get_type()==geq) {
        _build.addRow(row2Index.size(), row2Ind, row2Val, c->get_rhs(),numeric_limits<double>::infinity());
      }
      else if(c->get_type()==leq) {
        _build.addRow(row2Index.size(), row2Ind, row2Val, -numeric_limits<double>::infinity(), c->get_rhs());

      }
      else if(c->get_type()==eq) {
        _build.addRow(row2Index.size(), row2Ind, row2Val, c->get_rhs(), c->get_rhs());

      }
      inst++;
    }
  }
  _clp->addRows(_build);
}
