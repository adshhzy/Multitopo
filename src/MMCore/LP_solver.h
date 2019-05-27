
/*
 * [ A Linear Programming solver ]
 * Given a range, detect whether it is feasible,
 * if feasible, find a feasible solution for s0~s_k-1 which lies
 * as far away from the range boundary as possible (to avoid
 * generating ambigurious labeling later on).
 * Method: maximize the minimum margin to each side of the boudary.
 *
 * Variables: s0 ~ s_k-1, z
 * Objective: maximize the margin z
 * Constraints:
 *  si - sj >= (V2F[vij][j] - V2F[vij][i]) + z
 *  <=> si - sj - z >= V2F[vij][j] - V2F[vij][i] (at most K*(K-1))
 *      z >= 0                                   (1)
 *      0 <= si <= 1                             (K)
 */

#ifndef _LP_SOLVER_H_
#define _LP_SOLVER_H_

#include "gurobi_c++.h"
using namespace std;

class LPSolver {
 public:
  // Initialize the solver.
  LPSolver() : _env(GRBEnv()), _model(GRBModel(_env)), _upbound(1.0) {
    // Quiet output.
    _model.getEnv().set(GRB_IntParam_OutputFlag, 0);
  }
  ~LPSolver() {}
  inline void addVariables(int K);
  inline void addConstraint(int i, int j, float val);
  inline void addEqualConstraint(int i, int j, float val);
  inline void clearConstraints();
  inline bool solve();
  inline void getResults(vector<float>& res);
  inline bool solveBoundary();
  inline void setVarUpbound(float upbound);

 private:
  int _K;
  float _upbound;
  GRBEnv _env;
  GRBModel _model;
  GRBVar _z;
  vector<GRBVar> _s;
};

// Set a upper bounding on the variables (si is in [0, _upbound]).
inline void LPSolver::setVarUpbound(float upbound) { _upbound = upbound; }

// Initialize the solver, set the variables and objective.
inline void LPSolver::addVariables(int K) {
  _K = K;

  // Create variables
  _z = _model.addVar(0.0, 10.0, 0.0, GRB_CONTINUOUS);  // z
  _s.resize(_K);
  for (int i = 0; i < _K; ++i) {
    _s[i] = _model.addVar(0.0, _upbound, 0.0, GRB_CONTINUOUS);  // s_i
  }
  _model.update();

  // Set objective: maximize z
  _model.setObjective(0.0 + _z, GRB_MAXIMIZE);
}

// Add one constraint.
inline void LPSolver::addConstraint(int i, int j, float val) {
  // [1] use operator overload, the easiest way to write code.
  // _model.addConstr(_s[i] - _s[j] - _z >= val);
  // [2] get rid of operator overload, improve efficiency
  _model.addConstr(_s[i] - _s[j] - _z, GRB_GREATER_EQUAL, (double)val);
  // // [3] try to get rid of LinearExpresion operator overload by using
  // addTerms to improve efficiency even more
  // // wrong! does not support this operation. addTerms takes an array of terms
  // and an array of coeffs.
  // GRBLinExpr exp(_s[i]);
  // exp.addTerms(_s[j], -1.0, 1);
  // exp.addTerms(_z, -1.0, 1);
  // _model.addConstr(exp, GRB_GREATER_EQUAL, (double)val);
}

// Add one equal constraint.
inline void LPSolver::addEqualConstraint(int i, int j, float val) {
  _model.addConstr(_s[i] - _s[j] - _z, GRB_EQUAL, (double)val);
}

// Clear all the constraints.
inline void LPSolver::clearConstraints() {
  GRBConstr* constrs = _model.getConstrs();
  int nConstr = _model.get(GRB_IntAttr_NumConstrs);
  for (int i = nConstr - 1; i >= 0; --i) {
    _model.remove(constrs[i]);
  }
  // update the model each time before solve, to update constraints
  _model.update();
}

// Solve.
inline bool LPSolver::solve() {
  try {
    _model.reset();
    _model.optimize();
    // Return true if the optimization is solved optimally && margin is not
    // zero.
    return (_model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) &&
           (_z.get(GRB_DoubleAttr_X) != 0.0);
  } catch (...) {
    return false;
  }
}

// Solve.
inline bool LPSolver::solveBoundary() {
  try {
    _model.reset();
    _model.optimize();
    // Return true if the optimization is solved optimally.
    return (_model.get(GRB_IntAttr_Status) == GRB_OPTIMAL);
  } catch (...) {
    return false;
  }
}

// Get the results s0 ~ s_k-1.
inline void LPSolver::getResults(vector<float>& res) {
  res.resize(_K);
  for (int i = 0; i < _K; ++i) {
    res[i] = _s[i].get(GRB_DoubleAttr_X);
  }
}

#endif  // _LP_SOLVER_H_
