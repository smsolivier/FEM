#ifndef _LINEAR_SOLVER_H_
#define _LINEAR_SOLVER_H_

#include "SparseMatrix.H" 
#include <vector>
using namespace std; 

/// base class for all linear system solvers
/** 
	Allows EigenvalueSolver to take a generalized linear solver object. 
	EigenvalueSolver only needs to be able to call solve on the solver object 
*/ 
class LinearSolver {
public: 
	/// common interface for solving Ax = b
	/** 
		solution is returned in a_rhs
	*/ 
	virtual void solve(vector<double>& a_rhs) = 0; 

	/// get the sparse matrix, must be provided in constructor of inherited classes 
	const SparseMatrix& getA() const {return m_A; }; 
protected:
	SparseMatrix m_A; /// store the sparse matrix 
}; 

#endif