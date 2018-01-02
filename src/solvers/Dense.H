#ifndef _DENSE_H_
#define _DENSE_H_

#include "LinearSolver.H"
#include "SparseMatrix.H"
#include <vector>

using namespace std; 

/// Dense linear solver. Runs Gauss Elimination on a dense matrix 
class Dense : public LinearSolver {
public:
	/// constructor 
	/** 
		the create dense matrix function is called from SparseMatrix 
	*/ 
	Dense(const SparseMatrix& a_A); 
	/// Gauss Elimination to solve Ax = a_rhs 
	/** returns solution in a_rhs */ 
	virtual void solve(vector<double>& a_rhs); 
private:
	/// the gauss elimination function 
	int gauss_elim(int dim, vector<vector<double>>& a, vector<double>& x, 
		vector<double>& b); 

	int m_m; 
	int m_n; 
	vector<vector<double>> m_A; 
}; 

#endif 