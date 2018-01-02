#ifndef _CG_H_
#define _CG_H_

#include "LinearSolver.H"
#include "SparseMatrix.H"
#include <vector>
using namespace std; 

/// Conjugate gradient linear solver 
class CG : public LinearSolver {
public:
	/// Constructor. Takes an absolute tolerance and a maximum number of iterations 
	CG(const SparseMatrix& A, const double a_tol, const int a_maxiter); 

	/// solve interface. solution is returned in a_rhs 
	virtual void solve(vector<double>& a_rhs); 

	/// prints number of iterations and final residual 
	void printStats() {m_print = true; }; 
	void plot() {m_plot = true; }; 
private:
	double m_tol; 
	int m_maxiter; 
	int m_N; 
	bool m_print; 
	bool m_plot; 
}; 

#endif