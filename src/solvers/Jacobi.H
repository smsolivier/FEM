#ifndef _JACOBI_H_
#define _JACOBI_H_

#include "LinearSolver.H"
#include "SparseMatrix.H"
#include "VectorMath.H"
#include <vector>

using namespace std; 

/// Point Jacobi iteration class 
class Jacobi : public LinearSolver {
public:
	/// Constructor 
	/** 
		set the absolute tolerance, maximum number of iterations, and relaxation parameter
		Relaxation parameter defaults to .85 
	*/ 
	Jacobi(const SparseMatrix& a_A, const double a_tol, 
		const int a_maxiter, const double a_alpha=.85); 

	/// solve interface 
	virtual void solve(vector<double> &a_rhs); 

	/// call to have iteration information output after the solve 
	void printStats() {m_print = true; }; 
private:
	int m_N; 
	double m_tol; 
	int m_maxiter; 
	double m_alpha; 
	bool m_print; 
}; 

#endif