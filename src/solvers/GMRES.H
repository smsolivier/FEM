#ifndef _GMRES_H_
#define _GMRES_H_

#include "LinearSolver.H"
#include "SparseMatrix.H"
#include <vector> 

using namespace std; 

/// Solve a linear system with GMRES 
class GMRES : public LinearSolver {
public:
	/// Constructor 
	GMRES(const SparseMatrix& a_A, const double a_tol, 
		const int a_innerIter, const int a_outerIter); 
	/// solve \f$Ax=b\f$\n Solution is returned in a_rhs 
	void solve(vector<double>& a_rhs); 

	/// output iteration information 
	void setVerbose(); 
	/// print iteration information after solve
	void printStats(); 
private:
	/// tolerance 
	double m_tol; 
	/// maximum number of outer iterations 
	int m_outerIter;  
	/// maximum number of inner iterations 
	int m_innerIter; 
	/// number of rows/columns 
	int m_N; 
	/// number of non zero entries 
	int m_nnz; 
	/// array of matrix entries 
	vector<double> m_vals; 
	/// row of the entries
	vector<int> m_rows; 
	/// columns of the entries 
	vector<int> m_cols;
	/// control output of GMRES solver 
	bool m_verbose;  
	/// if true, prints stats 
	bool m_print; 
}; 

#endif 