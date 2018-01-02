#ifndef _CHOLESKY_H_
#define _CHOLESKY_H_

#include "LinearSolver.H"
#include "SparseMatrix.H"
#include "Sparse"
#include <vector>

using namespace std; 

/// solve with Eigen's Cholesky factorization 
class Cholesky : public LinearSolver {
public:
	/// constructor 
	Cholesky(const SparseMatrix& a_A); 
	/// solve 
	void solve(vector<double>& a_rhs); 
private:
	Eigen::SparseMatrix<double> m_eigen; 
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> m_chol; 
	int m_m, m_n; 
}; 

#endif 