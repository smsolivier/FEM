#include "Cholesky.H"
#include "CH_Timer.H"
#include <iostream>

Cholesky::Cholesky(const SparseMatrix& a_A) {

	CH_TIMERS("Cholesky Factorization"); 

	m_m = a_A.getM(); 
	m_n = a_A.getN(); 

	if (!a_A.symmetric()) {
		cout << "ERROR: matrix is not symmetric in Cholesky solver" << endl; 
		exit(0); 
	}

	m_eigen = a_A.getEigenMatrix(); 
	m_chol.compute(m_eigen); 

}

void Cholesky::solve(vector<double>& a_rhs) {

	CH_TIMERS("Cholesky Solve"); 

	// convert to Eigen format
	Eigen::VectorXd rhs = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(a_rhs.data(), a_rhs.size());

	Eigen::VectorXd sol = m_chol.solve(rhs); 

	for (int i=0; i<m_m; i++) {
		a_rhs[i] = sol[i]; 
	}

}