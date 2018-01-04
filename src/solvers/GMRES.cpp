#include "GMRES.H"
#include "CH_Timer.H"
#include "gmres.hpp"
#include <cassert>

GMRES::GMRES(const SparseMatrix& a_A, const double a_tol, 
	const int a_innerIter, const int a_outerIter) {

	m_tol = a_tol; 
	m_innerIter = a_innerIter; 
	m_outerIter = a_outerIter; 

	m_verbose = false; 
	m_print = false; 

	m_N = a_A.getN(); 
	m_nnz = a_A.getNNZ(); 

	// convert SpMat to GMRES format 
	a_A.convertToTriplet(m_vals, m_rows, m_cols); 
}

void GMRES::solve(vector<double>& a_rhs) {

	CH_TIMERS("GMRES solve"); 

	assert(a_rhs.size() == m_N); 

	vector<double> rhs = a_rhs; 
	for (int i=0; i<a_rhs.size(); i++) {a_rhs[i] = 0; }

	double atol = m_tol; 

	int iter = mgmres_st(m_N, m_nnz, &m_rows[0], &m_cols[0], &m_vals[0], &a_rhs[0], 
		&rhs[0], m_outerIter, m_innerIter, atol, m_tol, m_verbose); 

	if (m_print) {
		cout << "number of iterations = " << iter << endl; 
		cout << "final residual = " << atol << endl; 
	}
}

void GMRES::setVerbose() {m_verbose = true; }
void GMRES::printStats() {m_print = true; }