#include "Jacobi.H"
#include "CH_Timer.H"
#include <cmath>
#include <iostream>

double max(const vector<double>& a_a) {

	int N = a_a.size(); 

	double max = 0.; 

	for (int i=0; i<N; i++) {

		if (abs(a_a[i]) > max) max = abs(a_a[i]); 

	}

	return max; 

}

Jacobi::Jacobi(const SparseMatrix& a_A, const double a_tol, 
	const int a_maxiter, const double a_alpha) {

	// store sparse matrix, tolerance and maxiter
	m_A = a_A; 
	m_N = m_A.getN(); 
	m_tol = a_tol; 
	m_maxiter = a_maxiter; 
	m_alpha = a_alpha; 
	m_print = false; 

}

void Jacobi::solve(vector<double> &a_rhs) {

	CH_TIMERS("Jacobi solve"); 

	// find largest diagonal 
	double max_diag = 0; 
	for (int i=0; i<a_rhs.size(); i++) {

		if (m_A(i,i) > max_diag) max_diag = m_A(i,i); 

	}

	// set relaxation parameter
	double relax = m_alpha/max_diag; 

	vector<double> phi_old(a_rhs.size(), 0); 
	vector<double> phi(a_rhs.size(), 0); 

	const double max_rhs = max(a_rhs); 

	double R; 

	int i; 
	for (i=0; i<m_maxiter; i++) {

		phi = phi_old + (a_rhs - (m_A*phi_old))*relax; 

		vector<double> resid = (m_A*phi) - a_rhs; 

		R = max(resid)/max_rhs; 

		if (R < m_tol) break; 

		phi_old = phi; 

	}

	if (i == m_maxiter) {
		cout << "WARNING: maximum number of iterations reached in Jacobi solver" << endl; 
	}

	if (m_print) {
		cout << "number of iterations = " << i << endl; 
		cout << "final residual = " << R << endl; 
	}

	// return solution in right hand side 
	a_rhs = phi; 

}