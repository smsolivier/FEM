#include "CG.H"
#include "VectorMath.H"
#include "CH_Timer.H"
#include <iostream>
#include <fstream>

CG::CG(const SparseMatrix& a_A, const double a_tol, const int a_maxiter) {

	m_A = a_A; 
	m_tol = a_tol; 
	m_maxiter = a_maxiter; 

	m_N = m_A.getN(); 

	m_print = false; 

}

void CG::solve(vector<double>& a_rhs) {

	CH_TIMERS("CG solve"); 

	vector<double> x(a_rhs.size(), 0);

	vector<double> r = a_rhs - m_A*x; 
	vector<double> s = r; 

	// temporary variables
	double denom, alpha, beta, norm; 
	vector<double> As; 

	vector<double> residual; 

	int i; 
	for (i=0; i<m_maxiter; i++) {

		As = m_A*s; 
		denom = dot(s, As); 
		if (denom == 0) {cout << "denom = 0" << endl; break; }
		alpha = dot(s, r/denom); 
		x = x + s*alpha; 
		r = a_rhs - m_A*x; 
		beta = -dot(r, As/denom); 
		s = r + s*beta; 
		norm = L2_norm(r);  
		if (m_plot) residual.push_back(norm); 
		if (norm < m_tol) break; 

	}

	if (i == m_maxiter) cout << "WARNING: maximum iterations reached in CG solver" << endl; 

	if (m_print) {
		cout << "number of iterations = " << i << endl; 
		cout << "final residual = " << norm << endl; 
	}

	if (m_plot) {
		ofstream out; 
		out.open("residual"); 
		for (int i=0; i<residual.size(); i++) {
			out << i << " " << residual[i] << endl; 
		}
		out.close(); 
	}

	// return solution in rhs 
	a_rhs = x; 

}