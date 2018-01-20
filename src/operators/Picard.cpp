#include "Picard.H"
#include "NSOp.H"
#include "VectorMath.H"
#include "CH_Timer.H"

void Picard::solve(vector<double>& a_sol, FEGrid& a_grid, 
	const Materials& a_materials, const Fields& a_fields, LinearSolver& a_solver, 
	double a_tol, int a_maxiter) {

	CH_TIMERS("picard iterations"); 

	// initialize solution 
	a_sol.resize(a_fields.size()); 
	for (int i=0; i<a_sol.size(); i++) {a_sol[i] = 0; }

	vector<double> psol = a_sol; 

	int i; 
	double norm; 
	for (i=0; i<a_maxiter; i++) {
		NSOp * ns = new NSOp(a_grid, a_materials, a_fields);
		ns->buildLHS(psol); 
		SparseMatrix& A = ns->matrix(); 
		ns->makeRHS(a_sol); 
		a_solver(A); 
		a_solver.solve(a_sol); 

		// convergence 
		norm = L2_norm(a_sol - psol); 
		if (norm < a_tol) break; 

		if (m_verbose > 0) {
			cout << i << " " << norm << endl; 
		}

		if (m_verbose > 1) {
			string fname = "solution"+to_string(i)+".vtk"; 
			a_grid.writeFields(a_sol, a_fields, fname.c_str()); 
		}

		psol = a_sol; 
		delete(ns); 
	}

	if (i == a_maxiter) {
		cout << "WARNING: maximum number of iterations reached in Picard solver" << endl; 
	}

	if (m_print) {
		cout << "final norm = " << norm << endl; 
		cout << "number of iterations = " << i << endl; 
	}
}

void Picard::printStats() {m_print = true; }
void Picard::setVerbose() {m_verbose = 1; }
void Picard::setVerbose(int v) {m_verbose = v; }