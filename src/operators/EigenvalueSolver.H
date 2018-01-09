#ifndef _EIGENVALUE_SOLVER_H_
#define _EIGENVALUE_SOLVER_H_

#include <vector>
#include "DiffusionOperator.H"
#include "LinearSolver.H" 

using namespace std; 

/// Computes the dominant eigenvalue and corresponding eigenvector with Power Iteration 
class EigenvalueSolver 
{
public:
	/// solve for eigenvalue/vector using the solver provided by LinearSolver object 
	void solve(
		vector<double>& a_phi, // dominant eigenvector
		double& a_k, // dominant eigenvalue 
		DiffusionOperator& a_op, 
		LinearSolver& a_solver,
		const double a_tolerance, 
		const int a_iter); 

	/// outputs residuals and iteration information 
	void printStats() { m_print = true; }; 
	/// output norms as it goes 
	void setVerbose() {m_verbose = true; }; 

private:
	bool m_print; /// stores whether to print iterative information 
	bool m_verbose; 
}; 

#endif 