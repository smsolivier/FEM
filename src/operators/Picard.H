#ifndef _PICARD_H_
#define _PICARD_H_

#include "FEGrid.H"
#include "Materials.H"
#include "LinearSolver.H"
#include <vector> 

using namespace std; 

/// uses picard (fixed point iteration) to solve the nonlinear Navier Stokes equations
class Picard {
public:
	/// default constructor 
	Picard() {}; 
	/// solve 
	void solve(vector<double>& a_sol, FEGrid& a_grid, const Materials& a_materials, 
		LinearSolver& a_solver, double a_tol, int a_maxiter); 

	/// output iteration information 
	void printStats(); 
	/// print norm at the end of each iteration 
	void setVerbose(); 
	/// set verbosity
	/** v = 0: no output \n
		v = 1: print norm at the end of each iteration \n
		v = 2: print norm and write solution to file 
	*/ 
	void setVerbose(int v); 
private:
	/// if true prints iterative information at end 
	bool m_print; 
	/// if true prints norm at the end of each iteration 
	int m_verbose; 
}; 

#endif