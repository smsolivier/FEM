#include "SparseMatrix.H"
#include "CG.H"
#include "LU.H"
#include "Jacobi.H"
#include "Dense.H"
#include "GMRES.H"
#include "VectorMath.H"
#include <string> 
#include <vector>
#include <cmath>
#include <iostream>

using namespace std; 

/*  tests all solvers for accuracy 
	solves Ax = b with b = A*ones(N) so that the solution will be all ones 
	if all entries of the resulting solution are within 1e-3 of 1 the solver passes 
*/ 

string testPass(LinearSolver& a_solver, vector<double> a_rhs) {

	a_solver.solve(a_rhs); 
	bool passed = true; 
	for (int i=0; i<a_rhs.size(); i++) {
		if (abs(a_rhs[i] - 1) > 1e-3) passed = false; 
	}

	if (passed) return "passed"; 
	else return "failed";  
}

int main() {

	// size of matrix 
	int N = 100000; 

	// setup 1D finite difference matrix (CG requires SPD) 
	SparseMatrix A(N,N); 
	A(0,0) = 2; 
	A(0,1) = -1; 
	A(N-1,N-1) = 2; 
	A(N-1,N-2) = -1; 
	for (int i=1; i<N-1; i++) {

		A(i,i) = 2; 
		A(i,i-1) = -1; 
		A(i,i+1) = -1; 

	}

	// force solution to ones 
	vector<double> x(N, 1); 
	vector<double> b = A*x; 

	// create solver objects 
	// CG cg(A, 1e-8, 10000); 
	// Jacobi jacobi(A, 1e-8, 100000); 
	LU lu(A); 
	// Dense dense(A); 
	// GMRES gmres(A, 1e-8, N/5, 100); 
	// gmres.setVerbose(); 
	// gmres.printStats(); 

	// make sure solution is ones 
	// cout << "CG " << testPass(cg, b) << endl; 
	// cout << "Jacobi " << testPass(jacobi, b) << endl; 
	cout << "LU " << testPass(lu, b) << endl; 
	// cout << "Dense " << testPass(dense, b) << endl; 
	// cout << "GMRES " << testPass(gmres, b) << endl; 

}