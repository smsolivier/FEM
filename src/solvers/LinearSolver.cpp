#include "LinearSolver.H"
#include <iostream>

void LinearSolver::operator()(const SparseMatrix& a_A) {
	cout << "operator() not defined in LinearSolver object" << endl; 
	exit(0); 
}