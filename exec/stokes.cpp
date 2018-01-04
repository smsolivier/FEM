#include "NSOp.H"
#include "FEGrid.H"
#include "CG.H"
#include "LU.H"
#include "Cholesky.H"
#include "GMRES.H"
#include <iostream>

using namespace std; 

int main() {

	FEGrid grid("../mesh/cavity"); 
	grid.meshInfo(); 
	grid.precomputeIntegrals(); 
	NSOp ns(grid); 

	SparseMatrix& A = ns.matrix(); 
	A.sparsity(); 

	if (!A.symmetric()) {
		cout << "A is not symmetric" << endl; 
	} 

	if (!A.positiveDiagonal()) {
		cout << "Diagonal < 0" << endl; 
	}

	// setup solvers 
	// CG linsol(A, 1e-6, 1000); 
	// linsol.printStats(); 

	// GMRES linsol(A, 1e-6, 50, 500); 
	// linsol.printStats(); 

	// LU linsol(A); 
	Cholesky linsol(A); 

	vector<double> sol; 
	ns.makeRHS(sol); 

	linsol.solve(sol); 

	grid.writeFields(sol, "solution.vtk"); 
}