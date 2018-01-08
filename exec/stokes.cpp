#include "NSOp.H"
#include "FEGrid.H"
#include "LU.H"
#include "Cholesky.H"
#include "Materials.H"
#include <iostream>

using namespace std; 

int main(int argc, char* argv[]) {

	if (argc < 2) {
		cout << "give mesh name" << endl; 
		exit(0); 
	}

	string prefix(argv[1]); 
	FEGrid grid("../mesh/"+prefix); 
	grid.meshInfo(); 
	grid.precomputeIntegrals(); 

	Materials mat; 
	mat("mat", "mu") = 1; 
	mat("mat", "rho") = 1; 
	NSOp ns(grid, mat); 
	ns.buildLHS(); 

	SparseMatrix& A = ns.matrix(); 
	A.sparsity(); 

	if (!A.symmetric()) {
		cout << "A is not symmetric" << endl; 
	} 

	if (!A.positiveDiagonal()) {
		cout << "Diagonal < 0" << endl; 
	}

	vector<double> sol; 
	ns.makeRHS(sol); 

	LU linsol(A); 
	// Cholesky linsol(A); 

	cout << "solving system" << endl; 
	linsol.solve(sol); 

	cout << "writing fields to vtk" << endl; 
	grid.writeFields(sol, "solution.vtk"); 
}
