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

	// setup field object for P2P1 ordering of unknowns 
	Fields fields; 
	for (int i=0; i<grid.getNumElements(); i++) {
		Element& el = grid.getElement(i); 
		for (int j=0; j<el.getNumNodes(); j++) {
			if (el[j].isInterior()) {
				int jid = el[j].interiorNodeID(); 
				fields.set(jid, "u"); 
				fields.set(jid, "v"); 
				if (j < 3) {
					fields.set(jid, "p"); 
				}
			}
		}
	}

	Materials mat; 
	double Re = 500; 
	double L = .01; 
	mat("mat", "rho") = 10;
	mat("mat", "mu") = mat("mat", "rho")*1*L/Re; 
	cout << "nu = " << mat("mat", "mu")/mat("mat", "rho") << endl; 
	cout << "Re = " << Re << endl;
	
	NSOp ns(grid, mat, fields); 
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

	// LU linsol(A); 
	Cholesky linsol(A); 

	cout << "solving system" << endl; 
	linsol.solve(sol); 

	cout << "writing fields to vtk" << endl; 
	grid.writeFields(sol, fields, "solution.vtk"); 
}
