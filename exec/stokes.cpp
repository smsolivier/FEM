#include "NSOp.H"
#include "FEGrid.H"
#include "CG.H"
#include "LU.H"
#include "Cholesky.H"
#include "GMRES.H"
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
	NSOp ns(grid); 

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

	// find a node
	// Node node; 
	// for (int i=0; i<grid.getNumNodes(); i++) {
	// 	node = grid.node(i); 
	// 	if (node.isInterior()) {
	// 		node.printPosition(); 
	// 		break; 
	// 	}
	// }

	// int id = node.interiorNodeID(); 

	// Fields& fields = grid.getFields(); 
	// int pid = fields[id]["p"]; 
	// sol[pid] = 0; 
	// for (int i=0; i<A.getM(); i++) {
	// 	if (i == pid) A(pid, pid) = 1; 
	// 	else {
	// 		if (A.at(pid, i) != 0) A(pid, i) = 0; 
	// 	}
	// }
	// setup solvers 
	// CG linsol(A, 1e-6, 1000); 
	// linsol.printStats(); 

	// GMRES linsol(A, 1e-6, 50, 500); 
	// linsol.printStats(); 

	// LU linsol(A); 
	Cholesky linsol(A); 

	linsol.solve(sol); 

	grid.writeFields(sol, "solution.vtk"); 
}
