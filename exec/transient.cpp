#include <vector>
#include <iostream>
#include <cmath> 

#include "FEGrid.H"
#include "TransientPoisson.H"
#include "CG.H"
#include "LU.H"
#include "Cholesky.H"
#include "Materials.H"
#include "VectorMath.H"

double sourceFunction(vector<double> pos, double t) {
	// return sin(M_PI*pos[0]*t)*sin(M_PI*pos[1]*t); 
	return 1; 
}

int main(int argc, char* argv[]) {

	if (argc < 2) {
		cout << "give mesh name" << endl; 
		exit(0); 
	}

	string prefix(argv[1]); 
	FEGrid grid("../mesh/"+prefix); 
	grid.meshInfo(); 

	int nIntNodes = grid.getNumInteriorNodes(); 
	int nEl = grid.getNumElements(); 

	Materials mat; 
	mat("mat", "k") = 1; 
	mat("mat", "alpha") = 1/mat("mat", "k"); 

	vector<double> F(nIntNodes); 
	vector<double> T_old(nIntNodes, 0); 
	vector<double> T(nIntNodes); 

	int Nt = 200; 
	double tend = .2;
	double h = tend/Nt; 
	cout << "h = " << h << endl; 

	TransientPoisson op(grid, mat, h); 

	SparseMatrix& A = op.matrix(); 

	// LU linsol(A); 
	Cholesky linsol(A); 

	for (int i=1; i<Nt+1; i++) {

		double t= i*h; 
		for (int j=0; j<grid.getNumNodes(); j++) {
			Node node = grid.node(j); 
			int id = node.getInteriorNodeID(); 
			if (id >= 0) {
				F[id] = sourceFunction(node.getPosition(), t/tend); 
			}
		}

		op.makeRHS(T, F, T_old, h); 
		linsol.solve(T); 
		string fname = "solution"+to_string(i)+".vtk"; 
		grid.write(T, fname.c_str()); 
		T_old = T; 
	} 

}