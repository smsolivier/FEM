#include <vector>
#include <iostream>
#include <cmath>

#include "FEGrid.H"
#include "PoissonOperator.H"
#include "CG.H"
#include "LU.H"
#include "GMRES.H"
#include "Materials.H"

using namespace std; 

double MMS(vector<double> a_position) {

	return sin(M_PI*a_position[0])*sin(M_PI*a_position[1]); 

}

double sourceFunction(vector<double> a_position) {

	// return 2*M_PI*M_PI*MMS(a_position); 
	return 1; 
	// return 0; 

}

int main(int argc, char* argv[]) {

	if (argc < 2) {
		cout << "give mesh name" << endl; 
		exit(0); 
	}

	string prefix(argv[1]); 
	FEGrid grid("../mesh/"+prefix); 
	grid.meshInfo(); 

	Materials mat; 
	mat("core", "k") = 1; 
	mat("reflector", "k") = 1; 
	mat("test", "k") = 1; 
	mat("test2", "k") = 1; 
	mat("test3", "k") = 1; 
	PoissonOperator op(grid, mat); 

	SparseMatrix A = op.matrix(); 
	A.sparsity(); 

	if (!A.symmetric()) cout << "A is not symmetric" << endl; 
	if (!A.positiveDiagonal()) cout << "A does not have positive diagonal" << endl; 

	int nEl = grid.getNumElements(); 
	int nNodes = grid.getNumInteriorNodes(); 

	vector<double> FNodes(nNodes); 
	for (int i=0; i<grid.getNumNodes(); i++) {
		Node node = grid.node(i); 
		int id = node.getInteriorNodeID(); 
		if (id >= 0) {
			FNodes[id] = sourceFunction(node.getPosition());
		}
	}

	vector<double> FCentroids(nEl, 0); 
	for (int i=0; i<nEl; i++) {
		Element& el = grid.getElement(i); 
		int group = el.getGroup(); 
		if (group > 0) {
			FCentroids[i] = 1; 
		}
	}

	vector<double> b; 

	CG * cg = new CG(A, 1e-6, 1000); 
	cg->printStats(); 
	cg->plot(); 
	op.makeRHS(b, FNodes); 
	// op.makeRHSAtCentroids(b, FCentroids); 
	cg->solve(b); 
	delete(cg); 

	// GMRES * gmres = new GMRES(A, 1e-6, 50, 1000); 
	// gmres->printStats(); 
	// op.makeRHS(b, FNodes); 
	// gmres->solve(b); 
	// delete(gmres); 

	// LU * lu = new LU(A); 
	// op.makeRHS(b, FCentroids); 
	// lu->solve(b); 
	// delete(lu); 

	grid.write(b, "solution.vtk"); 

}