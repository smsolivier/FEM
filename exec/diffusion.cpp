#include <vector>
#include <iostream>
#include <cmath>

#include "FEGrid.H"
#include "DiffusionOperator.H"
#include "CG.H"
#include "LU.H"
#include "Material.H"

using namespace std; 

double MMS(vector<double> a_position) {

	return sin(M_PI*a_position[0])*sin(M_PI*a_position[1]); 

}

double sourceFunction(vector<double> a_position) {

	// return 2./3*M_PI*M_PI*MMS(a_position) + .1*MMS(a_position); 
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

	// setup materials 
	Materials materials; 
	materials("core", "Sigma_t") = 1; 
	materials("core", "Sigma_a") = .1; 
	materials("reflector", "Sigma_t") = 10; 
	materials("reflector", "Sigma_a") = 0; 

	DiffusionOperator op(grid, materials); 

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

	vector<double> b; 

	CG * cg = new CG(A, 1e-6, 1000); 
	cg->printStats(); 
	cg->plot(); 
	op.makeRHS(b, FNodes); 
	cg->solve(b); 
	delete(cg); 

	// LU * lu = new LU(A); 
	// op.makeRHS(b, FCentroids); 
	// lu->solve(b); 
	// delete(lu); 

	grid.write(b, "solution.vtk"); 

}