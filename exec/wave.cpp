#include <vector> 
#include <iostream>
#include <cmath> 

#include "WaveOperator.H"
#include "FEGrid.H"
#include "LU.H"
#include "CG.H"
#include "Cholesky.H"
#include "Materials.H"
#include "VectorMath.H"

using namespace std; 

double sourceFunction(vector<double> a_pos, double t) {
	double cx = 1; 
	double cy = 3; 
	double width = .3; 
	double depth = .5; 
	double rad = sqrt(pow(a_pos[0]-cx,2) + pow(a_pos[1]-cy,2)); 
	if (rad <= width) {
		return -cos(rad/width/2*M_PI)*depth;  
		// return -pow(M_E, 30*abs(width - rad))/pow(M_E, 30*width); 
	}
	else return 0; 
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
	mat("mat", "a") = .5; 
	mat("mat", "c") = 2; 

	int Nt = 3000; 
	double tend = 3; 
	double h = tend/Nt; 
	cout << "dt = " << h << endl; 

	WaveOperator op(grid, mat, h); 
	SparseMatrix& A = op.matrix(); 
	// LU linsol(A); 
	// CG linsol(A, 1e-6, 1000); 
	Cholesky linsol(A); 

	vector<double> u2(nIntNodes, 0); 
	vector<double> u1(nIntNodes, 0); 
	vector<double> u0(nIntNodes, 0); 

	double subdiv = 1; 
	int n = 1; 

	// initialize u0 
	for (int i=0; i<grid.getNumNodes(); i++) {
		Node node = grid.node(i); 
		int id = node.getInteriorNodeID(); 
		if (id >= 0) {
			vector<double> pos = node.getPosition(); 
			u0[id] = sourceFunction(pos, 0); 
		}
	}

	// initialize u1 
	for (int i=0; i<grid.getNumNodes(); i++) {
		Node node = grid.node(i); 
		int id = node.getInteriorNodeID(); 
		if (id >= 0) {
			vector<double> pos = node.getPosition(); 
			u1[id] = u0[id]; 
		}
	}

	grid.write(u0, "solution0.vtk"); 

	for (int i=1; i<Nt+1; i++) {

		cout << (double)i/Nt << "\r"; 
		cout.flush(); 
		double t = i*h; 

		vector<double> F(nIntNodes, 0); 
		// for (int j=0; j<grid.getNumNodes(); j++) {
		// 	Node node = grid.node(j); 
		// 	int id = node.getInteriorNodeID(); 
		// 	if (id >= 0) {
		// 		vector<double> pos = node.getPosition(); 
		// 		F[id] = pow(M_E, -10*t)*sourceFunction(pos, t); 
		// 	}
		// }

		op.makeRHS(u2, u1, u0, F, h); 
		linsol.solve(u2); 
		string fname = "solution"+to_string(i)+".vtk"; 
		grid.write(u2, fname.c_str()); 
		u0 = u1; 
		u1 = u2; 
	}

}