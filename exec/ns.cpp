#include "Picard.H"
#include "FEGrid.H"
#include "LU.H"
#include "Materials.H"
#include "LinearSolver.H"

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
	double Re = 500; 
	mat("mat", "rho") = 10;
	mat("mat", "mu") = mat("mat", "rho")*1*.1/Re; 

	Picard picard; 
	picard.printStats(); 
	picard.setVerbose(2); 
	LU lu; 

	vector<double> sol;  
	double tol = 1e-3; 
	int maxiter = 100; 
	cout << "starting picard iterations" << endl; 
	picard.solve(sol, grid, mat, lu, tol, maxiter); 

	grid.writeFields(sol, "solution.vtk"); 
}