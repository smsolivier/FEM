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
	double Re = 100; 
	mat("mat", "rho") = 10;
	mat("mat", "mu") = mat("mat", "rho")*1*.1/Re; 
	cout << "Re = " << Re << endl; 

	Picard picard; 
	picard.printStats(); 
	picard.setVerbose(2); 
	LU lu; 

	vector<double> sol;  
	double tol = 1e-2; 
	int maxiter = 100; 
	cout << "starting picard iterations" << endl; 
	picard.solve(sol, grid, mat, fields, lu, tol, maxiter); 

	grid.writeFields(sol, fields, "solution.vtk"); 
}