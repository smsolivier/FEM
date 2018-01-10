#include "FEGrid.H"
#include "DiffusionOperator.H"
#include "EigenvalueSolver.H"
#include "LU.H"
#include "CG.H"
#include "Cholesky.H"
#include "CH_Timer.H"
#include <iostream>
using namespace std;

// runs power iteration on a command line provided mesh 

int main(int argc, char** argv) {

	if (argc != 2) {
		cout << "this program takes mesh arguments" << endl; 
		exit(0); 
	}

	string prefix(argv[1]);

	FEGrid grid("../mesh/"+prefix);
	grid.meshInfo(); 

	// setup materials 
	Materials mat; 
	mat("core", "Sigma_t") = 6; 
	mat("core", "Sigma_a") = 5; 
	mat("core", "nuSigma_f") = 2.43*2; 

	mat("gap", "Sigma_t") = .001; 
	mat("gap", "Sigma_a") = 0; 
	mat("gap", "nuSigma_f") = 0; 

	mat("clad", "Sigma_t") = 1; 
	mat("clad", "Sigma_a") = .1; 
	mat("clad", "nuSigma_f") = 0; 

	mat("reflector", "Sigma_t") = 6; 
	mat("reflector", "Sigma_a") = .1; 
	mat("reflector", "nuSigma_f") = 0;

	DiffusionOperator op(grid, mat);

	vector<double> phi; 
	double k; 
	EigenvalueSolver eig; 
	eig.printStats(); 
	// eig.setVerbose(); 
	double tol = 1e-5; 
	int maxiter = 1000; 
	LU lu(op.matrix()); 
	CG cg(op.matrix(), 1e-6, 10000); 
	Cholesky chol(op.matrix()); 
	eig.solve(phi, k, op, chol, tol, maxiter); 

	cout << "Numerical k = " << k << endl; 
	// double D = 1./(3*core("Sigma_t")); 
	// double kexact = core("nuSigma_f")/(core("Sigma_a") + 2*D*M_PI*M_PI); 
	// cout << "Exact k = " << kexact << endl; 

	grid.write(phi, "solution.vtk"); 
}