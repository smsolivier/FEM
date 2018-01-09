#include "EigenvalueSolver.H"
#include "FEGrid.H"
#include "Material.H"
#include "SparseMatrix.H"
#include "VectorMath.H"
#include <iostream>

void EigenvalueSolver::solve(vector<double>& a_phi, double& a_k, 
	DiffusionOperator& a_op, LinearSolver& a_solver, 
	const double a_tol, const int a_iter) {

	FEGrid grid = a_op.getFEGrid(); 
	int N = grid.getNumInteriorNodes(); 
	a_phi.resize(N, 0); 
	vector<double> phi_old(N,1); 
	double norm, knorm; 
	double k_old; 

	const Materials& mat = a_op.getMaterials(); 

	k_old = grid.integrateBi(phi_old); 

	// renormalize 
	phi_old /= k_old; 

	int i; 
	for (i=0; i<a_iter; i++) {

		// apply nuSigmaf and integrate 
		a_op.applyNuSigmaf(a_phi, phi_old);

		// solve (solution returned in a_phi) 
		a_solver.solve(a_phi); 

		// compute k by integrating a_phi
		a_k = grid.integrateBi(a_phi); 

		// renormalize 
		a_phi /= a_k; 

		// get L2 norm of difference 
		norm = L2_norm(a_phi - phi_old); 

		// get difference between k 
		knorm = abs(a_k - k_old);  

		if (m_verbose) {
			cout << i << " " << norm << " " << knorm << "\r"; 
			cout.flush(); 
		}

		if (knorm < a_tol && norm < a_tol) break; 

		phi_old = a_phi; 
		k_old = a_k; 

	}
	if (i == a_iter) {
		cout << "WARNING: maximum number of iterations reached in EigenvalueSolver" << endl; 
	}

	// normalize to 1 
	double max = 0; 
	for (int i=0; i<N; i++) {
		if (a_phi[i] > max) max = a_phi[i]; 
	}
	for (int i=0; i<N; i++) {
		a_phi[i] /= max; 
	}

	if (m_print) {
		cout << "Number of Iterations = " << i << endl; 
		cout << "Final Phi Norm = " << norm << endl; 
		cout << "Final k Norm = " << knorm << endl; 
	}

}