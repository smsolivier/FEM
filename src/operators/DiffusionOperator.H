#ifndef _DIFFUSIONOPERATOR_H_
#define _DIFFUSIONOPERATOR_H_

#include "Operator.H"
#include "FEGrid.H"
#include "Materials.H"
#include <vector>

using namespace std; 

/// builds LHS and RHS for neutron diffusion 
class DiffusionOperator : public Operator {
public:
	/// builds LHS 
	DiffusionOperator(FEGrid& a_grid, const Materials& a_materials); 
	/// builds RHS 
	void makeRHS(vector<double>& a_rhs, const vector<double>& a_FNodes); 
	/// build the RHS and multiply by \f$\nu\Sigma_f\f$ 
	void applyNuSigmaf(vector<double>& a_rhs, const vector<double>& a_FNodes); 
}; 

#endif