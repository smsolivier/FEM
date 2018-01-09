#ifndef _POISSON_OPERATOR_H_
#define _POISSON_OPERATOR_H_

#include <vector>

#include "Operator.H"
#include "FEGrid.H"
#include "Material.H"

using namespace std; 

/// Sets up the LHS and RHS for the poisson equation 
class PoissonOperator : public Operator {
public:
	/// Builds LHS 
	PoissonOperator(FEGrid& a_grid, const Materials& a_materials); 
	/// Builds RHS 
	void makeRHS(vector<double>& a_rhs, 
		const vector<double>& a_FNodes); 
	void makeRHSAtCentroids(vector<double>& a_rhs, 
		const vector<double>& a_FCentroids); 
private:

}; 

#endif 