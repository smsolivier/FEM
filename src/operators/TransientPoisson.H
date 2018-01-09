#ifndef _TRANSIENTPOISSON_H_
#define _TRANSIENTPOISSON_H_

#include <vector>
#include "Operator.H"

using namespace std; 

/// Generates LHS and RHS for transient heat conduction 
class TransientPoisson : public Operator {
public:
	/// constructor 
	TransientPoisson(FEGrid& a_grid, 
		const Materials& a_materials, const double a_h); 
	/// makes RHS 
	void makeRHS(vector<double>& a_rhs, 
		const vector<double>& FNodes, const vector<double>& a_Told, const double h); 

}; 

#endif 