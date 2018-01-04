#ifndef _NSOP_H_
#define _NSOP_H_

#include "Operator.H"
#include "FEGrid.H"

/// set up the LHS and RHS for the navier stokes equations 
class NSOp : public Operator {
public:
	/// builds lhs 
	NSOp(FEGrid& a_grid); 
	/// builds rhs 
	void makeRHS(vector<double>& a_rhs); 
}; 

#endif