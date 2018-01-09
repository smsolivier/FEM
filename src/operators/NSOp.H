#ifndef _NSOP_H_
#define _NSOP_H_

#include "Operator.H"
#include "FEGrid.H"
#include "Materials.H"

/// set up the LHS and RHS for the navier stokes equations 
class NSOp : public Operator {
public:
	/// builds lhs for Stokes equation 
	NSOp(FEGrid& a_grid, const Materials& a_materials, const Fields& a_fields); 
	/// setup LHS for Navier Stokes 
	void buildLHS(vector<double>& a_sol); 
	/// setup LHS for Stokes 
	void buildLHS(); 
	/// builds rhs 
	void makeRHS(vector<double>& a_rhs); 
private:
	/// controls whether to include nonlinear terms in makeRHS 
	bool m_stokes; 
}; 

#endif