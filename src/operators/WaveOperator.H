#ifndef _WAVEOPERATOR_H_
#define _WAVEOPERATOR_H_

#include "Operator.H"
#include "FEGrid.H"
#include "Material.H"
#include <vector> 

using namespace std; 

/// Sets up LHS and RHS for the wave equation 
/** 
	\f$\frac{\partial^2 u}{\partial t^2} + a\frac{\partial u}{\partial t} = c\nabla^2 u \f$ 
*/  
class WaveOperator : public Operator {
public:
	WaveOperator(FEGrid& a_grid, 
		const Materials& a_materials, const double a_h); 
	void makeRHS(vector<double>& a_rhs, 
		const vector<double>& a_u, const vector<double>& a_v,
		const vector<double>& a_F, const double h);  
}; 

#endif 