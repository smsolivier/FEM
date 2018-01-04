#ifndef _OPERATOR_H_
#define _OPERATOR_H_

#include "SparseMatrix.H"
#include "FEGrid.H"

using namespace std; 

/// Base class for all operators 
/** Operators provide functions to build the left and right hand sides */ 
class Operator {
public:
	/// Base constructor 
	/** stores a_grid and initializes the SparseMatrix \n
		derived classes must construct the LHS in the constructor
	*/ 
	Operator(FEGrid& a_grid); 
	/// Generates the right hand side vector 
	// virtual void makeRHS(vector<double>& a_rhs, 
		// const vector<double>& a_FCentroids); 

	/// get the left hand side matrix 
	SparseMatrix& matrix(); 
	/// returns the FEGrid object 
	FEGrid& getFEGrid(); 
	
protected:
	/// stores the LHS matrix
	SparseMatrix m_matrix; 
	/// stores the FEGrid object
	FEGrid m_grid;  
	/// number of nodes 
	int m_nInteriorNodes; 
	/// number of unknowns 
	int m_nUnknowns; 
}; 

#endif