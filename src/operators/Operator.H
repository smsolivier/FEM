#ifndef _OPERATOR_H_
#define _OPERATOR_H_

#include "SparseMatrix.H"
#include "Materials.H"
#include "FEGrid.H"
#include "Fields.H"

using namespace std; 

/// Base class for all operators 
/** Operators provide functions to build the left and right hand sides */ 
class Operator {
public:
	/// Base constructor 
	/** stores a_grid and initializes the SparseMatrix \n
		derived classes must construct the LHS in the constructor
	*/ 
	Operator(FEGrid& a_grid, const Materials& a_materials, const Fields& a_fields); 
	/// base constructor for programs that don't require the fields class 
	Operator(FEGrid& a_grid, const Materials& a_materials); 

	/// get the left hand side matrix 
	SparseMatrix& matrix(); 
	/// returns the FEGrid object 
	FEGrid& getFEGrid(); 
	/// return the Materials object 
	Materials& getMaterials(); 
	
protected:
	/// stores the LHS matrix
	SparseMatrix m_matrix; 
	/// stores the FEGrid object
	FEGrid m_grid;  
	/// number of nodes 
	int m_nInteriorNodes; 
	/// number of unknowns 
	int m_nUnknowns; 
	/// store collection of materials 
	Materials m_materials; 
	/// store the fields object 
	Fields m_fields; 
}; 

#endif