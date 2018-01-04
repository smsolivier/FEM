#include "Operator.H"

Operator::Operator(FEGrid& a_grid) : m_grid(a_grid) {

	// store number of interior nodes 
	m_nInteriorNodes = m_grid.getNumInteriorNodes(); 

	// store number of unknowns 
	m_nUnknowns = m_grid.getFields().size(); 

	// initialize the sparse matrix 
	m_matrix = SparseMatrix(m_nUnknowns, m_nUnknowns);

}

SparseMatrix& Operator::matrix() {return m_matrix; } 
FEGrid& Operator::getFEGrid() {return m_grid; }