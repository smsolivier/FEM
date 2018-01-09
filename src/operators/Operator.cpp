#include "Operator.H"

Operator::Operator(FEGrid& a_grid, const Materials& a_materials, 
	const Fields& a_fields) 
	: m_grid(a_grid), m_materials(a_materials), m_fields(a_fields) 
{
	// store number of interior nodes 
	m_nInteriorNodes = m_grid.getNumInteriorNodes(); 

	// store number of unknowns 
	m_nUnknowns = m_fields.size(); 

	// initialize the sparse matrix 
	m_matrix = SparseMatrix(m_nUnknowns, m_nUnknowns);
}

Operator::Operator(FEGrid& a_grid, const Materials& a_materials)
	: m_grid(a_grid), m_materials(a_materials) {

	// store number of interior nodes 
	m_nInteriorNodes = m_grid.getNumInteriorNodes(); 

	// store number of unknowns 
	m_nUnknowns = a_grid.getNumInteriorNodes(); 

	// initialize the sparse matrix 
	m_matrix = SparseMatrix(m_nUnknowns, m_nUnknowns);
}

SparseMatrix& Operator::matrix() {return m_matrix; } 
FEGrid& Operator::getFEGrid() {return m_grid; }
Materials& Operator::getMaterials() {return m_materials; }