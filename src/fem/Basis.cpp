#include "Basis.H"
#include "CH_Timer.H"

Basis::Basis() {} 
Basis::Basis(vector<Node>& a_nodeList) {

	m_nodeList = a_nodeList; 
	m_nNodes = m_nodeList.size(); 

	// initialize basis functions 
	if (m_nNodes == 3) {
		m_basis.push_back(Poly2D({1, -1, -1})); 
		m_basis.push_back(Poly2D({0, 1, 0})); 
		m_basis.push_back(Poly2D({0, 0, 1})); 
	} else if (m_nNodes == 6) {
		m_basis.push_back(Poly2D({1, -3, -3, 2, 4, 2})); 
		m_basis.push_back(Poly2D({0, -1, 0, 2, 0, 0})); 
		m_basis.push_back(Poly2D({0, 0, -1, 0, 0, 2})); 
		m_basis.push_back(Poly2D({0, 4, 0, -4, -4, 0})); 
		m_basis.push_back(Poly2D({0, 0, 0, 0, 4, 0})); 
		m_basis.push_back(Poly2D({0, 0, 4, 0, -4, -4})); 
	} else if (m_nNodes == 10) {
		// 0: 9/2(1-x-y)(1/3-x-y)(2/3-x-y)
		m_basis.push_back(Poly2D({
			1, -11./2, -11./2, 9, 18, 9, -9./2, -27./2, -27./2, -9./2
		})); 
		// 1: 9/2x(x-1/3)(x-2/3)
		m_basis.push_back(Poly2D({
			0, 1, 0, -9./2, 0, 0, 9./2, 0, 0, 0
		})); 
		// 2: 9/2y(y-1/3)(y-2/3)
		m_basis.push_back(Poly2D({
			0, 0, 1, 0, 0, -9./2, 0, 0, 0, 9./2
		})); 
		// 3: 27/2(1-x-y)x(2/3-x-y)
		m_basis.push_back(Poly2D({
			0, 9, 0, -45./2, -45./2, 0, 27./2, 27., 27./2, 0
		})); 
		// 4: 27/2(1 -x-y)x(x - 1/3)
		m_basis.push_back(Poly2D({
			0, -9./2, 0, 18., 9./2, 0, -27./2, -27./2, 0, 0
		})); 
		// 5: 27/2xy(x-1/3)
		m_basis.push_back(Poly2D({
			0, 0, 0, 0, -9./2, 0, 0, 27./2, 0, 0
		})); 
		// 6: 27/2xy(y-1/3)
		m_basis.push_back(Poly2D({
			0, 0, 0, 0, -9./2, 0, 0, 0, 27./2, 0
		})); 
		// 7: 27/2(1-x-y)y(y-1/3)
		m_basis.push_back(Poly2D({
			0, 0, -9./2, 0, 9./2, 18, 0, 0, -27./2, -27./2
		})); 
		// 8: 27/2(1-x-y)y(2/3 - x- y))
		m_basis.push_back(Poly2D({
			0, 0, 9, 0, -45./2, -45./2, 0, 27./2, 27, 27./2
		})); 
		// 9: 27xy(1 - x - y)
		m_basis.push_back(Poly2D({
			0, 0, 0, 0, 27, 0, 0, -27, -27, 0
		})); 

	}

	// initialize derivatives of basis functions 
	for (int i=0; i<m_nNodes; i++) {
		m_dbasis.push_back(m_basis[i].gradient()); 
	}
}