#include "GQ1D.H"
#include <cassert>
#include <iostream>

GQ1D::GQ1D(int a_p) {

	m_p = a_p; 

	// weights and points scaled to integrate between 0 and 1 
	if (m_p == 2) {
		m_w = {.5, .5}; 
		m_x = {0.21132487, 0.78867513}; 
	} else if (m_p == 3) {
		m_w = {0.27777778, 0.44444444, 0.27777778}; 
		m_x = {0.11270167, 0.5, 0.88729833}; 
	} else if (m_p == 4) {
		m_w = {0.17392742,  0.32607258,  0.32607258,  0.17392742}; 
		m_x = {0.06943184,  0.33000948,  0.66999052,  0.93056816}; 
	} else {
		cout << "ERROR: GQ1D order not defined" << endl; 
		exit(0); 
	}

	assert(m_w.size() == m_x.size()); 
	m_Np = m_w.size(); 

}