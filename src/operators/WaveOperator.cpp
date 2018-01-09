#include "WaveOperator.H"

WaveOperator::WaveOperator(FEGrid& a_grid, 
	const Materials& a_materials, const double a_h) 
	: Operator(a_grid, a_materials) {

	const int nEl = m_grid.getNumElements(); 

	for (int i=0; i<nEl; i++) {
		Element& el = m_grid.getElement(i); 
		int gid = el.getGroup(); 
		Material mat = m_materials[gid]; 
		double a = mat("a"); 
		double c= mat("c"); 
		for (int j=0; j<el.getNumNodes(); j++) {
			int jid = el.interiorNodeID(j); 
			if (jid >= 0) {
				for (int k=0; k<el.getNumNodes(); k++) {
					int kid = el.interiorNodeID(k); 
					if (kid >= 0) {
						m_matrix(jid,kid) += el.BiBj(j,k)
							+ c*.5*a_h*a_h*el.gradBiGradBj(j,k) 
							+ a*a_h*el.BiBj(j,k); 
					}
				}
			}
		}
	}
}

void WaveOperator::makeRHS(vector<double>& a_rhs, 
	const vector<double>& a_u1, const vector<double>& a_u0, 
	const vector<double>& a_F, const double a_h) {

	const int nNodes = m_grid.getNumInteriorNodes(); 
	const int nEl = m_grid.getNumElements(); 

	a_rhs.resize(nNodes); 
	for (int i=0; i<nNodes; i++) a_rhs[i] = 0; 

	for (int i=0; i<nEl; i++) {
		Element& el = m_grid.getElement(i); 
		int gid = el.getGroup(); 
		Material mat = m_materials[gid]; 
		double a = mat("a"); 
		double c= mat("c"); 
		for (int j=0; j<el.getNumNodes(); j++) {
			int jid = el.interiorNodeID(j); 
			if (jid >= 0) {
				for (int k=0; k<el.getNumNodes(); k++) {
					int kid = el.interiorNodeID(k); 
					if (kid >= 0) {
						a_rhs[jid] += 2*el.BiBj(j,k)*a_u1[kid] 
							- el.BiBj(j,k)*a_u0[kid] 
							+ a_h*a_h*el.BiBj(j,k)*a_F[kid]
							- c*.5*a_h*a_h*el.gradBiGradBj(j,k)*a_u0[kid]
							+ a*a_h*el.BiBj(j,k)*a_u1[kid]; 
					}
				}
			}
		}
	}
}