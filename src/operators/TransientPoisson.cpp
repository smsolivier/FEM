#include "TransientPoisson.H"
#include "VectorMath.H"

TransientPoisson::TransientPoisson(FEGrid& a_grid, 
	const Materials& a_materials, 
	const double a_h) : Operator(a_grid, a_materials) {

	const int nEl = m_grid.getNumElements(); 

	for (int i=0; i<nEl; i++) {
		Element& el = m_grid.getElement(i); 
		int gid = el.getGroup(); 
		double alpha = m_materials[gid]("alpha"); 
		double kappa = m_materials[gid]("k"); 
		for (int j=0; j<el.getNumNodes(); j++) {
			int jid = el.interiorNodeID(j); 
			if (jid >= 0) {
				for (int k=0; k<el.getNumNodes(); k++) {
					int kid = el.interiorNodeID(k); 
					if (kid >= 0) {
						m_matrix(jid,kid) += el.BiBj(j,k);
						m_matrix(jid,kid) += .5*a_h*alpha*kappa*el.gradBiGradBj(j,k); 
					}
				}
			}
		}
	}
}

void TransientPoisson::makeRHS(vector<double>& a_rhs, 
	const vector<double>& a_FNodes, const vector<double>& a_Told, 
	const double a_h) {

	const int nNodes = m_grid.getNumInteriorNodes(); 
	const int nEl = m_grid.getNumElements(); 

	a_rhs.resize(nNodes); 
	for (int i=0; i<nNodes; i++) a_rhs[i] = 0; 

	for (int i=0; i<nEl; i++) {
		Element& el = m_grid.getElement(i); 
		int gid = el.getGroup(); 
		double kappa = m_materials[gid]("k"); 
		double alpha = m_materials[gid]("alpha"); 
		for (int j=0; j<el.getNumNodes(); j++) {
			int jid = el.interiorNodeID(j); 
			if (jid >= 0) {
				for (int k=0; k<el.getNumNodes(); k++) {
					int kid = el.interiorNodeID(k); 
					if (kid >= 0) {
						a_rhs[jid] += el.BiBj(j,k)*a_Told[kid]; 
						a_rhs[jid] += alpha*a_h*el.BiBj(j,k)*a_FNodes[kid]; 
						a_rhs[jid] -= .5*alpha*a_h*kappa*el.gradBiGradBj(j,k)*a_Told[kid]; 
					}
				}
				// boundary integral 
				vector<double> bound = el.boundaryIntegral(j); 
				double val = 0; 
				if (el.boundaryFace(j)) {
					val -= dot(el.normal(j), bound); 
				} else if (el.boundaryFace(j-1)) {
					val -= dot(el.normal(j-1), bound); 
				}
				a_rhs[jid] -= kappa*val*a_h; 
			}
		}
	}
}