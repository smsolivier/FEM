#include "PoissonOperator.H"
#include "CH_Timer.H"
#include "VectorMath.H"
#include <iostream>

PoissonOperator::PoissonOperator(FEGrid& a_grid, const Materials& a_materials) 
	: Operator(a_grid, a_materials) {

	CH_TIMERS("setup LHS"); 

	const int nEl = m_grid.getNumElements(); 

	// construct element by element 
	for (int i=0; i<nEl; i++) {
		Element& el = m_grid.getElement(i); 
		double kappa = m_materials[el.getGroup()]("k"); 
		for (int j=0; j<el.getNumNodes(); j++) {
			int jid = el.interiorNodeID(j); 
			if (jid >= 0) {
				for (int k=0; k<el.getNumNodes(); k++) {
					int kid = el.interiorNodeID(k); 
					if (kid >= 0) {
						m_matrix(jid, kid) += kappa*el.gradBiGradBj(j, k); 
					}
				}
			}
		}
	}
}

void PoissonOperator::makeRHS(vector<double>& a_rhs, 
	const vector<double>& a_FNodes) {

	CH_TIMERS("setup RHS"); 

	const int nNodes = m_grid.getNumInteriorNodes(); 
	const int nEl = m_grid.getNumElements(); 

	a_rhs.resize(nNodes); 
	for (int i=0; i<nNodes; i++) a_rhs[i] = 0; 

	vector<double> Q = {1, 0}; 

	for (int i=0; i<nEl; i++) {
		Element& el = m_grid.getElement(i); 
		double kappa = m_materials[el.getGroup()]("k"); 
		for (int j=0; j<el.getNumNodes(); j++) {
			int jid = el.interiorNodeID(j); 
			if (jid >= 0) {
				for (int k=0; k<el.getNumNodes(); k++) {
					int kid = el.interiorNodeID(k); 
					if (kid >= 0) {
						a_rhs[jid] += el.BiBj(j,k) * a_FNodes[kid]; 
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
				a_rhs[jid] -= kappa*val; 
			}
		}
	}
}

void PoissonOperator::makeRHSAtCentroids(vector<double>& a_rhs, 
	const vector<double>& a_FCentroids) {

	CH_TIMERS("setup RHS"); 

	const int nNodes = m_grid.getNumInteriorNodes(); 
	const int nEl = m_grid.getNumElements(); 

	a_rhs.resize(nNodes); 

	for (int i=0; i<nEl; i++) {
		Element el = m_grid.getElement(i); 
		for (int j=0; j<el.getNumNodes(); j++) {
			int jid = el.interiorNodeID(j); 
			if (jid >= 0) {
				a_rhs[jid] += el.Bi(j) * a_FCentroids[i]; 
			}
		}
	}
}