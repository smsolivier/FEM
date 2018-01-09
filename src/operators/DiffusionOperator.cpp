#include "DiffusionOperator.H"
#include "CH_Timer.H"
#include <cassert>
#include <iostream>

DiffusionOperator::DiffusionOperator(FEGrid& a_grid, 
	const Materials& a_materials) : Operator(a_grid, a_materials) {

	CH_TIMERS("setup LHS"); 

	const int nEl = m_grid.getNumElements(); 

	// construct element by element 
	for (int i=0; i<nEl; i++) {
		Element& el = m_grid.getElement(i); 
		int gid = el.getGroup();
		const Material& mat = m_materials[gid]; 
		double D = 1./(3*mat("Sigma_t")); 
		double Sigma_a = mat("Sigma_a"); 
		for (int j=0; j<el.getNumNodes(); j++) {
			int jid = el.interiorNodeID(j); 
			if (jid >= 0) {
				for (int k=0; k<el.getNumNodes(); k++) {
					int kid = el.interiorNodeID(k); 
					if (kid >= 0) {
						m_matrix(jid, kid) += D*el.gradBiGradBj(j, k) 
							+ Sigma_a * el.BiBj(j, k); 
					}
				}
			}
		}
	}
}

void DiffusionOperator::makeRHS(vector<double>& a_rhs, 
	const vector<double>& a_FNodes) {

	CH_TIMERS("setup RHS"); 

	const int nNodes = m_grid.getNumInteriorNodes(); 
	const int nEl = m_grid.getNumElements(); 

	a_rhs.resize(nNodes); 
	for (int i=0; i<nNodes; i++) a_rhs[i] = 0; 

	for (int i=0; i<nEl; i++) {
		Element& el = m_grid.getElement(i); 
		for (int j=0; j<el.getNumNodes(); j++) {
			int jid = el.interiorNodeID(j); 
			if (jid >= 0) {
				for (int k=0; k<el.getNumNodes(); k++) {
					int kid = el.interiorNodeID(k); 
					if (kid >= 0) {
						a_rhs[jid] += el.BiBj(j,k) * a_FNodes[kid]; 
					}
				}
			}
		}
	}
}

void DiffusionOperator::applyNuSigmaf(vector<double>& a_rhs, 
	const vector<double>& a_FNodes) {

	CH_TIMERS("apply nu sigmaf"); 

	const int nNodes = m_grid.getNumInteriorNodes(); 
	const int nEl = m_grid.getNumElements(); 

	a_rhs.resize(nNodes); 
	for (int i=0; i<nNodes; i++) a_rhs[i] = 0; 

	for (int i=0; i<nEl; i++) {
		Element& el = m_grid.getElement(i); 
		int gid = el.getGroup(); 
		double nuSigmaf = m_materials[gid]("nuSigma_f"); 
		for (int j=0; j<el.getNumNodes(); j++) {
			int jid = el.interiorNodeID(j); 
			if (jid >= 0) {
				for (int k=0; k<el.getNumNodes(); k++) {
					int kid = el.interiorNodeID(k); 
					if (kid >= 0) {
						a_rhs[jid] += nuSigmaf*el.BiBj(j,k) * a_FNodes[kid]; 
					}
				}
			}
		}
	}
}