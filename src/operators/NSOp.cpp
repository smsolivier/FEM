#include "NSOp.H"
#include "CH_Timer.H"
#include "VectorMath.H"
#include <iostream>

NSOp::NSOp(FEGrid& a_grid, const Materials& a_materials) 
	: Operator(a_grid, a_materials) {

	m_stokes = false; 
}

void NSOp::buildLHS() {

	vector<double> sol(m_grid.getFields().size()); 

	m_stokes = true; 

	buildLHS(sol); 
}

void NSOp::buildLHS(vector<double>& a_sol) {

	CH_TIMERS("setup LHS"); 

	const int nEl = m_grid.getNumElements(); 
	Fields& fields = m_grid.getFields(); 

	for (int i=0; i<nEl; i++) {
		Element& el = m_grid.getElement(i); 
		int group = el.getGroup(); 
		double mu = m_materials[group]("mu"); 
		double rho = m_materials[group]("rho"); 

		// x momentum equation 
		for (int j=0; j<el.getNumNodes(); j++) {
			if (el[j].isInterior()) {
				int jid = el[j].interiorNodeID(); 
				int gjid = fields[jid]["u"]; // global id for x velocity at node jid 

				// u and v terms 
				for (int k=0; k<el.getNumNodes(); k++) { 
					if (el[k].isInterior()) {
						int kid = el[k].interiorNodeID(); 
						double u0 = a_sol[fields[kid]["u"]]; 
						double v0 = a_sol[fields[kid]["v"]]; 
						// nonlinear term 
						m_matrix(gjid, fields[kid]["u"]) += 
							rho*(u0*el.BidBjdx(j,k) + v0*el.BidBjdy(j,k)); 
						// u term 
						m_matrix(gjid, fields[kid]["u"]) += mu*el.NSx(j,k);  
						// v term 
						m_matrix(gjid, fields[kid]["v"]) += mu*el.dBidy_dBjdx(j,k); 
					}
				}

				// p term 
				for (int k=0; k<3; k++) {
					if (el[k].isInterior()) {
						int kid = el[k].interiorNodeID(); 
						m_matrix(gjid, fields[kid]["p"]) -= el.dBidxBhatj(j,k); 
					}
				}
			} 
		}

		// y momentum equation 
		for (int j=0; j<el.getNumNodes(); j++) {
			if (el[j].isInterior()) {
				int jid = el[j].interiorNodeID();  
				int gjid = fields[jid]["v"]; // global id for y velocity at node jid 

				// u and v terms 
				for (int k=0; k<el.getNumNodes(); k++) {
					if (el[k].isInterior()) {
						int kid = el[k].interiorNodeID(); 
						double u0 = a_sol[fields[kid]["u"]]; 
						double v0 = a_sol[fields[kid]["v"]]; 
						// nonlinear term 
						m_matrix(gjid, fields[kid]["v"]) += 
							rho*(u0*el.BidBjdx(j,k) + v0*el.BidBjdy(j,k)); 
						// u term 
						m_matrix(gjid, fields[kid]["u"]) += mu*el.dBidx_dBjdy(j,k); 
						// v term 
						m_matrix(gjid, fields[kid]["v"]) += mu*el.NSy(j,k); 
					}
				}

				// p term 
				for (int k=0; k<3; k++) {
					if (el[k].isInterior()) {
						int kid = el[k].interiorNodeID(); 
						m_matrix(gjid, fields[kid]["p"]) -= el.dBidyBhatj(j,k); 
					}
				}
			}
		}

		// continuity 
		for (int j=0; j<3; j++) {
			if (el[j].isInterior()) {
				int jid = el.interiorNodeID(j); 
				int gjid = fields[jid]["p"]; // global id for p at node jid 
				for (int k=0; k<el.getNumNodes(); k++) {
					if (el[k].isInterior()) {
						int kid = el[k].interiorNodeID();  
						m_matrix(gjid, fields[kid]["u"]) -= el.dBidxBhatj(k,j); 
						m_matrix(gjid, fields[kid]["v"]) -= el.dBidyBhatj(k,j); 
					}
				}
			}
		}
	}
}

void NSOp::makeRHS(vector<double>& a_rhs) {

	const int nEl = m_grid.getNumElements(); 
	Fields& fields = m_grid.getFields(); 

	a_rhs.resize(m_nUnknowns); 
	for (int i=0; i<m_nUnknowns; i++) a_rhs[i] = 0; 

	double u0 = 1; 
	double v0 = 0; 

	for (int i=0; i<nEl; i++) {
		Element& el = m_grid.getElement(i); 
		int group = el.getGroup(); 
		double mu = m_materials[group]("mu"); 
		double rho = m_materials[group]("rho"); 

		// x momentum equation 
		for (int j=0; j<el.getNumNodes(); j++) {
			if (el[j].isInterior()) {
				int jid = el[j].interiorNodeID(); 
				int gjid = fields[jid]["u"]; // global id for x velocity at node jid 

				// u and v terms 
				for (int k=0; k<el.getNumNodes(); k++) { 
					int kid = el[k].interiorNodeID(); 
					if (kid == -2) {
						// nonlinear term 
						if (m_stokes) {
							a_rhs[gjid] -= rho*u0*(u0*el.BidBjdx(j,k) + 
								v0*el.BidBjdy(j,k)); 
						}
						// u term 
						a_rhs[gjid] -= mu*el.NSx(j,k)*u0; 
						// v term 
						a_rhs[gjid] -= mu*el.dBidy_dBjdx(j,k)*v0; 
					}
				}
			} 
		}

		// y momentum equation 
		for (int j=0; j<el.getNumNodes(); j++) {
			if (el[j].isInterior()) {
				int jid = el[j].interiorNodeID();  
				int gjid = fields[jid]["v"]; // global id for y velocity at node jid 

				// u and v terms 
				for (int k=0; k<el.getNumNodes(); k++) {
					int kid = el[k].interiorNodeID(); 
					if (kid == -2) {
						// nonlinear term 
						if (m_stokes) {
							a_rhs[gjid] -= rho*v0*(u0*el.BidBjdx(j,k) + 
								v0*el.BidBjdy(j,k)); 
						}
						// u term 
						a_rhs[gjid] -= mu*el.dBidx_dBjdy(j,k)*u0; 
						// v term 
						a_rhs[gjid] -= mu*el.NSy(j,k)*v0; 
					}
				}
			}
		}

		// continuity 
		for (int j=0; j<3; j++) {
			if (el[j].isInterior()) {
				int jid = el.interiorNodeID(j); 
				int gjid = fields[jid]["p"]; // global id for p at node jid 
				for (int k=0; k<el.getNumNodes(); k++) {
					int kid = el[k].interiorNodeID(); 
					if (kid == -2) {
						a_rhs[gjid] += el.dBidxBhatj(k,j)*u0; 
						a_rhs[gjid] += el.dBidyBhatj(k,j)*v0; 
					}
				}
			}
		}
	}
}