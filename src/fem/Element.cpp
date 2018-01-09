#include "Element.H"
// #include "CH_Timer.H"
#include "VectorMath.H"
#include <iostream>
#include <cassert>
#include <cmath>

Element::Element(vector<Node>& a_nodeList, int a_group) {

	// CH_TIMERS("element initialization"); 

	m_group = a_group; 

	m_nodeList = a_nodeList; 
	m_nNodes = m_nodeList.size(); 

	// set the polynomial order of the element 
	// solves n_p = (p+1)(p+2)/2
	m_pOrder = -3/2 + sqrt(9/4 - 2*(1 - m_nNodes)); 
	assert(m_pOrder > 0);

	// initialize basis functions 
	m_vbasis = Basis(m_nodeList); 
	vector<Node> pnodes = {m_nodeList[0], m_nodeList[1], m_nodeList[2]}; 
	m_pbasis = Basis(pnodes); 

	// compute centroid 
	m_centroid.resize(2,0); 
	for (int i=0; i<2; i++) {
		for (int j=0; j<m_nNodes; j++) {
			m_centroid[i] += m_vbasis[j](1./3, 1./3) 
				* m_nodeList[j].getPosition()[i]; 
		}
	}

	// set up integrators
	int p = 1; 
	if (m_nNodes == 6) p = 4; 
	else if (m_nNodes == 10) p = 8; 
	m_gq = GQ(p); 
	m_gq1d = GQ1D(4); 

	// determine boundary faces 
	m_boundaryFaces.resize(3, false); 
	if (!m_nodeList[0].isInterior() && !m_nodeList[1].isInterior()) {
		m_boundaryFaces[0] = true; 
	} else if (!m_nodeList[1].isInterior() && !m_nodeList[2].isInterior()) {
		m_boundaryFaces[1] = true; 
	} else if (!m_nodeList[2].isInterior() && !m_nodeList[0].isInterior()) {
		m_boundaryFaces[2] = true; 
	}

	// side lengths 
	m_h.resize(3); 
	vector<double> vec = m_nodeList[1].getPosition()
		- m_nodeList[0].getPosition(); 
	m_h[0] = sqrt(dot(vec, vec)); 
	vec = m_nodeList[2].getPosition() - m_nodeList[1].getPosition(); 
	m_h[1] = sqrt(dot(vec, vec)); 
	vec = m_nodeList[0].getPosition() - m_nodeList[2].getPosition(); 
	m_h[2] = sqrt(dot(vec, vec)); 

	// initialize storage 
	m_NSx.resize(m_nNodes); 
	m_NSy.resize(m_nNodes); 
	m_dBidy_dBjdx.resize(m_nNodes); 
	m_dBidx_dBjdy.resize(m_nNodes); 
	m_dBidxBhatj.resize(m_nNodes); 
	m_dBidyBhatj.resize(m_nNodes);
	m_BidBjdx.resize(m_nNodes); 
	m_BidBjdy.resize(m_nNodes); 
	for (int i=0; i<m_nNodes; i++) {
		m_NSx[i].resize(m_nNodes); 
		m_NSy[i].resize(m_nNodes); 
		m_dBidy_dBjdx[i].resize(m_nNodes); 
		m_dBidx_dBjdy[i].resize(m_nNodes); 
		m_dBidxBhatj[i].resize(3); 
		m_dBidyBhatj[i].resize(3);
		m_BidBjdx[i].resize(m_nNodes); 
		m_BidBjdy[i].resize(m_nNodes); 
	}

	// initialize booleans 
	m_fNSx = true; 
	m_fNSy = true; 
	m_fdBidy_dBjdx = true; 
	m_fdBidx_dBjdy = true; 
	m_fdBidxBhatj = true; 
	m_fdBidyBhatj = true; 
	m_fBidBjdx = true; 
	m_fBidBjdy = true; 
}

double Element::gradBiGradBj(int a_i, int a_j) {

	// CH_TIMERS("gradBiGradBj"); 

	if (m_nNodes == 3) {
		if (m_gradBiGradBj.size() == 0) {
			m_gradBiGradBj.resize(m_nNodes); 
			for (int i=0; i<m_nNodes; i++) m_gradBiGradBj[i].resize(m_nNodes); 
			for (int i=0; i<m_nNodes; i++) {
				for (int j=i; j<m_nNodes; j++) {
					vector<double> grad_i = m_vbasis.gradient(i, 0, 0); 
					vector<double> grad_j = m_vbasis.gradient(j, 0, 0); 
					// symmetry 
					double integral = dot(grad_i, grad_j)*m_vbasis.detJ(0,0)/2;
					m_gradBiGradBj[i][j] = integral;
					m_gradBiGradBj[j][i] = integral;  
				}
			}
		}
		return m_gradBiGradBj[a_i][a_j]; 
	} else {
		if (m_gradBiGradBj.size() == 0) {
			m_gradBiGradBj.resize(m_nNodes); 
			for (int i=0; i<m_nNodes; i++) m_gradBiGradBj[i].resize(m_nNodes); 
			for (int i=0; i<m_nNodes; i++) {
				for (int j=i; j<m_nNodes; j++) {
					auto f = [this, i, j] (double xi, double eta) {
						vector<double> grad_i = m_vbasis.gradient(i, xi, eta); 
						vector<double> grad_j = m_vbasis.gradient(j, xi, eta); 
						return dot(grad_i, grad_j) * m_vbasis.detJ(xi, eta); 
					}; 
					// take advantage of symmetry 
					// grad B_i grad B_j = grad B_j grad B_i 
					m_gradBiGradBj[i][j] = m_gq.integrate(f);;
					m_gradBiGradBj[j][i] = m_gradBiGradBj[i][j]; 
				}
			}
		} 
		return m_gradBiGradBj[a_i][a_j]; 
	}

}

double Element::BiBj(int a_i, int a_j) {

	// CH_TIMERS("BiBj"); 

	if (m_nNodes == 3) {
		return 1./9*m_vbasis.detJ(0,0)/2; 
	} else {
		if (m_BiBj.size() == 0) {
			m_BiBj.resize(m_nNodes); 
			for (int i=0; i<m_nNodes; i++) m_BiBj[i].resize(m_nNodes); 
			for (int i=0; i<m_nNodes; i++) {
				for (int j=i; j<m_nNodes; j++) {
					auto f = [this, i, j] (double xi, double eta) {
						return m_vbasis[i](xi, eta) * 
							m_vbasis[j](xi, eta) * m_vbasis.detJ(xi, eta); 
					}; 
					// symmetry 
					m_BiBj[i][j] = m_gq.integrate(f); 
					m_BiBj[j][i] = m_BiBj[i][j]; 
				}
			}
		}

		return m_BiBj[a_i][a_j]; 
	}

}

double Element::Bi(int a_i) {

	// CH_TIMERS("Bi"); 

	if (m_nNodes == 3) {
		return 1./3*m_vbasis.detJ(0,0)/2; 
	} else {
		if (m_Bi.size() == 0) {
			m_Bi.resize(m_nNodes);
			for (int i=0; i<m_nNodes; i++) {
				auto f = [this, i] (double xi, double eta) {
					return m_vbasis[i](xi, eta) * m_vbasis.detJ(xi, eta); 
				}; 
				m_Bi[i] = m_gq.integrate(f); 
			}
		}
		return m_Bi[a_i]; 
	}
}

vector<double> Element::boundaryIntegral(int a_i) {

	// CH_TIMERS("boundary integral"); 

	vector<double> ret(2, 0); 

	int f1 = a_i; 
	// if (a_i == 3) f1 = 0; 
	// else if (a_i == 4) f1 = 1; 
	// else if (a_i == 5) f1 = 2; 

	if (boundaryFace(f1)) {

		double integral = 0; 
		if (f1 == 0) {
			auto f = [this, a_i] (double xi) {
				return m_vbasis[a_i](xi, 0)*
					sqrt(pow(m_vbasis.J(0,0,xi,0),2)+pow(m_vbasis.J(0,1,xi,0),2)); 
			}; 
			integral = m_gq1d(f);
		} else if (f1 == 1) {
			auto f = [this, a_i] (double xi) {
				return m_vbasis[a_i](xi, 1 - xi)*
					sqrt(pow(m_vbasis.J(0,1,xi,1-xi),2) + pow(m_vbasis.J(1,0,xi,1-xi),2)); 
			}; 
			integral = m_gq1d(f); 
		} else if (f1 == 2) {
			auto f = [this, a_i] (double xi) {
				return m_vbasis[a_i](0, xi)*
					sqrt(pow(m_vbasis.J(1,0,0,xi),2) + pow(m_vbasis.J(1,1,0,xi),2)); 
			}; 
			integral = m_gq1d(f); 
		} else {
			cout << "a_i not defined" << endl; 
			exit(0); 
		}

		ret += normal(f1)*integral; 
	} 

	int f2 = f1 - 1; 
	if (f2 == -1) f2 = 2; 
	if (boundaryFace(f2)) {

		double integral = 0; 
		if (f2 == 0) {
			auto f = [this, a_i] (double xi) {
				return m_vbasis[a_i](xi, 0)*
					sqrt(pow(m_vbasis.J(0,0,xi,0),2)+pow(m_vbasis.J(0,1,xi,0),2)); 
			}; 
			integral = m_gq1d(f);
		} else if (f2 == 1) {
			auto f = [this, a_i] (double xi) {
				return m_vbasis[a_i](xi, 1 - xi)*
					sqrt(pow(m_vbasis.J(0,1,xi,1-xi),2) + pow(m_vbasis.J(1,0,xi,1-xi),2)); 
			}; 
			integral = m_gq1d(f); 
		} else if (f2 == 2) {
			auto f = [this, a_i] (double xi) {
				return m_vbasis[a_i](0, xi)*
					sqrt(pow(m_vbasis.J(1,0,0,xi),2) + pow(m_vbasis.J(1,1,0,xi),2)); 
			}; 
			integral = m_gq1d(f); 
		} else {
			cout << "a_i not defined" << endl; 
			exit(0); 
		}

		ret += normal(f2)*integral; 
	}

	return ret; 

}

double Element::NSx(int a_i, int a_j) {

	// CH_TIMERS("NSx"); 

	// if first time through 
	if (m_fNSx) {
		m_NSx.resize(m_nNodes); 

		for (int i=0; i<m_nNodes; i++) {
			for (int j=0; j<m_nNodes; j++) {
				auto f = [this, i, j] (double xi, double eta) {
					vector<double> grad_i = m_vbasis.gradient(i, xi, eta); 
					vector<double> grad_j = m_vbasis.gradient(j, xi, eta); 
					return (2*grad_i[0]*grad_j[0] + grad_i[1]*grad_j[1]) * m_vbasis.detJ(xi, eta); 
				};
				m_NSx[i][j] = m_gq(f); 
			}
		}
		m_fNSx = false; 
	}

	return m_NSx[a_i][a_j]; 
}

double Element::NSy(int a_i, int a_j) {

	// CH_TIMERS("NSy"); 

	// if first time through 
	if (m_fNSy) {
		m_NSy.resize(m_nNodes); 

		for (int i=0; i<m_nNodes; i++) {
			for (int j=0; j<m_nNodes; j++) {
				auto f = [this, i, j] (double xi, double eta) {
					vector<double> grad_i = m_vbasis.gradient(i, xi, eta); 
					vector<double> grad_j = m_vbasis.gradient(j, xi, eta); 
					return (2*grad_i[1]*grad_j[1] + grad_i[0]*grad_j[0]) * m_vbasis.detJ(xi, eta); 
				};
				m_NSy[i][j] = m_gq(f); 
			}
		}
		m_fNSy = false; 
	}

	return m_NSy[a_i][a_j]; 
}

double Element::dBidy_dBjdx(int a_i, int a_j) {

	// CH_TIMERS("dBidy_dBjdx"); 

	if (m_fdBidy_dBjdx) {

		for (int i=0; i<m_nNodes; i++) {
			for (int j=0; j<m_nNodes; j++) {
				auto f = [this, i, j] (double xi, double eta) {
					return m_vbasis.dy(i,xi,eta) * m_vbasis.dx(j,xi,eta) * m_vbasis.detJ(xi, eta); 
				}; 
				m_dBidy_dBjdx[i][j] = m_gq(f); 
			}
		}
		m_fdBidy_dBjdx = false; 
	}

	return m_dBidy_dBjdx[a_i][a_j]; 
}

double Element::dBidx_dBjdy(int a_i, int a_j) {

	// CH_TIMERS("dBidx_dBjdy"); 

	if (m_fdBidx_dBjdy) {

		for (int i=0; i<m_nNodes; i++) {
			for (int j=0; j<m_nNodes; j++) {
				auto f = [this, i, j] (double xi, double eta) {
					return m_vbasis.dx(i,xi,eta)*m_vbasis.dy(j,xi,eta)*m_vbasis.detJ(xi, eta); 
				}; 
				m_dBidx_dBjdy[i][j] = m_gq(f); 
			}
		}
		m_fdBidx_dBjdy = false; 
	}

	return m_dBidx_dBjdy[a_i][a_j]; 
}

double Element::dBidxBhatj(int a_i, int a_j) {

	// CH_TIMERS("dBidxBhatj"); 

	if (m_fdBidxBhatj) {

		for (int i=0; i<m_nNodes; i++) {
			for (int j=0; j<3; j++) {
				auto f = [this, i, j] (double xi, double eta) {
					return m_vbasis.dx(i,xi,eta)*m_pbasis[j](xi,eta)*m_vbasis.detJ(xi,eta); 
				}; 
				m_dBidxBhatj[i][j] = m_gq(f); 
			}
		}
		m_fdBidxBhatj = false; 
	}

	return m_dBidxBhatj[a_i][a_j]; 
}

double Element::dBidyBhatj(int a_i, int a_j) {

	// CH_TIMERS("dBidyBhatj"); 

	if (m_fdBidyBhatj) {

		for (int i=0; i<m_nNodes; i++) {
			for (int j=0; j<3; j++) {
				auto f = [this, i, j] (double xi, double eta) {
					return m_vbasis.dy(i,xi,eta)*m_pbasis[j](xi,eta)*m_vbasis.detJ(xi,eta); 
				}; 
				m_dBidyBhatj[i][j] = m_gq(f); 
			}
		}
		m_fdBidyBhatj = false; 
	}

	return m_dBidyBhatj[a_i][a_j]; 
}

double Element::BidBjdx(int a_i, int a_j) {

	if (m_fBidBjdx) {

		for (int i=0; i<m_nNodes; i++) {
			for (int j=0; j<m_nNodes; j++) {
				auto f = [this, i, j] (double xi, double eta) {
					return m_vbasis[i](xi, eta)*m_vbasis.dx(j, xi, eta)*m_vbasis.detJ(xi, eta); 
				}; 
				m_BidBjdx[i][j] = m_gq(f); 
			}
		}
		m_fBidBjdx = false; 
	}

	return m_BidBjdx[a_i][a_j]; 
}

double Element::BidBjdy(int a_i, int a_j) {

	if (m_fBidBjdy) {

		for (int i=0; i<m_nNodes; i++) {
			for (int j=0; j<m_nNodes; j++) {
				auto f = [this, i, j] (double xi, double eta) {
					return m_vbasis[i](xi, eta)*m_vbasis.dy(j, xi, eta)*m_vbasis.detJ(xi, eta); 
				}; 
				m_BidBjdy[i][j] = m_gq(f); 
			}
		}
		m_fBidBjdy = false; 
	}

	return m_BidBjdy[a_i][a_j]; 
}

bool Element::boundaryFace(int a_i) {

	// loop face numbering back around 
	if (a_i == -1) a_i = 2;

	return m_boundaryFaces[a_i]; 

}

vector<double> Element::normal(int a_i) {

	if (a_i == -1) a_i = 2; 

	vector<double> normal(2); 
	int i, j; 
	if (a_i == 0) {
		i = 1; 
		j = 0; 
	} else if (a_i == 1) {
		i = 2; 
		j = 1; 
	} else if (a_i == 2) {
		i = 0; 
		j = 2; 
	} else {
		cout << "ERROR: face number not defined" << endl; 
		exit(0); 
	}

	// P_i - P_j -> flip and multiply y by -1 
	vector<double> Pi = m_nodeList[i].getPosition(); 
	vector<double> Pj = m_nodeList[j].getPosition(); 

	normal = Pi - Pj; // difference 
	// swap x and y
	double tmp = normal[0]; 
	normal[0] = normal[1]; 
	normal[1] = tmp; 
	// negative y 
	normal[1] *= -1; 

	// normalize the vector 
	normal /= sqrt(dot(normal, normal)); 

	return normal; 
}

vector<Node> Element::nodeList() const {return m_nodeList; }
int Element::interiorNodeID(int i) {return m_nodeList[i].getInteriorNodeID(); }
int Element::getNodeID(int i) {return m_nodeList[i].getNodeID(); }
Node& Element::getNode(int i) {return m_nodeList[i]; }
Node& Element::operator[](int i) {return m_nodeList[i]; }
vector<double> Element::getCentroid() {return m_centroid; }
int Element::getNumNodes() const {return m_nNodes; }

double Element::getArea() {

	auto f = [this] (double xi, double eta) {return m_vbasis.detJ(xi, eta); }; 

	return m_gq.integrate(f); 

}
int Element::getOrder() {return m_pOrder; } 
int Element::getGroup() {return m_group; }