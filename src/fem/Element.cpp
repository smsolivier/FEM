#include "Element.H"
#include "CH_Timer.H"
#include "VectorMath.H"
#include <iostream>
#include <cassert>
#include <cmath>

Element::Element(vector<Node>& a_nodeList, int a_group) {

	CH_TIMERS("element initialization"); 

	m_group = a_group; 

	m_nodeList = a_nodeList; 
	m_nNodes = m_nodeList.size(); 

	// set the polynomial order of the element 
	// solves n_p = (p+1)(p+2)/2
	m_pOrder = -3/2 + sqrt(9/4 - 2*(1 - m_nNodes)); 
	assert(m_pOrder > 0); 

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

	// compute centroid 
	m_centroid.resize(2,0); 
	for (int i=0; i<2; i++) {
		for (int j=0; j<m_nNodes; j++) {
			m_centroid[i] += m_basis[j](1./3, 1./3) 
				* m_nodeList[j].getPosition()[i]; 
		}
	}

	// set up integrators
	int p = 1; 
	if (m_nNodes == 6) p = 6; 
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
}

double Element::gradBiGradBj(int a_i, int a_j) {

	CH_TIMERS("gradBiGradBj"); 

	if (m_nNodes == 3) {
		if (m_gradBiGradBj.size() == 0) {
			m_gradBiGradBj.resize(m_nNodes); 
			for (int i=0; i<m_nNodes; i++) m_gradBiGradBj[i].resize(m_nNodes); 
			for (int i=0; i<m_nNodes; i++) {
				for (int j=i; j<m_nNodes; j++) {
					vector<double> grad_i = gradient(i, 0, 0); 
					vector<double> grad_j = gradient(j, 0, 0); 
					// symmetry 
					double integral = dot(grad_i, grad_j)*detJ(0,0)/2;
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
						vector<double> grad_i = gradient(i, xi, eta); 
						vector<double> grad_j = gradient(j, xi, eta); 
						return dot(grad_i, grad_j) * detJ(xi, eta); 
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

	CH_TIMERS("BiBj"); 

	if (m_nNodes == 3) {
		return 1./9*detJ(0,0)/2; 
	} else {
		if (m_BiBj.size() == 0) {
			m_BiBj.resize(m_nNodes); 
			for (int i=0; i<m_nNodes; i++) m_BiBj[i].resize(m_nNodes); 
			for (int i=0; i<m_nNodes; i++) {
				for (int j=i; j<m_nNodes; j++) {
					auto f = [this, i, j] (double xi, double eta) {
						return m_basis[i](xi, eta) * 
							m_basis[j](xi, eta) * detJ(xi, eta); 
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

	CH_TIMERS("Bi"); 

	if (m_nNodes == 3) {
		return 1./3*detJ(0,0)/2; 
	} else {
		if (m_Bi.size() == 0) {
			m_Bi.resize(m_nNodes);
			for (int i=0; i<m_nNodes; i++) {
				auto f = [this, i] (double xi, double eta) {
					return m_basis[i](xi, eta) * detJ(xi, eta); 
				}; 
				m_Bi[i] = m_gq.integrate(f); 
			}
		}
		return m_Bi[a_i]; 
	}
}

vector<double> Element::boundaryIntegral(int a_i) {

	CH_TIMERS("boundary integral"); 

	vector<double> ret(2, 0); 

	int f1 = a_i; 
	// if (a_i == 3) f1 = 0; 
	// else if (a_i == 4) f1 = 1; 
	// else if (a_i == 5) f1 = 2; 

	if (boundaryFace(f1)) {

		double integral = 0; 
		if (f1 == 0) {
			auto f = [this, a_i] (double xi) {
				return m_basis[a_i](xi, 0)*
					sqrt(pow(J(0,0,xi,0),2)+pow(J(0,1,xi,0),2)); 
			}; 
			integral = m_gq1d(f);
		} else if (f1 == 1) {
			auto f = [this, a_i] (double xi) {
				return m_basis[a_i](xi, 1 - xi)*
					sqrt(pow(J(0,1,xi,1-xi),2) + pow(J(1,0,xi,1-xi),2)); 
			}; 
			integral = m_gq1d(f); 
		} else if (f1 == 2) {
			auto f = [this, a_i] (double xi) {
				return m_basis[a_i](0, xi)*
					sqrt(pow(J(1,0,0,xi),2) + pow(J(1,1,0,xi),2)); 
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
				return m_basis[a_i](xi, 0)*
					sqrt(pow(J(0,0,xi,0),2)+pow(J(0,1,xi,0),2)); 
			}; 
			integral = m_gq1d(f);
		} else if (f2 == 1) {
			auto f = [this, a_i] (double xi) {
				return m_basis[a_i](xi, 1 - xi)*
					sqrt(pow(J(0,1,xi,1-xi),2) + pow(J(1,0,xi,1-xi),2)); 
			}; 
			integral = m_gq1d(f); 
		} else if (f2 == 2) {
			auto f = [this, a_i] (double xi) {
				return m_basis[a_i](0, xi)*
					sqrt(pow(J(1,0,0,xi),2) + pow(J(1,1,0,xi),2)); 
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
vector<double> Element::getCentroid() {return m_centroid; }
int Element::getNumNodes() const {return m_nNodes; }

double Element::getArea() {

	auto f = [this] (double xi, double eta) {return detJ(xi, eta); }; 

	return m_gq.integrate(f); 

}
int Element::getOrder() {return m_pOrder; } 
int Element::getGroup() {return m_group; }

vector<double> Element::gradient(int a_i, double xi, double eta) {

	CH_TIMERS("gradient"); 

	vector<double> ret(2); 
	ret[0] = Jinv(0,0,xi,eta)*m_dbasis[a_i][0](xi,eta)
		 + Jinv(0,1,xi,eta)*m_dbasis[a_i][1](xi,eta); 
	ret[1] = Jinv(1,0,xi,eta)*m_dbasis[a_i][0](xi,eta)
		 + Jinv(1,1,xi,eta)*m_dbasis[a_i][1](xi,eta); 

	return ret; 
}

double Element::J(int a_i, int a_j, double xi, double eta) {

	double val = 0; 
	for (int i=0; i<m_nNodes; i++) {
		val += m_dbasis[i][a_i](xi, eta) * m_nodeList[i].getPosition()[a_j]; 
	}

	return val; 
}

double Element::Jinv(int i, int j, double xi, double eta) {

	double val = 0; 
	if (i==0 && j==0) {
		val = 1/detJ(xi,eta)*J(1,1,xi,eta); 
	} else if (i==0 && j==1) {
		val = -1/detJ(xi,eta)*J(0,1,xi,eta); 
	} else if (i==1 && j==0) {
		val = -1/detJ(xi,eta)*J(1,0,xi,eta); 
	} else if (i==1 && j==1) {
		val = 1/detJ(xi,eta)*J(0,0,xi,eta); 
	} else {
		cout << "ERROR: jacobian indices out of range in Element.cpp" << endl; 
		exit(0); 
	}

	return val; 

}

double Element::detJ(double xi, double eta) {

	double det = J(0,0,xi,eta)*J(1,1,xi,eta) - J(0,1,xi,eta)*J(1,0,xi,eta); 
	if (det == 0) {
		cout << "ERROR: |Jacobian| zero in Element" << endl; 
		exit(0); 
	}

	return abs(det); 
}