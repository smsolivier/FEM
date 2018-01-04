#include "CH_Timer.H"
#include <iostream>
#include <cmath>

using namespace std; 

Poly2D& Basis::operator[](int a_i) {
	return m_basis[a_i]; 
}

Poly2D& Basis::operator()(int a_i, int a_j) {
	return m_dbasis[a_i][a_j]; 
}

vector<double> Basis::gradient(int a_i, double xi, double eta) {

	// CH_TIMERS("gradient"); 

	vector<double> ret(2); 
	ret[0] = dx(a_i, xi, eta); 
	ret[1] = dy(a_i, xi, eta); 

	return ret; 
}

double Basis::dx(int a_i, double xi, double eta) {

	// CH_TIMERS("dx"); 
	return Jinv(0,0,xi,eta)*m_dbasis[a_i][0](xi,eta)
		+ Jinv(0,1,xi,eta)*m_dbasis[a_i][1](xi,eta); 
}

double Basis::dy(int a_i, double xi, double eta) {

	// CH_TIMERS("dy"); 
	return Jinv(1,0,xi,eta)*m_dbasis[a_i][0](xi,eta)
		+ Jinv(1,1,xi,eta)*m_dbasis[a_i][1](xi,eta);
}

double Basis::J(int a_i, int a_j, double xi, double eta) {

	double val = 0; 
	for (int i=0; i<m_nNodes; i++) {
		val += m_dbasis[i][a_i](xi, eta) * m_nodeList[i].getPosition()[a_j]; 
	}

	return val; 
}

double Basis::Jinv(int i, int j, double xi, double eta) {

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
		cout << "ERROR: jacobian indices out of range in Basis.cpp" << endl; 
		exit(0); 
	}

	return val; 

}

double Basis::detJ(double xi, double eta) {

	double det = J(0,0,xi,eta)*J(1,1,xi,eta) - J(0,1,xi,eta)*J(1,0,xi,eta); 
	if (det == 0) {
		cout << "ERROR: |Jacobian| zero in Basis" << endl; 
		exit(0); 
	}

	return abs(det); 
}