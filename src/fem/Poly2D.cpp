#include "Poly2D.H"
#include <cmath>
#include <cassert>
#include <iostream>

Poly2D::Poly2D() {
	
}

Poly2D::Poly2D(vector<double> a_coef) {

	m_coef = a_coef;
	int n = m_coef.size(); 

	assert(n == 1 || n == 3 || n == 6 || n == 10); 
}

double Poly2D::operator()(double x, double y) {

	double sum = m_coef[0]; 
	int ind = 1; 
	int lev = 1; 
	while (ind < m_coef.size()) {
		for (int i=0; i<lev + 1; i++) {
			int x_pow = lev - i; 
			sum += m_coef[ind] * pow(x, x_pow) * pow(y, i); 
			ind++; 
		}
		lev++; 
	}

	return sum; 

}

vector<Poly2D> Poly2D::gradient() {

	vector<Poly2D> ret; 
	ret.push_back((*this).dx()); 
	ret.push_back((*this).dy()); 

	return ret; 
}

Poly2D Poly2D::dx() {

	int p = m_coef.size(); 
	vector<double> c; 

	if (p == 1) {
		c = {0}; 
	} else if (p == 3) {
		c = {m_coef[1]}; 
	} else if (p == 6) {
		c = {m_coef[1], 2*m_coef[3], m_coef[4]}; 
	} else if (p == 10) {
		c = {m_coef[1], 2*m_coef[3], m_coef[4], 3*m_coef[6], 2*m_coef[7], m_coef[8]}; 
	} else {
		cout << "ERROR: polynomial order not defined in Poly2D derivative" << endl; 
		exit(0); 
	}

	return Poly2D(c); 

}

Poly2D Poly2D::dy() {

	int p = m_coef.size(); 
	vector<double> c; 

	if (p == 1) {
		c = {0}; 
	} else if (p == 3) {
		c = {m_coef[2]}; 
	} else if (p == 6) {
		c = {m_coef[2], m_coef[4], 2*m_coef[5]}; 
	} else if (p == 10) {
		c = {m_coef[2], m_coef[4], 2*m_coef[5], m_coef[7], 2*m_coef[8], 3*m_coef[9]}; 
	} else {
		cout << "ERROR: polynomial order not defined in Poly2D derivative" << endl; 
		exit(0); 
	}

	return Poly2D(c); 

}

