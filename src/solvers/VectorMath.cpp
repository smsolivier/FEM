#include "VectorMath.H"
#include <cassert>
#include <cmath>
#include <iostream>
using namespace std; 

vector<double> operator+(const vector<double> &a, const vector<double> &b) {

	int N = a.size(); 

	// assert same size 
	assert(N == b.size()); 

	vector<double> c(N); 

	for (int i=0; i<N; i++) {

		c[i] = a[i] + b[i]; 

	}

	return c; 

}

vector<double> operator-(const vector<double> &a, const vector<double> &b) {

	int N = a.size(); 

	// assert same size 
	assert(N == b.size()); 

	vector<double> c(N); 

	for (int i=0; i<N; i++) {

		c[i] = a[i] - b[i]; 

	}

	return c; 

}

vector<double> operator*(const vector<double> &a, const double alpha) {

	vector<double> b(a.size()); 

	for (int i=0; i<a.size(); i++) {

		b[i] = a[i] * alpha; 

	}

	return b; 

}

vector<double> operator/(const vector<double>& a, const double alpha) {

	vector<double> b(a.size()); 
	for (int i=0; i<a.size(); i++) {
		b[i] = a[i] / alpha; 
	}

	return b; 
}

void operator/=(vector<double>& a, const double alpha) {

	for (int i=0; i<a.size(); i++) {
		a[i] /= alpha; 
	}
	
}

void operator*=(vector<double>& a, const double alpha) {

	for (int i=0; i<a.size(); i++) {
		a[i] *= alpha; 
	}
	
}

void operator+=(vector<double>& a, const vector<double>& b) {

	for (int i=0; i<a.size(); i++) {
		a[i] += b[i]; 
	}
}

double dot(const vector<double>& a, const vector<double>& b) {

	int N = a.size(); 
	assert(N == b.size()); 

	double sum = 0; 
	for (int i=0; i<N; i++) {

		sum += a[i] * b[i]; 

	}

	return sum; 

}

double L2_norm(const vector<double>& a) {

	double sum = 0; 
	for (int i=0; i<a.size(); i++) {
		sum += a[i]*a[i]; 
	}

	return sqrt(sum); 

}

vector<double> copy(const vector<double>& a) {

	vector<double> b(a.size()); 
	for (int i=0; i<a.size(); i++) {
		b[i] = a[i]; 
	}

	return b; 
}

void printVector(const vector<double> a) {

	for (int i=0; i<a.size(); i++) {
		cout << a[i] << " "; 
	}
	cout << endl; 
}