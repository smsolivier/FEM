#ifndef _VECTORMATH_H_
#define _VECTORMATH_H_

#include <vector>
using namespace std; 

vector<double> operator+(const vector<double> &a, const vector<double> &b); 
vector<double> operator-(const vector<double> &a, const vector<double> &b); 
vector<double> operator*(const vector<double> &a, const double alpha); 
vector<double> operator/(const vector<double> &a, const double alpha); 
void operator/=(vector<double>& a, const double alpha); 
void operator*=(vector<double>& a, const double alpha); 
void operator+=(vector<double>& a, const vector<double>& b); 
double dot(const vector<double>& a, const vector<double>& b); 
double L2_norm(const vector<double>& a); 
vector<double> copy(const vector<double>& a); 
void printVector(vector<double> a); 

#endif