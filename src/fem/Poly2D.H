#ifndef _POLY2D_H_
#define _POLY2D_H_

#include <vector>

using namespace std; 

/// Stores a complete 2D polynomial\n Coefficients are stored in order of Pascal's triangle 
class Poly2D {
public: 
	/// default constructor 
	Poly2D(); 
	/// Constructor 
	/** provide a vector of the coefficients */ 
	Poly2D(vector<double> coef); 
	/// Evaluates the polynomial at \f$(x,y)\f$ 
	double operator()(double x, double y); 

	/// returns the gradient of the polynomial 
	vector<Poly2D> gradient(); 

private:
	/// compute x derivative 
	Poly2D dx(); 
	/// computes y derivative 
	Poly2D dy(); 

	/// stores the coefficients 
	vector<double> m_coef; 
}; 

#endif 