#ifndef _GQ1D_H_
#define _GQ1D_H_

#include <vector> 
#include "CH_Timer.H"

using namespace std; 

/// Gauss Quadrature between 0 and 1 
class GQ1D {
public:
	/// default constructor 
	GQ1D() {}; 
	/// constructor 
	/** a_p: integrate polynomials up to a_p order exactly */ 
	GQ1D(int a_p); 
	/// integrate the lambda function f between 0 and 1 
	double integrate(auto& f) {
		// CH_TIMERS("GQ1D integrate"); 
		double val = 0; 
		for (int i=0; i<m_Np; i++) {
			val += m_w[i] * f(m_x[i]); 
		}

		return val; 
	}

	double operator()(auto& f) {return integrate(f); } 

private:
	/// quadrature weights 
	vector<double> m_w; 
	/// quadrature points 
	vector<double> m_x; 
	/// number of quad points 
	int m_Np; 
	/// polynomial order that can be integrated exactly 
	int m_p; 
}; 

#endif 