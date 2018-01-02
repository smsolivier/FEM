#ifndef _GQ_H_
#define	_GQ_H_

#include <vector>
#include "CH_Timer.H"
using namespace std; 

/// performs Gauss Quadrature integration on a reference triangle 
class GQ {
public:
	/// default constructor 
	GQ() {}; 
	/// constructor 
	/** a_p: integrate up to a_p order polynomials exactly */ 
	GQ(int a_p); 
	/// integrate the lambda function f 
	double integrate(auto& f) {
		CH_TIMERS("GQ integrate"); 
		double val = 0; 
		for (int i=0; i<m_Np; i++) {
			val += m_w[i] * f(m_xy[i][0], m_xy[i][1]); 
		}	

		return val/2; 
	}

private:
	/// store the quadrature weights 
	vector<double> m_w; 
	/// store the quadrature abscissae 
	vector<vector<double>> m_xy; 
	/// number of quad points 
	int m_Np; 
	/// polynomial order that can be integrated exactly 
	int m_p; 
}; 

#endif