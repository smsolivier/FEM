#ifndef _BASIS_H_
#define _BASIS_H_

#include <vector>
#include "Node.H"
#include "Poly2D.H"

using namespace std; 

/// Evaluate basis functions 
class Basis {
public:
	/// default constructor
	Basis(); 
	/// constructor 
	Basis(vector<Node>& a_nodeList); 
	/// get the basis function for node a_i
	inline Poly2D& operator[](int a_i); 
	/// get the a_jth component of the a_ith basis function 
	inline Poly2D& operator()(int a_i, int a_j); 
	/// computes \f$ \nabla B_i(\xi, \eta) \f$ 
	inline vector<double> gradient(int i, double xi, double eta);
	/// computes \f$ \frac{\partial B_i}{\partial x} \f$ 
	inline double dx(int i, double xi, double eta); 
	/// computes \f$ \frac{\partial B_i}{\partial y} \f$ 
	inline double dy(int i, double xi, double eta);  
	/// computes the \f$(i,j)\f$ entry of the Jacobian of the transformation 
	inline double J(int i, int j, double xi, double eta); 
	/// computes the \f$(i,j)\f$ entry of the inverse of the Jacobian 
	inline double Jinv(int i, int j, double xi, double eta); 
	/// computes the determinant of the transformation 
	inline double detJ(double xi, double eta); 
private:
	/// store the nodes 
	vector<Node> m_nodeList; 
	/// basis functions 
	vector<Poly2D> m_basis; 
	/// gradients of basis function 
	vector<vector<Poly2D>> m_dbasis; 
	/// number of nodes 
	int m_nNodes; 
}; 
#include "BasisImplem.H"

#endif