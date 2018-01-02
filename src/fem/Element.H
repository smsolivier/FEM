#ifndef _ELEMENT_H_
#define _ELEMENT_H_

#include <vector>
#include "Node.H"
#include "Poly2D.H"
#include "GQ.H"
#include "GQ1D.H"

using namespace std; 

/// Stores a collection of nodes and basis functions for an element. 
/** provides integration routines for commonly used combinations of basis functions */ 
class Element {
public: 
	/// default constructor 
	Element() {}; 
	/// constructor 
	Element(vector<Node>& a_nodeList, int a_group=0); 

	/// get a vector of the nodes in the element 
	vector<Node> nodeList() const; 
	/// get the interior node id of node i 
	int interiorNodeID(int i); 
	/// get total node id (includes boundary nodes) of node id 
	int getNodeID(int i); 
	/// get access to node i 
	Node& getNode(int i); 
	/// get the centroid of the element 
	vector<double> getCentroid(); 
	/// get the number of nodes per element 
	int getNumNodes() const; 
	/// get the element's area 
	/** found by integrating the triangle */ 
	double getArea(); 
	/// get the element's polynomial order 
	int getOrder(); 
	/// get the material group of the element 
	int getGroup(); 

	/// computes \f$\int \nabla B_i \cdot \nabla B_j d\Omega \f$
	double gradBiGradBj(int i, int j); 
	/// computes \f$\int B_i B_j d\Omega \f$
	double BiBj(int i, int j); 
	/// computes \f$\int B_i d\Omega \f$ 
	double Bi(int i); 
	/// computes \f$ \hat{n} \oint B_i d\Gamma \f$ on for node a_i
	vector<double> boundaryIntegral(int a_i); 

	/// returns true if face a_i is a boundary face 
	bool boundaryFace(int a_i); 
	/// get outward unit normal for face i 
	vector<double> normal(int a_i); 

private:
	/// computes \f$ \nabla B_i(\xi, \eta) \f$ 
	vector<double> gradient(int i, double xi, double eta); 
	/// computes the \f$(i,j)\f$ entry of the Jacobian of the transformation 
	double J(int i, int j, double xi, double eta); 
	/// computes the \f$(i,j)\f$ entry of the inverse of the Jacobian 
	double Jinv(int i, int j, double xi, double eta); 
	/// computes the determinant of the transformation 
	double detJ(double xi, double eta); 

	/// stores basis functions
	vector<Poly2D> m_basis;  
	/// stores derivative of  basis functions
	vector<vector<Poly2D>> m_dbasis;  
	/// vector of nodes
	vector<Node> m_nodeList;  
	/// number of nodes per element
	int m_nNodes;  
	/// location of centroid
	vector<double> m_centroid;  
	/// polynomial order of basis functions 
	int m_pOrder; 
	/// Gauss Quadrature object 
	GQ m_gq; 
	/// 1D GQ object for boundary integrals 
	GQ1D m_gq1d; 
	/// store which material group 
	int m_group; 

	/// store \f$\int B_i d\Omega \f$ 
	vector<double> m_Bi; 
	/// store \f$\int B_i B_j d\Omega \f$ 
	vector<vector<double>> m_BiBj; 
	/// store \f$\int \nabla B_i \cdot \nabla B_j d\Omega \f$ 
	vector<vector<double>> m_gradBiGradBj; 

	/// true if face is a boundary face 
	vector<bool> m_boundaryFaces; 
	/// lengths of each face 
	vector<double> m_h; 

}; 

#endif