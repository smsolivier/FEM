#ifndef _ELEMENT_H_
#define _ELEMENT_H_

#include <vector>
#include "Basis.H"
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
	/// access to node i by indexing into m_nodeList 
	/** can call any node member function through:\n 
		<code> Element[i].function() </code> 
	*/ 
	Node& operator[](int i); 
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
	/// computes \f$\int 2\frac{\partial B_i}{\partial x} \frac{\partial B_j}{\partial x} + \frac{\partial B_i}{\partial y} \frac{\partial B_j}{\partial y} d\Omega \f$ 
	double NSx(int i, int j); 
	/// computes \f$\int 2\frac{\partial B_i}{\partial y} \frac{\partial B_j}{\partial y} + \frac{\partial B_i}{\partial x} \frac{\partial B_j}{\partial x} d\Omega \f$
	double NSy(int i, int j); 
	/// computes \f$\int \frac{\partial B_i}{\partial y} \frac{\partial B_j}{\partial x} d\Omega \f$ 
	double dBidy_dBjdx(int i, int j); 
	/// computes \f$\int \frac{\partial B_i}{\partial x} \frac{\partial B_j}{\partial y} d\Omega \f$ 
	double dBidx_dBjdy(int i, int j); 
	/// computes \f$\int \frac{\partial B_i}{\partial x}\hat{B_j} d\Omega \f$
	double dBidxBhatj(int i, int j); 
	/// computes \f$\int \frac{\partial B_i}{\partial y}\hat{B_j} d\Omega \f$ 
	double dBidyBhatj(int i, int j); 

	/// returns true if face a_i is a boundary face 
	bool boundaryFace(int a_i); 
	/// get outward unit normal for face i 
	vector<double> normal(int a_i); 

private:
	/// quadratic velocity basis functions 
	Basis m_vbasis; 
	/// linear pressure basis functions 
	Basis m_pbasis; 
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
	/// store \f$\int 2\frac{\partial B_i}{\partial x} \frac{\partial B_j}{\partial x} + \frac{\partial B_i}{\partial y} \frac{\partial B_j}{\partial y} d\Omega \f$ 
	vector<vector<double>> m_NSx; 
	/// store \f$\int 2\frac{\partial B_i}{\partial y} \frac{\partial B_j}{\partial y} + \frac{\partial B_i}{\partial x} \frac{\partial B_j}{\partial x} d\Omega \f$
	vector<vector<double>> m_NSy; 
	/// store \f$\int \frac{\partial B_i}{\partial y} \frac{\partial B_j}{\partial x} d\Omega \f$
	vector<vector<double>> m_dBidy_dBjdx; 
	/// stores \f$\int \frac{\partial B_i}{\partial x} \frac{\partial B_j}{\partial y} d\Omega \f$
	vector<vector<double>> m_dBidx_dBjdy; 
	/// store \f$\int \frac{\partial B_i}{\partial x}\hat{B_j} d\Omega \f$
	vector<vector<double>> m_dBidxBhatj; 
	/// store \f$\int \frac{\partial B_i}{\partial y}\hat{B_j} d\Omega \f$
	vector<vector<double>> m_dBidyBhatj; 

	// first time booleans 
	bool m_fNSx; 
	bool m_fNSy; 
	bool m_fdBidy_dBjdx; 
	bool m_fdBidx_dBjdy; 
	bool m_fdBidxBhatj; 
	bool m_fdBidyBhatj; 

	/// true if face is a boundary face 
	vector<bool> m_boundaryFaces; 
	/// lengths of each face 
	vector<double> m_h; 

}; 

#endif