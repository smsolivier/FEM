#ifndef _FEGRID_H_ 
#define _FEGRID_H_

#include <string>
#include <vector> 

#include "Node.H"
#include "Element.H"
#include "Fields.H"

using namespace std; 

/// Loads the mesh and initializes elements and nodes 
class FEGrid {
public:
	/// Loads .ele and .node file and initializes elements and nodes 
	FEGrid(const string& a_mshFile); 
	/// get access to node i. i is the global index number 
	const Node node(int i) const; 
	/// get total number of nodes (includes boundary nodes) 
	int getNumNodes() const; 
	/// get number of interior nodes 
	int getNumInteriorNodes() const; 
	/// get number of elements 
	int getNumElements() const;
	/// compute the average element area in the mesh  
	double getAverageElementArea(); 
	/// access to element i 
	Element& getElement(int i); 

	/// \f$\int B_i \phi_i d\Omega \f$ 
	double integrateBi(const vector<double>& a_phi); 

	/// write a .vtk file 
	void write(vector<double> a_scalarField, const char* a_filename); 

	/// prints node locations in each element 
	void printNodes(); 

	/// print number of elements and unknowns and the average element area 
	void meshInfo(); 

	/// return field object 
	Fields& getFields(); 
private:
	/// stores the nodes 
	vector<Node> m_nodes;
	/// stores the elements  
	vector<Element> m_elements; 
	/// stores the number of interior nodes 
	int m_nInteriorNodes; 
	/// stores the total number of nodes 
	int m_nTotalNodes; 
	/// store unknown ordering 
	Fields m_fields; 
}; 

#endif 