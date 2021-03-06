#ifndef _NODE_H_
#define _NODE_H_

#include <vector>

using namespace std; 

/// Stores information for an FEM node 
class Node {
public:
	/// default constructor 
	Node() {}; 
	/// Constructor 
	/** a_position: x, y position of the node \n
		a_interiorNodeID: global ID of interior nodes only \n
		a_globalID: global ID of nodes including the boundary nodes \n
		a_isInterior: stores whether the node is on the interior
	*/ 
	Node(vector<double> a_position, 
		const int a_interiorNodeID, 
		const int a_globalID, 
		const bool a_isInterior); 

	/// returns the x, y position of the node 
	inline vector<double> getPosition() {return m_position; };  
	/// returns the nodes interior ID
	inline int getInteriorNodeID() {return m_intID; }; 
	/// same as getInteriorNodeID
	inline int interiorNodeID() {return m_intID; }; 
	/// returns the node's global ID 
	inline int getNodeID() {return m_gID; };
	/// returns true if the node is on the interior  
	inline bool isInterior() {return m_interior; }; 

	/// print the x,y,z position of the node to screen 
	void printPosition() const; 
	
private:
	/// store x, y, z position 
	vector<double> m_position; 
	/// interior id 
	int m_intID; 
	/// global id
	int m_gID; 
	/// store if interior or not 
	bool m_interior; 
}; 

#endif 