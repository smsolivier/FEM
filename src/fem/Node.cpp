#include "Node.H"
#include <iostream>

Node::Node(vector<double> a_position, const int a_intID, const int a_ID, 
	const bool a_isInterior) {

	m_position = a_position; 
	m_intID = a_intID; 
	m_gID = a_ID; 
	m_interior = a_isInterior; 

}

void Node::printPosition() const {
	for (int i=0; i<m_position.size(); i++) {
		cout << m_position[i] << " "; 
	}
	cout << endl; 
}