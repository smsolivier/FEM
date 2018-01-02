#include "FEGrid.H"
#include "VisitWriter.H"
#include "CH_Timer.H"
#include <fstream>
#include <iostream>
#include <sstream>
#include <cassert>
#include <cmath>

FEGrid::FEGrid(const string& a_mshFile) {

	ifstream nodes(a_mshFile+".node"); 
	int n_nodes; 
	nodes >> n_nodes; 

	m_nodes.resize(n_nodes); 
	m_nInteriorNodes = 0; 
	m_nTotalNodes = 0; 
	vector<double> position(2);
	double aa;  
	int ID, boundary; 
	for (int i=0; i<n_nodes; i++) {
		nodes >> ID >> boundary >> position[0] >> position[1] >> aa;
		if (boundary == 0) {
			m_nodes[ID] = Node(position, m_nInteriorNodes, m_nTotalNodes, true); 
			m_nInteriorNodes++; 
		} else if (boundary == 2) {
			m_nodes[ID] = Node(position, m_nInteriorNodes, m_nTotalNodes, false); 
			m_nInteriorNodes++; 
		} else {
			m_nodes[ID] = Node(position, -1, m_nTotalNodes, false); 
		}
		m_nTotalNodes++; 
	}

	nodes.close(); 

	int max_group = 0; 
	int min_group = 10000; 
	ifstream elements(a_mshFile+".ele"); 
	int n_el; 
	elements >> n_el; 

	m_elements.resize(n_el); 
	string line; 
	getline(elements, line); 
	int val; 
	for (int i=0; i<n_el; i++) {
		getline(elements, line); 
		stringstream ss(line); 
		int tmp, group; 
		ss >> tmp >> group; 
		if (group > max_group) max_group = group; 
		if (group < min_group) min_group = group; 
		vector<Node> node_list; 
		while (ss >> val) {
			node_list.push_back(m_nodes[val]); 
		}
		m_elements[i] = Element(node_list, group); 
	}

	if (min_group < 0) {
		cout << "WARNING: material group less than zero" << endl; 
	}

	// setup field object for P2P1 ordering of unknowns 
	for (int i=0; i<m_elements.size(); i++) {
		Element& el = m_elements[i]; 
		for (int j=0; j<el.getNumNodes(); j++) {
			int jid = el.interiorNodeID(j); 
			if (jid >= 0) {
				m_fields.set(jid, "u"); 
				m_fields.set(jid, "v"); 
				if (j < 3) {
					m_fields.set(jid, "p"); 
				}
			}
		}
	}
}

const Node FEGrid::node(int i) const {return m_nodes[i]; }

int FEGrid::getNumNodes() const { return m_nTotalNodes; }
int FEGrid::getNumInteriorNodes() const { return m_nInteriorNodes; }
int FEGrid::getNumElements() const {return m_elements.size(); }
double FEGrid::getAverageElementArea() {
	double area = 0; 
	for (int i=0; i<m_elements.size(); i++) {
		area += m_elements[i].getArea(); 
	}

	return area/m_elements.size(); 
}

Element& FEGrid::getElement(int i) {

	if (i > m_elements.size()) {
		cout << "ERROR: indexed past m_elements" << endl; 
		exit(0); 
	}
	return m_elements[i]; 
}

double FEGrid::integrateBi(const vector<double>& a_phi) {

	double sum = 0; 
	int nEl = this->getNumElements(); 

	for (int i=0; i<nEl; i++) {
		Element& el = this->getElement(i); 
		for (int j=0; j<el.getNumNodes(); j++) {
			int jid = el.interiorNodeID(j); 
			if (jid >= 0) {
				sum += el.Bi(j) * a_phi[jid]; 
			}
		}
	}

	return sum; 

}

void FEGrid::printNodes() {

	for (int i=0; i<m_elements.size(); i++) {
		cout << i << " "; 
		vector<Node> nodes = m_elements[i].nodeList(); 
		for (int j=0; j<nodes.size(); j++) {
			cout << nodes[j].getNodeID() << " "; 
		}
		cout << endl; 
	}
}

void FEGrid::meshInfo() {

	cout << "Number of Elements = " << m_elements.size() << endl; 
	cout << "Number of Unknowns = " << m_nInteriorNodes << endl; 
	cout << "Average Element Area = " << getAverageElementArea() << endl; 
	cout << "Average Element Length = " << sqrt(getAverageElementArea()) << endl; 
	cout << "Element Order = " << m_elements[0].getOrder() << endl; 
}

Fields& FEGrid::getFields() {return m_fields; }

void FEGrid::write(vector<double> a_scalarField, const char* a_filename) {

	CH_TIMERS("write to vtk"); 

	// insert zero for boundary nodes 
	vector<double> tmp(m_nTotalNodes); 
	for (int i=0; i<m_nTotalNodes; i++) {
		if (m_nodes[i].getInteriorNodeID() >= 0) {
			tmp[i] = a_scalarField[m_nodes[i].getInteriorNodeID()]; 
		}
		else tmp[i] = 0.; 
	}

	int nNodes = tmp.size(); 

	assert(m_nodes.size() == nNodes);
	vector<float> scalarField(nNodes); 
	for (int i=0; i<nNodes; i++) {
		scalarField[i] = tmp[i]; 
	}

	vector<float> pts(3*nNodes); 
	for (int i=0; i<nNodes; i++) {
		int p = 3*i; 
		vector<double> position = m_nodes[i].getPosition(); 
		pts[p] = position[0]; 
		pts[p+1] = position[1]; 
		pts[p+2] = position[2]; 
	}

	float* vars[1];
	vars[0] = &(scalarField[0]);
	int vardim[1] = {1};
	int centering[1] = {1};
	const char * const varnames[] = { "nodeData" };

	int ncell = m_elements.size(); 
	if (m_elements[0].getNumNodes() == 6) ncell *= 4; 
	else if (m_elements[0].getNumNodes() == 10) ncell *= 9; 
	vector<int> cellType;  
	vector<int> conns; 
	int ii = 0; 
	for (int i=0; i<m_elements.size(); i++) {
		Element& el = m_elements[i]; 
		vector<Node> nodes = el.nodeList(); 
		int nNodes = el.getNumNodes(); 

		if (nNodes == 3) {
			conns.push_back(nodes[0].getNodeID()); 
			conns.push_back(nodes[1].getNodeID()); 
			conns.push_back(nodes[2].getNodeID()); 
			cellType.push_back(VISIT_TRIANGLE); 
		}

		else if (nNodes == 6) {
			conns.push_back(nodes[0].getNodeID()); 
			conns.push_back(nodes[3].getNodeID()); 
			conns.push_back(nodes[5].getNodeID()); 
			cellType.push_back(VISIT_TRIANGLE); 

			conns.push_back(nodes[3].getNodeID()); 
			conns.push_back(nodes[4].getNodeID()); 
			conns.push_back(nodes[5].getNodeID()); 
			cellType.push_back(VISIT_TRIANGLE); 

			conns.push_back(nodes[3].getNodeID()); 
			conns.push_back(nodes[1].getNodeID()); 
			conns.push_back(nodes[4].getNodeID()); 
			cellType.push_back(VISIT_TRIANGLE); 

			conns.push_back(nodes[5].getNodeID()); 
			conns.push_back(nodes[4].getNodeID()); 
			conns.push_back(nodes[2].getNodeID()); 
			cellType.push_back(VISIT_TRIANGLE); 
		}

		else if (nNodes == 10) {
			conns.push_back(nodes[0].getNodeID()); 
			conns.push_back(nodes[3].getNodeID()); 
			conns.push_back(nodes[8].getNodeID()); 
			cellType.push_back(VISIT_TRIANGLE); 

			conns.push_back(nodes[3].getNodeID()); 
			conns.push_back(nodes[4].getNodeID()); 
			conns.push_back(nodes[9].getNodeID()); 
			cellType.push_back(VISIT_TRIANGLE); 

			conns.push_back(nodes[4].getNodeID()); 
			conns.push_back(nodes[1].getNodeID()); 
			conns.push_back(nodes[5].getNodeID()); 
			cellType.push_back(VISIT_TRIANGLE); 
			
			conns.push_back(nodes[5].getNodeID()); 
			conns.push_back(nodes[9].getNodeID()); 
			conns.push_back(nodes[4].getNodeID()); 
			cellType.push_back(VISIT_TRIANGLE); 

			conns.push_back(nodes[9].getNodeID()); 
			conns.push_back(nodes[8].getNodeID()); 
			conns.push_back(nodes[3].getNodeID()); 
			cellType.push_back(VISIT_TRIANGLE); 

			conns.push_back(nodes[9].getNodeID()); 
			conns.push_back(nodes[5].getNodeID()); 
			conns.push_back(nodes[6].getNodeID()); 
			cellType.push_back(VISIT_TRIANGLE); 

			conns.push_back(nodes[8].getNodeID()); 
			conns.push_back(nodes[9].getNodeID()); 
			conns.push_back(nodes[7].getNodeID()); 
			cellType.push_back(VISIT_TRIANGLE); 

			conns.push_back(nodes[9].getNodeID()); 
			conns.push_back(nodes[6].getNodeID()); 
			conns.push_back(nodes[7].getNodeID()); 
			cellType.push_back(VISIT_TRIANGLE); 

			conns.push_back(nodes[7].getNodeID()); 
			conns.push_back(nodes[6].getNodeID()); 
			conns.push_back(nodes[2].getNodeID()); 
			cellType.push_back(VISIT_TRIANGLE); 
		}
	}

	write_unstructured_mesh(a_filename, 0, nNodes, &(pts[0]), cellType.size(), 
		&(cellType[0]), &(conns[0]), 1, vardim, centering, varnames, vars); 

}