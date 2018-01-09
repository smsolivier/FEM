#include "Fields.H"
#include <iostream>
#include <cassert>

Fields::Fields() {
	m_nUnknowns = 0; 
}

void Fields::set(int a_nodeNum, string a_name) {

	// extend m_points to the required size 
	if (a_nodeNum >= m_points.size()) {
		int size = m_points.size(); 
		for (int i=0; i<a_nodeNum-size+1; i++) {
			m_points.push_back(Point()); 
		}
	}

	// make sure size is correct 
	assert(m_points.size() > a_nodeNum); 

	// add the field to the point 
	Point& p = m_points[a_nodeNum]; 
	if (!p.inPoint(a_name)) {
		p.addField(a_name, m_nUnknowns); 
		m_nUnknowns++; 
	}
}

int Fields::operator()(int a_pointNum, string a_name) {
	return m_points[a_pointNum][a_name]; 
}

const Point& Fields::operator[](int a_pointNum) const {
	if (a_pointNum > m_points.size()) {
		cout << "ERROR (Fields.cpp): index out of m_points" << endl; 
		exit(0); 
	}
	return m_points[a_pointNum]; 
}

int Fields::size() const {return m_nUnknowns; }

vector<double> Fields::extract(vector<double>& a_vals, string a_name) const {
	assert(a_vals.size() == m_nUnknowns); 

	int nPoints = m_points.size(); 
	vector<double> ret(nPoints); 
	for (int i=0; i<nPoints; i++) {
		if (m_points[i].inPoint(a_name)) {
			ret[i] = a_vals[m_points[i][a_name]]; 
		}
	}

	return ret; 
}

void Fields::printUnknownOrdering() const {
	vector<string> ret(m_nUnknowns); 
	for (int i=0; i<m_points.size(); i++) {
		Point p = m_points[i]; 
		for (int j=0; j<m_points[i].numFields(); j++) {
			ret[p[j]] = p.name(j); 
		}
	}

	for (int i=0; i<ret.size(); i++) {
		cout << i << " " << ret[i] << endl; 
	}
}

void Fields::printFieldsPerPoint() const {
	for (int i=0; i<m_points.size(); i++) {
		cout << i << ": "; 
		for (int j=0; j<m_points[i].numFields(); j++) {
			cout << m_points[i].name(j) << " "; 
		}
		cout << endl; 
	}
}