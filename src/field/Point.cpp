#include "Point.H"
#include <iostream>

Point::Point() {}

void Point::addField(string a_name, int a_globalIndex) {
	if (!inPoint(a_name)) {
		m_names.push_back(a_name); 
		m_ind.push_back(a_globalIndex); 
	}
}

int Point::operator[](string a_name) {
	int id = -1; 
	for (int i=0; i<m_names.size(); i++) {
		if (m_names[i].compare(a_name) == 0) {
			id = i; 
		}
	}

	if (id == -1) {
		cout << "ERROR: field " << a_name << " not in Point" << endl; 
		exit(0); 
	}

	return m_ind[id]; 
}

int Point::operator[](int a_i) {

	if (a_i >= m_ind.size()) {
		cout << "ERROR (Point.cpp): field index out of bounds in operator[]" << endl; 
		exit(0); 
	}
	return m_ind[a_i];

}

bool Point::inPoint(string a_name) {

	for (int i=0; i<m_names.size(); i++) {
		if (m_names[i].compare(a_name) == 0) return true; 
	}

	return false; 

}

int Point::numFields() {return m_names.size(); }
string Point::name(int a_i) {return m_names[a_i]; }