#include "Materials.H"
#include <iostream>

Materials::Materials() {
	m_zero = 0; 
}

double& Materials::operator()(string a_name, string a_var) {

	// search for material name in m_materials 
	int id = -1; 
	for (int i=0; i<m_materials.size(); i++) {
		if (a_name.compare(m_materials[i].getName()) == 0) id = i; 
	}

	// if not found initialize a Material object and add to m_materials 
	if (id == -1) {
		m_materials.push_back(Material(a_name)); 
		id = m_materials.size() - 1; 
	}

	// return the variable 
	return m_materials[id](a_var); 
}

double Materials::operator()(string a_name, string a_var) const {

	// search for material name in m_materials 
	int id = -1; 
	for (int i=0; i<m_materials.size(); i++) {
		if (a_name.compare(m_materials[i].getName()) == 0) {
			id = i; 
		}
	}

	if (id == -1) {
		cout << "ERROR: material " << a_name 
			<< " not found" << endl; 
		exit(0); 
	} 

	return m_materials[id](a_var); 
}

// double Materials::operator()(int a_id, string a_var) const {

// 	// check if index out of bounds 
// 	if (a_id < 0) {
// 		cout << "ERROR: material id < 0" << endl; 
// 		exit(0); 
// 	} else if (a_id >= m_materials.size()) {
// 		cout << "ERROR: materials index out of bounds" << endl; 
// 		exit(0); 
// 	}

// 	return m_materials[a_id](a_var); 
// }

const Material& Materials::operator[](int a_id) const {

	if (a_id >= m_materials.size()) {
		cout << "ERROR (Materials.cpp): index out of bounds in operator[]" << endl; 
		exit(0); 
	}

	return m_materials[a_id]; 

}