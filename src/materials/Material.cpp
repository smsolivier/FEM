#include "Material.H"
#include <iostream>

Material::Material(string a_name) {

	m_name = a_name; 
	m_zero = 0.; 

}

double& Material::operator()(string a_var) {

	// search for a_var in m_vars 
	int id = -1; 
	for (int i=0; i<m_vars.size(); i++) {
		if (a_var.compare(m_vars[i]) == 0) id = i; 
	}

	if (id == -1) {
		m_vars.push_back(a_var); 
		m_vals.push_back(m_zero); 
		return m_vals[m_vals.size() - 1]; 
	} else {
		return m_vals[id]; 
	}
} 

double Material::operator()(string a_var) const {

	// search for a_var in m_vars 
	int id = -1; 
	for (int i=0; i<m_vars.size(); i++) {
		if (a_var.compare(m_vars[i]) == 0) id = i; 
	}

	if (id == -1) {
		cout << "ERROR: variable " << a_var << 
			" not found in Material " << m_name << endl; 
		exit(0); 
	} else {
		return m_vals[id]; 
	}

}

string Material::getName() const {return m_name; }
