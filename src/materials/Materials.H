#ifndef _MATERIALS_H_
#define _MATERIALS_H_

#include <string> 
#include <vector>
#include "Material.H"

using namespace std; 

/// Stores a collection of Material objects 
class Materials {
public:
	/// default constructor 
	Materials(); 
	/// access and change a material variable 
	/** if material a_name does not exist, a_name is added and a_var is initialized to zero 
	*/ 
	double& operator()(string a_name, string a_var);  
	/// access only to variable a_var from material a_name 
	double operator()(string a_name, string a_var) const; 
	/// access only to a_var from the a_id'th material defined 
	/** requires materials to be added in order of use */ 
	const Material& operator[](int a_id) const; 
private:
	vector<Material> m_materials; 
	double m_zero; 
}; 

#endif 