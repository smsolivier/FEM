#ifndef _MATERIAL_H_
#define _MATERIAL_H_
#include <string> 
#include <vector>

using namespace std; 

/// interface for storing and accessing material information 
class Material {
public:
	/// default constructor 
	Material() {}; 
	/// constructor 
	Material(string a_name); 
	/// access a variable name a_var 
	/** can fetch and set values if *this is not const\n
		adds zero value if a_var is not already present
	 */ 
	double& operator()(string a_var); 
	/// retrieve variable a_var 
	/** can only retrieve if *this is const \n
		throws an error if a_var is not in m_vars 
	*/ 
	double operator()(string a_var) const; 
	/// return the name of the material 
	string getName() const; 
private:
	/// store the name of the material 
	string m_name; 
	/// store the variable names 
	vector<string> m_vars; 
	/// store the corresponding values 
	vector<double> m_vals; 	
	/// zero for returning in operator() 
	double m_zero; 
}; 

#endif 