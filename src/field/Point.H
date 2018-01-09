#ifndef _POINT_H_
#define _POINT_H_
#include <string>
#include <vector> 

using namespace std; 

/// stores variable names and global indices for a solution point 
class Point {
public:
	/// constructor 
	Point(); 
	/// add the field name a_name with global index a_globalIndex
	void addField(string a_name, int a_globalIndex); 
	/// get the global index corresponding to a_name 
	int operator[](string a_name) const; 
	/// get the global index corresponding to the a_ith entry in m_ind
	int operator[](int a_i) const; 
	/// check if a variable is in m_names 
	bool inPoint(string a_name) const; 
	/// get number of fields 
	int numFields() const; 
	/// get name of variable a_i 
	string name(int a_i) const; 
private:
	/// store global indices corresponding to m_names 
	vector<int> m_ind; 
	/// store the field names 
	vector<string> m_names; 
}; 

#endif