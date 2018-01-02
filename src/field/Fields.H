#ifndef _FIELDS_H_
#define _FIELDS_H_
#include "Point.H"
#include <vector>
#include <string>

using namespace std; 

/// Maps a node number and variable name to a global index 
class Fields {
public:
	/// constructor 
	Fields(); 
	/// assigns a global index to the variable a_name at point a_pointNum
	void set(int a_pointNum, string a_name); 
	/// get global index from a point number and variable name 
	int operator()(int a_pointNum, string a_name); 
	/// access to the Point object 
	/** can get the global index for a variable with:\n
		<code> int gid = Fields[a_pointNum][a_name]; </code>\n 
		This syntax is equivalent to calling operator()
	*/ 
	Point& operator[](int a_pointNum); 
	/// get number of unknowns 
	int size(); 
	/// pull out a_name values from a_vals 
	vector<double> extract(vector<double>& a_vals, string a_name);
	/// print out the ordering of the variables in the global vector 
	void printUnknownOrdering();  
	/// print the fields per point 
	void printFieldsPerPoint(); 
private:
	/// store the solution points 
	vector<Point> m_points; 
	/// store the number of unknowns
	int m_nUnknowns; 
}; 

#endif