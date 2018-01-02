#include "Fields.H"
#include <iostream>

using namespace std; 

int main() {

	Fields x; 

	x.set(1, "u"); 
	x.set(1, "v"); 
	x.set(3, "u"); 

	x.printUnknownOrdering(); 
	// x.printFieldsPerPoint(); 

	cout << x[1][1] << endl; 

}