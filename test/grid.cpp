#include "FEGrid.H"

using namespace std; 
int main() {
	FEGrid grid("../mesh/hole"); 
	grid.meshInfo(); 

	Fields& f = grid.getFields(); 
	f.printUnknownOrdering(); 
	// f.printFieldsPerPoint(); 
}