#ifndef _SPARSEMATRIX_H_
#define _SPARSEMATRIX_H_

#ifdef SUPERLU
#include "slu_ddefs.h"
#endif
#include "Sparse"
#include <vector> 
#include <string> 

using namespace std; 

/// stores a matrix in a 2D compressed column format
class SparseMatrix 
{
public:
	/// default constructor
	SparseMatrix() {}; 
	/// construct an a_M by a_N sparse matrix 
	SparseMatrix(int a_M, int a_N); 

	/// matrix vector product 
	vector<double> operator*(const vector<double>& a_v) const; 
	/// multipy all non-zero entrys by a double 
	void scale(const double alpha); 
	/// add elements to the matrix
	/**
		returns a reference to the entry 
		if the entry does not exist one is created and assigned a value of zero 
		only use if adding non-zeros in otherwise it will clutter the sparsematrix with zeros
	*/
	double& operator()(int a_row, int a_col);
	/// get values from the matrix 
	/** 
		if a value does not exist zero is returned 
		a new element is not created 
		does not fill the matrix in if the index does not exist like operator()
	*/ 
	double at(int a_row, int a_col) const; // extract values from the matrix 

	/// set the diagonal to 1 
	void eye(); // create identity matrix 

	/// prints the dense matrix to screen 
	void print(); 
	/// prints the dense matrix to file 
	void print(string a_outname); 
	/// print m_rowIndex to screen  
	void printRowIndex();
	/// output sparsity information to screen 
	void sparsity() const;

	/// converts 2D to 1D compressed column format and creates a SuperMatrix (SuperLU formatted sparse matrix)
#ifdef SUPERLU
	const SuperMatrix getSuperMatrix() const; 
#endif
	/// creates a dense matrix 
	const vector<vector<double>> getDenseMatrix() const; 
	/// convert to Eigen matrix format 
	const Eigen::SparseMatrix<double> getEigenMatrix() const; 
	/// convert to triplet format 
	void convertToTriplet(vector<double>& a_vals, 
		vector<int>& a_rows, vector<int>& a_cols) const; 

	/// return number of rows 
	int getM() const; 
	/// returns number of columns 
	int getN() const;
	/// return number of nonzeros 
	int getNNZ() const; 

	/// returns a SparseMatrix that is the transpose of *this
	SparseMatrix transpose() const;
	/// tests if *this is symmetric 
	bool symmetric() const;  
	/// tests for positive diagonal 
	bool positiveDiagonal() const; 

private:
	/// store number of rows, cols 
	int m_m, m_n;
	/// store number of nonzeros
	int m_nnz;  
	/// store a zero to be returned as reference 
	double m_zero; 
	/// stores the non-zero data 
	vector<vector<double>> m_data; 
	/// stores where the nonzeros go in the matrix 
	vector<vector<int>> m_rowIndex; 
};

#endif 