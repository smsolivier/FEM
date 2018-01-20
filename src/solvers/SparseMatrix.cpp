#include "SparseMatrix.H"
#include <fstream>
#include <iostream>

#include "CH_Timer.H"

SparseMatrix::SparseMatrix(int a_M, int a_N) {

	// store matrix size 
	m_m = a_M; 
	m_n = a_N; 

	// set zero 
	m_zero = 0.0; 

	// set initial size to number of columns 
	m_rowIndex.resize(m_n); 
	m_data.resize(m_n); 

	// initialize m_nnz 
	m_nnz = 0; 

}

vector<double> SparseMatrix::operator*(const vector<double>& a_v) const {

	CH_TIMERS("MVP"); 

	// check sizing 
	int N = a_v.size(); 
	if (N != m_m) {
		cout << "ERROR: sizes don't agree in SparseMatrix * vector<double>" << endl; 
		exit(0); 
	}

	vector<double> sol(m_m, 0); 

	for (int i=0; i<m_n; i++) {

		for (int j=0; j<m_rowIndex[i].size(); j++) {

			int col = m_rowIndex[i][j]; 

			sol[col] += m_data[i][j] * a_v[i]; 

		}

	}

	return sol; 

}

void SparseMatrix::scale(const double alpha) {

	for (int i=0; i<m_n; i++) {
		for (int j=0; j<m_rowIndex[i].size(); j++) {
			m_data[i][j] *= alpha; 
		}
	}
}

double& SparseMatrix::operator()(int a_row, int a_col) {

	if (a_row >= m_m || a_col >= m_n) {
		cout << "ERROR: sparse matrix indices out of bounds" << endl; 
		exit(0); 
	}

	// extract the row information for column a_col
	vector<int> row = m_rowIndex[a_col]; 

	// search if (a_row, a_col) already stored 
	int index = -1; 
	for (int i=0; i<row.size(); i++) {

		if (row[i] == a_row) index = i; 

	}

	// if not found, add into m_data 
	if (index == -1) {

		// make it zero in data 
		m_data[a_col].push_back(m_zero); 

		// add the index in to rowIndex 
		m_rowIndex[a_col].push_back(a_row); 

		// position into m_data[a_col] 
		index = row.size(); 

		m_nnz += 1; 

	}

	return m_data[a_col][index]; 

}

double SparseMatrix::at(int a_row, int a_col) const {

	if (a_row >= m_m || a_col >= m_n) {
		cout << "ERROR: sparse matrix indices out of bounding" << endl; 
		exit(0); 
	}

	// extract the row information from column a_col 
	vector<int> row = m_rowIndex[a_col]; 

	// search for a_row 
	int index = -1; 
	for (int i=0; i<row.size(); i++) {

		if (row[i] == a_row) {

			// cout << "found " << a_row << " " << a_col << endl; 
			return m_data[a_col][i]; 

		}

	}

	return m_zero; 

}

#if defined SUPERLU || defined PSUPERLU
SuperMatrix SparseMatrix::getSuperMatrix() const {

	CH_TIMERS("SPMAT: convert to SuperMatrix"); 

	double * val = new double[m_nnz]; 
	int * row = new int[m_nnz]; 
	int * colptr = new int[m_n+1]; 

	colptr[0] = 0; 

	int index = 0; 
	for (int i=0; i<m_n; i++) {

		int nvals = m_rowIndex[i].size(); // number of values in column i 

		for (int j=0; j<m_rowIndex[i].size(); j++) {

			// copy row over 
			row[index] = m_rowIndex[i][j]; 

			// copy data over
			val[index] = m_data[i][j]; 

			// update index 
			index++; 

		}

		// update colptr 
		colptr[i+1] = colptr[i] + nvals; 

	}

	SuperMatrix A; 

	dCreate_CompCol_Matrix(&A, m_m, m_n, m_nnz, val, row, colptr, 
		SLU_NC, SLU_D, SLU_GE); 

	return A; 

}
#endif

const vector<vector<double>> SparseMatrix::getDenseMatrix() const {

	CH_TIMERS("SPMAT: convert to dense format"); 

	vector<vector<double>> A(m_m, vector<double>(m_n)); 

	for (int i=0; i<m_m; i++) {

		for (int j=0; j<m_n; j++) {

			A[i][j] = (*this).at(i,j); 

		}
	}

	return A; 

}

const Eigen::SparseMatrix<double> SparseMatrix::getEigenMatrix() const {

	CH_TIMERS("SpMat: convert to Eigen"); 

	Eigen::SparseMatrix<double> eigen(m_m, m_n); 

	// transfer values over 
	for (int i=0; i<m_m; i++) {
		for (int j=0; j<m_rowIndex[i].size(); j++) {
			eigen.insert(m_rowIndex[i][j], i) = m_data[i][j]; 
		}
	}

	return eigen; 

}

void SparseMatrix::convertToTriplet(vector<double>& a_vals, 
	vector<int>& a_rows, vector<int>& a_cols) const {

	CH_TIMERS("SpMat: convert to triplet"); 

	// make sure vectors are correct size 
	a_vals.resize(m_nnz); 
	a_rows.resize(m_nnz); 
	a_cols.resize(m_nnz); 

	int index = 0; 
	for (int i=0; i<m_m; i++) {
		for (int j=0; j<m_rowIndex[i].size(); j++) {
			a_vals[index] = m_data[i][j]; 
			a_rows[index] = m_rowIndex[i][j]; 
			a_cols[index] = i; 
			index++; 
		}
	}

}

int SparseMatrix::getM() const { return m_m; }
int SparseMatrix::getN() const { return m_n; } 
int SparseMatrix::getNNZ() const {return m_nnz; }

SparseMatrix SparseMatrix::transpose() const {

	// create new sparse matrix 
	SparseMatrix sp(m_m, m_n); 

	for (int i=0; i<m_n; i++) {

		for (int j=0; j<m_rowIndex[i].size(); j++) {

			sp(i, m_rowIndex[i][j]) = m_data[i][j]; 

		}

	}

	return sp; 

}

bool SparseMatrix::symmetric() const {

	CH_TIMERS("SPMAT: test matrix symmetry"); 

	const SparseMatrix T = this->transpose(); 

	// check transpose == this 
	for (int i=0; i<m_n; i++) {

		for (int j=0; j<m_rowIndex[i].size(); j++) {

			if (T.at(m_rowIndex[i][j],i) != m_data[i][j]) return false; 

		}

	}

	return true; 

}

bool SparseMatrix::positiveDiagonal() const {

	CH_TIMERS("SPMAT: check for positive diagonal"); 

	for (int i=0; i<m_n; i++) {
		for (int j=0; j<m_rowIndex[i].size(); j++) {
			if (m_rowIndex[i][j] == i) {
				if (m_data[i][j] <= 0) {
					cout << "SparseMatrix(" << m_rowIndex[i][j] << "," << i
						<< ") = " << m_data[i][j] << endl; 
					return false; 
				}
			}
		}
	}

	return true; 

}

void SparseMatrix::eye() {

	for (int i=0; i<m_n; i++) {
		(*this)(i,i) = 1; 
	}

}

void SparseMatrix::print() {

	CH_TIMERS("SPMAT: print to screen"); 

	for (int i=0; i<m_m; i++) {

		for (int j=0; j<m_n; j++) {

			cout << (*this).at(i,j) << " "; 

		}

		cout << endl; 

	}

}

void SparseMatrix::print(string a_outname) {

	CH_TIMERS("SPMAT print to file"); 

	ofstream out; 
	out.open(a_outname); 

	for (int i=0; i<m_m; i++) {

		for (int j=0; j<m_n; j++) {

			out << (*this).at(i,j) << " "; 

		}

		out << endl; 

	}

	out.close(); 

}

void SparseMatrix::printRowIndex() {

	for (int i=0; i<m_rowIndex.size(); i++) {

		for (int j=0; j<m_rowIndex[i].size(); j++) {

			cout << m_rowIndex[i][j] << " "; 

		}

		cout << endl; 

	}

}

void SparseMatrix::sparsity() const {

	// compute number of non zeros 
	int nnz = 0; 
	for (int i=0; i<m_n; i++) {

		nnz += m_data[i].size(); 

	}

	cout << "number of unknowns = " << m_n << endl; 
	cout << "number of non zeros = " << nnz << endl; 
	cout << "sparsity factor = " << (double)nnz/(m_n*m_m) << endl; 

}