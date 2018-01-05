#ifndef _LU_H_
#define _LU_H_

#include "LinearSolver.H"
#ifdef SUPERLU
#include "slu_ddefs.h"
#endif
#include <vector> 
#include "SparseMatrix.H"

using namespace std; 

/// uses SuperLU to solve a linear system. Converts a SparseMatrix to a SuperMatrix 
class LU : public LinearSolver {
public:
	LU(); 
	/// constructor 
	LU(const SparseMatrix& a_A); 
	/// set the matrix to solve
	void operator()(const SparseMatrix& a_A); 

	/// solve interface 
	/** 
		performs a factorization and solve if the first solve. LU factorization is stored for subsequent solves. 
	*/ 
	virtual void solve(vector<double>& a_rhs); 
private:
	// store size of matrix 
	int m_m; 
	int m_n; 

	#ifdef SUPERLU
	SuperMatrix m_super; /// SuperLU formatted matrix 
	SuperMatrix m_L, m_U; /// store L, U factorization 
	int * m_perm_r; 
	int * m_perm_c; 
	superlu_options_t m_options; 
	SuperLUStat_t m_stat; 
	int m_info; 
	#endif

	bool m_first; /// if true: do factorization and solve, false: just solve 
}; 

#endif 