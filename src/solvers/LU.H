#ifndef _LU_H_
#define _LU_H_

#include "LinearSolver.H"
#ifdef SUPERLU
#include "slu_ddefs.h"
#endif
#ifdef PSUPERLU
#include "slu_mt_ddefs.h"
#endif
#include <vector> 
#include "SparseMatrix.H"

using namespace std; 

/// uses SuperLU to solve a linear system. Converts a SparseMatrix to a SuperMatrix 
/** OpenMP parallel SuperLU is used if linked correctly */ 
class LU : public LinearSolver {
public:
	LU(); 
	/// constructor 
	LU(const SparseMatrix& a_A); 
	/// set the matrix to solve. Used for multiple solves where the LHS changes 
	void operator()(const SparseMatrix& a_A); 

	/// solve interface 
	/** 
		performs a factorization and solve if the first solve. LU factorization is stored for subsequent solves. 
	*/ 
	void solve(vector<double>& a_rhs); 
private:
	// store size of matrix 
	int m_m; 
	int m_n; 

	#if defined SUPERLU || defined PSUPERLU
	/// superlu formatted matrix 
	SuperMatrix m_super; 
	/// store L, U factorization for repeated solves 
	SuperMatrix m_L, m_U; 
	#endif
	#ifdef SUPERLU 
	superlu_options_t m_options; 
	SuperLUStat_t m_stat; 
	#endif
	#ifdef PSUPERLU
	superlumt_options_t m_options; 
	#endif
	vector<int> m_perm_r; 
	vector<int> m_perm_c; 
	int m_info; 

	/// if true: do factorization and solve, false: just solve
	bool m_first;  
}; 

#endif 