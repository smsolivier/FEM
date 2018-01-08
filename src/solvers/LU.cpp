#include "LU.H"
#include "CH_Timer.H"

LU::LU() {

#ifdef SUPERLU 
	// setup solver options 
	set_default_options(&m_options); 
	m_options.ColPerm = NATURAL; 
	StatInit(&m_stat); 
#elif defined PSUPERLU
	cout << "define default constructor for parallel LU" << endl; 
#else 
	cout << "ERROR: SuperLU not defined for LU solver" << endl; 
	exit(0); 
#endif
}

LU::LU(const SparseMatrix& a_A) {

#ifdef SUPERLU 
	m_A = a_A; 

	m_m = a_A.getM(); 
	m_n = a_A.getN(); 

	// get the super matrix form from a_A 
	m_super = a_A.getSuperMatrix(); 

	// setup solver options 
	set_default_options(&m_options); 
	m_options.ColPerm = NATURAL; 
	if (a_A.symmetric()) {
		m_options.SymmetricMode = YES; 
	}
	StatInit(&m_stat); 

	// initialize permutation arrays 
	m_perm_r.resize(m_m);  
	m_perm_c.resize(m_m);  

	// set to true so factorization will occur 
	m_first = true; 
#elif defined PSUPERLU 
	m_A = a_A; 
	m_n = m_A.getN(); 
	m_m = m_A.getM(); 

	m_super = m_A.getSuperMatrix();

	m_perm_c.resize(m_m); 
	m_perm_r.resize(m_m); 
#else 
	cout << "ERROR: SuperLU not defined for LU solver" << endl; 
	exit(0); 
#endif
}

void LU::operator()(const SparseMatrix& a_A) {

#ifdef SUPERLU
	m_A = a_A; 
	m_m = m_A.getM(); 
	m_n = m_A.getN(); 

	m_super = m_A.getSuperMatrix();

	// initialize permutation arrays 
	m_perm_r.resize(m_m);  
	m_perm_c.resize(m_m);  

	// set to true so factorization will occur 
	m_first = true; 
#endif
}

void LU::solve(vector<double>& a_rhs) {
#ifdef SUPERLU

	// setup rhs 
	SuperMatrix b; 
	int nrhs = 1; 
	dCreate_Dense_Matrix(&b, m_m, nrhs, &a_rhs[0], m_m, SLU_DN, SLU_D, SLU_GE); 

	// factor and solve if first solve 
	if (m_first) {

		CH_TIMERS("LU factorize and solve"); 

		dgssv(&m_options, &m_super, &m_perm_c[0], &m_perm_r[0], &m_L, &m_U, 
			&b, &m_stat, &m_info); 

		m_first = false; // factorization is done now (stored in m_L and m_U)

	}

	// just solve if factorization is already done 
	else {

		CH_TIMERS("LU solve"); 

		dgstrs(NOTRANS, &m_L, &m_U, &m_perm_c[0], &m_perm_r[0], &b, &m_stat, &m_info); 

	}
#elif defined PSUPERLU
	CH_TIMERS("PLU solve"); 

	SuperMatrix b; 
	int nrhs = 1; 
	dCreate_Dense_Matrix(&b, m_m, nrhs, &a_rhs[0], m_m, SLU_DN, SLU_D, SLU_GE); 

	int_t nproc = 4; 

	get_perm_c(0, &m_super, &m_perm_c[0]); 
	pdgssv(nproc, &m_super, &m_perm_c[0], &m_perm_r[0], &m_L, &m_U, 
		&b, &m_info); 
#endif
}