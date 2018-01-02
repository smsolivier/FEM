#include "Dense.H"
#include "CH_Timer.H"

Dense::Dense(const SparseMatrix& a_A) {

	// get dense version of matrix from SparseMatrix 
	m_A = a_A.getDenseMatrix(); 

	m_m = a_A.getM(); 
	m_n = a_A.getN(); 

}

void Dense::solve(vector<double>& a_rhs) {

    CH_TIMERS("Dense solve"); 

	// make a copy of the rhs (will return solution in a_rhs) 
	vector<double> b = a_rhs; 

	int err = gauss_elim(m_m, m_A, a_rhs, b); 

}

int Dense::gauss_elim(int dim, vector<vector<double> > &a, vector<double> &x, vector <double> &b) 
{
	vector<int> mGaussLastRow; 
	vector<int> mGaussLastCol; 
    int ierr = 0;
    if(mGaussLastCol.size() != dim) {
       mGaussLastCol.clear();
       mGaussLastRow.clear();
       mGaussLastCol.resize(dim, dim-1);
       mGaussLastRow.resize(dim, dim-1);
       // if(mGaussPatternOn) {
       //    setGaussPattern(dim, a);
       // }
    }

    // We first do dim reduction steps.
    for(int k = 0; k < dim; k++) {
        double aSubK = a[k][k];
        double bSubK = b[k];
        if(aSubK == 0.0) {              /* check */
            //fprintf(stderr, "gauss: a[%d][%d] = 0\n", k, k);
            //cout << "gauss: a[" << k << "][" << k << "] = 0\n";
            //exit(1);
            return -1;
        }

        for(int i = k+1; i <= mGaussLastCol[k]; i++) {
            double xtemp = a[i][k];
            if(xtemp == 0.0) continue;
            a[i][k] = 0.0;                      /* eq (4.1a) */
            xtemp /= aSubK;
            for(int j = k+1; j <= mGaussLastRow[k]; j++) {
                a[i][j] -= a[k][j]*xtemp;
            }
            b[i] -= bSubK*xtemp;              /* eq (4.1c) */
        }
    }
 
    // Now we perform dim back substitutions.
    for(int i = dim-1; i >= 0; i--) {
        b[i] /= a[i][i];
        for(int k = i - 1; k >= 0; k -= 1) {
            b[k] -= a[k][i] * b[i];
        }
    }

    // Return solution in x and b vectors
    for(int i = 0; i < dim; i++) x[i] = b[i];
    return ierr;
}