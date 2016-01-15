#include "eigenvalues.h"
#include "general.h"

void eigs(double * eigenvalues, double * eigenvectors, double n_eigs, bool getEigenvectors, double * a, double * b, int siz){
  char jobz, range, uplo;
  int itype, n, lda, ldb, il, iu, m, ldz, lwork, *iwork, *ifail, info;
  double *vl, *vu, abstol,/* *w = eigenvalues, *z = eigenvectors,*/ *work;

  if (getEigenvectors){ 
    jobz = 'V';
    ldz = siz;
  }
  else{
    jobz = 'N';
    ldz = 1;
  }
  range = 'I';
  uplo = 'L';
  
  n = siz;
  itype = 1;
  lda = siz;
  ldb = siz;
  il = 1;
  iu = (n_eigs > n) ? n : n_eigs;

  iwork = new int[8*n];
  ifail = new int[5*n];
  
  lwork = 15 * n;
  work = new double[lwork];
  vu = new double[1];
  vl = new double[1];

  abstol = 1E-8;

  dsygvx_(&itype, &jobz, &range, &uplo, &n, a, &lda, b, &ldb, vl, vu, &il, &iu, &abstol, &m, eigenvalues, eigenvectors, &ldz, work, &lwork, iwork, ifail, &info);

  delete [] iwork;
  delete [] ifail;
  delete [] work;
  delete [] vu;
  delete [] vl;
  if (info != 0){
    cout << "Diag problem? info = " << info << endl;
  }
}
