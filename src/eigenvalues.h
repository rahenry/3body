#ifndef EIGENVALUES_H
#define EIGENVALUES_H

extern "C" void dsygvx_(int *itype, char *jobz, char *range, char *uplo, int *n,
                           double *a, int *lda, double *b, int *ldb,
                                      double *vl, double *vu, int *il, int *iu, double *abstol, int *m,
                                                             double *w, double *z, int *ldz, double *work,
                                                                        int *lwork, int *iwork, int *ifail, int *info);

void eigs(double * eigenvalues, double * eigenvectors, double n_eigs, bool getEigenvectors, double * a, double * b, int siz);

#endif

