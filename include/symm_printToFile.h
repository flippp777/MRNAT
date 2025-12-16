#ifndef SYMM_PRINTTOFILE
#define SYMM_PRINTTOFILE
#include <gmp.h>
#include <stdio.h>

void
symm_printToFile_2p2h(FILE* fptr,
                       int smax,
                       int **phi,
                       mpz_t phib[],
                       int kmax,
                       int Mref);

void
symm_printToFile_3p3h(FILE* fptr,
                      int smax,
                      int** phi,
                      mpz_t phib[],
                      int kmax,
                      int Mref);

#endif // !SYMM_PRINTTOFILE
