#ifndef SYMM_PRINTTOFILE
#define SYMM_PRINTTOFILE
#include <gmp.h>
#include <stdio.h>

void
symm_printToFile_2p2h(FILE* fptr,
                       mp_bitcnt_t smax,
                       int **phi,
                       mpz_t phib[],
                       size_t kmax,
                       size_t Mref);

void
symm_printToFile_3p3h(FILE* fptr,
                      mp_bitcnt_t smax,
                      int** phi,
                      mpz_t phib[],
                      size_t kmax,
                      size_t Mref);

#endif // !SYMM_PRINTTOFILE
