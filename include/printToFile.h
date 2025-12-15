#ifndef PRINTTOFILE
#define PRINTTOFILE
#include <gmp.h>
#include <stdio.h>

void
printToFile(FILE* fptr,
            mpz_t phib[],
            int smax,
            int** phi,
            int kmax,
            int mu,
            int nu,
            int b,
            int a,
            int s,
            int r,
            mpz_t phi_nub_exc,
            int sign);

#endif // !PRINTTOFILE
