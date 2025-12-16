#ifndef PRINTTOFILE
#define PRINTTOFILE
#include <gmp.h>
#include <stdio.h>

void
printToFile_2p2h(FILE* fptr,
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
            int perm,
            int sign);

void
printToFile_3p3h(FILE* fptr,
                 mpz_t phib[],
                 int smax,
                 int** phi,
                 int kmax,
                 int mu,
                 int nu,
                 int b,
                 int a,
                 int d,
                 int c,
                 int s,
                 int r,
                 mpz_t phi_nub_exc,
                 int perm,
                 int sign);

#endif // !PRINTTOFILE
