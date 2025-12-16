#ifndef BRUTEFORCE_PRINTTOFILE
#define BRUTEFORCE_PRINTTOFILE
#include <gmp.h>
#include <stdio.h>

void
bruteForce_printToFile_2p2h(FILE* fptr,
                            int smax,
                            int** phi,
                            mpz_t phib[],
                            int kmax,
                            int Mref);

void
bruteForce_printToFile_3p3h(FILE* fptr,
                            int smax,
                            int** phi,
                            mpz_t phib[],
                            int kmax,
                            int Mref);

#endif // !BRUTEFORCE_PRINTTOFILE
