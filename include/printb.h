#ifndef PRINTB
#define PRINTB
#include <gmp.h>
#include <stdio.h>

void printb(mpz_t num, int smax);

void fprintb(FILE *fptr, mpz_t num, int smax);

#endif // !PRINTB
