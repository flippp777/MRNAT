#ifndef SYMM_PRINTTOFILE
#define SYMM_PRINTTOFILE
#include <stdio.h>

void
symm_printToFile_2p2h(FILE* fptr,
                       int states[],
                       unsigned long long statesb[],
                       int smax,
                       int **phi,
                       int phib[],
                       int kmax,
                       int Mref);

#endif // !SYMM_PRINTTOFILE
