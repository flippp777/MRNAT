#ifndef BRUTEFORCE_PRINTTOFILE
#define BRUTEFORCE_PRINTTOFILE
#include <stdio.h>

void bruteForce_printToFile_2p2h(FILE *fptr, int states[],
                                 unsigned long long statesb[], int smax,
                                 int **phi, int phib[], int kmax, int Mref);

void bruteForce_printToFile_3p3h(FILE *fptr, int states[],
                                 unsigned long long statesb[], int smax,
                                 int **phi, int phib[], int kmax, int Mref);

#endif // !BRUTEFORCE_PRINTTOFILE
