#ifndef BRUTEFORCE_EVAL
#define BRUTEFORCE_EVAL

#include <gmp.h>
void
bruteForce_eval_2p2h(int states[],
                     unsigned long long statesb[],
                     int smax,
                     int** phi,
                     int phib[],
                     int kmax,
                     int Mref);

void
bruteForce_eval_3p3h(int states[],
                     unsigned long long statesb[],
                     int smax,
                     int** phi,
                     int phib[],
                     int kmax,
                     int Mref);
void
bruteForce_eval_4p4h(int smax, mpz_t phib[], int kmax, int Mref);

#endif // !BRUTEFORCE_EVAL
