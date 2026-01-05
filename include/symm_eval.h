#ifndef SYMM_EVAL
#define SYMM_EVAL

#include <gmp.h>
void
symm_eval_2p2h(int states[],
               unsigned long long statesb[],
               int smax,
               int** phi,
               int phib[],
               int kmax,
               int Mref);

void
symm_eval_3p3h(int states[],
               unsigned long long statesb[],
               int smax,
               int** phi,
               int phib[],
               int kmax,
               int Mref);

void
symm_eval_4p4h(mp_bitcnt_t smax, int** phi, mpz_t phib[], size_t kmax, size_t Mref);

void
symm_eval_5p5h(mp_bitcnt_t smax, int** phi, mpz_t phib[], size_t kmax, size_t Mref);

#endif // !SYMM_EVAL
