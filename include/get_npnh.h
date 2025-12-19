#ifndef GET_NPNH
#define GET_NPNH
#include <gmp.h>

void
get_npnh(mp_bitcnt_t particles[],
         mp_bitcnt_t holes[],
         size_t two_n,
         mpz_t phi_ph,
         mpz_t phi_nub);

#endif // !GET_NPNH
