#ifndef GET_NPNH
#define GET_NPNH
#include <gmp.h>

void
get_npnh(mpz_t particles[],
         mpz_t holes[],
         int two_n,
         int smax,
         mpz_t phi_ph,
         mpz_t phi_nub);

#endif // !GET_NPNH
