#ifndef GET_NPNH
#define GET_NPNH
#include <gmp.h>

void
get_npnh(int particles[],
         int holes[],
         int two_n,
         int smax,
         mpz_t phi_ph,
         mpz_t phi_nub);

#endif // !GET_NPNH
