#include "../include/get_npnh.h"
#include <gmp.h>
void
get_npnh(mpz_t particles[],
         mpz_t holes[],
         int two_n,
         int smax,
         mpz_t phi_ph,
         mpz_t phi_nub)
{
  int state_idx = 0;
  int hole_idx = 0;
  int particle_idx = 0;

  for (int i; i < two_n; i++) {
    state_idx = mpz_scan1(phi_ph, state_idx);

    if (mpz_tstbit(phi_nub, state_idx) == 1) {
      mpz_setbit(holes[hole_idx], state_idx);
      hole_idx++;
    } else {
      mpz_setbit(particles[particle_idx], state_idx);
      particle_idx++;
    }
  }
}
