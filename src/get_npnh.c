#include "../include/get_npnh.h"
#include <gmp.h>
void
get_npnh(mp_bitcnt_t particles[],
         mp_bitcnt_t holes[],
         size_t two_n,
         mpz_t phi_ph,
         mpz_t phi_nub)
{
  mp_bitcnt_t state_idx = 0;
  size_t hole_idx = 0;
  size_t particle_idx = 0;

  for (size_t i = 0; i < two_n; i++) {
    state_idx = mpz_scan1(phi_ph, state_idx);

    if (mpz_tstbit(phi_nub, state_idx) == 1) {
      holes[hole_idx] = state_idx;
      hole_idx++;
    } else {
      particles[particle_idx] = state_idx;
      particle_idx++;
    }
    state_idx++;
  }
}
