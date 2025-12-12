#include "../include/count_bits_right.h"
#include <gmp.h>

mp_bitcnt_t count_bits_right(mpz_t x, mp_bitcnt_t p) {
    mpz_t tmp;
    mpz_init(tmp);

    mpz_tdiv_r_2exp(tmp, x, p);   // tmp = x mod 2^p
    mp_bitcnt_t count = mpz_popcount(tmp);

    mpz_clear(tmp);
    return count;
}

