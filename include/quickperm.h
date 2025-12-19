#ifndef QUICKPERM
#define QUICKPERM

#include <gmp.h>
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

int
quickperm(mp_bitcnt_t* input, size_t n, mp_bitcnt_t** out, size_t* count);

int
quickperm_ptr(mp_bitcnt_t* const* input, size_t n, mp_bitcnt_t*** out, size_t* count);

#endif // !QUICKPERM
