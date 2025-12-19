#ifndef HASH
#define HASH
#include <gmp.h>
#include <stddef.h>
#include <stdint.h>

uint64_t
hash_func(const mp_bitcnt_t* a, size_t block, size_t n);

uint64_t
hash_func_ptr(mp_bitcnt_t** a, size_t block, size_t n);

typedef struct
{
  uint64_t* keys; // saved hashes
  uint8_t* used;  // 0/1 found
  size_t size;    // power of two (faster)
  size_t mask;    // size - 1
} HashSet;

void
hashset_init(HashSet* s, uint64_t* keys, uint8_t* used, size_t size_pow2);

int
hashset_check(HashSet* s, uint64_t h);

#endif // !HASH
