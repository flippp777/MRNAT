#include "../include/hash_primitive.h"
#include <gmp.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>


uint64_t
hash_func_prim(const mp_bitcnt_t* a, size_t block, size_t n)
{
  uint64_t h = 0;
  for (size_t i = 0; i < n; i++) {
    // Position einmischen, damit ähnliche Wertefolgen stärker auseinanderlaufen
    h += a[block * n + i] * pow(10, n-i+1);
  }
  return h;
}

uint64_t
hash_func_ptr_prim(mp_bitcnt_t** a, size_t block, size_t n)
{
  uint64_t h = 0;
  for (size_t i = 0; i < n; i++) {
    // Position einmischen, damit ähnliche Wertefolgen stärker auseinanderlaufen
    h += *a[block * n + i] * pow(10, n-i+1);
  }
  return h;
}

// void
// hashset_init(HashSet* s, uint64_t* keys, uint8_t* used, size_t size_pow2)
// {
//   s->keys = keys;
//   s->used = used;
//   s->size = size_pow2;
//   s->mask = size_pow2 - 1;
//   memset(s->used, 0, s->size * sizeof(s->used[0]));
// }
//
// int
// hashset_check(HashSet* s, uint64_t h)
// {
//   size_t idx = (size_t)h & s->mask;
//
//   while (s->used[idx]) {
//     if (s->keys[idx] == h)
//       return 1; // already used
//     idx = (idx + 1) & s->mask;
//   }
//
//   s->used[idx] = 1;
//   s->keys[idx] = h;
//   return 0; // new
// }
