#include "../include/hash.h"
#include <gmp.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

static inline uint64_t
splitmix64(uint64_t x)
{
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

uint64_t
hash_func(const mp_bitcnt_t* a, size_t block, size_t n)
{
  uint64_t h = 0x6a09e667f3bcc909ULL ^ (uint64_t)n; // Seed + Länge
  for (size_t i = 0; i < n; i++) {
    // Position einmischen, damit ähnliche Wertefolgen stärker auseinanderlaufen
    uint64_t v = (uint64_t)(uint32_t)a[block * n + i] ^
                 (0x9e3779b97f4a7c15ULL * (uint64_t)(i + 1));
    h ^= splitmix64(v);
    h *= 0x9ddfea08eb382d69ULL; // schnell, gute Durchmischung
  }
  return h;
}

uint64_t
hash_func_ptr(mp_bitcnt_t** a, size_t block, size_t n)
{
  uint64_t h = 0x6a09e667f3bcc909ULL ^ (uint64_t)n; // Seed + Länge
  for (size_t i = 0; i < n; i++) {
    // Position einmischen, damit ähnliche Wertefolgen stärker auseinanderlaufen
    uint64_t v = (uint64_t)(uint32_t)*a[block * n + i] ^
                 (0x9e3779b97f4a7c15ULL * (uint64_t)(i + 1));
    h ^= splitmix64(v);
    h *= 0x9ddfea08eb382d69ULL; // schnell, gute Durchmischung
  }
  return h;
}

void
hashset_init(HashSet* s, uint64_t* keys, uint8_t* used, size_t size_pow2)
{
  s->keys = keys;
  s->used = used;
  s->size = size_pow2;
  s->mask = size_pow2 - 1;
  memset(s->used, 0, s->size * sizeof(s->used[0]));
}

int
hashset_check(HashSet* s, uint64_t h)
{
  size_t idx = (size_t)h & s->mask;

  while (s->used[idx]) {
    if (s->keys[idx] == h)
      return 1; // already used
    idx = (idx + 1) & s->mask;
  }

  s->used[idx] = 1;
  s->keys[idx] = h;
  return 0; // new
}
