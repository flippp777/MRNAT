#include "hash.h"
#include <stddef.h>
#include <stdint.h>

static inline uint64_t
splitmix64(uint64_t x)
{
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

uint64_t
hash_func(const int* a, size_t block, size_t n)
{
  uint64_t h = 0x6a09e667f3bcc909ULL ^ (uint64_t)n; // Seed + Länge
  for (size_t i = 0; i < n; i++) {
    // Position einmischen, damit ähnliche Wertefolgen stärker auseinanderlaufen
    uint64_t v =
      (uint64_t)(uint32_t)a[block * n + i] ^ (0x9e3779b97f4a7c15ULL * (uint64_t)(i + 1));
    h ^= splitmix64(v);
    h *= 0x9ddfea08eb382d69ULL; // schnell, gute Durchmischung
  }
  return h;
}

uint64_t
hash_func_ptr(int** a, size_t block, size_t n)
{
  uint64_t h = 0x6a09e667f3bcc909ULL ^ (uint64_t)n; // Seed + Länge
  for (size_t i = 0; i < n; i++) {
    // Position einmischen, damit ähnliche Wertefolgen stärker auseinanderlaufen
    uint64_t v =
      (uint64_t)(uint32_t)*a[block * n + i] ^ (0x9e3779b97f4a7c15ULL * (uint64_t)(i + 1));
    h ^= splitmix64(v);
    h *= 0x9ddfea08eb382d69ULL; // schnell, gute Durchmischung
  }
  return h;
}
