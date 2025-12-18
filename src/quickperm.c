#include "../include/quickperm.h"
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static inline void
pswap(int** a, int** b)
{
  int* t = *a;
  *a = *b;
  *b = t;
}

static int
factorial_size_t(size_t n, size_t* out)
{
  size_t r = 1;
  for (size_t i = 2; i <= n; ++i) {
    if (r > SIZE_MAX / i)
      return 0; // overflow
    r *= i;
  }
  *out = r;
  return 1;
}

static void
iswap(int* a, int* b)
{
  int t = *a;
  *a = *b;
  *b = t;
}

// QuickPerm: speichert alle Permutationen in *out (Layout: [perm0][perm1]...)
// Rückgabe: 0 ok, -1 overflow/zu groß, -2 malloc fehlgeschlagen
int
quickperm(int* input, size_t n, int** out, size_t* count)
{
  if (!out || !count)
    return -1;
  *out = NULL;
  *count = 0;

  if (n == 0)
    return 0;

  size_t perms = 0;
  if (!factorial_size_t(n, &perms))
    return -1;

  // Speicherbedarf prüfen: perms * n * sizeof(int)
  if (perms != 0 && n > SIZE_MAX / perms)
    return -1;
  size_t total_elems = perms * n;
  if (total_elems > SIZE_MAX / sizeof(int))
    return -1;

  int* buf = (int*)malloc(total_elems * sizeof(int));
  if (!buf)
    return -2;

  int* a = (int*)malloc(n * sizeof(int));
  if (!a) {
    free(buf);
    return -2;
  }
  memcpy(a, input, n * sizeof(int));

  // Controller-Array p
  size_t* p = (size_t*)malloc(n * sizeof(size_t));
  if (!p) {
    free(a);
    free(buf);
    return -2;
  }
  for (size_t i = 0; i < n; ++i)
    p[i] = i;

  // 1) erste Permutation speichern
  size_t w = 0;
  memcpy(buf + (w++ * n), a, n * sizeof(int));

  // 2) QuickPerm iterativ
  size_t i = 1;
  while (i < n) {
    if (p[i] == 0) {
      p[i] = i;
      ++i;
      continue;
    }

    p[i] -= 1;
    size_t j = (i & 1) ? p[i] : 0; // odd i -> j=p[i], even i -> j=0
    iswap(&a[j], &a[i]);

    memcpy(buf + (w++ * n), a, n * sizeof(int));

    i = 1;
  }

  free(p);
  free(a);

  *out = buf;
  *count = perms;
  return 0;
}

// QuickPerm: speichert alle Permutationen in *out (Layout: [perm0][perm1]...)
// Rückgabe: 0 ok, -1 overflow/zu groß, -2 malloc fehlgeschlagen
int
quickperm_ptr(int* const* input, size_t n, int*** out, size_t* count)
{
  if (!out || !count)
    return -1;
  *out = NULL;
  *count = 0;
  if (n == 0)
    return 0;

  size_t perms = 0;
  if (!factorial_size_t(n, &perms))
    return -1;

  if (perms != 0 && n > SIZE_MAX / perms)
    return -1;
  size_t total = perms * n;
  if (total > SIZE_MAX / sizeof(int*))
    return -1;

  int** buf = (int**)malloc(total * sizeof(int*));
  if (!buf)
    return -2;

  int** a = (int**)malloc(n * sizeof(int*));
  if (!a) {
    free(buf);
    return -2;
  }
  memcpy(a, input, n * sizeof(int*));

  size_t* p = (size_t*)malloc(n * sizeof(size_t));
  if (!p) {
    free(a);
    free(buf);
    return -2;
  }
  for (size_t i = 0; i < n; ++i)
    p[i] = i;

  size_t w = 0;
  memcpy(buf + (w++ * n), a, n * sizeof(int*));

  size_t i = 1;
  while (i < n) {
    if (p[i] == 0) {
      p[i] = i;
      ++i;
      continue;
    }
    p[i] -= 1;
    size_t j = (i & 1) ? p[i] : 0;
    pswap(&a[j], &a[i]);

    memcpy(buf + (w++ * n), a, n * sizeof(int*));
    i = 1;
  }

  free(p);
  free(a);

  *out = buf;
  *count = perms;
  return 0;
}
