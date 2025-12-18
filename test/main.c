#include "hash.h"
#include "quickperm.h"
#include <stdint.h>
#include <stdio.h>

int
main()
{
  size_t n = 4;
  int arr[] = { 1, 2, 3, 3 };

  size_t count;
  int* out;
  quickperm(arr, n, &out, &count);
  printf("count: %zu\n", count);
  printf("\n");

  size_t h_size = 32;
  int used[32];
  memset(used, 0, h_size * sizeof(int));
  uint64_t hashes[32];

  for (int i = 0; i < 24; i++) {
    printf("Permutation: ");
    for (int j = 0; j < 4; j++) {
      printf("%d", out[i * 4 + j]);
    }
    printf("\n");

    uint64_t hash = hash_func(out, i, n);
    printf("hash: %lu\n", hash);

    int already = 0;
    int idx = hash & (h_size - 1);

    while (used[idx]) {
      printf("used[idx] == %d\n", used[idx]);
      if (hashes[idx] == hash) {
        already = 1;
        printf("Already used\n");
        printf("\n");
      }
      idx = (idx + 1) & (h_size - 1);
    }
    if (already)
      continue;
    used[idx] = 1;
    hashes[idx] = hash;
    printf("Added\n");

    printf("\n");
  }
}
