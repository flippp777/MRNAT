#include "../include/permute.h"

void permute(int *arr, int start, int n, int **perms, int *count) {
  if (start == n) {
    // Permutation speichern
    for (int i = 0; i < n; i++)
      perms[*count][i] = arr[i];
    (*count)++;
    return;
  }

  for (int i = start; i < n; i++) {
    int tmp = arr[start];
    arr[start] = arr[i];
    arr[i] = tmp;

    permute(arr, start + 1, n, perms, count);

    tmp = arr[start];
    arr[start] = arr[i];
    arr[i] = tmp;
  }
}

int permute_sign(int *arr, int n) {
  int inv = 0;
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      if (arr[i] > arr[j])
        inv++;
    }
  }
  return (inv % 2 == 0) ? 1 : -1;
}
