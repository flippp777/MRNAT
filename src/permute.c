#include "../include/permute.h"
#include <stdlib.h>
void
permute_rek(int** perms, int* arr, int start, int n, int perm_count, int* count)
{
  if (start == n) {
    for (int i = 0; i < n; i++)
      perms[*count][i] = arr[i];
    (*count)++;
    return;
  }

  for (int i = start; i < n; i++) {
    int tmp = arr[start];
    arr[start] = arr[i];
    arr[i] = tmp;

    permute_rek(perms, arr, start + 1, n, perm_count, count);

    tmp = arr[start];
    arr[start] = arr[i];
    arr[i] = tmp;
  }
}

int**
permute(int* arr, int start, int n, int perm_count)
{
  int** perms = malloc(perm_count * sizeof(int*));
  for (int i = 0; i < perm_count; i++) {
    perms[i] = calloc(n, sizeof(int));
  }
  int count = 0;

  if (start == n) {
    for (int i = 0; i < n; i++)
      perms[count][i] = arr[i];
    count++;
    return perms;
  }

  for (int i = start; i < n; i++) {
    int tmp = arr[start];
    arr[start] = arr[i];
    arr[i] = tmp;

    permute_rek(perms, arr, start + 1, n, perm_count, &count);

    tmp = arr[start];
    arr[start] = arr[i];
    arr[i] = tmp;
  }

  return perms;
}

void
permute_ptr_rek(int*** perms, int** arr, int start, int n, int perm_count, int* count)
{
  if (start == n) {
    for (int i = 0; i < n; i++)
      perms[*count][i] = arr[i];
    (*count)++;
    return;
  }

  for (int i = start; i < n; i++) {
    int* tmp = arr[start];
    arr[start] = arr[i];
    arr[i] = tmp;

    permute_ptr_rek(perms, arr, start + 1, n, perm_count, count);

    tmp = arr[start];
    arr[start] = arr[i];
    arr[i] = tmp;
  }
}

int***
permute_ptr(int** arr, int start, int n, int perm_count)
{
  int*** perms = malloc(perm_count * sizeof(int**));
  for (int i = 0; i < perm_count; i++) {
    perms[i] = calloc(n, sizeof(int*));
  }
  int count = 0;

  if (start == n) {
    for (int i = 0; i < n; i++)
      perms[count][i] = arr[i];
    count++;
    return perms;
  }

  for (int i = start; i < n; i++) {
    int* tmp = arr[start];
    arr[start] = arr[i];
    arr[i] = tmp;

    permute_ptr_rek(perms, arr, start + 1, n, perm_count, &count);

    tmp = arr[start];
    arr[start] = arr[i];
    arr[i] = tmp;
  }

  return perms;
}
