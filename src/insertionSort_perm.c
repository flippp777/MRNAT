#include "../include/insertionSort_perm.h"

void
insertionSort_perm(int arr[], int n, int* perm)
{
  for (int i = 1; i < n; ++i) {
    int key = arr[i];
    int j = i - 1;

    /* Move elements of arr[0..i-1], that are
       greater than key, to one position ahead
       of their current position */
    while (j >= 0 && arr[j] > key) {
      arr[j + 1] = arr[j];
      j = j - 1;
      *perm += 1;
    }
    arr[j + 1] = key;
  }
}

void
insertionSort_perm_steps(int arr[], int n, int* perm, int* steps)
{
  for (int i = 1; i < n; ++i) {
    *steps += 1;
    int key = arr[i];
    int j = i - 1;

    /* Move elements of arr[0..i-1], that are
       greater than key, to one position ahead
       of their current position */
    while (j >= 0 && arr[j] > key) {
      *steps += 1;
      arr[j + 1] = arr[j];
      j = j - 1;
      *perm += 1;
    }
    arr[j + 1] = key;
  }
}
