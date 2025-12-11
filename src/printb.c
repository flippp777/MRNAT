#include "../include/printb.h"
#include <stdio.h>

void
printb(int n, int smax)
{
  if (n == 0) {
    for (int i = 0; i < smax; i++)
      printf("0");
    return;
  }

  int binary[smax];
  for (int i = 0; i < smax; i++)
    binary[i] = 0;
  int index = 0;

  while (n > 0) {
    binary[index] = n % 2;
    n = n / 2;
    index++;
  }

  for (int i = smax - 1; i >= 0; i--) {
    printf("%d", binary[i]);
  }
}

void
fprintb(FILE* fptr, int n, int smax)
{
  if (n == 0) {
    for (int i = 0; i < smax; i++)
      fprintf(fptr, "0");
    return;
  }

  int binary[smax];
  for (int i = 0; i < smax; i++)
    binary[i] = 0;
  int index = 0;

  while (n > 0) {
    binary[index] = n % 2;
    n = n / 2;
    index++;
  }

  for (int i = smax - 1; i >= 0; i--) {
    fprintf(fptr, "%d", binary[i]);
  }
}
